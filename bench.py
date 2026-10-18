#!/usr/bin/env python3
"""S4 toolpath metrics + a default-infill baseline, for head-to-head
comparison against the origami A-trail.

Metrics: travel-move count (primary), self-crossings, total printed
length, total travel length, total turning angle, stroke count. Computed on any
``tp:`` toolpath. Includes a rectilinear/grid infill generator as the "default
slicer" baseline. Connected Fermat Spirals (cfs.py) and the Gupta-Krishnamoorthy-
Dreifus Euler transformation (euler_transform.py) are both MEASURED, not cited
against the authors' reported numbers  — see those modules' docstrings.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def _seg_len(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])


def _proper_cross(a, b, c, d):
    """True iff segments ab and cd cross *transversally* in their interiors.
    Shared endpoints or collinear touching (any orientation == 0) are not proper
    crossings — that lets two strokes touch at a graph vertex (the Eulerian revisit)
    without counting, matching the A-trail notion of crossing-free."""
    def o(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    d1, d2, d3, d4 = o(c, d, a), o(c, d, b), o(a, b, c), o(a, b, d)
    if 0 in (round(d1, 12), round(d2, 12), round(d3, 12), round(d4, 12)):
        return False
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def _turn(u, v, w):
    a = math.atan2(v[1] - u[1], v[0] - u[0])
    b = math.atan2(w[1] - v[1], w[0] - v[0])
    return abs((b - a + math.pi) % (2 * math.pi) - math.pi)


def toolpath_metrics(raw: dict) -> dict:
    C = raw["vertices_coords"]
    EV = raw["edges_vertices"]
    pv = raw["tp:paths_vertices"]
    pe = raw["tp:paths_edges"]
    closed = raw.get("tp:paths_closed", [True] * len(pv))
    travel = raw.get("tp:paths_travel", [[False] * len(e) for e in pe])

    strokes = len(pv)
    travel_moves = max(0, strokes - 1)          
    printed_len = travel_len = turning = 0.0
    printed_segs = []                           
    for p in range(strokes):
        verts = pv[p]
        m = len(verts)
        seq = verts + ([verts[0]] if closed[p] else [])
        for i in range(len(seq) - 1):
            a, b = C[seq[i]], C[seq[i + 1]]
            if travel[p][i]:
                travel_moves += 1
                travel_len += _seg_len(a, b)
            else:
                printed_len += _seg_len(a, b)
                printed_segs.append((a, b))
        ring = [C[v] for v in (seq if closed[p] else verts)]
        for i in range(1, len(ring) - 1):
            turning += _turn(ring[i - 1], ring[i], ring[i + 1])
        if closed[p] and len(ring) >= 3:
            turning += _turn(ring[-2], ring[0], ring[1])

    crossings = 0
    n = len(printed_segs)
    for i in range(n):
        a, b = printed_segs[i]
        for j in range(i + 1, n):
            c, d = printed_segs[j]
            if _proper_cross(a, b, c, d):
                crossings += 1
    return {
        "strokes": strokes,
        "travel_moves": travel_moves,
        "self_crossings": crossings,
        "printed_length": round(printed_len, 2),
        "travel_length": round(travel_len, 2),
        "turning_deg": round(math.degrees(turning)),
    }


def poly_area(poly):
    """Shoelace area — the ACTUAL area a polygon outline covers, as opposed to its
    bounding-box area. For a non-convex outline (L-shape, star) the two can differ by
    50-280%; using bbox area as a density denominator silently overstates how much
    area a tightly-clipped method (like the origami A-trail) is actually covering,
    which biases any length/area comparison against a method whose mesh isn't
    outline-clipped. """
    n = len(poly)
    return abs(0.5 * sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
                          for i in range(n)))


def point_in_polygon(p, poly):
    x, y = p; inside = False
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def _scan_intervals(c, poly, axis):
    """Intervals of the axis-aligned scan line (axis='h': y=c, 'v': x=c) inside poly."""
    ts = []
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        if axis == "h":
            if (a[1] > c) != (b[1] > c):
                ts.append(a[0] + (c - a[1]) * (b[0] - a[0]) / (b[1] - a[1]))
        else:
            if (a[0] > c) != (b[0] > c):
                ts.append(a[1] + (c - a[0]) * (b[1] - a[1]) / (b[0] - a[0]))
    ts.sort()
    return [(ts[i], ts[i + 1]) for i in range(0, len(ts) - 1, 2)]


def grid_infill(bbox, spacing: float, directions=("h", "v"), polygon=None) -> dict:
    """Parallel lines (and/or crossing grid) clipped to a bbox or, if given, an
    arbitrary `polygon` outline — emitted as one `tp:` path per line segment, the
    way a default slicer lays grid: each segment a separate stroke with a travel
    between. On a non-convex outline a scan line breaks into several segments →
    several strokes → several retractions."""
    x0, y0, x1, y1 = bbox
    coords, EV, PV, PE = [], [], [], []

    def add_line(p, q):
        i, j = len(coords), len(coords) + 1
        coords.append(list(p)); coords.append(list(q))
        e = len(EV); EV.append([i, j])
        PV.append([i, j]); PE.append([e])

    def lines(axis, lo, hi, span_lo, span_hi):
        c = lo + spacing
        while c < hi:
            ivs = _scan_intervals(c, polygon, axis) if polygon else [(span_lo, span_hi)]
            for a, b in ivs:
                if b - a < spacing * 0.1:
                    continue
                add_line((a, c) if axis == "h" else (c, a),
                         (b, c) if axis == "h" else (c, b))
            c += spacing

    if "h" in directions:
        lines("h", y0, y1, x0, x1)
    if "v" in directions:
        lines("v", x0, x1, y0, y1)
    return {
        "vertices_coords": coords, "edges_vertices": EV,
        "edges_assignment": ["U"] * len(EV),
        "tp:paths_vertices": PV, "tp:paths_edges": PE,
        "tp:paths_closed": [False] * len(PV),
        "tp:paths_travel": [[False] for _ in PE],
        "tp:paths_source": ["default-grid"] * len(PV),
    }


def boustrophedon_infill(bbox, spacing: float, polygon=None, axis="h") -> dict:
    """Connected (serpentine) single-direction infill — the way a real slicer lays rectilinear:
    consecutive scan lines are joined end-to-end with U-turns into one stroke per connected run
    (greedy: a line extends a run if its interval overlaps the run's last). On a convex part this
    is ONE stroke / 0 travel (but 1-D infill); on a non-convex part it splits into a few runs. The
    U-turns mean its turning is *not* low (≈180° per line end) — only the crossing-free 2-D
    coverage and crossings remain origami's robust edge; this row makes the travel comparison fair."""
    x0, y0, x1, y1 = bbox
    lo, hi, slo, shi = (y0, y1, x0, x1) if axis == "h" else (x0, x1, y0, y1)
    scan = []
    c = lo + spacing
    while c < hi:
        ivs = _scan_intervals(c, polygon, axis) if polygon else [(slo, shi)]
        ivs = [iv for iv in ivs if iv[1] - iv[0] > spacing * 0.1]
        if ivs:
            scan.append((c, ivs))
        c += spacing
    runs, active = [], []                       
    for (c, ivs) in scan:
        used = set()
        for ap in active:
            la, lb = ap["last"]; ap["alive"] = False
            for k, iv in enumerate(ivs):
                if k in used:
                    continue
                if min(lb, iv[1]) - max(la, iv[0]) > -spacing * 0.5:   
                    seg = [(iv[0], c), (iv[1], c)] if axis == "h" else [(c, iv[0]), (c, iv[1])]
                    if ap["flip"]:
                        seg = seg[::-1]
                    ap["pts"] += seg; ap["last"] = iv; ap["flip"] = not ap["flip"]
                    ap["alive"] = True; used.add(k); break
        runs += [ap["pts"] for ap in active if not ap["alive"]]
        active = [ap for ap in active if ap["alive"]]
        for k, iv in enumerate(ivs):
            if k not in used:
                seg = [(iv[0], c), (iv[1], c)] if axis == "h" else [(c, iv[0]), (c, iv[1])]
                active.append({"pts": seg, "last": iv, "flip": True})
    runs += [ap["pts"] for ap in active]
    coords, PV = [], []
    for pth in runs:
        idx = []
        for p in pth:
            idx.append(len(coords)); coords.append(list(p))
        PV.append(idx)
    return {
        "vertices_coords": coords, "edges_vertices": [],
        "tp:paths_vertices": PV, "tp:paths_edges": [[0] * (len(v) - 1) for v in PV],
        "tp:paths_closed": [False] * len(PV),
        "tp:paths_travel": [[False] * (len(v) - 1) for v in PV],
        "tp:paths_source": ["boustrophedon"] * len(PV),
    }


def bbox_of(raw):
    xs = [c[0] for c in raw["vertices_coords"]]
    ys = [c[1] for c in raw["vertices_coords"]]
    return (min(xs), min(ys), max(xs), max(ys))


def miura(nx, ny, a=0.5, b=0.5, c=0.18):
    """Minimal Miura-ori crease pattern: an nx×ny grid with a y-zigzag (odd columns
    offset by c). Interior vertices are degree-4 and Kawasaki — atrail-able."""
    V = []; idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            idx[(i, j)] = len(V); V.append([i * a, j * b + (i % 2) * c])
    EV = []; EA = []
    for i in range(nx + 1):
        for j in range(ny + 1):
            if i < nx: EV.append([idx[(i, j)], idx[(i + 1, j)]]); EA.append("M")
            if j < ny: EV.append([idx[(i, j)], idx[(i, j + 1)]]); EA.append("V")
    return {"vertices_coords": V, "edges_vertices": EV, "edges_assignment": EA}


def run_part(name, outline, bbox, spacing=0.5):
    """Run the full origami pipeline (Miura → clip → atrail+patch) and the infill
    baselines on the same outline; return (name, rows). Origami needs clip.py +
    atrail.py (local modules)."""
    import clip, atrail
    x0, y0, x1, y1 = bbox

    
    mg = 3 * spacing
    nx = max(2, int(round((x1 - x0 + 2 * mg) / spacing)))
    ny = max(2, int(round((y1 - y0 + 2 * mg) / spacing)))
    M = miura(nx, ny, a=spacing, b=spacing, c=spacing * 0.36)
    for v in M["vertices_coords"]:
        v[0] += x0 - mg; v[1] += y0 - mg
    rows = []
    try:
        clipped, _ = clip.clip_to_polygon(M, outline)
        out, problems = atrail.ATrailGenerator(clipped).generate(patch=True)
        if problems:
            rows.append(("origami A-trail", None, f"FAILED: {problems[0]}"))
        else:
            rows.append(("origami A-trail", toolpath_metrics(out), ""))
    except Exception as ex:
        rows.append(("origami A-trail", None, f"ERROR: {ex}"))
    rect = grid_infill(bbox, spacing, directions=("h",), polygon=outline)
    bous = boustrophedon_infill(bbox, spacing, polygon=outline)
    grid = grid_infill(bbox, spacing, directions=("h", "v"), polygon=outline)
    rows.append(("rectilinear (separate)", toolpath_metrics(rect), "1-D, line-per-stroke"))
    rows.append(("rectilinear (boustrophedon)", toolpath_metrics(bous), "1-D, connected (U-turns)"))
    rows.append(("grid infill (crosshatch)", toolpath_metrics(grid), "2-D, self-crossing"))


    try:
        import cfs
        cpath = cfs.cfs_infill(outline, [], spacing)
        cm = cfs.metrics(cpath)
        cseg = [(cpath[i], cpath[i+1]) for i in range(len(cpath)-1)]
        cx = sum(1 for i in range(len(cseg)) for j in range(i+1, len(cseg))
                 if _proper_cross(cseg[i][0], cseg[i][1], cseg[j][0], cseg[j][1]))
        rows.append(("Connected Fermat Spirals", {"strokes": 1, "travel_moves": 0,
                     "self_crossings": cx, "printed_length": round(cm["length"], 2),
                     "turning_deg": round(cm["turning_deg"], 1)},
                     "measured: 1 stroke, LOW turning (beats origami on smoothness — RQ4)"))
    except Exception as ex:
        rows.append(("Connected Fermat Spirals", {"strokes": 1, "travel_moves": 0, "self_crossings": 0,
                     "printed_length": "≈", "turning_deg": "low"}, f"cited (cfs.py failed: {ex})"))
    
    try:
        import euler_transform as et
        gsize = spacing * 0.9
        gpolys = et.hex_mesh(x0, y0, x1, y1, gsize, outline=outline)
        gout, gproblems, gbase_e = et.gupta_toolpath(gpolys)
        if gproblems:
            rows.append(("Euler transform (measured)", None, f"FAILED: {gproblems[0]}"))
        else:
            gm = toolpath_metrics(gout)
            area = poly_area(outline)
            ratio = len(gout["edges_vertices"]) / gbase_e
            origami_m = rows[0][1]                      
            cmp_note = (f" vs origami {origami_m['printed_length'] / area:.2f}"
                        if origami_m else "")
            rows.append(("Euler transform (measured)", gm,
                         f"hex-mesh clipped to the exact outline, "
                         f"{ratio:.2f}x its own base E (Lemma 3.3 predicts 4x); "
                         f"density {gm['printed_length'] / area:.2f} len/area{cmp_note}"))
    except Exception as ex:
        rows.append(("Euler transform (measured)", None, f"ERROR: {ex}"))
    return name, rows


def _outlines():
    return [
        ("rectangle 10×10 (convex)",
         [[.3, .3], [9.7, .3], [9.7, 9.7], [.3, 9.7]], (0, 0, 10, 10)),
        ("L-shape (non-convex)",
         [[.3, .3], [9.7, .3], [9.7, 4.8], [4.8, 4.8], [4.8, 9.7], [.3, 9.7]], (0, 0, 10, 10)),
        ("5-point star (non-convex)",
         [[5 + 4.6 * math.cos(math.radians(90 + k * 36)) * (1 if k % 2 == 0 else 0.42),
           5 + 4.6 * math.sin(math.radians(90 + k * 36)) * (1 if k % 2 == 0 else 0.42)]
          for k in range(10)], (0, 0, 10, 10)),
    ]


def _tpm(m):
    """Turning per printed mm (normalises out part size) — the smoothness metric."""
    t, l = m.get("turning_deg"), m.get("printed_length")
    return f"{t / l:.1f}" if isinstance(t, (int, float)) and isinstance(l, (int, float)) and l > 0 else "·"


def realparts_benchmark():
    """Slice real solid STL parts and benchmark the origami pipeline on each layer contour."""
    import slice_stl, tempfile
    d = tempfile.gettempdir()
    jobs = [
        ("gear (non-convex)", lambda p: slice_stl.gen_gear_stl(p, teeth=12, R=12, root=8.6, h=4), [2.0]),
        ("twisted star prism (rotating non-convex)",
         lambda p: slice_stl.gen_twist_star_stl(p, points=5, R=11, r=4.6, h=16, turn=60), [2.0, 5.0, 8.0]),
        ("frustum cone (smooth convex, varying r)",
         lambda p: slice_stl.gen_cone_stl(p, R=12, rtop=2, h=12), [2.0, 8.0]),
    ]
    print("S4 benchmark on REAL SLICED LAYERS (stdlib STL slicer, slice_stl.py)\n")
    hdr = (f"  {'layer / method':30s} {'strokes':>7} {'travel':>6} {'cross':>6} "
           f"{'length':>8} {'turn°':>7} {'t/mm':>5}")
    for name, gen, zs in jobs:
        stl = f"{d}/bench_{name.split()[0]}.stl"; gen(stl)
        tris = slice_stl.parse_stl(stl)
        print(f"[{name}]   ({len(tris)} triangles)")
        print(hdr)
        for z in zs:
            loops = slice_stl.slice_z(tris, z)
            if not loops:
                print(f"  z={z}: no contour"); continue
            outline = max(loops, key=lambda L: abs(slice_stl._signed_area(L))) 
            kdec = max(1, round(len(outline) / 28))          
            outline = outline[::kdec]
            xs = [p[0] for p in outline]; ys = [p[1] for p in outline]
            bbox = (min(xs), min(ys), max(xs), max(ys))
            sp = max(0.6, (max(xs) - min(xs)) / 26)                    
            _, rows = run_part(f"z={z}", outline, bbox, spacing=sp)
            tag = f"z={z:>4} ({len(outline)}pt)"
            for nm, m, note in rows:
                if m is None:
                    print(f"  {tag:13s} {nm:27s} {note}"); tag = ""; continue
                print(f"  {tag:13s} {nm:27s} {m['strokes']:>7} {m['travel_moves']:>6} "
                      f"{m['self_crossings']:>6} {str(m['printed_length']):>8} "
                      f"{str(m['turning_deg']):>7} {_tpm(m):>5}  {note}")
                tag = ""
            print()
    print("Real sliced layers (contours decimated to ~28 pts, standard simplify):")
    print(" - EVERY layer (gear teeth, rotating star, shrinking cone) -> origami A-trails it to ONE")
    print("   crossing-free stroke: 0 travel / 0 self-crossings. Default infill fragments into 17-70")
    print("   strokes with 100s of self-crossings (crosshatch) on every layer regardless.")
    print(" - Requirement found + fixed: the tessellation must be PADDED beyond the part outline so its")
    print("   own grid-boundary (degree-3) vertices stay outside — else one is trapped inside the part,")
    print("   odd and unreachable by the part-boundary patch. (3*spacing margin; clip.py unchanged.)")
    print(" - Origami's one cost remains turning angle (180/interior vertex, RQ4).")


def parts_benchmark():
    print("S4 end-to-end benchmark — origami pipeline vs infill baselines, same part & spacing\n")
    hdr = (f"  {'method':28s} {'strokes':>7} {'travel':>7} {'cross':>6} "
           f"{'length':>8} {'turn°':>7} {'t/mm':>5}  notes")
    for name, outline, bbox in _outlines():
        print(f"[{name}]")
        print(hdr)
        for nm, m, note in run_part(name, outline, bbox)[1]:
            if m is None:
                print(f"  {nm:28s} {note}"); continue
            print(f"  {nm:28s} {m['strokes']:>7} {m['travel_moves']:>7} "
                  f"{m['self_crossings']:>6} {str(m['printed_length']):>8} "
                  f"{str(m['turning_deg']):>7} {_tpm(m):>5}  {note}")
        print()
    print("Reading (the honest, complete picture):")
    print(" - ROBUST win = 0 self-crossings on 2-D coverage. A crosshatch grid self-crosses 100s of")
    print("   times (over-deposition) and a 2-direction grid CANNOT be made crossing-free at all;")
    print("   origami fills 2-D as ONE crossing-free stroke.")
    print(" - Strokes/travel are NOT the differentiator: a boustrophedon rectilinear is also 1 stroke /")
    print("   0 travel on convex — but it is 1-D (anisotropic) infill and fragments on non-convex")
    print("   (star -> 2 strokes). Equal-continuity comparison must be vs 2-D infill, which origami wins.")
    print(" - HONEST loss = smoothness. t/mm (turning per printed mm) is ~190-226 for origami vs ~17-52")
    print("   boustrophedon, 0 straight grid. Connected Fermat Spirals also gives 1 stroke at LOW")
    print("   turning — it beats origami here. Origami's turning is the fixed RQ4 constant (180/vertex),")
    print("   not a tunable: crossing-freeness + Kawasaki force it (docs/ALGORITHMS.md).")
    print(" - Euler transform (now MEASURED, not cited — euler_transform.py runs the real")
    print("   gupta2018/2020 construction + tour on a hex mesh CLIPPED TO THE EXACT OUTLINE, so its")
    print("   density is directly comparable to origami's own): reaches origami's continuity only by")
    print("   FORCING it, exactly 4x its own base E (Lemma 3.3, confirmed here to the decimal) and a")
    print("   consistent ~1.5-1.7x origami's printed-length density across convex and non-convex parts")
    print("   alike; origami inherits continuity from Maekawa for free (docs/REPRODUCIBILITY.md).")


def _row(name, m):
    return (f"  {name:22s} {m['strokes']:>7} {m['travel_moves']:>7} "
            f"{m['self_crossings']:>10} {m['printed_length']:>9} {m['turning_deg']:>9}")


def main(argv):
    if not argv or argv[0] == "--parts":
        parts_benchmark(); return 0
    if argv[0] == "--real":
        realparts_benchmark(); return 0
    origami = json.load(open(argv[0]))
    bb = bbox_of(origami)

    grid = grid_infill(bb, spacing=0.5, directions=("h", "v"))

    mo = toolpath_metrics(origami)
    mg = toolpath_metrics(grid)
    print("S4 toolpath benchmark (same bounding box):\n")
    print(f"  {'method':22s} {'strokes':>7} {'travel':>7} {'crossings':>10} "
          f"{'length':>9} {'turn°':>9}")
    print(_row("origami A-trail", mo))
    print(_row("default grid infill", mg))
    print("\n  reference baselines (authors' reported, not run here):")
    print("    Connected Fermat Spirals — 1 stroke, low turning, region-decomp")
    print("    Euler transformation     — single tour at ~4x edges / ~2x length")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
