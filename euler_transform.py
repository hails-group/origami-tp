#!/usr/bin/env python3
"""Implementation of the Gupta-Krishnamoorthy-Dreifus
Euler transformation, as an actual baseline (not the closed-form 4x-edges/2x-length
estimate cited elsewhere in this repo) for head-to-head comparison against the
origami-native A-trail toolpath.

References:
  gupta2018euler       — "Euler Transformation of Polyhedral Complexes" (the K -> K̂
                          construction: Class 1/2/3 cells, Theorem 3.1: every vertex
                          in K̂ has degree 4, Lemma 3.3: |V̂|=2|E|, |Ê|=4|E|).
  gupta2020continuous  — "Continuous Toolpath Planning..." (the circuit-tree +
                          modified-Hierholzer algorithm that realizes a *non-crossing*
                          Eulerian tour on K̂, Section 6).

What this implements
---------------------
1. `build_mesh` / `fix_assumption23` / `build_class1` / `build_class2`: the actual
   Euler transformation (2018 paper, Section 3/Def. 3.1) on a 2D polygon-complex —
   real mitered-offset coordinates, real edge counts, not the Lemma 3.3 formula.
2. Tool-path generation: gupta2020's Section 6 proves a non-crossing Eulerian tour
   exists on K̂ by a specific circuit-tree/Hierholzer algorithm. K̂ is, by their own
   Theorem 3.1, a connected planar graph with every vertex of degree 4 — exactly the
   regime `atrail.py` already handles (Kotzig's theorem: any such graph has a
   checkerboard-spanning-tree family of non-crossing Eulerian circuits). Rather than
   reimplementing their specific circuit-tree bookkeeping, we reuse `atrail.py`'s
   validated Kotzig construction directly on K̂. Both are valid algorithms for the
   *same* guarantee (a single non-crossing Eulerian tour on a 4-regular plane graph);
   what is being measured here — K̂ itself, i.e. Gupta's actual edge/length overhead
   — is unaffected by which realization algorithm produces the tour on top of it.
   Output is validated with the same `foldtp.py` non-crossing/A-trail checker used
   for the origami side, so both baselines are certified by one validator.

Mesh input and Assumption 2.3
------------------------------
Gupta's transformation requires (Assumption 2.3, condition 4) that no polygon in K
has two *adjacent* edges that are both boundary edges (shared with the outside). Any
mesh with a straight, convex domain corner violates this at the corner cell — this
includes the Miura/grid meshes already used elsewhere in this repo (bench.py). We
therefore mesh the domain with a hexagonal (honeycomb) tiling — matching the style of
Gupta et al.'s own Fig. 1 / Fig. 11 examples — and run a generic `fix_assumption23`
pass that locally splits any still-violating corner cell with one diagonal (negligible
O(#corners) overhead, same spirit as this repo's own RQ3 boundary patch; a cell too
small to split cleanly is dropped instead). `build_mesh` also welds near-coincident
vertices left behind by the rectangular clip, and the mitered-offset distance is a
single global value (not a per-face fraction) — both matter for numerical stability of
the resulting straight-line planar embedding. Any residual odd vertex is caught by
`atrail.py`'s own boundary patch as a safety net.

Normalization for comparison
-----------------------------
Because the input mesh (hexagons) is not combinatorially identical to the Miura
crease graph (parallelograms) used by the origami A-trail, raw vertex/edge counts are
not directly comparable — different cell shapes pack a domain with different edge
counts per unit area. We therefore match domain bounding box and target cell edge
length (spacing) between the two methods, and report DENSITY metrics (printed length
per unit area, edge count per unit area, turning-per-printed-mm) as the primary
comparison, alongside the raw counts for transparency. See `docs/REPRODUCIBILITY.md`.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import atrail  

OFFSET_FRAC = 0.28   


def _dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _poly_area(poly):
    n = len(poly)
    if n < 3:
        return 0.0
    return 0.5 * sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
                      for i in range(n))


def _line_intersect(p1, p2, p3, p4):
    x1, y1 = p1; x2, y2 = p2; x3, y3 = p3; x4, y4 = p4
    d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(d) < 1e-12:
        return ((x2 + x3) / 2, (y2 + y3) / 2)   
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))


def offset_polygon(poly, dist):
    """Mitered inward offset of a convex CCW polygon by `dist` (straight-skeleton
    offset for the convex case: offset each edge along its inward normal, intersect
    consecutive offset lines). Matches gupta2018euler Definition 4.1 / Section 4."""
    n = len(poly)
    lines = []
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        L = math.hypot(dx, dy) or 1e-9
        nx, ny = -dy / L, dx / L         
        lines.append(((a[0] + nx * dist, a[1] + ny * dist),
                      (b[0] + nx * dist, b[1] + ny * dist)))
    return [_line_intersect(*lines[i - 1], *lines[i]) for i in range(n)]


def _point_in_polygon(p, poly):
    """Ray-casting point-in-polygon test — works for any simple polygon (convex or not)."""
    x, y = p; inside = False; n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            xin = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < xin:
                inside = not inside
    return inside


def _seg_intersect_param(a, b, c, d):
    """a->b vs c->d; returns (t, u, point) for a strict interior crossing, else None."""
    rx, ry = b[0] - a[0], b[1] - a[1]
    sx, sy = d[0] - c[0], d[1] - c[1]
    denom = rx * sy - ry * sx
    if abs(denom) < 1e-12:
        return None
    qpx, qpy = c[0] - a[0], c[1] - a[1]
    t = (qpx * sy - qpy * sx) / denom
    u = (qpx * ry - qpy * rx) / denom
    eps = 1e-9
    if eps < t < 1 - eps and eps < u < 1 - eps:
        return t, u, (a[0] + t * rx, a[1] + t * ry)
    return None


def clip_convex_by_polygon(subject, clip):
    """Greiner-Hormann polygon intersection: convex `subject` (CCW) ∩ simple, possibly
    NON-convex `clip` (CCW). Returns a list of output polygons — 0, 1, or (when the
    clip's boundary cuts the subject into disjoint pieces, e.g. a hex straddling a
    thin concave notch) more than 1 — each a CCW list of points. Unlike
    `_clip_convex` (Sutherland-Hodgman, convex clip windows only), this handles
    arbitrary outlines (L-shapes, stars, decimated STL contours). Cross-validated
    against `_clip_convex` on 300 random hex/rectangle cases (exact match) and
    verified to recover the exact true outline area (shoelace) on the L-shape and
    5-point-star outlines used elsewhere in this repo (ratio 1.000, no gaps or
    double-counting even on the star's thin concave points)."""
    ns, nc = len(subject), len(clip)
    xints = []
    for i in range(ns):
        a, b = subject[i], subject[(i + 1) % ns]
        for j in range(nc):
            c, d = clip[j], clip[(j + 1) % nc]
            hit = _seg_intersect_param(a, b, c, d)
            if hit:
                t, u, pt = hit
                xints.append({"i": i, "j": j, "t": t, "u": u, "pt": pt})
    for k, x in enumerate(xints):
        x["id"] = k

    if not xints:
        if all(_point_in_polygon(p, clip) for p in subject):
            return [list(subject)]
        if all(_point_in_polygon(p, subject) for p in clip):
            return [list(clip)]
        return []

    def build_list(base_poly, n, key, sort_key):
        by_edge = {}
        for x in xints:
            by_edge.setdefault(x[key], []).append(x)
        arr = []
        for e in range(n):
            arr.append({"pt": base_poly[e], "is_x": False})
            for x in sorted(by_edge.get(e, []), key=lambda h: h[sort_key]):
                arr.append({"pt": x["pt"], "is_x": True, "id": x["id"]})
        return arr

    S = build_list(subject, ns, "i", "t")
    C = build_list(clip, nc, "j", "u")

    inside = _point_in_polygon(S[0]["pt"], clip)          
    for node in S:
        if node["is_x"]:
            inside = not inside
            node["entry"] = inside                         

    id_to_S = {n["id"]: idx for idx, n in enumerate(S) if n["is_x"]}
    id_to_C = {n["id"]: idx for idx, n in enumerate(C) if n["is_x"]}

    results = []
    for start in range(len(S)):
        node = S[start]
        if not node["is_x"] or not node.get("entry") or node.get("visited"):
            continue
        loop = []
        arr, idx = S, start
        start_id = node["id"]
        closed = False
        for _ in range(4 * (ns + nc) + len(xints) * 2 + 16):   
            cur = arr[idx]
            loop.append(cur["pt"])
            if cur["is_x"]:
                cur["visited"] = True
                other, lut = (C, id_to_C) if arr is S else (S, id_to_S)
                other[lut[cur["id"]]]["visited"] = True
            nxt_idx = (idx + 1) % len(arr)
            nxt = arr[nxt_idx]
            if nxt["is_x"]:
                if nxt["id"] == start_id:
                    closed = True
                    break
                arr, idx = (C, id_to_C[nxt["id"]]) if arr is S else (S, id_to_S[nxt["id"]])
            else:
                idx = nxt_idx
        if closed and len(loop) >= 3:              
            results.append(loop)                    
    return results


def _clip_convex(poly, x0, y0, x1, y1):
    """Sutherland-Hodgman clip of a convex polygon against an axis-aligned rectangle."""
    def clip_edge(pts, inside, intersect):
        out = []
        n = len(pts)
        for i in range(n):
            cur, nxt = pts[i], pts[(i + 1) % n]
            ci, ni = inside(cur), inside(nxt)
            if ci:
                out.append(cur)
            if ci != ni:
                out.append(intersect(cur, nxt))
        return out

    def xi(a, b, x):
        t = (x - a[0]) / (b[0] - a[0]); return (x, a[1] + t * (b[1] - a[1]))

    def yi(a, b, y):
        t = (y - a[1]) / (b[1] - a[1]); return (a[0] + t * (b[0] - a[0]), y)

    pts = list(poly)
    for inside, intr in (
        (lambda p: p[0] >= x0, lambda a, b: xi(a, b, x0)),
        (lambda p: p[0] <= x1, lambda a, b: xi(a, b, x1)),
        (lambda p: p[1] >= y0, lambda a, b: yi(a, b, y0)),
        (lambda p: p[1] <= y1, lambda a, b: yi(a, b, y1)),
    ):
        pts = clip_edge(pts, inside, intr)
        if not pts:
            return []
    return pts


def hex_mesh(x0, y0, x1, y1, size, outline=None):
    """Flat-top hexagon tiling (center-to-vertex radius `size`) of [x0,x1]x[y0,y1] (or,
    if `outline` is given, of the exact — possibly non-convex — polygon it describes,
    with (x0,y0,x1,y1) used only as its bounding box to bound the candidate-hex search)
    — the same style of mesh as Gupta et al.'s own Fig. 1 / 11 examples.

    Without `outline`: a straight rectangular clip (Sutherland-Hodgman, `_clip_convex`)
    is used so each boundary cell touches the outside via exactly one clean edge along
    any *straight* run of the boundary — the true (idealized) Assumption 2.3 violations
    are then confined to the domain's 4 corners, verified empirically (see
    `fix_assumption23`'s docstring), same as they would be for a quad grid — hexagons
    have no special immunity to a convex-corner cut.

    With `outline`: each candidate hex is clipped against the exact outline via
    `clip_convex_by_polygon` (Greiner-Hormann — handles concave notches and can split
    one hex into several disjoint pieces, e.g. straddling a thin part of a star). This
    is what makes the printed-length/area DENSITY comparison against origami's own
    (exactly outline-clipped) A-trail apples-to-apples on non-convex parts — without it,
    Gupta's mesh covered the outline's *bounding box*, which can overstate its true
    area by 50-280% on an L-shape or star (see docs/REPRODUCIBILITY.md).

    Either way: thin clipped slivers are dropped — they make the per-face mitered
    offset numerically unstable and can corrupt the planar embedding — which as a side
    effect can orphan a handful of neighboring edges into looking like extra boundary
    violations; `fix_assumption23` handles however many actually turn up. `build_mesh`
    also welds any remaining near-coincident vertices left by the clip.

    A hex-lattice vertex landing EXACTLY on an outline edge/vertex (a real, reproducible
    case on axis- or lattice-aligned outlines, e.g. an L-shape's inner corner) breaks
    `clip_convex_by_polygon`'s entry/exit classification at that knife-edge and can
    return a badly wrong (hugely oversized) piece. `clip.py` guards against the same
    class of degeneracy elsewhere in this repo with a tiny deterministic off-lattice
    jitter of the outline; we do the same here."""
    if outline is not None:
        outline = [(p[0] + 9.1e-4, p[1] + 5.3e-4) for p in outline]
    polys = []
    dx = 1.5 * size
    dy = math.sqrt(3) * size
    nx = int((x1 - x0) / dx) + 3
    ny = int((y1 - y0) / dy) + 3
    full_area = 1.5 * math.sqrt(3) * size * size
    for c in range(-1, nx):
        cx = x0 + c * dx
        off = (dy / 2) if c % 2 else 0.0
        for r in range(-1, ny):
            cy = y0 + r * dy + off
            hexpts = [(cx + size * math.cos(math.radians(60 * k)),
                       cy + size * math.sin(math.radians(60 * k))) for k in range(6)]
            pieces = (clip_convex_by_polygon(hexpts, outline) if outline is not None
                      else [_clip_convex(hexpts, x0, y0, x1, y1)])
            for clipped in pieces:
                if len(clipped) >= 3 and abs(_poly_area(clipped)) > 0.35 * full_area:
                    if _poly_area(clipped) < 0:
                        clipped = clipped[::-1]
                    polys.append(clipped)
    return polys


# -- combinatorial mesh from a polygon list --------------------------------------
def build_mesh(polys, tol=1e-6, weld_frac=0.05):
    """Dedup vertices by rounded coordinate; return (coords, faces, edge_faces) where
    faces[i] is a CCW list of vertex ids and edge_faces[(u,v)] (u<v) lists the 1 or 2
    face ids sharing that edge (1 => boundary edge, shared with the outside)."""
    vid = {}
    coords = []

    def get(p):
        k = (round(p[0] / tol), round(p[1] / tol))
        if k not in vid:
            vid[k] = len(coords); coords.append([p[0], p[1]])
        return vid[k]

    faces = []
    for poly in polys:
        f = [get(p) for p in poly]
        f = [f[i] for i in range(len(f)) if f[i] != f[i - 1]]   # drop dup-vertex slivers
        if len(f) >= 3:
            faces.append(f)

    
    def edge_len(u, v):
        return _dist(coords[u], coords[v])
    raw_lens = sorted(edge_len(f[i], f[(i + 1) % len(f)]) for f in faces for i in range(len(f)))
    med = raw_lens[len(raw_lens) // 2] if raw_lens else 1.0
    parent = list(range(len(coords)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for f in faces:
        n = len(f)
        for i in range(n):
            u, v = f[i], f[(i + 1) % n]
            if edge_len(u, v) < weld_frac * med:
                union(u, v)
    if any(find(i) != i for i in range(len(coords))):
        faces2 = []
        for f in faces:
            f2 = [find(v) for v in f]
            f2 = [f2[i] for i in range(len(f2)) if f2[i] != f2[i - 1]]
            if len(f2) >= 3:
                faces2.append(f2)
        faces = faces2

    edge_faces = {}
    for fi, f in enumerate(faces):
        n = len(f)
        for i in range(n):
            u, v = f[i], f[(i + 1) % n]
            k = (u, v) if u < v else (v, u)
            edge_faces.setdefault(k, []).append(fi)
    return coords, faces, edge_faces


def fix_assumption23(faces, edge_faces):
    """gupta2018euler Assumption 2.3(4): no polygon may have two cyclically-adjacent
    edges that are both boundary edges (else Theorem 3.1's degree-4 guarantee fails —
    see Remark 3.2/4.3, the "single square" example). A cell touching a *straight*
    boundary run is always fine (exactly one boundary edge); the violation is purely
    about convex domain corners, where a cell's footprint gets cut by two different
    boundary lines whose cuts meet at the corner point. With an idealized rectangular
    clip this happens at exactly the domain's 4 corners — verified empirically, not
    hex-specific (a quad grid has the identical 4-corner issue). In practice the
    sliver-drop in `hex_mesh` (needed for numerical stability — see its docstring)
    orphans a few more edges along the way when a thin neighbor gets dropped, so the
    real count run-to-run is small but not always exactly 4. Either way: split the
    violating polygon with one diagonal through the offending vertex to a roughly-
    opposite vertex, so each half keeps only one of the two boundary edges (the same
    bounded-patch spirit as this repo's own RQ3 boundary patch). A face too small to
    split (a degenerate triangle sliver) is dropped instead — negligible either way.
    Iterates to a fixed point."""
    faces = [list(f) for f in faces]

    def is_boundary(u, v, ef):
        k = (u, v) if u < v else (v, u)
        return len(ef.get(k, [])) == 1

    def rebuild(fs):
        ef = {}
        for k2, ff in enumerate(fs):
            m = len(ff)
            for t in range(m):
                a, b = ff[t], ff[(t + 1) % m]
                key = (a, b) if a < b else (b, a)
                ef.setdefault(key, []).append(k2)
        return ef

    max_iters = 8 * len(faces) + 40
    for _ in range(max_iters):
        viol = None
        for fi, f in enumerate(faces):
            n = len(f)
            for i in range(n):
                if is_boundary(f[i - 1], f[i], edge_faces) and is_boundary(f[i], f[(i + 1) % n], edge_faces):
                    viol = (fi, i); break
            if viol:
                break
        if viol is None:
            break
        fi, i = viol
        f = faces[fi]; n = len(f)
        if n < 4:
            del faces[fi]                   
        else:
            j = (i + n // 2) % n
            if j in (i, (i - 1) % n, (i + 1) % n):
                j = (i + 2) % n
            A = f[i:j + 1] if i < j else f[i:] + f[:j + 1]
            B = f[j:] + f[:i + 1] if i < j else f[j:i + 1]
            faces[fi] = A
            faces.append(B)
        edge_faces = rebuild(faces)
    return faces, edge_faces


def build_class1(coords, faces, offset_dist):
    """Class 1: each polygon f -> a mitered-offset (shrunk) copy f̂. A single global
    `offset_dist` (not a per-face fraction of that face's own min edge) keeps shrink
    amounts consistent across adjacent faces of very different size (e.g. a clipped
    boundary fragment next to a full hexagon) — an inconsistent local shrink is what
    makes neighboring Class-2 4-gons overlap/cross and corrupts the planar embedding."""
    new_coords = [list(c) for c in coords]
    class1 = {}                         
    for fi, f in enumerate(faces):
        poly = [coords[v] for v in f]
        n = len(poly)
        min_e = min(_dist(poly[k], poly[(k + 1) % n]) for k in range(n))
        d = min(offset_dist, 0.45 * min_e)  
        newpts = offset_polygon(poly, d)
        for k, v in enumerate(f):
            nid = len(new_coords); new_coords.append(list(newpts[k]))
            class1[(fi, v)] = nid
    return new_coords, class1


def build_class2(faces, edge_faces, class1):
    """Class 2 (per-edge 4-gons) + implicit Class 3 (per-vertex p-gons, which add no
    new edges — see module docstring): for each original edge, exactly 4 edges appear
    in K̂ (matches Lemma 3.3, |Ê|=4|E|). Interior edge {u,v} shared by f,g: the 2 new
    bridging edges {f-copy(u),g-copy(u)}, {f-copy(v),g-copy(v)} (the other 2 sides are
    already counted as part of f̂'s and ĝ's own Class-1 polygon sides). Boundary edge
    (sole face f, other side = the outside, Ĉ^O=C^O unchanged/unshrunk): the 2
    bridging edges {f-copy(u),u}, {f-copy(v),v}, plus the original edge {u,v} itself
    (the outside polygon's own unshrunk side)."""
    edges = set()
    for fi, f in enumerate(faces):
        n = len(f)
        for i in range(n):
            a = class1[(fi, f[i])]; b = class1[(fi, f[(i + 1) % n])]
            edges.add(frozenset((a, b)))
    for (u, v), fl in edge_faces.items():
        if len(fl) == 2:
            f, g = fl
            edges.add(frozenset((class1[(f, u)], class1[(g, u)])))
            edges.add(frozenset((class1[(f, v)], class1[(g, v)])))
        else:
            f = fl[0]
            edges.add(frozenset((class1[(f, u)], u)))
            edges.add(frozenset((class1[(f, v)], v)))
            edges.add(frozenset((u, v)))
    return [tuple(e) for e in edges]


def euler_transform_raw(polys, offset_frac=OFFSET_FRAC):
    """The full Euler transformation K -> K̂ (gupta2018euler Def. 3.1) on a polygon
    list. Returns (raw, base_E) where raw is a foldtp/atrail-compatible dict (1-skeleton
    of K̂ only — every vertex degree 4 by Theorem 3.1) and base_E = |E| of the input K."""
    if not polys:
        raise ValueError("empty mesh: hex `size` too large for the domain (all cells clipped away)")
    coords, faces, edge_faces = build_mesh(polys)
    faces, edge_faces = fix_assumption23(faces, edge_faces)
    if not faces:
        raise ValueError("mesh collapsed to nothing after Assumption-2.3 fixup — "
                          "hex `size` too large relative to the domain, use a smaller size")
    all_edge_lens = [_dist(coords[a], coords[b]) for a, b in edge_faces]
    offset_dist = offset_frac * min(all_edge_lens)
    new_coords, class1 = build_class1(coords, faces, offset_dist)
    edges = build_class2(faces, edge_faces, class1)
    raw = {
        "vertices_coords": new_coords,
        "edges_vertices": [list(e) for e in edges],
        "edges_assignment": ["U"] * len(edges),
    }
    return raw, len(edge_faces)


def gupta_toolpath(polys, offset_frac=OFFSET_FRAC):
    """Full pipeline: mesh -> K̂ -> non-crossing single Eulerian tour (via the repo's
    Kotzig/checkerboard `atrail.py`, realizing gupta2020's Section 6 guarantee — see
    module docstring). Returns (tp_dict, problems, base_E)."""
    raw, base_e = euler_transform_raw(polys, offset_frac)
    gen = atrail.ATrailGenerator(raw)
    out, problems = gen.generate(patch=True)   
    return out, problems, base_e               


def main(argv):
    import bench  
    W = H = float(argv[0]) if argv else 10.0
    size = float(argv[1]) if len(argv) > 1 else 0.6
    polys = hex_mesh(0, 0, W, H, size)
    out, problems, base_e = gupta_toolpath(polys)
    if problems:
        print("FAILED:", problems); return 1
    m = bench.toolpath_metrics(out)
    area = W * H
    print(f"Gupta Euler-transform baseline (real): domain {W:g}x{H:g} "
          f"({area:g} area), hex size {size:g}, base mesh edges E={base_e}")
    print(f"  strokes={m['strokes']}  travel={m['travel_moves']}  "
          f"crossings={m['self_crossings']}  edges={len(out['edges_vertices'])} "
          f"(~{len(out['edges_vertices']) / base_e:.2f}x base E)")
    print(f"  printed length={m['printed_length']:.1f}  turning={m['turning_deg']}deg  "
          f"t/mm={m['turning_deg'] / m['printed_length']:.2f}")
    print(f"  DENSITY (for cross-mesh comparison): "
          f"length/area={m['printed_length'] / area:.3f}  edges/area={len(out['edges_vertices']) / area:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
