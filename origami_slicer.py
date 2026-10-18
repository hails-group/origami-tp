#!/usr/bin/env python3
"""STL/OBJ mesh -> sliced layers -> origami crossing-free
A-trail toolpath per layer -> G-code. This is the fabrication layer; it orchestrates the core tools (slice_stl, clip, atrail) and
emits standard G-code following the single-stroke origami toolpath.

The infill is a continuous crossing-free stroke when generation succeeds. Perimeter
walls are separate strokes. This is a research prototype; see README.md for limits.
"""
from __future__ import annotations
import math, sys, struct, time
from pathlib import Path

import slice_stl                                   
import clip                                        
import atrail                                      


def parse_obj(data: bytes):
    verts, tris = [], []
    for line in data.decode("utf-8", "replace").splitlines():
        if line.startswith("v "):
            p = line.split()
            verts.append((float(p[1]), float(p[2]), float(p[3])))
        elif line.startswith("f "):
            idx = []
            for tok in line.split()[1:]:
                k = int(tok.split("/")[0])
                idx.append(k - 1 if k > 0 else len(verts) + k)
            for i in range(1, len(idx) - 1):       
                tris.append((verts[idx[0]], verts[idx[i]], verts[idx[i + 1]]))
    return tris


def parse_stl_bytes(data: bytes):
    """STL (binary or ASCII) from raw bytes — no filesystem, so it works under Pyodide."""
    if data[:5] == b"solid" and b"facet" in data[:2000]:
        nums = []
        for line in data.decode("ascii", "replace").splitlines():
            s = line.strip()
            if s.startswith("vertex"):
                p = s.split()
                nums.append((float(p[1]), float(p[2]), float(p[3])))
        return [(nums[i], nums[i + 1], nums[i + 2]) for i in range(0, len(nums) - 2, 3)]
    n = struct.unpack_from("<I", data, 80)[0]; off = 84; tris = []
    for _ in range(n):
        v = struct.unpack_from("<9f", data, off + 12)         
        tris.append(((v[0], v[1], v[2]), (v[3], v[4], v[5]), (v[6], v[7], v[8]))); off += 50
    return tris


def parse_mesh(data: bytes, name: str):
    if isinstance(data, str):
        data = data.encode("latin-1")
    return parse_obj(data) if name.lower().endswith(".obj") else parse_stl_bytes(data)


def mesh_bounds(tris):
    xs = [p[0] for t in tris for p in t]; ys = [p[1] for t in tris for p in t]
    zs = [p[2] for t in tris for p in t]
    return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def prep_mesh(tris, *, up="z", fit_mm=60.0, bed=(110, 110)):
    """Orient (up='y' rotates Y-up models to Z-up), uniformly scale so the largest dimension is
    `fit_mm`, set min-Z to 0, and center XY on the bed. Handles normalized/Y-up graphics meshes
    (e.g. the Stanford bunny) as well as already-mm Z-up CAD parts."""
    if up == "y":                                    
        tris = [tuple((p[0], p[2], p[1]) for p in t) for t in tris]
    x0, y0, z0, x1, y1, z1 = mesh_bounds(tris)
    span = max(x1 - x0, y1 - y0, z1 - z0) or 1.0
    scale = fit_mm / span
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    bx, by = bed
    return [tuple(((p[0] - cx) * scale + bx,
                   (p[1] - cy) * scale + by,
                   (p[2] - z0) * scale) for p in t) for t in tris]


def miura(nx, ny, a, b, c):
    V, idx = [], {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            idx[(i, j)] = len(V); V.append([i * a, j * b + (i % 2) * c])
    EV, EA = [], []
    for i in range(nx + 1):
        for j in range(ny + 1):
            if i < nx: EV.append([idx[(i, j)], idx[(i + 1, j)]]); EA.append("M")
            if j < ny: EV.append([idx[(i, j)], idx[(i, j + 1)]]); EA.append("V")
    return {"vertices_coords": V, "edges_vertices": EV, "edges_assignment": EA}


def _line_isect(e0, e1):
    """Intersect two lines given as (px,py,ux,uy) point+unit-direction."""
    px0, py0, ux0, uy0 = e0; px1, py1, ux1, uy1 = e1
    den = ux0 * uy1 - uy0 * ux1
    if abs(den) < 1e-9:
        return None
    t = ((px1 - px0) * uy1 - (py1 - py0) * ux1) / den
    return (px0 + t * ux0, py0 + t * uy0)


def inset_polygon(poly, d):
    """Inset a CCW polygon inward by d (miter joins, clamped at sharp corners). Returns a new
    polygon or None if it collapses. Approximate — adequate for FDM wall offsets."""
    n = len(poly)
    if n < 3 or d <= 0:
        return [list(p) for p in poly] if d == 0 else None
    edges = []
    for i in range(n):
        ax, ay = poly[i]; bx, by = poly[(i + 1) % n]
        ex, ey = bx - ax, by - ay; L = math.hypot(ex, ey)
        if L < 1e-9:
            edges.append(None); continue
        ux, uy = ex / L, ey / L
        edges.append((ax - uy * d, ay + ux * d, ux, uy))   
    out = []
    for i in range(n):
        e0, e1 = edges[(i - 1) % n], edges[i]
        if e0 is None or e1 is None:
            out.append(list(poly[i])); continue
        p = _line_isect(e0, e1) or (e1[0], e1[1])
        ox, oy = poly[i]
        if math.hypot(p[0] - ox, p[1] - oy) > 3 * d:      
            p = ((e0[0] + e1[0]) / 2, (e0[1] + e1[1]) / 2)
        out.append([p[0], p[1]])
    if abs(slice_stl._signed_area(out)) < (2 * d) ** 2:     
        return None
    return out if slice_stl._signed_area(out) > 0 else out[::-1]


def layer_walls(outline, n_walls, line_width):
    """n_walls closed perimeter loops, inset (k+0.5)*line_width from the outline."""
    walls = []
    for k in range(n_walls):
        w = inset_polygon(outline, (k + 0.5) * line_width)
        if not w:
            break
        walls.append([list(p) for p in w] + [list(w[0])])   
    return walls


def layer_toolpath(outline, spacing):
    """Return a list of strokes (each a list of (x,y)) covering `outline` with the origami
    A-trail, or None if the clip/A-trail fails on this layer."""
    xs = [p[0] for p in outline]; ys = [p[1] for p in outline]
    x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
    mg = 3 * spacing                                  
    nx = max(2, int(round((x1 - x0 + 2 * mg) / spacing)))
    ny = max(2, int(round((y1 - y0 + 2 * mg) / spacing)))
    M = miura(nx, ny, spacing, spacing, spacing * 0.36)
    for v in M["vertices_coords"]:
        v[0] += x0 - mg; v[1] += y0 - mg
    try:
        clipped, _ = clip.clip_to_polygon(M, outline)
        tp, problems = atrail.ATrailGenerator(clipped).generate(patch=True)
        if problems:
            return None
    except Exception:
        return None
    C = tp["vertices_coords"]; pv = tp["tp:paths_vertices"]
    closed = tp.get("tp:paths_closed", [True] * len(pv))
    strokes = []
    for p, verts in enumerate(pv):
        pts = [(C[v][0], C[v][1]) for v in verts]
        if closed[p] and pts:
            pts.append(pts[0])
        strokes.append(pts)
    return strokes


def slice_to_gcode(tris, *, layer_height=0.3, spacing=1.2, line_width=0.45, n_walls=2,
                   filament_d=1.75, feed=1200, travel_feed=3000, max_layers=60):
    """Full pipeline: prepped mesh -> per-layer origami toolpath -> G-code string + stats."""
    _, _, zmin, _, _, zmax = mesh_bounds(tris)
    height = zmax - zmin
    n_layers = min(max_layers, max(1, int(height / layer_height)))
    fil_area = math.pi * (filament_d / 2) ** 2
    e_per_mm = (layer_height * line_width) / fil_area

    g = ["; origami single-stroke toolpath  (generated by origami_slicer.py)",
         "; infill uses an A-trail; perimeter walls are separate strokes",
         "G21 ; mm", "G90 ; absolute XYZ", "M83 ; relative E",
         "G92 E0", "; --- begin ---"]
    E = 0.0
    stats = {"layers": 0, "strokes": 0, "travels": 0, "print_mm": 0.0,
             "failed_layers": 0, "layer_strokes": [], "n_layers": n_layers}
    layers = []                                  
    px = py = None

    def emit_travel(x, y):
        nonlocal px, py
        g.append("G1 E-1.0 F1800 ; retract")
        g.append(f"G0 X{x:.3f} Y{y:.3f} F{travel_feed}")
        g.append("G1 E1.0 F1800 ; prime")
        px, py = x, y
        stats["travels"] += 1

    for li in range(n_layers):
        z = zmin + (li + 0.5) * layer_height
        loops = slice_stl.slice_z(tris, z)
        if not loops:
            continue
        outline = max(loops, key=lambda L: abs(slice_stl._signed_area(L)))
        kdec = max(1, round(len(outline) / 40))        
        outline = outline[::kdec]
        if len(outline) < 3:
            continue
        walls = layer_walls(outline, n_walls, line_width)
        infreg = inset_polygon(outline, n_walls * line_width) if n_walls > 0 else outline
        infill = None
        if infreg and abs(slice_stl._signed_area(infreg)) > (2 * spacing) ** 2:
            infill = layer_toolpath(infreg, spacing)
        strokes = list(walls) + (infill or [])
        if not strokes:
            stats["failed_layers"] += 1
            continue
        g.append(f"; layer {li}  z={z:.3f}  ({len(walls)} walls + {len(infill or [])} infill)")
        g.append(f"G1 Z{z + layer_height:.3f} F600")
        for pts in strokes:
            if len(pts) < 2:
                continue
            emit_travel(pts[0][0], pts[0][1])
            for (x, y) in pts[1:]:
                d = math.hypot(x - px, y - py)
                g.append(f"G1 X{x:.3f} Y{y:.3f} E{d * e_per_mm:.5f} F{feed}")
                stats["print_mm"] += d; px, py = x, y
        nwall = len(walls); ninf = len(infill or [])
        stats["layers"] += 1; stats["strokes"] += nwall + ninf
        stats["layer_strokes"].append(ninf)           
        layers.append({
            "z": round(z, 3), "nwall": nwall,
            "outline": [[round(x, 2), round(y, 2)] for (x, y) in outline],
            "strokes": [[[round(x, 2), round(y, 2)] for (x, y) in s] for s in strokes],
        })

    g += ["; --- end ---", "G1 E-2 F1800", "G91", "G0 Z5", "G90",
          "M104 S0", "M140 S0", "M84 ; motors off"]
    stats["print_mm"] = round(stats["print_mm"], 1)
    return "\n".join(g) + "\n", stats, layers


def run_slice(mesh_bytes, name="part.stl", *, layer_height=0.3, spacing=1.2,
              fit_mm=60.0, up="z", walls=2, max_layers=40):
    """Single entry for the web app (called from JS via Pyodide): raw mesh bytes + params ->
    {gcode, stats, layers} (JSON-ready)."""
    tris = parse_mesh(mesh_bytes, name)
    if not tris:
        return {"error": "no triangles parsed from file"}
    tris = prep_mesh(tris, up=up, fit_mm=float(fit_mm), bed=(110, 110))
    gcode, stats, layers = slice_to_gcode(
        tris, layer_height=float(layer_height), spacing=float(spacing),
        n_walls=int(walls), max_layers=int(max_layers))
    return {"gcode": gcode, "stats": stats, "layers": layers}


def main(argv=None):
    """Slice a mesh into an explicitly named G-code file."""
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="STL or OBJ input mesh")
    parser.add_argument("-o", "--output", required=True, type=Path, help="output G-code path")
    parser.add_argument("--layer-height", type=float, default=0.3)
    parser.add_argument("--spacing", type=float, default=1.2)
    parser.add_argument("--fit-mm", type=float, default=60.0)
    parser.add_argument("--up", choices=("x", "y", "z"), default="z")
    parser.add_argument("--walls", type=int, default=2)
    parser.add_argument("--max-layers", type=int, default=40)
    args = parser.parse_args(argv)
    if any(not math.isfinite(v) or v <= 0 for v in
           (args.layer_height, args.spacing, args.fit_mm)):
        parser.error("layer height, spacing and fit size must be positive finite numbers")
    if args.walls < 0 or args.max_layers < 1:
        parser.error("walls must be non-negative and max-layers must be positive")
    if args.mesh.resolve() == args.output.resolve():
        parser.error("input mesh and output must be different files")
    try:
        result = run_slice(args.mesh.read_bytes(), args.mesh.name,
                           layer_height=args.layer_height, spacing=args.spacing,
                           fit_mm=args.fit_mm, up=args.up, walls=args.walls,
                           max_layers=args.max_layers)
        if "error" in result:
            parser.error(result["error"])
        if not result["stats"]["layers"]:
            parser.error("no printable layers generated")
        args.output.write_text(result["gcode"])
    except (OSError, ValueError, IndexError, struct.error) as exc:
        parser.error(str(exc))
    print(json.dumps(result["stats"], indent=2))
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
