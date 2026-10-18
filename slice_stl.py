#!/usr/bin/env python3
"""Minimal stdlib STL slicer (mesh → 2D cross-section contours), plus a few
parametric solid-part generators, so the S4 benchmark can run on *real sliced layers* rather
than hand-drawn outlines. 

  parse_stl(path)        -> list of triangles [(p0,p1,p2), ...], each p=(x,y,z)
  slice_z(tris, z)       -> list of contour polygons (CCW outer) at height z
  gen_gear_stl / gen_cone_stl / gen_vase_stl  -> write an ASCII STL of a solid part

Slicing: intersect each triangle with the plane z=z0 (edges whose endpoints straddle z0 give a
segment), then stitch segments end-to-end into closed loops (endpoints quantized). Single outer
contour per loop; holes (nested loops) are returned separately and currently unused downstream
(clip.py doesn't yet handle interior boundaries)."""
from __future__ import annotations
import math, struct
from pathlib import Path

def parse_stl(path):
    data = Path(path).read_bytes()
    tris = []
    if data[:5] == b"solid" and b"facet" in data[:2000]:
        # ASCII
        nums = []
        for line in data.decode("ascii", "replace").splitlines():
            line = line.strip()
            if line.startswith("vertex"):
                _, x, y, z = line.split()[:4]
                nums.append((float(x), float(y), float(z)))
        for i in range(0, len(nums) - 2, 3):
            tris.append((nums[i], nums[i + 1], nums[i + 2]))
        return tris
    n = struct.unpack_from("<I", data, 80)[0]
    off = 84
    for _ in range(n):
        vals = struct.unpack_from("<9f", data, off + 12)   
        p0 = (vals[0], vals[1], vals[2]); p1 = (vals[3], vals[4], vals[5])
        p2 = (vals[6], vals[7], vals[8]); tris.append((p0, p1, p2)); off += 50
    return tris


def write_ascii_stl(path, tris):
    out = ["solid part"]
    for p0, p1, p2 in tris:
        out.append(" facet normal 0 0 0\n  outer loop")
        for p in (p0, p1, p2):
            out.append(f"   vertex {p[0]:.6g} {p[1]:.6g} {p[2]:.6g}")
        out.append("  endloop\n endfacet")
    out.append("endsolid part\n")
    Path(path).write_text("\n".join(out))


def _plane_seg(tri, z):
    """Segment where triangle crosses plane z=z0, or None."""
    pts = []
    for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
        za, zb = a[2], b[2]
        if (za > z) != (zb > z):                      
            t = (z - za) / (zb - za)
            pts.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return (pts[0], pts[1]) if len(pts) == 2 else None


def slice_z(tris, z, q=1e4):
    segs = []
    for t in tris:
        s = _plane_seg(t, z)
        if s:
            segs.append(s)
    def key(p): return (round(p[0] * q) / q, round(p[1] * q) / q)
    adj = {}
    for a, b in segs:
        adj.setdefault(key(a), []).append((key(b), b))
        adj.setdefault(key(b), []).append((key(a), a))
    seen_edge = set(); loops = []
    for a, b in segs:
        ka, kb = key(a), key(b)
        if (ka, kb) in seen_edge or (kb, ka) in seen_edge:
            continue
        loop = [a]; cur, prev = kb, ka; seen_edge.add((ka, kb))
        pt = b
        while cur != ka and len(loop) < len(segs) + 2:
            loop.append(pt)
            nxts = [(k, p) for (k, p) in adj.get(cur, []) if k != prev]
            if not nxts:
                break
            nk, npt = nxts[0]
            seen_edge.add((cur, nk)); prev, cur, pt = cur, nk, npt
        if len(loop) >= 3:
            loops.append(_ccw(loop))
    return loops


def _signed_area(poly):
    s = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]; s += x1 * y2 - x2 * y1
    return s / 2


def _ccw(poly):
    return poly if _signed_area(poly) > 0 else poly[::-1]


def _extrude(profile, h):
    """Solid prism from a CCW 2D profile, z in [0,h]. Returns triangles (caps + walls)."""
    n = len(profile); tris = []
    cx = sum(p[0] for p in profile) / n; cy = sum(p[1] for p in profile) / n
    bot = [(x, y, 0.0) for x, y in profile]; top = [(x, y, h) for x, y in profile]
    cb = (cx, cy, 0.0); ct = (cx, cy, h)
    for i in range(n):
        j = (i + 1) % n
        tris.append((cb, bot[j], bot[i]))             
        tris.append((ct, top[i], top[j]))             
        tris.append((bot[i], bot[j], top[j]))      
        tris.append((bot[i], top[j], top[i]))
    return tris


def gear_profile(teeth=12, R=10.0, root=7.5, tip_frac=0.42):
    """Simple involute-free cog: alternating tip/root arcs with radial flanks."""
    pts = []; step = 2 * math.pi / teeth
    for k in range(teeth):
        a0 = k * step
        tw = step * tip_frac
        for ang, rad in ((a0, root), (a0 + (step - tw) / 2, R),
                         (a0 + (step + tw) / 2, R), (a0 + step - 1e-6, root)):
            pts.append((rad * math.cos(ang), rad * math.sin(ang)))
    return pts


def gen_gear_stl(path, teeth=12, R=10.0, root=7.5, h=4.0):
    write_ascii_stl(path, _extrude(gear_profile(teeth, R, root), h)); return path


def gen_cone_stl(path, R=10.0, rtop=1.5, h=12.0, n=48, layers=16):
    """Truncated cone (frustum), ring-stacked: circular cross-sections shrinking linearly with z."""
    def rad(z): return R + (rtop - R) * (z / h)
    zs = [h * k / layers for k in range(layers + 1)]
    rings = [[(rad(z) * math.cos(2 * math.pi * i / n), rad(z) * math.sin(2 * math.pi * i / n), z)
              for i in range(n)] for z in zs]
    tris = []; cb = (0, 0, zs[0]); ct = (0, 0, zs[-1])
    for i in range(n):
        j = (i + 1) % n
        tris.append((cb, rings[0][j], rings[0][i]))
        tris.append((ct, rings[-1][i], rings[-1][j]))
        for k in range(layers):
            a, b = rings[k][i], rings[k][j]; c, d = rings[k + 1][i], rings[k + 1][j]
            tris.append((a, b, d)); tris.append((a, d, c))
    write_ascii_stl(path, tris); return path


def gen_vase_stl(path, h=20.0, n=48, layers=24):
    """Revolved 'vase': radius varies sinusoidally with z — non-trivial varying convex layers."""
    def rad(z): return 6.0 + 2.5 * math.sin(z / h * 2 * math.pi) + 0.04 * z
    tris = []
    zs = [h * k / layers for k in range(layers + 1)]
    rings = [[(rad(z) * math.cos(2 * math.pi * i / n), rad(z) * math.sin(2 * math.pi * i / n), z)
              for i in range(n)] for z in zs]
    cb = (0, 0, zs[0]); ct = (0, 0, zs[-1])
    for i in range(n):
        j = (i + 1) % n
        tris.append((cb, rings[0][j], rings[0][i]))
        tris.append((ct, rings[-1][i], rings[-1][j]))
        for k in range(layers):
            a, b = rings[k][i], rings[k][j]; c, d = rings[k + 1][i], rings[k + 1][j]
            tris.append((a, b, d)); tris.append((a, d, c))
    write_ascii_stl(path, tris); return path


def gen_twist_star_stl(path, points=5, R=11.0, r=4.6, h=16.0, turn=60.0, layers=20):
    """Star prism twisted about z (rotates `turn`° over the height): each layer is a rotated
    non-convex star — varying non-convex cross-sections, the case grid infill handles worst."""
    def ring(z):
        rot = math.radians(turn) * (z / h)
        pts = []
        for k in range(2 * points):
            a = rot + math.pi / 2 + k * math.pi / points
            rad = R if k % 2 == 0 else r
            pts.append((rad * math.cos(a), rad * math.sin(a), z))
        return pts
    zs = [h * k / layers for k in range(layers + 1)]
    rings = [ring(z) for z in zs]; n = 2 * points
    tris = []; cb = (0, 0, zs[0]); ct = (0, 0, zs[-1])
    for i in range(n):
        j = (i + 1) % n
        tris.append((cb, rings[0][j], rings[0][i]))
        tris.append((ct, rings[-1][i], rings[-1][j]))
        for k in range(layers):
            a, b = rings[k][i], rings[k][j]; c, d = rings[k + 1][i], rings[k + 1][j]
            tris.append((a, b, d)); tris.append((a, d, c))
    write_ascii_stl(path, tris); return path


if __name__ == "__main__":
    import tempfile
    d = tempfile.gettempdir()
    for name, fn in [("gear", lambda p: gen_gear_stl(p, teeth=12, R=10, root=7.2, h=4)),
                     ("cone", lambda p: gen_cone_stl(p, R=10, h=12)),
                     ("vase", lambda p: gen_vase_stl(p))]:
        p = f"{d}/{name}.stl"; fn(p); tris = parse_stl(p)
        zc = 2.0 if name != "vase" else 7.0
        loops = slice_z(tris, zc)
        sizes = [len(l) for l in loops]
        print(f"{name}.stl: {len(tris)} tris -> slice z={zc}: {len(loops)} contour(s), pts={sizes}")
