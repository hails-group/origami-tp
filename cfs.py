#!/usr/bin/env python3
"""Connected-Fermat-Spiral-style infill, so
the S4 benchmark can report CFS's smoothness from a generated path instead of citing
Zhao et al. (2016).

Approach (stdlib only, robust on convex, non-convex, and holed regions):
  1. signed distance-to-boundary field on a grid (positive inside the region);
  2. inward offset contours = iso-level loops at depths s/2, 3s/2, ... (marching
     squares), i.e. the contour-parallel offsets spaced by the line spacing s;
  3. connect the nested loops into ONE low-curvature path by nearest-point bridges
     (a connected contour-parallel path; the Fermat double-spiral would additionally
     return both ends to the boundary, removing one end move — not needed to measure
     the turning, which is the point of the CFS comparison).

This is an approximation of Zhao et al.'s full region-decomposition CFS, sufficient to
*measure* the low turning that origami trades away (RQ4). Returns polylines (mm).
"""
from __future__ import annotations
import math


def _pip(x, y, poly):
    inside = False; n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i+1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside

def _dist_seg(px, py, ax, ay, bx, by):
    dx, dy = bx-ax, by-ay
    L2 = dx*dx + dy*dy
    t = 0.0 if L2 < 1e-12 else max(0.0, min(1.0, ((px-ax)*dx + (py-ay)*dy)/L2))
    qx, qy = ax + t*dx, ay + t*dy
    return math.hypot(px-qx, py-qy)

def _loops(outer, holes):
    return [outer] + list(holes)

def signed_distance_field(outer, holes, res):
    """Sample a grid (cell size `res`) of signed distance to the region boundary,
    positive inside (inside outer and outside every hole)."""
    loops = _loops(outer, holes)
    xs = [p[0] for lp in loops for p in lp]; ys = [p[1] for lp in loops for p in lp]
    x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
    nx = max(2, int(math.ceil((x1-x0)/res)))
    ny = max(2, int(math.ceil((y1-y0)/res)))
    edges = []
    for lp in loops:
        for i in range(len(lp)):
            a = lp[i]; b = lp[(i+1) % len(lp)]
            edges.append((a[0], a[1], b[0], b[1]))
    F = [[0.0]*(nx+1) for _ in range(ny+1)]
    def inside(px, py):
        if not _pip(px, py, outer): return False
        return not any(_pip(px, py, h) for h in holes)
    for j in range(ny+1):
        py = y0 + j*res
        for i in range(nx+1):
            px = x0 + i*res
            d = min(_dist_seg(px, py, *e) for e in edges)
            F[j][i] = d if inside(px, py) else -d
    return F, x0, y0, res, nx, ny

def _march(F, x0, y0, res, nx, ny, level):
    """Marching squares: return segments [(p,q)] of the iso-contour F==level."""
    def xy(i, j): return (x0 + i*res, y0 + j*res)
    def interp(i1, j1, i2, j2):
        a = F[j1][i1]; b = F[j2][i2]
        t = 0.5 if abs(b-a) < 1e-12 else (level - a) / (b - a)
        x1, y1 = xy(i1, j1); x2, y2 = xy(i2, j2)
        return (x1 + t*(x2-x1), y1 + t*(y2-y1))
    segs = []
    for j in range(ny):
        for i in range(nx):
            tl = F[j+1][i] >= level; tr = F[j+1][i+1] >= level
            br = F[j][i+1] >= level;  bl = F[j][i] >= level
            code = (tl<<3) | (tr<<2) | (br<<1) | bl
            if code in (0, 15): continue
            top    = lambda: interp(i, j+1, i+1, j+1)
            bottom = lambda: interp(i, j,   i+1, j)
            left   = lambda: interp(i, j,   i,   j+1)
            right  = lambda: interp(i+1, j, i+1, j+1)
            E = {1:[(left,bottom)],2:[(bottom,right)],3:[(left,right)],4:[(top,right)],
                 6:[(top,bottom)],7:[(left,top)],8:[(left,top)],9:[(top,bottom)],
                 11:[(top,right)],12:[(left,right)],13:[(bottom,right)],14:[(left,bottom)],
                 5:[(left,top),(bottom,right)],10:[(left,bottom),(top,right)]}
            for fa, fb in E.get(code, []):
                segs.append((fa(), fb()))
    return segs

def _stitch(segs, tol):
    """Stitch unordered segments into closed loops by endpoint matching."""
    used = [False]*len(segs)
    def close(a, b): return abs(a[0]-b[0]) < tol and abs(a[1]-b[1]) < tol
    loops = []
    for s in range(len(segs)):
        if used[s]: continue
        used[s] = True
        loop = [segs[s][0], segs[s][1]]
        extended = True
        while extended:
            extended = False
            for k in range(len(segs)):
                if used[k]: continue
                p, q = segs[k]
                if close(loop[-1], p): loop.append(q); used[k]=True; extended=True
                elif close(loop[-1], q): loop.append(p); used[k]=True; extended=True
                if extended: break
        loops.append(loop)
    return loops

def _simplify(poly, eps):
    """Ramer-Douglas-Peucker: drop the marching-squares staircase, recovering the
    true (smooth) contour so turning reflects geometry, not grid quantization."""
    if len(poly) < 3: return poly
    ax, ay = poly[0]; bx, by = poly[-1]
    dx, dy = bx-ax, by-ay; L = math.hypot(dx, dy)
    dmax, idx = 0.0, 0
    for i in range(1, len(poly)-1):
        px, py = poly[i]
        d = (abs(dy*px - dx*py + bx*ay - by*ax) / L) if L > 1e-9 else math.hypot(px-ax, py-ay)
        if d > dmax: dmax, idx = d, i
    if dmax > eps:
        left = _simplify(poly[:idx+1], eps); right = _simplify(poly[idx:], eps)
        return left[:-1] + right
    return [poly[0], poly[-1]]

def offset_contours(outer, holes, spacing, res=None):
    res = res or spacing/6.0
    F, x0, y0, r, nx, ny = signed_distance_field(outer, holes, res)
    maxd = max(F[j][i] for j in range(ny+1) for i in range(nx+1))
    loops = []
    level = spacing/2.0
    while level < maxd:
        segs = _march(F, x0, y0, r, nx, ny, level)
        for lp in _stitch(segs, tol=res*1.5):
            if len(lp) >= 4:
                loops.append((level, _simplify(lp, eps=res*1.2)))
        level += spacing
    return loops

def _path_len(p):
    return sum(math.hypot(p[i+1][0]-p[i][0], p[i+1][1]-p[i][1]) for i in range(len(p)-1))

def connect(loops):
    """Connect nested offset loops into one low-curvature path by nearest-point bridges
    (outer to inner). A connected contour-parallel path: low turning (the CFS property
    we measure). It is not guaranteed crossing-free here -- the raster marching-squares
    contours and the simple bridges leave a few crossings that a faithful Fermat reroute
    on exact polygon offsets would avoid; we therefore report turning, the measured
    quantity, and treat crossing-freeness as the cited CFS property (Zhao et al.)."""
    if not loops: return []
    items = [lp for _, lp in sorted(loops, key=lambda t: t[0])]  
    used = [False]*len(items)
    path = list(items[0]); used[0] = True
    for _ in range(len(items)-1):
        ex, ey = path[-1]
        best, bk, bidx = 1e30, -1, 0
        for k in range(len(items)):
            if used[k]: continue
            for idx, (qx, qy) in enumerate(items[k]):
                d = (qx-ex)**2 + (qy-ey)**2
                if d < best: best, bk, bidx = d, k, idx
        if bk < 0: break
        path += items[bk][bidx:] + items[bk][:bidx]             
        used[bk] = True
    return path

def cfs_infill(outer, holes, spacing, res=None):
    """One connected contour-parallel (Fermat-style) path filling the region."""
    return connect(offset_contours(outer, holes, spacing, res))


def metrics(path):
    """Total turning (deg) and t/mm of a polyline."""
    if len(path) < 3: return {"length": _path_len(path), "turning_deg": 0.0, "t_per_mm": 0.0}
    tot = 0.0; prev = None
    for i in range(len(path)-1):
        a = math.atan2(path[i+1][1]-path[i][1], path[i+1][0]-path[i][0])
        if prev is not None:
            tot += math.degrees(abs((a - prev + math.pi) % (2*math.pi) - math.pi))
        prev = a
    L = _path_len(path)
    return {"length": L, "turning_deg": tot, "t_per_mm": (tot/L if L > 1e-9 else 0.0)}


if __name__ == "__main__":
    star = [[5 + 4.6*math.cos(math.radians(90+k*36))*(1 if k%2==0 else 0.42),
             5 + 4.6*math.sin(math.radians(90+k*36))*(1 if k%2==0 else 0.42)] for k in range(10)]
    cases = [("square 20", [[0,0],[20,0],[20,20],[0,20]], []),
             ("L-shape",   [[0,0],[20,0],[20,8],[8,8],[8,20],[0,20]], []),
             ("5-star",    star, []),
             ("square+hole", [[0,0],[20,0],[20,20],[0,20]], [[[7,7],[13,7],[13,13],[7,13]]])]
    def _cross(a, b, c, d):
        def o(p, q, r): return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])
        for p in (a, b):
            for q in (c, d):
                if abs(p[0]-q[0]) < 1e-7 and abs(p[1]-q[1]) < 1e-7: return False
        return ((o(c,d,a) > 0) != (o(c,d,b) > 0)) and ((o(a,b,c) > 0) != (o(a,b,d) > 0))
    def self_cross(p):
        s = [(p[i], p[i+1]) for i in range(len(p)-1)]
        return sum(1 for i in range(len(s)) for j in range(i+2, len(s))
                   if _cross(s[i][0], s[i][1], s[j][0], s[j][1]))
    print(f"{'part':14}{'loops':>6}{'len(mm)':>10}{'turn(deg)':>11}{'t/mm':>8}{'self-x':>8}")
    for name, o, h in cases:
        loops = offset_contours(o, h, spacing=1.0)
        path = connect(loops)
        m = metrics(path)
        print(f"{name:14}{len(loops):>6}{m['length']:>10.1f}{m['turning_deg']:>11.0f}"
              f"{m['t_per_mm']:>8.1f}{self_cross(path):>8}")
