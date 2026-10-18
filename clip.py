#!/usr/bin/env python3
"""S2 boundary-tracing clip: cut a crease pattern to an arbitrary 2D
outline so the result is printable (atrail.py --patch then yields an A-trail).

Unlike a naive keep-inside clip, this inserts a vertex at each crease∩outline
crossing and builds a real **border cycle** through {polygon vertices} ∪ {crossings}
in perimeter order. The crossings are degree-1 crease-ends *on* the boundary, so
the boundary patch pairs them *along* the outline (short, planarity-preserving)
instead of with chords that cut across a non-convex shape.

Handles any number of crossings per crease: the crease is split at every
crease∩outline crossing, and an inside sub-segment is kept only if it touches an
original interior vertex — a segment between two crossings (a disconnected sliver
in a concave notch) is dropped, which keeps the tour a single stroke. Crossings
are deduped at shared points, and the outline is given a tiny deterministic
off-lattice perturbation so cuts that pass exactly through tessellation vertices
(e.g. an axis-aligned cut on a grid) don't hit degeneracies. The boundary patch
(`atrail.py --patch`) then connects odd crease-ends *along* the border cycle.
Verified single crossing-free A-trail on star / rectangle / L-shape / hexagon.

Outlines with **holes** (interior boundaries) are handled by `clip_to_region`: the
outer outline and each hole become their own Border cycle, and the patch pairs odd
crease-ends along each loop (with a direct-chord fallback for the rare cross-loop
leftover). Validated to a single crossing-free A-trail on a square with one and two
holes and an L-shape with a hole.
"""
from __future__ import annotations

import math


def point_in_polygon(p, poly) -> bool:
    x, y = p[0], p[1]; inside = False
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def _seg_intersect(a, b, c, d):
    """Segment a→b vs c→d; return (t_on_ab, point) or None."""
    rx, ry = b[0] - a[0], b[1] - a[1]
    sx, sy = d[0] - c[0], d[1] - c[1]
    rxs = rx * sy - ry * sx
    if abs(rxs) < 1e-12:
        return None
    qpx, qpy = c[0] - a[0], c[1] - a[1]
    t = (qpx * sy - qpy * sx) / rxs
    u = (qpx * ry - qpy * rx) / rxs
    if -1e-9 <= t <= 1 + 1e-9 and -1e-9 <= u <= 1 + 1e-9:
        return t, (a[0] + t * rx, a[1] + t * ry)
    return None


def clip_to_polygon(raw: dict, poly) -> tuple[dict, int]:
    """Clip the M/V creases of ``raw`` to a simple outline ``poly`` (CCW (x,y))."""
    return clip_to_region(raw, poly, [])


def clip_to_region(raw: dict, outer, holes=()) -> tuple[dict, int]:
    """Clip the M/V creases of ``raw`` to the region inside ``outer`` (CCW) and
    OUTSIDE every loop in ``holes`` (interior boundaries). Each boundary loop
    (outer + each hole) becomes its own Border cycle, so the RQ3 patch
    (``atrail.py --patch``) pairs odd crease-ends *along each* loop. Returns
    (clipped_fold, n_anomalies)."""
    
    loops = [[(p[0] + 9.1e-4, p[1] + 5.3e-4) for p in lp]
             for lp in ([list(outer)] + [list(h) for h in holes])]

    C = [list(c[:2]) for c in raw["vertices_coords"]]
    EV = raw["edges_vertices"]
    EA = raw.get("edges_assignment", ["M"] * len(EV))

    def in_region(p):
        if not point_in_polygon(p, loops[0]):
            return False
        return not any(point_in_polygon(p, h) for h in loops[1:])
    inside = [in_region(c) for c in C]

    coords = [c[:] for c in C]
    out_edges, out_assign = [], []
    boundary = []                      
    for li, lp in enumerate(loops):     
        for i, p in enumerate(lp):
            boundary.append((li, float(i), len(coords))); coords.append([p[0], p[1]])

    def crossings_of(u, v):
        """All crease∩(any loop) crossings, sorted along the crease, deduped."""
        cs = []
        for li, lp in enumerate(loops):
            n = len(lp)
            for i in range(n):
                hit = _seg_intersect(C[u], C[v], lp[i], lp[(i + 1) % n])
                if hit:
                    t, pt = hit
                    e = lp[i]; f = lp[(i + 1) % n]
                    fl = (math.hypot(pt[0] - e[0], pt[1] - e[1]) /
                          max(1e-12, math.hypot(f[0] - e[0], f[1] - e[1])))
                    cs.append((t, pt, li, i + fl))
        cs.sort(key=lambda c: c[0])
        ded = []
        for c in cs:
            if not ded or abs(c[0] - ded[-1][0]) > 1e-9:
                ded.append(c)
        return ded

    anomalies = 0
    for (u, v), asg in zip(EV, EA):
        if asg not in ("M", "V", "U"):
            continue
        cs = crossings_of(u, v)
        if not cs:
            if inside[u] and inside[v]:          
                out_edges.append([u, v]); out_assign.append(asg)
            continue
        
        bps = [(0.0, u)]
        for t, pt, li, perim in cs:
            pid = len(coords); coords.append([pt[0], pt[1]])
            boundary.append((li, perim, pid)); bps.append((t, pid))
        bps.append((1.0, v))
        cur = inside[u]
        last = len(bps) - 2
        for k in range(len(bps) - 1):
            if cur and (k == 0 or k == last):
                out_edges.append([bps[k][1], bps[k + 1][1]]); out_assign.append(asg)
            cur = not cur
        if (len(cs) % 2) != (inside[u] != inside[v]):
            anomalies += 1                         

    boundary.sort()
    i = 0
    while i < len(boundary):
        li = boundary[i][0]; j = i
        while j < len(boundary) and boundary[j][0] == li:
            j += 1
        comp = [boundary[k][2] for k in range(i, j)]
        m = len(comp)
        for k in range(m):
            out_edges.append([comp[k], comp[(k + 1) % m]]); out_assign.append("B")
        i = j

    return {
        "vertices_coords": coords,
        "edges_vertices": out_edges,
        "edges_assignment": out_assign,
        "frame_classes": ["creasePattern"],
    }, anomalies
