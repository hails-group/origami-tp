#!/usr/bin/env python3
"""Generate a crossing-free Eulerian toolpath (A-trail) for a crease
pattern, via the checkerboard spanning-tree construction.

This is the *generator* counterpart to ``foldtp.py`` (the validator): given a FOLD
crease pattern it computes the ``tp:paths`` ordered walk and the
``tp:vertices_transition`` system, so the result can be round-tripped through the
validator and confirmed Eulerian-complete / crossing-free / single-stroke.

Method (see docs/ALGORITHMS.md): the printed crease graph is
planar with even degree; 2-color its faces (front/back); form the black
checkerboard graph G_B (black faces joined at each degree-4 vertex's black
corners); any spanning tree of G_B picks a merge/separate transition at every
vertex, whose boundary is a single non-crossing Eulerian circuit — an A-trail.

Scope: printed subgraph must be connected with all even degrees (Kotzig regime).
Odd-degree vertices (clipped boundaries) need RQ3 patching first and are reported.
Degree-2 and degree-4 vertices are handled; higher even degrees raise NotImplemented.
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict, deque
from pathlib import Path

PRINTED = {"M", "V", "F", "U"}   


class ATrailGenerator:
    def __init__(self, raw: dict):
        self.raw = raw
        self.coords = raw["vertices_coords"]
        self.edges = raw["edges_vertices"]
        self.assign = raw.get("edges_assignment", ["U"] * len(self.edges))
        self.nv = len(self.coords)

        self.printed = [e for e, a in enumerate(self.assign) if a in PRINTED]
        self.pset = set(self.printed)
        self.vinc = defaultdict(list)
        for e in self.printed:
            u, v = self.edges[e]
            self.vinc[u].append(e); self.vinc[v].append(e)

        self.badj = defaultdict(set)
        for e, a in enumerate(self.assign):
            if a == "B":
                u, v = self.edges[e]
                self.badj[u].add(v); self.badj[v].add(u)


    def _border_loops(self) -> list[list[int]]:
        """Walk each border component once, returning its vertices in cyclic order."""
        loops, visited = [], set()
        for s in list(self.badj):
            if s in visited:
                continue
            order, prev, cur = [], None, s
            while cur is not None and cur not in visited:
                order.append(cur); visited.add(cur)
                nxts = [w for w in self.badj[cur] if w != prev and w not in visited]
                prev, cur = cur, (nxts[0] if nxts else None)
            loops.append(order)
        return loops

    def _pairing_plan(self, loops, phase: int):
        """Edges the patch would add for one of the two alternating matchings.

        On a loop with 2m severed ends there are exactly two perfect matchings whose
        boundary arcs have disjoint interiors: pair (q1,q2),(q3,q4),... or the shift
        (q2,q3),(q4,q5),...,(q2m,q1). Both are valid Q-joins and both stay planar;
        they differ in which components they merge, so we try both (RQ3)."""
        plan, leftovers = [], []
        for order in loops:
            op = [i for i, v in enumerate(order) if len(self.vinc.get(v, ())) % 2]
            if not op:
                continue
            if phase and len(op) >= 2:
                op = op[1:] + [op[0]]           
            k = 0
            while k + 1 < len(op):
                a, b = op[k], op[k + 1]
                if a < b:
                    plan += [(order[j], order[j + 1]) for j in range(a, b)]
                else:                           
                    plan += [(order[j], order[j + 1]) for j in range(a, len(order) - 1)]
                    plan.append((order[-1], order[0]))
                    plan += [(order[j], order[j + 1]) for j in range(0, b)]
                k += 2
            if len(op) % 2:
                leftovers.append(order[op[-1]])
        for k in range(0, len(leftovers) - 1, 2):
            plan.append((leftovers[k], leftovers[k + 1]))
        return plan

    def _components(self, extra=()) -> list[set[int]]:
        """Connected components of the printed graph, optionally with planned edges."""
        g = defaultdict(set)
        for e in self.printed:
            u, v = self.edges[e]; g[u].add(v); g[v].add(u)
        for u, v in extra:
            g[u].add(v); g[v].add(u)
        seen, out = set(), []
        for v in g:
            if v in seen:
                continue
            st, cc = [v], {v}; seen.add(v)
            while st:
                u = st.pop()
                for w in g[u]:
                    if w not in seen:
                        seen.add(w); cc.add(w); st.append(w)
            out.append(cc)
        out.sort(key=len, reverse=True)
        return out

    def _add_printed(self, a, b) -> None:
        e = len(self.edges)
        self.edges.append([a, b]); self.assign.append("U")
        if "edges_foldAngle" in self.raw:
            self.raw["edges_foldAngle"].append(0)
        self.printed.append(e); self.pset.add(e)
        self.vinc[a].append(e); self.vinc[b].append(e)

    def patch_boundary(self) -> int:
        """Restore even printed-degree by pairing odd vertices (clip-induced
        crease-ends) consecutively around the border cycle and adding one printed
        edge per pair, then repair connectivity if the clip stranded a component.

        Parity: of the two alternating matchings on each border loop we take the one
        leaving fewer components; both are Q-joins of the same order, so the choice
        is free. Connectivity: any component still separate is joined to the main one
        by DOUBLING the border arc between them, which adds 2 to every degree on the
        arc and so preserves parity, and which crosses no crease because it runs
        along the boundary. Overhead is O(boundary) ~ O(1/n) of the crease count for
        the parity patch, plus one short doubled arc per extra component, vs. the
        Euler transformation's ~4x. Returns the number of edges added."""
        if not any(len(self.vinc[v]) % 2 for v in self.vinc):
            return self._bridge_components()
        added = 0
        if self.badj:                             
            loops = self._border_loops()
            best = None
            for phase in (0, 1):
                plan = self._pairing_plan(loops, phase)
                ncomp = len(self._components(plan))
                if best is None or ncomp < best[0]:
                    best = (ncomp, plan)
            for a, b in best[1]:
                self._add_printed(a, b); added += 1
        else:                                
            odd = [v for v in self.vinc if len(self.vinc[v]) % 2]
            cx = sum(self.coords[v][0] for v in odd) / len(odd)
            cy = sum(self.coords[v][1] for v in odd) / len(odd)
            oc = sorted(odd, key=lambda v: math.atan2(self.coords[v][1] - cy,
                                                      self.coords[v][0] - cx))
            for k in range(0, len(oc) - 1, 2):
                self._add_printed(oc[k], oc[k + 1]); added += 1
        return added + self._bridge_components()

    def _bridge_components(self) -> int:
        """Join stranded components by doubling the border arc between them.

        A clip can disconnect the lattice (a part feature thinner than one cell
        strands an island); parity patching does not fix that, and a disconnected
        graph needs one stroke per component. Doubling a border arc preserves every
        degree's parity and stays clear of the interior creases."""
        if not self.badj:
            return 0
        added = 0
        while True:
            cc = self._components()
            if len(cc) < 2:
                return added
            main, other = cc[0], cc[1]
            starts = [v for v in other if v in self.badj]
            if not starts:
                return added                     
            prev, seen, dq, hit = {}, set(starts), deque(starts), None
            while dq:
                u = dq.popleft()
                if u in main:
                    hit = u; break
                for w in self.badj.get(u, ()):
                    if w not in seen:
                        seen.add(w); prev[w] = u; dq.append(w)
            if hit is None:
                return added
            path = [hit]
            while path[-1] in prev:
                path.append(prev[path[-1]])
            for i in range(len(path) - 1):
                self._add_printed(path[i], path[i + 1])
                self._add_printed(path[i], path[i + 1])    
                added += 2


    def other(self, e: int, v: int) -> int:
        a, b = self.edges[e]
        return b if a == v else a

    def rotation(self, v: int) -> list[int]:
        """CCW incident printed edges by angle (matches foldtp.rotation)."""
        x0, y0 = self.coords[v][0], self.coords[v][1]
        return sorted(self.vinc[v], key=lambda e:
                      math.atan2(self.coords[self.other(e, v)][1] - y0,
                                 self.coords[self.other(e, v)][0] - x0))

    def check_even_connected(self) -> list[str]:
        problems = []
        odd = [v for v in self.vinc if len(self.vinc[v]) % 2]
        if odd:
            problems.append(f"{len(odd)} odd-degree printed vertices (need RQ3 patch): "
                            f"{odd[:8]}{'...' if len(odd) > 8 else ''}")
        hi = [v for v in self.vinc if len(self.vinc[v]) not in (2, 4)]
        if hi:
            problems.append(f"{len(hi)} vertices with degree not in (2,4) — unsupported v1")
        return problems


    def trace_faces(self):
        """Return (faces, dart_face) where a dart is (v, e) = leave v along e, and
        each face is the orbit under next-dart (predecessor in ccw rotation at the
        arrival vertex). Faces border to the dart's left."""
        rot_idx = {v: {e: i for i, e in enumerate(self.rotation(v))} for v in self.vinc}
        rot = {v: self.rotation(v) for v in self.vinc}

        def nxt(dart):
            v, e = dart
            w = self.other(e, v)
            r = rot[w]; i = rot_idx[w][e]
            return (w, r[(i - 1) % len(r)])          

        darts = [(v, e) for e in self.printed for v in self.edges[e]]
        seen = set(); faces = []; dart_face = {}
        for d in darts:
            if d in seen:
                continue
            fid = len(faces); face = []
            cur = d
            while cur not in seen:
                seen.add(cur); dart_face[cur] = fid; face.append(cur)
                cur = nxt(cur)
            faces.append(face)
        return faces, dart_face


    def two_color(self, faces, dart_face):
        adj = defaultdict(set)
        for e in self.printed:
            u, v = self.edges[e]
            f1, f2 = dart_face[(u, e)], dart_face[(v, e)]
            if f1 != f2:
                adj[f1].add(f2); adj[f2].add(f1)
        color = {}
        for s in range(len(faces)):              
            if s in color:
                continue
            color[s] = 0; q = deque([s])
            while q:
                u = q.popleft()
                for w in adj[u]:
                    if w not in color:
                        color[w] = 1 - color[u]; q.append(w)
        ok = all(color[u] != color[w] for u in adj for w in adj[u])
        return color, ok


    def build_gb(self, faces, dart_face, color):
        gb = defaultdict(set); gb_vertex = {}
        for v in self.vinc:
            if len(self.vinc[v]) != 4:
                continue
            r = self.rotation(v)
            corners = [dart_face[(v, e)] for e in r]      
            blacks = [c for c in corners if color[c] == 0]
            bset = list(dict.fromkeys(blacks))
            if len(bset) == 2:
                a, b = bset
                gb[a].add(b); gb[b].add(a); gb_vertex[frozenset((a, b))] = v
        return gb, gb_vertex

    def spanning_tree(self, gb):
        """Spanning FOREST of G_B — one tree per connected component, so a clipped
        (disconnected) pattern yields one A-trail stroke per component."""
        seen = set(); tree = set()
        for root in gb:
            if root in seen:
                continue
            seen.add(root); q = deque([root])
            while q:
                u = q.popleft()
                for w in gb[u]:
                    if w not in seen:
                        seen.add(w); tree.add(frozenset((u, w))); q.append(w)
        return tree


    def transitions(self, faces, dart_face, color, gb_vertex, tree):
        """Per vertex: a perfect matching (list of [eA,eB]) of incident printed edges."""
        merge_v = {gb_vertex[e] for e in tree}      
        trans = {}
        for v in self.vinc:
            r = self.rotation(v)
            d = len(r)
            if d == 2:
                trans[v] = [[r[0], r[1]]]
                continue
            if d != 4:
                raise NotImplementedError(f"degree {d} vertex {v}")
            corners = [dart_face[(v, r[k])] for k in range(4)]
            black_corner = [k for k in range(4) if color[corners[k]] == 0]  
            def pairing(corner_ks):
                return [[r[k], r[(k + 1) % 4]] for k in corner_ks]
            white_corner = [k for k in range(4) if k not in black_corner]
            merge = v in merge_v
            ks = white_corner if merge else black_corner  
            trans[v] = pairing(ks)
        return trans


    def trace_paths(self, trans):
        partner = {}
        for v, prs in trans.items():
            partner[v] = {}
            for a, b in prs:
                partner[v][a] = b; partner[v][b] = a
        used = set(); paths = []
        for start in self.printed:
            if start in used:
                continue
            e = start; v = self.edges[e][1]           
            verts = [self.edges[e][0]]; eseq = []
            while e not in used:
                used.add(e); eseq.append(e); verts.append(v)
                out = partner[v][e]
                v = self.other(out, v); e = out
            paths.append((verts, eseq))
        return paths


    def generate(self, patch: bool = False) -> tuple[dict, list[str]]:
        if patch:
            self.patch_boundary()
        problems = self.check_even_connected()
        if problems:
            return {}, problems
        faces, dart_face = self.trace_faces()
        color, ok = self.two_color(faces, dart_face)
        if not ok:
            return {}, ["face 2-coloring failed (graph not Eulerian-planar?)"]
        gb, gb_vertex = self.build_gb(faces, dart_face, color)
        tree = self.spanning_tree(gb)
        trans = self.transitions(faces, dart_face, color, gb_vertex, tree)
        paths = self.trace_paths(trans)
        out = dict(self.raw)
        out["frame_classes"] = list(dict.fromkeys(
            (self.raw.get("frame_classes") or []) + ["tp:toolpath"]))
        out["tp:paths_vertices"] = [p[0][:-1] for p in paths]   
        out["tp:paths_edges"] = [p[1] for p in paths]
        out["tp:paths_closed"] = [True] * len(paths)
        out["tp:paths_travel"] = [[False] * len(p[1]) for p in paths]
        out["tp:paths_role"] = ["infill"] * len(paths)
        out["tp:paths_source"] = ["checkerboard-tree"] * len(paths)
        vt = [[] for _ in range(self.nv)]
        for v, prs in trans.items():
            vt[v] = [[int(a), int(b)] for a, b in prs]
        out["tp:vertices_transition"] = vt
        return out, []


def main(argv):
    patch = "--patch" in argv
    argv = [a for a in argv if a != "--patch"]
    if not argv:
        print("usage: atrail.py [--patch] <crease.fold> [out.fold]"); return 2
    raw = json.load(open(argv[0]))
    gen = ATrailGenerator(raw)
    n_crease = len(gen.printed)
    out, problems = gen.generate(patch=patch)
    if problems:
        print(f"cannot generate single A-trail for {argv[0]}:")
        for p in problems:
            print(f"  - {p}")
        if not patch:
            print("  (try --patch to add an RQ3 boundary patch)")
        return 1
    npaths = len(out["tp:paths_edges"])
    added = len(gen.printed) - n_crease
    extra = f" (+{added} RQ3 patch edges, {added/n_crease*100:.1f}% overhead)" if added else ""
    print(f"generated A-trail for {argv[0]}: {npaths} path(s), "
          f"{sum(len(e) for e in out['tp:paths_edges'])} printed edges{extra}")
    if len(argv) > 1:
        json.dump(out, open(argv[1], "w"), indent=1)
        print(f"  wrote {argv[1]}")
    # self-validate if foldtp is importable
    try:
        from foldtp import FoldTP
        rep = FoldTP(out).validate(argv[0])
        print(f"  validate: eulerian-complete={rep.is_eulerian_complete} "
              f"a-trail={rep.is_a_trail} single-stroke={rep.is_single_stroke}")
    except Exception as ex:                       # pragma: no cover
        print(f"  (self-validation skipped: {ex})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
