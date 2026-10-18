#!/usr/bin/env python3
"""
CLI:
    python3 foldtp.py FILE.fold [FILE2.fold ...]      # validate + report
    python3 foldtp.py examples/*.fold

Library:
    from foldtp import FoldTP
    f = FoldTP.load("examples/figure8_tour.fold")
    report = f.validate()
    print(report.ok, report.is_a_trail, report.is_single_stroke)
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

K_PATHS_V = "tp:paths_vertices"
K_PATHS_E = "tp:paths_edges"
K_PATHS_CLOSED = "tp:paths_closed"
K_PATHS_TRAVEL = "tp:paths_travel"
K_PATHS_ROLE = "tp:paths_role"
K_PATHS_SRC = "tp:paths_source"
K_VTRANS = "tp:vertices_transition"

CREASE_ASSIGN = frozenset("MVFU")


@dataclass
class Report:
    path: str
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    has_toolpath: bool = False
    is_eulerian_complete: bool = False
    is_a_trail: bool = False
    is_single_stroke: bool = False

    @property
    def ok(self) -> bool:
        return not self.errors


class FoldTP:
    """A parsed FOLD file with optional ``tp:`` toolpath data."""

    def __init__(self, raw: dict):
        self.raw = raw
        self.coords = raw.get("vertices_coords", [])
        self.edges = raw.get("edges_vertices", [])
        self.assignment = raw.get("edges_assignment")
        self.n_v = len(self.coords)
        self.n_e = len(self.edges)
        self._vinc: list[list[int]] = [[] for _ in range(self.n_v)]
        for ei, (u, v) in enumerate(self.edges):
            if 0 <= u < self.n_v:
                self._vinc[u].append(ei)
            if 0 <= v < self.n_v:
                self._vinc[v].append(ei)

    @classmethod
    def load(cls, path: str | Path) -> "FoldTP":
        with open(path) as fh:
            return cls(json.load(fh))

    def other(self, edge: int, vertex: int) -> int:
        """The endpoint of ``edge`` that is not ``vertex``."""
        u, v = self.edges[edge]
        return v if u == vertex else u

    def rotation(self, vertex: int) -> list[int]:
        """CCW-ordered incident edge IDs at ``vertex``.

        Uses ``vertices_edges`` if present (spec mandates ccw for orientable
        manifolds); otherwise computes the order from 2D coordinates by angle.
        """
        ve = self.raw.get("vertices_edges")
        if ve is not None and vertex < len(ve):
            return list(ve[vertex])
        x0, y0 = self.coords[vertex][0], self.coords[vertex][1]

        def ang(ei: int) -> float:
            o = self.other(ei, vertex)
            return math.atan2(self.coords[o][1] - y0, self.coords[o][0] - x0)

        return sorted(self._vinc[vertex], key=ang)

    def validate(self, src_path: str = "<dict>") -> Report:
        r = Report(path=src_path)
        self._check_structure(r)
        if not (K_PATHS_V in self.raw):
            r.notes.append("no tp: toolpath present — standard FOLD only.")
            self._report_crease_pattern(r)
            return r
        r.has_toolpath = True
        paths = self._check_paths(r)
        if r.errors:
            return r
        printed = self._collect_printed(paths, r)
        self._check_eulerian_complete(printed, r)
        self._check_transition_system(paths, r)
        self._report_single_stroke(paths, r)
        return r

    def _check_structure(self, r: Report) -> None:
        if not self.coords:
            r.errors.append("missing/empty vertices_coords")
        if not self.edges:
            r.errors.append("missing/empty edges_vertices")
        for ei, e in enumerate(self.edges):
            if len(e) != 2 or any(not (0 <= x < self.n_v) for x in e):
                r.errors.append(f"edge {ei}={e} has bad vertex index")

    def _check_paths(self, r: Report) -> list[dict]:
        pv = self.raw.get(K_PATHS_V, [])
        pe = self.raw.get(K_PATHS_E, [])
        pc = self.raw.get(K_PATHS_CLOSED, [])
        pt = self.raw.get(K_PATHS_TRAVEL)
        k = len(pv)
        if not (len(pe) == len(pc) == k):
            r.errors.append(
                f"tp:paths_* length mismatch: vertices={len(pv)} edges={len(pe)} "
                f"closed={len(pc)}")
            return []
        if pt is not None and len(pt) != k:
            r.errors.append(f"tp:paths_travel length {len(pt)} != {k}")
            return []
        out = []
        for p in range(k):
            verts, edgs, closed = list(pv[p]), list(pe[p]), bool(pc[p])
            m = len(verts)
            s = m if closed else m - 1
            travel = list(pt[p]) if pt is not None else [False] * s
            if len(edgs) != s:
                r.errors.append(
                    f"path {p}: edges len {len(edgs)} != expected segments {s} "
                    f"(closed={closed}, m={m})")
                continue
            if len(travel) != s:
                r.errors.append(f"path {p}: travel len {len(travel)} != {s}")
                continue
            for i in range(s):
                a = verts[i]
                b = verts[(i + 1) % m] if closed else verts[i + 1]
                ei = edgs[i]
                if travel[i]:
                    if ei is not None:
                        r.warnings.append(
                            f"path {p} seg {i}: travel segment should have null edge")
                    continue
                if ei is None or not (0 <= ei < self.n_e):
                    r.errors.append(f"path {p} seg {i}: bad printed edge id {ei}")
                    continue
                if {a, b} != set(self.edges[ei]):
                    r.errors.append(
                        f"path {p} seg {i}: edge {ei}={self.edges[ei]} does not "
                        f"connect vertices {a}->{b}")
            out.append(dict(verts=verts, edges=edgs, closed=closed, travel=travel))
        return out

    def _collect_printed(self, paths: list[dict], r: Report) -> set[int]:
        seen: dict[int, int] = {}
        for p in paths:
            for ei, tv in zip(p["edges"], p["travel"]):
                if tv or ei is None:
                    continue
                seen[ei] = seen.get(ei, 0) + 1
        dupes = sorted(e for e, c in seen.items() if c > 1)
        if dupes:
            r.errors.append(f"edges printed more than once (not edge-disjoint): {dupes}")
        return set(seen)

    def _check_eulerian_complete(self, printed: set[int], r: Report) -> None:
        if self.assignment is None:
            creases = set(range(self.n_e))  
        else:
            creases = {i for i, a in enumerate(self.assignment) if a in CREASE_ASSIGN}
        missing = sorted(creases - printed)
        extra = sorted(printed - creases)
        r.is_eulerian_complete = not missing
        if missing:
            r.warnings.append(f"crease edges not covered (not Eulerian-complete): {missing}")
        if extra:
            r.notes.append(f"printed non-crease (boundary) edges: {extra}")

    def _implied_transitions(self, paths: list[dict]) -> dict[int, list[frozenset]]:
        """For each vertex, the (incoming,outgoing) printed-edge pairs from the walks."""
        trans: dict[int, list[frozenset]] = {}
        for p in paths:
            verts, edgs, closed, travel = (
                p["verts"], p["edges"], p["closed"], p["travel"])
            m = len(verts)
            s = m if closed else m - 1
            for idx, v in enumerate(verts):
                out_i = idx if idx < s else None
                in_i = (idx - 1) if idx > 0 else (s - 1 if closed else None)
                e_in = edgs[in_i] if (in_i is not None and not travel[in_i]) else None
                e_out = edgs[out_i] if (out_i is not None and not travel[out_i]) else None
                if e_in is not None and e_out is not None:
                    trans.setdefault(v, []).append(frozenset((e_in, e_out)))
        return trans

    def _check_transition_system(self, paths: list[dict], r: Report) -> None:
        implied = self._implied_transitions(paths)
        declared_raw = self.raw.get(K_VTRANS)
        declared: dict[int, list[frozenset]] = {}
        if declared_raw is not None:
            for v, pairs in enumerate(declared_raw):
                if pairs:
                    declared[v] = [frozenset(pr) for pr in pairs]
        if declared_raw is None:
            r.notes.append("no tp:vertices_transition given — using path-implied transitions.")
            trans = implied
        else:
            trans = declared
            for v in set(implied) | set(declared):
                if sorted(map(sorted, implied.get(v, []))) != \
                   sorted(map(sorted, declared.get(v, []))):
                    r.warnings.append(
                        f"vertex {v}: declared transition {[sorted(x) for x in declared.get(v, [])]} "
                        f"!= path-implied {[sorted(x) for x in implied.get(v, [])]}")
        all_non_crossing = True
        for v, pairs in trans.items():
            inc = set(self._vinc[v])
            matched = [e for pr in pairs for e in pr]
            
            if len(matched) != len(set(matched)):
                r.errors.append(f"vertex {v}: transition pairs reuse an edge: {[sorted(x) for x in pairs]}")
            if not set(matched) <= inc:
                r.errors.append(f"vertex {v}: transition references non-incident edge")
            if not self._non_crossing(v, pairs):
                all_non_crossing = False
                r.notes.append(f"vertex {v}: transition CROSSES under rotation system.")
        r.is_a_trail = all_non_crossing and r.is_eulerian_complete and not r.errors

    def _non_crossing(self, vertex: int, pairs: list[frozenset]) -> bool:
        """True iff the matching's chords don't interleave in the ccw edge order."""
        order = self.rotation(vertex)
        pos = {e: i for i, e in enumerate(order)}
        chords = []
        for pr in pairs:
            es = [e for e in pr if e in pos]
            if len(es) != 2:
                continue
            i, j = sorted((pos[es[0]], pos[es[1]]))
            chords.append((i, j))
        for a in range(len(chords)):
            for b in range(a + 1, len(chords)):
                i1, j1 = chords[a]
                i2, j2 = chords[b]
                inside2 = (i1 < i2 < j1) + (i1 < j2 < j1)
                if inside2 == 1:
                    return False
        return True

    def _report_single_stroke(self, paths: list[dict], r: Report) -> None:
        single = (len(paths) == 1 and paths[0]["closed"]
                  and not any(paths[0]["travel"]))
        r.is_single_stroke = single

    def _report_crease_pattern(self, r: Report) -> None:
        """Light crease-pattern report for standard FOLD files (Maekawa check)."""
        if self.assignment is None:
            return
        counts: dict[str, int] = {}
        for a in self.assignment:
            counts[a] = counts.get(a, 0) + 1
        r.notes.append(f"edge assignments: {counts}")
        for v in range(self.n_v):
            inc = self._vinc[v]
            if any(self.assignment[e] == "B" for e in inc):
                continue  
            if len(inc) < 3:
                continue
            m = sum(1 for e in inc if self.assignment[e] == "M")
            val = sum(1 for e in inc if self.assignment[e] == "V")
            ok = abs(m - val) == 2
            r.notes.append(
                f"  interior vertex {v}: deg {len(inc)}, M={m} V={val} -> "
                f"Maekawa {'OK' if ok else 'VIOLATED'}")


def _fmt(r: Report) -> str:
    lines = [f"== {r.path} =="]
    if r.has_toolpath:
        verdict = "VALID" if r.ok else "INVALID"
        lines.append(
            f"  toolpath: {verdict} | eulerian-complete={r.is_eulerian_complete} "
            f"| crossing-free(A-trail)={r.is_a_trail} | single-stroke={r.is_single_stroke}")
    for e in r.errors:
        lines.append(f"  ERROR:   {e}")
    for w in r.warnings:
        lines.append(f"  WARN:    {w}")
    for n in r.notes:
        lines.append(f"  note:    {n}")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    rc = 0
    for arg in argv[1:]:
        try:
            f = FoldTP.load(arg)
            r = f.validate(arg)
        except Exception as exc: 
            print(f"== {arg} ==\n  ERROR:   could not parse: {exc}")
            rc = 1
            continue
        print(_fmt(r))
        if not r.ok:
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
