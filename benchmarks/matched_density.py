"""Submission benchmark with nominal density matching.

The Euler base lattice is tuned before refinement; CFS contour spacing is derived
from the target but its final density is not matched. See docs/REPRODUCIBILITY.md.
"""
import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bench, clip, atrail, cfs
import euler_transform as et

def origami(outline, bbox, s):
    x0,y0,x1,y1 = bbox; mg = 3*s
    nx = max(2, round((x1-x0+2*mg)/s)); ny = max(2, round((y1-y0+2*mg)/s))
    M = bench.miura(nx, ny, a=s, b=s, c=s*0.36)
    for v in M["vertices_coords"]: v[0]+=x0-mg; v[1]+=y0-mg
    clipped,_ = clip.clip_to_polygon(M, outline)
    out, problems = atrail.ATrailGenerator(clipped).generate(patch=True)
    assert not problems, problems
    return bench.toolpath_metrics(out)

def hex_base_density(outline, bbox, size):
    x0,y0,x1,y1 = bbox
    polys = et.hex_mesh(x0,y0,x1,y1,size,outline=outline)
    coords, faces, ef = et.build_mesh(polys)
    faces, ef = et.fix_assumption23(faces, ef)
    L = sum(et._dist(coords[a],coords[b]) for a,b in ef)
    return L, polys

def run(outline, bbox, s0):
    A = bench.poly_area(outline)
    rows = {}
    om = origami(outline, bbox, s0)
    lam = om["printed_length"]/A
    rows["origami"] = om
    rows["rect-sep"] = bench.toolpath_metrics(bench.grid_infill(bbox, 1/lam, directions=("h",), polygon=outline))
    rows["boustro"]  = bench.toolpath_metrics(bench.boustrophedon_infill(bbox, 1/lam, polygon=outline))
    rows["grid"]     = bench.toolpath_metrics(bench.grid_infill(bbox, 2/lam, directions=("h","v"), polygon=outline))
    # CFS at matched density
    cp = cfs.cfs_infill(outline, [], 1/lam)
    cm = cfs.metrics(cp)
    rows["cfs"] = {"strokes":1,"travel_moves":0,"self_crossings":"-","printed_length":round(cm["length"],2),
                   "turning_deg":round(cm["turning_deg"],1)}
    # ET: bisect hex size so BASE mesh density == lam
    lo, hi = 0.4*s0, 4*s0
    for _ in range(22):
        mid = (lo+hi)/2
        L,_p = hex_base_density(outline, bbox, mid)
        if L/A > lam: lo = mid
        else: hi = mid
    size = (lo+hi)/2
    Lb, polys = hex_base_density(outline, bbox, size)
    gout, gproblems, gbase_e = et.gupta_toolpath(polys)
    assert not gproblems, gproblems
    gm = bench.toolpath_metrics(gout)
    gm["_base_density"] = Lb/A; gm["_edge_ratio"] = len(gout["edges_vertices"])/gbase_e
    gm["_len_ratio"] = gm["printed_length"]/Lb; gm["_hex_size"]=size
    rows["hex+ET"] = gm
    print(f"area={A:.2f}  lambda(origami)={lam:.3f} len/area")
    for k,m in rows.items():
        t = m.get("turning_deg"); l = m.get("printed_length")
        tpm = f"{t/l:.1f}" if isinstance(t,(int,float)) and isinstance(l,(int,float)) else "-"
        print(f"  {k:9s} strokes={m['strokes']:>3} travel={m['travel_moves']:>3} "
              f"cross={m['self_crossings']:>4} len={m['printed_length']:>7} dens={float(m['printed_length'])/A:.3f} "
              f"turn={t} t/mm={tpm}", {x:round(v,3) for x,v in m.items() if x.startswith('_')})
    return rows

RECT = ([[.3,.3],[9.7,.3],[9.7,9.7],[.3,9.7]], (0,0,10,10))
STAR = ([[5+4.6*math.cos(math.radians(90+k*36))*(1 if k%2==0 else 0.42),
          5+4.6*math.sin(math.radians(90+k*36))*(1 if k%2==0 else 0.42)] for k in range(10)], (0,0,10,10))
def main():
    print("== rectangle 10x10, s0=1.0 =="); run(*RECT, 1.0)
    print("\n== star (non-convex), s0=1.0 =="); run(*STAR, 1.0)


if __name__ == "__main__":
    main()
