"""Generate a small, deterministic test solid without distributing mesh assets."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from slice_stl import gen_gear_stl

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    gen_gear_stl(str(args.output), teeth=10, R=12, root=8.4, h=6)
    print(f"Wrote {args.output}")
