"""Build the static browser bundle from this repository's core modules."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent
MODULES = ("origami_slicer.py", "slice_stl.py", "clip.py", "atrail.py", "foldtp.py")


def main():
    destination = ROOT / "web" / "py"
    destination.mkdir(parents=True, exist_ok=True)
    for name in MODULES:
        shutil.copyfile(ROOT / name, destination / name)
    print("Built web/py. Run python3 server.py to preview the slicer.")


if __name__ == "__main__":
    main()
