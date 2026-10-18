# Origami toolpaths

Research implementation accompanying **Flat-Foldable Origami Tessellations as Native Continuous Toolpaths for Fused Deposition Modeling**, submitted to ICRA 2027.

Generate continuous, crossing-free infill paths from structured degree-4 origami crease patterns. This standalone release brings together the working graph algorithms, validation, benchmarks, and experimental STL/OBJ slicer. No submodules or parent checkout are required.

## Requirements

Python 3.10 or newer. All Python code uses the standard library; no pip packages are required. The optional browser slicer downloads Pyodide 0.27.2 from jsDelivr on first load. Slicing itself happens locally in the browser.

## Quickstart

Run these commands from this repository's root:

```sh
python3 -m unittest discover -s tests -v
python3 foldtp.py examples/square_tour.fold examples/figure8_tour.fold
mkdir -p outputs
python3 atrail.py --patch examples/figure8_tour.fold outputs/toolpath.fold
python3 foldtp.py outputs/toolpath.fold
```

Open `toolpath_viewer.html` in a browser to inspect FOLD toolpaths.

## Benchmarks

```sh
python3 benchmarks/matched_density.py
python3 bench.py --parts
python3 bench.py --real
```

The first command runs the paper's available matched-*base*-density script on a square and non-convex star. It does **not** reproduce every table in the manuscript: the final-density-matched Euler row, anisotropy analysis, and full theory experiments are not included in this release. See [reproduction notes](docs/REPRODUCIBILITY.md).

## Browser slicer

```sh
python3 build_web.py
python3 server.py
```

Open http://localhost:8765 and select an STL or OBJ. The app previews per-layer toolpaths and exports experimental G-code. To host the slicer, build once and publish the contents of `web/` to a static host. `web/py/` is generated, not a second source tree.

## Command-line slicer

```sh
python3 examples/make_gear.py outputs/gear.stl
python3 origami_slicer.py outputs/gear.stl -o outputs/gear.gcode --fit-mm 24 --layer-height 0.4 --spacing 1.5 --max-layers 3
```

Use `python3 origami_slicer.py --help` for options. Output is written only to the requested path. `--max-layers` caps the output; the example deliberately produces only three layers.

### Prototype limits

- The mesh slicer selects the largest contour on each layer and simplifies it. Holes and disconnected components are not preserved: use simple, single-component meshes without internal voids.
- The single-stroke result concerns successful **infill**, not the entire layer or part. Perimeter walls are separate strokes; failed infill generation can leave walls only. Inspect `layer_strokes` and `failed_layers` plus the preview.
- G-code uses a fixed bed center at (110, 110) mm and generic extrusion/retraction settings. It has no printer-specific homing, heating, or startup sequence. Review it in a printer-aware slicer and adapt it to your machine before fabrication.
- The crossing-free construction targets the supported degree-4 families. Even degree alone does not guarantee an A-trail in arbitrary higher-degree patterns.
- Physical surface quality, stiffness, and comparative print times are not quantitatively validated here.

## Layout

| File | Purpose |
| --- | --- |
| `atrail.py`, `clip.py` | Boundary repair and checkerboard spanning-tree traversal |
| `foldtp.py`, `MICROSPEC.md` | FOLD toolpath representation and validator |
| `bench.py` | Miura generation, infill baselines, structural metrics |
| `euler_transform.py` | Euler-transformation baseline |
| `cfs.py` | Connected-contour approximation used for smoothness comparison |
| `slice_stl.py`, `origami_slicer.py` | Mesh slicing and experimental G-code output |
| `benchmarks/` | Matched-base-density evaluation |
| `examples/`, `tests/` | Small fixtures, mesh generator, regression checks |
| `web/`, `build_web.py`, `server.py` | Browser slicer and local preview tooling |


