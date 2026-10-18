# Reproduction scope

`benchmarks/matched_density.py` is adapted from the submission's `gen_eval_table.py`. It derives target length/area from the generated Miura path, tunes the Euler baseline's **base** hex lattice, then measures the refined path. Compare the printed density column as well as length: final Euler density is higher after refinement. The target is computed, not forced to the manuscript's rounded 2.20.

The script runs square and star outlines. `bench.py --parts` compares additional synthetic outlines at fixed spacing; `bench.py --real` procedurally generates and slices test solids. These protocols differ and their raw values should not be mixed.

The CFS-style row intentionally leaves self-crossings unreported. It measures an approximate connected-contour path's turning; it is not the original authors' CFS implementation.

The manuscript's final-density-matched Euler row and other theory/anisotropy tables are not reproduced by these commands. Results can differ from rounded manuscript tables; use this release's actual output and disclose the protocol. No mechanical coupon dataset is included.

Tests verify graph invariants, reject a deliberately crossing fixture, exercise clipping on square and star outlines, and run the slicer and browser-bundle construction without parent-repository dependencies.

## Observed release check

On the square, this checkout produces origami length 194.72 and turning/length 99.8; the Euler baseline produces length 391.94 and turning/length 236.9. These agree with the corresponding rounded manuscript rows. The CFS-style implementation produces length **298.85**, rather than the manuscript table's **207.5**, with turning/length 25.0. Its resulting density is 3.382 rather than the target 2.204. This is a known reproduction discrepancy, not a matched-final-density result. The available script sets contour spacing from the target but does not retune the finished CFS-style path.
