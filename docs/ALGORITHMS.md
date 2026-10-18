# Algorithm overview

The supported flat-foldable degree-4 patterns have even-degree interior vertices. Clipping cuts edges; the boundary patch restores parity. Face checkerboarding and a spanning tree choose local non-crossing transitions that join circuits into an A-trail. `foldtp.py` validates edge coverage, transitions, and single-stroke status independently.

`bench.toolpath_metrics` counts transversal interior segment intersections, not shared endpoints or collinear overlaps. An A-trail revisits graph vertices. This combinatorial crossing-free definition does not model finite extrusion width.

The Euler baseline implements refinement followed by traversal. `cfs.py` is a connected contour-parallel approximation, not the complete Connected Fermat Spirals algorithm; its bridge crossings are not a verified CFS result.

`MICROSPEC.md` specifies the `tp:` FOLD extension and its validation semantics.
