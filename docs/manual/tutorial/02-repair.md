# 2. Repairing the surface

Chapter 1 found that the dragon is neither watertight nor manifold. Before it goes into a grid, `sanitize` repairs
what can be repaired, and the validity predicates say what is left.

<<< @/examples/snippets/tut02_repair-load.f90

<<< @/examples/snippets/tut02_repair-sanitize.f90

- `sanitize` runs the repair passes in order: `remove_degenerate_facets` (zero-area and sliver facets),
  `connect_nearby_vertices` (only when facets have disconnected edges), `analyze`, `remove_duplicate_facets` (the
  same triangle twice, in any winding), then `sanitize_normals` (consistent winding, outward by the sign of the
  volume). The surface is analyzed again after every pass that changes it.
- `fill_holes=.true.` adds `fill_holes` after the orientation pass: every closed loop of boundary edges with at most
  `max_hole_edges` edges (default 200) and a bounding-box diagonal at most `max_hole_size` times the surface's
  (default 0.1) is closed by a minimum-weight patch. Larger loops stay open, so a genuine opening (an inlet, a symmetry
  cut) survives. `refine_holes=.true.` and `fair_holes=.true.` refine the patches and fair them to a smooth surface.
  The default is `fill_holes=.false.`.

<div class="fossil-pair">

![a round hole cut into the repaired dragon, its boundary loop in ember: the facet edges without a neighbour](../../examples/images/tut02-hole.svg)

![the same hole closed by fill_holes with refine and fair: the patch in ember, refined to the density of the loop and faired to a smooth surface](../../examples/images/tut02-patch.svg)

</div>

A hole is a loop of boundary edges, the facet edges without a neighbour (`fcon_edge == 0`); `fill_holes` triangulates
it, `refine` matches the density of the loop and `fair` smooths the patch into the surface around it.
- Each defect class left or removed is reported once on stderr, as a `WARNING: sanitize:` line.

<<< @/examples/snippets/tut02_repair-check.f90

- `get_degenerate_facets_removed` and `get_duplicate_facets_removed` count what the last passes dropped.
- `is_volume()` is true when the surface is manifold, its signed volume positive and its centroid finite: the
  condition for the pseudo-normal sign of chapter 4 to be exact.
- `save_into_file` without `is_ascii` writes ASCII.

## Running it

<<< @/examples/output/tut02.ansi{ansi}

The repair removed 62 duplicate facets and filled 54 holes, but the dragon is still not watertight: 54 edges carry
three or more facets, and 136 boundary edges do not close into a loop, so `fill_holes` cannot patch them. `sanitize`
fixes local defects; it does not rebuild a surface. The rest of the tutorial uses this cleaned but open dragon on
purpose: it shows where the choice of sign algorithm matters (chapters 4 and 5) and which cut cells cannot be trusted
(chapter 10). For a watertight surrogate of a broken surface see [Alpha wrapping](/guide/advanced/alpha-wrap).

::: tip What you learned
`sanitize` (with `fill_holes=.true.` for small holes) repairs local defects and reports them on stderr;
`is_watertight`, `is_manifold` and `is_volume` tell whether the result is a closed solid. Reference:
[Hole filling](/guide/advanced/hole-filling).
:::

Next: [3. Placing the body](./03-manipulation).
