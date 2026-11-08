---
title: The tutorial
---

# The tutorial

The tutorial teaches FOSSIL by following one body, from its STL file into a simulation grid. It is for the author of
a solver, or of a grid generator, who needs geometry from a triangulated surface: distances, signs, wall points, tags,
cut cells. It assumes Fortran 2003 and the basics of immersed-boundary or cut-cell methods; it does not assume
any knowledge of FOSSIL. The [guide](/guide/) explains the algorithms and the [reference](/guide/api/) lists every
method and argument.

## The chapters

The body is `dragon.stl`, a real mesh with real defects: it is not watertight, and repairing it does not make it so.
The tutorial keeps it that way on purpose, so that it shows where an open surface changes the answer and which
queries tolerate it. Chapters 1 to 3 load, repair and place the dragon in the domain [-0.6, 0.6]³; chapters 4 to 6
compute its signed distance, point by point, in batches and as a field; chapter 7 goes to 2-D with a section of it;
chapters 8 to 11 answer what an immersed-boundary or cut-cell solver asks, for a fixed and then a moving body.
Chapters 12 to 15 widen the scope: processing and analysing a surface, structured blocks and overset grids, contact
between bodies, and body-fitted grids; they build their own geometry and need no file from the earlier chapters.
Chapter 16 moves the batch to the GPU.

| Chapter | You learn |
|---|---|
| [1. A first surface](./tutorial/01-first-surface) | `load_from_file`, what `analyze` computes, accessors, `facet_at` and `fcon_edge`, `is_watertight`, `is_manifold`, `feature_edges` (boundary, non-manifold), shells and bodies, `save_into_file`, `set_patch` and multi-solid ASCII files |
| [2. Repairing the surface](./tutorial/02-repair) | `sanitize` and its passes one by one, `fill_holes` (`max_hole_size`, refine, fair), the stderr report, `is_volume`, what is left, `deviation` |
| [3. Placing the body](./tutorial/03-manipulation) | `translate`, `rotate`, `resize` about the world origin, `analyze` after a transformation, `clip` with a remainder, `mirror`, `merge_solids` |
| [4. Signed distance](./tutorial/04-distance) | `distance`, `is_point_inside`, the four sign algorithms and where the pseudo-normal fails, `compute_distance` (closest facet, feature, point, normal) |
| [5. Batched queries](./tutorial/05-batched) | `distance_many`, SAH BVH and octree (what each splits), pseudo-normal and fast winding number on a grid |
| [6. A distance field](./tutorial/06-distance-field) | `distance_field_object` on a stretched grid: `compute`, and `compute_band`, `compute_signs`, `sweep` one by one; checking the field |
| [7. Curves in 2-D](./tutorial/07-curves) | `slice_surface` at many levels, `curve_object`, `self_intersections`, 2-D distance and batch, 2-D field |
| [8. Immersed-boundary queries](./tutorial/08-ib-queries) | ghost points, `query_many`: wall point, feature, normal, body, ambiguity; `image_points`, `probe_many`; `thickness`, `thin_regions` |
| [9. AMR tagging](./tutorial/09-amr-tagging) | `box_classify_many` (cut, near, inside, outside) and an octree refined level by level, `tag_block` on a block and on a section, feature-edge and thin-wall refinement (`build_features`, `build_thin`) |
| [10. Cut cells](./tutorial/10-cut-cells) | exact `line_crossings`, `curve%cut_cells` and `cut_cells_3d` (fractions, apertures, fluid pieces), then the solver support: piece links, SRD neighbourhoods, merge targets, moments, quadrature, coarsening |
| [11. Moving bodies](./tutorial/11-moving-bodies) | poses, `pose=` queries and `move`, `state_changes`, `swept_cells_2d` / `swept_cells_3d` for a translation and a rotation, the inertia of a closed surrogate, the gap between two bodies |
| [12. Geometry processing and analysis](./tutorial/12-geometry) | `boolean`, `set_patch`, `isotropic_remesh`, `deviation`, `feature_edges`, `thickness`, `thin_regions`, `mass_properties`, on a bracket built step by step |
| [13. Structured blocks and overset grids](./tutorial/13-blocks-overset) | an O-grid's topology (`detect_topology_2d`), `validity`, `quality`, `locate_many`, `stencils`; its overset assembly (`wall_cutter_2d`, `hole_flags_2d`, `donor_candidates_2d`, `classify_2d`, fringes, orphans, `save_connectivity`) and a moving near-body grid |
| [14. Contact: distance and collision](./tutorial/14-contact) | why sampling misses contacts; `distance_to` with poses, `time_of_impact` (plain and `directional=`), `first_impact`, `proximity_pairs_of`; the 2-D twins |
| [15. Body-fitted grids](./tutorial/15-body-fitted) | `projection_target`, free versus constrained projection (`classify_points`, `project_points`), `project_block_face` without collapsing cells, `geometric_spacing` from a y₁⁺ target, `hyperbolic_grid_2d`, `block_report_2d`, `save_block_report` |
| [16. On the GPU](./tutorial/16-gpu) | `enter_device` and the `*_device` batches (distances, full and wall queries, box classification, curves), equal to the host bit for bit; one batch versus one slice per call; what is refused; `exit_device`; the nvfortran build |

The chapters pass files to each other, so run them in order:

```mermaid
flowchart LR
  dragon[(dragon.stl)] --> c1["1: load"]
  dragon --> c2["2: repair"] --> clean[(dragon-clean.stl)]
  clean --> c3["3: place"] --> placed[(dragon-placed.stl)]
  placed --> c4to6["4-6, 8, 9, 16"]
  placed --> c7["7: slice"] --> section[(dragon-section.dat)]
  section --> c10["10: cut cells"]
  placed --> c10
  section --> c11["11: moving"]
  placed --> c11
```

## How the examples are made

Every snippet in these pages is part of a real program in
[`docs/examples/src`](https://github.com/szaghi/FOSSIL/tree/master/docs/examples/src), one per chapter
(`tutNN_<name>.f90`), and every output shown is what that program printed: captured, not typed. Comment markers in
the sources (`!region NAME` ... `!endregion NAME`) cut the programs into the parts the pages show. The programs print
with explicit formats only, so their output does not depend on the compiler.

`bash scripts/docs_examples.sh` rebuilds the library (`fobis build --mode static-gnu`), builds every program, runs
them in order in one scratch directory holding `dragon.stl`, `cube-01.stl` and `naca0012-ascii.stl`, and regenerates
the snippets and outputs these pages include. To build one chapter by hand, with the library built in `static/`:

```bash
gfortran -I static/mod docs/examples/src/tut01_first_surface.f90 static/libfossil.a -o tut01
cp src/tests/dragon.stl . && ./tut01
```
