---
title: Upgrading to 3.0
---

# Upgrading to 3.0

FOSSIL 3.0 adds roadmap #68 (cut-cell solver support, curvilinear blocks, overset assembly, continuous collision
detection, body-fitted grids), the documentation rebuild (#77, #85), and the fixes the documentation work exposed (#90,
#91). No public name was removed or renamed: code that compiled against 2.2 compiles against 3.0. Some routines now
**behave differently**, which is why the major version changes. This page lists every change in behaviour since 2.2.0,
with its opt-out where one exists, then the public names that are new.

The [changelog](./changelog) lists every commit.

## Changes in behaviour

| Routine | 2.2 | 3.0 | Opt-out | Issue |
|---|---|---|---|---|
| `intersect_ray_first`, `_all`, `_any` | Moller-Trumbore on the closed range: a ray through a shared edge or vertex hit every facet there (12 hits for 2 crossings through cube-01's face centres); the octree could lose hits | one hit per crossing, by exact tests with a symbolic tie-break; hits at t = 0 are kept; equal t sort by facet id, so `_first` equals `_all(1)` | `_all(exclude=)` skips a facet (the ray's own, for instance) | [#82](https://github.com/szaghi/FOSSIL/issues/82) |
| `isotropic_remesh` | almost never collapsed an edge: it did not coarsen | coarsens (link condition), tangential relaxation, and every collapse, flip and move is refused when its facets leave the input by more than `REMESH_ERROR_TOL` = 0.07 L. Facet counts and deviations differ from 2.2 | none (`preserve_features` unchanged) | [#79](https://github.com/szaghi/FOSSIL/issues/79), [#91](https://github.com/szaghi/FOSSIL/issues/91) |
| `decimate` | collapses could leave zero-area facets and non-manifold edges with `DEC_STATUS_OK`; `target_facets <= 0` returned an empty surface; patch rims never coarsened | link condition and zero-area refusal; `target_facets <= 0` is `DEC_STATUS_BAD_INPUT`; patch rims coarsen along themselves onto input rim vertices (corners fixed); an unreachable target gives `DEC_STATUS_NO_PROGRESS` with a valid surface | `rim_collapses=.false.` locks the rims; `max_error=` (new, opt-in) is a cost ceiling | [#81](https://github.com/szaghi/FOSSIL/issues/81), [#91](https://github.com/szaghi/FOSSIL/issues/91) |
| `resolve_self_intersections` | returned the arrangement unselected (volume counted twice, not closed) with `BOOL_STATUS_OK` | the self-union (patches of winding number 0 in front, Zhou et al.); `BOOL_STATUS_NOT_CLOSED` when the selection is open (the surface is left unchanged); triple crossings get a three-plane implicit point | none | [#78](https://github.com/szaghi/FOSSIL/issues/78) |
| `jet_curvature`, `query_many(curvatures=)` | a rank-deficient fit silently kept the cotangent value (H = 7.4e9 on near-coincident rings) | near-coincident ring vertices merged, the ring grown; a fit that stays rank-deficient gives NaN and `CURV_STATUS_FIT_FAILED` | none | [#91](https://github.com/szaghi/FOSSIL/issues/91) |
| `extract_isosurface` | node values missing the level by round-off produced zero-area slivers | such nodes are snapped to the level (within 1e-12 of their edge crossings): no slivers; other cells bit-identical | none | [#87](https://github.com/szaghi/FOSSIL/issues/87) |
| `fill_holes`, `sanitize` | `open_edges` and the "boundary edges" warning counted non-manifold edges too | `open_edges` counts true boundary edges; `nonmanifold_edges=` reports the others; the filled surfaces are unchanged | `manifold_only=.true.` (new, opt-in; `sanitize(manifold_holes=)`) | [#88](https://github.com/szaghi/FOSSIL/issues/88), [#91](https://github.com/szaghi/FOSSIL/issues/91) |
| `save_into_file` | a binary save over a larger file kept its tail | the file is replaced | none | [#80](https://github.com/szaghi/FOSSIL/issues/80) |
| `patch_name` (binary STL) | the 80-byte header with its NUL padding | the header cut at the first NUL; `get_header` keeps the raw bytes | none | [#83](https://github.com/szaghi/FOSSIL/issues/83) |
| `swept_cells_2d`, `swept_cells_3d` | `exact = .false.` on some rotations whose residual was round-off | `exact = .true.` there; the values are unchanged | none | [#84](https://github.com/szaghi/FOSSIL/issues/84), [#86](https://github.com/szaghi/FOSSIL/issues/86) |
| `proximity_pairs_of` | grew its result by self-referencing reshapes | allocates once; results bit-identical | none | [#76](https://github.com/szaghi/FOSSIL/issues/76) |

What to check in your code:

- **Remesh and decimate outputs** are different meshes: re-measure anything that depends on their facet counts or
  deviations. The [isotropic remeshing](./advanced/isotropic-remesh) and [decimation](./advanced/decimate) pages
  give the new numbers.
- **Code that counted ray hits** to decide inside/outside, or merged equal-t hits by hand, can drop that workaround.
- **Code that tested `decimate`'s status** gets `DEC_STATUS_BAD_INPUT` for a non-positive target and
  `DEC_STATUS_NO_PROGRESS` where 2.2 returned a broken surface with OK.
- **Curvature consumers** should accept NaN with `CURV_STATUS_FIT_FAILED` where 2.2 returned a meaningless number.
- **`resolve_self_intersections` callers** should handle `BOOL_STATUS_NOT_CLOSED` (open inputs).

## New public names

All exported by `use fossil` unless noted; no name was removed or renamed.

- **Cut-cell solver support** ([page](./advanced/cut-solver)): `cut_links_t`, `cut_links_2d`, `srd_neighbourhoods`,
  `srd_t`, `SRD_*`, `frd_neighbourhoods`, `frd_t`, `merge_targets`, `merge_t`, `MERGE_*`, `cut_coarsen_2d`,
  `cut_coarsen_3d`, `MOMENTS_MAX`, `mono_count`, `mono_exponents`; modules `fossil_cut_graph`, `fossil_cut_moments`,
  `fossil_cut_coarsen`; new arguments of `cut_cells_3d` (`links=`, `edges=`, `moments_order=`, `quadrature=`) and
  `curve%cut_cells`.
- **Curvilinear blocks** ([page](./advanced/curvilinear-blocks)): `curvilinear_block_t`, `curvilinear_block_2d_t`,
  `topo_record_t`, `TOPO_*`, `record_map`, `detect_topology(_2d)`, `CELL_*`, `LOC_*`, `locate_in_blocks(_2d)`,
  `STENCIL_*`, `block_real_t`, `block_logical_t`, `BLOCK_STATUS_*`; modules `fossil_curvilinear_block`,
  `fossil_block_geometry`, `fossil_block_locate`, `fossil_block_stencil`.
- **Overset assembly** ([page](./advanced/overset)): `wall_cutter(_2d)` (strict trimming by default, `strict=.false.` for
  the centroid rule), `hole_flags(_2d)`, `hole_map(_2d)`,
  `xray_flags(_2d)`, `donor_candidates(_2d)`, `donor_set_t`, `overlap_matrix(_2d)`, `PAIR_*`, `classify(_2d)`,
  `SELECT_*`, `overset_mask_t`, `ORPHAN_*`, `fringe_levels(_2d)`, `FRINGE_*`, `overlap_width(_2d)`,
  `wall_projection(_2d)`, `exposed_points`, `cell_centre_block(_2d)`, `connectivity_t`, `build_connectivity`,
  `save_connectivity`; module `fossil_overset`.
- **Collision** ([page](./advanced/collision)): `surface%time_of_impact` (adaptive default; `directional=` forces the
  plain or the directional step), `curve%time_of_impact`, `first_impact(_2d)`,
  `motion_bound_3d/_2d`, `ccd_advance`, `ccd_evaluator_t`, `CCD_*` (statuses, `CCD_MAX_ITERATIONS`,
  `CCD_STALL_STEPS`, `CCD_STALL_FRACTION`), `safe_step`; module `fossil_ccd`.
- **Body-fitted grids** ([page](./advanced/body-fitted)): `surface%projection_target`, `curve%projection_target`,
  `projection_target_t`, `PROJ_*`, `project_points(_2d)`, `classify_points(_2d)`, `project_block_face`,
  `project_block_side`, `block_report(_2d)`, `block_report_t`, `save_block_report`, `utau_flat_plate`,
  `hyperbolic_grid_2d`, `geometric_spacing`; modules `fossil_projection`, `fossil_block_report`, `fossil_hyperbolic`.
- **Self-intersections and booleans**: `SELFINT_*` kinds, `BOOL_STATUS_NOT_CLOSED`; the implicit point kind `IPK_TPI`
  (module `fossil_implicit_points`).
- **Queries and trees**: `surface%tree_boxes`, `curve%tree_boxes`; `intersect_ray_all(exclude=)`.
- **Processing**: `decimate(max_error=, rim_collapses=)`, `REMESH_ERROR_TOL` (module `fossil_remesh`),
  `CURV_STATUS_FIT_FAILED`, `fill_holes(nonmanifold_edges=, manifold_only=)`, `sanitize(manifold_holes=)`.

## The documentation

The site was rebuilt around runnable examples: 16 tutorial chapters, 85 cookbook recipes (one for every entry point of
the [feature map](./features)), and a generated figure for every explaining section of the reference pages, 373
figures in all, every one computed by a program that runs when the docs are built.
