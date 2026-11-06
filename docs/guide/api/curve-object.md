---
title: curve_object
---

# `curve_object`

Closed (or open) polylines in an axis-aligned plane: the 2-D geometry of
FOSSIL, for quadtree solvers. See
[2-D geometry](/guide/advanced/2d-geometry) for the conventions and the
algorithms.

```fortran
program curve_tour
use fossil, only : curve_object, CURVE_VALID, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
implicit none
type(curve_object) :: c
real(R8P)          :: d, cp(2)
integer(I4P)       :: seg, vtx

! unit square, counter-clockwise (solid on the left)
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
if (c%check_validity() /= CURVE_VALID) error stop
call c%compute_distance(p=[2._R8P, 2._R8P], distance=d, closest_point=cp, closest_segment=seg, closest_vertex=vtx)
print *, d, cp, vtx                                        ! sqrt(2), (1,1), vertex 3 (seg: 2 or 3, both meet there)
print *, c%winding_number([0.5_R8P, 0.5_R8P])              ! 1
print *, c%is_point_inside([0.5_R8P, 0.5_R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER)
endprogram curve_tour
```

## Building and I/O

| Method | Purpose |
|---|---|
| `adopt_loops(u, v, loop_sizes, normal_axis, closed, status)` | Build from vertex arrays: `loop_sizes(l)` vertices per loop, closed implicitly unless `closed(l)` is false. `STATUS_INVALID_INPUT` on inconsistent sizes, a bad axis, or non-finite coordinates. |
| `load_from_file(file_name, normal_axis, status)` | ASCII loop format: `u v` per line, blank lines between loops, `#` comments; `# open` before a loop marks it open; a repeated closing vertex is dropped. |
| `save_into_file(file_name, status)` | Write the ASCII loop format, bit-exact on reload. |
| `destroy` | Reset to an empty curve. |

## Metrics and validity

| Method | Purpose |
|---|---|
| `get_vertices_number`, `get_loops_number`, `get_loop_range(l, first, last)`, `get_vertex(i)`, `get_next(i)` | Topology and coordinates; `get_next(i)` is the end vertex of segment i (0 at the end of an open loop). |
| `get_area`, `get_loop_area(l)`, `get_perimeter`, `get_centroid`, `get_bmin`, `get_bmax`, `get_normal_axis` | Metrics; areas are signed (solid area for a valid curve; open loops have none). |
| `is_loop_closed(l)` | False for an open loop. |
| `check_validity(loop, check_self_intersections)` | First defect: `CURVE_OPEN_LOOP`, `CURVE_TOO_FEW_VERTICES`, `CURVE_ZERO_LENGTH_SEGMENT`, `CURVE_ZERO_AREA_LOOP`, `CURVE_WRONG_ORIENTATION` (orientation not matching the exact nesting depth), or `CURVE_VALID`. With `check_self_intersections=.true.` (default `.false.`), also `CURVE_SELF_INTERSECTING` when two segments cross, touch or overlap. |
| `nesting_depth(l)` | Number of closed loops containing loop l (exact). |
| `mass_properties(mass, center, second_moments, polar_moment, principal_moments, principal_angle, density, about, area, status)` | Per unit depth: mass, centroid, second moments [I_uu, I_vv, I_uv] about the centroid (or `about`), polar moment, principal values (ascending) and the angle of the long axis; exact Green's-theorem sums. Closed loops and a positive area required, else `STATUS_INVALID_INPUT`. See [mass properties](/guide/advanced/mass-properties#2-d-curves) (issue #59). |
| `deviation(other, lower, upper, mean, rms, per_vertex, tolerance, samples_per_length, max_segments, samples, created, device, status)`, `hausdorff(other, lower, upper, ...)` | One-sided deviation h(self, other) and the symmetric Hausdorff distance, bracketed lower ≤ h ≤ upper (certified by halving segments); mean, rms, per-vertex distances; open loops allowed; `device=.true.` through `query_many_device`. See [surface deviation](/guide/advanced/surface-deviation#2-d-curves) (issue #60). |
| `feature_corners(corners, angle, status)`; `tag_block(..., feature_band, feature_tags, feature_angle)`; `query_many(..., feature_distances, feature_angle)` | Corners (turning above `angle`, default 30 degrees; ends of open loops), cell distances to them, point distances to them (CPU and device). See [feature edges](/guide/advanced/feature-edges#2-d-corners) (issue #61). |
| `thickness(thickness, gap, method, cone_angle, rays, invalid, status, separation_angle)`; `thin_regions(thickness, gap, h, factor, region, region_length, status)`; `build_thin(thickness, gap, status)`, `box_thin_many(_device)(lo, hi, factor, touched, status)`, `tag_block(..., thin_factor, thin_band, thin_tags)` | Thickness and gap per segment from its midpoint (a ray fan or the medial disc), runs of segments thinner than `factor * h` along each loop, and the AMR criterion (exact box-segment test, CPU and device). See [thickness and gaps](/guide/advanced/thickness#2-d-curves) (issue #62). |
| `self_intersections(pairs, n)` | Segment pairs that cross, touch or fold back, decided exactly. |
| `repair(tolerance, gaps_closed, vertices_merged, loops_removed, loops_reversed, loops_open)` | Join and close open chains within the tolerance, merge vertices, drop degenerate loops, orient closed loops by nesting depth. |

## Queries

| Method | Purpose |
|---|---|
| `distance(p, is_signed, sign_algorithm)` | (Signed) distance of p = [u, v]; negative inside. |
| `compute_distance(p, distance, closest_point, closest_segment, closest_vertex, is_signed, use_brute_force, sign_algorithm, normal, body)` | Distance with the closest point, segment and vertex (0 in a segment interior), normal and body. |
| `winding_number(p, use_brute_force, beta)` | Generalized winding number: fast hierarchical, or the exact angle sum. |
| `winding_number_exact(p)` | Exact integer winding number over the closed loops. |
| `is_point_inside(p, sign_algorithm)` | Inside test: pseudo-normal sign (default) or winding number > 0.5. |
| `distance_many(points, distances, is_signed, sign_algorithm, status)` | Batch over `points(2, n)`, OpenMP; bit-exact with `compute_distance`. |
| `enter_device`, `exit_device`, `distance_many_device(...)`, `get_is_on_device` | The same batch on the device (OpenACC). |
| `query_many(points, distances, closest_points, segments, vertices, normals, bodies, is_signed, sign_algorithm, ambiguous, second_foot_points, rel_tolerance, foot_tolerance, status, pose, curvatures)` | Rich batch (issue #37); `query_many_device` on the device. `curvatures`: curvature at the foot point (issue #42, CPU). |
| `probe_many(points, lengths, probe_points, visible, foot_points, normals, status)` | Probe points along the normal; exact visibility. |
| `image_points(points, images, visible, foot_points, status)` | Mirror through the closest point; exact visibility. |
| `get_bodies_number`, `loop_body(l)` | Bodies: outer loops with their holes; islands and open loops are their own. |
| `box_classify_many(lo, hi, codes, radius, sign_algorithm, inside, status)` | Classify boxes `lo(2, n)`, `hi(2, n)`: `BOX_CUT` (exact), `BOX_NEAR`, `BOX_INSIDE`, `BOX_OUTSIDE` (issue #38); `box_classify_many_device` on the device. |
| `box_intersects_curve(lo, hi)` | True if the closed box meets the curve (exact). |
| `line_crossings(axis, a, t0, t1, t, segments, entering, n, status)` | Exact crossings of an axis-aligned line (issue #40); `line_crossings_many(_device)` for many lines into fixed buffers. |
| `cell_wall_fractions(x, y, frac, inside, bad_lines, status)` | Wall fractions toward the four neighbours of every cell centre, parity inside flags. |
| `move(pose)`, `state_changes(pose0, pose1, points, fresh, dead, sign_algorithm, checked, status)` | Rigid in-plane move with a tree refit; fresh and dead points between two poses (issue #41). `distance_many`, `query_many`, `is_point_inside` and the device twins take an optional `pose` (`rigid_pose_2d_t`). The swept geometry between two poses is the free routine `swept_cells_2d(curve, ...)`. See [Moving bodies](/guide/advanced/moving-bodies). |
| `wall_distance_many(points, distances, bodies, cutoff, status)` | Unsigned wall distance and body, (cutoff, 0) beyond the cutoff (issue #42). |
| `curve_markers(dx, centroids, lengths, normals, spacing_fraction, spacings, status)` | Exact arc-length resampling into nint(L/dx) pieces per loop: markers on the polyline, lengths, normals (issue #42). |
| `cut_cells(x, y, cut, status)` | Exact cut cells (`cut_cells_2d_t`): fractions, apertures, centroids, wetted boundary, pieces. See [Cut-cell geometry](/guide/advanced/cut-cells). |
| `tag_block(x, y, band, tags, sign_algorithm, status)` | Signed tag bands on a block of cells given its edges. See [AMR tagging](/guide/advanced/amr-tagging). |

`sign_algorithm` is `SIGN_PSEUDO_NORMAL` or `SIGN_FAST_WINDING_NUMBER`; the
latter is robust on open loops and gaps. The point queries default to the
pseudo-normal; `box_classify_many(_device)` and `tag_block` default to the fast
winding number (issue #39).

## Related

- `slice_surface(surface, normal_axis, level, curve, merge_collinear, open_loops, status)`
  builds a curve from the section of a `surface_stl_object`.
- `distance_field_object%initialize(x, y)` (no `z`) makes a 2-D field;
  `compute(curve=...)` fills it.
