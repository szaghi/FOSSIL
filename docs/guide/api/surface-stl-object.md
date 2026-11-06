---
title: surface_stl_object
---

# `surface_stl_object`

The triangulated surface — the library's centrepiece. A `surface_stl_object`
owns:

- a private array of [`facet_object`](/guide/api/facet-object) holding the
  geometry,
- a private [`aabb_tree_object`](/guide/api/aabb-tree-object) accelerating
  distance and inside queries,
- a private unique-vertex pool deduplicating coincident vertices,
- cached scalar properties (volume, centroid, bounding box, manifold counts)
  populated by [`analyze`](#analyze).

Every public component is `private`. All state is read through getters
([`get_facets_number`](#getters-and-predicates),
[`get_bmin`](#getters-and-predicates),
[`get_volume`](#getters-and-predicates), …) and mutated through type-bound
procedures ([`load_from_file`](#load_from_file),
[`sanitize`](#sanitize),
[`translate`](#translate), …).

The TBPs below are grouped by user-flow: I/O, repair, analysis, queries,
manipulation, and a final catch-all section for the read-only accessors and
validity predicates.

[[toc]]

---

## I/O

### `load_from_file`

```fortran
call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
```

**Purpose.** Read an STL file from disk into the surface. On return the surface
holds the parsed facets and a full set of cached scalars: bounding box,
connectivity, volume, centroid, AABB tree, and per-vertex / per-edge
pseudo-normals. In other words, `load_from_file` is `read + analyze` fused
into one call — you can immediately query [`distance`](#distance) without an
intermediate [`analyze`](#analyze).

The format is detected automatically when `guess_format=.true.`: the file size
is compared against the size implied by the binary header's facet count, and
the format that matches wins. Pass `is_ascii=.true./.false.` explicitly if you
already know the format and want to skip detection.

**Arguments.**

| Argument                  | Intent       | Type               | Notes                                                                                              |
| ------------------------- | ------------ | ------------------ | -------------------------------------------------------------------------------------------------- |
| `file_name`               | `in`         | `character(*)`     | (required) Path to the STL file.                                                                   |
| `is_ascii`                | `in`, opt.   | `logical`          | Force ASCII or binary parsing. Ignored when `guess_format=.true.`. Defaults to ASCII otherwise.     |
| `guess_format`            | `in`, opt.   | `logical`          | Auto-detect ASCII vs binary from the file size. Recommended for end-user-supplied paths.            |
| `clip_min`, `clip_max`    | `in`, opt.   | `type(vector_R8P)` | Both must be present together. Facets entirely outside the AABB `[clip_min, clip_max]` are dropped during load. Useful for opening a single subregion of a huge file without paying for full I/O. |
| `aabb_refinement_levels`  | `in`, opt.   | `integer(I4P)`     | Octree depth, or `AABB_AUTO_REFINEMENT`. Ignored when the kind is `AABB_TREE_SAH_BVH`.              |
| `aabb_tree_kind`          | `in`, opt.   | `integer(I4P)`     | `AABB_TREE_SAH_BVH` (default) or `AABB_TREE_OCTREE`.                                                |
| `status`                  | `out`, opt.  | `integer(I4P)`     | `STATUS_OK`, `STATUS_FILE_NOT_FOUND`, `STATUS_INVALID_INPUT`, or `STATUS_ALLOC_FAIL`.                |

**Example.** Load the bundled cube, print summary, and check it's a valid solid.

```fortran
program ex_load_from_file
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: stat

call surface%load_from_file(file_name='cube.stl', guess_format=.true., status=stat)
if (stat /= STATUS_OK) error stop 'load failed'
print '(A, I0)',    'facets : ', surface%get_facets_number()
print '(A, ES12.5)','volume : ', surface%get_volume()
print '(A, L1)',    'is_volume? ', surface%is_volume()
end program ex_load_from_file
```

**Example.** Load with an octree (overriding the default BVH).

```fortran
program ex_load_octree
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.,        &
                            aabb_tree_kind=AABB_TREE_OCTREE,                  &
                            aabb_refinement_levels=AABB_AUTO_REFINEMENT)
print '(A, I0)', 'octree levels = ', surface%aabb%get_refinement_levels()
end program ex_load_octree
```

**Notes.**

- NaN / Inf coordinate values cause an immediate `STATUS_INVALID_INPUT`; the
  surface is left untouched. This is policy: load refuses corrupt input rather
  than passing it downstream to be mis-repaired.
- A successful load **does not** imply a clean mesh. The file may still have
  degenerate, duplicate, or inconsistently-oriented facets. Call
  [`sanitize`](#sanitize) immediately after loading user-supplied input.
- Repeated calls reuse the surface — the previous content is destroyed first.

**See also.** [`save_into_file`](#save_into_file),
[`sanitize`](#sanitize),
[`analyze`](#analyze),
[`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh),
[`STATUS_*`](/guide/api/constants#status-codes).

---

### `save_into_file`

```fortran
call surface%save_into_file(file_name='out.stl', is_ascii=.false.)
```

**Purpose.** Write the surface to disk in STL format. The current header (set
via [`set_header`](#getters-and-predicates) or carried over from
[`load_from_file`](#load_from_file)) is written first.

**Arguments.**

| Argument    | Intent       | Type             | Notes                                                              |
| ----------- | ------------ | ---------------- | ------------------------------------------------------------------ |
| `file_name` | `in`         | `character(*)`   | (required) Output path. The file is truncated.                     |
| `is_ascii`  | `in`, opt.   | `logical`        | Default `.true.` (ASCII output). Set `.false.` for binary STL.     |
| `status`    | `out`, opt.  | `integer(I4P)`   | `STATUS_OK` or `STATUS_FILE_OPEN_FAIL`.                            |

**Example.** Round-trip through binary STL.

```fortran
program ex_save_into_file
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%save_into_file(file_name='cube-binary.stl', is_ascii=.false.)
call surface%save_into_file(file_name='cube-ascii.stl',  is_ascii=.true.)
end program ex_save_into_file
```

**Notes.**

- Binary is ~6× smaller and ~20× faster to read back than ASCII. Use it for
  intermediate data; reserve ASCII for human inspection.
- The STL format is lossy w.r.t. connectivity — the file stores only triangles.
  After a save / load round trip you must call [`analyze`](#analyze) (or rely
  on `load_from_file`'s internal analyse) to rebuild connectivity.

**See also.** [`load_from_file`](#load_from_file),
[`set_header`](#getters-and-predicates).

---

## Repair pipeline

Real-world STL files are routinely broken — zero-area slivers, duplicate
triangles, unwelded vertices, mixed winding. FOSSIL's repair API is layered:
a single entry point ([`sanitize`](#sanitize)) that orchestrates every pass in
the right order, plus the individual passes exposed so you can run them
manually when the orchestrator's policy is wrong for your data.

### `sanitize`

```fortran
call surface%sanitize
```

**Purpose.** Full repair pipeline. Run this immediately after every
user-supplied load. The pipeline is the result of a long-running iteration on
the order in which repair passes must compose; running them in a different
order can leave subtle corruptions.

The fixed pipeline is:

1. `remove_degenerate_facets` — drop zero-area / sliver triangles first so
   their NaN normals do not poison later steps.
2. `connect_nearby_vertices` — only if disconnected edges remain. Snaps
   coincident vertices via union-find on a tolerance.
3. `analyze` — rebuild bounding box, connectivity, volume, centroid, AABB tree.
4. `remove_duplicate_facets` — drop literal duplicates (orientation-agnostic).
5. `analyze` again, if step 4 removed anything, so winding fixup starts from
   a clean state.
6. `sanitize_normals` — BFS-propagate winding consistency from facet 1, then
   flip globally if the signed volume came out negative. If anything was
   flipped, `analyze` runs once more so the pseudo-normals and the AABB
   distance payload match the new winding (without it, an inside-out input
   reported its interior as outside).
7. **Warnings** — one line to `stderr` for every defect class with non-zero
   count. Suppress them by redirecting unit `error_unit`; preserve them in
   production logs so you notice when your inputs degrade.

**Arguments.**

| Argument       | Intent       | Type           | Notes                                                                                |
| -------------- | ------------ | -------------- | ------------------------------------------------------------------------------------ |
| `do_analysis`  | `in`, opt.   | `logical`      | If `.true.`, run an initial `analyze` before step 1. Default `.false.` — the caller normally arrives here straight from `load_from_file`, which has already analysed. |
| `status`       | `out`, opt.  | `integer(I4P)` | `STATUS_OK`; reserved for future status codes.                                       |
| `fill_holes`   | `in`, opt.   | `logical`      | Default `.false.`. Close small holes after orientation, see [`fill_holes`](#fill_holes); orientation runs once more if anything was filled. |
| `max_hole_edges` | `in`, opt. | `integer(I4P)` | Passed to `fill_holes` (default 200).                                                 |
| `max_hole_size`  | `in`, opt. | `real(R8P)`    | Passed to `fill_holes` (default 0.1 of the bounding-box diagonal).                    |
| `refine_holes`   | `in`, opt. | `logical`      | Passed to `fill_holes` as `refine` (default `.false.`): patches refined to the loop's edge density. |
| `fair_holes`     | `in`, opt. | `logical`      | Passed to `fill_holes` as `fair` (default `.false.`, implies refine): thin-plate fairing of the patches. |

**Example.** Load a deliberately broken mesh and verify repair.

```fortran
program ex_sanitize
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube-inconsistent.stl', guess_format=.true.)
print '(A, L1)', 'is_volume before sanitize? ', surface%is_volume()
call surface%sanitize
print '(A, L1)', 'is_volume after  sanitize? ', surface%is_volume()
print '(A, I0)', 'degenerate facets removed = ', surface%get_degenerate_facets_removed()
print '(A, I0)', 'duplicate  facets removed = ', surface%get_duplicate_facets_removed()
end program ex_sanitize
```

**Notes.**

- `sanitize` is **idempotent on a clean mesh**: a second call is a no-op
  (every pass detects zero defects). Use this property to assert cleanliness
  in tests.
- The connect / duplicate-removal tolerances are derived from the bounding-box
  diagonal. If your mesh has wildly different feature scales, run the
  individual passes manually with explicit tolerances.

**See also.** [`remove_degenerate_facets`](#remove_degenerate_facets),
[`remove_duplicate_facets`](#remove_duplicate_facets),
[`connect_nearby_vertices`](#connect_nearby_vertices),
[`sanitize_normals`](#sanitize_normals),
[`is_volume`](#getters-and-predicates).

---

### `remove_degenerate_facets`

```fortran
call surface%remove_degenerate_facets
```

**Purpose.** Drop facets whose area is below a tolerance relative to the
mesh's bounding-box diagonal. Degenerate facets have undefined normals (often
NaN) and contaminate every downstream computation; removing them is the first
step of any repair pipeline.

**Arguments.** None.

**Example.**

```fortran
program ex_remove_degenerate
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%remove_degenerate_facets
print '(A, I0)', 'degenerate facets removed = ', surface%get_degenerate_facets_removed()
end program ex_remove_degenerate
```

**Notes.**

- The count is queryable through
  [`get_degenerate_facets_removed`](#getters-and-predicates) for post-run
  inspection.
- Removal does not re-run [`analyze`](#analyze) — call it yourself if you need
  fresh connectivity afterwards. The orchestrator [`sanitize`](#sanitize) does
  this automatically.

**See also.** [`remove_duplicate_facets`](#remove_duplicate_facets),
[`sanitize`](#sanitize).

---

### `remove_duplicate_facets`

```fortran
call surface%remove_duplicate_facets
```

**Purpose.** Drop facets that duplicate another facet up to vertex permutation
(including reversed winding). Duplicates skew volume, double-count surface
area, and bias the winding-consistency BFS.

**Arguments.** None.

**Example.**

```fortran
program ex_remove_duplicate
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%remove_duplicate_facets
print '(A, I0)', 'duplicate facets removed = ', surface%get_duplicate_facets_removed()
end program ex_remove_duplicate
```

**Notes.**

- Should run *after* [`remove_degenerate_facets`](#remove_degenerate_facets),
  because degenerate slivers tend to come in duplicate pairs.
- The detection key is the *sorted triple of unique vertex ids*; this catches
  reversed-winding duplicates that a naïve vertex-coordinate comparison would
  miss.

**See also.** [`sanitize`](#sanitize),
[`remove_degenerate_facets`](#remove_degenerate_facets).

---

### `connect_nearby_vertices`

```fortran
call surface%connect_nearby_vertices
```

**Purpose.** Snap coincident-but-not-shared vertices via union-find. STL is
soup-of-triangles by spec: nothing in the format enforces that two triangles
sharing an edge use the same vertex objects. After this pass, the
unique-vertex pool reflects geometric coincidence, and connectivity (computed
by [`analyze`](#analyze)) can resolve all the edges that visually look shared
but were previously disconnected.

**Arguments.** None.

**Example.**

```fortran
program ex_connect_nearby
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube_disconnect.stl', guess_format=.true.)
call surface%connect_nearby_vertices
call surface%analyze
print '(A, L1)', 'manifold after connect? ', surface%is_manifold()
end program ex_connect_nearby
```

**Notes.**

- Uses a bounding-box-diagonal-relative tolerance. The tolerance is fixed
  internally; for adversarial inputs with feature-scale separation you should
  instead pre-snap your vertices.
- Always re-run [`analyze`](#analyze) afterwards (the orchestrator does this
  automatically).

**See also.** [`sanitize`](#sanitize),
[`analyze`](#analyze).

---

### `fill_holes`

```fortran
call surface%fill_holes(max_hole_size=0.05_R8P, holes_filled=n, holes_skipped=m)
```

**Purpose.** Close small holes: every closed loop of boundary edges with at
most `max_hole_edges` edges and a bounding-box diagonal of at most
`max_hole_size` times the surface's is triangulated by Liepa's minimum-weight
dynamic programming, oriented like the surrounding facets. The surface is
re-analysed afterwards. See [Hole filling](/guide/advanced/hole-filling).

**Arguments.**

| Argument         | Intent       | Type           | Notes                                                              |
| ---------------- | ------------ | -------------- | ------------------------------------------------------------------ |
| `max_hole_edges` | `in`, opt.   | `integer(I4P)` | Largest loop filled, in edges (default 200; the DP is O(n³)).      |
| `max_hole_size`  | `in`, opt.   | `real(R8P)`    | Largest loop bbox diagonal / surface bbox diagonal (default 0.1).  |
| `holes_filled`   | `out`, opt.  | `integer(I4P)` | Loops filled.                                                      |
| `holes_skipped`  | `out`, opt.  | `integer(I4P)` | Loops left open by the thresholds.                                 |
| `open_edges`     | `out`, opt.  | `integer(I4P)` | Boundary edges that do not close into a loop.                      |
| `status`         | `out`, opt.  | `integer(I4P)` | `STATUS_OK`.                                                       |

**See also.** [`sanitize`](#sanitize), [`is_watertight`](#getters-and-predicates).

---

### `sanitize_normals`

```fortran
call surface%sanitize_normals
call surface%sanitize_normals(flipped_count=n)
```

**Purpose.** Make every facet's winding consistent with its neighbours, then
flip the whole surface if the signed volume is negative — yielding outward
normals on a closed body. This is the precondition for the default
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal) sign algorithm.

**Arguments.**

| Argument        | Intent       | Type           | Notes                                                                 |
| --------------- | ------------ | -------------- | --------------------------------------------------------------------- |
| `flipped_count` | `out`, opt.  | `integer(I4P)` | Facets whose winding was reversed; a global flip counts every facet.  |

::: warning Call analyze after a flip
A flip reverses facet windings, but the pseudo-normals and the AABB distance
payload built by [`analyze`](#analyze) still describe the old winding, so
signed-distance and inside queries are wrong until `analyze` runs again.
[`sanitize`](#sanitize) does this for you. Called on its own,
`sanitize_normals` leaves it to you: re-run `analyze` whenever `flipped_count`
is positive.
:::

**Example.**

```fortran
program ex_sanitize_normals
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: n_flipped

call surface%load_from_file(file_name='cube-inconsistent.stl', guess_format=.true.)
call surface%sanitize_normals(flipped_count=n_flipped)
if (n_flipped > 0) call surface%analyze   ! rebuild pseudo-normals for the new winding
print '(A, I0, A, ES12.5)', 'flipped ', n_flipped, ' facet(s); volume = ', surface%get_volume()
end program ex_sanitize_normals
```

**Notes.**

- The BFS starts at facet 1 and walks every connected component in turn, so
  each component ends up internally consistent. Each component is then
  oriented as a whole:
  - a **closed** component (every facet edge has a neighbour) by nesting
    parity: inside an even number of other closed shells it must enclose
    positive volume (an outer boundary), inside an odd number negative (a
    cavity wall). Disjoint bodies are oriented independently, and a hollow
    body keeps its cavity wall facing inward. Containment is tested by
    generalized winding number, which does not depend on the other shell's
    orientation;
  - the **open** components (a body with holes, a sheet) together, by the sign
    of their summed volume. With no closed component this is the whole-surface
    volume rule, and on open shells the volume is not meaningful, so the flip
    may or may not happen.
- Cost: for every pair of closed shells whose bounding boxes nest, one
  winding-number evaluation over the outer shell's facets.

**See also.** [`reverse_normals`](#reverse_normals),
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal),
[`sanitize`](#sanitize).

---

### `reverse_normals`

```fortran
call surface%reverse_normals
```

**Purpose.** Flip every facet's normal and reverse its vertex winding. Useful
when you load a mesh that is correctly *consistent* but oriented inward
(volume negative) and you want to fix it without running the full BFS in
[`sanitize_normals`](#sanitize_normals).

**Arguments.** None.

**Example.**

```fortran
program ex_reverse_normals
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: v_before, v_after

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
v_before = surface%get_volume()
call surface%reverse_normals
call surface%analyze
v_after  = surface%get_volume()
print '(A, ES12.5, A, ES12.5)', 'volume: before=', v_before, ' after=', v_after
end program ex_reverse_normals
```

**Notes.**

- Re-run [`analyze`](#analyze) afterwards to refresh the signed volume.

**See also.** [`sanitize_normals`](#sanitize_normals).

---

## Analysis

### `analyze`

```fortran
call surface%analyze
```

**Purpose.** Compute, in one pass, every derived quantity that the surface
caches: bounding box, edge connectivity, signed volume, centroid, AABB tree,
and per-edge / per-vertex pseudo-normals. Called automatically by
[`load_from_file`](#load_from_file), [`sanitize`](#sanitize), and
[`clip`](#clip); you only call it manually after a transformation that you
chose **not** to follow with `recompute_metrix=.true.`.

**Arguments.**

| Argument                 | Intent       | Type           | Notes                                                                              |
| ------------------------ | ------------ | -------------- | ---------------------------------------------------------------------------------- |
| `aabb_refinement_levels` | `in`, opt.   | `integer(I4P)` | Octree depth or `AABB_AUTO_REFINEMENT`. Ignored by `AABB_TREE_SAH_BVH`.             |
| `aabb_tree_kind`         | `in`, opt.   | `integer(I4P)` | `AABB_TREE_SAH_BVH` or `AABB_TREE_OCTREE`. Defaults to the surface's current kind.  |
| `status`                 | `out`, opt.  | `integer(I4P)` | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                                |

**Example.** Translate, then re-analyse to refresh the AABB tree and centroid.

```fortran
program ex_analyze
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

type(vector_R8P)         :: c

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%translate(delta=10._R8P * ex_R8P)
call surface%analyze
c = surface%get_centroid()
print '(A, 3ES12.5)', 'new centroid = ', c%x, c%y, c%z
end program ex_analyze
```

**Notes.**

- Switching `aabb_tree_kind` between calls *does* rebuild the tree.
- Pure rotations / translations / mirrors leave connectivity intact; if
  performance matters you can call `recompute_metrix=.true.` on the
  transformation directly and avoid the full analyse.

**See also.** [`build_connectivity`](#build_connectivity),
[`statistics`](#statistics),
[`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh).

---

### `build_connectivity`

```fortran
call surface%build_connectivity
```

**Purpose.** Rebuild only the facet edge-adjacency (which facet shares which
edge with which). Cheaper than [`analyze`](#analyze) when you don't need the
full bounding box / volume / AABB refresh.

**Arguments.** None.

**Example.**

```fortran
program ex_build_connectivity
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%build_connectivity
print '(A, I0)', 'non-manifold edges = ', surface%get_non_manifold_edges_number()
end program ex_build_connectivity
```

**See also.** [`analyze`](#analyze).

---

### `statistics`

```fortran
print '(A)', surface%statistics()
```

**Purpose.** Return a multi-line formatted string describing the surface:
bounding box, volume, centroid, facet count, manifold / degenerate /
duplicate counts, and AABB-tree configuration. The format is human-readable;
do not parse it.

**Arguments.** None.

**Example.**

```fortran
program ex_statistics
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
end program ex_statistics
```

**See also.** [`analyze`](#analyze).

---

### `mass_properties`

```fortran
call surface%mass_properties(mass, center, inertia, principal_moments, principal_axes, density, about, volume, status)
```

**Purpose.** Mass properties of the closed body bounded by the surface (issue #59): mass, centre of mass, inertia
tensor about the centre (or about `about`), principal moments (ascending) and principal axes (columns,
right-handed). Exact polyhedral integrals, compensated sums about the bounding-box centre: the only error is
round-off.

**Arguments.** All optional.

| Argument | Intent | Meaning |
|---|---|---|
| `mass` | out | `density` × enclosed volume |
| `center(3)` | out | centre of mass |
| `inertia(3,3)` | out | inertia tensor about `center`, or about `about` |
| `principal_moments(3)` | out | eigenvalues of the tensor about the centre, ascending |
| `principal_axes(3,3)` | out | the matching eigenvectors as columns, det = +1 |
| `density` | in | uniform density, default 1 |
| `about(3)` | in | the point of `inertia` |
| `volume` | out | enclosed volume |
| `status` | out | `STATUS_OK`; `STATUS_INVALID_INPUT` if not watertight or volume ≤ 0 (outputs zero) |

**Example.**

```fortran
program ex_mass_properties
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: body
real(R8P)                :: mass, center(3), moments(3), axes(3,3)
integer(I4P)             :: status

call body%load_from_file(file_name='cube.stl', guess_format=.true.)
call body%sanitize
call body%mass_properties(mass=mass, center=center, principal_moments=moments, principal_axes=axes, &
                          density=7850._R8P, status=status)
if (status == STATUS_OK) print '(A,ES12.5,A,3F9.4)', 'mass ', mass, ' centre ', center
end program ex_mass_properties
```

**See also.** [Mass properties](/guide/advanced/mass-properties) (method, validation), [`move`](#move).

---

### `deviation`

```fortran
call surface%deviation(other, lower, upper, mean, rms, per_vertex, tolerance, samples_per_area, max_triangles, &
                       samples, triangles, device, status)
call surface%hausdorff(other, lower, upper, tolerance, samples_per_area, max_triangles, device, status)
```

**Purpose.** How far this surface lies from `other` (issue #60): `lower <= h <= upper` brackets the one-sided deviation
h(self, other) = sup over x in self of d(x, other). `lower` comes from deterministic samples (vertices, edge midpoints,
sub-triangle centroids); `upper` is certified by a branch and bound over the facets (Lipschitz and convexity bounds of
d), to `tolerance` (default 1e-3 of the larger bounding-box diagonal) or until `max_triangles`. `mean` and `rms` are
area-weighted estimates of d over self; `per_vertex` holds d at every vertex of the vertex pool. `hausdorff` is the
symmetric max(h(self, other), h(other, self)). With `device=.true.` the queries run on the device (`other%enter_device`
first), bit-identical to the CPU.

**Arguments.** All optional except `other`; `status` is `DEV_STATUS_OK`, `DEV_STATUS_BAD_INPUT` (an empty surface, or
`device` without `enter_device`) or
`DEV_STATUS_NOT_CONVERGED` (the work cap stopped it; the bounds are valid).

**Example.**

```fortran
program ex_deviation
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: before, after
real(R8P)                :: lower, upper
integer(I4P)             :: status

call before%load_from_file(file_name='cube.stl', guess_format=.true.)
after = before
call after%decimate(target_facets=before%get_facets_number() / 2)
call after%deviation(before, lower=lower, upper=upper, status=status)
if (status == DEV_STATUS_OK) print '(ES10.3,A,ES10.3)', lower, ' <= h <= ', upper
end program ex_deviation
```

**See also.** [Surface deviation](/guide/advanced/surface-deviation) (method, validation).

---

### `feature_edges`

```fortran
call surface%feature_edges(edges, kinds, corners, polyline_start, polyline_vertices, angle, corner_angle, &
                           include_boundary, include_nonmanifold, degenerate, status)
```

**Purpose.** The feature edges of the surface (issue #61): sharp edges (dihedral above `angle`, default 30 degrees),
boundary edges and non-manifold edges, their kinds (`FEATURE_EDGE_SHARP`, `_BOUNDARY`, `_NONMANIFOLD`), the corners and
the polylines between them (CSR), as vertex-pool ids.

**Arguments.** All optional; `status` is `STATUS_OK` or `STATUS_INVALID_INPUT` (empty surface, angle outside [0, 180]).

**Example.**

```fortran
program ex_feature_edges
use fossil
use penf, only : I4P
implicit none
type(surface_stl_object)  :: part
integer(I4P), allocatable :: kinds(:), corners(:)

call part%load_from_file(file_name='cube.stl', guess_format=.true.)
call part%feature_edges(kinds=kinds, corners=corners)
print '(A,I0,A,I0)', 'sharp edges ', count(kinds == FEATURE_EDGE_SHARP), ', corners ', size(corners)
end program ex_feature_edges
```

**See also.** [Feature edges](/guide/advanced/feature-edges) (definitions, method, validation).

---

### `build_features`, `box_feature_many`

```fortran
call surface%build_features(angle, corner_angle, include_boundary, include_nonmanifold, status)
call surface%box_feature_many(lo, hi, touched, status)         ! and box_feature_many_device
```

**Purpose.** `build_features` stores the feature edges (same arguments as `feature_edges`) in a segment tree for the
feature queries (issue #61); `analyze` and `move` drop it. `box_feature_many` marks the boxes `lo(3, n)`, `hi(3, n)`
that touch a feature edge (conservative, never misses a contact); the device twin needs `build_features` before
`enter_device`. `tag_block` takes `feature_band=` and `feature_tags=` for the Chebyshev cell distance to them.

**See also.** [Feature edges](/guide/advanced/feature-edges#refinement-along-the-features),
[AMR tagging](/guide/advanced/amr-tagging).

---

### `thickness`

```fortran
call surface%thickness(thickness, gap, method, cone_angle, rays, vertex_thickness, vertex_gap, invalid, status, &
                       separation_angle)
```

**Purpose.** The local thickness of the body behind every facet and the gap in front of it (issue #62), `huge` where
undefined or where nothing is in front. `method=THICKNESS_RAY` (default) is the shape diameter: the lower median of
`rays` ray lengths (default 31) in a cone of full aperture `cone_angle` (default 60°) around -n (thickness) or +n
(gap) from the facet centroid. `method=THICKNESS_MEDIAL` is the shrinking ball: the diameter of the largest ball
tangent to the facet at its centroid with no surface point inside (the local feature size), never above the axial
chord; `separation_angle` (default 30°) is its denoising angle. Vertex fields are area-weighted averages of the finite
facet values. OpenMP over the facets, independent of the thread count.

**Arguments.** All optional, at least one output; `invalid` counts the facets whose rays did not give a value (an open
or misoriented surface); `status` is `STATUS_OK` or `STATUS_INVALID_INPUT` (empty surface, `rays < 1`, `cone_angle`
or `separation_angle` outside [0, 180), no output requested).

**Example.**

```fortran
program ex_thickness
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: part
real(R8P), allocatable   :: thickness(:), gap(:)
integer(I4P)             :: invalid

call part%load_from_file(file_name='cube.stl', guess_format=.true.)
call part%thickness(thickness=thickness, gap=gap, invalid=invalid)
print '(A,ES10.3,A,I0)', 'thinnest wall ', minval(thickness), ', facets without a value ', invalid
end program ex_thickness
```

**See also.** [Thickness and gaps](/guide/advanced/thickness).

---

### `thin_regions`, `build_thin`, `box_thin_many`

```fortran
call surface%thin_regions(thickness, gap, h, factor, region, region_area, status)
call surface%build_thin(thickness, gap, status)
call surface%box_thin_many(lo, hi, factor, touched, status)     ! and box_thin_many_device
```

**Purpose.** The under-resolution report and the AMR criterion (issue #62); the local size of a facet is
min(thickness, gap) of the arrays given (from `thickness`). `thin_regions` marks the facets whose local size is below
`factor * h` and groups them into regions connected through shared edges, sorted by decreasing area: `region(nf)` is 0
or the region id, `region_area(nr)` the areas. `build_thin` stores the local sizes in the SAH BVH (per triangle, with a
minimum per node) for `box_thin_many`, which marks the boxes `lo(3, n)`, `hi(3, n)` meeting a facet thinner than
`factor` times the largest side of the box (conservative); `analyze` drops the stored sizes, `move` keeps them. The
device twin needs `build_thin` before `enter_device`. `tag_block` takes `thin_factor=`, `thin_band=` and
`thin_tags=` for the Chebyshev cell distance to those cells.

**See also.** [Thickness and gaps](/guide/advanced/thickness#under-resolution-report),
[AMR tagging](/guide/advanced/amr-tagging).

---

## Distance and inside queries

### `distance`

```fortran
d = surface%distance(point=p)
```

**Purpose.** Return the (minimum) distance from a point to the triangulated
surface. By default the result is **unsigned and squared** — useful for
nearest-facet ranking without a square-root cost. Pass `is_square_root=.true.`
for Euclidean distance and `is_signed=.true.` to apply one of the
[`SIGN_*`](/guide/api/constants#sign-algorithms) sign algorithms.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                                                          |
| ----------------- | ------------ | ------------------ | ---------------------------------------------------------------------------------------------- |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                                                        |
| `is_signed`       | `in`, opt.   | `logical`          | Default `.false.` (unsigned).                                                                  |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | `SIGN_PSEUDO_NORMAL` (default), `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`, `SIGN_FAST_WINDING_NUMBER`. |
| `is_square_root`  | `in`, opt.   | `logical`          | Default `.false.` — returns d². Set `.true.` for Euclidean distance.                          |

**Returns.** `real(R8P)` — the distance.

**Example.** All four common combinations from a single query point.

```fortran
program ex_distance
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
p = 2._R8P * ex_R8P
print '(A, ES12.5)', 'unsigned d^2     = ', surface%distance(point=p)
print '(A, ES12.5)', 'unsigned d       = ', surface%distance(point=p, is_square_root=.true.)
print '(A, ES12.5)', 'signed   d (pn)  = ', surface%distance(point=p, is_signed=.true.,         &
                                                              is_square_root=.true.)
print '(A, ES12.5)', 'signed   d (ray) = ', surface%distance(point=p, is_signed=.true.,         &
                                                              sign_algorithm=SIGN_RAY_INTERSECTIONS,&
                                                              is_square_root=.true.)
end program ex_distance
```

**Notes.**

- This is a thin wrapper over [`compute_distance`](#compute_distance) — use
  that one when you want the closest facet / edge / vertex indices in
  addition to the distance.
- Tree-accelerated: O(log N) expected on uniform meshes. The traversal is
  best-first with d² pruning, so it never visits a node whose AABB is farther
  than the current best.
- Bit-exact against brute-force on the regression suite, by construction —
  the pruning predicate is conservative.

**See also.** [`compute_distance`](#compute_distance),
[`is_point_inside`](#is_point_inside),
[`SIGN_*`](/guide/api/constants#sign-algorithms).

---

### `compute_distance`

```fortran
call surface%compute_distance(point=p, distance=d)
```

**Purpose.** Same as [`distance`](#distance) but with the distance returned
through an `intent(out)` argument plus four optional outputs that locate the
closest point on the mesh: facet index, local edge index, and local vertex
index. Use this when you need to know *which* feature of the mesh produced the
closest distance — for visualisation, classification, or as input to a
projection.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                                                                 |
| ----------------- | ------------ | ------------------ | ----------------------------------------------------------------------------------------------------- |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                                                               |
| `distance`        | `out`        | `real(R8P)`        | (required) Distance from `point` to surface.                                                          |
| `is_signed`       | `in`, opt.   | `logical`          | Default `.false.`.                                                                                    |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | `SIGN_*`; default `SIGN_PSEUDO_NORMAL`.                                                               |
| `is_square_root`  | `in`, opt.   | `logical`          | Default `.false.` — d² is returned.                                                                  |
| `facet_index`     | `out`, opt.  | `integer(I4P)`     | 1-based index of the facet holding the closest point. 0 if the surface has no facets.                 |
| `edge_index`      | `out`, opt.  | `integer(I4P)`     | Local edge id (1, 2, or 3) when the closest point lies on an edge interior; 0 otherwise.              |
| `vertex_index`    | `out`, opt.  | `integer(I4P)`     | Local vertex id (1, 2, or 3) when the closest point coincides with a vertex; 0 otherwise.             |
| `closest_point`   | `out`, opt.  | `type(vector_R8P)` | Closest point on the surface (issue #37). Does not change `distance`.                                 |
| `normal`          | `out`, opt.  | `type(vector_R8P)` | Pseudo-normal of the closest feature (facet, edge or vertex; issue #37).                              |

**Example.** Identify the closest feature class on every query.

```fortran
program ex_compute_distance
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p
real(R8P)                :: d
integer(I4P)             :: fid, eid, vid
character(len=8)         :: feature

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
p = 2._R8P * ex_R8P
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=fid, edge_index=eid, vertex_index=vid)
if (vid /= 0) then
   feature = 'vertex'
else if (eid /= 0) then
   feature = 'edge'
else
   feature = 'face'
end if
print '(A,ES12.5,A,A,A,I0)', 'd=', d, ' closest feature=', trim(feature), &
                              ' on facet ', fid
end program ex_compute_distance
```

**Notes.**

- The `edge_index` / `vertex_index` outputs are mutually exclusive with
  "interior face": at most one of them is non-zero for a given query.
- When `sign_algorithm = SIGN_PSEUDO_NORMAL`, sign and closest-feature are
  computed in a single tree traversal. With `SIGN_RAY_INTERSECTIONS` or
  `SIGN_SOLID_ANGLE` a second O(N) pass classifies the sign; with
  `SIGN_FAST_WINDING_NUMBER` a second, hierarchical (O(log N) amortized)
  winding-number traversal does.

**See also.** [`distance`](#distance),
[`is_point_inside`](#is_point_inside),
[`facet%compute_distance_with_region`](/guide/api/facet-object#compute_distance_with_region).

---

### `query_many`

```fortran
call surface%query_many(points=p, distances=d, closest_points=cp, facets=f, features=ft, normals=n, bodies=b)
```

**Purpose.** Batch rich query for immersed-boundary closures (issue #37):
per point the distance, closest point, closest facet, feature, pseudo-normal
and body. OpenMP over the points; bit-exact with `compute_distance`. See
[Immersed-boundary queries](/guide/advanced/ib-queries).

| Argument             | Intent      | Type                  | Notes                                                                 |
| -------------------- | ----------- | --------------------- | --------------------------------------------------------------------- |
| `points`             | `in`        | `type(vector_R8P)(n)` | Query points.                                                         |
| `distances`          | `out`       | `real(R8P)(n)`        | Distances, signed by default (negative inside), not squared.          |
| `closest_points`     | `out`, opt. | `real(R8P)(3, n)`     | Closest points.                                                       |
| `facets`             | `out`, opt. | `integer(I4P)(n)`     | Closest facets.                                                       |
| `features`           | `out`, opt. | `integer(I4P)(n)`     | `FEATURE_FACE`, `FEATURE_EDGE`, `FEATURE_VERTEX`.                     |
| `feature_indexes`    | `out`, opt. | `integer(I4P)(n)`     | Local edge or vertex index (0 for a face).                            |
| `normals`            | `out`, opt. | `real(R8P)(3, n)`     | Pseudo-normals of the closest features.                               |
| `bodies`             | `out`, opt. | `integer(I4P)(n)`     | Bodies (`facet_body`).                                                |
| `is_signed`          | `in`, opt.  | `logical`             | Default `.true.`.                                                     |
| `sign_algorithm`     | `in`, opt.  | `integer(I4P)`        | Default `SIGN_PSEUDO_NORMAL`.                                         |
| `ambiguous`          | `out`, opt. | `logical(n)`          | Another foot point within the tolerances (CPU only).                  |
| `second_foot_points` | `out`, opt. | `real(R8P)(3, n)`     | The nearest other foot point (the foot point if none).                |
| `rel_tolerance`      | `in`, opt.  | `real(R8P)`           | Default 1e-6.                                                         |
| `foot_tolerance`     | `in`, opt.  | `real(R8P)`           | Default 1e-9 of the box diagonal.                                     |
| `status`             | `out`, opt. | `integer(I4P)`        | `STATUS_INVALID_INPUT` on a size mismatch.                            |
| `pose`               | `in`, opt.  | `rigid_pose_t`        | Body pose: world-frame points and outputs (issue #41).                |
| `curvatures`         | `out`, opt. | `real(R8P)(2, n)`     | Mean and Gaussian curvature at the foot point (issue #42, CPU only).  |

**Related.** `query_point` (one point), `query_many_device` (same, on the
device after `enter_device`; pseudo-normal or fast-winding-number sign),
`probe_many(points, lengths, probe_points, visible, foot_points, normals)`,
`image_points(points, images, visible, foot_points)`,
`get_shells_number`, `get_bodies_number`, `facet_shell(i)`, `facet_body(i)`.

---


### `surface_markers`

```fortran
call surface%surface_markers(dx=0.05_R8P, centroids=c, areas=a, normals=n, spacing_fraction=f)
```

**Purpose.** Grid-adapted surface markers for wall models and force
integration (issue #42). A copy of the surface is remeshed so that the
centroids of neighbouring facets lie about `dx` apart (target edge √3 dx,
plus one calibration pass). The surface itself is not changed. See
[Immersed-boundary queries](/guide/advanced/ib-queries).

| Argument            | Intent      | Type                         | Notes                                                        |
| ------------------- | ----------- | ---------------------------- | ------------------------------------------------------------ |
| `dx`                | `in`        | `real(R8P)`                  | Grid spacing (> 0).                                          |
| `centroids`         | `out`       | `real(R8P), allocatable(3, m)` | Marker positions (facet centroids).                        |
| `areas`             | `out`       | `real(R8P), allocatable(m)`  | Marker areas.                                                |
| `normals`           | `out`       | `real(R8P), allocatable(3, m)` | Unit outward normals.                                      |
| `spacing_fraction`  | `out`, opt. | `real(R8P)`                  | Fraction of nearest-neighbour spacings in [0.8, 1.2] dx.     |
| `spacings`          | `out`, opt. | `real(R8P), allocatable(m)`  | Distance to the nearest edge-neighbour centroid.             |
| `iterations`        | `in`, opt.  | `integer(I4P)`               | Remesh iterations (default 5).                               |
| `preserve_features` | `in`, opt.  | `logical`                    | Keep sharp edges (default `.false.`).                        |
| `status`            | `out`, opt. | `integer(I4P)`               | `STATUS_INVALID_INPUT` for dx <= 0, an empty surface or a failed remesh. |

**Related.** `curve%curve_markers` (exact 2-D resampling), `isotropic_remesh`.
### `move`

```fortran
call surface%move(pose)
```

**Purpose.** Move the surface rigidly by a `rigid_pose_t` and refit the tree
(O(N)) instead of rebuilding it (issue #41): queries afterwards equal those
on the mesh moved and analyzed from scratch.

**Related.** The optional `pose` argument of `distance_many`, `query_many`,
`is_point_inside`, `distance_many_device`, `query_many_device` (body-frame
queries, no move); `state_changes(pose0, pose1, points, fresh, dead,
sign_algorithm, checked, status)` (fresh and dead points between two poses).
See [Moving bodies](/guide/advanced/moving-bodies).

---

### `line_crossings`

```fortran
call surface%line_crossings(axis=1_I4P, a=y0, b=z0, t0=xmin, t1=xmax, t=t, facets=f, entering=e, n=n)
```

**Purpose.** Every crossing of the line along `axis` at (a, b) (the other two
coordinates in cyclic order) with the surface, with parameter in [t0, t1],
sorted, with its facet and whether the line, going in +axis, enters the solid
(issue #40). Exact test (exact orient2d on the projected triangle with a
symbolic tie-break): on a watertight surface every crossing is counted once.

**Related.** `line_crossings_many(axis, a, b, t0, t1, t, facets, entering, counts, overflow)`
and its device twin (fixed buffers, bit-exact);
`cell_wall_fractions(x, y, z, frac, inside, bad_lines)`: the fraction toward each of
the six neighbours of every cell centre and the parity inside flags;
`cut_cells_3d(surface, x, y, z, cut)` (module procedure): face fractions and cell
volumes. See [Cut-cell geometry](/guide/advanced/cut-cells).

---

### `box_classify_many`

```fortran
call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=r)
```

**Purpose.** Classify boxes (AMR blocks or cells) against the surface
(issue #38): `BOX_CUT`, `BOX_NEAR` (within `radius`), `BOX_INSIDE` or
`BOX_OUTSIDE`. SAH BVH traversal with box pruning, OpenMP over the boxes;
the cut test is the separating-axis test on the box inflated by
`box_tolerance`, so it never misses a contact. See
[AMR tagging](/guide/advanced/amr-tagging).

| Argument         | Intent      | Type              | Notes                                                                  |
| ---------------- | ----------- | ----------------- | ---------------------------------------------------------------------- |
| `lo`, `hi`       | `in`        | `real(R8P)(3, n)` | Box corners.                                                           |
| `codes`          | `out`       | `integer(I4P)(n)` | `BOX_*` codes.                                                         |
| `radius`         | `in`, opt.  | `real(R8P)`       | Near radius, default 0 (no `BOX_NEAR`). Up to √3·radius at a corner.   |
| `sign_algorithm` | `in`, opt.  | `integer(I4P)`    | Centre query of uncut boxes, default `SIGN_FAST_WINDING_NUMBER`.       |
| `inside`         | `out`, opt. | `logical(n)`      | Centre sign of every box not cut (also near boxes).                    |
| `status`         | `out`, opt. | `integer(I4P)`    | `STATUS_INVALID_INPUT` on a size mismatch or a negative radius.        |

**Related.** `box_classify_many_device` (same, on the device after
`enter_device`; pseudo-normal or fast-winding-number sign),
`box_intersects_surface(lo, hi)` (one box, cut only),
`tag_block(x, y, z, band, tags, sign_algorithm, status)`: signed tag bands on a
block of cells given its edges `x(nx+1)`, `y(ny+1)`, `z(nz+1)`; `tags(nx, ny, nz)` is
0 on cut cells, else the Chebyshev cell distance to the nearest cut cell
capped at band + 1, negative inside. Optional extra fields: `feature_band`, `feature_tags` (after `build_features`)
and `thin_factor`, `thin_band`, `thin_tags` (after `build_thin`). Box classification and tag bands default
to `SIGN_FAST_WINDING_NUMBER` (issue #39); the point queries keep
`SIGN_PSEUDO_NORMAL`.

---

### `is_point_inside`

```fortran
inside = surface%is_point_inside(point=p)
```

**Purpose.** Return `.true.` if the query point lies strictly inside the
closed surface. A wrapper over [`distance`](#distance) with `is_signed=.true.`
that returns only the sign — slightly more readable at call sites that don't
care about the magnitude.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                              |
| ----------------- | ------------ | ------------------ | ------------------------------------------------------------------ |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                            |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | `SIGN_*`; default `SIGN_PSEUDO_NORMAL`.                            |

**Returns.** `logical` — `.true.` inside, `.false.` outside.

**Example.**

```fortran
program ex_is_point_inside
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
print '(A, L1)', 'centre inside? ', surface%is_point_inside(point=surface%get_centroid())
end program ex_is_point_inside
```

**Notes.**

- Behaviour on the boundary is **algorithm-dependent**: pseudo-normal returns
  whichever side wins the dot product (numerically unstable on the surface
  itself); ray-cast is technically undefined on the surface. Treat the result
  as meaningless when the query point is closer to the surface than the
  mesh's floating-point resolution.
- On open shells, `SIGN_PSEUDO_NORMAL` returns a well-defined value but the
  notion of "inside" is mesh-dependent. Prefer `SIGN_FAST_WINDING_NUMBER` on
  open or holed shells: the winding number degrades smoothly near a gap
  instead of flipping.

**See also.** [`distance`](#distance),
[`compute_distance`](#compute_distance),
[`is_watertight`](#getters-and-predicates).

---

### `winding_number`

```fortran
w = surface%winding_number(point=p, beta=2._R8P)
```

**Purpose.** Generalized / fast winding number (Jacobson 2013 + Barill 2018).
Returns a continuous scalar `w(q)` that measures how many times the surface
wraps around `q`: ≈ 1.0 strictly inside a closed outward-oriented surface,
≈ 0.0 strictly outside, intermediate on the boundary or on **non-watertight
inputs** — the headline "graceful degradation on dirty STL" property that
the legacy sign algorithms cannot deliver. Conceptual narrative and algorithm
overview on the [§1.4 feature page](/guide/advanced/winding-number).

**Arguments.**

| Argument | Intent     | Type               | Notes                                                                                                                                                                |
| -------- | ---------- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `point`  | `in`       | `type(vector_R8P)` | (required) Query point.                                                                                                                                              |
| `beta`   | `in`, opt. | `real(R8P)`        | Barnes-Hut admissibility ratio. Default `2.0`. `beta ≤ 0` forces the exact O(N) per-facet sum (useful for ground truth). Larger β means stricter admissibility, more accuracy, slower. |

**Returns.** `real(R8P)` — the unrounded scalar `w(q)`. Threshold at `0.5`
to use as a binary inside/outside classifier on dirty input.

**Example.**

```fortran
program ex_winding_number_api
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: w_inside, w_outside

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
w_inside  = surface%winding_number(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
w_outside = surface%winding_number(point=vector_R8P(2._R8P, 2._R8P, 2._R8P))
print '(A,F8.4,A,F8.4)', 'inside w = ', w_inside, '   outside w = ', w_outside
endprogram ex_winding_number_api
```

**Notes.**

- Hierarchical traversal uses the existing AABB tree (SAH BVH or octree) —
  no separate data structure built. Per-node dipole moments are computed
  lazily on first query.
- At surface vertices and edges, `w` returns the local solid-angle
  fraction (e.g. `0.125` at a cube corner) — this is the correct
  mathematical answer, not a bug. Threshold at 0.5 for a clean classifier
  or query a small offset away from the surface.
- For typical CAD STL with dirty topology, prefer `winding_number > 0.5`
  over `is_point_inside` — the latter uses sign algorithms that can fail
  silently on non-watertight inputs.

**See also.** [`distance`](#distance), [`is_point_inside`](#is_point_inside),
[§1.4 feature page](/guide/advanced/winding-number).

---

### `intersect_ray_first`

```fortran
call surface%intersect_ray_first(ray_origin, ray_direction, hit, has_hit, status, exclude)
```

**Purpose.** AABB-tree-accelerated **closest hit** ray query (Möller-
Trumbore + best-first traversal with `t_best` pruning). For picking
and "what facet did this ray hit first?" queries. Conceptual overview
plus all three variants on the [§2.5 feature page](/guide/advanced/ray-queries).

**Arguments.**

| Argument        | Intent       | Type               | Notes                                                                                          |
| --------------- | ------------ | ------------------ | ---------------------------------------------------------------------------------------------- |
| `ray_origin`    | `in`         | `type(vector_R8P)` | (required) Ray origin.                                                                         |
| `ray_direction` | `in`         | `type(vector_R8P)` | (required) Ray direction. Need not be unit; `t` is parametric.                                 |
| `hit`           | `out`        | `type(ray_hit_t)`  | (required) Hit record (`%facet_id`, `%t`, `%point`). **Valid only when `has_hit = .true.`**.   |
| `has_hit`       | `out`        | `logical`          | (required) `.true.` if any facet was hit.                                                      |
| `status`        | `out`, opt.  | `integer(I4P)`     | `RAY_STATUS_OK` or `RAY_STATUS_BAD_INPUT`.                                                     |
| `exclude`       | `in`, opt.   | `integer(I4P)`     | A facet to skip: a ray cast from a facet never hits it at t = 0 (issue #62).                   |

**Example.**

```fortran
program ex_intersect_ray_first_api
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: cube
type(ray_hit_t)          :: hit
logical                  :: has_hit
integer(I4P)             :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%intersect_ray_first(ray_origin=vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P), &
                              ray_direction=ex_R8P, hit=hit, has_hit=has_hit, &
                              status=status)
if (has_hit) print '(A,I0,A,F8.4)', 'hit facet ', hit%facet_id, ' at t=', hit%t
endprogram ex_intersect_ray_first_api
```

**See also.**
[§2.5 feature page](/guide/advanced/ray-queries),
[`intersect_ray_all`](#intersect_ray_all),
[`intersect_ray_any`](#intersect_ray_any),
[`facet%intersect_ray`](/guide/api/facet-object).

---

### `intersect_ray_all`

```fortran
call surface%intersect_ray_all(ray_origin, ray_direction, hits, status)
```

**Purpose.** Return **every** ray-facet intersection, sorted ascending
by `t`. For visibility / volume-counting queries (e.g. ray-cast inside
test = odd number of hits).

**Arguments.**

| Argument        | Intent       | Type                          | Notes                                                                                  |
| --------------- | ------------ | ----------------------------- | -------------------------------------------------------------------------------------- |
| `ray_origin`    | `in`         | `type(vector_R8P)`            | (required) Ray origin.                                                                 |
| `ray_direction` | `in`         | `type(vector_R8P)`            | (required) Ray direction.                                                              |
| `hits`          | `out`        | `type(ray_hit_t), allocatable`| (required) Always allocated (zero-length if no hits). Sorted ascending by `t`.         |
| `status`        | `out`, opt.  | `integer(I4P)`                | `RAY_STATUS_OK` or `RAY_STATUS_BAD_INPUT`.                                             |

**Example.**

```fortran
program ex_intersect_ray_all_api
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: cube
type(ray_hit_t), allocatable :: hits(:)
integer(I4P)                 :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%intersect_ray_all(ray_origin=vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P), &
                            ray_direction=ex_R8P, hits=hits, status=status)
print '(A,I0)', 'n_hits = ', size(hits)
endprogram ex_intersect_ray_all_api
```

**See also.**
[§2.5 feature page](/guide/advanced/ray-queries),
[`intersect_ray_first`](#intersect_ray_first),
[`intersect_ray_any`](#intersect_ray_any).

---

### `intersect_ray_any`

```fortran
call surface%intersect_ray_any(ray_origin, ray_direction, max_t, found, status)
```

**Purpose.** Early-exit "is there *any* hit within distance `max_t`?".
For shadow rays, occlusion tests, and conservative inside-test
rejection.

**Arguments.**

| Argument        | Intent       | Type               | Notes                                                                                          |
| --------------- | ------------ | ------------------ | ---------------------------------------------------------------------------------------------- |
| `ray_origin`    | `in`         | `type(vector_R8P)` | (required) Ray origin.                                                                         |
| `ray_direction` | `in`         | `type(vector_R8P)` | (required) Ray direction.                                                                      |
| `max_t`         | `in`, opt.   | `real(R8P)`        | Upper t bound. Default `MaxR8P` (any forward hit). Pass Euclidean distance when direction unit.|
| `found`         | `out`        | `logical`          | (required) `.true.` if at least one hit in `[0, max_t]`.                                       |
| `status`        | `out`, opt.  | `integer(I4P)`     | `RAY_STATUS_OK` or `RAY_STATUS_BAD_INPUT`.                                                     |

**Example.**

```fortran
program ex_intersect_ray_any_api
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: cube
logical                  :: blocked
integer(I4P)             :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%intersect_ray_any(ray_origin=vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P), &
                            ray_direction=ex_R8P, max_t=1.5_R8P, &
                            found=blocked, status=status)
print '(A,L1)', 'blocked? ', blocked
endprogram ex_intersect_ray_any_api
```

**See also.**
[§2.5 feature page](/guide/advanced/ray-queries),
[`intersect_ray_first`](#intersect_ray_first),
[`intersect_ray_all`](#intersect_ray_all).

---

## Manipulation

### After a transform

`translate`, `rotate`, `mirror` and `resize` move the facet vertices and
nothing else. `recompute_metrix=.true.` refreshes the per-facet plane data,
plus the surface bounding box for `translate` and `resize` only. **No
transform updates** the centroid, volume, area, connectivity, pseudo-normals
or AABB tree: getters and queries keep reading the values from before the
transform until [`analyze`](#analyze) runs. Measured on the unit cube:

| Transform | Read before `analyze` | Correct value |
|---|---|---|
| `translate` by 10 along x | signed distance at the new centre −9.5 | −0.5 |
| `rotate` by 90° about z | signed distance at the new centre +0.5 (wrong sign) | −0.5 |
| `resize` by 3 | volume 1, centroid (0.5, 0.5, 0.5) | 27, (1.5, 1.5, 1.5) |

**Rule:** call `analyze` after the last transform of a batch, before any getter
or query. `analyze` recomputes the facet metrix itself, so `recompute_metrix`
is then unnecessary. Pivots read the stored centroid too
(`center=surface%get_centroid()` in `rotate`, `respect_centroid=.true.` in
`resize`): after an earlier transform, call `analyze` first or the pivot is
the old centroid.

**Reflections keep the orientation.** `mirror`, and `resize` with an odd
number of negative factors, are reflections: they reverse the handedness of
each facet's vertex order. FOSSIL swaps two vertices of every facet in that
case (`mirror` by `matrix` does it when `det(matrix) < 0`), so an outward
surface stays outward and no `reverse_normals` is needed.

### `translate`

```fortran
call surface%translate(delta=vec)
call surface%translate(x=1._R8P, y=2._R8P, z=3._R8P)
```

**Purpose.** Rigid-body translation. Pass **either** a `vector_R8P` `delta`
**or** any subset of scalar `x`, `y`, `z` — mixing the two forms returns
`STATUS_AMBIGUOUS_ARGS`.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                              |
| ------------------- | ------------ | ------------------- | ---------------------------------------------------------------------------------- |
| `x`, `y`, `z`       | `in`, opt.   | `real(R8P)`         | Per-axis increments. Missing axes are not translated.                              |
| `delta`             | `in`, opt.   | `type(vector_R8P)`  | Vectorial increment.                                                               |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`. Refreshes facet plane data and bounding box only; see notes.    |
| `status`            | `out`, opt.  | `integer(I4P)`      | `STATUS_OK` or `STATUS_AMBIGUOUS_ARGS`.                                            |

**Example.**

```fortran
program ex_translate
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

type(vector_R8P)         :: c

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%translate(delta=10._R8P * ex_R8P)
call surface%analyze   ! without it, get_centroid returns the pre-translation centroid
c = surface%get_centroid()
print '(A, 3ES12.5)', 'centroid = ', c%x, c%y, c%z
end program ex_translate
```

**Notes.**

- `recompute_metrix=.true.` refreshes the facet plane data and the bounding
  box, not the centroid or the AABB tree: call [`analyze`](#analyze) before
  any getter or query (see [After a transform](#after-a-transform)). For a
  batch of transforms, one `analyze` after the last one is enough.

**See also.** [`rotate`](#rotate),
[`mirror`](#mirror),
[`resize`](#resize).

---

### `rotate`

```fortran
call surface%rotate(axis=ax, angle=theta)
call surface%rotate(matrix=R)
call surface%rotate(axis=ax, angle=theta, center=surface%get_centroid())
```

**Purpose.** Rigid-body rotation. Either give an axis (a `vector_R8P` not
required to be normalised) and an angle in radians, or pass a 3×3 rotation
matrix directly. Pass `center=` to rotate about an arbitrary pivot — by
default the rotation pivots about the world origin.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                                |
| ------------------- | ------------ | ------------------- | ------------------------------------------------------------------------------------ |
| `axis`              | `in`         | `type(vector_R8P)`  | (axis-angle form) Rotation axis.                                                     |
| `angle`             | `in`         | `real(R8P)`         | (axis-angle form) Angle in radians.                                                  |
| `matrix`            | `in`         | `real(R8P)(3,3)`    | (matrix form) Rotation matrix. An improper one (det < 0) keeps the orientation.      |
| `center`            | `in`, opt.   | `type(vector_R8P)`  | Pivot. Defaults to the world origin. Pass `surface%get_centroid()` for body-frame.   |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                                   |

**Example.** Body-frame rotation around the centroid (issue #6).

```fortran
program ex_rotate
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ez_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P),    parameter  :: pi = acos(-1._R8P)

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%rotate(axis=ez_R8P, angle=0.25_R8P * pi, center=surface%get_centroid())
call surface%analyze   ! refresh bounding box, AABB tree and pseudo-normals before queries
end program ex_rotate
```

**Notes.**

- The axis is normalised internally; magnitude is irrelevant.
- World-origin rotation moves the body; body-frame rotation does not. If you
  want "rotate in place", pass `center=surface%get_centroid()`, which reads the
  centroid stored by the last [`analyze`](#analyze).
- Even with `recompute_metrix=.true.` the bounding box, centroid and AABB tree
  stay stale: call `analyze` before any getter or query (see
  [After a transform](#after-a-transform)).

**See also.** [`translate`](#translate),
[`mirror`](#mirror).

---

### `mirror`

```fortran
call surface%mirror(normal=plane_normal)
call surface%mirror(matrix=M)
```

**Purpose.** Mirror the surface through a plane. The plane is specified
either by its unit normal (passing through the origin) or by an explicit 3×3
reflection matrix. The facet winding is reversed along with the reflection, so
an outward-oriented surface stays outward; call [`analyze`](#analyze) before
queries (see [After a transform](#after-a-transform)). A proper (det > 0)
matrix is applied without reversing the winding.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                |
| ------------------- | ------------ | ------------------- | -------------------------------------------------------------------- |
| `normal`            | `in`         | `type(vector_R8P)`  | (normal form) Plane normal. Plane passes through the origin.         |
| `matrix`            | `in`         | `real(R8P)(3,3)`    | (matrix form) Reflection matrix.                                     |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                   |

**Example.**

```fortran
program ex_mirror
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%mirror(normal=ex_R8P)
call surface%analyze   ! rebuild volume, centroid, AABB tree, pseudo-normals
end program ex_mirror
```

**Notes.**

- To mirror through a plane *not* through the origin, translate the surface
  so the plane passes through the origin, mirror, translate back.

**See also.** [`rotate`](#rotate),
[`translate`](#translate).

---

### `resize`

```fortran
call surface%resize(factor=2._R8P * ex_R8P + ey_R8P + ez_R8P)
call surface%resize(x=2._R8P, y=2._R8P, z=2._R8P)
call surface%resize(factor=2._R8P*ex_R8P+ey_R8P+ez_R8P, respect_centroid=.true.)
```

**Purpose.** Anisotropic scaling. As with [`translate`](#translate), pass
either a vector `factor` or scalar `x`/`y`/`z` — never both. By default the
scale is around the world origin; pass `respect_centroid=.true.` to pivot
about the centroid stored by the last [`analyze`](#analyze) (the scale is then
"in place").

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                            |
| ------------------- | ------------ | ------------------- | -------------------------------------------------------------------------------- |
| `x`, `y`, `z`       | `in`, opt.   | `real(R8P)`         | Per-axis factors. Missing axes are not scaled (effective factor 1).              |
| `factor`            | `in`, opt.   | `type(vector_R8P)`  | Vectorial factor.                                                                |
| `respect_centroid`  | `in`, opt.   | `logical`           | Default `.false.`. Pivot at centroid when `.true.`.                              |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                               |
| `status`            | `out`, opt.  | `integer(I4P)`      | `STATUS_OK` or `STATUS_AMBIGUOUS_ARGS`.                                          |

**Example.**

```fortran
program ex_resize
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%resize(factor=2._R8P*ex_R8P + 2._R8P*ey_R8P + 2._R8P*ez_R8P, &
                    respect_centroid=.true.)
call surface%analyze   ! volume, centroid and AABB tree are stale until this call
end program ex_resize
```

**Notes.**

- The method is called `resize` and not `scale` because `scale` is a Fortran
  intrinsic.
- Negative factors mirror across the corresponding axis. With an odd number of
  them the facet winding is reversed as for [`mirror`](#mirror), so the
  orientation is preserved.
- The volume and centroid change with the scale but are only recomputed by
  [`analyze`](#analyze).

**See also.** [`translate`](#translate),
[`mirror`](#mirror).

---

### `clip`

```fortran
call surface%clip(bmin=lo, bmax=hi)
call surface%clip(bmin=lo, bmax=hi, remainder=outside)
```

**Purpose.** Drop facets that do not overlap the AABB `[bmin, bmax]`,
keeping only the part of the surface inside (or intersecting) the box. With
the optional `remainder` argument the discarded facets are moved into a fresh
surface so nothing is lost.

**Arguments.**

| Argument     | Intent       | Type                       | Notes                                                                |
| ------------ | ------------ | -------------------------- | -------------------------------------------------------------------- |
| `bmin`       | `in`         | `type(vector_R8P)`         | (required) Clip box minimum corner.                                  |
| `bmax`       | `in`         | `type(vector_R8P)`         | (required) Clip box maximum corner.                                  |
| `remainder`  | `out`, opt.  | `type(surface_stl_object)` | Receives the facets that were dropped.                               |
| `status`     | `out`, opt.  | `integer(I4P)`             | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                  |

**Example.** Split a surface into two halves on the xy-plane.

```fortran
program ex_clip
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P
implicit none
type(surface_stl_object) :: surface, top
type(vector_R8P)         :: lo, hi

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
lo = surface%get_bmin()
hi = surface%get_bmax()
lo%z = 0.5_R8P * (lo%z + hi%z)
call surface%clip(bmin=lo, bmax=hi, remainder=top)
print '(A,I0,A,I0)', 'lower half facets=', surface%get_facets_number(), &
                     '  upper half facets=', top%get_facets_number()
end program ex_clip
```

**Notes.**

- Clip is **inclusion-by-overlap**, not by-containment: facets that straddle
  the box boundary go into the *inside* surface. There is no triangle
  subdivision along the clipping planes; if you need precise cuts, use a CSG
  boolean operation (planned — see roadmap issue #18).
- Both the kept and remainder surfaces are re-analysed automatically.

**See also.** [`merge_solids`](#merge_solids).

---

### `merge_solids`

```fortran
call surface%merge_solids(other=second_surface)
```

**Purpose.** Append every facet from `other` to `self`. The merged surface is
re-analysed afterwards. **Geometry is not de-duplicated** — overlapping
regions remain as overlapping triangles; the merge is concatenation, not a
boolean union.

**Arguments.**

| Argument | Intent       | Type                       | Notes                                                              |
| -------- | ------------ | -------------------------- | ------------------------------------------------------------------ |
| `other`  | `in`         | `type(surface_stl_object)` | (required) The donor surface. Unchanged on return.                 |
| `status` | `out`, opt.  | `integer(I4P)`             | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                |

**Example.** Reassemble a dragon from two halves.

```fortran
program ex_merge_solids
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: dragon, half2

call dragon%load_from_file(file_name='dragon_part_1.stl', guess_format=.true.)
call half2 %load_from_file(file_name='dragon_part_2.stl', guess_format=.true.)
call dragon%merge_solids(other=half2)
print '(A, I0)', 'merged facets = ', dragon%get_facets_number()
end program ex_merge_solids
```

**Notes.**

- A subsequent [`sanitize`](#sanitize) is recommended: the seam between the
  two donor meshes typically produces near-coincident vertices and duplicate
  facets.

**See also.** [`clip`](#clip),
[`sanitize`](#sanitize).

---

### `boolean`

```fortran
call surface%boolean(other=other_surface, op=BOOL_DIFFERENCE, status=status)
```

**Purpose.** Mesh boolean (Zhou et al. 2016). Replaces `self` in-place
with the result of `op(self, other)`, where `op` is one of `BOOL_UNION`,
`BOOL_INTERSECT`, `BOOL_DIFFERENCE`, `BOOL_SYMDIFF`. Conceptual overview
and full pipeline on the [§1.1 feature page](/guide/advanced/booleans).

**Arguments.**

| Argument | Intent       | Type                        | Notes                                                                                                       |
| -------- | ------------ | --------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `other`  | `in`         | `class(surface_stl_object)` | (required) The right-hand-side solid. Unchanged on return.                                                  |
| `op`     | `in`         | `integer(I4P)`              | (required) `BOOL_UNION`, `BOOL_INTERSECT`, `BOOL_DIFFERENCE`, or `BOOL_SYMDIFF`.                            |
| `status` | `out`, opt.  | `integer(I4P)`              | `BOOL_STATUS_OK`, `BOOL_STATUS_CDT_FAILED`, `BOOL_STATUS_NOT_IMPLEMENTED`, or `BOOL_STATUS_EMPTY_INPUT`.    |

Pre-conditions for clean output: both inputs must be watertight
manifold solids with outward-oriented normals (run `sanitize_normals`
if you're not sure). Both must have their AABB tree built (true after
`load_from_file` / `analyze` / `adopt_facets`).

**Example.**

```fortran
program ex_boolean_api
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: a, b
integer(I4P)             :: status

call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call b%analyze

call a%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
print '(A,I0,A,F8.4)', 'A \ B status=', status, '  vol=', a%get_volume()
endprogram ex_boolean_api
```

**Notes.**

- **Axis-aligned coplanar face overlaps silently produce wrong volumes**
  even with `status == BOOL_STATUS_OK`. Workaround: perturb one input
  by `~1e-6 * bbox_diagonal` along each axis before the boolean. Full
  discussion + the underlying failure modes in the
  [§1.1 known limitations](/guide/advanced/booleans#known-limitations).
- The result replaces `self` in place — if you need the original, copy
  it first (`saved = self`) before calling.
- `BOOL_SYMDIFF` is implemented as
  `(A ∪ B) \ (A ∩ B)` per the standard truth table; performance is
  comparable to UNION + INTERSECT.

**See also.**
[§1.1 feature page](/guide/advanced/booleans),
[`resolve_self_intersections`](#resolve_self_intersections),
[Constants](/guide/api/constants) (BOOL_*, BOOL_STATUS_*).

---

### `find_self_intersections`

```fortran
call surface%find_self_intersections(pairs=pairs, status=status)
```

**Purpose.** Locate every facet pair whose triangles geometrically
cross (Möller tri-tri intersection in the narrow phase, AABB tree-vs-
tree broad phase). Adjacent facets sharing a vertex or edge are
filtered out — only **real geometric crossings** remain. Conceptual
overview on the [§1.2 feature page](/guide/advanced/self-intersection).

**Arguments.**

| Argument | Intent       | Type                                     | Notes                                                                                       |
| -------- | ------------ | ---------------------------------------- | ------------------------------------------------------------------------------------------- |
| `pairs`  | `out`        | `type(intersection_pair_t), allocatable` | (required) Returned list of intersection records. Always allocated (zero-length if none).   |
| `status` | `out`, opt.  | `integer(I4P)`                           | `STATUS_OK` on success.                                                                     |

Each `pairs(i)` record has fields `%a` and `%b` (facet IDs with
`a < b`), and `%p` and `%q` (segment endpoints in world coordinates).

**Example.**

```fortran
program ex_find_self_intersections_api
use fossil
use penf, only : I4P
implicit none
type(surface_stl_object)               :: dragon
type(intersection_pair_t), allocatable :: pairs(:)
integer(I4P)                           :: status

call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%find_self_intersections(pairs=pairs, status=status)
print '(A,I0)', 'self-intersections: ', size(pairs)
endprogram ex_find_self_intersections_api
```

**Notes.**

- Run `sanitize` **first** to filter degenerate-facet noise out of the
  intersection list.
- On clean closed meshes (`cube.stl`, `bunny.stl`), `size(pairs) == 0`.

**See also.**
[§1.2 feature page](/guide/advanced/self-intersection),
[`resolve_self_intersections`](#resolve_self_intersections).

---

### `resolve_self_intersections`

```fortran
call surface%resolve_self_intersections(status=status)
```

**Purpose.** Clean up every self-intersection in-place by running
`boolean(self, self, BOOL_UNION)`. The §1.1 arrangement engine
collects every self-crossing, retriangulates the affected facets, and
the union selection rule keeps only the outer manifold. Result is a
self-intersection-free version of the input.

**Arguments.**

| Argument | Intent       | Type           | Notes                                                              |
| -------- | ------------ | -------------- | ------------------------------------------------------------------ |
| `status` | `out`, opt.  | `integer(I4P)` | `BOOL_STATUS_*` (inherits status codes from the boolean engine).   |

Mutates `self` in place.

**Example.**

```fortran
program ex_resolve_self_intersections_api
use fossil
use penf, only : I4P
implicit none
type(surface_stl_object) :: dragon
integer(I4P)             :: status

call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%resolve_self_intersections(status=status)
print '(A,I0)', 'resolve status = ', status
endprogram ex_resolve_self_intersections_api
```

**Notes.**

- **Inherits the §1.1 coplanar limitation.** Self-intersections that
  involve coplanar facet pairs may produce a result with silently
  wrong volume; perturb by `~1e-6 * bbox_diag` if you suspect
  coplanar self-overlaps.
- Interior cavities formed by self-intersections are **dropped**, not
  re-meshed. If the input had meaningful interior structure, it's
  gone after resolve.
- For inputs too dirty for resolve to handle (large holes, missing
  facets), reach for [`alpha_wrap`](#alpha_wrap) instead.

**See also.**
[§1.2 feature page](/guide/advanced/self-intersection),
[`find_self_intersections`](#find_self_intersections),
[`boolean`](#boolean),
[`alpha_wrap`](#alpha_wrap).

---

### `alpha_wrap`

```fortran
call surface%alpha_wrap(alpha=0.1_R8P, offset=0.033_R8P, wrapped=wrapped, status=status)
```

**Purpose.** Alpha wrapping (Portaneri 2022) — a **watertight,
2-manifold** surface around any input, however broken (holes,
self-intersections, non-manifold patches, missing facets): the boundary
of the space a ball of diameter `alpha` carves from outside, laid at
distance `offset` from the input. Gaps wider than about
`alpha + 2 offset` are kept (both sides wrapped), narrower ones are
capped. Pipeline, measured behaviour and limitations on the
[§1.6 feature page](/guide/advanced/alpha-wrap).

**Arguments.**

| Argument         | Intent       | Type                       | Notes                                                                                                                                  |
| ---------------- | ------------ | -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `alpha`          | `in`, opt.   | `real(R8P)`                | Diameter of the carving ball. Default `bbox_diag / 50` (CGAL's recommendation). The wrap is coarsened to edges of about `alpha/2`.             |
| `offset`         | `in`, opt.   | `real(R8P)`                | Distance of the wrap from the input. Default `alpha / 30`.                                                                             |
| `max_iterations` | `in`, opt.   | `integer(I4P)`             | Grid levels of the refinement loop (cell size `alpha/4`, halved per level). Default 5.                                                 |
| `wrapped`        | `out`        | `type(surface_stl_object)` | (required) Fresh output surface. `self` is unchanged.                                                                                  |
| `status`         | `out`, opt.  | `integer(I4P)`             | `AWRAP_STATUS_OK / _NOT_CONVERGED / _BAD_INPUT / _DEGENERATE`.                                                                         |

**Example.**

```fortran
program ex_alpha_wrap_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: source, wrapped
integer(I4P)             :: status

call source%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call source%alpha_wrap(wrapped=wrapped, status=status)
print '(A,I0,A,L1)', 'wrap status=', status, '  watertight? ', wrapped%is_watertight()
endprogram ex_alpha_wrap_api
```

**Notes.**

- **`AWRAP_STATUS_OK`** means the wrap encloses the input without
  touching it and does not cross itself; **`AWRAP_STATUS_NOT_CONVERGED`**
  that the finest grid level allowed still had a defect there. The wrap
  is closed and 2-manifold in both cases.
- A gap wider than the ball is a feature: the wrap of an open box lines
  both sides of its walls (a thin hollow shell). Check `get_volume()`,
  not only `is_watertight()`.
- **For dirty-input use cases** (the whole point of the algorithm),
  prefer `alpha_wrap` over running `sanitize` followed by `boolean`
  — the latter assumes a watertight input which `sanitize` cannot
  guarantee on truly broken topology.
- The output is a **fresh surface**. `self` is unchanged; copy the
  result into `self` if needed (`self = wrapped`).
- See the [§1.6 known limitations](/guide/advanced/alpha-wrap#known-limitations)
  before using this on inputs with large defects or microscale features.

**See also.**
[§1.6 feature page](/guide/advanced/alpha-wrap),
[`sanitize`](#sanitize) (cheaper local-defect repair),
[Constants](/guide/api/constants) (AWRAP_STATUS_*).

---

### `decimate`

```fortran
call surface%decimate(target_facets=n0/4_I4P, status=status)
```

**Purpose.** Reduce facet count to ≤ `target_facets` via QEM
(Garland & Heckbert 1997) edge-collapse decimation. Quadric error
metric per vertex; min-heap of edge collapse costs; three safety
checks (normal-flip, duplicate-vertex, non-manifold-edge).
Conceptual overview on the [§1.3 feature page](/guide/advanced/decimate).

**Arguments.**

| Argument        | Intent       | Type           | Notes                                                                                            |
| --------------- | ------------ | -------------- | ------------------------------------------------------------------------------------------------ |
| `target_facets` | `in`         | `integer(I4P)` | (required) Upper bound on the output facet count.                                                |
| `status`        | `out`, opt.  | `integer(I4P)` | `DEC_STATUS_OK` (reached target), `_NO_PROGRESS` (couldn't reach), or `_BAD_INPUT`.              |

Mutates `self` in place. After return, the AABB tree, vertex pool,
connectivity, and pseudo-normals are all rebuilt — no further user
action needed before downstream queries.

**Example.**

```fortran
program ex_decimate_api
use fossil
use penf, only : I4P
implicit none
type(surface_stl_object) :: bunny
integer(I4P)             :: n0, status

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
n0 = bunny%get_facets_number()
call bunny%decimate(target_facets=n0 / 4_I4P, status=status)
print '(A,I0,A,I0)', 'after 25%: n_facets = ', bunny%get_facets_number(), &
                     '  status=', status
endprogram ex_decimate_api
```

**Notes.**

- **Volume drifts inward** — see the
  [§1.3 known limitations](/guide/advanced/decimate#known-limitations)
  for the missing volume-preservation term and per-percentage drift
  numbers on a sphere fixture. For volume-critical use cases, run a
  post-decimate volume sanity check via `get_volume()`.
- On perfectly-flat inputs (cube), the no-tie-breaking limitation can
  cause unexpectedly small effective decimation; for general meshes
  this isn't an issue.
- For tessellation-regularity rather than count reduction, reach for
  [`isotropic_remesh`](#isotropic_remesh) instead.

**See also.**
[§1.3 feature page](/guide/advanced/decimate),
[`isotropic_remesh`](#isotropic_remesh),
[Constants](/guide/api/constants) (DEC_STATUS_*).

---

### `resample_via_distance_field`

```fortran
call surface%resample_via_distance_field(resolution=64_I4P, &
                                          surface_out=remeshed, status=status)
```

**Purpose.** "Repair via level set" — sample `self`'s signed distance
field on a regular grid, then extract the iso=0 surface via
[Marching Cubes](/guide/advanced/marching-cubes). Output is a watertight
re-tessellation whose triangle distribution is determined by the grid
spacing rather than by the input's tessellation.

**Arguments.**

| Argument      | Intent       | Type                       | Notes                                                                        |
| ------------- | ------------ | -------------------------- | ---------------------------------------------------------------------------- |
| `resolution`  | `in`         | `integer(I4P)`             | (required) Grid corners along the **longest** bbox axis; cells are cubic.    |
| `surface_out` | `out`        | `type(surface_stl_object)` | (required) Result. Must be **distinct from `self`** — aliasing breaks it.    |
| `status`      | `out`, opt.  | `integer(I4P)`             | `MC_STATUS_*`.                                                               |

`self` is unchanged; the result lives in `surface_out`.

**Example.**

```fortran
program ex_resample_via_distance_field_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: cube_in, cube_out
integer(I4P)             :: status

call cube_in%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_in%resample_via_distance_field(resolution=64_I4P, &
                                          surface_out=cube_out, status=status)
print '(A,F8.4)', 'resampled cube vol = ', cube_out%get_volume()  ! ≈ 1.0
endprogram ex_resample_via_distance_field_api
```

**Notes.**

- **Output triangle count grows as `resolution²`** — 32³ ≈ thousands,
  128³ ≈ hundreds of thousands. Pair with [`decimate`](#decimate) if
  count matters downstream.
- **`surface_out` cannot alias `self`** (the routine consumes
  `self%facet` mid-loop on alias). Always use a distinct destination.
- A 5 % bbox-diagonal margin is added around the source to prevent
  geometry from touching the grid boundary.
- See the [§1.5 limitations](/guide/advanced/marching-cubes#known-limitations)
  for the Lorensen-Cline ambiguity (cases 105 / 150).

**See also.**
[§1.5 feature page](/guide/advanced/marching-cubes),
[`distance`](#distance) (the signed-distance query sampled internally),
[`decimate`](#decimate) (count reduction follow-up).

---

### `segment_sdf`

```fortran
call surface%segment_sdf(facet_labels=labels, sdf=sdf, &
                         num_clusters=4_I4P, status=status)
```

**Purpose.** Per-facet semantic labels via the Shape Diameter Function
(Shapira, Shamir & Cohen-Or 2008). Three-stage pipeline: cone-of-rays
SDF → dual-graph Laplacian smoothing → 1D Gaussian-mixture clustering
(k-means++ init + EM). Conceptual overview on the
[§1.9 feature page](/guide/advanced/sdf-segmentation).

**Arguments.**

| Argument                | Intent       | Type                       | Notes                                                                            |
| ----------------------- | ------------ | -------------------------- | -------------------------------------------------------------------------------- |
| `facet_labels`          | `out`        | `integer(I4P), allocatable`| (required) Per-facet labels in `[0, num_clusters]`; `0` = `SDF_LABEL_UNASSIGNED`.|
| `sdf`                   | `out`, opt.  | `real(R8P), allocatable`   | Smoothed per-facet SDF for inspection.                                           |
| `num_clusters`          | `in`, opt.   | `integer(I4P)`             | Default `4`. Mechanical-CAD-friendly; CGAL uses 5.                               |
| `smoothing_lambda`      | `in`, opt.   | `real(R8P)`                | Default `0.5`. Laplacian blend weight ∈ `[0, 1]`.                                |
| `smoothing_iterations`  | `in`, opt.   | `integer(I4P)`             | Default `2` (Shapira's value).                                                   |
| `num_rays`              | `in`, opt.   | `integer(I4P)`             | Default `30`. Cone rays per facet.                                               |
| `cone_angle_deg`        | `in`, opt.   | `real(R8P)`                | Default `120.0`. Full cone aperture, degrees.                                    |
| `status`                | `out`, opt.  | `integer(I4P)`             | `SDF_STATUS_OK` or `SDF_STATUS_BAD_INPUT`.                                       |

**Example.**

```fortran
program ex_segment_sdf_api
use fossil
use penf, only : I4P
implicit none
type(surface_stl_object)  :: bunny
integer(I4P), allocatable :: labels(:)
integer(I4P)              :: status, k

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny%segment_sdf(facet_labels=labels, num_clusters=4_I4P, status=status)
print '(A,I0,A,I0,A,I0)', 'labels nf=', size(labels), &
      '  range=[', minval(labels), ',', maxval(labels)
do k = 1_I4P, 4_I4P
   print '(A,I0,A,I0)', '  cluster ', k, ': nf = ', count(labels == k)
enddo
endprogram ex_segment_sdf_api
```

**Notes.**

- Cluster ordering is **ascending by SDF mean** — label 1 = thinnest
  features, label `num_clusters` = thickest. Stable across runs
  (k-means++ init is deterministic) for reproducible downstream
  pipelines.
- **GMM is a mode-finder, not a mode-merger** — see the
  [§1.9 known limitations](/guide/advanced/sdf-segmentation#known-limitations)
  for why uniform inputs with `k > 1` may split into near-identical
  clusters. The cube + `k = 2` case is a misleading fixture; use
  sphere + `k = 1` to test trivial cluster behaviour.
- **No graph-cut post-pass** in MVP — labels can flip across
  low-dihedral edges. Deferred to §1.9b.

**See also.**
[§1.9 feature page](/guide/advanced/sdf-segmentation),
[Constants](/guide/api/constants) (SDF_STATUS_*, SDF_LABEL_UNASSIGNED, SDF_SENTINEL).

---

### `isotropic_remesh`

```fortran
call surface%isotropic_remesh(target_length=0.05_R8P, iterations=3_I4P, &
                               preserve_features=.true., status=status)
```

**Purpose.** Botsch & Kobbelt 2004 isotropic remeshing — equalizes edge
length across the surface while preserving its geometry. Each outer
iteration runs four passes: split (edges > 4L/3), collapse (< 4L/5),
flip (valence balancing), tangential relaxation + projection onto the
original surface. Optional sharp-feature preservation via dihedral-angle
vertex lock. Conceptual overview on the
[§1.7 feature page](/guide/advanced/isotropic-remesh).

**Arguments.**

| Argument             | Intent       | Type           | Notes                                                                                                                                  |
| -------------------- | ------------ | -------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `target_length`      | `in`, opt.   | `real(R8P)`    | Target uniform edge length. `≤ 0` (default) → use the median input edge length.                                                        |
| `iterations`         | `in`, opt.   | `integer(I4P)` | Outer-loop count. Default `5`. Convergence is typical at 3-5; more iterations smooth distribution further but also accumulate volume drift. |
| `preserve_features`  | `in`, opt.   | `logical`      | Lock vertices on edges with dihedral > 30°. Default `.false.`. Set to `.true.` for mechanical CAD inputs.                              |
| `status`             | `out`, opt.  | `integer(I4P)` | `REM_STATUS_OK` or `REM_STATUS_BAD_INPUT`.                                                                                             |

Mutates `self` in place.

**Example.**

```fortran
program ex_isotropic_remesh_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
integer(I4P)             :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%isotropic_remesh(target_length=0.2_R8P, iterations=3_I4P, &
                            preserve_features=.true., status=status)
print '(A,I0)', 'after remesh: n_facets = ', cube%get_facets_number()
endprogram ex_isotropic_remesh_api
```

**Notes.**

- **Set `preserve_features=.true.`** for mechanical CAD inputs with
  sharp edges (cubes, prisms, brackets). Without it, the 12 edges of a
  cube round out into a ball after a few iterations.
- **Volume drifts inward on convex surfaces** — see the
  [§1.7 known limitations](/guide/advanced/isotropic-remesh#known-limitations)
  for the inward Laplacian bias and the ~5% / iteration sphere drift.
  Use 2-3 iterations on convex inputs where volume preservation matters.
- A subsequent [`sanitize`](#sanitize) is rarely needed — remeshing
  produces clean topology by construction.

**See also.**
[§1.7 feature page](/guide/advanced/isotropic-remesh),
[`statistics`](#statistics) (verify edge-length distribution before/after).

---

### `smooth`

```fortran
call surface%smooth(method=SMOOTH_METHOD_TAUBIN, status=status)
call surface%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.25_R8P, iterations=1_I4P)
```

**Purpose.** In-place mesh smoothing via either an explicit
uniform-weight Laplacian or Taubin's `λ|μ` non-shrinking alternation.
Defaults to Taubin (band-pass, volume-preserving). The implicit
variant `(M − τL) V_new = M V` is deferred until a sparse SPD solver
is bound. Conceptual overview on the
[§2.3 feature page](/guide/advanced/smoothing).

**Arguments.**

| Argument     | Intent       | Type           | Notes                                                                                                  |
| ------------ | ------------ | -------------- | ------------------------------------------------------------------------------------------------------ |
| `method`     | `in`, opt.   | `integer(I4P)` | `SMOOTH_METHOD_EXPLICIT` or `SMOOTH_METHOD_TAUBIN` (default).                                          |
| `lambda`     | `in`, opt.   | `real(R8P)`    | Positive step. Default `SMOOTH_DEFAULT_LAMBDA = 0.5`.                                                  |
| `mu`         | `in`, opt.   | `real(R8P)`    | Negative counter-step (Taubin only). Default `SMOOTH_DEFAULT_MU = -0.53`. **Must satisfy `mu < -lambda`** or `SMOOTH_STATUS_BAD_INPUT` is returned. |
| `iterations` | `in`, opt.   | `integer(I4P)` | Explicit: step count. Taubin: pair count (each pair = `λ` + `μ`). Default `SMOOTH_DEFAULT_ITERATIONS = 5`. |
| `status`     | `out`, opt.  | `integer(I4P)` | `SMOOTH_STATUS_OK`, `SMOOTH_STATUS_BAD_INPUT`, or `SMOOTH_STATUS_DEGENERATE` (some triangles skipped). |

Mutates `self` in place; AABB tree, vertex pool, connectivity, and
pseudo-normals are all rebuilt before return.

**Example.**

```fortran
program ex_smooth_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: bunny
real(R8P)                :: v_before, v_after
integer(I4P)             :: status

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
v_before = bunny%get_volume()
call bunny%smooth(method=SMOOTH_METHOD_TAUBIN, status=status)
v_after = bunny%get_volume()
print '(A,F6.3,A)', 'volume drift = ', &
                    100._R8P * abs(v_after - v_before) / v_before, ' %'
endprogram ex_smooth_api
```

**Notes.**

- **Prefer Taubin** for production denoising of closed solids — it
  preserves volume to within ~0.5 % on the standard bunny while
  attenuating high-frequency surface noise.
- The uniform-weight Laplacian (not the area-normalised cotangent
  `M⁻¹ L`) is used deliberately; see
  [§2.3 page — Why uniform weights, not cotangent](/guide/advanced/smoothing#why-uniform-weights-not-cotangent).
- After return the surface re-adopts the smoothed facets; downstream
  queries (`distance`, `winding_number`, `mean_curvature`) work
  without further user action.

**See also.**
[§2.3 feature page](/guide/advanced/smoothing),
[`isotropic_remesh`](#isotropic_remesh) (the right primitive when sharp
features must be preserved),
[`mean_curvature`](#mean_curvature) (cleaner curvature fields after
a few Taubin passes).

---

## Discrete differential geometry

Operators and per-vertex scalar fields derived from the surface mesh.
Conceptual coverage on the
[Cotangent Laplacian](/guide/advanced/cotangent-laplacian) and
[Curvature](/guide/advanced/curvature) feature pages.

### `cotangent_laplacian`

```fortran
call surface%cotangent_laplacian(L=L, M=M, status=status)
```

**Purpose.** Builds the sparse cotangent Laplacian `L` and diagonal
barycentric mass `M` on the surface mesh. Both returned as
`csr_matrix_t`. The diagonal post-pass enforces `L_ii = -Σ_{j ≠ i} L_ij`
so the constant-vector kernel `L · 1 = 0` is floating-point exact.
Conceptual overview on the
[§2.1 feature page](/guide/advanced/cotangent-laplacian).

**Arguments.**

| Argument | Intent       | Type                 | Notes                                                                  |
| -------- | ------------ | -------------------- | ---------------------------------------------------------------------- |
| `L`      | `out`        | `type(csr_matrix_t)` | Sparse SPD cotangent Laplacian, `n_vertices × n_vertices`.             |
| `M`      | `out`        | `type(csr_matrix_t)` | Diagonal barycentric mass, `n_vertices × n_vertices`.                  |
| `status` | `out`, opt.  | `integer(I4P)`       | `LAPL_STATUS_OK`, `LAPL_STATUS_BAD_INPUT`, or `LAPL_STATUS_DEGENERATE_TRIANGLE`. |

Does **not** mutate `self`.

**Example.**

```fortran
program ex_cotangent_laplacian_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: sphere
type(csr_matrix_t)       :: L, M
real(R8P), allocatable   :: ones(:), Lones(:)
integer(I4P)             :: n, status

call sphere%load_from_file(file_name='sphere.stl', guess_format=.true.)
call sphere%cotangent_laplacian(L=L, M=M, status=status)
n = L%get_nrows()
allocate(ones(n), Lones(n))
ones = 1._R8P
call L%multiply_vector(x=ones, y=Lones)
print '(A,ES12.5)', 'max |L * 1| = ', maxval(abs(Lones))   ! ~ 0
endprogram ex_cotangent_laplacian_api
```

**See also.**
[§2.1 feature page](/guide/advanced/cotangent-laplacian),
[`csr_matrix_t`](/guide/api/constants),
[`mean_curvature`](#mean_curvature),
[`gaussian_curvature`](#gaussian_curvature).

---

### `gaussian_curvature`

```fortran
call surface%gaussian_curvature(K=K, status=status)
```

**Purpose.** Per-vertex Gaussian curvature `K_i` from the angle defect:
`(2π − Σ θ) / A_i` for interior vertices, `(π − Σ θ) / A_i` for
boundary vertices. `K > 0` on convex caps, `K < 0` on saddles, `K = 0`
on developable surfaces. Integrated `K · A_i` recovers `4π χ` for
closed orientable manifolds (Gauss–Bonnet). Conceptual overview on the
[§2.4 feature page](/guide/advanced/curvature).

**Arguments.**

| Argument | Intent       | Type                   | Notes                                                                                 |
| -------- | ------------ | ---------------------- | ------------------------------------------------------------------------------------- |
| `K`      | `out`        | `real(R8P), alloc(:)`  | Length `n_vertices`. Indexed by pool vertex id.                                       |
| `status` | `out`, opt.  | `integer(I4P)`         | `CURV_STATUS_OK`, `CURV_STATUS_BAD_INPUT`, or `CURV_STATUS_DEGENERATE_TRIANGLE`.      |

Does **not** mutate `self`. Pre-condition: `analyze` has run (every
public construction path guarantees this).

**Example.**

```fortran
program ex_gaussian_curvature_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: sphere
real(R8P), allocatable   :: K(:)
integer(I4P)             :: status

call sphere%load_from_file(file_name='sphere.stl', guess_format=.true.)
call sphere%gaussian_curvature(K=K, status=status)
print '(A,ES12.5)', 'mean K = ', sum(K) / real(size(K, kind=I4P), R8P)
endprogram ex_gaussian_curvature_api
```

**See also.**
[§2.4 feature page](/guide/advanced/curvature),
[`mean_curvature`](#mean_curvature),
[`cotangent_laplacian`](#cotangent_laplacian).

---

### `mean_curvature`

```fortran
call surface%mean_curvature(H=H, status=status)
```

**Purpose.** Per-vertex **signed** mean curvature `H_i` from the
mean-curvature normal identity `H_i n_i = (1/2) M⁻¹ (L V)_i`. Magnitude
is `||H_i n_i||`; sign comes from projecting onto the per-vertex
angle-weighted pseudo-normal. `H > 0` outward-convex, `H < 0`
saddle / inward-concave, `H ≈ 0` flat. Conceptual overview on the
[§2.4 feature page](/guide/advanced/curvature).

**Arguments.**

| Argument | Intent       | Type                   | Notes                                                                                 |
| -------- | ------------ | ---------------------- | ------------------------------------------------------------------------------------- |
| `H`      | `out`        | `real(R8P), alloc(:)`  | Length `n_vertices`. Indexed by pool vertex id.                                       |
| `status` | `out`, opt.  | `integer(I4P)`         | `CURV_STATUS_OK`, `CURV_STATUS_BAD_INPUT`, or `CURV_STATUS_DEGENERATE_TRIANGLE`.      |

Does **not** mutate `self`. Pre-condition: `analyze` has run (every
public construction path guarantees this).

**Example.**

```fortran
program ex_mean_curvature_api
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
real(R8P), allocatable   :: H(:)
integer(I4P)             :: status

call surface%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call surface%mean_curvature(H=H, status=status)
print '(A,I0)',     'curv status = ', status
print '(A,ES12.5)', 'mean H      = ', sum(H) / real(size(H, kind=I4P), R8P)
endprogram ex_mean_curvature_api
```

**Notes.**

- Sign relies on outward-consistent pseudo-normals — run
  [`sanitize_normals`](#sanitize_normals) on dirty inputs first.
- A few Taubin [`smooth`](#smooth) passes before this call produce
  noticeably cleaner `H` fields without volume loss.

**See also.**
[§2.4 feature page](/guide/advanced/curvature),
[`gaussian_curvature`](#gaussian_curvature),
[`jet_curvature`](#jet_curvature),
[`cotangent_laplacian`](#cotangent_laplacian),
[`smooth`](#smooth).

---

### `jet_curvature`

```fortran
call surface%jet_curvature(H=H, K=K, status=status)
```

**Purpose.** Per-vertex mean and Gaussian curvature by a local polynomial
fit (issue #42).
- In the frame of each vertex's pseudo-normal, a cubic height function is
  fitted to the 2-ring by Householder least squares (Cazals & Pouget 2005).
  A quadratic fit is used below 12 neighbours, and the cotangent /
  angle-defect value below 6.
- The error is O(h²) pointwise on any sampling. `mean_curvature` and
  `gaussian_curvature` stall near 15 % at the valence-5 vertices of an
  icosphere.
- Sign as `mean_curvature`: a sphere gives (1/R, 1/R²). Flat
  neighbourhoods give 0 exactly. `query_many(curvatures=)` uses it.

| Argument | Intent      | Type                  | Notes                                               |
| -------- | ----------- | --------------------- | --------------------------------------------------- |
| `H`      | `out`       | `real(R8P), alloc(:)` | Mean curvature, indexed by pool vertex id.          |
| `K`      | `out`       | `real(R8P), alloc(:)` | Gaussian curvature, indexed by pool vertex id.      |
| `status` | `out`, opt. | `integer(I4P)`        | `CURV_STATUS_*` (as `mean_curvature`).              |

---

## Getters and predicates

The read-only accessors are listed here as a single table — they are
self-explanatory and do not benefit from the full per-method template. Every
one of them is `pure`.

| Method                                | Returns                                                                                |
| ------------------------------------- | -------------------------------------------------------------------------------------- |
| `get_facets_number()`                 | `integer(I4P)` — number of facets.                                                     |
| `get_bmin()`, `get_bmax()`            | `type(vector_R8P)` — min / max corners of the axis-aligned bounding box.               |
| `get_volume()`                        | `real(R8P)` — signed volume. Positive for closed bodies with outward normals.          |
| `get_area()`                          | `real(R8P)` — total surface area (sum of facet areas). Added in issue #7.              |
| `get_centroid()`                      | `type(vector_R8P)` — centre of the enclosed volume (from the vertices, as the volume; meaningful on a closed body). Fixed in issue #59: before, it scaled with the facet size. |
| `get_header()`                        | `character(FRLEN)` — STL header string.                                                |
| `set_header(header)`                  | Set the header for the next `save_into_file`.                                          |
| `get_non_manifold_edges_number()`     | `integer(I4P)` — edges shared by 3+ facets.                                            |
| `get_degenerate_facets_removed()`     | `integer(I4P)` — count from the most recent `remove_degenerate_facets` pass.           |
| `get_duplicate_facets_removed()`      | `integer(I4P)` — count from the most recent `remove_duplicate_facets` pass.            |
| `get_vertex_pool()`                   | Read-only pointer to the unique-vertex pool.                                           |
| `facet_at(i)`                         | Pointer to facet `i`; `null()` if `i` is out of range.                                 |
| `facets_ref()`                        | Pointer to the entire facet array.                                                     |

And the three composite validity predicates:

| Predicate          | True iff                                                                                          |
| ------------------ | ------------------------------------------------------------------------------------------------- |
| `is_watertight()`  | Every edge has exactly two incident facets (no boundary, no non-manifold edges).                  |
| `is_manifold()`    | No edge is shared by three or more facets (boundary edges are allowed — open shells qualify).     |
| `is_volume()`      | `is_watertight()` **and** signed volume > 0 **and** centroid is finite.                           |

**Idiomatic usage.** Read the predicates after [`sanitize`](#sanitize) to
classify the input:

```fortran
program ex_predicates
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
if (.not. surface%is_manifold()) then
   print '(A)', 'mesh has non-manifold edges; distance queries may be unreliable'
elseif (.not. surface%is_watertight()) then
   print '(A)', 'mesh is manifold but open; signed distance is mesh-dependent'
elseif (.not. surface%is_volume()) then
   print '(A)', 'mesh is closed but oriented inward; sanitize_normals should have fixed this'
else
   print '(A)', 'mesh is a valid solid'
end if
end program ex_predicates
```

**See also.** [`sanitize`](#sanitize),
[`statistics`](#statistics).
