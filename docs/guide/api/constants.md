---
title: Constants
---

# Constants

This page documents the named integer constants re-exported by `use fossil`.
They are passed as optional arguments to control the behaviour of three things:

1. **Signed-distance queries** — the `sign_algorithm` argument on `distance`,
   `compute_distance`, and `is_point_inside` picks one of the `SIGN_*` codes.
2. **The AABB acceleration structure** — the `aabb_tree_kind` argument on
   `load_from_file` and `analyze` picks one of `AABB_TREE_SAH_BVH` /
   `AABB_TREE_OCTREE`; `aabb_refinement_levels` accepts `AABB_AUTO_REFINEMENT`
   to delegate the octree depth to the auto-tune heuristic.
3. **Error reporting** — the optional `status` argument on every mutating
   procedure returns one of the `STATUS_*` codes.

Each constant is a small `integer(I4P)`. You never need to know its numeric
value — refer to it by name and you stay portable.

---

## Sign algorithms

The sign of a signed-distance query (whether the query point lies inside or
outside the closed surface) is decided by one of three algorithms. The default
is `SIGN_PSEUDO_NORMAL`; the other two are kept for benchmarking and as
fallbacks for meshes whose orientation cannot be sanitised.

### `SIGN_PSEUDO_NORMAL` *(default)*

**Purpose.** Bærentzen–Aanæs angle-weighted pseudo-normal test. At the closest
point on the surface, the sign of `dot(point - closest, pseudo_normal)`
classifies the query point as inside or outside. The pseudo-normal is fused
into the existing distance traversal so signed distance costs the same as
unsigned distance.

**When to use.** Almost always. This is the default for a reason: it is the
fastest of the three, robust on closed meshes, and produces a continuous sign
field (no flickering on flat regions, no ray-grazing artefacts).

**Precondition.** Normals must be consistently outward-oriented. `sanitize`
produces this state; `sanitize_normals` alone is enough if you already trust
the rest of the mesh.

**Example.**

```fortran
program ex_sign_pseudo_normal
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: d

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
d = surface%distance(point=2._R8P * ex_R8P, is_signed=.true., &
                     sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
print '(A, ES12.5)', 'signed distance (pseudo-normal) = ', d
end program ex_sign_pseudo_normal
```

**See also.** [`distance`](/guide/api/surface-stl-object#distance),
[`is_point_inside`](/guide/api/surface-stl-object#is_point_inside),
[`sanitize_normals`](/guide/api/surface-stl-object#sanitize_normals).

**Reference.** Bærentzen & Aanæs, *Signed Distance Computation Using the Angle
Weighted Pseudo-normal*, IEEE TVCG 11(3), 2005.

---

### `SIGN_RAY_INTERSECTIONS`

**Purpose.** Axis-aligned ray cast from the query point: count surface
intersections along the +x axis; odd means inside, even means outside. The sign
is decided *after* the unsigned distance has been computed (two passes total).

**When to use.** Meshes with mixed normal orientation that you cannot or do
not want to repair (open shells, scan fragments). Also useful as a regression
oracle against `SIGN_PSEUDO_NORMAL` on hand-crafted test surfaces.

**Cost.** ~3× slower than `SIGN_PSEUDO_NORMAL` because of the second pass and
the per-facet ray-triangle test. Ray-grazing tie-breaks are handled but
introduce a small failure rate near edge-on geometry.

**Example.**

```fortran
program ex_sign_ray
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: d

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
d = surface%distance(point=2._R8P * ex_R8P, is_signed=.true., &
                     sign_algorithm=SIGN_RAY_INTERSECTIONS, is_square_root=.true.)
print '(A, ES12.5)', 'signed distance (ray) = ', d
end program ex_sign_ray
```

**See also.** [`SIGN_PSEUDO_NORMAL`](#sign_pseudo_normal),
[`SIGN_SOLID_ANGLE`](#sign_solid_angle).

---

### `SIGN_SOLID_ANGLE`

**Purpose.** Sum the projected solid angles subtended by each facet from the
query point; the sum approaches ±4π inside and 0 outside.

**When to use.** Reference / regression. **Not accelerated** — every query is
O(N) in the facet count regardless of the AABB tree state. Use on small
surfaces or low-frequency queries only.

**Example.**

```fortran
program ex_sign_solid_angle
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: d

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
d = surface%distance(point=2._R8P * ex_R8P, is_signed=.true., &
                     sign_algorithm=SIGN_SOLID_ANGLE, is_square_root=.true.)
print '(A, ES12.5)', 'signed distance (solid angle) = ', d
end program ex_sign_solid_angle
```

**See also.** [`SIGN_PSEUDO_NORMAL`](#sign_pseudo_normal),
[`facet%solid_angle`](/guide/api/facet-object#solid_angle).

---

### `sign_algorithm_from_string`

Helper for CLI / config parsing. Maps a lower-case name to the matching
`SIGN_*` constant, or raises `error stop` on an unknown string.

```fortran
algo = sign_algorithm_from_string('pseudo_normal')      ! == SIGN_PSEUDO_NORMAL
algo = sign_algorithm_from_string('ray_intersections')  ! == SIGN_RAY_INTERSECTIONS
algo = sign_algorithm_from_string('solid_angle')        ! == SIGN_SOLID_ANGLE
```

---

## AABB tree kind

Each `surface_stl_object` owns one acceleration structure for distance and
inside-test queries. The kind is chosen at load / analyse time and cannot be
changed without rebuilding (which is what `set_tree_kind` followed by
`analyze` does for you).

### `AABB_TREE_SAH_BVH` *(default)*

**Purpose.** Binary bounding-volume hierarchy built top-down by partitioning
triangles along the longest axis of their centroid bounding box, using the
bucketed surface-area heuristic (16 buckets per axis). The tree adapts to
local triangle density: leaves cluster tightly where the surface is dense and
spread thinly where it is not.

**When to use.** Almost always. On `dragon-fine.stl` (24 k facets, 32³ query
grid) it is ~110× faster than the octree. Build cost is O(N log N).

**Note.** Ignores the `aabb_refinement_levels` argument — depth is a function
of the data, not of a user knob.

**See also.** [`AABB_TREE_OCTREE`](#aabb_tree_octree),
[`AABB_AUTO_REFINEMENT`](#aabb_auto_refinement),
[`aabb_tree_object`](/guide/api/aabb-tree-object).

---

### `AABB_TREE_OCTREE`

**Purpose.** Eight-way space-partitioning octree. Each non-leaf node is split
into eight axis-aligned children regardless of triangle distribution. Depth is
fixed up front via `aabb_refinement_levels` (or `AABB_AUTO_REFINEMENT`).

**When to use.** Benchmarking the SAH BVH; rare meshes where the BVH's
top-down split happens to perform poorly; teaching / debugging because the
geometry is easier to visualise. **Not the default** — the BVH wins on every
realistic input we have measured.

**Example.** (switching kind on a freshly-loaded surface)

```fortran
program ex_octree
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true., &
                            aabb_tree_kind=AABB_TREE_OCTREE,            &
                            aabb_refinement_levels=AABB_AUTO_REFINEMENT)
print '(A, I0)', 'octree refinement levels = ', surface%aabb%get_refinement_levels()
end program ex_octree
```

**See also.** [`AABB_TREE_SAH_BVH`](#aabb_tree_sah_bvh),
[`aabb%set_tree_kind`](/guide/api/aabb-tree-object#set_tree_kind).

---

### `AABB_AUTO_REFINEMENT`

**Purpose.** Sentinel passed in place of an integer depth: ask the octree to
pick a depth from the facet count, `ceil(log8(N/64))`, clamped to `[1, 6]`.

**When to use.** Whenever you choose `AABB_TREE_OCTREE` and have no specific
depth in mind. Has no effect when `AABB_TREE_SAH_BVH` is active (the BVH
ignores it).

---

## Status codes

Every mutating procedure on `surface_stl_object` and a few on `facet_object`
accepts an optional `intent(out)` `status` argument. When you supply it, the
procedure **never** calls `error stop` — failures are reported through the
status code instead. When you omit it, the procedure aborts on unrecoverable
errors.

| Constant                 | Meaning                                                                                              |
| ------------------------ | ---------------------------------------------------------------------------------------------------- |
| `STATUS_OK`              | Success. Numeric value `0`.                                                                          |
| `STATUS_ALLOC_FAIL`      | An internal `allocate` failed (usually out of memory).                                               |
| `STATUS_AMBIGUOUS_ARGS`  | Conflicting optionals — e.g. both `delta=` and `x=/y=/z=` passed to `translate`.                     |
| `STATUS_FILE_NOT_FOUND`  | `load_from_file` could not open the input path for reading.                                          |
| `STATUS_FILE_OPEN_FAIL`  | `save_into_file` could not open the output path for writing.                                         |
| `STATUS_INVALID_INPUT`   | The input file contains NaN / Inf vertex coordinates and was refused. The mesh is left untouched.    |

**Idiom for fail-soft loading.**

```fortran
program ex_status
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: stat

call surface%load_from_file(file_name='maybe-missing.stl', guess_format=.true., status=stat)
select case (stat)
case (STATUS_OK)
   print '(A,I0,A)', 'loaded ', surface%get_facets_number(), ' facets'
case (STATUS_FILE_NOT_FOUND)
   print '(A)', 'file not found'
case (STATUS_INVALID_INPUT)
   print '(A)', 'file contains NaN or Inf coords; load refused'
case default
   print '(A,I0)', 'load failed with status ', stat
end select
end program ex_status
```

**See also.** [`load_from_file`](/guide/api/surface-stl-object#load_from_file),
[`save_into_file`](/guide/api/surface-stl-object#save_into_file),
[`translate`](/guide/api/surface-stl-object#translate),
[`resize`](/guide/api/surface-stl-object#resize).

---

## Smoothing — methods, defaults, status codes

Constants consumed by [`surface%smooth`](/guide/api/surface-stl-object#smooth)
and the underlying `fossil_smoothing` module. Conceptual overview on the
[§2.3 feature page](/guide/advanced/smoothing).

### Methods

| Constant                  | Meaning                                                                                                                  |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `SMOOTH_METHOD_EXPLICIT`  | One `V_new = V + λ ΔV` step per iteration. Shrinks. Use only as a building block.                                        |
| `SMOOTH_METHOD_TAUBIN`    | One `(λ, μ)` pair per iteration (Taubin 1995 band-pass). Volume-preserving. **The default**; pick this for production.   |

### Defaults

| Constant                       | Value      | Used by                            |
| ------------------------------ | ---------- | ---------------------------------- |
| `SMOOTH_DEFAULT_LAMBDA`        | `0.5`      | `lambda` default for both methods. |
| `SMOOTH_DEFAULT_MU`            | `-0.53`    | `mu` default for Taubin.           |
| `SMOOTH_DEFAULT_ITERATIONS`    | `5`        | Step count (explicit) or pair count (Taubin). |

### Status codes

| Constant                    | Meaning                                                                                                                   |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `SMOOTH_STATUS_OK`          | Smoothing completed normally.                                                                                             |
| `SMOOTH_STATUS_BAD_INPUT`   | Empty surface, `iterations ≤ 0`, unknown `method`, or Taubin called with `μ ≥ −λ` (band-pass property violated).          |
| `SMOOTH_STATUS_DEGENERATE`  | At least one triangle had area below tolerance during Laplacian assembly; iteration proceeded on the rest.                |

**See also.** [`surface%smooth`](/guide/api/surface-stl-object#smooth),
[§2.3 feature page](/guide/advanced/smoothing).

---

## Curvature / Laplacian status codes

Constants consumed by
[`surface%cotangent_laplacian`](/guide/api/surface-stl-object#cotangent_laplacian),
[`surface%gaussian_curvature`](/guide/api/surface-stl-object#gaussian_curvature),
and [`surface%mean_curvature`](/guide/api/surface-stl-object#mean_curvature).
Conceptual overview on the
[Cotangent Laplacian](/guide/advanced/cotangent-laplacian) and
[Curvature](/guide/advanced/curvature) feature pages.

### Cotangent Laplacian

| Constant                              | Meaning                                                                                            |
| ------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `LAPL_STATUS_OK`                      | Operator built successfully.                                                                       |
| `LAPL_STATUS_BAD_INPUT`               | Empty surface or uninitialised vertex pool.                                                        |
| `LAPL_STATUS_DEGENERATE_TRIANGLE`     | At least one triangle had area below tolerance and was skipped; operator returned on the survivors.|

### Curvature (Gaussian + mean, shared code set)

| Constant                              | Meaning                                                                                                |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `CURV_STATUS_OK`                      | Curvature field built successfully.                                                                    |
| `CURV_STATUS_BAD_INPUT`               | Empty surface or uninitialised pool.                                                                   |
| `CURV_STATUS_DEGENERATE_TRIANGLE`     | At least one triangle had area below tolerance and was skipped; degenerate vertices receive 0.         |

---

## `csr_matrix_t` — sparse-matrix container

A minimal compressed-sparse-row (CSR) container exported by
`use fossil`. Returned by
[`surface%cotangent_laplacian`](/guide/api/surface-stl-object#cotangent_laplacian).
Designed for read-only use after construction (matrix-vector product +
element / structure inspection); not a general linear-algebra library.

**Member procedures.**

| Procedure                  | Signature                                                                                          | Notes                                                                |
| -------------------------- | -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| `get_nrows()`              | `pure function get_nrows() result(n)` — `integer(I4P)`                                             | Row count.                                                           |
| `get_ncols()`              | `pure function get_ncols() result(n)` — `integer(I4P)`                                             | Column count. For the operators returned by §2.1 this equals `get_nrows()`. |
| `get_nnz()`                | `pure function get_nnz() result(nz)` — `integer(I4P)`                                              | Number of structurally stored non-zeros.                             |
| `get_value(i, j, status)`  | `function get_value(i, j, status) result(v)` — `real(R8P)`                                         | Element access; returns `0` for non-stored entries. Optional `status` for bounds violations. |
| `multiply_vector(x, y)`    | `subroutine multiply_vector(x, y)` — `real(R8P), intent(in) :: x(:)` / `real(R8P), intent(out) :: y(:)` | Computes `y = A · x`.                                                |
| `row_sum(i)`               | `function row_sum(i) result(s)` — `real(R8P)`                                                      | Sum of all stored entries on row `i`. Used internally for the `L_ii = -Σ_j L_ij` post-pass. |
| `is_symmetric()`           | `function is_symmetric() result(sym)` — `logical`                                                  | True iff `A_ij == A_ji` to within a small tolerance — useful as a sanity check on the cotangent Laplacian. |

**Status codes** (currently informational — the construction paths
shipped today always succeed):

| Constant                       | Meaning                                                       |
| ------------------------------ | ------------------------------------------------------------- |
| `CSR_STATUS_OK`                | Operation completed.                                          |
| `CSR_STATUS_BAD_INPUT`         | Constructor / mutator received an inconsistent argument.      |
| `CSR_STATUS_OUT_OF_RANGE`      | `get(i, j)` / `matvec` called with out-of-range index.        |

**See also.**
[`surface%cotangent_laplacian`](/guide/api/surface-stl-object#cotangent_laplacian),
[§2.1 feature page](/guide/advanced/cotangent-laplacian),
[§2.4 feature page](/guide/advanced/curvature).
