---
title: Usage
---

# Usage

All examples use `use fossil` which re-exports `surface_stl_object` and `facet_object`, plus the constants for sign algorithms (`SIGN_PSEUDO_NORMAL`, `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`), tree kinds (`AABB_TREE_SAH_BVH`, `AABB_TREE_OCTREE`), and status codes (`STATUS_OK`, `STATUS_INVALID_INPUT`, …).

Numeric kinds (`I4P`, `R8P`) come from [PENF](https://github.com/szaghi/PENF); 3D vectors and unit versors (`ex_R8P`, `ey_R8P`, `ez_R8P`) come from [VecFor](https://github.com/szaghi/VecFor).

## Loading an STL file

### Auto-detect format

Pass `guess_format=.true.` to let FOSSIL determine whether the file is ASCII or binary automatically.

```fortran
use fossil
use penf, only: R8P

type(surface_stl_object) :: surface

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
```

`load_from_file` runs `analyze` internally — bounding box, connectivity, volume, centroid, AABB tree, and pseudo-normals are all populated before it returns. You typically call `sanitize` next to repair common STL defects.

### Explicit format

```fortran
! Explicitly ASCII
call surface%load_from_file(file_name='naca0012.stl', is_ascii=.true.)

! Explicitly binary
call surface%load_from_file(file_name='part.stl', is_ascii=.false.)
```

### Load with on-the-fly clipping

Facets whose vertices all lie outside the bounding box are discarded during the load itself, without a separate clip pass:

```fortran
use fossil
use penf, only: R8P
use vecfor, only: vector_R8P

type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax

bmin%x = -15.0_R8P ; bmin%y = -5.0_R8P ; bmin%z = 0.0_R8P
bmax%x =   0.0_R8P ; bmax%y =  5.0_R8P ; bmax%z = 20.0_R8P

call surface%load_from_file(file_name='dragon.stl', guess_format=.true., &
                            clip_min=bmin, clip_max=bmax)
```

### Rejecting bad input

If any vertex coordinate is `NaN` or `Inf`, `load_from_file` refuses the file via the optional `status` argument rather than letting the bad value propagate:

```fortran
use fossil
use penf, only: I4P
integer(I4P) :: stat

call surface%load_from_file(file_name='maybe_bad.stl', guess_format=.true., status=stat)
if (stat == STATUS_INVALID_INPUT) then
   write(*, '(A)') 'STL contains non-finite coordinates; aborting.'
   stop
endif
```

---

## Printing statistics

```fortran
print '(A)', surface%statistics()
```

Example output (cube fixture, post-sanitize):

```
X extents: [+0.000000000000000E+000, +0.100000000000000E+001]
Y extents: [+0.000000000000000E+000, +0.100000000000000E+001]
Z extents: [+0.000000000000000E+000, +0.100000000000000E+001]
volume: +0.100000000000000E+001
centroid: [+0.500000000000000E+000, +0.500000000000000E+000, +0.500000000000000E+000]
number of facets: +12
number of facets with 1 edges disconnected: +0
number of facets with 2 edges disconnected: +0
number of facets with 3 edges disconnected: +0
number of non-manifold edges (3+ incident facets): +0
degenerate facets removed (last pass): +0
duplicate facets removed (last pass): +0
number of AABB refinement levels: +1
```

::: tip
Volume is positive for a closed body with outward-pointing normals (the standard divergence-theorem convention; `sanitize` produces this state). If you see a negative volume after sanitize, the mesh is probably non-watertight — check the non-manifold-edge and disconnected-edge counts.
:::

---

## Saving an STL file

```fortran
! Save as ASCII
call surface%save_into_file(file_name='output.stl', is_ascii=.true.)

! Save as binary (default)
call surface%save_into_file(file_name='output.stl')
```

---

## Sanitize (the full repair pipeline)

`sanitize` is the orchestrator. It runs every repair pass in the right order and prints one stderr warning per defect class detected.

```fortran
use fossil

type(surface_stl_object) :: surface

call surface%load_from_file(file_name='messy.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
```

The passes that run inside `sanitize`:

1. `remove_degenerate_facets` — drop zero-area / sliver triangles before their NaN normals contaminate downstream queries.
2. `connect_nearby_vertices` (if disconnected edges were detected) — union-find vertex deduplication on a tolerance.
3. `analyze` — rebuild metrix, connectivity, vertex occurrences.
4. `remove_duplicate_facets` — drop literal duplicates, including reversed-orientation ones.
5. `analyze` — rebuild after duplicate removal (only if any duplicates were dropped).
6. `sanitize_normals` — BFS-propagate winding consistency, then global flip so volume > 0.

Each individual pass is also callable directly if you want fine-grained control.

---

## Validity predicates

After `sanitize`, three named predicates summarise the mesh state:

```fortran
if (.not. surface%is_watertight()) print '(A)', 'surface has boundary or non-manifold edges'
if (.not. surface%is_manifold())   print '(A)', 'surface has non-manifold edges (3+ incidences)'
if (.not. surface%is_volume()) then
   print '(A)', 'volume is not physically meaningful — do not trust get_volume()/get_centroid()'
endif
```

| Predicate | True iff |
|---|---|
| `is_watertight` | Every edge has exactly 2 incident facets |
| `is_manifold`   | No non-manifold edges (boundary edges allowed) |
| `is_volume`     | watertight AND volume > 0 AND finite centroid |

---

## Translate

Move the surface by a vectorial delta or by individual axis components:

```fortran
use fossil
use vecfor, only: vector_R8P

type(vector_R8P) :: delta

delta%x = 1.0_R8P ; delta%y = 0.0_R8P ; delta%z = 0.0_R8P

! Translate by a 3D vector
call surface%translate(delta=delta)

! Translate by scalar components (any combination)
call surface%translate(x=2.0_R8P)
call surface%translate(y=1.5_R8P)
call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
```

---

## Rotate

Rotate around an arbitrary axis by an angle in radians:

```fortran
use fossil
use vecfor, only: ex_R8P, vector_R8P
use penf, only: R8P

type(vector_R8P) :: axis
real(R8P)        :: angle

axis  = ex_R8P        ! rotate around x-axis
angle = 1.5707963_R8P ! π/2 radians

call surface%rotate(axis=axis, angle=angle)
```

Alternatively, provide a pre-built rotation matrix:

```fortran
use vecfor, only: rotation_matrix_R8P, vector_R8P

real(R8P) :: matrix(3,3)
matrix = rotation_matrix_R8P(axis=ex_R8P, angle=1.5707963_R8P)
call surface%rotate(matrix=matrix)
```

---

## Mirror

Mirror with respect to a plane defined by its outward normal:

```fortran
use fossil
use vecfor, only: ex_R8P, ey_R8P, vector_R8P

type(vector_R8P) :: normal

! Mirror with respect to the YZ-plane (normal = x-axis)
normal = ex_R8P
call surface%mirror(normal=normal)

! Mirror with respect to an oblique plane
normal = ex_R8P + ey_R8P
call surface%mirror(normal=normal)
```

Or by a pre-built mirror matrix:

```fortran
use vecfor, only: mirror_matrix_R8P

real(R8P) :: matrix(3,3)
matrix = mirror_matrix_R8P(normal=ex_R8P)
call surface%mirror(matrix=matrix)
```

---

## Resize (scale)

Scale all facets by a uniform or per-axis factor. The default pivot is the origin; pass `respect_centroid=.true.` to scale about the surface centroid instead:

```fortran
use fossil
use vecfor, only: ex_R8P, ey_R8P, ez_R8P, vector_R8P
use penf, only: R8P

type(vector_R8P) :: factor

! Uniform scale by a 3D vector factor
factor = 2.0_R8P * ex_R8P + 2.0_R8P * ey_R8P + 2.0_R8P * ez_R8P
call surface%resize(factor=factor)

! Per-axis scalar scale
call surface%resize(x=0.5_R8P, z=1.2_R8P)

! Scale about the centroid
call surface%resize(factor=factor, respect_centroid=.true.)
```

---

## Clip

Discard all facets outside a given axis-aligned bounding box. The removed facets can be captured as a separate surface:

```fortran
use fossil
use vecfor, only: vector_R8P

type(surface_stl_object) :: surface, remainder
type(vector_R8P)         :: bmin, bmax

bmin%x = -15.0_R8P ; bmin%y = -5.0_R8P ; bmin%z = 0.0_R8P
bmax%x =   0.0_R8P ; bmax%y =  5.0_R8P ; bmax%z = 20.0_R8P

call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)
! `clip` re-runs analyze internally on both `surface` and `remainder`.
```

---

## Merge

Combine two surfaces into one. The result is stored in the first surface:

```fortran
use fossil

type(surface_stl_object) :: surface_1, surface_2

call surface_1%load_from_file(file_name='dragon_part_1.stl', guess_format=.true.)
call surface_2%load_from_file(file_name='dragon_part_2.stl', guess_format=.true.)

call surface_1%merge_solids(other=surface_2)
call surface_1%analyze        ! rebuild connectivity / tree after the merge
call surface_1%save_into_file(file_name='dragon_merged.stl')
```

---

## Distance queries

### Unsigned distance (squared, default)

```fortran
use fossil
use vecfor, only: ex_R8P, ey_R8P, ez_R8P
use penf, only: R8P

real(R8P) :: d2

d2 = surface%distance(point=2.0_R8P * ex_R8P + 0.0_R8P * ey_R8P + 0.0_R8P * ez_R8P)
```

### Euclidean (non-squared) distance

```fortran
d = surface%distance(point=p, is_square_root=.true.)
```

### Signed distance

Negative inside, positive outside. Default sign algorithm is `SIGN_PSEUDO_NORMAL` (Bærentzen–Aanæs); pass `sign_algorithm=` to override.

```fortran
real(R8P) :: sd
sd = surface%distance(point=p, is_signed=.true., is_square_root=.true.)

! Override the sign algorithm for benchmarking or for meshes whose orientation
! cannot be sanitised (e.g. open shells where outward is ill-defined).
sd = surface%distance(point=p, is_signed=.true., is_square_root=.true., &
                      sign_algorithm=SIGN_RAY_INTERSECTIONS)
```

### Point-in-polyhedron test

```fortran
logical :: inside

inside = surface%is_point_inside(point=0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.5_R8P * ez_R8P)
if (inside) print *, 'point is inside the surface'
```

### Choosing the AABB tree kind

The default `AABB_TREE_SAH_BVH` is faster on every measured mesh (see [Features](./features.md#aabb-acceleration)). The legacy octree is still selectable:

```fortran
call surface%load_from_file(file_name='dragon.stl', guess_format=.true., &
                            aabb_tree_kind=AABB_TREE_OCTREE)
```

Both kinds produce bit-exact distance results — the difference is purely query speed.

---

## Connect nearby vertices (repair)

`sanitize` calls this internally when needed; you only invoke it directly if you want to repair a surface without running the full pipeline:

```fortran
use fossil

type(surface_stl_object) :: surface

call surface%load_from_file(file_name='disconnected.stl', guess_format=.true.)
print '(A)', surface%statistics()   ! shows disconnected edge counts

call surface%connect_nearby_vertices
call surface%analyze
print '(A)', surface%statistics()   ! disconnected edges should be reduced
```

---

## Full workflow example

```fortran
use fossil
use penf, only: R8P
use vecfor, only: ex_R8P, ey_R8P, ez_R8P

type(surface_stl_object) :: surface

! 1. Load (analyze runs internally).
call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)

! 2. Run the full repair pipeline.
call surface%sanitize
print '(A)', surface%statistics()

! 3. Center at origin and rescale to fit in a unit box.
associate(c => surface%get_centroid(), bmin => surface%get_bmin(), bmax => surface%get_bmax())
   call surface%translate(delta=-c)
   call surface%resize(factor=(1.0_R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])) * &
                              (ex_R8P + ey_R8P + ez_R8P))
end associate
call surface%analyze

! 4. Save.
call surface%save_into_file(file_name='dragon-normalized.stl')
```
