---
title: Usage
---

# Usage

All examples use `use fossil` which re-exports `file_stl_object`, `surface_stl_object`, and `facet_object`. Numeric kinds (`I4P`, `R8P`) come from [PENF](https://github.com/szaghi/PENF); 3D vectors and unit versors (`ex_R8P`, `ey_R8P`, `ez_R8P`) come from [VecFor](https://github.com/szaghi/VecFor).

## Loading an STL file

### Auto-detect format

Pass `guess_format=.true.` to let FOSSIL determine whether the file is ASCII or binary automatically.

```fortran
use fossil
use penf, only: R8P

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

call file_stl%load_from_file(facet=surface%facet, file_name='dragon.stl', guess_format=.true.)
call surface%analize
```

### Explicit format

```fortran
! Explicitly ASCII
call file_stl%load_from_file(facet=surface%facet, file_name='naca0012.stl', is_ascii=.true.)

! Explicitly binary
call file_stl%load_from_file(facet=surface%facet, file_name='part.stl', is_ascii=.false.)
```

### Load with on-the-fly clipping

Facets whose centroid lies outside the bounding box are discarded during the load itself, without a separate clip pass:

```fortran
use fossil
use vecfor, only: vector_R8P

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax

bmin%x = -15.0_R8P ; bmin%y = -5.0_R8P ; bmin%z = 0.0_R8P
bmax%x =   0.0_R8P ; bmax%y =  5.0_R8P ; bmax%z = 20.0_R8P

call file_stl%load_from_file(facet=surface%facet, file_name='dragon.stl', &
                             guess_format=.true., clip_min=bmin, clip_max=bmax)
call surface%analize
```

---

## Printing statistics

After `analize`, both the file handler and the surface can print their statistics:

```fortran
print '(A)', file_stl%statistics()
print '(A)', surface%statistics()
```

Example output:

```
Mesh_1
file name:   src/tests/cube.stl
file format: ascii
X extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
Y extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
Z extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
volume: -0.100000000000000E+001
centroid: [+0.500000000000000E+000, +0.500000000000000E+000, +0.500000000000000E+000]
number of facets: +12
number of facets with 1 edges disconnected: +0
number of facets with 2 edges disconnected: +0
number of facets with 3 edges disconnected: +0
number of AABB refinement levels: +2
```

::: tip
Volume is negative when normals point inward. Use `sanitize_normals` to make them consistent, then recompute.
:::

---

## Saving an STL file

```fortran
! Save as ASCII
call file_stl%save_into_file(facet=surface%facet, file_name='output.stl', is_ascii=.true.)

! Save as binary (default)
call file_stl%save_into_file(facet=surface%facet, file_name='output.stl')
```

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

axis  = ex_R8P       ! rotate around x-axis
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
call surface%analize
call remainder%analize
```

---

## Merge

Combine two surfaces into one. The result is stored in the first surface:

```fortran
use fossil

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface_1, surface_2

call file_stl%load_from_file(facet=surface_1%facet, file_name='dragon_part_1.stl', guess_format=.true.)
call file_stl%load_from_file(facet=surface_2%facet, file_name='dragon_part_2.stl', guess_format=.true.)
call surface_1%analize
call surface_2%analize

call surface_1%merge_solids(other=surface_2)

call file_stl%save_into_file(facet=surface_1%facet, file_name='dragon_merged.stl')
```

---

## Sanitize normals

Make all facet normals consistent (all outward). This is typically needed after loading STL files that have mixed or inward-pointing normals:

```fortran
use fossil

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

call file_stl%load_from_file(facet=surface%facet, file_name='cube-inconsistent.stl', guess_format=.true.)
call surface%analize

print *, 'volume before:', surface%volume   ! likely negative or wrong

call surface%sanitize_normals
call surface%compute_volume

print *, 'volume after: ', surface%volume   ! positive and correct
```

---

## Distance queries

### Minimum distance from a point

```fortran
use fossil
use vecfor, only: ex_R8P, ey_R8P, ez_R8P
use penf, only: R8P

real(R8P) :: d

call surface%analize
d = surface%distance(point=2.0_R8P * ex_R8P + 0.0_R8P * ey_R8P + 0.0_R8P * ez_R8P)
print *, 'distance =', d
```

### Point-in-polyhedron test

```fortran
logical :: inside

inside = surface%is_point_inside(point=0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.5_R8P * ez_R8P)
if (inside) print *, 'point is inside the surface'
```

---

## Connect nearby vertices (repair)

Repair a surface with disconnected edges by snapping vertices that are geometrically close together:

```fortran
use fossil

type(surface_stl_object) :: surface

call file_stl%load_from_file(facet=surface%facet, file_name='disconnected.stl', guess_format=.true.)
call surface%analize
print '(A)', surface%statistics()   ! shows disconnected edge counts

call surface%connect_nearby_vertices
call surface%analize
print '(A)', surface%statistics()   ! disconnected edges should be reduced
```

---

## Full workflow example

```fortran
use fossil
use penf, only: R8P
use vecfor, only: ex_R8P, ey_R8P, ez_R8P, vector_R8P

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

! 1. Load
call file_stl%load_from_file(facet=surface%facet, file_name='src/tests/dragon.stl', guess_format=.true.)
call surface%analize
print '(A)', surface%statistics()

! 2. Fix normals
call surface%sanitize_normals

! 3. Center the dragon at the origin
call surface%translate(delta=-surface%centroid)

! 4. Scale to fit in a unit box
call surface%resize(factor=(1.0_R8P / surface%bmax) * ex_R8P + &
                           (1.0_R8P / surface%bmax) * ey_R8P + &
                           (1.0_R8P / surface%bmax) * ez_R8P)

! 5. Save
call file_stl%save_into_file(facet=surface%facet, file_name='dragon-normalized.stl')
```
