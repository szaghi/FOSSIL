---
title: API Reference
---

# API Reference

FOSSIL exposes a single module:

```fortran
use fossil
```

The module re-exports three public types from internal modules. All geometry operations use [PENF](https://github.com/szaghi/PENF) numeric kinds (`I4P`, `R8P`) and [VecFor](https://github.com/szaghi/VecFor) vector types (`vector_R8P`).

---

## `file_stl_object`

Handles STL file I/O. A single instance can be reused to load or save multiple surfaces sequentially.

### Public members

| Member | Type | Description |
|--------|------|-------------|
| `file_name` | `character(len=:), allocatable` | Path of the STL file |
| `is_ascii` | `logical` | `.true.` for ASCII format, `.false.` for binary |
| `header` | `character(FRLEN)` | STL file header string |

### Public methods

| Method | Description |
|--------|-------------|
| [`initialize`](#initialize) | Initialize the file handler |
| [`load_from_file`](#load_from_file) | Load facets from an STL file into a surface |
| [`save_into_file`](#save_into_file) | Save facets to an STL file |
| [`save_aabb_into_file`](#save_aabb_into_file) | Save the AABB tree as an STL file |
| [`open_file`](#open_file) | Open the file manually (after `initialize`) |
| [`close_file`](#close_file) | Close the file |
| [`destroy`](#destroy-file) | Reset to default state |
| [`statistics`](#statistics-file) | Return a formatted string with file statistics |

---

### `initialize` {#initialize}

Set up the file handler before loading or saving.

```fortran
use fossil
type(file_stl_object) :: file_stl

call file_stl%initialize(file_name='cube.stl', is_ascii=.true.)
```

| Argument | Intent | Type | Description |
|----------|--------|------|-------------|
| `file_name` | `in`, optional | `character(*)` | File path |
| `is_ascii` | `in`, optional | `logical` | Format flag (default `.true.`) |
| `skip_destroy` | `in`, optional | `logical` | Skip reset before init (default `.false.`) |

---

### `load_from_file` {#load_from_file}

Load STL facets from a file directly into `surface%facet`. The surface must call `analize` afterwards to build connectivity and the AABB tree.

```fortran
use fossil
type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

! Auto-detect format
call file_stl%load_from_file(facet=surface%facet, file_name='part.stl', guess_format=.true.)

! With on-the-fly clipping
call file_stl%load_from_file(facet=surface%facet, file_name='part.stl', &
                             guess_format=.true., clip_min=bmin, clip_max=bmax)
call surface%analize
```

| Argument | Intent | Type | Description |
|----------|--------|------|-------------|
| `facet` | `out` | `type(facet_object), allocatable(:)` | Receives the loaded facets |
| `file_name` | `in`, optional | `character(*)` | File path |
| `is_ascii` | `in`, optional | `logical` | Format flag |
| `guess_format` | `in`, optional | `logical` | Auto-detect ASCII vs binary |
| `clip_min` | `in`, optional | `type(vector_R8P)` | Discard facets below this bounding corner |
| `clip_max` | `in`, optional | `type(vector_R8P)` | Discard facets above this bounding corner |

---

### `save_into_file` {#save_into_file}

Write facets to an STL file.

```fortran
call file_stl%save_into_file(facet=surface%facet, file_name='output.stl')
```

| Argument | Intent | Type | Description |
|----------|--------|------|-------------|
| `facet` | `in` | `type(facet_object), allocatable(:)` | Facets to write |
| `file_name` | `in`, optional | `character(*)` | Output file path |
| `is_ascii` | `in`, optional | `logical` | Format flag (default: binary) |

---

### `statistics` (file) {#statistics-file}

Returns a multi-line formatted string describing the file (name, format, bounding extents).

```fortran
print '(A)', file_stl%statistics()
```

---

## `surface_stl_object`

Holds the triangulated surface geometry and provides all analysis and manipulation methods.

### Public members

| Member | Type | Description |
|--------|------|-------------|
| `facets_number` | `integer(I4P)` | Number of triangular facets |
| `facet` | `type(facet_object), allocatable(:)` | Array of facets |
| `bmin`, `bmax` | `type(vector_R8P)` | Axis-aligned bounding box corners |
| `volume` | `real(R8P)` | Enclosed volume (negative if normals point inward) |
| `centroid` | `type(vector_R8P)` | Surface centroid |
| `aabb` | `type(aabb_tree_object)` | AABB octree for acceleration |
| `facet_1_de` | `type(list_id_object)` | IDs of facets with one disconnected edge |
| `facet_2_de` | `type(list_id_object)` | IDs of facets with two disconnected edges |
| `facet_3_de` | `type(list_id_object)` | IDs of facets with three disconnected edges |

### Public methods

| Method | Description |
|--------|-------------|
| [`analize`](#analize) | Compute bounding box, connectivity, volume, centroid, and AABB tree |
| [`translate`](#translate) | Translate facets |
| [`rotate`](#rotate) | Rotate facets |
| [`mirror`](#mirror) | Mirror facets |
| [`resize`](#resize) | Scale facets |
| [`clip`](#clip) | Clip surface to an axis-aligned bounding box |
| [`merge_solids`](#merge_solids) | Merge another surface into this one |
| [`sanitize_normals`](#sanitize_normals) | Make all normals consistent |
| [`reverse_normals`](#reverse_normals) | Flip all normals |
| [`connect_nearby_vertices`](#connect_nearby_vertices) | Repair disconnected edges |
| [`distance`](#distance) | Minimum distance from a point to the surface |
| [`is_point_inside`](#is_point_inside) | Point-in-polyhedron test |
| [`compute_mesh_distance`](#compute_mesh_distance) | Distance field over a structured mesh |
| [`compute_volume`](#compute_volume) | Recompute volume from current facets |
| [`compute_centroid`](#compute_centroid) | Recompute centroid from current facets |
| [`build_connectivity`](#build_connectivity) | Rebuild facet connectivity only |
| [`statistics`](#statistics-surface) | Return a formatted string with surface statistics |
| [`destroy`](#destroy-surface) | Reset to default state |
| [`initialize`](#initialize-surface) | Initialize (optionally skip destroy) |

---

### `analize` {#analize}

Full surface analysis: bounding box, volume, centroid, facet connectivity, disconnected-edge detection, and AABB octree. Call this after every load or geometry modification.

```fortran
call surface%analize
```

---

### `translate` {#translate}

```fortran
! By 3D vector
call surface%translate(delta=delta_vector)

! By scalar components (any subset)
call surface%translate(x=1.0_R8P)
call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
```

---

### `rotate` {#rotate}

```fortran
! By axis + angle (radians)
call surface%rotate(axis=axis_vector, angle=angle_radians)

! By rotation matrix
call surface%rotate(matrix=rotation_matrix)
```

---

### `mirror` {#mirror}

```fortran
! By plane normal
call surface%mirror(normal=normal_vector)

! By mirror matrix
call surface%mirror(matrix=mirror_matrix)
```

---

### `resize` {#resize}

```fortran
! Vectorial scale factor
call surface%resize(factor=factor_vector)

! Per-axis scalars
call surface%resize(x=0.5_R8P, z=1.2_R8P)

! Scale about centroid
call surface%resize(factor=factor_vector, respect_centroid=.true.)
```

---

### `clip` {#clip}

```fortran
type(surface_stl_object) :: remainder

call surface%clip(bmin=bmin_vector, bmax=bmax_vector, remainder=remainder)
call surface%analize
call remainder%analize
```

The `remainder` argument is optional — omit it when the removed facets are not needed.

---

### `merge_solids` {#merge_solids}

```fortran
call surface_1%merge_solids(other=surface_2)
```

Appends `surface_2` facets into `surface_1`. Call `analize` afterwards to rebuild connectivity.

---

### `sanitize_normals` {#sanitize_normals}

Make all normals consistent (propagates orientation from each facet to its neighbors via connectivity).

```fortran
call surface%sanitize_normals
call surface%compute_volume   ! re-check volume after sanitization
```

---

### `reverse_normals` {#reverse_normals}

Flip every facet normal. Useful after `mirror` or when normals point inward and need to be flipped uniformly.

```fortran
call surface%reverse_normals
```

---

### `connect_nearby_vertices` {#connect_nearby_vertices}

Snap vertices that are geometrically close to repair disconnected edges.

```fortran
call surface%connect_nearby_vertices
call surface%analize
```

---

### `distance` {#distance}

Returns the minimum distance (Euclidean, unsigned) from a query point to the surface. Uses the AABB tree for acceleration.

```fortran
use vecfor, only: ex_R8P, ey_R8P, ez_R8P
use penf, only: R8P

real(R8P) :: d

d = surface%distance(point=1.0_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.0_R8P * ez_R8P)
```

---

### `is_point_inside` {#is_point_inside}

Returns `.true.` if the point is inside the closed surface. Internally calls `is_point_inside_polyhedron_sa` (solid angle method).

```fortran
logical :: inside
inside = surface%is_point_inside(point=query_point)
```

---

### `compute_mesh_distance` {#compute_mesh_distance}

Compute the signed distance from every grid node of a structured block to the surface.

```fortran
call surface%compute_mesh_distance(block=my_block)
```

---

### `statistics` (surface) {#statistics-surface}

Returns a multi-line string with: bounding extents, volume, centroid, facet count, disconnected-edge counts, and AABB refinement levels.

```fortran
print '(A)', surface%statistics()
```

---

## `facet_object`

Represents a single triangular facet. Most user code accesses facets through `surface%facet(:)` rather than directly.

### Public members

| Member | Type | Description |
|--------|------|-------------|
| `normal` | `type(vector_R8P)` | Outward unit normal |
| `vertex(3)` | `type(vector_R8P)` | Three vertices |
| `centroid` | `type(vector_R8P)` | Facet centroid |
| `id` | `integer(I4P)` | Global facet ID |
| `fcon_edge_12` | `integer(I4P)` | ID of connected facet along edge 1-2 (0 = disconnected) |
| `fcon_edge_23` | `integer(I4P)` | ID of connected facet along edge 2-3 |
| `fcon_edge_31` | `integer(I4P)` | ID of connected facet along edge 3-1 |
| `bb(2)` | `type(vector_R8P)` | Facet AABB: `bb(1)` = min corner, `bb(2)` = max corner |

### Selected methods

| Method | Description |
|--------|-------------|
| `compute_normal` | Compute normal from vertices |
| `compute_metrix` | Compute plane equation coefficients used for distance queries |
| `compute_distance` | Unsigned squared distance from a point to this facet |
| `solid_angle` | Projected solid angle subtended by this facet from a given point |
| `do_ray_intersect` | Test whether a ray intersects this facet |
| `translate`, `rotate`, `mirror`, `resize` | Per-facet geometry transforms |
| `load_from_file_ascii` / `load_from_file_binary` | Low-level facet I/O |
| `save_into_file_ascii` / `save_into_file_binary` | Low-level facet I/O |
