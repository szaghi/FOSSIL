---
title: facet_object
---

# `facet_object`

A single triangular facet — three vertices, an outward unit normal, a
centroid, an id, the three connected-facet ids by edge, and a pair of
angle-weighted pseudo-normal fields needed by the Bærentzen–Aanæs sign test.

You **rarely instantiate a `facet_object` yourself**. Facets are created and
owned by [`surface_stl_object`](/guide/api/surface-stl-object); you walk over
them via [`surface%facet_at(i)`](/guide/api/surface-stl-object#getters-and-predicates).
This page documents the public surface a caller would use *through that
pointer* — typically for diagnostics, custom traversals, or per-facet metrics
that the surface API does not provide directly.

The page is in two parts: a list of **public members** (data fields you read
from a facet pointer), then the **public methods**.

[[toc]]

---

## Public members

| Member               | Type                       | Description                                                                                                  |
| -------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `normal`             | `type(vector_R8P)`         | Outward unit normal.                                                                                         |
| `vertex(3)`          | `type(vector_R8P)`         | The three vertices, in winding order.                                                                        |
| `centroid`           | `type(vector_R8P)`         | Facet centroid (average of the three vertices).                                                              |
| `id`                 | `integer(I4P)`             | Global 1-based facet id within the surface.                                                                  |
| `fcon_edge(3)`       | `integer(I4P)`             | Connected-facet ids by edge. Index 1 = edge `v1→v2`, 2 = `v2→v3`, 3 = `v3→v1`. Value 0 means disconnected.   |
| `bb(2)`              | `type(vector_R8P)`         | Facet AABB: `bb(1)` = min corner, `bb(2)` = max corner.                                                      |
| `edge_pnormal(3)`    | `type(vector_R8P)`         | Angle-weighted pseudo-normals at each edge (populated by `surface%analyze`).                                 |
| `vertex_pnormal(3)`  | `type(vector_R8P)`         | Angle-weighted pseudo-normals at each vertex (populated by `surface%analyze`).                               |

**Example.** Walk every facet, print summary.

```fortran
program ex_facet_members
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
integer(I4P)                 :: i

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   print '(A,I0,A,3F8.3,A,3F8.3)', 'facet ', f%id,                            &
                                   ' normal=(', f%normal%x, f%normal%y, f%normal%z, ')', &
                                   ' centroid=(', f%centroid%x, f%centroid%y, f%centroid%z, ')'
end do
end program ex_facet_members
```

**Notes.**

- `facet_at(i)` returns `null()` for `i` outside `[1, get_facets_number()]`;
  always check.
- `fcon_edge` and the pseudo-normal arrays are only meaningful after
  `surface%analyze` has run (which `load_from_file` and `sanitize` do for you).

---

## Methods

### `compute_distance`

```fortran
call facet%compute_distance(point=p, distance=d2)
```

**Purpose.** Squared unsigned distance from a point to this facet. The thin
wrapper around [`compute_distance_with_region`](#compute_distance_with_region)
when you only want the scalar — closest-point and Voronoi-region outputs are
discarded internally.

**Arguments.**

| Argument   | Intent | Type               | Notes                                       |
| ---------- | ------ | ------------------ | ------------------------------------------- |
| `point`    | `in`   | `type(vector_R8P)` | (required) Query point.                     |
| `distance` | `out`  | `real(R8P)`        | (required) Squared distance point–facet.    |

**Example.** Find the closest facet by brute force (the surface API does this
faster — this is for illustration).

```fortran
program ex_facet_distance
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
type(vector_R8P)             :: p
real(R8P)                    :: d2, best
integer(I4P)                 :: i, best_i

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
p = 2._R8P * ex_R8P
best   = huge(0._R8P)
best_i = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   call f%compute_distance(point=p, distance=d2)
   if (d2 < best) then ; best = d2 ; best_i = i ; end if
end do
print '(A, I0, A, ES12.5)', 'closest facet=', best_i, ' d=', sqrt(best)
end program ex_facet_distance
```

**Notes.**

- Use [`surface%distance`](/guide/api/surface-stl-object#distance) instead for
  real workloads — it is tree-accelerated.
- The facet's "metrix" (plane-equation coefficients) must already be
  computed; `surface%analyze` (called by `load_from_file` and `sanitize`)
  takes care of this.

**See also.** [`compute_distance_with_region`](#compute_distance_with_region),
[`surface%distance`](/guide/api/surface-stl-object#distance).

---

### `compute_distance_with_region`

```fortran
call facet%compute_distance_with_region(point=p, distance=d2, closest=c, region=r)
```

**Purpose.** Same squared distance as [`compute_distance`](#compute_distance),
plus the **closest point** on the facet and a tag identifying *which Voronoi
region of the triangle* that point lies in (interior face, one of three
edges, or one of three vertices). The region tag is the input to
[`pseudo_normal_for_region`](#pseudo_normal_for_region), and together they
form the building block of the default
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal) sign test.

**Arguments.**

| Argument   | Intent | Type               | Notes                                                                  |
| ---------- | ------ | ------------------ | ---------------------------------------------------------------------- |
| `point`    | `in`   | `type(vector_R8P)` | (required) Query point.                                                |
| `distance` | `out`  | `real(R8P)`        | (required) Squared distance.                                           |
| `closest`  | `out`  | `type(vector_R8P)` | (required) Closest point on the facet.                                 |
| `region`   | `out`  | `integer(I4P)`     | (required) Voronoi region tag — see encoding below.                    |

**Region encoding.**

| Value         | Meaning                                                       |
| ------------- | ------------------------------------------------------------- |
|  0            | Interior of the face. Use `facet%normal`.                     |
|  1, 2, 3      | On edge 1 (`v1→v2`), 2 (`v2→v3`), 3 (`v3→v1`) interior.       |
| -1, -2, -3    | At vertex 1, 2, or 3.                                         |

**Example.** Classify the closest feature.

```fortran
program ex_compute_distance_with_region
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
type(vector_R8P)             :: p, closest
real(R8P)                    :: d2
integer(I4P)                 :: region

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
f => surface%facet_at(1)
p = 2._R8P * ex_R8P
call f%compute_distance_with_region(point=p, distance=d2, closest=closest, region=region)
print '(A, I0)', 'region tag = ', region
print '(A, 3ES12.5)', 'closest pt = ', closest%x, closest%y, closest%z
end program ex_compute_distance_with_region
```

**Notes.**

- Algorithm: David Eberly's *Distance Between Point and Triangle in 3D*
  (Geometric Tools).
- A region tag of `0` (interior face) is by far the most common; ties go to
  the highest-dimensional region (face beats edge beats vertex).

**See also.** [`pseudo_normal_for_region`](#pseudo_normal_for_region),
[`surface%compute_distance`](/guide/api/surface-stl-object#compute_distance).

---

### `pseudo_normal_for_region`

```fortran
n = facet%pseudo_normal_for_region(region=region)
```

**Purpose.** Return the angle-weighted pseudo-normal of the facet at the
Voronoi region selected by `region`. Combined with the closest-point output
of [`compute_distance_with_region`](#compute_distance_with_region) this gives
the Bærentzen–Aanæs sign for signed-distance queries:

```
sign = sign( dot(point - closest, pseudo_normal_for_region(region)) )
```

**Arguments.**

| Argument   | Intent | Type           | Notes                                                                                   |
| ---------- | ------ | -------------- | --------------------------------------------------------------------------------------- |
| `region`   | `in`   | `integer(I4P)` | (required) Region tag from `compute_distance_with_region`. See encoding on that page.   |

**Returns.** `type(vector_R8P)` — pseudo-normal at the requested feature.

**Example.** Custom signed-distance loop, equivalent to
[`surface%distance(..., is_signed=.true.)`](/guide/api/surface-stl-object#distance):

```fortran
program ex_pseudo_normal_for_region
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
type(vector_R8P)             :: p, closest, n, delta
real(R8P)                    :: d2, signed_d
integer(I4P)                 :: region

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
f => surface%facet_at(1)
p = 0.5_R8P * ex_R8P
call f%compute_distance_with_region(point=p, distance=d2, closest=closest, region=region)
n         = f%pseudo_normal_for_region(region=region)
delta     = p - closest
signed_d  = sign(sqrt(d2), delta%dotproduct(rhs=n))
print '(A, ES12.5)', 'signed distance = ', signed_d
end program ex_pseudo_normal_for_region
```

**Notes.**

- The pseudo-normal fields (`facet%edge_pnormal`, `facet%vertex_pnormal`) are
  populated by `surface%analyze`. If you mutate facets behind the surface's
  back and don't re-analyse, this function returns stale data.

**See also.** [`compute_distance_with_region`](#compute_distance_with_region),
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal).

---

### `solid_angle`

```fortran
omega = facet%solid_angle(point=p)
```

**Purpose.** Return the *projected* solid angle subtended by the facet from
`point`, in steradians. Used by
[`SIGN_SOLID_ANGLE`](/guide/api/constants#sign_solid_angle): summing
`solid_angle` over every facet of a closed surface gives ±4π inside and 0
outside, with intermediate values for non-closed meshes (which is the
foundation of the *generalized winding number*).

**Arguments.**

| Argument   | Intent | Type               | Notes                            |
| ---------- | ------ | ------------------ | -------------------------------- |
| `point`    | `in`   | `type(vector_R8P)` | (required) Query point.          |

**Returns.** `real(R8P)` — solid angle in steradians.

**Example.** Manual generalized winding number (slow — for didactic use).

```fortran
program ex_solid_angle
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
type(vector_R8P)             :: p
real(R8P)                    :: w
integer(I4P)                 :: i

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
p = 0.0_R8P * ex_R8P
w = 0._R8P
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   w = w + f%solid_angle(point=p)
end do
print '(A, ES12.5)', 'winding number = ', w / (4._R8P * acos(-1._R8P))
end program ex_solid_angle
```

**Notes.**

- Implementation: Van Oosterom–Strackee formula. Robust at the cost of two
  square roots and one `atan2` per facet — hence `SIGN_SOLID_ANGLE` is O(N)
  per query without acceleration.
- The fast *hierarchical* generalized winding number (Barill et al., 2018)
  is shipped — see [`surface%winding_number`](/guide/advanced/winding-number).

**See also.** [`SIGN_SOLID_ANGLE`](/guide/api/constants#sign_solid_angle).

---

### `do_ray_intersect`

```fortran
hit = facet%do_ray_intersect(ray_origin=o, ray_direction=d)
```

**Purpose.** Boolean ray–triangle intersection test (Möller–Trumbore).
Returns `.true.` if the ray from `o` along `d` hits the facet at `t > 0`.

**Arguments.**

| Argument        | Intent | Type               | Notes                                       |
| --------------- | ------ | ------------------ | ------------------------------------------- |
| `ray_origin`    | `in`   | `type(vector_R8P)` | (required) Ray origin.                      |
| `ray_direction` | `in`   | `type(vector_R8P)` | (required) Ray direction. Need not be unit. |

**Returns.** `logical` — true on intersection.

**Example.** Count intersections of a +x ray (the classic odd-parity inside
test).

```fortran
program ex_do_ray_intersect
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: surface
type(facet_object), pointer  :: f
type(vector_R8P)             :: o
integer(I4P)                 :: i, hits

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
o    = surface%get_centroid()
hits = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   if (f%do_ray_intersect(ray_origin=o, ray_direction=ex_R8P)) hits = hits + 1
end do
print '(A, I0, A, L1)', 'hits=', hits, '  inside? ', mod(hits, 2) == 1
end program ex_do_ray_intersect
```

**Notes.**

- The bare test is **boolean**; the parameter `t` along the ray is computed
  internally but not returned. A future ray-API surface on `surface_stl_object`
  (issue #18 §2.5) will expose `t`, hit point, and hit facet ids.
- Returns `.false.` when the ray is parallel to the facet (back-face hits are
  not distinguished from misses).

**See also.** [`SIGN_RAY_INTERSECTIONS`](/guide/api/constants#sign_ray_intersections).

---

### `compute_normal` *(advanced)*

```fortran
call facet%compute_normal
```

**Purpose.** Recompute the outward unit normal from the current winding
(`(v2-v1) × (v3-v1)`, normalised). Called automatically by
`surface%analyze` and `surface%sanitize_normals`; you only invoke it
manually after rewriting vertices behind the surface's back.

**Arguments.** None.

**See also.** [`compute_metrix`](#compute_metrix),
[`surface%sanitize_normals`](/guide/api/surface-stl-object#sanitize_normals).

---

### `compute_metrix` *(advanced)*

```fortran
call facet%compute_metrix
```

**Purpose.** Recompute the plane-equation coefficients used by
[`compute_distance_with_region`](#compute_distance_with_region) (edge vectors
`E12`, `E13`, dot-product cache `a`, `b`, `c`, determinant `det`, centroid,
normal). Called automatically by `surface%analyze` and by the manipulation
TBPs when `recompute_metrix=.true.`.

**Arguments.** None.

**See also.** [`compute_normal`](#compute_normal),
[`surface%analyze`](/guide/api/surface-stl-object#analyze).
