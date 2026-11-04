---
title: Boolean operations
---

# Boolean operations

## What it does

Combines two surfaces using one of the four standard set operations:
union (`A ∪ B`), intersection (`A ∩ B`), difference (`A \ B`), or
symmetric difference (`A △ B`). The result is a fresh watertight
solid that replaces `self` in place.

This is the single biggest missing primitive for an STL toolkit —
**every meshing pipeline that combines parts needs booleans**. CSG
modelling, cutting holes for cooling channels, removing fixturing from
a part for CFD analysis, building wings from intersecting NACA
sections. FOSSIL's boolean is on par with libigl's `igl::mesh_boolean`
and CGAL's nef-polyhedra-free boolean for the *generic* case.

## Pipeline

The implementation follows Zhou et al. (SIGGRAPH 2016) "Mesh
Arrangements for Solid Geometry":

1. **Cross-mesh tree-vs-tree intersection collection**
   (`fossil_arrangement`). Both surfaces' AABB trees walked in
   lockstep; every triangle pair whose AABBs overlap is tested for
   actual triangle-triangle intersection (Möller). The collected
   intersection segments form the *arrangement*.
2. **CDT-based facet retriangulation**. Every input triangle that
   participates in an intersection is split into sub-triangles along
   its share of the arrangement segments. FOSSIL's own `fossil_dt`
   module (Bowyer-Watson + Sloan constraint recovery) does the
   constrained Delaunay triangulation in the facet's plane.
3. **Per-sub-triangle classification** (`fossil_boolean`). For each
   sub-triangle's centroid, evaluate the [winding
   number](/guide/advanced/winding-number) against both input solids:
   `(in_A?, in_B?)` is a 2-bit tag.
4. **Truth-table selection** per operation:

   | Tag       | UNION | INTERSECT | DIFFERENCE (A\B) | SYMDIFF |
   |-----------|:-----:|:---------:|:----------------:|:-------:|
   | (in, in)  | keep  | keep      | drop             | drop    |
   | (in, out) | keep  | drop      | keep             | keep    |
   | (out, in) | keep  | drop      | drop             | keep    |
   | (out, out)| drop  | drop      | drop             | drop    |

   Special-case row for sub-triangles **on the shared boundary** of A
   and B: their orientation determines which operation wants them.
5. **Adopt** the surviving sub-triangles into `self` via
   `adopt_facets` (rebuilds AABB tree, vertex pool, connectivity,
   pseudo-normals).

## API

```fortran
call surface%boolean(other, op, status)
   class(surface_stl_object),  intent(inout)         :: self
   class(surface_stl_object),  intent(in)            :: other
   integer(I4P),               intent(in)            :: op      ! BOOL_*
   integer(I4P),               intent(out), optional :: status  ! BOOL_STATUS_*
```

- **`op`** is one of:
  - `BOOL_UNION` — everything inside either body
  - `BOOL_INTERSECT` — only what's inside both
  - `BOOL_DIFFERENCE` — inside `self` (A) but not inside `other` (B)
  - `BOOL_SYMDIFF` — inside exactly one of `self`, `other`
- **`status`** ∈ `{BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED,
  BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT}`. Most
  failures are degenerate-input-related (CDT couldn't recover a
  near-collinear constraint). The CDT failure mode is documented in
  `fossil_dt`'s docstring.

**Pre-conditions for clean output** (enforced by example, not by
runtime check — set them up yourself):

- Both inputs are watertight, manifold solids (`is_volume() ==
  .true.`).
- Both have outward-oriented normals — run `surface%sanitize_normals`
  if you're not sure.
- The intersection between A and B doesn't trip the CDT's recovery
  loop. Coplanar faces, shared edges and touching contacts are handled
  (issue #33); see [Known limitations](#known-limitations) for the
  exactness statement.

## Example

Two unit cubes offset diagonally by `(0.5, 0.5, 0.5)` — the canonical
mesh-arrangement test case:

```fortran
program ex_boolean
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

! In-place boolean: a is replaced by (a OP b).
call a%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
print '(A,I0,A,F8.4)', 'A \ B  status=', status, '  vol=', a%get_volume()
! Expected analytic volume: 1 - 0.5^3 = 0.875
endprogram ex_boolean
```

## Visual reference

Two unit cubes offset diagonally by `(0.5, 0.5, 0.5)`. From the same
input pair, the four boolean ops produce:

**`A ∪ B` — UNION** (volume 1.875 = 2 − 0.125):

![Boolean UNION](/pictures/advanced/boolean-union.png)

**`A ∩ B` — INTERSECT** (volume 0.125 = 0.5³):

![Boolean INTERSECT](/pictures/advanced/boolean-intersect.png)

**`A \ B` — DIFFERENCE** (volume 0.875 = 1 − 0.125):

![Boolean DIFFERENCE](/pictures/advanced/boolean-difference.png)

All three results match their analytic volume to FP precision.

### Coplanar faces

Two unit cubes offset along x only: their y and z faces are coplanar and
overlap. This is the box ∖ body configuration of cut-cell meshing, and it
needs no perturbation:

```fortran
program ex_boolean_coplanar
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: a, b
integer(I4P)             :: status

call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P))
call b%analyze
call a%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
print '(A,I0,A,F8.5,A,L1,A,L1)', 'status ', status, '  volume ', a%get_volume(), &
      ' (analytic 0.5)  watertight ', a%is_watertight(), '  manifold ', a%is_manifold()
endprogram ex_boolean_coplanar
```

How the coplanar overlap is handled:

1. Each facet pair is classified exactly (separated, crossing, touching,
   coplanar) from `orient3d` signs, not from floating-point distances.
2. A coplanar pair cuts each facet along the parts of the other facet's
   edges that lie inside it, so the overlap becomes sub-triangles of its
   own in both facets.
3. Cut endpoints on a facet edge are also inserted in the facet across
   that edge, so no crack (T-junction) is left.
4. A sub-triangle lying on the other surface is classified by the
   shared-face table (same or opposite orientation), not by the winding
   number, which is ill-defined on the surface.

## Known limitations

- **Exactness.** The decisions (how two facets meet, which side a vertex
  is on, whether points are collinear in a shared plane) use the exact
  `orient3d` (Shewchuk's filter, then an exact expansion), with one snap: a vertex closer to the other facet's
  plane, or to the other edge's line, than 10⁻¹⁰ of the local edge length
  counts as lying on it. Faces coplanar by construction are rarely exactly
  coplanar after a rotation or a scaling in floating point, and the snap
  keeps them coplanar. Features smaller than that relative distance are
  merged by design. The intersection points are floating-point
  constructions, but canonical ones (each built by one formula from
  lexicographically ordered inputs), so every facet that shares a point
  gets the same bits and the result welds exactly.
- **What is validated.** Coplanar overlaps, shared edges and contacts are
  tested on box pairs: a 6-case matrix of all four ops, axis-aligned and
  rotated (`fossil_test_boolean_coplanar`), and 400 booleans on snapped
  random box pairs, half of them rotated (`fossil_test_boolean_stress`),
  all with analytic volumes. Curved meshes with coplanar patches have not
  been stress-tested.
- **Non-manifold results by definition.** Solids that touch only along an
  edge or at a point give a union (or difference) that is non-manifold
  there, and a symmetric difference whose two pieces `A \ B` and `B \ A`
  touch along the curves where the surfaces cross has four facets at each
  such edge. The volume is right; `is_manifold()` is false, as it should be.

- **Self-intersection in inputs** (each surface intersecting itself)
  produces undefined results. Run
  [`surface%resolve_self_intersections`](/guide/api/surface-stl-object)
  first if you suspect this is your case.

- **CDT failures** (`BOOL_STATUS_CDT_FAILED`) signal that
  near-collinear arrangement segments couldn't be recovered into the
  retriangulation. Rare in practice; usually accompanies pathological
  input geometry.

## See also

- [`resolve_self_intersections`](/guide/api/surface-stl-object) —
  self-intersection cleanup, a pre-requisite for clean booleans.
- [winding number](/guide/advanced/winding-number) — used internally
  to classify each sub-triangle.
- [Constants](/guide/api/constants) — `BOOL_*` and `BOOL_STATUS_*`.
- [Roadmap issue](https://github.com/szaghi/FOSSIL/issues/18) — full
  multi-commit history of the boolean pipeline.
- [Issue #33](https://github.com/szaghi/FOSSIL/issues/33) — coplanar
  overlaps, exact classification and the validation matrix.

## References

- Zhou, Grinspun, Zorin & Jacobson, *Mesh Arrangements for Solid
  Geometry*, SIGGRAPH 2016. The reference algorithm.
- Bernstein & Fussell, *Fast, Exact, Linear Booleans*, Symposium on
  Geometry Processing 2009. The "Cork" precursor — historical
  context.
- Möller, *A Fast Triangle-Triangle Intersection Test*, Journal of
  Graphics Tools 2(2), 1997. Used for the tree-vs-tree narrow phase.
- Sloan, *A fast algorithm for generating constrained Delaunay
  triangulations*, Computers & Structures 47(3), 1993. The constraint
  recovery used in `fossil_dt` Phase 2.
