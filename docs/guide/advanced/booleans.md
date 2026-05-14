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
  loop (true for any pair of solids in *generic position*; coplanar
  faces are the only known trouble).

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

## Known limitations

- **Axis-aligned coplanar face overlaps produce wrong volumes** even
  though `status == BOOL_STATUS_OK`. Two failure modes contribute:
  - The Möller-style tri-tri test produces degenerate or NaN
    intervals when one triangle's edge lies on the other triangle's
    plane.
  - Adjacent A-facets and B-facets sharing an edge along a coplanar
    region produce sub-triangles whose centroids are equidistant
    from both surfaces; the winding-number classifier evaluates to
    ~1 from both (rather than the textbook ~0.5), defeating the
    boundary detector that the shared-face truth table relies on.

  This is a **silent** failure (no status code flags it) — the
  boolean returns successfully but the result volume is wrong.

- **Workaround for coplanar inputs**: perturb one of the inputs by a
  small offset along **each** axis before the boolean, breaking the
  coplanarity and routing through the bit-exact generic codepath:

  ```fortran
  ! Avoid axis-aligned coplanar faces by 3D-offsetting B.
  call b%translate(delta=vector_R8P(eps, eps, eps))
  call b%analyze
  call a%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
  ```

  `eps = 1e-6 * bbox_diagonal` is a safe value: small enough to be
  numerically negligible in the output volume, large enough to
  defeat the coplanar pathology. The full discussion (with the two
  contributing failure modes) is in the
  `surface%boolean` source docstring at
  `src/lib/fossil_surface_stl.f90`.

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
  multi-commit history including the coplanar-aware shared-boundary
  tagging that handles the well-behaved subset of degenerate cases.

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
