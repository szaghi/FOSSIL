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

Two unit cubes offset diagonally by `(0.5, 0.5, 0.5)` (`cube-01.stl`, the
unit cube with 8 facets per face). From the same input pair, three of the
four boolean ops produce:

**`A ∪ B` — UNION** (volume 1.875 = 2 − 0.125):

![the union of two unit cubes offset by (0.5, 0.5, 0.5): one solid, the two cubes joined at their overlap](../../examples/images/ref-boolean-union.svg){.fossil-figure}

**`A ∩ B` — INTERSECT** (volume 0.125 = 0.5³):

![the intersection of the two cubes: the overlap cube of side 0.5, drawn inside the two inputs (translucent)](../../examples/images/ref-boolean-intersect.svg){.fossil-figure}

**`A \ B` — DIFFERENCE** (volume 0.875 = 1 − 0.125):

![the difference A minus B: the unit cube with its corner of side 0.5 cut away, next to the subtracted cube B (translucent)](../../examples/images/ref-boolean-difference.svg){.fossil-figure}

All three results match their analytic volume to FP precision (errors
below 5e-16, in the captions).

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

## Validation identities

A boolean result can be checked without knowing the answer: each identity below holds for any pair of closed,
outward-oriented inputs (issue #63, roadmap #58). The test-side module `src/tests/fossil_boolean_harness.f90` runs the
four ops on one pair and checks them all.

| Identity | Check | Tolerance |
|---|---|---|
| Volume | V(A ∪ B) + V(A ∩ B) = V(A) + V(B); V(A − B) = V(A) − V(A ∩ B); V(A △ B) = V(A ∪ B) − V(A ∩ B) | 1e-9 of V(A) + V(B) |
| Second moments | the same identities on $C = \int x\,x^T$ about the centre of the joint box (so on the inertia tensor $\mathrm{tr}(C)\,I - C$): where the volume is, not only how much | 1e-9 of $\lvert C(A)\rvert + \lvert C(B)\rvert$ |
| Area | area(A ∪ B) + area(A ∩ B) = area(A) + area(B), without coplanar overlaps | 1e-9 relative |
| Surface | the result lies on the inputs: its deviation from A and B merged ([surface deviation](/guide/advanced/surface-deviation), sampled lower bound) is 0 to round-off | 1e-9 of the joint diagonal |
| Membership | at fixed-seed points farther than 1e-6 of the diagonal from both inputs, inside(result) = op(inside(A), inside(B)), solid-angle sign | no mismatch |
| Validity | ∪, ∩, − without contact: watertight, manifold, V > 0, no self-intersecting facet pair | exact |

The volume and the moments come from signed tetrahedra about a fixed point, so they are defined on the non-manifold
results of contacts and of △ too. The pass criterion is #33's: a result with status OK satisfies every identity that
applies; a non-OK status (`BOOL_STATUS_CDT_FAILED`) is allowed and counted.

`src/tests/validation/fossil_val_boolean_identities.f90` runs the harness on the box cases validated so far against
analytic volumes: the 6 coplanar cases, axis-aligned and rotated, and the 18 regression stress configurations. Its
regression subset, `src/tests/fossil_test_boolean_identities.f90`, runs 3 rotated coplanar cases and one stress
configuration (the boolean dominates the cost: 0.25 s per pair in the debug coverage build). All 30 pairs
pass. Every residual is below 2e-16 and the deviation below 1e-15, with no membership mismatch. As a negative control,
the union with one facet dropped fails: its volume is off by 6 %, its moments by 14 %, its area by 4 %.

### Curved and dense meshes

The bodies are closed marching-cubes surfaces (test-side module `src/tests/fossil_boolean_bodies.f90`): a sphere of
radius 0.5, a torus of radii 0.4 and 0.15, a capped cylinder of radius 0.3 and height 1, each sampled on n³ points of
[-0.65, 0.65]³ and closed by `sanitize(fill_holes=.true.)`, plus the unit cube. The second body of each pair is moved
by a fixed-seed proper rigid motion (random axis and angle, each translation component within ±0.2). `bunny.stl` and `dragon.stl` are
not closed, so they enter through [`alpha_wrap`](/guide/advanced/alpha-wrap) (α = d/50, offset d/500 of their
diagonal d), scaled to a unit diagonal.

| Configurations | Pairs | Status OK | Pass | max volume residual | max moment residual | max deviation |
|---|---|---|---|---|---|---|
| 24³ bodies, all 10 kind pairs × 4 motions | 40 | 40 | 40 | 2.5e-16 | 1.6e-16 | 6.7e-15 |
| 40³ bodies, all 10 kind pairs × 4 motions | 40 | 40 | 40 | 2.1e-16 | 1.8e-16 | 3.5e-14 |
| wrapped bunny (16 280 facets) × wrapped dragon (9 678) × 4 motions | 4 | 4 | 4 | 1.1e-16 | 1.4e-16 | 5.2e-15 |

No membership mismatch on 200 points per pair (`src/tests/validation/fossil_val_boolean_curved.f90`, 244 s at -O2).
The regression tier, `src/tests/fossil_test_boolean_curved.f90`, runs two 10³ pairs (sphere × torus, cylinder ×
cube), 1.7 s in the debug coverage build. The first run of this sweep found the defect below.

### Near-degenerate configurations

`src/tests/validation/fossil_val_boolean_degenerate.f90` sets a body against itself, coincident or moved by a tiny
amount, plus nested bodies and contacts. The bodies are the unit cube and a marching-cubes sphere (24³, 2 972 facets).
Translations follow a generic direction; rotations turn about a generic axis through the centre. Outcome per case (✓:
status OK and every identity holds; CDT: `BOOL_STATUS_CDT_FAILED`, allowed; ✗: status OK, identities within 2e-10, but
the result is not watertight):

| Configuration | 1e-4 | 1e-6 | 1e-8 | 1e-10 | 1e-12 | 1e-14 |
|---|---|---|---|---|---|---|
| cube, translated (of the size) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cube, rotated (rad) | ✓ | ✓ | CDT | ✓ | ✓ | |
| sphere, translated | ✓ | ✓ | ✗ | CDT | ✓ | ✓ |
| sphere, rotated | ✓ | ✗ | CDT | ✗ | ✓ | |

Also ✓: the cube and the sphere against themselves, each against itself scaled by 0.5 (nested), the cube against a
half cube sharing part of a face, contacts at a vertex, along an edge and on a face patch, and the sphere against the
same sphere re-facetted (25³). 28 of the 31 cases pass; before the fixes below, 23 did, and five of the failures were
silently wrong volumes, up to 2.3 %.

Where the snap decides: a copy moved by less than 1e-10 of the local edge length is coplanar with the original
(columns 1e-12 and 1e-14, and the cube at 1e-10), and the result is the original body. A copy moved by more than
about 1e-6 is an ordinary crossing. In between, the cuts are 1e-10 to 1e-8 apart, and the retriangulation decides
with fixed tolerances (points merged within 1e-10 per coordinate in the facet's plane, points on a segment within 1e-9 of its
length) that do not match the exact classification. On the flat cube that only costs a `CDT_FAILED`. On the curved
sphere three cases come out with holes a few 1e-9 wide, with the volume, the moments and the membership still right. The
band is relative to the facets, not to the body: a coarser sphere (10³) against itself translated by 1e-6 cracks
too.
This remains open: matching the retriangulation's decisions to the classification's (symbolic points with their
provenance, rather than coordinates compared within tolerances) is a design change. Aligning the tolerances to 1e-10
was tried, and it brought back a silent 2.3 % volume error on the cube rotated by 1e-8.

### CSG chains

`src/tests/validation/fossil_val_boolean_csg.f90` runs chains of booleans whose exact volume and second moments
follow from independent booleans (test-side module `src/tests/fossil_boolean_csg.f90`). Every step is checked, so
drift along the chain would show:

- **Drilled block:** the unit cube minus k² disjoint vertical capped cylinders on a k × k grid, each through the
  block. Step j must give V(R_j) = V(R_j−1) − V(cube ∩ C_j), with the intersection computed separately against the
  original cube.
- **Sphere lattice:** the union of k² spheres on a grid, radius 0.6 of the spacing. Neighbours along the rows and
  columns overlap and diagonal ones do not, so step j must add V(S_j) minus the overlaps with its earlier neighbours.

Both use 12³ marching-cubes bodies.

| Chain | Operations | Final facets | Worst step residual (volume, moments) | Final residual (volume, moments) | Time at -O2 |
|---|---|---|---|---|---|
| block − 9 cylinders | 9 | 3 696 | 2.2e-16, 2.5e-16 | 2.2e-16, 1.8e-16 | 0.14 s |
| block − 25 cylinders | 25 | 10 688 | 1.1e-16, 1.1e-15 | 1.1e-16, 1.1e-15 | 1.1 s |
| block − 49 cylinders | 49 | 21 944 | 1.1e-16, 7.0e-16 | 1.1e-16, 7.0e-16 | 4.6 s |
| ∪ 9 spheres | 8 | 5 820 | 1.1e-16, 8.7e-17 | 1.1e-16, 8.7e-17 | 0.37 s |
| ∪ 25 spheres | 24 | 16 300 | 6.3e-16, 7.4e-17 | 6.3e-16, 4.8e-17 | 1.9 s |
| ∪ 49 spheres | 48 | 32 060 | 4.1e-16, 9.0e-17 | 4.1e-16, 4.5e-17 | 6.3 s |

Every status is OK and every intermediate body is watertight and manifold. The residuals stay at round-off and do
not grow with the length: there is no drift. One step costs about N^1.2 in the facets of the running body (0.027 s
at 3 700 facets, 0.22 s at 22 000). Before the fix of *Performance* it was N^1.8 (2.0 s at 22 000), and the 49-step
chains took 35 s and 40 s. The
regression tier, `src/tests/fossil_test_boolean_csg.f90`, runs 4-operation chains of both kinds.

### Performance

`src/tests/validation/fossil_val_boolean_perf.f90` times a boolean (best of 3, -O2) against the facet count:

| Case | Facets (A + B) | Time | µs per facet | Before (#63) |
|---|---|---|---|---|
| sphere 16³ − small cylinder | 1 656 | 0.021 s | 12.9 | 0.049 s |
| sphere 32³ − small cylinder | 5 688 | 0.053 s | 9.4 | 0.18 s |
| sphere 64³ − small cylinder | 22 440 | 0.19 s | 8.4 | 0.76 s |
| sphere 96³ − small cylinder | 50 808 | 0.43 s | 8.4 | 1.81 s |
| sphere 16³ − moved copy | 2 680 | 0.064 s | 23.8 | 0.12 s |
| sphere 32³ − moved copy | 10 744 | 0.26 s | 23.8 | 0.56 s |
| sphere 64³ − moved copy | 44 248 | 0.80 s | 18.1 | 2.26 s |
| sphere 96³ − moved copy | 100 984 | 1.65 s | 16.4 | 5.51 s |
| wrapped bunny − wrapped dragon (each op) | 25 958 | 1.1–1.2 s | 42–46 | 2.7–2.8 s |

The cost per facet is flat in N: the boolean is linear in practice. Selection dominated (about 95 % of the time).
For every sub-triangle it computed the winding number of both surfaces, although only the other surface's decides,
and the one of its own surface, at a point lying on that surface, is the expensive query. It also ran an unbounded
closest-point search for the geometric shared test. Selection now computes only the other surface's winding number,
none for a shared sub-triangle, and the search stops at the test's own bound, 1e-6 of the longest edge
(`distance_tree_with_region(max_distance2=)`): same decisions, 2.4–4.2× faster. What remains is mostly the winding
number queries, one per sub-triangle. Flood-filling the side over the uncut facets would skip most of them; it was not
done, because one wrong query would then misplace a whole patch.

### Defects found by the harness

- **Delaunay store overflow** (fixed, issue #63). The constrained triangulation of a cut facet
  (`fossil_dt`) allocated its triangle store once, at 6 (n + 3) entries. Bowyer–Watson keeps the triangles it
  replaces, so the store fills with the sum of the cavity sizes. When the cut points arrive along a curve, as a torus
  cutting a large cube face does, the cavities are large, and the store was overrun: heap corruption, an abort, or a
  silently wrong result. The store now grows by doubling. `src/tests/fossil_test_dt_capacity.f90` covers both the
  triangulation (600 points on an elliptic arc create 27 087 triangles against 3 618 entries) and the boolean that
  failed.
- **Facet edge lost by the retriangulation** (fixed, issue #63). Each cut facet was triangulated with its cuts as
  constraints, but not its own three edges. Bowyer–Watson with a finite super-triangle does not keep a hull edge when
  a point lies close inside it, and the facet across that edge keeps it whole: a crack. The cube against itself
  rotated by 1e-4 rad (cuts 3.5e-5 from the unit edges) lost 28 half-edges, with status OK and a volume off by 9e-6.
  The facet edges are now constraints too.
- **Folded Bowyer–Watson cavity** (fixed, issue #63). `incircle` is not exact, and on near-cocircular points (the
  cuts of near-coincident surfaces) the cavity of an inserted point was not star-shaped from it, so the fan of new
  triangles folded. The cube against itself rotated by 1e-8 rad lost 2.3 % of its volume with status OK. A neighbour
  across an edge that the point does not see strictly (exact `orient2d`) now always joins the cavity. That case now
  returns `BOOL_STATUS_CDT_FAILED`, which is allowed.
- **Shared sub-triangles misdetected** (fixed, issue #63). A sub-triangle was "shared" with the other surface (a
  coplanar overlap, kept once) when its centroid lay within 1e-9 L of the closest facet of the other surface and the
  normals were parallel to 0.99, that is within 8°. Sub-triangles a few 1e-9 wide beside a crossing of two facets a
  few degrees apart passed both tests: a marching-cubes sphere against itself translated by 1e-6 cracked. A
  sub-triangle is now shared when its centroid lies inside a facet of the other surface that the exact pair
  classification found coplanar with its parent. The geometric test remains only for faces within the snap whose
  boxes are disjoint, so that they never met the classifier (a face and its copy moved by 1e-12), with normals
  parallel to 1e-15.

`src/tests/fossil_test_boolean_degenerate.f90` covers the first two with the cube cases, which fail on the code
before the fixes. The third needs a curved body and is too slow for the regression tier, so the validation sweep
gates its case, the sphere translated by 1e-6.

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
