---
title: Alpha wrapping
---

# Alpha wrapping

## What it does

Takes any triangle soup — broken, holed, self-intersecting,
non-manifold, or all of the above — and wraps it in a closed,
2-manifold surface at the scale `α` and the offset `ε` (Portaneri et
al. 2022), the "make any CAD junk usable for CFD" primitive. The wrap
is the boundary of the space a ball of diameter `α` carves from
outside, laid on the surface at distance `ε` from the input:

- a gap wider than about `α + 2ε` lets the ball in, and the wrap lines
  both of its sides (a gap is a feature to keep);
- a narrower gap, a hole or a crack, is closed by a cap where the ball
  stops;
- everywhere the ball reaches the input, the wrap lies on the
  `ε`-offset surface.

The existing `surface%sanitize` pipeline only fixes *local* defects:
degenerate facets, duplicate facets, near-coincident vertices, normal
orientation. It cannot recover from missing facets, large
self-intersections, internal walls, or non-manifold patches. For
AMR cut-cell tagging on CAD-derived STL, an irrecoverably broken input
silently produces wrong inside/outside labels — `is_watertight()`
returning `.true.` is **necessary but not sufficient** for a usable
solid. Alpha wrapping is what makes that assumption true.

## Pipeline

Portaneri et al. carve a 3-D Delaunay triangulation. FOSSIL carves a
uniform grid of cells of size `h` (`α/4` at first) instead, with the
same rule, and puts the vertices on the offset surface (issue #50):

1. **Distance.** The unsigned distance `u` to the input at every cell
   centre, by the input's AABB tree. Blocks of cells far from the
   input are culled.
2. **Flood fill.** Free cells are those where the ball fits outside
   the offset surface, `u > α/2 + ε`. A flood fill from the grid
   corner over free cells gives the reached set: where the ball's
   centre can go from outside. Cavities it cannot reach are filled.
3. **Carving.** A cell is carved if a ball whose centre is reachable
   covers it, which holds within `u(r) − ε` of a reached centre `r`.
   That is an exact power-distance transform, separable like the
   squared distance transform. A cell that touches the input is never
   carved, so the input always lies in solid cells, even a sheet of
   zero thickness carved on both sides. This is the offset-surface
   barrier: the two sides of a wall are separate surfaces and never
   collapse onto each other.
4. **Well-composed.** Carved cells are turned solid until no 2×2
   square has a checkerboard and no 2×2×2 block has an isolated
   diagonal pair. Marching Cubes on such a lattice gives a closed
   2-manifold.
5. **Extraction.** Marching Cubes on the cell centres. Each
   carved–solid edge gets one vertex, computed once:
   - the first point from the carved centre where `u = ε` (sphere
     tracing, which never passes the first root, so it stays on the
     carved side of any wall);
   - if the edge never comes within `ε` of the input, its midpoint
     projected onto the offset surface;
   - on a cap, the point where the carving balls end.
6. **Repair.** A triangle that meets the input, or crosses another
   wrap triangle, has its vertices moved halfway back towards their
   cell faces, where the Marching Cubes surface is embedded, for up to
   6 rounds. A chord across a convex edge of the input, with its ends
   on the offset surface, cuts the input when the cells are larger
   than `ε`. The contact test is exact (`orient3d`); the crossing test
   is in floating point, with the exact test inside a `1e-8` band.
7. **Refinement loop.** A level that still has a defect is redone with
   `h/2`, up to `max_iterations` levels (default 5) and `2^23` cells.

## API

```fortran
call surface%alpha_wrap(alpha, offset, max_iterations, wrapped, status)
   class(surface_stl_object), intent(in)            :: self
   real(R8P),                 intent(in),  optional :: alpha          ! default bbox_diag/50
   real(R8P),                 intent(in),  optional :: offset         ! default alpha/30
   integer(I4P),              intent(in),  optional :: max_iterations ! grid levels, default 5
   type(surface_stl_object),  intent(out)           :: wrapped
   integer(I4P),              intent(out), optional :: status         ! AWRAP_STATUS_*
```

- **`alpha`** is the diameter of the carving ball: the scale below
  which gaps are closed. The defaults follow CGAL's `Alpha_wrap_3`.
  The output has about `(4/α)²` triangles per unit of area, so halving
  `alpha` quadruples it.
- **`offset`** is the distance from the input at which the wrap lies.
- **`max_iterations`** is the number of grid levels the refinement
  loop may use.
- **`wrapped`** is a fresh output surface; `self` is unchanged.
- **`status`**:
  - `AWRAP_STATUS_OK`: the wrap encloses the input without touching
    it and does not cross itself.
  - `AWRAP_STATUS_NOT_CONVERGED`: the finest level allowed still has a
    triangle that touches the input or crosses the wrap. The wrap of
    that level is returned, still closed and 2-manifold.
  - `AWRAP_STATUS_BAD_INPUT`: empty surface, `alpha ≤ 0`, or
    `offset ≤ 0`.
  - `AWRAP_STATUS_DEGENERATE`: `alpha` at least the input's bounding
    box diagonal. The wrap is empty.

## Example

A "broken" cube, `cube.stl` minus one facet (a triangular hole of
inradius 0.29), wrapped twice:

```fortran
program ex_alpha_wrap
use fossil
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
implicit none

type(surface_stl_object) :: source, holed, wrapped
type(facet_object), pointer :: facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
integer(I4P)            :: status, i

! Build a "broken" cube by dropping one facet from cube.stl.
call source%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call holed%adopt_facets(facets=holed_facets)
print '(A,L1)', 'input watertight? ', holed%is_watertight()  ! .false. — has a hole

! A ball larger than the hole: the wrap caps it.
call holed%alpha_wrap(alpha=0.6_R8P, wrapped=wrapped, status=status)
print '(A,I0,A,L1,A,F8.4)', 'output status=', status, &
      '  watertight? ', wrapped%is_watertight(), &
      '  vol = ', wrapped%get_volume()
endprogram ex_alpha_wrap
```

With `alpha = 0.6` the ball cannot pass the hole: `status` is
`AWRAP_STATUS_OK`, the wrap is watertight and manifold, and its volume
is 1.109. The cube grown by the offset (0.02) has volume 1.124; the
difference is the dimple the cap makes into the hole. With the default
`alpha` (0.035) the hole is far wider than the ball. The wrap is then
one closed shell lining both sides of every wall, of volume 0.0135 and
not a solid. Check `get_volume()`, not only `is_watertight()`.

## Visual reference

A unit cube minus one facet (`cube-01.stl`, 8 facets per face, minus
its first facet: a triangle of legs 0.5 on the face x = 0, its
boundary in ember):

![the unit cube minus one facet: a triangular hole of legs 0.5 on the face x = 0, its boundary in ember](../../examples/images/ref-alpha-wrap-input.svg){.fossil-figure}

After `alpha_wrap` with `alpha` 0.15 and `offset` 0.033, cut open by
the plane y = 0.83 through the hole:

![the wrap cut open through the hole: one closed shell lining both sides of every wall at the offset, two layers 0.066 apart, entering through the hole](../../examples/images/ref-alpha-wrap-output.svg){.fossil-figure}

The ball, of diameter 0.15, fits through the hole (inradius 0.146 >
0.075 + 0.033). The wrap is one closed shell: one layer outside each
wall, one inside, 2ε apart, joined around the rim of the hole. The
barrier keeps the two layers apart: a cell touching the wall is never
carved, so each layer's vertices are traced from its own side.

With `alpha` 0.6 (and `offset` 0.02) the ball does not fit, and the
wrap closes the hole with a cap that dips into it:

![the wrap of alpha 0.6 cut open the same way: the solid cube grown by the offset, the hole closed by a cap that dips into it](../../examples/images/ref-alpha-wrap-capped.svg){.fossil-figure}

## Known limitations

- **Measured behaviour (issue #50).** `cube-01.stl` (closed) and the
  same cube without one facet, default `offset` (`α/30`). Every wrap
  is status OK, watertight and manifold:

  | `α` | closed: volume, facets | holed: volume, facets |
  |---|---|---|
  | 0.035 | 1.007, 161468 | 0.014, 309188 |
  | 0.1 | 1.020, 20168 | 0.048, 37928 |
  | 0.3 | 1.060, 2696 | 1.059, 2728 |
  | 0.7 | 1.143, 584 | 1.143, 584 |

  The closed cube's volume is the cube grown by the offset,
  `(1 + 2ε)³` less the rounding of the edges and corners: at
  `α = 0.1` the Minkowski sum is 1.0201. The hole (inradius 0.146)
  lets the ball in for `α/2 + ε < 0.146`, that is `α < 0.27`, but the
  flood fill needs a cell centre where the ball fits, so on the grid
  the hole closes earlier: open at `α = 0.1` and 0.15, closed at 0.2
  (cells of `α/4`). Below the threshold the wrap is a hollow shell. `dragon.stl` with the defaults is status OK, watertight and
  manifold, with 126916 facets and volume 484.2 (the dragon's is
  465.9); 97% of its vertices lie on the offset surface.

- **Output size and time.** Cells of `α/4` give many more triangles
  than Portaneri's (whose triangles are of size `α`): 161k for the
  closed `cube-01` at the default `α`. The time grows with the output
  (its `analyze`, and the repair pass's checks): 20 s for the closed
  `cube-01` at the default `α`, 200 s for the holed one (309k facets),
  single core.
  Decimate the wrap (`surface%decimate`) if its size matters.

- **Vertices off the offset surface.** On caps the wrap lies where
  the ball stops, up to about `α/2 + ε` from the input. A vertex the
  repair pass moved lies between the offset surface and its cell
  face, which can be closer than `ε`: 0.85 `ε` at the closest on a
  coarse marching-cubes sphere, with 98% of its vertices on the offset
  surface. Portaneri's wrap never comes closer than `ε`.

- **Concave detail below `α`** is not carved, by design. The ball
  does not enter a groove narrower than `α + 2ε`, so the wrap bridges
  it.

## See also

- [`surface%alpha_wrap` API entry](/guide/api/surface-stl-object#alpha_wrap) —
  the formal signature with intent declarations.
- [Marching Cubes](/guide/advanced/marching-cubes) — the extraction
  stage, here on a binary lattice.
- [`surface%sanitize`](/guide/api/surface-stl-object#sanitize) — for
  inputs that are *already* manifold-modulo-local-defects (degenerate
  / duplicate facets, near-coincident vertices). `sanitize` is much
  cheaper than `alpha_wrap`; reach for `alpha_wrap` only when
  sanitize cannot recover the topology.
- [`surface%decimate`](/guide/api/surface-stl-object) — to coarsen the
  wrap.
- Issue [#50](https://github.com/szaghi/FOSSIL/issues/50): the
  diagnosis of the first, voxel-shell implementation and the design of
  this one.

## References

- Portaneri, Hemmer, Birdal, Mandad & Alliez, *Alpha Wrapping with
  an Offset*, SIGGRAPH 2022. The reference algorithm.
- CGAL `Alpha_wrap_3` package documentation. The reference
  implementation that defines the canonical defaults (α =
  bbox_diag/50, offset = α/30) FOSSIL inherits.
- Felzenszwalb & Huttenlocher, *Distance Transforms of Sampled
  Functions*, Theory of Computing 8 (2012). The separable transform
  used for the carving.
- Lorensen & Cline, *Marching Cubes*, SIGGRAPH 1987; Latecki,
  *Discrete Representation of Spatial Objects in Computer Vision*
  (1998), on well-composed sets.
