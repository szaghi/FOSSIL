---
title: 2-D geometry (curves)
---

# 2-D geometry (curves)

## What it does

A 2-D solver on a quadtree needs the same geometry queries as a 3-D one on an
octree: distance, inside/outside, a signed-distance field. FOSSIL answers them
natively in 2-D with `curve_object`, a set of closed polylines (loops) in an
axis-aligned plane.

It does not extrude the curve into an STL. That would make every 2-D query
pay for 3-D triangles. It would also put the query plane along the
extrusion's facet edges, where sign and closest-feature decisions are most
fragile.

What it covers:

- **curves**: loops, metrics, validity, an ASCII format;
- **queries**: signed distance with closest point, segment and vertex;
  generalized winding number; inside test; exact self-intersections. They run
  on a segment BVH;
- **slicing** an STL by an axis-aligned plane, so a 2-D and a 3-D run share
  one body;
- **repair**: joining and closing open chains, merging vertices, orienting
  loops;
- **batches** of distances on the CPU (OpenMP) or the device (OpenACC);
- **fields**: 2-D signed-distance fields with
  [`distance_field_object`](/guide/advanced/distance-field), one quadtree block
  at a time.

## Conventions

- **Plane.** A curve lives in the plane normal to `normal_axis` (1, 2 or 3;
  default 3). Its in-plane coordinates (u, v) are the other two axes in
  cyclic order: (y, z), (z, x) or (x, y), so u × v points along
  +`normal_axis`. Points are passed as `p(2) = [u, v]`.
- **Orientation.** The solid is on the left of every segment. Outer loops run
  counter-clockwise, holes clockwise, and islands inside holes
  counter-clockwise again. The signed distance is negative inside, as in 3-D.
- **Open loops.** A loop may be open: a polyline, such as a slice of a mesh
  with holes. An open loop has no closing segment and no area, and
  `check_validity` reports it. The distance and the generalized winding
  number still handle it, and the fast winding number sign is meant for it.
- **Validity.** `check_validity` returns the first defect loop by loop (open
  loop, too few vertices, zero-length segment, zero area, wrong orientation for
  the nesting depth). By default it does not look for crossings;
  `check_validity(check_self_intersections=.true.)` also returns
  `CURVE_SELF_INTERSECTING` when two segments cross, touch or overlap
  (`self_intersections` lists them, decided exactly). A curve can pass the
  default check and still cross itself, a section of a non-manifold mesh for
  example; the cut cells and the swept geometry assume it does not.
- **Exactness.** Decisions that must agree from one query to the next use the
  exact `orient2d`: nesting depth, self-intersections, slicing, repair. The
  point queries are plain floating point, as in 3-D, so they run unchanged
  on the device.

## Usage

```fortran
program naca_field
use fossil, only : curve_object, distance_field_object
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter        :: PI = acos(-1._R8P)
type(curve_object)          :: airfoil
type(distance_field_object) :: field
real(R8P)                   :: u(400), v(400), x(257), y(129), xi, t
integer(I4P)                :: k, j, st

! NACA 0012, closed trailing edge, cosine spacing, counter-clockwise
j = 0
do k = 200, 0, -1
   xi = 0.5_R8P * (1._R8P - cos(PI * k / 200)) ; t = thickness(xi)
   j = j + 1 ; u(j) = xi ; v(j) = t
enddo
do k = 1, 199
   xi = 0.5_R8P * (1._R8P - cos(PI * k / 200)) ; t = thickness(xi)
   j = j + 1 ; u(j) = xi ; v(j) = -t
enddo
call airfoil%adopt_loops(u=u, v=v, loop_sizes=[400_I4P])

print *, 'signed distance at mid-chord: ', airfoil%distance([0.5_R8P, 0._R8P])
print *, 'inside: ', airfoil%is_point_inside([0.5_R8P, 0._R8P])

! 2-D field on a stretched grid: initialize without z
do k = 1, 257
   x(k) = -1._R8P + 3._R8P * (1._R8P - cos(PI * (k - 1) / 256._R8P)) / 2._R8P
enddo
do k = 1, 129
   y(k) = tanh(2._R8P * (-1._R8P + 2._R8P * (k - 1) / 128._R8P)) / tanh(2._R8P)
enddo
call field%initialize(x=x, y=y)
call field%compute(curve=airfoil, status=st)
print *, 'band cells: ', field%get_band_cells(), '  phi(1,1): ', field%phi(1, 1, 1)

contains
   pure function thickness(s) result(h)
   real(R8P), intent(in) :: s
   real(R8P)             :: h
   h = 0.6_R8P * (0.2969_R8P * sqrt(s) - 0.1260_R8P * s - 0.3516_R8P * s**2 + 0.2843_R8P * s**3 - 0.1036_R8P * s**4)
   endfunction thickness
endprogram naca_field
```

Slicing an STL and repairing the section:

```fortran
program section
use fossil, only : surface_stl_object, curve_object, slice_surface, CURVE_VALID
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: body
type(curve_object)       :: cut
integer(I4P)             :: n_open, gaps

call body%load_from_file(file_name='body.stl', guess_format=.true.)
call slice_surface(body, 3_I4P, 0.4_R8P, cut, merge_collinear=.true., open_loops=n_open)
if (n_open > 0) call cut%repair(tolerance=1.e-6_R8P, gaps_closed=gaps)
if (cut%check_validity() == CURVE_VALID) call cut%save_into_file(file_name='section.dat')
endprogram section
```

## How it works

- **Segment BVH** (`fossil_segment_bvh`). The 2-D twin of the SAH BVH: binary,
  16 buckets per axis on box half-perimeters, with median splits below depth
  48 so the fixed traversal stack cannot overflow. Everything is stored in
  flat arrays. Queries are iterative and deterministic: the nearer child
  first, pruning on squared distance, and a candidate wins only if strictly
  closer.
- **Sign.** The sign comes from the closest feature: the side of the segment,
  or the vertex pseudo-normal. In 2-D the angle-weighted pseudo-normal of
  Bærentzen and Aanæs reduces to the bisector of the two segment normals.
  At a vertex the distance is computed from the vertex itself, so the two
  segments sharing it give the same bits.
- **Winding number.** w = (1/2π) Σ (angle subtended by each segment): 1
  inside, 0 outside or in a hole, fractional near an open end.
  `winding_number` evaluates it hierarchically, the 2-D analogue of the fast
  winding number (Barill et al. 2018). A node far enough away contributes
  the second-order expansion of (1/2π) ∫ n·(x − p)/|x − p|² ds from its
  dipole and second moments, which are exact for straight segments.
- **Slicing.**
  - A vertex on the plane counts as above it. A facet is therefore crossed
    along exactly two edges or not at all, and a facet lying in the plane is
    not crossed.
  - Crossing points are built from the edge endpoints in lexicographic
    order, so the two facets sharing an edge produce the same bits. Loops
    close by exact equality, without tolerances.
  - Each segment's direction comes from the facet's vertex order, so
    outward-oriented facets put the solid on the left.
- **Device.** The traversals are free `_oac` kernels (`!$acc routine seq`) on
  the flat arrays, and the CPU queries call them too. `distance_many_device`
  follows the [GPU rules of #20](/guide/advanced/device-distance).

## Results

| Case | Result |
|---|---|
| near-collinear `orient2d` grid (2048 points) | naive formula wrong on 482, `orient2d` on 0 |
| 2000 random points (NACA 0012, cusp, 64-gon) | BVH distance = brute force bit for bit; sign = exact winding number |
| fast winding number, NACA 0012 | inside decisions exact; max error 1.4e-2, 2.0e-3, 2.3e-4 at β = 2, 4, 8 |
| slice of the unit cube at z = 0.4, 0.5, 0.25, 1 | the exact unit square after merging, area exactly 1 |
| slice of an icosphere at z = 0.3, levels 2-4 | area error 6.9e-2, 1.8e-2, 4.7e-3 (O(h²)) |
| `distance_many` / device path vs `compute_distance` | bit for bit (4000 points, both signs) |
| 2-D field, circle, 32², 64², 128² | band bit-exact; far-field error 2.8e-2, 2.0e-2, 1.3e-2 (O(h \|log h\|)) |

## Limitations

- **GPU not validated yet.** The device path is written and runs bit-exact on
  the host, through the same kernel, in the `tests-gnu` build. The GPU check
  (`fobis build --mode tests-nvidia-curve --varset local_nvf`, then
  `fossil_test_curve_bench`) is pending.
- **Coarse fast winding number near thin bodies.** It is about 1e-2 at the
  default β = 2. The inside decision has a margin of 0.5; raise β where the
  value itself matters.
- **No curved segments.** Curves are polylines: arcs and splines must be
  sampled first.
- **Slicing with planes along axes only.** An arbitrary plane needs a
  rotation of the surface first.
- **A coarse join rule in repair.** It joins open chains greedily by nearest
  end within the tolerance. Ambiguous gaps (several ends within the
  tolerance) take the nearest.

## API

See [`curve_object`](/guide/api/curve-object) for every method and argument.
