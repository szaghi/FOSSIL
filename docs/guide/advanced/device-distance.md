---
title: GPU distance queries
---

# GPU-offloaded distance queries

::: warning Scope: SAH BVH only, pseudo-normal sign only
The device path runs **only on surfaces built with the SAH BVH**
(`AABB_TREE_SAH_BVH`, the default tree kind). An octree-backed surface
(`AABB_TREE_OCTREE`) stays CPU-only: `enter_device` refuses it with a
stderr diagnostic and leaves the surface on the host. The sign is always
the Bærentzen–Aanæs pseudo-normal (`SIGN_PSEUDO_NORMAL`); the ray-cast and
solid-angle sign algorithms have no device port. The reasons, and what an
octree port would cost, are under [Known limitations](#known-limitations).
:::

## What it does

`distance_many_device` is the GPU version of `distance_many`, the batch
query that returns one distance per point against a static surface. It is
the access pattern of an immersed-boundary (IB) marker pass on a Cartesian
grid, where every cell centre needs its signed distance to the body to
decide whether the cell is inside, outside, or cut by the surface. The loop
over query points is an OpenACC `parallel loop`, one GPU thread per point,
and each thread runs the same best-first BVH traversal and pseudo-normal
sign test as the CPU path.

The surface's search structures (the flat BVH node array, the flat
triangle payload, and the per-facet pseudo-normals) are copied to the GPU
**once**, by `enter_device`, and stay there for as many query calls as the
caller makes. Each `distance_many_device` call moves only the query points
in and the distances out.

The results are **bit-exact** against the serial CPU path on every case
validated so far (nvfortran 26.9, RTX 4070, see [Performance](#performance)).

## Pipeline

1. **Residency** (`enter_device`, once per surface). `!$acc enter data
   copyin` of three flat arrays: the BVH nodes (`aabb%node`), the triangle
   metrics (`aabb%payload`), and the per-facet pseudo-normals
   (`pnormal_payload`, 7 vectors per facet, built by `analyze`). Only those
   member arrays are mapped, never the `surface_stl_object` itself.
2. **Launch** (`distance_many_device`, per batch). A host-side check refuses
   the call if the surface is not resident; then one `parallel loop` over
   the points, with `present(...)` on the three arrays, `copyin(points)` and
   `copyout(distances)`.
3. **Traversal** (per thread). `bvh_distance_node_with_region_oac`: the
   same best-first descent with d² pruning as on the CPU, written with an
   explicit stack instead of recursion (recursion inside an `!$acc routine`
   costs per-thread call stacks). It returns the best d², the winning facet
   and its payload index.
4. **Sign** (per thread, only with `is_signed=.true.`). Recompute the closest
   point and its Voronoi region (face, edge or vertex) on the winning
   triangle, take the matching pseudo-normal from `pnormal_payload`, and
   sign the distance by the dot product with `point - closest`.
5. **Release** (`exit_device`, once). `!$acc exit data delete` of the three
   arrays.

## API

### Make the surface resident

```fortran
call surface%enter_device()
   class(surface_stl_object), intent(inout) :: self
```

Call it after the surface is final: after `load_from_file`, `sanitize` and
any transform. The call is idempotent (a second call is a no-op) and returns
without acting on an empty surface. On an octree surface it writes
`enter_device: refused -- device path requires AABB_TREE_SAH_BVH.` to
stderr; a later `distance_many_device` then refuses too (see below).

### Query

```fortran
call surface%distance_many_device(points, distances, is_signed, is_square_root, status)
   class(surface_stl_object), intent(inout), target   :: self
   type(vector_R8P),          intent(in)              :: points(:)
   real(R8P),                 intent(out)             :: distances(:)  ! caller-allocated, size(points)
   logical,                   intent(in),  optional   :: is_signed       ! default .false.
   logical,                   intent(in),  optional   :: is_square_root  ! default .false. (d²)
   integer(I4P),              intent(out), optional   :: status          ! STATUS_OK / STATUS_INVALID_INPUT
```

The arguments mirror `distance_many`, minus `sign_algorithm`, so a caller
switches between host and device by changing the call name. As for
`distance_many`, the default result is the **squared** distance; pass
`is_square_root=.true.` for the Euclidean one.

`status = STATUS_INVALID_INPUT` when `size(distances) /= size(points)`, or
when the surface is not resident (the call also writes `surface not on
device -- call enter_device first.` to stderr). `distances` is left
untouched in both cases.

### Release

```fortran
call surface%exit_device()
   class(surface_stl_object), intent(inout) :: self
```

Idempotent. **`destroy` does not call it**: `destroy` is `elemental`, and
`!$acc exit data` is not allowed in an `elemental` procedure. Call
`exit_device` yourself before the surface is destroyed or goes out of scope,
otherwise the device allocations stay registered against host addresses that
no longer exist.

### Building

The device code is compiled only by nvfortran with `-acc`, through the
`static-nvidia` mode (library) or `tests-nvidia` mode (the distance
benchmark):

```bash
fobis build --mode static-nvidia --varset local_nvf
```

The `local_nvf` varset sets the GPU compute capability (`$NVF_CC`, default
`cc89`; use `cc80` for A100). nvfortran **26.9 or newer** is required. See
[Installation](/guide/installation).

Under every other build (gfortran, Intel) the OpenACC directives are inert
comments: the three procedures still exist, `enter_device` only sets the
residency flag, and `distance_many_device` runs its loop **serially on the
host** with the same results. Code written against the device API therefore
builds and runs everywhere, and only goes faster where nvfortran is used.

## Example

A minimal IB marker: classify the cell centres of a Cartesian block around
the body as inside, outside or cut, with one device call for the whole block.
A cell counts as cut when the surface passes within half a cell diagonal of
its centre.

```fortran
program ex_device_ib_marker
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none

integer(I4P), parameter :: NC = 32_I4P            ! cells per direction
type(surface_stl_object)      :: body
type(vector_R8P)              :: bmin, bmax, h
type(vector_R8P), allocatable :: centre(:)
real(R8P),        allocatable :: d(:)
real(R8P)                     :: half_diag
integer(I4P)                  :: i, j, k, c, status
integer(I4P)                  :: n_inside, n_outside, n_cut

call body%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call body%sanitize  ! the pseudo-normal sign needs outward normals

! Block = body bounding box grown by 25% on each side.
bmin = body%get_bmin()
bmax = body%get_bmax()
h    = bmax - bmin
bmin = bmin - 0.25_R8P * h
bmax = bmax + 0.25_R8P * h
h    = (bmax - bmin) / real(NC, R8P)
half_diag = 0.5_R8P * sqrt(h%x**2 + h%y**2 + h%z**2)

allocate(centre(NC**3), d(NC**3))
c = 0
do k = 1, NC
   do j = 1, NC
      do i = 1, NC
         c = c + 1
         centre(c)%x = bmin%x + (real(i, R8P) - 0.5_R8P) * h%x
         centre(c)%y = bmin%y + (real(j, R8P) - 0.5_R8P) * h%y
         centre(c)%z = bmin%z + (real(k, R8P) - 0.5_R8P) * h%z
      enddo
   enddo
enddo

call body%enter_device()                         ! copy the BVH to the GPU once
call body%distance_many_device(points=centre, distances=d, is_signed=.true., &
                               is_square_root=.true., status=status)
call body%exit_device()                          ! release before the surface dies
if (status /= STATUS_OK) error stop 'distance_many_device failed'

n_cut     = count(abs(d) <= half_diag)
n_inside  = count(d < -half_diag)
n_outside = count(d >  half_diag)
print '(A,3(I0,A))', 'inside: ', n_inside, '  outside: ', n_outside, '  cut: ', n_cut, ''
endprogram ex_device_ib_marker
```

In a real marker pass the surface stays resident across many calls (one per
AMR block or per time step on a moving grid), and `exit_device` comes at the
end of the run. The `enter_device` cost is paid once.

## Performance

Measured with `fossil_test_distance_bench` (`tests-nvidia`, nvfortran 26.9,
NVIDIA RTX 4070, `cc89`), signed Euclidean distance, query points on an
N³ grid. Baseline: the serial CPU loop on **one core**.

| Mesh | Query grid | Speedup vs 1 CPU core | max \|d_device − d_cpu\| |
|---|---|---|---|
| bunny | 16³ | 3.05× | 0 |
| bunny | 32³ | ~8.8× (7.1–8.9× over 5 runs) | 0 |
| bunny | 48³ | 12.6× | 0 |
| cube-01 | 32³ | 4.77× | 0 |
| dragon-fine (24k facets) | 32³ | 8.6× | 0 |

How to read these numbers:

- The device time is **one call** with no warm-up: it includes the
  point/distance transfers and the launch, and excludes `enter_device`.
- The baseline is one core. Against `distance_many` built with OpenMP on a
  multi-core CPU the GPU margin is smaller, and on small batches it can
  vanish: the speedup grows with the number of query points because the
  fixed launch and transfer cost is amortized.
- The RTX 4070 is a consumer GPU, with an FP64 rate of 1/64 of FP32. The
  kernel is FP64 throughout, so a data-centre GPU (A100, H100) should do
  better. That has not been measured.

## Known limitations

- **SAH BVH only.** The octree is not supported on the device. The original
  obstacle is gone: the octree allocates all 8 child slots per level but
  populates only the refined ones, and the "empty slot, skip it" signal used
  to be `allocated(node%aabb) == .false.`. Making the node device-mappable
  stored the AABB by value, so that signal is now an explicit
  `is_present` flag on `aabb_node_object` (`is_allocated()` returns it, and
  every octree call site reads it unchanged).

  **What an octree port still needs:**
  1. A flat payload for octree leaves. BVH leaves own a contiguous slice of
     the `facet_distance_payload` array (`payload_first` / `payload_count`);
     octree leaves have `payload_count = 0` and keep their facet ids in an
     allocatable `list_id_object`, which is not mapped to the device.
  2. A device traversal over 8 implicit children
     (`first_child = 8 * parent + 1`, skipping slots with `is_present`
     false) next to the 2-child BVH one.
  3. `enter_device` accepting `AABB_TREE_OCTREE`, and the benchmark gate
     run on an octree surface.

  The work is left out on purpose: the octree is opt-in, the default is the
  BVH, and on dragon-fine the octree is about 110× slower than the BVH even
  on the CPU, so a device octree would have no user.
- **Pseudo-normal sign only.** `SIGN_RAY_INTERSECTIONS` and
  `SIGN_SOLID_ANGLE` need their own traversals, which have not been
  ported. Like the CPU default, the pseudo-normal sign needs consistently
  outward normals: call `sanitize` before `enter_device`.
- **Distances only.** `compute_distance`'s richer outputs (closest facet,
  edge, vertex), ray queries and winding numbers have no device version.
- **No automatic coherence.** The device copy is a snapshot taken by
  `enter_device`, and nothing checks that it is still current. Operations
  that rebuild the tree (`analyze` and the procedures that call it:
  `sanitize`, `clip`, `merge_solids`, `adopt_facets`, `load_from_file`)
  reallocate the host arrays behind the device mapping. The transforms
  (`translate`, `rotate`, `mirror`, `resize`) do not rebuild the tree at
  all, so after them the tree is stale on the host as well, until `analyze`
  runs. In both cases: `exit_device`, change the surface, `analyze` if
  needed, `enter_device` again.
- **`destroy` does not release the device copy** (see
  [Release](#release)).
- **Traversal stack depth.** The explicit traversal stack holds 128
  entries (`MAX_TRAVERSAL_STACK`). A binary BVH needs about one entry per
  tree level. Measured depths are small: 12 levels on bunny (69451 facets),
  13 on dragon-fine (24318 facets). But the builder does not cap the depth
  and the push is not bounds-checked, so a pathologically skewed tree
  deeper than about 127 levels would write past the stack. This applies to
  the CPU path too, which uses the same traversal.
- **Validated on one GPU.** Bit-exactness and the timings above come from
  one consumer GPU (RTX 4070) with nvfortran 26.9. The device build is not
  part of CI (no GPU runner): it is validated manually with
  `fossil_test_distance_bench`, which fails unless the device result
  matches the serial CPU result bit for bit. Bit-exactness is not
  guaranteed by construction: a different GPU or compiler version could
  contract a multiply-add differently and move the last bits.

## See also

- [Distance and point-in-polyhedron queries](/guide/features#distance-and-point-in-polyhedron-queries):
  the CPU `distance` / `distance_many` API and the sign algorithms.
- [AABB acceleration](/guide/features#aabb-acceleration): the BVH and octree tree kinds.
- [Installation](/guide/installation): the `static-nvidia` / `tests-nvidia`
  modes and the nvfortran version requirement.
- [Issue #20](https://github.com/szaghi/FOSSIL/issues/20): the
  implementation plan, the device-code rules and the validation history.

## References

- Bærentzen & Aanæs, *Signed Distance Computation Using the Angle Weighted
  Pseudonormal*, IEEE TVCG 11(3), 2005. The sign test.
- Wald, *On Fast Construction of SAH-based Bounding Volume Hierarchies*,
  IEEE Symposium on Interactive Ray Tracing, 2007. The bucketed SAH build.
- OpenACC-Standard.org, *The OpenACC Application Programming Interface*,
  version 3.4, 2025. Data clauses (§2.6–2.7) and `routine` (§2.15).
