---
title: Ray queries (§2.5)
---

# §2.5 — Ray-mesh intersection queries

## What it does

Three TBPs that wrap the AABB-tree-accelerated **ray-versus-surface**
intersection in flavours that match the canonical use cases:

- **`intersect_ray_first`** — closest hit (best-first traversal with
  `t_best` pruning). For *picking* (which facet did the ray hit?).
- **`intersect_ray_all`** — every hit, sorted ascending by `t`. For
  *visibility* counting and inside / outside ray-cast tests.
- **`intersect_ray_any`** — early-exit "is there *any* hit within
  `max_t`?" For *shadow rays* and conservative inside-test
  rejection.

All three share the same Möller-Trumbore facet-test
(`facet%intersect_ray`) and the same AABB tree traversal pattern (the
`enumerate_children` dispatcher used by §1.4 and §1.2). Both SAH BVH
and octree tree kinds are supported transparently.

## Pipeline

For each query:

1. **Tree traversal**. Starting from the root, walk down to leaves
   whose AABB the ray pierces. The slab-test produces a `[t_near,
   t_far]` interval per box; subtree pruning uses `t_near` for
   ordering.
2. **Per-leaf facet sweep**. For each owned facet, run Möller-
   Trumbore (`facet%intersect_ray`) which returns the parametric `t`
   plus barycentrics `(u, v)`.
3. **Per-flavour assembly**:
   - `_first` — keep the smallest `t ≥ 0`; prune any subtree whose
     `t_near ≥ t_best`.
   - `_all` — accumulate every hit with `t ≥ 0`; sort ascending.
   - `_any` — return on the first hit with `t ∈ [0, max_t]`; prune any
     subtree whose `t_near > max_t`.

A flat-O(N) fallback exists for surfaces without a built AABB tree;
it operates on `facet(:)` directly. Standard surfaces from
`load_from_file` always have a tree, so the fallback rarely fires.

## API

### Closest hit

```fortran
call surface%intersect_ray_first(ray_origin, ray_direction, hit, has_hit, status)
   class(surface_stl_object), intent(in)            :: self
   type(vector_R8P),          intent(in)            :: ray_origin
   type(vector_R8P),          intent(in)            :: ray_direction  ! need not be unit
   type(ray_hit_t),           intent(out)           :: hit            ! valid iff has_hit
   logical,                   intent(out)           :: has_hit
   integer(I4P),              intent(out), optional :: status         ! RAY_STATUS_*
```

`hit` carries `%facet_id`, `%t`, `%point` (= `ray_origin + t·ray_direction`).
**Always check `has_hit` first** — `hit` is undefined when `has_hit
= .false.`.

### All hits

```fortran
call surface%intersect_ray_all(ray_origin, ray_direction, hits, status)
   type(ray_hit_t), allocatable, intent(out) :: hits(:)
```

`hits` is **always allocated** (zero-length if no hits found), and
**sorted ascending by `t`**.

### Any hit (shadow ray)

```fortran
call surface%intersect_ray_any(ray_origin, ray_direction, max_t, found, status)
   real(R8P), intent(in), optional :: max_t      ! default MaxR8P
   logical,   intent(out)          :: found
```

Returns immediately on the first hit with `0 ≤ t ≤ max_t`. Default
`max_t` accepts any forward hit. Pass the Euclidean distance to a
target point (when direction is unit) to test occlusion of that
target.

## Examples

### Picking with `_first`

```fortran
program ex_intersect_ray_first
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none

type(surface_stl_object) :: cube
type(ray_hit_t)          :: hit
type(vector_R8P)         :: origin, direction
logical                  :: has_hit
integer(I4P)             :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
origin    = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
direction = ex_R8P
call cube%intersect_ray_first(ray_origin=origin, ray_direction=direction, &
                              hit=hit, has_hit=has_hit, status=status)
if (has_hit) then
   print '(A,I0,A,F8.4)', 'first hit on facet ', hit%facet_id, ' at t = ', hit%t
else
   print '(A)', 'ray missed the cube'
endif
endprogram ex_intersect_ray_first
```

### Counting hits with `_all`

```fortran
program ex_intersect_ray_all
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none

type(surface_stl_object)     :: cube
type(ray_hit_t), allocatable :: hits(:)
type(vector_R8P)             :: origin, direction
integer(I4P)                 :: status, i

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
origin    = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
direction = ex_R8P
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, &
                            hits=hits, status=status)
print '(A,I0)', 'n_hits = ', size(hits)
do i = 1_I4P, size(hits, kind=I4P)
   print '(A,I0,A,I0,A,F8.4)', '  hit ', i, ': facet=', hits(i)%facet_id, '  t=', hits(i)%t
enddo
endprogram ex_intersect_ray_all
```

### Shadow-ray with `_any`

```fortran
program ex_intersect_ray_any
use fossil
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P
implicit none

type(surface_stl_object) :: cube
type(vector_R8P)         :: origin, direction
logical                  :: blocked
integer(I4P)             :: status

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
origin    = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
direction = ex_R8P

! Is anything blocking the ray within distance 0.5? (cube is at x=0, so no)
call cube%intersect_ray_any(ray_origin=origin, ray_direction=direction, &
                            max_t=0.5_R8P, found=blocked, status=status)
print '(A,L1)', 'blocked within 0.5? ', blocked

! Within distance 1.5? (cube enters at t=1.0, so yes)
call cube%intersect_ray_any(ray_origin=origin, ray_direction=direction, &
                            max_t=1.5_R8P, found=blocked, status=status)
print '(A,L1)', 'blocked within 1.5? ', blocked
endprogram ex_intersect_ray_any
```

## Known limitations

- **Direction need not be unit.** `t` is the *parametric* position
  (`origin + t·direction`), not the Euclidean distance, when
  `direction` isn't unit. Normalize the direction yourself if you
  want `t` to mean Euclidean distance.
- **Hits exactly at `t = 0`** (origin lies on the surface) are
  filtered out of `_first` and `_all` (they'd produce false self-
  hits). For ray casts whose origin must be exactly on the surface,
  offset the origin by a small `ε` along the direction first.
- **Boundary-touching cases at `t = max_t`** are *included* by
  `_any` (the contract is `t ≤ max_t`, inclusive).
- **Two latent slab-test bugs** in the underlying AABB ray test were
  fixed during the §2.5 implementation: a parallel-ray classifier
  that misread `EPS = 0.0` as "negative directions are parallel," and
  a slab-interval accumulator that took `max` of far-plane crossings
  instead of `min`. The pre-§2.5 ray-cast inside-test sign algorithm
  worked anyway — the bugs only inflated AABB-hit counts; the per-
  facet test downstream filtered correctly. The fix made the inside-
  test sign algorithm strictly faster and unblocked the new
  intersect_ray_* TBPs.

## See also

- [`facet%intersect_ray`](/guide/api/facet-object) — the per-facet
  Möller-Trumbore that the AABB traversal calls into. Useful if you
  want to test against a single triangle without building a surface.
- [§1.9 SDF segmentation](/guide/advanced/sdf-segmentation) — the
  cone-of-rays SDF query uses `intersect_ray_first` for every cone
  ray.
- [`ray_hit_t`](/guide/api/constants) — the record type returned by
  the queries.
- [Constants](/guide/api/constants) — `RAY_STATUS_*`.

## References

- Möller & Trumbore, *Fast, Minimum Storage Ray-Triangle
  Intersection*, Journal of Graphics Tools 2(1), 1997. The per-facet
  test.
- Smits, *Efficient Bounding Box Intersection*, Ray Tracing News
  15(1), 2002. The slab-based AABB ray test.
- Williams, Barrus, Morley & Shirley, *An Efficient and Robust
  Ray-Box Intersection Algorithm*, J. of Graphics Tools 10(1), 2005.
  Refinement of Smits' slab test that handles the parallel-ray
  case correctly (the formulation FOSSIL adopts after the §2.5
  fix).
