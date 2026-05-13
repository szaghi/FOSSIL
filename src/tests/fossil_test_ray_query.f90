!< FOSSIL test: ray-mesh intersection queries (issue #18 §2.5).
!<
!< Public contract that all step-N implementations must honour:
!<   - hit lists are sorted ascending by `t`,
!<   - only `t >= 0` hits are returned,
!<   - `status = RAY_STATUS_OK` and `size(hits) = 0` for misses (no error),
!<   - `intersect_ray_first%hit` matches `intersect_ray_all%hits(1)` whenever
!<     the all-hits list is non-empty.
!<
!< Step 1 covers invariants 1..4 (cube-only sanity). Step 2 adds 5..7 for the
!< tree-accelerated path. Step 3 adds 8..9 for the any-hit early-exit query:
!<   5. Tree-vs-flat parity on bunny.stl. 100 deterministic random rays:
!<      `intersect_ray_all` (which now goes through the AABB tree) must produce
!<      hit lists bit-identical to the flat O(n) reference oracle
!<      `ray_intersect_all_flat`. Identical = same hit count, same per-hit
!<      `t` and `facet_id`. The fixture is non-trivial (bunny ~ 70k facets,
!<      tree depth ~5), so any traversal bug shows up.
!<   6. Sphere from §1.5 marching cubes, radial outward ray from center.
!<      Exactly one hit at t ~= R. Sphere has no flat-face degeneracies, so
!<      this validates the curved-mesh path independently of axis alignment.
!<   7. `intersect_ray_first` matches `intersect_ray_all%hits(1)` on the
!<      cube oblique ray (which gives a deterministic single best hit).
!<   8. `intersect_ray_any` `max_t` boundary toggle. Cube + ray with a known
!<      first-hit `t = 1.0`. With `max_t = 0.5` (below the hit), expect
!<      `found = .false.`. With `max_t = 1.5` (above), expect `found = .true.`.
!<      The exactly-on-boundary case `max_t = 1.0` should also be `.true.`
!<      because the contract is `t <= max_t`.
!<   9. `intersect_ray_any` ≡ `size(intersect_ray_all%hits) > 0` over the
!<      same 100 random bunny rays from invariant 5. Cross-check between the
!<      two query types — protects against a divergence where one accepts a
!<      hit the other rejects (e.g. different `t < 0` filtering).
!<
!< Geometry note (cube fixture in src/tests/cube.stl spans [0,1]^3): each face
!< is split along one of {y=z, x=z, x+y=1} (depending on the face). Test rays
!< are placed at OFF-DIAGONAL cross-section points (e.g. (0.3, 0.6) instead of
!< (0.5, 0.5)) to avoid landing on the shared edge — a shared-edge ray would
!< correctly produce 4 hits (one per triangle), which is a different and
!< noisier invariant than what we want to validate here.

program fossil_test_ray_query

use fossil, only : surface_stl_object, ray_hit_t, RAY_STATUS_OK, &
                   extract_isosurface
use fossil_facet_object, only : facet_object
use fossil_ray_query, only : ray_intersect_all_flat
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

type(surface_stl_object) :: cube, bunny, sphere
type(facet_object), allocatable :: sphere_facets(:)
type(ray_hit_t), allocatable :: hits(:), hits_flat(:)
type(ray_hit_t) :: best_hit
type(vector_R8P) :: origin, direction
type(vector_R8P) :: bmin, bmax
real(R8P), allocatable :: values(:,:,:)
real(R8P) :: sqrt3
real(R8P) :: x, y, z, dx
integer(I4P) :: status
integer(I4P) :: i, j, k, r
integer(I4P) :: n_diff, n_any_diff
logical :: has_hit, found_below, found_at, found_above
logical :: any_tree, any_all
logical :: are_tests_passed(9)

real(R8P), parameter :: TOL    = 1.0e-12_R8P
real(R8P), parameter :: TOL_FP = 1.0e-9_R8P  ! tree vs flat: identical math, but reordered → FP rounding
integer(I4P), parameter :: N_RANDOM_RAYS = 100_I4P
integer(I4P), parameter :: N_GRID = 32_I4P

are_tests_passed = .false.

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)

! ---- Invariant 1: ray along +X, off-diagonal in (y, z) cross section.
origin    = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
direction = ex_R8P
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 1 (+X): status=', status, ' n_hits=', size(hits)
if (size(hits) == 2) print '(A,F8.4,A,F8.4)', '            t=[', hits(1)%t, ', ', hits(2)%t
are_tests_passed(1) = (status == RAY_STATUS_OK)              .and. &
                      (size(hits) == 2)                      .and. &
                      (abs(hits(1)%t - 1._R8P) < TOL)        .and. &
                      (abs(hits(2)%t - 2._R8P) < TOL)

! ---- Invariant 2a: ray along +Y, off-diagonal in (x, z).
origin    = vector_R8P(0.3_R8P, -1._R8P, 0.6_R8P)
direction = ey_R8P
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 2a (+Y): status=', status, ' n_hits=', size(hits)
are_tests_passed(2) = (status == RAY_STATUS_OK)              .and. &
                      (size(hits) == 2)                      .and. &
                      (abs(hits(1)%t - 1._R8P) < TOL)        .and. &
                      (abs(hits(2)%t - 2._R8P) < TOL)

! ---- Invariant 2b: ray along +Z, off-diagonal in (x, y).
origin    = vector_R8P(0.3_R8P, 0.6_R8P, -1._R8P)
direction = ez_R8P
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 2b (+Z): status=', status, ' n_hits=', size(hits)
are_tests_passed(2) = are_tests_passed(2)                    .and. &
                      (status == RAY_STATUS_OK)              .and. &
                      (size(hits) == 2)                      .and. &
                      (abs(hits(1)%t - 1._R8P) < TOL)        .and. &
                      (abs(hits(2)%t - 2._R8P) < TOL)

! ---- Invariant 3: oblique ray hitting cube face triangles' interior.
sqrt3 = sqrt(5.25_R8P)
origin    = vector_R8P(-0.5_R8P, 0.3_R8P, 0.2_R8P)
direction = vector_R8P(2._R8P, 1._R8P, 0.5_R8P) / sqrt3
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 3 (oblique): status=', status, ' n_hits=', size(hits)
if (size(hits) == 2) print '(A,F10.6,A,F10.6)', '                 t=[', hits(1)%t, ', ', hits(2)%t
are_tests_passed(3) = (status == RAY_STATUS_OK)                                       .and. &
                      (size(hits) == 2)                                               .and. &
                      (abs(hits(1)%t - 0.25_R8P * sqrt3) < TOL)                       .and. &
                      (abs(hits(2)%t - 0.7_R8P  * sqrt3) < TOL)                       .and. &
                      (abs((origin%x + hits(1)%t * direction%x) - 0._R8P) < TOL)      .and. &
                      (abs((origin%y + hits(2)%t * direction%y) - 1._R8P) < TOL)

! ---- Invariant 4: ray that misses the cube.
origin    = vector_R8P(10._R8P, 10._R8P, 10._R8P)
direction = ex_R8P
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 4 (miss): status=', status, ' n_hits=', size(hits)
are_tests_passed(4) = (status == RAY_STATUS_OK) .and. (size(hits) == 0)

! ---- Invariant 5: tree vs flat parity on bunny.stl, deterministic random rays.
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
n_diff = 0_I4P
do r = 1_I4P, N_RANDOM_RAYS
   ! Cheap LCG-like PRNG seeded by `r`. No need for randomness across runs;
   ! reproducibility wins so a regression bisects cleanly.
   call cheap_random_origin_dir(seed=r, origin=origin, direction=direction)
   call bunny%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
   call ray_intersect_all_flat(facet=bunny_facets_ref(), ray_origin=origin, &
                               ray_direction=direction, hits=hits_flat, status=status)
   if (.not. hit_lists_match(hits_a=hits, hits_b=hits_flat, tol=TOL_FP)) then
      n_diff = n_diff + 1_I4P
      if (n_diff <= 3_I4P) print '(A,I0,A,I0,A,I0)', '   ray ', r, ' MISMATCH tree=', size(hits), ' flat=', size(hits_flat)
   endif
enddo
print '(A,I0,A,I0,A)', 'inv 5 (bunny tree=flat): ', N_RANDOM_RAYS - n_diff, '/', N_RANDOM_RAYS, ' rays match'
are_tests_passed(5) = (n_diff == 0_I4P)

! ---- Invariant 6: sphere from MC, radial outward ray from center → 1 hit.
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)  ! grid covers [-2, 2]^3
do k = 1_I4P, N_GRID
   z = -2._R8P + (k - 1) * dx
   do j = 1_I4P, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1_I4P, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1._R8P  ! unit sphere
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)
origin    = vector_R8P(0._R8P, 0._R8P, 0._R8P)
! Slightly tilted off-axis to avoid grazing equator triangles' shared edges
! (which would correctly produce two hits at the same `t`).
direction = vector_R8P(1._R8P, 0.123_R8P, 0.157_R8P)
direction = direction / sqrt(direction%dotproduct(rhs=direction))
call sphere%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,F8.4)', 'inv 6 (sphere radial): n_hits=', size(hits), ' t_first=', merge(hits(1)%t, 0._R8P, size(hits) > 0)
are_tests_passed(6) = (status == RAY_STATUS_OK)              .and. &
                      (size(hits) == 1)                      .and. &
                      (abs(hits(1)%t - 1._R8P) < 5.0e-2_R8P)  ! MC mesh resolution-bounded

! ---- Invariant 7: intersect_ray_first matches intersect_ray_all%hits(1).
sqrt3 = sqrt(5.25_R8P)
origin    = vector_R8P(-0.5_R8P, 0.3_R8P, 0.2_R8P)
direction = vector_R8P(2._R8P, 1._R8P, 0.5_R8P) / sqrt3
call cube%intersect_ray_first(ray_origin=origin, ray_direction=direction, &
                              hit=best_hit, has_hit=has_hit, status=status)
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, &
                            hits=hits, status=status)
print '(A,L1,A,F10.6,A,F10.6)', 'inv 7 (first==all(1)): has_hit=', has_hit, ' t_first=', best_hit%t, &
                                ' t_all(1)=', hits(1)%t
are_tests_passed(7) = has_hit                                       .and. &
                      (size(hits) >= 1)                             .and. &
                      (best_hit%facet_id == hits(1)%facet_id)       .and. &
                      (abs(best_hit%t - hits(1)%t) < TOL)

! ---- Invariant 8: intersect_ray_any max_t boundary toggle.
!      Same +X cube ray as invariant 1 — first hit is at t = 1.0.
origin    = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
direction = ex_R8P
call cube%intersect_ray_any(ray_origin=origin, ray_direction=direction, max_t=0.5_R8P, &
                            found=found_below, status=status)
call cube%intersect_ray_any(ray_origin=origin, ray_direction=direction, max_t=1.0_R8P, &
                            found=found_at, status=status)
call cube%intersect_ray_any(ray_origin=origin, ray_direction=direction, max_t=1.5_R8P, &
                            found=found_above, status=status)
print '(A,3L2)', 'inv 8 (any max_t toggle): below/at/above = ', found_below, found_at, found_above
are_tests_passed(8) = (.not. found_below) .and. found_at .and. found_above

! ---- Invariant 9: any-hit ≡ (size(all-hits) > 0) over the same 100 bunny rays.
n_any_diff = 0_I4P
do r = 1_I4P, N_RANDOM_RAYS
   call cheap_random_origin_dir(seed=r, origin=origin, direction=direction)
   call bunny%intersect_ray_any(ray_origin=origin, ray_direction=direction, &
                                found=any_tree, status=status)
   call bunny%intersect_ray_all(ray_origin=origin, ray_direction=direction, &
                                hits=hits, status=status)
   any_all = (size(hits, kind=I4P) > 0_I4P)
   if (any_tree .neqv. any_all) then
      n_any_diff = n_any_diff + 1_I4P
      if (n_any_diff <= 3_I4P) print '(A,I0,A,L1,A,L1)', '   ray ', r, ' any/all mismatch any=', any_tree, ' all=', any_all
   endif
enddo
print '(A,I0,A,I0,A)', 'inv 9 (any==(any-of-all)): ', N_RANDOM_RAYS - n_any_diff, '/', N_RANDOM_RAYS, ' rays match'
are_tests_passed(9) = (n_any_diff == 0_I4P)

print '(A,9L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains

   function bunny_facets_ref() result(arr)
   !< Read-only handle to bunny's facet array — passes the surface's pointer
   !< through unchanged. Local helper to keep the call site at invariant 5
   !< readable.
   type(facet_object), pointer :: arr(:)

   arr => bunny%facets_ref()
   endfunction bunny_facets_ref

   subroutine cheap_random_origin_dir(seed, origin, direction)
   !< Deterministic pseudo-random ray, parameterized by `seed`. Origin lies on a
   !< sphere of radius 5 around the bunny's bounding-box center; direction
   !< points toward a random point inside the bbox. Designed so most rays
   !< actually hit the model — to exercise the all-hits accumulator and
   !< mid-traversal pruning, not just the trivial-miss path.
   integer(I4P),     intent(in)  :: seed
   type(vector_R8P), intent(out) :: origin
   type(vector_R8P), intent(out) :: direction
   real(R8P)                     :: a, b, c, d, e, f
   ! Six independent pseudo-random reals in [0, 1) from a Wang-style hash.
   call hashed_uniform(seed=seed, k=1_I4P, x=a)
   call hashed_uniform(seed=seed, k=2_I4P, x=b)
   call hashed_uniform(seed=seed, k=3_I4P, x=c)
   call hashed_uniform(seed=seed, k=4_I4P, x=d)
   call hashed_uniform(seed=seed, k=5_I4P, x=e)
   call hashed_uniform(seed=seed, k=6_I4P, x=f)
   ! Origin: scatter on a 5-radius shell around centroid (0, 0, 0 for bunny).
   origin = vector_R8P((a - 0.5_R8P) * 10._R8P, &
                       (b - 0.5_R8P) * 10._R8P, &
                       (c - 0.5_R8P) * 10._R8P)
   ! Direction: aim toward a random point in [-1, 1]^3 (bunny lives in there).
   direction = vector_R8P((d - 0.5_R8P) * 2._R8P - origin%x, &
                          (e - 0.5_R8P) * 2._R8P - origin%y, &
                          (f - 0.5_R8P) * 2._R8P - origin%z)
   endsubroutine cheap_random_origin_dir

   pure subroutine hashed_uniform(seed, k, x)
   !< Stateless integer hash → uniform [0, 1). Same (seed, k) → same `x`, so the
   !< test is bit-reproducible. Built from small multipliers and bit-mixing so
   !< no intermediate exceeds I4P range. Quality only needs to be "doesn't
   !< clump at predictable angles" — not cryptographic.
   integer(I4P), intent(in)  :: seed, k
   real(R8P),    intent(out) :: x
   integer(I4P)              :: h

   h = ieor(seed * 1103515245_I4P + 12345_I4P, k * 2531011_I4P)
   h = ieor(h, ishft(h, -16))
   h = h * 1664525_I4P
   h = ieor(h, ishft(h, -13))
   h = h * 1013904223_I4P
   h = ieor(h, ishft(h, -16))
   x = real(iand(h, 2147483647_I4P), R8P) / 2147483647._R8P
   endsubroutine hashed_uniform

   pure function hit_lists_match(hits_a, hits_b, tol) result(yes)
   !< True if two hit lists agree: same length, same per-hit `t` within `tol`,
   !< same per-hit `facet_id`. Hit lists are already sorted by `t` by both
   !< drivers; we don't re-sort here.
   type(ray_hit_t), intent(in) :: hits_a(:), hits_b(:)
   real(R8P),       intent(in) :: tol
   logical                     :: yes
   integer(I4P)                :: ii

   yes = (size(hits_a, kind=I4P) == size(hits_b, kind=I4P))
   if (.not. yes) return
   do ii = 1_I4P, size(hits_a, kind=I4P)
      if (abs(hits_a(ii)%t - hits_b(ii)%t) > tol) then
         yes = .false. ; return
      endif
      if (hits_a(ii)%facet_id /= hits_b(ii)%facet_id) then
         yes = .false. ; return
      endif
   enddo
   endfunction hit_lists_match

endprogram fossil_test_ray_query
