!< FOSSIL test: ray-mesh intersection queries (issue #18 §2.5, step 1).
!<
!< Step-1 scope: validate the public contract of `surface%intersect_ray_all`
!< on the flat O(n) implementation. The same contract — sorted ascending by
!< `t`, `t >= 0` only, deterministic ordering — must survive the tree-
!< accelerated rewrite in step 2 unchanged.
!<
!< Four invariants:
!<   The cube fixture in src/tests/cube.stl spans [0,1]^3 (12 facets, 2 per
!<   axis-aligned face). Each face is split along its (y=z, x=z, x=y) diagonal,
!<   so test rays are placed at OFF-DIAGONAL cross-section points (e.g. (0.3,
!<   0.7) instead of (0.5, 0.5)) to avoid landing on the shared edge — a
!<   shared-edge ray would correctly produce 4 hits (one per triangle), which
!<   is a different and noisier invariant than what we want to validate here.
!<
!<   1. Cube + axis-aligned ray (+X). Ray from (-1, 0.3, 0.7) along +X. Exactly
!<      two hits at t = 1.0 (entry face x = 0) and t = 2.0 (exit face x = 1),
!<      in that order. Validates: hit count, t values, sort ordering.
!<   2. Same axis test along +Y and +Z (independent rays, off-diagonal cross
!<      section). Validates that the result is not axis-biased — a cube
!<      symmetry test.
!<   3. Oblique ray through cube interior. From (-0.5, -0.5, -0.5) along
!<      (1, 1, 1)/sqrt(3) — chosen so the entry point is strictly interior to
!<      a face triangle (not on a vertex or shared edge). Tests non-axis-
!<      aligned direction and hit-point reconstruction.
!<   4. Miss-the-bbox: ray far from cube, direction not aimed at the cube.
!<      Status = RAY_STATUS_OK and size(hits) = 0 — the empty-result case
!<      must not be misreported as an error.

program fossil_test_ray_query

use fossil, only : surface_stl_object, ray_hit_t, RAY_STATUS_OK
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

type(surface_stl_object) :: cube
type(ray_hit_t), allocatable :: hits(:)
type(vector_R8P) :: origin, direction
real(R8P) :: sqrt3
integer(I4P) :: status
logical :: are_tests_passed(4)

real(R8P), parameter :: TOL = 1.0e-12_R8P

are_tests_passed = .false.

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)

! ---- Invariant 1: ray along +X, off-diagonal in (y, z) cross section.
!      Each face of the unit cube is split along one of {y=z, x=z, x+y=1}
!      (depending on the face). The point (0.3, 0.6) avoids all three
!      diagonals when used as a (y,z), (x,z), or (x,y) cross section.
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
!      Origin (-0.5, 0.3, 0.2), direction (2, 1, 0.5)/sqrt(5.25). Geometry:
!        - Enters -X face (x=0) at parameter s=0.25, position (0, 0.55, 0.325).
!          Off-diagonal (y+z=0.875, not 1), so strictly interior to one triangle.
!        - Exits +Y face (y=1) at s=0.7,  position (0.9, 1.0, 0.55).
!          Off-diagonal (x=0.9, z=0.55, not equal), so interior to one triangle.
!      In ray-parameter units: t_entry = 0.25*sqrt(5.25), t_exit = 0.7*sqrt(5.25).
sqrt3 = sqrt(5.25_R8P)  ! local: |direction| before normalization
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
!      Origin far away, direction also pointing away from the cube.
origin    = vector_R8P(10._R8P, 10._R8P, 10._R8P)
direction = ex_R8P  ! still points along +X; the cube is at -10 in X, never reached
call cube%intersect_ray_all(ray_origin=origin, ray_direction=direction, hits=hits, status=status)
print '(A,I0,A,I0)', 'inv 4 (miss): status=', status, ' n_hits=', size(hits)
are_tests_passed(4) = (status == RAY_STATUS_OK) .and. (size(hits) == 0)

print '(A,4L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram fossil_test_ray_query
