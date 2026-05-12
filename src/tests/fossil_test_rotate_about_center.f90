!< FOSSIL test: rotate around an arbitrary centre of rotation (issue #6).
!<
!< Verifies three properties of the new optional `center` argument on
!< surface_stl_object%rotate_by_*:
!<
!<   1. Centroid invariance: rotating a non-origin-centred mesh by an arbitrary
!<      angle about its own centroid must leave the centroid fixed.
!<   2. Round-trip: a +theta followed by a -theta rotation about the same centre
!<      must return every facet vertex to its starting position (mod EPS).
!<   3. Backward compatibility: omitting `center` is equivalent to passing
!<      `center=vector_R8P(0,0,0)` (rotation about the world origin).

program fossil_test_rotate_about_center

use fossil,   only : surface_stl_object
use penf,     only : I4P, R8P
use vecfor,   only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface, surface_default, surface_explicit
type(vector_R8P)         :: offset, centroid_before, centroid_after, origin
type(vector_R8P)         :: bmin_before, bmax_before, bmin_after, bmax_after
real(R8P), parameter     :: HALF_PI = 1.5707963267948966_R8P
real(R8P), parameter     :: TOL = 1.0e-9_R8P
logical                  :: are_tests_passed(3)
real(R8P)                :: dx, dy, dz

are_tests_passed = .false.
origin = vector_R8P(0._R8P, 0._R8P, 0._R8P)

! --- 1. Centroid invariance under rotation about its own centroid. ---

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
! Shove the surface away from the origin so the centroid is non-zero -- this is
! the regime where the old API (rotation about origin) silently translated the
! mesh, which is the surprise the issue #6 reporter hit.
offset = vector_R8P(100._R8P, -50._R8P, 25._R8P)
call surface%translate(delta=offset)
centroid_before = surface%get_centroid()

call surface%rotate(axis=ey_R8P, angle=HALF_PI, center=centroid_before)
centroid_after = surface%get_centroid()

dx = abs(centroid_after%x - centroid_before%x)
dy = abs(centroid_after%y - centroid_before%y)
dz = abs(centroid_after%z - centroid_before%z)
are_tests_passed(1) = (dx <= TOL) .and. (dy <= TOL) .and. (dz <= TOL)

! --- 2. Round trip: +theta then -theta about the same centre restores bbox. ---

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call surface%translate(delta=offset)
bmin_before = surface%get_bmin()
bmax_before = surface%get_bmax()

call surface%rotate(axis=ez_R8P, angle= HALF_PI, center=surface%get_centroid())
call surface%rotate(axis=ez_R8P, angle=-HALF_PI, center=surface%get_centroid())

bmin_after = surface%get_bmin()
bmax_after = surface%get_bmax()
are_tests_passed(2) = (abs(bmin_after%x - bmin_before%x) <= TOL) .and. &
                      (abs(bmin_after%y - bmin_before%y) <= TOL) .and. &
                      (abs(bmin_after%z - bmin_before%z) <= TOL) .and. &
                      (abs(bmax_after%x - bmax_before%x) <= TOL) .and. &
                      (abs(bmax_after%y - bmax_before%y) <= TOL) .and. &
                      (abs(bmax_after%z - bmax_before%z) <= TOL)

! --- 3. Backward compat: no `center` arg == `center=origin`. ---

call surface_default %load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call surface_explicit%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call surface_default %rotate(axis=ex_R8P, angle=HALF_PI)
call surface_explicit%rotate(axis=ex_R8P, angle=HALF_PI, center=origin)

centroid_before = surface_default %get_centroid()
centroid_after  = surface_explicit%get_centroid()
are_tests_passed(3) = (abs(centroid_after%x - centroid_before%x) <= TOL) .and. &
                      (abs(centroid_after%y - centroid_before%y) <= TOL) .and. &
                      (abs(centroid_after%z - centroid_before%z) <= TOL)

print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_rotate_about_center
