!< VALIDATION version (issue #45): the full, slow test, run by scripts/run_validation.sh at -O2.
!< The fast regression variant run per push is src/tests/fossil_test_alpha_wrap.f90.
!<
!< FOSSIL test: alpha wrapping, `surface%alpha_wrap` (issue #18 §1.6, issue #50). Every wrap must be status OK,
!< watertight, manifold and free of self-crossings (find_self_intersections).
!<   1. The cube, alpha 0.1: volume within 1.5% of the offset cube (Minkowski sum; chords of the coarsened wrap).
!<   2. A marching-cubes sphere (32^3): >= 90% of the vertices on the offset surface, none closer than the offset.
!<   3. The cube minus one facet, hole wider than the ball: a hollow shell lining both sides.
!<   4. The same, hole narrower than the ball: capped.
!<   5. One triangle: a closed shell around both sides.
!<   6. The issue #50 table: cube-01 closed and minus one facet, alpha 0.1, 0.3, 0.7.
!<   7. dragon.stl with the defaults; the wrap encloses it, vertices as in 2.
!<   8. cube-01 with the defaults (alpha = bbox/50): under 20000 facets (161468 before the collapse stage, #52).

program fossil_val_alpha_wrap

use fossil, only : surface_stl_object, extract_isosurface, intersection_pair_t, AWRAP_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_utils, only : PI
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere, cube_holed, wrapped_tbp, plate, cube01, cube01_holed, dragon
type(facet_object), allocatable :: plate_facets(:)
real(R8P), parameter :: TABLE_ALPHA(3) = [0.1_R8P, 0.3_R8P, 0.7_R8P]
type(facet_object), allocatable :: sphere_facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
type(facet_object), pointer :: facets(:)
real(R8P), allocatable :: values(:,:,:)
type(vector_R8P) :: bmin, bmax
real(R8P) :: x, y, z, dx, alpha
real(R8P) :: wrap_volume, expected_volume, offset
integer(I4P) :: status, i, j, k
logical :: are_tests_passed(8)

integer(I4P), parameter :: N_GRID = 32_I4P

are_tests_passed = .false.

! fixtures: the cube, the cube minus its first facet, the marching-cubes unit sphere
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
facets => cube%facets_ref()
allocate(holed_facets(size(facets, kind=I4P) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call cube_holed%adopt_facets(facets=holed_facets)
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)
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
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)

! ---- Invariant 1: capstone TBP on the cube (issue #50): status OK, closed, manifold, no self-crossing, and its
!      volume is the cube grown by the offset (Minkowski sum with a ball: 1 + 6e + 3 pi e^2 + 4/3 pi e^3).
offset = 0.033_R8P
call cube%alpha_wrap(alpha=0.1_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
are_tests_passed(1) = wrap_ok('inv 1 (cube TBP)', wrapped_tbp, status) .and. &
                       (abs(wrap_volume - expected_volume) <= 0.015_R8P * expected_volume)
print '(A,F8.5,A,F8.5)', '        volume=', wrap_volume, ' offset cube=', expected_volume

! ---- Invariant 2: capstone TBP on the marching-cubes sphere: at least 90% of the wrap vertices on the offset
!      surface, none closer to the input than the offset (issue #53).
call sphere%alpha_wrap(alpha=0.1_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
are_tests_passed(2) = wrap_ok('inv 2 (sphere TBP)', wrapped_tbp, status) .and. offset_ok(sphere, wrapped_tbp, offset)

! ---- Invariant 3: the holed cube (a hole of inradius 0.29) at alpha = 0.15: the ball gets through the hole, the
!      wrap lines both sides of every wall, two offset layers the barrier keeps apart.
call cube_holed%alpha_wrap(alpha=0.15_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
are_tests_passed(3) = wrap_ok('inv 3 (holed cube, hole open)', wrapped_tbp, status) .and. &
                       (wrap_volume >= 11._R8P * offset) .and. (wrap_volume <= 0.7_R8P)
print '(A,F8.5)', '        volume=', wrap_volume

! ---- Invariant 4: the same holed cube at alpha = 0.6: the hole is capped, volume just below the offset cube's.
offset = 0.02_R8P
call cube_holed%alpha_wrap(alpha=0.6_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
are_tests_passed(4) = wrap_ok('inv 4 (holed cube, hole capped)', wrapped_tbp, status) .and. &
                       (wrap_volume < expected_volume) .and. (wrap_volume > 0.95_R8P * expected_volume)
print '(A,F8.5,A,F8.5)', '        volume=', wrap_volume, ' offset cube=', expected_volume

! ---- Invariant 5: one triangle (zero thickness, area 1/2): a closed shell around both sides, 2 offset x area + rim.
allocate(plate_facets(1))
plate_facets(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P)
call plate_facets(1)%compute_metrix
call plate%adopt_facets(facets=plate_facets)
call plate%alpha_wrap(alpha=0.2_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
are_tests_passed(5) = wrap_ok('inv 5 (one triangle)', wrapped_tbp, status) .and. &
                       (wrap_volume >= 0.019_R8P) .and. (wrap_volume <= 0.0225_R8P)
print '(A,F8.5)', '        volume=', wrap_volume

! ---- Invariant 6: the measured table of issue #50: cube-01.stl closed and minus one facet, alpha 0.1, 0.3, 0.7,
!      default offset (alpha / 30): every wrap OK, closed, manifold, not self-crossing.
call cube01%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
facets => cube01%facets_ref()
if (allocated(holed_facets)) deallocate(holed_facets) ! adopt_facets may have taken it
allocate(holed_facets(size(facets, kind=I4P) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call cube01_holed%adopt_facets(facets=holed_facets)
are_tests_passed(6) = .true.
do i = 1_I4P, 3_I4P
   alpha = TABLE_ALPHA(i)
   call cube01%alpha_wrap(alpha=alpha, wrapped=wrapped_tbp, status=status)
   are_tests_passed(6) = wrap_ok('inv 6 (cube-01, closed)', wrapped_tbp, status) .and. are_tests_passed(6)
   print '(A,F5.2,A,F8.5)', '        alpha=', alpha, ' volume=', wrapped_tbp%get_volume()
   call cube01_holed%alpha_wrap(alpha=alpha, wrapped=wrapped_tbp, status=status)
   are_tests_passed(6) = wrap_ok('inv 6 (cube-01, holed)', wrapped_tbp, status) .and. are_tests_passed(6)
   print '(A,F5.2,A,F8.5)', '        alpha=', alpha, ' volume=', wrapped_tbp%get_volume()
enddo

! ---- Invariant 7: dragon.stl with the defaults (it was neither watertight nor manifold before issue #50): OK, closed,
!      manifold, not self-crossing, and it encloses the dragon (larger volume).
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%alpha_wrap(wrapped=wrapped_tbp, status=status)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
offset = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2) / 50._R8P / 30._R8P
are_tests_passed(7) = wrap_ok('inv 7 (dragon, defaults)', wrapped_tbp, status) .and. &
                       offset_ok(dragon, wrapped_tbp, offset) .and. (wrapped_tbp%get_volume() > dragon%get_volume())
print '(A,F10.3,A,F10.3)', '        volume=', wrapped_tbp%get_volume(), ' dragon=', dragon%get_volume()

! ---- Invariant 8: cube-01 with the defaults, coarsened (issue #52): under 20000 facets, and a valid wrap.
call cube01%alpha_wrap(wrapped=wrapped_tbp, status=status)
are_tests_passed(8) = wrap_ok('inv 8 (cube-01, defaults)', wrapped_tbp, status) .and. &
                      (wrapped_tbp%get_facets_number() < 20000_I4P)

print '(A,8L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. (all(are_tests_passed))) error stop 1

contains

   function wrap_ok(label, w, st) result(ok)
   !< Status OK, watertight, manifold, no self-crossing (find_self_intersections); prints the line.
   character(*),             intent(in) :: label
   type(surface_stl_object), intent(in) :: w
   integer(I4P),             intent(in) :: st
   logical                              :: ok
   type(intersection_pair_t), allocatable :: pairs(:)

   call w%find_self_intersections(pairs=pairs)
   print '(A,A,I0,A,I0,2(A,L1),A,I0)', label, ': status=', st, ' n_facets=', w%get_facets_number(), &
         ' watertight=', w%is_watertight(), ' manifold=', w%is_manifold(), ' self-crossings=', size(pairs)
   ok = (st == AWRAP_STATUS_OK) .and. w%is_watertight() .and. w%is_manifold() .and. (size(pairs) == 0)
   endfunction wrap_ok

   function offset_ok(input, w, eps) result(ok)
   !< At least 90% of the vertices of w at eps from input (relative 1e-6), none closer than eps (relative 1e-9; the
   !< repair pass projects a vertex it moves back onto the offset surface, issue #53).
   type(surface_stl_object), intent(in) :: input, w
   real(R8P),                intent(in) :: eps
   logical                              :: ok
   type(facet_object), pointer          :: f
   real(R8P)                            :: d, dmin
   integer(I4P)                         :: i, v, n_on, n_all

   dmin = huge(1._R8P) ; n_on = 0 ; n_all = 0
   do i = 1, w%get_facets_number()
      f => w%facet_at(i)
      do v = 1, 3
         d = input%distance(point=f%vertex(v), is_signed=.false., is_square_root=.true.)
         dmin = min(dmin, d) ; n_all = n_all + 1
         if (abs(d - eps) <= 1.e-6_R8P * eps) n_on = n_on + 1
      enddo
   enddo
   print '(A,F12.9,A,F6.3)', '        min vertex distance / offset=', dmin / eps, ' on the offset surface: ', &
         real(n_on, R8P) / real(n_all, R8P)
   ok = (dmin >= (1._R8P - 1.e-9_R8P) * eps) .and. (real(n_on, R8P) >= 0.9_R8P * real(n_all, R8P))
   endfunction offset_ok
endprogram fossil_val_alpha_wrap
