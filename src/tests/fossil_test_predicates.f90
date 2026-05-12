!< FOSSIL, assert the composite predicates is_watertight / is_manifold / is_volume
!< produce sensible answers across known-good and known-bad inputs.
!<
!< Cases:
!<   1. Clean cube after sanitize:                watertight = T, manifold = T, volume = T
!<   2. Dragon (62 non-manifold edges):           watertight = F, manifold = F, volume = F
!<   3. Open shell (single triangle):             watertight = F, manifold = T, volume = F
!<
!< The dragon and open-shell results depend on the corresponding counters being
!< populated correctly by `build_connectivity` and `sanitize_normals`.

program fossil_test_predicates

use fossil, only : surface_stl_object, facet_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P

implicit none

type(surface_stl_object)        :: surface
type(facet_object), allocatable :: facets(:)
logical                         :: all_passed

all_passed = .true.

call check_cube(all_passed)
call check_dragon(all_passed)
call check_open_shell(all_passed)

print '(A,L1)', 'Are all tests passed? ', all_passed
if (.not. all_passed) error stop 1

contains

   subroutine check_cube(passed)
   logical, intent(inout) :: passed
   logical                :: w, m, v
   logical                :: this_ok

   call surface%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call surface%sanitize
   w = surface%is_watertight()
   m = surface%is_manifold()
   v = surface%is_volume()
   this_ok = w .and. m .and. v
   passed = passed .and. this_ok
   print '(A,3L2,A,L1)', 'cube       : (W,M,V)=', w, m, v, ' passed=', this_ok
   endsubroutine check_cube

   subroutine check_dragon(passed)
   logical, intent(inout) :: passed
   logical                :: w, m, v
   logical                :: this_ok

   call surface%destroy
   call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
   call surface%sanitize
   w = surface%is_watertight()
   m = surface%is_manifold()
   v = surface%is_volume()
   ! Dragon has 62 non-manifold edges and disconnected edges; none of the three
   ! predicates should hold.
   this_ok = (.not. w) .and. (.not. m) .and. (.not. v)
   passed = passed .and. this_ok
   print '(A,3L2,A,L1)', 'dragon     : (W,M,V)=', w, m, v, ' passed=', this_ok
   endsubroutine check_dragon

   subroutine check_open_shell(passed)
   logical, intent(inout) :: passed
   logical                :: w, m, v
   logical                :: this_ok

   call surface%destroy
   allocate(facets(1))
   facets(1)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
   facets(1)%vertex(2) = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
   facets(1)%vertex(3) = 0._R8P * ex_R8P + 1._R8P * ey_R8P + 0._R8P * ez_R8P
   call surface%adopt_facets(facets)
   call surface%sanitize
   w = surface%is_watertight()
   m = surface%is_manifold()
   v = surface%is_volume()
   ! A lone triangle: 3 boundary edges, but no non-manifold edges. So manifold = T,
   ! watertight = F (has boundaries), is_volume = F (not watertight).
   this_ok = (.not. w) .and. m .and. (.not. v)
   passed = passed .and. this_ok
   print '(A,3L2,A,L1)', 'open shell : (W,M,V)=', w, m, v, ' passed=', this_ok
   endsubroutine check_open_shell

endprogram fossil_test_predicates
