!< FOSSIL, test distance queries on an empty surface (regression test for audit S5).
!<
!< An empty `surface_stl_object` (`facets_number == 0`) used to leave the
!< `intent(out) :: distance` argument undefined — undefined behaviour per the Fortran
!< standard. After the fix in audit #14 (S5) the routine initializes `distance` to
!< `MaxR8P` so callers always get a defined value.

program fossil_test_empty_surface

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface             !< Empty surface (never analized, never loaded).
real(R8P)                :: distance            !< Distance output.
integer(I4P)             :: facet_index         !< Facet-index output.
logical                  :: are_tests_passed(4) !< Result of tests check.

are_tests_passed = .false.

! T1: distance() function on empty surface returns MaxR8P (not UB).
distance = surface%distance(point = 0._R8P * ex_R8P)
are_tests_passed(1) = distance == MaxR8P

! T2: compute_distance() subroutine on empty surface defines distance.
call surface%compute_distance(point = 0._R8P * ex_R8P, distance = distance)
are_tests_passed(2) = distance == MaxR8P

! T3: compute_distance with optional facet_index — index must be the sentinel 0.
call surface%compute_distance(point = 0._R8P * ex_R8P, distance = distance, facet_index = facet_index)
are_tests_passed(3) = (distance == MaxR8P) .and. (facet_index == 0_I4P)

! T4: get_facets_number on default-initialized surface is 0.
are_tests_passed(4) = surface%get_facets_number() == 0_I4P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1
endprogram fossil_test_empty_surface
