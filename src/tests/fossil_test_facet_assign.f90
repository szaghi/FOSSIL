!< FOSSIL, test facet assignment preserves pseudo-normals (regression test for audit D6).
!<
!< Before the fix, `facet_assign_facet` was silently dropping the four pseudo-normal
!< fields (edge_pnormal(3) and vertex_pnormal(3)) on `=`. After the D7 collapse these
!< are now indexed arrays — this test populates them on the rhs, assigns, and checks
!< the lhs has the same values.

program fossil_test_facet_assign

use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(facet_object) :: rhs                 !< Source facet.
type(facet_object) :: lhs                 !< Destination facet.
integer(I4P)       :: e                   !< Edge index.
logical            :: are_tests_passed(8) !< Result of tests check.

are_tests_passed = .false.

! Populate rhs pseudo-normals with exactly-representable values (powers of 2 scaled
! integers) so the verification can use `==` without floating-point traps.
do e=1, 3
   rhs%edge_pnormal(e)   = real(e,        R8P) * ex_R8P + real(2 * e,  R8P) * ey_R8P
   rhs%vertex_pnormal(e) = real(4 * e,    R8P) * ey_R8P + real(8 * e,  R8P) * ez_R8P
enddo
rhs%fcon_edge = [7_I4P, 11_I4P, 13_I4P]

! Assignment.
lhs = rhs

! Verify: each component copied.
are_tests_passed(1) = all(lhs%fcon_edge == [7_I4P, 11_I4P, 13_I4P])

are_tests_passed(2) = (lhs%edge_pnormal(1)%x == 1._R8P) .and. (lhs%edge_pnormal(1)%y == 2._R8P)
are_tests_passed(3) = (lhs%edge_pnormal(2)%x == 2._R8P) .and. (lhs%edge_pnormal(2)%y == 4._R8P)
are_tests_passed(4) = (lhs%edge_pnormal(3)%x == 3._R8P) .and. (lhs%edge_pnormal(3)%y == 6._R8P)

are_tests_passed(5) = (lhs%vertex_pnormal(1)%y ==  4._R8P) .and. (lhs%vertex_pnormal(1)%z ==  8._R8P)
are_tests_passed(6) = (lhs%vertex_pnormal(2)%y ==  8._R8P) .and. (lhs%vertex_pnormal(2)%z == 16._R8P)
are_tests_passed(7) = (lhs%vertex_pnormal(3)%y == 12._R8P) .and. (lhs%vertex_pnormal(3)%z == 24._R8P)

! T8: assignment is independent — mutating rhs does not affect lhs (deep copy, not aliasing).
rhs%fcon_edge = [0_I4P, 0_I4P, 0_I4P]
are_tests_passed(8) = all(lhs%fcon_edge == [7_I4P, 11_I4P, 13_I4P])

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1
endprogram fossil_test_facet_assign
