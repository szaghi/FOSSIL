!< FOSSIL, assert remove_duplicate_facets drops literal duplicate triangles.
!<
!< STL files emitted from CAD tools occasionally contain literal duplicate
!< triangles (same three vertices, possibly reversed winding) from overlapping
!< exported shells. Without removal these double-count area, volume, and signed
!< distance contributions.
!<
!< Test mesh: one healthy triangle plus three duplicates of it:
!<   - same orientation (v1, v2, v3)
!<   - reversed orientation (v1, v3, v2)
!<   - cyclic permutation (v2, v3, v1) — same orientation, different starting vertex
!< All three must be detected as duplicates of the first; survivor count must be 1.

program fossil_test_duplicate_removal

use fossil, only : surface_stl_object, facet_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(surface_stl_object)        :: surface
type(facet_object), allocatable :: facets(:)
type(vector_R8P)                :: a, b, c
integer(I4P)                    :: nf, removed
logical                         :: tests_passed(2)

! Healthy triangle vertices.
a = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
b = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
c = 0._R8P * ex_R8P + 1._R8P * ey_R8P + 0._R8P * ez_R8P

allocate(facets(4))
facets(1)%vertex = [a, b, c]   ! original
facets(2)%vertex = [a, b, c]   ! exact duplicate
facets(3)%vertex = [a, c, b]   ! reversed orientation
facets(4)%vertex = [b, c, a]   ! cyclic permutation

call surface%adopt_facets(facets)
call surface%sanitize

nf      = surface%get_facets_number()
removed = surface%get_duplicate_facets_removed()

tests_passed(1) = (nf == 1)
tests_passed(2) = (removed == 3)

print '(A,I0)', 'facets_after_sanitize:    ', nf
print '(A,I0)', 'duplicate_facets_removed: ', removed
print '(A,L1)', 'Are all tests passed? ', all(tests_passed)
if (.not. all(tests_passed)) error stop 1

endprogram fossil_test_duplicate_removal
