!< FOSSIL, assert degenerate facets are removed before they corrupt downstream queries.
!<
!< Latent bug being guarded against: `compute_normal` calls `face_normal3(..., norm='y')`,
!< which divides by the cross-product magnitude to produce a unit normal. For a
!< zero-area or near-zero-area triangle this yields NaN/Inf. The cached `self%normal`
!< field then propagates NaN into:
!<   - `compute_pseudo_normals` (edge/vertex pseudo-normals on adjacent facets)
!<   - signed-distance via SIGN_PSEUDO_NORMAL
!<   - and every downstream query touching the affected vertex.
!<
!< This test builds an in-memory mesh with one healthy triangle plus three degenerate
!< triangles (zero-area, sliver, and repeated-vertex variants). After sanitize, the
!< degenerate facets must be gone, the survivor must have a finite normal, and a
!< signed-distance query must return a finite value.
!<
!< Without the degenerate-removal pass, the NaN propagation makes the signed distance
!< at a query point near the survivor return NaN.

program fossil_test_degenerate_removal

use fossil, only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none

type(surface_stl_object)        :: surface
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: query
real(R8P)                       :: d
integer(I4P)                    :: nf, dropped
logical                         :: normal_finite
logical                         :: tests_passed(3)

! Build 4 facets:
!   1) healthy: a unit-sized triangle on the z=0 plane
!   2) degenerate: three coincident vertices (collapsed)
!   3) degenerate: collinear vertices (zero-area sliver)
!   4) degenerate: near-collinear sliver (cross-product magnitude well below tol)
allocate(facets(4))

facets(1)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(1)%vertex(2) = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(1)%vertex(3) = 0._R8P * ex_R8P + 1._R8P * ey_R8P + 0._R8P * ez_R8P

facets(2)%vertex(1) = 0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0._R8P * ez_R8P
facets(2)%vertex(2) = 0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0._R8P * ez_R8P
facets(2)%vertex(3) = 0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0._R8P * ez_R8P

facets(3)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(3)%vertex(2) = 0.5_R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(3)%vertex(3) = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P

facets(4)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(4)%vertex(2) = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
facets(4)%vertex(3) = 0.5_R8P * ex_R8P + 1.0e-14_R8P * ey_R8P + 0._R8P * ez_R8P

call surface%adopt_facets(facets)
call surface%sanitize

nf      = surface%get_facets_number()
dropped = surface%get_degenerate_facets_removed()

! 1. Three degenerate facets removed, one survivor remains.
tests_passed(1) = (nf == 1) .and. (dropped == 3)

! 2. Survivor's cached normal is finite (would be NaN without the removal pass).
normal_finite = .false.
f => surface%facet_at(1)
if (associated(f)) then
   normal_finite = ieee_is_finite(f%normal%x) .and. &
                   ieee_is_finite(f%normal%y) .and. &
                   ieee_is_finite(f%normal%z)
endif
tests_passed(2) = normal_finite

! 3. A signed-distance query above the survivor returns a finite value (would be NaN
!    if any sliver leaked through and corrupted pseudo-normals on shared vertices).
query = 0.25_R8P * ex_R8P + 0.25_R8P * ey_R8P + 0.5_R8P * ez_R8P
d = surface%distance(point=query, is_signed=.true., is_square_root=.true., &
                     sign_algorithm=SIGN_PSEUDO_NORMAL)
tests_passed(3) = ieee_is_finite(d)

print '(A,I0)',     'facets_after_sanitize:      ', nf
print '(A,I0)',     'degenerate_facets_removed:  ', dropped
print '(A,L1)',     'survivor_normal_finite:     ', normal_finite
print '(A,ES12.5)', 'signed_distance_at_query:   ', d

print '(A,L1)', 'Are all tests passed? ', all(tests_passed)
if (.not. all(tests_passed)) error stop 1

endprogram fossil_test_degenerate_removal
