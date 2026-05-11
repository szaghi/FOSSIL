!< FOSSIL, strict invariants on facets after sanitize_normals.
!<
!< Whereas fossil_test_sanitize_normals only verifies the bulk-volume sign,
!< this test asserts the per-facet invariants the bulk volume can mask:
!<
!<   1. winding consistency: the cross product `(v2-v1) x (v3-v1)` agrees with
!<      the cached `self%normal` (same direction). After any winding flip,
!<      self%normal must be refreshed; mismatch here means a flip path forgot.
!<
!<   2. outward orientation on a convex body: for the unit cube, every facet's
!<      stored normal points away from the surface centroid (positive dot
!<      product with `centroid_facet - centroid_surface`). Mixed orientations
!<      can still average to a positive volume on symmetric shapes, so this
!<      check is what bulk-volume cannot do.
!<
!< Failure on cube.stl historically came from a depth-first chain in
!< sanitize_normals that did not visit every facet, plus flip_edge leaving
!< the cached normal stale. Both must be fixed for pseudo-normal sign to work.

program fossil_test_sanitize_normals_strict

use fossil, only : surface_stl_object, facet_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter :: DIRECTION_TOL = 1.0e-10_R8P

type(surface_stl_object)    :: surface
type(facet_object), pointer :: f
type(vector_R8P)            :: surface_centroid
type(vector_R8P)            :: cross_n
real(R8P)                   :: dot_outward
real(R8P)                   :: cos_winding
integer(I4P)                :: i, n_facets
integer(I4P)                :: bad_winding, bad_outward
logical                     :: tests_passed(2)

! Pipeline: load -> sanitize. NO extra analyze afterwards: sanitize must leave the
! cached self%normal in agreement with the (now-consistent) winding. Failure to do so
! is bug 2 (flip_edge/reverse_normal forget to refresh the cache); detecting it
! requires not running another compute_metrix that would mask it.
call surface%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call surface%sanitize

n_facets = surface%get_facets_number()
surface_centroid = surface%get_centroid()

bad_winding = 0
bad_outward = 0
do i = 1, n_facets
   f => surface%facet_at(i)
   if (.not. associated(f)) cycle

   ! Invariant 1: stored normal direction matches winding cross product.
   ! Use a unit-cosine check rather than equality so a (legal) magnitude
   ! difference does not register as a winding mismatch.
   cross_n = (f%vertex(2) - f%vertex(1)) .cross. (f%vertex(3) - f%vertex(1))
   cos_winding = (cross_n%x * f%normal%x + cross_n%y * f%normal%y + cross_n%z * f%normal%z)
   if (cos_winding < DIRECTION_TOL) bad_winding = bad_winding + 1

   ! Invariant 2 (convex body): outward = points away from surface centroid.
   associate(d => f%centroid - surface_centroid)
      dot_outward = d%x * f%normal%x + d%y * f%normal%y + d%z * f%normal%z
   end associate
   if (dot_outward < DIRECTION_TOL) bad_outward = bad_outward + 1
enddo

tests_passed(1) = bad_winding == 0
tests_passed(2) = bad_outward == 0

print '(A,I0,A,I0)',  'facets:                       ', n_facets, ' total'
print '(A,I0)',       'winding/normal mismatches:    ', bad_winding
print '(A,I0)',       'inward-pointing normals:      ', bad_outward
print '(A,F12.6)',    'surface volume:               ', surface%get_volume()
print '(A,L1)',       'Are all tests passed? ', all(tests_passed)
if (.not. all(tests_passed)) error stop 1

endprogram fossil_test_sanitize_normals_strict
