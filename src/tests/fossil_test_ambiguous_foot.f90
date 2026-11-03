!< FOSSIL test: multiple foot points (issue #37 step 7).
program fossil_test_ambiguous_foot
!< FOSSIL test: multiple foot points (issue #37 step 7).
!<
!< 3-D, unit cube (cube-02): the centre (equidistant from six faces) and (0.3, 0.3, 0.5) (on the
!< diagonal plane between the faces x = 0 and y = 0) are flagged, the second foot point on another
!< face; a generic point (0.2, 0.35, 0.6) is not, and its second foot point is its foot point.
!< Two cubes 0.2 apart: the mid-gap point (1.1, 0.55, 0.4) is flagged, with feet on both cubes.
!< 2-D, L-shape with a reflex vertex at (1, 1): the medial points (1.5, 1.5) and (1.5, 0.5) are
!< flagged; (0.3, 0.2) is not.

use fossil, only : surface_stl_object, facet_object, curve_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)        :: cube, two
type(curve_object)              :: c
type(facet_object), allocatable :: f(:), g(:)
type(facet_object), pointer     :: fp
real(R8P)                       :: d(3), cp(3, 3), sf(3, 3), d2(3), cp2(2, 3), sf2(2, 3)
logical                         :: amb(3), ok(3)
integer(I4P)                    :: i, k, n

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call cube%query_many(points=[vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(0.3_R8P, 0.3_R8P, 0.5_R8P), &
                             vector_R8P(0.2_R8P, 0.35_R8P, 0.6_R8P)], distances=d, closest_points=cp, ambiguous=amb, &
                     second_foot_points=sf)
ok(1) = amb(1) .and. amb(2) .and. .not. amb(3) .and. norm2(sf(:, 2) - cp(:, 2)) > 0.1_R8P .and. all(sf(:, 3) == cp(:, 3)) .and. &
        abs(norm2(sf(:, 2) - [0.3_R8P, 0.3_R8P, 0.5_R8P]) - 0.3_R8P) < 1.e-12_R8P
print '(A,L1,A,3L2)', 'cube: centre, diagonal plane flagged, generic not: ', ok(1), '  flags:', amb

n = cube%get_facets_number()
allocate(f(n), g(n))
do i = 1, n
   fp => cube%facet_at(i)
   f(i) = fp ; f(i)%vertex_id = 0_I4P ; g(i) = f(i)
   do k = 1, 3
      g(i)%vertex(k)%x = g(i)%vertex(k)%x + 1.2_R8P
   enddo
   call g(i)%compute_metrix
enddo
f = [f, g]
call two%adopt_facets(facets=f)
call two%query_many(points=[vector_R8P(1.1_R8P, 0.55_R8P, 0.4_R8P)], distances=d(1:1), closest_points=cp(:, 1:1), &
                    ambiguous=amb(1:1), second_foot_points=sf(:, 1:1))
ok(2) = amb(1) .and. minval([abs(cp(1, 1) - 1._R8P), abs(cp(1, 1) - 1.2_R8P)]) < 1.e-12_R8P .and. &
        abs(abs(sf(1, 1) - cp(1, 1)) - 0.2_R8P) < 1.e-12_R8P
print '(A,L1)', 'two cubes, mid-gap point flagged, feet on both: ', ok(2)

call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P], &
                   loop_sizes=[6_I4P])
call c%query_many(points=reshape([1.5_R8P, 1.5_R8P, 1.5_R8P, 0.5_R8P, 0.3_R8P, 0.2_R8P], [2, 3]), distances=d2, &
                  closest_points=cp2, ambiguous=amb, second_foot_points=sf2)
ok(3) = amb(1) .and. amb(2) .and. .not. amb(3) .and. all(sf2(:, 3) == cp2(:, 3))
print '(A,L1,A,3L2)', '2-D L-shape: medial points flagged, generic not: ', ok(3), '  flags:', amb

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_ambiguous_foot
