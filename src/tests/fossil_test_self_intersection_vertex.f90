!< FOSSIL test: coplanar facets sharing one vertex, decided exactly (issue #69, step 2).
!<
!< Two-facet surfaces sharing the vertex s = (0,0,0) with A = s, (1,0,0), (0,1,0):
!< 1. coplanar, B's wedge inside A's: one SELFINT_COPLANAR (skipped before: unit normals parallel to 1e-12 rad).
!< 2. coplanar, wedges apart (meeting at s only): none, also with include_touching (adjacency).
!< 3. coplanar, B's edge along A's edge from s, interiors disjoint: none by default, one SELFINT_TOUCHING with
!<    include_touching.
!< 4. B through s with its free vertices at +h and -h off A's plane, crossing A along s-(0.5,0.5,0): one crossing for
!<    h = 1e-6, 1e-9, 1e-12, 1e-13, 1e-15 (below 1e-12 rad they were taken as coplanar and skipped).
!< 5. the flat-faced cube-05 (fans of coplanar facets around every vertex): none.
!< Negative control (base 04e5bbf): 1 and 4 (h = 1e-13, 1e-15) report no pair.
program fossil_test_self_intersection_vertex

use fossil,              only : surface_stl_object, intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, &
                                SELFINT_COPLANAR
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

real(R8P), parameter                   :: A(9) = [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P]
real(R8P), parameter                   :: H(5) = [1e-6_R8P, 1e-9_R8P, 1e-12_R8P, 1e-13_R8P, 1e-15_R8P]
type(intersection_pair_t), allocatable :: found(:)
type(surface_stl_object)               :: cube
logical                                :: ok(5)
integer(I4P)                           :: k

call pairs(A, [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0.2_R8P, 0._R8P, 0.2_R8P, 1._R8P, 0._R8P], .false., found)
ok(1) = size(found) == 1
if (ok(1)) ok(1) = found(1)%kind == SELFINT_COPLANAR
print '(A,L1,A,I0)', '1. coplanar overlap at a shared vertex: ', ok(1), '  pairs ', size(found)

call pairs(A, [0._R8P, 0._R8P, 0._R8P, -1._R8P, 0.5_R8P, 0._R8P, -0.5_R8P, -1._R8P, 0._R8P], .true., found)
ok(2) = size(found) == 0
print '(A,L1,A,I0)', '2. coplanar wedges apart: ', ok(2), '  pairs ', size(found)

call pairs(A, [0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P, 0._R8P, 1._R8P, -1._R8P, 0._R8P], .false., found)
ok(3) = size(found) == 0
call pairs(A, [0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P, 0._R8P, 1._R8P, -1._R8P, 0._R8P], .true., found)
ok(3) = ok(3) .and. size(found) == 1
if (ok(3)) ok(3) = found(1)%kind == SELFINT_TOUCHING
print '(A,L1)', '3. coplanar common edge segment: ', ok(3)

ok(4) = .true.
do k = 1, size(H)
   call pairs(A, [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0.3_R8P, H(k), 0.3_R8P, 1._R8P, -H(k)], .false., found)
   if (size(found) /= 1) then
      ok(4) = .false.
   elseif (found(1)%kind /= SELFINT_CROSSING) then
      ok(4) = .false.
   endif
   print '(A,ES8.1,A,I0)', '   h ', H(k), ': pairs ', size(found)
enddo
print '(A,L1)', '4. near-coplanar crossings: ', ok(4)

call cube%load_from_file(file_name='src/tests/cube-05.stl', guess_format=.true.)
call cube%find_self_intersections(pairs=found)
ok(5) = size(found) == 0
print '(A,L1,A,I0,A,I0)', '5. cube-05: ', ok(5), '  pairs ', size(found), '  facets ', cube%get_facets_number()

print '(A,L1)', 'Are all tests passed? ', all(ok)
contains
   subroutine pairs(va, vb, touching, found)
   !< Self-intersections of the two-facet surface (va, vb): 9 coordinates each.
   real(R8P),                              intent(in)  :: va(9), vb(9) !< Vertices of the two facets.
   logical,                                intent(in)  :: touching     !< include_touching.
   type(intersection_pair_t), allocatable, intent(out) :: found(:)     !< Pairs.
   type(facet_object),        allocatable              :: two(:)
   type(surface_stl_object)                            :: s2
   integer(I4P)                                        :: v

   allocate(two(2))
   do v = 1, 3
      two(1)%vertex(v) = vector_R8P(va(3 * v - 2), va(3 * v - 1), va(3 * v))
      two(2)%vertex(v) = vector_R8P(vb(3 * v - 2), vb(3 * v - 1), vb(3 * v))
   enddo
   two(1)%id = 1_I4P ; two(2)%id = 2_I4P
   call s2%adopt_facets(facets=two)
   call s2%find_self_intersections(pairs=found, include_touching=touching)
   endsubroutine pairs
endprogram fossil_test_self_intersection_vertex
