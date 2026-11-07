!< FOSSIL test: exact self-intersection narrow phase for facets with no vertex in common (issue #69, step 1).
!<
!< Two-facet surfaces against the triangle A = (0,0,0), (1,0,0), (0,1,0):
!< 1. a crossing 2e-9 long (B pierces A, one vertex on A's interior): one SELFINT_CROSSING with that segment.
!< 2. a coplanar overlap with positive area: one SELFINT_COPLANAR. Möller never reported coplanar pairs.
!< 3. coplanar, touching along a boundary segment (B across A's hypotenuse): none by default, one SELFINT_TOUCHING with
!<    include_touching.
!< 4. B touching A's interior at one vertex: none by default, one SELFINT_TOUCHING with include_touching.
!< 5. B parallel to A, 1e-13 above it and overlapping it in projection: none (separated, decided exactly).
!< Negative control (base 04e5bbf, the floating Möller test): 2 reports no pair, 3 and 4 have no include_touching.
program fossil_test_self_intersection_exact

use fossil,              only : surface_stl_object, intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, &
                                SELFINT_COPLANAR
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

real(R8P), parameter                   :: A(9) = [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P]
type(intersection_pair_t), allocatable :: found(:)
logical                                :: ok(5)

! 1. crossing 2e-9 long
call pairs(A, [0.3_R8P, 0.3_R8P, -1._R8P, 0.3_R8P, 0.3_R8P, 1._R8P, 0.3_R8P + 2e-9_R8P, 0.3_R8P, 0._R8P], .false., found)
ok(1) = size(found) == 1
if (ok(1)) ok(1) = found(1)%kind == SELFINT_CROSSING .and. &
                   abs(abs(found(1)%p%x - found(1)%q%x) - 2e-9_R8P) < 1e-15_R8P .and. &
                   abs(found(1)%p%y - 0.3_R8P) < 1e-15_R8P .and. found(1)%p%z == 0._R8P .and. found(1)%q%z == 0._R8P
print '(A,L1,A,I0)', '1. crossing 2e-9 long: ', ok(1), '  pairs ', size(found)

! 2. coplanar overlap
call pairs(A, [0.2_R8P, 0.2_R8P, 0._R8P, 1.2_R8P, 0.2_R8P, 0._R8P, 0.2_R8P, 1.2_R8P, 0._R8P], .false., found)
ok(2) = size(found) == 1
if (ok(2)) ok(2) = found(1)%kind == SELFINT_COPLANAR
print '(A,L1,A,I0)', '2. coplanar overlap: ', ok(2), '  pairs ', size(found)

! 3. coplanar, touching along a boundary segment
call pairs(A, [0.25_R8P, 0.75_R8P, 0._R8P, 0.75_R8P, 0.25_R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P], .false., found)
ok(3) = size(found) == 0
call pairs(A, [0.25_R8P, 0.75_R8P, 0._R8P, 0.75_R8P, 0.25_R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P], .true., found)
ok(3) = ok(3) .and. size(found) == 1
if (ok(3)) ok(3) = found(1)%kind == SELFINT_TOUCHING
print '(A,L1)', '3. coplanar touching: ', ok(3)

! 4. touching at one vertex
call pairs(A, [0.3_R8P, 0.3_R8P, 0._R8P, 1._R8P, 0.3_R8P, 1._R8P, 0.3_R8P, 1._R8P, 1._R8P], .false., found)
ok(4) = size(found) == 0
call pairs(A, [0.3_R8P, 0.3_R8P, 0._R8P, 1._R8P, 0.3_R8P, 1._R8P, 0.3_R8P, 1._R8P, 1._R8P], .true., found)
ok(4) = ok(4) .and. size(found) == 1
if (ok(4)) ok(4) = found(1)%kind == SELFINT_TOUCHING .and. found(1)%p%x == 0.3_R8P .and. found(1)%p%z == 0._R8P
print '(A,L1)', '4. touching at a vertex: ', ok(4)

! 5. parallel, 1e-13 above
call pairs(A, [0.2_R8P, 0.2_R8P, 1e-13_R8P, 1.2_R8P, 0.2_R8P, 1e-13_R8P, 0.2_R8P, 1.2_R8P, 2e-13_R8P], .true., found)
ok(5) = size(found) == 0
print '(A,L1,A,I0)', '5. parallel 1e-13 above: ', ok(5), '  pairs ', size(found)

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
endprogram fossil_test_self_intersection_exact
