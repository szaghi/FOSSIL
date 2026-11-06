!< FOSSIL test: the opt-in weld of the exact boolean (issue #65, step 5).
!<
!< 1. `weld_vertices`: a point of B near one of A takes it; two points of B near one of A are both left alone (the weld
!<    never merges B's vertices); a point equidistant from two of A is left alone; nothing beyond the distance.
!< 2. A sphere (10^3 samples) and its copy translated by 1e-12 along (1, 2, 3):
!<    - with `weld=1e-9` B becomes A bit for bit: the union and the intersection are A (same facets, same volume),
!<      the difference is empty;
!<    - without it the result is the exact one: status OK, watertight, manifold, volume within 1e-11 of A's (the
!<      slivers between the two copies).
program fossil_test_boolean_weld

use fossil,                only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_boolean,        only : weld_vertices
use fossil_boolean_bodies, only : BODY_SPHERE, mc_body, move_body
use penf,                  only : I4P, R8P

implicit none

integer(I4P), parameter :: OPS3(3) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE]

real(R8P)                 :: xa(3, 2), xb(3, 4), va, d(3)
integer(I4P), allocatable :: map(:)
type(surface_stl_object)  :: a, b, r
integer(I4P)              :: st, k
logical                   :: ok1, ok2, ok3

! ---- 1. weld_vertices
xa(:, 1) = [0._R8P, 0._R8P, 0._R8P] ; xa(:, 2) = [1._R8P, 0._R8P, 0._R8P]
xb(:, 1) = [1._R8P + 1.e-12_R8P, 0._R8P, 0._R8P]   ! -> 2
xb(:, 2) = [5.e-13_R8P, 0._R8P, 0._R8P]            ! two of B near A(1): left alone
xb(:, 3) = [0._R8P, 7.e-13_R8P, 0._R8P]
xb(:, 4) = [0._R8P, 5._R8P, 0._R8P]                ! far
call weld_vertices(xa, xb, 1.e-9_R8P, map)
ok1 = all(map == [2, 0, 0, 0])
xb(:, 1) = [0.5_R8P, 0._R8P, 0._R8P]               ! equidistant from A(1) and A(2)
call weld_vertices(xa, xb(:, 1:1), 0.6_R8P, map)
ok1 = ok1 .and. all(map == [0])
print '(A,L1)', '1. weld_vertices: ', ok1

! ---- 2. sphere and its copy moved by 1e-12
call mc_body(BODY_SPHERE, 10_I4P, a)
d = 1.e-12_R8P * [1._R8P, 2._R8P, 3._R8P] / sqrt(14._R8P)
b = a ; call move_body(b, [0._R8P, 0._R8P, 1._R8P], 0._R8P, d)
va = a%get_volume()
ok2 = .true.
do k = 1, 3
   r = a
   call r%boolean(other=b, op=OPS3(k), status=st, weld=1.e-9_R8P)
   if (k < 3) then
      ok2 = ok2 .and. st == BOOL_STATUS_OK .and. r%get_facets_number() == a%get_facets_number() .and. &
            abs(r%get_volume() - va) <= 1.e-14_R8P * va .and. r%is_watertight() .and. r%is_manifold()
   else
      ok2 = ok2 .and. st == BOOL_STATUS_OK .and. r%get_facets_number() == 0
   endif
   print '(A,I0,A,I0,A,I0,A,ES10.2)', '   weld op ', k, ': status ', st, ', facets ', r%get_facets_number(), &
         ', dV/V ', (r%get_volume() - merge(va, 0._R8P, k < 3)) / va
enddo
print '(A,L1)', '2. welded copy gives the clean body: ', ok2
r = a
call r%boolean(other=b, op=BOOL_UNION, status=st)
ok3 = st == BOOL_STATUS_OK .and. r%is_watertight() .and. r%is_manifold() .and. abs(r%get_volume() - va) <= 1.e-11_R8P * va
print '(A,I0,A,I0,A,ES10.2,A,L1)', '3. without the weld: status ', st, ', facets ', r%get_facets_number(), &
      ', dV/V ', (r%get_volume() - va) / va, '  ok ', ok3

print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3
endprogram fossil_test_boolean_weld
