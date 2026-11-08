program ref_booleans
use fossil, only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: a, b, c
character(10), parameter :: NAMES(4) = ['union     ', 'intersect ', 'difference', 'symdiff   ']
integer(I4P),  parameter :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
real(R8P),     parameter :: EXACT(4) = [1.875_R8P, 0.125_R8P, 0.875_R8P, 1.75_R8P]
integer(I4P)             :: status, k
call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze ! B overlaps A's corner
do k = 1, 4
   c = a                                                        ! boolean replaces self: work on a copy
   call c%boolean(other=b, op=OPS(k), status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   print '(A,A,I4,A,F6.4,A,ES8.1,A,L1,A,L1)', NAMES(k), ': facets', c%get_facets_number(), ', volume ', &
         c%get_volume(), ', error ', abs(c%get_volume() - EXACT(k)), ', watertight ', c%is_watertight(), &
         ', manifold ', c%is_manifold()
enddo
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P)) ; call b%analyze   ! 4 faces coplanar with A's
c = a
call c%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
print '(A,I3,A,F6.4,A,L1,A,L1)', 'coplanar difference: facets', c%get_facets_number(), ', volume ', c%get_volume(), &
      ', watertight ', c%is_watertight(), ', manifold ', c%is_manifold()
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(1._R8P + 1.e-13_R8P, 0._R8P, 0._R8P)) ; call b%analyze ! a gap of 1e-13
c = a
call c%boolean(other=b, op=BOOL_UNION, status=status)
print '(A,I0,A,F6.4)', 'side by side, no weld: shells ', c%get_shells_number(), ', volume ', c%get_volume()
c = a
call c%boolean(other=b, op=BOOL_UNION, status=status, weld=1.e-9_R8P)
print '(A,I0,A,F6.4,A,L1)', 'weld 1e-9:             shells ', c%get_shells_number(), ', volume ', c%get_volume(), &
      ', watertight ', c%is_watertight()
endprogram ref_booleans
