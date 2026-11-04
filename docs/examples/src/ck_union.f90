!run ck_union ck_union
program ck_union
use fossil, only : surface_stl_object, BOOL_UNION, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: a, b
integer(I4P)             :: status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze
call a%boolean(other=b, op=BOOL_UNION, status=status) ! a is replaced by a + b
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
print '(A,F8.5,A,L1)', 'union volume ', a%get_volume(), ', watertight ', a%is_watertight()
endprogram ck_union
