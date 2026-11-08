program ck_tag_block
use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
real(R8P)                :: x(21)
integer(I4P)             :: tags(20, 20, 20), i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
x = [(-0.55_R8P + 0.105_R8P * i, i = 0, 20)]                 ! cell edges: 20^3 cells on [-0.55, 1.55]
call cube%tag_block(x=x, y=x, z=x, band=3, tags=tags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A)', 'tag     -4    -3    -2    -1     0     1     2     3     4'
print '(A,9I6)', 'cells', (count(tags == i), i = -4, 4)
endprogram ck_tag_block
