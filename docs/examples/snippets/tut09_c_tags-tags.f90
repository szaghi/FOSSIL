use fossil, only : surface_stl_object, curve_object, slice_surface, STATUS_OK
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nc = 32                     ! cells per direction of the tagged block
integer(I4P), parameter  :: kc = 16                     ! the layer compared with the 2-D twin
type(surface_stl_object) :: surface
type(curve_object)       :: section
real(R8P)                :: x(nc + 1), zc
integer(I4P)             :: tags(nc, nc, nc), tags2(nc, nc), i, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]          ! cell edges of one 32^3 block over the domain
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A)', 'tag        -5     -4     -3     -2     -1      0      1      2      3      4      5'
print '(A,11I7)', '3-D   ', (count(tags == i), i = -5, 5)
