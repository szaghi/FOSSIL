!run ck_rotate ck_rotate
program ck_rotate
use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : vector_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: wing
type(vector_R8P)         :: bmin, bmax

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
bmin = wing%get_bmin() ; bmax = wing%get_bmax()
print '(A,2F8.4,A,2F8.4)', 'x range ', bmin%x, bmax%x, ', y range ', bmin%y, bmax%y
! 10 degrees nose up: a clockwise turn about the span axis z, through the quarter chord (the origin)
call wing%rotate(axis=ez_R8P, angle=-10._R8P * PI / 180._R8P, center=vector_R8P(0._R8P, 0._R8P, 0._R8P))
call wing%analyze                                             ! box, normals and tree of the new position
bmin = wing%get_bmin() ; bmax = wing%get_bmax()
print '(A,2F8.4,A,2F8.4)', 'x range ', bmin%x, bmax%x, ', y range ', bmin%y, bmax%y
print '(A,F9.6)', 'volume  ', wing%get_volume()
endprogram ck_rotate
