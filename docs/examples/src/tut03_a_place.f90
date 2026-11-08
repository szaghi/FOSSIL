!run tut03 tut03_a_place
!< Tutorial 3, first run: place the repaired dragon in the domain (translate, rotate, resize) and save it.
program tut03_a_place
!region place
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: scale

call surface%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
print '(A,3F9.4)', 'CAD box min ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'CAD box max ', bmax%x, bmax%y, bmax%z
call surface%translate(delta=-0.5_R8P * (bmin + bmax))             ! box centre to the origin
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)                 ! a quarter turn about z
call surface%analyze                                                ! rebuild box, volume, tree
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
print '(A,F9.6)', 'scale       ', scale
call surface%resize(x=scale, y=scale, z=scale)                     ! longest side 1, about the origin
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
print '(A,3F9.4)', 'box min     ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'box max     ', bmax%x, bmax%y, bmax%z
print '(A,F9.6)',  'volume      ', surface%get_volume()
call surface%save_into_file(file_name='dragon-placed.stl')
!endregion place
endprogram tut03_a_place
