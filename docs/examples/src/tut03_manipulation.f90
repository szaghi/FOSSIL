!run tut03 tut03_manipulation
!< Tutorial 3: place the dragon in the domain (translate, rotate, resize), clip it, mirror and merge.
program tut03_manipulation
!region place
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, half, rest, twin
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: scale

call surface%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))             ! box centre to the origin
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)                 ! a quarter turn about z
call surface%analyze                                                ! rebuild box, volume, tree
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)                     ! longest side 1, about the origin
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
print '(A,3F9.4)', 'box min   ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'box max   ', bmax%x, bmax%y, bmax%z
print '(A,F9.6)',  'volume    ', surface%get_volume()
call surface%save_into_file(file_name='dragon-placed.stl')
!endregion place
!region clip
half = surface
call half%clip(bmin=bmin, bmax=bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P, remainder=rest)
print '(A,I0,A,I0)', 'clip y<=0 ', half%get_facets_number(), ' facets, remainder ', rest%get_facets_number()
!endregion clip
!region merge
twin = half
call twin%mirror(normal=ey_R8P)                                     ! reflect through the plane y = 0
call twin%analyze
call half%merge_solids(other=twin)                                  ! merge_solids analyzes
print '(A,I0,A)', 'merged    ', half%get_facets_number(), ' facets (half + mirrored half)'
!endregion merge
endprogram tut03_manipulation
