!run tut03-clip tut03_b_clip
!< Tutorial 3, second run: split the placed dragon by a box (clip with a remainder), mirror the half and merge it.
program tut03_b_clip
!region clip
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface, half, rest, twin
type(vector_R8P)         :: bmin, bmax

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
half = surface
call half%clip(bmin=bmin, bmax=bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P, remainder=rest)  ! the box y <= 0
print '(A,I0,A,I0,A,I0)', 'clip y<=0 ', half%get_facets_number(), ' facets, remainder ', rest%get_facets_number(), &
                          ', sum ', half%get_facets_number() + rest%get_facets_number()
!endregion clip
!region merge
twin = half
call twin%mirror(normal=ey_R8P)                                     ! reflect through the plane y = 0
call twin%analyze
call half%merge_solids(other=twin)                                  ! merge_solids analyzes
print '(A,I0,A)', 'merged    ', half%get_facets_number(), ' facets (half + mirrored half)'
bmin = half%get_bmin() ; bmax = half%get_bmax()
print '(A,3F9.4)', 'box min   ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'box max   ', bmax%x, bmax%y, bmax%z
!endregion merge
endprogram tut03_b_clip
