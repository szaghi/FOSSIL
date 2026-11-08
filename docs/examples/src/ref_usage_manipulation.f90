!run ref-usage-manipulation ref_usage_manipulation
!< Usage reference, manipulation: translate, rotate, mirror, resize, clip and merge the dragon; each operation starts
!< from the dragon as loaded and prints the box and the volume that analyze computes afterwards.
program ref_usage_manipulation
!region declarations
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P, rotation_matrix_R8P, mirror_matrix_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: dragon, surface, remainder, surface_1, surface_2
type(vector_R8P)         :: delta, axis, normal, factor, bmin, bmax
real(R8P)                :: angle, matrix(3,3)

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
!endregion declarations
call report('as loaded', dragon)
!region translate
surface = dragon
delta%x = 1._R8P ; delta%y = 0._R8P ; delta%z = 0._R8P
call surface%translate(delta=delta)                       ! by a 3-D vector
call surface%translate(x=2._R8P)                          ! by scalar components, any combination
call surface%translate(y=1.5_R8P)
call surface%translate(x=1._R8P, y=2._R8P, z=0.5_R8P)
call surface%analyze                                      ! box, volume, tree of the moved surface
!endregion translate
call report('translate', surface)
!region rotate
surface = dragon
axis  = ex_R8P                                            ! about the x axis
angle = pi / 2._R8P                                       ! a quarter turn, in radians
call surface%rotate(axis=axis, angle=angle)               ! about the origin
call surface%analyze
!endregion rotate
call report('rotate', surface)
!region rotate_matrix
surface = dragon
matrix = rotation_matrix_R8P(axis=ex_R8P, angle=pi / 2._R8P)
call surface%rotate(matrix=matrix)                        ! the same quarter turn
call surface%analyze
!endregion rotate_matrix
call report('rotate matrix', surface)
!region rotate_center
surface = dragon
call surface%rotate(axis=ex_R8P, angle=pi / 2._R8P, center=surface%get_centroid())  ! about the centroid
call surface%analyze
!endregion rotate_center
call report('rotate centre', surface)
!region mirror
surface = dragon
normal = ex_R8P                                           ! through the yz plane (normal = x axis)
call surface%mirror(normal=normal)
call surface%analyze
!endregion mirror
call report('mirror x', surface)
!region mirror_oblique
surface = dragon
normal = ex_R8P + ey_R8P                                  ! through an oblique plane (normalised inside)
call surface%mirror(normal=normal)
call surface%analyze
!endregion mirror_oblique
call report('mirror x+y', surface)
!region mirror_matrix
surface = dragon
matrix = mirror_matrix_R8P(normal=ex_R8P)
call surface%mirror(matrix=matrix)                        ! the same as normal=ex_R8P
call surface%analyze
!endregion mirror_matrix
call report('mirror matrix', surface)
!region resize
surface = dragon
factor = 2._R8P * ex_R8P + 2._R8P * ey_R8P + 2._R8P * ez_R8P
call surface%resize(factor=factor)                        ! uniform, about the origin
call surface%analyze
!endregion resize
call report('resize x2', surface)
!region resize_axes
surface = dragon
call surface%resize(x=0.5_R8P, z=1.2_R8P)                 ! per axis; y keeps its size
call surface%analyze
!endregion resize_axes
call report('resize x,z', surface)
!region resize_centroid
surface = dragon
call surface%resize(factor=factor, respect_centroid=.true.)  ! about the centroid of the last analyze
call surface%analyze
!endregion resize_centroid
call report('resize centre', surface)
!region clip
surface = dragon
bmin%x = -15._R8P ; bmin%y = -5._R8P ; bmin%z =  0._R8P
bmax%x =   0._R8P ; bmax%y =  5._R8P ; bmax%z = 20._R8P
call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)  ! clip re-runs analyze on both
print '(A,I0,A,I0,A,I0)', 'clip           kept ', surface%get_facets_number(), ', remainder ', &
                          remainder%get_facets_number(), ', of ', dragon%get_facets_number()
!endregion clip
call report('clip', surface)
!region merge
call surface%save_into_file(file_name='dragon_part_1.stl')
call remainder%save_into_file(file_name='dragon_part_2.stl')
call surface_1%load_from_file(file_name='dragon_part_1.stl', guess_format=.true.)
call surface_2%load_from_file(file_name='dragon_part_2.stl', guess_format=.true.)
call surface_1%merge_solids(other=surface_2)              ! merge_solids re-runs analyze
call surface_1%save_into_file(file_name='dragon_merged.stl')
print '(A,I0,A)', 'merge          ', surface_1%get_facets_number(), ' facets'
!endregion merge
call report('merge', surface_1)
contains
   subroutine report(label, surface)
   !< The box and the volume of a surface.
   character(*),             intent(in) :: label   !< What was done.
   type(surface_stl_object), intent(in) :: surface !< Surface.
   type(vector_R8P)                     :: bmin, bmax
   character(15)                        :: l

   l = label ; bmin = surface%get_bmin() ; bmax = surface%get_bmax()
   print '(A15,A,3F9.3,A,3F9.3,A,F10.3)', l, 'box ', bmin%x, bmin%y, bmin%z, ' ..', bmax%x, bmax%y, bmax%z, &
                                         '  volume', surface%get_volume()
   endsubroutine report
endprogram ref_usage_manipulation
