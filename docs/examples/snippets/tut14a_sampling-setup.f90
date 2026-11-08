use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_interpolate
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: dragon, pebble
type(rigid_pose_t)       :: p0, p1
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: scale, d, s, first, last
integer(I4P)             :: k, touching
logical                  :: hit, was

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! the dragon, placed as in chapter 3
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze

call pebble%load_from_file(file_name='cube-01.stl', guess_format=.true.)         ! a pebble: a cube of side 0.04
call pebble%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)                       ! centred on its own origin
call pebble%resize(x=0.04_R8P, y=0.04_R8P, z=0.04_R8P)
call pebble%analyze
! one step: thrown along x through the dragon's body, turning a quarter turn about y
p0 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
