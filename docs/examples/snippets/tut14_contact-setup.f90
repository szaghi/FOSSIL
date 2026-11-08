use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, pose_2d, &
                   pose_interpolate, first_impact, first_impact_2d, slice_surface, CCD_CONTACT, CCD_SEPARATE, &
                   CCD_NOT_CONVERGED
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter          :: PI = acos(-1._R8P)
type(surface_stl_object)      :: dragon, pebble, naca, bodies(4)
type(curve_object)            :: wing, wall, disc, curves(4)
type(rigid_pose_t)            :: p0, p1, poses0(4), poses1(4)
type(rigid_pose_2d_t)         :: q0, q1, curve_poses0(4), curve_poses1(4)
type(vector_R8P)              :: bmin, bmax
real(R8P)                     :: scale, d, toi, toi_all, pa(3), pb(3), qa(2), s
integer(I4P)                  :: i, j, k, fa, fb, steps, pairs, status, pair(2), tested, sa, sb
logical                       :: hit

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
