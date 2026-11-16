program ck_loads_bodies
use fossil, only : curve_object, surface_stl_object, load_samples_t, load_samples_2d_bodies, load_samples_3d_bodies, &
                   integrate_loads, load_samples_reset, pose_2d, pose_from_axis_angle
use penf,   only : I4P, R8P
implicit none
type(curve_object)       :: bodies(2)
type(surface_stl_object) :: cubes(2)
type(load_samples_t)     :: ls
real(R8P), allocatable   :: fb(:,:), mb(:)
real(R8P)                :: u(64), v(64), force(2), moment
integer(I4P)             :: i

u = [(0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]    ! a circle of radius 0.3, counter-clockwise
v = [(0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call bodies(1)%adopt_loops(u=u, v=v, loop_sizes=[64_I4P]) ; call bodies(2)%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
call load_samples_2d_bodies(bodies, samples=ls, lengths=[0.1_R8P], &   ! one probe 0.1 out of every sample
                            poses=[pose_2d(0._R8P, [0.42_R8P, 0.5_R8P]), pose_2d(0._R8P, [1.08_R8P, 0.5_R8P])]) ! 0.06 apart
print '(A,I0,A,I0,A,F7.4)', 'two circles: samples ', ls%n, ', probes blocked by the other body ', &
      count(.not. ls%visible), ', least clearance ', minval(ls%clearance)
call integrate_loads(ls, force, moment, pressure=1._R8P - ls%points(2, :), &  ! a fluid at rest, gravity -v
                     per_body_force=fb, per_body_moment=mb, about=[0.75_R8P, 0.5_R8P])
print '(A,F9.5,A,F9.5,A,F9.5)', '  buoyancy (v) of body 1', fb(2, 1), ', of body 2', fb(2, 2), ', area', &
      bodies(1)%get_area()
call cubes(1)%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube, twice
call cubes(2)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call load_samples_3d_bodies(cubes, samples=ls, lengths=[0.1_R8P, 0.5_R8P], &
                            poses=[pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P), &
                                   pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [1.25_R8P, 0._R8P, 0._R8P])])
print '(A,I0,A,I0,A,F6.3)', 'two cubes 0.25 apart: samples ', ls%n, ', probes blocked ', count(.not. ls%visible), &
      ', least clearance ', minval(ls%clearance)
call load_samples_reset(ls)
endprogram ck_loads_bodies
