!run ck_closure_bodies ck_closure_bodies
program ck_closure_bodies
use fossil, only : surface_stl_object, curve_object, ib_closure_t, fd_closure_2d_bodies, fd_closure_3d_bodies, &
                   ib_closure_reset, rigid_pose_2d_t, pose_2d, rigid_pose_t, pose_from_axis_angle
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: h = 0.05_R8P
type(curve_object)       :: bodies(2)
type(surface_stl_object) :: cubes(2)
type(ib_closure_t)       :: table
real(R8P)                :: u(64), v(64), xc(30), yc(20), xn(13), yn(9)
integer(I4P)             :: i

u = [(0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]    ! a circle of radius 0.3, counter-clockwise
v = [(0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call bodies(1)%adopt_loops(u=u, v=v, loop_sizes=[64_I4P]) ; call bodies(2)%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xc = [((i - 0.5_R8P) * h, i = 1, 30)] ; yc = [((i - 0.5_R8P) * h, i = 1, 20)]   ! the solver's cell centres
call fd_closure_2d_bodies(bodies, xc, yc, [2_I4P, 2_I4P], table, &
                          poses=[pose_2d(0._R8P, [0.42_R8P, 0.5_R8P]), pose_2d(0._R8P, [1.08_R8P, 0.5_R8P])], & ! 0.06 apart
                          lengths=[2 * h], angular_velocities=[1._R8P, -1._R8P])  ! each body spins about its centre
print '(A,I0,A,I0,A,I0,A,I0)', 'two circles: ghosts ', table%n, ', shared ', count(table%shared), &
      ', probes blocked by the other body ', count(.not. table%visible), ', overlapping ', count(table%overlap)
print '(A,F6.3,A,F7.4)', '  shared below ', table%shared_distance, ', least clearance ', minval(table%clearance)
call cubes(1)%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube, twice
call cubes(2)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
xn = [(-0.5_R8P + 0.25_R8P * i, i = 0, 12)] ; yn = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call fd_closure_3d_bodies(cubes, xn, yn, yn, [2_I4P, 2_I4P, 2_I4P], table, ilw=.true., &
                          poses=[pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P), &
                                 pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [1.25_R8P, 0._R8P, 0._R8P])])
print '(A,I0,A,I0,A,I0,A,I0)', 'two cubes 0.25 apart: ghosts ', table%n, ' (', count(table%body == 1), ' + ', &
      count(table%body == 2), '), shared ', count(table%shared)
call ib_closure_reset(table)
endprogram ck_closure_bodies
