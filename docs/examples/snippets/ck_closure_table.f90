program ck_closure_table
use fossil, only : surface_stl_object, curve_object, ib_closure_t, fd_closure_2d, fd_closure_3d, ib_closure_reset
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: h = 0.05_R8P
type(curve_object)       :: body
type(surface_stl_object) :: cube
type(ib_closure_t)       :: table
real(R8P)                :: u(64), v(64), xc(20), xn(9)
integer(I4P)             :: i

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]   ! an ellipse, counter-clockwise
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xc = [((i - 0.5_R8P) * h, i = 1, 20)]                     ! the solver's cell centres
call fd_closure_2d(body, xc, xc, [2_I4P, 2_I4P], table, lengths=[sqrt(2._R8P) * h, 2 * sqrt(2._R8P) * h], &
                   angular_velocity=1._R8P, centre=[0.51_R8P, 0.52_R8P], ilw=.true.)
print '(A,I0,A,I0,A,I0)', 'ellipse: ghosts ', table%n, ', visible probes ', count(table%visible), ', ILW crossings ', &
      count(table%ilw_element > 0)
print '(A,2F8.4,A,2F8.4,A,2F8.4)', 'ghost 1: foot', table%foot(:, 1), '  normal', table%normal(:, 1), '  wall velocity', &
      table%velocity(:, 1)
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube
call cube%enter_device
xn = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]                ! node lines through its faces
call fd_closure_3d(cube, xn, xn, xn, [1_I4P, 1_I4P, 1_I4P], table, ilw=.true., resident=.true.) ! left on the device
call table%update_host                                    ! only to print it here
print '(A,I0,A,L1,A,F6.3,A,F6.3)', 'cube: ghosts ', table%n, ', resident ', table%on_device, ', ILW x from ', &
      minval(table%ilw_t(1:2, :), mask=table%ilw_element(1:2, :) > 0) + 0._R8P, ' to ', &
      maxval(table%ilw_t(1:2, :), mask=table%ilw_element(1:2, :) > 0)
call ib_closure_reset(table)                              ! releases the device copies
call cube%exit_device
endprogram ck_closure_table
