program ck_ghost_points
use fossil, only : surface_stl_object, curve_object, ghost_points_t, ghost_points_2d, ghost_points_3d, GHOST_XM
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: h = 0.05_R8P
type(curve_object)       :: body
type(surface_stl_object) :: cube
type(ghost_points_t)     :: gp
real(R8P)                :: u(64), v(64), xc(20), xn(9)
integer(I4P)             :: i

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]   ! an ellipse, counter-clockwise
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xc = [((i - 0.5_R8P) * h, i = 1, 20)]                      ! the solver's cell centres, h = 0.05
call ghost_points_2d(body, xc, xc, [2_I4P, 2_I4P], gp)    ! stencil half-width 2: two ghost layers
print '(A,I0,A,I0)', 'ellipse: ghost points ', gp%n, ', of them seeing fluid at -u ', count(iand(gp%mask, GHOST_XM) /= 0)
print '(A,2I3,A,4I2)', 'ghost 1 at ij', gp%ijk(:, 1), ', steps to fluid (-u +u -v +v)', gp%steps(:, 1)
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube
xn = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]                ! node lines through its faces
call ghost_points_3d(cube, xn, xn, xn, [1_I4P, 1_I4P, 1_I4P], gp, device=.true.)
print '(A,I0,A,F5.2,A,F5.2)', 'cube: ghost points ', gp%n, ', x from ', minval(gp%points(1, :)), ' to ', maxval(gp%points(1, :))
endprogram ck_ghost_points
