!run ck_cut_loads ck_cut_loads
program ck_cut_loads
use fossil, only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t, cut_walls_gyration, cut_loads_3d, &
                   curve_object, cut_cells_2d_t, cut_walls_2d_t, cut_walls_2d, cut_loads_2d
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube
type(cut_cells_3d_t)     :: cut
type(cut_walls_3d_t)     :: walls
type(curve_object)       :: disc
type(cut_cells_2d_t)     :: cut2
type(cut_walls_2d_t)     :: walls2
real(R8P), allocatable   :: p(:), g(:,:), gyration(:,:)
real(R8P)                :: x(7), force(3), moment(3), f2(2), m2
integer(I4P)             :: i

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube [0, 1]^3
x = [(-0.13_R8P + 0.21_R8P * i, i = 0, 6)]                ! 6^3 cells around it
call cut_cells_3d(cube, x, x, x, cut, walls=walls)          ! the cut cells and their wall terms
p = 1._R8P - walls%centroid(3, :)                           ! the solver's pressure at the walls: p = 1 - z
g = spread([0._R8P, 0._R8P, -1._R8P], 2, walls%n)           ! and its gradient (for the exact moment)
call cut_walls_gyration(walls, cube, x, x, x, gyration)
call cut_loads_3d(walls, force, moment, pressure=p, gradient=g, gyration=gyration)   ! about the origin
print '(A,I0,A,3F9.5,A,3F9.5)', 'cube, ', walls%n, ' wall terms: force', r5(force), '  moment', r5(moment)
call disc%adopt_loops(u=[(0.5_R8P + 0.4_R8P * cos(2 * PI * i / 64), i = 0, 63)], &
                      v=[(0.5_R8P + 0.4_R8P * sin(2 * PI * i / 64), i = 0, 63)], loop_sizes=[64_I4P])
call disc%cut_cells(x, x, cut2, boundaries=.true.)          ! 2-D: the pieces' oriented boundaries
call cut_walls_2d(cut2, walls2)                             ! their wall segments
call cut_loads_2d(walls2, f2, m2, pressure=1._R8P - walls2%centroid(2, :), &   ! p = 1 - v, per unit depth
                  gradient=spread([0._R8P, -1._R8P], 2, walls2%n))
print '(A,I0,A,2F9.5,A,F9.5)', 'disc, ', walls2%n, ' wall segments: force', r5(f2), '  moment', r5([m2])

contains
   pure function r5(x) result(y)
   !< Rounded to 5 decimals, no negative zero.
   real(R8P), intent(in) :: x(:)       !< Values.
   real(R8P)             :: y(size(x)) !< Rounded.

   y = anint(x * 1.e5_R8P) / 1.e5_R8P + 0._R8P
   endfunction r5
endprogram ck_cut_loads
