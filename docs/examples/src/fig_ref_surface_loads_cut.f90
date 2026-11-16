!run figref-surface-loads-cut fig_ref_surface_loads_cut
!fig ref-surface-loads-cut
!< Figure of the reference page "Surface loads", section "Loads from the cut cells" (`cut_walls_2d`, `cut_loads_2d`): a
!< NACA 0012 of unit chord (60 vertices per side) at 20 degrees nose up about its leading edge, cut by 14 x 8 cells of 0.1
!< (grey). The wall segments of the cut pieces in gold, each with its contribution -p A_w in the pressure p = -v of a fluid
!< at rest (ember, drawn ending on the wall at the segment's midpoint), the total from the centroid (teal): the buoyancy
!< (0, A). The output compares the force and the moment about the leading edge with the buoyancy's, with and without
!< the pressure gradient on each wall.
program fig_ref_surface_loads_cut
use fossil,     only : curve_object, cut_cells_2d_t, cut_walls_2d_t, cut_walls_2d, cut_loads_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter         :: M = 60, NX = 14, NY = 8
real(R8P),    parameter         :: PI = acos(-1._R8P), H = 0.1_R8P, LIFT = 1.e-3_R8P, SF = 2._R8P, SW = 5._R8P
type(curve_object)              :: naca
type(cut_cells_2d_t)            :: cut
type(cut_walls_2d_t)            :: walls
type(vector_R8P), allocatable   :: a(:), b(:)
real(R8P),        allocatable   :: p(:), g(:,:)
real(R8P)                       :: u(2 * M), v(2 * M), x(NX + 1), y(NY + 1), force(2), m0, m1, area, xc(2), xi, ca, sa, s
integer(I4P)                    :: i, k

k = 0
do i = M, 0, -1                                              ! trailing edge, upper side, leading edge
   xi = 0.5_R8P * (1._R8P - cos(PI * i / M)) ; k = k + 1 ; u(k) = xi ; v(k) = yt(xi)
enddo
do i = 1, M - 1                                              ! lower side
   xi = 0.5_R8P * (1._R8P - cos(PI * i / M)) ; k = k + 1 ; u(k) = xi ; v(k) = -yt(xi)
enddo
ca = cos(20._R8P * PI / 180._R8P) ; sa = sin(20._R8P * PI / 180._R8P)
do k = 1, 2 * M                                              ! nose up about the leading edge
   s = u(k) ; u(k) = ca * s + sa * v(k) ; v(k) = -sa * s + ca * v(k)
enddo
call naca%adopt_loops(u=u, v=v, loop_sizes=[2 * M])
x = [(-0.13_R8P + H * i, i = 0, NX)] ; y = [(-0.47_R8P + H * i, i = 0, NY)]
call naca%cut_cells(x, y, cut, boundaries=.true.)            ! the cut pieces' oriented boundaries
call cut_walls_2d(cut, walls)                                ! their wall segments
p = -walls%centroid(2, :)                                    ! the pressure at the midpoints: a fluid at rest, p = -v
g = spread([0._R8P, -1._R8P], 2, walls%n)                    ! its gradient on every wall
call cut_loads_2d(walls, force, m0, pressure=p, about=[0._R8P, 0._R8P])
call cut_loads_2d(walls, force, m1, pressure=p, gradient=g, about=[0._R8P, 0._R8P])
call naca%mass_properties(area=area, center=xc)

call fig_open(k, 'ref-surface-loads-cut', &
              'NACA 0012 on 14 x 8 cells, p = -v: walls (gold), -p A_w (ember), F (teal)', &
              -90._R8P, 90._R8P)
allocate(a(NX + NY + 2), b(NX + NY + 2))
do i = 1, NX + 1
   a(i) = vec(x(i), y(1)) ; b(i) = vec(x(i), y(NY + 1))
enddo
do i = 1, NY + 1
   a(NX + 1 + i) = vec(x(1), y(i)) ; b(NX + 1 + i) = vec(x(NX + 1), y(i))
enddo
call fig_segments(k, 'box', a, b)
deallocate(a, b) ; allocate(a(walls%n), b(walls%n))
do i = 1, walls%n
   a(i) = vec(walls%ends(1, i), walls%ends(2, i), LIFT) ; b(i) = vec(walls%ends(3, i), walls%ends(4, i), LIFT)
enddo
call fig_segments(k, 'foot', a, b)
do i = 1, walls%n                                            ! -p A_w, drawn ending on the wall
   a(i) = vec(walls%centroid(1, i) + SW * p(i) * walls%area(1, i), walls%centroid(2, i) + SW * p(i) * walls%area(2, i), &
              2 * LIFT)
   b(i) = vec(walls%centroid(1, i), walls%centroid(2, i), 2 * LIFT)
enddo
call fig_segments(k, 'patch', a, b)
call fig_segments(k, 'probe', [vec(xc(1), xc(2), 3 * LIFT)], [vec(xc(1) + SF * force(1), xc(2) + SF * force(2), 3 * LIFT)])
call fig_points(k, 'probe', [vec(xc(1), xc(2), 3 * LIFT)])
call fig_label(k, 'probe', vec(xc(1) + SF * force(1) + 0.02_R8P, xc(2) + SF * force(2), 3 * LIFT), 'F')
call fig_close(k)
print '(A,I0,A,I0,A,F9.6)', 'cut pieces ', cut%np, ', wall segments ', walls%n, ', area A ', area
print '(A,2F11.7,A,2F11.7)', 'force', clean(force), '   buoyancy (0, A)', clean([0._R8P, area])
print '(A,F11.7,A,F11.7,A,F11.7)', 'moment about the leading edge: centroid values', m0, ', with the gradient', m1, &
      ', A x_c', area * xc(1)
print '(A,ES8.1,A,L1)', 'moment error: centroid values', abs(m0 - area * xc(1)), ';  with the gradient below 1e-15: ', &
      abs(m1 - area * xc(1)) < 1.e-15_R8P .and. abs(force(1)) < 1.e-15_R8P .and. abs(force(2) - area) < 1.e-15_R8P

contains
   function vec(xx, yy, z) result(pt)
   !< A point of the drawing plane.
   real(R8P), intent(in)           :: xx, yy !< Coordinates.
   real(R8P), intent(in), optional :: z      !< Lift.
   type(vector_R8P)                :: pt     !< Point.

   pt = vector_R8P(xx, yy, 0._R8P)
   if (present(z)) pt%z = z
   endfunction vec

   pure function clean(q) result(r)
   !< Rounded to 7 decimals, no negative zero.
   real(R8P), intent(in) :: q(:)       !< Values.
   real(R8P)             :: r(size(q)) !< Rounded.

   r = anint(q * 1.e7_R8P) / 1.e7_R8P + 0._R8P
   endfunction clean

   pure real(R8P) function yt(xx)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: xx !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.1260_R8P * xx - 0.3516_R8P * xx ** 2 + 0.2843_R8P * xx ** 3 - &
                   0.1036_R8P * xx ** 4)
   endfunction yt
endprogram fig_ref_surface_loads_cut
