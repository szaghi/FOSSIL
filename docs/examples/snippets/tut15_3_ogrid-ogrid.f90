use fossil,     only : surface_stl_object, curve_object, curvilinear_block_2d_t, slice_surface, hyperbolic_grid_2d, &
                       geometric_spacing, utau_flat_plate
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_quad, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)     :: wing
type(curve_object)           :: section
type(curvilinear_block_2d_t) :: grid
real(R8P),        allocatable :: spacing(:), u(:), v(:)
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                    :: nu, u_inf, u_tau, h1, p(2), q(2)
integer(I4P)                 :: i, j, k, m, first, last, invalid, status, un
character(96)                :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)  ! the mid-span section, (u, v) = (x, y)
call section%get_loop_range(1_I4P, first, last)
allocate(u(2 * (last - first + 1)), v(2 * (last - first + 1)))
do i = first, last                                                     ! each segment split in 2: the same polygon
   p = section%get_vertex(i) ; q = section%get_vertex(merge(first, i + 1_I4P, i == last))
   do j = 0, 1
      u(2 * (i - first) + j + 1) = p(1) + (q(1) - p(1)) * j / 2._R8P
      v(2 * (i - first) + j + 1) = p(2) + (q(2) - p(2)) * j / 2._R8P
   enddo
enddo
print '(A,I0,A,I0,A)', 'section: ', last - first + 1, ' vertices, ', size(u), ' wall points'
nu = 1.e-6_R8P ; u_inf = 1._R8P                                        ! chord 1, Re = u_inf c / nu = 1e6
u_tau = utau_flat_plate(u_inf, nu, 0.5_R8P)                            ! flat-plate estimate at mid-chord
h1 = 1._R8P * nu / u_tau                                               ! y1+ = 1 by design: y1 = y1+ nu / u_tau
spacing = geometric_spacing(h1, 1.2_R8P, 40_I4P)                       ! 40 layers, ratio 1.2
print '(A,ES10.3,A,ES10.3,A,ES10.3,A,F7.4)', 'u_tau ', u_tau, ', h1 ', h1, ', last layer ', spacing(40), &
                                             ', outer height ', sum(spacing)
call hyperbolic_grid_2d(u, v, spacing, grid, invalid=invalid, status=status)
print '(A,I0,A,I0,A,I0,A,I0)', 'O-grid ', grid%ni, ' x ', grid%nj, ', status ', status, &
                               ', cells not certified valid ', invalid
