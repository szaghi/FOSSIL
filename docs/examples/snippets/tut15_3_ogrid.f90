!< Tutorial 15, step 3: a hyperbolic O-grid around the mid-span section of the wing, its first height set from a y1+
!< target. Figures: the whole O-grid; an animation of the march at the blunt trailing edge, two layers per frame.
program tut15_3_ogrid
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
do k = 10, 40, 10                                                      ! layer k: its height, the distance marched
   print '(A,I2,A,ES10.3,A,ES10.3)', 'layer ', k, ': height ', spacing(k), ', from the wall ', sum(spacing(1:k))
enddo

allocate(a(2 * grid%ni * grid%nj), b(2 * grid%ni * grid%nj))
call fig_open(un, 'tut15-ogrid', '121 x 41 O-grid: 40 layers from h1 = 2.2e-5 out to 0.16, the wall in yellow', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.25_R8P, 0._R8P, 0._R8P), 0.52_R8P)
call lines(41_I4P, m)
call fig_segments(un, 'box', a(1:m), b(1:m))
call fig_segments(un, 'accent', [(pt(i, 1), i = 1, grid%ni - 1)], [(pt(i + 1, 1), i = 1, grid%ni - 1)])
call fig_close(un)

call fig_open(un, 'tut15-march', '', -90._R8P, 90._R8P)               ! the march at the blunt trailing edge
call fig_zoom(un, vector_R8P(0.746_R8P, 0._R8P, 0._R8P), 0.014_R8P)
call fig_segments(un, 'accent', [(pt(i, 1), i = 1, grid%ni - 1)], [(pt(i + 1, 1), i = 1, grid%ni - 1)])
call fig_anim(un, 0.6_R8P)
do k = 2, 30, 2
   write(caption, '(A,I2,A,ES9.2,A,F7.4)') 'layer ', k, ': height ', spacing(k), ', from the wall ', sum(spacing(1:k))
   call fig_frame(un, caption=trim(caption), seconds=merge(2._R8P, 0.6_R8P, k == 30))
   do j = k - 1, k                                                     ! the two new layers, filled
      do i = 1, grid%ni - 1
         call fig_quad(un, 'frac', pt(i, j), pt(i + 1, j), pt(i + 1, j + 1), pt(i, j + 1), real(j, R8P) / 30._R8P)
      enddo
   enddo
   call lines(k + 1_I4P, m)
   call fig_segments(un, 'box', a(1:m), b(1:m))
enddo
call fig_close(un)
contains
   function pt(i, j) result(r)
   !< Vertex (i, j) of the O-grid, in the plane z = 0.
   integer(I4P), intent(in) :: i, j !< Vertex.
   type(vector_R8P)         :: r    !< Point.

   r = vector_R8P(grid%u(i, j), grid%v(i, j), 0._R8P)
   endfunction pt

   subroutine lines(nj, n)
   !< The grid lines of the layers 1..nj into a, b.
   integer(I4P), intent(in)  :: nj !< Layers drawn.
   integer(I4P), intent(out) :: n  !< Segments.
   integer(I4P)              :: i, j

   n = 0
   do j = 1, nj
      do i = 1, grid%ni - 1
         n = n + 1 ; a(n) = pt(i, j) ; b(n) = pt(i + 1, j)
      enddo
   enddo
   do i = 1, grid%ni - 1
      do j = 1, nj - 1
         n = n + 1 ; a(n) = pt(i, j) ; b(n) = pt(i, j + 1)
      enddo
   enddo
   endsubroutine lines
endprogram tut15_3_ogrid
