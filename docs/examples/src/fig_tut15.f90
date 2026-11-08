!run fig15 fig_tut15
!fig tut15-face
!fig tut15-ogrid
!< Figures of tutorial 15. First the NACA0012 wing near its tip, the chapter's 41 x 21 block face before the projection
!< (grey lines, 0.005 above the upper side) and after project_block_face (yellow lines): the
!< trailing-edge row on the upper trailing-edge line, the tip row on z = 1, the corner on the wing's corner. Then the
!< hyperbolic O-grid around the mid-span section (h1 from y1+ = 1, ratio 1.2, 40 layers): its grid lines, the wall in
!< yellow.
program fig_tut15
use fossil,     only : surface_stl_object, curve_object, projection_target_t, curvilinear_block_t, &
                       curvilinear_block_2d_t, project_block_face, slice_surface, hyperbolic_grid_2d, geometric_spacing, &
                       utau_flat_plate
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(surface_stl_object)     :: wing, shown
type(projection_target_t)    :: target
type(curvilinear_block_t)    :: block, block0
type(curve_object)           :: section
type(curvilinear_block_2d_t) :: g
type(vector_R8P)             :: w0(120), w1(120), l0(10000), l1(10000)
real(R8P)                    :: x(41, 21, 2), y(41, 21, 2), z(41, 21, 2), xc, yt, a(2), b(2), h1
real(R8P), allocatable       :: u(:), v(:)
integer(I4P)                 :: i, j, un, first, last, n

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%projection_target(target)
do j = 1, 21
   do i = 1, 41
      xc = 0.5_R8P * (1._R8P - cos(PI * (i - 1) / 40._R8P))
      yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
      x(i, j, :) = xc - 0.25_R8P ; y(i, j, :) = [yt + 0.005_R8P, yt + 0.1_R8P] ; z(i, j, :) = (j - 1) / 20._R8P
   enddo
enddo
call block%initialize(x, y, z)
block0 = block
call project_block_face(block, face=5_I4P, target=target, tolerance=0.01_R8P)
shown = wing
call shown%rotate(axis=ex_R8P, angle=PI / 2._R8P)           ! drawn with the span horizontal: (x, y, z) -> (x, -z, y)
call fig_open(un, 'tut15-face', 'project_block_face: the face grid on the wing near its tip', -35._R8P, 35._R8P)
call fig_zoom(un, vector_R8P(0.45_R8P, -0.8_R8P, 0._R8P), 0.33_R8P)
call fig_surface(un, 'ghost', shown)
call grid_lines(un, 'box', block0)
call grid_lines(un, 'foot', block)
call fig_close(un)

call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)
call section%get_loop_range(1_I4P, first, last)
allocate(u(2 * (last - first + 1)), v(2 * (last - first + 1)))
do i = first, last
   a = section%get_vertex(i) ; b = section%get_vertex(merge(first, i + 1_I4P, i == last))
   do j = 0, 1
      u(2 * (i - first) + j + 1) = a(1) + (b(1) - a(1)) * j / 2._R8P
      v(2 * (i - first) + j + 1) = a(2) + (b(2) - a(2)) * j / 2._R8P
   enddo
enddo
h1 = 1.e-6_R8P / utau_flat_plate(1._R8P, 1.e-6_R8P, 0.5_R8P)
call hyperbolic_grid_2d(u, v, geometric_spacing(h1, 1.2_R8P, 40_I4P), g)
call fig_open(un, 'tut15-ogrid', 'hyperbolic_grid_2d: 121 x 41 O-grid, y1+ = 1 by design', -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.25_R8P, 0._R8P, 0._R8P), 0.52_R8P)
n = 0
do j = 1, g%nj
   do i = 1, g%ni - 1
      n = n + 1 ; l0(n) = p(i, j) ; l1(n) = p(i + 1, j)
   enddo
enddo
do i = 1, g%ni - 1
   do j = 1, g%nj - 1
      n = n + 1 ; l0(n) = p(i, j) ; l1(n) = p(i, j + 1)
   enddo
enddo
call fig_segments(un, 'box', l0(1:n), l1(1:n))
do i = 1, 120
   w0(i) = p(i, 1) ; w1(i) = p(i + 1, 1)
enddo
call fig_segments(un, 'accent', w0, w1)
call fig_close(un)
contains
   function p(i, j)
   !< Vertex (i, j) of the O-grid.
   integer(I4P), intent(in) :: i, j
   type(vector_R8P)         :: p

   p = vector_R8P(g%u(i, j), g%v(i, j), 0._R8P)
   endfunction p

   subroutine grid_lines(un, style, bl)
   !< The lines of the face k = 1 of a block, drawn with the span horizontal: (x, y, z) -> (x, -z, y).
   integer(I4P),              intent(in) :: un
   character(*),              intent(in) :: style
   type(curvilinear_block_t), intent(in) :: bl
   type(vector_R8P)                      :: a(1660), b(1660)
   integer(I4P)                          :: i, j, n

   n = 0
   do j = 1, 21
      do i = 1, 40
         n = n + 1 ; a(n) = vector_R8P(bl%x(i, j, 1), -bl%z(i, j, 1), bl%y(i, j, 1))
         b(n) = vector_R8P(bl%x(i + 1, j, 1), -bl%z(i + 1, j, 1), bl%y(i + 1, j, 1))
      enddo
   enddo
   do i = 1, 41
      do j = 1, 20
         n = n + 1 ; a(n) = vector_R8P(bl%x(i, j, 1), -bl%z(i, j, 1), bl%y(i, j, 1))
         b(n) = vector_R8P(bl%x(i, j + 1, 1), -bl%z(i, j + 1, 1), bl%y(i, j + 1, 1))
      enddo
   enddo
   call fig_segments(un, style, a(1:n), b(1:n))
   endsubroutine grid_lines
endprogram fig_tut15
