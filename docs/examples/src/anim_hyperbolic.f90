!run anim-hyperbolic anim_hyperbolic
!fig anim-hyperbolic
!< Animation of the body-fitted reference page and of the landing page: an O-grid marched out of a NACA 0012 section
!< by hyperbolic_grid_2d (24 layers, first height 4e-3 chord, growth 1.15). Each frame adds three layers of the
!< computed block: the new layers are filled (darker outwards), the grid already marched is drawn in lines; the caption
!< gives the layer, its height, the distance marched from the wall and the cells not certified valid.
program anim_hyperbolic
use fossil,     only : curvilinear_block_2d_t, hyperbolic_grid_2d, geometric_spacing
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_segments, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter      :: N = 48, NL = 24, STEP = 3
real(R8P),    parameter      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: g
type(vector_R8P)             :: p((N + 1) * (NL + 1) * 2), q((N + 1) * (NL + 1) * 2)
real(R8P)                    :: u(N), v(N), xc, yt, h(NL)
integer(I4P)                 :: i, j, k, m, un, invalid, status
character(96)                :: caption

do i = 1, N                                                    ! the section, closed, cosine-clustered at the ends
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / real(N, R8P)))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= N / 2 + 1)
enddo
h = geometric_spacing(4e-3_R8P, 1.15_R8P, NL)
call hyperbolic_grid_2d(u, v, h, g, invalid=invalid, status=status)
print '(A,I0,A,I0,A,I0,A,I0)', 'block ', g%ni, ' x ', g%nj, ', invalid cells ', invalid, ', status ', status

call fig_open(un, 'anim-hyperbolic', '', -90._R8P, 90._R8P, width=720_I4P, height=460_I4P)
call fig_zoom(un, vector_R8P(0.83_R8P, 0._R8P, 0._R8P), 0.95_R8P)
call fig_anim(un, 0.5_R8P)
do k = STEP, NL, STEP
   write(caption, '(A,I2,A,F6.4,A,F6.4,A,I0)') 'layer ', k, ': height ', h(k), ', marched ', sum(h(1:k)), &
                                              ', invalid cells ', invalid
   print '(A)', trim(caption)
   call fig_frame(un, caption=trim(caption), seconds=merge(1.5_R8P, 0.5_R8P, k == NL))
   do j = k - STEP + 1, k                                      ! the new layers, filled
      do i = 1, g%ni - 1
         call fig_quad(un, 'frac', pt(i, j), pt(i + 1, j), pt(i + 1, j + 1), pt(i, j + 1), real(j, R8P) / NL)
      enddo
   enddo
   m = 0                                                       ! the grid marched so far, in lines
   do j = 1, k + 1
      do i = 1, g%ni - 1
         m = m + 1 ; p(m) = pt(i, j) ; q(m) = pt(i + 1, j)
      enddo
   enddo
   do j = 1, k
      do i = 1, g%ni - 1
         m = m + 1 ; p(m) = pt(i, j) ; q(m) = pt(i, j + 1)
      enddo
   enddo
   call fig_segments(un, 'body', p(1:m), q(1:m))
   m = 0
   do i = 1, g%ni - 1
      m = m + 1 ; p(m) = pt(i, 1) ; q(m) = pt(i + 1, 1)
   enddo
   call fig_segments(un, 'accent', p(1:m), q(1:m))
enddo
call fig_close(un)
contains
   function pt(i, j) result(r)
   !< Vertex (i, j) of the block, in the plane z = 0.
   integer(I4P), intent(in) :: i, j !< Vertex.
   type(vector_R8P)         :: r    !< Point.

   r = vector_R8P(g%u(i, j), g%v(i, j), 0._R8P)
   endfunction pt
endprogram anim_hyperbolic
