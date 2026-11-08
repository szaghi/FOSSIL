!run figref-body-fitted-notch fig_ref_body_fitted_notch
!fig ref-body-fitted-notch
!< Figure of the reference page "Body-fitted grid support": hyperbolic extrusion (25 layers, 0.2 in all) from three
!< unit squares with a notch in the top side: a V notch 0.15 deep without smoothing (epsilon 0), the same with the
!< default smoothing, and a right-angled notch 0.6 wide and 0.2 deep with it. Cells not certified valid are filled.
program fig_ref_body_fitted_notch
use fossil,     only : curvilinear_block_2d_t, hyperbolic_grid_2d, geometric_spacing, CELL_VALID
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: su(80), sv(80), h(25), shift
integer(I4P), allocatable    :: cs(:,:)
integer(I4P)                 :: n, invalid, status, u

h = geometric_spacing(2e-3_R8P, 1.1_R8P, 25_I4P)
call fig_open(u, 'ref-body-fitted-notch', 'V notch: no smoothing, smoothing; right-angled notch: folded cells (ember)', &
              -90._R8P, 90._R8P, height=330_I4P)
call fig_zoom(u, vector_R8P(2.27_R8P, 0.85_R8P, 0._R8P), 1.25_R8P)
n = 0
call side(0._R8P, 0._R8P, 1._R8P, 0._R8P) ; call side(1._R8P, 0._R8P, 1._R8P, 1._R8P)
call side(1._R8P, 1._R8P, 0.5_R8P, 0.85_R8P) ; call side(0.5_R8P, 0.85_R8P, 0._R8P, 1._R8P)
call side(0._R8P, 1._R8P, 0._R8P, 0._R8P)
call hyperbolic_grid_2d(su(1:n), sv(1:n), h, b, epsilon=0._R8P, invalid=invalid, status=status)
print '(A,I0,A,I0)', 'V notch, epsilon 0:      status ', status, ', cells not certified valid ', invalid
call draw(0._R8P, 'V notch, epsilon 0')
call hyperbolic_grid_2d(su(1:n), sv(1:n), h, b, invalid=invalid, status=status)
print '(A,I0,A,I0)', 'V notch, epsilon 0.5:    status ', status, ', cells not certified valid ', invalid
call draw(1.8_R8P, 'V notch, epsilon 0.5')
n = 0
call side(0._R8P, 0._R8P, 1._R8P, 0._R8P) ; call side(1._R8P, 0._R8P, 1._R8P, 1._R8P)
call side(1._R8P, 1._R8P, 0.8_R8P, 1._R8P) ; call side(0.8_R8P, 1._R8P, 0.8_R8P, 0.8_R8P)
call side(0.8_R8P, 0.8_R8P, 0.2_R8P, 0.8_R8P) ; call side(0.2_R8P, 0.8_R8P, 0.2_R8P, 1._R8P)
call side(0.2_R8P, 1._R8P, 0._R8P, 1._R8P) ; call side(0._R8P, 1._R8P, 0._R8P, 0._R8P)
call hyperbolic_grid_2d(su(1:n), sv(1:n), h, b, invalid=invalid, status=status)
print '(A,I0,A,I0)', 'right-angled notch, 0.5: status ', status, ', cells not certified valid ', invalid
call draw(3.6_R8P, 'square notch, 0.5')
call fig_close(u)
contains
   subroutine side(x0, y0, x1, y1)
   !< Append 10 points of the side from (x0, y0) toward (x1, y1) (its end is the next side's start).
   real(R8P), intent(in) :: x0, y0, x1, y1 !< Ends.
   integer(I4P)          :: j

   do j = 0, 9
      n = n + 1 ; su(n) = x0 + (x1 - x0) * j / 10._R8P ; sv(n) = y0 + (y1 - y0) * j / 10._R8P
   enddo
   endsubroutine side

   subroutine draw(dx, title)
   !< The grid b shifted by dx: its lines, the wall, the cells not certified valid.
   real(R8P),    intent(in)      :: dx    !< Shift.
   character(*), intent(in)      :: title !< Label.
   type(vector_R8P), allocatable :: a(:), c(:)
   integer(I4P)                  :: i, j, m

   shift = dx
   call b%validity(cs)
   allocate(a(2 * b%ni * b%nj), c(2 * b%ni * b%nj)) ; m = 0
   do j = 2, b%nj ! the lines above v = 0.5
      do i = 1, b%ni - 1
         if (min(b%v(i, j), b%v(i + 1, j)) < 0.5_R8P) cycle
         m = m + 1 ; a(m) = p(i, j) ; c(m) = p(i + 1_I4P, j)
      enddo
   enddo
   do j = 1, b%nj - 1
      do i = 1, b%ni
         if (min(b%v(i, j), b%v(i, j + 1)) < 0.5_R8P) cycle
         m = m + 1 ; a(m) = p(i, j) ; c(m) = p(i, j + 1_I4P)
      enddo
   enddo
   call fig_segments(u, 'box', a(1:m), c(1:m))
   m = 0
   do j = 1, b%nj - 1 ! the edges of the cells not certified valid
      do i = 1, b%ni - 1
         if (cs(i, j) == CELL_VALID) cycle
         a(m + 1:m + 4) = [p(i, j), p(i + 1_I4P, j), p(i + 1_I4P, j + 1_I4P), p(i, j + 1_I4P)]
         c(m + 1:m + 4) = [p(i + 1_I4P, j), p(i + 1_I4P, j + 1_I4P), p(i, j + 1_I4P), p(i, j)]
         m = m + 4
      enddo
   enddo
   if (m > 0) call fig_segments(u, 'hole', a(1:m), c(1:m))
   call fig_segments(u, 'accent', [(p(i, 1_I4P), i = 1, b%ni - 1)], [(p(i, 1_I4P), i = 2, b%ni)])
   call fig_label(u, 'accent', vector_R8P(dx + 0.05_R8P, 0.2_R8P, 0._R8P), title)
   endsubroutine draw

   function p(i, j) result(q)
   !< Vertex (i, j) of b, shifted.
   integer(I4P), intent(in) :: i, j !< Indices.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(b%u(i, j) + shift, b%v(i, j), 0._R8P)
   endfunction p
endprogram fig_ref_body_fitted_notch
