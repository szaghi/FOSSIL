!< Figures of the cookbook recipe "Coarsen cut cells": the disc of ck_cut_coarsen, seen from above, its cells coloured
!< by the fluid fraction alpha: the fine 32 x 32 cut cells, then the 16 x 16 ones cut_coarsen_2d builds from them.
program fig_ck_cut_coarsen
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_coarsen_2d
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curve_object)        :: disc
type(cut_cells_2d_t)      :: fine, coarse
type(cut_links_t)         :: links, coarse_links
real(R8P)                 :: u(48), v(48), xf(33), xc(17)
integer(I4P), allocatable :: fine_to_coarse(:)
integer(I4P)              :: i, status

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
xf = [(i / 32._R8P, i = 0, 32)] ; xc = [(i / 16._R8P, i = 0, 16)]
call disc%cut_cells(x=xf, y=xf, cut=fine)
call cut_links_2d(fine, links)
call cut_coarsen_2d(fine, links, xf, xf, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
call draw('ck-cut-coarsen-fine', 'the fine cut cells, 32 x 32, by fluid fraction', xf, fine)
call draw('ck-cut-coarsen-coarse', 'cut_coarsen_2d, ratio 2: 16 x 16', xc, coarse)
contains
   subroutine draw(name, caption, x, cut)
   !< One figure: the cells coloured by alpha, the grid lines and the disc's polygon.
   character(*),         intent(in) :: name, caption !< Figure name, caption.
   real(R8P),            intent(in) :: x(:)          !< Cell edges (both directions).
   type(cut_cells_2d_t), intent(in) :: cut           !< Cut cells.
   type(vector_R8P), allocatable    :: a(:), b(:)    !< Segments.
   integer(I4P)                     :: i, j, n, un   !< Work.

   n = size(x) - 1
   call fig_open(un, name, caption, -90._R8P, 90._R8P, width=420_I4P, height=420_I4P)
   do j = 1, n ; do i = 1, n
      call fig_quad(un, 'frac', vector_R8P(x(i), x(j), 0._R8P), vector_R8P(x(i + 1), x(j), 0._R8P), &
                    vector_R8P(x(i + 1), x(j + 1), 0._R8P), vector_R8P(x(i), x(j + 1), 0._R8P), value=max(cut%alpha(i, j), 0._R8P))
   enddo ; enddo
   a = [(vector_R8P(x(i), 0._R8P, 1.e-3_R8P), i = 1, n + 1), (vector_R8P(0._R8P, x(i), 1.e-3_R8P), i = 1, n + 1)]
   b = [(vector_R8P(x(i), 1._R8P, 1.e-3_R8P), i = 1, n + 1), (vector_R8P(1._R8P, x(i), 1.e-3_R8P), i = 1, n + 1)]
   call fig_segments(un, 'box', a, b)
   a = [(vector_R8P(u(i), v(i), 2.e-3_R8P), i = 1, 48)]
   b = [(vector_R8P(u(mod(i, 48) + 1), v(mod(i, 48) + 1), 2.e-3_R8P), i = 1, 48)]
   call fig_segments(un, 'accent', a, b)
   call fig_close(un)
   endsubroutine draw
endprogram fig_ck_cut_coarsen
