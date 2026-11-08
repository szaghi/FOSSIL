!< Figure of the cookbook recipe "Flux-redistribution weights of a small cell": the disc of ck_frd on 16 x 16 cells,
!< zoomed on its smallest fluid piece; the cut cells coloured by fluid fraction, the cells of the piece's FRD
!< neighbourhood outlined, each labelled with its weight.
program fig_ck_frd
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, frd_t, frd_neighbourhoods
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P), H = 0.0625_R8P
type(curve_object)   :: disc
type(cut_cells_2d_t) :: cut
type(cut_links_t)    :: links
type(frd_t)          :: frd
type(vector_R8P)     :: a(48), b(48), c(4)
real(R8P)            :: u(48), v(48), x(17), x0, y0
integer(I4P)         :: i, j, k, p, m, un
character(10)        :: text

u = [(0.5_R8P + 0.3_R8P * cos(2 * PI * i / 48), i = 0, 47)]
v = [(0.5_R8P + 0.3_R8P * sin(2 * PI * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(H * i, i = 0, 16)]
call disc%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)
call frd_neighbourhoods(links, cut%piece_cell, cut%piece_area, x, x, frd)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)
call fig_open(un, 'ck-frd', 'FRD neighbourhood of the smallest piece (ember), weight per cell', &
              -90._R8P, 90._R8P, height=440_I4P)
call fig_zoom(un, vector_R8P(x(cut%piece_cell(1, p)) + H / 2, x(cut%piece_cell(2, p)) + H / 2, 0._R8P), 2.2_R8P * H)
do j = 1, 16 ; do i = 1, 16
   if (cut%alpha(i, j) <= 0._R8P) cycle
   call fig_quad(un, 'frac', vector_R8P(x(i), x(j), 0._R8P), vector_R8P(x(i + 1), x(j), 0._R8P), &
                 vector_R8P(x(i + 1), x(j + 1), 0._R8P), vector_R8P(x(i), x(j + 1), 0._R8P), value=1._R8P - cut%alpha(i, j))
enddo ; enddo
a = [(vector_R8P(u(i), v(i), 1e-3_R8P), i = 1, 48)] ; b = [a(2:48), a(1)]
call fig_segments(un, 'body', a, b)
do k = frd%start(p), frd%start(p + 1) - 1
   m = frd%member(k) ; x0 = x(cut%piece_cell(1, m)) ; y0 = x(cut%piece_cell(2, m))
   c = [vector_R8P(x0, y0, 2e-3_R8P), vector_R8P(x0 + H, y0, 2e-3_R8P), vector_R8P(x0 + H, y0 + H, 2e-3_R8P), &
        vector_R8P(x0, y0 + H, 2e-3_R8P)]
   call fig_segments(un, merge('hole  ', 'accent', m == p), c, [c(2:4), c(1)])
   write(text, '(F6.4)') frd%weight(k)
   call fig_label(un, merge('hole  ', 'accent', m == p), vector_R8P(x0 + 0.12_R8P * H, y0 + 0.4_R8P * H, 3e-3_R8P), &
                  trim(text))
enddo
call fig_close(un)
endprogram fig_ck_frd
