!run ck_frd ck_frd
program ck_frd
use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, frd_t, frd_neighbourhoods
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: disc
type(cut_cells_2d_t) :: cut
type(cut_links_t)    :: links
type(frd_t)          :: frd
real(R8P)            :: u(48), v(48), x(17)
integer(I4P)         :: i, p, m

u = [(0.5_R8P + 0.3_R8P * cos(2 * PI * i / 48), i = 0, 47)]   ! a disc, 48 segments
v = [(0.5_R8P + 0.3_R8P * sin(2 * PI * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(0.0625_R8P * i, i = 0, 16)]                             ! 16 x 16 cells
call disc%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)
call frd_neighbourhoods(links, cut%piece_cell, cut%piece_area, x, x, frd)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)   ! the smallest piece
print '(A,I0,A,2I3,A,F7.4)', 'smallest piece ', p, ' in cell', cut%piece_cell(:, p), ', area fraction ', &
                             cut%piece_area(p) / 0.0625_R8P**2
do i = frd%start(p), frd%start(p + 1) - 1                     ! the piece first, then its neighbours
   m = frd%member(i)
   print '(A,I4,A,2I3,A,F7.4,A,F7.4)', '  member ', m, ' in cell', cut%piece_cell(:, m), ', kappa ', &
                                       cut%piece_area(m) / 0.0625_R8P**2, ', weight ', frd%weight(i)
enddo
print '(A,F12.9)', '  sum of the weights ', sum(frd%weight(frd%start(p):frd%start(p + 1) - 1))
endprogram ck_frd
