!run ck_cut_solver ck_cut_solver
program ck_cut_solver
use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, merge_t, &
                   merge_targets, SRD_REGULAR, MERGE_OK
use penf,   only : I4P, R8P
implicit none
type(curve_object)   :: disc
type(cut_cells_2d_t) :: cut
type(cut_links_t)    :: links
type(srd_t)          :: srd
type(merge_t)        :: merge
real(R8P)            :: u(48), v(48), x(17)
integer(I4P)         :: i, p, small

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]   ! a disc, 48 segments
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(0.0625_R8P * i, i = 0, 16)]                                               ! 16 x 16 cells
call disc%cut_cells(x=x, y=x, cut=cut, moments_order=2_I4P)
call cut_links_2d(cut, links)
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, x, srd)
call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, x, merge)
small = count(srd%status /= SRD_REGULAR .and. cut%piece_area > 0._R8P)
print '(A,I0,A,I0,A,I0)', 'pieces ', cut%np, ', links ', size(links%piece), ', small pieces ', small
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)                    ! the smallest piece
print '(A,I0,A,2I3,A,F7.4)', 'smallest piece ', p, ' in cell', cut%piece_cell(:, p), ', area fraction ', &
                             cut%piece_area(p) / 0.0625_R8P**2
do i = srd%start(p), srd%start(p + 1) - 1
   print '(A,I4,A,2I3)', '  SRD member ', srd%member(i), ' in cell', cut%piece_cell(:, srd%member(i))
enddo
print '(A,2F8.4)', '  wall normal ', cut%piece_wet_normal(:, p) / norm2(cut%piece_wet_normal(:, p))
print '(A,F7.4,A,F7.4)', '  weighted volume / cell ', srd%vhat(p) / 0.0625_R8P**2, ', eta_nb ', srd%eta_nb(p)
if (merge%status(p) == MERGE_OK) print '(A,I0,A,2I3)', '  merge target ', merge%target(p), ' in cell', &
                                                      cut%piece_cell(:, merge%target(p))
print '(A,3F9.5)', '  moments 1, X, Y   ', cut%piece_moments(1:3, p) / 0.0625_R8P**2
endprogram ck_cut_solver
