program ck_cut_coarsen
use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_coarsen_2d
use penf,   only : I4P, R8P
implicit none
type(curve_object)        :: disc
type(cut_cells_2d_t)      :: fine, coarse, direct
type(cut_links_t)         :: links, coarse_links, direct_links
real(R8P)                 :: u(48), v(48), xf(33), xc(17)
integer(I4P), allocatable :: fine_to_coarse(:)
integer(I4P)              :: i, status

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]   ! a disc, 48 segments
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
xf = [(i / 32._R8P, i = 0, 32)] ; xc = [(i / 16._R8P, i = 0, 16)]
call disc%cut_cells(x=xf, y=xf, cut=fine)                                       ! 32 x 32
call cut_links_2d(fine, links)
call cut_coarsen_2d(fine, links, xf, xf, 2_I4P, coarse, coarse_links, fine_to_coarse, status)   ! to 16 x 16
call disc%cut_cells(x=xc, y=xc, cut=direct)                                     ! 16 x 16 directly
call cut_links_2d(direct, direct_links)
print '(A,I0,A,I0,A,I0)', 'fine pieces ', fine%np, ' -> coarse ', coarse%np, ', direct ', direct%np
print '(A,ES9.2)', 'max |alpha coarse - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha))
print '(A,ES9.2)', 'max |aperture difference|         ', maxval(abs(coarse%aperture - direct%aperture))
print '(A,I0,A,I0)', 'links: coarsened ', size(coarse_links%piece), ', direct ', size(direct_links%piece)
endprogram ck_cut_coarsen
