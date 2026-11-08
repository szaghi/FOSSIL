program ck_cut_coarsen_3d
use fossil, only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_links_t, cut_coarsen_3d
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, ez_R8P, vector_R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: cube
type(cut_cells_3d_t)      :: fine, coarse, direct
type(cut_links_t)         :: fine_links, coarse_links
integer(I4P), allocatable :: fine_to_coarse(:)
real(R8P)                 :: xf(17), xc(9)
integer(I4P)              :: i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%rotate(axis=ez_R8P, angle=PI / 6._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%rotate(axis=ex_R8P, angle=PI / 9._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%analyze
xf = [(-0.5_R8P + 0.125_R8P * i, i = 0, 16)]                  ! fine: 16^3 cells
call cut_cells_3d(cube, xf, xf, xf, fine, links=fine_links)
call cut_coarsen_3d(fine, fine_links, xf, xf, xf, 2_I4P, coarse, coarse_links, fine_to_coarse, status)  ! to 8^3
if (status /= 0) error stop 'cut_coarsen_3d failed'
xc = xf(1::2)                                                 ! check: cut the coarse grid directly
call cut_cells_3d(cube, xc, xc, xc, direct)
print '(A,I0,A,I0)', 'pieces: fine ', fine%np, ', coarse ', coarse%np
print '(A,ES9.2)', 'largest |alpha coarsened - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha))
print '(A,ES9.2)', 'largest |wet area coarsened - direct|    ', maxval(abs(coarse%wet_area - direct%wet_area))
print '(A,I0,A,I0)', 'coarse cells cut ', count(coarse%alpha > 0._R8P .and. coarse%alpha < 1._R8P), &
                     ', fine pieces mapped ', count(fine_to_coarse > 0)
endprogram ck_cut_coarsen_3d
