!run figck-cut-coarsen-3d fig_ck_cut_coarsen_3d
!fig ck-cut-coarsen-3d-fine
!fig ck-cut-coarsen-3d-coarse
!< Figures of the cookbook recipe "Coarsen 3-D cut cells": the fluid fraction of the turned cube's cut cells, seen from
!< above, on the fine 16^3 grid (layer 0.5 < z < 0.625) and on the 8^3 grid cut_coarsen_3d builds from it (layer
!< 0.5 < z < 0.75); cream fluid, sepia solid.
program fig_ck_cut_coarsen_3d
use fossil,     only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_links_t, cut_coarsen_3d
use fossil_fig, only : fig_open, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : ex_R8P, ez_R8P, vector_R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: cube
type(cut_cells_3d_t)      :: fine, coarse
type(cut_links_t)         :: fine_links, coarse_links
integer(I4P), allocatable :: fine_to_coarse(:)
real(R8P)                 :: xf(17)
integer(I4P)              :: i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%rotate(axis=ez_R8P, angle=PI / 6._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%rotate(axis=ex_R8P, angle=PI / 9._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%analyze
xf = [(-0.5_R8P + 0.125_R8P * i, i = 0, 16)]
call cut_cells_3d(cube, xf, xf, xf, fine, links=fine_links)
call cut_coarsen_3d(fine, fine_links, xf, xf, xf, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
call layer('ck-cut-coarsen-3d-fine', 'fine, 16^3 cells: layer 0.5 < z < 0.625', fine, xf, 9_I4P)
call layer('ck-cut-coarsen-3d-coarse', 'cut_coarsen_3d, 8^3: layer 0.5 < z < 0.75', coarse, &
           xf(1::2), 5_I4P)
contains
   subroutine layer(name, caption, cut, x, k)
   !< One layer of cells, from above, coloured by solid fraction.
   character(*),         intent(in) :: name, caption
   type(cut_cells_3d_t), intent(in) :: cut
   real(R8P),            intent(in) :: x(:)
   integer(I4P),         intent(in) :: k
   integer(I4P)                     :: i, j, u

   call fig_open(u, name, caption, -90._R8P, 90._R8P, width=420_I4P, height=420_I4P)
   do j = 1, cut%ny ; do i = 1, cut%nx
      call fig_quad(u, 'frac', vector_R8P(x(i), x(j), 0._R8P), vector_R8P(x(i + 1), x(j), 0._R8P), &
                    vector_R8P(x(i + 1), x(j + 1), 0._R8P), vector_R8P(x(i), x(j + 1), 0._R8P), &
                    value=1._R8P - cut%alpha(i, j, k))
   enddo ; enddo
   call fig_close(u)
   endsubroutine layer
endprogram fig_ck_cut_coarsen_3d
