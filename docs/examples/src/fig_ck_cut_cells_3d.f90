!run figck-cut-cells-3d fig_ck_cut_cells_3d
!fig ck-cut-cells-3d
!< Figure of the cookbook recipe "Cut a 3-D grid by a body": the turned unit cube of ck_cut_cells_3d (translucent) and
!< one layer of its 16^3 cells, 0.5 <= z <= 0.625, each drawn at the layer's mid-plane and coloured by its fluid
!< fraction alpha (cream fluid, sepia solid; the cut cells in between).
program fig_ck_cut_cells_3d
use fossil,     only : surface_stl_object, cut_cells_3d_t, cut_cells_3d
use fossil_fig, only : fig_open, fig_surface, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : ex_R8P, ez_R8P, vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
integer(I4P), parameter  :: K = 9
type(surface_stl_object) :: cube
type(cut_cells_3d_t)     :: cut
real(R8P)                :: x(17), z
integer(I4P)             :: i, j, u
character(90)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%rotate(axis=ez_R8P, angle=PI / 6._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%rotate(axis=ex_R8P, angle=PI / 9._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube%analyze
x = [(-0.5_R8P + 0.125_R8P * i, i = 0, 16)]
call cut_cells_3d(cube, x, x, x, cut)
write(caption, '(A,I0,A)') 'fluid fraction alpha of the layer k = ', K, ' (0.5 < z < 0.625): 1 cream, 0 sepia'
call fig_open(u, 'ck-cut-cells-3d', trim(caption), -60._R8P, 35._R8P)
z = 0.5_R8P * (x(K) + x(K + 1))
do j = 1, 16 ; do i = 1, 16
   call fig_quad(u, 'frac', vector_R8P(x(i), x(j), z), vector_R8P(x(i + 1), x(j), z), vector_R8P(x(i + 1), x(j + 1), z), &
                 vector_R8P(x(i), x(j + 1), z), value=1._R8P - cut%alpha(i, j, K))
enddo ; enddo
call fig_surface(u, 'ghost', cube)
call fig_close(u)
endprogram fig_ck_cut_cells_3d
