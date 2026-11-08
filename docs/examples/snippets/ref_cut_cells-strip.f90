use fossil, only : curve_object, cut_cells_2d_t, surface_stl_object, cut_cells_3d_t, cut_cells_3d, rigid_pose_t, &
                   pose_from_axis_angle
use penf,   only : I4P, R8P
implicit none
type(curve_object)       :: strip
type(cut_cells_2d_t)     :: cut2
type(surface_stl_object) :: cube
type(cut_cells_3d_t)     :: cut3
real(R8P)                :: e(5), g(7), c(6), h, res, frac(6, 6, 6, 6)
logical                  :: inside(6, 6, 6)
integer(I4P)             :: i, j, k, p

e = [(0.25_R8P * i, i = 0, 4)]                                         ! 4 x 4 cells of side 1/4
call strip%adopt_loops(u=[-0.1_R8P, 1.1_R8P, 1.1_R8P, -0.1_R8P], &      ! a solid strip 0.05 thick across row 2
                       v=[0.35_R8P, 0.35_R8P, 0.40_R8P, 0.40_R8P], loop_sizes=[4_I4P])
call strip%cut_cells(x=e, y=e, cut=cut2)
print '(A,4I2)', 'strip: pieces of row 2 ', cut2%pieces(:, 2)
p = cut2%first_piece(2, 2)
print '(A,2F8.4,A,F8.4)', '  cell (2, 2): piece areas / cell ', cut2%piece_area(p:p + 1) / 0.0625_R8P, &
                          ', alpha ', cut2%alpha(2, 2)
print '(A,4F6.2,A,4F6.2)', '  apertures (-u, +u, -v, +v) ', cut2%piece_aperture(:, p), ' and ', cut2%piece_aperture(:, p + 1)
