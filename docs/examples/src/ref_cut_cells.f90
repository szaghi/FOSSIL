!run ref_cut_cells ref_cut_cells
!< Reference page "Cut-cell geometry": the checks of the page on small exact cases. 2-D: a solid strip thinner than
!< the cell across a row of cells (two fluid pieces per cell). 3-D: the unit cube turned about (1, 1, 1) and shifted,
!< on a 6^3 block: covering sum, per-cell divergence residual, small cells, wall fractions along the centre lines.
program ref_cut_cells
!region strip
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
!endregion strip
!region cube
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
h = 0.25_R8P
g = [(-0.1_R8P + h * i, i = 0, 6)]                                     ! 6^3 cells, the block covers the body
call cut_cells_3d(surface=cube, x=g, y=g, z=g, cut=cut3)
print '(A,I0,A,I0,A,I0)', 'cube: cut cells ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1._R8P), &
                          ', multi-piece ', count(cut3%pieces > 1), ', flagged ', count(cut3%pieces < 0)
print '(A,ES9.2)', '  covering sum: solid volume - body volume ', sum(1._R8P - cut3%alpha) * h**3 - cube%get_volume()
res = 0._R8P
do k = 1, 6 ; do j = 1, 6 ; do i = 1, 6                                ! divergence theorem, cell by cell
   res = max(res, maxval(abs([cut3%face_x(i + 1, j, k) - cut3%face_x(i, j, k), &
                              cut3%face_y(i, j + 1, k) - cut3%face_y(i, j, k), &
                              cut3%face_z(i, j, k + 1) - cut3%face_z(i, j, k)] * h**2 - cut3%wet_normal(:, i, j, k))))
enddo ; enddo ; enddo
print '(A,ES9.2)', '  largest divergence residual ', res
print '(A,I0,A,ES10.3)', '  small cells (0 < alpha < 1e-2) ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1.e-2_R8P), &
                         ', smallest alpha ', minval(cut3%alpha, cut3%alpha > 0._R8P)
!endregion cube
!region fractions
c = 0.5_R8P * (g(1:6) + g(2:7))                                       ! the cell centres
call cube%cell_wall_fractions(x=c, y=c, z=c, frac=frac, inside=inside)
print '(A,I0,A,I0)', 'centres inside ', count(inside), ' of ', size(inside)
print '(A,L1,A,6F7.3)', '  centre (2, 3, 3): inside ', inside(2, 3, 3), ', fractions (-x .. +z)', frac(:, 2, 3, 3)
!endregion fractions
endprogram ref_cut_cells
