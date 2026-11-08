!run figref-cut-cells-checks fig_ref_cut_cells_checks
!fig ref-cut-cells-checks
!< Figure of the reference page "Cut-cell geometry", section "Checks you can make": cube-01 turned by 0.3 rad about
!< (1, 1, 1) and shifted, on a 6^3 block of cells of 0.25 (the case of ref_cut_cells). The covering sum and the
!< per-cell divergence residual are printed; the cells whose fluid fraction is below 1e-2, the small cells a solver
!< must treat, are outlined in red over the opaque cube, the smallest labelled: almost all solid, only a corner of each
!< sticks out of a face of the cube into the fluid.
program fig_ref_cut_cells_checks
use fossil,     only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_box, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter     :: H = 0.25_R8P
type(surface_stl_object) :: cube
type(cut_cells_3d_t)     :: cut
real(R8P)                :: g(7), res
integer(I4P)             :: i, j, k, w, small(3)
character(32)            :: lab

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
g = [(-0.1_R8P + H * i, i = 0, 6)]
call cut_cells_3d(surface=cube, x=g, y=g, z=g, cut=cut)
res = 0._R8P
do k = 1, 6 ; do j = 1, 6 ; do i = 1, 6
   res = max(res, maxval(abs([cut%face_x(i + 1, j, k) - cut%face_x(i, j, k), cut%face_y(i, j + 1, k) - cut%face_y(i, j, k), &
                              cut%face_z(i, j, k + 1) - cut%face_z(i, j, k)] * H**2 - cut%wet_normal(:, i, j, k))))
enddo ; enddo ; enddo
print '(A,ES9.2,A,ES9.2)', 'covering sum - body volume ', sum(1._R8P - cut%alpha) * H**3 - cube%get_volume(), &
                           ', largest divergence residual ', res
small = minloc(cut%alpha, cut%alpha > 0._R8P)
print '(A,I0,A,3I2,A,ES10.3)', 'small cells (alpha < 1e-2) ', count(cut%alpha > 0._R8P .and. cut%alpha < 1.e-2_R8P), &
                               ', smallest in cell', small, ', alpha ', cut%alpha(small(1), small(2), small(3))

call fig_open(w, 'ref-cut-cells-checks', 'cells with fluid fraction below 1e-2 (red): only a corner sticks out of '// &
              'the cube', -55._R8P, 25._R8P, width=720_I4P, height=520_I4P)
call fig_zoom(w, vector_R8P(g(small(1)), g(small(2)), g(small(3))), 0.75_R8P)
call fig_surface(w, 'body', cube)
do k = 1, 6 ; do j = 1, 6 ; do i = 1, 6
   if (cut%alpha(i, j, k) <= 0._R8P .or. cut%alpha(i, j, k) >= 1.e-2_R8P) cycle
   call fig_box(w, 'hole', vector_R8P(g(i), g(j), g(k)), vector_R8P(g(i + 1), g(j + 1), g(k + 1)))
enddo ; enddo ; enddo
write(lab, '(A,ES9.2)') 'alpha ', cut%alpha(small(1), small(2), small(3))
call fig_label(w, 'accent', vector_R8P(g(small(1) + 1), g(small(2) + 1), g(small(3) + 1)), trim(adjustl(lab)))
call fig_close(w)
endprogram fig_ref_cut_cells_checks
