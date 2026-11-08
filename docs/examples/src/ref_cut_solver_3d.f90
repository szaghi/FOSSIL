!run ref_cut_solver_3d ref_cut_solver_3d
!< Reference page "Cut cells: solver support", the 3-D calls: cube-01 turned by 0.3 rad about (1, 1, 1) and shifted, on
!< 8^3 cells of 0.2. One cut_cells_3d call with links, edge fractions, moments to order 2 and quadrature rules; the
!< SRD neighbourhoods and merge targets in 3-D; then cut_coarsen_3d at ratio 2, compared with the cut cells computed
!< directly on the coarse grid.
program ref_cut_solver_3d
use fossil, only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_links_t, srd_t, srd_neighbourhoods, &
                   SRD_COMPLETE, SRD_INCOMPLETE, merge_t, merge_targets, MERGE_NONE, MERGE_OK, cut_coarsen_3d, &
                   pose_from_axis_angle
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: cube
type(cut_cells_3d_t)      :: cut, coarse, direct
type(cut_links_t)         :: links, coarse_links
type(srd_t)               :: srd
type(merge_t)             :: mt
real(R8P)                 :: g(9), gc(5), q
integer(I4P), allocatable :: f2c(:)
integer(I4P)              :: p, status

!region cut
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
g = [(-0.15_R8P + 0.2_R8P * p, p = 0, 8)]
call cut_cells_3d(cube, g, g, g, cut, status=status, links=links, edges=.true., moments_order=2_I4P, quadrature=.true.)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'status ', status, ', pieces ', cut%np, ', multi-piece cells ', &
                                    count(cut%pieces > 1), ', links ', size(links%piece), ', unmatched ', links%unmatched
print '(A,I0,A,I0)', 'cut edges (x, y, z) ', count(cut%edge_x > 0._R8P .and. cut%edge_x < 1._R8P) + &
                     count(cut%edge_y > 0._R8P .and. cut%edge_y < 1._R8P) + &
                     count(cut%edge_z > 0._R8P .and. cut%edge_z < 1._R8P), ', nodes inside ', count(cut%node_in)
q = 0._R8P
do p = 1, cut%np                                       ! the zeroth moment is the volume
   q = max(q, abs(cut%piece_moments(1, p) - cut%piece_volume(p)))
enddo
print '(A,I0,A,ES9.2,A,I0)', 'moments per piece ', size(cut%piece_moments, 1), ', |m0 - volume| ', q, &
                             ', pieces with quad_status 1 ', count(cut%quad_status /= 0)
!endregion cut
!region graph
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_volume, cut%piece_centroid, cut%piece_wet_normal, g, g, srd, z=g)
call merge_targets(links, cut%piece_cell, cut%piece_volume, cut%piece_wet_normal, g, g, mt, z=g)
print '(A,I0,A,I0,A,ES9.2)', 'SRD: small ', count(srd%status == SRD_COMPLETE), ', incomplete ', &
                             count(srd%status == SRD_INCOMPLETE), ', sum vhat - fluid volume ', &
                             sum(srd%vhat) - sum(cut%piece_volume)
print '(A,I0,A,I0)', 'merge: served ', count(mt%status /= MERGE_NONE), ', direct ', count(mt%status == MERGE_OK)
!endregion graph
!region coarsen
call cut_coarsen_3d(cut, links, g, g, g, 2_I4P, coarse, coarse_links, f2c, status)
gc = g(1:9:2)
call cut_cells_3d(cube, gc, gc, gc, direct)
print '(A,I0,A,I0,A,I0)', 'coarse: status ', status, ', pieces ', coarse%np, ', direct ', direct%np
print '(A,ES9.2,A,L1)', '  largest |alpha coarsened - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha)), &
                        ', same piece counts ', all(coarse%pieces == direct%pieces)
!endregion coarsen
endprogram ref_cut_solver_3d
