program ref_cut_solver
use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, cut_links_t, &
                   cut_links_2d, srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, frd_t, frd_neighbourhoods, &
                   merge_t, merge_targets, MERGE_OK, cut_coarsen_2d, mono_count
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: cube
type(curve_object)        :: half
type(cut_cells_3d_t)      :: cut3
type(cut_cells_2d_t)      :: cut, coarse
type(cut_links_t)         :: links3, links, coarse_links
type(srd_t)               :: srd
type(frd_t)               :: frd
type(merge_t)             :: mt
real(R8P)                 :: x(9), h(11)
integer(I4P), allocatable :: fine_to_coarse(:)
integer(I4P)              :: i, k, p, status

x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)           ! [0, 1]^3 on the cell faces
call cut_cells_3d(cube, x, x, x, cut3, links=links3, edges=.true.)             ! 3-D: links with the cut cells
print '(A,I0,A,I0,A,I0,A,I0)', '3-D: pieces ', cut3%np, ' (', count(cut3%piece_volume == 0._R8P), &
                               ' of zero volume), links ', size(links3%piece), ', unmatched regions ', links3%unmatched
p = 1
do k = links3%start(p), links3%start(p + 1) - 1 ! piece p meets piece(k) over its face face(k), area area(k)
   print '(A,I0,A,I4,A,I0,A,F7.4)', '  piece ', p, ' -> ', links3%piece(k), ' over face ', links3%face(k), &
                                    ', area ', links3%area(k)
enddo
print '(A,I0,A,I0)', '  z edges with fluid ', count(cut3%edge_z > 0._R8P), ' of ', size(cut3%edge_z)

h = [(0.1_R8P * i, i = 0, 10)] ! 10 x 10 cells of 0.1; the solid is x < 0.58: fluid fraction 0.2 in the cut column
call half%adopt_loops(u=[-1._R8P, 0.58_R8P, 0.58_R8P, -1._R8P], v=[-1._R8P, -1._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
call half%cut_cells(x=h, y=h, cut=cut, moments_order=2_I4P)
call cut_links_2d(cut, links)                                                  ! 2-D: from the cut cells
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, h, h, srd)
print '(A,I0,A,I0,A,I0)', '2-D: pieces ', cut%np, ', small ', count(srd%status == SRD_COMPLETE), ', regular ', &
                          count(srd%status == SRD_REGULAR)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P .and. cut%piece_cell(2, :) == 5)
print '(A,2I3,A,F5.2,A,10I4)', '  piece in cell', cut%piece_cell(:, p), ', fraction ', cut%piece_area(p) / 0.01_R8P, &
                               ', SRD members', srd%member(srd%start(p):srd%start(p + 1) - 1)
print '(A,2I3,A,F6.3,A,F6.3)', '  member cells', cut%piece_cell(:, srd%member(srd%start(p) + 1)), &
                               ', eta_nb ', srd%eta_nb(p), ', eta_self ', srd%eta_self(p)

call frd_neighbourhoods(links, cut%piece_cell, cut%piece_area, h, h, frd)
print '(A,I0,A,F6.3)', '  FRD set size ', frd%start(p + 1) - frd%start(p), ', own weight ', frd%weight(frd%start(p))

call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, h, h, mt)
if (mt%status(p) == MERGE_OK) print '(A,2I3,A,I0,A,F5.2)', '  merge target in cell', cut%piece_cell(:, mt%target(p)), &
   ', normal-directed targets ', mt%mix_start(p + 1) - mt%mix_start(p), ', weight ', mt%mix_weight(mt%mix_start(p))

print '(A,I0,A,6F9.5)', '  moments up to order 2 (', mono_count(2_I4P, 2_I4P), '), / cell area: ', &
                        cut%piece_moments(:, p) / 0.01_R8P

call cut_coarsen_2d(cut, links, h, h, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
print '(A,I0,A,I0,A,I0,A,F5.2)', 'coarse (ratio 2): status ', status, ', pieces ', coarse%np, ', links ', &
                                 size(coarse_links%piece), ', fraction of the cut column ', coarse%alpha(3, 3)
print '(A,I0,A,I0)', '  the fine piece lands in coarse piece ', fine_to_coarse(p), ', in coarse cell ', &
                     coarse%piece_cell(1, fine_to_coarse(p))
endprogram ref_cut_solver
