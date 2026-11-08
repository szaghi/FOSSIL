use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, cut_links_t, &
                   cut_links_2d, srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, merge_t, &
                   merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)       :: section
type(surface_stl_object) :: surface
type(cut_cells_2d_t)     :: cut
type(cut_cells_3d_t)     :: cut3
type(cut_links_t)        :: links, links3
type(srd_t)              :: srd
type(merge_t)            :: mt
real(R8P)                :: x(49), h
integer(I4P)             :: i, k, p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]
call section%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)                        ! piece p meets links%piece(k) over face links%face(k)
print '(A,I0,A,I0)', 'pieces ', cut%np, ', links (both directions) ', size(links%piece)
p = cut%first_piece(22, 20) + 1                      ! the small second piece of cell (22, 20), step 2
do k = links%start(p), links%start(p + 1) - 1
   print '(A,I0,A,I0,A,2I3,A,I0,A,F7.4)', '  piece ', p, ' -> ', links%piece(k), ' in cell', &
         cut%piece_cell(:, links%piece(k)), ' over face ', links%face(k), ', length / h ', links%area(k) / h
enddo
