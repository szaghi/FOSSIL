!< Example of the reference page "Hole filling": three holes in the unit cube cube-01.stl (48 facets, each face a 2 x 2
!< grid of quads split in two). The face z = 0 is removed whole (one loop of 8 edges); on the face z = 1 two facets are
!< removed that touch only at the centre of the face: one walk through that vertex twice, split into two loops of 3.
!< The size threshold decides which loops are filled.
program ref_hole_filling
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
implicit none
type(surface_stl_object)        :: cube, holed
type(facet_object), allocatable :: kept(:), every(:)
type(facet_object), pointer     :: f
integer(I4P)                    :: i, m, filled, skipped, open_edges

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(kept(cube%get_facets_number()))
m = 0
do i = 1, cube%get_facets_number()
   f => cube%facet_at(i)
   if (f%normal%z < -0.5_R8P) cycle                                               ! the face z = 0, whole
   if (f%normal%z > 0.5_R8P .and. abs(f%centroid%x + f%centroid%y - 7._R8P / 6._R8P) < 1.e-9_R8P .and. &
       abs(abs(f%centroid%x - f%centroid%y) - 0.5_R8P) < 1.e-9_R8P) cycle       ! two facets of z = 1, a pinch
   m = m + 1 ; kept(m) = f
enddo
call move_alloc(kept, every)
allocate(kept, source=every(1:m))
call holed%adopt_facets(facets=kept)
print '(A,I0,A,L1)', 'cut cube: ', holed%get_facets_number(), ' facets, watertight ', holed%is_watertight()
call holed%fill_holes(holes_filled=filled, holes_skipped=skipped, open_edges=open_edges)
print '(A,I0,A,I0,A,I0)', 'default max_hole_size 0.1: filled ', filled, ', skipped ', skipped, ', open edges ', open_edges
call holed%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled, holes_skipped=skipped)
print '(A,I0,A,I0)', 'max_hole_size 0.5:         filled ', filled, ', skipped ', skipped
call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, holes_skipped=skipped)
print '(A,I0,A,I0)', 'max_hole_size 1:           filled ', filled, ', skipped ', skipped
print '(A,I0,A,L1,A,F8.5)', 'result: ', holed%get_facets_number(), ' facets, watertight ', holed%is_watertight(), &
                            ', volume ', holed%get_volume()
endprogram ref_hole_filling
