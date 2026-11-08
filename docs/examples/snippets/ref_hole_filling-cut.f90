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
