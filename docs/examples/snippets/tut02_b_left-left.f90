use fossil,                    only : surface_stl_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_facet_object,       only : facet_object
use penf,                      only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface
type(facet_object), pointer :: f
integer(I4P), allocatable   :: edges(:,:), kinds(:), size_of(:)
integer(I4P)                :: i, n_open, u, main
character(96)               :: caption

call surface%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
call surface%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
print '(A,I0,A,I0)', 'boundary edges ', count(kinds == FEATURE_EDGE_BOUNDARY), &
                     ', non-manifold edges ', count(kinds == FEATURE_EDGE_NONMANIFOLD)
n_open = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i) ; n_open = n_open + count(f%fcon_edge == 0_I4P)
enddo
print '(A,I0)', 'facet edges without a single neighbour ', n_open
print '(A,I0,A,I0)', 'shells ', surface%get_shells_number(), ', bodies ', surface%get_bodies_number()
