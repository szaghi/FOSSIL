!< Figure of the cookbook recipe "Find the sharp edges of a part": the dragon with its sharp edges above 60 degrees
!< (gold) and its non-manifold edges (ember), from feature_edges, zoomed on the head.
program fig_ck_feature_edges
use fossil,                    only : surface_stl_object, FEATURE_EDGE_SHARP, FEATURE_EDGE_NONMANIFOLD
use fossil_fig,                only : fig_open, fig_zoom, fig_surface, fig_segments, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P
implicit none
type(surface_stl_object)          :: dragon
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:)
integer(I4P), allocatable         :: edges(:,:), kinds(:)
integer(I4P)                      :: u, i
character(96)                     :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%feature_edges(edges=edges, kinds=kinds, angle=60._R8P, include_boundary=.false.)
pool => dragon%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
write(caption, '(I0,A,I0,A)') count(kinds == FEATURE_EDGE_SHARP), ' sharp edges above 60 degrees (gold), ', &
                              count(kinds == FEATURE_EDGE_NONMANIFOLD), ' non-manifold (ember)'
call fig_open(u, 'ck-feature-edges', trim(caption), -60._R8P, 22._R8P)
call fig_zoom(u, vector_R8P(-6._R8P, 0._R8P, 16._R8P), 5._R8P)                ! the head
call fig_surface(u, 'body', dragon)
call fig_segments(u, 'accent', pack(a, kinds == FEATURE_EDGE_SHARP), pack(b, kinds == FEATURE_EDGE_SHARP))
call fig_segments(u, 'hole', pack(a, kinds == FEATURE_EDGE_NONMANIFOLD), pack(b, kinds == FEATURE_EDGE_NONMANIFOLD))
call fig_close(u)
endprogram fig_ck_feature_edges
