!run tut01-open tut01_open_edges
!fig tut01-open
!< Tutorial 1, second run: why the dragon is not a closed surface. One facet and its neighbours, the boundary edges
!< (an edge with one facet) and the non-manifold edges (three or more), the shells and bodies analyze found.
program tut01_open_edges
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_fig,                only : fig_open, fig_surface, fig_triangle, fig_segments, fig_close
!region facet
use fossil,                    only : surface_stl_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_facet_object,       only : facet_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
type(surface_stl_object)          :: surface
type(facet_object), pointer       :: f
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:)
integer(I4P), allocatable         :: edges(:,:), kinds(:)
integer(I4P)                      :: i, v, open_edges, u
character(96)                     :: caption

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
f => surface%facet_at(1_I4P)                                            ! a pointer into the surface, read-only use
do v = 1, 3
   print '(A,I0,A,3F9.4)', 'facet 1 vertex ', v, '   ', f%vertex(v)%x, f%vertex(v)%y, f%vertex(v)%z
enddo
print '(A,3F9.4)', 'facet 1 normal     ', f%normal%x, f%normal%y, f%normal%z
print '(A,3I6)',   'facet 1 neighbours ', f%fcon_edge                    ! across edges 1-2, 2-3, 3-1; 0: none
!endregion facet
!region edges
call surface%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)    ! no sharp edge at 180 degrees
print '(A,I0)', 'boundary edges     ', count(kinds == FEATURE_EDGE_BOUNDARY)      ! one facet
print '(A,I0)', 'non-manifold edges ', count(kinds == FEATURE_EDGE_NONMANIFOLD)   ! three or more
open_edges = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   if (any(f%fcon_edge == 0_I4P)) open_edges = open_edges + 1             ! no single neighbour across an edge
enddo
print '(A,I0)', 'facets on them     ', open_edges
print '(A,I0,A,I0)', 'shells             ', surface%get_shells_number(), ', bodies ', surface%get_bodies_number()
!endregion edges
pool => surface%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
write(caption, '(I0,A,I0,A)') count(kinds == FEATURE_EDGE_NONMANIFOLD), ' non-manifold edges: the ', &
                              open_edges, ' facets on them (ember)'
call fig_open(u, 'tut01-open', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', surface)
do i = 1, surface%get_facets_number()                                   ! the facets on those edges
   f => surface%facet_at(i)
   if (any(f%fcon_edge == 0_I4P)) call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_segments(u, 'accent', pack(a, kinds == FEATURE_EDGE_NONMANIFOLD), pack(b, kinds == FEATURE_EDGE_NONMANIFOLD))
call fig_close(u)
endprogram tut01_open_edges
