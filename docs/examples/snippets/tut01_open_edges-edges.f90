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
