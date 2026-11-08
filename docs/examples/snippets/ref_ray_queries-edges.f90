origin = vector_R8P(0._R8P, 0.3_R8P, 0.6_R8P)                                ! on the face x = 0
call cube%intersect_ray_first(ray_origin=origin, ray_direction=ex_R8P, hit=hit, has_hit=has_hit)
print '(A,I0,A,F6.3)', 'origin on the surface: first facet ', hit%facet_id, ', t ', hit%t
own = hit%facet_id
call cube%intersect_ray_first(ray_origin=origin, ray_direction=ex_R8P, hit=hit, has_hit=has_hit, exclude=own)
print '(A,I0,A,F6.3)', '  the same, that facet excluded: facet ', hit%facet_id, ', t ', hit%t
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A,2F6.3)', '  all: ', size(hits), ' hits, t', hits(:)%t
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits, exclude=own)
print '(A,I0,A,F6.3)', '  all, that facet excluded: ', size(hits), ' hit, t', hits(1)%t
origin = vector_R8P(-1._R8P, 0.3_R8P, 0.3_R8P)                               ! through the edge of two facets
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A,2I4)', 'through a shared edge: ', size(hits), ' hits, facets', hits(:)%facet_id
origin = vector_R8P(-1._R8P, 0.5_R8P, 0.5_R8P)                               ! through the vertex of four facets
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A,2F6.3)', 'through the face centres: ', size(hits), ' hits, t', hits(:)%t
