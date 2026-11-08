origin = vector_R8P(0._R8P, 0.3_R8P, 0.6_R8P)                                ! on the face x = 0
call cube%intersect_ray_first(ray_origin=origin, ray_direction=ex_R8P, hit=hit, has_hit=has_hit)
print '(A,I0,A,F6.3)', 'origin on the surface: first facet ', hit%facet_id, ', t ', hit%t
own = hit%facet_id
call cube%intersect_ray_first(ray_origin=origin, ray_direction=ex_R8P, hit=hit, has_hit=has_hit, exclude=own)
print '(A,I0,A,F6.3)', '  the same, that facet excluded: facet ', hit%facet_id, ', t ', hit%t
origin = vector_R8P(-1._R8P, 0.3_R8P, 0.3_R8P)                               ! through the edge of two facets
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A,I0,A)', 'through a shared edge: ', size(hits), ' hits, ', count(hits(2:)%t /= hits(:size(hits) - 1)%t) + 1, &
                       ' distinct t'
print '(A,8I4)', '  facets', hits(:)%facet_id
