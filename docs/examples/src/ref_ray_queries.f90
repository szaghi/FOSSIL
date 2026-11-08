!run ref_ray_queries ref_ray_queries
!< Example of the reference page "Ray queries": the three queries on the unit cube cube-01.stl (48 facets, four per
!< face around its centre), then the cases at its edges: an origin on the surface, and a ray through a shared edge.
program ref_ray_queries
!region queries
use fossil,              only : surface_stl_object, ray_hit_t, RAY_STATUS_OK
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object)     :: cube
type(ray_hit_t)              :: hit
type(ray_hit_t), allocatable :: hits(:)
type(facet_object), pointer  :: f
type(vector_R8P)             :: origin
logical                      :: has_hit, blocked
integer(I4P)                 :: status, k, own

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
origin = vector_R8P(-1._R8P, 0.3_R8P, 0.6_R8P)
call cube%intersect_ray_first(ray_origin=origin, ray_direction=ex_R8P, hit=hit, has_hit=has_hit, status=status)
print '(A,L1,A,I0,A,F6.3,A,3F6.3)', 'first: ', has_hit, ', facet ', hit%facet_id, ', t ', hit%t, ', point ', &
                                    hit%point%x, hit%point%y, hit%point%z
call cube%intersect_ray_first(ray_origin=origin, ray_direction=2._R8P * ex_R8P, hit=hit, has_hit=has_hit)
print '(A,F6.3)', 'first, direction 2 ex: t ', hit%t                       ! parametric: half the distance
call cube%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
do k = 1, size(hits)
   f => cube%facet_at(hits(k)%facet_id)
   print '(A,I0,A,I0,A,F6.3,A,A)', 'all: hit ', k, ', facet ', hits(k)%facet_id, ', t ', hits(k)%t, ', ', &
                                   merge('enters', 'leaves', f%normal%dotproduct(ex_R8P) < 0._R8P)
enddo
do k = 1, 3                                                                   ! max_t = 0.5, 1, 1.5
   call cube%intersect_ray_any(ray_origin=origin, ray_direction=ex_R8P, max_t=0.5_R8P * k, found=blocked)
   print '(A,F4.1,A,L1)', 'any: within t <= ', 0.5_R8P * k, ': ', blocked
enddo
!endregion queries
!region edges
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
!endregion edges
endprogram ref_ray_queries
