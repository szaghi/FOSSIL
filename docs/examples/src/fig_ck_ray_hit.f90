!run figck-ray_hit fig_ck_ray_hit
!fig ck-ray_hit
!< Figure of the cookbook recipe "Cast a ray": the dragon, the ray of the recipe from its origin to the first hit
!< intersect_ray_first returns, and the facet hit (ember). The ray is drawn as 64 pieces so the painter sort hides
!< the parts behind the body.
program fig_ck_ray_hit
use fossil,     only : surface_stl_object, ray_hit_t
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: dragon
type(ray_hit_t)          :: hit
type(vector_R8P)         :: origin, a(64), b(64)
logical                  :: has_hit
integer(I4P)             :: u, k
character(80)            :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
origin = vector_R8P(0._R8P, -20._R8P, 12._R8P)
call dragon%intersect_ray_first(ray_origin=origin, ray_direction=vector_R8P(0._R8P, 1._R8P, 0._R8P), hit=hit, &
                                has_hit=has_hit)
if (.not. has_hit) error stop 'fig_ck_ray_hit: no hit'
write(caption, '(A,I0,A,F7.4)') 'intersect_ray_first: facet ', hit%facet_id, ' at t = ', hit%t
call fig_open(u, 'ck-ray_hit', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_surface(u, 'patch', dragon, first=hit%facet_id, last=hit%facet_id)
do k = 1, 64
   a(k) = origin + ((k - 1) / 64._R8P) * (hit%point - origin) ; b(k) = origin + (k / 64._R8P) * (hit%point - origin)
enddo
call fig_segments(u, 'probe', a, b)
call fig_points(u, 'probe', [origin])
call fig_points(u, 'foot', [hit%point])
call fig_label(u, 'probe', origin, 'origin')
call fig_close(u)
endprogram fig_ck_ray_hit
