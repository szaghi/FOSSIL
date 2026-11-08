!< Figures of the reference page "Ray queries", What it does and Pipeline, on dragon.stl:
!< - picking: 7 x 5 rays from an eye point in front of the dragon, aimed at a grid over its box; each ray that hits is
!<   drawn up to its first hit (intersect_ray_first, gold) and the facet it picked is drawn in ember; the rays that miss
!<   are drawn faint, past the dragon;
!< - the per-facet test of one of those hits (even weights, a facet facing the eye), zoomed: the picked facet split at
!<   the hit point into three sub-triangles, each shaded by twice the barycentric weight of the vertex opposite it (its area
!<   fraction), (1 - u - v, u, v) from facet%intersect_ray (Moller-Trumbore), labelled at those vertices; the ray.
program fig_ref_ray_queries_picking
use fossil,     only : surface_stl_object, facet_object, ray_hit_t
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segment, fig_points, fig_label, fig_triangle, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter     :: NX = 7, NZ = 5
type(surface_stl_object)    :: dragon
type(facet_object), pointer :: fp
type(ray_hit_t)             :: hit, hits(NX * NZ)
type(vector_R8P)            :: eye, bmin, bmax, target(NX * NZ), far(NX * NZ), dir, p, lift
real(R8P)                   :: t, uu, vv, best, score, az, el
logical                     :: has_hit(NX * NZ), ok
integer(I4P)                :: i, j, n, u, pick
character(90)               :: caption
character(32)               :: text

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
eye = vector_R8P(bmin%x - 0.4_R8P * (bmax%x - bmin%x), bmin%y - 1.6_R8P * (bmax%y - bmin%y), &
                 0.5_R8P * (bmin%z + bmax%z))
n = 0
do j = 1, NZ ; do i = 1, NX
   n = n + 1
   target(n) = vector_R8P(bmin%x + (bmax%x - bmin%x) * (i - 0.5_R8P) / NX, 0.5_R8P * (bmin%y + bmax%y), &
                          bmin%z + (bmax%z - bmin%z) * (j - 0.5_R8P) / NZ)
   dir = target(n) - eye
   call dragon%intersect_ray_first(ray_origin=eye, ray_direction=dir, hit=hits(n), has_hit=has_hit(n))
   far(n) = eye + 1.6_R8P * dir
enddo ; enddo
print '(A,I0,A,I0,A)', 'intersect_ray_first: ', NX * NZ, ' rays, ', count(has_hit), ' hits'
write(caption, '(I0,A,I0,A)') NX * NZ, ' rays, ', count(has_hit), ' first hits (gold) and the facets picked (ember)'
call fig_open(u, 'ref-ray-queries-picking', trim(caption), -130._R8P, 25._R8P)
call fig_zoom(u, 0.5_R8P * (bmin + bmax), 0.45_R8P * dist(bmin, bmax))
call fig_surface(u, 'body', dragon)
do n = 1, NX * NZ
   if (has_hit(n)) then
      call fig_surface(u, 'patch', dragon, first=hits(n)%facet_id, last=hits(n)%facet_id)
      call fig_segment(u, 'accent', eye, hits(n)%point)
   else
      call fig_segment(u, 'box', eye, far(n))
   endif
enddo
call fig_points(u, 'foot', pack(hits(:)%point, has_hit))
call fig_points(u, 'probe', [eye])
call fig_close(u)

! ---- the per-facet test of a hit with even barycentric weights on a facet facing the eye
best = -1._R8P ; pick = 0
do n = 1, NX * NZ
   if (.not. has_hit(n)) cycle
   fp => dragon%facet_at(hits(n)%facet_id)
   call fp%intersect_ray(ray_origin=eye, ray_direction=target(n) - eye, t=t, u=uu, v=vv, hit=ok)
   dir = target(n) - eye
   score = min(1._R8P - uu - vv, uu, vv) * max(0._R8P, -(fp%normal%x * dir%x + fp%normal%y * dir%y + fp%normal%z * dir%z)) &
           / dist(target(n), eye)                                   ! even weights, facing the eye
   if (score > best) then
      best = score ; pick = n
   endif
enddo
hit = hits(pick)
fp => dragon%facet_at(hit%facet_id)
dir = target(pick) - eye
call fp%intersect_ray(ray_origin=eye, ray_direction=dir, t=t, u=uu, v=vv, hit=ok)
p = (1._R8P - uu - vv) * fp%vertex(1) + uu * fp%vertex(2) + vv * fp%vertex(3)
print '(A,I0,A,F7.4,A,F6.4,A,F6.4,A,F6.4)', 'facet ', hit%facet_id, ': t ', t, ', weights ', 1._R8P - uu - vv, ' ', uu, &
                                            ' ', vv
print '(A,L1)', 'eye + t dir = hit%point = barycentric point, to 1e-12: ', &
                max(dist(eye + t * dir, hit%point), dist(p, hit%point)) < 1.e-12_R8P
write(caption, '(A,I0,A,F6.4,A,F5.3,A,F5.3,A,F5.3,A)') 'facet ', hit%facet_id, ': t = ', t, ', weights (', &
                                                         1._R8P - uu - vv, ', ', uu, ', ', vv, ')'
az = atan2(eye%y - hit%point%y, eye%x - hit%point%x) * 180._R8P / acos(-1._R8P) + 20._R8P   ! from near the eye
el = atan2(eye%z - hit%point%z, hypot(eye%x - hit%point%x, eye%y - hit%point%y)) * 180._R8P / acos(-1._R8P) + 10._R8P
call fig_open(u, 'ref-ray-queries-facet', trim(caption), az, el)
call fig_zoom(u, hit%point, 1.3_R8P * maxval([dist(fp%vertex(1), p), dist(fp%vertex(2), p), dist(fp%vertex(3), p)]))
call fig_surface(u, 'body', dragon)
lift = 0.002_R8P * (eye - hit%point)                                ! just in front of the facet
call fig_triangle(u, 'frac', p + lift, fp%vertex(2) + lift, fp%vertex(3) + lift, value=min(1._R8P, 2._R8P * (1._R8P - uu - vv)))
call fig_triangle(u, 'frac', p + lift, fp%vertex(3) + lift, fp%vertex(1) + lift, value=min(1._R8P, 2._R8P * uu))
call fig_triangle(u, 'frac', p + lift, fp%vertex(1) + lift, fp%vertex(2) + lift, value=min(1._R8P, 2._R8P * vv))
call fig_segment(u, 'accent', hit%point - 0.3_R8P * (hit%point - eye), hit%point)
call fig_points(u, 'foot', [hit%point])
write(text, '(A,F5.3)') '1 - u - v = ', 1._R8P - uu - vv ; call fig_label(u, 'foot', fp%vertex(1), trim(text))
write(text, '(A,F5.3)') 'u = ', uu ; call fig_label(u, 'foot', fp%vertex(2), trim(text))
write(text, '(A,F5.3)') 'v = ', vv ; call fig_label(u, 'foot', fp%vertex(3), trim(text))
call fig_close(u)
contains
   real(R8P) function dist(a, b)
   !< Distance between two points.
   type(vector_R8P), intent(in) :: a, b

   dist = sqrt((a%x - b%x)**2 + (a%y - b%y)**2 + (a%z - b%z)**2)
   endfunction dist
endprogram fig_ref_ray_queries_picking
