!run figref-ray-queries fig_ref_ray_queries
!fig ref-ray-queries-all
!fig ref-ray-queries-shadow
!< Figures of the reference page "Ray queries", on dragon.stl.
!< - intersect_ray_all along +x from (-14, -1.5, 7.857): the dragon translucent, every hit marked (gold), the stretches
!<   of the ray inside the body (between an odd hit and the next) in ember, the stretches outside in teal.
!< - intersect_ray_any as a shadow ray: from the centre of each cell of a ground plane under the dragon towards a
!<   distant light; the cells whose ray is blocked are in shadow (sepia), the others lit (cream).
program fig_ref_ray_queries
use fossil,     only : surface_stl_object, ray_hit_t
use fossil_facet_object, only : facet_object
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: NP = 24                    ! pieces per stretch of the ray (painter's sort)
type(surface_stl_object)      :: dragon
type(facet_object), pointer   :: f
type(ray_hit_t), allocatable  :: hits(:)
type(vector_R8P)              :: origin, light, bmin, bmax, p
type(vector_R8P), allocatable :: a(:), b(:), pts(:)
type(vector_R8P)              :: s0, s1
real(R8P)                     :: z0, h, x, y
logical                       :: blocked
integer(I4P)                  :: u, i, j, k, n, shadowed
character(80)                 :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
origin = vector_R8P(-14._R8P, -1.5_R8P, 7.857_R8P)
call dragon%intersect_ray_all(ray_origin=origin, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A)', 'intersect_ray_all: ', size(hits), ' hits at t ='
print '(10F7.3)', hits(:)%t
do k = 1, size(hits)                                         ! enters (E) or leaves (L): the facet normal against +x
   f => dragon%facet_at(hits(k)%facet_id)
   write(*, '(A6)', advance='no') merge('E', 'L', f%normal%x < 0._R8P)
enddo
write(*, '(A)') ''
write(caption, '(A,I0,A)') 'intersect_ray_all: ', size(hits), ' hits; inside the body in ember, outside in teal'
call fig_open(u, 'ref-ray-queries-all', trim(caption), -60._R8P, 20._R8P)
call fig_surface(u, 'ghost', dragon)
allocate(pts(size(hits)))
do k = 1, size(hits)
   pts(k) = hits(k)%point
enddo
do k = 0, size(hits)                                         ! stretch k: from hit k (or the origin) to hit k + 1
   s0 = origin ; if (k > 0) s0 = hits(k)%point
   s1 = origin + (hits(size(hits))%t + 4._R8P) * ex_R8P ; if (k < size(hits)) s1 = hits(k + 1)%point
   allocate(a(NP), b(NP))
   do i = 1, NP
      a(i) = s0 + ((i - 1) / real(NP, R8P)) * (s1 - s0) ; b(i) = s0 + (i / real(NP, R8P)) * (s1 - s0)
   enddo
   call fig_segments(u, merge('patch', 'probe', mod(k, 2) == 1), a, b)
   deallocate(a, b)
enddo
call fig_points(u, 'foot', pts)
call fig_label(u, 'probe', origin, 'origin')
call fig_close(u)

bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
z0 = bmin%z - 0.5_R8P                                        ! the ground
light = vector_R8P(-0.4_R8P, 0.7_R8P, 0.8_R8P)                   ! towards the light (need not be unit)
h = 0.3_R8P ; shadowed = 0 ; n = 0
call fig_open(u, 'ref-ray-queries-shadow', 'intersect_ray_any towards the light: cells in shadow (sepia) and lit (cream)', &
              -60._R8P, 35._R8P)
call fig_surface(u, 'body', dragon)
do j = 0, nint((bmax%y - bmin%y + 12._R8P) / h) - 1               ! the shadow falls towards +x, -y
   do i = 0, nint((bmax%x - bmin%x + 9._R8P) / h) - 1
      x = bmin%x - 1._R8P + (i + 0.5_R8P) * h ; y = bmin%y - 11._R8P + (j + 0.5_R8P) * h
      p = x * ex_R8P + y * ey_R8P + z0 * ez_R8P
      call dragon%intersect_ray_any(ray_origin=p, ray_direction=light, found=blocked)
      n = n + 1 ; if (blocked) shadowed = shadowed + 1
      call fig_quad(u, 'frac', p - 0.5_R8P * h * (ex_R8P + ey_R8P), p + 0.5_R8P * h * (ex_R8P - ey_R8P), &
                    p + 0.5_R8P * h * (ex_R8P + ey_R8P), p - 0.5_R8P * h * (ex_R8P - ey_R8P), &
                    value=merge(0.85_R8P, 0.05_R8P, blocked))
   enddo
enddo
call fig_close(u)
print '(A,I0,A,I0,A)', 'intersect_ray_any: ', shadowed, ' of ', n, ' ground cells in shadow'
endprogram fig_ref_ray_queries
