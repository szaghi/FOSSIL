!< FOSSIL test: ray thickness and gap (issue #62, roadmap #58 WP4, step 1).
!<
!< 1. A plate 10 x 10 x t (t = 0.1): on the facets of the large faces, the thickness is t / cos(theta_16) (the lower
!<    median of t / cos(theta_k) over the 31 rays of the default 60-degree cone) to 1e-13; t exactly with rays=1; the gap
!<    is huge everywhere (every outward ray misses) and no facet is invalid.
!< 2. A faceted sphere of radius 1 (icosahedron subdivided 3 times, 1280 facets): thickness 2 cos(theta_16) to the
!<    faceting: 1e-2 (measured 6.1e-3 and 4.5e-3); with rays=1, 2 to the faceting.
!< 3. Two unit boxes at distance g = 0.2 (one surface, two shells): gap g / cos(theta_16) on the facing faces, huge on
!<    the others.
!< 4. An open surface (a box without the face x = 0): rays escape, `invalid` > 0, thickness huge there.
!< 5. intersect_ray_first(exclude=f) equals the first hit of intersect_ray_all other than f (facet and t, bit for bit).
!< 6. Bit-identical thickness and gap with 1 and 4 OpenMP threads; vertex fields inside the facet range.
program fossil_test_thickness_ray

use fossil,           only : surface_stl_object, facet_object, STATUS_OK, THICKNESS_RAY
use fossil_ray_query, only : ray_hit_t
use fossil_thickness, only : cone_cosines
use fossil_thickness_bodies, only : make_boxes, make_sphere
use penf,             only : I4P, R8P
use vecfor,           only : vector_R8P
!$ use omp_lib,       only : omp_set_num_threads, omp_get_max_threads

implicit none

real(R8P), parameter            :: T = 0.1_R8P, G = 0.2_R8P
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(ray_hit_t)                 :: hit
type(ray_hit_t), allocatable    :: hits(:)
type(vector_R8P)                :: d
real(R8P), allocatable          :: th(:), gp(:), th2(:), gp2(:), vt(:)
real(R8P)                       :: c(31), err, err1, nz
integer(I4P)                    :: i, k, status, invalid, nt
logical                         :: has_hit, ok, all_ok

all_ok = .true.
c = cone_cosines(60._R8P, 31_I4P)

! ---- 1. plate
call make_boxes([0._R8P, 0._R8P, 0._R8P], [10._R8P, 10._R8P, T], 1_I4P, G, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, gap=gp, invalid=invalid, status=status)
ok = status == STATUS_OK .and. invalid == 0 .and. all(gp == huge(1._R8P))
err = 0._R8P
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if (abs(fp%normal%z) > 0.5_R8P) err = max(err, abs(th(i) - T / c(16)) / T)
enddo
call s%thickness(thickness=th2, rays=1_I4P)
err1 = 0._R8P
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if (abs(fp%normal%z) > 0.5_R8P) err1 = max(err1, abs(th2(i) - T) / T)
enddo
ok = ok .and. err <= 1e-13_R8P .and. err1 <= 1e-14_R8P
print '(A,2ES10.2,A,L1)', '1. plate: rel err cone, axial ', err, err1, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. sphere
call make_sphere(3_I4P, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, gap=gp, invalid=invalid, vertex_thickness=vt)
err = maxval(abs(th - 2._R8P * c(16))) / 2._R8P
call s%thickness(thickness=th2, rays=1_I4P)
err1 = maxval(abs(th2 - 2._R8P)) / 2._R8P
ok = invalid == 0 .and. all(gp == huge(1._R8P)) .and. err <= 1e-2_R8P .and. err1 <= 1e-2_R8P
ok = ok .and. minval(vt) >= minval(th) .and. maxval(vt) <= maxval(th)
print '(A,I0,A,2ES10.2,A,L1)', '2. sphere (', s%get_facets_number(), ' facets): rel err cone, axial ', err, err1, &
      '  ok ', ok
all_ok = all_ok .and. ok

! ---- 6. threads (on the sphere)
nt = 1
!$ nt = omp_get_max_threads()
!$ call omp_set_num_threads(1)
call s%thickness(thickness=th, gap=gp)
!$ call omp_set_num_threads(4)
call s%thickness(thickness=th2, gap=gp2)
!$ call omp_set_num_threads(nt)
ok = all(th == th2) .and. all(gp == gp2)
print '(A,L1)', '6. 1 and 4 threads bit-identical: ', ok
all_ok = all_ok .and. ok

! ---- 3. two boxes
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, G, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(gap=gp, invalid=invalid)
ok = invalid == 0
err = 0._R8P
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if ((fp%normal%x > 0.5_R8P .and. fp%centroid%x < 1.5_R8P) .or. (fp%normal%x < -0.5_R8P .and. fp%centroid%x > 0.5_R8P)) &
      then
      err = max(err, abs(gp(i) - G / c(16)) / G)
   else
      ok = ok .and. gp(i) == huge(1._R8P)
   endif
enddo
ok = ok .and. err <= 1e-13_R8P
print '(A,ES10.2,A,L1)', '3. two boxes: facing-gap rel err ', err, ', others huge  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. open surface
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 1_I4P, G, f)
f = f(3:)                                   ! the face x = 0 removed
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, invalid=invalid)
ok = invalid > 0 .and. count(th == huge(1._R8P)) == invalid
print '(A,I0,A,L1)', '4. open surface: invalid facets ', invalid, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. exclude
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, G, f)
call s%adopt_facets(facets=f) ; call s%analyze
ok = .true.
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   do k = 1, 3
      nz = 0.3_R8P * k
      d = vector_R8P(-fp%normal%x + nz * fp%normal%y, -fp%normal%y + nz * fp%normal%z, -fp%normal%z + nz * fp%normal%x)
      call s%intersect_ray_first(ray_origin=fp%centroid, ray_direction=d, hit=hit, has_hit=has_hit, exclude=i)
      call s%intersect_ray_all(ray_origin=fp%centroid, ray_direction=d, hits=hits)
      hits = pack(hits, hits%facet_id /= i)
      if (size(hits) == 0) then
         ok = ok .and. .not. has_hit
      else
         ok = ok .and. has_hit .and. hit%facet_id == hits(1)%facet_id .and. hit%t == hits(1)%t
      endif
   enddo
enddo
print '(A,L1)', '5. exclude equals the first other hit of intersect_ray_all: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

endprogram fossil_test_thickness_ray
