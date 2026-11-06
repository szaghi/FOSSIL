!< FOSSIL test: medial thickness and gap, the shrinking ball (issue #62, roadmap #58 WP4, step 2).
!<
!< 1. A plate 10 x 10 x t (t = 0.1): medial thickness t on the facets of the large faces to 1e-12, equal to the axial ray
!<    thickness (rays=1) to 1e-14; gap huge everywhere.
!< 2. A faceted unit sphere (1280 facets): medial thickness 2 to the faceting, from below: 6e-2 (measured 4.8e-2; 1.8e-2
!<    on 5120 facets, O(h): the ball is tangent at the facet centroid, off the circumcentre, so its centre misses the
!<    sphere's by O(h)); never above the axial ray chord, its starting bound.
!< 3. An extruded NACA 0012 (100 points per side, span 4): on the side facets with 0.1 < x < 0.9 and 1.5 < z < 2.5, the
!<    closed form 2 y_t sqrt(1 + y_t'^2), the diameter of the circle centred on the chord and tangent at (x, +-y_t), to
!<    the facet resolution: 2e-3 (measured 7.3e-4).
!< 4. Two unit boxes 0.2 apart: medial gap 0.2 on the facing faces to 1e-12, huge on the others.
!< 5. Bit-identical thickness and gap with 1 and 4 OpenMP threads (the sphere).
program fossil_test_thickness_medial

use fossil,                  only : surface_stl_object, facet_object, STATUS_OK, THICKNESS_MEDIAL
use fossil_thickness_bodies, only : make_boxes, make_sphere, make_wing, naca_half_thickness
use penf,                    only : I4P, R8P
!$ use omp_lib,              only : omp_set_num_threads, omp_get_max_threads

implicit none

real(R8P), parameter            :: T = 0.1_R8P, G = 0.2_R8P
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
real(R8P), allocatable          :: th(:), gp(:), th2(:), gp2(:)
real(R8P)                       :: err, err1, yt, dyt
integer(I4P)                    :: i, status, invalid, nt, n
logical                         :: ok, all_ok

all_ok = .true.

! ---- 1. plate
call make_boxes([0._R8P, 0._R8P, 0._R8P], [10._R8P, 10._R8P, T], 1_I4P, G, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL, invalid=invalid, status=status)
call s%thickness(thickness=th2, rays=1_I4P)
ok = status == STATUS_OK .and. invalid == 0 .and. all(gp == huge(1._R8P))
err = 0._R8P ; err1 = 0._R8P
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if (abs(fp%normal%z) < 0.5_R8P) cycle
   err = max(err, abs(th(i) - T) / T) ; err1 = max(err1, abs(th(i) - th2(i)) / T)
enddo
ok = ok .and. err <= 1e-12_R8P .and. err1 <= 1e-14_R8P
print '(A,2ES10.2,A,L1)', '1. plate: medial rel err, medial vs axial ray ', err, err1, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. sphere
call make_sphere(3_I4P, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL, invalid=invalid)
call s%thickness(thickness=th2, rays=1_I4P)
err = maxval(abs(th - 2._R8P)) / 2._R8P
ok = invalid == 0 .and. all(gp == huge(1._R8P)) .and. err <= 6e-2_R8P .and. maxval(th) <= 2._R8P .and. &
     all(th <= th2 * (1._R8P + 1e-12_R8P))
print '(A,ES10.2,A,L1)', '2. sphere: medial rel err (below 2, below the axial chord) ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. threads
nt = 1
!$ nt = omp_get_max_threads()
!$ call omp_set_num_threads(1)
call s%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL)
!$ call omp_set_num_threads(4)
call s%thickness(thickness=th2, gap=gp2, method=THICKNESS_MEDIAL)
!$ call omp_set_num_threads(nt)
ok = all(th == th2) .and. all(gp == gp2)
print '(A,L1)', '5. 1 and 4 threads bit-identical: ', ok
all_ok = all_ok .and. ok

! ---- 3. NACA 0012 wing
call make_wing(0.12_R8P, 100_I4P, 4._R8P, 8_I4P, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, method=THICKNESS_MEDIAL, invalid=invalid)
err = 0._R8P ; n = 0
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if (abs(fp%normal%z) > 0.5_R8P .or. fp%centroid%x < 0.1_R8P .or. fp%centroid%x > 0.9_R8P .or. &
       fp%centroid%z < 1.5_R8P .or. fp%centroid%z > 2.5_R8P) cycle
   call naca_half_thickness(0.12_R8P, fp%centroid%x, yt, dyt)
   err = max(err, abs(th(i) - 2._R8P * yt * sqrt(1._R8P + dyt ** 2)) / (2._R8P * yt)) ; n = n + 1
enddo
ok = invalid == 0 .and. n > 0 .and. err <= 2e-3_R8P
print '(A,I0,A,ES10.2,A,L1)', '3. NACA 0012 wing: ', n, ' facets, rel err vs closed form ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. two boxes
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, G, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(gap=gp, method=THICKNESS_MEDIAL, invalid=invalid)
ok = invalid == 0
err = 0._R8P
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if ((fp%normal%x > 0.5_R8P .and. fp%centroid%x < 1.5_R8P) .or. (fp%normal%x < -0.5_R8P .and. fp%centroid%x > 0.5_R8P)) &
      then
      err = max(err, abs(gp(i) - G) / G)
   else
      ok = ok .and. gp(i) == huge(1._R8P)
   endif
enddo
ok = ok .and. err <= 1e-12_R8P
print '(A,ES10.2,A,L1)', '4. two boxes: facing medial gap rel err ', err, ', others huge  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok
endprogram fossil_test_thickness_medial
