!< FOSSIL test: coplanar overlap cuts and retriangulation (issue #33, step 3).
program fossil_test_coplanar_cuts
!< FOSSIL test: coplanar overlap cuts and retriangulation (issue #33, step 3).
!<
!<  - `coplanar_portion`: a segment crossing a triangle gives the analytic chord; a segment along
!<    an edge gives the shared part; a segment outside gives nothing.
!<  - Arrangement of two coplanar triangles T1 = (0,0)-(2,0)-(0,2) and T2 = T1 + (0.5, 0) (bottom
!<    edges collinear, so the cuts overlap the facets' edges): the CDT succeeds; each facet's
!<    sub-triangles cover it exactly (area sum = facet area); the sub-triangles of T1 inside T2 add
!<    up to the analytic overlap 1.125, and so do those of T2 inside T1.

use fossil_arrangement,  only : arrangement_t, arrangement_initialize, arrangement_collect_intersections, &
                                arrangement_retriangulate, ARR_STATUS_OK
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_tritri_exact, only : coplanar_portion
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

type(vector_R8P)                :: t1(3), t2(3), p, q
type(facet_object)              :: fa(1), fb(1)
type(aabb_tree_object)          :: tree_a, tree_b
type(arrangement_t)             :: arr
type(facet_object), allocatable :: sub(:)
integer(I4P),       allocatable :: owner(:), source(:)
integer(I4P)                    :: st, k
real(R8P)                       :: area_a, area_b, over_a, over_b
logical                         :: ok, found

ok = .true.
t1 = [v(0._R8P, 0._R8P), v(2._R8P, 0._R8P), v(0._R8P, 2._R8P)]
t2 = [v(0.5_R8P, 0._R8P), v(2.5_R8P, 0._R8P), v(0.5_R8P, 2._R8P)]

! coplanar_portion
call coplanar_portion(u=v(-1._R8P, 0.5_R8P), v=v(3._R8P, 0.5_R8P), t=t1, found=found, p=p, q=q)
ok = ok .and. found .and. ends(p, q, v(0._R8P, 0.5_R8P), v(1.5_R8P, 0.5_R8P))
print '(A,L1)', 'chord through T1: ', found .and. ends(p, q, v(0._R8P, 0.5_R8P), v(1.5_R8P, 0.5_R8P))
call coplanar_portion(u=t2(1), v=t2(2), t=t1, found=found, p=p, q=q)
ok = ok .and. found .and. ends(p, q, v(0.5_R8P, 0._R8P), v(2._R8P, 0._R8P))
print '(A,L1)', 'collinear edge part: ', found .and. ends(p, q, v(0.5_R8P, 0._R8P), v(2._R8P, 0._R8P))
call coplanar_portion(u=v(3._R8P, 3._R8P), v=v(4._R8P, 3._R8P), t=t1, found=found, p=p, q=q)
ok = ok .and. .not. found
print '(A,L1)', 'segment outside: none ', .not. found

! arrangement of the two triangles (empty trees: brute-force pair loop)
fa(1)%vertex = t1 ; call fa(1)%compute_metrix
fb(1)%vertex = t2 ; call fb(1)%compute_metrix
call arrangement_initialize(arr=arr, facet_a=fa, tree_a=tree_a, facet_b=fb, tree_b=tree_b)
call arrangement_collect_intersections(arr=arr)
call arrangement_retriangulate(arr=arr, sub_facet=sub, sub_owner=owner, sub_source=source, status=st)
ok = ok .and. st == ARR_STATUS_OK
area_a = 0._R8P ; area_b = 0._R8P ; over_a = 0._R8P ; over_b = 0._R8P
do k = 1, size(sub)
   if (owner(k) == 1) then
      area_a = area_a + sub(k)%area()
      if (inside(sub(k)%centroid, t2)) over_a = over_a + sub(k)%area()
   else
      area_b = area_b + sub(k)%area()
      if (inside(sub(k)%centroid, t1)) over_b = over_b + sub(k)%area()
   endif
enddo
print '(A,I0,A,I0,A,4F12.8)', 'arrangement status ', st, '  sub-triangles ', size(sub), &
      '  areas A, B, overlap in A, in B: ', area_a, area_b, over_a, over_b
ok = ok .and. abs(area_a - 2._R8P) < 1.e-14_R8P .and. abs(area_b - 2._R8P) < 1.e-14_R8P .and. &
     abs(over_a - 1.125_R8P) < 1.e-14_R8P .and. abs(over_b - 1.125_R8P) < 1.e-14_R8P

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   pure function v(x, y) result(r)
   real(R8P), intent(in) :: x, y
   type(vector_R8P)      :: r

   r%x = x ; r%y = y ; r%z = 0._R8P
   endfunction v

   pure function ends(a, b, e1, e2) result(s)
   type(vector_R8P), intent(in) :: a, b, e1, e2
   logical                      :: s

   s = (near(a, e1) .and. near(b, e2)) .or. (near(a, e2) .and. near(b, e1))
   endfunction ends

   pure function near(a, b) result(s)
   type(vector_R8P), intent(in) :: a, b
   logical                      :: s

   s = abs(a%x - b%x) < 1.e-14_R8P .and. abs(a%y - b%y) < 1.e-14_R8P .and. abs(a%z - b%z) < 1.e-14_R8P
   endfunction near

   pure function inside(c, t) result(s)
   !< Point strictly inside triangle t (z = 0 plane), by the signs of the three edge cross products.
   type(vector_R8P), intent(in) :: c, t(3)
   logical                      :: s
   real(R8P)                    :: d1, d2, d3

   d1 = (t(2)%x - t(1)%x) * (c%y - t(1)%y) - (t(2)%y - t(1)%y) * (c%x - t(1)%x)
   d2 = (t(3)%x - t(2)%x) * (c%y - t(2)%y) - (t(3)%y - t(2)%y) * (c%x - t(2)%x)
   d3 = (t(1)%x - t(3)%x) * (c%y - t(3)%y) - (t(1)%y - t(3)%y) * (c%x - t(3)%x)
   s = (d1 > 0._R8P .and. d2 > 0._R8P .and. d3 > 0._R8P) .or. (d1 < 0._R8P .and. d2 < 0._R8P .and. d3 < 0._R8P)
   endfunction inside
endprogram fossil_test_coplanar_cuts
