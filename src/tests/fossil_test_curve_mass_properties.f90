!< FOSSIL test: mass properties of the region bounded by a curve, per unit depth (issue #59, roadmap #58 WP1, step 3).
!<
!< Second moments are [I_uu, I_vv, I_uv] = [int (u-uc)^2, int (v-vc)^2, int (u-uc)(v-vc)] dA about the centroid.
!< 1. Rectangle a x b away from the origin: area, centroid, A a^2/12, A b^2/12, 0 (1e-15).
!< 2. Right triangle (0,0), (a,0), (0,b): centroid (a/3, b/3), A a^2/18, A b^2/18, -A a b/36 (1e-15).
!< 3. L-shape with integer coordinates: the two rectangles moved to the common centroid (1e-15).
!< 4. Annulus (outer regular N-gon CCW, inner CW): the polygon formula J = N R^4 sin(2 pi/N)(2 + cos(2 pi/N))/12 of each
!<    loop (1e-14), and J -> pi (R^4 - r^4)/2 at second order in 1/N.
!< 5. Ellipse polygon (vertices on the ellipse): I_vv -> pi a b^3/4 at second order.
!< 6. The L-shape rotated by 0.6 rad: R M R^T, polar moment unchanged, principal values and angle (1e-14); then moved
!<    by 1e6: the same moments as the rounded polygon moved back exactly (1e-12, the conditioning of the sums).
!< 7. Rejections: an open loop and a clockwise outer loop give STATUS_INVALID_INPUT.
!< 8. Extrusion: the 3-D prism of the L-shape, depth d, has mass A d, I_zz = d J and I_xx = d I_vv + A d^3/12 (1e-14):
!<    the 2-D and 3-D paths agree.
program fossil_test_curve_mass_properties

use fossil, only : curve_object, surface_stl_object, facet_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
real(R8P), parameter :: LU(6) = [0._R8P, 2._R8P, 2._R8P, 1._R8P, 1._R8P, 0._R8P] !< L-shape, u.
real(R8P), parameter :: LV(6) = [0._R8P, 0._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P] !< L-shape, v.
type(curve_object)       :: c
type(surface_stl_object) :: s
real(R8P)                :: mass, area, center(2), sm(3), ref(3), polar, pmom(2), angle, e(4), order, j(3)
real(R8P)                :: a, b, m1, m2, c1(2), c2(2), cc(2), rot(2,2), mm(2,2), inertia(3,3), d
real(R8P), allocatable   :: u(:), v(:)
integer(I4P)             :: status, k, n, i
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. rectangle a x b at (3, -2)
a = 2._R8P ; b = 0.5_R8P
call c%adopt_loops(u=[3._R8P, 3._R8P + a, 3._R8P + a, 3._R8P], v=[-2._R8P, -2._R8P, -2._R8P + b, -2._R8P + b], &
                   loop_sizes=[4_I4P])
call c%mass_properties(mass=mass, center=center, second_moments=sm, density=3._R8P, area=area, status=status)
ref = 3._R8P * a * b * [a * a / 12._R8P, b * b / 12._R8P, 0._R8P]
ok = status == STATUS_OK .and. abs(area - a * b) <= 1e-15_R8P .and. abs(mass - 3._R8P * a * b) <= 1e-15_R8P * 3._R8P .and. &
     maxval(abs(center - [3._R8P + a / 2, -2._R8P + b / 2])) <= 1e-15_R8P * 4._R8P .and. relerr(sm, ref) <= 1e-15_R8P
print '(A,ES10.3,A,L1)', '1. rectangle: moments rel err ', relerr(sm, ref), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. right triangle
a = 3._R8P ; b = 2._R8P
call c%adopt_loops(u=[0._R8P, a, 0._R8P], v=[0._R8P, 0._R8P, b], loop_sizes=[3_I4P])
call c%mass_properties(center=center, second_moments=sm, area=area, status=status)
ref = area * [a * a / 18._R8P, b * b / 18._R8P, -a * b / 36._R8P]
ok = status == STATUS_OK .and. maxval(abs(center - [a / 3, b / 3])) <= 1e-15_R8P * a .and. relerr(sm, ref) <= 1e-15_R8P
print '(A,ES10.3,A,L1)', '2. triangle: moments rel err ', relerr(sm, ref), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. L-shape: [0,2]x[0,1] and [0,1]x[1,2]
call c%adopt_loops(u=LU, v=LV, loop_sizes=[6_I4P])
call c%mass_properties(center=center, second_moments=sm, area=area, status=status)
m1 = 2._R8P ; c1 = [1._R8P, 0.5_R8P] ; m2 = 1._R8P ; c2 = [0.5_R8P, 1.5_R8P]
cc = (m1 * c1 + m2 * c2) / (m1 + m2)
ref = m1 * [4._R8P / 12, 1._R8P / 12, 0._R8P] + m1 * shift(c1 - cc) + m2 * [1._R8P / 12, 1._R8P / 12, 0._R8P] + &
      m2 * shift(c2 - cc)
ok = status == STATUS_OK .and. abs(area - 3._R8P) <= 1e-15_R8P * 3 .and. maxval(abs(center - cc)) <= 1e-15_R8P * 2 .and. &
     relerr(sm, ref) <= 1e-15_R8P
print '(A,ES10.3,A,L1)', '3. L-shape: moments rel err ', relerr(sm, ref), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. annulus R = 1, r = 0.4: the polygon formula, then convergence
do k = 1, 3
   n = 32 * 2 ** k
   call make_annulus(c, n, 1._R8P, 0.4_R8P)
   call c%mass_properties(polar_moment=polar, center=center, status=status)
   j(k) = polar
   e(k) = abs(polar - jpoly(n, 1._R8P) + jpoly(n, 0.4_R8P)) / polar
enddo
order = log(abs(j(1) - PI * (1._R8P - 0.4_R8P ** 4) / 2) / abs(j(2) - PI * (1._R8P - 0.4_R8P ** 4) / 2)) / log(2._R8P)
ok = status == STATUS_OK .and. maxval(e(1:3)) <= 1e-14_R8P .and. order > 1.9_R8P .and. order < 2.1_R8P .and. &
     maxval(abs(center)) <= 1e-15_R8P
print '(A,ES10.3,A,F5.2,A,L1)', '4. annulus: polygon formula err ', maxval(e(1:3)), '  order ', order, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. ellipse a = 2, b = 1: I_vv -> pi a b^3/4
a = 2._R8P ; b = 1._R8P
do k = 1, 3
   n = 32 * 2 ** k
   u = [(a * cos(2 * PI * (i - 1) / n), i = 1, n)] ; v = [(b * sin(2 * PI * (i - 1) / n), i = 1, n)]
   call c%adopt_loops(u=u, v=v, loop_sizes=[n])
   call c%mass_properties(second_moments=sm, status=status)
   e(k) = abs(sm(2) - PI * a * b ** 3 / 4) / (PI * a * b ** 3 / 4)
enddo
order = log(e(2) / e(3)) / log(2._R8P)
ok = status == STATUS_OK .and. order > 1.9_R8P .and. order < 2.1_R8P
print '(A,3ES10.3,A,F5.2,A,L1)', '5. ellipse: I_vv errs ', e(1:3), '  order ', order, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 6. rotation and translation
call c%adopt_loops(u=LU, v=LV, loop_sizes=[6_I4P])
call c%mass_properties(second_moments=ref, polar_moment=j(1), status=status)
rot = reshape([cos(0.6_R8P), sin(0.6_R8P), -sin(0.6_R8P), cos(0.6_R8P)], [2, 2])
u = rot(1,1) * LU + rot(1,2) * LV ; v = rot(2,1) * LU + rot(2,2) * LV          ! rotated about the origin
call c%adopt_loops(u=u, v=v, loop_sizes=[6_I4P])
call c%mass_properties(second_moments=sm, polar_moment=polar, principal_moments=pmom, principal_angle=angle, status=k)
mm = matmul(rot, matmul(reshape([ref(1), ref(3), ref(3), ref(2)], [2, 2]), transpose(rot)))
e(1) = relerr(sm, [mm(1,1), mm(2,2), mm(1,2)])
e(2) = abs(polar - j(1)) / j(1)
e(3) = abs(sum(pmom) - polar) / polar                             ! the principal values add up to the polar moment
e(4) = abs(cos(angle) * (sm(1) * cos(angle) + sm(3) * sin(angle)) + sin(angle) * (sm(3) * cos(angle) + &
       sm(2) * sin(angle)) - pmom(2)) / pmom(2)                   ! the angle's axis carries the larger value
ok = status == STATUS_OK .and. k == STATUS_OK .and. e(1) <= 1e-14_R8P .and. e(2) <= 1e-14_R8P .and. &
     e(3) <= 1e-14_R8P .and. e(4) <= 1e-14_R8P
print '(A,4ES10.3,A,L1)', '6. rotated L: moments, polar, principal errs ', e, '  ok ', ok
all_ok = all_ok .and. ok
! Moved by 1e6: the coordinates round to the spacing of doubles there (1.2e-10), so the polygon is compared with the
! same rounded polygon moved back exactly (Sterbenz): what is checked is the conditioning of the sums, not the input.
u = 1.e6_R8P + u ; v = -1.e6_R8P + v
call c%adopt_loops(u=u, v=v, loop_sizes=[6_I4P])
call c%mass_properties(second_moments=sm, polar_moment=polar, status=k)
call c%adopt_loops(u=u - 1.e6_R8P, v=v + 1.e6_R8P, loop_sizes=[6_I4P])
call c%mass_properties(second_moments=ref, polar_moment=j(1), status=status)
e(1) = relerr(sm, ref) ; e(2) = abs(polar - j(1)) / j(1)
ok = status == STATUS_OK .and. k == STATUS_OK .and. e(1) <= 1e-12_R8P .and. e(2) <= 1e-12_R8P
print '(A,2ES10.3,A,L1)', '6. the rotated L at 1e6: moments, polar errs ', e(1:2), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 7. rejections
call c%adopt_loops(u=LU, v=LV, loop_sizes=[6_I4P], closed=[.false.])
call c%mass_properties(mass=mass, status=status)
ok = status == STATUS_INVALID_INPUT .and. mass == 0._R8P
call c%adopt_loops(u=LU(6:1:-1), v=LV(6:1:-1), loop_sizes=[6_I4P])
call c%mass_properties(mass=mass, status=status)
ok = ok .and. status == STATUS_INVALID_INPUT .and. mass == 0._R8P
print '(A,L1)', '7. open and clockwise loops rejected: ', ok
all_ok = all_ok .and. ok

! ---- 8. extrusion: the 3-D L prism of depth d against the 2-D moments
d = 0.5_R8P
call c%adopt_loops(u=LU, v=LV, loop_sizes=[6_I4P])
call c%mass_properties(second_moments=sm, polar_moment=polar, area=area, status=status)
call make_l_prism(s, d)
call s%mass_properties(mass=mass, inertia=inertia, status=k)
e(1) = abs(mass - area * d) / (area * d)
e(2) = abs(inertia(3,3) - d * polar) / (d * polar)
e(3) = abs(inertia(1,1) - (d * sm(2) + area * d ** 3 / 12)) / inertia(1,1)
e(4) = abs(inertia(1,2) + d * sm(3)) / inertia(1,1)
ok = status == STATUS_OK .and. k == STATUS_OK .and. maxval(e) <= 1e-14_R8P
print '(A,4ES10.3,A,L1)', '8. extruded L: mass, I_zz, I_xx, I_xy errs ', e, '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   pure function relerr(x, y) result(r)
   !< Max error relative to the largest entry of y.
   real(R8P), intent(in) :: x(:), y(:)
   real(R8P)             :: r

   r = maxval(abs(x - y)) / maxval(abs(y))
   endfunction relerr

   pure function shift(dd) result(t)
   !< Parallel-axis terms of the second moments per unit area.
   real(R8P), intent(in) :: dd(2)
   real(R8P)             :: t(3)

   t = [dd(1) * dd(1), dd(2) * dd(2), dd(1) * dd(2)]
   endfunction shift

   pure function jpoly(nn, r) result(jj)
   !< Polar moment of the regular nn-gon of circumradius r about its centre.
   integer(I4P), intent(in) :: nn
   real(R8P),    intent(in) :: r
   real(R8P)                :: jj

   jj = nn * r ** 4 * sin(2 * PI / nn) * (2._R8P + cos(2 * PI / nn)) / 12._R8P
   endfunction jpoly

   subroutine make_annulus(cv, nn, ro, ri)
   !< Outer regular nn-gon of radius ro counter-clockwise, inner one of radius ri clockwise.
   type(curve_object), intent(inout) :: cv
   integer(I4P),       intent(in)    :: nn
   real(R8P),          intent(in)    :: ro, ri
   real(R8P)                         :: t(nn)
   integer(I4P)                      :: i

   t = [(2 * PI * (i - 1) / nn, i = 1, nn)]
   call cv%adopt_loops(u=[ro * cos(t), ri * cos(t(nn:1:-1))], v=[ro * sin(t), ri * sin(t(nn:1:-1))], loop_sizes=[nn, nn])
   endsubroutine make_annulus

   subroutine make_l_prism(surf, depth)
   !< The L-shape extruded over z in [0, depth]: caps fanned from the reflex vertex 4, sides outward.
   type(surface_stl_object), intent(inout) :: surf
   real(R8P),                intent(in)    :: depth
   integer,                  parameter     :: T(3,4) = reshape([4, 5, 6, 4, 6, 1, 4, 1, 2, 4, 2, 3], [3, 4])
   type(facet_object), allocatable         :: f(:)
   integer                                 :: i, jj, m

   allocate(f(2 * 4 + 2 * 6))
   m = 0
   do i = 1, 4
      m = m + 1 ; call tri(f(m), pt(T(1,i), depth), pt(T(2,i), depth), pt(T(3,i), depth))
      m = m + 1 ; call tri(f(m), pt(T(1,i), 0._R8P), pt(T(3,i), 0._R8P), pt(T(2,i), 0._R8P))
   enddo
   do i = 1, 6
      jj = mod(i, 6) + 1
      m = m + 1 ; call tri(f(m), pt(i, 0._R8P), pt(jj, 0._R8P), pt(jj, depth))
      m = m + 1 ; call tri(f(m), pt(i, 0._R8P), pt(jj, depth), pt(i, depth))
   enddo
   call surf%adopt_facets(facets=f)
   endsubroutine make_l_prism

   pure function pt(i, z) result(p)
   !< Vertex i of the L-shape at height z.
   integer,   intent(in) :: i
   real(R8P), intent(in) :: z
   type(vector_R8P)      :: p

   p = vector_R8P(LU(i), LV(i), z)
   endfunction pt

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   type(vector_R8P),   intent(in)    :: p1, p2, p3

   f%vertex(1) = p1 ; f%vertex(2) = p2 ; f%vertex(3) = p3
   call f%compute_metrix
   endsubroutine tri
endprogram fossil_test_curve_mass_properties
