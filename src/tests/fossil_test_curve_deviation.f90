!< FOSSIL test: deviation and Hausdorff distance of curves (issue #60, roadmap #58 WP2, step 4).
!<
!< 1. A unit square against itself translated by (0.3, 0): h = 0.3 both ways (at vertices); certified bounds bracket it
!<    within the tolerance 1e-4.
!< 2. A regular 64-gon on the unit circle against itself scaled by 1.05: every vertex is on the circle, so
!<    0.98 * 0.05 <= h <= 0.05.
!< 3. The maximum inside a segment: the open polyline (0,0)-(3,0) against two tiny segments at x = 0 and x = 2.2: d is
!<    min(x, |x - 2.2|), largest at x = 1.1 (value 1.1), between the samples: the sampled bound misses it by more than
!<    1e-3, the certified bounds bracket it.
!< 4. The work cap: max_segments = 3 gives DEV_STATUS_NOT_CONVERGED with valid bounds.
!< 5. The extrusion link: the 64-gon and its scaled copy extruded to prisms (z in [0, 1]); the 3-D certified bracket of
!<    surface%deviation (tolerance 1e-3: the coplanar caps, fanned differently, refine along their crossing edges at a
!<    cost linear in 1 / tolerance) and the 2-D one overlap (the side distances are the 2-D ones).
!< 6. Device path: device=.true. without enter_device is refused; entered (host build: device code, directives inert),
!<    lower, upper and per_vertex equal the CPU path bit for bit.
!< 7. hausdorff of the two squares: both directions, 0.3.
program fossil_test_curve_deviation

use fossil, only : curve_object, surface_stl_object, facet_object, DEV_STATUS_OK, DEV_STATUS_BAD_INPUT, &
                   DEV_STATUS_NOT_CONVERGED
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P), TOL = 1.e-4_R8P
type(curve_object)       :: a, b
type(surface_stl_object) :: pa, pb
real(R8P)                :: lower, upper, lower2, upper2, e(4)
real(R8P), allocatable   :: t(:), pv(:), pv2(:)
integer(I4P)             :: status, status2, i, created
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. translated square
call a%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call b%adopt_loops(u=[0.3_R8P, 1.3_R8P, 1.3_R8P, 0.3_R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call a%deviation(b, lower=lower, upper=upper, tolerance=TOL, status=status)
call b%deviation(a, lower=lower2, upper=upper2, tolerance=TOL, status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. abs(lower - 0.3_R8P) <= 1e-15_R8P .and. &
     0.3_R8P <= upper .and. upper - lower <= TOL .and. abs(lower2 - 0.3_R8P) <= 1e-15_R8P .and. upper2 - lower2 <= TOL
print '(A,4ES12.4,A,L1)', '1. translated square: bounds both ways ', lower, upper, lower2, upper2, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 7. hausdorff
call a%hausdorff(b, lower=lower, upper=upper, tolerance=TOL, status=status)
ok = status == DEV_STATUS_OK .and. lower <= 0.3_R8P + 1e-15_R8P .and. 0.3_R8P <= upper .and. upper - lower <= TOL
print '(A,2ES12.4,A,L1)', '7. hausdorff of the squares ', lower, upper, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. scaled 64-gon
t = [(2 * PI * (i - 1) / 64, i = 1, 64)]
call a%adopt_loops(u=cos(t), v=sin(t), loop_sizes=[64_I4P])
call b%adopt_loops(u=1.05_R8P * cos(t), v=1.05_R8P * sin(t), loop_sizes=[64_I4P])
call a%deviation(b, lower=lower, upper=upper, tolerance=TOL, status=status)
ok = status == DEV_STATUS_OK .and. lower >= 0.98_R8P * 0.05_R8P .and. upper <= 0.05_R8P + TOL .and. upper - lower <= TOL
print '(A,2ES23.15,A,L1)', '2. scaled 64-gon: lower, upper ', lower, upper, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. maximum inside a segment
call a%adopt_loops(u=[0._R8P, 3._R8P], v=[0._R8P, 0._R8P], loop_sizes=[2_I4P], closed=[.false.])
call b%adopt_loops(u=[0._R8P, 0._R8P, 2.2_R8P, 2.2_R8P], v=[0._R8P, 1.e-7_R8P, 0._R8P, 1.e-7_R8P], &
                   loop_sizes=[2_I4P, 2_I4P], closed=[.false., .false.])
call a%deviation(b, lower=lower2, status=status2)
call a%deviation(b, lower=lower, upper=upper, tolerance=TOL, created=created, status=status)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. lower2 < 1.1_R8P - 1e-3_R8P .and. &
     lower <= 1.1_R8P + 1e-14_R8P .and. 1.1_R8P - 1e-7_R8P <= upper .and. upper - lower <= TOL
print '(A,ES12.4,A,2ES23.15,A,I0,A,L1)', '3. interior maximum: sampled ', lower2, '  certified ', lower, upper, &
                                         '  segments ', created, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. work cap
call a%deviation(b, lower=lower, upper=upper, tolerance=TOL, max_segments=3_I4P, status=status)
ok = status == DEV_STATUS_NOT_CONVERGED .and. lower <= 1.1_R8P + 1e-14_R8P .and. 1.1_R8P - 1e-7_R8P <= upper .and. &
     upper - lower > TOL
print '(A,2ES12.4,A,L1)', '4. capped: lower, upper ', lower, upper, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. extrusion link
call a%adopt_loops(u=cos(t), v=sin(t), loop_sizes=[64_I4P])
call b%adopt_loops(u=1.05_R8P * cos(t), v=1.05_R8P * sin(t), loop_sizes=[64_I4P])
call a%deviation(b, lower=lower, upper=upper, tolerance=TOL, status=status)
call make_prism(pa, cos(t), sin(t))
call make_prism(pb, 1.05_R8P * cos(t), 1.05_R8P * sin(t))
call pa%deviation(pb, lower=lower2, upper=upper2, tolerance=1.e-3_R8P, status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. max(lower, lower2) <= min(upper, upper2)
print '(A,4ES12.4,A,L1)', '5. extruded 64-gons: 2-D and 3-D brackets ', lower, upper, lower2, upper2, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 6. device path
call a%deviation(b, lower=lower, upper=upper, per_vertex=pv, tolerance=TOL, status=status)
call a%deviation(b, lower=lower2, upper=upper2, per_vertex=pv2, tolerance=TOL, device=.true., status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_BAD_INPUT
call b%enter_device
call a%deviation(b, lower=lower2, upper=upper2, per_vertex=pv2, tolerance=TOL, device=.true., status=status2)
call b%exit_device
ok = ok .and. status2 == DEV_STATUS_OK .and. lower2 == lower .and. upper2 == upper .and. all(pv2 == pv)
print '(A,L1)', '6. device path bit for bit: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine make_prism(surf, u, v)
   !< The convex polygon (u, v), counter-clockwise, extruded over z in [0, 1]; caps fanned from vertex 1.
   type(surface_stl_object), intent(inout) :: surf
   real(R8P),                intent(in)    :: u(:), v(:)
   type(facet_object), allocatable         :: f(:)
   integer                                 :: n, k, j, m

   n = size(u)
   allocate(f(2 * (n - 2) + 2 * n))
   m = 0
   do k = 2, n - 1
      m = m + 1 ; call tri(f(m), [u(1), v(1), 1._R8P], [u(k), v(k), 1._R8P], [u(k+1), v(k+1), 1._R8P])
      m = m + 1 ; call tri(f(m), [u(1), v(1), 0._R8P], [u(k+1), v(k+1), 0._R8P], [u(k), v(k), 0._R8P])
   enddo
   do k = 1, n
      j = mod(k, n) + 1
      m = m + 1 ; call tri(f(m), [u(k), v(k), 0._R8P], [u(j), v(j), 0._R8P], [u(j), v(j), 1._R8P])
      m = m + 1 ; call tri(f(m), [u(k), v(k), 0._R8P], [u(j), v(j), 1._R8P], [u(k), v(k), 1._R8P])
   enddo
   call surf%destroy
   call surf%adopt_facets(facets=f)
   endsubroutine make_prism

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_metrix
   endsubroutine tri
endprogram fossil_test_curve_deviation
