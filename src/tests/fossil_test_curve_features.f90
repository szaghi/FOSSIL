!< FOSSIL test: corners of curves, corner tags and distances (issue #61, roadmap #58 WP3, step 4).
!<
!< 1. A unit square: 4 corners. A regular 64-gon: none (5.6 degrees per vertex). With angle = 1 degree: all 64.
!< 2. A NACA 0012 section (200 points, cosine spacing, closed trailing edge): the trailing edge only, at 30 degrees.
!< 3. An open polyline of 3 points with a 10-degree bend: its two ends.
!< 4. tag_block(feature_band=1, feature_tags=) on the square, 10 x 10 cells over [-0.25, 1.25]^2: 0 on the cells holding
!<    a corner (a corner on a cell edge marks both cells), the Chebyshev distance elsewhere, against a brute force.
!< 5. query_many(feature_distances=): closed forms on the square ((0.5, 0.5): sqrt(0.5); (2, 0): 1), and query_many_device
!<    equal to the CPU (host build: the device code path).
program fossil_test_curve_features

use fossil, only : curve_object, STATUS_OK
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(curve_object)        :: c
integer(I4P), allocatable :: corners(:)
real(R8P), allocatable    :: u(:), v(:), t(:)
real(R8P)                 :: x(11), p(2, 3), d(3), fd(3), fd2(3), xc(4), yc(4)
integer(I4P)              :: tags(10, 10), ft(10, 10), ref(10, 10), status, i, j, k, m
logical                   :: ok, all_ok

all_ok = .true.

! ---- 1. square and polygon
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%feature_corners(corners, status=status)
ok = status == STATUS_OK .and. size(corners) == 4
t = [(2 * PI * (i - 1) / 64, i = 1, 64)]
call c%adopt_loops(u=cos(t), v=sin(t), loop_sizes=[64_I4P])
call c%feature_corners(corners)
ok = ok .and. size(corners) == 0
call c%feature_corners(corners, angle=1._R8P)
ok = ok .and. size(corners) == 64
print '(A,L1)', '1. square 4 corners, 64-gon none (all at 1 degree): ', ok
all_ok = all_ok .and. ok

! ---- 2. NACA 0012
call naca(100_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[size(u, kind=I4P)])
call c%feature_corners(corners)
ok = size(corners) == 1
if (ok) ok = abs(u(corners(1)) - 1._R8P) <= 1e-15_R8P
print '(A,I0,A,L1)', '2. NACA 0012: corners ', size(corners), ' (the trailing edge)  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. open polyline
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P + cos(PI / 18)], v=[0._R8P, 0._R8P, sin(PI / 18)], loop_sizes=[3_I4P], &
                   closed=[.false.])
call c%feature_corners(corners)
ok = size(corners) == 2
if (ok) ok = corners(1) == 1 .and. corners(2) == 3
print '(A,L1)', '3. open polyline: its two ends: ', ok
all_ok = all_ok .and. ok

! ---- 4. corner tags on the square
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
x = [(-0.25_R8P + 0.15_R8P * (i - 1), i = 1, 11)]
call c%tag_block(x, x, 1_I4P, tags, status=status, feature_band=1_I4P, feature_tags=ft)
xc = [0._R8P, 1._R8P, 1._R8P, 0._R8P] ; yc = [0._R8P, 0._R8P, 1._R8P, 1._R8P]
do j = 1, 10
   do i = 1, 10
      ref(i, j) = 2
      do k = 1, 10            ! brute force: Chebyshev distance to the nearest cell holding a corner, capped at 2
         do m = 1, 10
            if (.not. holds(k, m)) cycle
            ref(i, j) = min(ref(i, j), max(abs(i - k), abs(j - m)))
         enddo
      enddo
   enddo
enddo
ok = status == STATUS_OK .and. all(ft == ref) .and. count(ft == 0) > 0
print '(A,I0,A,L1)', '4. corner tags: corner cells ', count(ft == 0), ', equal to brute force: ', ok
all_ok = all_ok .and. ok

! ---- 5. corner distances
p = reshape([0.5_R8P, 0.5_R8P, 2._R8P, 0._R8P, 1._R8P, 1._R8P], [2, 3])
call c%query_many(points=p, distances=d, feature_distances=fd, status=status)
ok = status == STATUS_OK .and. abs(fd(1) - sqrt(0.5_R8P)) <= 1e-15_R8P .and. abs(fd(2) - 1._R8P) <= 1e-15_R8P .and. &
     fd(3) == 0._R8P
call c%enter_device
call c%query_many_device(points=p, distances=d, feature_distances=fd2, status=status)
call c%exit_device
ok = ok .and. status == STATUS_OK .and. all(fd2 == fd)
print '(A,3ES12.4,A,L1)', '5. corner distances ', fd, ', device equal: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   logical function holds(k, m)
   !< Cell (k, m) holds a corner of the square (closed box).
   integer(I4P), intent(in) :: k, m
   integer(I4P)             :: q

   holds = .false.
   do q = 1, 4
      holds = holds .or. (xc(q) >= x(k) .and. xc(q) <= x(k + 1) .and. yc(q) >= x(m) .and. yc(q) <= x(m + 1))
   enddo
   endfunction holds

   subroutine naca(np, uu, vv)
   !< NACA 0012 of unit chord, np points per side (cosine spacing, closed trailing edge), counter-clockwise.
   integer(I4P),           intent(in)  :: np
   real(R8P), allocatable, intent(out) :: uu(:), vv(:)
   real(R8P)                           :: xs(np + 1), ts(np + 1)
   integer(I4P)                        :: q, n

   xs = [(0.5_R8P * (1._R8P - cos(PI * q / np)), q = 0, np)]
   ts = 0.6_R8P * (0.2969_R8P * sqrt(xs) - 0.1260_R8P * xs - 0.3516_R8P * xs ** 2 + 0.2843_R8P * xs ** 3 - &
                   0.1036_R8P * xs ** 4)
   allocate(uu(2 * np), vv(2 * np))
   n = 0
   do q = np + 1, 2, -1
      n = n + 1 ; uu(n) = xs(q) ; vv(n) = ts(q)
   enddo
   do q = 1, np
      n = n + 1 ; uu(n) = xs(q) ; vv(n) = -ts(q)
   enddo
   endsubroutine naca
endprogram fossil_test_curve_features
