!< FOSSIL test: 2-D hyperbolic extrusion (issue #74, step 6).
!<
!< 1. From a regular 64-gon on the unit circle (counter-clockwise), 20 geometric layers (h1 = 1e-3, ratio 1.2): its
!<    exact offsets, vertices on the rays at radii 1 + (h1 + ... + h_k) / cos(pi / 64), to 1e-13, every cell valid.
!< 2. From a NACA0012 section (129 points, cosine spacing, blunt trailing edge closed), 30 layers from h1 = 1e-4 with
!<    ratio 1.15: every cell certified valid, the first grid line normal to the wall (|cos| of the angle with the
!<    wall tangent below 1e-2 away from the trailing-edge corner).
!< 3. Concave corners, 25 layers (0.2 in all): a square with a V notch 0.15 deep folds without smoothing (epsilon 0,
!<    status 2) and not with the default smoothing; a right-angled notch (0.6 wide, 0.2 deep) folds with it too, once
!<    the height passes the point spacing at the corners: the status reports it (a limit, measured).
!< 4. From an open straight line (a flat plate, 11 points), 5 layers: the Cartesian grid exactly.
!< 5. Bad input (2 points, a non-positive height): status 1.
program fossil_test_hyperbolic_2d

use fossil, only : curvilinear_block_2d_t, hyperbolic_grid_2d, geometric_spacing
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: cu(64), cv(64), r, err, h(30), au(129), av(129), xc, yt, t(2), g(2), cmax, th
real(R8P), allocatable       :: su(:), sv(:)
integer(I4P)                 :: i, k, st, bad, bad0, n
logical                      :: all_ok

all_ok = .true.

! 1. circle -> polar grid
do i = 1, 64
   cu(i) = cos(2 * PI * (i - 1) / 64) ; cv(i) = sin(2 * PI * (i - 1) / 64)
enddo
h(1:20) = geometric_spacing(1e-3_R8P, 1.2_R8P, 20_I4P)
call hyperbolic_grid_2d(cu, cv, h(1:20), b, invalid=bad, status=st)
err = 0._R8P
do k = 1, 21
   r = 1._R8P + sum(h(1:k - 1)) / cos(PI / 64)
   do i = 1, 65
      err = max(err, abs(b%u(i, k) - r * cu(mod(i - 1, 64) + 1)), abs(b%v(i, k) - r * cv(mod(i - 1, 64) + 1)))
   enddo
enddo
call check(st == 0 .and. bad == 0 .and. b%ni == 65 .and. b%nj == 21 .and. err <= 1e-13_R8P, &
           '64-gon: its exact offsets on the rays, every cell valid (max error)', err)

! 2. NACA0012
do i = 1, 129 ! from the trailing edge along the upper side, round the nose, back along the lower side (CCW)
   th = 2 * PI * (i - 1) / 129._R8P
   xc = 0.5_R8P * (1._R8P + cos(th))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   au(i) = xc ; av(i) = merge(yt, -yt, i <= 65)
enddo
h = geometric_spacing(1e-4_R8P, 1.15_R8P, 30_I4P)
call hyperbolic_grid_2d(au, av, h, b, invalid=bad, status=st)
cmax = 0._R8P
do i = 4, 126
   t = [b%u(i + 1, 1) - b%u(i - 1, 1), b%v(i + 1, 1) - b%v(i - 1, 1)]
   g = [b%u(i, 2) - b%u(i, 1), b%v(i, 2) - b%v(i, 1)]
   cmax = max(cmax, abs(dot_product(t, g)) / (norm2(t) * norm2(g)))
enddo
err = 0._R8P
call check(st == 0 .and. bad == 0 .and. cmax <= 1e-2_R8P .and. err <= 1e-15_R8P, &
           'NACA0012: every cell valid, first line normal to the wall (max |cos|)', cmax)

! 3. concave corners
n = 0 ; allocate(su(400), sv(400))
call side(0._R8P, 0._R8P, 1._R8P, 0._R8P) ; call side(1._R8P, 0._R8P, 1._R8P, 1._R8P)
call side(1._R8P, 1._R8P, 0.5_R8P, 0.85_R8P) ; call side(0.5_R8P, 0.85_R8P, 0._R8P, 1._R8P)
call side(0._R8P, 1._R8P, 0._R8P, 0._R8P)
h(1:25) = geometric_spacing(2e-3_R8P, 1.1_R8P, 25_I4P)
call hyperbolic_grid_2d(su(1:n), sv(1:n), h(1:25), b, epsilon=0._R8P, invalid=bad0, status=st)
call check(st == 2 .and. bad0 > 0, 'V notch 0.15 deep without smoothing: folds (cells not certified valid)', real(bad0, R8P))
call hyperbolic_grid_2d(su(1:n), sv(1:n), h(1:25), b, invalid=bad, status=st)
call check(st == 0 .and. bad == 0, 'V notch 0.15 deep with the default smoothing: every cell valid', real(bad, R8P))
n = 0
call side(0._R8P, 0._R8P, 1._R8P, 0._R8P) ; call side(1._R8P, 0._R8P, 1._R8P, 1._R8P)
call side(1._R8P, 1._R8P, 0.8_R8P, 1._R8P) ; call side(0.8_R8P, 1._R8P, 0.8_R8P, 0.8_R8P)
call side(0.8_R8P, 0.8_R8P, 0.2_R8P, 0.8_R8P) ; call side(0.2_R8P, 0.8_R8P, 0.2_R8P, 1._R8P)
call side(0.2_R8P, 1._R8P, 0._R8P, 1._R8P) ; call side(0._R8P, 1._R8P, 0._R8P, 0._R8P)
call hyperbolic_grid_2d(su(1:n), sv(1:n), h(1:25), b, invalid=bad, status=st)
call check(st == 2 .and. bad > 0, 'right-angled notch: the folds are reported (the smoothing does not prevent them)', &
           real(bad, R8P))

! 4. open flat plate
call hyperbolic_grid_2d([(0.1_R8P * i, i = 0, 10)], [(0._R8P, i = 0, 10)], [0.1_R8P, 0.1_R8P, 0.1_R8P, 0.1_R8P, 0.1_R8P], &
                        b, closed=.false., status=st)
err = 0._R8P
do k = 1, 6
   do i = 1, 11
      err = max(err, abs(b%u(i, k) - 0.1_R8P * (i - 1)), abs(b%v(i, k) + 0.1_R8P * (k - 1)))
   enddo
enddo
call check(st == 0 .and. err <= 1e-15_R8P, 'open flat plate: the Cartesian grid (to the right of +x: -y)', err)

! 5. bad input
call hyperbolic_grid_2d([0._R8P, 1._R8P], [0._R8P, 0._R8P], [0.1_R8P], b, status=st)
call check(st == 1, 'two points: bad input', real(st, R8P))
call hyperbolic_grid_2d(cu, cv, [0.1_R8P, 0._R8P], b, status=st)
call check(st == 1, 'a zero height: bad input', real(st, R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES23.15)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   subroutine side(x0, y0, x1, y1)
   !< Append 10 points of the side from (x0, y0) toward (x1, y1) (its end is the next side's start).
   real(R8P), intent(in) :: x0, y0, x1, y1 !< Ends.
   integer(I4P)          :: j

   do j = 0, 9
      n = n + 1 ; su(n) = x0 + (x1 - x0) * j / 10._R8P ; sv(n) = y0 + (y1 - y0) * j / 10._R8P
   enddo
   endsubroutine side
endprogram fossil_test_hyperbolic_2d
