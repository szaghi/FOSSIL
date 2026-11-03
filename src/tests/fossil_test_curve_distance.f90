!< FOSSIL test: curve distance queries and self-intersections (issue #36 step 3).
program fossil_test_curve_distance
!< FOSSIL test: curve distance queries and self-intersections (issue #36 step 3).
!<
!< - Square [0,2]^2: exact signed distances inside, outside, in a vertex region; closest point,
!<   segment and vertex.
!< - L-shape with a reflex vertex at (1,1): the reflex-vertex region inside (negative), a medial
!<   point outside equidistant from two segments (+0.5).
!< - 64-gon, procedural NACA 0012 (cosine spacing, 200 vertices) and a cusp (tip angle ~1e-4 rad):
!<   2000 random points each; the BVH distance equals the brute-force one bit for bit, and the
!<   sign agrees with the exact winding number (|d| > 1e-12).
!< - Self-intersections (exact): a figure eight has one crossing; a loop folding back on itself
!<   has a fold-back and a touching pair; the square and the NACA 0012 have none.

use fossil, only : curve_object
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(curve_object)        :: c
real(R8P)                 :: d, q(2), cp(2), u(200), v(200)
integer(I4P)              :: i, n, sg, vx
integer(I4P), allocatable :: pairs(:,:)
logical                   :: ok(6)

! square
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
ok(1) = c%distance([1._R8P, 1._R8P]) == -1._R8P .and. c%distance([3._R8P, 1._R8P]) == 1._R8P .and. &
        c%distance([1._R8P, 0.5_R8P]) == -0.5_R8P .and. c%distance([3._R8P, 1._R8P], is_signed=.false.) == 1._R8P
call c%compute_distance(p=[3._R8P, 3._R8P], distance=d, closest_point=cp, closest_segment=sg, closest_vertex=vx)
ok(1) = ok(1) .and. d == sqrt(2._R8P) .and. all(cp == [2._R8P, 2._R8P]) .and. vx == 3
call c%compute_distance(p=[1._R8P, -1._R8P], distance=d, closest_point=cp, closest_segment=sg, closest_vertex=vx)
ok(1) = ok(1) .and. d == 1._R8P .and. all(cp == [1._R8P, 0._R8P]) .and. sg == 1 .and. vx == 0
print '(A,L1)', 'square: ', ok(1)

! L-shape, reflex vertex (1,1)
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P], &
                   loop_sizes=[6_I4P])
call c%compute_distance(p=[0.75_R8P, 0.75_R8P], distance=d, closest_vertex=vx)
ok(2) = d == -sqrt(0.125_R8P) .and. vx == 4 .and. c%distance([1.5_R8P, 1.5_R8P]) == 0.5_R8P .and. &
        c%distance([1.5_R8P, 0.5_R8P]) == -0.5_R8P
print '(A,L1)', 'L-shape (reflex vertex, medial point): ', ok(2)

! 64-gon
n = 64
do i = 1, n
   u(i) = cos(2._R8P * PI * (i - 1) / n) ; v(i) = sin(2._R8P * PI * (i - 1) / n)
enddo
call c%adopt_loops(u=u(1:n), v=v(1:n), loop_sizes=[n])
ok(3) = random_check(c, [-1.5_R8P, -1.5_R8P], [1.5_R8P, 1.5_R8P], '64-gon')

! NACA 0012
call naca0012(100_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(4) = random_check(c, [-0.2_R8P, -0.2_R8P], [1.2_R8P, 0.2_R8P], 'NACA 0012')

! cusp: tip at (10, 0.001), interior angle ~2e-4 rad
call c%adopt_loops(u=[0._R8P, 10._R8P, 0._R8P], v=[0._R8P, 0.001_R8P, 0.002_R8P], loop_sizes=[3_I4P])
ok(5) = random_check(c, [9._R8P, -0.01_R8P], [11._R8P, 0.012_R8P], 'cusp') .and. &
        c%distance([10.5_R8P, 0.001_R8P]) == 0.5_R8P .and. c%distance([9.9_R8P, 0.001_R8P]) < 0._R8P

! self-intersections
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P], v=[0._R8P, 2._R8P, 0._R8P, 2._R8P], loop_sizes=[4_I4P])
call c%self_intersections(pairs, n)
ok(6) = n == 1
if (ok(6)) ok(6) = all(pairs(:, 1) == [1, 3])
call c%adopt_loops(u=[0._R8P, 2._R8P, 1._R8P, 1._R8P], v=[0._R8P, 0._R8P, 0._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%self_intersections(pairs, n)
ok(6) = ok(6) .and. n == 2
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
call c%self_intersections(pairs, n)
ok(6) = ok(6) .and. n == 0
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
call c%self_intersections(pairs, n)
ok(6) = ok(6) .and. n == 0
print '(A,L1)', 'self-intersections: ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function random_check(cur, lo, hi, label) result(good)
   !< 2000 random points in [lo, hi]: BVH == brute force bit for bit; sign == exact winding.
   type(curve_object), intent(in) :: cur    !< Curve.
   real(R8P),          intent(in) :: lo(2)  !< Box minimum.
   real(R8P),          intent(in) :: hi(2)  !< Box maximum.
   character(*),       intent(in) :: label  !< Label.
   logical                        :: good   !< Result.
   real(R8P)                      :: p(2), db, dt !< Work.
   integer(I4P)                   :: k, mism, wrong, seed_size !< Work.

   call random_seed(size=seed_size)
   call random_seed(put=[(4242 + k, k = 1, seed_size)])
   mism = 0 ; wrong = 0
   do k = 1, 2000
      call random_number(p)
      p = lo + p * (hi - lo)
      call cur%compute_distance(p=p, distance=dt)
      call cur%compute_distance(p=p, distance=db, use_brute_force=.true.)
      if (abs(dt) /= abs(db)) mism = mism + 1
      if (abs(dt) > 1.e-12_R8P .and. ((dt < 0._R8P) .neqv. (cur%winding_number_exact(p) /= 0))) wrong = wrong + 1
   enddo
   good = mism == 0 .and. wrong == 0
   print '(A,A,I0,A,I0)', label, ': BVH/brute mismatches ', mism, '  wrong signs ', wrong
   endfunction random_check

   subroutine naca0012(m, x, y)
   !< NACA 0012 with a closed trailing edge, 2m vertices, cosine spacing, counter-clockwise.
   integer(I4P), intent(in)  :: m       !< Points per side.
   real(R8P),    intent(out) :: x(:), y(:) !< Vertices [2m].
   real(R8P)                 :: xi      !< Abscissa.
   integer(I4P)              :: k, j    !< Counters.

   j = 0
   do k = m, 0, -1 ! upper surface, trailing edge to leading edge
      xi = 0.5_R8P * (1._R8P - cos(PI * k / m))
      j = j + 1 ; x(j) = xi ; y(j) = thickness(xi)
   enddo
   do k = 1, m - 1 ! lower surface, leading edge to trailing edge
      xi = 0.5_R8P * (1._R8P - cos(PI * k / m))
      j = j + 1 ; x(j) = xi ; y(j) = -thickness(xi)
   enddo
   endsubroutine naca0012

   pure function thickness(x) result(y)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge (-0.1036).
   real(R8P), intent(in) :: x !< Abscissa.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_curve_distance
