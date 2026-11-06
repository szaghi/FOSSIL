!< FOSSIL test: thickness and gaps of curves, thin regions and tags (issue #62, roadmap #58 WP4, step 4).
!<
!< 1. A NACA 0012 section (100 points per side, cosine spacing, closed trailing edge): the medial thickness of the
!<    segments with 0.1 < x < 0.9 equals the closed form 2 y_t sqrt(1 + y_t'^2) at the midpoint to the discretisation
!<    (measured, stated with a margin); never above the axial ray chord (rays=1), its starting bound.
!< 2. An annulus, radii 1 (256 sides) and 0.5 (128 sides, clockwise): medial thickness 0.5 on the outer segments, gap
!<    2 x 0.5 on the hole's segments, huge on the outer ones (convex); the ray thickness (default fan) is the lower median
!<    over the fan of the exact chord to the inner circle, R cos(theta) - sqrt(r^2 - R^2 sin(theta)^2), to the
!<    discretisation.
!< 3. h decreasing on the section: the thin length never grows, and the trailing-edge segments stay thin to the end.
!< 4. Before build_thin: box_thin_many and tag_block(thin_tags=) give STATUS_INVALID_INPUT. After it, box_thin_many on
!<    32 x 32 cells equals a brute force over every segment (box_segment_overlap_2d, exact), thin_tags equals
!<    cell_bands of it, and box_thin_many_device equals the CPU (host build: the device code path).
program fossil_test_curve_thickness

use fossil,                only : curve_object, STATUS_OK, STATUS_INVALID_INPUT, THICKNESS_MEDIAL
use fossil_box_predicates, only : box_segment_overlap_2d
use fossil_box_query,      only : cell_bands
use penf,                  only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
integer(I4P), parameter   :: NB = 32
type(curve_object)        :: c
real(R8P), allocatable    :: u(:), v(:), th(:), gp(:), tr(:), rl(:), lo(:,:), hi(:,:)
integer(I4P), allocatable :: region(:)
logical                   :: touched(NB * NB), brute(NB * NB), dev(NB * NB)
integer(I4P)              :: tags(NB, NB), ttags(NB, NB), ref(NB, NB, 1)
real(R8P)                 :: err, yt, dyt, xm, h, prev, x(NB + 1), y(NB + 1), cmed
integer(I4P)              :: i, j, k, n, status, invalid, s
logical                   :: ok, all_ok

all_ok = .true.

! ---- 1. NACA 0012
call naca(100_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[size(u, kind=I4P)])
call c%thickness(thickness=th, method=THICKNESS_MEDIAL, invalid=invalid, status=status)
call c%thickness(thickness=tr, rays=1_I4P)
err = 0._R8P ; n = 0
do i = 1, size(u)
   j = mod(i, size(u)) + 1
   xm = 0.5_R8P * (u(i) + u(j))
   if (xm < 0.1_R8P .or. xm > 0.9_R8P) cycle
   call half_thickness(xm, yt, dyt)
   err = max(err, abs(th(i) - 2._R8P * yt * sqrt(1._R8P + dyt ** 2)) / (2._R8P * yt)) ; n = n + 1
enddo
ok = status == STATUS_OK .and. invalid == 0 .and. n > 0 .and. err <= 2e-3_R8P .and. all(th <= tr * (1._R8P + 1e-12_R8P))
print '(A,I0,A,ES10.2,A,L1)', '1. NACA 0012: ', n, ' segments, rel err vs closed form ', err, ', below the chord  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. h sweep
prev = huge(1._R8P) ; h = 0.2_R8P
do while (h > 1e-3_R8P)
   call c%thin_regions(thickness=th, h=h, factor=1._R8P, region=region, region_length=rl)
   ok = ok .and. sum(rl) <= prev .and. region(1) > 0 .and. region(size(u)) > 0   ! the two trailing-edge segments
   prev = sum(rl) ; h = 0.7_R8P * h
enddo
print '(A,I0,A,L1)', '3. h sweep: thin length never grows, trailing edge thin to the end (', size(rl), ' regions left)  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. tagging
x = [(-0.1_R8P + 1.2_R8P * (i - 1) / NB, i = 1, NB + 1)] ; y = [(-0.3_R8P + 0.6_R8P * (i - 1) / NB, i = 1, NB + 1)]
allocate(lo(2, NB * NB), hi(2, NB * NB))
k = 0
do j = 1, NB ; do i = 1, NB
   k = k + 1 ; lo(:, k) = [x(i), y(j)] ; hi(:, k) = [x(i + 1), y(j + 1)]
enddo ; enddo
call c%box_thin_many(lo, hi, 2._R8P, touched, status=status)
ok = status == STATUS_INVALID_INPUT
call c%tag_block(x, y, 1_I4P, tags, status=status, thin_factor=2._R8P, thin_band=1_I4P, thin_tags=ttags)
ok = ok .and. status == STATUS_INVALID_INPUT
call c%build_thin(thickness=th, status=status)
call c%box_thin_many(lo, hi, 2._R8P, touched)
brute = .false.
do k = 1, NB * NB
   do s = 1, size(u)
      if (th(s) >= 2._R8P * maxval(hi(:, k) - lo(:, k))) cycle
      j = mod(s, size(u)) + 1
      if (box_segment_overlap_2d(lo(:, k), hi(:, k), [u(s), v(s)], [u(j), v(j)])) then
         brute(k) = .true. ; exit
      endif
   enddo
enddo
call c%tag_block(x, y, 1_I4P, tags, status=k, thin_factor=2._R8P, thin_band=1_I4P, thin_tags=ttags)
call cell_bands(reshape(brute, [NB, NB, 1]), 1_I4P, ref)
call c%enter_device
call c%box_thin_many_device(lo, hi, 2._R8P, dev)
call c%exit_device
ok = ok .and. status == STATUS_OK .and. k == STATUS_OK .and. all(touched .eqv. brute) .and. count(touched) > 0 .and. &
     all(ttags == ref(:, :, 1)) .and. all(dev .eqv. touched)
print '(A,I0,A,I0,A,L1)', '4. thin cells ', count(touched), ' of ', NB * NB, &
                          ': refused before build_thin, equal to brute force, tags, device  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. annulus
deallocate(u, v) ; allocate(u(384), v(384))
do i = 1, 256
   u(i) = cos(2 * PI * (i - 1) / 256) ; v(i) = sin(2 * PI * (i - 1) / 256)
enddo
do i = 1, 128
   u(256 + i) = 0.5_R8P * cos(-2 * PI * (i - 1) / 128) ; v(256 + i) = 0.5_R8P * sin(-2 * PI * (i - 1) / 128)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[256_I4P, 128_I4P])
call c%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL, invalid=invalid)
call c%thickness(thickness=tr)
cmed = fan_median_chord(1._R8P, 0.5_R8P)
ok = invalid == 0 .and. maxval(abs(th(1:256) - 0.5_R8P)) <= 1e-3_R8P .and. all(gp(1:256) == huge(1._R8P)) .and. &
     maxval(abs(gp(257:) - 1._R8P)) <= 1e-3_R8P .and. maxval(abs(tr(1:256) - cmed)) <= 2e-3_R8P
print '(A,3ES10.2,A,L1)', '2. annulus: thickness, hole gap, ray errors ', maxval(abs(th(1:256) - 0.5_R8P)), &
      maxval(abs(gp(257:) - 1._R8P)), maxval(abs(tr(1:256) - cmed)), '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   pure subroutine half_thickness(xx, yy, dy)
   !< NACA 0012 half thickness (closed trailing edge) and its slope.
   real(R8P), intent(in)  :: xx
   real(R8P), intent(out) :: yy, dy

   yy = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.1260_R8P * xx - 0.3516_R8P * xx ** 2 + 0.2843_R8P * xx ** 3 - 0.1036_R8P * xx ** 4)
   dy = 0.6_R8P * (0.14845_R8P / sqrt(xx) - 0.1260_R8P - 0.7032_R8P * xx + 0.8529_R8P * xx ** 2 - 0.4144_R8P * xx ** 3)
   endsubroutine half_thickness

   function fan_median_chord(ro, ri) result(cm)
   !< Lower median, over the default fan (60 degrees, 31 rays), of the chord from the circle of radius ro inwards to the
   !< circle of radius ri (the ray grazes it at 30 degrees when ro sin 30 = ri).
   real(R8P), intent(in) :: ro, ri
   real(R8P)             :: cm, a(31), key, t
   integer(I4P)          :: kk, jj

   do kk = 1, 31
      t = 0._R8P ; if (kk > 1) t = 0.5_R8P * 60._R8P * PI / 180._R8P * real(kk / 2, R8P) / 15._R8P
      a(kk) = ro * cos(t) - sqrt(max(0._R8P, ri ** 2 - (ro * sin(t)) ** 2))
   enddo
   do kk = 2, 31
      key = a(kk) ; jj = kk - 1
      do while (jj >= 1)
         if (a(jj) <= key) exit
         a(jj + 1) = a(jj) ; jj = jj - 1
      enddo
      a(jj + 1) = key
   enddo
   cm = a(16)
   endfunction fan_median_chord

   subroutine naca(np, uu, vv)
   !< NACA 0012 of unit chord, np points per side (cosine spacing, closed trailing edge), counter-clockwise from the
   !< trailing edge (vertex 1).
   integer(I4P),           intent(in)  :: np
   real(R8P), allocatable, intent(out) :: uu(:), vv(:)
   real(R8P)                           :: xs(np + 1), ts(np + 1), dd
   integer(I4P)                        :: q, nn

   xs = [(0.5_R8P * (1._R8P - cos(PI * q / np)), q = 0, np)]
   do q = 1, np + 1
      call half_thickness(xs(q), ts(q), dd)
   enddo
   ts(np + 1) = 0._R8P
   allocate(uu(2 * np), vv(2 * np))
   nn = 0
   do q = np + 1, 2, -1
      nn = nn + 1 ; uu(nn) = xs(q) ; vv(nn) = ts(q)
   enddo
   do q = 1, np
      nn = nn + 1 ; uu(nn) = xs(q) ; vv(nn) = -ts(q)
   enddo
   endsubroutine naca
endprogram fossil_test_curve_thickness
