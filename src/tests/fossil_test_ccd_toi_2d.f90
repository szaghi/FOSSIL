!< FOSSIL test: time of first contact of two moving curves (issue #73, step 4).
!<
!< Unit squares ([0, 1]^2, body origin at a corner), tolerance 1e-6:
!< 1. a square translating by 3 along u toward a square 1 away: toi in [(1 - tol) / 3, 1/3], one step;
!< 2. a square turning a quarter turn about its origin toward a post above it (v >= 1.3): first contact when
!<    sin(phi) + cos(phi) = 1.3, toi within 1e-5 below it;
!< 3. tunnelling: a needle 1e-3 wide crossing a square within the step (apart at s = 0 and 1): contact at
!<    (1 - 1e-3) / 4 in one step;
!< 4. a near miss: CCD_SEPARATE; overlapping at s = 0: CCD_INITIAL_CONTACT;
!< 5. a circle (64 segments) turning and translating against a square: d >= tol / 2 at 400 samples of [0, toi], toi
!<    within 1e-4 of the brute-force first contact; a missing step gives CCD_SEPARATE.
program fossil_test_ccd_toi_2d

use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, pose_interpolate, CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter  :: TOL = 1e-6_R8P, PI = acos(-1._R8P)
type(curve_object)    :: a, b
type(rigid_pose_2d_t) :: a0, a1, b0
real(R8P)             :: toi, d, pa(2), pb(2), sc, d0, d1, cu(64), cv(64)
integer(I4P)          :: st, it, k
logical               :: all_ok

all_ok = .true.
b0 = pose_2d(0._R8P)

! 1. translation
call rect(a, 0._R8P, 0._R8P, 1._R8P, 1._R8P) ; call rect(b, 2._R8P, 0._R8P, 1._R8P, 1._R8P)
a0 = pose_2d(0._R8P) ; a1 = pose_2d(0._R8P, [3._R8P, 0._R8P])
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, point_self=pa, point_other=pb, iterations=it)
call check(st == CCD_CONTACT .and. toi <= 1._R8P / 3 .and. toi >= (1._R8P - TOL) / 3 .and. it == 1 .and. d <= TOL .and. &
           abs(norm2(pb - pa) - d) <= 1e-14_R8P, 'translation: toi in [(1 - tol) / 3, 1/3], one step (toi)', toi)

! 2. rotation onto a post
call rect(b, -0.2_R8P, 1.3_R8P, 1._R8P, 1._R8P)
a1 = pose_2d(PI / 2)
sc = (asin(1.3_R8P / sqrt(2._R8P)) - PI / 4) / (PI / 2)
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. toi <= sc .and. sc - toi <= 1e-5_R8P .and. d <= TOL, &
           'rotation onto a post: toi below the analytic angle by < 1e-5 (sc - toi)', sc - toi)
print '(A,I0)', '     steps: ', it

! 3. tunnelling: a needle 1e-3 wide
call rect(a, 0._R8P, 0._R8P, 1e-3_R8P, 1._R8P) ; call rect(b, 0._R8P, 0._R8P, 1._R8P, 1._R8P)
a0 = pose_2d(0._R8P, [-1._R8P, 0._R8P]) ; a1 = pose_2d(0._R8P, [3._R8P, 0._R8P])
call a%distance_to(b, d0, pose_self=a0) ; call a%distance_to(b, d1, pose_self=a1)
call check(d0 > 0.9_R8P .and. d1 > 0.9_R8P, 'tunnelling: the end poses are apart (min end distance)', min(d0, d1))
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, iterations=it)
call check(st == CCD_CONTACT .and. toi <= (1._R8P - 1e-3_R8P) / 4 .and. toi >= (1._R8P - 1e-3_R8P - TOL) / 4 .and. &
           it == 1, 'tunnelling: contact at (1 - 1e-3) / 4 in one step (toi)', toi)

! 4. near miss and initial contact
call rect(a, 0._R8P, 0._R8P, 1._R8P, 1._R8P) ; call rect(b, 2._R8P, 1.1_R8P, 1._R8P, 1._R8P)
a0 = pose_2d(0._R8P) ; a1 = pose_2d(0._R8P, [4._R8P, 0._R8P])
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d)
call check(st == CCD_SEPARATE .and. toi == 1._R8P, 'near miss (gap 0.1): separate, toi = 1 (d(1))', d)
call rect(b, 0.5_R8P, 0.5_R8P, 1._R8P, 1._R8P)
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, iterations=it)
call check(st == CCD_INITIAL_CONTACT .and. toi == 0._R8P .and. it == 0, 'overlapping at s = 0: initial contact', toi)

! 5. circle vs square against a brute force
do k = 1, 64
   cu(k) = 0.5_R8P * cos(2 * PI * (k - 1) / 64._R8P) + 0.3_R8P ; cv(k) = 0.5_R8P * sin(2 * PI * (k - 1) / 64._R8P)
enddo
call a%adopt_loops(u=cu, v=cv, loop_sizes=[64_I4P]) ; call rect(b, 0.8_R8P, -0.5_R8P, 1._R8P, 1._R8P)
do k = 1, 2
   a0 = pose_2d(0.4_R8P * k, [-2._R8P, 0.1_R8P * k]) ; a1 = pose_2d(-1.5_R8P * k, [1.5_R8P, -0.2_R8P * k])
   call brute(k)
enddo
a0 = pose_2d(0._R8P, [-2._R8P, 2._R8P]) ; a1 = pose_2d(2._R8P, [1.5_R8P, 2._R8P])
call brute(3)

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

   subroutine brute(case_)
   !< time_of_impact of a (a0 -> a1) and b (fixed) against 400 samples and a bisection of the first contact.
   integer(I4P), intent(in) :: case_ !< Case number.
   integer(I4P), parameter  :: NS = 400_I4P
   real(R8P)                :: t, dt, lo, hi, mid, dm, first
   integer(I4P)             :: j, nlow, n
   character(64)            :: label

   call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, iterations=it)
   nlow = 0 ; first = 2._R8P
   do j = 0, NS
      t = real(j, R8P) / NS
      call a%distance_to(b, dt, pose_self=pose_interpolate(a0, a1, t))
      if (t < toi .and. dt < 0.5_R8P * TOL * (1._R8P - 1e-9_R8P)) nlow = nlow + 1
      if (dt <= TOL .and. first > 1._R8P) then
         lo = max(0._R8P, t - 1._R8P / NS) ; hi = t
         do n = 1, 60
            mid = 0.5_R8P * (lo + hi)
            call a%distance_to(b, dm, pose_self=pose_interpolate(a0, a1, mid))
            if (dm <= TOL) then
               hi = mid
            else
               lo = mid
            endif
         enddo
         first = hi
      endif
   enddo
   write(label, '(A,I0,A,I0,A)') 'circle vs square, step ', case_, ' (', it, ' steps): toi vs brute first contact'
   if (first > 1._R8P) then
      call check(st == CCD_SEPARATE .and. nlow == 0, trim(label) // ' (separate)', d)
   else
      call check(st == CCD_CONTACT .and. nlow == 0 .and. d <= TOL .and. abs(toi - first) <= 1e-4_R8P, trim(label), &
                 first - toi)
   endif
   endsubroutine brute

   subroutine rect(c, u0, v0, w, h)
   !< The rectangle [u0, u0 + w] x [v0, v0 + h], counter-clockwise.
   type(curve_object), intent(inout) :: c          !< Curve.
   real(R8P),          intent(in)    :: u0, v0, w, h !< Corner and sides.
   type(curve_object)                :: fresh

   c = fresh
   call c%adopt_loops(u=[u0, u0 + w, u0 + w, u0], v=[v0, v0, v0 + h, v0 + h], loop_sizes=[4_I4P])
   endsubroutine rect
endprogram fossil_test_ccd_toi_2d
