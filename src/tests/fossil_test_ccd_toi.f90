!< FOSSIL test: time of first contact of two moving surfaces (issue #73, step 4).
!<
!< Unit cubes (src/tests/cube.stl, [0, 1]^3, body origin at a corner), tolerance 1e-6:
!< 1. a cube translating by 3 along x toward a cube 1 away: contact at s = 1/3, toi in [(1 - tol) / 3, 1/3], one step (the
!<    bound is the speed), witnesses tol apart at most;
!< 2. a cube turning a quarter turn about z through its origin toward a post above it (y >= 1.3): first contact when
!<    sin(phi) + cos(phi) = 1.3 (vertex (1, 1) on the post's face), toi within 1e-5 below it;
!< 3. tunnelling: a cube crossing a plate 1e-3 thick within the step; the end-pose check (distance at s = 0 and 1) sees
!<    two separate bodies; time_of_impact reports the contact at (2 - 1e-3) / 5;
!< 4. a near miss (gap 0.1): CCD_SEPARATE, toi = 1; overlapping at s = 0: CCD_INITIAL_CONTACT, toi = 0;
!< 5. a cube turning and translating against another (2 fixed steps; the torus and the larger bodies are in
!<    fossil_val_ccd): d >= tol / 2 at 200 samples of [0, toi] (no earlier contact), d(toi) <= tol, toi within 1e-4 of
!<    the brute-force first contact (first sample with d <= tol, bisected); a step that misses gives CCD_SEPARATE and
!<    the samples agree.
program fossil_test_ccd_toi

use fossil,                only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_identity, pose_interpolate, &
                                  CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT
use penf,                  only : I4P, R8P

implicit none

real(R8P), parameter     :: TOL = 1e-6_R8P, PI = acos(-1._R8P)
type(surface_stl_object) :: a, b
type(rigid_pose_t)       :: a0, a1, b0
real(R8P)                :: toi, d, pa(3), pb(3), sc, d0, d1
integer(I4P)             :: st, it, k
logical                  :: all_ok

all_ok = .true.
b0 = pose_identity()

! 1. translation
call cube(a, [0._R8P, 0._R8P, 0._R8P]) ; call cube(b, [2._R8P, 0._R8P, 0._R8P])
a0 = pose_identity() ; a1 = a0 ; a1%t = [3._R8P, 0._R8P, 0._R8P]
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, point_self=pa, point_other=pb, iterations=it)
call check(st == CCD_CONTACT .and. toi <= 1._R8P / 3 .and. toi >= (1._R8P - TOL) / 3 .and. it == 1 .and. d <= TOL .and. &
           abs(norm2(pb - pa) - d) <= 1e-14_R8P, 'translation: toi in [(1 - tol) / 3, 1/3], one step (toi)', toi)

! 2. rotation onto a post
call cube(b, [-0.2_R8P, 1.3_R8P, 0._R8P])
a1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 2)
sc = (asin(1.3_R8P / sqrt(2._R8P)) - PI / 4) / (PI / 2)
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. toi <= sc .and. sc - toi <= 1e-5_R8P .and. d <= TOL, &
           'rotation onto a post: toi below the analytic angle by < 1e-5 (sc - toi)', sc - toi)
print '(A,I0)', '     steps: ', it

! 3. tunnelling through a thin plate
call cube(b, [0._R8P, 0._R8P, 0._R8P])
call b%resize(x=4._R8P, y=4._R8P, z=1e-3_R8P) ; call b%translate(x=-2._R8P, y=-2._R8P, z=0._R8P) ; call b%analyze
a0 = pose_identity() ; a0%t = [0._R8P, 0._R8P, 2._R8P] ; a1 = a0 ; a1%t = [0._R8P, 0._R8P, -3._R8P]
call a%distance_to(b, d0, pose_self=a0) ; call a%distance_to(b, d1, pose_self=a1)
call check(d0 > 1._R8P .and. d1 > 1._R8P, 'tunnelling: the end poses are apart (min end distance)', min(d0, d1))
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. toi <= (2._R8P - 1e-3_R8P) / 5 .and. toi >= (2._R8P - 1e-3_R8P - TOL) / 5 .and. &
           it == 1, 'tunnelling: contact at (2 - 1e-3) / 5 in one step (toi)', toi)

! 4. near miss and initial contact
call cube(b, [2._R8P, 1.1_R8P, 0._R8P])
a0 = pose_identity() ; a1 = a0 ; a1%t = [4._R8P, 0._R8P, 0._R8P]
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d)
call check(st == CCD_SEPARATE .and. toi == 1._R8P, 'near miss (gap 0.1): separate, toi = 1 (d(1))', d)
call cube(b, [0.5_R8P, 0.5_R8P, 0.5_R8P])
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, iterations=it)
call check(st == CCD_INITIAL_CONTACT .and. toi == 0._R8P .and. it == 0, 'overlapping at s = 0: initial contact', toi)

! 5. cube vs cube against a brute force
call cube(b, [0.6_R8P, -0.5_R8P, -0.5_R8P])
do k = 1, 2
   a0 = pose_from_axis_angle([1._R8P, 0.3_R8P * k, 0.2_R8P], 0.4_R8P * k, [-1.5_R8P, 0.1_R8P * k, 0._R8P])
   a1 = pose_from_axis_angle([0.2_R8P, 1._R8P, -0.3_R8P * k], 1.1_R8P + 0.3_R8P * k, [1.5_R8P, -0.2_R8P * k, 0.1_R8P])
   call brute(k)
enddo
a0 = pose_identity() ; a0%t = [-1.5_R8P, 2._R8P, 0._R8P] ; a1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 2._R8P, &
                                                                                       [1.5_R8P, 2._R8P, 0._R8P])
call brute(4)

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
   !< time_of_impact of a (a0 -> a1) and b (fixed) against 2000 samples and a bisection of the first contact.
   integer(I4P), intent(in) :: case_ !< Case number.
   integer(I4P), parameter  :: NS = 200_I4P
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
   write(label, '(A,I0,A,I0,A)') 'cube vs cube, step ', case_, ' (', it, ' steps): toi vs brute first contact'
   if (first > 1._R8P) then
      call check(st == CCD_SEPARATE .and. nlow == 0, trim(label) // ' (separate)', d)
   else
      call check(st == CCD_CONTACT .and. nlow == 0 .and. d <= TOL .and. abs(toi - first) <= 1e-4_R8P, trim(label), &
                 first - toi)
   endif
   endsubroutine brute

   subroutine cube(s, o)
   !< The unit cube [0, 1]^3 moved by o.
   type(surface_stl_object), intent(inout) :: s    !< Cube.
   real(R8P),                intent(in)    :: o(3) !< Offset.
   type(surface_stl_object)                :: fresh

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%translate(x=o(1), y=o(2), z=o(3)) ; call s%analyze
   endsubroutine cube
endprogram fossil_test_ccd_toi
