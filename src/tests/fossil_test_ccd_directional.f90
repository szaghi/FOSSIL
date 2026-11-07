!< FOSSIL test: the directional (C2A) step of the time of first contact (issue #73, step 5).
!<
!< 1. 3-D (unit cubes): head-on translation, an oblique approach (45 degrees to a face), a quarter turn onto a post and
!<    two turning-and-translating steps; with directional=.true. every toi keeps d >= tol / 2 at 200 samples of
!<    [0, toi] and d(toi) <= tol, never takes more steps than the plain bound, and the plain / directional steps and the
!<    exact pair evaluations are printed (measured, not asserted);
!< 2. 2-D: the same cases with squares;
!< 3. negative control (2-D): an L-shaped body (a bar 0.5 above a square, an arm 0.7 to its left) sliding right: the
!<    whole-body separating-plane bound (witness direction vertical, normal to the motion) never closes and misses the
!<    arm's contact at s = 0.35; the per-pair bound finds it.
program fossil_test_ccd_directional

use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, &
                   pose_2d, pose_interpolate, CCD_CONTACT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter     :: TOL = 1e-6_R8P, PI = acos(-1._R8P)
type(surface_stl_object) :: a, b
type(curve_object)       :: ca, cb
type(rigid_pose_t)       :: a0, a1, b0
type(rigid_pose_2d_t)    :: q0, q1, r0
real(R8P)                :: toi, d, pw(2), qw(2), s, dd, n(2), stp
integer(I4P)             :: st, k
logical                  :: all_ok

all_ok = .true.
b0 = pose_identity() ; r0 = pose_2d(0._R8P)
print '(A)', '  case                         plain steps  directional steps  pairs'

! 1. 3-D
call cube(a, [0._R8P, 0._R8P, 0._R8P]) ; call cube(b, [2._R8P, 0._R8P, 0._R8P])
a0 = pose_identity() ; a1 = a0 ; a1%t = [3._R8P, 0._R8P, 0._R8P]
call case_3d('3-D head-on')
a1%t = [3._R8P, 3._R8P, 0._R8P]
call cube(b, [0._R8P, 0._R8P, 0._R8P]) ; call b%resize(y=11._R8P) ; call b%translate(x=2._R8P, y=-5._R8P) ; call b%analyze
call case_3d('3-D oblique 45')
call cube(b, [-0.2_R8P, 1.3_R8P, 0._R8P]) ; a1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 2)
call case_3d('3-D quarter turn onto post')
call cube(b, [0.6_R8P, -0.5_R8P, -0.5_R8P])
do k = 1, 2
   a0 = pose_from_axis_angle([1._R8P, 0.3_R8P * k, 0.2_R8P], 0.4_R8P * k, [-1.5_R8P, 0.1_R8P * k, 0._R8P])
   a1 = pose_from_axis_angle([0.2_R8P, 1._R8P, -0.3_R8P * k], 1.1_R8P + 0.3_R8P * k, [1.5_R8P, -0.2_R8P * k, 0.1_R8P])
   call case_3d('3-D turn and translate')
enddo

! 2. 2-D
call rect(ca, 0._R8P, 0._R8P, 1._R8P, 1._R8P) ; call rect(cb, 2._R8P, 0._R8P, 1._R8P, 1._R8P)
q0 = pose_2d(0._R8P) ; q1 = pose_2d(0._R8P, [3._R8P, 0._R8P])
call case_2d('2-D head-on')
call rect(cb, 2._R8P, -5._R8P, 1._R8P, 11._R8P) ; q1 = pose_2d(0._R8P, [3._R8P, 3._R8P])
call case_2d('2-D oblique 45')
call rect(cb, -0.2_R8P, 1.3_R8P, 1._R8P, 1._R8P) ; q1 = pose_2d(PI / 2)
call case_2d('2-D quarter turn onto post')

! 3. negative control: the L-shaped body
call rect(ca, 0._R8P, 0._R8P, 1._R8P, 1._R8P)
call cb%adopt_loops(u=[-0.8_R8P, -0.7_R8P, -0.7_R8P, 2._R8P, 2._R8P, -0.8_R8P], &
                    v=[-0.5_R8P, -0.5_R8P, 1.5_R8P, 1.5_R8P, 1.6_R8P, 1.6_R8P], loop_sizes=[6_I4P])
q0 = pose_2d(0._R8P) ; q1 = pose_2d(0._R8P, [2._R8P, 0._R8P])
call ca%distance_to(cb, dd, point_self=pw, point_other=qw, pose_other=q0)
n = (qw - pw) / dd
stp = huge(1._R8P)
if (dot_product(n, [0._R8P, 0._R8P] - [2._R8P, 0._R8P]) > 0._R8P) &
   stp = (dd - 0.5_R8P * TOL) / dot_product(n, [0._R8P, 0._R8P] - [2._R8P, 0._R8P])
call check(abs(n(2) - 1._R8P) <= 1e-14_R8P .and. stp > 0.35_R8P, &
           'negative control: the whole-body plane bound steps past the contact at 0.35 (its step)', min(stp, 9._R8P))
call ca%time_of_impact(cb, r0, r0, q0, q1, toi, st, tolerance=TOL, directional=.true.)
call check(st == CCD_CONTACT .and. toi <= 0.35_R8P .and. toi >= 0.35_R8P - 1e-6_R8P, &
           'per-pair bound: contact of the arm at 0.35 (toi)', toi)

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

   subroutine case_3d(label)
   !< Plain vs directional on a (a0 -> a1) and b (fixed); the directional result certified by 200 samples.
   character(*), intent(in) :: label !< Case.
   real(R8P)                :: t1, t2, dt
   integer(I4P)             :: s1, s2, i1, i2, np, j, nlow

   call a%time_of_impact(b, a0, a1, b0, b0, t1, s1, tolerance=TOL, iterations=i1)
   call a%time_of_impact(b, a0, a1, b0, b0, t2, s2, tolerance=TOL, distance=d, iterations=i2, directional=.true., pairs=np)
   nlow = 0
   do j = 0, 199
      call a%distance_to(b, dt, pose_self=pose_interpolate(a0, a1, t2 * j / 200._R8P))
      if (dt < 0.5_R8P * TOL * (1._R8P - 1e-9_R8P)) nlow = nlow + 1
   enddo
   print '(2X,A28,I12,I19,I7)', label, i1, i2, np
   call check(s1 == s2 .and. i2 <= i1 .and. nlow == 0 .and. (s2 /= CCD_CONTACT .or. d <= TOL), &
              label // ': certified, no more steps (plain toi - directional toi)', t1 - t2)
   endsubroutine case_3d

   subroutine case_2d(label)
   !< 2-D twin of case_3d.
   character(*), intent(in) :: label !< Case.
   real(R8P)                :: t1, t2, dt
   integer(I4P)             :: s1, s2, i1, i2, np, j, nlow

   call ca%time_of_impact(cb, q0, q1, r0, r0, t1, s1, tolerance=TOL, iterations=i1)
   call ca%time_of_impact(cb, q0, q1, r0, r0, t2, s2, tolerance=TOL, distance=d, iterations=i2, directional=.true., &
                          pairs=np)
   nlow = 0
   do j = 0, 199
      call ca%distance_to(cb, dt, pose_self=pose_interpolate(q0, q1, t2 * j / 200._R8P))
      if (dt < 0.5_R8P * TOL * (1._R8P - 1e-9_R8P)) nlow = nlow + 1
   enddo
   print '(2X,A28,I12,I19,I7)', label, i1, i2, np
   call check(s1 == s2 .and. i2 <= i1 .and. nlow == 0 .and. (s2 /= CCD_CONTACT .or. d <= TOL), &
              label // ': certified, no more steps (plain toi - directional toi)', t1 - t2)
   endsubroutine case_2d

   subroutine cube(sf, o)
   !< The unit cube [0, 1]^3 moved by o.
   type(surface_stl_object), intent(inout) :: sf   !< Cube.
   real(R8P),                intent(in)    :: o(3) !< Offset.
   type(surface_stl_object)                :: fresh

   sf = fresh
   call sf%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call sf%translate(x=o(1), y=o(2), z=o(3)) ; call sf%analyze
   endsubroutine cube

   subroutine rect(c, u0, v0, w, h)
   !< The rectangle [u0, u0 + w] x [v0, v0 + h], counter-clockwise.
   type(curve_object), intent(inout) :: c            !< Curve.
   real(R8P),          intent(in)    :: u0, v0, w, h !< Corner and sides.
   type(curve_object)                :: fresh

   c = fresh
   call c%adopt_loops(u=[u0, u0 + w, u0 + w, u0], v=[v0, v0, v0 + h, v0 + h], loop_sizes=[4_I4P])
   endsubroutine rect
endprogram fossil_test_ccd_directional
