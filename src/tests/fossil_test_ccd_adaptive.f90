!< FOSSIL test: the adaptive default of the time of first contact (issue #75).
program fossil_test_ccd_adaptive
!< FOSSIL test: the adaptive default of the time of first contact (issue #75).
!<
!< Absent `directional`, conservative advancement takes plain steps and switches to max(plain, directional) after 8
!< consecutive stalled steps (each advancing s by less than 1 % of the remaining interval and closing less than 1 % of
!< the margin d - tol / 2); `.true.` / `.false.` force either.
!< 1. The rule on an analytic evaluator (constant d = 2 tol, mu = 4, a directional step of 1): default CCD_SEPARATE in
!<    9 steps (8 plain, then one directional), one directional call; `.true.` in 1 step; `.false.` hits the cap
!<    (CCD_NOT_CONVERGED, 1000 steps, no directional call). A linear approach at the bound's speed converges in one
!<    plain step, and with mu = 20 times the speed in about 270 steps, most of them under 1 % of the remaining
!<    interval but each closing 5 % of the margin d - tol / 2: the default never calls the directional step and gives
!<    the plain toi bit for bit.
!< 2. 3-D: a unit cube sliding 4 along x past another, face parallel to face, at a gap of 32 tol: the default reaches
!<    CCD_SEPARATE in at most 8 more steps than `directional=.true.`; a head-on approach (plain converges in one step)
!<    gives the plain toi bit for bit, no pair evaluation.
!< 3. 2-D: the same with unit squares, plain alone (`directional=.false.`) hitting the cap.
!< 4. first_impact_2d: the sliding pair plus a far body: CCD_SEPARATE by default, CCD_NOT_CONVERGED with
!<    `directional=.false.`.
!< The gap sweep, the plain cap in 3-D and the 3-D first_impact are in the validation (fossil_val_ccd).

use fossil,                    only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                                      pose_2d, ccd_advance, first_impact_2d, CCD_CONTACT, CCD_SEPARATE, &
                                      CCD_NOT_CONVERGED
use fossil_ccd_test_distances, only : stepping_distance_t
use penf,                      only : I4P, R8P

implicit none

real(R8P), parameter      :: TOL = 1e-6_R8P
type(stepping_distance_t) :: f
type(surface_stl_object)  :: a, b
type(curve_object)        :: ca, cb, c(3)
type(rigid_pose_t)        :: a0, a1, b0
type(rigid_pose_2d_t)     :: q0, q1, r0, w0(3), w1(3)
real(R8P)                 :: toi, t1, t2
integer(I4P)              :: st, s1, s2, it, i1, i2, np, k, pair(2)
logical                   :: all_ok

all_ok = .true.
b0 = pose_identity() ; r0 = pose_2d(0._R8P)

! 1. the rule
f = stepping_distance_t(e=2._R8P * TOL, dir_step=1._R8P)
call ccd_advance(f, 4._R8P, TOL, toi, st, iterations=it)
call check(st == CCD_SEPARATE .and. toi == 1._R8P .and. it == 9 .and. f%safe_calls == 1, &
           'stall: default separates after 8 plain steps and one directional (steps)', real(it, R8P))
f = stepping_distance_t(e=2._R8P * TOL, dir_step=1._R8P)
call ccd_advance(f, 4._R8P, TOL, toi, st, iterations=it, directional=.true.)
call check(st == CCD_SEPARATE .and. it == 1 .and. f%safe_calls == 1, 'stall: directional=.true. in one step (steps)', &
           real(it, R8P))
f = stepping_distance_t(e=2._R8P * TOL, dir_step=1._R8P)
call ccd_advance(f, 4._R8P, TOL, toi, st, iterations=it, directional=.false.)
call check(st == CCD_NOT_CONVERGED .and. it == 1000 .and. f%safe_calls == 0, &
           'stall: directional=.false. hits the cap, no directional call (toi)', toi)
f = stepping_distance_t(a=2._R8P, c=0.5_R8P, dir_step=1._R8P)
call ccd_advance(f, 2._R8P, TOL, t1, s1, directional=.false.)
f = stepping_distance_t(a=2._R8P, c=0.5_R8P, dir_step=1._R8P)
call ccd_advance(f, 2._R8P, TOL, t2, s2)
call check(s2 == CCD_CONTACT .and. s1 == s2 .and. t1 == t2 .and. f%safe_calls == 0, &
           'no stall: default = plain bit for bit, no directional call (toi)', t2)
f = stepping_distance_t(a=2._R8P, c=0.5_R8P, dir_step=1._R8P)
call ccd_advance(f, 40._R8P, TOL, t1, s1, iterations=i1, directional=.false.)
f = stepping_distance_t(a=2._R8P, c=0.5_R8P, dir_step=1._R8P)
call ccd_advance(f, 40._R8P, TOL, t2, s2, iterations=i2)
call check(s2 == CCD_CONTACT .and. s1 == s2 .and. t1 == t2 .and. i1 == i2 .and. i2 > 100 .and. f%safe_calls == 0, &
           'slow approach (mu = 20 speed): short steps closing the margin, no switch (steps)', real(i2, R8P))

! 2. 3-D sliding past a cube
call cube(a, [0._R8P, 0._R8P, 0._R8P])
a0 = pose_identity() ; a1 = a0 ; a1%t = [4._R8P, 0._R8P, 0._R8P]
call cube(b, [2._R8P, 1._R8P + 32._R8P * TOL, 0._R8P])
call a%time_of_impact(b, a0, a1, b0, b0, t2, s2, tolerance=TOL, iterations=i2, directional=.true.)
call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, iterations=it)
call check(s2 == CCD_SEPARATE .and. st == CCD_SEPARATE .and. it <= i2 + 8, &
           '3-D sliding, gap 32 tol: default separates (default steps)', real(it, R8P))
call cube(b, [2._R8P, 0._R8P, 0._R8P]) ; a1%t = [3._R8P, 0._R8P, 0._R8P]
call a%time_of_impact(b, a0, a1, b0, b0, t1, s1, tolerance=TOL, directional=.false.)
call a%time_of_impact(b, a0, a1, b0, b0, t2, s2, tolerance=TOL, pairs=np)
call check(s2 == CCD_CONTACT .and. s1 == s2 .and. t1 == t2 .and. np == 0, &
           '3-D head-on: default = plain bit for bit, no pair evaluation (toi)', t2)

! 3. 2-D sliding past a square
call rect(ca, 0._R8P, 0._R8P)
q0 = pose_2d(0._R8P) ; q1 = pose_2d(0._R8P, [4._R8P, 0._R8P])
call rect(cb, 2._R8P, 1._R8P + 32._R8P * TOL)
call ca%time_of_impact(cb, q0, q1, r0, r0, t1, s1, tolerance=TOL, iterations=i1, directional=.false.)
call ca%time_of_impact(cb, q0, q1, r0, r0, t2, s2, tolerance=TOL, iterations=i2, directional=.true.)
call ca%time_of_impact(cb, q0, q1, r0, r0, toi, st, tolerance=TOL, iterations=it)
call check(s1 == CCD_NOT_CONVERGED .and. s2 == CCD_SEPARATE .and. st == CCD_SEPARATE .and. it <= i2 + 8, &
           '2-D sliding, gap 32 tol: plain caps, default separates (default steps)', real(it, R8P))
call rect(cb, 2._R8P, 0._R8P) ; q1 = pose_2d(0._R8P, [3._R8P, 0._R8P])
call ca%time_of_impact(cb, q0, q1, r0, r0, t1, s1, tolerance=TOL, directional=.false.)
call ca%time_of_impact(cb, q0, q1, r0, r0, t2, s2, tolerance=TOL, pairs=np)
call check(s2 == CCD_CONTACT .and. s1 == s2 .and. t1 == t2 .and. np == 0, &
           '2-D head-on: default = plain bit for bit, no pair evaluation (toi)', t2)

! 4. several curves
call rect(c(1), 0._R8P, 0._R8P) ; call rect(c(2), 2._R8P, 1._R8P + 16._R8P * TOL) ; call rect(c(3), 50._R8P, 50._R8P)
do k = 1, 3
   w0(k) = pose_2d(0._R8P) ; w1(k) = w0(k)
enddo
w1(1) = pose_2d(0._R8P, [4._R8P, 0._R8P])
call first_impact_2d(c, w0, w1, toi, pair, st, tolerance=TOL, directional=.false.)
call first_impact_2d(c, w0, w1, t2, pair, s2, tolerance=TOL)
call check(st == CCD_NOT_CONVERGED .and. s2 == CCD_SEPARATE .and. t2 == 1._R8P, &
           'first_impact_2d, sliding pair: plain caps, default separates (plain toi)', toi)

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

   subroutine cube(sf, o)
   !< The unit cube [0, 1]^3 moved by o.
   type(surface_stl_object), intent(inout) :: sf   !< Cube.
   real(R8P),                intent(in)    :: o(3) !< Offset.
   type(surface_stl_object)                :: fresh

   sf = fresh
   call sf%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call sf%translate(x=o(1), y=o(2), z=o(3)) ; call sf%analyze
   endsubroutine cube

   subroutine rect(cv, u0, v0)
   !< The unit square [u0, u0 + 1] x [v0, v0 + 1], counter-clockwise.
   type(curve_object), intent(inout) :: cv     !< Curve.
   real(R8P),          intent(in)    :: u0, v0 !< Corner.
   type(curve_object)                :: fresh

   cv = fresh
   call cv%adopt_loops(u=[u0, u0 + 1._R8P, u0 + 1._R8P, u0], v=[v0, v0, v0 + 1._R8P, v0 + 1._R8P], loop_sizes=[4_I4P])
   endsubroutine rect
endprogram fossil_test_ccd_adaptive
