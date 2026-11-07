!< FOSSIL test: the earliest contact among several moving bodies (issue #73, step 6).
!<
!< Four unit cubes (squares in 2-D), tolerance 1e-6: cube 1 slides +x by 3 toward cube 2 (contact at 1/3), cube 3 drops
!< -y by 4 onto cube 2 (at 1/2, and onto cube 1's path at 1/2), cube 4 is far away:
!< 1. first_impact gives pair (1, 2) at toi in [(1 - tol) / 3, 1/3], equal bit for bit (toi, pair, status) to the loop of
!<    time_of_impact over all 6 pairs with the same tie rule; the far cube is never tested;
!< 2. cube 1 at rest: the earliest contact is (2, 3) at 1/2, as the loop;
!< 3. nothing moves toward anything: CCD_SEPARATE, toi = 1, pair 0; two overlapping cubes: CCD_INITIAL_CONTACT, toi = 0;
!<    a size mismatch: CCD_INVALID_INPUT;
!< 4. the 2-D twin with squares.
program fossil_test_ccd_bodies

use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_2d, first_impact, &
                   first_impact_2d, CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT, CCD_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter     :: TOL = 1e-6_R8P
type(surface_stl_object) :: s(4)
type(curve_object)       :: c(4)
type(rigid_pose_t)       :: p0(4), p1(4)
type(rigid_pose_2d_t)    :: q0(4), q1(4)
real(R8P)                :: toi, tb
integer(I4P)             :: pair(2), pb(2), st, sb, nt, k
logical                  :: all_ok

all_ok = .true.

! 1. four cubes
call cube(s(1), [0._R8P, 0._R8P, 0._R8P]) ; call cube(s(2), [2._R8P, 0._R8P, 0._R8P])
call cube(s(3), [2._R8P, 3._R8P, 0._R8P]) ; call cube(s(4), [50._R8P, 50._R8P, 50._R8P])
do k = 1, 4
   p0(k) = pose_identity() ; p1(k) = p0(k)
enddo
p1(1)%t = [3._R8P, 0._R8P, 0._R8P] ; p1(3)%t = [0._R8P, -4._R8P, 0._R8P]
call first_impact(s, p0, p1, toi, pair, st, tolerance=TOL, tested=nt)
call brute_3d(4_I4P)
call check(st == CCD_CONTACT .and. all(pair == [1, 2]) .and. toi <= 1._R8P / 3 .and. toi >= (1._R8P - TOL) / 3, &
           'four cubes: pair (1, 2) at 1/3 (toi)', toi)
call check(toi == tb .and. all(pair == pb) .and. st == sb, 'four cubes: equal to the loop over all pairs, bit for bit', toi - tb)
call check(nt < 6, 'four cubes: pairs tested (of 6)', real(nt, R8P))

! 2. cube 1 at rest
p1(1) = p0(1)
call first_impact(s, p0, p1, toi, pair, st, tolerance=TOL, tested=nt)
call brute_3d(4_I4P)
call check(st == CCD_CONTACT .and. all(pair == [2, 3]) .and. all(pair == pb) .and. toi == tb .and. st == sb .and. &
           abs(toi - 0.5_R8P) <= 1e-6_R8P, 'cube 1 at rest: (2, 3) at 1/2, as the loop (toi)', toi)

! 3. separate, initial contact, invalid input
p1(3) = p0(3)
call first_impact(s, p0, p1, toi, pair, st, tolerance=TOL)
call check(st == CCD_SEPARATE .and. toi == 1._R8P .and. all(pair == 0), 'nothing moves: separate, toi = 1', toi)
call cube(s(4), [2.5_R8P, 0.5_R8P, 0.5_R8P])
call first_impact(s, p0, p1, toi, pair, st, tolerance=TOL)
call check(st == CCD_INITIAL_CONTACT .and. toi == 0._R8P .and. all(pair == [2, 4]), 'overlapping: initial contact (2, 4)', &
           toi)
call first_impact(s, p0(1:3), p1, toi, pair, st)
call check(st == CCD_INVALID_INPUT, 'size mismatch: invalid input', real(st, R8P))

! 4. 2-D
call rect(c(1), 0._R8P, 0._R8P) ; call rect(c(2), 2._R8P, 0._R8P) ; call rect(c(3), 2._R8P, 3._R8P)
call rect(c(4), 50._R8P, 50._R8P)
do k = 1, 4
   q0(k) = pose_2d(0._R8P) ; q1(k) = q0(k)
enddo
q1(1) = pose_2d(0._R8P, [3._R8P, 0._R8P]) ; q1(3) = pose_2d(0._R8P, [0._R8P, -4._R8P])
call first_impact_2d(c, q0, q1, toi, pair, st, tolerance=TOL, tested=nt)
call brute_2d(4_I4P)
call check(st == CCD_CONTACT .and. all(pair == [1, 2]) .and. toi <= 1._R8P / 3 .and. toi >= (1._R8P - TOL) / 3 .and. &
           toi == tb .and. all(pair == pb) .and. st == sb .and. nt < 6, '2-D: pair (1, 2) at 1/3, equal to the loop (tested)', &
           real(nt, R8P))

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

   subroutine brute_3d(n)
   !< The loop over all pairs: the earliest contact, ties to the smallest pair.
   integer(I4P), intent(in) :: n !< Bodies.
   integer(I4P)             :: i, j, stt
   real(R8P)                :: t

   tb = 1._R8P ; pb = 0 ; sb = CCD_SEPARATE
   do i = 1, n
      do j = i + 1, n
         call s(i)%time_of_impact(s(j), p0(i), p1(i), p0(j), p1(j), t, stt, tolerance=TOL)
         if (stt == CCD_SEPARATE) cycle
         if (t < tb .or. (t == tb .and. sb == CCD_SEPARATE)) then
            tb = t ; pb = [i, j] ; sb = stt
         endif
      enddo
   enddo
   endsubroutine brute_3d

   subroutine brute_2d(n)
   !< 2-D twin of brute_3d.
   integer(I4P), intent(in) :: n !< Bodies.
   integer(I4P)             :: i, j, stt
   real(R8P)                :: t

   tb = 1._R8P ; pb = 0 ; sb = CCD_SEPARATE
   do i = 1, n
      do j = i + 1, n
         call c(i)%time_of_impact(c(j), q0(i), q1(i), q0(j), q1(j), t, stt, tolerance=TOL)
         if (stt == CCD_SEPARATE) cycle
         if (t < tb .or. (t == tb .and. sb == CCD_SEPARATE)) then
            tb = t ; pb = [i, j] ; sb = stt
         endif
      enddo
   enddo
   endsubroutine brute_2d

   subroutine cube(sf, o)
   !< The unit cube [0, 1]^3 moved by o.
   type(surface_stl_object), intent(inout) :: sf   !< Cube.
   real(R8P),                intent(in)    :: o(3) !< Offset.
   type(surface_stl_object)                :: fresh

   sf = fresh
   call sf%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call sf%translate(x=o(1), y=o(2), z=o(3)) ; call sf%analyze
   endsubroutine cube

   subroutine rect(cc, u0, v0)
   !< The unit square at (u0, v0), counter-clockwise.
   type(curve_object), intent(inout) :: cc     !< Curve.
   real(R8P),          intent(in)    :: u0, v0 !< Corner.
   type(curve_object)                :: fresh

   cc = fresh
   call cc%adopt_loops(u=[u0, u0 + 1._R8P, u0 + 1._R8P, u0], v=[v0, v0, v0 + 1._R8P, v0 + 1._R8P], loop_sizes=[4_I4P])
   endsubroutine rect
endprogram fossil_test_ccd_bodies
