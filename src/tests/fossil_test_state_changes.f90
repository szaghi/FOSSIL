!< FOSSIL test: fresh and dead points between two poses (issue #41 step 4).
program fossil_test_state_changes
!< FOSSIL test: fresh and dead points between two poses (issue #41 step 4).
!<
!< - The unit cube translated, and rotated about an off-centre axis: on a 21^3 lattice of points the
!<   flags equal the analytic change of membership (body point in (0, 1)^3 at each pose), at every
!<   point farther than 1e-9 from the faces at both poses; pseudo-normal and fast winding number.
!< - The band filter never misses a change: the flags equal those from both queries at every point,
!<   and fewer than all points need the second query.
!< - A square translated and a NACA 0012 pitching about its quarter chord: the flags equal the change
!<   of the exact winding number at the body points.

use fossil,                   only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                                     pose_from_axis_angle, pose_2d, to_body, SIGN_FAST_WINDING_NUMBER, SIGN_PSEUDO_NORMAL
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
type(rigid_pose_t)        :: p0, p1
real(R8P)                 :: u(200), v(200), r(3, 3), ctr(3)
integer(I4P)              :: k
logical                   :: ok(3)

call bench_bodies(b)
call to_surface(b(1)%soup, s)
p0 = pose_identity()
p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0.13_R8P, 0.07_R8P, -0.05_R8P])
ok(1) = cube(p0, p1, 'cube translated')
ctr = [2._R8P, 0.5_R8P, 0.5_R8P]
p1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.3_R8P)
r = p1%r ; p1%t = ctr - matmul(r, ctr)
ok(1) = ok(1) .and. cube(p0, p1, 'cube rotated about an off-centre axis')
p0 = p1 ; p1 = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], -0.2_R8P, [0.05_R8P, 0._R8P, 0.1_R8P])
ok(1) = ok(1) .and. cube(p0, p1, 'cube between two general poses')
! curves
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
ok(2) = curve(pose_2d(0._R8P), pose_2d(0._R8P, [0.21_R8P, -0.13_R8P]), 'square translated')
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(3) = curve(pose_2d(0._R8P), pitch(0.15_R8P), 'NACA 0012 pitching')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function cube(q0, q1, label) result(good)
   !< Flags on a lattice against the analytic membership change of the unit cube.
   type(rigid_pose_t), intent(in) :: q0, q1 !< Poses.
   character(*),       intent(in) :: label  !< Label.
   logical                        :: good   !< Result.
   integer, parameter :: M = 21 !< Lattice.
   type(vector_R8P)   :: q(M**3) !< Points.
   logical            :: fr(M**3), de(M**3), fr2(M**3), de2(M**3), a0, a1, near !< Work.
   real(R8P)          :: x(3), c0(3), c1(3) !< Work.
   integer(I4P)       :: i, j, kk, m_, nbad, nchk, nchg, a !< Work.
   integer(I4P), parameter :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.

   m_ = 0
   do kk = 1, M
      do j = 1, M
         do i = 1, M
            m_ = m_ + 1
            x = -0.3_R8P + 1.6_R8P * ([i, j, kk] - 1) / real(M - 1, R8P) + 0.0137_R8P
            q(m_) = vector_R8P(x(1), x(2), x(3))
         enddo
      enddo
   enddo
   good = .true.
   do a = 1, 2
      call s%state_changes(pose0=q0, pose1=q1, points=q, fresh=fr, dead=de, sign_algorithm=ALGO(a), checked=nchk)
      nbad = 0 ; nchg = 0
      do i = 1, M**3
         c0 = to_body(q0, [q(i)%x, q(i)%y, q(i)%z]) ; c1 = to_body(q1, [q(i)%x, q(i)%y, q(i)%z])
         near = minval(min(abs(c0), abs(c0 - 1._R8P))) < 1.e-9_R8P .or. minval(min(abs(c1), abs(c1 - 1._R8P))) < 1.e-9_R8P
         if (near) cycle
         a0 = all(c0 > 0._R8P .and. c0 < 1._R8P) ; a1 = all(c1 > 0._R8P .and. c1 < 1._R8P)
         if ((fr(i) .neqv. (a0 .and. .not. a1)) .or. (de(i) .neqv. (a1 .and. .not. a0))) nbad = nbad + 1
         if (a0 .neqv. a1) nchg = nchg + 1
         ! both queries at every point
         fr2(i) = s%is_point_inside(q(i), sign_algorithm=ALGO(a), pose=q0) .and. &
                  .not. s%is_point_inside(q(i), sign_algorithm=ALGO(a), pose=q1)
         de2(i) = s%is_point_inside(q(i), sign_algorithm=ALGO(a), pose=q1) .and. &
                  .not. s%is_point_inside(q(i), sign_algorithm=ALGO(a), pose=q0)
         if ((fr2(i) .neqv. fr(i)) .or. (de2(i) .neqv. de(i))) nbad = nbad + 1
      enddo
      good = good .and. nbad == 0 .and. nchk < M**3 .and. nchg > 0
      print '(A,A,I0,A,I0,A,I0,A,I0)', label, ': sign ', ALGO(a), '  changed ', nchg, '  second queries ', nchk, '  errors ', nbad
   enddo
   endfunction cube

   function curve(q0, q1, label) result(good)
   !< Flags on a lattice against the change of the exact winding number at the body points.
   type(rigid_pose_2d_t), intent(in) :: q0, q1 !< Poses.
   character(*),          intent(in) :: label  !< Label.
   logical                           :: good   !< Result.
   integer, parameter :: M = 61 !< Lattice.
   real(R8P)          :: q(2, M * M), lo(2), hi(2), b0(2), b1(2) !< Work.
   logical            :: fr(M * M), de(M * M), a0, a1 !< Work.
   integer(I4P)       :: i, j, m_, nbad, nchk, nchg !< Work.

   lo = c%get_bmin() - 0.3_R8P ; hi = c%get_bmax() + 0.3_R8P
   m_ = 0
   do j = 1, M
      do i = 1, M
         m_ = m_ + 1 ; q(:, m_) = lo + (hi - lo) * ([i, j] - 1) / real(M - 1, R8P) + 0.0013_R8P
      enddo
   enddo
   call c%state_changes(pose0=q0, pose1=q1, points=q, fresh=fr, dead=de, checked=nchk)
   nbad = 0 ; nchg = 0
   do i = 1, M * M
      b0 = to_body(q0, q(:, i)) ; b1 = to_body(q1, q(:, i))
      if (abs(c%distance(b0)) < 1.e-9_R8P .or. abs(c%distance(b1)) < 1.e-9_R8P) cycle
      a0 = c%winding_number_exact(b0) /= 0 ; a1 = c%winding_number_exact(b1) /= 0
      if ((fr(i) .neqv. (a0 .and. .not. a1)) .or. (de(i) .neqv. (a1 .and. .not. a0))) nbad = nbad + 1
      if (a0 .neqv. a1) nchg = nchg + 1
   enddo
   good = nbad == 0 .and. nchk < M * M .and. nchg > 0
   print '(A,A,I0,A,I0,A,I0)', label, ': changed ', nchg, '  second queries ', nchk, '  errors ', nbad
   endfunction curve

   pure function pitch(angle) result(p)
   !< Rotation by angle about the quarter chord (0.25, 0).
   real(R8P), intent(in) :: angle !< Angle.
   type(rigid_pose_2d_t) :: p     !< Pose.

   p = pose_2d(angle)
   p%t = [0.25_R8P, 0._R8P] - matmul(p%r, [0.25_R8P, 0._R8P])
   endfunction pitch

   pure function thickness(xx) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xx !< Chord position.
   real(R8P)             :: y  !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.126_R8P * xx - 0.3516_R8P * xx**2 + 0.2843_R8P * xx**3 - 0.1036_R8P * xx**4)
   endfunction thickness
endprogram fossil_test_state_changes
