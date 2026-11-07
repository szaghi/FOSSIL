!< FOSSIL test: motion bounds of `pose_interpolate` (issue #73, step 2).
program fossil_test_ccd_bounds
!< FOSSIL test: motion bounds of `pose_interpolate` (issue #73, step 2).
!<
!< A box of 8 vertices away from its frame origin (2-D: a quadrilateral):
!< 1. a pure translation: mu = |t1 - t0|;
!< 2. a pure rotation about z: mu = theta max rho, rho the distance to the axis through the body origin;
!< 3. rotation + translation with general axes: over 10^4 sub-intervals no vertex moves farther than mu |s2 - s1|;
!< 4. the bound is attained: in a pure rotation the farthest vertex moves on a circle of radius mu / theta (chord
!<    2 rho sin(theta ds / 2), to 1e-11);
!< 5. negative control: the bound without its translation term fails the check of 3;
!< 6. 2-D: the same checks, the angle difference across +-pi taken as `pose_interpolate` does.

use fossil, only : rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, pose_2d, pose_interpolate, to_world, &
                   motion_bound_3d, motion_bound_2d
use penf,   only : I4P, R8P

implicit none

integer(I4P), parameter :: NS = 10000_I4P
real(R8P)               :: x(8), y(8), z(8), u(4), v(4), mu, mu0, ratio, worst
type(rigid_pose_t)      :: p0, p1
type(rigid_pose_2d_t)   :: q0, q1
integer(I4P)            :: i
logical                 :: all_ok

all_ok = .true.
x = [1._R8P, 2._R8P, 1._R8P, 2._R8P, 1._R8P, 2._R8P, 1._R8P, 2._R8P]
y = [-0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P, -0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P]
z = [3._R8P, 3._R8P, 3._R8P, 3._R8P, 4._R8P, 4._R8P, 4._R8P, 4._R8P]

! 1. pure translation
p0 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.3_R8P, [1._R8P, 2._R8P, 3._R8P])
p1 = p0 ; p1%t = p0%t + [0.3_R8P, -0.4_R8P, 1.2_R8P]
mu = motion_bound_3d(p0, p1, x, y, z)
call check(abs(mu - 1.3_R8P) <= 1e-14_R8P, 'pure translation: mu = |t1 - t0| = 1.3', mu)

! 2. pure rotation about z through the body origin
p0 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P)
p1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.7_R8P)
mu = motion_bound_3d(p0, p1, x, y, z)
call check(abs(mu - 0.7_R8P * sqrt(4.25_R8P)) <= 1e-14_R8P, 'pure rotation: mu = theta max rho = 0.7 sqrt(4.25)', mu)

! 3. rotation + translation, general axes
p0 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.5_R8P, -1._R8P, 2._R8P])
p1 = pose_from_axis_angle([-1._R8P, 0.5_R8P, 2._R8P], 1.3_R8P, [1.5_R8P, 0.2_R8P, 1._R8P])
mu = motion_bound_3d(p0, p1, x, y, z)
worst = sampled_3d(p0, p1, mu)
call check(worst <= 1._R8P + 1e-12_R8P, 'rotation + translation: displacement <= mu ds on 10^4 sub-intervals (max ratio)', &
           worst)

! 4. attained: pure rotation, general axis, chord of the farthest vertex
p1 = pose_from_axis_angle([2._R8P, -1._R8P, 0.5_R8P], 0.9_R8P)
p1%r = matmul(p1%r, p0%r) ; p1%t = p0%t
mu = motion_bound_3d(p0, p1, x, y, z)
ratio = attained_3d(p0, p1, mu)
call check(abs(ratio - 1._R8P) <= 1e-11_R8P, 'pure rotation: the farthest vertex attains the bound (chord ratio)', ratio)

! 5. negative control: the translation term dropped
p1 = pose_from_axis_angle([-1._R8P, 0.5_R8P, 2._R8P], 1.3_R8P, [1.5_R8P, 0.2_R8P, 1._R8P])
mu = motion_bound_3d(p0, p1, x, y, z)
mu0 = mu - norm2(p1%t - p0%t)
worst = sampled_3d(p0, p1, mu0)
call check(worst > 1._R8P + 1e-3_R8P, 'negative control: without the translation term the check fails (max ratio)', worst)

! 6. 2-D
u = [1._R8P, 2._R8P, 2._R8P, 1._R8P] ; v = [-0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P]
q0 = pose_2d(0._R8P, [1._R8P, 2._R8P]) ; q1 = pose_2d(0._R8P, [1.3_R8P, 1.6_R8P])
mu = motion_bound_2d(q0, q1, u, v)
call check(abs(mu - 0.5_R8P) <= 1e-14_R8P, '2-D pure translation: mu = 0.5', mu)
q0 = pose_2d(3._R8P, [0._R8P, 0._R8P]) ; q1 = pose_2d(-3._R8P, [0._R8P, 0._R8P])
mu = motion_bound_2d(q0, q1, u, v)
call check(abs(mu - (2._R8P * acos(-1._R8P) - 6._R8P) * sqrt(4.25_R8P)) <= 1e-13_R8P, &
           '2-D rotation across pi: theta = 2 pi - 6, mu = theta max |x_b|', mu)
ratio = attained_2d(q0, q1, mu)
call check(abs(ratio - 1._R8P) <= 1e-11_R8P, '2-D pure rotation: the bound is attained (chord ratio)', ratio)
q0 = pose_2d(0.4_R8P, [0.5_R8P, -1._R8P]) ; q1 = pose_2d(-1.1_R8P, [1.5_R8P, 0.2_R8P])
mu = motion_bound_2d(q0, q1, u, v)
worst = sampled_2d(q0, q1, mu)
call check(worst <= 1._R8P + 1e-12_R8P, '2-D rotation + translation: displacement <= mu ds (max ratio)', worst)
worst = sampled_2d(q0, q1, mu - norm2(q1%t - q0%t))
call check(worst > 1._R8P + 1e-3_R8P, '2-D negative control: without the translation term the check fails', worst)

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

   function sampled_3d(a, b, bound) result(w)
   !< Largest displacement / (bound ds) over NS sub-intervals and the 8 vertices.
   type(rigid_pose_t), intent(in) :: a, b  !< End poses.
   real(R8P),          intent(in) :: bound !< Bound tested.
   real(R8P)                      :: w     !< Worst ratio.
   real(R8P)                      :: s1, s2, d
   integer(I4P)                   :: j, k

   w = 0._R8P
   do j = 1, NS
      s1 = real(j - 1, R8P) / NS ; s2 = real(j, R8P) / NS
      do k = 1, 8
         d = norm2(to_world(pose_interpolate(a, b, s2), [x(k), y(k), z(k)]) - &
                   to_world(pose_interpolate(a, b, s1), [x(k), y(k), z(k)]))
         w = max(w, d / (bound * (s2 - s1)))
      enddo
   enddo
   endfunction sampled_3d

   function attained_3d(a, b, bound) result(r)
   !< Largest chord over 1/100 of the step, over the chord of the bound's circle: 2 (mu / theta) sin(theta ds / 2).
   type(rigid_pose_t), intent(in) :: a, b  !< End poses (pure rotation).
   real(R8P),          intent(in) :: bound !< Bound.
   real(R8P)                      :: r     !< Ratio.
   real(R8P)                      :: c, theta, ds
   type(rigid_pose_t)             :: rel
   integer(I4P)                   :: k

   ds = 0.01_R8P
   rel%r = matmul(b%r, transpose(a%r))
   theta = acos(max(-1._R8P, min(1._R8P, 0.5_R8P * (rel%r(1, 1) + rel%r(2, 2) + rel%r(3, 3) - 1._R8P))))
   c = 0._R8P
   do k = 1, 8
      c = max(c, norm2(to_world(pose_interpolate(a, b, 0.5_R8P + ds), [x(k), y(k), z(k)]) - &
                       to_world(pose_interpolate(a, b, 0.5_R8P), [x(k), y(k), z(k)])))
   enddo
   r = c / (2._R8P * bound / theta * sin(theta * ds / 2._R8P))
   endfunction attained_3d

   function sampled_2d(a, b, bound) result(w)
   !< 2-D twin of sampled_3d.
   type(rigid_pose_2d_t), intent(in) :: a, b  !< End poses.
   real(R8P),             intent(in) :: bound !< Bound tested.
   real(R8P)                         :: w     !< Worst ratio.
   real(R8P)                         :: s1, s2, d
   integer(I4P)                      :: j, k

   w = 0._R8P
   do j = 1, NS
      s1 = real(j - 1, R8P) / NS ; s2 = real(j, R8P) / NS
      do k = 1, 4
         d = norm2(to_world(pose_interpolate(a, b, s2), [u(k), v(k)]) - to_world(pose_interpolate(a, b, s1), [u(k), v(k)]))
         w = max(w, d / (bound * (s2 - s1)))
      enddo
   enddo
   endfunction sampled_2d

   function attained_2d(a, b, bound) result(r)
   !< 2-D twin of attained_3d (theta = 2 pi - 6 here).
   type(rigid_pose_2d_t), intent(in) :: a, b  !< End poses (pure rotation).
   real(R8P),             intent(in) :: bound !< Bound.
   real(R8P)                         :: r     !< Ratio.
   real(R8P)                         :: c, theta, ds
   integer(I4P)                      :: k

   ds = 0.01_R8P ; theta = 2._R8P * acos(-1._R8P) - 6._R8P
   c = 0._R8P
   do k = 1, 4
      c = max(c, norm2(to_world(pose_interpolate(a, b, 0.5_R8P + ds), [u(k), v(k)]) - &
                       to_world(pose_interpolate(a, b, 0.5_R8P), [u(k), v(k)])))
   enddo
   r = c / (2._R8P * bound / theta * sin(theta * ds / 2._R8P))
   endfunction attained_2d
endprogram fossil_test_ccd_bounds
