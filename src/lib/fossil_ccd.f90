!< FOSSIL, continuous collision detection: time of first contact by conservative advancement (issue #73).
module fossil_ccd
!< FOSSIL, continuous collision detection: time of first contact by conservative advancement (issue #73).
!<
!< Design: docs/guide/advanced/collision.md. The motion over a step is the one of `pose_interpolate`: x(s) =
!< Q(s) R0 x_b + t0 + s (t1 - t0), Q(s) turning about the fixed unit axis k of R1 R0^T by s theta (the axis through
!< the body-frame origin, which moves with t(s)). The speed of a body point is at most |t1 - t0| + |theta| rho(x_b),
!< rho = |k x (R0 x_b)| (constant in s); rho is convex, so its maximum over a body is at a vertex: the motion bound
!< mu = |t1 - t0| + |theta| max_v rho_v holds for every point of the body and every s (`motion_bound_3d`,
!< `motion_bound_2d`; in 2-D rho = |x_b|).
use penf,              only : I4P, R8P
use fossil_rigid_pose, only : rigid_pose_t, rigid_pose_2d_t, pose_log_3d, pose_angle_2d

implicit none
private
public :: motion_bound_3d, motion_bound_2d

contains
   pure function motion_bound_3d(p0, p1, x, y, z) result(mu)
   !< Bound on the speed |dx/ds| of the vertices (x, y, z) (body frame) moved by `pose_interpolate(p0, p1, s)`, s in
   !< [0, 1]: no point of their hull moves farther than mu |s2 - s1| between s1 and s2.
   type(rigid_pose_t), intent(in) :: p0, p1         !< End poses.
   real(R8P),          intent(in) :: x(:), y(:), z(:) !< Vertices (body frame).
   real(R8P)                      :: mu             !< Bound.
   real(R8P)                      :: k(3), theta, w(3), c(3), rho !< Axis, angle, rotated vertex, cross, distance.
   integer(I4P)                   :: i              !< Counter.

   call pose_log_3d(p0, p1, k, theta)
   rho = 0._R8P
   if (theta > 0._R8P) then
      do i = 1, size(x, kind=I4P)
         w = p0%r(:, 1) * x(i) + p0%r(:, 2) * y(i) + p0%r(:, 3) * z(i)
         c = [k(2) * w(3) - k(3) * w(2), k(3) * w(1) - k(1) * w(3), k(1) * w(2) - k(2) * w(1)]
         rho = max(rho, norm2(c))
      enddo
   endif
   mu = norm2(p1%t - p0%t) + abs(theta) * rho
   endfunction motion_bound_3d

   pure function motion_bound_2d(p0, p1, u, v) result(mu)
   !< 2-D twin of `motion_bound_3d`: the angle difference in (-pi, pi] as `pose_interpolate`, rho = |(u, v)|.
   type(rigid_pose_2d_t), intent(in) :: p0, p1     !< End poses.
   real(R8P),             intent(in) :: u(:), v(:) !< Vertices (body frame).
   real(R8P)                         :: mu         !< Bound.
   real(R8P)                         :: theta, rho !< Angle, distance.
   integer(I4P)                      :: i          !< Counter.

   theta = pose_angle_2d(p1) - pose_angle_2d(p0)
   theta = atan2(sin(theta), cos(theta))
   rho = 0._R8P
   if (theta /= 0._R8P) then
      do i = 1, size(u, kind=I4P)
         rho = max(rho, hypot(u(i), v(i)))
      enddo
   endif
   mu = norm2(p1%t - p0%t) + abs(theta) * rho
   endfunction motion_bound_2d
endmodule fossil_ccd
