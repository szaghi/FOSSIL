!< FOSSIL, rigid poses for moving bodies (issue #41).
module fossil_rigid_pose
!< FOSSIL, rigid poses for moving bodies (issue #41).
!<
!< A pose maps the body frame to the world frame, x_world = R x_body + t, with R a rotation. The
!< geometry stays in the body frame; a query at a pose maps its points to the body frame,
!< x_body = R^T (x_world - t), and maps closest points and normals back (decision 2 of #41).
!< `rigid_pose_t` is the 3-D pose (R(3, 3), t(3)), `rigid_pose_2d_t` the in-plane pose of a curve
!< (R(2, 2), t(2), a rotation about the curve's normal axis). Plain arrays, so the device kernels can
!< use them: the point maps are `!$acc routine seq`.
use penf, only : I4P, R8P

implicit none
private
public :: rigid_pose_t, rigid_pose_2d_t
public :: pose_identity, pose_from_axis_angle, pose_from_quaternion, pose_2d, pose_compose, pose_inverse
public :: to_body, to_world, rotate, pose_interpolate, pose_angle_2d, pose_log_3d

type :: rigid_pose_t
   !< 3-D rigid pose: x_world = r x_body + t.
   real(R8P) :: r(3, 3) = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) !< Rotation.
   real(R8P) :: t(3)    = 0._R8P !< Translation.
endtype rigid_pose_t

type :: rigid_pose_2d_t
   !< In-plane rigid pose of a curve: (u, v)_world = r (u, v)_body + t.
   real(R8P) :: r(2, 2) = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) !< Rotation.
   real(R8P) :: t(2)    = 0._R8P !< Translation.
endtype rigid_pose_2d_t

interface pose_compose
   module procedure compose_3d, compose_2d
endinterface pose_compose

interface pose_inverse
   module procedure inverse_3d, inverse_2d
endinterface pose_inverse

interface to_body
   module procedure to_body_3d, to_body_2d
endinterface to_body

interface to_world
   module procedure to_world_3d, to_world_2d
endinterface to_world

interface rotate
   module procedure rotate_3d, rotate_2d
endinterface rotate

interface pose_interpolate
   module procedure interpolate_3d, interpolate_2d
endinterface pose_interpolate

contains
   pure function pose_identity() result(p)
   !< The identity pose.
   type(rigid_pose_t) :: p !< Pose.

   p = rigid_pose_t()
   endfunction pose_identity

   pure function pose_from_axis_angle(axis, angle, translation) result(p)
   !< Rotation by `angle` (radians, right-handed) about `axis` (any length > 0), then translation.
   real(R8P), intent(in)           :: axis(3)        !< Rotation axis.
   real(R8P), intent(in)           :: angle          !< Angle [rad].
   real(R8P), intent(in), optional :: translation(3) !< Translation (default 0).
   type(rigid_pose_t)              :: p              !< Pose.
   real(R8P)                       :: k(3), c, s     !< Unit axis, cosine, sine.

   k = axis / norm2(axis) ; c = cos(angle) ; s = sin(angle)
   p%r(1, :) = [c + k(1) * k(1) * (1 - c), k(1) * k(2) * (1 - c) - k(3) * s, k(1) * k(3) * (1 - c) + k(2) * s]
   p%r(2, :) = [k(2) * k(1) * (1 - c) + k(3) * s, c + k(2) * k(2) * (1 - c), k(2) * k(3) * (1 - c) - k(1) * s]
   p%r(3, :) = [k(3) * k(1) * (1 - c) - k(2) * s, k(3) * k(2) * (1 - c) + k(1) * s, c + k(3) * k(3) * (1 - c)]
   p%t = 0._R8P ; if (present(translation)) p%t = translation
   endfunction pose_from_axis_angle

   pure function pose_from_quaternion(q, translation) result(p)
   !< Rotation of the quaternion q = (w, x, y, z) (normalized here), then translation.
   real(R8P), intent(in)           :: q(4)           !< Quaternion (w, x, y, z).
   real(R8P), intent(in), optional :: translation(3) !< Translation (default 0).
   type(rigid_pose_t)              :: p              !< Pose.
   real(R8P)                       :: w, x, y, z, n  !< Work.

   n = norm2(q) ; w = q(1) / n ; x = q(2) / n ; y = q(3) / n ; z = q(4) / n
   p%r(1, :) = [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)]
   p%r(2, :) = [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)]
   p%r(3, :) = [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]
   p%t = 0._R8P ; if (present(translation)) p%t = translation
   endfunction pose_from_quaternion

   pure function pose_2d(angle, translation) result(p)
   !< In-plane rotation by `angle` (counter-clockwise in (u, v)), then translation.
   real(R8P), intent(in)           :: angle          !< Angle [rad].
   real(R8P), intent(in), optional :: translation(2) !< Translation (default 0).
   type(rigid_pose_2d_t)           :: p              !< Pose.

   p%r = reshape([cos(angle), sin(angle), -sin(angle), cos(angle)], [2, 2])
   p%t = 0._R8P ; if (present(translation)) p%t = translation
   endfunction pose_2d

   pure function compose_3d(a, b) result(p)
   !< a o b: apply b, then a.
   type(rigid_pose_t), intent(in) :: a, b !< Poses.
   type(rigid_pose_t)             :: p    !< Composition.

   p%r = matmul(a%r, b%r) ; p%t = matmul(a%r, b%t) + a%t
   endfunction compose_3d

   pure function compose_2d(a, b) result(p)
   !< a o b: apply b, then a.
   type(rigid_pose_2d_t), intent(in) :: a, b !< Poses.
   type(rigid_pose_2d_t)             :: p    !< Composition.

   p%r = matmul(a%r, b%r) ; p%t = matmul(a%r, b%t) + a%t
   endfunction compose_2d

   pure function inverse_3d(a) result(p)
   !< Inverse pose: R^T, -R^T t.
   type(rigid_pose_t), intent(in) :: a !< Pose.
   type(rigid_pose_t)             :: p !< Inverse.

   p%r = transpose(a%r) ; p%t = -matmul(p%r, a%t)
   endfunction inverse_3d

   pure function inverse_2d(a) result(p)
   !< Inverse pose: R^T, -R^T t.
   type(rigid_pose_2d_t), intent(in) :: a !< Pose.
   type(rigid_pose_2d_t)             :: p !< Inverse.

   p%r = transpose(a%r) ; p%t = -matmul(p%r, a%t)
   endfunction inverse_2d

   pure function to_body_3d(p, x) result(y)
   !< World point to body frame: R^T (x - t).
   !$acc routine seq
   type(rigid_pose_t), intent(in) :: p    !< Pose.
   real(R8P),          intent(in) :: x(3) !< World point.
   real(R8P)                      :: y(3) !< Body point.
   real(R8P)                      :: d(3) !< x - t.

   d = x - p%t
   y(1) = p%r(1, 1) * d(1) + p%r(2, 1) * d(2) + p%r(3, 1) * d(3)
   y(2) = p%r(1, 2) * d(1) + p%r(2, 2) * d(2) + p%r(3, 2) * d(3)
   y(3) = p%r(1, 3) * d(1) + p%r(2, 3) * d(2) + p%r(3, 3) * d(3)
   endfunction to_body_3d

   pure function to_world_3d(p, x) result(y)
   !< Body point to world frame: R x + t.
   !$acc routine seq
   type(rigid_pose_t), intent(in) :: p    !< Pose.
   real(R8P),          intent(in) :: x(3) !< Body point.
   real(R8P)                      :: y(3) !< World point.

   y = rotate_3d(p, x) + p%t
   endfunction to_world_3d

   pure function rotate_3d(p, v) result(w)
   !< Rotate a vector: R v.
   !$acc routine seq
   type(rigid_pose_t), intent(in) :: p    !< Pose.
   real(R8P),          intent(in) :: v(3) !< Vector.
   real(R8P)                      :: w(3) !< Rotated vector.

   w(1) = p%r(1, 1) * v(1) + p%r(1, 2) * v(2) + p%r(1, 3) * v(3)
   w(2) = p%r(2, 1) * v(1) + p%r(2, 2) * v(2) + p%r(2, 3) * v(3)
   w(3) = p%r(3, 1) * v(1) + p%r(3, 2) * v(2) + p%r(3, 3) * v(3)
   endfunction rotate_3d

   pure function to_body_2d(p, x) result(y)
   !< World point to body frame: R^T (x - t).
   !$acc routine seq
   type(rigid_pose_2d_t), intent(in) :: p    !< Pose.
   real(R8P),             intent(in) :: x(2) !< World point.
   real(R8P)                         :: y(2) !< Body point.
   real(R8P)                         :: d(2) !< x - t.

   d = x - p%t
   y(1) = p%r(1, 1) * d(1) + p%r(2, 1) * d(2)
   y(2) = p%r(1, 2) * d(1) + p%r(2, 2) * d(2)
   endfunction to_body_2d

   pure function to_world_2d(p, x) result(y)
   !< Body point to world frame: R x + t.
   !$acc routine seq
   type(rigid_pose_2d_t), intent(in) :: p    !< Pose.
   real(R8P),             intent(in) :: x(2) !< Body point.
   real(R8P)                         :: y(2) !< World point.

   y = rotate_2d(p, x) + p%t
   endfunction to_world_2d

   pure function rotate_2d(p, v) result(w)
   !< Rotate a vector: R v.
   !$acc routine seq
   type(rigid_pose_2d_t), intent(in) :: p    !< Pose.
   real(R8P),             intent(in) :: v(2) !< Vector.
   real(R8P)                         :: w(2) !< Rotated vector.

   w(1) = p%r(1, 1) * v(1) + p%r(1, 2) * v(2)
   w(2) = p%r(2, 1) * v(1) + p%r(2, 2) * v(2)
   endfunction rotate_2d

   pure function pose_angle_2d(p) result(theta)
   !< Rotation angle of an in-plane pose, in (-pi, pi].
   type(rigid_pose_2d_t), intent(in) :: p     !< Pose.
   real(R8P)                         :: theta !< Angle [rad].

   theta = atan2(p%r(2, 1), p%r(1, 1))
   endfunction pose_angle_2d

   pure subroutine pose_log_3d(p0, p1, axis, angle)
   !< Relative rotation R1 R0^T of two poses as a unit axis and an angle in [0, pi] (axis = e_z for no
   !< rotation).
   type(rigid_pose_t), intent(in)  :: p0, p1  !< Poses.
   real(R8P),          intent(out) :: axis(3) !< Unit axis.
   real(R8P),          intent(out) :: angle   !< Angle [rad].
   real(R8P)                       :: r(3, 3), w(3), c, s !< Relative rotation, work.
   integer                         :: k       !< Column.

   r = matmul(p1%r, transpose(p0%r))
   w = [r(3, 2) - r(2, 3), r(1, 3) - r(3, 1), r(2, 1) - r(1, 2)] ! 2 sin(angle) axis
   s = 0.5_R8P * norm2(w) ; c = 0.5_R8P * (r(1, 1) + r(2, 2) + r(3, 3) - 1._R8P)
   angle = atan2(s, c)
   axis = [0._R8P, 0._R8P, 1._R8P]
   if (angle <= 0._R8P) return
   if (s > 1.e-6_R8P) then
      axis = w / norm2(w)
   else ! angle near pi: axis from R + I = 2 k k^T (+ O(s)), its largest column
      r = r + reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3])
      k = maxloc([r(1, 1), r(2, 2), r(3, 3)], dim=1)
      axis = r(:, k) / norm2(r(:, k))
      if (dot_product(axis, w) < 0._R8P) axis = -axis
   endif
   endsubroutine pose_log_3d

   pure function interpolate_2d(p0, p1, s) result(p)
   !< Screw interpolation of in-plane poses (issue #41): the angle and the translation linear in s in
   !< [0, 1], the angle difference taken in (-pi, pi]; the end poses are returned exactly.
   type(rigid_pose_2d_t), intent(in) :: p0, p1 !< End poses.
   real(R8P),             intent(in) :: s      !< Parameter in [0, 1].
   type(rigid_pose_2d_t)             :: p      !< Pose at s.
   real(R8P)                         :: a0, da !< Angles.

   if (s == 0._R8P) then
      p = p0 ; return
   elseif (s == 1._R8P) then
      p = p1 ; return
   endif
   a0 = pose_angle_2d(p0) ; da = pose_angle_2d(p1) - a0
   da = atan2(sin(da), cos(da))
   p = pose_2d(a0 + s * da, p0%t + s * (p1%t - p0%t))
   endfunction interpolate_2d

   pure function interpolate_3d(p0, p1, s) result(p)
   !< Screw interpolation of 3-D poses (issue #41): the rotation R(s) = exp(s log(R1 R0^T)) R0 (fixed
   !< axis, angle linear in s) and the translation linear in s; the end poses are returned exactly.
   type(rigid_pose_t), intent(in) :: p0, p1  !< End poses.
   real(R8P),          intent(in) :: s       !< Parameter in [0, 1].
   type(rigid_pose_t)             :: p       !< Pose at s.
   type(rigid_pose_t)             :: q       !< Partial relative rotation.
   real(R8P)                      :: axis(3), angle !< Relative rotation.

   if (s == 0._R8P) then
      p = p0 ; return
   elseif (s == 1._R8P) then
      p = p1 ; return
   endif
   call pose_log_3d(p0, p1, axis, angle)
   q = pose_from_axis_angle(axis, s * angle)
   p%r = matmul(q%r, p0%r)
   p%t = p0%t + s * (p1%t - p0%t)
   endfunction interpolate_3d
endmodule fossil_rigid_pose
