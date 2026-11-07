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
!<
!< Conservative advancement (`ccd_advance`, Mirtich 1996): the distance d(s) between two bodies with bounds mu_A, mu_B
!< is Lipschitz with constant mu = mu_A + mu_B. From s_k with d_k > tol the step s_{k+1} = min(1, s_k + (d_k - tol / 2)
!< / mu) keeps d >= tol / 2 on [s_k, s_{k+1}]; the loop stops at the first iterate with d <= tol (CCD_CONTACT), at
!< s = 1 with d(1) > tol (CCD_SEPARATE), with d(0) <= tol (CCD_INITIAL_CONTACT, toi = 0) or after `max_iterations`
!< steps (CCD_NOT_CONVERGED, toi = the last iterate: d >= tol / 2 on [0, toi], a certified lower bound on the contact
!< time). Each step with d_k > tol advances at least tol / (2 mu): at most ceiling(2 mu / tol) steps. The distance
!< comes from an extension of `ccd_evaluator_t` (no body type here; an internal procedure as actual argument made
!< gfortran put trampolines on an executable stack, issue #60).
!<
!< `directional=.true.`: the evaluator's `safe_step` (C2A's per-pair directional bound, e.g. `body_safe_step`) gives a
!< second certified step from s_k; the loop takes the larger of the two, so it never takes more steps than the plain
!< bound. The default `safe_step` gives none.
use penf,              only : I4P, R8P
use fossil_rigid_pose, only : rigid_pose_t, rigid_pose_2d_t, pose_log_3d, pose_angle_2d

implicit none
private
public :: motion_bound_3d, motion_bound_2d
public :: ccd_evaluator_t, ccd_advance
public :: CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT, CCD_NOT_CONVERGED, CCD_INVALID_INPUT, CCD_MAX_ITERATIONS

integer(I4P), parameter :: CCD_CONTACT         = 1_I4P    !< The bodies come within the tolerance.
integer(I4P), parameter :: CCD_SEPARATE        = 2_I4P    !< No contact over [0, 1].
integer(I4P), parameter :: CCD_INITIAL_CONTACT = 3_I4P    !< Within the tolerance at s = 0.
integer(I4P), parameter :: CCD_NOT_CONVERGED   = 4_I4P    !< Iteration cap hit: toi is a certified lower bound.
integer(I4P), parameter :: CCD_INVALID_INPUT   = 5_I4P    !< tolerance <= 0, mu < 0 or max_iterations < 1.
integer(I4P), parameter :: CCD_MAX_ITERATIONS  = 1000_I4P !< Default iteration cap.

type, abstract :: ccd_evaluator_t
   !< Distance between two moving bodies at the step parameter s (extended by the caller).
   contains
      procedure(ccd_eval_interface), deferred, pass(self) :: distance  !< d(s).
      procedure,                           pass(self) :: safe_step => ccd_safe_step_none !< Directional step (none).
endtype ccd_evaluator_t

abstract interface
   subroutine ccd_eval_interface(self, s, d)
   !< Distance at s in [0, 1].
   import :: ccd_evaluator_t, R8P
   class(ccd_evaluator_t), intent(inout) :: self !< Evaluator (may keep the witnesses of the last call).
   real(R8P),              intent(in)    :: s    !< Step parameter.
   real(R8P),              intent(out)   :: d    !< Distance.
   endsubroutine ccd_eval_interface
endinterface

contains
   subroutine ccd_advance(evaluator, mu, tolerance, toi, status, distance, iterations, max_iterations, directional)
   !< Time of first contact by conservative advancement (see the module notes). The evaluator's last call is at toi
   !< (with mu = 0, at s = 0: d is constant).
   class(ccd_evaluator_t), intent(inout)         :: evaluator      !< Distance d(s).
   real(R8P),              intent(in)            :: mu             !< Lipschitz bound of d (mu_A + mu_B).
   real(R8P),              intent(in)            :: tolerance      !< Contact tolerance (> 0).
   real(R8P),              intent(out)           :: toi            !< Time of first contact (or 1, or the bound).
   integer(I4P),           intent(out)           :: status         !< CCD_* status.
   real(R8P),              intent(out), optional :: distance       !< d(toi).
   integer(I4P),           intent(out), optional :: iterations     !< Advancement steps (evaluations after s = 0).
   integer(I4P),           intent(in),  optional :: max_iterations !< Cap (default CCD_MAX_ITERATIONS).
   logical,                intent(in),  optional :: directional    !< Also the evaluator's directional step.
   real(R8P)                                     :: s, d, s_next, d_next, dt, dt_dir !< Iterates, steps.
   logical                                       :: dir            !< Directional.
   integer(I4P)                                  :: it, cap        !< Counters.

   toi = 0._R8P ; d = huge(1._R8P) ; it = 0_I4P
   cap = CCD_MAX_ITERATIONS ; if (present(max_iterations)) cap = max_iterations
   dir = .false. ; if (present(directional)) dir = directional
   s = 0._R8P
   if (.not. (tolerance > 0._R8P) .or. .not. (mu >= 0._R8P) .or. cap < 1_I4P) then
      status = CCD_INVALID_INPUT
   else
      call evaluator%distance(s, d)
      if (d <= tolerance) then
         status = CCD_INITIAL_CONTACT
      elseif (mu == 0._R8P) then ! no relative motion: d is constant
         s = 1._R8P ; status = CCD_SEPARATE
      else
         status = CCD_NOT_CONVERGED
         do while (it < cap)
            dt = (d - 0.5_R8P * tolerance) / mu
            if (dir) then
               call evaluator%safe_step(s, tolerance, 1._R8P - s, dt_dir)
               dt = max(dt, dt_dir)
            endif
            s_next = min(1._R8P, s + dt)
            if (.not. (s_next > s)) exit ! the step is below the resolution of s: not converged
            it = it + 1_I4P
            call evaluator%distance(s_next, d_next)
            s = s_next ; d = d_next
            if (d <= tolerance) then
               status = CCD_CONTACT
               exit
            elseif (s == 1._R8P) then
               status = CCD_SEPARATE
               exit
            endif
         enddo
      endif
   endif
   toi = s
   if (present(distance)) distance = d
   if (present(iterations)) iterations = it
   endsubroutine ccd_advance

   subroutine ccd_safe_step_none(self, s, tolerance, cap, step)
   !< No directional step (the default): step = -1.
   class(ccd_evaluator_t), intent(inout) :: self      !< Evaluator.
   real(R8P),              intent(in)    :: s         !< Step parameter.
   real(R8P),              intent(in)    :: tolerance !< Contact tolerance.
   real(R8P),              intent(in)    :: cap       !< Largest step wanted (1 - s).
   real(R8P),              intent(out)   :: step      !< Safe step.

   step = -1._R8P
   endsubroutine ccd_safe_step_none

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
