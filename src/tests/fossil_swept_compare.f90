!< FOSSIL test support: bitwise comparison of two 2-D swept-geometry results (issue #100 step 4).
module fossil_swept_compare
!< FOSSIL test support: bitwise comparison of two 2-D swept-geometry results (issue #100 step 4).
!<
!< `swept_cells_2d_mismatches(a, b)` counts the differing elements of every output of `swept_cells_2d`: the scalars nx,
!< ny, exact, error_estimate (bits), evaluations, self_intersecting, and the arrays alpha0, alpha1, swept, face_u,
!< face_v, wall, residual (reals by their bits: `transfer` to a 64-bit integer), cleared, covered. An array of another
!< shape, or allocated in one result only, counts as max(1, its size). `swept_twin_mismatches(c, p0, p1, x, y, ...)` runs
!< `swept_cells_2d` and `swept_cells_2d_device` with the same arguments and compares them. Used by
!< `fossil_test_swept_2d_device_twin` (CI, host) and by the curve bench (device vs serial).
use fossil, only : curve_object, rigid_pose_2d_t, swept_cells_2d_t, swept_cells_2d, swept_cells_2d_device
use penf,   only : I4P, I8P, R8P

implicit none
private
public :: swept_cells_2d_mismatches, swept_twin_mismatches

contains
   function swept_twin_mismatches(c, p0, p1, x, y, nq, poses, exact, refine) result(n)
   !< Mismatches of `swept_cells_2d_device` against `swept_cells_2d` (same arguments), statuses included.
   type(curve_object),    intent(in)           :: c          !< Curve.
   type(rigid_pose_2d_t), intent(in)           :: p0, p1     !< End poses.
   real(R8P),             intent(in)           :: x(:), y(:) !< Cell edges.
   integer(I4P),          intent(in), optional :: nq         !< Gauss nodes.
   type(rigid_pose_2d_t), intent(in), optional :: poses(:)   !< Node poses.
   logical,               intent(in), optional :: exact      !< Exact paths.
   logical,               intent(in), optional :: refine     !< Gauss path refined per event interval (issue #104).
   integer(I4P)                                :: n          !< Mismatches.
   type(swept_cells_2d_t)                      :: h, d       !< Host, device.
   integer(I4P)                                :: sh, sd     !< Statuses.

   call swept_cells_2d(c, p0, p1, x, y, h, nq=nq, poses=poses, exact=exact, status=sh, refine=refine)
   call swept_cells_2d_device(c, p0, p1, x, y, d, nq=nq, poses=poses, exact=exact, status=sd, refine=refine)
   n = swept_cells_2d_mismatches(h, d) + abs(sh - sd)
   endfunction swept_twin_mismatches

   function swept_cells_2d_mismatches(a, b) result(n)
   !< Differing elements of two 2-D swept results (see the module notes).
   type(swept_cells_2d_t), intent(in) :: a, b !< Results.
   integer(I4P)                       :: n    !< Mismatches.

   n = 0_I4P
   if (a%nx /= b%nx) n = n + 1_I4P
   if (a%ny /= b%ny) n = n + 1_I4P
   if (a%exact .neqv. b%exact) n = n + 1_I4P
   if (transfer(a%error_estimate, 0_I8P) /= transfer(b%error_estimate, 0_I8P)) n = n + 1_I4P
   if (a%evaluations /= b%evaluations) n = n + 1_I4P
   if (a%self_intersecting .neqv. b%self_intersecting) n = n + 1_I4P
   n = n + r2(a%alpha0, b%alpha0) + r2(a%alpha1, b%alpha1) + r2(a%swept, b%swept) + r2(a%face_u, b%face_u) + &
       r2(a%face_v, b%face_v) + r2(a%residual, b%residual) + r3(a%wall, b%wall) + l2(a%cleared, b%cleared) + &
       l2(a%covered, b%covered)
   endfunction swept_cells_2d_mismatches

   function r2(a, b) result(n)
   !< Bitwise mismatches of two real rank-2 arrays.
   real(R8P), allocatable, intent(in) :: a(:,:), b(:,:) !< Arrays.
   integer(I4P)                       :: n              !< Mismatches.

   n = 0_I4P
   if (allocated(a) .neqv. allocated(b)) then
      if (allocated(a)) n = max(1_I4P, size(a, kind=I4P))
      if (allocated(b)) n = max(1_I4P, size(b, kind=I4P))
      return
   endif
   if (.not. allocated(a)) return
   if (any(shape(a) /= shape(b))) then
      n = max(1_I4P, size(a, kind=I4P)) ; return
   endif
   n = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   endfunction r2

   function r3(a, b) result(n)
   !< Bitwise mismatches of two real rank-3 arrays.
   real(R8P), allocatable, intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.
   integer(I4P)                       :: n                  !< Mismatches.

   n = 0_I4P
   if (allocated(a) .neqv. allocated(b)) then
      if (allocated(a)) n = max(1_I4P, size(a, kind=I4P))
      if (allocated(b)) n = max(1_I4P, size(b, kind=I4P))
      return
   endif
   if (.not. allocated(a)) return
   if (any(shape(a) /= shape(b))) then
      n = max(1_I4P, size(a, kind=I4P)) ; return
   endif
   n = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   endfunction r3

   function l2(a, b) result(n)
   !< Mismatches of two logical rank-2 arrays.
   logical, allocatable, intent(in) :: a(:,:), b(:,:) !< Arrays.
   integer(I4P)                     :: n              !< Mismatches.

   n = 0_I4P
   if (allocated(a) .neqv. allocated(b)) then
      if (allocated(a)) n = max(1_I4P, size(a, kind=I4P))
      if (allocated(b)) n = max(1_I4P, size(b, kind=I4P))
      return
   endif
   if (.not. allocated(a)) return
   if (any(shape(a) /= shape(b))) then
      n = max(1_I4P, size(a, kind=I4P)) ; return
   endif
   n = count(a .neqv. b, kind=I4P)
   endfunction l2
endmodule fossil_swept_compare
