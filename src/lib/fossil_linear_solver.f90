!< FOSSIL, iterative linear solver for sparse symmetric positive-definite systems (issue #31).
module fossil_linear_solver
!< FOSSIL, iterative linear solver for sparse symmetric positive-definite systems (issue #31).
!<
!< `cg_solve` is the conjugate gradient method with a Jacobi (diagonal) preconditioner on a
!< finalized `csr_matrix_t`. It is meant for the small SPD systems of the geometry operators:
!< the bi-Laplacian fairing of hole patches, and any future use of the cotangent Laplacian. The
!< matrix must be symmetric positive definite; the method detects a non-positive curvature
!< p . A p (an indefinite or singular matrix along the search direction) and stops with
!< CG_STATUS_NOT_SPD instead of returning garbage.
use fossil_csr_matrix, only : csr_matrix_t
use penf,              only : I4P, R8P

implicit none
private
public :: cg_solve
public :: CG_STATUS_OK, CG_STATUS_BAD_INPUT, CG_STATUS_NOT_SPD, CG_STATUS_MAX_ITER

integer(I4P), parameter :: CG_STATUS_OK        = 0_I4P !< Converged: ||r|| <= tol ||b||.
integer(I4P), parameter :: CG_STATUS_BAD_INPUT = 1_I4P !< Matrix not square / not finalized, size mismatch, zero diagonal.
integer(I4P), parameter :: CG_STATUS_NOT_SPD   = 2_I4P !< p . A p <= 0 met: the matrix is not SPD.
integer(I4P), parameter :: CG_STATUS_MAX_ITER  = 3_I4P !< Iteration cap reached before convergence.

contains
   subroutine cg_solve(a, b, x, tol, max_iter, iterations, residual, status)
   !< Solve A x = b by Jacobi-preconditioned conjugate gradient.
   !<
   !< `x` is the initial guess on entry and the solution on return. Converged when the residual
   !< norm is at most `tol` times the norm of b (default 1e-12); a zero b gives x = 0 at once.
   !< `residual` returns the final relative residual ||b - A x|| / ||b|| (0 for b = 0).
   type(csr_matrix_t), intent(in)            :: a          !< SPD matrix (finalized, square).
   real(R8P),          intent(in)            :: b(:)       !< Right-hand side.
   real(R8P),          intent(inout)         :: x(:)       !< Initial guess in, solution out.
   real(R8P),          intent(in),  optional :: tol        !< Relative residual tolerance (default 1e-12).
   integer(I4P),       intent(in),  optional :: max_iter   !< Iteration cap (default 10 n).
   integer(I4P),       intent(out), optional :: iterations !< Iterations done.
   real(R8P),          intent(out), optional :: residual   !< Final relative residual.
   integer(I4P),       intent(out), optional :: status     !< CG_STATUS_*.
   real(R8P),          allocatable           :: r(:), z(:), p(:), ap(:), dinv(:)
   real(R8P)                                 :: tol_, bnorm, rz, rz_new, pap, alpha, rnorm
   integer(I4P)                              :: n, it, max_it, i

   if (present(iterations)) iterations = 0_I4P
   if (present(residual)) residual = 0._R8P
   if (present(status)) status = CG_STATUS_OK
   n = a%get_nrows()
   if (.not. a%is_finalized .or. a%get_ncols() /= n .or. size(b, kind=I4P) /= n .or. size(x, kind=I4P) /= n) then
      if (present(status)) status = CG_STATUS_BAD_INPUT
      return
   endif
   if (n == 0_I4P) return
   tol_   = 1.e-12_R8P ; if (present(tol)) tol_ = tol
   max_it = 10_I4P * n ; if (present(max_iter)) max_it = max_iter

   allocate(r(n), z(n), p(n), ap(n), dinv(n))
   do i = 1, n
      dinv(i) = a%get_value(i, i)
      if (.not. (dinv(i) > 0._R8P)) then ! an SPD matrix has a positive diagonal
         if (present(status)) status = CG_STATUS_BAD_INPUT
         return
      endif
      dinv(i) = 1._R8P / dinv(i)
   enddo
   bnorm = sqrt(sum(b * b))
   if (bnorm == 0._R8P) then
      x = 0._R8P
      return
   endif

   call a%multiply_vector(x, ap)
   r  = b - ap
   z  = dinv * r
   p  = z
   rz = sum(r * z)
   rnorm = sqrt(sum(r * r))
   it = 0_I4P
   do while (rnorm > tol_ * bnorm)
      if (it >= max_it) then
         if (present(status)) status = CG_STATUS_MAX_ITER
         exit
      endif
      call a%multiply_vector(p, ap)
      pap = sum(p * ap)
      if (.not. (pap > 0._R8P)) then
         if (present(status)) status = CG_STATUS_NOT_SPD
         exit
      endif
      alpha = rz / pap
      x = x + alpha * p
      r = r - alpha * ap
      z = dinv * r
      rz_new = sum(r * z)
      p = z + (rz_new / rz) * p
      rz = rz_new
      rnorm = sqrt(sum(r * r))
      it = it + 1_I4P
   enddo
   if (present(iterations)) iterations = it
   if (present(residual)) then
      call a%multiply_vector(x, ap)
      residual = sqrt(sum((b - ap)**2)) / bnorm
   endif
   endsubroutine cg_solve
endmodule fossil_linear_solver
