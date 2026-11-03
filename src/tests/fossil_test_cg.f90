!< FOSSIL test: conjugate gradient on CSR matrices (issue #31, step 5).
program fossil_test_cg
!< FOSSIL test: conjugate gradient on CSR matrices (issue #31, step 5).
!<
!<  - 1D Poisson (tridiagonal 2, -1; n = 50) with the right-hand side of a known solution: the
!<    solution recovered to 1e-10, converged within n iterations (exact arithmetic would take n);
!<  - 2D Poisson (5-point, 12 x 12): residual below the tolerance, status OK;
!<  - b = 0 gives x = 0 at once;
!<  - an indefinite matrix (diag(1, -1)) is refused with CG_STATUS_NOT_SPD;
!<  - an iteration cap of 2 on the 2D system gives CG_STATUS_MAX_ITER;
!<  - a size mismatch gives CG_STATUS_BAD_INPUT.
!< All systems are tiny.

use fossil_csr_matrix,    only : csr_matrix_t
use fossil_linear_solver, only : cg_solve, CG_STATUS_OK, CG_STATUS_BAD_INPUT, CG_STATUS_NOT_SPD, CG_STATUS_MAX_ITER
use penf,                 only : I4P, R8P

implicit none

type(csr_matrix_t)     :: a1, a2, ai
real(R8P), allocatable :: x(:), b(:), xs(:)
real(R8P)              :: res
integer(I4P)           :: n, i, j, k, it, st
logical                :: ok

ok = .true.

! 1D Poisson, known solution xs(i) = sin(i)
n = 50
call a1%initialize(n, n)
do i = 1, n
   call a1%append(i, i, 2._R8P)
   if (i > 1) call a1%append(i, i - 1, -1._R8P)
   if (i < n) call a1%append(i, i + 1, -1._R8P)
enddo
call a1%finalize
allocate(xs(n), b(n), x(n))
do i = 1, n
   xs(i) = sin(real(i, R8P))
enddo
call a1%multiply_vector(xs, b)
x = 0._R8P
call cg_solve(a1, b, x, iterations=it, residual=res, status=st)
print '(A,I0,A,I0,A,ES10.3,A,ES10.3)', '1D Poisson: status ', st, '  iterations ', it, '  residual ', res, &
      '  max error ', maxval(abs(x - xs))
ok = ok .and. st == CG_STATUS_OK .and. it <= n .and. maxval(abs(x - xs)) < 1.e-10_R8P

! 2D Poisson 12 x 12
n = 12
call a2%initialize(n * n, n * n)
do j = 1, n
   do i = 1, n
      k = i + n * (j - 1)
      call a2%append(k, k, 4._R8P)
      if (i > 1) call a2%append(k, k - 1, -1._R8P)
      if (i < n) call a2%append(k, k + 1, -1._R8P)
      if (j > 1) call a2%append(k, k - n, -1._R8P)
      if (j < n) call a2%append(k, k + n, -1._R8P)
   enddo
enddo
call a2%finalize
deallocate(b, x) ; allocate(b(n * n), x(n * n))
b = 1._R8P ; x = 0._R8P
call cg_solve(a2, b, x, tol=1.e-12_R8P, iterations=it, residual=res, status=st)
print '(A,I0,A,I0,A,ES10.3)', '2D Poisson: status ', st, '  iterations ', it, '  residual ', res
ok = ok .and. st == CG_STATUS_OK .and. res <= 1.e-12_R8P

! zero right-hand side
b = 0._R8P ; x = 3._R8P
call cg_solve(a2, b, x, iterations=it, status=st)
ok = ok .and. st == CG_STATUS_OK .and. it == 0 .and. all(x == 0._R8P)
print '(A,L1)', 'b = 0 -> x = 0: ', st == CG_STATUS_OK .and. all(x == 0._R8P)

! iteration cap
b = 1._R8P ; x = 0._R8P
call cg_solve(a2, b, x, max_iter=2_I4P, iterations=it, status=st)
ok = ok .and. st == CG_STATUS_MAX_ITER .and. it == 2
print '(A,I0)', 'cap of 2 iterations: status ', st

! indefinite matrix: diag(1, -1) has a negative diagonal -> refused before iterating;
! [[1, 2], [2, 1]] has a positive diagonal but is indefinite -> NOT_SPD during the iteration
call ai%initialize(2, 2)
call ai%append(1, 1, 1._R8P) ; call ai%append(1, 2, 2._R8P)
call ai%append(2, 1, 2._R8P) ; call ai%append(2, 2, 1._R8P)
call ai%finalize
deallocate(b, x) ; allocate(b(2), x(2))
b = [1._R8P, -1._R8P] ; x = 0._R8P
call cg_solve(ai, b, x, status=st)
print '(A,I0)', 'indefinite matrix: status ', st
ok = ok .and. st == CG_STATUS_NOT_SPD

! size mismatch
deallocate(b) ; allocate(b(3)) ; b = 1._R8P
call cg_solve(ai, b, x, status=st)
ok = ok .and. st == CG_STATUS_BAD_INPUT
print '(A,I0)', 'size mismatch: status ', st

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_cg
