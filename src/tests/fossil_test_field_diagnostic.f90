!< FOSSIL test: eikonal diagnostic of a distance field (issue #42 step 1).
program fossil_test_field_diagnostic
!< FOSSIL test: eikonal diagnostic of a distance field (issue #42 step 1).
!<
!< - A plane field phi = n . x - d on a stretched 3-D grid: q = 0 to round-off (the three-point
!<   non-uniform difference is exact on linear data), nothing flagged.
!< - The same field scaled by 1.2: q = 0.2 everywhere, every node flagged.
!< - Unknown nodes (FAR_VALUE) and their neighbours: q = -1.
!< - A 2-D field computed from a 64-gon circle: outside the polygon, away from it, q is below 0.05
!<   (no medial axis outside a convex body).

use fossil, only : distance_field_object, curve_object, FAR_VALUE
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter        :: PI = acos(-1._R8P)
type(distance_field_object) :: f
type(curve_object)          :: c
real(R8P)                   :: x(21), y(17), z(13), n(3), u(64), v(64), q(21, 17, 13), qq(41, 41, 1), g(41), xc(3)
integer(I4P)                :: i, j, k, nf, nq
logical                     :: ok(4)

do i = 1, 21
   x(i) = -1._R8P + 2._R8P * ((i - 1) / 20._R8P)**1.3_R8P
enddo
do j = 1, 17
   y(j) = -0.5_R8P + ((j - 1) / 16._R8P)**0.8_R8P
enddo
do k = 1, 13
   z(k) = 0.3_R8P * (k - 1) / 12._R8P
enddo
call f%initialize(x=x, y=y, z=z)
n = [1._R8P, 2._R8P, -2._R8P] / 3._R8P
do k = 1, 13
   do j = 1, 17
      do i = 1, 21
         f%phi(i, j, k) = dot_product(n, [x(i), y(j), z(k)]) - 0.1_R8P
      enddo
   enddo
enddo
call f%gradient_diagnostic(q=q, flagged=nf)
ok(1) = maxval(q) < 1.e-13_R8P .and. nf == 0
print '(A,ES9.1,A,I0,A,L1)', 'plane, stretched grid: max q ', maxval(q), '  flagged ', nf, '  ok ', ok(1)
f%phi = 1.2_R8P * f%phi
call f%gradient_diagnostic(q=q, flagged=nf)
ok(2) = maxval(abs(q - 0.2_R8P)) < 1.e-13_R8P .and. nf == size(q)
print '(A,ES9.1,A,I0,A,L1)', 'plane x 1.2: max |q - 0.2| ', maxval(abs(q - 0.2_R8P)), '  flagged ', nf, '  ok ', ok(2)
f%phi(10, 8, 6) = FAR_VALUE
call f%gradient_diagnostic(q=q, flagged=nf)
ok(3) = q(10, 8, 6) == -1._R8P .and. q(11, 8, 6) == -1._R8P .and. q(10, 9, 6) == -1._R8P .and. q(10, 8, 7) == -1._R8P .and. &
        q(12, 8, 6) >= 0._R8P .and. count(q == -1._R8P) == 7
print '(A,I0,A,L1)', 'unknown node: q = -1 at it and its 6 neighbours, count ', count(q == -1._R8P), '  ok ', ok(3)
! 2-D field of a 64-gon circle of radius 0.5
do k = 1, 64
   u(k) = 0.5_R8P * cos(2 * PI * (k - 1) / 64) ; v(k) = 0.5_R8P * sin(2 * PI * (k - 1) / 64)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
g = [(-1._R8P + 2._R8P * (i - 1) / 40._R8P, i = 1, 41)]
call f%initialize(x=g, y=g)
call f%compute(curve=c)
call f%gradient_diagnostic(q=qq, flagged=nf)
nq = 0
do j = 2, 40
   do i = 2, 40
      xc = [g(i), g(j), 0._R8P]
      if (norm2(xc(1:2)) > 0.5_R8P + 0.1_R8P .and. qq(i, j, 1) > 0.05_R8P) nq = nq + 1
   enddo
enddo
ok(4) = nq == 0
print '(A,I0,A,I0,A,L1)', '64-gon field: outside nodes with q > 0.05: ', nq, '  (flagged overall ', nf, ')  ok ', ok(4)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_field_diagnostic
