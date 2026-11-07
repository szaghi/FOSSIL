!< FOSSIL test: quality metrics of curvilinear blocks (issue #71, step 9).
program fossil_test_block_quality
!< FOSSIL test: quality metrics of curvilinear blocks (issue #71, step 9).
!<
!< Analytic values (1e-12):
!< 1. sheared grid x = X + s Y (shear angle a = atan s), 2-D and 3-D: skewness a / 90 degrees, orthogonality cos a,
!<    aspect sqrt(1 + s^2), warp 0, min det J / mean det J = 1;
!< 2. stretched grid (spacings 1, 0.1, 0.5): aspect 10, skewness 0, orthogonality 1, warp 0;
!< 3. a unit cube with the vertex (1, 1, 1) lifted by d along z: warp acos(1 / (1 + d^2)) on its top face (the diagonal
!<    through the lifted vertex; the other gives acos(1 / sqrt(1 + 2 d^2))),
!<    the other metrics finite; min det J below 1.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter         :: RAD2DEG = 180._R8P / acos(-1._R8P)
type(curvilinear_block_t)    :: b
type(curvilinear_block_2d_t) :: q
real(R8P), allocatable       :: x(:,:,:), y(:,:,:), z(:,:,:), u(:,:), v(:,:)
real(R8P), allocatable       :: sk(:,:,:), ort(:,:,:), asp(:,:,:), wp(:,:,:), jm(:,:,:)
real(R8P), allocatable       :: sk2(:,:), ort2(:,:), asp2(:,:), jm2(:,:)
real(R8P)                    :: s, a, d, err
integer(I4P)                 :: i, j, k
logical                      :: all_ok

all_ok = .true.

! 1. sheared
s = 0.4_R8P ; a = atan(s) * RAD2DEG
allocate(x(6, 5, 4), y(6, 5, 4), z(6, 5, 4), u(6, 5), v(6, 5))
do k = 1, 4
   do j = 1, 5
      do i = 1, 6
         x(i, j, k) = 0.2_R8P * (i - 1) + s * 0.2_R8P * (j - 1) ; y(i, j, k) = 0.2_R8P * (j - 1) ; z(i, j, k) = 0.2_R8P * (k - 1)
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%quality(sk, ort, asp, wp, jm)
err = max(maxval(abs(sk - a / 90)), maxval(abs(ort - cos(a / RAD2DEG))), maxval(abs(asp - sqrt(1 + s * s))), &
          maxval(abs(wp)), maxval(abs(jm - 1)))
call check(err <= 1e-12_R8P, '3-D sheared: skewness, orthogonality, aspect, warp, min J', err)
u = x(:,:,1) ; v = y(:,:,1)
call q%initialize(u, v) ; call q%quality(sk2, ort2, asp2, jm2)
err = max(maxval(abs(sk2 - a / 90)), maxval(abs(ort2 - cos(a / RAD2DEG))), maxval(abs(asp2 - sqrt(1 + s * s))), &
          maxval(abs(jm2 - 1)))
call check(err <= 1e-12_R8P, '2-D sheared: skewness, orthogonality, aspect, min J', err)

! 2. stretched
do k = 1, 4
   do j = 1, 5
      do i = 1, 6
         x(i, j, k) = 1._R8P * (i - 1) ; y(i, j, k) = 0.1_R8P * (j - 1) ; z(i, j, k) = 0.5_R8P * (k - 1)
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%quality(sk, ort, asp, wp)
err = max(maxval(abs(sk)), maxval(abs(ort - 1)), maxval(abs(asp - 10)), maxval(abs(wp)))
call check(err <= 1e-12_R8P, '3-D stretched: aspect 10, orthogonal, unwarped', err)
u = x(:,:,1) ; v = y(:,:,1)
call q%initialize(u, v) ; call q%quality(sk2, ort2, asp2)
err = max(maxval(abs(sk2)), maxval(abs(ort2 - 1)), maxval(abs(asp2 - 10)))
call check(err <= 1e-12_R8P, '2-D stretched: aspect 10, orthogonal', err)

! 3. warped face
deallocate(x, y, z) ; allocate(x(2, 2, 2), y(2, 2, 2), z(2, 2, 2))
d = 0.3_R8P
do k = 1, 2
   do j = 1, 2
      do i = 1, 2
         x(i, j, k) = i - 1._R8P ; y(i, j, k) = j - 1._R8P ; z(i, j, k) = k - 1._R8P
      enddo
   enddo
enddo
z(2, 2, 2) = 1 + d
call b%initialize(x, y, z) ; call b%quality(sk, ort, asp, wp, jm)
err = abs(wp(1, 1, 1) - acos(1 / (1 + d * d)) * RAD2DEG)
call check(err <= 1e-12_R8P .and. jm(1, 1, 1) < 1 .and. jm(1, 1, 1) > 0, 'lifted vertex: warp angle, 0 < min J < 1', err)

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES11.3)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check
endprogram fossil_test_block_quality
