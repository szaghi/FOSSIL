!< FOSSIL test: `fossil_ib_generator` output against the analytic box distance (issue #32, step 7).
program fossil_test_ib_generator
!< FOSSIL test: `fossil_ib_generator` output against the analytic box distance (issue #32, step 7).
!<
!< Runs `exe/fossil_ib_generator` on cube-02 (unit cube) in the box [-1, 2]^3, 30 cells per
!< axis plus a 2-cell frame, reads the Xall `.grd` / `.ib` files back and compares the field
!< (positive inside, Xall convention) with the analytic signed distance of the cube at the
!< 34^3 cell centres:
!<  - every sign is right (cells off the surface);
!<  - the max error is below 2 h (h = 0.1): the far field comes from fast sweeping, whose
!<    error measured on this case at h = 0.125 / 0.0625 / 0.03125 is 1.10 h / 1.44 h / 1.80 h
!<    (fossil_test_distance_field_fsm); the former generator only reached cells it computed
!<    exactly, and copied boundary values beyond.
!< A second run with `--squared` must write sign(d) d**2 of the first run, bit for bit (the
!< former generator's output quantity).
!< The test needs the generator built in the same mode (it is, by every tests-* mode that
!< builds src/app) and is run from the repository root, like every FOSSIL test.

use penf, only : I4P, R8P

implicit none

character(*), parameter :: BASE = 'fossil_test_ibgen'
integer(I4P), parameter :: N = 30
real(R8P),    parameter :: H = 0.1_R8P
real(R8P), allocatable  :: xn(:,:,:), yn(:,:,:), zn(:,:,:), ib(:,:,:), ib2(:,:,:)
real(R8P)               :: t, c(3), d, err
integer(I4P)            :: exitstat, exitstat2, u, nb, n1, n2, n3, i, j, k, n_sign
logical                 :: ok

call execute_command_line('exe/fossil_ib_generator --stl src/tests/cube-02.stl --out '//BASE// &
                          ' --bmin -1 -1 -1 --bmax 2 2 2 --ni 30 --nj 30 --nk 30 --gi 2 --gj 2 2 --gk 2'// &
                          ' --sign_algorithm pseudo_normal > /dev/null', exitstat=exitstat)
ok = exitstat == 0

open(newunit=u, file=BASE//'.grd', form='unformatted', action='read', status='old')
read(u) nb
read(u) n1, n2, n3
allocate(xn(0:n1, 0:n2, 0:n3), yn(0:n1, 0:n2, 0:n3), zn(0:n1, 0:n2, 0:n3))
read(u) xn
read(u) yn
read(u) zn
close(u, status='delete')
open(newunit=u, file=BASE//'.ib', form='unformatted', action='read', status='old')
read(u) t
read(u) n1, n2, n3
allocate(ib(-1:n1 + 2, -1:n2 + 2, -1:n3 + 2))
read(u) ib
close(u, status='delete')
open(newunit=u, file=BASE//'.vts', status='old')
close(u, status='delete')
ok = ok .and. nb == 1 .and. n1 == N .and. n2 == N .and. n3 == N .and. &
     abs(xn(0, 0, 0) + 1._R8P) < 1.e-14_R8P .and. abs(xn(N, 0, 0) - 2._R8P) < 1.e-14_R8P

err = 0._R8P ; n_sign = 0
do k = -1, N + 2
   do j = -1, N + 2
      do i = -1, N + 2
         c = -1._R8P + (real([i, j, k], R8P) - 0.5_R8P) * H
         d = box_sdf(c)
         err = max(err, abs(-ib(i, j, k) - d))
         if (abs(d) > 1.e-12_R8P .and. ((d < 0._R8P) .neqv. (ib(i, j, k) > 0._R8P))) n_sign = n_sign + 1
      enddo
   enddo
enddo
ok = ok .and. n_sign == 0 .and. err < 2._R8P * H
print '(A,I0,A,I0,A,ES10.3,A,F6.3,A)', 'generator exit ', exitstat, '  sign errors ', n_sign, '  max |err| ', err, &
      '  (', err / H, ' h)'

! --squared: sign(d) d**2 of the same field
call execute_command_line('exe/fossil_ib_generator --stl src/tests/cube-02.stl --out '//BASE// &
                          ' --bmin -1 -1 -1 --bmax 2 2 2 --ni 30 --nj 30 --nk 30 --gi 2 --gj 2 2 --gk 2'// &
                          ' --sign_algorithm pseudo_normal --squared > /dev/null', exitstat=exitstat2)
allocate(ib2, mold=ib)
open(newunit=u, file=BASE//'.ib', form='unformatted', action='read', status='old')
read(u) t
read(u) n1, n2, n3
read(u) ib2
close(u, status='delete')
open(newunit=u, file=BASE//'.grd', status='old')
close(u, status='delete')
open(newunit=u, file=BASE//'.vts', status='old')
close(u, status='delete')
ok = ok .and. exitstat2 == 0 .and. all(ib2 == sign(ib * ib, ib))
print '(A,I0,A,L1)', '--squared: exit ', exitstat2, '  ib == sign(d) d**2 of the first run: ', all(ib2 == sign(ib * ib, ib))
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   pure function box_sdf(p) result(dist)
   !< Signed distance to the unit cube [0, 1]^3 (negative inside).
   real(R8P), intent(in) :: p(3) !< Point.
   real(R8P)             :: dist !< Distance.
   real(R8P)             :: q(3) !< |p - centre| - half side.

   q = abs(p - 0.5_R8P) - 0.5_R8P
   dist = sqrt(sum(max(q, 0._R8P)**2)) + min(maxval(q), 0._R8P)
   endfunction box_sdf
endprogram fossil_test_ib_generator
