!< FOSSIL test: device signs on the sign-benchmark points (issue #39 step 4).
program fossil_test_sign_bench_device
!< FOSSIL test: device signs on the sign-benchmark points (issue #39 step 4).
!<
!< `distance_many_device` and `query_many_device` against `distance_many` and `query_many` on every
!< point of the 3-D benchmark (9 bodies, delta 1e-1 ... 1e-12 L) and near the mid-height sections
!< of the bodies (2-D), pseudo-normal and fast-winding-number signs: the signed distances must be
!< equal bit for bit. In tests-gnu the device loops run on the host; the GPU benches print the same
!< count (`sign bench device vs serial mismatches`) for the manual gate.

use fossil_sign_bench_bodies, only : body_t, bench_bodies, device_sign_mismatches_3d, device_sign_mismatches_2d

implicit none

type(body_t), allocatable :: b(:)
integer                   :: n3, n2

call bench_bodies(b)
n3 = device_sign_mismatches_3d(b)
n2 = device_sign_mismatches_2d(b)
print '(A,I0)', '3-D sign bench device vs serial mismatches: ', n3
print '(A,I0)', '2-D sign bench device vs serial mismatches: ', n2
print '(A,L1)', 'Are all tests passed? ', n3 == 0 .and. n2 == 0
if (n3 /= 0 .or. n2 /= 0) error stop 1
endprogram fossil_test_sign_bench_device
