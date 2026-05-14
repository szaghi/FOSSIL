!< FOSSIL test: per-vertex signed mean curvature (issue #18 §2.4 mean
!< variant — Gaussian variant tested separately in
!< `fossil_test_curvature`).
!<
!< Three invariants on the public contract of `surface%mean_curvature`:
!<
!<   1. **Sphere H ~ 1/R = 1**: unit sphere via §1.5 MC, `H_i = 1`
!<      analytically (positive sign by convention — sphere surface is
!<      convex outward). Asserts `median(H)` within 15% of 1.0; loose
!<      tolerance because MC tessellation introduces vertex-position
!<      noise that the discrete formula amplifies more than the
!<      angle-defect Gaussian. A tighter check is invariant 2.
!<   2. **Sphere H sign**: at least 90% of valid vertices have H > 0
!<      (convex-outward sphere). The positive sign comes from the
!<      H*n vector projecting onto the per-vertex pseudo-normal;
!<      this catches a regression where the sign convention flips.
!<   3. **Bunny finiteness**: 70k facets, 35k vertices, no NaN /
!<      Inf in any H value. Same large-real-world sanity check as
!<      the Gaussian test.
!<
!< Cube is a poor fixture for mean-curvature testing — its 8 corners
!< are sharp features where H is dominated by the discretization
!< artifact rather than by smooth curvature, and cube has no
!< flat-region interior vertices to test "H ≈ 0 on flat regions".
!< Sphere + bunny are sufficient.

program fossil_test_mean_curvature

use fossil, only : surface_stl_object, extract_isosurface, CURV_STATUS_OK
use fossil_facet_object, only : facet_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: sphere, bunny
type(facet_object), allocatable :: sphere_facets(:)
real(R8P), allocatable :: H(:), values(:,:,:)
type(vector_R8P) :: bmin, bmax
real(R8P) :: x, y, z, dx
real(R8P) :: H_median
integer(I4P), parameter :: N_GRID = 32_I4P
integer(I4P) :: status, i, j, kk, n_pos, n_total, n_bad
real(R8P), parameter :: TOL_SPHERE_REL = 0.15_R8P
real(R8P), parameter :: SIGN_FRACTION_MIN = 0.90_R8P
logical :: are_tests_passed(3)

are_tests_passed = .false.

! ---- Build sphere via MC.
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)
do kk = 1_I4P, N_GRID
   z = -2._R8P + (kk - 1) * dx
   do j = 1_I4P, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1_I4P, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, kk) = sqrt(x**2 + y**2 + z**2) - 1._R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)

! ---- Inv 1: sphere median H ~ 1/R = 1.
call sphere%mean_curvature(H=H, status=status)
H_median = median_inplace(H)
print '(A,I0,A,F8.4)', 'inv 1 (sphere H median): nv=', size(H), '  median(H)=', H_median
are_tests_passed(1) = (status == CURV_STATUS_OK) .and. &
                      (abs(H_median - 1._R8P) < TOL_SPHERE_REL)

! ---- Inv 2: at least 90% of sphere vertices have H > 0 (convex outward).
n_total = size(H, kind=I4P)
n_pos = count(H > 0._R8P)
print '(A,I0,A,I0,A,F6.3)', 'inv 2 (sphere sign): n_positive=', n_pos, ' / ', n_total, &
      '  fraction=', real(n_pos, R8P) / real(n_total, R8P)
are_tests_passed(2) = (real(n_pos, R8P) / real(n_total, R8P) >= SIGN_FRACTION_MIN)
deallocate(H)

! ---- Inv 3: bunny finiteness.
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny%mean_curvature(H=H, status=status)
n_bad = 0_I4P
do i = 1_I4P, size(H, kind=I4P)
   if (H(i) /= H(i)) n_bad = n_bad + 1_I4P                   ! NaN
   if (abs(H(i)) > huge(1._R8P) * 0.5_R8P) n_bad = n_bad + 1_I4P  ! Inf
enddo
print '(A,I0,A,I0,A,I0)', 'inv 3 (bunny finite): nf=', bunny%get_facets_number(), &
      '  nv=', size(H), '  n_bad=', n_bad
are_tests_passed(3) = (status == CURV_STATUS_OK) .and. (n_bad == 0_I4P)

print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains

   function median_inplace(arr) result(m)
   real(R8P), intent(in) :: arr(:)
   real(R8P)             :: m
   real(R8P), allocatable :: tmp(:)
   real(R8P)              :: key
   integer(I4P)           :: n, ii, jj

   n = size(arr, kind=I4P)
   if (n == 0_I4P) then
      m = 0._R8P
      return
   endif
   allocate(tmp(n))
   tmp = arr
   do ii = 2_I4P, n
      key = tmp(ii)
      jj = ii - 1_I4P
      do while (jj >= 1_I4P)
         if (tmp(jj) <= key) exit
         tmp(jj + 1_I4P) = tmp(jj)
         jj = jj - 1_I4P
      enddo
      tmp(jj + 1_I4P) = key
   enddo
   if (mod(n, 2_I4P) == 1_I4P) then
      m = tmp((n + 1_I4P) / 2_I4P)
   else
      m = 0.5_R8P * (tmp(n / 2_I4P) + tmp(n / 2_I4P + 1_I4P))
   endif
   deallocate(tmp)
   endfunction median_inplace

endprogram fossil_test_mean_curvature
