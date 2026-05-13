!< FOSSIL test: cotangent Laplacian + barycentric mass matrix (issue #18 §2.1).
!<
!< Five invariants on the public contract of `surface%cotangent_laplacian`:
!<
!<   1. **Constant kernel**: `L * 1 = 0` on every vertex. This is the
!<      defining property of the Laplacian and is automatic if the
!<      diagonal `L_ii = -sum_{j ≠ i} L_ij` was assembled correctly.
!<      Asserted to within `1e-9` (tight FP tolerance) on the cube.
!<
!<   2. **Symmetry**: `L_ij = L_ji` for every nonzero entry. The cotangent
!<      Laplacian is symmetric by construction, but a bookkeeping bug in
!<      the per-triangle assembly would break this. Asserted via
!<      `L%is_symmetric(tol=1e-12)` on the cube.
!<
!<   3. **Mass conservation**: `sum(M_diag) = total_surface_area`.
!<      Barycentric mass distributes one third of each triangle's area to
!<      each of its three vertices; summing the diagonal must reproduce
!<      `sum(triangle areas)`. Asserted to within `1e-9` of the
!<      surface's `get_area()` result.
!<
!<   4. **Row sums on the cube**: every row sum of `L` is zero (a
!<      consequence of invariant 1, but a stronger statement: it holds
!<      per-row, not just on the global sum). Asserted via the CSR
!<      `row_sum` accessor.
!<
!<   5. **Sphere mean-curvature check**: for a unit sphere built via
!<      §1.5 marching cubes, the mean-curvature-normal magnitude
!<      `||L V|| / (2 M_ii)` should be approximately `1` (= 1/R) at
!<      every vertex. Tessellation noise from MC means we can't expect
!<      tight tolerance — assert that the MEDIAN magnitude is within
!<      30% of 1 on the unit sphere. This is the "is the Laplacian
!<      actually doing what a Laplacian does" sanity test.

program fossil_test_laplacian

use fossil, only : surface_stl_object, csr_matrix_t, extract_isosurface, &
                   LAPL_STATUS_OK
use fossil_facet_object, only : facet_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere
type(facet_object), allocatable :: sphere_facets(:)
type(csr_matrix_t) :: L, M
real(R8P), allocatable :: ones(:), Lx(:), Vx(:), Vy(:), Vz(:), LVx(:), LVy(:), LVz(:), Hmag(:)
real(R8P), allocatable :: values(:,:,:)
type(vector_R8P) :: bmin, bmax
real(R8P) :: x, y, z, dx, mass_sum, row_sum_max, area_total
real(R8P) :: H_median, H_mean, M_ii
integer(I4P), parameter :: N_GRID = 32_I4P
integer(I4P) :: status, i, j, k, n_vertices, n_valid
real(R8P), parameter :: TOL_KERNEL = 1.0e-9_R8P
real(R8P), parameter :: TOL_SYMMETRY = 1.0e-12_R8P
real(R8P), parameter :: TOL_MASS = 1.0e-9_R8P
real(R8P), parameter :: TOL_CURVATURE_REL = 0.30_R8P
logical :: are_tests_passed(5)

are_tests_passed = .false.

! ---- Build the matrices on the cube fixture.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%cotangent_laplacian(L=L, M=M, status=status)
n_vertices = L%get_nrows()
print '(A,I0,A,I0,A,I0)', 'cube: n_vertices=', n_vertices, &
      '  L_nnz=', L%get_nnz(), '  M_nnz=', M%get_nnz()

! ---- Invariant 1: L * 1 = 0.
allocate(ones(n_vertices))
allocate(Lx(n_vertices))
ones = 1._R8P
call L%multiply_vector(x=ones, y=Lx)
print '(A,ES12.4)', 'inv 1 (L*1=0):     max|Lx| = ', maxval(abs(Lx))
are_tests_passed(1) = (status == LAPL_STATUS_OK) .and. (maxval(abs(Lx)) < TOL_KERNEL)

! ---- Invariant 2: symmetry.
print '(A,L1)', 'inv 2 (symmetric):  ', L%is_symmetric(tol=TOL_SYMMETRY)
are_tests_passed(2) = L%is_symmetric(tol=TOL_SYMMETRY)

! ---- Invariant 3: sum(M_diag) = total_area.
mass_sum = 0._R8P
do i = 1_I4P, n_vertices
   mass_sum = mass_sum + M%row_sum(row=i)
enddo
area_total = cube%get_area()
print '(A,F12.6,A,F12.6)', 'inv 3 (mass=area):  mass_sum=', mass_sum, '  area=', area_total
are_tests_passed(3) = (abs(mass_sum - area_total) < TOL_MASS)

! ---- Invariant 4: every row of L sums to zero.
row_sum_max = 0._R8P
do i = 1_I4P, n_vertices
   row_sum_max = max(row_sum_max, abs(L%row_sum(row=i)))
enddo
print '(A,ES12.4)', 'inv 4 (per-row=0):  max|row_sum| = ', row_sum_max
are_tests_passed(4) = (row_sum_max < TOL_KERNEL)

deallocate(ones, Lx)

! ---- Invariant 5: sphere mean-curvature ~= 1/R.
! Build sphere via MC.
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)
do k = 1_I4P, N_GRID
   z = -2._R8P + (k - 1) * dx
   do j = 1_I4P, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1_I4P, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1._R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)
call sphere%cotangent_laplacian(L=L, M=M, status=status)
n_vertices = L%get_nrows()

! Pack per-vertex coordinates into three arrays.
allocate(Vx(n_vertices), Vy(n_vertices), Vz(n_vertices))
allocate(LVx(n_vertices), LVy(n_vertices), LVz(n_vertices))
allocate(Hmag(n_vertices))
block
   use fossil_vertex_pool_object, only : vertex_pool_object
   type(vertex_pool_object), pointer :: pool
   type(vector_R8P) :: p
   pool => sphere%get_vertex_pool()
   do i = 1_I4P, n_vertices
      p = pool%coord(vid=i)
      Vx(i) = p%x
      Vy(i) = p%y
      Vz(i) = p%z
   enddo
endblock

call L%multiply_vector(x=Vx, y=LVx)
call L%multiply_vector(x=Vy, y=LVy)
call L%multiply_vector(x=Vz, y=LVz)

! Per-vertex mean-curvature magnitude: ||LV|| / (2 M_ii). Skip vertices
! with vanishing mass (defensive — shouldn't happen on a clean sphere).
n_valid = 0_I4P
do i = 1_I4P, n_vertices
   M_ii = M%row_sum(row=i)
   if (M_ii <= 0._R8P) cycle
   n_valid = n_valid + 1_I4P
   Hmag(n_valid) = sqrt(LVx(i)**2 + LVy(i)**2 + LVz(i)**2) / (2._R8P * M_ii)
enddo

! Median via copy-sort.
H_median = median_inplace(Hmag(1:n_valid))
H_mean = sum(Hmag(1:n_valid)) / real(n_valid, R8P)
print '(A,F8.4,A,F8.4,A,I0,A,I0)', 'inv 5 (sphere H):   median=', H_median, &
      '  mean=', H_mean, '  valid=', n_valid, '/', n_vertices
are_tests_passed(5) = (status == LAPL_STATUS_OK) .and. &
                      (abs(H_median - 1._R8P) < TOL_CURVATURE_REL)

deallocate(Vx, Vy, Vz, LVx, LVy, LVz, Hmag)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains

   function median_inplace(arr) result(m)
   !< Median via insertion sort on a copy. n is small (~few thousand), fine.
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

endprogram fossil_test_laplacian
