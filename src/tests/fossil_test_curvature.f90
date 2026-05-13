!< FOSSIL test: per-vertex Gaussian curvature (issue #18 §2.4).
!<
!< Four invariants on the public contract of `surface%gaussian_curvature`:
!<
!<   1. **Cube Gauss-Bonnet**: for a closed topological sphere,
!<      `sum_i K_i * A_i = 4 * pi`. The cube has 8 corner vertices each
!<      contributing a `pi / 2` angle defect; total = `4 * pi`. This is
!<      the topological invariant — independent of mesh density.
!<      Asserted to within `1e-9` (FP precision on the exact cube
!<      angles).
!<
!<   2. **Cube vertex sign**: every cube corner is convex, so K > 0
!<      everywhere on the cube. Asserts minval(K) > 0. The cube has only
!<      8 corner vertices (no interior face vertices, because each face
!<      is one diagonal-split quad).
!<
!<   3. **Sphere K ~ 1/R^2**: for a unit sphere, K = 1 everywhere
!<      analytically. Asserts `median(K)` is within 10% of 1.0 on the
!<      MC-generated unit sphere. Tessellation noise inflates outliers,
!<      so use median + loose tolerance.
!<
!<   4. **Bunny finiteness**: on the bunny fixture (~70k facets), every
!<      computed `K(i)` is finite (no NaN, no Inf). This is the basic
!<      "did the algorithm run without arithmetic disasters on a
!<      large real-world mesh" sanity test.
!<      (Gauss-Bonnet on the bunny would be the natural sphere-genus
!<      check but the bunny fixture turns out to NOT be watertight —
!<      `is_watertight() = .false.` — even though it passes the
!<      §1.4 winding-number test on dirty input. Same situation for
!<      the MC sphere via §1.5: documented Lorensen-Cline 105/150
!<      ambiguity leaves ~200 disconnected edges. Neither closed-
!<      manifold fixture is available in the test suite, so the cube
!<      test alone covers Gauss-Bonnet; the bunny test instead
!<      validates "didn't blow up on a large dense input.")

program fossil_test_curvature

use fossil, only : surface_stl_object, extract_isosurface, CURV_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere, bunny
type(facet_object), allocatable :: sphere_facets(:)
type(vertex_pool_object), pointer :: pool
real(R8P), allocatable :: K(:), area_third(:), values(:,:,:)
type(vector_R8P) :: bmin, bmax, p1, p2, p3
real(R8P) :: x, y, z, dx
real(R8P) :: gb_sum, K_median, K_min, K_max, area_2x, area_tri
real(R8P) :: cross_x, cross_y, cross_z
integer(I4P), parameter :: N_GRID = 32_I4P
integer(I4P) :: status, i, j, kk, f, v1, v2, v3, n_vertices
real(R8P), parameter :: PI = 4._R8P * atan(1._R8P)
real(R8P), parameter :: TOL_CUBE_GB = 1.0e-9_R8P
real(R8P), parameter :: TOL_SPHERE_REL = 0.10_R8P  ! median K within 10% of 1.0
real(R8P), parameter :: TOL_SPHERE_GB = 0.05_R8P   ! Gauss-Bonnet within 5% of 4*pi
logical :: are_tests_passed(4)

are_tests_passed = .false.

! ---- Invariant 1 & 2: cube.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%gaussian_curvature(K=K, status=status)
n_vertices = size(K, kind=I4P)
print '(A,I0,A,I0)', 'cube: status=', status, '  n_vertices=', n_vertices

! For Gauss-Bonnet, we need barycentric A_i to multiply by K_i. Compute
! it here from the facets directly — this is a TEST-internal sanity, not
! something the public API exposes.
allocate(area_third(n_vertices))
area_third = 0._R8P
pool => cube%get_vertex_pool()
do f = 1_I4P, cube%get_facets_number()
   v1 = pool%facet_vid(facet_id=f, local_v=1_I4P)
   v2 = pool%facet_vid(facet_id=f, local_v=2_I4P)
   v3 = pool%facet_vid(facet_id=f, local_v=3_I4P)
   p1 = pool%coord(vid=v1)
   p2 = pool%coord(vid=v2)
   p3 = pool%coord(vid=v3)
   ! |cross((p2-p1), (p3-p1))| = 2 * area
   block
      type(vector_R8P) :: e1, e2, c
      e1 = p2 - p1
      e2 = p3 - p1
      c = e1%crossproduct(rhs=e2)
      area_2x = sqrt(c%dotproduct(rhs=c))
   endblock
   area_tri = 0.5_R8P * area_2x
   area_third(v1) = area_third(v1) + area_tri / 3._R8P
   area_third(v2) = area_third(v2) + area_tri / 3._R8P
   area_third(v3) = area_third(v3) + area_tri / 3._R8P
enddo
gb_sum = 0._R8P
do i = 1_I4P, n_vertices
   gb_sum = gb_sum + K(i) * area_third(i)
enddo
print '(A,F12.6,A,F12.6)', 'inv 1 (cube GB): sum(K*A)=', gb_sum, '  target=', 4._R8P * PI
are_tests_passed(1) = (status == CURV_STATUS_OK) .and. (abs(gb_sum - 4._R8P * PI) < TOL_CUBE_GB)

K_min = minval(K)
K_max = maxval(K)
print '(A,F12.6,A,F12.6)', 'inv 2 (cube sign): K_min=', K_min, '  K_max=', K_max
are_tests_passed(2) = (K_min > 0._R8P)
deallocate(K, area_third)

! ---- Invariant 3: sphere from MC, median(K) ~ 1/R^2.
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
call sphere%gaussian_curvature(K=K, status=status)
n_vertices = size(K, kind=I4P)
print '(A,I0,A,I0)', 'sphere: status=', status, '  n_vertices=', n_vertices

K_median = median_inplace(K)
print '(A,F12.6,A,F12.6)', 'inv 3 (sphere K=1): median(K)=', K_median, '  target=', 1._R8P
are_tests_passed(3) = (status == CURV_STATUS_OK) .and. &
                      (abs(K_median - 1._R8P) < TOL_SPHERE_REL)
deallocate(K)

! ---- Invariant 4: bunny finiteness — every K(i) is a finite real.
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny%gaussian_curvature(K=K, status=status)
n_vertices = size(K, kind=I4P)
block
   integer(I4P) :: n_bad
   n_bad = 0_I4P
   do i = 1_I4P, n_vertices
      if (K(i) /= K(i)) n_bad = n_bad + 1_I4P  ! NaN check
      if (abs(K(i)) > huge(1._R8P) * 0.5_R8P) n_bad = n_bad + 1_I4P  ! Inf check
   enddo
   print '(A,I0,A,I0,A,I0)', 'inv 4 (bunny finite): nf=', bunny%get_facets_number(), &
         '  nv=', n_vertices, '  n_bad=', n_bad
   are_tests_passed(4) = (status == CURV_STATUS_OK) .and. (n_bad == 0_I4P)
endblock

print '(A,4L2)', 'per-case results: ', are_tests_passed
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

endprogram fossil_test_curvature
