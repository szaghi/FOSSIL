!< FOSSIL test: surface deviation, the one-sided and symmetric Hausdorff distances (issue #60, roadmap #58 WP2).
!<
!< 1. Identical surfaces (cube-02 against itself): 0 at the vertices, round-off (1e-15) at the interior samples.
!< 2. cube-02 against the same cube translated by t = (0.3, 0, 0): h = 0.3 both ways, attained on the face x = 0 (at its
!<    vertices), so the sampled bound is exact.
!< 3. A marching-cubes unit sphere against itself scaled by 1.05 about the centre (resize): every vertex lies on the
!<    unit sphere, so d(x, B) <= 0.05 |x| <= 0.05, and h is above 0.98 * 0.05 (the facet planes of B).
!< 4. Asymmetric pair: the cube [0.4, 0.6]^3 inside the unit cube: h(small, unit) = 0.4, h(unit, small) = 0.4 sqrt(3)
!<    (the unit cube's corners), both exact on vertices; H = 0.4 sqrt(3).
!< 5. Quadrature: d is constant (0.3) on the face x = 0 and smaller elsewhere; on two parallel planes at distance 0.3
!<    (a unit square against the same square lifted by 0.3) the mean and the rms are 0.3 to round-off.
!< 6. The result does not depend on the thread count (1 and 4 OpenMP threads, bit for bit; only with OpenMP).
!< Certified upper bound (step 2):
!< 8. On cases 2, 4 and 5: lower <= h <= upper (up to the round-off of the distances, 1e-14) and upper - lower <= tolerance.
!< 9. The maximum inside a facet: the acute triangle (0,0,0), (3,0,0), (1,2.5,0) against three tiny triangles (size
!<    1e-7) at its corners. d is the distance to the nearest corner, largest at the circumcentre (1.5, 0.85, 0), where it
!<    is the circumradius R = sqrt(2.9725), away from every vertex, midpoint and sample: the sampled lower bound misses
!<    it by more than 1e-3, the certified bounds bracket it (within the 1e-7 of the tiny triangles).
!< 10. The work cap: max_triangles = 8 gives DEV_STATUS_NOT_CONVERGED with valid, wider bounds.
!< Device path (step 3):
!< 12. device=.true. without enter_device: DEV_STATUS_BAD_INPUT. After enter_device: lower, upper and per_vertex equal
!<     to the CPU path bit for bit (tolerance 1e-3: the comparison is bitwise, a tight tolerance only costs time) (in a host build the device code runs with its directives inert; the GPU bench
!<     fossil_test_distance_bench checks the device itself).
!< 11. triangle_rho on 2000 random triangles (acute and obtuse, fixed seed): the largest distance to the nearest corner
!<     over a 60 x 60 lattice of each triangle is <= rho and >= 0.98 rho (the lattice's resolution).
program fossil_test_deviation

use fossil, only : surface_stl_object, facet_object, extract_isosurface, DEV_STATUS_OK, DEV_STATUS_BAD_INPUT, &
                   DEV_STATUS_NOT_CONVERGED
use fossil_deviation, only : triangle_rho
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
!$ use omp_lib, only : omp_set_num_threads

implicit none

type(surface_stl_object) :: a, b
real(R8P)                :: lower, lower2, mean, rms, e(4), upper, tol, exact
real(R8P), allocatable   :: pv(:), pv2(:)
integer(I4P)             :: status, status2, samples, triangles
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. identical
call a%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call a%deviation(a, lower=lower, mean=mean, rms=rms, per_vertex=pv, samples=samples, status=status)
ok = status == DEV_STATUS_OK .and. lower <= 1e-15_R8P .and. mean <= 1e-15_R8P .and. rms <= 1e-15_R8P .and. all(pv == 0._R8P)
print '(A,ES10.3,A,I0,A,L1)', '1. identical: lower ', lower, '  samples ', samples, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. translated cube
b = a
call b%translate(delta=0.3_R8P * ex_R8P)
call b%analyze
call a%deviation(b, lower=lower, status=status)
call b%deviation(a, lower=lower2, status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. abs(lower - 0.3_R8P) <= 1e-15_R8P .and. &
     abs(lower2 - 0.3_R8P) <= 1e-15_R8P
print '(A,2ES23.15,A,L1)', '2. translated cube: h(a,b), h(b,a) ', lower, lower2, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. sphere against the sphere scaled by 1.05
call make_sphere(a, 33_I4P, 1._R8P)
b = a
call b%resize(factor=1.05_R8P * (ex_R8P + ey_R8P + ez_R8P))
call b%analyze
call a%deviation(b, lower=lower, status=status)
ok = status == DEV_STATUS_OK .and. lower <= 0.05_R8P + 1e-15_R8P .and. lower >= 0.98_R8P * 0.05_R8P
print '(A,ES23.15,A,L1)', '3. scaled sphere: h ', lower, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. asymmetric: a small cube inside the unit cube
call a%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)       ! [0, 1]^3
b = a
call b%resize(factor=0.2_R8P * (ex_R8P + ey_R8P + ez_R8P))
call b%translate(delta=0.4_R8P * (ex_R8P + ey_R8P + ez_R8P))
call b%analyze                                                                      ! [0.4, 0.6]^3
call b%deviation(a, lower=lower, status=status)
call a%deviation(b, lower=lower2, status=status2)
e(1) = abs(lower - 0.4_R8P) ; e(2) = abs(lower2 - 0.4_R8P * sqrt(3._R8P))
call a%hausdorff(b, lower=e(3), status=status)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. e(1) <= 1e-15_R8P .and. e(2) <= 1e-15_R8P .and. &
     e(3) == lower2
print '(A,2ES23.15,A,L1)', '4. small in unit: h(small,unit), h(unit,small) ', lower, lower2, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. quadrature: two parallel squares 0.3 apart
call make_square(a, 0._R8P)
call make_square(b, 0.3_R8P)
call a%deviation(b, lower=lower, mean=mean, rms=rms, status=status)
ok = status == DEV_STATUS_OK .and. abs(lower - 0.3_R8P) <= 1e-15_R8P .and. abs(mean - 0.3_R8P) <= 1e-15_R8P .and. &
     abs(rms - 0.3_R8P) <= 1e-15_R8P
print '(A,3ES23.15,A,L1)', '5. parallel squares: lower, mean, rms ', lower, mean, rms, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 6. thread count
call make_sphere(a, 25_I4P, 1._R8P)
call make_sphere(b, 21_I4P, 1.02_R8P)
!$ call omp_set_num_threads(1)
call a%deviation(b, lower=lower, mean=mean, per_vertex=pv, status=status)
!$ call omp_set_num_threads(4)
call a%deviation(b, lower=lower2, mean=e(1), per_vertex=pv2, status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. lower == lower2 .and. mean == e(1) .and. all(pv == pv2)
print '(A,L1)', '6. 1 and 4 threads bit for bit: ', ok
all_ok = all_ok .and. ok

! ---- empty input
call b%destroy
call a%deviation(b, lower=lower, status=status)
ok = status == DEV_STATUS_BAD_INPUT
print '(A,L1)', '7. empty other rejected: ', ok
all_ok = all_ok .and. ok

! ---- 8. certified bounds bracket the exact values
tol = 1.e-4_R8P
call a%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
b = a
call b%translate(delta=0.3_R8P * ex_R8P)
call b%analyze
call a%deviation(b, lower=lower, upper=upper, tolerance=tol, triangles=triangles, status=status)
ok = status == DEV_STATUS_OK .and. lower <= 0.3_R8P + 1e-14_R8P .and. 0.3_R8P <= upper .and. upper - lower <= tol
print '(A,2ES23.15,A,I0,A,L1)', '8. translated cube: lower, upper ', lower, upper, '  triangles ', triangles, '  ok ', ok
all_ok = all_ok .and. ok
b = a
call b%resize(factor=0.2_R8P * (ex_R8P + ey_R8P + ez_R8P))
call b%translate(delta=0.4_R8P * (ex_R8P + ey_R8P + ez_R8P))
call b%analyze
call a%deviation(b, lower=lower, upper=upper, tolerance=tol, status=status)
exact = 0.4_R8P * sqrt(3._R8P)
ok = status == DEV_STATUS_OK .and. lower <= exact + 1e-14_R8P .and. exact <= upper .and. upper - lower <= tol
call b%deviation(a, lower=lower2, upper=e(1), tolerance=tol, status=status2)
ok = ok .and. status2 == DEV_STATUS_OK .and. lower2 <= 0.4_R8P + 1e-14_R8P .and. 0.4_R8P <= e(1) .and. e(1) - lower2 <= tol
call a%hausdorff(b, lower=e(2), upper=e(3), tolerance=tol, status=status)
ok = ok .and. status == DEV_STATUS_OK .and. e(2) <= exact + 1e-14_R8P .and. exact <= e(3)
print '(A,4ES12.4,A,L1)', '8. small in unit: bounds both ways ', lower, upper, lower2, e(1), '  ok ', ok
all_ok = all_ok .and. ok
call make_square(a, 0._R8P)
call make_square(b, 0.3_R8P)
call a%deviation(b, lower=lower, upper=upper, tolerance=tol, status=status)
ok = status == DEV_STATUS_OK .and. lower <= 0.3_R8P + 1e-15_R8P .and. 0.3_R8P <= upper .and. upper - lower <= tol
print '(A,2ES23.15,A,L1)', '8. parallel squares: lower, upper ', lower, upper, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 9. maximum inside a facet
call make_acute(a)
call make_tiny_corners(b)
exact = sqrt(2.9725_R8P)
call a%deviation(b, lower=lower2, status=status2)
call a%deviation(b, lower=lower, upper=upper, tolerance=tol, triangles=triangles, status=status)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_OK .and. lower2 < exact - 1e-3_R8P .and. &
     lower <= exact + 1e-14_R8P .and. exact - 1e-7_R8P <= upper .and. upper - lower <= tol
print '(A,ES12.4,A,2ES23.15,A,I0,A,L1)', '9. interior maximum: sampled ', lower2, '  certified ', lower, upper, &
                                         '  triangles ', triangles, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 10. work cap
call a%deviation(b, lower=lower, upper=upper, tolerance=tol, max_triangles=8_I4P, status=status)
ok = status == DEV_STATUS_NOT_CONVERGED .and. lower <= exact + 1e-14_R8P .and. exact - 1e-7_R8P <= upper .and. &
     upper - lower > tol
print '(A,2ES12.4,A,L1)', '10. capped: lower, upper ', lower, upper, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 11. triangle_rho
ok = check_rho()
print '(A,L1)', '11. triangle_rho on random triangles: ', ok
all_ok = all_ok .and. ok

! ---- 12. device path
call make_sphere(a, 25_I4P, 1._R8P)
call make_sphere(b, 21_I4P, 1.02_R8P)
call a%deviation(b, lower=lower, upper=upper, per_vertex=pv, tolerance=1e-3_R8P, status=status)
call a%deviation(b, lower=lower2, upper=e(1), per_vertex=pv2, tolerance=1e-3_R8P, device=.true., status=status2)
ok = status == DEV_STATUS_OK .and. status2 == DEV_STATUS_BAD_INPUT
call b%enter_device
call a%deviation(b, lower=lower2, upper=e(1), per_vertex=pv2, tolerance=1e-3_R8P, device=.true., status=status2)
call b%exit_device
ok = ok .and. status2 == DEV_STATUS_OK .and. lower2 == lower .and. e(1) == upper .and. all(pv2 == pv)
print '(A,2ES23.15,A,L1)', '12. device path: lower, upper ', lower2, e(1), '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine make_acute(surf)
   !< The acute triangle (0,0,0), (3,0,0), (1,2.5,0): circumcentre (1.5, 0.85, 0), circumradius sqrt(2.9725).
   type(surface_stl_object), intent(inout) :: surf
   type(facet_object), allocatable         :: f(:)

   allocate(f(1))
   f(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; f(1)%vertex(2) = vector_R8P(3._R8P, 0._R8P, 0._R8P)
   f(1)%vertex(3) = vector_R8P(1._R8P, 2.5_R8P, 0._R8P)
   call f(1)%compute_metrix
   call surf%destroy
   call surf%adopt_facets(facets=f)
   endsubroutine make_acute

   subroutine make_tiny_corners(surf)
   !< Three triangles of size 1e-7 containing the corners of make_acute's triangle (each corner is one of their vertices):
   !< d(x, surf) is within 1e-7 of the distance from x to the nearest corner.
   type(surface_stl_object), intent(inout) :: surf
   type(facet_object), allocatable         :: f(:)
   real(R8P), parameter                    :: S = 1.e-7_R8P
   real(R8P)                               :: p(3,3)
   integer                                 :: k

   p = reshape([0._R8P, 0._R8P, 0._R8P, 3._R8P, 0._R8P, 0._R8P, 1._R8P, 2.5_R8P, 0._R8P], [3, 3])
   allocate(f(3))
   do k = 1, 3
      f(k)%vertex(1) = vector_R8P(p(1,k), p(2,k), p(3,k))
      f(k)%vertex(2) = vector_R8P(p(1,k) + S, p(2,k), p(3,k) + S)
      f(k)%vertex(3) = vector_R8P(p(1,k), p(2,k) + S, p(3,k) + S)
      call f(k)%compute_metrix
   enddo
   call surf%destroy
   call surf%adopt_facets(facets=f)
   endsubroutine make_tiny_corners

   function check_rho() result(good)
   !< triangle_rho against a fine lattice of random triangles (Park-Miller, fixed seed).
   logical    :: good
   real(R8P)  :: a(3), b(3), c(3), x(3), rho, m, r(9)
   integer    :: t, i, j, k, n
   integer(8) :: seed

   good = .true.
   seed = 4242_8
   n = 60
   do t = 1, 2000
      do k = 1, 9
         seed = mod(16807_8 * seed, 2147483647_8) ; r(k) = real(seed, R8P) / 2147483647._R8P
      enddo
      a = r(1:3) ; b = r(4:6) ; c = r(7:9)
      rho = triangle_rho(a, b, c)
      m = 0._R8P
      do i = 0, n
         do j = 0, n - i
            x = a + (i * (b - a) + j * (c - a)) / real(n, R8P)
            m = max(m, min(norm2(x - a), norm2(x - b), norm2(x - c)))
         enddo
      enddo
      good = good .and. m <= rho * (1._R8P + 1e-12_R8P) .and. m >= 0.98_R8P * rho
   enddo
   endfunction check_rho

   subroutine make_sphere(surf, n, r)
   !< Marching-cubes sphere of radius r on an n^3 grid over [-1.5, 1.5]^3.
   type(surface_stl_object), intent(inout) :: surf
   integer(I4P),             intent(in)    :: n
   real(R8P),                intent(in)    :: r
   real(R8P), allocatable                  :: values(:,:,:)
   type(facet_object), allocatable         :: facets(:)
   real(R8P)                               :: h
   integer(I4P)                            :: i, j, k, st

   allocate(values(n, n, n))
   h = 3._R8P / (n - 1)
   do k = 1, n ; do j = 1, n ; do i = 1, n
      values(i,j,k) = norm2([-1.5_R8P + (i - 1) * h, -1.5_R8P + (j - 1) * h, -1.5_R8P + (k - 1) * h]) - r
   enddo ; enddo ; enddo
   call extract_isosurface(values=values, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                           bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=st)
   call surf%destroy
   call surf%adopt_facets(facets=facets)
   endsubroutine make_sphere

   subroutine make_square(surf, z)
   !< The unit square [0, 1]^2 at height z, two facets.
   type(surface_stl_object), intent(inout) :: surf
   real(R8P),                intent(in)    :: z
   type(facet_object), allocatable         :: f(:)

   allocate(f(2))
   f(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, z) ; f(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, z)
   f(1)%vertex(3) = vector_R8P(1._R8P, 1._R8P, z)
   f(2)%vertex(1) = vector_R8P(0._R8P, 0._R8P, z) ; f(2)%vertex(2) = vector_R8P(1._R8P, 1._R8P, z)
   f(2)%vertex(3) = vector_R8P(0._R8P, 1._R8P, z)
   call f(1)%compute_metrix ; call f(2)%compute_metrix
   call surf%destroy
   call surf%adopt_facets(facets=f)
   endsubroutine make_square
endprogram fossil_test_deviation
