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
program fossil_test_deviation

use fossil, only : surface_stl_object, facet_object, extract_isosurface, DEV_STATUS_OK, DEV_STATUS_BAD_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
!$ use omp_lib, only : omp_set_num_threads

implicit none

type(surface_stl_object) :: a, b
real(R8P)                :: lower, lower2, mean, rms, e(4)
real(R8P), allocatable   :: pv(:), pv2(:)
integer(I4P)             :: status, status2, samples
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

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
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
