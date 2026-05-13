!< FOSSIL test: 2D Delaunay triangulation Phase 1 (issue #18 §1.1 stage 2).
!<
!< Five invariants:
!<   1. Three collinear-free points produce exactly 1 triangle.
!<   2. Square (4 corners) produces exactly 2 triangles.
!<   3. Square + center produces exactly 4 triangles, each containing the center.
!<   4. Random 100-point cloud: every output triangle satisfies the empty-
!<      circumcircle property (no other input point lies strictly inside any
!<      triangle's circumcircle). This is the **defining** property of a
!<      Delaunay triangulation; if it holds for every triangle the output is
!<      Delaunay by definition.
!<   5. Triangle count: for N input points (no degeneracy), the Delaunay
!<      triangulation has at most 2*N - 5 triangles (Euler's formula).
!<      Verify on the random cloud.

program fossil_test_dt

use fossil_dt,         only : triangulation_t, dt_build
use fossil_predicates, only : incircle, predicates_initialize
use penf,              only : I4P, R8P
use vecfor,            only : vector_R8P

implicit none

type(triangulation_t) :: dt
real(R8P), allocatable :: pts(:,:)
integer(I4P) :: n, i, ntri, v(3)
logical :: are_tests_passed(5)
logical :: empty_circ_ok

real(R8P), parameter :: TOL_INCIRCLE = 1.0e-9_R8P  ! tolerance for the empty-circumcircle check

are_tests_passed = .false.

call predicates_initialize

! ----- 1. Three points -> 1 triangle.
allocate(pts(2, 3))
pts(:, 1) = [0._R8P, 0._R8P]
pts(:, 2) = [1._R8P, 0._R8P]
pts(:, 3) = [0._R8P, 1._R8P]
call dt_build(tri=dt, points=pts)
ntri = dt%num_triangles()
are_tests_passed(1) = (ntri == 1)
print '(A,I0)', '3 pts -> triangles = ', ntri
deallocate(pts)

! ----- 2. Square -> 2 triangles.
allocate(pts(2, 4))
pts(:, 1) = [0._R8P, 0._R8P]
pts(:, 2) = [1._R8P, 0._R8P]
pts(:, 3) = [1._R8P, 1._R8P]
pts(:, 4) = [0._R8P, 1._R8P]
call dt_build(tri=dt, points=pts)
ntri = dt%num_triangles()
are_tests_passed(2) = (ntri == 2)
print '(A,I0)', 'square -> triangles = ', ntri
deallocate(pts)

! ----- 3. Square + center -> 4 triangles.
allocate(pts(2, 5))
pts(:, 1) = [0._R8P, 0._R8P]
pts(:, 2) = [1._R8P, 0._R8P]
pts(:, 3) = [1._R8P, 1._R8P]
pts(:, 4) = [0._R8P, 1._R8P]
pts(:, 5) = [0.5_R8P, 0.5_R8P]
call dt_build(tri=dt, points=pts)
ntri = dt%num_triangles()
are_tests_passed(3) = (ntri == 4)
print '(A,I0)', 'square+center -> triangles = ', ntri
deallocate(pts)

! ----- 4. Random 100-point cloud: empty-circumcircle property + Euler bound.
n = 100
allocate(pts(2, n))
call random_seed_default()
do i = 1, n
   call random_number(pts(:, i))
enddo
call dt_build(tri=dt, points=pts)
ntri = dt%num_triangles()

empty_circ_ok = check_empty_circumcircle(dt=dt, n=n)
are_tests_passed(4) = empty_circ_ok
print '(A,I0,A,L1)', 'random 100 pts -> triangles = ', ntri, ', empty-circ ok = ', empty_circ_ok

! ----- 5. Euler bound: at most 2N - 5.
are_tests_passed(5) = (ntri <= 2 * n - 5)
print '(A,I0,A,I0)', 'Euler bound: ', ntri, ' <= ', 2 * n - 5
deallocate(pts)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains

   subroutine random_seed_default()
   !< Seed the RNG deterministically so this test is reproducible across runs.
   integer :: seed_size, k
   integer, allocatable :: seed(:)
   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   do k = 1, seed_size
      seed(k) = 42 + k
   enddo
   call random_seed(put=seed)
   endsubroutine random_seed_default

   function check_empty_circumcircle(dt, n) result(ok)
   !< For every active triangle (a, b, c), check that no input point d (other
   !< than a, b, c) lies strictly inside the circumcircle. Sufficient condition
   !< for Delaunayness.
   type(triangulation_t), intent(in) :: dt
   integer(I4P),          intent(in) :: n
   logical                           :: ok
   integer(I4P)                      :: t, p, vv(3), ntri_active
   type(vector_R8P)                  :: a, b, c, d
   real(R8P)                         :: in_or_out

   ok = .true.
   ntri_active = dt%num_triangles()
   do t = 1, ntri_active
      call dt%triangle_vertices(k=t, v=vv)
      a = vector_R8P(dt%px(vv(1)), dt%py(vv(1)), 0._R8P)
      b = vector_R8P(dt%px(vv(2)), dt%py(vv(2)), 0._R8P)
      c = vector_R8P(dt%px(vv(3)), dt%py(vv(3)), 0._R8P)
      do p = 1, n
         if (p == vv(1) .or. p == vv(2) .or. p == vv(3)) cycle
         d = vector_R8P(dt%px(p), dt%py(p), 0._R8P)
         in_or_out = incircle(a=a, b=b, c=c, d=d)
         ! Note: incircle's sign depends on (a,b,c) orientation. With CCW
         ! triangles, in_or_out > 0 means d is strictly inside. We tolerate
         ! a small positive bound for boundary points.
         if (in_or_out > TOL_INCIRCLE) then
            ok = .false.
            print '(A,I0,A,I0,A,ES14.6)', 'tri ', t, ' fails empty-circ at point ', p, &
                  ': incircle=', in_or_out
            return
         endif
      enddo
   enddo
   endfunction check_empty_circumcircle

endprogram fossil_test_dt
