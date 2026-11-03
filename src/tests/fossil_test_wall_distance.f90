!< FOSSIL test: wall distance with body id and cutoff (issue #42 step 2).
program fossil_test_wall_distance
!< FOSSIL test: wall distance with body id and cutoff (issue #42 step 2).
!<
!< - Two unit cubes 1.5 apart (two bodies), 2000 seeded points: without a cutoff the distances and
!<   bodies equal those of query_many bit for bit; with cutoff 0.3 they equal them below the cutoff
!<   and are (0.3, 0) at or beyond it; the device twin (on the host) equals the CPU; an octree surface
!<   (compute_distance fallback) gives the same distances and bodies.
!< - An annulus with an island (two bodies): the curve twin, same checks.

use fossil,                   only : surface_stl_object, curve_object, AABB_TREE_OCTREE
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

integer, parameter        :: N = 2000
type(body_t), allocatable :: b(:)
type(soup_t)              :: two
type(surface_stl_object)  :: s, so
type(curve_object)        :: c
type(vector_R8P)          :: q(N)
real(R8P)                 :: r(3), d0(N), d1(N), d2(N), d3(N), p2(2, N), co
integer(I4P)              :: b0(N), b1(N), b2(N), b3(N), i, k, nf, nbad
logical                   :: ok(2)

call random_seed(size=k)
call random_seed(put=[(99 + i, i = 1, k)])
call bench_bodies(b)
nf = size(b(1)%soup%v, 3)
allocate(two%v(3, 3, 2 * nf))
two%v(:, :, 1:nf) = b(1)%soup%v ; two%v(:, :, nf + 1:) = b(1)%soup%v ; two%v(1, :, nf + 1:) = two%v(1, :, nf + 1:) + 1.5_R8P
call to_surface(two, s) ; call to_surface(two, so) ; call so%analyze(aabb_tree_kind=AABB_TREE_OCTREE)
do i = 1, N
   call random_number(r) ; q(i) = vector_R8P(-0.5_R8P + 3.5_R8P * r(1), -0.5_R8P + 2 * r(2), -0.5_R8P + 2 * r(3))
enddo
co = 0.3_R8P
call s%query_many(points=q, distances=d0, bodies=b0, is_signed=.false.)
call s%wall_distance_many(points=q, distances=d1, bodies=b1)
nbad = count(d1 /= d0) + count(b1 /= b0)
call s%wall_distance_many(points=q, distances=d1, bodies=b1, cutoff=co)
do i = 1, N
   if (d0(i) < co) then
      if (d1(i) /= d0(i) .or. b1(i) /= b0(i)) nbad = nbad + 1
   else
      if (d1(i) /= co .or. b1(i) /= 0) nbad = nbad + 1
   endif
enddo
call s%enter_device
call s%wall_distance_many_device(points=q, distances=d2, bodies=b2, cutoff=co)
call s%exit_device
nbad = nbad + count(d2 /= d1) + count(b2 /= b1)
call so%wall_distance_many(points=q, distances=d3, bodies=b3, cutoff=co)
nbad = nbad + count(abs(d3 - d1) > 1.e-15_R8P) + count(b3 /= b1)
ok(1) = nbad == 0 .and. count(b1 == 1) > 0 .and. count(b1 == 2) > 0 .and. count(b1 == 0) > 0
print '(A,I0,A,3I6,A,L1)', 'two cubes: mismatches ', nbad, '  bodies 1/2/none ', count(b1 == 1), count(b1 == 2), count(b1 == 0), &
      '  ok ', ok(1)
! curve: annulus with an island
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P, 2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P, 2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
call random_number(p2) ; p2 = -1._R8P + 8._R8P * p2
call c%query_many(points=p2, distances=d0, bodies=b0, is_signed=.false.)
call c%wall_distance_many(points=p2, distances=d1, bodies=b1, cutoff=0.4_R8P)
nbad = 0
do i = 1, N
   if (d0(i) < 0.4_R8P) then
      if (d1(i) /= d0(i) .or. b1(i) /= b0(i)) nbad = nbad + 1
   else
      if (d1(i) /= 0.4_R8P .or. b1(i) /= 0) nbad = nbad + 1
   endif
enddo
ok(2) = nbad == 0 .and. count(b1 == 1) > 0 .and. count(b1 == 2) > 0
print '(A,I0,A,L1)', 'annulus with island: mismatches ', nbad, '  ok ', ok(2)
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_wall_distance
