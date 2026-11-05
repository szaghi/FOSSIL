!< FOSSIL test: distance to the nearest feature edge (issue #61, roadmap #58 WP3, step 3).
!<
!< 1. Before build_features: query_many(feature_distances=) gives STATUS_INVALID_INPUT.
!< 2. cube-01, closed-form points: 0.5 from the edge y = z = 0 at (0.5, -0.3, -0.4); 0.5 from the face centre (0.5, 0.5,
!<    1); sqrt(0.29) from the corner (1, 1, 1) at (1.3, 1.4, 1.2); 0 on an edge (1e-15).
!< 3. 500 points (Park-Miller, fixed seed) in [-0.5, 1.5]^3: equal bit for bit to a brute force over every segment with
!<    the same clamped-projection formula.
!< 4. With a pose: the world-frame query equals the query of the body-frame points.
!< 5. query_many_device(feature_distances=) after enter_device equals the CPU (host build: device code path).
program fossil_test_feature_distance

use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, to_world, STATUS_OK, STATUS_INVALID_INPUT
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

integer(I4P), parameter           :: NP = 500
type(surface_stl_object)          :: s
type(vertex_pool_object), pointer :: pool
type(rigid_pose_t)                :: pose
type(vector_R8P)                  :: q(NP), qw(NP), pa, pb, p4(4)
real(R8P)                         :: d(NP), fd(NP), fd2(NP), ref(NP), d4(4), f4(4), x(3), e(3), w(3), t, best
integer(I4P), allocatable         :: edges(:,:)
integer(I4P)                      :: status, i, k
integer(8)                        :: seed
logical                           :: ok, all_ok

all_ok = .true.
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)

! ---- 1. not built
call s%query_many(points=q(1:1), distances=d(1:1), feature_distances=fd(1:1), status=status)
ok = status == STATUS_INVALID_INPUT
print '(A,L1)', '1. refused before build_features: ', ok
all_ok = all_ok .and. ok

! ---- 2. closed forms
call s%build_features()
p4 = [vector_R8P(0.5_R8P, -0.3_R8P, -0.4_R8P), vector_R8P(0.5_R8P, 0.5_R8P, 1._R8P), &
      vector_R8P(1.3_R8P, 1.4_R8P, 1.2_R8P), vector_R8P(1._R8P, 0.3_R8P, 0._R8P)]
call s%query_many(points=p4, distances=d4, feature_distances=f4, status=status)
ok = status == STATUS_OK .and. abs(f4(1) - 0.5_R8P) <= 1e-15_R8P .and. abs(f4(2) - 0.5_R8P) <= 1e-15_R8P .and. &
     abs(f4(3) - sqrt(0.29_R8P)) <= 1e-15_R8P .and. f4(4) <= 1e-15_R8P
print '(A,4ES12.4,A,L1)', '2. closed forms: ', f4, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. brute force
seed = 777_8
do i = 1, NP
   do k = 1, 3
      seed = mod(16807_8 * seed, 2147483647_8) ; x(k) = -0.5_R8P + 2._R8P * real(seed, R8P) / 2147483647._R8P
   enddo
   q(i) = vector_R8P(x(1), x(2), x(3))
enddo
call s%query_many(points=q, distances=d, feature_distances=fd, status=status)
call s%feature_edges(edges=edges)
pool => s%get_vertex_pool()
do i = 1, NP
   best = huge(1._R8P)
   do k = 1, size(edges, 2)
      pa = pool%coord(edges(1, k)) ; pb = pool%coord(edges(2, k))
      e = [pb%x - pa%x, pb%y - pa%y, pb%z - pa%z] ; w = [q(i)%x - pa%x, q(i)%y - pa%y, q(i)%z - pa%z]
      t = e(1) * e(1) + e(2) * e(2) + e(3) * e(3)
      t = max(0._R8P, min(1._R8P, (w(1) * e(1) + w(2) * e(2) + w(3) * e(3)) / t))
      w = w - t * e
      best = min(best, w(1) * w(1) + w(2) * w(2) + w(3) * w(3))
   enddo
   ref(i) = sqrt(best)
enddo
ok = status == STATUS_OK .and. all(fd == ref)
print '(A,L1)', '3. 500 points equal to brute force bit for bit: ', ok
all_ok = all_ok .and. ok

! ---- 4. pose
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.2_R8P, -0.1_R8P, 0.3_R8P])
do i = 1, NP
   x = to_world(pose, [q(i)%x, q(i)%y, q(i)%z]) ; qw(i) = vector_R8P(x(1), x(2), x(3))
enddo
call s%query_many(points=qw, distances=d, feature_distances=fd2, pose=pose, status=status)
ok = status == STATUS_OK .and. maxval(abs(fd2 - fd)) <= 1e-14_R8P
print '(A,ES10.3,A,L1)', '4. posed query against the body-frame one: ', maxval(abs(fd2 - fd)), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. device path
call s%enter_device
call s%query_many_device(points=q, distances=d, feature_distances=fd2, status=status)
call s%exit_device
ok = status == STATUS_OK .and. all(fd2 == fd)
print '(A,L1)', '5. device path equals the CPU: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok
endprogram fossil_test_feature_distance
