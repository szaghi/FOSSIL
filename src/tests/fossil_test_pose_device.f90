!< FOSSIL test: body-frame queries on the device (issue #41 step 2).
program fossil_test_pose_device
!< FOSSIL test: body-frame queries on the device (issue #41 step 2).
!<
!< `distance_many_device` and `query_many_device` with a pose (after `enter_device`; on the host in
!< tests-gnu) against `distance_many` and `query_many` with the same pose, every output bit for bit,
!< pseudo-normal and fast-winding-number signs: the cube, the cone and the L-prism of the sign
!< benchmark under a rotation with a large translation; a NACA 0012 and an annulus under in-plane
!< poses.

use fossil,                   only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, &
                                     pose_2d, to_world, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
type(rigid_pose_t)        :: p
real(R8P)                 :: u(200), v(200)
integer(I4P)              :: k, nbad3, nbad2

call bench_bodies(b)
p = pose_from_axis_angle([1._R8P, -2._R8P, 0.5_R8P], 1.1_R8P, [20._R8P, 3._R8P, -7._R8P])
nbad3 = 0
do k = 1, 3
   call to_surface(b(merge(1, merge(7, 8, k == 2), k == 1))%soup, s)
   nbad3 = nbad3 + check3(s, p)
enddo
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
nbad2 = check2(c, pose_2d(0.9_R8P, [5._R8P, -2._R8P]))
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
nbad2 = nbad2 + check2(c, pose_2d(-0.4_R8P, [30._R8P, 7._R8P]))
print '(A,I0,A,I0)', 'posed device vs serial mismatches: 3-D ', nbad3, '  2-D ', nbad2
print '(A,L1)', 'Are all tests passed? ', nbad3 == 0 .and. nbad2 == 0
if (nbad3 /= 0 .or. nbad2 /= 0) error stop 1

contains
   function check3(sf, pose) result(nbad)
   !< Posed device vs posed CPU, every output, both signs.
   type(surface_stl_object), intent(inout) :: sf   !< Surface.
   type(rigid_pose_t),       intent(in)    :: pose !< Pose.
   integer                                 :: nbad !< Mismatches.
   integer, parameter        :: N = 500 !< Points.
   integer(I4P), parameter   :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   type(vector_R8P) :: q(N) !< World points.
   real(R8P)        :: r(3), w(3), d1(N), d2(N), c1(3, N), c2(3, N), n1(3, N), n2(3, N) !< Work.
   integer(I4P)     :: f1(N), f2(N), t1(N), t2(N), i1(N), i2(N), b1(N), b2(N), i, a !< Work.

   do i = 1, N
      call random_number(r) ; w = to_world(pose, -0.5_R8P + 2._R8P * r) ; q(i) = vector_R8P(w(1), w(2), w(3))
   enddo
   nbad = 0
   call sf%enter_device
   do a = 1, 2
      call sf%distance_many(points=q, distances=d1, is_signed=.true., sign_algorithm=ALGO(a), pose=pose)
      call sf%distance_many_device(points=q, distances=d2, is_signed=.true., sign_algorithm=ALGO(a), pose=pose)
      nbad = nbad + count(d1 /= d2)
      call sf%query_many(points=q, distances=d1, closest_points=c1, facets=f1, features=t1, feature_indexes=i1, normals=n1, &
                         bodies=b1, sign_algorithm=ALGO(a), pose=pose)
      call sf%query_many_device(points=q, distances=d2, closest_points=c2, facets=f2, features=t2, feature_indexes=i2, &
                                normals=n2, bodies=b2, sign_algorithm=ALGO(a), pose=pose)
      nbad = nbad + count(d1 /= d2) + count(c1 /= c2) + count(n1 /= n2) + count(f1 /= f2 .or. t1 /= t2 .or. i1 /= i2 .or. b1 /= b2)
   enddo
   call sf%exit_device
   endfunction check3

   function check2(cur, pose) result(nbad)
   !< Posed device vs posed CPU for a curve, every output, both signs.
   type(curve_object),    intent(inout) :: cur  !< Curve.
   type(rigid_pose_2d_t), intent(in)    :: pose !< Pose.
   integer                              :: nbad !< Mismatches.
   integer, parameter        :: N = 500 !< Points.
   integer(I4P), parameter   :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   real(R8P)    :: lo(2), hi(2), r(2), q(2, N), d1(N), d2(N), c1(2, N), c2(2, N), n1(2, N), n2(2, N) !< Work.
   integer(I4P) :: s1(N), s2(N), x1(N), x2(N), b1(N), b2(N), i, a !< Work.

   lo = cur%get_bmin() ; hi = cur%get_bmax()
   do i = 1, N
      call random_number(r) ; q(:, i) = to_world(pose, lo - 0.2_R8P * (hi - lo) + 1.4_R8P * r * (hi - lo))
   enddo
   nbad = 0
   call cur%enter_device
   do a = 1, 2
      call cur%distance_many(points=q, distances=d1, sign_algorithm=ALGO(a), pose=pose)
      call cur%distance_many_device(points=q, distances=d2, sign_algorithm=ALGO(a), pose=pose)
      nbad = nbad + count(d1 /= d2)
      call cur%query_many(points=q, distances=d1, closest_points=c1, segments=s1, vertices=x1, normals=n1, bodies=b1, &
                          sign_algorithm=ALGO(a), pose=pose)
      call cur%query_many_device(points=q, distances=d2, closest_points=c2, segments=s2, vertices=x2, normals=n2, bodies=b2, &
                                 sign_algorithm=ALGO(a), pose=pose)
      nbad = nbad + count(d1 /= d2) + count(c1 /= c2) + count(n1 /= n2) + count(s1 /= s2 .or. x1 /= x2 .or. b1 /= b2)
   enddo
   call cur%exit_device
   endfunction check2

   pure function thickness(xx) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xx !< Chord position.
   real(R8P)             :: y  !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.126_R8P * xx - 0.3516_R8P * xx**2 + 0.2843_R8P * xx**3 - 0.1036_R8P * xx**4)
   endfunction thickness
endprogram fossil_test_pose_device
