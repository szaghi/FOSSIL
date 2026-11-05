!< FOSSIL test: rigid poses and body-frame point queries (issue #41 step 1).
program fossil_test_pose
!< FOSSIL test: rigid poses and body-frame point queries (issue #41 step 1).
!<
!< - Pose algebra: pose o inverse = identity, to_world(to_body(x)) = x, axis-angle = quaternion of the
!<   same rotation, to round-off.
!< - Identity pose = no pose, bit for bit (surface and curve, every output).
!< - Against the mesh (curve) moved explicitly and rebuilt: cube, cone, L-prism of the sign
!<   benchmark and cube-04.stl (3072 facets, clean and closed; it replaced the bunny, issue #45), under
!<   rotations about several axes with large translations, and a 256-sided faceted cone (needle facets of aspect
!<   163); a NACA 0012 and an annulus under in-plane poses: signed distances to 1e-14 (L + |t|), closest points to
!<   max(1e-14, 0.1 eps a^2) (L + |t|) with a the largest facet aspect (longest edge over height; issue #56: the foot
!<   point on a needle is conditioned by a^2, measured 0.020-0.023 eps a^2 (L + |t|) at a = 10 to 163), normals to
!<   1e-14 (L + |t|) / h (h the shortest edge: the reference normals are rebuilt from moved, rounded
!<   vertices), inside flags equal at points farther than 1e-12 L. Closest
!<   points and normals are compared where the foot point is unique (`ambiguous` with a 1e-9 relative
!<   tolerance): at a near-tie either foot is a correct answer. The points are seeded (deterministic).

use fossil,                   only : surface_stl_object, curve_object, facet_object, rigid_pose_t, rigid_pose_2d_t, &
                                     pose_identity, pose_from_axis_angle, pose_from_quaternion, pose_2d, pose_compose, &
                                     pose_inverse, to_body, to_world, SIGN_FAST_WINDING_NUMBER
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface, make_apex_body
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(body_t), allocatable :: b(:)
type(soup_t)              :: fine
type(body_t)              :: needles
type(surface_stl_object)  :: s, body
type(facet_object), pointer :: fp
type(curve_object)        :: c
type(rigid_pose_t)        :: p1, p2, p3, id
type(rigid_pose_2d_t)     :: q1
real(R8P)                 :: u(200), v(200), x(3), e
integer(I4P)              :: k, f, m
logical                   :: ok(4)

call random_seed(size=k)
call random_seed(put=[(4141 + m, m = 1, k)]) ! deterministic points
! 1. algebra
p1 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.7_R8P, [10._R8P, -3._R8P, 5._R8P])
p2 = pose_from_quaternion([cos(0.35_R8P), sin(0.35_R8P) * [1._R8P, 2._R8P, 3._R8P] / sqrt(14._R8P)], [10._R8P, -3._R8P, 5._R8P])
id = pose_compose(p1, pose_inverse(p1))
x = [0.3_R8P, -1.7_R8P, 2.2_R8P]
p3 = pose_identity()
e = max(maxval(abs(id%r - p3%r)), maxval(abs(id%t)), maxval(abs(to_world(p1, to_body(p1, x)) - x)) / 20._R8P, &
        maxval(abs(p1%r - p2%r)), maxval(abs(p1%t - p2%t)))
ok(1) = e < 4.e-15_R8P
print '(A,ES9.1,A,L1)', 'pose algebra: max err ', e, '  ok ', ok(1)

call bench_bodies(b)
! cube-04.stl as a soup: a clean closed body, 3072 facets (issue #45 replaced the bunny; dragon.stl fails: a facet
! tied at a non-manifold edge flips one pseudo-normal sign between posed and moved queries); a 256-sided faceted cone,
! needle facets of aspect 163, under the aspect-scaled closest-point bound (issue #56)
call make_apex_body(256, 15._R8P, 'cone 15 (256)', needles)
call body%load_from_file(file_name='src/tests/cube-04.stl', guess_format=.true.)
allocate(fine%v(3, 3, body%get_facets_number()))
do f = 1, body%get_facets_number()
   fp => body%facet_at(f)
   do m = 1, 3
      fine%v(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
   enddo
enddo
! 2. identity = no pose; 3. against moved and rebuilt meshes
p3 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], -2.1_R8P, [-100._R8P, 40._R8P, 0.25_R8P])
call to_surface(b(1)%soup, s)
ok(2) = same_identity(s)
ok(3) = moved3(b(1)%soup, p1, 'cube') .and. moved3(b(1)%soup, p3, 'cube, far') .and. moved3(b(7)%soup, p2, 'cone') .and. &
        moved3(b(8)%soup, p1, 'L-prism') .and. moved3(fine, p1, 'cube-04') .and. &
        moved3(needles%soup, p1, 'cone 256') .and. moved3(needles%soup, p2, 'cone 256, p2')
! 4. curves
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
q1 = pose_2d(0.9_R8P, [5._R8P, -2._R8P])
ok(4) = moved2(u, v, [200_I4P], q1, 'NACA 0012') .and. &
        moved2([0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
               [0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], [4_I4P, 4_I4P], &
               pose_2d(-0.4_R8P, [30._R8P, 7._R8P]), 'annulus')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function same_identity(sf) result(same)
   !< Identity pose gives the unposed outputs bit for bit.
   type(surface_stl_object), intent(in) :: sf   !< Surface.
   logical                              :: same !< Result.
   type(vector_R8P) :: q(50)  !< Points.
   real(R8P)        :: d1(50), d2(50), c1(3, 50), c2(3, 50), n1(3, 50), n2(3, 50), r(3) !< Work.
   integer(I4P)     :: i, f1(50), f2(50) !< Work.

   do i = 1, 50
      call random_number(r) ; q(i) = vector_R8P(-0.5_R8P + 2 * r(1), -0.5_R8P + 2 * r(2), -0.5_R8P + 2 * r(3))
   enddo
   call sf%query_many(points=q, distances=d1, closest_points=c1, facets=f1, normals=n1)
   call sf%query_many(points=q, distances=d2, closest_points=c2, facets=f2, normals=n2, pose=pose_identity())
   same = all(d1 == d2) .and. all(c1 == c2) .and. all(f1 == f2) .and. all(n1 == n2)
   call sf%distance_many(points=q, distances=d1, is_signed=.true.)
   call sf%distance_many(points=q, distances=d2, is_signed=.true., pose=pose_identity())
   same = same .and. all(d1 == d2)
   print '(A,L1)', 'identity pose = no pose, bit for bit: ', same
   endfunction same_identity

   function moved3(sp, pose, label) result(good)
   !< Posed queries against the mesh moved explicitly and rebuilt.
   type(soup_t),       intent(in) :: sp    !< Body-frame mesh.
   type(rigid_pose_t), intent(in) :: pose  !< Pose.
   character(*),       intent(in) :: label !< Label.
   logical                        :: good  !< Result.
   integer, parameter :: N = 400 !< Points.
   type(soup_t)             :: mv      !< Moved mesh.
   type(surface_stl_object) :: sb, sm  !< Body and moved surfaces.
   type(vector_R8P)         :: q(N)    !< World points.
   real(R8P) :: lo(3), hi(3), r(3), w(3), l, scale, d1(N), d2(N), c1(3, N), c2(3, N), n1(3, N), n2(3, N), ed, ec, en, hmin !< Work.
   real(R8P) :: amax, cbound !< Largest facet aspect, closest-point bound.
   integer(I4P) :: i, ff, mm, nsign, ninc !< Work.
   logical      :: in1, in2, amb(N) !< Work.

   mv = sp
   do ff = 1, size(sp%v, 3)
      do mm = 1, 3
         mv%v(:, mm, ff) = to_world(pose, sp%v(:, mm, ff))
      enddo
   enddo
   call to_surface(sp, sb) ; call to_surface(mv, sm)
   lo = [minval(mv%v(1, :, :)), minval(mv%v(2, :, :)), minval(mv%v(3, :, :))]
   hi = [maxval(mv%v(1, :, :)), maxval(mv%v(2, :, :)), maxval(mv%v(3, :, :))]
   l = norm2(hi - lo) ; scale = l + norm2(pose%t)
   do i = 1, N
      call random_number(r) ; w = lo - 0.2_R8P * (hi - lo) + 1.4_R8P * r * (hi - lo)
      q(i) = vector_R8P(w(1), w(2), w(3))
   enddo
   call sb%query_many(points=q, distances=d1, closest_points=c1, normals=n1, pose=pose)
   call sm%query_many(points=q, distances=d2, closest_points=c2, normals=n2, ambiguous=amb, rel_tolerance=1.e-9_R8P)
   ed = 0._R8P ; ec = 0._R8P ; en = 0._R8P ; nsign = 0 ; ninc = 0
   do i = 1, N
      ed = max(ed, abs(abs(d1(i)) - abs(d2(i))))
      if (.not. amb(i)) ec = max(ec, norm2(c1(:, i) - c2(:, i))) ! a near-tie may pick either foot
      if (abs(d2(i)) > 1.e-12_R8P * l) then
         if ((d1(i) < 0._R8P) .neqv. (d2(i) < 0._R8P)) nsign = nsign + 1
         if (.not. amb(i)) en = max(en, norm2(n1(:, i) - n2(:, i)))
         in1 = sb%is_point_inside(q(i), sign_algorithm=SIGN_FAST_WINDING_NUMBER, pose=pose)
         in2 = sm%is_point_inside(q(i), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
         if (in1 .neqv. in2) ninc = ninc + 1
      endif
   enddo
   call sb%distance_many(points=q, distances=d1, is_signed=.true., is_square_root=.true., pose=pose)
   ed = max(ed, maxval(abs(abs(d1) - abs(d2))))
   ! the reference normals come from moved vertices: their own error is ~eps (L + |t|) / (shortest edge)
   hmin = huge(1._R8P)
   do ff = 1, size(sp%v, 3)
      hmin = min(hmin, norm2(sp%v(:, 2, ff) - sp%v(:, 1, ff)), norm2(sp%v(:, 3, ff) - sp%v(:, 2, ff)), &
                 norm2(sp%v(:, 1, ff) - sp%v(:, 3, ff)))
   enddo
   ! the foot point on a needle is conditioned by its aspect squared (issue #56): largest aspect, longest edge over height
   amax = 0._R8P
   do ff = 1, size(sp%v, 3)
      amax = max(amax, max(norm2(sp%v(:, 2, ff) - sp%v(:, 1, ff)), norm2(sp%v(:, 3, ff) - sp%v(:, 2, ff)), &
                           norm2(sp%v(:, 1, ff) - sp%v(:, 3, ff)))**2 / norm2(cross3(sp%v(:, 2, ff) - sp%v(:, 1, ff), &
                                                                                   sp%v(:, 3, ff) - sp%v(:, 1, ff))))
   enddo
   cbound = max(1.e-14_R8P, 0.1_R8P * epsilon(1._R8P) * amax**2) * scale
   good = ed <= 1.e-14_R8P * scale .and. ec <= cbound .and. en <= 1.e-14_R8P * scale / hmin .and. nsign == 0 &
          .and. ninc == 0
   print '(A,A,ES9.1,A,ES9.1,A,ES9.1,A,F6.1,A,ES9.1,A,I0,A,I0,A,L1)', label, ': err/(L+|t|) distance ', ed / scale, &
         '  closest ', ec / scale, ' (bound ', cbound / scale, ', aspect ', amax, ')  normal*h/(L+|t|) ', en * hmin / scale, &
         '  sign mismatches ', nsign, '  inside mismatches ', ninc, '  ok ', good
   endfunction moved3

   pure function cross3(a, c) result(x)
   !< Cross product of two arrays.
   real(R8P), intent(in) :: a(3), c(3) !< Vectors.
   real(R8P)             :: x(3)       !< Result.

   x = [a(2) * c(3) - a(3) * c(2), a(3) * c(1) - a(1) * c(3), a(1) * c(2) - a(2) * c(1)]
   endfunction cross3

   function moved2(uu, vv, sizes, pose, label) result(good)
   !< Posed curve queries against the curve moved explicitly and rebuilt.
   real(R8P),             intent(in) :: uu(:), vv(:) !< Body-frame vertices.
   integer(I4P),          intent(in) :: sizes(:)     !< Loop sizes.
   type(rigid_pose_2d_t), intent(in) :: pose         !< Pose.
   character(*),          intent(in) :: label        !< Label.
   logical                           :: good         !< Result.
   integer, parameter :: N = 400 !< Points.
   type(curve_object) :: cb, cm !< Body and moved curves.
   real(R8P) :: mu(size(uu)), mvv(size(uu)), w(2), lo(2), hi(2), r(2), q(2, N), d1(N), d2(N), c1(2, N), c2(2, N), n1(2, N), &
                n2(2, N), l, scale, ed, ec, en, hmin !< Work.
   integer(I4P) :: i, nsign, ninc !< Work.
   logical      :: amb(N) !< Ambiguous foot points of the reference.

   do i = 1, size(uu)
      w = to_world(pose, [uu(i), vv(i)]) ; mu(i) = w(1) ; mvv(i) = w(2)
   enddo
   call cb%adopt_loops(u=uu, v=vv, loop_sizes=sizes) ; call cm%adopt_loops(u=mu, v=mvv, loop_sizes=sizes)
   lo = cm%get_bmin() ; hi = cm%get_bmax() ; l = norm2(hi - lo) ; scale = l + norm2(pose%t)
   do i = 1, N
      call random_number(r) ; q(:, i) = lo - 0.2_R8P * (hi - lo) + 1.4_R8P * r * (hi - lo)
   enddo
   call cb%query_many(points=q, distances=d1, closest_points=c1, normals=n1, pose=pose)
   call cm%query_many(points=q, distances=d2, closest_points=c2, normals=n2, ambiguous=amb, rel_tolerance=1.e-9_R8P)
   ed = maxval(abs(abs(d1) - abs(d2))) ; ec = 0._R8P ; en = 0._R8P ; nsign = 0 ; ninc = 0
   do i = 1, N
      if (.not. amb(i)) ec = max(ec, norm2(c1(:, i) - c2(:, i))) ! a near-tie may pick either foot
      if (abs(d2(i)) > 1.e-12_R8P * l) then
         if ((d1(i) < 0._R8P) .neqv. (d2(i) < 0._R8P)) nsign = nsign + 1
         if (.not. amb(i)) en = max(en, norm2(n1(:, i) - n2(:, i)))
         if (cb%is_point_inside(q(:, i), pose=pose) .neqv. cm%is_point_inside(q(:, i))) ninc = ninc + 1
      endif
   enddo
   call cb%distance_many(points=q, distances=d1, pose=pose)
   ed = max(ed, maxval(abs(abs(d1) - abs(d2))))
   hmin = minval([(hypot(uu(merge(i + 1, i - sizes(1) + 1, i < sizes(1))) - uu(i), &
                         vv(merge(i + 1, i - sizes(1) + 1, i < sizes(1))) - vv(i)), i = 1, sizes(1))])
   good = ed <= 1.e-14_R8P * scale .and. ec <= 1.e-14_R8P * scale .and. en <= 1.e-14_R8P * scale / hmin .and. nsign == 0 &
          .and. ninc == 0
   print '(A,A,ES9.1,A,ES9.1,A,ES9.1,A,I0,A,I0,A,L1)', label, ': err/(L+|t|) distance ', ed / scale, '  closest ', ec / scale, &
         '  normal*h/(L+|t|) ', en * hmin / scale, '  sign mismatches ', nsign, '  inside mismatches ', ninc, '  ok ', good
   endfunction moved2

   pure function thickness(xx) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xx !< Chord position.
   real(R8P)             :: y  !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.126_R8P * xx - 0.3516_R8P * xx**2 + 0.2843_R8P * xx**3 - 0.1036_R8P * xx**4)
   endfunction thickness
endprogram fossil_test_pose
