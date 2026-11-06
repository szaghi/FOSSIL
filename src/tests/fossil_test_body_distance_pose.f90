!< FOSSIL test: body-body distance with poses (issue #67, step 2).
!<
!< A 10^3 sphere and a 10^3 torus, each at 8 fixed-seed rigid poses: distance_to(pose_self=, pose_other=) on the bodies
!< as they are (trees not refitted) equals distance_to on copies moved by `move`, bit for bit (distance, facets,
!< witness points, contact), apart, touching or crossing.
program fossil_test_body_distance_pose

use fossil,                only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use fossil_boolean_bodies, only : BODY_SPHERE, BODY_TORUS, mc_body
use penf,                  only : I4P, I8P, R8P

implicit none

type(surface_stl_object) :: a, b, am, bm
type(rigid_pose_t)       :: pa_, pb_
real(R8P)                :: d1, d2, p1(3), q1(3), p2(3), q2(3)
integer(I4P)             :: f1, g1, f2, g2, k, nbad, nhit
logical                  :: h1, h2
integer(I8P)             :: seed

call mc_body(BODY_SPHERE, 10_I4P, a) ; call mc_body(BODY_TORUS, 10_I4P, b)
seed = 20261107_I8P ; nbad = 0 ; nhit = 0
do k = 1, 8
   pa_ = pose_from_axis_angle([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P], 6.28_R8P * rnd(), &
                              [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   pb_ = pose_from_axis_angle([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P], 6.28_R8P * rnd(), &
                              2._R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   call a%distance_to(b, d1, point_self=p1, point_other=q1, facet_self=f1, facet_other=g1, intersecting=h1, &
                      pose_self=pa_, pose_other=pb_)
   am = a ; bm = b
   call am%move(pa_) ; call bm%move(pb_)
   call am%distance_to(bm, d2, point_self=p2, point_other=q2, facet_self=f2, facet_other=g2, intersecting=h2)
   if (d1 /= d2 .or. f1 /= f2 .or. g1 /= g2 .or. any(p1 /= p2) .or. any(q1 /= q2) .or. (h1 .neqv. h2)) nbad = nbad + 1
   if (h1) nhit = nhit + 1
enddo
print '(A,I0,A,I0)', 'pose vs move mismatches: ', nbad, ', intersecting cases ', nhit
print '(A,L1)', 'Are all tests passed? ', nbad == 0 .and. nhit > 0 .and. nhit < 8

contains
   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_body_distance_pose
