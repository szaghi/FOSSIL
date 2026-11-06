!< FOSSIL test: distance between two curves (issue #67, step 4).
!<
!< 1. Two unit squares apart by 1 along u, and by d on both axes: distance 1 and sqrt(2) d, not intersecting.
!< 2. Touching at a corner, along a side, crossing: 0 and intersecting.
!< 3. A circle (64 segments) against a square at 6 fixed-seed poses: distance_to(pose_self=, pose_other=) equals the
!<    distance of the curves moved by `move`, bit for bit (distance, segments, points).
program fossil_test_curve_distance_to

use fossil, only : curve_object, rigid_pose_2d_t, pose_2d
use penf,   only : I4P, I8P, R8P

implicit none

type(curve_object)    :: a, b, am, bm
type(rigid_pose_2d_t) :: p1, p2
real(R8P)             :: d, d2, x1(2), y1(2), x2(2), y2(2), cu(64), cv(64)
integer(I4P)          :: s1, t1, s2, t2, k, nbad
logical               :: hit, ok1, ok2, ok3
integer(I8P)          :: seed

call square(a, 0._R8P, 0._R8P) ; call square(b, 2._R8P, 0._R8P)
call a%distance_to(b, d, intersecting=hit)
ok1 = abs(d - 1._R8P) < 1e-15_R8P .and. .not. hit
call square(b, 1.5_R8P, 1.5_R8P)
call a%distance_to(b, d)
ok1 = ok1 .and. abs(d - sqrt(2._R8P) * 0.5_R8P) < 1e-15_R8P
print '(A,L1)', '1. apart: ', ok1

ok2 = .true.
call square(b, 1._R8P, 1._R8P) ; call a%distance_to(b, d, intersecting=hit) ; ok2 = ok2 .and. d == 0._R8P .and. hit
call square(b, 1._R8P, 0.3_R8P) ; call a%distance_to(b, d, intersecting=hit) ; ok2 = ok2 .and. d == 0._R8P .and. hit
call square(b, 0.5_R8P, 0.5_R8P) ; call a%distance_to(b, d, intersecting=hit) ; ok2 = ok2 .and. d == 0._R8P .and. hit
print '(A,L1)', '2. touching and crossing: ', ok2

do k = 1, 64
   cu(k) = 0.5_R8P * cos(6.283185307179586_R8P * (k - 1) / 64._R8P)
   cv(k) = 0.5_R8P * sin(6.283185307179586_R8P * (k - 1) / 64._R8P)
enddo
call a%adopt_loops(u=cu, v=cv, loop_sizes=[64_I4P])
call square(b, 0._R8P, 0._R8P)
seed = 20261108_I8P ; nbad = 0
do k = 1, 6
   p1 = pose_2d(6.28_R8P * rnd(), [rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   p2 = pose_2d(6.28_R8P * rnd(), 2._R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   call a%distance_to(b, d, point_self=x1, point_other=y1, segment_self=s1, segment_other=t1, pose_self=p1, pose_other=p2)
   am = a ; bm = b
   call am%move(p1) ; call bm%move(p2)
   call am%distance_to(bm, d2, point_self=x2, point_other=y2, segment_self=s2, segment_other=t2)
   if (d /= d2 .or. s1 /= s2 .or. t1 /= t2 .or. any(x1 /= x2) .or. any(y1 /= y2)) nbad = nbad + 1
enddo
ok3 = nbad == 0
print '(A,L1,A,I0)', '3. poses: ', ok3, '  mismatches ', nbad
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3

contains
   subroutine square(c, u0, v0)
   type(curve_object), intent(inout) :: c
   real(R8P),          intent(in)    :: u0, v0

   call c%adopt_loops(u=[u0, u0 + 1._R8P, u0 + 1._R8P, u0], v=[v0, v0, v0 + 1._R8P, v0 + 1._R8P], loop_sizes=[4_I4P])
   endsubroutine square

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_curve_distance_to
