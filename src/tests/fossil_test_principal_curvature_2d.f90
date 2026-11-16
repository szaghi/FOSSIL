!< FOSSIL test: curvature at the foot point of a curve for a wall closure (roadmap #109 WP1 step 03).
program fossil_test_principal_curvature_2d
!< FOSSIL test: curvature at the foot point of a curve for a wall closure (roadmap #109 WP1 step 03).
!<
!< `curve%build_curvature` stores a signed curvature per vertex (cubic jet over five vertices, O(h^2) on any spacing);
!< `query_many(_device)(principal_curvatures=, principal_directions=)` interpolate it along the closest segment and return
!< the unit tangent.
!< - Ellipses (2, 1) with alternating spacing (angles 2 pi (k + (-1)^k / 4) / N, N = 64, 128, 256): the maximum vertex
!<   error falls by more than 3.5 per halving (second order), the circle through three vertices (`menger_curvature`) by
!<   less than 2.5 (its O(h2 - h1) term). A uniform circle: exact at the vertices and the feet (round-off).
!< - A clockwise hole of radius 1/2: -2; an ellipse (2, 1) with vertices on its axes: a / b^2 = 2 and b / a^2 = 1/4 at
!<   the tips, second order; an open arc of a circle: exact, the one-sided circles at its ends included.
!< - The tangent is unit, normal to the pseudo-normal, along the segment (the solid on its left).
!< - Device path as host twin: query_many_device = query_many bit for bit, posed too; posed = the body-frame query with
!<   the tangent rotated, bit for bit; `move` keeps the curvatures.
!< - Refusals: not built, not mapped, a wrong size; `adopt_loops` drops the curvatures.

use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, to_world, to_body, rotate, STATUS_INVALID_INPUT, STATUS_OK, &
                   normal_curvature
use fossil_curvature, only : menger_curvature
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter   :: PI = acos(-1._R8P)
integer, parameter     :: NP = 300
type(curve_object)     :: c, cm
type(rigid_pose_2d_t)  :: p
real(R8P), allocatable :: u(:), v(:), kx(:)
real(R8P)              :: pts(2, NP), pw(2, NP), pb(2, NP), d(NP), nm(2, NP), pk(1, NP), pd(2, 1, NP), pk2(1, NP), &
                          pd2(2, 1, NP), kt(NP), ej(3), em(3), et(3), ef(3), r(2), th, e
integer(I4P)           :: i, lev, n, st, nn
logical                :: ok(6)

call random_seed(size=nn)
call random_seed(put=[(5 + i, i = 1, nn)])
do i = 1, NP
   call random_number(r) ; th = 2 * PI * r(1)
   pts(:, i) = (1.3_R8P + 0.4_R8P * r(2)) * [cos(th), sin(th)]
enddo
! ellipses (2, 1) with alternating spacing: second order where the circle through three vertices is first order
do lev = 1, 3
   n = 32 * 2**lev
   if (allocated(u)) deallocate(u, v, kx)
   allocate(u(n), v(n), kx(n))
   do i = 1, n
      th = 2 * PI * ((i - 1) + 0.25_R8P * (-1)**i) / n
      u(i) = 2 * cos(th) ; v(i) = sin(th) ; kx(i) = 2._R8P / (4._R8P * sin(th)**2 + cos(th)**2)**1.5_R8P
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[n]) ; call c%build_curvature(status=st)
   ej(lev) = maxval(abs(c%get_vertex_curvature([(i, i = 1, n)]) - kx))
   em(lev) = 0._R8P
   do i = 1, n
      em(lev) = max(em(lev), abs(menger_curvature(u(modulo(i - 2, n) + 1), v(modulo(i - 2, n) + 1), u(i), v(i), &
                                                  u(modulo(i, n) + 1), v(modulo(i, n) + 1)) - kx(i)))
   enddo
enddo
! a uniform circle: exact at the vertices, the feet within round-off
n = 64
deallocate(u, v) ; allocate(u(n), v(n))
do i = 1, n
   u(i) = cos(2 * PI * (i - 1) / n) ; v(i) = sin(2 * PI * (i - 1) / n)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[n]) ; call c%build_curvature
call c%query_many(points=pts, distances=d, normals=nm, curvatures=kt, principal_curvatures=pk, principal_directions=pd)
et(1) = maxval(abs(c%get_vertex_curvature([(i, i = 1, n)]) - 1._R8P)) ; ef(1) = maxval(abs(pk(1, :) - 1._R8P))
et(2) = maxval(abs(kt - 1._R8P))
e = 0._R8P
do i = 1, NP
   e = max(e, abs(norm2(pd(:, 1, i)) - 1._R8P), abs(dot_product(pd(:, 1, i), nm(:, i))), &
           abs(pd(1, 1, i) * nm(2, i) - pd(2, 1, i) * nm(1, i) + 1._R8P), &
           abs(normal_curvature(pk(:, i), pd(:, :, i), 3._R8P * pd(:, 1, i)) - pk(1, i)))
enddo
ok(1) = st == STATUS_OK .and. ej(1) / ej(2) > 3.5_R8P .and. ej(2) / ej(3) > 3.5_R8P .and. em(2) / em(3) < 2.5_R8P .and. &
        et(1) < 1.e-12_R8P .and. ef(1) < 1.e-12_R8P .and. e < 1.e-14_R8P
print '(A,3ES10.2,A,3ES10.2,A,ES9.1,A,ES9.1,A,ES9.1,A,ES9.1,A,L1)', &
      'ellipses (2, 1), alternating spacing, max err: vertex ', ej, '  circle through 3 ', em, &
      '  | uniform circle: vertices ', et(1), ' feet ', ef(1), ' (turning angle ', et(2), ')  tangent ', e, '  ok ', ok(1)
! hole, ellipse tips, open arc
deallocate(u, v) ; allocate(u(96), v(96))
do i = 1, 64
   u(i) = 2 * cos(2 * PI * (i - 1) / 64) ; v(i) = 2 * sin(2 * PI * (i - 1) / 64)
enddo
do i = 1, 32
   u(64 + i) = 0.5_R8P * cos(-2 * PI * (i - 1) / 32) ; v(64 + i) = 0.5_R8P * sin(-2 * PI * (i - 1) / 32)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 32_I4P]) ; call c%build_curvature
e = maxval(abs(c%get_vertex_curvature([(i, i = 65, 96)]) + 2._R8P))
ok(2) = e < 1.e-12_R8P
do lev = 1, 3
   n = 16 * 2**lev
   deallocate(u, v) ; allocate(u(n), v(n))
   do i = 1, n
      u(i) = 2 * cos(2 * PI * (i - 1) / n) ; v(i) = sin(2 * PI * (i - 1) / n)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[n]) ; call c%build_curvature
   call c%query_many(points=reshape([2.5_R8P, 0._R8P, 0._R8P, 1.5_R8P], [2, 2]), distances=d(1:2), &
                     principal_curvatures=pk(:, 1:2))
   em(lev) = max(abs(pk(1, 1) - 2._R8P) / 2._R8P, abs(pk(1, 2) - 0.25_R8P) / 0.25_R8P)
enddo
ok(2) = ok(2) .and. em(1) / em(2) > 3.5_R8P .and. em(2) / em(3) > 3.5_R8P .and. em(3) < 3.e-3_R8P
deallocate(u, v) ; allocate(u(20), v(20))
do i = 1, 20
   th = PI * (i - 1) / 19 ; u(i) = cos(th) ; v(i) = sin(th)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[20_I4P], closed=[.false.]) ; call c%build_curvature
ok(2) = ok(2) .and. maxval(abs(c%get_vertex_curvature([(i, i = 1, 20)]) - 1._R8P)) < 1.e-12_R8P
print '(A,ES9.2,A,3ES10.2,A,ES9.2,A,L1)', 'hole max |k + 2| ', e, '  ellipse tips relative err ', em, &
      '  open arc max |k - 1| ', maxval(abs(c%get_vertex_curvature([(i, i = 1, 20)]) - 1._R8P)), '  ok ', ok(2)
! device twin, posed, move
deallocate(u, v) ; allocate(u(80), v(80))
do i = 1, 80
   th = 2 * PI * (i - 1) / 80 ; u(i) = 1.2_R8P * cos(th) + 0.1_R8P * cos(3 * th) ; v(i) = 0.8_R8P * sin(th)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[80_I4P]) ; call c%build_curvature ; call c%enter_device
call c%query_many(points=pts, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
call c%query_many_device(points=pts, distances=d, principal_curvatures=pk2, principal_directions=pd2, status=st)
ok(3) = st == STATUS_OK .and. all(pk2 == pk) .and. all(pd2 == pd)
p = pose_2d(0.9_R8P, [5._R8P, -2._R8P])
do i = 1, NP
   pw(:, i) = to_world(p, pts(:, i)) ; pb(:, i) = to_body(p, pw(:, i))
enddo
call c%query_many(points=pb, distances=d, principal_curvatures=pk, principal_directions=pd)
call c%query_many(points=pw, distances=d, principal_curvatures=pk2, principal_directions=pd2, pose=p)
do i = 1, NP
   pd(:, 1, i) = rotate(p, pd(:, 1, i))
enddo
ok(3) = ok(3) .and. all(pk2 == pk) .and. all(pd2 == pd)
call c%query_many_device(points=pw, distances=d, principal_curvatures=pk, principal_directions=pd, pose=p)
ok(3) = ok(3) .and. all(pk2 == pk) .and. all(pd2 == pd)
print '(A,L1)', 'device = host bit for bit, posed = rotated body-frame query bit for bit: ', ok(3)
call c%exit_device
cm = c ; call cm%move(p)
call cm%query_many(points=pw, distances=d, principal_curvatures=pk, principal_directions=pd)
ok(4) = maxval(abs(pk - pk2)) < 1.e-12_R8P .and. maxval(abs(pd - pd2)) < 1.e-12_R8P
print '(A,ES9.2,A,L1)', 'move keeps the curvatures, max |k - k(posed)| ', maxval(abs(pk - pk2)), '  ok ', ok(4)
! refusals
call c%query_many(points=pts, distances=d, principal_curvatures=pk(:, 1:NP - 1), status=st)
ok(5) = st == STATUS_INVALID_INPUT
call c%query_many(points=pts, distances=d, principal_directions=pd(:, :, 1:NP - 1), status=st)
ok(5) = ok(5) .and. st == STATUS_INVALID_INPUT
call c%query_many_device(points=pts, distances=d, principal_curvatures=pk, status=st) ! not on the device
ok(5) = ok(5) .and. st == STATUS_INVALID_INPUT
call c%adopt_loops(u=u, v=v, loop_sizes=[80_I4P]) ! drops the curvatures
call c%query_many(points=pts, distances=d, principal_curvatures=pk, status=st)
ok(5) = ok(5) .and. st == STATUS_INVALID_INPUT .and. c%get_vertex_curvature(1_I4P) == 0._R8P
call c%enter_device
call c%query_many_device(points=pts, distances=d, principal_directions=pd, status=st)
ok(5) = ok(5) .and. st == STATUS_INVALID_INPUT
call c%exit_device
print '(A,L1)', 'not built, wrong size, re-adopted, not mapped: refused: ', ok(5)
! a straight open line: 0
call c%adopt_loops(u=[0._R8P, 1._R8P, 2._R8P, 3._R8P, 4._R8P], v=[1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
call c%build_curvature(status=st)
call c%query_many(points=pts, distances=d, principal_curvatures=pk)
ok(6) = st == STATUS_OK .and. maxval(abs(pk)) == 0._R8P
call c%adopt_loops(u=u(1:0), v=v(1:0), loop_sizes=[integer(I4P) ::])
call c%build_curvature(status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'straight open line: 0; empty curve refused: ', ok(6)
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_principal_curvature_2d
