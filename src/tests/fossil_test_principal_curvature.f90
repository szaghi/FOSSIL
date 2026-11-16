!< FOSSIL test: principal curvatures and directions at the foot point (roadmap #109 WP1 step 03).
program fossil_test_principal_curvature
!< FOSSIL test: principal curvatures and directions at the foot point (roadmap #109 WP1 step 03).
!<
!< `build_curvature` stores the per-vertex shape operator of the jet fit; `query_many(_device)(principal_curvatures=,
!< principal_directions=)` interpolate it at the foot and eigen-decompose its tangent restriction.
!< - Icospheres of radius 2 (2, 3, 4 subdivisions), feet of points at radius 3: k1, k2 -> 1/2 at second order (the
!<   maximum relative error falls by more than 3 per level); d1, d2 unit, orthogonal, tangent to the foot normal,
!<   (d1, d2, n) right-handed (1e-12); the normal curvature (`normal_curvature`) of any tangent is 1/R within the same
!<   error. Per vertex, k1 + k2 = 2 H and k1 k2 = K of `jet_curvature` to round-off.
!< - Closed cylinders of radius 1/2 (16, 32, 64 sectors, feet away from the caps): k1 -> 2, k2 -> 0 at second order,
!<   d2 along the axis (|d2 . e_z| -> 1), d1 around it.
!< - An ellipsoid (the level-4 icosphere scaled by (2, 1.5, 1)): at the feet near the tip (2, 0, 0) the analytic
!<   curvatures a / c^2 = 2 (d1 along z) and a / b^2 = 8/9 (d2 along y).
!< - Device path as host twin (tests-gnu; the GPU in tests-nvidia): query_many_device = query_many bit for bit, posed
!<   too; the posed query = the body-frame query with the directions rotated, bit for bit; `move` keeps the shape
!<   operator with rotated directions (equal to the posed query to round-off).
!< - Refusals: not built (host and device), a wrong output size, a surface re-analyzed after the build.
!< - `eig_sym2_oac`: an exact umbilic gives the first basis vector; a diagonal matrix its sorted eigenpairs.

use fossil, only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle, to_world, to_body, rotate, &
                   STATUS_INVALID_INPUT, STATUS_OK, normal_curvature
use fossil_shape_operator, only : eig_sym2_oac
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
integer, parameter              :: NP = 200
type(surface_stl_object)        :: s, sm
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: q(NP), qw(NP), qb(NP)
type(rigid_pose_t)              :: p
real(R8P)                       :: d(NP), pk(2, NP), pd(3, 2, NP), nm(3, NP), pk2(2, NP), pd2(3, 2, NP), nm2(3, NP)
real(R8P)                       :: dir(3, NP), r(3), mx(3), e2(3), ax(3), frame, kn, t(3), c, sn, k1, k2
real(R8P), allocatable          :: hv(:), kv(:)
integer(I4P)                    :: i, lev, st, nn
logical                         :: ok(8)

call random_seed(size=nn)
call random_seed(put=[(11 + i, i = 1, nn)])
do i = 1, NP
   call random_number(r) ; dir(:, i) = 2._R8P * r - 1._R8P ; dir(:, i) = dir(:, i) / norm2(dir(:, i))
   q(i) = vector_R8P(3._R8P * dir(1, i), 3._R8P * dir(2, i), 3._R8P * dir(3, i))
enddo
! spheres: accuracy, frame, normal curvature, per-vertex consistency with the jet's H and K
frame = 0._R8P ; kn = 0._R8P
do lev = 2, 4
   call icosphere(lev, 2._R8P, f) ; call s%adopt_facets(facets=f) ; call s%build_curvature(status=st)
   call s%query_many(points=q, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
   mx(lev - 1) = maxval(abs(pk - 0.5_R8P)) / 0.5_R8P
   if (lev == 4) then
      frame = frame_error(nm, pd)
      do i = 1, NP
         call random_number(t) ; t = t - 0.5_R8P ; t = t - dot_product(t, nm(:, i)) * nm(:, i)
         kn = max(kn, abs(normal_curvature(pk(:, i), pd(:, :, i), t) - 0.5_R8P) / 0.5_R8P)
      enddo
   endif
enddo
call s%jet_curvature(H=hv, K=kv)
ok(1) = st == STATUS_OK .and. mx(1) / mx(2) > 3._R8P .and. mx(2) / mx(3) > 3._R8P .and. mx(3) < 0.02_R8P .and. &
        frame < 1.e-12_R8P .and. kn <= mx(3) * (1._R8P + 1.e-9_R8P) .and. &
        maxval(abs(s%curv_k(1, :) + s%curv_k(2, :) - 2._R8P * hv)) < 1.e-12_R8P .and. &
        maxval(abs(s%curv_k(1, :) * s%curv_k(2, :) - kv)) < 1.e-12_R8P
print '(A,3F8.4,A,ES9.1,A,F8.4,A,ES9.1,A,L1)', 'sphere R=2, max relative |k - 1/2| ', mx, '  frame ', frame, &
      '  normal curvature ', kn, '  vertex |k1 + k2 - 2H| ', maxval(abs(s%curv_k(1, :) + s%curv_k(2, :) - 2._R8P * hv)), &
      '  ok ', ok(1)
! cylinders: k1 -> 1/R around, k2 -> 0 along the axis
do lev = 1, 3
   call cylinder(0.5_R8P, 8 * 2**lev, 4 * 2**lev, f) ; call s%adopt_facets(facets=f) ; call s%build_curvature
   do i = 1, NP
      call random_number(r)
      q(i) = vector_R8P(0.9_R8P * cos(2 * PI * r(1)), 0.9_R8P * sin(2 * PI * r(1)), 0.3_R8P + 0.4_R8P * r(2))
   enddo
   call s%query_many(points=q, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
   e2(lev) = max(maxval(abs(pk(1, :) - 2._R8P)), maxval(abs(pk(2, :))))
   ax(lev) = maxval(1._R8P - abs(pd(3, 2, :)))
enddo
kn = 0._R8P
do i = 1, NP
   kn = max(kn, abs(normal_curvature(pk(:, i), pd(:, :, i), [0._R8P, 0._R8P, 1._R8P])), &
            abs(normal_curvature(pk(:, i), pd(:, :, i), [-nm(2, i), nm(1, i), 0._R8P]) - 2._R8P))
enddo
ok(2) = e2(1) / e2(2) > 3._R8P .and. e2(2) / e2(3) > 3._R8P .and. e2(3) < 0.05_R8P .and. ax(3) < 1.e-6_R8P .and. &
        kn < 0.05_R8P .and. frame_error(nm, pd) < 1.e-12_R8P
print '(A,3ES10.2,A,3ES10.2,A,ES9.2,A,L1)', 'cylinder R=1/2, max err (k1 - 2, k2) ', e2, '  1 - |d2.e_z| ', ax, &
      '  normal curvature err ', kn, '  ok ', ok(2)
! ellipsoid near the tip (2, 0, 0)
call icosphere(4_I4P, 1._R8P, f)
do i = 1, size(f)
   do nn = 1, 3
      f(i)%vertex(nn)%x = 2._R8P * f(i)%vertex(nn)%x ; f(i)%vertex(nn)%y = 1.5_R8P * f(i)%vertex(nn)%y
   enddo
   call f(i)%compute_metrix
enddo
call s%adopt_facets(facets=f) ; call s%build_curvature
do i = 1, NP
   call random_number(r) ; q(i) = vector_R8P(2.5_R8P, 0.04_R8P * (r(1) - 0.5_R8P), 0.04_R8P * (r(2) - 0.5_R8P))
enddo
call s%query_many(points=q, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
mx(1) = maxval(abs(pk(1, :) - 2._R8P)) / 2._R8P ; mx(2) = maxval(abs(pk(2, :) - 8._R8P / 9._R8P)) / (8._R8P / 9._R8P)
mx(3) = max(maxval(1._R8P - abs(pd(3, 1, :))), maxval(1._R8P - abs(pd(2, 2, :))))
ok(3) = mx(1) < 0.05_R8P .and. mx(2) < 0.05_R8P .and. mx(3) < 0.01_R8P
print '(A,2F8.4,A,ES9.2,A,L1)', 'ellipsoid (2, 1.5, 1) at the tip, relative err (k1, k2) ', mx(1:2), &
      '  direction err ', mx(3), '  ok ', ok(3)
! device path (as host twin under gfortran), posed, move
call s%enter_device
do i = 1, NP
   q(i) = vector_R8P(2.6_R8P * dir(1, i), 2.1_R8P * dir(2, i), 1.6_R8P * dir(3, i))
enddo
call s%query_many(points=q, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
call s%query_many_device(points=q, distances=d, normals=nm2, principal_curvatures=pk2, principal_directions=pd2, status=st)
ok(4) = st == STATUS_OK .and. all(pk2 == pk) .and. all(pd2 == pd) .and. all(nm2 == nm)
p = pose_from_axis_angle([1._R8P, 2._R8P, 0.5_R8P], 0.7_R8P, [0.3_R8P, -1._R8P, 2._R8P])
do i = 1, NP
   r = to_world(p, [q(i)%x, q(i)%y, q(i)%z]) ; qw(i) = vector_R8P(r(1), r(2), r(3))
   r = to_body(p, r) ; qb(i) = vector_R8P(r(1), r(2), r(3)) ! the body-frame points the posed query uses
enddo
call s%query_many(points=qb, distances=d, principal_curvatures=pk, principal_directions=pd)
call s%query_many(points=qw, distances=d, principal_curvatures=pk2, principal_directions=pd2, pose=p)
do i = 1, NP
   pd(:, 1, i) = rotate(p, pd(:, 1, i)) ; pd(:, 2, i) = rotate(p, pd(:, 2, i))
enddo
ok(4) = ok(4) .and. all(pk2 == pk) .and. all(pd2 == pd)
call s%query_many_device(points=qw, distances=d, principal_curvatures=pk, principal_directions=pd, pose=p)
ok(4) = ok(4) .and. all(pk2 == pk) .and. all(pd2 == pd)
print '(A,L1)', 'device = host bit for bit, posed = rotated body-frame query bit for bit: ', ok(4)
call s%exit_device
sm = s ; call sm%move(p)
call sm%query_many(points=qw, distances=d, principal_curvatures=pk, principal_directions=pd)
ok(5) = maxval(abs(pk - pk2)) < 1.e-9_R8P .and. maxval(abs(abs(sum(pd(:, 1, :) * pd2(:, 1, :), 1)) - 1._R8P)) < 1.e-6_R8P
print '(A,ES9.2,A,L1)', 'move keeps the shape operator, max |k - k(posed)| ', maxval(abs(pk - pk2)), '  ok ', ok(5)
! refusals
call s%query_many(points=q, distances=d, principal_curvatures=pk(:, 1:NP - 1), status=st)
ok(6) = st == STATUS_INVALID_INPUT
call s%query_many(points=q, distances=d, principal_directions=pd(:, :, 1:NP - 1), status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call s%query_many_device(points=q, distances=d, principal_curvatures=pk, status=st) ! not on the device
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call s%analyze ! drops the shape operator
call s%query_many(points=q, distances=d, principal_curvatures=pk, status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call s%enter_device
call s%query_many_device(points=q, distances=d, principal_directions=pd, status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call s%query_many(points=q, distances=d, status=st)
ok(6) = ok(6) .and. st == STATUS_OK
call s%exit_device
print '(A,L1)', 'not built, wrong size, re-analyzed, not mapped: refused: ', ok(6)
! closed-form 2 x 2 eigen-decomposition
call eig_sym2_oac(0.7_R8P, 0._R8P, 0.7_R8P, k1, k2, c, sn)
ok(7) = k1 == 0.7_R8P .and. k2 == 0.7_R8P .and. c == 1._R8P .and. sn == 0._R8P
call eig_sym2_oac(-1._R8P, 0._R8P, 3._R8P, k1, k2, c, sn)
ok(7) = ok(7) .and. k1 == 3._R8P .and. k2 == -1._R8P .and. c == 0._R8P .and. sn == 1._R8P
call eig_sym2_oac(2._R8P, 1._R8P, 2._R8P, k1, k2, c, sn)
ok(7) = ok(7) .and. abs(k1 - 3._R8P) < 1.e-15_R8P .and. abs(k2 - 1._R8P) < 1.e-15_R8P .and. abs(c - sn) < 1.e-15_R8P
print '(A,L1)', 'eig_sym2_oac: umbilic -> (1, 0), sorted eigenpairs: ', ok(7)
! a flat face: k = 0 exactly
call box(8_I4P, f) ; call s%adopt_facets(facets=f) ; call s%build_curvature
do i = 1, NP
   call random_number(r) ; q(i) = vector_R8P(0.3_R8P + 0.4_R8P * r(1), 0.3_R8P + 0.4_R8P * r(2), 1.2_R8P)
enddo
call s%query_many(points=q, distances=d, normals=nm, principal_curvatures=pk, principal_directions=pd)
ok(8) = maxval(abs(pk)) == 0._R8P .and. frame_error(nm, pd) < 1.e-12_R8P
print '(A,ES9.1,A,L1)', 'flat face max |k| ', maxval(abs(pk)), '  ok ', ok(8)
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   pure function frame_error(n, dd) result(e)
   !< Largest defect of (d1, d2, n) as a right-handed orthonormal frame.
   real(R8P), intent(in) :: n(:,:), dd(:,:,:) !< Normals [3, m], directions [3, 2, m].
   real(R8P)             :: e                 !< Defect.
   real(R8P)             :: x(3)              !< d1 x d2.
   integer               :: k                 !< Counter.

   e = 0._R8P
   do k = 1, size(n, 2)
      x = [dd(2, 1, k) * dd(3, 2, k) - dd(3, 1, k) * dd(2, 2, k), dd(3, 1, k) * dd(1, 2, k) - dd(1, 1, k) * dd(3, 2, k), &
           dd(1, 1, k) * dd(2, 2, k) - dd(2, 1, k) * dd(1, 2, k)]
      e = max(e, abs(norm2(dd(:, 1, k)) - 1._R8P), abs(norm2(dd(:, 2, k)) - 1._R8P), abs(dot_product(dd(:, 1, k), dd(:, 2, k))), &
              abs(dot_product(dd(:, 1, k), n(:, k))), abs(dot_product(dd(:, 2, k), n(:, k))), maxval(abs(x - n(:, k))))
   enddo
   endfunction frame_error

   subroutine put(fo, k, a, b, cc)
   !< Facet k from three points.
   type(facet_object), intent(inout) :: fo(:)     !< Facets.
   integer(I4P),       intent(in)    :: k         !< Index.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo(k)%vertex(1) = vector_R8P(a(1), a(2), a(3))
   fo(k)%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo(k)%vertex(3) = vector_R8P(cc(1), cc(2), cc(3))
   call fo(k)%compute_metrix
   endsubroutine put

   subroutine icosphere(levels, radius, fo)
   !< Icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere.
   integer(I4P),                    intent(in)  :: levels   !< Subdivisions.
   real(R8P),                       intent(in)  :: radius   !< Radius.
   type(facet_object), allocatable, intent(out) :: fo(:)    !< Facets.
   real(R8P), allocatable  :: tri(:,:,:), nxt(:,:,:)
   real(R8P)               :: vv(3, 12), g
   integer(I4P), parameter :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                    2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                    4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                    5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   integer(I4P)            :: k, l, m

   g = (1._R8P + sqrt(5._R8P)) / 2._R8P
   vv = reshape([-1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, 0._R8P, &
                 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, &
                 g, 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, -g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P], [3, 12])
   do k = 1, 12
      vv(:, k) = vv(:, k) / norm2(vv(:, k))
   enddo
   allocate(tri(3, 3, 20))
   do k = 1, 20
      tri(:, :, k) = vv(:, IDX(:, k))
   enddo
   do l = 1, levels
      m = size(tri, dim=3)
      allocate(nxt(3, 3, 4 * m))
      do k = 1, m
         nxt(:, :, 4 * k - 3) = reshape([tri(:, 1, k), mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
         nxt(:, :, 4 * k - 2) = reshape([mid(tri(:, 1, k), tri(:, 2, k)), tri(:, 2, k), mid(tri(:, 2, k), tri(:, 3, k))], [3, 3])
         nxt(:, :, 4 * k - 1) = reshape([mid(tri(:, 3, k), tri(:, 1, k)), mid(tri(:, 2, k), tri(:, 3, k)), tri(:, 3, k)], [3, 3])
         nxt(:, :, 4 * k)     = reshape([mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 2, k), tri(:, 3, k)), &
                                         mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   m = size(tri, dim=3)
   allocate(fo(m))
   do k = 1, m
      call put(fo, k, radius * tri(:, 1, k), radius * tri(:, 2, k), radius * tri(:, 3, k))
   enddo
   endsubroutine icosphere

   pure function mid(a, b) result(m)
   !< Unit midpoint.
   real(R8P), intent(in) :: a(3), b(3) !< Points.
   real(R8P)             :: m(3)       !< Midpoint on the unit sphere.

   m = (a + b) / norm2(a + b)
   endfunction mid

   subroutine cylinder(radius, nt, nz, fo)
   !< Closed cylinder of height 1 on the z axis: nt sectors, nz bands, caps as fans.
   real(R8P),                       intent(in)  :: radius !< Radius.
   integer,                         intent(in)  :: nt, nz !< Sectors, bands.
   type(facet_object), allocatable, intent(out) :: fo(:)  !< Facets.
   real(R8P)    :: a(3), b(3), cc(3), dd(3), t0, t1, z0, z1
   integer(I4P) :: k, j, m

   allocate(fo(2 * nt * nz + 2 * nt))
   m = 0
   do k = 0, nt - 1
      t0 = 2 * PI * k / nt ; t1 = 2 * PI * mod(k + 1, nt) / nt ! mod: the seam vertices weld
      do j = 0, nz - 1
         z0 = real(j, R8P) / nz ; z1 = real(j + 1, R8P) / nz
         a = [radius * cos(t0), radius * sin(t0), z0] ; b = [radius * cos(t1), radius * sin(t1), z0]
         cc = [radius * cos(t1), radius * sin(t1), z1] ; dd = [radius * cos(t0), radius * sin(t0), z1]
         m = m + 1 ; call put(fo, m, a, b, cc)
         m = m + 1 ; call put(fo, m, a, cc, dd)
      enddo
      a = [radius * cos(t0), radius * sin(t0), 0._R8P] ; b = [radius * cos(t1), radius * sin(t1), 0._R8P]
      m = m + 1 ; call put(fo, m, [0._R8P, 0._R8P, 0._R8P], b, a)
      a(3) = 1._R8P ; b(3) = 1._R8P
      m = m + 1 ; call put(fo, m, [0._R8P, 0._R8P, 1._R8P], a, b)
   enddo
   endsubroutine cylinder

   subroutine box(n, fo)
   !< Unit cube, every face an n x n grid of squares split in two triangles (outward normals).
   integer(I4P),                    intent(in)  :: n  !< Squares per side.
   type(facet_object), allocatable, intent(out) :: fo(:) !< Facets.
   real(R8P)    :: a(3), b(3), cc(3), dd(3), s0, s1, t0, t1
   integer(I4P) :: ax_, side, i_, j, m, ia, ib

   allocate(fo(12 * n * n))
   m = 0
   do ax_ = 1, 3
      ia = mod(ax_, 3) + 1 ; ib = mod(ax_ + 1, 3) + 1
      do side = 0, 1
         do j = 0, n - 1
            do i_ = 0, n - 1
               s0 = real(i_, R8P) / n ; s1 = real(i_ + 1, R8P) / n ; t0 = real(j, R8P) / n ; t1 = real(j + 1, R8P) / n
               a(ax_) = side ; b(ax_) = side ; cc(ax_) = side ; dd(ax_) = side
               a(ia) = s0 ; a(ib) = t0 ; b(ia) = s1 ; b(ib) = t0 ; cc(ia) = s1 ; cc(ib) = t1 ; dd(ia) = s0 ; dd(ib) = t1
               if (side == 1) then
                  m = m + 1 ; call put(fo, m, a, b, cc)
                  m = m + 1 ; call put(fo, m, a, cc, dd)
               else
                  m = m + 1 ; call put(fo, m, a, cc, b)
                  m = m + 1 ; call put(fo, m, a, dd, cc)
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine box
endprogram fossil_test_principal_curvature
