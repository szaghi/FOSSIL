!< FOSSIL test: curvature at the foot point (issue #42 step 3).
program fossil_test_foot_curvature
!< FOSSIL test: curvature at the foot point (issue #42 step 3).
!<
!< - Icospheres of radius 2 (2, 3, 4 subdivisions): (H, K) at the feet of points at radius 3 tend to
!<   (1/2, 1/4) in the RMS over the feet; the maximum relative error stays near 0.1-0.15 (the
!<   barycentric-area operators of #18 do not converge pointwise at the 12 valence-5 vertices). The
!<   same sphere with the facets reversed (a cavity seen from the fluid): H changes sign, K does not.
!< - A closed cylinder of radius 1/2 (feet away from the caps): H = 1/(2R) = 1 and K = 0 to round-off
!<   (each ring vertex sees a circle through three points).
!< - A box with subdivided faces: feet in the face interiors give (0, 0) to round-off.
!< - The posed query returns the curvatures of the body-frame query bit for bit; a wrong output size
!<   is refused (negative control).
!< - Regular N-gons of circumradius 1 (N = 16, 32, 64): every foot gives (pi/N)/sin(pi/N), the error
!<   against 1/R falls by 4 per halving of h; a hole gives the opposite sign; an open straight line 0.

use fossil, only : surface_stl_object, facet_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, &
                   pose_2d, to_world, to_body, STATUS_INVALID_INPUT, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
integer, parameter              :: NP = 200
type(surface_stl_object)        :: s
type(curve_object)              :: c
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: q(NP), qw(NP)
type(rigid_pose_t)              :: p
real(R8P)                       :: mx(3), d(NP), cv(2, NP), cv2(2, NP), dir(3, NP), eh(3), ek(3), r(3), e2(3), &
                                   pts(2, NP), k1(NP), k2(NP), u(64), v(64), t
integer(I4P)                    :: i, lev, st, nn
logical                         :: ok(6), exact

call random_seed(size=nn)
call random_seed(put=[(7 + i, i = 1, nn)])
do i = 1, NP
   call random_number(r) ; dir(:, i) = 2._R8P * r - 1._R8P ; dir(:, i) = dir(:, i) / norm2(dir(:, i))
   q(i) = vector_R8P(3._R8P * dir(1, i), 3._R8P * dir(2, i), 3._R8P * dir(3, i))
enddo
! spheres
do lev = 2, 4
   call icosphere(lev, 2._R8P, .false., f) ; call s%adopt_facets(facets=f)
   call s%query_many(points=q, distances=d, curvatures=cv)
   eh(lev - 1) = sqrt(sum((cv(1, :) - 0.5_R8P)**2) / NP) ; ek(lev - 1) = sqrt(sum((cv(2, :) - 0.25_R8P)**2) / NP)
   mx(lev - 1) = max(maxval(abs(cv(1, :) - 0.5_R8P)) / 0.5_R8P, maxval(abs(cv(2, :) - 0.25_R8P)) / 0.25_R8P)
enddo
ok(1) = eh(3) < eh(2) .and. eh(2) < eh(1) .and. ek(3) < ek(2) .and. ek(2) < ek(1) .and. maxval(mx) < 0.2_R8P
print '(A,3ES10.2,A,3ES10.2,A,3F6.3,A,L1)', 'sphere R=2, rms |H - 1/2| ', eh, '  rms |K - 1/4| ', ek, &
      '  max relative ', mx, '  ok ', ok(1)
call icosphere(4_I4P, 2._R8P, .true., f) ; call s%adopt_facets(facets=f)
call s%query_many(points=q, distances=d, curvatures=cv2)
ok(2) = all(cv2(1, :) < 0._R8P) .and. all(cv2(2, :) > 0._R8P) .and. maxval(abs(cv2(1, :) + cv(1, :))) < 1.e-12_R8P .and. &
        maxval(abs(cv2(2, :) - cv(2, :))) < 1.e-12_R8P
print '(A,2F9.5,A,L1)', 'reversed sphere: H = -H, K = K of the sphere; mean H, mean K ', sum(cv2(1, :)) / NP, &
      sum(cv2(2, :)) / NP, '  ok ', ok(2)
! cylinder and flat faces
do lev = 1, 3
   call cylinder(0.5_R8P, 8 * 2**lev, 4 * 2**lev, f) ; call s%adopt_facets(facets=f)
   do i = 1, NP
      call random_number(r) ; t = 2 * PI * r(1)
      q(i) = vector_R8P(0.9_R8P * cos(t), 0.9_R8P * sin(t), 0.3_R8P + 0.4_R8P * r(2))
   enddo
   call s%query_many(points=q, distances=d, curvatures=cv)
   e2(lev) = max(maxval(abs(cv(1, :) - 1._R8P)), maxval(abs(cv(2, :))))
enddo
call box(4_I4P, f) ; call s%adopt_facets(facets=f)
do i = 1, NP
   call random_number(r) ; r(1:2) = 0.3_R8P + 0.4_R8P * r(1:2)
   q(i) = vector_R8P(r(1), r(2), 1.2_R8P) ; if (i > NP / 2) q(i) = vector_R8P(-0.2_R8P, r(1), r(2))
enddo
call s%query_many(points=q, distances=d, curvatures=cv)
ok(3) = maxval(e2) < 1.e-12_R8P .and. maxval(abs(cv)) < 1.e-12_R8P
print '(A,3ES10.2,A,ES9.1,A,L1)', 'cylinder R=1/2, max err (H, K) ', e2, '  flat faces max |(H, K)| ', maxval(abs(cv)), &
      '  ok ', ok(3)
! pose and refusal
call icosphere(3_I4P, 2._R8P, .false., f) ; call s%adopt_facets(facets=f)
p = pose_from_axis_angle([1._R8P, 2._R8P, 0.5_R8P], 0.7_R8P, [0.3_R8P, -1._R8P, 2._R8P])
do i = 1, NP
   r = to_world(p, 3._R8P * dir(:, i)) ; qw(i) = vector_R8P(r(1), r(2), r(3))
   r = to_body(p, r) ; q(i) = vector_R8P(r(1), r(2), r(3)) ! the body-frame points the posed query uses
enddo
call s%query_many(points=q, distances=d, curvatures=cv)
call s%query_many(points=q, distances=d, curvatures=cv2(:, 1:NP - 1), status=st)
ok(4) = st == STATUS_INVALID_INPUT
call s%query_many(points=qw, distances=d, curvatures=cv2, pose=p, status=st)
ok(4) = ok(4) .and. st == STATUS_OK .and. all(cv2 == cv)
print '(A,L1)', 'posed = body-frame bit for bit, wrong size refused: ', ok(4)
! curves: regular N-gons
do lev = 1, 3
   nn = 8 * 2**lev
   do i = 1, nn
      u(i) = cos(2 * PI * (i - 1) / nn) ; v(i) = sin(2 * PI * (i - 1) / nn)
   enddo
   call c%adopt_loops(u=u(1:nn), v=v(1:nn), loop_sizes=[nn])
   do i = 1, NP
      pts(:, i) = 1.7_R8P * [cos(dir(1, i) * PI), sin(dir(1, i) * PI)] * sign(1._R8P, dir(2, i))
   enddo
   call c%query_many(points=pts, distances=d, curvatures=k1)
   t = (PI / nn) / sin(PI / nn)
   e2(lev) = maxval(abs(k1 - 1._R8P))
   if (lev == 1) exact = .true.
   exact = exact .and. maxval(abs(k1 - t)) < 1.e-12_R8P
enddo
ok(5) = exact .and. e2(1) / e2(2) > 3.9_R8P .and. e2(2) / e2(3) > 3.9_R8P
print '(A,3ES10.2,A,2F6.2,A,L1)', 'N-gons R=1, max |k - 1| ', e2, '  ratios ', e2(1) / e2(2), e2(2) / e2(3), '  ok ', ok(5)
! hole (an annulus, the inner loop clockwise of radius 1/2) and an open straight line
do i = 1, 32
   u(i) = 2 * cos(2 * PI * (i - 1) / 32) ; v(i) = 2 * sin(2 * PI * (i - 1) / 32)
   u(32 + i) = 0.5_R8P * cos(-2 * PI * (i - 1) / 32) ; v(32 + i) = 0.5_R8P * sin(-2 * PI * (i - 1) / 32)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[32_I4P, 32_I4P])
do i = 1, NP
   pts(:, i) = 0.3_R8P * [cos(dir(1, i) * PI), sin(dir(1, i) * PI)] * sign(1._R8P, dir(2, i))
enddo
call c%query_many(points=pts, distances=d, curvatures=k1)
call c%adopt_loops(u=[0._R8P, 1._R8P, 2._R8P, 3._R8P], v=[1._R8P, 1._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P], &
                   closed=[.false.])
call c%query_many(points=pts, distances=d, curvatures=k2)
t = 2 * (PI / 32) / sin(PI / 32)
ok(6) = maxval(abs(k1 + t)) < 1.e-12_R8P .and. maxval(abs(k2)) == 0._R8P
print '(A,F9.5,A,ES9.1,A,L1)', 'hole r=1/2: k ', sum(k1) / NP, '  open line max |k| ', maxval(abs(k2)), '  ok ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
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

   subroutine icosphere(levels, radius, reversed, fo)
   !< Icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere.
   integer(I4P),                    intent(in)  :: levels   !< Subdivisions.
   real(R8P),                       intent(in)  :: radius   !< Radius.
   logical,                         intent(in)  :: reversed !< Inward normals.
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
      if (reversed) then
         call put(fo, k, radius * tri(:, 1, k), radius * tri(:, 3, k), radius * tri(:, 2, k))
      else
         call put(fo, k, radius * tri(:, 1, k), radius * tri(:, 2, k), radius * tri(:, 3, k))
      endif
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
   integer(I4P) :: ax, side, i, j, m, ia, ib

   allocate(fo(12 * n * n))
   m = 0
   do ax = 1, 3
      ia = mod(ax, 3) + 1 ; ib = mod(ax + 1, 3) + 1
      do side = 0, 1
         do j = 0, n - 1
            do i = 0, n - 1
               s0 = real(i, R8P) / n ; s1 = real(i + 1, R8P) / n ; t0 = real(j, R8P) / n ; t1 = real(j + 1, R8P) / n
               a(ax) = side ; b(ax) = side ; cc(ax) = side ; dd(ax) = side
               a(ia) = s0 ; a(ib) = t0 ; b(ia) = s1 ; b(ib) = t0 ; cc(ia) = s1 ; cc(ib) = t1 ; dd(ia) = s0 ; dd(ib) = t1
               if (side == 1) then ! (ia, ib, ax) is right-handed: CCW in (ia, ib) points along +ax
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
endprogram fossil_test_foot_curvature
