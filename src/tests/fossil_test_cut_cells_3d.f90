!< FOSSIL test: 3-D cut cells, volumes by the divergence theorem (issue #40 step 5).
program fossil_test_cut_cells_3d
!< FOSSIL test: 3-D cut cells, volumes by the divergence theorem (issue #40 step 5).
!<
!< - An axis-aligned box inside a 6^3 block: every fluid fraction and fluid centroid exact (to
!<   round-off) against the box overlaps.
!< - A half-space (a large tetrahedron whose slanted face crosses the block, the others far
!<   away): every fluid fraction against the closed-form volume of a box cut by a plane.
!< - A rotated cube, the 10-degree pyramid and the 32-face cone on uniform and stretched blocks:
!<   the solid volumes add up to the body volume, the solid moments to volume times centroid (cube),
!<   the wetted areas to the surface area, the wall normals to zero; the per-cell divergence
!<   residual (fluid face areas times outward normals minus the wall normal integral) at round-off;
!<   all cells closed.
!< - The unit cube on grid planes: every cell 0 or 1 to round-off (the divergence sums round), 64 solid
!<   cells, and the sums above (volume 1, area 6, centroid, wall normals, divergence residual).

use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface, make_prism
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(body_t)              :: box
type(soup_t)              :: rc, tet
type(surface_stl_object)  :: s
type(cut_cells_3d_t)      :: cut
real(R8P)                 :: u(7), w(7), g(9), q(13), rot(3, 3), cth, sth, cph, sph, n(3), d, p0(3), t1(3), t2(3)
integer(I4P)              :: k, f, m
logical                   :: ok(6)
integer(I4P), parameter   :: OPEN_UNTRUSTED = 233 !< Untrusted cells of the open box (measured, unchanged by #103).

call bench_bodies(b)
do k = 1, 7
   u(k) = -0.1_R8P + 1.2_R8P * (k - 1) / 6._R8P
   w(k) = -0.1_R8P + 1.2_R8P * ((k - 1) / 6._R8P)**1.3_R8P
enddo
! 1. axis-aligned box [0.3, 0.7] x [0.2, 0.65] x [0.35, 0.8]
call make_prism(reshape([0.3_R8P, 0.2_R8P, 0.7_R8P, 0.2_R8P, 0.7_R8P, 0.65_R8P, 0.3_R8P, 0.65_R8P], [2, 4]), 0.45_R8P, 1, &
                'box', box)
box%soup%v(3, :, :) = box%soup%v(3, :, :) + 0.35_R8P
call to_surface(box%soup, s)
ok(1) = aabox(s, u, w, [0.3_R8P, 0.2_R8P, 0.35_R8P], [0.7_R8P, 0.65_R8P, 0.8_R8P])
! 2. half-space n . x <= d (n > 0) inside the block
n = [1._R8P, 2._R8P, 3._R8P] / sqrt(14._R8P) ; d = 0.9_R8P
p0 = d * n ; t1 = [2._R8P, -1._R8P, 0._R8P] / sqrt(5._R8P) ; t2 = [n(2) * t1(3) - n(3) * t1(2), n(3) * t1(1) - n(1) * t1(3), &
                                                                     n(1) * t1(2) - n(2) * t1(1)]
allocate(tet%v(3, 3, 4))
tet%v(:, :, 1) = reshape([p0 + 50._R8P * t1, p0 - 25._R8P * t1 + 43.3_R8P * t2, p0 - 25._R8P * t1 - 43.3_R8P * t2], [3, 3])
call orient(tet%v(:, :, 1), n) ! outward normal +n
block
   real(R8P) :: apex(3), a1(3), a2(3), a3(3)
   apex = p0 - 60._R8P * n
   a1 = tet%v(:, 1, 1) ; a2 = tet%v(:, 2, 1) ; a3 = tet%v(:, 3, 1)
   tet%v(:, :, 2) = reshape([a2, a1, apex], [3, 3])
   tet%v(:, :, 3) = reshape([a3, a2, apex], [3, 3])
   tet%v(:, :, 4) = reshape([a1, a3, apex], [3, 3])
end block
call to_surface(tet, s)
ok(2) = halfspace(s, u, w, n, d)
! 3. rotated cube, pyramid, cone: sums and divergence residual
cth = cos(0.3_R8P) ; sth = sin(0.3_R8P) ; cph = cos(0.7_R8P) ; sph = sin(0.7_R8P)
rot = matmul(reshape([cth, sth, 0._R8P, -sth, cth, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]), &
             reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, cph, sph, 0._R8P, -sph, cph], [3, 3]))
rc = b(1)%soup
do f = 1, size(rc%v, 3)
   do m = 1, 3
      rc%v(:, m, f) = 0.5_R8P + 0.6_R8P * matmul(rot, rc%v(:, m, f) - 0.5_R8P)
   enddo
enddo
call to_surface(rc, s)
ok(3) = sums(s, u, u, u, 0.216_R8P, 2.16_R8P, [0.5_R8P, 0.5_R8P, 0.5_R8P], 'rotated cube, uniform') .and. &
        sums(s, w, u, w, 0.216_R8P, 2.16_R8P, [0.5_R8P, 0.5_R8P, 0.5_R8P], 'rotated cube, stretched')
call to_surface(b(6)%soup, s)
ok(3) = ok(3) .and. sums(s, [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                         [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                         [(-0.1_R8P + 3.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], b(6)%volume, -1._R8P, [0._R8P, 0._R8P, -1._R8P], &
                         'pyramid 10')
call to_surface(b(7)%soup, s)
ok(3) = ok(3) .and. sums(s, [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                         [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                         [(-0.1_R8P + 2.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], b(7)%volume, -1._R8P, [0._R8P, 0._R8P, -1._R8P], &
                         'cone 15 (32)')
! 4. unit cube on grid planes
g = [(-0.25_R8P + (k - 1) * 0.25_R8P, k = 1, 9)]
call to_surface(b(1)%soup, s)
call cut_cells_3d(surface=s, x=g, y=g, z=g, cut=cut)
ok(4) = all(cut%closed) .and. all(abs(cut%alpha) < 1.e-14_R8P .or. abs(cut%alpha - 1._R8P) < 1.e-14_R8P) .and. &
        count(cut%alpha < 0.5_R8P) == 64 .and. abs(sum(1._R8P - cut%alpha) / 64._R8P - 1._R8P) < 1.e-13_R8P
print '(A,L1)', 'unit cube on grid planes: cells 0 or 1 to round-off, 64 solid, total exact: ', ok(4)
ok(4) = ok(4) .and. sums(s, g, g, g, 1._R8P, 6._R8P, [0.5_R8P, 0.5_R8P, 0.5_R8P], 'unit cube on grid planes')

! 5. solid cells without a wall get alpha exactly 0 (issue #103): the unit cube on a grid of spacing 0.1 (not
!    representable), whose interior cells have no facet; before the fix the face-sum volume left +-1.3e-16
q = [(-0.1_R8P + 0.1_R8P * (k - 1), k = 1, 13)]
call to_surface(b(1)%soup, s)
call cut_cells_3d(surface=s, x=q, y=q, z=q, cut=cut)
ok(5) = count(cut%pieces == 0) > 0 .and. count(cut%pieces == 0 .and. cut%alpha /= 0._R8P) == 0 .and. &
        count(cut%pieces == 1 .and. cut%alpha == 1._R8P) > 0
print '(A,I0,A,I0,A,L1)', 'unit cube on a 0.1 grid: solid cells ', count(cut%pieces == 0), ', with alpha /= 0 ', &
      count(cut%pieces == 0 .and. cut%alpha /= 0._R8P), ': ', ok(5)
! 6. an open box (the unit cube without its +z facets) on the same grid: cells whose faces disagree keep their
!    fraction and stay untrusted (#103 snaps only cells whose three faces agree); solid cells still exactly 0
m = 0
do f = 1, size(b(1)%soup%v, 3)
   if (any(b(1)%soup%v(3, :, f) < maxval(b(1)%soup%v(3, :, :)))) m = m + 1
enddo
if (allocated(rc%v)) deallocate(rc%v)
allocate(rc%v(3, 3, m)) ; m = 0
do f = 1, size(b(1)%soup%v, 3)
   if (any(b(1)%soup%v(3, :, f) < maxval(b(1)%soup%v(3, :, :)))) then
      m = m + 1 ; rc%v(:, :, m) = b(1)%soup%v(:, :, f)
   endif
enddo
call to_surface(rc, s)
call cut_cells_3d(surface=s, x=q, y=q, z=q, cut=cut)
ok(6) = count(cut%pieces == -1) == OPEN_UNTRUSTED .and. count(cut%pieces == 0 .and. cut%alpha /= 0._R8P) == 0
print '(A,I0,A,I0,A,L1)', 'open box on a 0.1 grid: untrusted cells ', count(cut%pieces == -1), ' (expected ', &
      OPEN_UNTRUSTED, '): ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine orient(t, nn)
   !< Make the triangle t outward along nn.
   real(R8P), intent(inout) :: t(3, 3) !< Triangle.
   real(R8P), intent(in)    :: nn(3)   !< Desired normal direction.
   real(R8P)                :: e1(3), e2(3), c(3), tmp(3) !< Work.

   e1 = t(:, 2) - t(:, 1) ; e2 = t(:, 3) - t(:, 1)
   c = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   if (dot_product(c, nn) < 0._R8P) then
      tmp = t(:, 2) ; t(:, 2) = t(:, 3) ; t(:, 3) = tmp
   endif
   endsubroutine orient

   function aabox(sf, x, z, lo, hi) result(good)
   !< Axis-aligned box: fractions and centroids against the overlaps.
   type(surface_stl_object), intent(in) :: sf             !< Box surface.
   real(R8P),                intent(in) :: x(:), z(:)     !< Edges (x and y use x).
   real(R8P),                intent(in) :: lo(3), hi(3)   !< Box.
   logical                              :: good           !< Result.
   real(R8P)    :: ea, ec, ov(3), ol(3), vs, vc, cs(3), cc(3), cf(3) !< Work.
   integer(I4P) :: i, j, k !< Work.

   call cut_cells_3d(surface=sf, x=x, y=x, z=z, cut=cut)
   ea = 0._R8P ; ec = 0._R8P
   do k = 1, size(z) - 1
      do j = 1, size(x) - 1
         do i = 1, size(x) - 1
            ol = max(lo, [x(i), x(j), z(k)]) ; ov = min(hi, [x(i + 1), x(j + 1), z(k + 1)]) - ol
            vs = product(max(ov, 0._R8P)) ; vc = (x(i + 1) - x(i)) * (x(j + 1) - x(j)) * (z(k + 1) - z(k))
            cs = ol + 0.5_R8P * max(ov, 0._R8P) ; cc = 0.5_R8P * [x(i) + x(i + 1), x(j) + x(j + 1), z(k) + z(k + 1)]
            ea = max(ea, abs(cut%alpha(i, j, k) - (vc - vs) / vc))
            if (vc - vs > 0._R8P) then
               cf = (cc * vc - cs * vs) / (vc - vs)
               ec = max(ec, cut%alpha(i, j, k) * norm2(cut%centroid(:, i, j, k) - cf))
            endif
         enddo
      enddo
   enddo
   good = all(cut%closed) .and. ea < 1.e-13_R8P .and. ec < 1.e-13_R8P
   print '(A,ES9.1,A,ES9.1,A,L1)', 'axis-aligned box: max err fraction ', ea, '  moment ', ec, '  ok ', good
   endfunction aabox

   function halfspace(sf, x, z, nn, dd) result(good)
   !< Half-space n . x <= d (n > 0): fractions against the closed-form box-plane volume.
   type(surface_stl_object), intent(in) :: sf         !< Tetrahedron.
   real(R8P),                intent(in) :: x(:), z(:) !< Edges.
   real(R8P),                intent(in) :: nn(3), dd  !< Plane.
   logical                              :: good       !< Result.
   real(R8P)    :: ea, vs, vc, lo(3), h(3), vtx(3) !< Work.
   integer(I4P) :: i, j, k, c1, c2, c3, ncut !< Work.

   call cut_cells_3d(surface=sf, x=x, y=x, z=z, cut=cut)
   ea = 0._R8P ; ncut = 0
   do k = 1, size(z) - 1
      do j = 1, size(x) - 1
         do i = 1, size(x) - 1
            lo = [x(i), x(j), z(k)] ; h = [x(i + 1) - x(i), x(j + 1) - x(j), z(k + 1) - z(k)] ; vc = product(h)
            vs = 0._R8P
            do c3 = 0, 1
               do c2 = 0, 1
                  do c1 = 0, 1
                     vtx = lo + [c1, c2, c3] * h
                     vs = vs + (-1)**(c1 + c2 + c3) * max(0._R8P, dd - dot_product(nn, vtx))**3
                  enddo
               enddo
            enddo
            vs = vs / (6._R8P * product(nn))
            if (vs > 0._R8P .and. vs < vc) ncut = ncut + 1
            ea = max(ea, abs(cut%alpha(i, j, k) - (vc - vs) / vc))
         enddo
      enddo
   enddo
   good = all(cut%closed) .and. ea < 1.e-12_R8P .and. ncut > 0
   print '(A,ES9.1,A,I0,A,L1)', 'half-space: max err fraction ', ea, '  cut cells ', ncut, '  ok ', good
   endfunction halfspace

   function sums(sf, x, y, z, vol, area, cen, label) result(good)
   !< Covering sums and per-cell divergence residual.
   type(surface_stl_object), intent(in) :: sf               !< Surface.
   real(R8P),                intent(in) :: x(:), y(:), z(:) !< Edges (covering the body).
   real(R8P),                intent(in) :: vol, area        !< Body volume, area (area < 0: not checked).
   real(R8P),                intent(in) :: cen(3)           !< Body centroid (cen(3) < 0: not checked).
   character(*),             intent(in) :: label            !< Label.
   logical                              :: good             !< Result.
   real(R8P)    :: vs, ms(3), wa, wn(3), res, h(3), vc, cc(3), r(3) !< Work.
   integer(I4P) :: i, j, k !< Work.

   call cut_cells_3d(surface=sf, x=x, y=y, z=z, cut=cut)
   vs = 0._R8P ; ms = 0._R8P ; wa = 0._R8P ; wn = 0._R8P ; res = 0._R8P
   do k = 1, size(z) - 1
      do j = 1, size(y) - 1
         do i = 1, size(x) - 1
            h = [x(i + 1) - x(i), y(j + 1) - y(j), z(k + 1) - z(k)] ; vc = product(h)
            cc = 0.5_R8P * [x(i) + x(i + 1), y(j) + y(j + 1), z(k) + z(k + 1)]
            vs = vs + (1._R8P - cut%alpha(i, j, k)) * vc
            ms = ms + cc * vc - cut%centroid(:, i, j, k) * cut%alpha(i, j, k) * vc
            wa = wa + cut%wet_area(i, j, k) ; wn = wn + cut%wet_normal(:, i, j, k)
            r = [(cut%face_x(i + 1, j, k) - cut%face_x(i, j, k)) * h(2) * h(3), &
                 (cut%face_y(i, j + 1, k) - cut%face_y(i, j, k)) * h(3) * h(1), &
                 (cut%face_z(i, j, k + 1) - cut%face_z(i, j, k)) * h(1) * h(2)] - cut%wet_normal(:, i, j, k)
            res = max(res, norm2(r))
         enddo
      enddo
   enddo
   good = all(cut%closed) .and. abs(vs - vol) < 1.e-13_R8P .and. res < 1.e-14_R8P .and. norm2(wn) < 1.e-14_R8P
   if (area > 0._R8P) good = good .and. abs(wa - area) < 1.e-13_R8P
   if (cen(3) >= 0._R8P) good = good .and. norm2(ms / vs - cen) < 1.e-13_R8P
   print '(A,A,ES9.1,A,ES9.1,A,ES9.1,A,ES9.1,A,L1)', label, ': volume err ', vs - vol, '  area err ', wa - abs(area), &
         '  sum of wall normals ', norm2(wn), '  max divergence residual ', res, '  ok ', good
   endfunction sums
endprogram fossil_test_cut_cells_3d
