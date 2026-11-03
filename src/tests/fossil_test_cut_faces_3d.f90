!< FOSSIL test: 3-D cut-cell faces (issue #40 step 4).
program fossil_test_cut_faces_3d
!< FOSSIL test: 3-D cut-cell faces (issue #40 step 4).
!<
!< Reference for a convex polyhedron: the face rectangle clipped by every facet half-space (area and
!< centroid of the planar polygon).
!< - A rotated cube, the 10-degree pyramid and the 32-face cone of the sign benchmark, on uniform and
!<   stretched 6^3 blocks: every face fluid fraction to round-off, the face fluid first moment to
!<   round-off, every face closed.
!< - The unit cube with its faces on grid planes: a face is solid exactly when its rectangle lies
!<   in the cube and its plane in [0, 1) (the plane at x_i + eps), fluid otherwise.
!< - Bad edges refused.

use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, STATUS_INVALID_INPUT
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(soup_t)              :: rc
type(surface_stl_object)  :: s
type(cut_cells_3d_t)      :: cut
real(R8P)                 :: g(9), rot(3, 3), u(7), w(7), cth, sth, cph, sph
integer(I4P)              :: k, f, m, st
logical                   :: ok(3)

call bench_bodies(b)
! rotated cube about the centre (0.5, 0.5, 0.5)
cth = cos(0.3_R8P) ; sth = sin(0.3_R8P) ; cph = cos(0.7_R8P) ; sph = sin(0.7_R8P)
rot = matmul(reshape([cth, sth, 0._R8P, -sth, cth, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]), &
             reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, cph, sph, 0._R8P, -sph, cph], [3, 3]))
rc = b(1)%soup
do f = 1, size(rc%v, 3)
   do m = 1, 3
      rc%v(:, m, f) = 0.5_R8P + 0.6_R8P * matmul(rot, rc%v(:, m, f) - 0.5_R8P)
   enddo
enddo
do k = 1, 7
   u(k) = -0.1_R8P + 1.2_R8P * (k - 1) / 6._R8P
   w(k) = -0.1_R8P + 1.2_R8P * ((k - 1) / 6._R8P)**1.3_R8P
enddo
ok(1) = .true.
call to_surface(rc, s)
ok(1) = ok(1) .and. convex(s, rc, u, u, u, 'rotated cube, uniform') .and. convex(s, rc, w, u, w, 'rotated cube, stretched')
call to_surface(b(6)%soup, s)
ok(1) = ok(1) .and. convex(s, b(6)%soup, [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                           [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                           [(-0.1_R8P + 3.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], 'pyramid 10')
call to_surface(b(7)%soup, s)
ok(1) = ok(1) .and. convex(s, b(7)%soup, [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                           [(-0.6_R8P + 1.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], &
                           [(-0.1_R8P + 2.2_R8P * (k - 1) / 6._R8P, k = 1, 7)], 'cone 15 (32)')
! unit cube on grid planes
g = [(-0.25_R8P + (k - 1) * 0.25_R8P, k = 1, 9)]
call to_surface(b(1)%soup, s)
call cut_cells_3d(surface=s, x=g, y=g, z=g, cut=cut)
ok(2) = aligned(cut, g)
print '(A,L1)', 'unit cube on grid planes: every face 0 or 1 as the eps rule says: ', ok(2)
call cut_cells_3d(surface=s, x=[0._R8P, 0._R8P], y=g, z=g, cut=cut, status=st)
ok(3) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'bad edges refused: ', ok(3)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function convex(sf, sp, x, y, z, label) result(good)
   !< Face fractions and moments of a convex body against the half-space clipping reference.
   type(surface_stl_object), intent(in) :: sf       !< Surface.
   type(soup_t),             intent(in) :: sp       !< Its soup (facet planes).
   real(R8P),                intent(in) :: x(:), y(:), z(:) !< Edges.
   character(*),             intent(in) :: label    !< Label.
   logical                              :: good     !< Result.
   real(R8P)    :: ea, em, rect(3, 4) !< Work.
   integer(I4P) :: i, j, k, ncut !< Work.
   logical      :: allc !< Work.

   call cut_cells_3d(surface=sf, x=x, y=y, z=z, cut=cut)
   ea = 0._R8P ; em = 0._R8P ; ncut = 0
   allc = all(cut%face_x_closed) .and. all(cut%face_y_closed) .and. all(cut%face_z_closed)
   do k = 1, size(z) - 1
      do j = 1, size(y) - 1
         do i = 1, size(x)
            rect = reshape([x(i), y(j), z(k), x(i), y(j + 1), z(k), x(i), y(j + 1), z(k + 1), x(i), y(j), z(k + 1)], [3, 4])
            call check(rect, cut%face_x(i, j, k), cut%face_x_centroid(:, i, j, k), sp, ea, em, ncut)
         enddo
      enddo
   enddo
   do k = 1, size(z) - 1
      do j = 1, size(y)
         do i = 1, size(x) - 1
            rect = reshape([x(i), y(j), z(k), x(i + 1), y(j), z(k), x(i + 1), y(j), z(k + 1), x(i), y(j), z(k + 1)], [3, 4])
            call check(rect, cut%face_y(i, j, k), cut%face_y_centroid(:, i, j, k), sp, ea, em, ncut)
         enddo
      enddo
   enddo
   do k = 1, size(z)
      do j = 1, size(y) - 1
         do i = 1, size(x) - 1
            rect = reshape([x(i), y(j), z(k), x(i + 1), y(j), z(k), x(i + 1), y(j + 1), z(k), x(i), y(j + 1), z(k)], [3, 4])
            call check(rect, cut%face_z(i, j, k), cut%face_z_centroid(:, i, j, k), sp, ea, em, ncut)
         enddo
      enddo
   enddo
   good = allc .and. ea < 1.e-13_R8P .and. em < 1.e-13_R8P .and. ncut > 0
   print '(A,A,ES9.1,A,ES9.1,A,I0,A,L1)', label, ': max err fraction ', ea, '  moment ', em, '  cut faces ', ncut, '  ok ', good
   endfunction convex

   subroutine check(r, frac, centroid, sp, ea, em, ncut)
   !< Compare one face with the half-space clipping reference.
   real(R8P),    intent(in)    :: r(3, 4), frac, centroid(3) !< Face rectangle, FOSSIL fraction and centroid.
   type(soup_t), intent(in)    :: sp       !< Facets.
   real(R8P),    intent(inout) :: ea, em   !< Max errors.
   integer(I4P), intent(inout) :: ncut     !< Cut faces.
   real(R8P) :: q(3, 64), rc_(3), fa, fc(3), ar, af, cen(3) !< Work.
   integer   :: nq !< Work.

   fa = poly_area3(r, 4) ; fc = poly_centroid3(r, 4)
   call clip_halfspaces(r, sp, q, nq)
   ar = 0._R8P ; rc_ = fc
   if (nq >= 3) then
      ar = poly_area3(q, nq) ; rc_ = poly_centroid3(q, nq)
   endif
   af = fa - ar ! fluid area
   if (ar > 0._R8P .and. af > 0._R8P) ncut = ncut + 1
   ea = max(ea, abs(frac - af / fa))
   if (af > 0._R8P) then
      cen = (fc * fa - rc_ * ar) / af
      em = max(em, frac * norm2(centroid - cen))
   endif
   endsubroutine check

   subroutine clip_halfspaces(r, sp, q, nq)
   !< Clip the planar polygon r by the half-spaces n.(p - v1) <= 0 of every facet of sp.
   real(R8P),    intent(in)  :: r(3, 4) !< Polygon.
   type(soup_t), intent(in)  :: sp      !< Facets.
   real(R8P),    intent(out) :: q(3, 64) !< Clipped polygon.
   integer,      intent(out) :: nq      !< Vertices.
   real(R8P) :: wv(3, 64), o(3, 64), n(3), sa, sb !< Work.
   integer   :: nw, no, f, k, kk !< Work.

   nw = 4 ; wv(:, 1:4) = r
   do f = 1, size(sp%v, 3)
      n = cross(sp%v(:, 2, f) - sp%v(:, 1, f), sp%v(:, 3, f) - sp%v(:, 1, f))
      no = 0
      do k = 1, nw
         kk = mod(k, nw) + 1
         sa = -dot_product(n, wv(:, k) - sp%v(:, 1, f)) ; sb = -dot_product(n, wv(:, kk) - sp%v(:, 1, f))
         if (sa >= 0._R8P) then
            no = no + 1 ; o(:, no) = wv(:, k)
         endif
         if ((sa >= 0._R8P) .neqv. (sb >= 0._R8P)) then
            no = no + 1 ; o(:, no) = wv(:, k) + sa / (sa - sb) * (wv(:, kk) - wv(:, k))
         endif
      enddo
      nw = no ; wv(:, 1:nw) = o(:, 1:nw)
      if (nw == 0) exit
   enddo
   nq = nw ; q(:, 1:nq) = wv(:, 1:nq)
   endsubroutine clip_halfspaces

   pure real(R8P) function poly_area3(p, n)
   !< Area of a planar polygon.
   real(R8P), intent(in) :: p(:,:) !< Polygon.
   integer,   intent(in) :: n      !< Vertices.
   real(R8P)             :: a(3)   !< Vector area.
   integer               :: k      !< Counter.

   a = 0._R8P
   do k = 2, n - 1
      a = a + cross(p(:, k) - p(:, 1), p(:, k + 1) - p(:, 1))
   enddo
   poly_area3 = 0.5_R8P * norm2(a)
   endfunction poly_area3

   pure function poly_centroid3(p, n) result(c)
   !< Centroid of a planar convex polygon (fan of triangles).
   real(R8P), intent(in) :: p(:,:) !< Polygon.
   integer,   intent(in) :: n      !< Vertices.
   real(R8P)             :: c(3)   !< Centroid.
   real(R8P)             :: a, at  !< Areas.
   integer               :: k      !< Counter.

   c = 0._R8P ; a = 0._R8P
   do k = 2, n - 1
      at = 0.5_R8P * norm2(cross(p(:, k) - p(:, 1), p(:, k + 1) - p(:, 1)))
      c = c + at * (p(:, 1) + p(:, k) + p(:, k + 1)) / 3._R8P ; a = a + at
   enddo
   if (a > 0._R8P) then
      c = c / a
   else
      c = p(:, 1)
   endif
   endfunction poly_centroid3

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   function aligned(cc, g) result(good)
   !< Unit cube on grid planes: face solid iff its plane is in [0, 1) and its rectangle in the cube.
   type(cut_cells_3d_t), intent(in) :: cc   !< Cut cells.
   real(R8P),            intent(in) :: g(:) !< Edges (all axes).
   logical                          :: good !< Result.
   integer :: i, j, k, nbad !< Work.
   logical :: solid !< Expected.

   nbad = 0
   do k = 1, size(g) - 1
      do j = 1, size(g) - 1
         do i = 1, size(g)
            solid = g(i) >= 0._R8P .and. g(i) < 1._R8P .and. g(j) >= 0._R8P .and. g(j + 1) <= 1._R8P .and. &
                    g(k) >= 0._R8P .and. g(k + 1) <= 1._R8P
            if (cc%face_x(i, j, k) /= merge(0._R8P, 1._R8P, solid)) nbad = nbad + 1
            if (cc%face_y(j, i, k) /= merge(0._R8P, 1._R8P, solid)) nbad = nbad + 1
            if (cc%face_z(j, k, i) /= merge(0._R8P, 1._R8P, solid)) nbad = nbad + 1
         enddo
      enddo
   enddo
   good = nbad == 0 .and. all(cc%face_x_closed) .and. all(cc%face_y_closed) .and. all(cc%face_z_closed)
   if (nbad > 0) print '(A,I0)', '   aligned face errors: ', nbad
   endfunction aligned
endprogram fossil_test_cut_faces_3d
