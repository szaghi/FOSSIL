!< FOSSIL test: exact cell geometry of curvilinear blocks (issue #71, step 3).
program fossil_test_block_geometry
!< FOSSIL test: exact cell geometry of curvilinear blocks (issue #71, step 3).
!<
!< 1. affine map x = A X: every cell volume is det(A) h^3 and every centroid A X_c (1e-14 relative);
!< 2. the trilinear map F(X) = (X + a Y Z, Y + b Z X, Z + c X Y) is its own trilinear interpolant: cell volumes equal
!<    the analytic integral of det DF (1e-14), centroids the 3-point Gauss integral (exact for degree 5) (1e-14);
!< 3. GCL: the outward face vectors of every cell of a wavy O-grid sum to zero (1e-14 of the largest face);
!< 4. divergence theorem: the cell volumes of the wavy O-grid sum to (1/3) oint x . n dA over the block boundary
!<    (2 x 2 Gauss per bilinear face, exact) (1e-13 relative);
!< 5. the (theta, r, z) O-grid is left-handed (every volume negative); the mirrored block has exactly the negated volumes;
!< 6. 2-D: affine areas and centroids; the bilinear map (X + a X Y, Y + b X Y) analytic areas; edge closure; the block
!<    area equals the shoelace area of its boundary.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t
use fossil_block_geometry,    only : trilinear_point, trilinear_jacobian, det3, GAUSS2
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
real(R8P), parameter            :: A(3, 3) = reshape([1.2_R8P, 0.1_R8P, -0.2_R8P, 0.3_R8P, 0.9_R8P, 0.15_R8P, &
                                                      -0.1_R8P, 0.25_R8P, 1.1_R8P], [3, 3])
type(curvilinear_block_t)       :: b, m
type(curvilinear_block_2d_t)    :: q
real(R8P), allocatable          :: x(:,:,:), y(:,:,:), z(:,:,:), u(:,:), v(:,:)
real(R8P), allocatable          :: vol(:,:,:), cen(:,:,:,:), fi(:,:,:,:), fj(:,:,:,:), fk(:,:,:,:), vol2(:,:,:)
real(R8P), allocatable          :: ar(:,:), cen2(:,:,:), ei(:,:,:), ej(:,:,:)
real(R8P)                       :: h, err, xc(3), s(3), ex, fmax, aa, bb, cc, lo(3), hi(3), tot, flux
integer(I4P)                    :: i, j, k, n, ni, nj, nk
logical                         :: all_ok

all_ok = .true.

! 1. affine
n = 6 ; h = 0.2_R8P
allocate(x(n, n, n), y(n, n, n), z(n, n, n))
do k = 1, n
   do j = 1, n
      do i = 1, n
         s = matmul(A, h * [i - 1, j - 1, k - 1]) + [0.3_R8P, -0.2_R8P, 0.1_R8P]
         x(i, j, k) = s(1) ; y(i, j, k) = s(2) ; z(i, j, k) = s(3)
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%geometry(volume=vol, centroid=cen)
ex = det3(A) * h ** 3 ; err = maxval(abs(vol - ex)) / abs(ex)
call check(err <= 1e-14_R8P, 'affine volumes', err)
err = 0._R8P
do k = 1, n - 1
   do j = 1, n - 1
      do i = 1, n - 1
         xc = matmul(A, h * [i - 0.5_R8P, j - 0.5_R8P, k - 0.5_R8P]) + [0.3_R8P, -0.2_R8P, 0.1_R8P]
         err = max(err, maxval(abs(cen(:, i, j, k) - xc)))
      enddo
   enddo
enddo
call check(err <= 1e-14_R8P, 'affine centroids', err)

! 2. trilinear map
aa = 0.4_R8P ; bb = -0.3_R8P ; cc = 0.5_R8P
do k = 1, n
   do j = 1, n
      do i = 1, n
         s = h * [i - 1, j - 1, k - 1] - 0.5_R8P
         x(i, j, k) = s(1) + aa * s(2) * s(3) ; y(i, j, k) = s(2) + bb * s(3) * s(1) ; z(i, j, k) = s(3) + cc * s(1) * s(2)
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%geometry(volume=vol, centroid=cen)
err = 0._R8P ; tot = 0._R8P
do k = 1, n - 1
   do j = 1, n - 1
      do i = 1, n - 1
         lo = h * [i - 1, j - 1, k - 1] - 0.5_R8P ; hi = lo + h
         ex = h ** 3 - bb * cc * (hi(1) ** 3 - lo(1) ** 3) / 3 * h * h - aa * cc * (hi(2) ** 3 - lo(2) ** 3) / 3 * h * h &
              - aa * bb * (hi(3) ** 3 - lo(3) ** 3) / 3 * h * h &
              + 2 * aa * bb * cc * (hi(1) ** 2 - lo(1) ** 2) * (hi(2) ** 2 - lo(2) ** 2) * (hi(3) ** 2 - lo(3) ** 2) / 8
         err = max(err, abs(vol(i, j, k) - ex) / ex)
         tot = max(tot, maxval(abs(cen(:, i, j, k) - gauss3_centroid(b%corners(i, j, k)))))
      enddo
   enddo
enddo
call check(err <= 1e-14_R8P, 'trilinear map: analytic volumes', err)
call check(tot <= 1e-14_R8P, 'trilinear map: centroids = 3-point Gauss', tot)
deallocate(x, y, z)

! 3-4. wavy O-grid
ni = 25 ; nj = 7 ; nk = 5
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         h = (1._R8P + 0.15_R8P * (j - 1)) * (1._R8P + 0.05_R8P * sin(3 * 2 * PI * (i - 1) / (ni - 1) + 0.4_R8P * k))
         x(i, j, k) = h * cos(2 * PI * (i - 1) / (ni - 1)) ; y(i, j, k) = h * sin(2 * PI * (i - 1) / (ni - 1))
         z(i, j, k) = 0.2_R8P * (k - 1) + 0.03_R8P * cos(2 * PI * (i - 1) / (ni - 1)) * j
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%geometry(volume=vol, face_i=fi, face_j=fj, face_k=fk)
err = 0._R8P ; fmax = max(maxval(abs(fi)), maxval(abs(fj)), maxval(abs(fk)))
do k = 1, nk - 1
   do j = 1, nj - 1
      do i = 1, ni - 1
         s = fi(:, i + 1, j, k) - fi(:, i, j, k) + fj(:, i, j + 1, k) - fj(:, i, j, k) + fk(:, i, j, k + 1) - fk(:, i, j, k)
         err = max(err, maxval(abs(s)) / fmax)
      enddo
   enddo
enddo
call check(err <= 1e-14_R8P, 'GCL: outward face vectors sum to zero', err)
flux = boundary_flux(b)
err = abs(sum(vol) - flux) / abs(flux)
call check(err <= 1e-13_R8P, 'divergence theorem: block volume', err)
call check(all(vol < 0._R8P), '(theta, r, z) O-grid is left-handed: every volume negative', 0._R8P)

! 5. mirrored
call m%initialize(-x, y, z) ; call m%geometry(volume=vol2)
call check(all(vol2 == -vol), 'mirrored block: volumes exactly negated', 0._R8P)

! 6. 2-D
ni = 7 ; nj = 6 ; h = 0.25_R8P
allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   do i = 1, ni
      u(i, j) = 1.1_R8P * h * (i - 1) + 0.3_R8P * h * (j - 1) ; v(i, j) = -0.2_R8P * h * (i - 1) + 0.9_R8P * h * (j - 1)
   enddo
enddo
call q%initialize(u, v) ; call q%geometry(area=ar, centroid=cen2, face_i=ei, face_j=ej)
ex = (1.1_R8P * 0.9_R8P + 0.2_R8P * 0.3_R8P) * h * h
err = maxval(abs(ar - ex)) / ex
call check(err <= 1e-14_R8P, '2-D affine areas', err)
err = 0._R8P
do j = 1, nj - 1
   do i = 1, ni - 1
      xc(1:2) = [1.1_R8P * h * (i - 0.5_R8P) + 0.3_R8P * h * (j - 0.5_R8P), -0.2_R8P * h * (i - 0.5_R8P) + &
                 0.9_R8P * h * (j - 0.5_R8P)]
      err = max(err, maxval(abs(cen2(:, i, j) - xc(1:2))))
   enddo
enddo
call check(err <= 1e-14_R8P, '2-D affine centroids', err)
aa = 0.6_R8P ; bb = -0.4_R8P
do j = 1, nj
   do i = 1, ni
      s(1:2) = [h * (i - 1), h * (j - 1)]
      u(i, j) = s(1) + aa * s(1) * s(2) ; v(i, j) = s(2) + bb * s(1) * s(2)
   enddo
enddo
call q%initialize(u, v) ; call q%geometry(area=ar, face_i=ei, face_j=ej)
err = 0._R8P ; tot = 0._R8P ; fmax = max(maxval(abs(ei)), maxval(abs(ej)))
do j = 1, nj - 1
   do i = 1, ni - 1
      lo(1:2) = [h * (i - 1), h * (j - 1)] ; hi(1:2) = lo(1:2) + h
      ex = h * h + aa * h * (hi(2) ** 2 - lo(2) ** 2) / 2 + bb * h * (hi(1) ** 2 - lo(1) ** 2) / 2
      err = max(err, abs(ar(i, j) - ex) / ex)
      s(1:2) = ei(:, i + 1, j) - ei(:, i, j) + ej(:, i, j + 1) - ej(:, i, j)
      tot = max(tot, maxval(abs(s(1:2))) / fmax)
   enddo
enddo
call check(err <= 1e-14_R8P, '2-D bilinear map: analytic areas', err)
call check(tot <= 1e-15_R8P, '2-D edge normals close', tot)
flux = 0._R8P
do i = 1, ni - 1
   flux = flux + cross2(u(i, 1), v(i, 1), u(i + 1, 1), v(i + 1, 1)) + cross2(u(i + 1, nj), v(i + 1, nj), u(i, nj), v(i, nj))
enddo
do j = 1, nj - 1
   flux = flux + cross2(u(ni, j), v(ni, j), u(ni, j + 1), v(ni, j + 1)) + cross2(u(1, j + 1), v(1, j + 1), u(1, j), v(1, j))
enddo
flux = flux / 2
err = abs(sum(ar) - flux) / flux
call check(err <= 1e-14_R8P, '2-D block area = shoelace of the boundary', err)

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES10.2)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   pure real(R8P) function cross2(ax, ay, bx, by) result(c)
   !< a x b.
   real(R8P), intent(in) :: ax, ay, bx, by !< Vectors.

   c = ax * by - ay * bx
   endfunction cross2

   function gauss3_centroid(c) result(xc)
   !< Centroid by 3-point Gauss per direction.
   real(R8P), intent(in) :: c(3, 8) !< Corners.
   real(R8P)             :: xc(3)   !< Centroid.
   real(R8P)             :: g(3), w(3), d, vv, xi(3) !< Rule, work.
   integer(I4P)          :: a1, a2, a3 !< Counters.

   g = [0.5_R8P - sqrt(0.15_R8P), 0.5_R8P, 0.5_R8P + sqrt(0.15_R8P)] ; w = [5._R8P, 8._R8P, 5._R8P] / 18._R8P
   xc = 0._R8P ; vv = 0._R8P
   do a3 = 1, 3
      do a2 = 1, 3
         do a1 = 1, 3
            xi = [g(a1), g(a2), g(a3)] ; d = det3(trilinear_jacobian(c, xi)) * w(a1) * w(a2) * w(a3)
            xc = xc + d * trilinear_point(c, xi) ; vv = vv + d
         enddo
      enddo
   enddo
   xc = xc / vv
   endfunction gauss3_centroid

   function boundary_flux(b) result(f)
   !< (1/3) oint x . n dA over the six faces of the block, 2 x 2 Gauss per bilinear face.
   type(curvilinear_block_t), intent(in) :: b !< Block.
   real(R8P)                             :: f !< Volume.
   integer(I4P)                          :: i, j, k !< Counters.

   f = 0._R8P
   do k = 1, b%nk - 1
      do j = 1, b%nj - 1
         f = f + face_flux(b, [b%ni, j, k], [0, 1, 0], [0, 0, 1]) - face_flux(b, [1, j, k], [0, 1, 0], [0, 0, 1])
      enddo
   enddo
   do k = 1, b%nk - 1
      do i = 1, b%ni - 1
         f = f + face_flux(b, [i, b%nj, k], [0, 0, 1], [1, 0, 0]) - face_flux(b, [i, 1, k], [0, 0, 1], [1, 0, 0])
      enddo
   enddo
   do j = 1, b%nj - 1
      do i = 1, b%ni - 1
         f = f + face_flux(b, [i, j, b%nk], [1, 0, 0], [0, 1, 0]) - face_flux(b, [i, j, 1], [1, 0, 0], [0, 1, 0])
      enddo
   enddo
   f = f / 3
   endfunction boundary_flux

   function face_flux(b, p, d1, d2) result(f)
   !< int x . (x_s x x_t) ds dt over the bilinear face at vertex p spanned by d1, d2.
   type(curvilinear_block_t), intent(in) :: b             !< Block.
   integer(I4P),              intent(in) :: p(3), d1(3), d2(3) !< Corner, directions.
   real(R8P)                             :: f             !< Flux.
   real(R8P)                             :: c00(3), c10(3), c01(3), c11(3), xs(3), xt(3), xp(3), n(3) !< Work.
   integer(I4P)                          :: a1, a2        !< Gauss.
   real(R8P)                             :: s, t          !< Gauss nodes.

   c00 = pt(b, p) ; c10 = pt(b, p + d1) ; c01 = pt(b, p + d2) ; c11 = pt(b, p + d1 + d2)
   f = 0._R8P
   do a2 = 1, 2
      do a1 = 1, 2
         s = GAUSS2(a1) ; t = GAUSS2(a2)
         xp = c00 * (1 - s) * (1 - t) + c10 * s * (1 - t) + c01 * (1 - s) * t + c11 * s * t
         xs = (c10 - c00) * (1 - t) + (c11 - c01) * t ; xt = (c01 - c00) * (1 - s) + (c11 - c10) * s
         n = [xs(2) * xt(3) - xs(3) * xt(2), xs(3) * xt(1) - xs(1) * xt(3), xs(1) * xt(2) - xs(2) * xt(1)]
         f = f + dot_product(xp, n) / 4
      enddo
   enddo
   endfunction face_flux

   function pt(b, p) result(c)
   !< Vertex p.
   type(curvilinear_block_t), intent(in) :: b    !< Block.
   integer(I4P),              intent(in) :: p(3) !< Index.
   real(R8P)                             :: c(3) !< Coordinates.

   c = [b%x(p(1), p(2), p(3)), b%y(p(1), p(2), p(3)), b%z(p(1), p(2), p(3))]
   endfunction pt
endprogram fossil_test_block_geometry
