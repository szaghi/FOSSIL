!< FOSSIL test: exact 2-D cut cells (issue #40 step 3).
program fossil_test_cut_cells_2d
!< FOSSIL test: exact 2-D cut cells (issue #40 step 3).
!<
!< Reference: Sutherland-Hodgman clipping of the convex solid polygon by the cell (area, centroid),
!< Liang-Barsky clipping of the cell edges and of the polygon sides (apertures, wetted length).
!< - An offset square, a rotated square, a 128-gon on 8 x 8 and stretched blocks: alpha, centroid,
!<   apertures and wetted length equal the reference to round-off (the centroid as the fluid first
!<   moment, alpha times the centroid error: a nearly solid cell amplifies round-off in both); the sum of the fluid areas equals
!<   the block area minus the polygon area; the per-cell divergence residual (aperture-weighted edge
!<   normals minus the wall-normal integral) is at round-off; piece areas add up to alpha.
!< - A diamond with vertices on grid nodes and sides through grid nodes, against the reference too.
!< - Every case: a shared edge has the same aperture seen from both cells.
!< - Grid-aligned square (sides on grid lines, where a closed-set reference is ambiguous): every cell
!<   is 0 or 1, the sum exact, all closed, the wetted lengths sum to the perimeter, shared apertures
!<   agree.
!< - A thin solid strip across a row of cells: those cells have 2 fluid pieces with the right areas.
!< - Loops inside one cell: an island (1 piece, alpha = 1 - island), a hole in a solid (1 piece,
!<   alpha = hole), an island inside that hole (1 piece, alpha = hole - island).
!< - An open polyline: the cells it crosses are flagged not closed. Bad edges refused.

use fossil, only : curve_object, cut_cells_2d_t, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: c
type(cut_cells_2d_t) :: cut
real(R8P)            :: e8(9), s8(9), poly(2, 128), u(12), v(12)
integer(I4P)         :: k, st
logical              :: ok(6)

do k = 1, 9
   e8(k) = (k - 1) / 8._R8P
   s8(k) = ((k - 1) / 8._R8P)**1.4_R8P
enddo
! 1. convex polygons against the clipping reference
ok(1) = .true.
call square(0.31_R8P, 0.29_R8P, 0.43_R8P, 0._R8P, poly(:, 1:4))
ok(1) = ok(1) .and. convex(poly(:, 1:4), e8, e8, 'offset square, uniform')
call square(0.5_R8P, 0.5_R8P, 0.6_R8P, 0.37_R8P, poly(:, 1:4))
ok(1) = ok(1) .and. convex(poly(:, 1:4), e8, e8, 'rotated square, uniform')
ok(1) = ok(1) .and. convex(poly(:, 1:4), s8, e8, 'rotated square, stretched')
poly(:, 1:4) = reshape([0.5_R8P, 0.25_R8P, 0.75_R8P, 0.5_R8P, 0.5_R8P, 0.75_R8P, 0.25_R8P, 0.5_R8P], [2, 4])
ok(1) = ok(1) .and. convex(poly(:, 1:4), e8, e8, 'diamond through grid nodes')
do k = 1, 128
   poly(:, k) = [0.5_R8P + 0.37_R8P * cos(2 * PI * (k - 1) / 128 + 0.1_R8P), &
                 0.5_R8P + 0.37_R8P * sin(2 * PI * (k - 1) / 128 + 0.1_R8P)]
enddo
ok(1) = ok(1) .and. convex(poly, e8, e8, '128-gon, uniform') .and. convex(poly, s8, s8, '128-gon, stretched')
! 2. grid-aligned square
call c%adopt_loops(u=[0.25_R8P, 0.75_R8P, 0.75_R8P, 0.25_R8P], v=[0.25_R8P, 0.25_R8P, 0.75_R8P, 0.75_R8P], loop_sizes=[4_I4P])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(2) = all(cut%closed) .and. all(cut%alpha == 0._R8P .or. cut%alpha == 1._R8P) .and. &
        abs(sum(cut%alpha) / 64._R8P - (1._R8P - 0.25_R8P)) < 1.e-15_R8P .and. abs(sum(cut%wet_length) - 2._R8P) < 1.e-15_R8P &
        .and. shared_ok()
print '(A,L1,A,I0)', 'grid-aligned square: cells 0 or 1, exact sum: ', ok(2), '  solid cells ', count(cut%alpha == 0._R8P)
! 3. thin strip across row 4 (v in [0.375, 0.5]): solid v in [0.42, 0.43]
call c%adopt_loops(u=[-1._R8P, 2._R8P, 2._R8P, -1._R8P], v=[0.42_R8P, 0.42_R8P, 0.43_R8P, 0.43_R8P], loop_sizes=[4_I4P])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(3) = all(cut%pieces(:, 4) == 2) .and. all(abs(cut%alpha(:, 4) - (1._R8P - 0.01_R8P * 8._R8P)) < 1.e-13_R8P)
do k = 1, 8
   ok(3) = ok(3) .and. abs(cut%piece_area(cut%first_piece(k, 4)) - 0.125_R8P * (0.42_R8P - 0.375_R8P)) < 1.e-15_R8P .and. &
           abs(cut%piece_area(cut%first_piece(k, 4) + 1) - 0.125_R8P * (0.5_R8P - 0.43_R8P)) < 1.e-15_R8P
enddo
print '(A,L1)', 'thin strip: 2 fluid pieces per cell of the row, areas right: ', ok(3)
! 4. loops inside one cell (cell (3, 3) = [0.25, 0.375]^2)
call c%adopt_loops(u=[0.28_R8P, 0.33_R8P, 0.33_R8P, 0.28_R8P], v=[0.28_R8P, 0.28_R8P, 0.33_R8P, 0.33_R8P], loop_sizes=[4_I4P])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(4) = cut%pieces(3, 3) == 1 .and. abs(cut%alpha(3, 3) - (1._R8P - 0.0025_R8P / 0.015625_R8P)) < 1.e-13_R8P .and. &
        count(cut%alpha == 1._R8P) == 63
u(1:8) = [-1._R8P, 2._R8P, 2._R8P, -1._R8P, 0.28_R8P, 0.28_R8P, 0.34_R8P, 0.34_R8P]
v(1:8) = [-1._R8P, -1._R8P, 2._R8P, 2._R8P, 0.28_R8P, 0.34_R8P, 0.34_R8P, 0.28_R8P]
call c%adopt_loops(u=u(1:8), v=v(1:8), loop_sizes=[4_I4P, 4_I4P])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(4) = ok(4) .and. cut%pieces(3, 3) == 1 .and. abs(cut%alpha(3, 3) - 0.0036_R8P / 0.015625_R8P) < 1.e-13_R8P .and. &
        count(cut%alpha == 0._R8P) == 63
u = [-1._R8P, 2._R8P, 2._R8P, -1._R8P, 0.28_R8P, 0.28_R8P, 0.34_R8P, 0.34_R8P, 0.30_R8P, 0.32_R8P, 0.32_R8P, 0.30_R8P]
v = [-1._R8P, -1._R8P, 2._R8P, 2._R8P, 0.28_R8P, 0.34_R8P, 0.34_R8P, 0.28_R8P, 0.30_R8P, 0.30_R8P, 0.32_R8P, 0.32_R8P]
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 4_I4P, 4_I4P])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(4) = ok(4) .and. cut%pieces(3, 3) == 1 .and. abs(cut%alpha(3, 3) - (0.0036_R8P - 0.0004_R8P) / 0.015625_R8P) < 1.e-13_R8P
print '(A,L1)', 'island, hole, island in hole inside one cell: ', ok(4)
! 5. open polyline
call c%adopt_loops(u=[0.1_R8P, 0.9_R8P, 0.9_R8P], v=[0.1_R8P, 0.1_R8P, 0.9_R8P], loop_sizes=[3_I4P], closed=[.false.])
call c%cut_cells(x=e8, y=e8, cut=cut)
ok(5) = count(.not. cut%closed) > 0 .and. all(cut%alpha(1:7, 8) >= 0._R8P)
print '(A,L1,A,I0)', 'open polyline: crossed cells flagged not closed: ', ok(5), '  flagged ', count(.not. cut%closed)
! 6. refusal
call c%cut_cells(x=[0._R8P, 0._R8P], y=e8, cut=cut, status=st)
ok(6) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'bad edges refused: ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine square(cu, cv, side, rot, p)
   !< CCW square of centre (cu, cv), side `side`, rotated by rot.
   real(R8P), intent(in)  :: cu, cv, side, rot !< Centre, side, rotation.
   real(R8P), intent(out) :: p(2, 4)           !< Vertices.
   integer                :: k                 !< Counter.
   real(R8P)              :: a                 !< Angle.

   do k = 1, 4
      a = rot + PI / 4._R8P + (k - 1) * PI / 2._R8P
      p(:, k) = [cu, cv] + side / sqrt(2._R8P) * [cos(a), sin(a)]
   enddo
   endsubroutine square

   function convex(p, x, y, label) result(good)
   !< Cut cells of the convex polygon p against the clipping reference.
   real(R8P),    intent(in) :: p(:,:)     !< Polygon (CCW).
   real(R8P),    intent(in) :: x(:), y(:) !< Cell edges.
   character(*), intent(in) :: label      !< Label.
   logical                  :: good       !< Result.
   real(R8P)    :: ar, cen(2), q(2, 64), ap(4), wl, a_poly, tot, res, emax(5), cl(2), hx, hy, sp !< Work.
   integer(I4P) :: i, j, nq, k, nbadp !< Work.

   call c%adopt_loops(u=p(1, :), v=p(2, :), loop_sizes=[size(p, 2, kind=I4P)])
   call c%cut_cells(x=x, y=y, cut=cut)
   a_poly = poly_area(p, size(p, 2)) ; tot = 0._R8P ; emax = 0._R8P ; nbadp = 0
   do j = 1, size(y) - 1
      do i = 1, size(x) - 1
         hx = x(i + 1) - x(i) ; hy = y(j + 1) - y(j)
         call clip(p, size(p, 2), x(i), x(i + 1), y(j), y(j + 1), q, nq)
         ar = 0._R8P ; cen = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (y(j) + y(j + 1))]
         if (nq >= 3) then
            ar = poly_area(q, nq) ; cl = poly_centroid(q, nq)
         endif
         ! fluid = cell - solid
         if (hx * hy - ar > 0._R8P) cen = (cen * hx * hy - cl * ar) / (hx * hy - ar)
         if (nq < 3) cen = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (y(j) + y(j + 1))]
         ap(1) = 1._R8P - seg_in(p, [x(i), y(j)], [x(i), y(j + 1)]) / hy
         ap(2) = 1._R8P - seg_in(p, [x(i + 1), y(j)], [x(i + 1), y(j + 1)]) / hy
         ap(3) = 1._R8P - seg_in(p, [x(i), y(j)], [x(i + 1), y(j)]) / hx
         ap(4) = 1._R8P - seg_in(p, [x(i), y(j + 1)], [x(i + 1), y(j + 1)]) / hx
         wl = wet_ref(p, x(i), x(i + 1), y(j), y(j + 1))
         emax(1) = max(emax(1), abs(cut%alpha(i, j) - (1._R8P - ar / (hx * hy))))
         if (cut%alpha(i, j) > 0._R8P) emax(2) = max(emax(2), cut%alpha(i, j) * norm2(cut%centroid(:, i, j) - cen))
         emax(3) = max(emax(3), maxval(abs(cut%aperture(:, i, j) - ap)))
         emax(4) = max(emax(4), abs(cut%wet_length(i, j) - wl))
         res = norm2([cut%aperture(2, i, j) * hy - cut%aperture(1, i, j) * hy, &
                      cut%aperture(4, i, j) * hx - cut%aperture(3, i, j) * hx] - cut%wet_normal(:, i, j))
         emax(5) = max(emax(5), res)
         tot = tot + cut%alpha(i, j) * hx * hy
         sp = 0._R8P
         do k = cut%first_piece(i, j), cut%first_piece(i, j) + cut%pieces(i, j) - 1
            sp = sp + cut%piece_area(k)
         enddo
         if (abs(sp - cut%alpha(i, j) * hx * hy) > 1.e-15_R8P) nbadp = nbadp + 1
      enddo
   enddo
   good = all(cut%closed) .and. emax(1) < 1.e-13_R8P .and. emax(2) < 1.e-13_R8P .and. emax(3) < 1.e-13_R8P .and. &
          emax(4) < 1.e-13_R8P .and. emax(5) < 1.e-14_R8P .and. nbadp == 0 .and. &
          abs(tot - ((x(size(x)) - x(1)) * (y(size(y)) - y(1)) - a_poly)) < 1.e-14_R8P .and. shared_ok()
   print '(A,A,5ES9.1,A,ES9.1,A,I0,A,L1)', label, ': max err alpha/moment/aperture/wet/divergence ', emax, &
         '  sum err ', tot - ((x(size(x)) - x(1)) * (y(size(y)) - y(1)) - a_poly), '  piece sums bad ', nbadp, '  ok ', good
   endfunction convex

   logical function shared_ok()
   !< Shared edges have the same aperture from both cells.
   integer :: i, j !< Counters.

   shared_ok = .true.
   do j = 1, cut%ny
      do i = 1, cut%nx
         if (i < cut%nx) shared_ok = shared_ok .and. abs(cut%aperture(2, i, j) - cut%aperture(1, i + 1, j)) < 1.e-15_R8P
         if (j < cut%ny) shared_ok = shared_ok .and. abs(cut%aperture(4, i, j) - cut%aperture(3, i, j + 1)) < 1.e-15_R8P
      enddo
   enddo
   endfunction shared_ok

   subroutine clip(p, n, x0, x1, y0, y1, q, nq)
   !< Sutherland-Hodgman clipping of polygon p by the box.
   real(R8P),    intent(in)  :: p(:,:)          !< Polygon.
   integer,      intent(in)  :: n               !< Vertices.
   real(R8P),    intent(in)  :: x0, x1, y0, y1  !< Box.
   real(R8P),    intent(out) :: q(:,:)          !< Clipped polygon.
   integer(I4P), intent(out) :: nq              !< Its vertices.
   real(R8P)    :: w(2, 300), o(2, 300), sa, sb !< Work.
   integer      :: nw, no, k, kk, side          !< Work.

   nw = n ; w(:, 1:n) = p(:, 1:n)
   do side = 1, 4
      no = 0
      do k = 1, nw
         kk = mod(k, nw) + 1
         sa = meas(w(:, k), side, x0, x1, y0, y1) ; sb = meas(w(:, kk), side, x0, x1, y0, y1)
         if (sa >= 0._R8P) then
            no = no + 1 ; o(:, no) = w(:, k)
         endif
         if ((sa >= 0._R8P) .neqv. (sb >= 0._R8P)) then
            no = no + 1 ; o(:, no) = w(:, k) + sa / (sa - sb) * (w(:, kk) - w(:, k))
         endif
      enddo
      nw = no ; w(:, 1:nw) = o(:, 1:nw)
      if (nw == 0) exit
   enddo
   nq = nw ; q(:, 1:nq) = w(:, 1:nq)
   endsubroutine clip

   pure real(R8P) function meas(pt, sd, x0, x1, y0, y1)
   !< Inside measure of pt for side sd of the box (>= 0 inside).
   real(R8P), intent(in) :: pt(2)          !< Point.
   integer,   intent(in) :: sd             !< Side.
   real(R8P), intent(in) :: x0, x1, y0, y1 !< Box.

   select case (sd)
   case (1) ; meas = pt(1) - x0
   case (2) ; meas = x1 - pt(1)
   case (3) ; meas = pt(2) - y0
   case default ; meas = y1 - pt(2)
   endselect
   endfunction meas

   pure real(R8P) function poly_area(p, n)
   !< Signed area.
   real(R8P), intent(in) :: p(:,:) !< Polygon.
   integer,   intent(in) :: n      !< Vertices.
   integer               :: k      !< Counter.

   poly_area = 0._R8P
   do k = 1, n
      poly_area = poly_area + 0.5_R8P * (p(1, k) * p(2, mod(k, n) + 1) - p(1, mod(k, n) + 1) * p(2, k))
   enddo
   endfunction poly_area

   pure function poly_centroid(p, n) result(cc)
   !< Centroid.
   real(R8P), intent(in) :: p(:,:) !< Polygon.
   integer,   intent(in) :: n      !< Vertices.
   real(R8P)             :: cc(2)  !< Centroid.
   real(R8P)             :: cr, a  !< Work.
   integer               :: k, kk  !< Counters.

   cc = 0._R8P ; a = 0._R8P
   do k = 1, n
      kk = mod(k, n) + 1
      cr = p(1, k) * p(2, kk) - p(1, kk) * p(2, k)
      a = a + cr ; cc = cc + cr * (p(:, k) + p(:, kk))
   enddo
   cc = cc / (3._R8P * a)
   endfunction poly_centroid

   pure real(R8P) function seg_in(p, a, b)
   !< Length of the segment ab inside the convex polygon p (Cyrus-Beck).
   real(R8P), intent(in) :: p(:,:), a(2), b(2) !< Polygon, segment.
   real(R8P)             :: t0, t1, nrm(2), num, den !< Work.
   integer               :: k, kk !< Counters.

   t0 = 0._R8P ; t1 = 1._R8P
   do k = 1, size(p, 2)
      kk = mod(k, size(p, 2)) + 1
      nrm = [p(2, kk) - p(2, k), -(p(1, kk) - p(1, k))] ! outward
      num = dot_product(nrm, p(:, k) - a) ; den = dot_product(nrm, b - a)
      if (den == 0._R8P) then
         if (num < 0._R8P) then
            seg_in = 0._R8P ; return
         endif
      elseif (den > 0._R8P) then
         t1 = min(t1, num / den)
      else
         t0 = max(t0, num / den)
      endif
   enddo
   seg_in = max(0._R8P, t1 - t0) * norm2(b - a)
   endfunction seg_in

   pure real(R8P) function wet_ref(p, x0, x1, y0, y1)
   !< Length of the polygon boundary inside the box (Liang-Barsky per side).
   real(R8P), intent(in) :: p(:,:), x0, x1, y0, y1 !< Polygon, box.
   integer               :: k, kk, m !< Counters.
   real(R8P)             :: t0, t1, d(2), pp(4), qq(4) !< Work.

   wet_ref = 0._R8P
   do k = 1, size(p, 2)
      kk = mod(k, size(p, 2)) + 1
      d = p(:, kk) - p(:, k) ; t0 = 0._R8P ; t1 = 1._R8P
      pp = [-d(1), d(1), -d(2), d(2)] ; qq = [p(1, k) - x0, x1 - p(1, k), p(2, k) - y0, y1 - p(2, k)]
      do m = 1, 4
         if (pp(m) == 0._R8P) then
            if (qq(m) < 0._R8P) t1 = -1._R8P
         elseif (pp(m) < 0._R8P) then
            t0 = max(t0, qq(m) / pp(m))
         else
            t1 = min(t1, qq(m) / pp(m))
         endif
      enddo
      if (t1 > t0) wet_ref = wet_ref + (t1 - t0) * norm2(d)
   enddo
   endfunction wet_ref
endprogram fossil_test_cut_cells_2d
