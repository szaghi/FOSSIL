!< FOSSIL test: the intersection weights of a curvilinear receptor block against a Cartesian donor grid in 3-D (issue
!< #106 step 3, Picard et al. 2021).
program fossil_test_overlap_weights_3d
!< FOSSIL test: the intersection weights of a curvilinear receptor block against a Cartesian donor grid in 3-D (issue
!< #106 step 3, Picard et al. 2021).
!<
!< `overlap_weights_3d` on the two outer layers of a sphere shell (r in [1, 2], poles collapsed, a left-handed block)
!< over the donor grid [-2.6, 2.6]^3, and the same shifted by 100 h; a receptor on donor grid planes (2 x 2 x 2 donor
!< cells), one leaving the donor box; refusals: a folded cell, an inverted cell, a cell outside the box, a cell out of
!< range. Checks per receptor (its polyhedron rebuilt here: 12 triangles, shorter diagonals, ties by the least global
!< vertex index):
!< - sum of the weights = 1 to 2 eps (compensated sum);
!< - receptor_volume = the polyhedron's volume to 1e-13 relative (inside the box);
!< - sum weight x centroid = the polyhedron's centroid to 1e-13 of its size;
!< - a constant field is exact (4 eps); a linear field is exact with Picard's gradient term to 1e-13 of its variation;
!< - covered = 1 exactly inside the donor box, < 1 when the receptor leaves it; every volume >= 0, weight in [0, 1].
use fossil, only : curvilinear_block_t, overlap_weights_t, overlap_weights_3d, OVERLAP_WEIGHTS_OK, &
                   OVERLAP_WEIGHTS_BAD_INPUT, OVERLAP_WEIGHTS_DEGENERATE, OVERLAP_WEIGHTS_OUTSIDE
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)    !< Pi.
real(R8P), parameter      :: EPS = epsilon(1._R8P) !< Machine epsilon.
type(curvilinear_block_t) :: blk                   !< Receptor block.
type(overlap_weights_t)   :: w                     !< Weights.
real(R8P), allocatable    :: x(:)                  !< Donor grid.
real(R8P)                 :: s, h, v(3, 8)         !< Shift, cell size, hex corners.
integer(I4P), allocatable :: cells(:,:)            !< Receptor cells.
integer(I4P)              :: k, st, i, j, n        !< Counters.
logical                   :: ok                    !< Checks.

ok = .true.
! 1. the two outer layers of a 12 x 6 x 4 shell, at the origin and shifted by 100 h
h = 5.2_R8P / 16
do k = 0, 1
   s = k * 100._R8P * h
   call shell(12_I4P, 4_I4P, s, blk)
   x = s + grid_1d(-2.6_R8P, 2.6_R8P, 16_I4P)
   allocate(cells(3, 2 * 12 * 6)) ; n = 0
   do st = 3, 4
      do j = 1, 6
         do i = 1, 12
            n = n + 1 ; cells(:, n) = [i, j, st]
         enddo
      enddo
   enddo
   call check_case('shell outer layers, shift '//trim(itoa(k * 100))//' h', blk, cells, x, .true., ok)
   deallocate(cells)
enddo
x = grid_1d(-2.6_R8P, 2.6_R8P, 16_I4P) ; h = x(2) - x(1)
! 2. a receptor on donor grid planes: 2 x 2 x 2 donor cells, vertices on nodes
call hex_corners([x(3), x(5)], [x(2), x(4)], [x(6), x(8)], v) ; call hex_block(v, blk)
call check_case('receptor = 2 x 2 x 2 donor cells', blk, reshape([1_I4P, 1_I4P, 1_I4P], [3, 1]), x, .true., ok, &
                ndonors=8_I4P)
ok = ok .and. all(abs(w%weight - 0.125_R8P) <= 4 * EPS)
! 3. a sheared receptor leaving the donor box across x = 2.6
call hex_corners([2.45_R8P, 2.8_R8P], [0.1_R8P, 0.4_R8P], [-0.2_R8P, 0.15_R8P], v)
v(2, 2) = v(2, 2) + 0.05_R8P ; v(3, 8) = v(3, 8) + 0.07_R8P
call hex_block(v, blk)
call check_case('receptor leaving the donor box', blk, reshape([1_I4P, 1_I4P, 1_I4P], [3, 1]), x, .false., ok)
! refusals: 4. a folded cell (two corners swapped), an inverted cell next to a valid one, a cell outside the box
call hex_corners([0.1_R8P, 0.5_R8P], [0.1_R8P, 0.5_R8P], [0.1_R8P, 0.5_R8P], v)
v(:, [7, 8]) = v(:, [8, 7]) ; call hex_block(v, blk)
call overlap_weights_3d(blk, reshape([1_I4P, 1_I4P, 1_I4P], [3, 1]), x, x, x, w, status=st)
print '(A,I0,A,I0)', '  folded cell: status ', st, ', flag ', w%flag(1)
ok = ok .and. st == OVERLAP_WEIGHTS_DEGENERATE .and. w%flag(1) == OVERLAP_WEIGHTS_DEGENERATE .and. w%nw == 0 .and. &
     w%receptor_volume(1) == 0._R8P .and. w%covered(1) == 0._R8P
call two_cells(blk)
call overlap_weights_3d(blk, reshape([1_I4P, 1_I4P, 1_I4P, 2_I4P, 1_I4P, 1_I4P], [3, 2]), x, x, x, w, status=st)
print '(A,I0,A,2(I0,1X),A,ES10.2)', '  inverted cell next to a valid one: status ', st, ', flags ', w%flag, &
      '; valid cell |sum w - 1| ', abs(sum(w%weight) - 1._R8P)
ok = ok .and. st == OVERLAP_WEIGHTS_DEGENERATE .and. w%flag(1) == OVERLAP_WEIGHTS_OK .and. &
     w%flag(2) == OVERLAP_WEIGHTS_DEGENERATE .and. w%start(3) == w%start(2) .and. w%start(2) > 1 .and. &
     abs(sum(w%weight) - 1._R8P) <= 4 * EPS .and. abs(w%receptor_volume(1) - 0.064_R8P) <= 1e-15_R8P
call hex_corners([2.7_R8P, 3._R8P], [0.1_R8P, 0.4_R8P], [0.1_R8P, 0.4_R8P], v) ; call hex_block(v, blk)
call overlap_weights_3d(blk, reshape([1_I4P, 1_I4P, 1_I4P], [3, 1]), x, x, x, w, status=st)
print '(A,I0,A,I0,A,ES10.2)', '  cell outside the donor box: status ', st, ', flag ', w%flag(1), ', covered ', w%covered(1)
ok = ok .and. st == OVERLAP_WEIGHTS_OK .and. w%flag(1) == OVERLAP_WEIGHTS_OUTSIDE .and. w%nw == 0 .and. &
     w%covered(1) == 0._R8P
call overlap_weights_3d(blk, reshape([1_I4P, 1_I4P, 2_I4P], [3, 1]), x, x, x, w, status=st)
print '(A,I0)', '  cell out of range: status ', st
ok = ok .and. st == OVERLAP_WEIGHTS_BAD_INPUT .and. w%nr == 0
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check_case(tag, blk, cells, x, inside, ok, ndonors)
   !< overlap_weights_3d on the receptor cells, with the checks of the program notes.
   character(*),              intent(in)    :: tag        !< Case.
   type(curvilinear_block_t), intent(in)    :: blk        !< Receptor block.
   integer(I4P),              intent(in)    :: cells(:,:) !< Receptor cells.
   real(R8P),                 intent(in)    :: x(:)       !< Donor grid edges (all directions).
   logical,                   intent(in)    :: inside     !< Every receptor lies in the donor box.
   logical,                   intent(inout) :: ok         !< Checks.
   integer(I4P),    optional, intent(in)    :: ndonors    !< Expected donors of the first receptor.
   real(R8P), parameter :: C0 = 3.7_R8P, B(3) = [0.7_R8P, -2.1_R8P, 1.3_R8P] !< Fields.
   real(R8P)    :: vol, cen(3), diam, ws, wc, cs, cc, qs, qc, mom(3), xc(3), ref, e(5), tri(3, 3, 12), o(3) !< Work.
   integer(I4P) :: st, r, k, ncov, nbad, nt, nmin                                        !< Work.

   call overlap_weights_3d(blk, cells, x, x, x, w, status=st)
   e = 0._R8P ; ncov = 0 ; nbad = 0 ; nmin = 0
   do r = 1, w%nr
      if (w%flag(r) /= OVERLAP_WEIGHTS_OK) then
         nbad = nbad + 1 ; cycle
      endif
      call receptor(blk, cells(:, r), tri, nt)
      o = tri(:, 1, 1)
      call polyhedron(tri, nt, o, vol, cen, diam)
      ws = 0._R8P ; wc = 0._R8P ; cs = 0._R8P ; cc = 0._R8P ; qs = 0._R8P ; qc = 0._R8P ; mom = 0._R8P
      do k = w%start(r), w%start(r + 1) - 1
         call neumaier(ws, wc, w%weight(k))
         call neumaier(cs, cc, w%weight(k) * C0)
         xc = 0.5_R8P * [x(w%donor(1, k)) + x(w%donor(1, k) + 1), x(w%donor(2, k)) + x(w%donor(2, k) + 1), &
                         x(w%donor(3, k)) + x(w%donor(3, k) + 1)]
         call neumaier(qs, qc, w%weight(k) * (lin(xc - o) + dot_product(B, w%centroid(:, k) - xc)))
         mom = mom + w%weight(k) * (w%centroid(:, k) - o)
      enddo
      ws = ws + wc ; cs = cs + cc ; qs = qs + qc
      ref = lin(cen)
      e(1) = max(e(1), abs(ws - 1._R8P))
      e(3) = max(e(3), abs(cs - C0) / C0)
      if (inside) then
         e(2) = max(e(2), abs(w%receptor_volume(r) - vol) / vol)
         e(4) = max(e(4), norm2(mom - cen) / diam)
         e(5) = max(e(5), abs(qs - ref) / (norm2(B) * diam + abs(ref)))
      endif
      if (inside .neqv. w%covered(r) == 1._R8P) ncov = ncov + 1
      if (any(.not. (abs(w%weight(w%start(r):w%start(r + 1) - 1)) <= 1._R8P))) nbad = nbad + 1
   enddo
   print '(3A,I0,A,I0,A,F6.2,A,I0,A)', '  ', tag, ': ', w%nr, ' receptors, ', w%nw, ' intersections (', &
         real(w%nw, R8P) / max(1, w%nr), ' each; ', count(w%volume <= 1e-12_R8P * (x(2) - x(1)) ** 3), ' below 1e-12 h^3)'
   nbad = nbad + count(.not. (w%volume >= 0._R8P .and. w%weight >= 0._R8P .and. w%weight <= 1._R8P))
   print '(A,5ES10.2)', '    |sum w - 1|, volume, constant, centroid / size, linear: ', e
   ok = ok .and. st == OVERLAP_WEIGHTS_OK .and. nbad == 0 .and. ncov == 0 .and. e(1) <= 2 * EPS .and. &
        e(2) <= 1e-13_R8P .and. e(3) <= 4 * EPS .and. e(4) <= 1e-13_R8P .and. e(5) <= 1e-13_R8P
   if (present(ndonors)) ok = ok .and. w%start(2) - w%start(1) == ndonors
   if (.not. inside) ok = ok .and. all(w%covered < 1._R8P) .and. all(w%covered > 0._R8P)
   endsubroutine check_case

   pure real(R8P) function lin(p)
   !< A linear field (of coordinates relative to the receptor's first vertex).
   real(R8P), intent(in) :: p(3) !< Point.

   lin = 3.7_R8P + 0.7_R8P * p(1) - 2.1_R8P * p(2) + 1.3_R8P * p(3)
   endfunction lin

   subroutine receptor(blk, c, tri, nt)
   !< The 12 triangles of receptor cell c (built here: faces out of the cell by the block's handedness, shorter diagonal,
   !< ties by the least global vertex index, collapsed triangles dropped).
   type(curvilinear_block_t), intent(in)  :: blk           !< Block.
   integer(I4P),              intent(in)  :: c(3)          !< Cell.
   real(R8P),                 intent(out) :: tri(3, 3, 12) !< Triangles.
   integer(I4P),              intent(out) :: nt            !< Kept.
   integer(I4P), parameter :: FQ(4, 6) = reshape([1, 5, 7, 3,  2, 4, 8, 6,  1, 2, 6, 5,  3, 7, 8, 4,  1, 3, 4, 2, &
                                                  5, 6, 8, 7], [4, 6]) !< Faces out of a right-handed cell.
   real(R8P)    :: v(3, 8), q(3, 4), t(3, 3, 2), d13, d24, n(3) !< Work.
   integer(I4P) :: g(8), gq(4), m, f, ii, jj, kk               !< Work.
   logical      :: use13                                       !< Diagonal.

   do m = 1, 8
      ii = c(1) + mod(m - 1, 2) ; jj = c(2) + mod((m - 1) / 2, 2) ; kk = c(3) + (m - 1) / 4
      v(:, m) = [blk%x(ii, jj, kk), blk%y(ii, jj, kk), blk%z(ii, jj, kk)] ; g(m) = ii + blk%ni * ((jj - 1) + blk%nj * (kk - 1))
   enddo
   nt = 0
   do f = 1, 6
      q = v(:, FQ(:, f)) ; gq = g(FQ(:, f))
      if (blk%orientation() < 0._R8P) then
         q = q(:, [1, 4, 3, 2]) ; gq = gq([1, 4, 3, 2])
      endif
      d13 = norm2(q(:, 3) - q(:, 1)) ; d24 = norm2(q(:, 4) - q(:, 2))
      use13 = d13 < d24 .or. (d13 == d24 .and. min(gq(1), gq(3)) <= min(gq(2), gq(4)))
      if (use13) then
         t(:, :, 1) = q(:, [1, 2, 3]) ; t(:, :, 2) = q(:, [1, 3, 4])
      else
         t(:, :, 1) = q(:, [1, 2, 4]) ; t(:, :, 2) = q(:, [2, 3, 4])
      endif
      do m = 1, 2
         n = cross(t(:, 2, m) - t(:, 1, m), t(:, 3, m) - t(:, 1, m))
         if (norm2(n) == 0._R8P) cycle
         nt = nt + 1 ; tri(:, :, nt) = t(:, :, m)
      enddo
   enddo
   endsubroutine receptor

   pure subroutine polyhedron(tri, nt, o, vol, cen, diam)
   !< Volume, centroid (relative to o) and diameter of the closed triangulated polyhedron.
   real(R8P),    intent(in)  :: tri(3, 3, 12) !< Triangles (outward).
   integer(I4P), intent(in)  :: nt            !< Triangles.
   real(R8P),    intent(in)  :: o(3)          !< Origin.
   real(R8P),    intent(out) :: vol, cen(3), diam !< Volume, centroid - o, diameter.
   real(R8P)                 :: a(3), b(3), c(3), dv, lo(3), hi(3) !< Work.
   integer(I4P)              :: m, k          !< Work.

   vol = 0._R8P ; cen = 0._R8P ; lo = huge(1._R8P) ; hi = -huge(1._R8P)
   do m = 1, nt
      a = tri(:, 1, m) - o ; b = tri(:, 2, m) - o ; c = tri(:, 3, m) - o
      dv = dot_product(a, cross(b, c)) / 6._R8P
      vol = vol + dv ; cen = cen + dv * (a + b + c) / 4._R8P
      do k = 1, 3
         lo = min(lo, tri(:, k, m)) ; hi = max(hi, tri(:, k, m))
      enddo
   enddo
   cen = cen / vol ; diam = norm2(hi - lo)
   endsubroutine polyhedron

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   pure subroutine neumaier(s, comp, a)
   !< s += a with Neumaier's compensation.
   real(R8P), intent(inout) :: s, comp !< Sum, compensation.
   real(R8P), intent(in)    :: a       !< Term.
   real(R8P)                :: t       !< Work.

   t = s + a
   if (abs(s) >= abs(a)) then
      comp = comp + ((s - t) + a)
   else
      comp = comp + ((a - t) + s)
   endif
   s = t
   endsubroutine neumaier

   subroutine shell(na, nl, s, blk)
   !< A sphere shell block centred at (s, s, s), r in [1, 2] over nl cells, na x na / 2 cells on each sphere (i around,
   !< j pole to pole, k outward: left-handed), poles collapsed.
   integer(I4P),              intent(in)    :: na, nl !< Cells around, radially.
   real(R8P),                 intent(in)    :: s      !< Shift.
   type(curvilinear_block_t), intent(inout) :: blk    !< Block.
   real(R8P), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:) !< Vertices.
   real(R8P)              :: th, ph, rr                     !< Work.
   integer(I4P)           :: np, i, j, k                    !< Counters.

   np = na / 2
   allocate(xx(na + 1, np + 1, nl + 1), yy(na + 1, np + 1, nl + 1), zz(na + 1, np + 1, nl + 1))
   do k = 1, nl + 1
      rr = 1._R8P + (k - 1._R8P) / nl
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            xx(i, j, k) = s + rr * sin(ph) * cos(th) ; yy(i, j, k) = s + rr * sin(ph) * sin(th) ; zz(i, j, k) = s + rr * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            xx(:, j, k) = s ; yy(:, j, k) = s ; zz(:, j, k) = s + merge(rr, -rr, j == 1)
         endif
      enddo
      xx(na + 1, :, k) = xx(1, :, k) ; yy(na + 1, :, k) = yy(1, :, k) ; zz(na + 1, :, k) = zz(1, :, k)
   enddo
   call blk%destroy
   call blk%initialize(xx, yy, zz)
   endsubroutine shell

   pure subroutine hex_corners(a, b, c, v)
   !< The corners of the box [a(1), a(2)] x [b(1), b(2)] x [c(1), c(2)], v(:, 1 + bits(x, y, z)).
   real(R8P), intent(in)  :: a(2), b(2), c(2) !< Bounds.
   real(R8P), intent(out) :: v(3, 8)          !< Corners.
   integer(I4P)           :: m                !< Corner.

   do m = 1, 8
      v(:, m) = [a(1 + mod(m - 1, 2)), b(1 + mod((m - 1) / 2, 2)), c(1 + (m - 1) / 4)]
   enddo
   endsubroutine hex_corners

   subroutine hex_block(v, blk)
   !< A one-cell block of corners v(:, 1 + bits(i, j, k)).
   real(R8P),                 intent(in)    :: v(3, 8) !< Corners.
   type(curvilinear_block_t), intent(inout) :: blk     !< Block.

   call blk%destroy
   call blk%initialize(reshape(v(1, :), [2, 2, 2]), reshape(v(2, :), [2, 2, 2]), reshape(v(3, :), [2, 2, 2]))
   endsubroutine hex_block

   subroutine two_cells(blk)
   !< Two cells along i: a valid cube [0.1, 0.5]^3 and an inverted one (its far face pulled back across the shared face).
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P) :: xx(3, 2, 2), yy(3, 2, 2), zz(3, 2, 2) !< Vertices.

   xx(1, :, :) = 0.1_R8P ; xx(2, :, :) = 0.5_R8P ; xx(3, :, :) = 0.3_R8P
   yy(:, 1, :) = 0.1_R8P ; yy(:, 2, :) = 0.5_R8P ; zz(:, :, 1) = 0.1_R8P ; zz(:, :, 2) = 0.5_R8P
   call blk%destroy
   call blk%initialize(xx, yy, zz)
   endsubroutine two_cells

   pure function grid_1d(a, b, n) result(x)
   !< n + 1 equally spaced edges on [a, b].
   real(R8P),    intent(in) :: a, b     !< Ends.
   integer(I4P), intent(in) :: n        !< Cells.
   real(R8P)                :: x(n + 1) !< Edges.
   integer(I4P)             :: i        !< Counter.

   do i = 1, n + 1
      x(i) = a + (b - a) * (i - 1) / n
   enddo
   x(n + 1) = b
   endfunction grid_1d

   pure function itoa(i) result(s)
   !< Integer to string.
   integer, intent(in) :: i !< Integer.
   character(len=16)   :: s !< String.

   write(s, '(I0)') i
   endfunction itoa
endprogram fossil_test_overlap_weights_3d
