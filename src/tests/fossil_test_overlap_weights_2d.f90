!< FOSSIL test: the intersection weights of a curvilinear receptor block against a Cartesian donor grid (issue #106
!< step 2, Picard et al. 2021).
program fossil_test_overlap_weights_2d
!< FOSSIL test: the intersection weights of a curvilinear receptor block against a Cartesian donor grid (issue #106
!< step 2, Picard et al. 2021).
!<
!< `overlap_weights_2d` on the two outer layers of a cylinder O-grid (r in [1, 2]) over the donor grid [-3, 3]^2 at two
!< sizes, and the same shifted far from the origin (|x| / h ~ 100); a receptor on donor grid lines (2 x 2 donor cells),
!< one through donor nodes (a diamond), a triangle (a collapsed edge), one leaving the donor box; refusals: a folded
!< cell, an inverted cell, a collapsed one, one outside the box, a cell out of range. Checks per receptor:
!< - sum of the weights = 1 to 2 eps (compensated sum);
!< - receptor_volume = the polygon's area to 1e-14 relative (inside the box), = its clipped area outside;
!< - sum weight x centroid = the receptor's centroid to 1e-13 of its size;
!< - a constant field is exact (4 eps); a linear field is exact with Picard's gradient term,
!<   q_c = sum w (q_c' + grad q . (x_cc' - x_c')), to 1e-13 of its variation over the receptor;
!< - covered = 1 exactly inside the donor box, < 1 exactly when the receptor leaves it;
!< - every volume >= 0 and every weight in [0, 1].
use fossil, only : curvilinear_block_2d_t, overlap_weights_t, overlap_weights_2d, OVERLAP_WEIGHTS_OK, &
                   OVERLAP_WEIGHTS_BAD_INPUT, OVERLAP_WEIGHTS_DEGENERATE, OVERLAP_WEIGHTS_OUTSIDE
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)  !< Pi.
real(R8P), parameter :: EPS = epsilon(1._R8P) !< Machine epsilon.
type(curvilinear_block_2d_t) :: blk         !< Receptor block.
type(overlap_weights_t)      :: w           !< Weights.
real(R8P), allocatable       :: x(:)        !< Donor grid.
real(R8P)                    :: s, h        !< Shift, cell size.
integer(I4P), allocatable    :: cells(:,:)  !< Receptor cells.
integer(I4P)                 :: m, k, st, i !< Counters.
logical                      :: ok          !< Checks.

ok = .true.
! 1. the two outer layers of the O-grid, at the origin and shifted by 100 h
do m = 1, 2
   h = 6._R8P / (40 * m)
   do k = 0, 1
      s = k * 100._R8P * h
      call ogrid(m, s, blk)
      x = s + grid_1d(-3._R8P, 3._R8P, 40_I4P * m)
      allocate(cells(2, 2 * (blk%ni - 1)))
      do i = 1, blk%ni - 1
         cells(:, i) = [i, blk%nj - 2_I4P] ; cells(:, blk%ni - 1 + i) = [i, blk%nj - 1_I4P]
      enddo
      call check_case('O-grid outer layers, m = '//trim(itoa(m))//', shift '//trim(itoa(k * 100))//' h', blk, cells, x, &
                      .true., ok)
      deallocate(cells)
   enddo
enddo
x = grid_1d(-3._R8P, 3._R8P, 40_I4P) ; h = x(2) - x(1)
! 2. a receptor on donor grid lines: 2 x 2 donor cells, vertices on nodes
call quad_block([x(3), x(5), x(5), x(3)], [x(2), x(2), x(4), x(4)], blk)
call check_case('receptor = 2 x 2 donor cells', blk, reshape([1_I4P, 1_I4P], [2, 1]), x, .true., ok, ndonors=4_I4P)
ok = ok .and. all(abs(w%weight - 0.25_R8P) <= 4 * EPS)
! 3. a diamond with its vertices on nodes and its sides through nodes
call quad_block([x(21), x(25), x(21), x(17)], [x(17), x(21), x(25), x(21)], blk)
call check_case('diamond through donor nodes', blk, reshape([1_I4P, 1_I4P], [2, 1]), x, .true., ok)
! 4. a triangle: a collapsed edge (a pole)
call quad_block([0.13_R8P, 0.52_R8P, 0.52_R8P, 0.31_R8P], [0.07_R8P, 0.11_R8P, 0.11_R8P, 0.48_R8P], blk)
call check_case('triangle (collapsed edge)', blk, reshape([1_I4P, 1_I4P], [2, 1]), x, .true., ok)
! 5. a receptor leaving the donor box across x = 3
call quad_block([2.9_R8P, 3.2_R8P, 3.15_R8P, 2.85_R8P], [0.1_R8P, 0.12_R8P, 0.4_R8P, 0.35_R8P], blk)
call check_case('receptor leaving the donor box', blk, reshape([1_I4P, 1_I4P], [2, 1]), x, .false., ok)
! refusals: 6. a folded cell (bow tie), an inverted cell next to a valid one, a collapsed cell, a cell outside the box
call quad_block([0.1_R8P, 0.5_R8P, 0.1_R8P, 0.5_R8P], [0.1_R8P, 0.1_R8P, 0.5_R8P, 0.5_R8P], blk)
call overlap_weights_2d(blk, reshape([1_I4P, 1_I4P], [2, 1]), x, x, w, status=st)
print '(A,I0,A,I0)', '  folded cell: status ', st, ', flag ', w%flag(1)
ok = ok .and. st == OVERLAP_WEIGHTS_DEGENERATE .and. w%flag(1) == OVERLAP_WEIGHTS_DEGENERATE .and. w%nw == 0 .and. &
     w%receptor_volume(1) == 0._R8P .and. w%covered(1) == 0._R8P
call blk%destroy
call blk%initialize(reshape([0.1_R8P, 0.5_R8P, 0.3_R8P, 0.1_R8P, 0.5_R8P, 0.3_R8P], [3, 2]), &
                    reshape([0.1_R8P, 0.1_R8P, 0.1_R8P, 0.6_R8P, 0.6_R8P, 0.6_R8P], [3, 2]))
call overlap_weights_2d(blk, reshape([1_I4P, 1_I4P, 2_I4P, 1_I4P], [2, 2]), x, x, w, status=st)
print '(A,I0,A,2(I0,1X),A,ES10.2)', '  inverted cell next to a valid one: status ', st, ', flags ', w%flag, &
      '; valid cell |sum w - 1| ', abs(sum(w%weight) - 1._R8P)
ok = ok .and. st == OVERLAP_WEIGHTS_DEGENERATE .and. w%flag(1) == OVERLAP_WEIGHTS_OK .and. &
     w%flag(2) == OVERLAP_WEIGHTS_DEGENERATE .and. w%start(3) == w%start(2) .and. w%start(2) > 1 .and. &
     abs(sum(w%weight) - 1._R8P) <= 4 * EPS .and. abs(w%receptor_volume(1) - 0.2_R8P) <= 1e-15_R8P
call quad_block([0.1_R8P, 0.3_R8P, 0.5_R8P, 0.4_R8P], [0.1_R8P, 0.2_R8P, 0.3_R8P, 0.25_R8P], blk)
call overlap_weights_2d(blk, reshape([1_I4P, 1_I4P], [2, 1]), x, x, w, status=st)
print '(A,I0,A,I0)', '  collapsed cell (collinear vertices): status ', st, ', flag ', w%flag(1)
ok = ok .and. st == OVERLAP_WEIGHTS_DEGENERATE .and. w%flag(1) == OVERLAP_WEIGHTS_DEGENERATE .and. w%nw == 0
call quad_block([3.1_R8P, 3.4_R8P, 3.4_R8P, 3.1_R8P], [0.1_R8P, 0.1_R8P, 0.4_R8P, 0.4_R8P], blk)
call overlap_weights_2d(blk, reshape([1_I4P, 1_I4P], [2, 1]), x, x, w, status=st)
print '(A,I0,A,I0,A,ES10.2)', '  cell outside the donor box: status ', st, ', flag ', w%flag(1), ', covered ', w%covered(1)
ok = ok .and. st == OVERLAP_WEIGHTS_OK .and. w%flag(1) == OVERLAP_WEIGHTS_OUTSIDE .and. w%nw == 0 .and. &
     w%covered(1) == 0._R8P
call overlap_weights_2d(blk, reshape([2_I4P, 1_I4P], [2, 1]), x, x, w, status=st)
print '(A,I0)', '  cell out of range: status ', st
ok = ok .and. st == OVERLAP_WEIGHTS_BAD_INPUT .and. w%nr == 0
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check_case(tag, blk, cells, x, inside, ok, ndonors)
   !< overlap_weights_2d on the receptor cells, with the checks of the program notes.
   character(*),                 intent(in)    :: tag        !< Case.
   type(curvilinear_block_2d_t), intent(in)    :: blk        !< Receptor block.
   integer(I4P),                 intent(in)    :: cells(:,:) !< Receptor cells.
   real(R8P),                    intent(in)    :: x(:)       !< Donor grid edges (both directions).
   logical,                      intent(in)    :: inside     !< Every receptor lies in the donor box.
   logical,                      intent(inout) :: ok         !< Checks.
   integer(I4P),       optional, intent(in)    :: ndonors    !< Expected donors of the first receptor.
   real(R8P), parameter :: A0 = 1.3_R8P, B(2) = [0.7_R8P, -2.1_R8P], C0 = 3.7_R8P !< Fields.
   real(R8P)    :: pu(4), pv(4), area, cen(2), diam, ws, wc, cs, cc, qs, qc, mom(2), xc(2), ref, e(6), ar, mi(2) !< Work.
   integer(I4P) :: st, r, k, ncov, nbad                                                    !< Work.

   call overlap_weights_2d(blk, cells, x, x, w, status=st)
   e = 0._R8P ; ncov = 0 ; nbad = 0
   do r = 1, w%nr
      if (w%flag(r) /= OVERLAP_WEIGHTS_OK) then
         nbad = nbad + 1 ; cycle
      endif
      call receptor(blk, cells(:, r), pu, pv)
      if (inside) then
         call polygon(pu, pv, area, cen)
      else
         call clipped(pu, pv, x(size(x)), area, cen)
      endif
      diam = norm2([maxval(pu) - minval(pu), maxval(pv) - minval(pv)])
      ws = 0._R8P ; wc = 0._R8P ; cs = 0._R8P ; cc = 0._R8P ; qs = 0._R8P ; qc = 0._R8P ; mom = 0._R8P
      do k = w%start(r), w%start(r + 1) - 1
         call neumaier(ws, wc, w%weight(k))
         call neumaier(cs, cc, w%weight(k) * C0)
         xc = [0.5_R8P * (x(w%donor(1, k)) + x(w%donor(1, k) + 1)), 0.5_R8P * (x(w%donor(2, k)) + x(w%donor(2, k) + 1))]
         call neumaier(qs, qc, w%weight(k) * (lin(xc) + dot_product(B, w%centroid(:, k) - xc)))
         mom = mom + w%weight(k) * (w%centroid(:, k) - [pu(1), pv(1)])
      enddo
      ws = ws + wc ; cs = cs + cc ; qs = qs + qc
      mi = mom + [pu(1), pv(1)]
      ref = lin(cen)
      e(1) = max(e(1), abs(ws - 1._R8P))
      e(2) = max(e(2), abs(w%receptor_volume(r) - area) / area)
      e(3) = max(e(3), norm2(mi - cen) / diam)
      e(4) = max(e(4), abs(cs - C0) / C0)
      e(5) = max(e(5), abs(qs - ref) / (norm2(B) * diam + abs(ref)))
      if (inside .neqv. w%covered(r) == 1._R8P) ncov = ncov + 1
      if (.not. inside) then
         call polygon(pu, pv, ar, cen)
         e(6) = max(e(6), abs(w%covered(r) - area / ar))
      endif
      if (any(.not. (abs(w%weight(w%start(r):w%start(r + 1) - 1)) <= 1._R8P))) nbad = nbad + 1
   enddo
   print '(3A,I0,A,I0,A,F5.2,A,I0,A,ES10.2,A)', '  ', tag, ': ', w%nr, ' receptors, ', w%nw, ' intersections (', &
         real(w%nw, R8P) / max(1, w%nr), ' each; ', count(w%volume <= 1e-12_R8P * (x(2) - x(1)) ** 2), &
         ' below 1e-12 h^2, the smallest ', minval([huge(1._R8P), w%volume]) / (x(2) - x(1)) ** 2, ' h^2)'
   nbad = nbad + count(.not. (w%volume >= 0._R8P .and. w%weight >= 0._R8P .and. w%weight <= 1._R8P))
   print '(A,6ES10.2)', '    |sum w - 1|, volume, centroid / size, constant, linear, covered: ', e
   ok = ok .and. st == OVERLAP_WEIGHTS_OK .and. nbad == 0 .and. ncov == 0 .and. e(1) <= 2 * EPS .and. &
        e(2) <= 1e-14_R8P .and. e(3) <= 1e-13_R8P .and. e(4) <= 4 * EPS .and. e(5) <= 1e-13_R8P .and. e(6) <= 1e-14_R8P
   if (present(ndonors)) ok = ok .and. w%start(2) - w%start(1) == ndonors
   if (.not. inside) ok = ok .and. all(w%covered < 1._R8P)
   endsubroutine check_case

   pure real(R8P) function lin(p)
   !< The linear field.
   real(R8P), intent(in) :: p(2) !< Point.

   lin = 1.3_R8P + 0.7_R8P * p(1) - 2.1_R8P * p(2)
   endfunction lin

   subroutine receptor(blk, c, pu, pv)
   !< The vertices of cell c, counter-clockwise by the block's handedness.
   type(curvilinear_block_2d_t), intent(in)  :: blk          !< Block.
   integer(I4P),                 intent(in)  :: c(2)         !< Cell.
   real(R8P),                    intent(out) :: pu(4), pv(4) !< Vertices.

   pu = [blk%u(c(1), c(2)), blk%u(c(1) + 1, c(2)), blk%u(c(1) + 1, c(2) + 1), blk%u(c(1), c(2) + 1)]
   pv = [blk%v(c(1), c(2)), blk%v(c(1) + 1, c(2)), blk%v(c(1) + 1, c(2) + 1), blk%v(c(1), c(2) + 1)]
   if (blk%orientation() < 0._R8P) then
      pu = pu([1, 4, 3, 2]) ; pv = pv([1, 4, 3, 2])
   endif
   endsubroutine receptor

   pure subroutine polygon(pu, pv, area, cen)
   !< Area and centroid of a polygon, about its first vertex.
   real(R8P), intent(in)  :: pu(:), pv(:) !< Vertices.
   real(R8P), intent(out) :: area, cen(2) !< Area, centroid.
   real(R8P)              :: a(2), b(2), cr, mom(2) !< Work.
   integer(I4P)           :: k, n         !< Counters.

   n = size(pu, kind=I4P) ; area = 0._R8P ; mom = 0._R8P
   do k = 2, n - 1
      a = [pu(k) - pu(1), pv(k) - pv(1)] ; b = [pu(k + 1) - pu(1), pv(k + 1) - pv(1)]
      cr = a(1) * b(2) - a(2) * b(1)
      area = area + cr ; mom = mom + cr * (a + b)
   enddo
   cen = mom / (3._R8P * area) + [pu(1), pv(1)]
   area = 0.5_R8P * area
   endsubroutine polygon

   pure subroutine clipped(pu, pv, xmax, area, cen)
   !< Area and centroid of the polygon clipped to u <= xmax (one Sutherland-Hodgman pass).
   real(R8P), intent(in)  :: pu(:), pv(:) !< Vertices.
   real(R8P), intent(in)  :: xmax         !< Clip line.
   real(R8P), intent(out) :: area, cen(2) !< Area, centroid.
   real(R8P)              :: qu(2 * size(pu)), qv(2 * size(pu)), t !< Clipped polygon.
   integer(I4P)           :: k, k2, n, nq !< Counters.

   n = size(pu, kind=I4P) ; nq = 0
   do k = 1, n
      k2 = mod(k, n) + 1_I4P
      if (pu(k) <= xmax) then
         nq = nq + 1 ; qu(nq) = pu(k) ; qv(nq) = pv(k)
      endif
      if ((pu(k) <= xmax) .neqv. (pu(k2) <= xmax)) then
         t = (xmax - pu(k)) / (pu(k2) - pu(k))
         nq = nq + 1 ; qu(nq) = xmax ; qv(nq) = pv(k) + t * (pv(k2) - pv(k))
      endif
   enddo
   call polygon(qu(1:nq), qv(1:nq), area, cen)
   endsubroutine clipped

   pure subroutine neumaier(s, comp, a)
   !< s += a, compensated (the sum is s + comp).
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

   subroutine ogrid(m, s, blk)
   !< The cylinder O-grid about (s, s): r in [1, 2], 32 m x 12 m cells (left-handed: i turns counter-clockwise, j out).
   integer(I4P),                 intent(in)    :: m   !< Refinement.
   real(R8P),                    intent(in)    :: s   !< Shift.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P), allocatable :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)           :: ni, nj, i, j   !< Counters.
   real(R8P)              :: r              !< Radius.

   ni = 32 * m + 1 ; nj = 12 * m + 1
   allocate(u(ni, nj), v(ni, nj))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = s + r * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = s + r * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni, :) = u(1, :) ; v(ni, :) = v(1, :)
   call blk%destroy
   call blk%initialize(u, v)
   endsubroutine ogrid

   subroutine quad_block(qu, qv, blk)
   !< A one-cell block of the vertices (i, j) = (1, 1), (2, 1), (2, 2), (1, 2).
   real(R8P),                    intent(in)    :: qu(4), qv(4) !< Vertices.
   type(curvilinear_block_2d_t), intent(inout) :: blk          !< Block.

   call blk%destroy
   call blk%initialize(reshape([qu(1), qu(2), qu(4), qu(3)], [2, 2]), reshape([qv(1), qv(2), qv(4), qv(3)], [2, 2]))
   endsubroutine quad_block

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
   !< An integer as a string.
   integer(I4P), intent(in) :: i !< Integer.
   character(16)            :: s !< String.

   write(s, '(I0)') i
   endfunction itoa
endprogram fossil_test_overlap_weights_2d
