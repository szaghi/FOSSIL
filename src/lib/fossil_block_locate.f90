!< FOSSIL, point location in curvilinear blocks: the cell-box hierarchy and the per-point kernels (issue #71).
module fossil_block_locate
!< FOSSIL, point location in curvilinear blocks: the cell-box hierarchy and the per-point kernels (issue #71).
!<
!< Robust path. The candidates of a point are the cells whose box (closed) contains it, from a bounding volume
!< hierarchy over the cell boxes (median splits on the longest axis of the centroids, leaves of up to 4 cells), all in
!< flat arrays. Each candidate is tested; among the cells that contain the point, the one with the lowest linear index
!< (i fastest, then j, then k) wins, which is the low-side rule: a point on the face between cells i - 1 and i belongs
!< to i - 1. A candidate that is not valid reports only if no valid cell contains the point.
!<
!< - 2-D: a valid or singular cell is convex with straight edges, so containment (closed) is decided exactly by
!<   orient2d against its four edges; xi by the closed-form inverse of the bilinear map (a quadratic in eta), then two
!<   Newton steps, clamped to [0, 1], and exactly 0 or 1 for a point exactly on an edge. An invalid cell contains the point if the inverse has a root in [0, 1]^2 (warped).
!< - 3-D: Newton on the trilinear map from the cell centre (at most NEWTON_MAX steps, xi kept in [-1, 2]), converged
!<   when the residual is below 1e-14 of the cell size plus the rounding of the coordinates; inside when xi is in
!<   [0, 1]^3 within XI_TOL.
!<
!< The kernels are free `_oac` procedures (`!$acc routine seq`) on plain arrays, so the CPU and the device run the same
!< code (issue #20 rules): no type-bound calls, no array-valued functions.
use penf,                  only : I4P, R8P
use fossil_predicates,     only : orient2d_xy
use fossil_block_geometry, only : CELL_VALID, CELL_SINGULAR, CELL_INVALID, CELL_UNDETERMINED

implicit none
private
public :: cell_bvh_t, build_cell_bvh, locate_point_3d_oac, locate_point_2d_oac, newton_3d_oac, inverse_2d_oac
public :: LOC_INSIDE, LOC_OUTSIDE, LOC_WARPED, LOC_DEGENERATE, LOC_NOT_CONVERGED, XI_TOL, NEWTON_MAX, LOC_STACK

integer(I4P), parameter :: LOC_INSIDE        = 1_I4P !< A valid cell contains the point.
integer(I4P), parameter :: LOC_OUTSIDE       = 2_I4P !< No cell contains it.
integer(I4P), parameter :: LOC_WARPED        = 3_I4P !< Only an invalid cell contains it.
integer(I4P), parameter :: LOC_DEGENERATE    = 4_I4P !< Only an undetermined cell may contain it.
integer(I4P), parameter :: LOC_NOT_CONVERGED = 5_I4P !< The 3-D inverse did not converge in a candidate.
real(R8P),    parameter :: XI_TOL            = 1e-12_R8P !< 3-D containment tolerance on xi.
integer(I4P), parameter :: NEWTON_MAX        = 40_I4P  !< Newton steps at most.
integer(I4P), parameter :: LOC_STACK         = 64_I4P  !< Traversal stack.
integer(I4P), parameter :: LEAF              = 4_I4P   !< Cells per leaf at most.

type :: cell_bvh_t
   !< Hierarchy over the cell boxes of one block (d = 2 or 3 rows).
   integer(I4P)              :: nnode = 0_I4P          !< Nodes.
   integer(I4P)              :: ncell = 0_I4P          !< Cells.
   real(R8P),    allocatable :: node_lo(:,:), node_hi(:,:) !< Node boxes [d, nnode].
   integer(I4P), allocatable :: left(:), right(:)      !< Children (0 for a leaf) [nnode].
   integer(I4P), allocatable :: first(:), last(:)      !< Leaf ranges in `order` [nnode].
   integer(I4P), allocatable :: order(:)               !< Cells in leaf order [ncell].
   real(R8P),    allocatable :: cell_lo(:,:), cell_hi(:,:) !< Cell boxes [d, ncell].
endtype cell_bvh_t

contains
   subroutine build_cell_bvh(cell_lo, cell_hi, bvh)
   !< Build the hierarchy over the given cell boxes [d, ncell].
   real(R8P),        intent(in)  :: cell_lo(:,:), cell_hi(:,:) !< Cell boxes.
   type(cell_bvh_t), intent(out) :: bvh                       !< Hierarchy.
   integer(I4P)                  :: n, i                      !< Counters.

   n = size(cell_lo, 2, I4P)
   bvh%ncell = n ; bvh%cell_lo = cell_lo ; bvh%cell_hi = cell_hi
   allocate(bvh%order(n)) ; bvh%order = [(i, i = 1, n)]
   allocate(bvh%node_lo(size(cell_lo, 1), max(2 * n, 1)), bvh%node_hi(size(cell_lo, 1), max(2 * n, 1)))
   allocate(bvh%left(max(2 * n, 1)), bvh%right(max(2 * n, 1)), bvh%first(max(2 * n, 1)), bvh%last(max(2 * n, 1)))
   bvh%nnode = 0_I4P
   if (n > 0) call build_node(bvh, 1_I4P, n, i)
   endsubroutine build_cell_bvh

   recursive subroutine build_node(bvh, lo, hi, node)
   !< Node over order(lo:hi): box, then a median split on the longest axis of the cell centres.
   type(cell_bvh_t), intent(inout) :: bvh     !< Hierarchy.
   integer(I4P),     intent(in)    :: lo, hi  !< Range in `order`.
   integer(I4P),     intent(out)   :: node    !< This node.
   real(R8P), allocatable          :: key(:)  !< Centres along the axis.
   real(R8P)                       :: clo(size(bvh%cell_lo, 1)), chi(size(bvh%cell_lo, 1)) !< Centre box.
   integer(I4P)                    :: k, ax, mid, l, r !< Work.

   bvh%nnode = bvh%nnode + 1_I4P ; node = bvh%nnode
   bvh%node_lo(:, node) = minval(bvh%cell_lo(:, bvh%order(lo:hi)), dim=2)
   bvh%node_hi(:, node) = maxval(bvh%cell_hi(:, bvh%order(lo:hi)), dim=2)
   bvh%first(node) = lo ; bvh%last(node) = hi ; bvh%left(node) = 0_I4P ; bvh%right(node) = 0_I4P
   if (hi - lo + 1_I4P <= LEAF) return
   clo = huge(1._R8P) ; chi = -huge(1._R8P)
   do k = lo, hi
      clo = min(clo, 0.5_R8P * (bvh%cell_lo(:, bvh%order(k)) + bvh%cell_hi(:, bvh%order(k))))
      chi = max(chi, 0.5_R8P * (bvh%cell_lo(:, bvh%order(k)) + bvh%cell_hi(:, bvh%order(k))))
   enddo
   ax = maxloc(chi - clo, dim=1)
   allocate(key(lo:hi))
   do k = lo, hi
      key(k) = bvh%cell_lo(ax, bvh%order(k)) + bvh%cell_hi(ax, bvh%order(k))
   enddo
   mid = (lo + hi) / 2_I4P
   call select_kth(key, bvh%order, lo, hi, mid)
   call build_node(bvh, lo, mid, l)
   call build_node(bvh, mid + 1_I4P, hi, r)
   bvh%left(node) = l ; bvh%right(node) = r
   endsubroutine build_node

   subroutine select_kth(key, order, lo, hi, kth)
   !< Partial sort of order(lo:hi) by key (ties by cell index): the kth in place, smaller before, larger after.
   integer(I4P), intent(in)    :: lo, hi, kth !< Range, rank.
   real(R8P),    intent(inout) :: key(lo:)   !< Keys.
   integer(I4P), intent(inout) :: order(:)   !< Cells.
   integer(I4P)                :: l, r, i, j, ti !< Work.
   real(R8P)                   :: pk, tk     !< Pivot, swap.
   integer(I4P)                :: pc         !< Pivot cell.

   l = lo ; r = hi
   do while (l < r)
      pk = key((l + r) / 2) ; pc = order((l + r) / 2)
      i = l ; j = r
      do while (i <= j)
         do while (key(i) < pk .or. (key(i) == pk .and. order(i) < pc))
            i = i + 1
         enddo
         do while (key(j) > pk .or. (key(j) == pk .and. order(j) > pc))
            j = j - 1
         enddo
         if (i <= j) then
            tk = key(i) ; key(i) = key(j) ; key(j) = tk
            ti = order(i) ; order(i) = order(j) ; order(j) = ti
            i = i + 1 ; j = j - 1
         endif
      enddo
      if (kth <= j) then
         r = j
      elseif (kth >= i) then
         l = i
      else
         exit
      endif
   enddo
   endsubroutine select_kth

   pure subroutine locate_point_3d_oac(p, ni, nj, x, y, z, cstat, nnode, node_lo, node_hi, left, right, first, last, &
                                       order, cell_lo, cell_hi, cell, xi, status)
   !$acc routine seq
   !< Robust location of p in a 3-D block (module notes). cell = linear index (0 if outside).
   real(R8P),    intent(in)  :: p(3)                    !< Point.
   integer(I4P), intent(in)  :: ni, nj                  !< Vertices in i, j.
   real(R8P),    intent(in)  :: x(*), y(*), z(*)        !< Vertices (flattened, i fastest).
   integer(I4P), intent(in)  :: cstat(*)                !< Cell validity.
   integer(I4P), intent(in)  :: nnode                   !< Nodes.
   real(R8P),    intent(in)  :: node_lo(3, *), node_hi(3, *) !< Node boxes.
   integer(I4P), intent(in)  :: left(*), right(*), first(*), last(*), order(*) !< Tree.
   real(R8P),    intent(in)  :: cell_lo(3, *), cell_hi(3, *) !< Cell boxes.
   integer(I4P), intent(out) :: cell                    !< Cell (linear index).
   real(R8P),    intent(out) :: xi(3)                   !< Local coordinates.
   integer(I4P), intent(out) :: status                  !< LOC_*.
   integer(I4P)              :: stack(LOC_STACK), top, nd, k, c, st, alt, alt_st !< Traversal.
   real(R8P)                 :: cc(3, 8), q(3), alt_xi(3) !< Candidate.
   logical                   :: conv                    !< Newton converged.

   cell = 0_I4P ; xi = 0._R8P ; status = LOC_OUTSIDE ; alt = 0_I4P ; alt_st = LOC_OUTSIDE ; alt_xi = 0._R8P
   if (nnode < 1) return
   top = 1_I4P ; stack(1) = 1_I4P
   do while (top > 0)
      nd = stack(top) ; top = top - 1_I4P
      if (any(p < node_lo(:, nd)) .or. any(p > node_hi(:, nd))) cycle
      if (left(nd) /= 0) then
         top = top + 1_I4P ; stack(top) = right(nd) ; top = top + 1_I4P ; stack(top) = left(nd)
         cycle
      endif
      do k = first(nd), last(nd)
         c = order(k)
         if (cell /= 0 .and. c > cell) cycle
         if (any(p < cell_lo(:, c)) .or. any(p > cell_hi(:, c))) cycle
         call corners_flat_3d(c, ni, nj, x, y, z, cc)
         call newton_3d_oac(cc, p, q, conv)
         st = cstat(c)
         if (conv .and. all(q >= -XI_TOL) .and. all(q <= 1._R8P + XI_TOL)) then
            if (st == CELL_VALID .or. st == CELL_SINGULAR) then
               cell = c ; xi = q
            elseif (alt == 0 .or. alt_st == LOC_NOT_CONVERGED .or. c < alt) then
               alt = c ; alt_xi = q ; alt_st = merge(LOC_WARPED, LOC_DEGENERATE, st == CELL_INVALID)
            endif
         elseif (.not. conv) then
            if (alt == 0 .or. (alt_st == LOC_NOT_CONVERGED .and. c < alt)) then
               alt = c ; alt_xi = q ; alt_st = LOC_NOT_CONVERGED
            endif
         endif
      enddo
   enddo
   if (cell /= 0) then
      status = LOC_INSIDE
   elseif (alt /= 0) then
      cell = alt ; xi = alt_xi ; status = alt_st
   endif
   endsubroutine locate_point_3d_oac

   pure subroutine locate_point_2d_oac(p, ni, u, v, cstat, orient, nnode, node_lo, node_hi, left, right, first, last, &
                                       order, cell_lo, cell_hi, cell, xi, status)
   !$acc routine seq
   !< Robust location of p in a 2-D block (module notes). cell = linear index (0 if outside).
   real(R8P),    intent(in)  :: p(2)                    !< Point.
   integer(I4P), intent(in)  :: ni                      !< Vertices in i.
   real(R8P),    intent(in)  :: u(*), v(*)              !< Vertices (flattened, i fastest).
   integer(I4P), intent(in)  :: cstat(*)                !< Cell validity.
   real(R8P),    intent(in)  :: orient                  !< Block orientation.
   integer(I4P), intent(in)  :: nnode                   !< Nodes.
   real(R8P),    intent(in)  :: node_lo(2, *), node_hi(2, *) !< Node boxes.
   integer(I4P), intent(in)  :: left(*), right(*), first(*), last(*), order(*) !< Tree.
   real(R8P),    intent(in)  :: cell_lo(2, *), cell_hi(2, *) !< Cell boxes.
   integer(I4P), intent(out) :: cell                    !< Cell (linear index).
   real(R8P),    intent(out) :: xi(2)                   !< Local coordinates.
   integer(I4P), intent(out) :: status                  !< LOC_*.
   integer(I4P)              :: stack(LOC_STACK), top, nd, k, c, alt !< Traversal.
   real(R8P)                 :: cc(2, 4), q(2), alt_xi(2), o(4) !< Candidate, edge orientations.
   logical                   :: found                   !< Root in the unit square.

   cell = 0_I4P ; xi = 0._R8P ; status = LOC_OUTSIDE ; alt = 0_I4P ; alt_xi = 0._R8P
   if (nnode < 1) return
   top = 1_I4P ; stack(1) = 1_I4P
   do while (top > 0)
      nd = stack(top) ; top = top - 1_I4P
      if (any(p < node_lo(:, nd)) .or. any(p > node_hi(:, nd))) cycle
      if (left(nd) /= 0) then
         top = top + 1_I4P ; stack(top) = right(nd) ; top = top + 1_I4P ; stack(top) = left(nd)
         cycle
      endif
      do k = first(nd), last(nd)
         c = order(k)
         if (cell /= 0 .and. c > cell) cycle
         if (any(p < cell_lo(:, c)) .or. any(p > cell_hi(:, c))) cycle
         call corners_flat_2d(c, ni, u, v, cc)
         if (cstat(c) == CELL_VALID .or. cstat(c) == CELL_SINGULAR) then
            o(1) = orient * orient2d_xy(cc(1, 1), cc(2, 1), cc(1, 2), cc(2, 2), p(1), p(2)) ! eta = 0
            if (o(1) < 0._R8P) cycle
            o(2) = orient * orient2d_xy(cc(1, 2), cc(2, 2), cc(1, 4), cc(2, 4), p(1), p(2)) ! xi = 1
            if (o(2) < 0._R8P) cycle
            o(3) = orient * orient2d_xy(cc(1, 4), cc(2, 4), cc(1, 3), cc(2, 3), p(1), p(2)) ! eta = 1
            if (o(3) < 0._R8P) cycle
            o(4) = orient * orient2d_xy(cc(1, 3), cc(2, 3), cc(1, 1), cc(2, 1), p(1), p(2)) ! xi = 0
            if (o(4) < 0._R8P) cycle
            call inverse_2d_oac(cc, p, q, found)
            cell = c ; xi = min(max(q, 0._R8P), 1._R8P)
            ! a point exactly on an edge (of positive length) gets that edge's coordinate exactly
            if (o(1) == 0._R8P .and. any(cc(:, 1) /= cc(:, 2))) xi(2) = 0._R8P
            if (o(3) == 0._R8P .and. any(cc(:, 3) /= cc(:, 4))) xi(2) = 1._R8P
            if (o(4) == 0._R8P .and. any(cc(:, 1) /= cc(:, 3))) xi(1) = 0._R8P
            if (o(2) == 0._R8P .and. any(cc(:, 2) /= cc(:, 4))) xi(1) = 1._R8P
         elseif (alt == 0 .or. c < alt) then
            call inverse_2d_oac(cc, p, q, found)
            if (found) then
               alt = c ; alt_xi = q
            endif
         endif
      enddo
   enddo
   if (cell /= 0) then
      status = LOC_INSIDE
   elseif (alt /= 0) then
      cell = alt ; xi = alt_xi ; status = LOC_WARPED
   endif
   endsubroutine locate_point_2d_oac

   pure subroutine corners_flat_3d(c, ni, nj, x, y, z, cc)
   !$acc routine seq
   !< Corners of the cell with linear index c (vertices flattened, i fastest).
   integer(I4P), intent(in)  :: c, ni, nj        !< Cell, vertices in i, j.
   real(R8P),    intent(in)  :: x(*), y(*), z(*) !< Vertices.
   real(R8P),    intent(out) :: cc(3, 8)         !< Corners.
   integer(I4P)              :: i, j, k, a, b, e, w !< Indices.

   i = mod(c - 1_I4P, ni - 1_I4P) + 1_I4P
   j = mod((c - 1_I4P) / (ni - 1_I4P), nj - 1_I4P) + 1_I4P
   k = (c - 1_I4P) / ((ni - 1_I4P) * (nj - 1_I4P)) + 1_I4P
   do e = 0, 1
      do b = 0, 1
         do a = 0, 1
            w = (i + a) + ni * ((j + b - 1_I4P) + nj * (k + e - 1_I4P))
            cc(1, 1 + a + 2 * b + 4 * e) = x(w) ; cc(2, 1 + a + 2 * b + 4 * e) = y(w) ; cc(3, 1 + a + 2 * b + 4 * e) = z(w)
         enddo
      enddo
   enddo
   endsubroutine corners_flat_3d

   pure subroutine corners_flat_2d(c, ni, u, v, cc)
   !$acc routine seq
   !< Corners of the cell with linear index c (vertices flattened, i fastest).
   integer(I4P), intent(in)  :: c, ni      !< Cell, vertices in i.
   real(R8P),    intent(in)  :: u(*), v(*) !< Vertices.
   real(R8P),    intent(out) :: cc(2, 4)   !< Corners.
   integer(I4P)              :: i, j, w    !< Indices.

   i = mod(c - 1_I4P, ni - 1_I4P) + 1_I4P ; j = (c - 1_I4P) / (ni - 1_I4P) + 1_I4P
   w = i + ni * (j - 1_I4P)
   cc(1, 1) = u(w) ; cc(2, 1) = v(w) ; cc(1, 2) = u(w + 1) ; cc(2, 2) = v(w + 1)
   cc(1, 3) = u(w + ni) ; cc(2, 3) = v(w + ni) ; cc(1, 4) = u(w + ni + 1) ; cc(2, 4) = v(w + ni + 1)
   endsubroutine corners_flat_2d

   pure subroutine newton_3d_oac(cc, p, xi, conv)
   !$acc routine seq
   !< Newton on the trilinear map from the cell centre (module notes).
   real(R8P), intent(in)  :: cc(3, 8) !< Corners.
   real(R8P), intent(in)  :: p(3)     !< Point.
   real(R8P), intent(out) :: xi(3)    !< Local coordinates.
   logical,   intent(out) :: conv     !< Converged.
   real(R8P)              :: xp(3), jac(3, 3), r(3), d(3), det, tol, scl !< Work.
   integer(I4P)           :: it, a   !< Counters.

   xi = 0.5_R8P ; conv = .false.
   scl = 0._R8P ; tol = 0._R8P
   do a = 1, 3
      scl = max(scl, maxval(cc(a, :)) - minval(cc(a, :)))
      tol = max(tol, maxval(abs(cc(a, :))))
   enddo
   tol = 1e-14_R8P * scl + 8._R8P * epsilon(1._R8P) * max(tol, maxval(abs(p)))
   do it = 1, NEWTON_MAX
      call trilinear_eval(cc, xi, xp, jac)
      r = xp - p
      if (maxval(abs(r)) <= tol) then
         conv = .true. ; return
      endif
      det = jac(1, 1) * (jac(2, 2) * jac(3, 3) - jac(3, 2) * jac(2, 3)) &
          + jac(2, 1) * (jac(3, 2) * jac(1, 3) - jac(1, 2) * jac(3, 3)) &
          + jac(3, 1) * (jac(1, 2) * jac(2, 3) - jac(2, 2) * jac(1, 3))
      if (det == 0._R8P) return
      ! d = J^-1 r by Cramer: d_a = det(J with column a replaced by r) / det
      d(1) = (r(1) * (jac(2, 2) * jac(3, 3) - jac(3, 2) * jac(2, 3)) + r(2) * (jac(3, 2) * jac(1, 3) - jac(1, 2) * jac(3, 3)) &
             + r(3) * (jac(1, 2) * jac(2, 3) - jac(2, 2) * jac(1, 3))) / det
      d(2) = (jac(1, 1) * (r(2) * jac(3, 3) - r(3) * jac(2, 3)) + jac(2, 1) * (r(3) * jac(1, 3) - r(1) * jac(3, 3)) &
             + jac(3, 1) * (r(1) * jac(2, 3) - r(2) * jac(1, 3))) / det
      d(3) = (jac(1, 1) * (jac(2, 2) * r(3) - jac(3, 2) * r(2)) + jac(2, 1) * (jac(3, 2) * r(1) - jac(1, 2) * r(3)) &
             + jac(3, 1) * (jac(1, 2) * r(2) - jac(2, 2) * r(1))) / det
      xi = min(max(xi - d, -1._R8P), 2._R8P)
   enddo
   call trilinear_eval(cc, xi, xp, jac)
   conv = maxval(abs(xp - p)) <= tol
   endsubroutine newton_3d_oac

   pure subroutine trilinear_eval(cc, xi, x, jac)
   !$acc routine seq
   !< Point and Jacobian of the trilinear map.
   real(R8P), intent(in)  :: cc(3, 8)  !< Corners.
   real(R8P), intent(in)  :: xi(3)     !< Local coordinates.
   real(R8P), intent(out) :: x(3), jac(3, 3) !< Point, Jacobian (columns: derivatives).
   real(R8P)              :: s0, s1, t0, t1, u0, u1 !< 1 - xi, xi.
   integer(I4P)           :: a        !< Component.

   s1 = xi(1) ; s0 = 1._R8P - s1 ; t1 = xi(2) ; t0 = 1._R8P - t1 ; u1 = xi(3) ; u0 = 1._R8P - u1
   do a = 1, 3
      x(a) = ((cc(a, 1) * s0 + cc(a, 2) * s1) * t0 + (cc(a, 3) * s0 + cc(a, 4) * s1) * t1) * u0 + &
             ((cc(a, 5) * s0 + cc(a, 6) * s1) * t0 + (cc(a, 7) * s0 + cc(a, 8) * s1) * t1) * u1
      jac(a, 1) = ((cc(a, 2) - cc(a, 1)) * t0 + (cc(a, 4) - cc(a, 3)) * t1) * u0 + &
                  ((cc(a, 6) - cc(a, 5)) * t0 + (cc(a, 8) - cc(a, 7)) * t1) * u1
      jac(a, 2) = ((cc(a, 3) - cc(a, 1)) * s0 + (cc(a, 4) - cc(a, 2)) * s1) * u0 + &
                  ((cc(a, 7) - cc(a, 5)) * s0 + (cc(a, 8) - cc(a, 6)) * s1) * u1
      jac(a, 3) = ((cc(a, 5) - cc(a, 1)) * s0 + (cc(a, 6) - cc(a, 2)) * s1) * t0 + &
                  ((cc(a, 7) - cc(a, 3)) * s0 + (cc(a, 8) - cc(a, 4)) * s1) * t1
   enddo
   endsubroutine trilinear_eval

   pure subroutine inverse_2d_oac(cc, p, xi, found)
   !$acc routine seq
   !< Inverse of the bilinear map: eta from (c x d) eta^2 + (c x b - e x d) eta - e x b = 0 (e = p - c1, b = c2 - c1,
   !< c = c3 - c1, d = c4 - c3 - c2 + c1), the root nearest to [0, 1]; xi = (e - c eta) . g / |g|^2, g = b + d eta;
   !< then two Newton steps. found: the result is in [0, 1]^2 within 1e-12.
   real(R8P), intent(in)  :: cc(2, 4) !< Corners.
   real(R8P), intent(in)  :: p(2)     !< Point.
   real(R8P), intent(out) :: xi(2)    !< Local coordinates.
   logical,   intent(out) :: found    !< Inside the unit square.
   real(R8P)              :: e(2), b(2), c(2), d(2), g(2), qa, qb, qc, disc, qq, r1, r2, x(2), j11, j12, j21, j22, det !< Work.
   integer(I4P)           :: it       !< Newton.

   e = p - cc(:, 1) ; b = cc(:, 2) - cc(:, 1) ; c = cc(:, 3) - cc(:, 1) ; d = cc(:, 4) - cc(:, 3) - cc(:, 2) + cc(:, 1)
   qa = c(1) * d(2) - c(2) * d(1)
   qb = (c(1) * b(2) - c(2) * b(1)) - (e(1) * d(2) - e(2) * d(1))
   qc = -(e(1) * b(2) - e(2) * b(1))
   if (abs(qa) <= 1e-14_R8P * (abs(qb) + abs(qc))) then
      r1 = 0.5_R8P
      if (qb /= 0._R8P) r1 = -qc / qb
   else
      disc = max(qb * qb - 4._R8P * qa * qc, 0._R8P)
      qq = -0.5_R8P * (qb + sign(sqrt(disc), qb))
      r1 = qq / qa ; r2 = r1
      if (qq /= 0._R8P) r2 = qc / qq
      if (dist01(r2) < dist01(r1)) r1 = r2
   endif
   xi(2) = r1
   g = b + d * r1
   xi(1) = 0.5_R8P
   if (g(1) * g(1) + g(2) * g(2) > 0._R8P) xi(1) = ((e(1) - c(1) * r1) * g(1) + (e(2) - c(2) * r1) * g(2)) / &
                                                   (g(1) * g(1) + g(2) * g(2))
   do it = 1, 2
      x = cc(:, 1) * (1._R8P - xi(1)) * (1._R8P - xi(2)) + cc(:, 2) * xi(1) * (1._R8P - xi(2)) + &
          cc(:, 3) * (1._R8P - xi(1)) * xi(2) + cc(:, 4) * xi(1) * xi(2) - p
      j11 = b(1) + d(1) * xi(2) ; j21 = b(2) + d(2) * xi(2) ; j12 = c(1) + d(1) * xi(1) ; j22 = c(2) + d(2) * xi(1)
      det = j11 * j22 - j12 * j21
      if (det == 0._R8P) exit
      xi = xi - [(j22 * x(1) - j12 * x(2)) / det, (-j21 * x(1) + j11 * x(2)) / det]
   enddo
   found = all(xi >= -1e-12_R8P) .and. all(xi <= 1._R8P + 1e-12_R8P)
   endsubroutine inverse_2d_oac

   pure real(R8P) function dist01(t) result(dd)
   !$acc routine seq
   !< Distance of t from [0, 1].
   real(R8P), intent(in) :: t !< Value.

   dd = max(-t, t - 1._R8P, 0._R8P)
   endfunction dist01
endmodule fossil_block_locate
