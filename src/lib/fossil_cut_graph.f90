!< FOSSIL, the graph of the cut-cell pieces and what a solver builds on it (issue #70).
module fossil_cut_graph
!< FOSSIL, the graph of the cut-cell pieces and what a solver builds on it (issue #70).
!<
!< Links. A link (p, q, f, A) says that piece p and piece q of the cell across face f of p's cell share a fluid region
!< of that face, of area A (a length in 2-D). Faces are numbered -x, +x, -y, +y, -z, +z (1-6; -u, +u, -v, +v in 2-D).
!< Links are stored from each side (CSR over the pieces, sorted by face, then piece), with the same area both ways.
!<
!< 2-D: the fluid intervals of a cell edge are bounded by the crossings of the walk and the cell corners, the same
!< points in the two cells that share the edge (`fossil_cut_cells_2d`), so the intervals of two pieces across an edge
!< are equal or disjoint: a link per pair of equal intervals (overlap, in case), of their length. A zero-length
!< interval gives no link.
use fossil_cut_cells_2d, only : cut_cells_2d_t
use penf,                only : I4P, R8P

implicit none
private
public :: cut_links_t, cut_links_2d, build_csr

type :: cut_links_t
   !< Links between the pieces of a block (see the module notes).
   integer(I4P)              :: np = 0_I4P     !< Pieces.
   integer(I4P), allocatable :: start(:)       !< Links of piece p: start(p) .. start(p + 1) - 1 [np + 1].
   integer(I4P), allocatable :: piece(:)       !< Piece across.
   integer(I4P), allocatable :: face(:)        !< Face of p's cell (1-6; 1-4 in 2-D).
   real(R8P),    allocatable :: area(:)        !< Shared fluid area (length in 2-D).
   integer(I4P)              :: unmatched = 0_I4P !< Fluid regions of a face that could not be matched to a piece.
endtype cut_links_t

contains
   subroutine cut_links_2d(cut, links)
   !< Links of the pieces of the 2-D cut cells `cut` (see the module notes).
   type(cut_cells_2d_t), intent(in)  :: cut   !< Cut cells.
   type(cut_links_t),    intent(out) :: links !< Links.
   integer(I4P), allocatable :: owner(:), ord(:), lp(:), lq(:), lf(:) !< Interval piece, order by edge; raw links.
   real(R8P),    allocatable :: la(:)            !< Raw link lengths.
   integer(I4P)              :: n, k, kk, m, p, q, nl, cap !< Work.
   real(R8P)                 :: ov                !< Overlap.

   links%np = cut%np
   n = size(cut%piece_iv_edge, kind=I4P)
   allocate(owner(n), ord(n))
   do p = 1, cut%np
      owner(cut%piece_iv_start(p):cut%piece_iv_start(p + 1) - 1) = p
   enddo
   call sort_by_key(cut%piece_iv_edge, cut%piece_iv_lo, ord)
   cap = max(16_I4P, n) ; allocate(lp(cap), lq(cap), lf(cap), la(cap)) ; nl = 0_I4P
   k = 1_I4P
   do while (k <= n)
      m = k ! intervals of one edge: ord(k:m)
      do while (m < n)
         if (cut%piece_iv_edge(ord(m + 1)) /= cut%piece_iv_edge(ord(k))) exit
         m = m + 1_I4P
      enddo
      do kk = k, m
         do q = kk + 1, m
            if (owner(ord(kk)) == owner(ord(q))) cycle
            ov = min(cut%piece_iv_hi(ord(kk)), cut%piece_iv_hi(ord(q))) - max(cut%piece_iv_lo(ord(kk)), cut%piece_iv_lo(ord(q)))
            if (ov <= 0._R8P) cycle
            call add_pair(owner(ord(kk)), owner(ord(q)), cut%piece_iv_edge(ord(kk)), ov)
         enddo
      enddo
      k = m + 1_I4P
   enddo
   call build_csr(nl, lp, lq, lf, la, links)
   contains
      subroutine add_pair(pa, pb, ed, len)
      !< The two directed links of pieces pa, pb sharing length len of edge ed.
      integer(I4P), intent(in) :: pa, pb, ed !< Pieces, edge id.
      real(R8P),    intent(in) :: len        !< Length.

      if (nl + 2_I4P > size(lp)) then
         lp = [lp, lp] ; lq = [lq, lq] ; lf = [lf, lf] ; la = [la, la]
      endif
      nl = nl + 1_I4P ; lp(nl) = pa ; lq(nl) = pb ; lf(nl) = face_2d(pa, ed) ; la(nl) = len
      nl = nl + 1_I4P ; lp(nl) = pb ; lq(nl) = pa ; lf(nl) = face_2d(pb, ed) ; la(nl) = len
      endsubroutine add_pair

      integer(I4P) function face_2d(pp, ed) result(fc)
      !< Face (1 -u, 2 +u, 3 -v, 4 +v) of piece pp's cell that is edge ed.
      integer(I4P), intent(in) :: pp, ed !< Piece, edge id.
      integer(I4P)             :: i, j, ll !< Work.

      i = cut%piece_cell(1, pp) ; j = cut%piece_cell(2, pp)
      if (ed <= (cut%nx + 1_I4P) * cut%ny) then ! vertical edge on line ll
         ll = mod(ed - 1_I4P, cut%nx + 1_I4P) + 1_I4P
         fc = merge(1_I4P, 2_I4P, ll == i)
      else
         ll = (ed - (cut%nx + 1_I4P) * cut%ny - 1_I4P) / cut%nx + 1_I4P
         fc = merge(3_I4P, 4_I4P, ll == j)
      endif
      endfunction face_2d
   endsubroutine cut_links_2d

   subroutine build_csr(nl, lp, lq, lf, la, links)
   !< CSR of the directed links (lp -> lq, face lf, area la), merged per (p, face, q) and sorted by face, then q.
   integer(I4P),      intent(in)    :: nl                 !< Links.
   integer(I4P),      intent(in)    :: lp(:), lq(:), lf(:) !< From, to, face.
   real(R8P),         intent(in)    :: la(:)              !< Areas.
   type(cut_links_t), intent(inout) :: links              !< Links (np set).
   integer(I4P), allocatable :: ord(:)          !< Order.
   integer(I4P)              :: k, m, n, p       !< Work.

   allocate(ord(nl))
   ord = [(k, k = 1, nl)] ! by piece, then face, then piece across
   call merge_sort_links(lp, lf, lq, ord)
   allocate(links%start(links%np + 1)) ; links%start = 0_I4P
   allocate(links%piece(nl), links%face(nl), links%area(nl))
   n = 0_I4P ; k = 1_I4P
   do while (k <= nl)
      m = k
      do while (m < nl)
         if (lp(ord(m + 1)) /= lp(ord(k)) .or. lf(ord(m + 1)) /= lf(ord(k)) .or. lq(ord(m + 1)) /= lq(ord(k))) exit
         m = m + 1_I4P
      enddo
      n = n + 1_I4P
      links%piece(n) = lq(ord(k)) ; links%face(n) = lf(ord(k)) ; links%area(n) = sum(la(ord(k:m)))
      links%start(lp(ord(k)) + 1) = links%start(lp(ord(k)) + 1) + 1_I4P
      k = m + 1_I4P
   enddo
   links%start(1) = 1_I4P
   do p = 1, links%np
      links%start(p + 1) = links%start(p + 1) + links%start(p)
   enddo
   links%piece = links%piece(1:n) ; links%face = links%face(1:n) ; links%area = links%area(1:n)
   endsubroutine build_csr

   subroutine merge_sort_links(k1, k2, k3, ord)
   !< Sort ord by (k1, k2, k3) of its entries (bottom-up merge sort, stable).
   integer(I4P), intent(in)    :: k1(:), k2(:), k3(:) !< Keys.
   integer(I4P), intent(inout) :: ord(:)              !< Order.
   integer(I4P), allocatable   :: tmp(:)              !< Work.
   integer(I4P)                :: n, w, lo, mid, hi, a, b, t !< Work.

   n = size(ord, kind=I4P) ; allocate(tmp(n)) ; w = 1_I4P
   do while (w < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + w - 1_I4P, n) ; hi = min(lo + 2_I4P * w - 1_I4P, n)
         a = lo ; b = mid + 1_I4P ; t = lo
         do while (t <= hi)
            if (b > hi) then
               tmp(t) = ord(a) ; a = a + 1_I4P
            elseif (a > mid) then
               tmp(t) = ord(b) ; b = b + 1_I4P
            elseif (less(ord(b), ord(a))) then
               tmp(t) = ord(b) ; b = b + 1_I4P
            else
               tmp(t) = ord(a) ; a = a + 1_I4P
            endif
            t = t + 1_I4P
         enddo
         lo = lo + 2_I4P * w
      enddo
      ord = tmp ; w = 2_I4P * w
   enddo
   contains
      logical function less(i, j)
      !< Entry i before entry j.
      integer(I4P), intent(in) :: i, j !< Entries.

      less = k1(i) < k1(j) .or. (k1(i) == k1(j) .and. (k2(i) < k2(j) .or. (k2(i) == k2(j) .and. k3(i) < k3(j))))
      endfunction less
   endsubroutine merge_sort_links

   subroutine sort_by_key(ik, rk, ord)
   !< Order by (ik, rk), stable (bottom-up merge sort).
   integer(I4P), intent(in)  :: ik(:)  !< Integer key.
   real(R8P),    intent(in)  :: rk(:)  !< Real key.
   integer(I4P), intent(out) :: ord(:) !< Order.
   integer(I4P), allocatable :: tmp(:) !< Work.
   integer(I4P)              :: n, w, lo, mid, hi, a, b, t !< Work.

   n = size(ord, kind=I4P) ; ord = [(a, a = 1, n)] ; allocate(tmp(n)) ; w = 1_I4P
   do while (w < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + w - 1_I4P, n) ; hi = min(lo + 2_I4P * w - 1_I4P, n)
         a = lo ; b = mid + 1_I4P ; t = lo
         do while (t <= hi)
            if (b > hi) then
               tmp(t) = ord(a) ; a = a + 1_I4P
            elseif (a > mid) then
               tmp(t) = ord(b) ; b = b + 1_I4P
            elseif (ik(ord(b)) < ik(ord(a)) .or. (ik(ord(b)) == ik(ord(a)) .and. rk(ord(b)) < rk(ord(a)))) then
               tmp(t) = ord(b) ; b = b + 1_I4P
            else
               tmp(t) = ord(a) ; a = a + 1_I4P
            endif
            t = t + 1_I4P
         enddo
         lo = lo + 2_I4P * w
      enddo
      ord = tmp ; w = 2_I4P * w
   enddo
   endsubroutine sort_by_key
endmodule fossil_cut_graph
