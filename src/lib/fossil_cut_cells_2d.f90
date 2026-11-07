!< FOSSIL, exact 2-D cut cells of a polyline curve on a rectilinear block (issue #40 step 3).
module fossil_cut_cells_2d
!< FOSSIL, exact 2-D cut cells of a polyline curve on a rectilinear block (issue #40 step 3).
!<
!< The solid is on the left of every segment (outer loops CCW, holes CW); the fluid is outside it.
!< For every cell of the block given by its edges x(nx+1), y(ny+1) the module computes the fluid
!< area fraction, the fluid fraction (aperture) of each of the four edges, the fluid centroid, the
!< wetted length, centroid and integral of the outward (solid to fluid) normal, and the connected
!< fluid pieces of the cell with their own measures.
!<
!< Exact topology (decision 1 of #40). Every decision is a comparison or the exact `orient2d_xy`:
!< - the vertical grid line u = x_i is taken at x_i + eps and the horizontal one v = y_j at
!<   y_j + eps^2 (the rule of `fossil_line_query`: a point exactly on a line is on its low side),
!<   so every point lies in exactly one half-open cell (x_i, x_i+1] x (y_j, y_j+1];
!< - a segment is walked through the grid cell by cell; whether it leaves a cell through the
!<   vertical or the horizontal line first is decided by the exact side of the cell corner
!<   (X + eps, Y + eps^2) with respect to the segment, the tie-break of `fossil_line_query`;
!< - the inside flag of every grid node is the parity of the crossings of its horizontal line
!<   left of it, decided by the same corner test.
!< The crossing points use the formulas of `fossil_line_query` bit for bit, so the wall pieces and
!< the edge intervals share their end points. The measures follow by Green's theorem: the solid
!< area of S n C is half the circulation of (u dv - v du) along the wall pieces (solid on the left)
!< and the solid parts of the cell edges (counter-clockwise); the fluid quantities are the cell's
!< minus the solid's, piece by piece along the fluid boundary arcs and the reversed wall chains.
!<
!< Pieces: union-find over the fluid parts of the cell boundary, joined by the wall chains that
!< connect them; a clockwise loop (a hole) entirely inside the cell is a piece of its own; a
!< counter-clockwise loop (an island) entirely inside the cell is subtracted from the piece that
!< contains it (exact parity tests). A cell where a chain ends inside it (an open curve) or where
!< the boundary states do not close is flagged not closed and gets no measures.
use fossil_predicates, only : orient2d_xy
use penf,              only : I4P, R8P

implicit none
private
public :: cut_cells_2d_t, compute_cut_cells_2d

type :: cut_cells_2d_t
   !< Cut-cell geometry of a block (fluid = outside the solid).
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P !< Cells.
   real(R8P),    allocatable :: alpha(:,:)          !< Fluid area fraction [nx, ny] (-1 if not closed).
   real(R8P),    allocatable :: aperture(:,:,:)     !< Edge fluid fractions [4 (-u, +u, -v, +v), nx, ny].
   real(R8P),    allocatable :: centroid(:,:,:)     !< Fluid centroid [2, nx, ny] (cell centre if no fluid).
   real(R8P),    allocatable :: wet_length(:,:)     !< Wetted length [nx, ny].
   real(R8P),    allocatable :: wet_centroid(:,:,:) !< Wetted centroid [2, nx, ny].
   real(R8P),    allocatable :: wet_normal(:,:,:)   !< Integral of the outward normal over the wall [2, nx, ny].
   integer(I4P), allocatable :: pieces(:,:)         !< Fluid pieces per cell [nx, ny].
   integer(I4P), allocatable :: first_piece(:,:)    !< Index of the first piece of the cell [nx, ny].
   logical,      allocatable :: closed(:,:)         !< The cell boundary closed [nx, ny].
   integer(I4P)              :: np = 0_I4P          !< Pieces in the block.
   integer(I4P), allocatable :: piece_cell(:,:)     !< Cell of each piece [2, np].
   real(R8P),    allocatable :: piece_area(:)       !< Fluid area of each piece [np].
   real(R8P),    allocatable :: piece_centroid(:,:) !< Fluid centroid of each piece [2, np].
   real(R8P),    allocatable :: piece_wet_length(:) !< Wetted length of each piece [np].
   real(R8P),    allocatable :: piece_wet_normal(:,:) !< Outward-normal integral of each piece [2, np].
   real(R8P),    allocatable :: piece_aperture(:,:) !< Fluid fraction of the cell edges per piece [4 (-u, +u, -v, +v), np].
   integer(I4P), allocatable :: piece_iv_start(:)   !< Fluid edge intervals of piece p: piece_iv_start(p) .. (p+1) - 1.
   integer(I4P), allocatable :: piece_iv_edge(:)    !< Edge of an interval: vertical edges (i, j) = i + (j - 1) (nx + 1),
                                                    !< then horizontal ones (i, j) = (nx + 1) ny + i + (j - 1) nx.
   real(R8P),    allocatable :: piece_iv_lo(:)      !< Interval start along its edge (v on a vertical edge, u on a horizontal).
   real(R8P),    allocatable :: piece_iv_hi(:)      !< Interval end.
endtype cut_cells_2d_t

type :: xing_t
   !< A crossing of a segment with a block edge.
   integer(I4P) :: edge = 0_I4P  !< Edge id (vertical edges first, then horizontal).
   real(R8P)    :: pos = 0._R8P  !< Position along the edge (v on a vertical edge, u on a horizontal one).
   integer(I4P) :: seg = 0_I4P   !< Segment.
   real(R8P)    :: p(2) = 0._R8P !< Point.
endtype xing_t

type :: piece_t
   !< A part of a segment inside one cell.
   integer(I4P) :: cell = 0_I4P   !< Cell id.
   integer(I4P) :: seg = 0_I4P    !< Segment.
   real(R8P)    :: a(2) = 0._R8P, b(2) = 0._R8P !< End points (curve orientation).
   integer(I4P) :: ta = 0_I4P, tb = 0_I4P !< End tags: crossing id (> 0) or -vertex (< 0).
endtype piece_t

type :: ctx_t
   !< Working data of one block, shared by the cell routines.
   integer(I4P)               :: nx = 0_I4P, ny = 0_I4P !< Cells.
   real(R8P),     allocatable :: x(:), y(:)       !< Edges.
   logical,       allocatable :: node_in(:,:)     !< Node inside flags.
   type(xing_t),  allocatable :: xs(:)            !< Crossings.
   type(piece_t), allocatable :: ps(:)            !< Pieces.
   integer(I4P),  allocatable :: xstart(:), xord(:), pstart(:), pord(:) !< Indexes.
   integer(I4P)               :: npb = 0_I4P      !< Pieces in the buffers.
   real(R8P),     allocatable :: pa(:), pcx(:), pcy(:), pwl(:), pwx(:), pwy(:) !< Piece buffers.
   integer(I4P),  allocatable :: pci(:), pcj(:)   !< Piece cells.
   real(R8P),     allocatable :: pap(:,:)         !< Piece apertures [4, capacity].
   integer(I4P)               :: niv = 0_I4P      !< Piece intervals in the buffers.
   integer(I4P),  allocatable :: ivp(:), ive(:)   !< Interval piece, edge.
   real(R8P),     allocatable :: ivl(:), ivh(:)   !< Interval range.
endtype ctx_t

contains
   subroutine compute_cut_cells_2d(u, v, next, x, y, cc, status)
   !< Cut cells of the curve (u, v, next) on the block of edges x(nx+1), y(ny+1).
   !<
   !< Segment i runs from vertex i to vertex next(i) (next(i) = 0: no segment). `status` = 0, or 1 on
   !< bad edges.
   real(R8P),            intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)  :: next(:)    !< Next vertex (0: none).
   real(R8P),            intent(in)  :: x(:), y(:) !< Cell edges (strictly increasing).
   type(cut_cells_2d_t), intent(out) :: cc         !< Cut cells.
   integer(I4P),         intent(out) :: status     !< 0 ok, 1 bad edges.
   type(ctx_t)  :: cx                  !< Working data.
   integer(I4P) :: nx, ny, nxs, nps, ne, ncell, c, npb !< Work.

   status = 0_I4P
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   if (nx < 1 .or. ny < 1) then
      status = 1_I4P ; return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny))) then
      status = 1_I4P ; return
   endif
   cc%nx = nx ; cc%ny = ny ; cx%nx = nx ; cx%ny = ny ; cx%x = x ; cx%y = y
   allocate(cc%alpha(nx, ny), cc%aperture(4, nx, ny), cc%centroid(2, nx, ny), cc%wet_length(nx, ny), &
            cc%wet_centroid(2, nx, ny), cc%wet_normal(2, nx, ny), cc%pieces(nx, ny), cc%first_piece(nx, ny), &
            cc%closed(nx, ny))
   call node_flags(u, v, next, x, y, cx%node_in)
   call walk_all(u, v, next, x, y, cx%xs, nxs, cx%ps, nps)
   ne = (nx + 1_I4P) * ny + nx * (ny + 1_I4P) ; ncell = nx * ny
   call index_crossings(cx%xs, nxs, ne, cx%xstart, cx%xord)
   call index_pieces(cx%ps, nps, ncell, cx%pstart, cx%pord)
   allocate(cx%pa(16), cx%pcx(16), cx%pcy(16), cx%pwl(16), cx%pwx(16), cx%pwy(16), cx%pci(16), cx%pcj(16), cx%pap(4, 16))
   allocate(cx%ivp(16), cx%ive(16), cx%ivl(16), cx%ivh(16))
   cx%npb = 0_I4P ; cx%niv = 0_I4P
   do c = 1_I4P, ncell
      call one_cell(c, cx, cc)
   enddo
   npb = cx%npb ; cc%np = npb
   allocate(cc%piece_cell(2, npb))
   cc%piece_cell(1, :) = cx%pci(1:npb) ; cc%piece_cell(2, :) = cx%pcj(1:npb)
   cc%piece_area = cx%pa(1:npb)
   allocate(cc%piece_centroid(2, npb), cc%piece_wet_length(npb), cc%piece_wet_normal(2, npb))
   cc%piece_centroid(1, :) = cx%pcx(1:npb) ; cc%piece_centroid(2, :) = cx%pcy(1:npb)
   cc%piece_wet_length = cx%pwl(1:npb)
   cc%piece_wet_normal(1, :) = cx%pwx(1:npb) ; cc%piece_wet_normal(2, :) = cx%pwy(1:npb)
   cc%piece_aperture = cx%pap(:, 1:npb)
   ! intervals per piece (CSR; appended piece by piece, so already grouped)
   allocate(cc%piece_iv_start(npb + 1)) ; cc%piece_iv_start = 0_I4P
   do c = 1_I4P, cx%niv
      cc%piece_iv_start(cx%ivp(c) + 1) = cc%piece_iv_start(cx%ivp(c) + 1) + 1_I4P
   enddo
   cc%piece_iv_start(1) = 1_I4P
   do c = 1_I4P, npb
      cc%piece_iv_start(c + 1) = cc%piece_iv_start(c + 1) + cc%piece_iv_start(c)
   enddo
   cc%piece_iv_edge = cx%ive(1:cx%niv) ; cc%piece_iv_lo = cx%ivl(1:cx%niv) ; cc%piece_iv_hi = cx%ivh(1:cx%niv)
   endsubroutine compute_cut_cells_2d

   subroutine one_cell(cid, cx, cc)
   !< Measures and pieces of the cell cid.
   integer(I4P),         intent(in)    :: cid !< Cell id.
   type(ctx_t),          intent(inout) :: cx  !< Working data.
   type(cut_cells_2d_t), intent(inout) :: cc  !< Cut cells.
   integer(I4P), parameter :: MAXS = 512 !< Sub-intervals per cell.
   integer(I4P), parameter :: ORDER(4) = [4_I4P, 2_I4P, 1_I4P, 3_I4P] !< Edges (counter-clockwise ids) as (-u, +u, -v, +v).
   integer(I4P) :: i, j, e, k, m, ns, q, r, nb, h, s0 !< Work.
   real(R8P)    :: sa(2, MAXS), sb(2, MAXS) !< Sub-interval end points (counter-clockwise).
   logical      :: sfl(MAXS)                !< Sub-interval is fluid.
   integer(I4P) :: sx(MAXS)                 !< Crossing at the end of sub-interval k (0: corner).
   integer(I4P) :: sed(MAXS)                !< Edge (1 bottom, 2 right, 3 top, 4 left) of sub-interval k.
   integer(I4P) :: par(MAXS)                !< Union-find parents.
   integer(I4P) :: xsub(2, MAXS)            !< For each crossing slot: (crossing id, sub-interval before it).
   real(R8P)    :: x0, x1, y0, y1, area_c, len(4), flen(4), gs(3) !< Work.
   logical      :: state, ok                !< Work.
   integer(I4P) :: nxc                      !< Crossing slots.
   integer(I4P) :: nsb                      !< Boundary sub-intervals.
   integer(I4P) :: eid(4)                   !< Edge ids of the cell.
   real(R8P)    :: acc(7, MAXS)             !< Per root: area, mx, my, wet len, wet mx, wet my, (unused).
   real(R8P)    :: wn(2, MAXS)              !< Per root: wet normal.
   logical,      allocatable :: pvis(:)     !< Visited pieces.
   integer(I4P) :: p0, p1                   !< Work.
   real(R8P)    :: lg(3), lw(4), lwn(2)     !< Loop accumulators.
   real(R8P)    :: hole_area(MAXS)          !< Areas of the hole loops.
   integer(I4P) :: hole_root(MAXS), hole_first(MAXS), hole_n(MAXS) !< Hole pieces.
   integer(I4P) :: nh, nisl, isl_first(MAXS), isl_n(MAXS), loopp(MAXS), nloopp !< Loops.
   real(R8P)    :: isl_g(3, MAXS), isl_w(4, MAXS), isl_wn(2, MAXS) !< Island contributions.

   i = mod(cid - 1_I4P, cx%nx) + 1_I4P ; j = (cid - 1_I4P) / cx%nx + 1_I4P
   x0 = cx%x(i) ; x1 = cx%x(i + 1) ; y0 = cx%y(j) ; y1 = cx%y(j + 1)
   area_c = (x1 - x0) * (y1 - y0)
   len = [y1 - y0, y1 - y0, x1 - x0, x1 - x0]
   cc%closed(i, j) = .true. ; cc%pieces(i, j) = 0_I4P ; cc%first_piece(i, j) = cx%npb + 1_I4P
   cc%wet_length(i, j) = 0._R8P ; cc%wet_centroid(:, i, j) = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
   cc%wet_normal(:, i, j) = 0._R8P
   ! edges counter-clockwise: bottom (+u), right (+v), top (-u), left (-v)
   eid = [hedge(cx, i, j), vedge(cx, i + 1_I4P, j), hedge(cx, i, j + 1_I4P), vedge(cx, i, j)]
   ns = 0 ; nxc = 0 ; ok = .true. ; flen = 0._R8P
   state = cx%node_in(i, j)
   do e = 1, 4
      call edge_subs(e, eid(e))
   enddo
   if (state .neqv. cx%node_in(i, j)) ok = .false.
   nsb = ns
   ! wall pieces of the cell
   nb = cx%pstart(cid + 1) - cx%pstart(cid)
   if (.not. ok .or. ns > MAXS - 4) then
      call mark_open
      return
   endif
   if (nb == 0 .and. nxc == 0) then
      if (cx%node_in(i, j)) then
         cc%alpha(i, j) = 0._R8P ; cc%aperture(:, i, j) = 0._R8P
         cc%centroid(:, i, j) = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
      else
         cc%alpha(i, j) = 1._R8P ; cc%aperture(:, i, j) = 1._R8P
         cc%centroid(:, i, j) = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
         call add_piece(cx, i, j, area_c, 0.5_R8P * (x0 + x1) * area_c, 0.5_R8P * (y0 + y1) * area_c, 0._R8P, [0._R8P, 0._R8P])
         cx%pap(:, cx%npb) = 1._R8P
         call add_iv(cx, cx%npb, eid(4), y0, y1) ; call add_iv(cx, cx%npb, eid(2), y0, y1)
         call add_iv(cx, cx%npb, eid(1), x0, x1) ; call add_iv(cx, cx%npb, eid(3), x0, x1)
         cc%pieces(i, j) = 1_I4P
      endif
      return
   endif
   ! union-find over the fluid sub-intervals; consecutive ones through corners join
   do k = 1, ns
      par(k) = k
   enddo
   do k = 1, nsb
      m = mod(k, nsb) + 1
      if (sx(k) == 0 .and. sfl(k) .and. sfl(m)) call unite(k, m)
   enddo
   ! wall chains: from a piece starting at a crossing, follow the vertices to a crossing
   allocate(pvis(nb)) ; pvis = .false.
   acc = 0._R8P ; wn = 0._R8P
   do q = 1, nb
      p0 = cx%pord(cx%pstart(cid) + q - 1)
      if (cx%ps(p0)%ta <= 0) cycle
      call follow(q, p0, p1, lg, lw, lwn, ok)
      if (.not. ok) then
         call mark_open
         return
      endif
      ! join the fluid sub-intervals next to both end crossings, charge the chain to that root
      r = fluid_sub_at(cx%ps(p0)%ta) ; h = fluid_sub_at(cx%ps(p1)%tb)
      if (r == 0 .or. h == 0) then
         call mark_open
         return
      endif
      call unite(r, h)
      call charge_later(find(r), lg, lw, lwn)
   enddo
   ! remaining pieces: closed loops entirely inside the cell
   nh = 0 ; nisl = 0 ; nloopp = 0
   do q = 1, nb
      if (pvis(q)) cycle
      p0 = cx%pord(cx%pstart(cid) + q - 1)
      s0 = nloopp + 1
      call follow_loop(q, lg, lw, lwn, ok)
      if (.not. ok) then
         call mark_open
         return
      endif
      if (lg(1) < 0._R8P) then ! clockwise: a hole, a piece of its own
         nh = nh + 1 ; hole_area(nh) = -lg(1) ; hole_first(nh) = s0 ; hole_n(nh) = nloopp - s0 + 1
         ns = ns + 1 ; par(ns) = ns ; sfl(ns) = .false. ; hole_root(nh) = ns
         call charge_later(ns, lg, lw, lwn)
      else
         nisl = nisl + 1 ; isl_first(nisl) = s0 ; isl_n(nisl) = nloopp - s0 + 1
         isl_g(:, nisl) = lg ; isl_w(:, nisl) = lw ; isl_wn(:, nisl) = lwn
      endif
   enddo
   ! islands: into the innermost containing hole, else the containing boundary piece
   do k = 1, nisl
      r = 0
      call island_host(k, r)
      if (r == 0) then
         call mark_open
         return
      endif
      call charge_later(find(r), isl_g(:, k), isl_w(:, k), isl_wn(:, k))
   enddo
   ! fluid sub-intervals: Green along them (counter-clockwise)
   do k = 1, ns - nh
      if (.not. sfl(k)) cycle
      call green(sa(:, k), sb(:, k), gs)
      r = find(k)
      acc(1:3, r) = acc(1:3, r) + gs
   enddo
   ! collect the pieces (roots)
   cc%alpha(i, j) = 0._R8P ; cc%centroid(:, i, j) = 0._R8P
   do k = 1, ns
      if (k <= ns - nh) then
         if (.not. sfl(k)) cycle
      endif
      if (find(k) /= k) cycle
      call add_piece(cx, i, j, acc(1, k), acc(2, k), acc(3, k), acc(4, k), wn(:, k))
      ! its fluid edge intervals (issue #70), in the order (-u, +u, -v, +v) = edges (4, 2, 1, 3)
      cx%pap(:, cx%npb) = 0._R8P
      do e = 1, 4
         m = ORDER(e)
         do q = 1, nsb
            if (.not. sfl(q) .or. sed(q) /= m .or. find(q) /= k) cycle
            cx%pap(e, cx%npb) = cx%pap(e, cx%npb) + norm2(sb(:, q) - sa(:, q)) / len(e)
            if (m == 2 .or. m == 4) then
               call add_iv(cx, cx%npb, eid(m), min(sa(2, q), sb(2, q)), max(sa(2, q), sb(2, q)))
            else
               call add_iv(cx, cx%npb, eid(m), min(sa(1, q), sb(1, q)), max(sa(1, q), sb(1, q)))
            endif
         enddo
      enddo
      cc%pieces(i, j) = cc%pieces(i, j) + 1_I4P
      cc%alpha(i, j) = cc%alpha(i, j) + acc(1, k)
      cc%centroid(:, i, j) = cc%centroid(:, i, j) + acc(2:3, k)
      cc%wet_length(i, j) = cc%wet_length(i, j) + acc(4, k)
      cc%wet_normal(:, i, j) = cc%wet_normal(:, i, j) + wn(:, k)
   enddo
   if (cc%alpha(i, j) > 0._R8P) then
      cc%centroid(:, i, j) = cc%centroid(:, i, j) / cc%alpha(i, j)
   else
      cc%centroid(:, i, j) = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
   endif
   ! wetted centroid from the per-root moments
   if (cc%wet_length(i, j) > 0._R8P) then
      cc%wet_centroid(:, i, j) = [sum(acc(5, 1:ns), mask=roots()), sum(acc(6, 1:ns), mask=roots())] / cc%wet_length(i, j)
   endif
   cc%alpha(i, j) = cc%alpha(i, j) / area_c
   cc%aperture(:, i, j) = [flen(4), flen(2), flen(1), flen(3)] / len ! counter-clockwise edges -> (-u, +u, -v, +v)
   contains
      subroutine edge_subs(ed, id)
      !< Sub-intervals of edge ed (1 bottom, 2 right, 3 top, 4 left) in counter-clockwise order.
      integer(I4P), intent(in) :: ed, id !< Edge, edge id.
      integer(I4P) :: kk, n0, n1, idx, st, en, dir !< Work.
      real(R8P)    :: pa_(2), pb_(2), cs(2), ce(2) !< Work.

      select case (ed)
      case (1) ; cs = [x0, y0] ; ce = [x1, y0]
      case (2) ; cs = [x1, y0] ; ce = [x1, y1]
      case (3) ; cs = [x1, y1] ; ce = [x0, y1]
      case (4) ; cs = [x0, y1] ; ce = [x0, y0]
      endselect
      n0 = cx%xstart(id) ; n1 = cx%xstart(id + 1) - 1
      if (ed <= 2) then
         st = n0 ; en = n1 ; dir = 1
      else
         st = n1 ; en = n0 ; dir = -1
      endif
      pa_ = cs
      if (n1 >= n0) then
         do kk = st, en, dir
            idx = cx%xord(kk)
            pb_ = cx%xs(idx)%p
            if (ed == 1 .or. ed == 3) pb_(2) = cs(2)
            if (ed == 2 .or. ed == 4) pb_(1) = cs(1)
            call push(ed, pa_, pb_, .not. state, idx)
            state = .not. state
            pa_ = pb_
         enddo
      endif
      call push(ed, pa_, ce, .not. state, 0_I4P)
      ! the state at the end corner must be the corner's flag
      select case (ed)
      case (1) ; if (state .neqv. cx%node_in(i + 1, j)) ok = .false.
      case (2) ; if (state .neqv. cx%node_in(i + 1, j + 1)) ok = .false.
      case (3) ; if (state .neqv. cx%node_in(i, j + 1)) ok = .false.
      endselect
      endsubroutine edge_subs

      subroutine push(ed, aa, bb, fluid, xid)
      !< Append a sub-interval; record the crossing at its end.
      integer(I4P), intent(in) :: ed           !< Edge.
      real(R8P),    intent(in) :: aa(2), bb(2) !< End points.
      logical,      intent(in) :: fluid        !< Fluid.
      integer(I4P), intent(in) :: xid          !< Crossing at its end (0: corner).

      if (ns >= MAXS - 4) then
         ok = .false. ; return
      endif
      ns = ns + 1 ; sa(:, ns) = aa ; sb(:, ns) = bb ; sfl(ns) = fluid ; sx(ns) = xid ; sed(ns) = ed
      if (fluid) flen(ed) = flen(ed) + norm2(bb - aa)
      if (xid > 0) then
         nxc = nxc + 1 ; xsub(:, nxc) = [xid, ns]
      endif
      endsubroutine push

      integer(I4P) function fluid_sub_at(xid)
      !< The fluid sub-interval next to crossing xid (before or after it, counter-clockwise).
      integer(I4P), intent(in) :: xid !< Crossing id.
      integer(I4P) :: kk, a_, b_      !< Work.

      fluid_sub_at = 0
      do kk = 1, nxc
         if (xsub(1, kk) /= xid) cycle
         a_ = xsub(2, kk) ; b_ = mod(a_, nsb) + 1
         if (sfl(a_)) then
            fluid_sub_at = a_
         elseif (sfl(b_)) then
            fluid_sub_at = b_
         endif
         return
      enddo
      endfunction fluid_sub_at

      subroutine follow(qq, pfirst, plast, g, wl, wnn, good)
      !< Follow the chain from the piece pfirst (local index qq) to a piece ending at a crossing.
      integer(I4P), intent(in)  :: qq, pfirst !< Local index, piece id.
      integer(I4P), intent(out) :: plast      !< Last piece.
      real(R8P),    intent(out) :: g(3), wl(4), wnn(2) !< Curve Green, wet (len, mx, my, -), normal.
      logical,      intent(out) :: good       !< The chain closed on a crossing.
      integer(I4P) :: cur, lq, guard          !< Work.

      g = 0._R8P ; wl = 0._R8P ; wnn = 0._R8P ; good = .true.
      cur = pfirst ; lq = qq ; guard = 0
      do
         pvis(lq) = .true.
         call add_curve(cur, g, wl, wnn)
         if (cx%ps(cur)%tb > 0) exit
         lq = local_start(-cx%ps(cur)%tb)
         guard = guard + 1
         if (lq == 0 .or. guard > nb) then
            good = .false. ; exit
         endif
         cur = cx%pord(cx%pstart(cid) + lq - 1)
      enddo
      plast = cur
      endsubroutine follow

      subroutine follow_loop(qq, g, wl, wnn, good)
      !< Follow a closed loop of pieces starting at local index qq.
      integer(I4P), intent(in)  :: qq !< Local index.
      real(R8P),    intent(out) :: g(3), wl(4), wnn(2) !< Accumulators.
      logical,      intent(out) :: good !< The loop closed.
      integer(I4P) :: cur, lq, guard !< Work.

      g = 0._R8P ; wl = 0._R8P ; wnn = 0._R8P ; good = .true.
      lq = qq ; guard = 0
      do
         if (pvis(lq)) exit
         pvis(lq) = .true.
         cur = cx%pord(cx%pstart(cid) + lq - 1)
         nloopp = nloopp + 1 ; loopp(nloopp) = cur
         call add_curve(cur, g, wl, wnn)
         if (cx%ps(cur)%tb >= 0) then
            good = .false. ; exit
         endif
         lq = local_start(-cx%ps(cur)%tb)
         guard = guard + 1
         if (lq == 0 .or. guard > nb) then
            good = .false. ; exit
         endif
      enddo
      endsubroutine follow_loop

      integer(I4P) function local_start(vtx)
      !< Local index of the piece of the cell that starts at vertex vtx (0: none).
      integer(I4P), intent(in) :: vtx !< Vertex.
      integer(I4P) :: kk !< Counter.

      local_start = 0
      do kk = 1, nb
         if (cx%ps(cx%pord(cx%pstart(cid) + kk - 1))%ta == -vtx) then
            local_start = kk ; return
         endif
      enddo
      endfunction local_start

      subroutine add_curve(pid, g, wl, wnn)
      !< Curve Green, wetted length, moments and normal of piece pid.
      integer(I4P), intent(in)    :: pid !< Piece.
      real(R8P),    intent(inout) :: g(3), wl(4), wnn(2) !< Accumulators.
      real(R8P) :: gg(3), d(2), l !< Work.

      call green(cx%ps(pid)%a, cx%ps(pid)%b, gg)
      g = g + gg
      d = cx%ps(pid)%b - cx%ps(pid)%a ; l = norm2(d)
      wl(1) = wl(1) + l ; wl(2:3) = wl(2:3) + l * 0.5_R8P * (cx%ps(pid)%a + cx%ps(pid)%b)
      wnn = wnn + [d(2), -d(1)]
      endsubroutine add_curve

      subroutine charge_later(rt, g, wl, wnn)
      !< Charge a chain or loop to the root of rt: fluid Green is minus the curve's.
      integer(I4P), intent(in) :: rt !< Sub-interval or hole slot.
      real(R8P),    intent(in) :: g(3), wl(4), wnn(2) !< Curve Green, wet, normal.
      integer(I4P) :: rr !< Root.

      rr = rt
      acc(1:3, rr) = acc(1:3, rr) - g
      acc(4:6, rr) = acc(4:6, rr) + wl(1:3)
      wn(:, rr) = wn(:, rr) + wnn
      endsubroutine charge_later

      subroutine island_host(kk, host)
      !< Slot that receives island kk: the smallest hole containing it, else the boundary piece.
      integer(I4P), intent(in)  :: kk   !< Island.
      integer(I4P), intent(out) :: host !< Slot (0: none found).
      real(R8P)    :: pt(2), best !< Island vertex, smallest hole area.
      integer(I4P) :: hh, kk2, cnt, nroots, rlast !< Work.

      host = 0
      pt = cx%ps(loopp(isl_first(kk)))%a
      best = huge(1._R8P)
      do hh = 1, nh
         if (inside_pieces(pt, loopp(hole_first(hh):hole_first(hh) + hole_n(hh) - 1)) .and. hole_area(hh) < best) then
            best = hole_area(hh) ; host = hole_root(hh)
         endif
      enddo
      if (host /= 0) return
      nroots = 0 ; rlast = 0
      do kk2 = 1, ns - nh
         if (.not. sfl(kk2)) cycle
         if (find(kk2) /= kk2) cycle
         nroots = nroots + 1 ; rlast = kk2
      enddo
      if (nroots == 1) then
         host = rlast ; return
      endif
      ! several boundary pieces: parity of a ray along +u against each piece's boundary
      do kk2 = 1, ns - nh
         if (.not. sfl(kk2)) cycle
         if (find(kk2) /= kk2) cycle
         cnt = ray_count(pt, kk2)
         if (mod(cnt, 2_I4P) == 1_I4P) then
            host = kk2 ; return
         endif
      enddo
      endsubroutine island_host

      integer(I4P) function ray_count(pt, rt)
      !< Crossings of the ray from pt along +u with the boundary of the piece of root rt (its fluid
      !< sub-intervals and the wall pieces charged to it).
      real(R8P),    intent(in) :: pt(2) !< Point.
      integer(I4P), intent(in) :: rt    !< Root.
      integer(I4P) :: kk, qq, pid, rr !< Work.

      ray_count = 0
      do kk = 1, ns - nh
         if (.not. sfl(kk)) cycle
         if (find(kk) /= rt) cycle
         if (ray_hits(pt, sa(:, kk), sb(:, kk))) ray_count = ray_count + 1
      enddo
      do qq = 1, nb
         pid = cx%pord(cx%pstart(cid) + qq - 1)
         if (cx%ps(pid)%ta <= 0) cycle
         rr = fluid_sub_at(cx%ps(pid)%ta)
         if (rr == 0) cycle
         if (find(rr) /= rt) cycle
         call ray_chain(pt, pid, ray_count)
      enddo
      endfunction ray_count

      subroutine ray_chain(pt, pfirst, cnt)
      !< Add the ray crossings of the chain starting at pfirst.
      real(R8P),    intent(in)    :: pt(2)  !< Point.
      integer(I4P), intent(in)    :: pfirst !< First piece.
      integer(I4P), intent(inout) :: cnt    !< Count.
      integer(I4P) :: cur, lq, guard !< Work.

      cur = pfirst ; guard = 0
      do
         if (ray_hits(pt, cx%ps(cur)%a, cx%ps(cur)%b)) cnt = cnt + 1
         if (cx%ps(cur)%tb > 0) exit
         lq = local_start(-cx%ps(cur)%tb)
         guard = guard + 1
         if (lq == 0 .or. guard > nb) exit
         cur = cx%pord(cx%pstart(cid) + lq - 1)
      enddo
      endsubroutine ray_chain

      logical function inside_pieces(pt, plist_)
      !< Parity of the ray from pt along +u against a closed loop of pieces.
      real(R8P),    intent(in) :: pt(2)     !< Point.
      integer(I4P), intent(in) :: plist_(:) !< Pieces of the loop.
      integer(I4P) :: kk, cnt !< Work.

      cnt = 0
      do kk = 1, size(plist_)
         if (ray_hits(pt, cx%ps(plist_(kk))%a, cx%ps(plist_(kk))%b)) cnt = cnt + 1
      enddo
      inside_pieces = mod(cnt, 2_I4P) == 1_I4P
      endfunction inside_pieces

      subroutine mark_open
      !< The cell boundary does not close: no measures.
      cc%closed(i, j) = .false. ; cc%alpha(i, j) = -1._R8P ; cc%aperture(:, i, j) = -1._R8P
      cc%centroid(:, i, j) = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)] ; cc%pieces(i, j) = 0_I4P
      cc%wet_length(i, j) = 0._R8P ; cc%wet_normal(:, i, j) = 0._R8P
      cx%npb = cc%first_piece(i, j) - 1_I4P
      endsubroutine mark_open

      recursive integer(I4P) function find(kk) result(rr)
      !< Union-find root.
      integer(I4P), intent(in) :: kk !< Element.

      rr = kk
      do while (par(rr) /= rr)
         rr = par(rr)
      enddo
      endfunction find

      subroutine unite(a_, b_)
      !< Union (smaller root wins, deterministic).
      integer(I4P), intent(in) :: a_, b_ !< Elements.
      integer(I4P) :: ra, rb !< Roots.

      ra = find(a_) ; rb = find(b_)
      if (ra == rb) return
      if (ra < rb) then
         par(rb) = ra ; acc(:, ra) = acc(:, ra) + acc(:, rb) ; acc(:, rb) = 0._R8P
         wn(:, ra) = wn(:, ra) + wn(:, rb) ; wn(:, rb) = 0._R8P
      else
         par(ra) = rb ; acc(:, rb) = acc(:, rb) + acc(:, ra) ; acc(:, ra) = 0._R8P
         wn(:, rb) = wn(:, rb) + wn(:, ra) ; wn(:, ra) = 0._R8P
      endif
      endsubroutine unite

      function roots() result(mask)
      !< Mask of the piece roots among 1..ns.
      logical :: mask(ns) !< Mask.
      integer(I4P) :: kk  !< Counter.

      do kk = 1, ns
         mask(kk) = find(kk) == kk .and. (kk > ns - nh .or. sfl(kk))
      enddo
      endfunction roots
   endsubroutine one_cell


   subroutine add_iv(cx, ip, ed, lo, hi)
   !< Append the fluid interval [lo, hi] of edge ed to piece ip.
   type(ctx_t),  intent(inout) :: cx     !< Working data.
   integer(I4P), intent(in)    :: ip, ed !< Piece, edge id.
   real(R8P),    intent(in)    :: lo, hi !< Range.

   if (cx%niv == size(cx%ivp)) then
      cx%ivp = [cx%ivp, cx%ivp] ; cx%ive = [cx%ive, cx%ive] ; cx%ivl = [cx%ivl, cx%ivl] ; cx%ivh = [cx%ivh, cx%ivh]
   endif
   cx%niv = cx%niv + 1_I4P
   cx%ivp(cx%niv) = ip ; cx%ive(cx%niv) = ed ; cx%ivl(cx%niv) = lo ; cx%ivh(cx%niv) = hi
   endsubroutine add_iv

   subroutine add_piece(cx, ii, jj, ar, mx, my, wl, wnn)
   !< Append a piece (area, first moments, wetted length, normal).
   type(ctx_t),  intent(inout) :: cx      !< Working data.
   integer(I4P), intent(in) :: ii, jj  !< Cell.
   real(R8P),    intent(in) :: ar, mx, my, wl, wnn(2) !< Measures.

   if (cx%npb == size(cx%pa)) then
      cx%pa = [cx%pa, cx%pa] ; cx%pcx = [cx%pcx, cx%pcx] ; cx%pcy = [cx%pcy, cx%pcy] ; cx%pwl = [cx%pwl, cx%pwl]
      cx%pwx = [cx%pwx, cx%pwx] ; cx%pwy = [cx%pwy, cx%pwy] ; cx%pci = [cx%pci, cx%pci] ; cx%pcj = [cx%pcj, cx%pcj]
      cx%pap = reshape([cx%pap, cx%pap], [4, 2 * size(cx%pap, 2)])
   endif
   cx%pap(:, cx%npb + 1) = 0._R8P
   cx%npb = cx%npb + 1_I4P
   cx%pa(cx%npb) = ar ; cx%pwl(cx%npb) = wl ; cx%pwx(cx%npb) = wnn(1) ; cx%pwy(cx%npb) = wnn(2)
   cx%pci(cx%npb) = ii ; cx%pcj(cx%npb) = jj
   if (ar > 0._R8P) then
      cx%pcx(cx%npb) = mx / ar ; cx%pcy(cx%npb) = my / ar
   else
      cx%pcx(cx%npb) = 0._R8P ; cx%pcy(cx%npb) = 0._R8P
   endif
   endsubroutine add_piece

   integer(I4P) function vedge(cx, ll, kk)
   !< Id of the vertical edge on line ll (1..nx+1), row kk (1..ny).
   type(ctx_t),  intent(in) :: cx     !< Working data.
   integer(I4P), intent(in) :: ll, kk !< Line, row.

   vedge = ll + (kk - 1_I4P) * (cx%nx + 1_I4P)
   endfunction vedge

   integer(I4P) function hedge(cx, kk, ll)
   !< Id of the horizontal edge in column kk (1..nx), on line ll (1..ny+1).
   type(ctx_t),  intent(in) :: cx     !< Working data.
   integer(I4P), intent(in) :: kk, ll !< Column, line.

   hedge = (cx%nx + 1_I4P) * cx%ny + kk + (ll - 1_I4P) * cx%nx
   endfunction hedge

   pure subroutine green(a, b, g)
   !< Green integrals along the straight path a -> b: g(1) = 1/2 (u dv - v du) (area), g(2) = 1/2 u^2 dv
   !< (moment in u), g(3) = -1/2 v^2 du (moment in v).
   real(R8P), intent(in)  :: a(2), b(2) !< End points.
   real(R8P), intent(out) :: g(3)       !< Integrals.

   g(1) = 0.5_R8P * (a(1) * b(2) - b(1) * a(2))
   g(2) = (b(2) - a(2)) * (a(1) * a(1) + a(1) * b(1) + b(1) * b(1)) / 6._R8P
   g(3) = -(b(1) - a(1)) * (a(2) * a(2) + a(2) * b(2) + b(2) * b(2)) / 6._R8P
   endsubroutine green

   pure logical function ray_hits(pt, a, b)
   !< The ray from pt along +u crosses the segment ab (v = pt(2) + eps rule; exact side test).
   real(R8P), intent(in) :: pt(2), a(2), b(2) !< Point, segment.
   real(R8P)             :: o                 !< Orientation.

   ray_hits = .false.
   if ((a(2) > pt(2)) .eqv. (b(2) > pt(2))) return
   o = orient2d_xy(a(1), a(2), b(1), b(2), pt(1), pt(2))
   if (o == 0._R8P) return
   ray_hits = (b(2) > a(2) .and. o > 0._R8P) .or. (b(2) < a(2) .and. o < 0._R8P)
   endfunction ray_hits

   pure integer(I4P) function corner_side(pu, pv, qu, qv, cu, cv)
   !< Exact side (+1 left, -1 right) of the corner (cu + eps, cv + eps^2) of the directed segment pq.
   real(R8P), intent(in) :: pu, pv, qu, qv, cu, cv !< Segment, corner.
   real(R8P)             :: o                      !< Orientation.

   o = orient2d_xy(pu, pv, qu, qv, cu, cv)
   if (o > 0._R8P) then
      corner_side = 1
   elseif (o < 0._R8P) then
      corner_side = -1
   elseif (pv > qv) then
      corner_side = 1
   elseif (pv < qv) then
      corner_side = -1
   elseif (qu > pu) then
      corner_side = 1
   else
      corner_side = -1
   endif
   endfunction corner_side

   pure integer(I4P) function col_of(w, e)
   !< Column of coordinate w: 0 if w <= e(1), k if e(k) < w <= e(k+1), n+1 beyond.
   real(R8P), intent(in) :: w    !< Coordinate.
   real(R8P), intent(in) :: e(:) !< Edges.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   if (w <= e(1)) then
      col_of = 0 ; return
   endif
   if (w > e(size(e))) then
      col_of = size(e, kind=I4P) ; return
   endif
   lo = 1 ; hi = size(e, kind=I4P)
   do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (w > e(mid)) then
         lo = mid
      else
         hi = mid
      endif
   enddo
   col_of = lo
   endfunction col_of

   subroutine node_flags(u, v, next, x, y, node_in)
   !< Inside flags of the grid nodes: parity of the crossings of their horizontal line left of them.
   real(R8P),            intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)  :: next(:)    !< Next vertex.
   real(R8P),            intent(in)  :: x(:), y(:) !< Edges.
   logical, allocatable, intent(out) :: node_in(:,:) !< Flags [nx+1, ny+1].
   integer(I4P), allocatable :: segs(:) !< Segments crossing the line.
   integer(I4P)              :: j, i, s, n, k, cnt, sd !< Work.

   allocate(node_in(size(x), size(y)), segs(size(u)))
   do j = 1, size(y)
      n = 0
      do s = 1, size(u)
         if (next(s) == 0) cycle
         if ((v(s) > y(j)) .eqv. (v(next(s)) > y(j))) cycle
         n = n + 1 ; segs(n) = s
      enddo
      do i = 1, size(x)
         cnt = 0
         do k = 1, n
            s = segs(k)
            sd = corner_side(u(s), v(s), u(next(s)), v(next(s)), x(i), y(j))
            ! the crossing is left of the node when the node is right of an upward segment, left of a downward one
            if ((v(next(s)) > v(s) .and. sd < 0) .or. (v(next(s)) < v(s) .and. sd > 0)) cnt = cnt + 1
         enddo
         node_in(i, j) = mod(cnt, 2_I4P) == 1_I4P
      enddo
   enddo
   endsubroutine node_flags

   subroutine walk_all(u, v, next, x, y, xs, nxs, ps, nps)
   !< Walk every segment through the grid: crossings with block edges, pieces in block cells.
   real(R8P),                  intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),               intent(in)  :: next(:)    !< Next vertex.
   real(R8P),                  intent(in)  :: x(:), y(:) !< Edges.
   type(xing_t),  allocatable, intent(out) :: xs(:)      !< Crossings.
   integer(I4P),               intent(out) :: nxs        !< Crossings count.
   type(piece_t), allocatable, intent(out) :: ps(:)      !< Pieces.
   integer(I4P),               intent(out) :: nps        !< Pieces count.
   integer(I4P) :: nx, ny, s, w, ci, cj, qi, qj, lv, lh, sd, ta, guard !< Work.
   real(R8P)    :: pu, pv, qu, qv, du, dv, a(2), cpt(2), xx, yy !< Work.
   logical      :: vfirst, cb, hasv, hash !< Work.
   integer(I4P) :: last_id                !< Id of the last crossing recorded (0: not on a block edge).

   nx = size(x, kind=I4P) - 1 ; ny = size(y, kind=I4P) - 1
   allocate(xs(64), ps(64)) ; nxs = 0 ; nps = 0
   do s = 1, size(u, kind=I4P)
      w = next(s)
      if (w == 0) cycle
      pu = u(s) ; pv = v(s) ; qu = u(w) ; qv = v(w) ; du = qu - pu ; dv = qv - pv
      if (du == 0._R8P .and. dv == 0._R8P) cycle
      if (max(pu, qu) < x(1) .or. min(pu, qu) > x(nx + 1) .or. max(pv, qv) < y(1) .or. min(pv, qv) > y(ny + 1)) cycle
      ci = col_of(pu, x) ; cj = col_of(pv, y) ; qi = col_of(qu, x) ; qj = col_of(qv, y)
      a = [pu, pv] ; ta = -s ; guard = 0
      do
         if (ci == qi .and. cj == qj) then
            call add_p(ci, cj, a, [qu, qv], ta, -w)
            exit
         endif
         guard = guard + 1
         if (guard > 4 * (nx + ny) + 8) exit
         ! candidate lines
         hasv = .false. ; hash = .false.
         if (du > 0._R8P .and. ci < qi) then
            hasv = .true. ; lv = ci + 1 ; xx = x(lv)
         elseif (du < 0._R8P .and. ci > qi) then
            hasv = .true. ; lv = ci ; xx = x(lv)
         endif
         if (dv > 0._R8P .and. cj < qj) then
            hash = .true. ; lh = cj + 1 ; yy = y(lh)
         elseif (dv < 0._R8P .and. cj > qj) then
            hash = .true. ; lh = cj ; yy = y(lh)
         endif
         if (hasv .and. hash) then
            sd = corner_side(pu, pv, qu, qv, xx, yy)
            cb = (du > 0._R8P .and. sd < 0) .or. (du < 0._R8P .and. sd > 0) ! corner below the segment line
            if (dv > 0._R8P) then
               vfirst = .not. cb
            else
               vfirst = cb
            endif
         else
            vfirst = hasv
         endif
         if (vfirst) then
            cpt = [xx, pv + (xx - pu) * (qv - pv) / (qu - pu)]  ! fossil_line_query, axis 2
            call add_x(1, lv, cj, cpt(2), cpt)
            call add_p(ci, cj, a, cpt, ta, last_id)
            if (du > 0._R8P) then
               ci = ci + 1
            else
               ci = ci - 1
            endif
         else
            cpt = [pu + (yy - pv) * (qu - pu) / (qv - pv), yy]  ! fossil_line_query, axis 1
            call add_x(2, lh, ci, cpt(1), cpt)
            call add_p(ci, cj, a, cpt, ta, last_id)
            if (dv > 0._R8P) then
               cj = cj + 1
            else
               cj = cj - 1
            endif
         endif
         a = cpt ; ta = last_id
      enddo
   enddo
   contains
      subroutine add_x(kind, ll, kk, pos, pt)
      !< Record the crossing with line ll (vertical: kind 1, interval kk = row; horizontal: kind 2,
      !< kk = column) if it lies on a block edge.
      integer(I4P), intent(in) :: kind, ll, kk !< Kind, line, interval.
      real(R8P),    intent(in) :: pos, pt(2)   !< Position, point.
      type(xing_t), allocatable :: tmp(:)      !< Growth.

      last_id = 0
      if (kind == 1) then
         if (ll < 1 .or. ll > nx + 1 .or. kk < 1 .or. kk > ny) return
      else
         if (ll < 1 .or. ll > ny + 1 .or. kk < 1 .or. kk > nx) return
      endif
      if (nxs == size(xs)) then
         allocate(tmp(2 * nxs)) ; tmp(1:nxs) = xs ; call move_alloc(tmp, xs)
      endif
      nxs = nxs + 1
      if (kind == 1) then
         xs(nxs)%edge = ll + (kk - 1) * (nx + 1)
      else
         xs(nxs)%edge = (nx + 1) * ny + kk + (ll - 1) * nx
      endif
      xs(nxs)%pos = pos ; xs(nxs)%seg = s ; xs(nxs)%p = pt
      last_id = nxs
      endsubroutine add_x

      subroutine add_p(ii, jj, aa, bb, t_a, t_b)
      !< Record the piece aa -> bb if it lies in a block cell.
      integer(I4P), intent(in) :: ii, jj     !< Cell indices (0 / n+1 outside).
      real(R8P),    intent(in) :: aa(2), bb(2) !< End points.
      integer(I4P), intent(in) :: t_a, t_b   !< Tags.
      type(piece_t), allocatable :: tmp(:)   !< Growth.

      if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) return
      if (nps == size(ps)) then
         allocate(tmp(2 * nps)) ; tmp(1:nps) = ps ; call move_alloc(tmp, ps)
      endif
      nps = nps + 1
      ps(nps)%cell = ii + (jj - 1) * nx ; ps(nps)%seg = s ; ps(nps)%a = aa ; ps(nps)%b = bb
      ps(nps)%ta = t_a ; ps(nps)%tb = t_b
      endsubroutine add_p
   endsubroutine walk_all

   subroutine index_crossings(xs, n, ne, start, ord)
   !< Order the crossings by edge, then position, then segment; start(e) = first index of edge e.
   type(xing_t),              intent(in)  :: xs(:)    !< Crossings.
   integer(I4P),              intent(in)  :: n, ne    !< Count, edges.
   integer(I4P), allocatable, intent(out) :: start(:) !< [ne + 1].
   integer(I4P), allocatable, intent(out) :: ord(:)   !< [n].
   integer(I4P) :: k, e, i, j, key !< Work.

   allocate(start(ne + 1), ord(n))
   start = 0
   do k = 1, n
      start(xs(k)%edge) = start(xs(k)%edge) + 1
   enddo
   key = 1
   do e = 1, ne
      i = start(e) ; start(e) = key ; key = key + i
   enddo
   start(ne + 1) = key
   block
      integer(I4P) :: fill(ne)
      fill = 0
      do k = 1, n
         e = xs(k)%edge
         ord(start(e) + fill(e)) = k ; fill(e) = fill(e) + 1
      enddo
   end block
   ! insertion sort inside each edge by position, then segment
   do e = 1, ne
      do i = start(e) + 1, start(e + 1) - 1
         key = ord(i) ; j = i - 1
         do while (j >= start(e))
            if (xs(ord(j))%pos < xs(key)%pos .or. (xs(ord(j))%pos == xs(key)%pos .and. xs(ord(j))%seg <= xs(key)%seg)) exit
            ord(j + 1) = ord(j) ; j = j - 1
         enddo
         ord(j + 1) = key
      enddo
   enddo
   endsubroutine index_crossings

   subroutine index_pieces(ps, n, ncell, start, ord)
   !< Order the pieces by cell; start(c) = first index of cell c.
   type(piece_t),             intent(in)  :: ps(:)    !< Pieces.
   integer(I4P),              intent(in)  :: n, ncell !< Count, cells.
   integer(I4P), allocatable, intent(out) :: start(:) !< [ncell + 1].
   integer(I4P), allocatable, intent(out) :: ord(:)   !< [n].
   integer(I4P) :: k, c, key, i !< Work.
   integer(I4P), allocatable :: fill(:) !< Work.

   allocate(start(ncell + 1), ord(n), fill(ncell))
   start = 0
   do k = 1, n
      start(ps(k)%cell) = start(ps(k)%cell) + 1
   enddo
   key = 1
   do c = 1, ncell
      i = start(c) ; start(c) = key ; key = key + i
   enddo
   start(ncell + 1) = key
   fill = 0
   do k = 1, n
      c = ps(k)%cell
      ord(start(c) + fill(c)) = k ; fill(c) = fill(c) + 1
   enddo
   endsubroutine index_pieces
endmodule fossil_cut_cells_2d
