!< FOSSIL, exact multilevel coarsening of the cut cells (issue #70).
module fossil_cut_coarsen
!< FOSSIL, exact multilevel coarsening of the cut cells (issue #70).
!<
!< The geometry of a coarse grid (every `ratio` fine cells along each axis, nested edges) is built from the fine one,
!< never recomputed: coarse-fine cut cells are where re-redistribution fails at tiny volume fractions (llm-wiki
!< `high-order-capturing-schemes/wiki/methods/high-order-amr.md`), and a coarse level that does not add up to its fine
!< cells breaks conservation across the levels.
!<
!< - Coarse pieces: the connected components of the fine pieces of a coarse cell, joined by the fine links between two
!<   fine cells of that coarse cell (`fossil_cut_graph`); `fine_to_coarse` maps every fine piece to its coarse one.
!< - Sums, in the order of the fine pieces: volumes (and so alpha), first moments (centroids), wetted areas, moments and
!<   normal integrals, the fluid areas of the coarse faces (from the fine faces on them) and the higher moments
!<   (re-centred and rescaled to the coarse cell: X_c = alpha + beta X_f per axis, binomial expansion). The coarse
!<   fractions are the fine ones added up.
!< - Links: the fine links across a coarse face, mapped to coarse pieces and merged.
!< - Node flags: the fine ones at the coarse nodes; edge fractions (if the fine cut has them): the fine edges added up.
!< - A coarse cell with an untrusted fine cell (pieces = -1) is untrusted; face piece counts are not derived (-1).
!< Quadrature rules are not coarsened (concatenating the fine rules is exact; a solver can do it).
use fossil_cut_cells_2d, only : cut_cells_2d_t, cut_cells_2d_reset
use fossil_cut_cells_3d, only : cut_cells_3d_t
use fossil_cut_graph,    only : cut_links_t, build_csr
use fossil_cut_moments,  only : mono_count, mono_exponents
use penf,                only : I4P, R8P

implicit none
private
public :: cut_coarsen_3d, cut_coarsen_2d

contains
   subroutine cut_coarsen_3d(fine, links, x, y, z, ratio, coarse, coarse_links, fine_to_coarse, status)
   !< Coarse cut cells of `fine` (block x, y, z) at `ratio` (see the module notes). status 0, or 1 if a cell count is not
   !< a multiple of ratio.
   type(cut_cells_3d_t),      intent(in)  :: fine            !< Fine cut cells.
   type(cut_links_t),         intent(in)  :: links           !< Fine links (`cut_cells_3d(..., links=)`).
   real(R8P),                 intent(in)  :: x(:), y(:), z(:) !< Fine cell edges.
   integer(I4P),              intent(in)  :: ratio           !< Refinement ratio (2: every other edge).
   type(cut_cells_3d_t),      intent(out) :: coarse          !< Coarse cut cells.
   type(cut_links_t),         intent(out) :: coarse_links    !< Coarse links.
   integer(I4P), allocatable, intent(out) :: fine_to_coarse(:) !< Coarse piece of every fine piece.
   integer(I4P),              intent(out) :: status          !< 0 or 1.
   integer(I4P), allocatable :: par(:), root_piece(:), lp(:), lq(:), lf(:), e3(:,:)
   real(R8P),    allocatable :: la(:), acc(:,:), pm(:,:), pw(:,:)
   real(R8P)                 :: hf(3), hc(3), cf(3), ccn(3), af, ac(3), w
   integer(I4P)              :: nf(3), n(3), p, q, k, ci(3), cj(3), ic(3), np, f, ax, nm, i, j, kk, nl, c, cell_count
   logical,      allocatable :: bad(:,:,:)

   status = 0_I4P
   nf = [fine%nx, fine%ny, fine%nz]
   if (any(mod(nf, ratio) /= 0_I4P) .or. ratio < 1_I4P) then
      status = 1_I4P ; return
   endif
   n = nf / ratio
   coarse%nx = n(1) ; coarse%ny = n(2) ; coarse%nz = n(3)
   ! components: union-find over the fine pieces through the links inside a coarse cell
   allocate(par(fine%np)) ; par = [(p, p = 1, fine%np)]
   do p = 1, fine%np
      do k = links%start(p), links%start(p + 1) - 1
         q = links%piece(k)
         if (all((fine%piece_cell(:, p) - 1) / ratio == (fine%piece_cell(:, q) - 1) / ratio)) call unite(p, q)
      enddo
   enddo
   ! untrusted coarse cells
   allocate(bad(n(1), n(2), n(3))) ; bad = .false.
   do k = 1, nf(3)
      do j = 1, nf(2)
         do i = 1, nf(1)
            if (fine%pieces(i, j, k) < 0_I4P) bad((i - 1) / ratio + 1, (j - 1) / ratio + 1, (k - 1) / ratio + 1) = .true.
         enddo
      enddo
   enddo
   ! coarse pieces: roots in the order of the coarse cells (x fastest), then of the roots
   allocate(root_piece(fine%np), fine_to_coarse(fine%np)) ; root_piece = 0_I4P ; fine_to_coarse = 0_I4P
   allocate(coarse%pieces(n(1), n(2), n(3)), coarse%first_piece(n(1), n(2), n(3))) ; coarse%pieces = 0_I4P
   cell_count = product(n)
   block
      integer(I4P) :: cstart(cell_count + 1), corder(fine%np), fill(cell_count)
      cstart = 0_I4P
      do p = 1, fine%np
         c = cell_id((fine%piece_cell(:, p) - 1) / ratio + 1)
         cstart(c + 1) = cstart(c + 1) + 1_I4P
      enddo
      cstart(1) = 1_I4P
      do c = 1, cell_count
         cstart(c + 1) = cstart(c + 1) + cstart(c)
      enddo
      fill = 0_I4P
      do p = 1, fine%np
         c = cell_id((fine%piece_cell(:, p) - 1) / ratio + 1)
         corder(cstart(c) + fill(c)) = p ; fill(c) = fill(c) + 1_I4P
      enddo
      np = 0_I4P
      do c = 1, cell_count
         ci = [mod(c - 1, n(1)) + 1, mod((c - 1) / n(1), n(2)) + 1, (c - 1) / (n(1) * n(2)) + 1]
         coarse%first_piece(ci(1), ci(2), ci(3)) = np + 1_I4P
         if (bad(ci(1), ci(2), ci(3))) then
            coarse%pieces(ci(1), ci(2), ci(3)) = -1_I4P ; cycle
         endif
         do k = cstart(c), cstart(c + 1) - 1
            p = corder(k) ; q = find(p)
            if (root_piece(q) == 0_I4P) then
               np = np + 1_I4P ; root_piece(q) = np
               coarse%pieces(ci(1), ci(2), ci(3)) = coarse%pieces(ci(1), ci(2), ci(3)) + 1_I4P
            endif
            fine_to_coarse(p) = root_piece(q)
         enddo
      enddo
   endblock
   coarse%np = np
   ! sums per coarse piece: volume, first moments, wet area, wet moments, wet normal, six coarse face areas
   allocate(acc(17, np)) ; acc = 0._R8P
   allocate(coarse%piece_cell(3, np))
   do p = 1, fine%np
      q = fine_to_coarse(p) ; if (q == 0_I4P) cycle
      ic = fine%piece_cell(:, p) ; ci = (ic - 1) / ratio + 1
      coarse%piece_cell(:, q) = ci
      hf = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2)), z(ic(3) + 1) - z(ic(3))]
      acc(1, q) = acc(1, q) + fine%piece_volume(p)
      acc(2:4, q) = acc(2:4, q) + fine%piece_volume(p) * fine%piece_centroid(:, p)
      acc(5, q) = acc(5, q) + fine%piece_wet_area(p)
      acc(6:8, q) = acc(6:8, q) + fine%piece_wet_area(p) * fine%piece_wet_centroid(:, p)
      acc(9:11, q) = acc(9:11, q) + fine%piece_wet_normal(:, p)
      do f = 1, 6 ! a fine face on the coarse face: low side at the coarse cell's first fine index, high at its last
         ax = (f + 1) / 2 ; af = hf(mod(ax, 3) + 1) * hf(mod(ax + 1, 3) + 1)
         if (mod(f, 2) == 1 .and. mod(ic(ax) - 1, ratio) == 0) acc(11 + f, q) = acc(11 + f, q) + fine%piece_face(f, p) * af
         if (mod(f, 2) == 0 .and. mod(ic(ax), ratio) == 0) acc(11 + f, q) = acc(11 + f, q) + fine%piece_face(f, p) * af
      enddo
   enddo
   allocate(coarse%piece_volume(np), coarse%piece_centroid(3, np), coarse%piece_wet_area(np), coarse%piece_wet_centroid(3, np), &
            coarse%piece_wet_normal(3, np), coarse%piece_face(6, np))
   do q = 1, np
      ci = coarse%piece_cell(:, q) ; call coarse_cell(ci, cf, hc)
      ccn = cf + 0.5_R8P * hc ; ac = [hc(2) * hc(3), hc(3) * hc(1), hc(1) * hc(2)]
      coarse%piece_volume(q) = acc(1, q)
      coarse%piece_centroid(:, q) = ccn ; if (acc(1, q) > 0._R8P) coarse%piece_centroid(:, q) = acc(2:4, q) / acc(1, q)
      coarse%piece_wet_area(q) = acc(5, q)
      coarse%piece_wet_centroid(:, q) = ccn ; if (acc(5, q) > 0._R8P) coarse%piece_wet_centroid(:, q) = acc(6:8, q) / acc(5, q)
      coarse%piece_wet_normal(:, q) = acc(9:11, q)
      coarse%piece_face(:, q) = acc(12:17, q) / [ac(1), ac(1), ac(2), ac(2), ac(3), ac(3)]
   enddo
   ! cell totals
   allocate(coarse%alpha(n(1), n(2), n(3)), coarse%centroid(3, n(1), n(2), n(3)), coarse%wet_area(n(1), n(2), n(3)), &
            coarse%wet_centroid(3, n(1), n(2), n(3)), coarse%wet_normal(3, n(1), n(2), n(3)), coarse%closed(n(1), n(2), n(3)))
   coarse%closed = .not. bad
   do k = 1, n(3)
      do j = 1, n(2)
         do i = 1, n(1)
            ci = [i, j, k] ; call coarse_cell(ci, cf, hc) ; ccn = cf + 0.5_R8P * hc
            coarse%alpha(i, j, k) = 0._R8P ; coarse%centroid(:, i, j, k) = 0._R8P ; coarse%wet_area(i, j, k) = 0._R8P
            coarse%wet_centroid(:, i, j, k) = 0._R8P ; coarse%wet_normal(:, i, j, k) = 0._R8P
            if (bad(i, j, k)) then
               coarse%alpha(i, j, k) = -1._R8P ; coarse%centroid(:, i, j, k) = ccn ; coarse%wet_centroid(:, i, j, k) = ccn
               cycle
            endif
            do q = coarse%first_piece(i, j, k), coarse%first_piece(i, j, k) + coarse%pieces(i, j, k) - 1
               coarse%alpha(i, j, k) = coarse%alpha(i, j, k) + acc(1, q)
               coarse%centroid(:, i, j, k) = coarse%centroid(:, i, j, k) + acc(2:4, q)
               coarse%wet_area(i, j, k) = coarse%wet_area(i, j, k) + acc(5, q)
               coarse%wet_centroid(:, i, j, k) = coarse%wet_centroid(:, i, j, k) + acc(6:8, q)
               coarse%wet_normal(:, i, j, k) = coarse%wet_normal(:, i, j, k) + acc(9:11, q)
            enddo
            if (coarse%alpha(i, j, k) > 0._R8P) then
               coarse%centroid(:, i, j, k) = coarse%centroid(:, i, j, k) / coarse%alpha(i, j, k)
            else
               coarse%centroid(:, i, j, k) = ccn
            endif
            if (coarse%wet_area(i, j, k) > 0._R8P) then
               coarse%wet_centroid(:, i, j, k) = coarse%wet_centroid(:, i, j, k) / coarse%wet_area(i, j, k)
            else
               coarse%wet_centroid(:, i, j, k) = ccn
            endif
            coarse%alpha(i, j, k) = coarse%alpha(i, j, k) / product(hc)
         enddo
      enddo
   enddo
   call coarse_faces
   ! node flags, edge fractions
   if (allocated(fine%node_in)) coarse%node_in = fine%node_in(1::ratio, 1::ratio, 1::ratio)
   if (allocated(fine%edge_x)) call coarse_edges
   ! moments
   if (fine%moments_order > 0_I4P .and. allocated(fine%piece_moments)) then
      coarse%moments_order = fine%moments_order
      nm = mono_count(3_I4P, fine%moments_order) ; e3 = mono_exponents(3_I4P, fine%moments_order)
      allocate(pm(nm, np), pw(nm, np)) ; pm = 0._R8P ; pw = 0._R8P
      do p = 1, fine%np
         q = fine_to_coarse(p) ; if (q == 0_I4P) cycle
         ic = fine%piece_cell(:, p)
         hf = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2)), z(ic(3) + 1) - z(ic(3))]
         call coarse_cell(coarse%piece_cell(:, q), cf, hc)
         call shift_moments(fine%piece_moments(:, p), [x(ic(1)), y(ic(2)), z(ic(3))] + 0.5_R8P * hf, hf, &
                            cf + 0.5_R8P * hc, hc, e3, pm(:, q))
         call shift_moments(fine%piece_wall_moments(:, p), [x(ic(1)), y(ic(2)), z(ic(3))] + 0.5_R8P * hf, hf, &
                            cf + 0.5_R8P * hc, hc, e3, pw(:, q))
      enddo
      call move_alloc(pm, coarse%piece_moments) ; call move_alloc(pw, coarse%piece_wall_moments)
   endif
   ! links across the coarse faces
   allocate(lp(64), lq(64), lf(64), la(64)) ; nl = 0_I4P
   do p = 1, fine%np
      if (fine_to_coarse(p) == 0_I4P) cycle
      do k = links%start(p), links%start(p + 1) - 1
         q = links%piece(k)
         if (fine_to_coarse(q) == 0_I4P) cycle
         ci = (fine%piece_cell(:, p) - 1) / ratio + 1 ; cj = (fine%piece_cell(:, q) - 1) / ratio + 1
         if (all(ci == cj)) cycle
         if (nl == size(lp)) then
            lp = [lp, lp] ; lq = [lq, lq] ; lf = [lf, lf] ; la = [la, la]
         endif
         nl = nl + 1_I4P ; lp(nl) = fine_to_coarse(p) ; lq(nl) = fine_to_coarse(q) ; lf(nl) = links%face(k)
         la(nl) = links%area(k)
      enddo
   enddo
   coarse_links%np = np
   call build_csr(nl, lp, lq, lf, la, coarse_links)
   contains
      integer(I4P) function find(a) result(r)
      !< Root of a (path halving).
      integer(I4P), intent(in) :: a !< Piece.

      r = a
      do while (par(r) /= r)
         par(r) = par(par(r)) ; r = par(r)
      enddo
      endfunction find

      subroutine unite(a, b)
      !< Join the sets of a and b (the smaller root wins: the order of the pieces).
      integer(I4P), intent(in) :: a, b !< Pieces.
      integer(I4P)             :: ra, rb !< Roots.

      ra = find(a) ; rb = find(b)
      if (ra /= rb) par(max(ra, rb)) = min(ra, rb)
      endsubroutine unite

      integer(I4P) function cell_id(cc) result(id)
      !< Linear index of the coarse cell cc (x fastest).
      integer(I4P), intent(in) :: cc(3) !< Cell.

      id = cc(1) + n(1) * ((cc(2) - 1) + n(2) * (cc(3) - 1))
      endfunction cell_id

      subroutine coarse_cell(cc, lo, h)
      !< Low corner and size of the coarse cell cc.
      integer(I4P), intent(in)  :: cc(3)    !< Cell.
      real(R8P),    intent(out) :: lo(3), h(3) !< Corner, size.

      lo = [x((cc(1) - 1) * ratio + 1), y((cc(2) - 1) * ratio + 1), z((cc(3) - 1) * ratio + 1)]
      h = [x(cc(1) * ratio + 1), y(cc(2) * ratio + 1), z(cc(3) * ratio + 1)] - lo
      endsubroutine coarse_cell

      subroutine coarse_faces
      !< Coarse face fractions, centroids and closed flags: the fine faces on each coarse face added up.
      integer(I4P) :: a1, a2, b1, b2 !< Work.
      real(R8P)    :: s, sc(3), fa  !< Work.

      allocate(coarse%face_x(n(1) + 1, n(2), n(3)), coarse%face_y(n(1), n(2) + 1, n(3)), coarse%face_z(n(1), n(2), n(3) + 1))
      allocate(coarse%face_x_centroid(3, n(1) + 1, n(2), n(3)), coarse%face_y_centroid(3, n(1), n(2) + 1, n(3)), &
               coarse%face_z_centroid(3, n(1), n(2), n(3) + 1))
      allocate(coarse%face_x_closed(n(1) + 1, n(2), n(3)), coarse%face_y_closed(n(1), n(2) + 1, n(3)), &
               coarse%face_z_closed(n(1), n(2), n(3) + 1))
      allocate(coarse%face_x_pieces(n(1) + 1, n(2), n(3)), coarse%face_y_pieces(n(1), n(2) + 1, n(3)), &
               coarse%face_z_pieces(n(1), n(2), n(3) + 1))
      coarse%face_x_pieces = -1_I4P ; coarse%face_y_pieces = -1_I4P ; coarse%face_z_pieces = -1_I4P
      do k = 1, n(3)
         do j = 1, n(2)
            do i = 1, n(1) + 1
               s = 0._R8P ; sc = 0._R8P ; coarse%face_x_closed(i, j, k) = .true.
               do b2 = (k - 1) * ratio + 1, k * ratio
                  do b1 = (j - 1) * ratio + 1, j * ratio
                     a1 = (i - 1) * ratio + 1 ; fa = (y(b1 + 1) - y(b1)) * (z(b2 + 1) - z(b2)) * fine%face_x(a1, b1, b2)
                     s = s + fa ; sc = sc + fa * fine%face_x_centroid(:, a1, b1, b2)
                     coarse%face_x_closed(i, j, k) = coarse%face_x_closed(i, j, k) .and. fine%face_x_closed(a1, b1, b2)
                  enddo
               enddo
               fa = (y(j * ratio + 1) - y((j - 1) * ratio + 1)) * (z(k * ratio + 1) - z((k - 1) * ratio + 1))
               coarse%face_x(i, j, k) = s / fa
               coarse%face_x_centroid(:, i, j, k) = [x((i - 1) * ratio + 1), &
                                                     0.5_R8P * (y(j * ratio + 1) + y((j - 1) * ratio + 1)), &
                                                     0.5_R8P * (z(k * ratio + 1) + z((k - 1) * ratio + 1))]
               if (s > 0._R8P) coarse%face_x_centroid(:, i, j, k) = sc / s
            enddo
         enddo
      enddo
      do k = 1, n(3)
         do j = 1, n(2) + 1
            do i = 1, n(1)
               s = 0._R8P ; sc = 0._R8P ; coarse%face_y_closed(i, j, k) = .true.
               do b2 = (k - 1) * ratio + 1, k * ratio
                  do b1 = (i - 1) * ratio + 1, i * ratio
                     a2 = (j - 1) * ratio + 1 ; fa = (x(b1 + 1) - x(b1)) * (z(b2 + 1) - z(b2)) * fine%face_y(b1, a2, b2)
                     s = s + fa ; sc = sc + fa * fine%face_y_centroid(:, b1, a2, b2)
                     coarse%face_y_closed(i, j, k) = coarse%face_y_closed(i, j, k) .and. fine%face_y_closed(b1, a2, b2)
                  enddo
               enddo
               fa = (x(i * ratio + 1) - x((i - 1) * ratio + 1)) * (z(k * ratio + 1) - z((k - 1) * ratio + 1))
               coarse%face_y(i, j, k) = s / fa
               coarse%face_y_centroid(:, i, j, k) = [0.5_R8P * (x(i * ratio + 1) + x((i - 1) * ratio + 1)), &
                                                     y((j - 1) * ratio + 1), &
                                                     0.5_R8P * (z(k * ratio + 1) + z((k - 1) * ratio + 1))]
               if (s > 0._R8P) coarse%face_y_centroid(:, i, j, k) = sc / s
            enddo
         enddo
      enddo
      do k = 1, n(3) + 1
         do j = 1, n(2)
            do i = 1, n(1)
               s = 0._R8P ; sc = 0._R8P ; coarse%face_z_closed(i, j, k) = .true.
               do b2 = (j - 1) * ratio + 1, j * ratio
                  do b1 = (i - 1) * ratio + 1, i * ratio
                     a2 = (k - 1) * ratio + 1 ; fa = (x(b1 + 1) - x(b1)) * (y(b2 + 1) - y(b2)) * fine%face_z(b1, b2, a2)
                     s = s + fa ; sc = sc + fa * fine%face_z_centroid(:, b1, b2, a2)
                     coarse%face_z_closed(i, j, k) = coarse%face_z_closed(i, j, k) .and. fine%face_z_closed(b1, b2, a2)
                  enddo
               enddo
               fa = (x(i * ratio + 1) - x((i - 1) * ratio + 1)) * (y(j * ratio + 1) - y((j - 1) * ratio + 1))
               coarse%face_z(i, j, k) = s / fa
               coarse%face_z_centroid(:, i, j, k) = [0.5_R8P * (x(i * ratio + 1) + x((i - 1) * ratio + 1)), &
                                                     0.5_R8P * (y(j * ratio + 1) + y((j - 1) * ratio + 1)), z((k - 1) * ratio + 1)]
               if (s > 0._R8P) coarse%face_z_centroid(:, i, j, k) = sc / s
            enddo
         enddo
      enddo
      endsubroutine coarse_faces

      subroutine coarse_edges
      !< Coarse edge fractions: the fine edges along each coarse edge added up (length-weighted).
      integer(I4P) :: a1 !< Work.
      real(R8P)    :: s  !< Work.

      allocate(coarse%edge_x(n(1), n(2) + 1, n(3) + 1), coarse%edge_y(n(1) + 1, n(2), n(3) + 1), &
               coarse%edge_z(n(1) + 1, n(2) + 1, n(3)))
      do k = 1, n(3) + 1
         do j = 1, n(2) + 1
            do i = 1, n(1)
               s = 0._R8P
               do a1 = (i - 1) * ratio + 1, i * ratio
                  s = s + fine%edge_x(a1, (j - 1) * ratio + 1, (k - 1) * ratio + 1) * (x(a1 + 1) - x(a1))
               enddo
               coarse%edge_x(i, j, k) = s / (x(i * ratio + 1) - x((i - 1) * ratio + 1))
            enddo
         enddo
      enddo
      do k = 1, n(3) + 1
         do j = 1, n(2)
            do i = 1, n(1) + 1
               s = 0._R8P
               do a1 = (j - 1) * ratio + 1, j * ratio
                  s = s + fine%edge_y((i - 1) * ratio + 1, a1, (k - 1) * ratio + 1) * (y(a1 + 1) - y(a1))
               enddo
               coarse%edge_y(i, j, k) = s / (y(j * ratio + 1) - y((j - 1) * ratio + 1))
            enddo
         enddo
      enddo
      do k = 1, n(3)
         do j = 1, n(2) + 1
            do i = 1, n(1) + 1
               s = 0._R8P
               do a1 = (k - 1) * ratio + 1, k * ratio
                  s = s + fine%edge_z((i - 1) * ratio + 1, (j - 1) * ratio + 1, a1) * (z(a1 + 1) - z(a1))
               enddo
               coarse%edge_z(i, j, k) = s / (z(k * ratio + 1) - z((k - 1) * ratio + 1))
            enddo
         enddo
      enddo
      endsubroutine coarse_edges
   endsubroutine cut_coarsen_3d

   subroutine cut_coarsen_2d(fine, links, x, y, ratio, coarse, coarse_links, fine_to_coarse, status)
   !< 2-D twin of `cut_coarsen_3d`: areas, centroids, wetted lengths, centroids and normals, apertures (per piece too),
   !< node flags, moments and links (see the module notes).
   type(cut_cells_2d_t),      intent(in)  :: fine            !< Fine cut cells.
   type(cut_links_t),         intent(in)  :: links           !< Fine links (`cut_links_2d`).
   real(R8P),                 intent(in)  :: x(:), y(:)      !< Fine cell edges.
   integer(I4P),              intent(in)  :: ratio           !< Refinement ratio.
   type(cut_cells_2d_t),      intent(inout) :: coarse        !< Coarse cut cells (released from the device, reset).
   type(cut_links_t),         intent(out) :: coarse_links    !< Coarse links.
   integer(I4P), allocatable, intent(out) :: fine_to_coarse(:) !< Coarse piece of every fine piece.
   integer(I4P),              intent(out) :: status          !< 0 or 1.
   integer(I4P), allocatable :: par(:), root_piece(:), lp(:), lq(:), lf(:), e2(:,:), corder(:), cstart(:), fill(:)
   real(R8P),    allocatable :: acc(:,:), pm(:,:), pw(:,:)
   real(R8P)                 :: hf(2), hc(2), cf(2), ccn(2), ln
   real(R8P),    allocatable :: la(:)
   integer(I4P)              :: nf(2), n(2), p, q, k, ci(2), cj(2), ic(2), np, f, ax, nm, i, j, nl, c
   logical,      allocatable :: bad(:,:)

   call cut_cells_2d_reset(coarse)
   status = 0_I4P
   nf = [fine%nx, fine%ny]
   if (any(mod(nf, ratio) /= 0_I4P) .or. ratio < 1_I4P) then
      status = 1_I4P ; return
   endif
   n = nf / ratio ; coarse%nx = n(1) ; coarse%ny = n(2)
   allocate(par(fine%np)) ; par = [(p, p = 1, fine%np)]
   do p = 1, fine%np
      do k = links%start(p), links%start(p + 1) - 1
         q = links%piece(k)
         if (all((fine%piece_cell(:, p) - 1) / ratio == (fine%piece_cell(:, q) - 1) / ratio)) call unite(p, q)
      enddo
   enddo
   allocate(bad(n(1), n(2))) ; bad = .false.
   do j = 1, nf(2)
      do i = 1, nf(1)
         if (.not. fine%closed(i, j)) bad((i - 1) / ratio + 1, (j - 1) / ratio + 1) = .true.
      enddo
   enddo
   allocate(root_piece(fine%np), fine_to_coarse(fine%np), cstart(product(n) + 1), corder(fine%np), fill(product(n)))
   root_piece = 0_I4P ; fine_to_coarse = 0_I4P ; cstart = 0_I4P
   allocate(coarse%pieces(n(1), n(2)), coarse%first_piece(n(1), n(2))) ; coarse%pieces = 0_I4P
   do p = 1, fine%np
      ci = (fine%piece_cell(:, p) - 1) / ratio + 1 ; c = ci(1) + n(1) * (ci(2) - 1) ; cstart(c + 1) = cstart(c + 1) + 1_I4P
   enddo
   cstart(1) = 1_I4P
   do c = 1, product(n)
      cstart(c + 1) = cstart(c + 1) + cstart(c)
   enddo
   fill = 0_I4P
   do p = 1, fine%np
      ci = (fine%piece_cell(:, p) - 1) / ratio + 1 ; c = ci(1) + n(1) * (ci(2) - 1)
      corder(cstart(c) + fill(c)) = p ; fill(c) = fill(c) + 1_I4P
   enddo
   np = 0_I4P
   do c = 1, product(n)
      ci = [mod(c - 1, n(1)) + 1, (c - 1) / n(1) + 1]
      coarse%first_piece(ci(1), ci(2)) = np + 1_I4P
      if (bad(ci(1), ci(2))) cycle
      do k = cstart(c), cstart(c + 1) - 1
         p = corder(k) ; q = find(p)
         if (root_piece(q) == 0_I4P) then
            np = np + 1_I4P ; root_piece(q) = np ; coarse%pieces(ci(1), ci(2)) = coarse%pieces(ci(1), ci(2)) + 1_I4P
         endif
         fine_to_coarse(p) = root_piece(q)
      enddo
   enddo
   coarse%np = np
   ! sums: area (1), first moments (2:3), wet length (4), wet normal (5:6), four coarse edge lengths (8:11)
   allocate(acc(11, np), coarse%piece_cell(2, np)) ; acc = 0._R8P
   do p = 1, fine%np
      q = fine_to_coarse(p) ; if (q == 0_I4P) cycle
      ic = fine%piece_cell(:, p) ; ci = (ic - 1) / ratio + 1 ; coarse%piece_cell(:, q) = ci
      hf = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2))]
      acc(1, q) = acc(1, q) + fine%piece_area(p)
      acc(2:3, q) = acc(2:3, q) + fine%piece_area(p) * fine%piece_centroid(:, p)
      acc(4, q) = acc(4, q) + fine%piece_wet_length(p)
      acc(5:6, q) = acc(5:6, q) + fine%piece_wet_normal(:, p)
      do f = 1, 4 ! edges -u, +u (length hf(2)), -v, +v (length hf(1)) on the coarse boundary
         ax = (f + 1) / 2 ; ln = hf(3 - ax)
         if (mod(f, 2) == 1 .and. mod(ic(ax) - 1, ratio) == 0) acc(7 + f, q) = acc(7 + f, q) + fine%piece_aperture(f, p) * ln
         if (mod(f, 2) == 0 .and. mod(ic(ax), ratio) == 0) acc(7 + f, q) = acc(7 + f, q) + fine%piece_aperture(f, p) * ln
      enddo
   enddo
   allocate(coarse%piece_area(np), coarse%piece_centroid(2, np), coarse%piece_wet_length(np), coarse%piece_wet_normal(2, np), &
            coarse%piece_aperture(4, np))
   do q = 1, np
      call coarse_cell(coarse%piece_cell(:, q), cf, hc) ; ccn = cf + 0.5_R8P * hc
      coarse%piece_area(q) = acc(1, q)
      coarse%piece_centroid(:, q) = 0._R8P ; if (acc(1, q) > 0._R8P) coarse%piece_centroid(:, q) = acc(2:3, q) / acc(1, q)
      coarse%piece_wet_length(q) = acc(4, q) ; coarse%piece_wet_normal(:, q) = acc(5:6, q)
      coarse%piece_aperture(:, q) = acc(8:11, q) / [hc(2), hc(2), hc(1), hc(1)]
   enddo
   allocate(coarse%alpha(n(1), n(2)), coarse%centroid(2, n(1), n(2)), coarse%aperture(4, n(1), n(2)), coarse%closed(n(1), n(2)))
   allocate(coarse%wet_length(n(1), n(2)), coarse%wet_normal(2, n(1), n(2)), coarse%wet_centroid(2, n(1), n(2)))
   coarse%closed = .not. bad ; coarse%wet_length = 0._R8P ; coarse%wet_normal = 0._R8P
   do j = 1, n(2)
      do i = 1, n(1)
         call coarse_cell([i, j], cf, hc) ; ccn = cf + 0.5_R8P * hc
         coarse%alpha(i, j) = 0._R8P ; coarse%centroid(:, i, j) = 0._R8P ; coarse%aperture(:, i, j) = 0._R8P
         coarse%wet_centroid(:, i, j) = ccn ! the wetted centroid is not derived (no per-piece wetted moments in 2-D)
         if (bad(i, j)) then
            coarse%alpha(i, j) = -1._R8P ; coarse%centroid(:, i, j) = ccn ; cycle
         endif
         do q = coarse%first_piece(i, j), coarse%first_piece(i, j) + coarse%pieces(i, j) - 1
            coarse%alpha(i, j) = coarse%alpha(i, j) + acc(1, q) ; coarse%centroid(:, i, j) = coarse%centroid(:, i, j) + acc(2:3, q)
            coarse%aperture(:, i, j) = coarse%aperture(:, i, j) + coarse%piece_aperture(:, q)
            coarse%wet_length(i, j) = coarse%wet_length(i, j) + acc(4, q)
            coarse%wet_normal(:, i, j) = coarse%wet_normal(:, i, j) + acc(5:6, q)
         enddo
         if (coarse%alpha(i, j) > 0._R8P) then
            coarse%centroid(:, i, j) = coarse%centroid(:, i, j) / coarse%alpha(i, j)
         else
            coarse%centroid(:, i, j) = ccn
         endif
         coarse%alpha(i, j) = coarse%alpha(i, j) / product(hc)
      enddo
   enddo
   if (allocated(fine%node_in)) coarse%node_in = fine%node_in(1::ratio, 1::ratio)
   if (fine%moments_order > 0_I4P .and. allocated(fine%piece_moments)) then
      coarse%moments_order = fine%moments_order
      nm = mono_count(2_I4P, fine%moments_order) ; e2 = mono_exponents(2_I4P, fine%moments_order)
      allocate(pm(nm, np), pw(nm, np)) ; pm = 0._R8P ; pw = 0._R8P
      do p = 1, fine%np
         q = fine_to_coarse(p) ; if (q == 0_I4P) cycle
         ic = fine%piece_cell(:, p) ; hf = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2))]
         call coarse_cell(coarse%piece_cell(:, q), cf, hc)
         call shift_moments(fine%piece_moments(:, p), [x(ic(1)), y(ic(2))] + 0.5_R8P * hf, hf, cf + 0.5_R8P * hc, hc, e2, pm(:, q))
         call shift_moments(fine%piece_wall_moments(:, p), [x(ic(1)), y(ic(2))] + 0.5_R8P * hf, hf, cf + 0.5_R8P * hc, hc, e2, &
                            pw(:, q))
      enddo
      call move_alloc(pm, coarse%piece_moments) ; call move_alloc(pw, coarse%piece_wall_moments)
   endif
   allocate(lp(64), lq(64), lf(64), la(64)) ; nl = 0_I4P
   do p = 1, fine%np
      if (fine_to_coarse(p) == 0_I4P) cycle
      do k = links%start(p), links%start(p + 1) - 1
         q = links%piece(k)
         if (fine_to_coarse(q) == 0_I4P) cycle
         ci = (fine%piece_cell(:, p) - 1) / ratio + 1 ; cj = (fine%piece_cell(:, q) - 1) / ratio + 1
         if (all(ci == cj)) cycle
         if (nl == size(lp)) then
            lp = [lp, lp] ; lq = [lq, lq] ; lf = [lf, lf] ; la = [la, la]
         endif
         nl = nl + 1_I4P ; lp(nl) = fine_to_coarse(p) ; lq(nl) = fine_to_coarse(q) ; lf(nl) = links%face(k)
         la(nl) = links%area(k)
      enddo
   enddo
   coarse_links%np = np
   call build_csr(nl, lp, lq, lf, la, coarse_links)
   contains
      integer(I4P) function find(a) result(r)
      !< Root of a (path halving).
      integer(I4P), intent(in) :: a !< Piece.

      r = a
      do while (par(r) /= r)
         par(r) = par(par(r)) ; r = par(r)
      enddo
      endfunction find

      subroutine unite(a, b)
      !< Join the sets of a and b.
      integer(I4P), intent(in) :: a, b !< Pieces.
      integer(I4P)             :: ra, rb !< Roots.

      ra = find(a) ; rb = find(b)
      if (ra /= rb) par(max(ra, rb)) = min(ra, rb)
      endsubroutine unite

      subroutine coarse_cell(cc, lo, h)
      !< Low corner and size of the coarse cell cc.
      integer(I4P), intent(in)  :: cc(2)       !< Cell.
      real(R8P),    intent(out) :: lo(2), h(2) !< Corner, size.

      lo = [x((cc(1) - 1) * ratio + 1), y((cc(2) - 1) * ratio + 1)]
      h = [x(cc(1) * ratio + 1), y(cc(2) * ratio + 1)] - lo
      endsubroutine coarse_cell
   endsubroutine cut_coarsen_2d

   pure subroutine shift_moments(mf, cf, hf, cc, hc, e, mc)
   !< Add the moments mf of a fine cell (centre cf, size hf) to mc, re-centred and rescaled to the coarse cell (centre
   !< cc, size hc): X_c = alpha + beta X_f per axis, alpha = (cf - cc) / hc, beta = hf / hc, binomial expansion.
   real(R8P),    intent(in)    :: mf(:)         !< Fine moments.
   real(R8P),    intent(in)    :: cf(:), hf(:)  !< Fine centre, size.
   real(R8P),    intent(in)    :: cc(:), hc(:)  !< Coarse centre, size.
   integer(I4P), intent(in)    :: e(:,:)        !< Exponents [d, n].
   real(R8P),    intent(inout) :: mc(:)         !< Coarse moments.
   real(R8P)                   :: al(size(cf)), be(size(cf)), t !< Work.
   integer(I4P)                :: k, l, a       !< Work.

   al = (cf - cc) / hc ; be = hf / hc
   do k = 1, size(e, 2)
      do l = 1, size(e, 2)
         if (any(e(:, l) > e(:, k))) cycle
         t = mf(l)
         do a = 1, size(cf)
            t = t * binom(e(a, k), e(a, l)) * al(a) ** (e(a, k) - e(a, l)) * be(a) ** e(a, l)
         enddo
         mc(k) = mc(k) + t
      enddo
   enddo
   endsubroutine shift_moments

   pure real(R8P) function binom(n_, k_)
   !< Binomial coefficient.
   integer(I4P), intent(in) :: n_, k_ !< Arguments.
   integer(I4P)             :: j      !< Work.

   binom = 1._R8P
   do j = 1, k_
      binom = binom * (n_ - k_ + j) / j
   enddo
   endfunction binom
endmodule fossil_cut_coarsen
