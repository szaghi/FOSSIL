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
!<
!< State redistribution (issue #70 step 5; Berger & Giuliani 2021 §4, weighted factors of Giuliani et al. 2022 §3.1):
!< `srd_neighbourhoods` gives every small piece (volume below tau times its cell's) a neighbourhood grown through the
!< links inside the 3^d tile around its cell, then the 5^d tile: candidates ordered by the alignment of their cell
!< offset o with the piece's wall normal (n . o / |o|, largest first), ties by smaller |o|_1, lexicographic o, piece
!< index; they join until the members' volume reaches tau times the cell's. Zero-volume pieces carry no state: never
!< members, never traversed. Each neighbourhood depends on its piece only, so the result does not depend on the order
!< of the pieces. Outputs: CSR members (the piece first), overlap counts N^ov, weighted volumes and centroids, the
!< weighted-SRD factors eta^nb and eta^self, and a status per piece.
!<
!< Flux redistribution (issue #70 step 6; Natarajan et al. 2022 §2, Barrio Sanchez et al. 2024 App. A): for every cut
!< piece (volume below its cell's), `frd_neighbourhoods` gives the pieces of the 3^d block around its cell reached by a
!< monotone index path: a sequence of links, each moving one index one step away from the piece's cell and never back
!< (every component of the offset grows in magnitude, keeping its sign), through pieces with volume. Weights
!< kappa_j / sum kappa over the set, the piece included (kappa = volume / cell volume). A piece filling its cell has
!< itself only, weight 1; a zero-volume piece has none.
!<
!< Merge targets (issue #70 step 7). The pieces served are the small ones (volume below tau times the cell's) and the
!< ones in `mask` (fresh or dead pieces of a moving body). A piece is suitable as a target if it has volume, is not
!< small and is not in `exclude`. One target (cell merging, Falcovitz et al. 1997): the linked piece with volume that
!< maximises n . o / |o| (n the unit wall normal, o the cell offset), suitable if any is, ties by the larger link area,
!< then the lower index; an unsuitable target (small, or excluded: a dead piece whose neighbours die too) is followed to
!< the first suitable piece along the chain (`final_target`), a chain that loops or ends is reported. Normal-directed targets
!< (conservative mixing, Bezgin et al. 2023 §3.6): the pieces at the offsets o with o_a in {0, sign n_a}, o /= 0, reached
!< by a monotone path; several in one cell: the largest. Weights beta_o ~ (prod_{a: o_a /= 0} |n_a|)^(2 / |supp o|)
!< (after Hu et al. 2006, general knowledge: the cited weights are not in the literature base), normalised.
use fossil_cut_cells_2d, only : cut_cells_2d_t
use penf,                only : I4P, R8P

implicit none
private
public :: cut_links_t, cut_links_2d, build_csr
public :: srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, SRD_NO_VOLUME
public :: frd_t, frd_neighbourhoods
public :: merge_t, merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, MERGE_CYCLE, MERGE_ALONE

integer(I4P), parameter :: MERGE_NONE  = 0_I4P !< Not served (neither small nor in the mask).
integer(I4P), parameter :: MERGE_OK    = 1_I4P !< Target suitable.
integer(I4P), parameter :: MERGE_CHAIN = 2_I4P !< Target unsuitable; `final_target` is the first suitable one along the chain.
integer(I4P), parameter :: MERGE_CYCLE = 3_I4P !< The chain loops or ends before a suitable piece: `final_target` = 0.
integer(I4P), parameter :: MERGE_ALONE = 4_I4P !< No linked piece has volume.

integer(I4P), parameter :: SRD_REGULAR    = 0_I4P !< Not small: its neighbourhood is itself.
integer(I4P), parameter :: SRD_COMPLETE   = 1_I4P !< Small, neighbourhood reaches tau |cell|.
integer(I4P), parameter :: SRD_INCOMPLETE = 2_I4P !< Small, the 5^d tile did not reach tau |cell| (all of it is kept).
integer(I4P), parameter :: SRD_NO_VOLUME  = 3_I4P !< Zero volume: no state, no neighbourhood.

type :: cut_links_t
   !< Links between the pieces of a block (see the module notes).
   integer(I4P)              :: np = 0_I4P     !< Pieces.
   integer(I4P), allocatable :: start(:)       !< Links of piece p: start(p) .. start(p + 1) - 1 [np + 1].
   integer(I4P), allocatable :: piece(:)       !< Piece across.
   integer(I4P), allocatable :: face(:)        !< Face of p's cell (1-6; 1-4 in 2-D).
   real(R8P),    allocatable :: area(:)        !< Shared fluid area (length in 2-D).
   integer(I4P)              :: unmatched = 0_I4P !< Fluid regions of a face that could not be matched to a piece.
endtype cut_links_t

type :: srd_t
   !< State-redistribution neighbourhoods of the pieces of a block (see the module notes).
   integer(I4P)              :: np = 0_I4P  !< Pieces.
   real(R8P)                 :: tau = 0.5_R8P !< Threshold (fraction of the cell volume).
   integer(I4P), allocatable :: start(:)    !< Members of N_p: member(start(p) .. start(p + 1) - 1), p first [np + 1].
   integer(I4P), allocatable :: member(:)   !< Members.
   integer(I4P), allocatable :: nov(:)      !< Overlap count N^ov_j: the neighbourhoods containing j [np].
   integer(I4P), allocatable :: status(:)   !< SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE or SRD_NO_VOLUME [np].
   real(R8P),    allocatable :: vhat(:)     !< Weighted volume sum_{j in N_p} |I_j| / N^ov_j [np].
   real(R8P),    allocatable :: xhat(:,:)   !< Weighted centroid [d, np].
   real(R8P),    allocatable :: eta_nb(:)   !< Weighted SRD: (tau |cell| - |I_p|) / sum_{j in N_p, j /= p} |I_j|, 0 if regular.
   real(R8P),    allocatable :: eta_self(:) !< Weighted SRD: 1 - (1 / N^ov_p) sum_{k /= p, p in N_k} eta^nb_k.
endtype srd_t

type :: frd_t
   !< Flux-redistribution neighbourhoods of the pieces of a block (see the module notes).
   integer(I4P)              :: np = 0_I4P !< Pieces.
   integer(I4P), allocatable :: start(:)   !< Members of piece p: member(start(p) .. start(p + 1) - 1), p first [np + 1].
   integer(I4P), allocatable :: member(:)  !< Members.
   real(R8P),    allocatable :: weight(:)  !< kappa_j / sum kappa over the set.
endtype frd_t

type :: merge_t
   !< Merge targets of the pieces of a block (see the module notes).
   integer(I4P)              :: np = 0_I4P !< Pieces.
   integer(I4P), allocatable :: status(:)  !< MERGE_NONE, _OK, _CHAIN, _CYCLE or _ALONE [np].
   integer(I4P), allocatable :: target(:)  !< The single target (0 if none) [np].
   integer(I4P), allocatable :: final_target(:) !< The first target along the chain that is not small (0 if none) [np].
   integer(I4P), allocatable :: mix_start(:) !< Normal-directed targets of p: mix_piece(mix_start(p) .. (p + 1) - 1).
   integer(I4P), allocatable :: mix_piece(:) !< Targets.
   real(R8P),    allocatable :: mix_weight(:) !< Weights beta_o, normalised over the targets found.
endtype merge_t

contains
   subroutine merge_targets(links, piece_cell, piece_volume, wet_normal, x, y, mt, z, tau, mask, exclude)
   !< Merge targets (see the module notes). 2-D without z, 3-D with it.
   type(cut_links_t), intent(in)           :: links            !< Links.
   integer(I4P),      intent(in)           :: piece_cell(:,:)  !< Cell of each piece [d, np].
   real(R8P),         intent(in)           :: piece_volume(:)  !< Volume (area in 2-D) of each piece [np].
   real(R8P),         intent(in)           :: wet_normal(:,:)  !< Wall normal integral (into the fluid) [d, np].
   real(R8P),         intent(in)           :: x(:), y(:)       !< Cell edges.
   type(merge_t),     intent(out)          :: mt               !< Targets.
   real(R8P),         intent(in), optional :: z(:)             !< Cell edges along z (3-D).
   real(R8P),         intent(in), optional :: tau              !< Small-piece threshold, default 1/2.
   logical,           intent(in), optional :: mask(:)          !< Pieces served whatever their size [np].
   logical,           intent(in), optional :: exclude(:)       !< Pieces that may not be targets [np].
   logical,      allocatable :: small(:), served(:), ok_t(:)  !< Small, served, suitable as a target.
   integer(I4P), allocatable :: mark(:), cand(:), mp(:)      !< Search stamps, queue, mixing targets.
   real(R8P),    allocatable :: mw(:)                        !< Mixing weights.
   integer(I4P) :: np, d, p, k, q, best, nc, head, kk, a, o(3), oq(3), stp(3), guard, nm, ns, sg(3), t !< Work.
   real(R8P)    :: cv, th, nh(3), sc, bsc, bar, w   !< Work.
   logical      :: bsmall                           !< Best so far is small.

   np = size(piece_volume, kind=I4P) ; d = size(piece_cell, 1, kind=I4P)
   th = 0.5_R8P ; if (present(tau)) th = tau
   mt%np = np
   allocate(mt%status(np), mt%target(np), mt%final_target(np), mt%mix_start(np + 1), small(np), served(np), ok_t(np), mark(np))
   allocate(cand(64), mp(max(16_I4P, np)), mw(max(16_I4P, np)))
   mark = 0_I4P ; nm = 0_I4P
   do p = 1, np
      cv = (x(piece_cell(1, p) + 1) - x(piece_cell(1, p))) * (y(piece_cell(2, p) + 1) - y(piece_cell(2, p)))
      if (d == 3_I4P .and. present(z)) cv = cv * (z(piece_cell(3, p) + 1) - z(piece_cell(3, p)))
      small(p) = piece_volume(p) > 0._R8P .and. piece_volume(p) < th * cv
      served(p) = small(p)
      if (present(mask)) served(p) = served(p) .or. (mask(p) .and. piece_volume(p) > 0._R8P)
      ok_t(p) = piece_volume(p) > 0._R8P .and. .not. small(p)
      if (present(exclude)) ok_t(p) = ok_t(p) .and. .not. exclude(p)
   enddo
   ! single targets
   do p = 1, np
      mt%target(p) = 0_I4P ; mt%status(p) = MERGE_NONE
      if (.not. served(p)) cycle
      call unit_normal(p)
      best = 0_I4P ; bsc = -huge(1._R8P) ; bar = 0._R8P ; bsmall = .true.
      do k = links%start(p), links%start(p + 1) - 1
         q = links%piece(k)
         if (piece_volume(q) <= 0._R8P) cycle
         o = 0_I4P ; o(1:d) = piece_cell(:, q) - piece_cell(:, p)
         sc = dot_product(nh, real(o, R8P))
         if (better(.not. ok_t(q), sc, links%area(k), q)) then
            best = q ; bsc = sc ; bar = links%area(k) ; bsmall = .not. ok_t(q)
         endif
      enddo
      mt%target(p) = best
      if (best == 0_I4P) then
         mt%status(p) = MERGE_ALONE
      else
         mt%status(p) = merge(MERGE_CHAIN, MERGE_OK, .not. ok_t(best))
      endif
   enddo
   do p = 1, np ! follow the chains of unsuitable targets
      mt%final_target(p) = mt%target(p)
      if (mt%status(p) /= MERGE_CHAIN) cycle
      q = mt%target(p) ; guard = 0_I4P
      do while (q /= 0_I4P .and. guard <= np)
         if (ok_t(q)) exit
         q = mt%target(q) ; guard = guard + 1_I4P
      enddo
      if (q == 0_I4P .or. .not. ok_t(q)) then
         mt%status(p) = MERGE_CYCLE ; mt%final_target(p) = 0_I4P
      else
         mt%final_target(p) = q
      endif
   enddo
   ! normal-directed targets
   mt%mix_start(1) = 1_I4P
   do p = 1, np
      mt%mix_start(p + 1) = mt%mix_start(p)
      if (.not. served(p)) cycle
      call unit_normal(p)
      sg = 0_I4P
      do a = 1, d
         if (nh(a) > 0._R8P) sg(a) = 1_I4P
         if (nh(a) < 0._R8P) sg(a) = -1_I4P
      enddo
      ! monotone search inside the cells with offsets o_a in {0, sg_a}
      nc = 1_I4P ; cand(1) = p ; mark(p) = p ; head = 0_I4P
      do while (head < nc)
         head = head + 1_I4P ; q = cand(head)
         o = 0_I4P ; o(1:d) = piece_cell(:, q) - piece_cell(:, p)
         do k = links%start(q), links%start(q + 1) - 1
            kk = links%piece(k)
            if (mark(kk) == p .or. piece_volume(kk) <= 0._R8P) cycle
            oq = 0_I4P ; oq(1:d) = piece_cell(:, kk) - piece_cell(:, p)
            stp = oq - o ; a = maxloc(abs(stp), 1)
            if (o(a) /= 0_I4P .or. oq(a) /= sg(a) .or. sg(a) == 0_I4P) cycle
            mark(kk) = p
            if (nc == size(cand)) cand = [cand, cand]
            nc = nc + 1_I4P ; cand(nc) = kk
         enddo
      enddo
      ns = 0_I4P ; w = 0._R8P
      do k = 2, nc ! one target per offset: the largest piece reached there
         q = cand(k)
         if (.not. ok_t(q)) cycle
         t = 0_I4P
         do kk = 2, nc
            if (kk == k .or. .not. ok_t(cand(kk))) cycle
            if (any(piece_cell(:, cand(kk)) /= piece_cell(:, q))) cycle
            if (piece_volume(cand(kk)) > piece_volume(q) .or. &
                (piece_volume(cand(kk)) == piece_volume(q) .and. cand(kk) < q)) t = 1_I4P
         enddo
         if (t == 1_I4P) cycle
         o = 0_I4P ; o(1:d) = piece_cell(:, q) - piece_cell(:, p)
         if (nm == size(mp)) then
            mp = [mp, mp] ; mw = [mw, mw]
         endif
         nm = nm + 1_I4P ; ns = ns + 1_I4P ; mp(nm) = q
         mw(nm) = product(abs(nh(1:d)), mask=o(1:d) /= 0_I4P) ** (2._R8P / count(o(1:d) /= 0_I4P))
         w = w + mw(nm)
      enddo
      do k = 1, nc
         mark(cand(k)) = 0_I4P
      enddo
      if (w > 0._R8P) mw(nm - ns + 1:nm) = mw(nm - ns + 1:nm) / w
      mt%mix_start(p + 1) = mt%mix_start(p) + ns
   enddo
   mt%mix_piece = mp(1:nm) ; mt%mix_weight = mw(1:nm)
   contains
      subroutine unit_normal(pp)
      !< nh: the unit wall normal of piece pp (zero if it has no wall).
      integer(I4P), intent(in) :: pp !< Piece.

      nh = 0._R8P ; nh(1:d) = wet_normal(:, pp)
      if (norm2(nh) > 0._R8P) nh = nh / norm2(nh)
      endsubroutine unit_normal

      logical function better(sq, s_, ar_, q_)
      !< Candidate (unsuitable sq, score s_, link area ar_, index q_) beats the best so far.
      logical,      intent(in) :: sq   !< Candidate is unsuitable.
      real(R8P),    intent(in) :: s_   !< Score.
      real(R8P),    intent(in) :: ar_  !< Link area.
      integer(I4P), intent(in) :: q_   !< Index.

      if (best == 0_I4P) then
         better = .true.
      elseif (sq .neqv. bsmall) then
         better = bsmall ! a suitable target beats an unsuitable one
      elseif (s_ /= bsc) then
         better = s_ > bsc
      elseif (ar_ /= bar) then
         better = ar_ > bar
      else
         better = q_ < best
      endif
      endfunction better
   endsubroutine merge_targets

   subroutine frd_neighbourhoods(links, piece_cell, piece_volume, x, y, frd, z)
   !< Flux-redistribution neighbourhoods (see the module notes). 2-D without z, 3-D with it.
   type(cut_links_t), intent(in)           :: links           !< Links.
   integer(I4P),      intent(in)           :: piece_cell(:,:) !< Cell of each piece [d, np].
   real(R8P),         intent(in)           :: piece_volume(:) !< Volume (area in 2-D) of each piece [np].
   real(R8P),         intent(in)           :: x(:), y(:)      !< Cell edges.
   type(frd_t),       intent(out)          :: frd             !< Neighbourhoods.
   real(R8P),         intent(in), optional :: z(:)            !< Cell edges along z (3-D).
   integer(I4P), allocatable :: mark(:), cand(:), mem(:) !< Visited stamps, queue, members.
   real(R8P),    allocatable :: kap(:), wt(:)            !< Volume fractions, weights.
   integer(I4P) :: np, d, p, nc, head, k, kk, q, nmem, a, o(3), oq(3), stp(3) !< Work.
   real(R8P)    :: cv, ks                                 !< Cell volume, sum of kappa.

   np = size(piece_volume, kind=I4P) ; d = size(piece_cell, 1, kind=I4P)
   frd%np = np
   allocate(frd%start(np + 1), kap(np), mark(np), cand(64), mem(max(16_I4P, 2_I4P * np)), wt(max(16_I4P, 2_I4P * np)))
   mark = 0_I4P ; nmem = 0_I4P
   do p = 1, np
      cv = (x(piece_cell(1, p) + 1) - x(piece_cell(1, p))) * (y(piece_cell(2, p) + 1) - y(piece_cell(2, p)))
      if (d == 3_I4P .and. present(z)) cv = cv * (z(piece_cell(3, p) + 1) - z(piece_cell(3, p)))
      kap(p) = piece_volume(p) / cv
   enddo
   frd%start(1) = 1_I4P
   do p = 1, np
      if (piece_volume(p) <= 0._R8P) then
         frd%start(p + 1) = frd%start(p) ; cycle
      endif
      nc = 1_I4P ; cand(1) = p ; mark(p) = p ; head = 0_I4P
      if (kap(p) < 1._R8P) then
         do while (head < nc) ! breadth first; every step adds one to |o|_1, so a piece is reached by its shortest paths
            head = head + 1_I4P ; q = cand(head)
            o = 0_I4P ; o(1:d) = piece_cell(:, q) - piece_cell(:, p)
            do k = links%start(q), links%start(q + 1) - 1
               kk = links%piece(k)
               if (mark(kk) == p .or. piece_volume(kk) <= 0._R8P) cycle
               oq = 0_I4P ; oq(1:d) = piece_cell(:, kk) - piece_cell(:, p)
               stp = oq - o ; a = maxloc(abs(stp), 1)
               if (abs(oq(a)) /= abs(o(a)) + 1_I4P .or. abs(oq(a)) > 1_I4P) cycle ! away from the centre, inside the block
               mark(kk) = p
               if (nc == size(cand)) cand = [cand, cand]
               nc = nc + 1_I4P ; cand(nc) = kk
            enddo
         enddo
      endif
      ks = sum(kap(cand(1:nc)))
      do k = 1, nc
         mark(cand(k)) = 0_I4P
         if (nmem == size(mem)) then
            mem = [mem, mem] ; wt = [wt, wt]
         endif
         nmem = nmem + 1_I4P ; mem(nmem) = cand(k) ; wt(nmem) = kap(cand(k)) / ks
      enddo
      frd%start(p + 1) = frd%start(p) + nc
   enddo
   frd%member = mem(1:nmem) ; frd%weight = wt(1:nmem)
   endsubroutine frd_neighbourhoods

   subroutine srd_neighbourhoods(links, piece_cell, piece_volume, piece_centroid, wet_normal, x, y, srd, z, tau)
   !< State-redistribution neighbourhoods (see the module notes). 2-D without z (piece_cell [2, np], ...), 3-D with it.
   type(cut_links_t), intent(in)           :: links              !< Links.
   integer(I4P),      intent(in)           :: piece_cell(:,:)    !< Cell of each piece [d, np].
   real(R8P),         intent(in)           :: piece_volume(:)    !< Volume (area in 2-D) of each piece [np].
   real(R8P),         intent(in)           :: piece_centroid(:,:) !< Centroid [d, np].
   real(R8P),         intent(in)           :: wet_normal(:,:)    !< Wall normal integral (into the fluid) [d, np].
   real(R8P),         intent(in)           :: x(:), y(:)         !< Cell edges.
   type(srd_t),       intent(out)          :: srd                !< Neighbourhoods.
   real(R8P),         intent(in), optional :: z(:)               !< Cell edges along z (3-D).
   real(R8P),         intent(in), optional :: tau                !< Threshold, default 1/2.
   integer(I4P), allocatable :: mark(:), cand(:), off(:,:), ord(:), mem(:), nm(:) !< Work.
   real(R8P),    allocatable :: sc(:), cv(:)    !< Scores, cell volumes.
   integer(I4P) :: np, d, p, nc, k, r, cap, nmem, j, m, head !< Work.
   real(R8P)    :: nh(3), vol, need, w      !< Unit normal, volume, target, weight.

   np = size(piece_volume, kind=I4P) ; d = size(piece_cell, 1, kind=I4P)
   srd%np = np ; if (present(tau)) srd%tau = tau
   allocate(srd%start(np + 1), srd%nov(np), srd%status(np), srd%vhat(np), srd%xhat(d, np), srd%eta_nb(np), srd%eta_self(np))
   allocate(cv(np), mark(np), nm(np)) ; mark = 0_I4P
   do p = 1, np
      cv(p) = (x(piece_cell(1, p) + 1) - x(piece_cell(1, p))) * (y(piece_cell(2, p) + 1) - y(piece_cell(2, p)))
      if (d == 3_I4P .and. present(z)) cv(p) = cv(p) * (z(piece_cell(3, p) + 1) - z(piece_cell(3, p)))
   enddo
   cap = max(16_I4P, 2_I4P * np) ; allocate(mem(cap)) ; nmem = 0_I4P
   allocate(cand(64), off(3, 64), sc(64), ord(64))
   srd%start(1) = 1_I4P
   do p = 1, np
      if (piece_volume(p) <= 0._R8P) then
         srd%status(p) = SRD_NO_VOLUME ; nm(p) = 0_I4P
      elseif (piece_volume(p) >= srd%tau * cv(p)) then
         srd%status(p) = SRD_REGULAR ; call push_member(p) ; nm(p) = 1_I4P
      else
         nh = 0._R8P ; nh(1:d) = wet_normal(:, p)
         if (norm2(nh) > 0._R8P) nh = nh / norm2(nh)
         need = srd%tau * cv(p)
         do r = 1, 2 ! the 3^d tile, then the 5^d tile
            call gather(p, r)
            call order_candidates
            vol = piece_volume(p) ; m = 1_I4P
            do k = 1, nc
               if (vol >= need) exit
               vol = vol + piece_volume(cand(ord(k))) ; m = m + 1_I4P
            enddo
            if (vol >= need) exit
         enddo
         srd%status(p) = merge(SRD_COMPLETE, SRD_INCOMPLETE, vol >= need)
         call push_member(p)
         do k = 1, m - 1_I4P
            call push_member(cand(ord(k)))
         enddo
         nm(p) = m
      endif
      srd%start(p + 1) = srd%start(p) + nm(p)
   enddo
   srd%member = mem(1:nmem)
   ! overlap counts, weighted volumes and centroids
   srd%nov = 0_I4P
   do k = 1, nmem
      srd%nov(srd%member(k)) = srd%nov(srd%member(k)) + 1_I4P
   enddo
   do p = 1, np
      srd%vhat(p) = 0._R8P ; srd%xhat(:, p) = 0._R8P
      do k = srd%start(p), srd%start(p + 1) - 1
         j = srd%member(k) ; w = piece_volume(j) / srd%nov(j)
         srd%vhat(p) = srd%vhat(p) + w ; srd%xhat(:, p) = srd%xhat(:, p) + w * piece_centroid(:, j)
      enddo
      if (srd%vhat(p) > 0._R8P) then
         srd%xhat(:, p) = srd%xhat(:, p) / srd%vhat(p)
      else
         srd%xhat(:, p) = piece_centroid(:, p)
      endif
   enddo
   ! weighted state redistribution
   do p = 1, np
      srd%eta_nb(p) = 0._R8P
      if (srd%status(p) /= SRD_COMPLETE .and. srd%status(p) /= SRD_INCOMPLETE) cycle
      w = sum(piece_volume(srd%member(srd%start(p) + 1:srd%start(p + 1) - 1)))
      if (w > 0._R8P) srd%eta_nb(p) = (srd%tau * cv(p) - piece_volume(p)) / w
   enddo
   srd%eta_self = 0._R8P
   do p = 1, np
      do k = srd%start(p) + 1, srd%start(p + 1) - 1
         j = srd%member(k) ; srd%eta_self(j) = srd%eta_self(j) + srd%eta_nb(p)
      enddo
   enddo
   do p = 1, np
      if (srd%nov(p) > 0_I4P) then
         srd%eta_self(p) = 1._R8P - srd%eta_self(p) / srd%nov(p)
      else
         srd%eta_self(p) = 0._R8P
      endif
   enddo
   contains
      subroutine push_member(j_)
      !< Append a member.
      integer(I4P), intent(in) :: j_ !< Piece.

      if (nmem == size(mem)) mem = [mem, mem]
      nmem = nmem + 1_I4P ; mem(nmem) = j_
      endsubroutine push_member

      subroutine gather(p0, rr)
      !< Candidates: the pieces with volume reached from p0 through links inside the tile of radius rr (p0 excluded).
      integer(I4P), intent(in) :: p0, rr !< Piece, tile radius.
      integer(I4P) :: kk                  !< Work.

      nc = 0_I4P ; head = 0_I4P
      mark(p0) = p0 ! visited stamp: the centre piece
      call visit_from(p0, p0, rr)
      do while (head < nc)
         head = head + 1_I4P
         call visit_from(cand(head), p0, rr)
      enddo
      do kk = 1, nc ! clear the stamps
         mark(cand(kk)) = 0_I4P
      enddo
      mark(p0) = 0_I4P
      endsubroutine gather

      subroutine visit_from(pf, p0, rr)
      !< Push the unvisited linked pieces of pf with volume inside the tile of radius rr around p0's cell.
      integer(I4P), intent(in) :: pf, p0, rr !< Piece, centre piece, tile radius.
      integer(I4P) :: kk, qq, o(3)            !< Work.

      do kk = links%start(pf), links%start(pf + 1) - 1
         qq = links%piece(kk)
         if (mark(qq) == p0 .or. piece_volume(qq) <= 0._R8P) cycle
         o = 0_I4P ; o(1:d) = piece_cell(:, qq) - piece_cell(:, p0)
         if (maxval(abs(o)) > rr) cycle
         mark(qq) = p0
         if (nc == size(cand)) call grow_candidates
         nc = nc + 1_I4P ; cand(nc) = qq ; off(:, nc) = o
      enddo
      endsubroutine visit_from

      subroutine grow_candidates
      !< Double the candidate buffers.
      integer(I4P), allocatable :: io(:,:) !< Work.

      allocate(io(3, 2 * size(cand))) ; io(:, 1:size(cand)) = off ; call move_alloc(io, off)
      cand = [cand, cand] ; sc = [sc, sc] ; ord = [ord, ord]
      endsubroutine grow_candidates

      subroutine order_candidates
      !< Order the candidates: n . o / |o| descending, then |o|_1, lexicographic o, piece index (insertion sort).
      integer(I4P) :: kk, ll, t !< Work.

      do kk = 1, nc
         sc(kk) = dot_product(nh, real(off(:, kk), R8P)) / norm2(real(off(:, kk), R8P))
         ord(kk) = kk
      enddo
      do kk = 2, nc
         t = ord(kk) ; ll = kk - 1_I4P
         do while (ll >= 1_I4P)
            if (.not. before(t, ord(ll))) exit
            ord(ll + 1) = ord(ll) ; ll = ll - 1_I4P
         enddo
         ord(ll + 1) = t
      enddo
      endsubroutine order_candidates

      logical function before(a, b)
      !< Candidate a goes before candidate b.
      integer(I4P), intent(in) :: a, b !< Candidates.
      integer(I4P)             :: aa   !< Work.

      if (sc(a) /= sc(b)) then
         before = sc(a) > sc(b) ; return
      endif
      if (sum(abs(off(:, a))) /= sum(abs(off(:, b)))) then
         before = sum(abs(off(:, a))) < sum(abs(off(:, b))) ; return
      endif
      do aa = 1, 3
         if (off(aa, a) /= off(aa, b)) then
            before = off(aa, a) < off(aa, b) ; return
         endif
      enddo
      before = cand(a) < cand(b)
      endfunction before
   endsubroutine srd_neighbourhoods

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
