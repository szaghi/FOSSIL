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
use, intrinsic :: iso_c_binding,   only : c_intptr_t, c_loc
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use fossil_cut_cells_2d, only : cut_cells_2d_t, chunk_offsets_dev
use penf,                only : I4P, R8P

implicit none
private
public :: cut_links_t, cut_links_2d, build_csr, cut_links_2d_device, cut_links_reset
! the link kernels of the 3-D blocks (issue #100 step 5), for `fossil_cut_cells_3d` (not re-exported by `fossil`)
public :: links_face_id, links_face_bases, links_piece_count_3d, links_piece_fill_3d, links_scan_dev, links_take
public :: srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, SRD_NO_VOLUME
public :: frd_t, frd_neighbourhoods
public :: merge_t, merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, MERGE_CYCLE, MERGE_ALONE

integer(I4P), parameter :: MERGE_NONE  = 0_I4P !< Not served (neither small nor in the mask).
integer(I4P), parameter :: MERGE_OK    = 1_I4P !< Target suitable.
integer(I4P), parameter :: MERGE_CHAIN = 2_I4P !< Target unsuitable; `final_target` is the first suitable one along the chain.
integer(I4P), parameter :: MERGE_CYCLE = 3_I4P !< The chain loops or ends before a suitable piece: `final_target` = 0.
integer(I4P), parameter :: MERGE_ALONE = 4_I4P !< No linked piece has volume.

integer(I4P), parameter :: LINKS_SCAN_CHUNK = 64_I4P !< Entries per chunk of the device prefix sums (`links_scan_dev`).

integer(I4P), parameter :: SRD_REGULAR    = 0_I4P !< Not small: its neighbourhood is itself.
integer(I4P), parameter :: SRD_COMPLETE   = 1_I4P !< Small, neighbourhood reaches tau |cell|.
integer(I4P), parameter :: SRD_INCOMPLETE = 2_I4P !< Small, the 5^d tile did not reach tau |cell| (all of it is kept).
integer(I4P), parameter :: SRD_NO_VOLUME  = 3_I4P !< Zero volume: no state, no neighbourhood.

type :: cut_links_t
   !< Links between the pieces of a block (see the module notes).
   !<
   !< Resident on the device (issue #100 step 5, the model of `cut_cells_2d_t`, step 6a). `cut_links_2d_device` and
   !< `cut_cells_3d_device(..., links=, resident=.true.)` leave start, piece, face and area mapped (OpenACC unstructured
   !< data) and never write their host copies: until `update_host` they hold garbage (np and unmatched are set).
   !< `on_device` says so. A device consumer passes the member arrays to its stage routine (the #20 rule: member arrays,
   !< never the object). `enter_device` maps a host result, `update_host` brings the device values back, `exit_device`
   !< releases the device copies (call it before resident links go out of scope or are deallocated). Every routine that
   !< writes links releases them first (`cut_links_reset`); a routine that reads them on the host (SRD, FRD, merge targets,
   !< the coarsening) needs `update_host` first. A copy `b = a` of resident links is not resident: the type-bound
   !< procedures recognise it (`device_key`, the host address of start when mapped). Built without OpenACC the
   !< procedures only set the flag.
   integer(I4P)              :: np = 0_I4P     !< Pieces.
   integer(I4P), allocatable :: start(:)       !< Links of piece p: start(p) .. start(p + 1) - 1 [np + 1].
   integer(I4P), allocatable :: piece(:)       !< Piece across.
   integer(I4P), allocatable :: face(:)        !< Face of p's cell (1-6; 1-4 in 2-D).
   real(R8P),    allocatable :: area(:)        !< Shared fluid area (length in 2-D).
   integer(I4P)              :: unmatched = 0_I4P !< Fluid regions of a face that could not be matched to a piece.
   logical                   :: on_device = .false. !< The arrays are mapped on the device (issue #100 step 5).
   integer(c_intptr_t)       :: device_key = 0  !< Host address of start when mapped (a copy has another one).
   contains
      procedure, pass(self) :: enter_device => links_enter_device !< Map the arrays to the device (copyin).
      procedure, pass(self) :: update_host  => links_update_host  !< Copy the mapped arrays back to the host.
      procedure, pass(self) :: exit_device  => links_exit_device  !< Release the device copies.
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
   !<
   !< Flat arrays (issue #100 step 5): the fluid edge intervals by edge (a counting sort in interval order, then a stable
   !< insertion sort by start within each edge, `links_edge_sort`: the order of the former global stable sort by (edge,
   !< start)); per edge the pairs of overlapping intervals of distinct pieces, counted (`links_edge_pairs`), a prefix
   !< sum, then written in the (interval, later interval) order of the edge (`links_edge_pairs_fill`); the CSR per
   !< piece (`links_piece_count_2d`, prefix sum, `links_piece_fill_2d`), with no global sort.
   !<
   !< Same bits as the sorted CSR of `build_csr` (which this replaced): build_csr stable-sorted the directed links (p,
   !< face f, q, length) by (p, f, q) and summed the lengths of equal keys in their input order. Face f of p's cell is
   !< one edge of the grid, and a pair recorded on an edge gave p one directed link (p is either interval of the pair),
   !< with f the face of p's cell that is that edge: the entries of a key are the pairs of one edge holding p and q,
   !< in the pair order of that edge, which is the slot order here; `links_piece_fill_2d` sums them in that order (from
   !< 0, one add at a time, as SUM does) and emits the keys of p by face 1..4, then q increasing, the sort order.
   type(cut_cells_2d_t), intent(in)    :: cut   !< Cut cells.
   type(cut_links_t),    intent(inout) :: links !< Links (released from the device and reset first).
   integer(I4P), allocatable :: owner(:)          !< Piece of each interval.
   integer(I4P), allocatable :: eoff(:), eord(:)  !< Intervals of edge e: eord(eoff(e) + 1 .. eoff(e + 1)).
   integer(I4P), allocatable :: fill(:)           !< Work.
   integer(I4P), allocatable :: poff(:)           !< Pairs of edge e: slots poff(e) + 1 .. poff(e + 1).
   integer(I4P), allocatable :: ppa(:), ppb(:)    !< Pieces of each pair (earlier, later interval in the edge order).
   real(R8P),    allocatable :: plen(:)           !< Shared length of each pair.
   integer(I4P)              :: n, ne, e, k, m, p, nk !< Work.

   call cut_links_reset(links)
   links%np = cut%np
   n = size(cut%piece_iv_edge, kind=I4P)
   ne = (cut%nx + 1_I4P) * cut%ny + cut%nx * (cut%ny + 1_I4P)
   allocate(owner(n))
   do p = 1, cut%np
      owner(cut%piece_iv_start(p):cut%piece_iv_start(p + 1) - 1) = p
   enddo
   ! intervals by edge, in interval order, then by start within the edge (stable)
   allocate(eoff(ne + 1), eord(n), fill(ne)) ; eoff = 0_I4P ; fill = 0_I4P
   do k = 1, n
      eoff(cut%piece_iv_edge(k) + 1) = eoff(cut%piece_iv_edge(k) + 1) + 1_I4P
   enddo
   do e = 1, ne
      eoff(e + 1) = eoff(e + 1) + eoff(e)
   enddo
   do k = 1, n
      e = cut%piece_iv_edge(k) ; fill(e) = fill(e) + 1_I4P ; eord(eoff(e) + fill(e)) = k
   enddo
   do e = 1, ne
      m = eoff(e + 1) - eoff(e)
      if (m > 1_I4P) call links_edge_sort(m, eord(eoff(e) + 1:eoff(e + 1)), n, cut%piece_iv_lo)
   enddo
   ! pairs per edge: counts, prefix sum, fill
   allocate(poff(ne + 1)) ; poff = 0_I4P
   do e = 1, ne
      m = eoff(e + 1) - eoff(e)
      if (m > 1_I4P) poff(e + 1) = links_edge_pairs(m, eord(eoff(e) + 1:eoff(e + 1)), n, owner, cut%piece_iv_lo, &
                                                    cut%piece_iv_hi)
   enddo
   do e = 1, ne
      poff(e + 1) = poff(e + 1) + poff(e)
   enddo
   allocate(ppa(poff(ne + 1)), ppb(poff(ne + 1)), plen(poff(ne + 1)))
   do e = 1, ne
      if (poff(e + 1) == poff(e)) cycle
      m = eoff(e + 1) - eoff(e)
      call links_edge_pairs_fill(m, eord(eoff(e) + 1:eoff(e + 1)), n, owner, cut%piece_iv_lo, cut%piece_iv_hi, &
                                 poff(e + 1) - poff(e), ppa(poff(e) + 1:poff(e + 1)), ppb(poff(e) + 1:poff(e + 1)), &
                                 plen(poff(e) + 1:poff(e + 1)))
   enddo
   ! CSR per piece
   allocate(links%start(cut%np + 1)) ; links%start = 0_I4P
   do p = 1, cut%np
      links%start(p + 1) = links_piece_count_2d(p, cut%piece_cell(1, p), cut%piece_cell(2, p), cut%nx, cut%ny, ne, poff, &
                                                size(ppa, kind=I4P), ppa, ppb)
   enddo
   links%start(1) = 1_I4P
   do p = 1, cut%np
      links%start(p + 1) = links%start(p + 1) + links%start(p)
   enddo
   nk = links%start(cut%np + 1) - 1_I4P
   allocate(links%piece(nk), links%face(nk), links%area(nk))
   do p = 1, cut%np
      call links_piece_fill_2d(p, cut%piece_cell(1, p), cut%piece_cell(2, p), cut%nx, cut%ny, ne, poff, &
                               size(ppa, kind=I4P), ppa, ppb, plen, links%start(p) - 1_I4P, nk, links%piece, links%face, &
                               links%area)
   enddo
   endsubroutine cut_links_2d

   subroutine cut_links_2d_device(cut, links, resident)
   !< `cut_links_2d` on the device (OpenACC, issue #100 step 5): the same links bit for bit (np, start, piece, face,
   !< area, unmatched = 0), from a resident cut (`curve%cut_cells_device(..., resident=.true.)`: its arrays are used where
   !< they are) or a host one (copied in for the call and released). The links come from the fluid edge intervals
   !< (`intervals` on, the default; with them off there are none, as on the host); a cut without the interval arrays
   !< (a default cut, a status error) leaves the links reset.
   !<
   !< The same steps as `cut_links_2d`, one thread per edge or per piece (`links_2d_dev`):
   !< 1. the owner piece of every interval (one thread per piece);
   !< 2. the intervals of every edge in interval order: a counting sort by edge without a scatter, each edge gathering
   !<    the intervals of the pieces of its two cells, the lower cell first (an interval lies on an edge of its piece's
   !<    cell, and the pieces are numbered cell by cell: the order of the host's counting sort); count, prefix sum, fill;
   !< 3. per edge `links_edge_sort` (stable by start), the pair count `links_edge_pairs`, a prefix sum, the pair fill
   !<    `links_edge_pairs_fill`;
   !< 4. per piece `links_piece_count_2d`, a prefix sum, `links_piece_fill_2d` (the sums of one key in pair order).
   !< The kernels are the host's, the prefix sums integer: the same bits. With `resident` the links stay mapped (see
   !< `cut_links_t`), else they are copied back. Built without OpenACC the loops run on the host in order (the CI twin).
   type(cut_cells_2d_t), intent(in), target  :: cut      !< Cut cells (resident or not; intervals needed).
   type(cut_links_t),    intent(inout)       :: links    !< Links (released from the device and reset first).
   logical,              intent(in), optional :: resident !< Leave the links on the device (default .false.).
   integer(I4P), allocatable :: st(:), pc(:), fc(:) !< Links: start, piece across, face.
   real(R8P),    allocatable :: ar(:)               !< Links: length.
   integer(I4P), target      :: phi(1)              !< Placeholder of an empty interval array.
   real(R8P),    target      :: phr(1)              !< Placeholder of an empty interval array.
   integer(I4P), pointer, contiguous :: ive(:)      !< Interval edges (or the placeholder).
   real(R8P),    pointer, contiguous :: ivl(:), ivh(:) !< Interval ranges (or the placeholders).
   integer(I4P)              :: n, st0              !< Intervals, status.
   logical                   :: res                 !< Resident.

   call cut_links_reset(links)
   res = .false. ; if (present(resident)) res = resident
   if (.not. (allocated(cut%piece_iv_start) .and. allocated(cut%piece_iv_edge))) return
   links%np = cut%np
   n = size(cut%piece_iv_edge, kind=I4P)
   if (cut%np < 1_I4P) then ! no pieces: start = [1] (as the host), mapped if resident
      allocate(st(1), pc(0), fc(0), ar(0)) ; st(1) = 1_I4P
      if (res) then
         !$acc enter data copyin(st, pc, fc, ar)
      endif
      call links_take(links, st, pc, fc, ar, res)
      return
   endif
   phi = 0_I4P ; phr = 0._R8P
   if (n > 0_I4P) then
      ive => cut%piece_iv_edge ; ivl => cut%piece_iv_lo ; ivh => cut%piece_iv_hi
   else
      ive => phi ; ivl => phr ; ivh => phr
   endif
   call links_2d_dev(cut%nx, cut%ny, cut%np, n, max(1_I4P, n), cut%pieces, cut%first_piece, cut%piece_cell, &
                     cut%piece_iv_start, ive, ivl, ivh, res, st, pc, fc, ar, st0)
   if (st0 /= 0_I4P) then
      write(stderr, '(A)') 'fossil cut_links_2d_device: an interval lies off its piece''s cell edges; links not computed'
      call cut_links_reset(links) ; return
   endif
   call links_take(links, st, pc, fc, ar, res)
   endsubroutine cut_links_2d_device

   subroutine links_2d_dev(nx, ny, np, n, nv, pcs, fpc, pcell, ivs, ive, ivl, ivh, res, lstart, lpiece, lface, larea, &
                           status)
   !< The device stages of `cut_links_2d_device` (see there) on the arrays of a cut of nx x ny cells, np >= 1 pieces and
   !< n fluid edge intervals (nv = max(1, n): the declared size). The inputs are copied in unless present. The outputs
   !< are built in locals mapped by unstructured data (#101 item 8), copied back unless `res`, then moved out (the host
   !< addresses the present table is keyed on are kept). status 1 if the edges did not gather all the intervals.
   integer(I4P),              intent(in)  :: nx, ny, np, n, nv !< Cells, pieces, intervals, their declared size.
   integer(I4P),              intent(in)  :: pcs(nx, ny), fpc(nx, ny) !< Pieces and first piece of every cell.
   integer(I4P),              intent(in)  :: pcell(2, np)  !< Cell of every piece.
   integer(I4P),              intent(in)  :: ivs(np + 1)   !< Intervals of piece p: ivs(p) .. ivs(p + 1) - 1.
   integer(I4P),              intent(in)  :: ive(nv)       !< Interval edges.
   real(R8P),                 intent(in)  :: ivl(nv), ivh(nv) !< Interval ranges.
   logical,                   intent(in)  :: res           !< Leave the outputs mapped.
   integer(I4P), allocatable, intent(out) :: lstart(:), lpiece(:), lface(:) !< Links: start, piece across, face.
   real(R8P),    allocatable, intent(out) :: larea(:)      !< Links: length.
   integer(I4P),              intent(out) :: status        !< 0 ok, 1 intervals off their cell.
   integer(I4P), allocatable :: owner(:), ecnt(:), eoff(:), eord(:), pcnt(:), poff(:), scnt(:), soff(:) !< Work.
   integer(I4P), allocatable :: ppa(:), ppb(:)    !< Pieces of each pair.
   real(R8P),    allocatable :: plen(:)           !< Shared length of each pair.
   integer(I4P), allocatable :: st(:), pc(:), fc(:) !< Links: start, piece across, face.
   real(R8P),    allocatable :: ar(:)             !< Links: length.
   integer(I4P) :: ne, p, t, e, m, tot, ns, nsv, nk !< Work.

   status = 0_I4P
   ne = (nx + 1_I4P) * ny + nx * (ny + 1_I4P)
   allocate(owner(nv), ecnt(ne), eoff(ne + 1), eord(nv), pcnt(ne), poff(ne + 1), scnt(np), soff(np))
   !$acc data copyin(pcs, fpc, pcell, ivs, ive, ivl, ivh) create(owner, ecnt, eoff, eord, pcnt, poff, scnt, soff)
   ! 1. owners
   !$acc parallel loop gang vector default(present) present(ivs, owner)
   do p = 1, np
      do t = ivs(p), ivs(p + 1) - 1
         owner(t) = p
      enddo
   enddo
   ! 2. intervals by edge: counts, prefix sum, gather
   !$acc parallel loop gang vector default(present) present(pcs, fpc, ivs, ive, ecnt)
   do e = 1, ne
      ecnt(e) = links_edge_count_2d(e, nx, ny, np, nv, pcs, fpc, ivs, ive)
   enddo
   call links_scan_dev(ne, ecnt, eoff, tot)
   !$acc serial default(present) present(eoff)
   eoff(ne + 1) = tot
   !$acc end serial
   if (tot /= n) status = 1_I4P
   if (status == 0_I4P) then
   ! 3. per edge: gather, stable sort by start, pair counts; prefix sum; pair fill
   !$acc parallel loop gang vector default(present) present(pcs, fpc, ivs, ive, ivl, ivh, owner, eoff, eord, pcnt)
   do e = 1, ne
      pcnt(e) = 0_I4P
      m = eoff(e + 1) - eoff(e)
      if (m < 1_I4P) cycle
      call links_edge_fill_2d(e, nx, ny, np, nv, pcs, fpc, ivs, ive, m, eord(eoff(e) + 1))
      if (m > 1_I4P) then
         call links_edge_sort(m, eord(eoff(e) + 1), nv, ivl)
         pcnt(e) = links_edge_pairs(m, eord(eoff(e) + 1), nv, owner, ivl, ivh)
      endif
   enddo
   call links_scan_dev(ne, pcnt, poff, ns)
   !$acc serial default(present) present(poff)
   poff(ne + 1) = ns
   !$acc end serial
   nsv = max(1_I4P, ns)
   allocate(ppa(nsv), ppb(nsv), plen(nsv))
   !$acc data create(ppa, ppb, plen)
   !$acc parallel loop gang vector default(present) present(ivl, ivh, owner, eoff, eord, poff, ppa, ppb, plen)
   do e = 1, ne
      if (poff(e + 1) == poff(e)) cycle
      call links_edge_pairs_fill(eoff(e + 1) - eoff(e), eord(eoff(e) + 1), nv, owner, ivl, ivh, poff(e + 1) - poff(e), &
                                 ppa(poff(e) + 1), ppb(poff(e) + 1), plen(poff(e) + 1))
   enddo
   ! 4. CSR per piece: counts, prefix sum, fill
   !$acc parallel loop gang vector default(present) present(pcell, poff, ppa, ppb, scnt)
   do p = 1, np
      scnt(p) = links_piece_count_2d(p, pcell(1, p), pcell(2, p), nx, ny, ne, poff, nsv, ppa, ppb)
   enddo
   call links_scan_dev(np, scnt, soff, nk)
   allocate(st(np + 1), pc(nk), fc(nk), ar(nk))
   !$acc enter data create(st, pc, fc, ar)
   !$acc parallel loop gang vector default(present) present(soff, st)
   do p = 1, np
      st(p) = soff(p) + 1_I4P
   enddo
   !$acc serial default(present) present(st)
   st(np + 1) = nk + 1_I4P
   !$acc end serial
   if (nk > 0_I4P) then
      !$acc parallel loop gang vector default(present) present(pcell, poff, ppa, ppb, plen, st, pc, fc, ar)
      do p = 1, np
         call links_piece_fill_2d(p, pcell(1, p), pcell(2, p), nx, ny, ne, poff, nsv, ppa, ppb, plen, st(p) - 1_I4P, nk, &
                                  pc, fc, ar)
      enddo
   endif
   !$acc end data
   endif ! status == 0
   !$acc end data
   if (status /= 0_I4P) return
   if (.not. res) then
      !$acc exit data copyout(st, pc, fc, ar)
   endif
   call move_alloc(st, lstart) ; call move_alloc(pc, lpiece) ; call move_alloc(fc, lface) ; call move_alloc(ar, larea)
   endsubroutine links_2d_dev

   pure subroutine links_edge_cells_2d(e, nx, ny, ia, ja, ib, jb)
   !< The cells next to edge e (numbering of `cut_cells_2d_t%piece_iv_edge`) in piece order: (ia, ja) the lower one (left
   !< of a vertical edge, below a horizontal one), (ib, jb) the upper one; ia = 0 or ib = 0 if outside the block.
   !$acc routine seq
   integer(I4P), intent(in)  :: e, nx, ny          !< Edge, cells.
   integer(I4P), intent(out) :: ia, ja, ib, jb     !< Cells.
   integer(I4P)              :: l                  !< Work.

   if (e <= (nx + 1_I4P) * ny) then ! vertical edge (u line i, row j)
      ib = mod(e - 1_I4P, nx + 1_I4P) + 1_I4P ; jb = (e - 1_I4P) / (nx + 1_I4P) + 1_I4P ; ia = ib - 1_I4P ; ja = jb
      if (ib > nx) ib = 0_I4P
   else ! horizontal edge (column i, v line j)
      l = e - (nx + 1_I4P) * ny
      ib = mod(l - 1_I4P, nx) + 1_I4P ; jb = (l - 1_I4P) / nx + 1_I4P ; ia = ib ; ja = jb - 1_I4P
      if (ja < 1_I4P) ia = 0_I4P
      if (jb > ny) ib = 0_I4P
   endif
   endsubroutine links_edge_cells_2d

   pure integer(I4P) function links_edge_count_2d(e, nx, ny, np, nv, pcs, fpc, ivs, ive) result(m)
   !< Kernel (issue #100 step 5): the fluid edge intervals on edge e, among those of the pieces of its two cells.
   !$acc routine seq
   integer(I4P), intent(in) :: e, nx, ny, np, nv          !< Edge, cells, pieces, intervals (declared size).
   integer(I4P), intent(in) :: pcs(nx, ny), fpc(nx, ny)   !< Pieces and first piece of every cell.
   integer(I4P), intent(in) :: ivs(np + 1)                !< Intervals by piece.
   integer(I4P), intent(in) :: ive(nv)                    !< Interval edges.
   integer(I4P)             :: ia, ja, ib, jb, sd, ic, jc, p, t !< Work.

   m = 0_I4P
   call links_edge_cells_2d(e, nx, ny, ia, ja, ib, jb)
   do sd = 1, 2
      if (sd == 1) then
         ic = ia ; jc = ja
      else
         ic = ib ; jc = jb
      endif
      if (ic < 1_I4P) cycle
      do p = fpc(ic, jc), fpc(ic, jc) + pcs(ic, jc) - 1_I4P ! none if pieces < 1
         do t = ivs(p), ivs(p + 1) - 1_I4P
            if (ive(t) == e) m = m + 1_I4P
         enddo
      enddo
   enddo
   endfunction links_edge_count_2d

   pure subroutine links_edge_fill_2d(e, nx, ny, np, nv, pcs, fpc, ivs, ive, m, ord)
   !< Kernel (issue #100 step 5): the m intervals on edge e in interval order (the pieces of the lower cell, then of the
   !< upper one, each in piece order, its intervals in order: increasing indexes).
   !$acc routine seq
   integer(I4P), intent(in)  :: e, nx, ny, np, nv         !< Edge, cells, pieces, intervals (declared size).
   integer(I4P), intent(in)  :: pcs(nx, ny), fpc(nx, ny)  !< Pieces and first piece of every cell.
   integer(I4P), intent(in)  :: ivs(np + 1)               !< Intervals by piece.
   integer(I4P), intent(in)  :: ive(nv)                   !< Interval edges.
   integer(I4P), intent(in)  :: m                         !< Intervals of the edge.
   integer(I4P), intent(out) :: ord(m)                    !< Their indexes.
   integer(I4P)              :: ia, ja, ib, jb, sd, ic, jc, p, t, w !< Work.

   w = 0_I4P
   call links_edge_cells_2d(e, nx, ny, ia, ja, ib, jb)
   do sd = 1, 2
      if (sd == 1) then
         ic = ia ; jc = ja
      else
         ic = ib ; jc = jb
      endif
      if (ic < 1_I4P) cycle
      do p = fpc(ic, jc), fpc(ic, jc) + pcs(ic, jc) - 1_I4P
         do t = ivs(p), ivs(p + 1) - 1_I4P
            if (ive(t) == e .and. w < m) then
               w = w + 1_I4P ; ord(w) = t
            endif
         enddo
      enddo
   enddo
   endsubroutine links_edge_fill_2d

   subroutine links_scan_dev(n, cnt, off, total)
   !< Exclusive prefix sums off(l) = cnt(1) + ... + cnt(l - 1) on the device (issue #100 step 5; `scan_dev` of the 3-D
   !< cut cells), total = their sum on the host: chunks of LINKS_SCAN_CHUNK entries (one thread each), the chunk offsets
   !< by `chunk_offsets_dev`, then one thread per chunk writes its entries. Integer sums: the integers of a sequential
   !< scan. cnt and off must be present.
   integer(I4P), intent(in)    :: n        !< Entries.
   integer(I4P), intent(in)    :: cnt(n)   !< Counts.
   integer(I4P), intent(inout) :: off(n)   !< Offsets (written on the device).
   integer(I4P), intent(out)   :: total    !< Sum.
   integer(I4P), allocatable   :: tot(:,:) !< Chunk totals, then offsets [3, nch].
   integer(I4P)                :: nch, ch, l, a1, sums(3) !< Work.

   nch = max(1_I4P, (n + LINKS_SCAN_CHUNK - 1_I4P) / LINKS_SCAN_CHUNK)
   allocate(tot(3, nch))
   !$acc data create(tot)
   !$acc parallel loop gang vector default(present) present(cnt, tot)
   do ch = 1_I4P, nch
      a1 = 0_I4P
      !$acc loop seq
      do l = (ch - 1_I4P) * LINKS_SCAN_CHUNK + 1_I4P, min(ch * LINKS_SCAN_CHUNK, n)
         a1 = a1 + cnt(l)
      enddo
      tot(1, ch) = a1 ; tot(2, ch) = 0_I4P ; tot(3, ch) = 0_I4P
   enddo
   call chunk_offsets_dev(nch, [0_I4P, 0_I4P, 0_I4P], tot, sums)
   total = sums(1)
   !$acc parallel loop gang vector default(present) present(cnt, off, tot)
   do ch = 1_I4P, nch
      a1 = tot(1, ch)
      !$acc loop seq
      do l = (ch - 1_I4P) * LINKS_SCAN_CHUNK + 1_I4P, min(ch * LINKS_SCAN_CHUNK, n)
         off(l) = a1 ; a1 = a1 + cnt(l)
      enddo
   enddo
   !$acc end data
   endsubroutine links_scan_dev

   subroutine links_take(links, st, pc, fc, ar, res)
   !< Move the arrays of a device link stage into `links` (issue #100 step 5; move_alloc keeps the host addresses the
   !< present table is keyed on) and flag them resident if `res`.
   type(cut_links_t),         intent(inout) :: links      !< Links (np and unmatched set).
   integer(I4P), allocatable, intent(inout) :: st(:), pc(:), fc(:) !< Start, piece across, face.
   real(R8P),    allocatable, intent(inout) :: ar(:)      !< Area.
   logical,                   intent(in)    :: res        !< Mapped.

   call move_alloc(st, links%start) ; call move_alloc(pc, links%piece) ; call move_alloc(fc, links%face)
   call move_alloc(ar, links%area)
   if (res) then
      links%on_device = .true. ; links%device_key = links_host_key(links)
   endif
   endsubroutine links_take

   subroutine cut_links_reset(links)
   !< Release resident links from the device (`exit_device`), then reset them to the default (issue #100 step 5): what
   !< `intent(out)` did, without freeing arrays that are still mapped.
   type(cut_links_t), intent(inout) :: links !< Links.
   type(cut_links_t)                :: fresh !< Default.

   if (links%on_device) call links%exit_device
   links = fresh
   endsubroutine cut_links_reset

   function links_host_key(links) result(key)
   !< Host address of start (the key of mapped links; 0 if start is not allocated).
   type(cut_links_t), intent(in), target :: links !< Links.
   integer(c_intptr_t)                   :: key   !< Address.

   key = 0
   if (allocated(links%start)) then
      if (size(links%start) > 0) key = transfer(c_loc(links%start), key)
   endif
   endfunction links_host_key

   logical function links_is_mapped(self)
   !< The flag is set and the arrays are the ones mapped (not a copy of resident links).
   type(cut_links_t), intent(in) :: self !< Links.

   links_is_mapped = .false.
   if (self%on_device) links_is_mapped = self%device_key == links_host_key(self)
   endfunction links_is_mapped

   subroutine links_enter_device(self)
   !< Map start, piece, face and area of host links to the device (copyin; issue #100 step 5). Member arrays only, never
   !< `self` (#20). No-op if already mapped or without links (start unallocated). A copy of resident links is mapped as
   !< new links.
   class(cut_links_t), intent(inout) :: self !< Links.

   if (links_is_mapped(self)) return
   self%on_device = .false. ; self%device_key = 0
   if (.not. (allocated(self%start) .and. allocated(self%piece) .and. allocated(self%face) .and. allocated(self%area))) return
   !$acc enter data copyin(self%start, self%piece, self%face, self%area)
   self%on_device = .true. ; self%device_key = links_host_key(self)
   endsubroutine links_enter_device

   subroutine links_update_host(self)
   !< Copy the mapped arrays back to the host (issue #100 step 5); they stay mapped. No-op if not mapped; on a copy of
   !< resident links a warning on stderr (its arrays are not on the device: update the original before copying it).
   class(cut_links_t), intent(inout) :: self !< Links.

   if (.not. self%on_device) return
   if (.not. links_is_mapped(self)) then
      write(stderr, '(A)') 'fossil cut_links_t%update_host: a copy of resident links is not on the device; '// &
                           'call update_host on the original before copying it'
      self%on_device = .false. ; self%device_key = 0
      return
   endif
   !$acc update host(self%start, self%piece, self%face, self%area)
   endsubroutine links_update_host

   subroutine links_exit_device(self)
   !< Release the device copies (issue #100 step 5); the host arrays are kept as they are (call `update_host` first to
   !< keep the device values). No-op if not mapped; a copy of resident links only drops the flag.
   class(cut_links_t), intent(inout) :: self !< Links.

   if (.not. self%on_device) return
   if (links_is_mapped(self)) then
      !$acc exit data delete(self%start, self%piece, self%face, self%area)
   endif
   self%on_device = .false. ; self%device_key = 0
   endsubroutine links_exit_device

   pure subroutine links_edge_sort(m, ord, n, lo)
   !< Kernel (issue #100 step 5): the intervals ord(1:m) of one edge ordered by start, stable (insertion sort: an
   !< interval moves before an earlier one only if its start is strictly smaller).
   !$acc routine seq
   integer(I4P), intent(in)    :: m      !< Intervals of the edge.
   integer(I4P), intent(inout) :: ord(m) !< Their indexes.
   integer(I4P), intent(in)    :: n      !< Intervals.
   real(R8P),    intent(in)    :: lo(n)  !< Starts.
   integer(I4P)                :: t, u, v !< Work.

   do t = 2, m
      v = ord(t) ; u = t - 1_I4P
      do while (u >= 1_I4P)
         if (.not. lo(v) < lo(ord(u))) exit
         ord(u + 1) = ord(u) ; u = u - 1_I4P
      enddo
      ord(u + 1) = v
   enddo
   endsubroutine links_edge_sort

   pure integer(I4P) function links_edge_pairs(m, ord, n, owner, lo, hi) result(np)
   !< Kernel, count pass (issue #100 step 5): the pairs of overlapping intervals of distinct pieces among the intervals
   !< ord(1:m) of one edge.
   !$acc routine seq
   integer(I4P), intent(in) :: m                !< Intervals of the edge.
   integer(I4P), intent(in) :: ord(m)           !< Their indexes, by start.
   integer(I4P), intent(in) :: n                !< Intervals.
   integer(I4P), intent(in) :: owner(n)         !< Piece of each interval.
   real(R8P),    intent(in) :: lo(n), hi(n)     !< Ranges.
   integer(I4P)             :: kk, q            !< Work.
   real(R8P)                :: ov               !< Overlap.

   np = 0_I4P
   do kk = 1, m
      do q = kk + 1, m
         if (owner(ord(kk)) == owner(ord(q))) cycle
         ov = min(hi(ord(kk)), hi(ord(q))) - max(lo(ord(kk)), lo(ord(q)))
         if (ov <= 0._R8P) cycle
         np = np + 1_I4P
      enddo
   enddo
   endfunction links_edge_pairs

   pure subroutine links_edge_pairs_fill(m, ord, n, owner, lo, hi, np, pa, pb, len)
   !< Kernel, fill pass (issue #100 step 5): the pairs counted by `links_edge_pairs`, in the order (earlier interval,
   !< later interval) of the edge.
   !$acc routine seq
   integer(I4P), intent(in)  :: m                !< Intervals of the edge.
   integer(I4P), intent(in)  :: ord(m)           !< Their indexes, by start.
   integer(I4P), intent(in)  :: n                !< Intervals.
   integer(I4P), intent(in)  :: owner(n)         !< Piece of each interval.
   real(R8P),    intent(in)  :: lo(n), hi(n)     !< Ranges.
   integer(I4P), intent(in)  :: np               !< Pairs.
   integer(I4P), intent(out) :: pa(np), pb(np)   !< Pieces.
   real(R8P),    intent(out) :: len(np)          !< Shared lengths.
   integer(I4P)              :: kk, q, s         !< Work.
   real(R8P)                 :: ov               !< Overlap.

   s = 0_I4P
   do kk = 1, m
      do q = kk + 1, m
         if (owner(ord(kk)) == owner(ord(q))) cycle
         ov = min(hi(ord(kk)), hi(ord(q))) - max(lo(ord(kk)), lo(ord(q)))
         if (ov <= 0._R8P) cycle
         s = s + 1_I4P ; pa(s) = owner(ord(kk)) ; pb(s) = owner(ord(q)) ; len(s) = ov
      enddo
   enddo
   endsubroutine links_edge_pairs_fill

   pure integer(I4P) function links_edge_of_face_2d(i, j, nx, ny, fc) result(e)
   !< Edge id of face fc (1 -u, 2 +u, 3 -v, 4 +v) of cell (i, j) (numbering of `cut_cells_2d_t%piece_iv_edge`).
   !$acc routine seq
   integer(I4P), intent(in) :: i, j, nx, ny, fc !< Cell, cells, face.

   select case (fc)
   case (1) ; e = i + (j - 1_I4P) * (nx + 1_I4P)
   case (2) ; e = i + 1_I4P + (j - 1_I4P) * (nx + 1_I4P)
   case (3) ; e = (nx + 1_I4P) * ny + i + (j - 1_I4P) * nx
   case default ; e = (nx + 1_I4P) * ny + i + j * nx
   endselect
   endfunction links_edge_of_face_2d

   pure integer(I4P) function links_next_q_2d(pw, last, s0, s1, ns, pa, pb) result(q)
   !< The smallest piece above `last` paired with piece pw in slots s0..s1 (0: none).
   !$acc routine seq
   integer(I4P), intent(in) :: pw, last, s0, s1, ns !< Piece, last piece across, slots, all slots.
   integer(I4P), intent(in) :: pa(ns), pb(ns)       !< Pieces of the slots.
   integer(I4P)             :: s, ot               !< Work.

   q = 0_I4P
   do s = s0, s1
      if (pa(s) == pw) then
         ot = pb(s)
      elseif (pb(s) == pw) then
         ot = pa(s)
      else
         cycle
      endif
      if (ot > last .and. (q == 0_I4P .or. ot < q)) q = ot
   enddo
   endfunction links_next_q_2d

   pure integer(I4P) function links_piece_count_2d(pw, i, j, nx, ny, ne, poff, ns, pa, pb) result(nk)
   !< Kernel, count pass (issue #100 step 5): the links of piece pw of cell (i, j): distinct pieces across, per face.
   !$acc routine seq
   integer(I4P), intent(in) :: pw, i, j, nx, ny, ne !< Piece, its cell, cells, edges.
   integer(I4P), intent(in) :: poff(ne + 1)         !< Pairs of edge e: poff(e) + 1 .. poff(e + 1).
   integer(I4P), intent(in) :: ns                   !< Pairs.
   integer(I4P), intent(in) :: pa(ns), pb(ns)       !< Pieces of the pairs.
   integer(I4P)             :: fc, e, q             !< Work.

   nk = 0_I4P
   do fc = 1, 4
      e = links_edge_of_face_2d(i, j, nx, ny, fc)
      q = 0_I4P
      do
         q = links_next_q_2d(pw, q, poff(e) + 1_I4P, poff(e + 1), ns, pa, pb)
         if (q == 0_I4P) exit
         nk = nk + 1_I4P
      enddo
   enddo
   endfunction links_piece_count_2d

   pure subroutine links_piece_fill_2d(pw, i, j, nx, ny, ne, poff, ns, pa, pb, len, k0, nl, piece, face, area)
   !< Kernel, fill pass (issue #100 step 5): the links of piece pw of cell (i, j) written at k0 + 1.., by face, then
   !< piece across; the lengths of one key summed in slot order.
   !$acc routine seq
   integer(I4P), intent(in)    :: pw, i, j, nx, ny, ne !< Piece, its cell, cells, edges.
   integer(I4P), intent(in)    :: poff(ne + 1)         !< Pairs of edge e: poff(e) + 1 .. poff(e + 1).
   integer(I4P), intent(in)    :: ns                   !< Pairs.
   integer(I4P), intent(in)    :: pa(ns), pb(ns)       !< Pieces of the pairs.
   real(R8P),    intent(in)    :: len(ns)              !< Shared lengths.
   integer(I4P), intent(in)    :: k0, nl               !< Links before pw's, all links.
   integer(I4P), intent(inout) :: piece(nl), face(nl)  !< Piece across, face.
   real(R8P),    intent(inout) :: area(nl)             !< Shared length.
   integer(I4P)                :: fc, e, q, s, k       !< Work.
   real(R8P)                   :: sm                   !< Sum.

   k = k0
   do fc = 1, 4
      e = links_edge_of_face_2d(i, j, nx, ny, fc)
      q = 0_I4P
      do
         q = links_next_q_2d(pw, q, poff(e) + 1_I4P, poff(e + 1), ns, pa, pb)
         if (q == 0_I4P) exit
         sm = 0._R8P
         do s = poff(e) + 1, poff(e + 1)
            if ((pa(s) == pw .and. pb(s) == q) .or. (pb(s) == pw .and. pa(s) == q)) sm = sm + len(s)
         enddo
         k = k + 1_I4P ; piece(k) = q ; face(k) = fc ; area(k) = sm
      enddo
   enddo
   endsubroutine links_piece_fill_2d

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

   pure integer(I4P) function links_face_id(a, c, n, base) result(f)
   !< Face (a, c) of a 3-D block of n cells, the face of axis a at the low side of cell c, numbered as the face arrays
   !< of `cut_cells_3d_t` (face_x [n1 + 1, n2, n3], then face_y, then face_z: base(a) faces before axis a).
   !$acc routine seq
   integer(I4P), intent(in) :: a, c(3), n(3), base(4) !< Axis, cell, cells, faces before each axis.
   integer(I4P)             :: m1, m2                  !< Face array extents.

   m1 = n(1) ; m2 = n(2)
   if (a == 1_I4P) m1 = m1 + 1_I4P
   if (a == 2_I4P) m2 = m2 + 1_I4P
   f = base(a) + c(1) + m1 * ((c(2) - 1_I4P) + m2 * (c(3) - 1_I4P))
   endfunction links_face_id

   pure subroutine links_face_bases(n, base)
   !< Faces before each axis in the numbering of `links_face_id`; base(4) = all faces.
   !$acc routine seq
   integer(I4P), intent(in)  :: n(3)    !< Cells.
   integer(I4P), intent(out) :: base(4) !< Faces before axis a.

   base(1) = 0_I4P
   base(2) = (n(1) + 1_I4P) * n(2) * n(3)
   base(3) = base(2) + n(1) * (n(2) + 1_I4P) * n(3)
   base(4) = base(3) + n(1) * n(2) * (n(3) + 1_I4P)
   endsubroutine links_face_bases

   pure integer(I4P) function links_next_q_3d(pw, low, last, s0, s1, ns, pa, pb) result(q)
   !< The smallest piece above `last` linked to piece pw in slots s0..s1 of a face, pw on the low side (pa) if low,
   !< else on the high side (pb); 0: none. A slot with pa = 0 holds no link.
   !$acc routine seq
   integer(I4P), intent(in) :: pw, last, s0, s1, ns !< Piece, last piece across, slots, all slots.
   logical,      intent(in) :: low                  !< pw in the low cell of the face.
   integer(I4P), intent(in) :: pa(ns), pb(ns)       !< Pieces of the slots (low, high cell).
   integer(I4P)             :: s, me, ot            !< Work.

   q = 0_I4P
   do s = s0, s1
      if (pa(s) == 0_I4P) cycle
      if (low) then
         me = pa(s) ; ot = pb(s)
      else
         me = pb(s) ; ot = pa(s)
      endif
      if (me /= pw) cycle
      if (ot > last .and. (q == 0_I4P .or. ot < q)) q = ot
   enddo
   endfunction links_next_q_3d

   pure integer(I4P) function links_piece_count_3d(pw, c, n, base, nfc, off, ns, pa, pb) result(nk)
   !< Kernel, count pass (issue #100 step 5): the links of piece pw of cell c: distinct pieces across, per face.
   !$acc routine seq
   integer(I4P), intent(in) :: pw, c(3), n(3), base(4) !< Piece, its cell, cells, faces before each axis.
   integer(I4P), intent(in) :: nfc                     !< Faces.
   integer(I4P), intent(in) :: off(nfc + 1)            !< Slots of face f: off(f) + 1 .. off(f + 1).
   integer(I4P), intent(in) :: ns                      !< Slots.
   integer(I4P), intent(in) :: pa(ns), pb(ns)          !< Pieces of the slots (low, high cell; pa = 0: no link).
   integer(I4P)             :: fc, a, f, q, cc(3)      !< Work.
   logical                  :: low                     !< pw in the low cell of the face.

   nk = 0_I4P
   do fc = 1, 6
      a = (fc + 1_I4P) / 2_I4P ; low = mod(fc, 2_I4P) == 0_I4P ; cc = c
      if (low) then
         if (c(a) >= n(a)) cycle
         cc(a) = c(a) + 1_I4P
      else
         if (c(a) <= 1_I4P) cycle
      endif
      f = links_face_id(a, cc, n, base)
      q = 0_I4P
      do
         q = links_next_q_3d(pw, low, q, off(f) + 1_I4P, off(f + 1), ns, pa, pb)
         if (q == 0_I4P) exit
         nk = nk + 1_I4P
      enddo
   enddo
   endfunction links_piece_count_3d

   pure subroutine links_piece_fill_3d(pw, c, n, base, nfc, off, ns, pa, pb, ar, k0, nl, piece, face, area)
   !< Kernel, fill pass (issue #100 step 5): the links of piece pw of cell c written at k0 + 1.., by face, then piece
   !< across; the areas of one key summed in slot order.
   !$acc routine seq
   integer(I4P), intent(in)    :: pw, c(3), n(3), base(4) !< Piece, its cell, cells, faces before each axis.
   integer(I4P), intent(in)    :: nfc                     !< Faces.
   integer(I4P), intent(in)    :: off(nfc + 1)            !< Slots of face f: off(f) + 1 .. off(f + 1).
   integer(I4P), intent(in)    :: ns                      !< Slots.
   integer(I4P), intent(in)    :: pa(ns), pb(ns)          !< Pieces of the slots (low, high cell; pa = 0: no link).
   real(R8P),    intent(in)    :: ar(ns)                  !< Areas of the slots.
   integer(I4P), intent(in)    :: k0, nl                  !< Links before pw's, all links.
   integer(I4P), intent(inout) :: piece(nl), face(nl)     !< Piece across, face.
   real(R8P),    intent(inout) :: area(nl)                !< Shared area.
   integer(I4P)                :: fc, a, f, q, s, k, cc(3), me, ot !< Work.
   logical                     :: low                     !< pw in the low cell of the face.
   real(R8P)                   :: sm                      !< Sum.

   k = k0
   do fc = 1, 6
      a = (fc + 1_I4P) / 2_I4P ; low = mod(fc, 2_I4P) == 0_I4P ; cc = c
      if (low) then
         if (c(a) >= n(a)) cycle
         cc(a) = c(a) + 1_I4P
      else
         if (c(a) <= 1_I4P) cycle
      endif
      f = links_face_id(a, cc, n, base)
      q = 0_I4P
      do
         q = links_next_q_3d(pw, low, q, off(f) + 1_I4P, off(f + 1), ns, pa, pb)
         if (q == 0_I4P) exit
         sm = 0._R8P
         do s = off(f) + 1, off(f + 1)
            if (pa(s) == 0_I4P) cycle
            if (low) then
               me = pa(s) ; ot = pb(s)
            else
               me = pb(s) ; ot = pa(s)
            endif
            if (me == pw .and. ot == q) sm = sm + ar(s)
         enddo
         k = k + 1_I4P ; piece(k) = q ; face(k) = fc ; area(k) = sm
      enddo
   enddo
   endsubroutine links_piece_fill_3d
endmodule fossil_cut_graph
