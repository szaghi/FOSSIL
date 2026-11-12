!< FOSSIL, the piece map of a time step of a moving body (issue #102, WP3 of #92).
module fossil_piece_map
!< FOSSIL, the piece map of a time step of a moving body (issue #102, WP3 of #92).
!<
!< In a cell c, the space-time fluid region of a step is {(p, s) : p in c, p fluid at pose s, s in [0, 1]}; its
!< connected components are the **tracks** of c. Every piece of the cut cells at pose0 lies in exactly one track, and so
!< does every piece at pose1. A track with no old piece is fresh, one with no new piece is dead; several pieces on one
!< side make a split or a merge; a track with no piece at either end (fluid that crosses the cell within the step) is
!< transient. A body passing through a cell within the step separates its tracks, which an overlap of the two end poses
!< cannot see; a gap moving by more than its width keeps one track, which the overlap would call dead + fresh.
!<
!< The map is built by `swept_cells_2d(..., map=)` (`fossil_swept_cells`): the cell's topology changes only at the
!< instants the swept faces already isolate (a wall vertex crossing a face line, a corner crossing a wall segment); the
!< pieces are evaluated once per interval between them and linked across each instant by the boundary features they
!< share (wall segments, fluid corners: a feature present on both sides of one event persists through it, so the fluid
!< next to it is one track). The kernels here are the device-ready parts (plain arrays, no allocation): the union of
!< the evaluations' pieces into tracks and the status of a track.
!<
!< In 3-D (`swept_cells_3d(..., map=)`, issue #102 step 2) the instants are a mesh vertex crossing a face plane of the
!< cell inside the face, a mesh edge crossing one of its 12 edges and one of its 8 corners crossing a facet; the
!< features are the facets (a planar facet clipped to a convex cell bounds exactly one piece) and the fluid corners. The
!< construction is the 2-D one.
use penf, only : I4P, R8P

implicit none
private
public :: piece_map_2d_t, piece_map_3d_t
public :: TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, TRACK_SPLIT, TRACK_MERGE, TRACK_SPLIT_MERGE, TRACK_TRANSIENT, &
          TRACK_UNTRUSTED
public :: track_components, track_status_of

integer(I4P), parameter :: TRACK_CONTINUING  = 1_I4P  !< One old piece, one new piece.
integer(I4P), parameter :: TRACK_FRESH       = 2_I4P  !< No old piece (fluid that appears in the cell).
integer(I4P), parameter :: TRACK_DEAD        = 3_I4P  !< No new piece (fluid that leaves the cell or is covered).
integer(I4P), parameter :: TRACK_SPLIT       = 4_I4P  !< One old piece, several new pieces.
integer(I4P), parameter :: TRACK_MERGE       = 5_I4P  !< Several old pieces, one new piece.
integer(I4P), parameter :: TRACK_SPLIT_MERGE = 6_I4P  !< Several old pieces and several new pieces.
integer(I4P), parameter :: TRACK_TRANSIENT   = 7_I4P  !< No piece at either end (fluid crossing the cell within the step).
integer(I4P), parameter :: TRACK_UNTRUSTED   = -1_I4P !< The cell's tracks could not be built (one track holds all its
                                                      !< pieces and the cell's measures; see `piece_map_2d_t`).

type :: piece_map_2d_t
   !< The piece map of one time step on a 2-D block (issue #102): the tracks of every cell, CSR like `cut_links_t`.
   !<
   !< Pieces are those of `curve%cut_cells` at pose0 (old) and at pose1 (new), numbered as the cut cells number them
   !< (`swept_cells_2d`'s end-pose cuts: a copy moved to the pose, then `cut_cells`). Tracks of cell (i, j):
   !< first_track(i, j) .. first_track(i, j) + tracks(i, j) - 1, ordered by their first piece (old pieces first, in
   !< piece order, then the pieces born within the step in time order). Old pieces of track t:
   !< old_piece(old_start(t) : old_start(t + 1) - 1), new pieces likewise; old_track / new_track invert them.
   !<
   !< Measures per track: volume0 / volume1 (sum of its pieces' areas), the time-averaged fluid length of the four faces
   !< of its cell taken by the track (faces(-u, +u, -v, +v)), the wall -sum_f A_f n_f, and the residual: the volume
   !< balance of the rigid field, volume1 - volume0 + int_0^1 (outflow through the track's face parts) ds, as
   !< `swept_cells_2d_t%residual` per cell (round-off on the exact paths, the quadrature error on the Gauss path). A cell
   !< with one track takes the cell's own values of `swept_cells_2d_t` (faces, wall, residual: bit for bit); a cell with
   !< several tracks integrates them per track on the cell's event intervals, on every path: two-point Gauss per interval
   !< without rotation, Gauss doubled to round-off with one (a track lives between events, so the step's nq nodes of the
   !< Gauss path would miss or overweight it). Its tracks' residuals are then round-off, and their sums over the cell are
   !< the cell's values to round-off on the exact paths; on the Gauss path they differ from the cell's Gauss values by
   !< that quadrature error (about the nq-vs-2nq estimate). An untrusted cell (status TRACK_UNTRUSTED: an evaluation did
   !< not close, its pieces did not match the end-pose cuts, a breakpoint was not certified, an event was a contact not
   !< certified by a sign change (`contact_cells`), or a face part could not be assigned) has one track with all its
   !< pieces and the cell's values.
   !<
   !< `outflow(f, t)`: the time integral over the step of the outward flux of the **rigid** field of the screw motion
   !< through the track's part of face f, so that residual = volume1 - volume0 + sum_f outflow(f, t). It is what makes a
   !< uniform flow moving with the body exact (the geometric conservation law of a moving cut cell): split each face's
   !< flux among the tracks in proportion to it. A solver with its own velocity field uses the apertures, and the rigid
   !< outflow is its consistency term. One track: the cell's face outflows; several: integrated with the apertures.
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P !< Cells.
   integer(I4P)              :: nt = 0_I4P            !< Tracks.
   integer(I4P)              :: np0 = 0_I4P           !< Pieces at pose0.
   integer(I4P)              :: np1 = 0_I4P           !< Pieces at pose1.
   integer(I4P), allocatable :: tracks(:,:)           !< Tracks per cell [nx, ny].
   integer(I4P), allocatable :: first_track(:,:)      !< First track of each cell [nx, ny].
   integer(I4P), allocatable :: track_cell(:,:)       !< Cell (i, j) of each track [2, nt].
   integer(I4P), allocatable :: status(:)             !< TRACK_* of each track [nt].
   integer(I4P), allocatable :: old_start(:)          !< Old pieces of each track [nt + 1].
   integer(I4P), allocatable :: old_piece(:)          !< Old pieces (pose0 piece indexes) [np0].
   integer(I4P), allocatable :: new_start(:)          !< New pieces of each track [nt + 1].
   integer(I4P), allocatable :: new_piece(:)          !< New pieces (pose1 piece indexes) [np1].
   integer(I4P), allocatable :: old_track(:)          !< Track of each pose0 piece [np0].
   integer(I4P), allocatable :: new_track(:)          !< Track of each pose1 piece [np1].
   real(R8P),    allocatable :: volume0(:)            !< Fluid area of the old pieces [nt].
   real(R8P),    allocatable :: volume1(:)            !< Fluid area of the new pieces [nt].
   real(R8P),    allocatable :: faces(:,:)            !< Time-averaged fluid length of the faces (-u, +u, -v, +v) [4, nt].
   real(R8P),    allocatable :: wall(:,:)             !< Time-averaged projected wall, -sum_f A_f n_f [2, nt].
   real(R8P),    allocatable :: residual(:)           !< Volume balance of the rigid field [nt].
   real(R8P),    allocatable :: outflow(:,:)          !< Rigid-field outflow through each face, outward [4, nt].
   logical                   :: exact = .true.        !< Every multi-track cell's integration converged (rotation).
   integer(I4P)              :: n_status(-1:7) = 0_I4P !< Tracks per status (index TRACK_*; 0 unused).
   integer(I4P)              :: untrusted_cells = 0_I4P !< Cells marked untrusted.
   integer(I4P)              :: multi_track_cells = 0_I4P !< Cells with more than one track.
   integer(I4P)              :: evaluated_cells = 0_I4P !< Cells whose pieces were evaluated within the step.
   integer(I4P)              :: evaluations = 0_I4P   !< One-cell piece evaluations (topology pass).
   integer(I4P)              :: skipped_nodes = 0_I4P !< Quadrature nodes or pieces skipped inside intervals between two
                                                      !< events closer than 1e-10 (their weight is below that).
   integer(I4P)              :: contact_cells = 0_I4P !< Untrusted cells with an event not certified by a sign change: a
                                                      !< contact within the round-off of its function under a rotation
                                                      !< (a corner or vertex grazing, where a piece can shrink to a point
                                                      !< and come back as another track; issue #102).
endtype piece_map_2d_t

type :: piece_map_3d_t
   !< The piece map of one time step on a 3-D block (issue #102 step 2): the tracks of every cell, CSR like
   !< `cut_links_t`; the 3-D twin of `piece_map_2d_t`, built by `swept_cells_3d(..., map=)`.
   !<
   !< Pieces are those of `cut_cells_3d` at pose0 (old) and at pose1 (new), numbered as the cut cells number them
   !< (`swept_cells_3d`'s end-pose cuts: a copy moved to the pose, then `cut_cells_3d`). Tracks of cell (i, j, k):
   !< first_track(i, j, k) .. first_track(i, j, k) + tracks(i, j, k) - 1, ordered by their first piece (old pieces first,
   !< in piece order, then the pieces born within the step in time order). Old pieces of track t:
   !< old_piece(old_start(t) : old_start(t + 1) - 1), new pieces likewise; old_track / new_track invert them.
   !<
   !< Measures per track: volume0 / volume1 (sum of its pieces' volumes), the time-averaged fluid area of the six faces
   !< of its cell taken by the track (faces(-x, +x, -y, +y, -z, +z)), the time integral over the step of the rigid
   !< field's outflow through each of them (outflow, outward from the cell: a low face counts the flux in -axis), the
   !< wall -sum_f A_f n_f, and the residual: the volume balance of the rigid field, volume1 - volume0 +
   !< sum_f outflow(f), as `swept_cells_3d_t%residual` per cell (round-off on the exact paths, the quadrature error on
   !< the Gauss path). A cell with one track takes the cell's own values of `swept_cells_3d_t` (faces, wall, residual:
   !< bit for bit; outflow: the face outflows of the call); a cell with several tracks integrates them per track on the
   !< cell's event intervals with the exact paths' rules on every path (two-point Gauss without rotation, Gauss doubled
   !< to round-off with one): their sum over the cell's tracks is the cell's value to round-off on the exact paths, and
   !< differs from it by the cell's quadrature error on the Gauss path (`exact=.false.`; of the order of
   !< `swept_cells_3d_t%error_estimate`), where the per-track values are the accurate ones. An untrusted cell
   !< (status TRACK_UNTRUSTED: an end cut or an evaluation not trusted (pieces -1), its end evaluations not equal to the
   !< end-pose cuts, a breakpoint not certified, a contact not certified by a sign change (`contact_cells`), or a face
   !< part that could not be assigned) has one track with all its pieces and the cell's values.
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P, nz = 0_I4P !< Cells.
   integer(I4P)              :: nt = 0_I4P            !< Tracks.
   integer(I4P)              :: np0 = 0_I4P           !< Pieces at pose0.
   integer(I4P)              :: np1 = 0_I4P           !< Pieces at pose1.
   integer(I4P), allocatable :: tracks(:,:,:)         !< Tracks per cell [nx, ny, nz].
   integer(I4P), allocatable :: first_track(:,:,:)    !< First track of each cell [nx, ny, nz].
   integer(I4P), allocatable :: track_cell(:,:)       !< Cell (i, j, k) of each track [3, nt].
   integer(I4P), allocatable :: status(:)             !< TRACK_* of each track [nt].
   integer(I4P), allocatable :: old_start(:)          !< Old pieces of each track [nt + 1].
   integer(I4P), allocatable :: old_piece(:)          !< Old pieces (pose0 piece indexes) [np0].
   integer(I4P), allocatable :: new_start(:)          !< New pieces of each track [nt + 1].
   integer(I4P), allocatable :: new_piece(:)          !< New pieces (pose1 piece indexes) [np1].
   integer(I4P), allocatable :: old_track(:)          !< Track of each pose0 piece [np0].
   integer(I4P), allocatable :: new_track(:)          !< Track of each pose1 piece [np1].
   real(R8P),    allocatable :: volume0(:)            !< Fluid volume of the old pieces [nt].
   real(R8P),    allocatable :: volume1(:)            !< Fluid volume of the new pieces [nt].
   real(R8P),    allocatable :: faces(:,:)            !< Time-averaged fluid area of the faces (-x, +x, -y, +y, -z, +z)
                                                      !< [6, nt].
   real(R8P),    allocatable :: outflow(:,:)          !< Time-integrated outflow of the rigid field per face [6, nt].
   real(R8P),    allocatable :: wall(:,:)             !< Time-averaged projected wall, -sum_f A_f n_f [3, nt].
   real(R8P),    allocatable :: residual(:)           !< Volume balance of the rigid field [nt].
   logical                   :: exact = .true.        !< Every multi-track cell's integration converged (rotation).
   integer(I4P)              :: n_status(-1:7) = 0_I4P !< Tracks per status (index TRACK_*; 0 unused).
   integer(I4P)              :: untrusted_cells = 0_I4P !< Cells marked untrusted.
   integer(I4P)              :: multi_track_cells = 0_I4P !< Cells with more than one track.
   integer(I4P)              :: evaluated_cells = 0_I4P !< Cells whose pieces were evaluated within the step.
   integer(I4P)              :: evaluations = 0_I4P   !< One-cell piece evaluations (topology pass).
   integer(I4P)              :: skipped_nodes = 0_I4P !< Quadrature nodes or pieces skipped inside intervals between two
                                                      !< events closer than 1e-10 (their weight is below that).
   integer(I4P)              :: contact_cells = 0_I4P !< Untrusted cells with an event not certified by a sign change: a
                                                      !< contact within the round-off of its function under a rotation
                                                      !< (a corner or vertex grazing, where a piece can shrink to a point
                                                      !< and come back as another track; issue #102).
endtype piece_map_3d_t

contains
   pure subroutine track_components(ne, nstart, lstart, nl, lkey, lpiece, par, comp, nt)
   !< Kernel (issue #102): the tracks of one cell from the pieces of its ne successive evaluations. Node (e, p), piece p
   !< of evaluation e, is nstart(e) + p - 1; evaluation e lists its boundary features lkey(lstart(e) : lstart(e+1) - 1)
   !< (increasing keys, each once) with the piece lpiece next to each. Pieces of consecutive evaluations that share a
   !< key are united (union-find, the smaller root wins); comp(node) is the track, numbered 1..nt in the order of each
   !< track's first node.
   !$acc routine seq
   integer(I4P), intent(in)    :: ne                 !< Evaluations.
   integer(I4P), intent(in)    :: nstart(ne + 1)     !< First node of each evaluation (nodes 1 .. nstart(ne+1) - 1).
   integer(I4P), intent(in)    :: lstart(ne + 1)     !< First key of each evaluation.
   integer(I4P), intent(in)    :: nl                 !< Keys.
   integer(I4P), intent(in)    :: lkey(nl)           !< Keys (sorted per evaluation).
   integer(I4P), intent(in)    :: lpiece(nl)         !< Piece next to each key.
   integer(I4P), intent(inout) :: par(nstart(ne + 1) - 1)  !< Work: parents [nodes].
   integer(I4P), intent(inout) :: comp(nstart(ne + 1) - 1) !< Track of each node [nodes].
   integer(I4P), intent(out)   :: nt                 !< Tracks.
   integer(I4P) :: nn, e, a, b, ae, be, ra, rb, k    !< Work.

   nn = nstart(ne + 1) - 1_I4P
   do k = 1_I4P, nn
      par(k) = k
   enddo
   do e = 1_I4P, ne - 1_I4P
      a = lstart(e) ; ae = lstart(e + 1) - 1_I4P ; b = lstart(e + 1) ; be = lstart(e + 2) - 1_I4P
      do while (a <= ae .and. b <= be) ! merge-join of the two sorted key lists
         if (lkey(a) < lkey(b)) then
            a = a + 1_I4P
         elseif (lkey(a) > lkey(b)) then
            b = b + 1_I4P
         else
            ra = uf_root(nstart(e) + lpiece(a) - 1_I4P, nn, par) ; rb = uf_root(nstart(e + 1) + lpiece(b) - 1_I4P, nn, par)
            if (ra < rb) then
               par(rb) = ra
            elseif (rb < ra) then
               par(ra) = rb
            endif
            a = a + 1_I4P ; b = b + 1_I4P
         endif
      enddo
   enddo
   nt = 0_I4P
   do k = 1_I4P, nn ! a root is the smallest node of its track: numbered in node order
      if (par(k) == k) then
         nt = nt + 1_I4P ; comp(k) = nt
      endif
   enddo
   do k = 1_I4P, nn
      comp(k) = comp(uf_root(k, nn, par))
   enddo
   endsubroutine track_components

   pure integer(I4P) function uf_root(kk, nn, par) result(r)
   !< Union-find root of node kk.
   !$acc routine seq
   integer(I4P), intent(in) :: kk, nn   !< Node, nodes.
   integer(I4P), intent(in) :: par(nn)  !< Parents.

   r = kk
   do while (par(r) /= r)
      r = par(r)
   enddo
   endfunction uf_root

   pure integer(I4P) function track_status_of(nold, nnew) result(st)
   !< Status of a track with nold old pieces and nnew new pieces (issue #102).
   !$acc routine seq
   integer(I4P), intent(in) :: nold, nnew !< Old and new pieces.

   if (nold == 0_I4P .and. nnew == 0_I4P) then
      st = TRACK_TRANSIENT
   elseif (nold == 0_I4P) then
      st = TRACK_FRESH
   elseif (nnew == 0_I4P) then
      st = TRACK_DEAD
   elseif (nold == 1_I4P .and. nnew == 1_I4P) then
      st = TRACK_CONTINUING
   elseif (nold == 1_I4P) then
      st = TRACK_SPLIT
   elseif (nnew == 1_I4P) then
      st = TRACK_MERGE
   else
      st = TRACK_SPLIT_MERGE
   endif
   endfunction track_status_of
endmodule fossil_piece_map
