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
use fossil_cut_moments, only : mono_count, mono_exponents, segment_green_moments, cell_moments_into, triangle_rule, tensor_rule
use fossil_predicates, only : orient2d_xy
use iso_c_binding,     only : c_intptr_t, c_loc
use iso_fortran_env,   only : stderr => error_unit
use penf,              only : I4P, R8P

implicit none
private
public :: cut_cells_2d_t, compute_cut_cells_2d, compute_cut_cells_2d_device, cut_cells_2d_reset
public :: cut_cells_2d_blocks_t, compute_cut_cells_2d_blocks_device, compute_cut_cells_2d_blocks_grids_device, &
          cut_cells_2d_blocks_reset
! the stages of the two paths, public for the twin test (issue #100 step 6a; not re-exported by `fossil`)
public :: node_flags, cut_cells_2d_prepare, cut_cells_2d_classify, cut_cells_2d_scan, cut_cells_2d_classify_dev, &
          cut_cells_2d_scan_dev
! the device chunk scan, shared with the 3-D device path (issue #100 step 3c)
public :: chunk_offsets_dev
! one cell with its pieces' boundary features, for the piece map of a time step (issue #102)
public :: cut_cell_2d_pieces, cut_cell_2d_labels
! one cell's fluid fraction and centroid, for the faces of the 3-D piece map (issue #102 step 2)
public :: cut_cell_2d_measures
! the one-cell kernels of those, for the piece map on the device (issue #105)
public :: cut_cell_2d_walk, cut_cell_2d_index, cut_cell_2d_one_sizes, cut_cell_2d_one_count, cut_cell_2d_one_pieces

type :: cut_cells_2d_t
   !< Cut-cell geometry of a block (fluid = outside the solid).
   !<
   !< Resident on the device (issue #100 step 6a). `curve%cut_cells_device(..., resident=.true.)` leaves the outputs
   !< mapped (OpenACC unstructured data) and never writes their host copies: until `update_host` the host arrays hold
   !< garbage (only nx, ny, np, moments_order are set). `on_device` says so. The mapped outputs are closed, alpha,
   !< aperture, centroid, wet_length, wet_centroid, wet_normal, pieces, first_piece, node_in, piece_cell, piece_area,
   !< piece_centroid, piece_wet_length, piece_wet_normal, piece_aperture and piece_iv_* (moments, boundaries and
   !< quadrature are host-only). A device consumer reads them through pointer aliases with `present` (the #20 rule:
   !< member arrays, never the object). `enter_device` maps a host result, `update_host` brings the device values back,
   !< `exit_device` releases the device copies (call it before a resident cut goes out of scope or is deallocated: the
   !< runtime would keep the stale mapping). Every routine that writes a cut (`curve%cut_cells`, `cut_cells_device`,
   !< `cut_coarsen_2d`) releases it first. A routine that reads one on the host (`cut_links_2d`, `cut_coarsen_2d`, ...)
   !< needs `update_host` first. A copy `b = a` of a resident cut is not resident (its arrays are new host arrays, not
   !< mapped, holding a's host values): the type-bound procedures recognise it (`device_key`, the host address of alpha
   !< when mapped) and treat it as not mapped; `update_host` on it warns and does nothing. Built without OpenACC the
   !< procedures only set the flag.
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
   logical,      allocatable :: node_in(:,:)        !< Grid node inside the solid [nx+1, ny+1] (perturbed, exact; #70).
   integer(I4P)              :: moments_order = 0_I4P !< Order of the moments computed (issue #70; 0: none).
   real(R8P),    allocatable :: piece_moments(:,:) !< int X^a Y^b dA of each piece [mono_count(2, p), np] (`fossil_cut_moments`).
   real(R8P),    allocatable :: piece_wall_moments(:,:) !< int X^a Y^b ds over the wall of each piece [mono_count(2, p), np].
   integer(I4P), allocatable :: piece_seg_start(:)  !< Oriented fluid boundary of piece p (issue #70, `boundaries`):
   real(R8P),    allocatable :: piece_seg(:,:)      !< segments piece_seg(1:4, start(p) .. start(p+1)-1) = (a, b), the fluid
   logical,      allocatable :: piece_seg_wall(:)   !< on the left (counter-clockwise); wall segments flagged.
   integer(I4P), allocatable :: quad_start(:)       !< Area rule of piece p (issue #70, `quadrature`): signed triangles
   real(R8P),    allocatable :: quad_point(:,:)     !< from the cell centre over its boundary, 7 points each (degree 5),
   real(R8P),    allocatable :: quad_weight(:)      !< a fluid cell the 3 x 3 Gauss rule.
   integer(I4P), allocatable :: wall_quad_start(:)  !< Wall rule of piece p: 3-point Gauss on each wall segment.
   real(R8P),    allocatable :: wall_quad_point(:,:) !< Its points [2, nw].
   real(R8P),    allocatable :: wall_quad_weight(:) !< Its weights (lengths).
   logical                   :: on_device = .false. !< The outputs are mapped on the device (issue #100 step 6a).
   integer(c_intptr_t)       :: device_key = 0      !< Host address of alpha when mapped (a copy has another one).
   contains
      procedure, pass(self) :: enter_device !< Map the outputs to the device (copyin).
      procedure, pass(self) :: update_host  !< Copy the mapped outputs back to the host.
      procedure, pass(self) :: exit_device  !< Release the device copies.
endtype cut_cells_2d_t

type :: cut_cells_2d_blocks_t
   !< Cut cells of nb blocks of nx x ny cells, each with its own curve, on one grid x(nx+1), y(ny+1) (issue #100 step 3b:
   !< the plane sections of one axis of `cut_cells_3d`; `compute_cut_cells_2d_blocks_device`) or on a grid per block
   !< x(nx+1, nb), y(ny+1, nb) (step 4.1: the translated planes of the swept geometry, the 2-D AMR batching of equal-size
   !< blocks; `compute_cut_cells_2d_blocks_grids_device`). Block b of every cell array
   !< is the cut of its curve, bit for bit (moments order 0, no boundaries, no quadrature): `cell(..., b)` is the
   !< `cut_cells_2d_t` array of the same name. With pieces (`with_pieces`), the pieces of all blocks are concatenated
   !< in block order: block b owns pieces block_piece(b) + 1 .. block_piece(b + 1) and fluid edge intervals
   !< block_iv(b) + 1 .. block_iv(b + 1); first_piece and piece_iv_start index the concatenated arrays (minus
   !< block_piece(b), block_iv(b): the block's own numbering); piece_cell holds the cell (i, j) in its block.
   !< Resident (`on_device`): every allocated array is mapped (OpenACC unstructured data) and its host copy is garbage
   !< until `update_host`; `exit_device` releases them (before the result goes out of scope). Do not copy a resident
   !< result (the copy's arrays are not mapped).
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P, nb = 0_I4P !< Cells per block, blocks.
   logical,      allocatable :: closed(:,:,:)        !< The cell boundary closed [nx, ny, nb].
   real(R8P),    allocatable :: alpha(:,:,:)         !< Fluid area fraction [nx, ny, nb].
   real(R8P),    allocatable :: aperture(:,:,:,:)    !< Edge fluid fractions [4, nx, ny, nb].
   real(R8P),    allocatable :: centroid(:,:,:,:)    !< Fluid centroid [2, nx, ny, nb].
   real(R8P),    allocatable :: wet_length(:,:,:)    !< Wetted length [nx, ny, nb].
   real(R8P),    allocatable :: wet_centroid(:,:,:,:) !< Wetted centroid [2, nx, ny, nb].
   real(R8P),    allocatable :: wet_normal(:,:,:,:)  !< Outward-normal integral over the wall [2, nx, ny, nb].
   integer(I4P), allocatable :: pieces(:,:,:)        !< Fluid pieces per cell [nx, ny, nb].
   logical,      allocatable :: node_in(:,:,:)       !< Grid node inside the solid [nx+1, ny+1, nb].
   logical                   :: with_pieces = .false. !< The piece arrays below were computed.
   integer(I4P), allocatable :: first_piece(:,:,:)   !< First piece of the cell in the concatenated arrays [nx, ny, nb].
   integer(I4P), allocatable :: block_piece(:)       !< Pieces before block b [nb + 1].
   integer(I4P), allocatable :: block_iv(:)          !< Fluid edge intervals before block b [nb + 1].
   integer(I4P)              :: np = 0_I4P           !< Pieces of all blocks.
   integer(I4P), allocatable :: piece_cell(:,:)      !< Cell (i, j) of each piece in its block [2, np].
   real(R8P),    allocatable :: piece_area(:)        !< Fluid area of each piece [np].
   real(R8P),    allocatable :: piece_centroid(:,:)  !< Fluid centroid of each piece [2, np].
   real(R8P),    allocatable :: piece_wet_length(:)  !< Wetted length of each piece [np].
   real(R8P),    allocatable :: piece_wet_normal(:,:) !< Outward-normal integral of each piece [2, np].
   real(R8P),    allocatable :: piece_aperture(:,:)  !< Edge fluid fractions per piece [4, np].
   integer(I4P), allocatable :: piece_iv_start(:)    !< First interval of each piece [np + 1] (concatenated).
   integer(I4P), allocatable :: piece_iv_edge(:)     !< Interval edges (numbering of `cut_cells_2d_t`).
   real(R8P),    allocatable :: piece_iv_lo(:)       !< Interval starts.
   real(R8P),    allocatable :: piece_iv_hi(:)       !< Interval ends.
   logical                   :: on_device = .false.  !< The arrays are mapped on the device.
   contains
      procedure, pass(self) :: update_host => blocks_update_host !< Copy the mapped arrays back to the host.
      procedure, pass(self) :: exit_device => blocks_exit_device !< Release the device copies.
endtype cut_cells_2d_blocks_t

type :: block_prep_t
   !< Host preparation of one block of `compute_cut_cells_2d_blocks_device` (`node_rows`, `cut_cells_2d_prepare`).
   integer(I4P)              :: nxs, nps                    !< Crossings, wall pieces.
   integer(I4P), allocatable :: rstart(:), rseg(:)          !< Segments crossing each node row.
   integer(I4P), allocatable :: xstart(:), xord(:)          !< Crossings by edge.
   integer(I4P), allocatable :: pstart(:), pord(:)          !< Wall pieces by cell.
   real(R8P),    allocatable :: xp(:,:), pab(:,:)           !< Crossing points, wall piece end points.
   integer(I4P), allocatable :: ptg(:,:)                    !< Wall piece end tags.
endtype block_prep_t

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

interface take
   !< Reuse an output array of the previous result when its shape is the new one (issue #100 step 6a).
   module procedure take_l2, take_i1, take_i2, take_r1, take_r2, take_r3
endinterface take

integer(I4P), parameter :: SCAN_SEG = 64_I4P !< Cells per chunk of the device classification and prefix sums (#100).
integer(I4P), parameter :: SCAN_SUP = 64_I4P !< Chunks per super-chunk of their offsets (`chunk_offsets_dev`).
integer(I4P), parameter :: NHDR = 7_I4P !< Header of a cut cell's integer workspace (issue #100): status (0 nothing to
                                        !< write, 1 pieces, 2 the whole cell fluid), ns, nh, nsb, chains and loops, crossings,
                                        !< wall pieces.

contains
   subroutine compute_cut_cells_2d(u, v, next, x, y, cc, status, with_iv, order, with_bnd, with_quad)
   !< Cut cells of the curve (u, v, next) on the block of edges x(nx+1), y(ny+1).
   !<
   !< Segment i runs from vertex i to vertex next(i) (next(i) = 0: no segment). `status` = 0, or 1 on
   !< bad edges.
   !<
   !< Per-cell kernels, count then fill (issue #100): the segments are walked through the grid on the host (`walk_all`,
   !< `index_crossings`, `index_pieces`, as flat arrays); then
   !< 1. every cell is classified quiet (no wall piece, no crossing, equal corner flags) or cut, and each cut cell gets
   !<    a workspace slice sized by `cut_cell_2d_sizes` from its crossing and wall-piece counts;
   !< 2. `cut_cell_2d_count` gives each cut cell its cell measures and the exact numbers of its pieces, fluid edge
   !<    intervals and boundary segments, keeping its pieces in its workspace;
   !< 3. a prefix sum over the cells in order gives every cell its first piece, interval and boundary segment;
   !< 4. `quiet_cell_2d_kernel` and `cut_cell_2d_fill` write each cell's pieces into its slices of the outputs.
   !< The kernels take plain arrays and scalars only (the `!$acc routine seq` form of #20), so the same code can run per
   !< cell on the device; the outputs are those of the cell-by-cell loop with growing buffers that they replace, bit for
   !< bit and in the same order. The preparation (`cut_cells_2d_grid`, `node_flags`, `cut_cells_2d_prepare`,
   !< `cut_cells_2d_classify`) and the prefix sums (`cut_cells_2d_scan`) are shared with the device twin
   !< `compute_cut_cells_2d_device`, which runs the same kernels in OpenACC loops.
   real(R8P),            intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)  :: next(:)    !< Next vertex (0: none).
   real(R8P),            intent(in)    :: x(:), y(:) !< Cell edges (strictly increasing).
   type(cut_cells_2d_t), intent(inout) :: cc         !< Cut cells (released from the device and reset first).
   integer(I4P),         intent(out)   :: status     !< 0 ok, 1 bad edges.
   logical,              intent(in)  :: with_iv    !< Record the fluid edge intervals of the pieces (issue #70).
   integer(I4P),         intent(in)  :: order      !< Order of the moments (0: none; issue #70).
   logical,              intent(in)  :: with_bnd   !< Record the oriented boundaries of the pieces (issue #70).
   logical,              intent(in)  :: with_quad  !< Quadrature rules of the pieces (issue #70; implies with_bnd).
   logical,       allocatable :: node_in(:,:)            !< Node inside flags.
   integer(I4P),  allocatable :: xstart(:), xord(:)      !< Crossings by edge.
   integer(I4P),  allocatable :: pstart(:), pord(:)      !< Wall pieces by cell.
   real(R8P),     allocatable :: xp(:,:)                 !< Crossing points [2, nxs].
   real(R8P),     allocatable :: pab(:,:)                !< Wall piece end points (a, b) [4, nps].
   integer(I4P),  allocatable :: ptg(:,:)                !< Wall piece end tags (ta, tb) [2, nps].
   integer(I4P),  allocatable :: me(:,:)                 !< Moment exponents [2, nm].
   integer(I4P),  allocatable :: ckind(:)                !< Per cell: -1 quiet solid, 0 quiet fluid, k > 0 cut cell k.
   integer(I4P),  allocatable :: cutc(:)                 !< Cell of cut cell k.
   integer(I4P),  allocatable :: rwo(:), iwo(:)          !< Workspace slices of the cut cells [ncut + 1].
   integer(I4P),  allocatable :: cnt(:,:)                !< Intervals, boundary segments, pieces of each cut cell [3, ncut].
   integer(I4P),  allocatable :: ivo(:), bdo(:)          !< Intervals, boundary segments before each cell [ncell].
   real(R8P),     allocatable :: rw(:)                   !< Real workspace.
   integer(I4P),  allocatable :: iw(:)                   !< Integer workspace.
   real(R8P),     allocatable :: pmom(:,:), pwmom(:,:)   !< Piece moments, wall moments [nm, np].
   integer(I4P),  allocatable :: bstart(:)               !< Boundary segments of the pieces [np + 1].
   real(R8P),     allocatable :: bseg(:,:)               !< Boundary segments [4, nbd].
   logical,       allocatable :: bwall(:)                !< Wall flags [nbd].
   integer(I4P) :: nx, ny, nxs, nps, ne, ncell, c, i, j, k, ncut, nm, np, niv, nbd, eid(4) !< Work.
   logical      :: wbnd                                  !< Record the boundaries.

   call cut_cells_2d_reset(cc)
   call cut_cells_2d_grid(x, y, nx, ny, status)
   if (status /= 0_I4P) return
   cc%nx = nx ; cc%ny = ny
   allocate(cc%alpha(nx, ny), cc%aperture(4, nx, ny), cc%centroid(2, nx, ny), cc%wet_length(nx, ny), &
            cc%wet_centroid(2, nx, ny), cc%wet_normal(2, nx, ny), cc%pieces(nx, ny), cc%first_piece(nx, ny), &
            cc%closed(nx, ny))
   call node_flags(u, v, next, x, y, node_in)
   cc%node_in = node_in
   call cut_cells_2d_prepare(u, v, next, x, y, nx, ny, ne, ncell, nxs, nps, xstart, xord, pstart, pord, xp, pab, ptg)
   nm = 0_I4P ; if (order > 0_I4P) nm = mono_count(2_I4P, order)
   cc%moments_order = max(0_I4P, order)
   if (order > 0_I4P) then
      me = mono_exponents(2_I4P, order)
   else
      allocate(me(2, 0))
   endif
   wbnd = with_bnd .or. with_quad
   ! 1. classify the cells; workspace slices of the cut cells
   call cut_cells_2d_classify(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ckind, ncut, cutc, rwo, iwo)
   allocate(cnt(3, ncut), rw(rwo(ncut + 1) - 1), iw(iwo(ncut + 1) - 1))
   ! 2. count: the cut cells' measures and output counts
   do k = 1_I4P, ncut
      c = cutc(k) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      call cut_cell_2d_count(x(i), x(i + 1), y(j), y(j + 1), node_in(i, j), node_in(i + 1, j), node_in(i + 1, j + 1), &
                             node_in(i, j + 1), eid, ne, nxs, xstart, xord, xp, pstart(c + 1) - pstart(c), pstart(c) - 1_I4P, &
                             nps, pord, pab, ptg, nm, me, with_iv, wbnd, rwo(k + 1) - rwo(k), rw(rwo(k)), &
                             iwo(k + 1) - iwo(k), iw(iwo(k)), cc%closed(i, j), cc%alpha(i, j), cc%aperture(1, i, j), &
                             cc%centroid(1, i, j), cc%wet_length(i, j), cc%wet_centroid(1, i, j), cc%wet_normal(1, i, j), &
                             cnt(3, k), cnt(1, k), cnt(2, k))
      cc%pieces(i, j) = cnt(3, k)
   enddo
   ! 3. prefix sums over the cells in order
   allocate(ivo(ncell), bdo(ncell))
   call cut_cells_2d_scan(nx, ny, ncell, ckind, ncut, cnt, with_iv, wbnd, cc%first_piece, ivo, bdo, np, niv, nbd)
   ! 4. fill
   cc%np = np
   allocate(cc%piece_cell(2, np), cc%piece_area(np), cc%piece_centroid(2, np), cc%piece_wet_length(np), &
            cc%piece_wet_normal(2, np), cc%piece_aperture(4, np), pmom(nm, np), pwmom(nm, np), cc%piece_iv_start(np + 1), &
            cc%piece_iv_edge(niv), cc%piece_iv_lo(niv), cc%piece_iv_hi(niv), bstart(np + 1), bseg(4, nbd), bwall(nbd))
   do c = 1_I4P, ncell
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      k = ckind(c)
      if (k <= 0_I4P) then
         call quiet_cell_2d_kernel(x(i), x(i + 1), y(j), y(j + 1), i, j, k < 0_I4P, eid, nm, me, with_iv, wbnd, &
                                   cc%closed(i, j), cc%alpha(i, j), cc%aperture(1, i, j), cc%centroid(1, i, j), &
                                   cc%wet_length(i, j), cc%wet_centroid(1, i, j), cc%wet_normal(1, i, j), cc%pieces(i, j), &
                                   cc%first_piece(i, j) - 1_I4P, ivo(c), bdo(c), np, cc%piece_cell, cc%piece_area, &
                                   cc%piece_centroid, cc%piece_wet_length, cc%piece_wet_normal, cc%piece_aperture, pmom, &
                                   pwmom, cc%piece_iv_start, niv, cc%piece_iv_edge, cc%piece_iv_lo, cc%piece_iv_hi, bstart, &
                                   nbd, bseg, bwall)
      else
         call cut_cell_2d_fill(x(i), x(i + 1), y(j), y(j + 1), i, j, eid, nm, me, with_iv, wbnd, rwo(k + 1) - rwo(k), &
                               rw(rwo(k)), iwo(k + 1) - iwo(k), iw(iwo(k)), cc%first_piece(i, j) - 1_I4P, ivo(c), bdo(c), &
                               np, cc%piece_cell, cc%piece_area, cc%piece_centroid, cc%piece_wet_length, &
                               cc%piece_wet_normal, cc%piece_aperture, pmom, pwmom, cc%piece_iv_start, niv, &
                               cc%piece_iv_edge, cc%piece_iv_lo, cc%piece_iv_hi, bstart, nbd, bseg, bwall)
      endif
   enddo
   cc%piece_iv_start(np + 1) = niv + 1_I4P ; bstart(np + 1) = nbd + 1_I4P
   if (nm > 0_I4P) then
      call move_alloc(pmom, cc%piece_moments) ; call move_alloc(pwmom, cc%piece_wall_moments)
   endif
   if (wbnd) then
      call move_alloc(bstart, cc%piece_seg_start) ; call move_alloc(bseg, cc%piece_seg) ; call move_alloc(bwall, cc%piece_seg_wall)
   endif
   if (with_quad) call quadrature_out(cc, x, y)
   endsubroutine compute_cut_cells_2d

   subroutine compute_cut_cells_2d_device(u, v, next, x, y, cc, status, with_iv, resident)
   !< `compute_cut_cells_2d` on the device (OpenACC, issue #100 step 2), moments order 0, no boundaries, no quadrature:
   !< the same outputs bit for bit.
   !<
   !< The host does what `compute_cut_cells_2d` does before the kernels, through the same routines: the grid check, the
   !< walk and the indexes (`cut_cells_2d_prepare`) and the per-row segment lists of the node flags (`node_rows`). The
   !< device runs the kernels of the host path, one thread per node or per cell, and the integer stages (step 6a):
   !< 1. `node_parity` per node (today's exact per-pair tests over the row's segments);
   !< 2. the cell classification and the workspace sizing (`cut_cells_2d_classify_dev`, the integers of
   !<    `cut_cells_2d_classify`; ncut and the pool sizes come back);
   !< 3. `cut_cell_2d_count` over the cut cells, into a workspace pool created on the device;
   !< 4. the prefix sums (`cut_cells_2d_scan_dev`, the integers of `cut_cells_2d_scan`; np, niv, nbd come back);
   !< 5. `quiet_cell_2d_kernel` over all cells (quiet ones act) and `cut_cell_2d_fill` over the cut cells;
   !< 6. the outputs are copied back, or with `resident` left on the device (step 6a; see `cut_cells_2d_t`).
   !< Only scalars come back to the host before the outputs.
   !< The outputs are built in local arrays mapped by unstructured data (`enter data`), then moved into `cc` by
   !< `move_alloc`, which keeps their host addresses: the runtime's present table, keyed on them, stays valid.
   !< The host arrays of `cc` are reused when their shape is the new one (`take`; a resident `cc` is released first): a
   !< copy back into freshly allocated pageable memory pays the first-touch page faults (about 3.7 GB/s against 13 GB/s
   !< into touched pages at 1024^2, 97 ms against 25 ms per call in a loop on the same `cc`).
   !< The same kernels, the same operations in the same order, `-gpu=nofma` (#20): the same bits. Built without OpenACC
   !< the directives are comments and the loops run on the host in order.
   real(R8P),            intent(in)           :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)           :: next(:)    !< Next vertex (0: none).
   real(R8P),            intent(in)           :: x(:), y(:) !< Cell edges (strictly increasing).
   type(cut_cells_2d_t), intent(inout)        :: cc         !< Cut cells (released from the device and reset first).
   integer(I4P),         intent(out)          :: status     !< 0 ok, 1 bad edges.
   logical,              intent(in)           :: with_iv    !< Record the fluid edge intervals of the pieces (issue #70).
   logical,              intent(in), optional :: resident   !< Leave the outputs on the device (default .false.).
   logical,       allocatable :: node_in(:,:)            !< Node inside flags [nx+1, ny+1].
   integer(I4P),  allocatable :: rstart(:), rseg(:)      !< Segments crossing each node row.
   integer(I4P),  allocatable :: xstart(:), xord(:)      !< Crossings by edge.
   integer(I4P),  allocatable :: pstart(:), pord(:)      !< Wall pieces by cell.
   real(R8P),     allocatable :: xp(:,:)                 !< Crossing points [2, nxs].
   real(R8P),     allocatable :: pab(:,:)                !< Wall piece end points (a, b) [4, nps].
   integer(I4P),  allocatable :: ptg(:,:)                !< Wall piece end tags (ta, tb) [2, nps].
   integer(I4P),  allocatable :: ckind(:)                !< Per cell: -1 quiet solid, 0 quiet fluid, k > 0 cut cell k.
   integer(I4P),  allocatable :: cutc(:)                 !< Cell of cut cell k.
   integer(I4P),  allocatable :: rwo(:), iwo(:)          !< Workspace slices of the cut cells [ncut + 1].
   integer(I4P),  allocatable :: cnt(:,:)                !< Intervals, boundary segments, pieces of each cut cell [3, ncut].
   integer(I4P),  allocatable :: ivo(:), bdo(:)          !< Intervals, boundary segments before each cell [ncell].
   real(R8P),     allocatable :: rw(:)                   !< Real workspace pool.
   integer(I4P),  allocatable :: iw(:)                   !< Integer workspace pool.
   logical,       allocatable :: closed(:,:)             !< Cell outputs.
   real(R8P),     allocatable :: alpha(:,:), aperture(:,:,:), centroid(:,:,:), wet_length(:,:) !< Cell outputs.
   real(R8P),     allocatable :: wet_centroid(:,:,:), wet_normal(:,:,:) !< Cell outputs.
   integer(I4P),  allocatable :: pieces(:,:), first_piece(:,:) !< Cell outputs.
   integer(I4P),  allocatable :: pcell(:,:), ivs(:), ive(:), bds(:) !< Piece outputs (bds: discarded).
   real(R8P),     allocatable :: parea(:), pcen(:,:), pwl(:), pwn(:,:), pap(:,:), ivl(:), ivh(:) !< Piece outputs.
   integer(I4P) :: zme(2, 1)                             !< No moment exponents (nm = 0).
   real(R8P)    :: zmom(1, 1), zwmom(1, 1), zbseg(4, 1)  !< No moments, no wall moments, no boundary segments.
   logical      :: zbwall(1)                             !< No wall flags.
   logical      :: res                                   !< Leave the outputs on the device.
   integer(I4P) :: nx, ny, nxs, nps, ne, ncell, nu, nr, c, i, j, k, ncut, np, niv, nbd, nrw, niw, eid(4) !< Work.

   if (cc%on_device) call cc%exit_device ! the host arrays of cc are kept for reuse (`take`), then cc is reset
   res = .false. ; if (present(resident)) res = resident
   call cut_cells_2d_grid(x, y, nx, ny, status)
   if (status /= 0_I4P) then
      call cut_cells_2d_reset(cc) ; return
   endif
   nu = size(u, kind=I4P)
   call node_rows(u, v, next, y, rstart, rseg)
   nr = size(rseg, kind=I4P)
   call cut_cells_2d_prepare(u, v, next, x, y, nx, ny, ne, ncell, nxs, nps, xstart, xord, pstart, pord, xp, pab, ptg)
   zme = 0_I4P ; zmom = 0._R8P ; zwmom = 0._R8P ; zbseg = 0._R8P ; zbwall = .false.
   call take(cc%node_in, node_in, nx + 1_I4P, ny + 1_I4P)
   !$acc data copyin(x, y, zme, zmom, zwmom, zbseg, zbwall)
   ! 1. node flags, one thread per node (an output: it stays on the device)
   !$acc enter data create(node_in)
   !$acc data copyin(u, v, next, rstart, rseg)
   !$acc parallel loop gang vector collapse(2) default(present) present(node_in)
   do j = 1_I4P, ny + 1_I4P
      do i = 1_I4P, nx + 1_I4P
         node_in(i, j) = node_parity(x(i), y(j), rstart(j), rstart(j + 1) - 1_I4P, nr, rseg, nu, u, v, next)
      enddo
   enddo
   !$acc end data
   !$acc data copyin(xstart, xord, xp, pstart, pord, pab, ptg)
   ! 2. classification and workspace sizing (ckind, cutc, rwo, iwo stay on the device; ncut and the totals come back)
   call cut_cells_2d_classify_dev(nx, ny, ne, ncell, xstart, pstart, node_in, 0_I4P, ckind, ncut, cutc, rwo, iwo, nrw, niw)
   allocate(cnt(3, max(1_I4P, ncut)), rw(max(1_I4P, nrw)), iw(max(1_I4P, niw)))
   call take(cc%closed, closed, nx, ny) ; call take(cc%alpha, alpha, nx, ny) ; call take(cc%aperture, aperture, 4_I4P, nx, ny)
   call take(cc%centroid, centroid, 2_I4P, nx, ny) ; call take(cc%wet_length, wet_length, nx, ny)
   call take(cc%wet_centroid, wet_centroid, 2_I4P, nx, ny) ; call take(cc%wet_normal, wet_normal, 2_I4P, nx, ny)
   call take(cc%pieces, pieces, nx, ny) ; call take(cc%first_piece, first_piece, nx, ny)
   allocate(ivo(ncell), bdo(ncell))
   !$acc enter data create(closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces, first_piece)
   !$acc data create(rw, iw, cnt, ivo, bdo)
   ! 3. count over the cut cells
   !$acc parallel loop gang vector default(present) private(eid) &
   !$acc& present(node_in, cutc, rwo, iwo, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal)
   do k = 1_I4P, ncut
      c = cutc(k) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      call cut_cell_2d_count(x(i), x(i + 1), y(j), y(j + 1), node_in(i, j), node_in(i + 1, j), node_in(i + 1, j + 1), &
                             node_in(i, j + 1), eid, ne, nxs, xstart, xord, xp, pstart(c + 1) - pstart(c), pstart(c) - 1_I4P, &
                             nps, pord, pab, ptg, 0_I4P, zme, with_iv, .false., rwo(k + 1) - rwo(k), rw(rwo(k)), &
                             iwo(k + 1) - iwo(k), iw(iwo(k)), closed(i, j), alpha(i, j), aperture(1, i, j), &
                             centroid(1, i, j), wet_length(i, j), wet_centroid(1, i, j), wet_normal(1, i, j), &
                             cnt(3, k), cnt(1, k), cnt(2, k))
   enddo
   ! 4. prefix sums (first_piece, ivo, bdo stay on the device; np, niv, nbd come back)
   call cut_cells_2d_scan_dev(nx, ny, ncell, ckind, ncut, cnt, with_iv, .false., first_piece, ivo, bdo, np, niv, nbd)
   call take(cc%piece_cell, pcell, 2_I4P, np) ; call take(cc%piece_area, parea, np)
   call take(cc%piece_centroid, pcen, 2_I4P, np) ; call take(cc%piece_wet_length, pwl, np)
   call take(cc%piece_wet_normal, pwn, 2_I4P, np) ; call take(cc%piece_aperture, pap, 4_I4P, np)
   call take(cc%piece_iv_start, ivs, np + 1_I4P)
   if (niv > 0_I4P) then
      call take(cc%piece_iv_edge, ive, niv) ; call take(cc%piece_iv_lo, ivl, niv) ; call take(cc%piece_iv_hi, ivh, niv)
   else
      allocate(ive(1), ivl(1), ivh(1))
   endif
   allocate(bds(np + 1))
   !$acc enter data create(pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
   !$acc data create(bds)
   ! 5. fill: the quiet cells, then the cut cells
   !$acc parallel loop gang vector default(present) private(eid) &
   !$acc& present(ckind, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces, first_piece, &
   !$acc&         pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
   do c = 1_I4P, ncell
      if (ckind(c) > 0_I4P) cycle
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      call quiet_cell_2d_kernel(x(i), x(i + 1), y(j), y(j + 1), i, j, ckind(c) < 0_I4P, eid, 0_I4P, zme, with_iv, .false., &
                                closed(i, j), alpha(i, j), aperture(1, i, j), centroid(1, i, j), wet_length(i, j), &
                                wet_centroid(1, i, j), wet_normal(1, i, j), pieces(i, j), first_piece(i, j) - 1_I4P, ivo(c), &
                                bdo(c), np, pcell, parea, pcen, pwl, pwn, pap, zmom, zwmom, ivs, niv, ive, ivl, ivh, bds, &
                                0_I4P, zbseg, zbwall)
   enddo
   !$acc parallel loop gang vector default(present) private(eid) &
   !$acc& present(cutc, rwo, iwo, pieces, first_piece, pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
   do k = 1_I4P, ncut
      c = cutc(k) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      pieces(i, j) = cnt(3, k)
      call cut_cell_2d_fill(x(i), x(i + 1), y(j), y(j + 1), i, j, eid, 0_I4P, zme, with_iv, .false., rwo(k + 1) - rwo(k), &
                            rw(rwo(k)), iwo(k + 1) - iwo(k), iw(iwo(k)), first_piece(i, j) - 1_I4P, ivo(c), bdo(c), np, &
                            pcell, parea, pcen, pwl, pwn, pap, zmom, zwmom, ivs, niv, ive, ivl, ivh, bds, 0_I4P, zbseg, zbwall)
   enddo
   !$acc serial default(present) present(ivs)
   ivs(np + 1) = niv + 1_I4P
   !$acc end serial
   !$acc end data
   !$acc end data
   !$acc exit data delete(ckind, cutc, rwo, iwo)
   !$acc end data
   !$acc end data
   if (.not. res) then ! copy back the outputs
      !$acc exit data copyout(node_in, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces, &
      !$acc&                  first_piece, pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
   endif
   call cut_cells_2d_reset(cc) ! frees what was not reused
   cc%nx = nx ; cc%ny = ny ; cc%moments_order = 0_I4P
   call move_alloc(node_in, cc%node_in)
   call move_alloc(closed, cc%closed) ; call move_alloc(alpha, cc%alpha) ; call move_alloc(aperture, cc%aperture)
   call move_alloc(centroid, cc%centroid) ; call move_alloc(wet_length, cc%wet_length)
   call move_alloc(wet_centroid, cc%wet_centroid) ; call move_alloc(wet_normal, cc%wet_normal)
   call move_alloc(pieces, cc%pieces) ; call move_alloc(first_piece, cc%first_piece)
   cc%np = np
   call move_alloc(pcell, cc%piece_cell) ; call move_alloc(parea, cc%piece_area) ; call move_alloc(pcen, cc%piece_centroid)
   call move_alloc(pwl, cc%piece_wet_length) ; call move_alloc(pwn, cc%piece_wet_normal)
   call move_alloc(pap, cc%piece_aperture) ; call move_alloc(ivs, cc%piece_iv_start)
   if (niv > 0_I4P) then
      call move_alloc(ive, cc%piece_iv_edge) ; call move_alloc(ivl, cc%piece_iv_lo) ; call move_alloc(ivh, cc%piece_iv_hi)
   else ! the one-element placeholders go; the empty outputs are mapped too, so every output is present when resident
      if (res) then
         !$acc exit data delete(ive, ivl, ivh)
      endif
      allocate(cc%piece_iv_edge(0), cc%piece_iv_lo(0), cc%piece_iv_hi(0))
      if (res) then
         !$acc enter data create(cc%piece_iv_edge, cc%piece_iv_lo, cc%piece_iv_hi)
      endif
   endif
   if (res) call mark_resident(cc)
   endsubroutine compute_cut_cells_2d_device

   subroutine compute_cut_cells_2d_blocks_device(nbk, vstart, u, v, next, x, y, cb, status, with_iv, with_pieces, resident)
   !< The cut cells of nbk blocks that share the grid x(nx+1), y(ny+1), each with its own curve, in one set of launches
   !< (issue #100 step 3b; the plane sections of one axis of `cut_cells_3d`): block b is `compute_cut_cells_2d` of its
   !< curve bit for bit, for the outputs of `compute_cut_cells_2d_device` (moments order 0, no boundaries). The work is
   !< `cut_cells_2d_blocks_core` on the one grid; see there. status: 0, or 1 on bad edges or nbk < 1.
   integer(I4P),                intent(in)           :: nbk         !< Blocks.
   integer(I4P),                intent(in)           :: vstart(nbk + 1) !< First vertex of each curve (vstart(1) = 1).
   real(R8P),                   intent(in)           :: u(:), v(:)  !< Vertices of all curves.
   integer(I4P),                intent(in)           :: next(:)     !< Next vertex within the block (0: none).
   real(R8P),                   intent(in)           :: x(:), y(:)  !< Cell edges (strictly increasing).
   type(cut_cells_2d_blocks_t), intent(inout)        :: cb          !< Cut cells (released from the device, reset first).
   integer(I4P),                intent(out)          :: status      !< 0 ok, 1 bad edges or no block.
   logical,                     intent(in)           :: with_iv     !< Record the fluid edge intervals (with_pieces).
   logical,                     intent(in)           :: with_pieces !< Compute the pieces.
   logical,                     intent(in), optional :: resident    !< Leave the outputs on the device (default .false.).
   integer(I4P)                                      :: nx, ny      !< Cells per block.

   call cut_cells_2d_blocks_reset(cb)
   status = 1_I4P
   if (nbk < 1_I4P) return
   call cut_cells_2d_grid(x, y, nx, ny, status)
   if (status /= 0_I4P) return
   call cut_cells_2d_blocks_core(nbk, 1_I4P, vstart, u, v, next, nx, ny, reshape(x, [nx + 1_I4P, 1_I4P]), &
                                 reshape(y, [ny + 1_I4P, 1_I4P]), cb, with_iv, with_pieces, resident)
   endsubroutine compute_cut_cells_2d_blocks_device

   subroutine compute_cut_cells_2d_blocks_grids_device(nbk, vstart, u, v, next, x, y, cb, status, with_iv, with_pieces, &
                                                       resident)
   !< The cut cells of nbk blocks, each with its own curve and its own grid x(:, b), y(:, b) of the same nx x ny cells
   !< (issue #100 step 4.1: the translated planes of the swept geometry; the 2-D AMR batching of equal-size blocks), in
   !< one set of launches: block b is `compute_cut_cells_2d` of its curve on x(:, b), y(:, b) bit for bit, for the
   !< outputs of `compute_cut_cells_2d_device`. The work is `cut_cells_2d_blocks_core`; see there. status: 0, or 1 on
   !< bad edges of any block, a grid count other than nbk, or nbk < 1.
   integer(I4P),                intent(in)           :: nbk         !< Blocks.
   integer(I4P),                intent(in)           :: vstart(nbk + 1) !< First vertex of each curve (vstart(1) = 1).
   real(R8P),                   intent(in)           :: u(:), v(:)  !< Vertices of all curves.
   integer(I4P),                intent(in)           :: next(:)     !< Next vertex within the block (0: none).
   real(R8P),                   intent(in)           :: x(:,:), y(:,:) !< Cell edges of each block [nx+1, nbk], [ny+1, nbk].
   type(cut_cells_2d_blocks_t), intent(inout)        :: cb          !< Cut cells (released from the device, reset first).
   integer(I4P),                intent(out)          :: status      !< 0 ok, 1 bad edges or no block.
   logical,                     intent(in)           :: with_iv     !< Record the fluid edge intervals (with_pieces).
   logical,                     intent(in)           :: with_pieces !< Compute the pieces.
   logical,                     intent(in), optional :: resident    !< Leave the outputs on the device (default .false.).
   integer(I4P)                                      :: nx, ny, b   !< Cells per block, block.

   call cut_cells_2d_blocks_reset(cb)
   status = 1_I4P
   if (nbk < 1_I4P .or. size(x, 2) /= nbk .or. size(y, 2) /= nbk) return
   do b = 1_I4P, nbk
      call cut_cells_2d_grid(x(:, b), y(:, b), nx, ny, status)
      if (status /= 0_I4P) return
   enddo
   call cut_cells_2d_blocks_core(nbk, nbk, vstart, u, v, next, nx, ny, x, y, cb, with_iv, with_pieces, resident)
   endsubroutine compute_cut_cells_2d_blocks_grids_device

   subroutine cut_cells_2d_blocks_core(nbk, ngr, vstart, u, v, next, nx, ny, xg, yg, cb, with_iv, with_pieces, resident)
   !< The cut cells of nbk blocks, each with its own curve, on ngr grids of nx x ny cells (ngr = 1: one grid for all
   !< blocks; ngr = nbk: block b on grid b), in one set of launches (issue #100 steps 3b, 4.1): block b is
   !< `compute_cut_cells_2d` of its curve on its grid bit for bit, for the outputs of `compute_cut_cells_2d_device`
   !< (moments order 0, no boundaries). The grids are checked by the callers.
   !<
   !< Curve b is u, v, next(vstart(b) : vstart(b + 1) - 1), next numbered within the block (1.., 0: no segment).
   !< Per block the host runs what the one-block twin runs (`node_rows`, `cut_cells_2d_prepare`: the walk and the
   !< indexes); their arrays are concatenated (offsets per block), and each kernel reads its block's part through the
   !< start element of the part (explicit-shape dummies), in block-local numbering, so the kernels see exactly the
   !< arrays of the one-block path. The cells of all blocks are numbered g = (b - 1) nx ny + c; then, as in
   !< `compute_cut_cells_2d_device`:
   !< 1. `node_parity` per node of every block;
   !< 2. the classification and the workspace sizing over all cells (`cut_cells_2d_classify_blocks_dev`, the chunked
   !<    integer scan of `cut_cells_2d_classify_dev`): every cut cell of every block gets its workspace slice;
   !< 3. `cut_cell_2d_count` over all cut cells;
   !< 4. with_pieces: the prefix sums over all cells (`cut_cells_2d_scan_dev` on the stacked grid nx x (ny nbk)), the
   !<    block offsets, then `quiet_cell_2d_kernel` and `cut_cell_2d_fill`; without, only the quiet cells' measures
   !<    (`quiet_cell_2d_measures`): no piece arrays.
   !< Integer offsets and counts only between the stages (exact in any order); the same kernels in the same order per
   !< cell, `-gpu=nofma`: the same bits. Built without OpenACC the loops run on the host in order. `resident` leaves
   !< the outputs mapped (see `cut_cells_2d_blocks_t`). Block b reads grid gb = min(b, ngr): xg(:, gb), yg(:, gb), in
   !< the host preparation and in every kernel (the per-block offsets come from the same count-then-scan either way).
   integer(I4P),                intent(in)           :: nbk         !< Blocks.
   integer(I4P),                intent(in)           :: ngr         !< Grids: 1 (shared) or nbk (one per block).
   integer(I4P),                intent(in)           :: vstart(nbk + 1) !< First vertex of each curve (vstart(1) = 1).
   real(R8P),                   intent(in)           :: u(:), v(:)  !< Vertices of all curves.
   integer(I4P),                intent(in)           :: next(:)     !< Next vertex within the block (0: none).
   integer(I4P),                intent(in)           :: nx, ny      !< Cells per block.
   real(R8P),                   intent(in)           :: xg(nx + 1, ngr), yg(ny + 1, ngr) !< Cell edges of the grids.
   type(cut_cells_2d_blocks_t), intent(inout)        :: cb          !< Cut cells (reset by the callers).
   logical,                     intent(in)           :: with_iv     !< Record the fluid edge intervals (with_pieces).
   logical,                     intent(in)           :: with_pieces !< Compute the pieces.
   logical,                     intent(in), optional :: resident    !< Leave the outputs on the device (default .false.).
   type(block_prep_t), allocatable :: pp(:)              !< Host preparation per block.
   real(R8P),     allocatable :: ua(:), va(:)            !< Vertices (one padding element).
   integer(I4P),  allocatable :: na(:)                   !< Next vertex (one padding element).
   integer(I4P),  allocatable :: rsa(:,:), rso(:), rsg(:) !< Node rows: starts [ny + 2, nbk], offsets, segments.
   integer(I4P),  allocatable :: xsa(:,:), xo(:), xoa(:) !< Crossings by edge [ne + 1, nbk], offsets, order.
   real(R8P),     allocatable :: xpa(:,:)                !< Crossing points.
   integer(I4P),  allocatable :: psa(:,:), po(:), poa(:) !< Wall pieces by cell [ncell + 1, nbk], offsets, order.
   real(R8P),     allocatable :: paba(:,:)               !< Wall piece end points.
   integer(I4P),  allocatable :: ptga(:,:)               !< Wall piece end tags.
   integer(I4P),  allocatable :: ckind(:), cutc(:), rwo(:), iwo(:) !< Cell kinds, cut cells, workspace slices.
   integer(I4P),  allocatable :: cnt(:,:), ivo(:), bdo(:) !< Counts of the cut cells, offsets per cell.
   real(R8P),     allocatable :: rw(:)                   !< Real workspace pool.
   integer(I4P),  allocatable :: iw(:)                   !< Integer workspace pool.
   logical,       allocatable :: node_in(:,:,:), closed(:,:,:) !< Cell outputs.
   real(R8P),     allocatable :: alpha(:,:,:), aperture(:,:,:,:), centroid(:,:,:,:), wet_length(:,:,:) !< Cell outputs.
   real(R8P),     allocatable :: wet_centroid(:,:,:,:), wet_normal(:,:,:,:) !< Cell outputs.
   integer(I4P),  allocatable :: pieces(:,:,:), first_piece(:,:,:), bpc(:), biv(:) !< Cell outputs, block offsets.
   integer(I4P),  allocatable :: pcell(:,:), ivs(:), ive(:), bds(:) !< Piece outputs (bds: discarded).
   real(R8P),     allocatable :: parea(:), pcen(:,:), pwl(:), pwn(:,:), pap(:,:), ivl(:), ivh(:) !< Piece outputs.
   integer(I4P) :: zme(2, 1)                             !< No moment exponents (nm = 0).
   real(R8P)    :: zmom(1, 1), zwmom(1, 1), zbseg(4, 1)  !< No moments, no wall moments, no boundary segments.
   logical      :: zbwall(1)                             !< No wall flags.
   logical      :: res                                   !< Leave the outputs on the device.
   integer(I4P) :: ne, ne1, ncell, nc1, ntot, nvt, b, v0, v1, g, c, i, j, k, ncut, nrw, niw, np, niv, nbd, gb !< Work.
   integer(I4P) :: eid(4)                                !< Edge ids.

   res = .false. ; if (present(resident)) res = resident
   ne = (nx + 1_I4P) * ny + nx * (ny + 1_I4P) ; ncell = nx * ny ; ntot = ncell * nbk
   ! host: node rows, walk and indexes per block, then concatenated (one padding element: a start element in bounds)
   allocate(pp(nbk), rso(nbk + 1), xo(nbk + 1), po(nbk + 1))
   rso(1) = 0_I4P ; xo(1) = 0_I4P ; po(1) = 0_I4P
   do b = 1_I4P, nbk
      v0 = vstart(b) ; v1 = vstart(b + 1) - 1_I4P ; gb = min(b, ngr)
      call node_rows(u(v0:v1), v(v0:v1), next(v0:v1), yg(:, gb), pp(b)%rstart, pp(b)%rseg)
      call cut_cells_2d_prepare(u(v0:v1), v(v0:v1), next(v0:v1), xg(:, gb), yg(:, gb), nx, ny, ne1, nc1, pp(b)%nxs, &
                                pp(b)%nps, pp(b)%xstart, pp(b)%xord, pp(b)%pstart, pp(b)%pord, pp(b)%xp, pp(b)%pab, pp(b)%ptg)
      rso(b + 1) = rso(b) + size(pp(b)%rseg, kind=I4P) ; xo(b + 1) = xo(b) + pp(b)%nxs ; po(b + 1) = po(b) + pp(b)%nps
   enddo
   nvt = vstart(nbk + 1) - 1_I4P
   allocate(ua(nvt + 1), va(nvt + 1), na(nvt + 1), rsa(ny + 2, nbk), rsg(rso(nbk + 1) + 1), xsa(ne + 1, nbk), &
            xoa(xo(nbk + 1) + 1), xpa(2, xo(nbk + 1) + 1), psa(ncell + 1, nbk), poa(po(nbk + 1) + 1), &
            paba(4, po(nbk + 1) + 1), ptga(2, po(nbk + 1) + 1))
   ua(1:nvt) = u(1:nvt) ; va(1:nvt) = v(1:nvt) ; na(1:nvt) = next(1:nvt)
   ua(nvt + 1) = 0._R8P ; va(nvt + 1) = 0._R8P ; na(nvt + 1) = 0_I4P
   rsg(rso(nbk + 1) + 1) = 0_I4P ; xoa(xo(nbk + 1) + 1) = 0_I4P ; xpa(:, xo(nbk + 1) + 1) = 0._R8P
   poa(po(nbk + 1) + 1) = 0_I4P ; paba(:, po(nbk + 1) + 1) = 0._R8P ; ptga(:, po(nbk + 1) + 1) = 0_I4P
   do b = 1_I4P, nbk
      rsa(:, b) = pp(b)%rstart ; rsg(rso(b) + 1:rso(b + 1)) = pp(b)%rseg
      xsa(:, b) = pp(b)%xstart ; psa(:, b) = pp(b)%pstart
      do k = 1_I4P, pp(b)%nxs
         xoa(xo(b) + k) = pp(b)%xord(k) ; xpa(:, xo(b) + k) = pp(b)%xp(:, k)
      enddo
      do k = 1_I4P, pp(b)%nps
         poa(po(b) + k) = pp(b)%pord(k) ; paba(:, po(b) + k) = pp(b)%pab(:, k) ; ptga(:, po(b) + k) = pp(b)%ptg(:, k)
      enddo
   enddo
   deallocate(pp)
   zme = 0_I4P ; zmom = 0._R8P ; zwmom = 0._R8P ; zbseg = 0._R8P ; zbwall = .false.
   allocate(node_in(nx + 1, ny + 1, nbk), closed(nx, ny, nbk), alpha(nx, ny, nbk), aperture(4, nx, ny, nbk), &
            centroid(2, nx, ny, nbk), wet_length(nx, ny, nbk), wet_centroid(2, nx, ny, nbk), wet_normal(2, nx, ny, nbk), &
            pieces(nx, ny, nbk))
   !$acc data copyin(xg, yg, vstart, ua, va, na, rsa, rso, rsg, xsa, xo, xoa, xpa, psa, po, poa, paba, ptga, zme, zmom, &
   !$acc&            zwmom, zbseg, zbwall)
   !$acc enter data create(node_in, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces)
   ! 1. node flags, one thread per node of every block
   !$acc parallel loop gang vector collapse(3) default(present) present(node_in)
   do b = 1_I4P, nbk
      do j = 1_I4P, ny + 1_I4P
         do i = 1_I4P, nx + 1_I4P
            node_in(i, j, b) = node_parity(xg(i, min(b, ngr)), yg(j, min(b, ngr)), rsa(j, b), rsa(j + 1, b) - 1_I4P, &
                                           rso(b + 1) - rso(b), rsg(rso(b) + 1), vstart(b + 1) - vstart(b), ua(vstart(b)), &
                                           va(vstart(b)), na(vstart(b)))
         enddo
      enddo
   enddo
   ! 2. classification and workspace sizing over the cells of all blocks
   call cut_cells_2d_classify_blocks_dev(nx, ny, ne, ncell, nbk, xsa, psa, node_in, ckind, ncut, cutc, rwo, iwo, nrw, niw)
   allocate(cnt(3, max(1_I4P, ncut)), rw(max(1_I4P, nrw)), iw(max(1_I4P, niw)))
   !$acc data create(rw, iw, cnt)
   ! 3. count over the cut cells of all blocks (each reads its block's part of the concatenated arrays)
   !$acc parallel loop gang vector default(present) private(eid) &
   !$acc& present(node_in, cutc, rwo, iwo, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces)
   do k = 1_I4P, ncut
      g = cutc(k) ; b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell ; gb = min(b, ngr)
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      call cut_cell_2d_count(xg(i, gb), xg(i + 1, gb), yg(j, gb), yg(j + 1, gb), node_in(i, j, b), node_in(i + 1, j, b), &
                             node_in(i + 1, j + 1, b), node_in(i, j + 1, b), eid, ne, xo(b + 1) - xo(b), xsa(1, b), &
                             xoa(xo(b) + 1), xpa(1, xo(b) + 1), psa(c + 1, b) - psa(c, b), psa(c, b) - 1_I4P, &
                             po(b + 1) - po(b), poa(po(b) + 1), paba(1, po(b) + 1), ptga(1, po(b) + 1), 0_I4P, zme, with_iv, &
                             .false., rwo(k + 1) - rwo(k), rw(rwo(k)), iwo(k + 1) - iwo(k), iw(iwo(k)), closed(i, j, b), &
                             alpha(i, j, b), aperture(1, i, j, b), centroid(1, i, j, b), wet_length(i, j, b), &
                             wet_centroid(1, i, j, b), wet_normal(1, i, j, b), cnt(3, k), cnt(1, k), cnt(2, k))
      pieces(i, j, b) = cnt(3, k)
   enddo
   if (.not. with_pieces) then
      ! 4. the quiet cells' measures only
      !$acc parallel loop gang vector default(present) &
      !$acc& present(ckind, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces)
      do g = 1_I4P, ntot
         if (ckind(g) > 0_I4P) cycle
         b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell ; gb = min(b, ngr)
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call quiet_cell_2d_measures(xg(i, gb), xg(i + 1, gb), yg(j, gb), yg(j + 1, gb), ckind(g) < 0_I4P, &
                                     closed(i, j, b), alpha(i, j, b), aperture(1, i, j, b), centroid(1, i, j, b), &
                                     wet_length(i, j, b), wet_centroid(1, i, j, b), wet_normal(1, i, j, b), pieces(i, j, b))
      enddo
      np = 0_I4P ; niv = 0_I4P
   else
      ! 4. prefix sums over the cells of all blocks (the stacked grid nx x (ny nbk): cell g of block b is its cell g),
      !    then the offsets of the blocks
      allocate(first_piece(nx, ny, nbk), ivo(ntot), bdo(ntot), bpc(nbk + 1), biv(nbk + 1))
      !$acc enter data create(first_piece)
      !$acc data create(ivo, bdo, bpc, biv)
      call cut_cells_2d_scan_dev(nx, ny * nbk, ntot, ckind, ncut, cnt, with_iv, .false., first_piece, ivo, bdo, np, niv, &
                                 nbd)
      !$acc parallel loop gang vector default(present) present(first_piece, ivo, bpc, biv)
      do b = 1_I4P, nbk
         bpc(b) = first_piece(1, 1, b) - 1_I4P ; biv(b) = ivo((b - 1_I4P) * ncell + 1_I4P)
      enddo
      !$acc update host(bpc, biv)
      bpc(nbk + 1) = np ; biv(nbk + 1) = niv
      allocate(pcell(2, np), parea(np), pcen(2, np), pwl(np), pwn(2, np), pap(4, np), ivs(np + 1), &
               ive(max(1_I4P, niv)), ivl(max(1_I4P, niv)), ivh(max(1_I4P, niv)), bds(np + 1))
      !$acc enter data create(pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
      !$acc data create(bds)
      ! 5. fill: the quiet cells, then the cut cells
      !$acc parallel loop gang vector default(present) private(eid) &
      !$acc& present(ckind, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces, first_piece, &
      !$acc&         ivo, bdo, pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh, bds)
      do g = 1_I4P, ntot
         if (ckind(g) > 0_I4P) cycle
         b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell ; gb = min(b, ngr)
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         call quiet_cell_2d_kernel(xg(i, gb), xg(i + 1, gb), yg(j, gb), yg(j + 1, gb), i, j, ckind(g) < 0_I4P, eid, &
                                   0_I4P, zme, with_iv, .false., closed(i, j, b), alpha(i, j, b), aperture(1, i, j, b), &
                                   centroid(1, i, j, b), wet_length(i, j, b), wet_centroid(1, i, j, b), &
                                   wet_normal(1, i, j, b), pieces(i, j, b), first_piece(i, j, b) - 1_I4P, ivo(g), bdo(g), np, &
                                   pcell, parea, pcen, pwl, pwn, pap, zmom, zwmom, ivs, niv, ive, ivl, ivh, bds, 0_I4P, zbseg, &
                                   zbwall)
      enddo
      !$acc parallel loop gang vector default(present) private(eid) &
      !$acc& present(cutc, rwo, iwo, first_piece, ivo, bdo, pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh, bds)
      do k = 1_I4P, ncut
         g = cutc(k) ; b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell ; gb = min(b, ngr)
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         call cut_cell_2d_fill(xg(i, gb), xg(i + 1, gb), yg(j, gb), yg(j + 1, gb), i, j, eid, 0_I4P, zme, with_iv, .false., &
                               rwo(k + 1) - rwo(k), rw(rwo(k)), iwo(k + 1) - iwo(k), iw(iwo(k)), first_piece(i, j, b) - 1_I4P, &
                               ivo(g), bdo(g), np, pcell, parea, pcen, pwl, pwn, pap, zmom, zwmom, ivs, niv, &
                               ive, ivl, ivh, bds, 0_I4P, zbseg, zbwall)
      enddo
      !$acc serial default(present) present(ivs)
      ivs(np + 1) = niv + 1_I4P
      !$acc end serial
      !$acc end data
      !$acc end data
   endif
   !$acc end data
   !$acc exit data delete(ckind, cutc, rwo, iwo)
   !$acc end data
   if (.not. res) then ! copy back the outputs
      !$acc exit data copyout(node_in, closed, alpha, aperture, centroid, wet_length, wet_centroid, wet_normal, pieces)
      if (with_pieces) then
         !$acc exit data copyout(first_piece, pcell, parea, pcen, pwl, pwn, pap, ivs, ive, ivl, ivh)
      endif
   endif
   cb%nx = nx ; cb%ny = ny ; cb%nb = nbk ; cb%with_pieces = with_pieces
   call move_alloc(node_in, cb%node_in) ; call move_alloc(closed, cb%closed) ; call move_alloc(alpha, cb%alpha)
   call move_alloc(aperture, cb%aperture) ; call move_alloc(centroid, cb%centroid) ; call move_alloc(wet_length, cb%wet_length)
   call move_alloc(wet_centroid, cb%wet_centroid) ; call move_alloc(wet_normal, cb%wet_normal)
   call move_alloc(pieces, cb%pieces)
   if (with_pieces) then
      cb%np = np
      call move_alloc(first_piece, cb%first_piece) ; call move_alloc(bpc, cb%block_piece) ; call move_alloc(biv, cb%block_iv)
      call move_alloc(pcell, cb%piece_cell) ; call move_alloc(parea, cb%piece_area) ; call move_alloc(pcen, cb%piece_centroid)
      call move_alloc(pwl, cb%piece_wet_length) ; call move_alloc(pwn, cb%piece_wet_normal)
      call move_alloc(pap, cb%piece_aperture) ; call move_alloc(ivs, cb%piece_iv_start)
      if (niv > 0_I4P) then
         call move_alloc(ive, cb%piece_iv_edge) ; call move_alloc(ivl, cb%piece_iv_lo) ; call move_alloc(ivh, cb%piece_iv_hi)
      else ! the one-element placeholders go; the empty outputs are mapped too, so every output is present when resident
         if (res) then
            !$acc exit data delete(ive, ivl, ivh)
         endif
         allocate(cb%piece_iv_edge(0), cb%piece_iv_lo(0), cb%piece_iv_hi(0))
         if (res) then
            !$acc enter data create(cb%piece_iv_edge, cb%piece_iv_lo, cb%piece_iv_hi)
         endif
      endif
      if (res) then
         !$acc enter data copyin(cb%block_piece, cb%block_iv)
      endif
   endif
   cb%on_device = res
   endsubroutine cut_cells_2d_blocks_core

   subroutine cut_cells_2d_classify_blocks_dev(nx, ny, ne, ncell, nbk, xstart, pstart, node_in, ckind, ncut, cutc, rwo, iwo, &
                                               nrw_tot, niw_tot)
   !< `cut_cells_2d_classify_dev` over the cells g = (b - 1) ncell + c of nbk blocks of one grid (issue #100 step 3b),
   !< moments order 0: cut cells numbered in the order of g, workspace slices in that order. xstart, pstart and node_in
   !< must be present on the device; ckind, cutc, rwo and iwo are left mapped (filled on the device only): the caller
   !< deletes them. Only ncut and the two workspace totals come back. Built in locals moved into the dummies (#100 6a).
   integer(I4P),              intent(in)  :: nx, ny, ne, ncell, nbk          !< Cells, edges per block, blocks.
   integer(I4P),              intent(in)  :: xstart(ne + 1, nbk)             !< Crossings by edge of each block.
   integer(I4P),              intent(in)  :: pstart(ncell + 1, nbk)          !< Wall pieces by cell of each block.
   logical,                   intent(in)  :: node_in(nx + 1, ny + 1, nbk)    !< Node inside flags.
   integer(I4P), allocatable, intent(out) :: ckind(:)                        !< Cell kinds [ncell nbk].
   integer(I4P),              intent(out) :: ncut                            !< Cut cells.
   integer(I4P), allocatable, intent(out) :: cutc(:)                         !< Cell g of each cut cell [ncut].
   integer(I4P), allocatable, intent(out) :: rwo(:), iwo(:)                  !< Workspace slices [ncut + 1].
   integer(I4P),              intent(out) :: nrw_tot, niw_tot                !< Workspace pool sizes.
   integer(I4P), allocatable :: tot(:,:)                                     !< Chunk totals, then offsets [3, nch].
   integer(I4P), allocatable :: ck(:), cu(:), ro(:), io(:)                   !< Kinds, cut cells, offsets (built here).
   integer(I4P) :: ntot, nch, ch, g, b, c, i, j, k, e, nxc, nrw, niw, a1, a2, a3, eid(4), sums(3) !< Work.

   ntot = ncell * nbk
   nch = (ntot + SCAN_SEG - 1_I4P) / SCAN_SEG
   allocate(ck(ntot), tot(3, nch))
   !$acc enter data create(ck)
   !$acc data create(tot)
   !$acc parallel loop gang vector default(present) private(eid) present(xstart, pstart, node_in, ck, tot)
   do ch = 1_I4P, nch
      a1 = 0_I4P ; a2 = 0_I4P ; a3 = 0_I4P
      !$acc loop seq
      do g = (ch - 1_I4P) * SCAN_SEG + 1_I4P, min(ch * SCAN_SEG, ntot)
         b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         if (cell_is_quiet(pstart(c + 1, b) - pstart(c, b), eid, ne, xstart(1, b), node_in(i, j, b), &
                           node_in(i + 1, j, b), node_in(i + 1, j + 1, b), node_in(i, j + 1, b))) then
            ck(g) = 0_I4P ; if (node_in(i, j, b)) ck(g) = -1_I4P
         else
            nxc = 0_I4P
            do e = 1_I4P, 4_I4P
               nxc = nxc + xstart(eid(e) + 1, b) - xstart(eid(e), b)
            enddo
            call cut_cell_2d_sizes(nxc, pstart(c + 1, b) - pstart(c, b), 0_I4P, nrw, niw)
            ck(g) = 1_I4P ; a1 = a1 + 1_I4P ; a2 = a2 + nrw ; a3 = a3 + niw
         endif
      enddo
      tot(1, ch) = a1 ; tot(2, ch) = a2 ; tot(3, ch) = a3
   enddo
   call chunk_offsets_dev(nch, [0_I4P, 1_I4P, 1_I4P], tot, sums)
   ncut = sums(1) ; nrw_tot = sums(2) - 1_I4P ; niw_tot = sums(3) - 1_I4P
   allocate(cu(ncut), ro(ncut + 1), io(ncut + 1))
   !$acc enter data create(cu, ro, io)
   !$acc parallel loop gang vector default(present) private(eid) present(xstart, pstart, ck, cu, ro, io, tot)
   do ch = 1_I4P, nch
      k = tot(1, ch) ; a2 = tot(2, ch) ; a3 = tot(3, ch)
      !$acc loop seq
      do g = (ch - 1_I4P) * SCAN_SEG + 1_I4P, min(ch * SCAN_SEG, ntot)
         if (ck(g) <= 0_I4P) cycle
         b = (g - 1_I4P) / ncell + 1_I4P ; c = g - (b - 1_I4P) * ncell
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         nxc = 0_I4P
         do e = 1_I4P, 4_I4P
            nxc = nxc + xstart(eid(e) + 1, b) - xstart(eid(e), b)
         enddo
         call cut_cell_2d_sizes(nxc, pstart(c + 1, b) - pstart(c, b), 0_I4P, nrw, niw)
         k = k + 1_I4P ; ck(g) = k ; cu(k) = g ; ro(k) = a2 ; io(k) = a3
         a2 = a2 + nrw ; a3 = a3 + niw
      enddo
   enddo
   a2 = sums(2) ; a3 = sums(3)
   !$acc serial default(present) present(ro, io)
   ro(ncut + 1) = a2 ; io(ncut + 1) = a3
   !$acc end serial
   !$acc end data
   call move_alloc(ck, ckind) ; call move_alloc(cu, cutc) ; call move_alloc(ro, rwo) ; call move_alloc(io, iwo)
   endsubroutine cut_cells_2d_classify_blocks_dev

   subroutine cut_cells_2d_blocks_reset(cb)
   !< Release a resident multi-block result from the device, then reset it to the default (issue #100 step 3b).
   type(cut_cells_2d_blocks_t), intent(inout) :: cb    !< Cut cells.
   type(cut_cells_2d_blocks_t)                :: fresh !< Default.

   if (cb%on_device) call cb%exit_device
   cb = fresh
   endsubroutine cut_cells_2d_blocks_reset

   subroutine blocks_update_host(self)
   !< Copy the mapped arrays of a resident multi-block result back to the host (they stay mapped). No-op if not mapped.
   class(cut_cells_2d_blocks_t), intent(inout) :: self !< Cut cells.

   if (.not. self%on_device) return
   !$acc update host(self%node_in, self%closed, self%alpha, self%aperture, self%centroid, self%wet_length, &
   !$acc&            self%wet_centroid, self%wet_normal, self%pieces)
   if (self%with_pieces) then
      !$acc update host(self%first_piece, self%piece_cell, self%piece_area, self%piece_centroid, self%piece_wet_length, &
      !$acc&            self%piece_wet_normal, self%piece_aperture, self%piece_iv_start, self%piece_iv_edge, &
      !$acc&            self%piece_iv_lo, self%piece_iv_hi)
   endif
   endsubroutine blocks_update_host

   subroutine blocks_exit_device(self)
   !< Release the device copies of a resident multi-block result (the host arrays are kept as they are).
   class(cut_cells_2d_blocks_t), intent(inout) :: self !< Cut cells.

   if (.not. self%on_device) return
   !$acc exit data delete(self%node_in, self%closed, self%alpha, self%aperture, self%centroid, self%wet_length, &
   !$acc&                 self%wet_centroid, self%wet_normal, self%pieces)
   if (self%with_pieces) then
      !$acc exit data delete(self%first_piece, self%block_piece, self%block_iv, self%piece_cell, self%piece_area, &
      !$acc&                 self%piece_centroid, self%piece_wet_length, self%piece_wet_normal, self%piece_aperture, &
      !$acc&                 self%piece_iv_start, self%piece_iv_edge, self%piece_iv_lo, self%piece_iv_hi)
   endif
   self%on_device = .false.
   endsubroutine blocks_exit_device

   subroutine cut_cells_2d_reset(cc)
   !< Release a resident cut from the device (`exit_device`), then reset it to the default (issue #100 step 6a): what
   !< `intent(out)` did, without freeing arrays that are still mapped.
   type(cut_cells_2d_t), intent(inout) :: cc    !< Cut cells.
   type(cut_cells_2d_t)                :: fresh !< Default.

   if (cc%on_device) call cc%exit_device
   cc = fresh
   endsubroutine cut_cells_2d_reset

   subroutine take_l2(src, dst, n1, n2)
   !< dst = src (moved) if src has shape [n1, n2], else a new array of that shape (issue #100 step 6a).
   logical, allocatable, intent(inout) :: src(:,:) !< Old array.
   logical, allocatable, intent(out)   :: dst(:,:) !< New array.
   integer(I4P),         intent(in)    :: n1, n2   !< Shape.

   if (allocated(src)) then
      if (size(src, 1) == n1 .and. size(src, 2) == n2) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1, n2))
   endsubroutine take_l2

   subroutine take_i1(src, dst, n1)
   !< dst = src (moved) if src has size n1, else a new array.
   integer(I4P), allocatable, intent(inout) :: src(:) !< Old array.
   integer(I4P), allocatable, intent(out)   :: dst(:) !< New array.
   integer(I4P),              intent(in)    :: n1     !< Size.

   if (allocated(src)) then
      if (size(src, 1) == n1) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1))
   endsubroutine take_i1

   subroutine take_i2(src, dst, n1, n2)
   !< dst = src (moved) if src has shape [n1, n2], else a new array.
   integer(I4P), allocatable, intent(inout) :: src(:,:) !< Old array.
   integer(I4P), allocatable, intent(out)   :: dst(:,:) !< New array.
   integer(I4P),              intent(in)    :: n1, n2   !< Shape.

   if (allocated(src)) then
      if (size(src, 1) == n1 .and. size(src, 2) == n2) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1, n2))
   endsubroutine take_i2

   subroutine take_r1(src, dst, n1)
   !< dst = src (moved) if src has size n1, else a new array.
   real(R8P), allocatable, intent(inout) :: src(:) !< Old array.
   real(R8P), allocatable, intent(out)   :: dst(:) !< New array.
   integer(I4P),           intent(in)    :: n1     !< Size.

   if (allocated(src)) then
      if (size(src, 1) == n1) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1))
   endsubroutine take_r1

   subroutine take_r2(src, dst, n1, n2)
   !< dst = src (moved) if src has shape [n1, n2], else a new array.
   real(R8P), allocatable, intent(inout) :: src(:,:) !< Old array.
   real(R8P), allocatable, intent(out)   :: dst(:,:) !< New array.
   integer(I4P),           intent(in)    :: n1, n2   !< Shape.

   if (allocated(src)) then
      if (size(src, 1) == n1 .and. size(src, 2) == n2) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1, n2))
   endsubroutine take_r2

   subroutine take_r3(src, dst, n1, n2, n3)
   !< dst = src (moved) if src has shape [n1, n2, n3], else a new array.
   real(R8P), allocatable, intent(inout) :: src(:,:,:) !< Old array.
   real(R8P), allocatable, intent(out)   :: dst(:,:,:) !< New array.
   integer(I4P),           intent(in)    :: n1, n2, n3 !< Shape.

   if (allocated(src)) then
      if (size(src, 1) == n1 .and. size(src, 2) == n2 .and. size(src, 3) == n3) then
         call move_alloc(src, dst) ; return
      endif
   endif
   allocate(dst(n1, n2, n3))
   endsubroutine take_r3

   function host_key(cc) result(key)
   !< Host address of alpha (the key of a mapped cut; 0 if alpha is not allocated).
   type(cut_cells_2d_t), intent(in), target :: cc  !< Cut cells.
   integer(c_intptr_t)                      :: key !< Address.

   key = 0
   if (allocated(cc%alpha)) then
      if (size(cc%alpha) > 0) key = transfer(c_loc(cc%alpha), key)
   endif
   endfunction host_key

   subroutine mark_resident(cc)
   !< Flag the outputs of `cc` as mapped, keyed by the host address of alpha.
   type(cut_cells_2d_t), intent(inout) :: cc !< Cut cells.

   cc%on_device = .true. ; cc%device_key = host_key(cc)
   endsubroutine mark_resident

   logical function is_mapped(self)
   !< The flag is set and the arrays are the ones mapped (not a copy of a resident cut).
   type(cut_cells_2d_t), intent(in) :: self !< Cut cells.

   is_mapped = .false.
   if (self%on_device) is_mapped = self%device_key == host_key(self)
   endfunction is_mapped

   subroutine enter_device(self)
   !< Map the outputs of a host result to the device (copyin; issue #100 step 6a). Member arrays only, never `self`
   !< (#20). No-op if already mapped or if there are no outputs (status error, default cut). A copy of a resident cut is
   !< mapped as a new cut.
   class(cut_cells_2d_t), intent(inout) :: self !< Cut cells.

   if (is_mapped(self)) return
   self%on_device = .false. ; self%device_key = 0
   if (.not. allocated(self%alpha)) return
   !$acc enter data copyin(self%closed, self%alpha, self%aperture, self%centroid, self%wet_length, self%wet_centroid, &
   !$acc&                  self%wet_normal, self%pieces, self%first_piece, self%node_in, self%piece_cell, self%piece_area, &
   !$acc&                  self%piece_centroid, self%piece_wet_length, self%piece_wet_normal, self%piece_aperture, &
   !$acc&                  self%piece_iv_start, self%piece_iv_edge, self%piece_iv_lo, self%piece_iv_hi)
   call mark_resident(self)
   endsubroutine enter_device

   subroutine update_host(self)
   !< Copy the mapped outputs back to the host (issue #100 step 6a); they stay mapped. No-op if not mapped; on a copy of
   !< a resident cut a warning on stderr (its arrays are not on the device: update the original before copying it).
   class(cut_cells_2d_t), intent(inout) :: self !< Cut cells.

   if (.not. self%on_device) return
   if (.not. is_mapped(self)) then
      write(stderr, '(A)') 'fossil cut_cells_2d_t%update_host: a copy of a resident cut is not on the device; '// &
                           'call update_host on the original before copying it'
      self%on_device = .false. ; self%device_key = 0
      return
   endif
   !$acc update host(self%closed, self%alpha, self%aperture, self%centroid, self%wet_length, self%wet_centroid, &
   !$acc&            self%wet_normal, self%pieces, self%first_piece, self%node_in, self%piece_cell, self%piece_area, &
   !$acc&            self%piece_centroid, self%piece_wet_length, self%piece_wet_normal, self%piece_aperture, &
   !$acc&            self%piece_iv_start, self%piece_iv_edge, self%piece_iv_lo, self%piece_iv_hi)
   endsubroutine update_host

   subroutine exit_device(self)
   !< Release the device copies of the outputs (issue #100 step 6a); the host arrays are kept as they are (call
   !< `update_host` first to keep the device values). No-op if not mapped; a copy of a resident cut only drops the flag.
   class(cut_cells_2d_t), intent(inout) :: self !< Cut cells.

   if (.not. self%on_device) return
   if (is_mapped(self)) then
      !$acc exit data delete(self%closed, self%alpha, self%aperture, self%centroid, self%wet_length, self%wet_centroid, &
      !$acc&                 self%wet_normal, self%pieces, self%first_piece, self%node_in, self%piece_cell, self%piece_area, &
      !$acc&                 self%piece_centroid, self%piece_wet_length, self%piece_wet_normal, self%piece_aperture, &
      !$acc&                 self%piece_iv_start, self%piece_iv_edge, self%piece_iv_lo, self%piece_iv_hi)
   endif
   self%on_device = .false. ; self%device_key = 0
   endsubroutine exit_device

   subroutine cut_cells_2d_grid(x, y, nx, ny, status)
   !< Cells of the block of edges x(nx+1), y(ny+1); status 1 if there is no cell or the edges do not increase strictly.
   real(R8P),    intent(in)  :: x(:), y(:) !< Cell edges.
   integer(I4P), intent(out) :: nx, ny     !< Cells.
   integer(I4P), intent(out) :: status     !< 0 ok, 1 bad edges.

   status = 0_I4P
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   if (nx < 1 .or. ny < 1) then
      status = 1_I4P ; return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny))) then
      status = 1_I4P ; return
   endif
   endsubroutine cut_cells_2d_grid

   subroutine cut_cells_2d_prepare(u, v, next, x, y, nx, ny, ne, ncell, nxs, nps, xstart, xord, pstart, pord, xp, pab, ptg)
   !< Host preparation of the kernels: the segments walked through the grid (`walk_all`), the crossings by edge and the
   !< wall pieces by cell (`index_crossings`, `index_pieces`), as flat arrays.
   real(R8P),                 intent(in)  :: u(:), v(:)         !< Vertices.
   integer(I4P),              intent(in)  :: next(:)            !< Next vertex.
   real(R8P),                 intent(in)  :: x(:), y(:)         !< Cell edges.
   integer(I4P),              intent(in)  :: nx, ny             !< Cells.
   integer(I4P),              intent(out) :: ne, ncell          !< Edges, cells.
   integer(I4P),              intent(out) :: nxs, nps           !< Crossings, wall pieces.
   integer(I4P), allocatable, intent(out) :: xstart(:), xord(:) !< Crossings by edge.
   integer(I4P), allocatable, intent(out) :: pstart(:), pord(:) !< Wall pieces by cell.
   real(R8P),    allocatable, intent(out) :: xp(:,:)            !< Crossing points [2, nxs].
   real(R8P),    allocatable, intent(out) :: pab(:,:)           !< Wall piece end points (a, b) [4, nps].
   integer(I4P), allocatable, intent(out) :: ptg(:,:)           !< Wall piece end tags (ta, tb) [2, nps].
   type(xing_t),  allocatable :: xs(:)                          !< Crossings.
   type(piece_t), allocatable :: ps(:)                          !< Wall pieces.
   integer(I4P) :: k                                            !< Counter.

   call walk_all(u, v, next, x, y, xs, nxs, ps, nps)
   ne = (nx + 1_I4P) * ny + nx * (ny + 1_I4P) ; ncell = nx * ny
   call index_crossings(xs, nxs, ne, xstart, xord)
   call index_pieces(ps, nps, ncell, pstart, pord)
   allocate(xp(2, nxs), pab(4, nps), ptg(2, nps))
   do k = 1_I4P, nxs
      xp(:, k) = xs(k)%p
   enddo
   do k = 1_I4P, nps
      pab(1:2, k) = ps(k)%a ; pab(3:4, k) = ps(k)%b ; ptg(1, k) = ps(k)%ta ; ptg(2, k) = ps(k)%tb
   enddo
   endsubroutine cut_cells_2d_prepare

   subroutine cut_cells_2d_classify(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ckind, ncut, cutc, rwo, iwo)
   !< Classify the cells (quiet: `cell_is_quiet`; -1 solid, 0 fluid, k > 0 the k-th cut cell) and give each cut cell its
   !< workspace slices rwo(k):rwo(k+1)-1, iwo(k):iwo(k+1)-1 (`cut_cell_2d_sizes`).
   integer(I4P),              intent(in)  :: nx, ny, ne, ncell  !< Cells, edges.
   integer(I4P),              intent(in)  :: xstart(ne + 1)     !< Crossings by edge.
   integer(I4P),              intent(in)  :: pstart(ncell + 1)  !< Wall pieces by cell.
   logical,                   intent(in)  :: node_in(nx + 1, ny + 1) !< Node inside flags.
   integer(I4P),              intent(in)  :: nm                 !< Moments.
   integer(I4P), allocatable, intent(out) :: ckind(:)           !< Cell kinds [ncell].
   integer(I4P),              intent(out) :: ncut               !< Cut cells.
   integer(I4P), allocatable, intent(out) :: cutc(:)            !< Cell of each cut cell [ncut].
   integer(I4P), allocatable, intent(out) :: rwo(:), iwo(:)     !< Workspace slices [ncut + 1].
   integer(I4P) :: c, i, j, k, e, nxc, nrw, niw, eid(4)         !< Work.

   allocate(ckind(ncell)) ; ncut = 0_I4P
   do c = 1_I4P, ncell ! a cell without wall pieces or crossings takes the short path (issue #70)
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      if (cell_is_quiet(pstart(c + 1) - pstart(c), eid, ne, xstart, node_in(i, j), node_in(i + 1, j), &
                        node_in(i + 1, j + 1), node_in(i, j + 1))) then
         ckind(c) = 0_I4P ; if (node_in(i, j)) ckind(c) = -1_I4P
      else
         ncut = ncut + 1_I4P ; ckind(c) = ncut
      endif
   enddo
   allocate(cutc(ncut), rwo(ncut + 1), iwo(ncut + 1))
   rwo(1) = 1_I4P ; iwo(1) = 1_I4P
   do c = 1_I4P, ncell
      k = ckind(c)
      if (k <= 0_I4P) cycle
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call cell_edge_ids(nx, ny, i, j, eid)
      nxc = 0_I4P
      do e = 1_I4P, 4_I4P
         nxc = nxc + xstart(eid(e) + 1) - xstart(eid(e))
      enddo
      call cut_cell_2d_sizes(nxc, pstart(c + 1) - pstart(c), nm, nrw, niw)
      cutc(k) = c ; rwo(k + 1) = rwo(k) + nrw ; iwo(k + 1) = iwo(k) + niw
   enddo
   endsubroutine cut_cells_2d_classify

   subroutine cut_cells_2d_scan(nx, ny, ncell, ckind, ncut, cnt, with_iv, with_bnd, first_piece, ivo, bdo, np, niv, nbd)
   !< Prefix sums over the cells in order (issue #100): first piece, intervals and boundary segments before each cell,
   !< and the totals. A quiet fluid cell has one piece, 4 intervals and 4 boundary segments; cut cell k has cnt(3, k)
   !< pieces, cnt(1, k) intervals and cnt(2, k) boundary segments.
   integer(I4P), intent(in)  :: nx, ny, ncell             !< Cells.
   integer(I4P), intent(in)  :: ckind(ncell)              !< Cell kinds.
   integer(I4P), intent(in)  :: ncut                      !< Cut cells.
   integer(I4P), intent(in)  :: cnt(:,:)                  !< Counts of the cut cells [3, >= ncut].
   logical,      intent(in)  :: with_iv, with_bnd         !< Intervals, boundary segments recorded.
   integer(I4P), intent(out) :: first_piece(nx, ny)       !< First piece of each cell.
   integer(I4P), intent(out) :: ivo(ncell), bdo(ncell)    !< Intervals, boundary segments before each cell.
   integer(I4P), intent(out) :: np, niv, nbd              !< Pieces, intervals, boundary segments.
   integer(I4P) :: c, i, j, k                             !< Work.

   np = 0_I4P ; niv = 0_I4P ; nbd = 0_I4P
   do c = 1_I4P, ncell
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      first_piece(i, j) = np + 1_I4P ; ivo(c) = niv ; bdo(c) = nbd
      k = ckind(c)
      if (k == 0_I4P) then
         np = np + 1_I4P
         if (with_iv) niv = niv + 4_I4P
         if (with_bnd) nbd = nbd + 4_I4P
      elseif (k > 0_I4P) then
         np = np + cnt(3, k) ; niv = niv + cnt(1, k) ; nbd = nbd + cnt(2, k)
      endif
   enddo
   endsubroutine cut_cells_2d_scan

   subroutine cut_cells_2d_classify_dev(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ckind, ncut, cutc, rwo, iwo, nrw_tot, &
                                        niw_tot, seg)
   !< `cut_cells_2d_classify` on the device (issue #100 step 6a): the same ckind, ncut, cutc, rwo, iwo, integer for
   !< integer. xstart, pstart and node_in must be present on the device; ckind, cutc, rwo and iwo are allocated here and
   !< left mapped (filled on the device; their host copies are not written): the caller deletes them. They are built in
   !< local arrays mapped by `enter data` and moved into the dummies (`move_alloc` keeps the host addresses): entering
   !< data on an allocatable dummy leaks a mapping of its descriptor with nvfortran 26.9 (later calls fail with
   !< 'partially present'). Only ncut and the two workspace totals come back (nrw_tot = rwo(ncut + 1) - 1, niw_tot =
   !< iwo(ncut + 1) - 1).
   !<
   !< Scan over chunks of `seg` cells in cell order (default SCAN_SEG): each chunk counts its cut cells and their
   !< workspace sizes (one thread per chunk), `chunk_offsets_dev` turns the chunk totals into exclusive offsets, and
   !< each chunk then numbers its cut cells from its offset. Integer sums: any chunk size gives the same integers, the
   !< numbering is the cell order of the host's. Built without OpenACC the loops run on the host in order.
   integer(I4P),              intent(in)           :: nx, ny, ne, ncell      !< Cells, edges.
   integer(I4P),              intent(in)           :: xstart(ne + 1)         !< Crossings by edge.
   integer(I4P),              intent(in)           :: pstart(ncell + 1)      !< Wall pieces by cell.
   logical,                   intent(in)           :: node_in(nx + 1, ny + 1) !< Node inside flags.
   integer(I4P),              intent(in)           :: nm                     !< Moments.
   integer(I4P), allocatable, intent(out)          :: ckind(:)               !< Cell kinds [ncell].
   integer(I4P),              intent(out)          :: ncut                   !< Cut cells.
   integer(I4P), allocatable, intent(out)          :: cutc(:)                !< Cell of each cut cell [ncut].
   integer(I4P), allocatable, intent(out)          :: rwo(:), iwo(:)         !< Workspace slices [ncut + 1].
   integer(I4P),              intent(out)          :: nrw_tot, niw_tot       !< Workspace pool sizes.
   integer(I4P),              intent(in), optional :: seg                    !< Cells per chunk (tests).
   integer(I4P), allocatable :: tot(:,:)                                     !< Chunk totals, then offsets [3, nch].
   integer(I4P), allocatable :: ck(:), cu(:), ro(:), io(:)                   !< Kinds, cut cells, offsets (built here).
   integer(I4P) :: sz, nch, ch, c, i, j, k, e, nxc, nrw, niw, a1, a2, a3, eid(4), sums(3) !< Work.

   sz = SCAN_SEG ; if (present(seg)) sz = max(1_I4P, seg)
   nch = (ncell + sz - 1_I4P) / sz
   allocate(ck(ncell), tot(3, nch))
   !$acc enter data create(ck)
   !$acc data create(tot)
   ! 1. per chunk: the quiet cells' kinds, the number of cut cells and their workspace sizes
   !$acc parallel loop gang vector default(present) private(eid) present(xstart, pstart, node_in, ck, tot)
   do ch = 1_I4P, nch
      a1 = 0_I4P ; a2 = 0_I4P ; a3 = 0_I4P
      !$acc loop seq
      do c = (ch - 1_I4P) * sz + 1_I4P, min(ch * sz, ncell)
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         if (cell_is_quiet(pstart(c + 1) - pstart(c), eid, ne, xstart, node_in(i, j), node_in(i + 1, j), &
                           node_in(i + 1, j + 1), node_in(i, j + 1))) then
            ck(c) = 0_I4P ; if (node_in(i, j)) ck(c) = -1_I4P
         else
            nxc = 0_I4P
            do e = 1_I4P, 4_I4P
               nxc = nxc + xstart(eid(e) + 1) - xstart(eid(e))
            enddo
            call cut_cell_2d_sizes(nxc, pstart(c + 1) - pstart(c), nm, nrw, niw)
            ck(c) = 1_I4P ; a1 = a1 + 1_I4P ; a2 = a2 + nrw ; a3 = a3 + niw
         endif
      enddo
      tot(1, ch) = a1 ; tot(2, ch) = a2 ; tot(3, ch) = a3
   enddo
   ! 2. exclusive offsets of the chunks
   call chunk_offsets_dev(nch, [0_I4P, 1_I4P, 1_I4P], tot, sums)
   ncut = sums(1) ; nrw_tot = sums(2) - 1_I4P ; niw_tot = sums(3) - 1_I4P
   allocate(cu(ncut), ro(ncut + 1), io(ncut + 1))
   !$acc enter data create(cu, ro, io)
   ! 3. per chunk: number the cut cells from the chunk's offsets
   !$acc parallel loop gang vector default(present) private(eid) present(xstart, pstart, ck, cu, ro, io, tot)
   do ch = 1_I4P, nch
      k = tot(1, ch) ; a2 = tot(2, ch) ; a3 = tot(3, ch)
      !$acc loop seq
      do c = (ch - 1_I4P) * sz + 1_I4P, min(ch * sz, ncell)
         if (ck(c) <= 0_I4P) cycle
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         call cell_edge_ids(nx, ny, i, j, eid)
         nxc = 0_I4P
         do e = 1_I4P, 4_I4P
            nxc = nxc + xstart(eid(e) + 1) - xstart(eid(e))
         enddo
         call cut_cell_2d_sizes(nxc, pstart(c + 1) - pstart(c), nm, nrw, niw)
         k = k + 1_I4P ; ck(c) = k ; cu(k) = c ; ro(k) = a2 ; io(k) = a3
         a2 = a2 + nrw ; a3 = a3 + niw
      enddo
   enddo
   a2 = sums(2) ; a3 = sums(3)
   !$acc serial default(present) present(ro, io)
   ro(ncut + 1) = a2 ; io(ncut + 1) = a3
   !$acc end serial
   !$acc end data
   call move_alloc(ck, ckind) ; call move_alloc(cu, cutc) ; call move_alloc(ro, rwo) ; call move_alloc(io, iwo)
   endsubroutine cut_cells_2d_classify_dev

   subroutine cut_cells_2d_scan_dev(nx, ny, ncell, ckind, ncut, cnt, with_iv, with_bnd, first_piece, ivo, bdo, np, niv, &
                                    nbd, seg)
   !< `cut_cells_2d_scan` on the device (issue #100 step 6a): the same first_piece, ivo, bdo and totals, integer for
   !< integer. Every array must be present on the device; the outputs are written there only (their host copies are not).
   !< Only the totals come back. The two-level chunk scan of `cut_cells_2d_classify_dev`.
   integer(I4P), intent(in)           :: nx, ny, ncell          !< Cells.
   integer(I4P), intent(in)           :: ckind(ncell)           !< Cell kinds.
   integer(I4P), intent(in)           :: ncut                   !< Cut cells.
   integer(I4P), intent(in)           :: cnt(:,:)               !< Counts of the cut cells [3, >= ncut].
   logical,      intent(in)           :: with_iv, with_bnd      !< Intervals, boundary segments recorded.
   integer(I4P), intent(out)          :: first_piece(nx, ny)    !< First piece of each cell.
   integer(I4P), intent(out)          :: ivo(ncell), bdo(ncell) !< Intervals, boundary segments before each cell.
   integer(I4P), intent(out)          :: np, niv, nbd           !< Pieces, intervals, boundary segments.
   integer(I4P), intent(in), optional :: seg                    !< Cells per chunk (tests).
   integer(I4P), allocatable :: tot(:,:)                        !< Chunk totals, then offsets [3, nch].
   integer(I4P) :: sz, nch, ch, c, i, j, k, a1, a2, a3, qi, qb, sums(3) !< Work.

   sz = SCAN_SEG ; if (present(seg)) sz = max(1_I4P, seg)
   nch = (ncell + sz - 1_I4P) / sz
   qi = 0_I4P ; if (with_iv) qi = 4_I4P
   qb = 0_I4P ; if (with_bnd) qb = 4_I4P
   allocate(tot(3, nch))
   !$acc data create(tot)
   !$acc parallel loop gang vector default(present) present(ckind, cnt, tot)
   do ch = 1_I4P, nch
      a1 = 0_I4P ; a2 = 0_I4P ; a3 = 0_I4P
      !$acc loop seq
      do c = (ch - 1_I4P) * sz + 1_I4P, min(ch * sz, ncell)
         k = ckind(c)
         if (k == 0_I4P) then
            a1 = a1 + 1_I4P ; a2 = a2 + qi ; a3 = a3 + qb
         elseif (k > 0_I4P) then
            a1 = a1 + cnt(3, k) ; a2 = a2 + cnt(1, k) ; a3 = a3 + cnt(2, k)
         endif
      enddo
      tot(1, ch) = a1 ; tot(2, ch) = a2 ; tot(3, ch) = a3
   enddo
   call chunk_offsets_dev(nch, [0_I4P, 0_I4P, 0_I4P], tot, sums)
   np = sums(1) ; niv = sums(2) ; nbd = sums(3)
   !$acc parallel loop gang vector default(present) present(ckind, cnt, tot, first_piece, ivo, bdo)
   do ch = 1_I4P, nch
      a1 = tot(1, ch) ; a2 = tot(2, ch) ; a3 = tot(3, ch)
      !$acc loop seq
      do c = (ch - 1_I4P) * sz + 1_I4P, min(ch * sz, ncell)
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         first_piece(i, j) = a1 + 1_I4P ; ivo(c) = a2 ; bdo(c) = a3
         k = ckind(c)
         if (k == 0_I4P) then
            a1 = a1 + 1_I4P ; a2 = a2 + qi ; a3 = a3 + qb
         elseif (k > 0_I4P) then
            a1 = a1 + cnt(3, k) ; a2 = a2 + cnt(1, k) ; a3 = a3 + cnt(2, k)
         endif
      enddo
   enddo
   !$acc end data
   endsubroutine cut_cells_2d_scan_dev

   subroutine chunk_offsets_dev(nch, base, tot, sums)
   !< Exclusive prefix sums of the chunk totals tot(1:3, 1:nch) in place, from `base`, on the device (issue #100 step
   !< 6a); sums = base + the totals, on the host. tot must be present. Two levels over super-chunks of SCAN_SUP chunks
   !< (one thread each), one thread over the super-chunk totals, then one thread per chunk adds its super-chunk's
   !< offset: integer sums, the integers of a sequential scan.
   integer(I4P), intent(in)    :: nch          !< Chunks.
   integer(I4P), intent(in)    :: base(3)      !< Start values.
   integer(I4P), intent(inout) :: tot(3, nch)  !< Chunk totals in, exclusive offsets out.
   integer(I4P), intent(out)   :: sums(3)      !< base + the totals.
   integer(I4P), allocatable   :: sup(:,:)     !< Super-chunk totals, then offsets.
   integer(I4P) :: nsu, q, ch, a1, a2, a3, b1, b2, b3, t1, t2, t3, sm(3) !< Work.

   nsu = (nch + SCAN_SUP - 1_I4P) / SCAN_SUP
   b1 = base(1) ; b2 = base(2) ; b3 = base(3)
   allocate(sup(3, nsu))
   !$acc data create(sup, sm)
   !$acc parallel loop gang vector default(present) present(tot, sup)
   do q = 1_I4P, nsu
      a1 = 0_I4P ; a2 = 0_I4P ; a3 = 0_I4P
      !$acc loop seq
      do ch = (q - 1_I4P) * SCAN_SUP + 1_I4P, min(q * SCAN_SUP, nch)
         t1 = tot(1, ch) ; t2 = tot(2, ch) ; t3 = tot(3, ch)
         tot(1, ch) = a1 ; tot(2, ch) = a2 ; tot(3, ch) = a3
         a1 = a1 + t1 ; a2 = a2 + t2 ; a3 = a3 + t3
      enddo
      sup(1, q) = a1 ; sup(2, q) = a2 ; sup(3, q) = a3
   enddo
   !$acc serial default(present) present(sup, sm)
   a1 = b1 ; a2 = b2 ; a3 = b3
   do q = 1_I4P, nsu
      t1 = sup(1, q) ; t2 = sup(2, q) ; t3 = sup(3, q)
      sup(1, q) = a1 ; sup(2, q) = a2 ; sup(3, q) = a3
      a1 = a1 + t1 ; a2 = a2 + t2 ; a3 = a3 + t3
   enddo
   sm(1) = a1 ; sm(2) = a2 ; sm(3) = a3
   !$acc end serial
   !$acc parallel loop gang vector default(present) present(tot, sup)
   do ch = 1_I4P, nch
      q = (ch - 1_I4P) / SCAN_SUP + 1_I4P
      tot(1, ch) = tot(1, ch) + sup(1, q) ; tot(2, ch) = tot(2, ch) + sup(2, q) ; tot(3, ch) = tot(3, ch) + sup(3, q)
   enddo
   !$acc update host(sm)
   !$acc end data
   sums = sm
   endsubroutine chunk_offsets_dev

   pure subroutine cell_edge_ids(nx, ny, i, j, eid)
   !< Edge ids of cell (i, j) counter-clockwise: bottom, right, top, left. Vertical edge (l, k) = l + (k - 1) (nx + 1),
   !< horizontal edge (k, l) = (nx + 1) ny + k + (l - 1) nx.
   !$acc routine seq
   integer(I4P), intent(in)  :: nx, ny, i, j !< Cells, cell.
   integer(I4P), intent(out) :: eid(4)       !< Edge ids.

   eid(1) = (nx + 1_I4P) * ny + i + (j - 1_I4P) * nx
   eid(2) = i + 1_I4P + (j - 1_I4P) * (nx + 1_I4P)
   eid(3) = (nx + 1_I4P) * ny + i + j * nx
   eid(4) = i + (j - 1_I4P) * (nx + 1_I4P)
   endsubroutine cell_edge_ids

   pure logical function cell_is_quiet(nb, eid, ne, xstart, f00, f10, f11, f01)
   !< The cell has no wall piece, no crossing on its edges and four equal corner flags.
   !$acc routine seq
   integer(I4P), intent(in) :: nb                 !< Wall pieces in the cell.
   integer(I4P), intent(in) :: eid(4)             !< Edge ids.
   integer(I4P), intent(in) :: ne                 !< Edges of the block.
   integer(I4P), intent(in) :: xstart(ne + 1)     !< Crossings by edge.
   logical,      intent(in) :: f00, f10, f11, f01 !< Corner flags (i, j), (i+1, j), (i+1, j+1), (i, j+1).
   integer(I4P)             :: e                  !< Counter.

   cell_is_quiet = nb == 0_I4P
   if (.not. cell_is_quiet) return
   do e = 1, 4
      if (xstart(eid(e) + 1) /= xstart(eid(e))) then
         cell_is_quiet = .false. ; return
      endif
   enddo
   cell_is_quiet = (f10 .eqv. f00) .and. (f01 .eqv. f00) .and. (f11 .eqv. f00)
   endfunction cell_is_quiet

   pure subroutine cut_cell_2d_sizes(nxc, nb, nm, nrw, niw)
   !< Workspace of one cut cell with nxc crossings on its four edges, nb wall pieces and nm moments (issue #100).
   !<
   !< With S = 4 + nxc + nb, the most sub-intervals the cell can hold (nxc + 4 on its boundary, then one hole slot per
   !< closed loop, at most one loop per wall piece):
   !<    nrw = S (12 + 2 nm) + nb (13 + 2 nm) + 2 nm + 1 reals,  niw = NHDR + 4 S + 2 nxc + 14 nb + 1 integers
   !< (layout in `cut_cell_2d_layout`; the + 1 keeps every slice start inside the slice when the last parts are empty).
   !< Every count of the cell is bounded by its part (sub-intervals <= S, crossing slots <= nxc, loops, holes, islands,
   !< chains, loop pieces <= nb), so the overflow paths of the fixed 512-column arrays (issue #95) cannot occur.
   !$acc routine seq
   integer(I4P), intent(in)  :: nxc, nb, nm !< Crossings, wall pieces, moments.
   integer(I4P), intent(out) :: nrw, niw    !< Reals, integers.
   integer(I4P)              :: s           !< Sub-interval slots.

   s = 4_I4P + nxc + nb
   nrw = s * (12_I4P + 2_I4P * nm) + nb * (13_I4P + 2_I4P * nm) + 2_I4P * nm + 1_I4P
   niw = NHDR + 4_I4P * s + 2_I4P * nxc + 14_I4P * nb + 1_I4P
   endsubroutine cut_cell_2d_sizes

   pure subroutine cut_cell_2d_layout(nxc, nb, nm, ro, io)
   !< Start of each part of a cut cell's workspace (sizes in `cut_cell_2d_sizes`; S = 4 + nxc + nb).
   !< Reals: sa(2, S), sb(2, S), acc(6, S), wn(2, S), accm(nm, S), wallm(nm, S), hole area(nb), island Green(3, nb),
   !< island wet(3, nb), island normal(2, nb), island moments(nm, nb), island wall moments(nm, nb), the cell's wall
   !< pieces (4, nb), chain moments(nm), chain wall moments(nm).
   !< Integers (after the NHDR header): fluid flag(S), end crossing(S), edge(S), union-find parent(S), crossing slots
   !< (2, nxc), visited(nb), loop pieces(nb), hole root(nb), hole first(nb), hole length(nb), island first(nb), island
   !< length(nb), island traversal(nb), traversal kind(nb), start(nb), length(nb), slot(nb), wall piece tags(2, nb).
   !$acc routine seq
   integer(I4P), intent(in)  :: nxc, nb, nm  !< Crossings, wall pieces, moments.
   integer(I4P), intent(out) :: ro(15)       !< Real parts.
   integer(I4P), intent(out) :: io(18)       !< Integer parts.
   integer(I4P)              :: s, k         !< Work.

   s = 4_I4P + nxc + nb
   ro(1) = 1_I4P
   ro(2) = ro(1) + 2_I4P * s ; ro(3) = ro(2) + 2_I4P * s ; ro(4) = ro(3) + 6_I4P * s ; ro(5) = ro(4) + 2_I4P * s
   ro(6) = ro(5) + nm * s ; ro(7) = ro(6) + nm * s ; ro(8) = ro(7) + nb ; ro(9) = ro(8) + 3_I4P * nb
   ro(10) = ro(9) + 3_I4P * nb ; ro(11) = ro(10) + 2_I4P * nb ; ro(12) = ro(11) + nm * nb ; ro(13) = ro(12) + nm * nb
   ro(14) = ro(13) + 4_I4P * nb ; ro(15) = ro(14) + nm
   io(1) = NHDR + 1_I4P
   do k = 2, 4
      io(k) = io(k - 1) + s
   enddo
   io(5) = io(4) + s ; io(6) = io(5) + 2_I4P * nxc
   do k = 7, 18
      io(k) = io(k - 1) + nb
   enddo
   endsubroutine cut_cell_2d_layout

   subroutine cut_cell_2d_count(x0, x1, y0, y1, f00, f10, f11, f01, eid, ne, nxs, xstart, xord, xp, nb, pofs, nps, pord, &
                                pab, ptg, nm, me, with_iv, with_bnd, nrw, rw, niw, iw, closed, alpha, aperture, centroid, &
                                wet_length, wet_centroid, wet_normal, npieces, nivc, nbdc)
   !< Kernel, count pass (issue #100): measures and pieces of the cut cell [x0, x1] x [y0, y1]; its cell outputs are
   !< written, its pieces stay in its workspace (rw, iw) for `cut_cell_2d_fill`; npieces, nivc, nbdc are the numbers
   !< of pieces, fluid edge intervals and boundary segments that the fill writes.
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1         !< Cell.
   logical,      intent(in)    :: f00, f10, f11, f01     !< Corner flags (i, j), (i+1, j), (i+1, j+1), (i, j+1).
   integer(I4P), intent(in)    :: eid(4)                 !< Edge ids (bottom, right, top, left).
   integer(I4P), intent(in)    :: ne, nxs                !< Edges of the block, crossings.
   integer(I4P), intent(in)    :: xstart(ne + 1)         !< Crossings by edge.
   integer(I4P), intent(in)    :: xord(nxs)              !< Crossings ordered along each edge.
   real(R8P),    intent(in)    :: xp(2, nxs)             !< Crossing points.
   integer(I4P), intent(in)    :: nb, pofs, nps          !< Wall pieces of the cell, first - 1 in pord, wall pieces.
   integer(I4P), intent(in)    :: pord(nps)              !< Wall pieces by cell.
   real(R8P),    intent(in)    :: pab(4, nps)            !< Wall piece end points.
   integer(I4P), intent(in)    :: ptg(2, nps)            !< Wall piece end tags.
   integer(I4P), intent(in)    :: nm                     !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)              !< Moment exponents.
   logical,      intent(in)    :: with_iv, with_bnd      !< Count the intervals, the boundary segments.
   integer(I4P), intent(in)    :: nrw, niw               !< Workspace sizes (`cut_cell_2d_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)                !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)                !< Integer workspace.
   logical,      intent(out)   :: closed                 !< The cell boundary closed.
   real(R8P),    intent(out)   :: alpha, aperture(4), centroid(2), wet_length, wet_centroid(2), wet_normal(2) !< Measures.
   integer(I4P), intent(out)   :: npieces, nivc, nbdc    !< Pieces, intervals, boundary segments.
   integer(I4P)                :: nxc, e, s, ro(15), io(18) !< Work.

   nxc = 0_I4P
   do e = 1, 4
      nxc = nxc + xstart(eid(e) + 1) - xstart(eid(e))
   enddo
   call cut_cell_2d_layout(nxc, nb, nm, ro, io)
   s = 4_I4P + nxc + nb
   call count_core(x0, x1, y0, y1, f00, f10, f11, f01, eid, ne, nxs, xstart, xord, xp, nb, pofs, nps, pord, pab, ptg, &
                   nm, me, with_iv, with_bnd, s, nxc, rw(ro(1)), rw(ro(2)), rw(ro(3)), rw(ro(4)), rw(ro(5)), rw(ro(6)), &
                   rw(ro(7)), rw(ro(8)), rw(ro(9)), rw(ro(10)), rw(ro(11)), rw(ro(12)), rw(ro(13)), rw(ro(14)), &
                   rw(ro(15)), iw(1), iw(io(1)), iw(io(2)), iw(io(3)), iw(io(4)), iw(io(5)), iw(io(6)), iw(io(7)), &
                   iw(io(8)), iw(io(9)), iw(io(10)), iw(io(11)), iw(io(12)), iw(io(13)), iw(io(14)), iw(io(15)), &
                   iw(io(16)), iw(io(17)), iw(io(18)), closed, alpha, aperture, centroid, wet_length, wet_centroid, &
                   wet_normal, npieces, nivc, nbdc)
   endsubroutine cut_cell_2d_count

   subroutine cut_cell_2d_fill(x0, x1, y0, y1, i, j, eid, nm, me, with_iv, with_bnd, nrw, rw, niw, iw, p0, iv0, bd0, &
                               np, pcell, parea, pcen, pwl, pwn, pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, &
                               bwall)
   !< Kernel, fill pass (issue #100): the pieces of cut cell (i, j) kept in its workspace by `cut_cell_2d_count`, written
   !< as pieces p0 + 1.., fluid edge intervals iv0 + 1.. and boundary segments bd0 + 1.. of the block.
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1       !< Cell.
   integer(I4P), intent(in)    :: i, j                 !< Cell indexes.
   integer(I4P), intent(in)    :: eid(4)               !< Edge ids (bottom, right, top, left).
   integer(I4P), intent(in)    :: nm                   !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)            !< Moment exponents.
   logical,      intent(in)    :: with_iv, with_bnd    !< Write the intervals, the boundary segments.
   integer(I4P), intent(in)    :: nrw, niw             !< Workspace sizes.
   real(R8P),    intent(in)    :: rw(nrw)              !< Real workspace.
   integer(I4P), intent(in)    :: iw(niw)              !< Integer workspace.
   integer(I4P), intent(in)    :: p0, iv0, bd0         !< Pieces, intervals, boundary segments before the cell's.
   integer(I4P), intent(in)    :: np, niv, nbd         !< Pieces, intervals, boundary segments of the block.
   integer(I4P), intent(inout) :: pcell(2, np)         !< Cell of each piece.
   real(R8P),    intent(inout) :: parea(np), pcen(2, np), pwl(np), pwn(2, np), pap(4, np) !< Piece measures.
   real(R8P),    intent(inout) :: pmom(nm, np), pwmom(nm, np) !< Piece moments, wall moments.
   integer(I4P), intent(inout) :: ivs(np + 1)          !< First interval of each piece.
   integer(I4P), intent(inout) :: ive(niv)             !< Interval edges.
   real(R8P),    intent(inout) :: ivl(niv), ivh(niv)   !< Interval ranges.
   integer(I4P), intent(inout) :: bds(np + 1)          !< First boundary segment of each piece.
   real(R8P),    intent(inout) :: bseg(4, nbd)         !< Boundary segments.
   logical,      intent(inout) :: bwall(nbd)           !< Wall flags.
   integer(I4P)                :: nxc, nb, s, ro(15), io(18) !< Work.

   if (iw(1) == 0_I4P) return
   if (iw(1) == 2_I4P) then
      call quiet_piece_2d(x0, x1, y0, y1, i, j, eid, nm, me, with_iv, with_bnd, p0, iv0, bd0, np, pcell, parea, pcen, &
                          pwl, pwn, pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
      return
   endif
   nxc = iw(6) ; nb = iw(7)
   call cut_cell_2d_layout(nxc, nb, nm, ro, io)
   s = 4_I4P + nxc + nb
   call fill_core(x0, x1, y0, y1, i, j, eid, nb, nm, with_iv, with_bnd, s, rw(ro(1)), rw(ro(2)), rw(ro(3)), rw(ro(4)), &
                  rw(ro(5)), rw(ro(6)), rw(ro(13)), iw(1), iw(io(1)), iw(io(3)), iw(io(4)), iw(io(7)), iw(io(14)), &
                  iw(io(15)), iw(io(16)), iw(io(17)), iw(io(18)), p0, iv0, bd0, np, pcell, parea, pcen, pwl, pwn, pap, &
                  pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
   endsubroutine cut_cell_2d_fill

   ! nvfortran 26.9 miscompiles this routine's device code at -O2 (issue #100 step 2: the loop `do k = 1, ns` ran once,
   ! the corner constructors passed to edge_subs lost their second element); opt=1 for this routine only, host and
   ! device alike (the two stay bit-identical). A comment to every other compiler.
!pgi$r opt=1
   subroutine count_core(x0, x1, y0, y1, f00, f10, f11, f01, eid, ne, nxs, xstart, xord, xp, nb, pofs, nps, pord, pab, &
                         ptg, nm, me, with_iv, with_bnd, s, nxm, sa, sb, acc, wn, accm, wallm, harea, ig, iwl, iwn, igm, &
                         iwm, lp, cgm, cwm, hdr, sfl, sx, sed, par, xsub, pvis, loopp, hroot, hfirst, hn, ifirst, inn, it, &
                         tkind, tstart, tn, tslot, lt, closed, alpha, aperture, centroid, wet_length, wet_centroid, &
                         wet_normal, npieces, nivc, nbdc)
   !< The count pass on the workspace parts (see `cut_cell_2d_layout`). Sub-intervals of the cell boundary
   !< counter-clockwise (1 bottom, 2 right, 3 top, 4 left), union-find over the fluid ones, the wall chains between
   !< crossings, the closed loops inside the cell (holes: pieces of their own; islands: into their host), Green sums.
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1         !< Cell.
   logical,      intent(in)    :: f00, f10, f11, f01     !< Corner flags.
   integer(I4P), intent(in)    :: eid(4)                 !< Edge ids.
   integer(I4P), intent(in)    :: ne, nxs                !< Edges, crossings.
   integer(I4P), intent(in)    :: xstart(ne + 1)         !< Crossings by edge.
   integer(I4P), intent(in)    :: xord(nxs)              !< Crossings ordered along each edge.
   real(R8P),    intent(in)    :: xp(2, nxs)             !< Crossing points.
   integer(I4P), intent(in)    :: nb, pofs, nps          !< Wall pieces of the cell, offset, wall pieces.
   integer(I4P), intent(in)    :: pord(nps)              !< Wall pieces by cell.
   real(R8P),    intent(in)    :: pab(4, nps)            !< Wall piece end points.
   integer(I4P), intent(in)    :: ptg(2, nps)            !< Wall piece tags.
   integer(I4P), intent(in)    :: nm                     !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)              !< Exponents.
   logical,      intent(in)    :: with_iv, with_bnd      !< Count intervals, boundary segments.
   integer(I4P), intent(in)    :: s, nxm                 !< Sub-interval slots, crossing slots.
   real(R8P),    intent(inout) :: sa(2, s), sb(2, s)     !< Sub-interval end points (counter-clockwise).
   real(R8P),    intent(inout) :: acc(6, s)              !< Per root: area, mx, my, wet length, wet mx, wet my.
   real(R8P),    intent(inout) :: wn(2, s)               !< Per root: wet normal.
   real(R8P),    intent(inout) :: accm(nm, s), wallm(nm, s) !< Per root: moments, wall moments.
   real(R8P),    intent(inout) :: harea(nb)              !< Hole areas.
   real(R8P),    intent(inout) :: ig(3, nb), iwl(3, nb), iwn(2, nb) !< Island Green, wet, normal.
   real(R8P),    intent(inout) :: igm(nm, nb), iwm(nm, nb) !< Island moments, wall moments.
   real(R8P),    intent(inout) :: lp(4, nb)              !< The cell's wall pieces (a, b).
   real(R8P),    intent(inout) :: cgm(nm), cwm(nm)       !< Chain moments, wall moments.
   integer(I4P), intent(inout) :: hdr(NHDR)              !< Header.
   integer(I4P), intent(inout) :: sfl(s), sx(s), sed(s), par(s) !< Fluid (1), end crossing (0 corner), edge, parent.
   integer(I4P), intent(inout) :: xsub(2, nxm)           !< Crossing slots: (crossing id, sub-interval before it).
   integer(I4P), intent(inout) :: pvis(nb), loopp(nb)    !< Visited (1), loop pieces.
   integer(I4P), intent(inout) :: hroot(nb), hfirst(nb), hn(nb) !< Hole slot, loop range.
   integer(I4P), intent(inout) :: ifirst(nb), inn(nb), it(nb) !< Island loop range, traversal.
   integer(I4P), intent(inout) :: tkind(nb), tstart(nb), tn(nb), tslot(nb) !< Traversals: 1 chain from wall piece tstart,
                                                         !< 2 loop of tn pieces from loopp(tstart); charged slot.
   integer(I4P), intent(inout) :: lt(2, nb)              !< The cell's wall piece tags (ta, tb).
   logical,      intent(out)   :: closed                 !< Closed.
   real(R8P),    intent(out)   :: alpha, aperture(4), centroid(2), wet_length, wet_centroid(2), wet_normal(2) !< Measures.
   integer(I4P), intent(out)   :: npieces, nivc, nbdc    !< Counts.
   real(R8P)    :: area_c, len(4), flen(4), gs(3), sm(2), hm(2), lg(3), lw(4), lwn(2), sw5, sw6 !< Work.
   integer(I4P) :: q, k, m, r, h, rr, ns, nsb, nxc, nh, nisl, nloopp, nt, ncs, s0, p1, pid, nfl !< Work.
   logical      :: state, ok !< Work.

   area_c = (x1 - x0) * (y1 - y0)
   len = [y1 - y0, y1 - y0, x1 - x0, x1 - x0]
   sm = 0.5_R8P * [x0 + x1, y0 + y1] ; hm = [x1 - x0, y1 - y0]
   closed = .true. ; npieces = 0_I4P ; nivc = 0_I4P ; nbdc = 0_I4P
   wet_length = 0._R8P ; wet_centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)] ; wet_normal = 0._R8P
   hdr = 0_I4P ; hdr(6) = nxm ; hdr(7) = nb
   do q = 1, nb
      pid = pord(pofs + q) ; lp(:, q) = pab(:, pid) ; lt(:, q) = ptg(:, pid)
   enddo
   ! sub-intervals of the edges counter-clockwise: bottom (+u), right (+v), top (-u), left (-v)
   ns = 0 ; nxc = 0 ; ok = .true. ; flen = 0._R8P
   state = f00
   call edge_subs(1_I4P, [x0, y0], [x1, y0], xstart(eid(1)), xstart(eid(1) + 1) - 1_I4P, nxs, xord, xp, .true., f10, &
                  state, ok, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   call edge_subs(2_I4P, [x1, y0], [x1, y1], xstart(eid(2)), xstart(eid(2) + 1) - 1_I4P, nxs, xord, xp, .true., f11, &
                  state, ok, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   call edge_subs(3_I4P, [x1, y1], [x0, y1], xstart(eid(3)), xstart(eid(3) + 1) - 1_I4P, nxs, xord, xp, .true., f01, &
                  state, ok, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   call edge_subs(4_I4P, [x0, y1], [x0, y0], xstart(eid(4)), xstart(eid(4) + 1) - 1_I4P, nxs, xord, xp, .false., f00, &
                  state, ok, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   if (state .neqv. f00) ok = .false.
   nsb = ns
   if (.not. ok) then
      call open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal) ; return
   endif
   if (nb == 0 .and. nxc == 0) then
      if (f00) then
         alpha = 0._R8P ; aperture = 0._R8P
         centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
      else
         alpha = 1._R8P ; aperture = 1._R8P
         centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
         npieces = 1_I4P ; hdr(1) = 2_I4P
         if (with_iv) nivc = 4_I4P
         if (with_bnd) nbdc = 4_I4P
      endif
      return
   endif
   ! union-find over the fluid sub-intervals; consecutive ones through corners join. The accumulators are cleared first
   ! (the unions move zeros; before #100 they moved stale values that the clear after them overwrote: same state)
   do k = 1, ns
      par(k) = k ; acc(:, k) = 0._R8P ; wn(:, k) = 0._R8P ; accm(:, k) = 0._R8P ; wallm(:, k) = 0._R8P
   enddo
   do k = 1, nsb
      m = mod(k, nsb) + 1
      if (sx(k) == 0 .and. sfl(k) == 1 .and. sfl(m) == 1) call unite(k, m, s, par, nm, acc, wn, accm, wallm)
   enddo
   ! wall chains: from a piece starting at a crossing, follow the vertices to a crossing
   pvis = 0_I4P ; ncs = 0_I4P ; nt = 0_I4P
   do q = 1, nb
      if (lt(1, q) <= 0) cycle
      call follow(q, nb, lp, lt, pvis, p1, lg, lw, lwn, ok, nm, me, sm, hm, area_c, cgm, cwm, ncs)
      if (.not. ok) then
         call open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal) ; return
      endif
      ! join the fluid sub-intervals next to both end crossings, charge the chain to that root
      r = fluid_sub_at(lt(1, q), nxc, nxm, xsub, nsb, s, sfl) ; h = fluid_sub_at(lt(2, p1), nxc, nxm, xsub, nsb, s, sfl)
      if (r == 0 .or. h == 0) then
         call open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal) ; return
      endif
      call unite(r, h, s, par, nm, acc, wn, accm, wallm)
      rr = find(r, s, par)
      call charge(rr, lg, lw, lwn, cgm, cwm, s, nm, acc, wn, accm, wallm)
      nt = nt + 1_I4P ; tkind(nt) = 1_I4P ; tstart(nt) = q ; tn(nt) = 0_I4P ; tslot(nt) = rr
   enddo
   ! remaining pieces: closed loops entirely inside the cell
   nh = 0 ; nisl = 0 ; nloopp = 0
   do q = 1, nb
      if (pvis(q) == 1) cycle
      s0 = nloopp + 1
      call follow_loop(q, nb, lp, lt, pvis, nloopp, loopp, lg, lw, lwn, ok, nm, me, sm, hm, area_c, cgm, cwm, ncs)
      if (.not. ok) then
         call open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal) ; return
      endif
      nt = nt + 1_I4P ; tkind(nt) = 2_I4P ; tstart(nt) = s0 ; tn(nt) = nloopp - s0 + 1
      if (lg(1) < 0._R8P) then ! clockwise: a hole, a piece of its own (a slot <= s: one loop per wall piece at most)
         nh = nh + 1 ; harea(nh) = -lg(1) ; hfirst(nh) = s0 ; hn(nh) = nloopp - s0 + 1
         ns = ns + 1 ; par(ns) = ns ; sfl(ns) = 0_I4P ; hroot(nh) = ns
         acc(:, ns) = 0._R8P ; wn(:, ns) = 0._R8P ; accm(:, ns) = 0._R8P ; wallm(:, ns) = 0._R8P
         call charge(ns, lg, lw, lwn, cgm, cwm, s, nm, acc, wn, accm, wallm)
         tslot(nt) = ns
      else
         nisl = nisl + 1 ; ifirst(nisl) = s0 ; inn(nisl) = nloopp - s0 + 1 ; it(nisl) = nt
         ig(:, nisl) = lg ; iwl(:, nisl) = lw(1:3) ; iwn(:, nisl) = lwn
         igm(:, nisl) = cgm ; iwm(:, nisl) = cwm ; tslot(nt) = 0_I4P
      endif
   enddo
   ! islands: into the innermost containing hole, else the containing boundary piece
   do k = 1, nisl
      call island_host(k, nb, lp, lt, loopp, ifirst, nh, hfirst, hn, harea, hroot, nsb, s, sfl, par, sa, sb, nxc, nxm, &
                       xsub, r)
      if (r == 0) then
         call open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal) ; return
      endif
      rr = find(r, s, par)
      call charge(rr, ig(:, k), iwl(:, k), iwn(:, k), igm(:, k), iwm(:, k), s, nm, acc, wn, accm, wallm)
      tslot(it(k)) = rr
   enddo
   ! fluid sub-intervals: Green along them (counter-clockwise)
   nfl = 0_I4P
   do k = 1, nsb
      if (sfl(k) == 0) cycle
      nfl = nfl + 1_I4P
      call green(sa(:, k), sb(:, k), gs)
      r = find(k, s, par)
      acc(1:3, r) = acc(1:3, r) + gs
      if (nm > 0_I4P) call segment_green_moments((sa(:, k) - sm) / hm, (sb(:, k) - sm) / hm, me, area_c, accm(:, r))
   enddo
   ! the pieces (roots): cell sums in root order; every fluid sub-interval and every chain or loop piece lies in one
   alpha = 0._R8P ; centroid = 0._R8P
   do k = 1, ns
      if (k <= nsb) then
         if (sfl(k) == 0) cycle
      endif
      if (find(k, s, par) /= k) cycle
      npieces = npieces + 1_I4P
      alpha = alpha + acc(1, k)
      centroid = centroid + acc(2:3, k)
      wet_length = wet_length + acc(4, k)
      wet_normal = wet_normal + wn(:, k)
   enddo
   if (with_iv) nivc = nfl
   if (with_bnd) nbdc = nfl + ncs
   if (alpha > 0._R8P) then
      centroid = centroid / alpha
   else
      centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
   endif
   ! wetted centroid from the per-root moments (summed in root order from zero, as the masked sum it replaces)
   if (wet_length > 0._R8P) then
      sw5 = 0._R8P ; sw6 = 0._R8P
      do k = 1, ns
         if (find(k, s, par) == k .and. (k > nsb .or. sfl(k) == 1)) sw5 = sw5 + acc(5, k)
      enddo
      do k = 1, ns
         if (find(k, s, par) == k .and. (k > nsb .or. sfl(k) == 1)) sw6 = sw6 + acc(6, k)
      enddo
      wet_centroid = [sw5, sw6] / wet_length
   endif
   alpha = alpha / area_c
   aperture = [flen(4), flen(2), flen(1), flen(3)] / len ! counter-clockwise edges -> (-u, +u, -v, +v)
   hdr(1) = 1_I4P ; hdr(2) = ns ; hdr(3) = nh ; hdr(4) = nsb ; hdr(5) = nt
   endsubroutine count_core

   subroutine fill_core(x0, x1, y0, y1, i, j, eid, nb, nm, with_iv, with_bnd, s, sa, sb, acc, wn, accm, wallm, lp, hdr, &
                        sfl, sed, par, loopp, tkind, tstart, tn, tslot, lt, p0, iv0, bd0, np, pcell, parea, pcen, pwl, pwn, &
                        pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
   !< The fill pass on the workspace parts left by `count_core`: per piece (root, in order) its measures, its fluid edge
   !< intervals in the order (-u, +u, -v, +v) and its oriented boundary (fluid sub-intervals, then the chains and loops
   !< charged to it, reversed so that the fluid is on the left, in the order they were followed).
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1       !< Cell.
   integer(I4P), intent(in)    :: i, j, eid(4)         !< Cell indexes, edge ids.
   integer(I4P), intent(in)    :: nb, nm               !< Wall pieces of the cell, moments.
   logical,      intent(in)    :: with_iv, with_bnd    !< Write intervals, boundary segments.
   integer(I4P), intent(in)    :: s                    !< Sub-interval slots.
   real(R8P),    intent(in)    :: sa(2, s), sb(2, s), acc(6, s), wn(2, s), accm(nm, s), wallm(nm, s) !< Sub-intervals, roots.
   real(R8P),    intent(in)    :: lp(4, nb)            !< Wall pieces.
   integer(I4P), intent(in)    :: hdr(NHDR)            !< Header.
   integer(I4P), intent(in)    :: sfl(s), sed(s), par(s) !< Fluid, edge, parent.
   integer(I4P), intent(in)    :: loopp(nb), tkind(nb), tstart(nb), tn(nb), tslot(nb) !< Loops, traversals.
   integer(I4P), intent(in)    :: lt(2, nb)            !< Wall piece tags.
   integer(I4P), intent(in)    :: p0, iv0, bd0         !< Offsets.
   integer(I4P), intent(in)    :: np, niv, nbd         !< Sizes.
   integer(I4P), intent(inout) :: pcell(2, np)         !< Cell of each piece.
   real(R8P),    intent(inout) :: parea(np), pcen(2, np), pwl(np), pwn(2, np), pap(4, np) !< Piece measures.
   real(R8P),    intent(inout) :: pmom(nm, np), pwmom(nm, np) !< Piece moments.
   integer(I4P), intent(inout) :: ivs(np + 1), ive(niv) !< Intervals.
   real(R8P),    intent(inout) :: ivl(niv), ivh(niv)   !< Interval ranges.
   integer(I4P), intent(inout) :: bds(np + 1)          !< Boundary segments of the pieces.
   real(R8P),    intent(inout) :: bseg(4, nbd)         !< Boundary segments.
   logical,      intent(inout) :: bwall(nbd)           !< Wall flags.
   real(R8P)    :: len(4)                              !< Edge lengths (-u, +u, -v, +v).
   integer(I4P) :: ns, nsb, nt, k, e, m, q, t, p, iv, bd, cur, guard !< Work.

   len = [y1 - y0, y1 - y0, x1 - x0, x1 - x0]
   ns = hdr(2) ; nsb = hdr(4) ; nt = hdr(5)
   p = p0 ; iv = iv0 ; bd = bd0
   do k = 1, ns
      if (k <= nsb) then
         if (sfl(k) == 0) cycle
      endif
      if (find(k, s, par) /= k) cycle
      p = p + 1_I4P
      pcell(1, p) = i ; pcell(2, p) = j
      parea(p) = acc(1, k) ; pwl(p) = acc(4, k) ; pwn(1, p) = wn(1, k) ; pwn(2, p) = wn(2, k)
      if (acc(1, k) > 0._R8P) then
         pcen(1, p) = acc(2, k) / acc(1, k) ; pcen(2, p) = acc(3, k) / acc(1, k)
      else
         pcen(1, p) = 0._R8P ; pcen(2, p) = 0._R8P
      endif
      if (nm > 0_I4P) then
         pmom(:, p) = accm(:, k) ; pwmom(:, p) = wallm(:, k)
      endif
      ! its fluid edge intervals (issue #70), in the order (-u, +u, -v, +v) = edges (4, 2, 1, 3)
      pap(:, p) = 0._R8P ; ivs(p) = iv + 1_I4P
      do e = 1, 4
         m = edge_of_side(e)
         do q = 1, nsb
            if (sfl(q) == 0 .or. sed(q) /= m) cycle
            if (find(q, s, par) /= k) cycle
            pap(e, p) = pap(e, p) + norm2_xy(sb(1, q) - sa(1, q), sb(2, q) - sa(2, q)) / len(e)
            if (.not. with_iv) cycle
            iv = iv + 1_I4P ; ive(iv) = eid(m)
            if (m == 2 .or. m == 4) then
               ivl(iv) = min(sa(2, q), sb(2, q)) ; ivh(iv) = max(sa(2, q), sb(2, q))
            else
               ivl(iv) = min(sa(1, q), sb(1, q)) ; ivh(iv) = max(sa(1, q), sb(1, q))
            endif
         enddo
      enddo
      bds(p) = bd + 1_I4P
      if (.not. with_bnd) cycle
      ! its oriented boundary: fluid sub-intervals, then the chains and loops charged to it (reversed)
      do q = 1, nsb
         if (sfl(q) == 0) cycle
         if (find(q, s, par) /= k) cycle
         bd = bd + 1_I4P ; bseg(1:2, bd) = sa(:, q) ; bseg(3:4, bd) = sb(:, q) ; bwall(bd) = .false.
      enddo
      do t = 1, nt
         if (find(tslot(t), s, par) /= k) cycle
         if (tkind(t) == 1_I4P) then ! the chain as `follow` walked it
            cur = tstart(t) ; guard = 0_I4P
            do
               bd = bd + 1_I4P ; bseg(1:2, bd) = lp(3:4, cur) ; bseg(3:4, bd) = lp(1:2, cur) ; bwall(bd) = .true.
               if (lt(2, cur) > 0) exit
               cur = local_start(-lt(2, cur), nb, lt)
               guard = guard + 1_I4P
               if (cur == 0 .or. guard > nb) exit ! cannot happen: the count pass closed this chain
            enddo
         else
            do q = tstart(t), tstart(t) + tn(t) - 1_I4P
               cur = loopp(q)
               bd = bd + 1_I4P ; bseg(1:2, bd) = lp(3:4, cur) ; bseg(3:4, bd) = lp(1:2, cur) ; bwall(bd) = .true.
            enddo
         endif
      enddo
   enddo
   endsubroutine fill_core

   subroutine quiet_cell_2d_kernel(x0, x1, y0, y1, i, j, solid, eid, nm, me, with_iv, with_bnd, closed, alpha, aperture, &
                                   centroid, wet_length, wet_centroid, wet_normal, npieces, p0, iv0, bd0, np, pcell, parea, &
                                   pcen, pwl, pwn, pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
   !< Kernel of a quiet cell (no wall pieces or crossings): all solid, or all fluid as one piece written at p0 + 1.
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1       !< Cell.
   integer(I4P), intent(in)    :: i, j                 !< Cell indexes.
   logical,      intent(in)    :: solid                !< The cell is inside the solid.
   integer(I4P), intent(in)    :: eid(4)               !< Edge ids (bottom, right, top, left).
   integer(I4P), intent(in)    :: nm                   !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)            !< Moment exponents.
   logical,      intent(in)    :: with_iv, with_bnd    !< Write intervals, boundary segments.
   logical,      intent(out)   :: closed               !< Closed.
   real(R8P),    intent(out)   :: alpha, aperture(4), centroid(2), wet_length, wet_centroid(2), wet_normal(2) !< Measures.
   integer(I4P), intent(out)   :: npieces              !< Pieces.
   integer(I4P), intent(in)    :: p0, iv0, bd0         !< Offsets.
   integer(I4P), intent(in)    :: np, niv, nbd         !< Sizes.
   integer(I4P), intent(inout) :: pcell(2, np)         !< Cell of each piece.
   real(R8P),    intent(inout) :: parea(np), pcen(2, np), pwl(np), pwn(2, np), pap(4, np) !< Piece measures.
   real(R8P),    intent(inout) :: pmom(nm, np), pwmom(nm, np) !< Piece moments.
   integer(I4P), intent(inout) :: ivs(np + 1), ive(niv) !< Intervals.
   real(R8P),    intent(inout) :: ivl(niv), ivh(niv)   !< Interval ranges.
   integer(I4P), intent(inout) :: bds(np + 1)          !< Boundary segments of the pieces.
   real(R8P),    intent(inout) :: bseg(4, nbd)         !< Boundary segments.
   logical,      intent(inout) :: bwall(nbd)           !< Wall flags.

   call quiet_cell_2d_measures(x0, x1, y0, y1, solid, closed, alpha, aperture, centroid, wet_length, wet_centroid, &
                               wet_normal, npieces)
   if (solid) return
   call quiet_piece_2d(x0, x1, y0, y1, i, j, eid, nm, me, with_iv, with_bnd, p0, iv0, bd0, np, pcell, parea, pcen, pwl, &
                       pwn, pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
   endsubroutine quiet_cell_2d_kernel

   pure subroutine quiet_cell_2d_measures(x0, x1, y0, y1, solid, closed, alpha, aperture, centroid, wet_length, &
                                          wet_centroid, wet_normal, npieces)
   !< Cell measures of a quiet cell (no wall pieces or crossings): all solid, or all fluid (one piece).
   !$acc routine seq
   real(R8P),    intent(in)  :: x0, x1, y0, y1       !< Cell.
   logical,      intent(in)  :: solid                !< The cell is inside the solid.
   logical,      intent(out) :: closed               !< Closed.
   real(R8P),    intent(out) :: alpha, aperture(4), centroid(2), wet_length, wet_centroid(2), wet_normal(2) !< Measures.
   integer(I4P), intent(out) :: npieces              !< Pieces.

   closed = .true. ; npieces = 0_I4P
   wet_length = 0._R8P ; wet_centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
   wet_normal = 0._R8P
   centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)]
   if (solid) then
      alpha = 0._R8P ; aperture = 0._R8P
      return
   endif
   alpha = 1._R8P ; aperture = 1._R8P ; npieces = 1_I4P
   endsubroutine quiet_cell_2d_measures

   subroutine quiet_piece_2d(x0, x1, y0, y1, i, j, eid, nm, me, with_iv, with_bnd, p0, iv0, bd0, np, pcell, parea, pcen, &
                             pwl, pwn, pap, pmom, pwmom, ivs, niv, ive, ivl, ivh, bds, nbd, bseg, bwall)
   !< The piece of a fluid cell without walls, written at p0 + 1 (intervals at iv0 + 1.., boundary at bd0 + 1..).
   !$acc routine seq
   real(R8P),    intent(in)    :: x0, x1, y0, y1       !< Cell.
   integer(I4P), intent(in)    :: i, j, eid(4)         !< Cell indexes, edge ids.
   integer(I4P), intent(in)    :: nm                   !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)            !< Moment exponents.
   logical,      intent(in)    :: with_iv, with_bnd    !< Write intervals, boundary segments.
   integer(I4P), intent(in)    :: p0, iv0, bd0         !< Offsets.
   integer(I4P), intent(in)    :: np, niv, nbd         !< Sizes.
   integer(I4P), intent(inout) :: pcell(2, np)         !< Cell of each piece.
   real(R8P),    intent(inout) :: parea(np), pcen(2, np), pwl(np), pwn(2, np), pap(4, np) !< Piece measures.
   real(R8P),    intent(inout) :: pmom(nm, np), pwmom(nm, np) !< Piece moments.
   integer(I4P), intent(inout) :: ivs(np + 1), ive(niv) !< Intervals.
   real(R8P),    intent(inout) :: ivl(niv), ivh(niv)   !< Interval ranges.
   integer(I4P), intent(inout) :: bds(np + 1)          !< Boundary segments of the pieces.
   real(R8P),    intent(inout) :: bseg(4, nbd)         !< Boundary segments.
   logical,      intent(inout) :: bwall(nbd)           !< Wall flags.
   real(R8P)    :: area_c, mx, my                      !< Cell area, first moments.
   integer(I4P) :: p, k                                !< Work.

   area_c = (x1 - x0) * (y1 - y0)
   mx = 0.5_R8P * (x0 + x1) * area_c ; my = 0.5_R8P * (y0 + y1) * area_c
   p = p0 + 1_I4P
   pcell(1, p) = i ; pcell(2, p) = j
   parea(p) = area_c ; pwl(p) = 0._R8P ; pwn(1, p) = 0._R8P ; pwn(2, p) = 0._R8P
   if (area_c > 0._R8P) then
      pcen(1, p) = mx / area_c ; pcen(2, p) = my / area_c
   else
      pcen(1, p) = 0._R8P ; pcen(2, p) = 0._R8P
   endif
   pap(:, p) = 1._R8P
   if (nm > 0_I4P) then
      call cell_moments_into(2_I4P, nm, me, area_c, pmom(1, p)) ; pwmom(:, p) = 0._R8P
   endif
   ivs(p) = iv0 + 1_I4P
   if (with_iv) then ! (-u, +u, -v, +v)
      k = iv0
      ive(k + 1) = eid(4) ; ivl(k + 1) = y0 ; ivh(k + 1) = y1
      ive(k + 2) = eid(2) ; ivl(k + 2) = y0 ; ivh(k + 2) = y1
      ive(k + 3) = eid(1) ; ivl(k + 3) = x0 ; ivh(k + 3) = x1
      ive(k + 4) = eid(3) ; ivl(k + 4) = x0 ; ivh(k + 4) = x1
   endif
   bds(p) = bd0 + 1_I4P
   if (with_bnd) then ! counter-clockwise: bottom, right, top, left
      k = bd0
      bseg(:, k + 1) = [x0, y0, x1, y0] ; bseg(:, k + 2) = [x1, y0, x1, y1]
      bseg(:, k + 3) = [x1, y1, x0, y1] ; bseg(:, k + 4) = [x0, y1, x0, y0]
      bwall(k + 1:k + 4) = .false.
   endif
   endsubroutine quiet_piece_2d

   subroutine cut_cell_2d_pieces(u, v, next, x0, x1, y0, y1, flags, closed, np, parea, nw, wseg, wpiece, cpiece, plen, pmom)
   !< The pieces of the single cell [x0, x1] x [y0, y1] of the curve (u, v, next) with given corner inside flags, and the
   !< boundary features of each piece (issue #102, the piece map of a time step): the wall pieces (one per segment
   !< meeting the cell: a straight segment clipped to a convex cell bounds exactly one piece), the fluid corners, and the
   !< fluid length and first moment of each of the four faces per piece.
   !<
   !< The curve may be any part of a curve that holds every segment meeting the cell (the corner flags replace the
   !< parity of the far segments); the segments must be listed in increasing order of the whole curve's indexes. The
   !< cell is walked as a one-cell block (`walk_all`, the crossing formulas of `fossil_line_query`) and counted by the
   !< kernel of `compute_cut_cells_2d` (`cut_cell_2d_count`), so with the flags of `node_flags` the pieces, their order
   !< and their areas are those of the block's cell, bit for bit. Then `cut_cell_2d_labels` reads the features off the
   !< count's workspace. `closed` false (an open chain, or flags that disagree with the crossings): no pieces. Since
   !< issue #105 the device kernels (`cut_cell_2d_walk`, `cut_cell_2d_one_pieces`) on a host workspace, bit for bit.
   real(R8P),                 intent(in)  :: u(:), v(:)     !< Vertices.
   integer(I4P),              intent(in)  :: next(:)        !< Next vertex (0: none).
   real(R8P),                 intent(in)  :: x0, x1, y0, y1 !< Cell.
   logical,                   intent(in)  :: flags(4)       !< Corners inside the solid: (x0,y0), (x1,y0), (x1,y1), (x0,y1).
   logical,                   intent(out) :: closed         !< The cell boundary closed.
   integer(I4P),              intent(out) :: np             !< Pieces.
   real(R8P),    allocatable, intent(out) :: parea(:)       !< Piece areas [np].
   integer(I4P),              intent(out) :: nw             !< Wall pieces.
   integer(I4P), allocatable, intent(out) :: wseg(:)        !< Segment (its start vertex in u, v) of each wall piece [nw].
   integer(I4P), allocatable, intent(out) :: wpiece(:)      !< Piece on the fluid side of each wall piece [nw].
   integer(I4P),              intent(out) :: cpiece(4)      !< Piece of each corner (0: solid), corners as `flags`.
   real(R8P),    allocatable, intent(out) :: plen(:,:)      !< Fluid length of the faces (-u, +u, -v, +v) per piece [4, np].
   real(R8P),    allocatable, intent(out) :: pmom(:,:)      !< int (coordinate along the face) over it [4, np].
   real(R8P),    allocatable :: rw(:), pa(:), pl(:,:), pm(:,:) !< Workspace, piece outputs.
   integer(I4P), allocatable :: iw(:), ws(:), wp(:)          !< Workspace, wall piece outputs.
   integer(I4P)              :: nu, nxs, nps, nrw, niw, mp, idum(2, 1) !< Work.
   real(R8P)                 :: rdum(4, 1)                   !< Work.

   ! the kernels of the device path (issue #105): the walk's count pass, then the evaluation on a workspace
   nu = size(u, kind=I4P)
   call cut_cell_2d_walk(nu, u, v, next, x0, x1, y0, y1, .false., 0_I4P, 0_I4P, nxs, nps, idum, rdum, idum, rdum, idum, &
                         rdum, idum)
   call cut_cell_2d_one_sizes(nxs, nps, nrw, niw)
   mp = 4_I4P + nxs + nps
   allocate(rw(nrw), iw(niw), pa(mp), pl(4, mp), pm(4, mp), ws(max(1_I4P, nps)), wp(max(1_I4P, nps)))
   call cut_cell_2d_one_pieces(nu, u, v, next, x0, x1, y0, y1, flags(1), flags(2), flags(3), flags(4), nxs, nps, nrw, rw, &
                               niw, iw, closed, np, mp, pa, ws, wp, cpiece, pl, pm)
   nw = nps
   allocate(wseg(nw), wpiece(nw), parea(np), plen(4, np), pmom(4, np))
   wseg = ws(1:nw) ; wpiece = wp(1:nw)
   if (.not. closed .or. np == 0_I4P) return
   parea = pa(1:np) ; plen = pl(:, 1:np) ; pmom = pm(:, 1:np)
   endsubroutine cut_cell_2d_pieces

   subroutine cut_cell_2d_labels(nb, s, np, hdr, sa, sb, acc, sfl, sx, sed, par, loopp, tkind, tstart, tn, tslot, lt, rix, &
                                 parea, wpiece, cpiece, plen, pmom)
   !< Kernel (issue #102): the boundary features of the pieces of a cut cell left by `count_core` in its workspace
   !< (status 1). Pieces in the root order of `fill_core` (the piece numbering of the cut cells): their areas, the piece
   !< of every wall piece (through the chains and loops charged to it), the piece of every fluid corner, and per piece the
   !< fluid length and first moment of each face (-u, +u, -v, +v) from its fluid sub-intervals.
   !$acc routine seq
   integer(I4P), intent(in)    :: nb, s, np           !< Wall pieces, sub-interval slots, pieces.
   integer(I4P), intent(in)    :: hdr(NHDR)           !< Header.
   real(R8P),    intent(in)    :: sa(2, s), sb(2, s)  !< Sub-intervals.
   real(R8P),    intent(in)    :: acc(6, s)           !< Root accumulators.
   integer(I4P), intent(in)    :: sfl(s), sx(s), sed(s), par(s) !< Fluid, end crossing, edge, parent.
   integer(I4P), intent(in)    :: loopp(nb), tkind(nb), tstart(nb), tn(nb), tslot(nb) !< Loops, traversals.
   integer(I4P), intent(in)    :: lt(2, nb)           !< Wall piece tags.
   integer(I4P), intent(inout) :: rix(s)              !< Work: piece of each root.
   real(R8P),    intent(out)   :: parea(np)           !< Piece areas.
   integer(I4P), intent(out)   :: wpiece(nb)          !< Piece of each wall piece.
   integer(I4P), intent(out)   :: cpiece(4)           !< Piece of each corner (0: solid).
   real(R8P),    intent(out)   :: plen(4, np), pmom(4, np) !< Face lengths and first moments per piece.
   integer(I4P) :: ns, nsb, nt, k, p, t, m, cur, guard, f, cn !< Work.
   real(R8P)    :: a, b                                !< Interval along the face.

   ns = hdr(2) ; nsb = hdr(4) ; nt = hdr(5)
   p = 0_I4P
   do k = 1_I4P, ns
      rix(k) = 0_I4P
      if (k <= nsb) then
         if (sfl(k) == 0_I4P) cycle
      endif
      if (find(k, s, par) /= k) cycle
      p = p + 1_I4P ; rix(k) = p
      parea(p) = acc(1, k)
      do f = 1_I4P, 4_I4P
         plen(f, p) = 0._R8P ; pmom(f, p) = 0._R8P
      enddo
   enddo
   cpiece = 0_I4P
   do k = 1_I4P, nsb
      if (sfl(k) == 0_I4P) cycle
      p = rix(find(k, s, par))
      select case (sed(k)) ! counter-clockwise edges: 1 bottom (-v), 2 right (+u), 3 top (+v), 4 left (-u)
      case (1_I4P)
         f = 3_I4P ; cn = 2_I4P ; a = min(sa(1, k), sb(1, k)) ; b = max(sa(1, k), sb(1, k))
      case (2_I4P)
         f = 2_I4P ; cn = 3_I4P ; a = min(sa(2, k), sb(2, k)) ; b = max(sa(2, k), sb(2, k))
      case (3_I4P)
         f = 4_I4P ; cn = 4_I4P ; a = min(sa(1, k), sb(1, k)) ; b = max(sa(1, k), sb(1, k))
      case default
         f = 1_I4P ; cn = 1_I4P ; a = min(sa(2, k), sb(2, k)) ; b = max(sa(2, k), sb(2, k))
      endselect
      plen(f, p) = plen(f, p) + (b - a) ; pmom(f, p) = pmom(f, p) + 0.5_R8P * (b - a) * (b + a)
      if (sx(k) == 0_I4P) cpiece(cn) = p ! the sub-interval ends at the corner after its edge
   enddo
   do k = 1_I4P, nb
      wpiece(k) = 0_I4P
   enddo
   do t = 1_I4P, nt
      p = rix(find(tslot(t), s, par))
      if (tkind(t) == 1_I4P) then ! a chain, as `follow` walked it
         cur = tstart(t) ; guard = 0_I4P
         do
            wpiece(cur) = p
            if (lt(2, cur) > 0_I4P) exit
            cur = local_start(-lt(2, cur), nb, lt)
            guard = guard + 1_I4P
            if (cur == 0_I4P .or. guard > nb) exit ! cannot happen: the count pass closed this chain
         enddo
      else
         do m = tstart(t), tstart(t) + tn(t) - 1_I4P
            wpiece(loopp(m)) = p
         enddo
      endif
   enddo
   endsubroutine cut_cell_2d_labels

   subroutine cut_cell_2d_measures(u, v, next, x0, x1, y0, y1, flags, closed, alpha, centroid)
   !< The fluid fraction and centroid of the single cell [x0, x1] x [y0, y1] of the curve (u, v, next) with given corner
   !< inside flags (issue #102 step 2: a face of a 3-D cell at an instant within the step, from the section of the facets
   !< around it). As `cut_cell_2d_pieces`, the curve may be any part of a curve that holds every segment meeting the cell
   !< (open chains included) and the flags replace the parity of the far segments; the cell is walked as a one-cell
   !< block and measured by `cut_cell_2d_count`. `closed` false (flags that disagree with the crossings): no measures.
   !< Since issue #105 the device kernels (`cut_cell_2d_walk`, `cut_cell_2d_one_count`) on a host workspace.
   real(R8P),    intent(in)  :: u(:), v(:)     !< Vertices.
   integer(I4P), intent(in)  :: next(:)        !< Next vertex (0: none).
   real(R8P),    intent(in)  :: x0, x1, y0, y1 !< Cell.
   logical,      intent(in)  :: flags(4)       !< Corners inside the solid: (x0,y0), (x1,y0), (x1,y1), (x0,y1).
   logical,      intent(out) :: closed         !< The cell boundary closed.
   real(R8P),    intent(out) :: alpha          !< Fluid fraction.
   real(R8P),    intent(out) :: centroid(2)    !< Fluid centroid (the cell centre without fluid).
   real(R8P),    allocatable :: rw(:)        !< Workspace.
   integer(I4P), allocatable :: iw(:)        !< Workspace.
   integer(I4P)              :: nu, nxs, nps, nrw, niw, np, idum(2, 1) !< Work.
   real(R8P)                 :: rdum(4, 1)   !< Work.

   ! the kernels of the device path (issue #105): the walk's count pass, then the evaluation on a workspace
   nu = size(u, kind=I4P)
   call cut_cell_2d_walk(nu, u, v, next, x0, x1, y0, y1, .false., 0_I4P, 0_I4P, nxs, nps, idum, rdum, idum, rdum, idum, &
                         rdum, idum)
   call cut_cell_2d_one_sizes(nxs, nps, nrw, niw)
   allocate(rw(nrw), iw(niw))
   call cut_cell_2d_one_count(nu, u, v, next, x0, x1, y0, y1, flags(1), flags(2), flags(3), flags(4), nxs, nps, nrw, rw, &
                              niw, iw, closed, alpha, centroid, np)
   endsubroutine cut_cell_2d_measures

   ! One-cell kernels (issue #105, the piece map on the device): the host wrappers above, as routine-seq kernels on a
   ! workspace. The walk of the one-cell block, its crossing index, then the count kernel; no allocation, plain arrays.

   pure integer(I4P) function col1(w, e0, e1)
   !< `col_of` on the two edges of a one-cell block: 0 if w <= e0, 1 if e0 < w <= e1, 2 beyond.
   !$acc routine seq
   real(R8P), intent(in) :: w      !< Coordinate.
   real(R8P), intent(in) :: e0, e1 !< Edges.

   if (w <= e0) then
      col1 = 0_I4P
   elseif (w > e1) then
      col1 = 2_I4P
   else
      col1 = 1_I4P
   endif
   endfunction col1

   pure subroutine cut_cell_2d_walk(nu, u, v, next, x0, x1, y0, y1, fill, mx, mp, nxs, nps, xedge, xpos, xseg, xp, pseg, &
                                    pab, ptg)
   !< Kernel (issue #105): `walk_all` on the one-cell block [x0, x1] x [y0, y1], the same crossing expressions and the
   !< same `corner_side` tie-break, so the same crossings and wall pieces in the same (walk) order with the same tags.
   !< Edges of the one-cell block: 1 the line x0, 2 the line x1, 3 the line y0, 4 the line y1 (`cell_edge_ids`: bottom 3,
   !< right 2, top 4, left 1). Count pass (fill false): nxs, nps only, the arrays (capacities mx, mp) are not touched.
   !< Fill pass: also writes them, as far as the capacities go (a straight segment crosses each line once and lies in the
   !< convex cell along one piece: nxs <= 4 and nps <= 1 per segment).
   !$acc routine seq
   integer(I4P), intent(in)    :: nu                 !< Vertices.
   real(R8P),    intent(in)    :: u(nu), v(nu)       !< Vertices.
   integer(I4P), intent(in)    :: next(nu)           !< Next vertex (0: none).
   real(R8P),    intent(in)    :: x0, x1, y0, y1     !< Cell.
   logical,      intent(in)    :: fill               !< Write the crossings and pieces (else count them).
   integer(I4P), intent(in)    :: mx, mp             !< Capacities.
   integer(I4P), intent(out)   :: nxs, nps           !< Crossings, wall pieces.
   integer(I4P), intent(inout) :: xedge(mx), xseg(mx) !< Crossings: edge, segment (its start vertex).
   real(R8P),    intent(inout) :: xpos(mx), xp(2, mx) !< Crossings: position along the edge, point.
   integer(I4P), intent(inout) :: pseg(mp)           !< Wall pieces: segment.
   real(R8P),    intent(inout) :: pab(4, mp)         !< Wall pieces: end points.
   integer(I4P), intent(inout) :: ptg(2, mp)         !< Wall pieces: end tags (> 0 crossing, < 0 minus the vertex).
   integer(I4P) :: s, w, ci, cj, qi, qj, lv, lh, sd, ta, guard, last_id !< Work.
   real(R8P)    :: pu, pv, qu, qv, du, dv, au, av, cu, cv, xx, yy     !< Work.
   logical      :: vfirst, cb, hasv, hash                             !< Work.

   nxs = 0_I4P ; nps = 0_I4P
   do s = 1_I4P, nu
      w = next(s)
      if (w == 0_I4P) cycle
      pu = u(s) ; pv = v(s) ; qu = u(w) ; qv = v(w) ; du = qu - pu ; dv = qv - pv
      if (du == 0._R8P .and. dv == 0._R8P) cycle
      if (max(pu, qu) < x0 .or. min(pu, qu) > x1 .or. max(pv, qv) < y0 .or. min(pv, qv) > y1) cycle
      ci = col1(pu, x0, x1) ; cj = col1(pv, y0, y1) ; qi = col1(qu, x0, x1) ; qj = col1(qv, y0, y1)
      au = pu ; av = pv ; ta = -s ; guard = 0_I4P ; last_id = 0_I4P
      do
         if (ci == qi .and. cj == qj) then
            if (ci == 1_I4P .and. cj == 1_I4P) then
               nps = nps + 1_I4P
               if (fill .and. nps <= mp) then
                  pseg(nps) = s ; pab(1, nps) = au ; pab(2, nps) = av ; pab(3, nps) = qu ; pab(4, nps) = qv
                  ptg(1, nps) = ta ; ptg(2, nps) = -w
               endif
            endif
            exit
         endif
         guard = guard + 1_I4P
         if (guard > 16_I4P) exit ! 4 (nx + ny) + 8
         hasv = .false. ; hash = .false. ; lv = 0_I4P ; lh = 0_I4P ; xx = 0._R8P ; yy = 0._R8P
         if (du > 0._R8P .and. ci < qi) then
            hasv = .true. ; lv = ci + 1_I4P
         elseif (du < 0._R8P .and. ci > qi) then
            hasv = .true. ; lv = ci
         endif
         if (hasv) xx = merge(x0, x1, lv == 1_I4P)
         if (dv > 0._R8P .and. cj < qj) then
            hash = .true. ; lh = cj + 1_I4P
         elseif (dv < 0._R8P .and. cj > qj) then
            hash = .true. ; lh = cj
         endif
         if (hash) yy = merge(y0, y1, lh == 1_I4P)
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
            cu = xx ; cv = pv + (xx - pu) * (qv - pv) / (qu - pu) ! fossil_line_query, axis 2
            last_id = 0_I4P
            if (cj == 1_I4P) then ! on the edge lv of the block
               nxs = nxs + 1_I4P ; last_id = nxs
               if (fill .and. nxs <= mx) then
                  xedge(nxs) = lv ; xpos(nxs) = cv ; xseg(nxs) = s ; xp(1, nxs) = cu ; xp(2, nxs) = cv
               endif
            endif
         else
            cu = pu + (yy - pv) * (qu - pu) / (qv - pv) ; cv = yy ! fossil_line_query, axis 1
            last_id = 0_I4P
            if (ci == 1_I4P) then ! on the edge 2 + lh of the block
               nxs = nxs + 1_I4P ; last_id = nxs
               if (fill .and. nxs <= mx) then
                  xedge(nxs) = 2_I4P + lh ; xpos(nxs) = cu ; xseg(nxs) = s ; xp(1, nxs) = cu ; xp(2, nxs) = cv
               endif
            endif
         endif
         if (ci == 1_I4P .and. cj == 1_I4P) then
            nps = nps + 1_I4P
            if (fill .and. nps <= mp) then
               pseg(nps) = s ; pab(1, nps) = au ; pab(2, nps) = av ; pab(3, nps) = cu ; pab(4, nps) = cv
               ptg(1, nps) = ta ; ptg(2, nps) = last_id
            endif
         endif
         if (vfirst) then
            if (du > 0._R8P) then
               ci = ci + 1_I4P
            else
               ci = ci - 1_I4P
            endif
         else
            if (dv > 0._R8P) then
               cj = cj + 1_I4P
            else
               cj = cj - 1_I4P
            endif
         endif
         au = cu ; av = cv ; ta = last_id
      enddo
   enddo
   endsubroutine cut_cell_2d_walk

   pure subroutine cut_cell_2d_index(nxs, xedge, xpos, xseg, xstart, xord)
   !< Kernel (issue #105): `index_crossings` of a one-cell block (4 edges): the crossings by edge, then position, then
   !< segment (a total order: a segment crosses a line once), xstart(e) the first of edge e.
   !$acc routine seq
   integer(I4P), intent(in)  :: nxs                 !< Crossings.
   integer(I4P), intent(in)  :: xedge(nxs), xseg(nxs) !< Edges, segments.
   real(R8P),    intent(in)  :: xpos(nxs)           !< Positions.
   integer(I4P), intent(out) :: xstart(5)           !< First crossing of each edge.
   integer(I4P), intent(out) :: xord(nxs)           !< Crossings in order.
   integer(I4P)              :: k, e, i, j, key, fl(4) !< Work.

   do e = 1_I4P, 5_I4P
      xstart(e) = 0_I4P
   enddo
   do k = 1_I4P, nxs
      xstart(xedge(k)) = xstart(xedge(k)) + 1_I4P
   enddo
   key = 1_I4P
   do e = 1_I4P, 4_I4P
      i = xstart(e) ; xstart(e) = key ; key = key + i
   enddo
   xstart(5) = key
   do e = 1_I4P, 4_I4P
      fl(e) = 0_I4P
   enddo
   do k = 1_I4P, nxs
      e = xedge(k)
      xord(xstart(e) + fl(e)) = k ; fl(e) = fl(e) + 1_I4P
   enddo
   do e = 1_I4P, 4_I4P ! insertion sort inside each edge by position, then segment
      do i = xstart(e) + 1_I4P, xstart(e + 1) - 1_I4P
         key = xord(i) ; j = i - 1_I4P
         do while (j >= xstart(e))
            if (xpos(xord(j)) < xpos(key) .or. (xpos(xord(j)) == xpos(key) .and. xseg(xord(j)) <= xseg(key))) exit
            xord(j + 1) = xord(j) ; j = j - 1_I4P
         enddo
         xord(j + 1) = key
      enddo
   enddo
   endsubroutine cut_cell_2d_index

   pure subroutine cut_cell_2d_one_layout(nxs, nps, ro, io)
   !< Start of each part of a one-cell evaluation's workspace (issue #105; sizes in `cut_cell_2d_one_sizes`).
   !< Reals: crossing points xp(2, nxs), positions(nxs), wall piece ends(4, nps), then the count workspace
   !< (`cut_cell_2d_sizes`). Integers: crossing edges(nxs), segments(nxs), wall piece segments(nps), tags(2, nps), edge
   !< starts(5), crossing order(nxs), piece order(nps), piece of each root(S), then the count workspace; S = 4 + nxs +
   !< nps. ro(5), io(10): the total sizes.
   !$acc routine seq
   integer(I4P), intent(in)  :: nxs, nps     !< Crossings, wall pieces (`cut_cell_2d_walk`).
   integer(I4P), intent(out) :: ro(5), io(10) !< Starts.
   integer(I4P)              :: nrc, nic     !< Count workspace.

   call cut_cell_2d_sizes(nxs, nps, 0_I4P, nrc, nic)
   ro(1) = 1_I4P ; ro(2) = ro(1) + 2_I4P * nxs ; ro(3) = ro(2) + nxs ; ro(4) = ro(3) + 4_I4P * nps ; ro(5) = ro(4) + nrc
   io(1) = 1_I4P ; io(2) = io(1) + nxs ; io(3) = io(2) + nxs ; io(4) = io(3) + nps ; io(5) = io(4) + 2_I4P * nps
   io(6) = io(5) + 5_I4P ; io(7) = io(6) + nxs ; io(8) = io(7) + nps ; io(9) = io(8) + 4_I4P + nxs + nps
   io(10) = io(9) + nic
   endsubroutine cut_cell_2d_one_layout

   pure subroutine cut_cell_2d_one_sizes(nxs, nps, nrw, niw)
   !< Workspace of a one-cell evaluation (`cut_cell_2d_one_count`, `cut_cell_2d_one_pieces`) with nxs crossings and nps
   !< wall pieces (the count pass of `cut_cell_2d_walk`), issue #105:
   !<    nrw = 3 nxs + 4 nps + nrc reals,  niw = 4 nxs + 4 nps + 9 + niw_c integers
   !< with nrc, niw_c of `cut_cell_2d_sizes(nxs, nps, 0)`. Pieces: at most 4 + nxs + nps.
   !$acc routine seq
   integer(I4P), intent(in)  :: nxs, nps     !< Crossings, wall pieces.
   integer(I4P), intent(out) :: nrw, niw     !< Reals, integers.
   integer(I4P)              :: ro(5), io(10) !< Layout.

   call cut_cell_2d_one_layout(nxs, nps, ro, io)
   nrw = ro(5) - 1_I4P ; niw = io(10) - 1_I4P
   endsubroutine cut_cell_2d_one_sizes

   subroutine cut_cell_2d_one_count(nu, u, v, next, x0, x1, y0, y1, f00, f10, f11, f01, nxs, nps, nrw, rw, niw, iw, &
                                    closed, alpha, centroid, np)
   !< Kernel (issue #105): the one-cell evaluation of `cut_cell_2d_measures` on a workspace: the walk of the cell
   !< (`cut_cell_2d_walk`, fill pass with the counts nxs, nps of its count pass), the crossing index, then the count
   !< kernel (`cut_cell_2d_count`, no moments, intervals or boundaries), whose workspace keeps the pieces for
   !< `cut_cell_2d_one_pieces`. Bit for bit the host path of issue #102. closed false: the cell boundary does not close
   !< (or the walk disagrees with the counts: cannot happen).
   !$acc routine seq
   integer(I4P), intent(in)    :: nu                  !< Vertices.
   real(R8P),    intent(in)    :: u(nu), v(nu)        !< Vertices.
   integer(I4P), intent(in)    :: next(nu)            !< Next vertex (0: none).
   real(R8P),    intent(in)    :: x0, x1, y0, y1      !< Cell.
   logical,      intent(in)    :: f00, f10, f11, f01  !< Corners inside the solid: (x0,y0), (x1,y0), (x1,y1), (x0,y1).
   integer(I4P), intent(in)    :: nxs, nps            !< Crossings, wall pieces (count pass of the walk).
   integer(I4P), intent(in)    :: nrw, niw            !< Workspace sizes (`cut_cell_2d_one_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)             !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)             !< Integer workspace.
   logical,      intent(out)   :: closed              !< The cell boundary closed.
   real(R8P),    intent(out)   :: alpha, centroid(2)  !< Fluid fraction and centroid.
   integer(I4P), intent(out)   :: np                  !< Pieces.
   integer(I4P) :: ro(5), io(10), k, cx, cp, nivc, nbdc, eid(4), me(2, 1) !< Work.
   real(R8P)    :: aperture(4), wl, wc(2), wn(2)       !< Other measures (unused).

   do k = 1_I4P, nrw
      rw(k) = 0._R8P
   enddo
   do k = 1_I4P, niw
      iw(k) = 0_I4P
   enddo
   call cut_cell_2d_one_layout(nxs, nps, ro, io)
   call cut_cell_2d_walk(nu, u, v, next, x0, x1, y0, y1, .true., nxs, nps, cx, cp, iw(io(1)), rw(ro(2)), iw(io(2)), &
                         rw(ro(1)), iw(io(3)), rw(ro(3)), iw(io(4)))
   if (cx /= nxs .or. cp /= nps) then
      closed = .false. ; alpha = -1._R8P ; centroid(1) = 0.5_R8P * (x0 + x1) ; centroid(2) = 0.5_R8P * (y0 + y1)
      np = 0_I4P
      return
   endif
   call cut_cell_2d_index(nxs, iw(io(1)), rw(ro(2)), iw(io(2)), iw(io(5)), iw(io(6)))
   do k = 1_I4P, nps
      iw(io(7) + k - 1_I4P) = k
   enddo
   eid(1) = 3_I4P ; eid(2) = 2_I4P ; eid(3) = 4_I4P ; eid(4) = 1_I4P ! cell_edge_ids(1, 1, 1, 1)
   me(1, 1) = 0_I4P ; me(2, 1) = 0_I4P
   call cut_cell_2d_count(x0, x1, y0, y1, f00, f10, f11, f01, eid, 4_I4P, nxs, iw(io(5)), iw(io(6)), rw(ro(1)), nps, &
                          0_I4P, nps, iw(io(7)), rw(ro(3)), iw(io(4)), 0_I4P, me, .false., .false., ro(5) - ro(4), &
                          rw(ro(4)), io(10) - io(9), iw(io(9)), closed, alpha, aperture, centroid, wl, wc, wn, np, nivc, nbdc)
   endsubroutine cut_cell_2d_one_count

   subroutine cut_cell_2d_one_pieces(nu, u, v, next, x0, x1, y0, y1, f00, f10, f11, f01, nxs, nps, nrw, rw, niw, iw, &
                                     closed, np, mp, parea, wseg, wpiece, cpiece, plen, pmom)
   !< Kernel (issue #105): `cut_cell_2d_pieces` on a workspace: `cut_cell_2d_one_count`, then the features of its pieces
   !< (`cut_cell_2d_labels`). Outputs as `cut_cell_2d_pieces`, in arrays of capacity mp >= 4 + nxs + nps pieces and nps
   !< wall pieces (wseg, wpiece: one per wall piece, in walk order); with closed false or np = 0 only wseg is written.
   !$acc routine seq
   integer(I4P), intent(in)    :: nu                  !< Vertices.
   real(R8P),    intent(in)    :: u(nu), v(nu)        !< Vertices.
   integer(I4P), intent(in)    :: next(nu)            !< Next vertex (0: none).
   real(R8P),    intent(in)    :: x0, x1, y0, y1      !< Cell.
   logical,      intent(in)    :: f00, f10, f11, f01  !< Corners inside the solid: (x0,y0), (x1,y0), (x1,y1), (x0,y1).
   integer(I4P), intent(in)    :: nxs, nps            !< Crossings, wall pieces (count pass of the walk).
   integer(I4P), intent(in)    :: nrw, niw            !< Workspace sizes (`cut_cell_2d_one_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)             !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)             !< Integer workspace.
   logical,      intent(out)   :: closed              !< The cell boundary closed.
   integer(I4P), intent(out)   :: np                  !< Pieces.
   integer(I4P), intent(in)    :: mp                  !< Capacity of the piece outputs.
   real(R8P),    intent(inout) :: parea(mp)           !< Piece areas.
   integer(I4P), intent(inout) :: wseg(nps)           !< Segment (start vertex) of each wall piece.
   integer(I4P), intent(inout) :: wpiece(nps)         !< Piece on the fluid side of each wall piece.
   integer(I4P), intent(out)   :: cpiece(4)           !< Piece of each corner (0: solid).
   real(R8P),    intent(inout) :: plen(4, mp), pmom(4, mp) !< Face lengths and first moments per piece.
   integer(I4P) :: ro(5), io(10), k, s, ic, rc, cro(15), cio(18) !< Work.
   real(R8P)    :: alpha, centroid(2), dx, dy                  !< Work.

   call cut_cell_2d_one_count(nu, u, v, next, x0, x1, y0, y1, f00, f10, f11, f01, nxs, nps, nrw, rw, niw, iw, closed, &
                              alpha, centroid, np)
   call cut_cell_2d_one_layout(nxs, nps, ro, io)
   do k = 1_I4P, nps
      wseg(k) = iw(io(3) + k - 1_I4P) ; wpiece(k) = 0_I4P
   enddo
   do k = 1_I4P, 4_I4P
      cpiece(k) = 0_I4P
   enddo
   if (.not. closed .or. np == 0_I4P) return
   ic = io(9) ; rc = ro(4)
   if (iw(ic) == 2_I4P) then ! the whole cell fluid, no wall piece and no crossing
      dx = x1 - x0 ; dy = y1 - y0
      parea(1) = dx * dy
      do k = 1_I4P, 4_I4P
         cpiece(k) = 1_I4P
      enddo
      plen(1, 1) = dy ; plen(2, 1) = dy ; plen(3, 1) = dx ; plen(4, 1) = dx
      pmom(1, 1) = 0.5_R8P * dy * (y1 + y0) ; pmom(2, 1) = pmom(1, 1)
      pmom(3, 1) = 0.5_R8P * dx * (x1 + x0) ; pmom(4, 1) = pmom(3, 1)
      return
   endif
   call cut_cell_2d_layout(nxs, nps, 0_I4P, cro, cio)
   s = 4_I4P + nxs + nps
   call cut_cell_2d_labels(nps, s, np, iw(ic), rw(rc + cro(1) - 1_I4P), rw(rc + cro(2) - 1_I4P), rw(rc + cro(3) - 1_I4P), &
                           iw(ic + cio(1) - 1_I4P), iw(ic + cio(2) - 1_I4P), iw(ic + cio(3) - 1_I4P), &
                           iw(ic + cio(4) - 1_I4P), iw(ic + cio(7) - 1_I4P), iw(ic + cio(14) - 1_I4P), &
                           iw(ic + cio(15) - 1_I4P), iw(ic + cio(16) - 1_I4P), iw(ic + cio(17) - 1_I4P), &
                           iw(ic + cio(18) - 1_I4P), iw(io(8)), parea, wpiece, cpiece, plen, pmom)
   endsubroutine cut_cell_2d_one_pieces

   pure subroutine open_cell(x0, x1, y0, y1, closed, alpha, aperture, centroid, npieces, wet_length, wet_normal)
   !< The cell boundary does not close: no measures, no pieces (the wetted centroid keeps the cell centre).
   !$acc routine seq
   real(R8P),    intent(in)  :: x0, x1, y0, y1   !< Cell.
   logical,      intent(out) :: closed           !< Closed.
   real(R8P),    intent(out) :: alpha, aperture(4), centroid(2), wet_length, wet_normal(2) !< Measures.
   integer(I4P), intent(out) :: npieces          !< Pieces.

   closed = .false. ; alpha = -1._R8P ; aperture = -1._R8P
   centroid = [0.5_R8P * (x0 + x1), 0.5_R8P * (y0 + y1)] ; npieces = 0_I4P
   wet_length = 0._R8P ; wet_normal = 0._R8P
   endsubroutine open_cell

   pure integer(I4P) function edge_of_side(e)
   !< Counter-clockwise edge (1 bottom, 2 right, 3 top, 4 left) of side e in the order (-u, +u, -v, +v).
   !$acc routine seq
   integer(I4P), intent(in) :: e !< Side.

   select case (e)
   case (1) ; edge_of_side = 4_I4P
   case (2) ; edge_of_side = 2_I4P
   case (3) ; edge_of_side = 1_I4P
   case default ; edge_of_side = 3_I4P
   endselect
   endfunction edge_of_side

   pure subroutine edge_subs(ed, cs, ce, n0, n1, nxs, xord, xp, check, fend, state, ok, ns, s, sa, sb, sfl, sx, sed, nxc, &
                             nxm, xsub, flen)
   !< Sub-intervals of edge ed (1 bottom, 2 right, 3 top, 4 left) from corner cs to corner ce, counter-clockwise, its
   !< crossings xord(n0:n1); with check, the state at the end corner must be its flag fend.
   !$acc routine seq
   integer(I4P), intent(in)    :: ed                 !< Edge.
   real(R8P),    intent(in)    :: cs(2), ce(2)       !< Start, end corner.
   integer(I4P), intent(in)    :: n0, n1, nxs        !< Crossings of the edge, crossings.
   integer(I4P), intent(in)    :: xord(nxs)          !< Crossings along the edges.
   real(R8P),    intent(in)    :: xp(2, nxs)         !< Crossing points.
   logical,      intent(in)    :: check, fend        !< Check the end corner, its flag.
   logical,      intent(inout) :: state, ok          !< Inside state, the states agree.
   integer(I4P), intent(inout) :: ns                 !< Sub-intervals.
   integer(I4P), intent(in)    :: s                  !< Slots.
   real(R8P),    intent(inout) :: sa(2, s), sb(2, s) !< Sub-intervals.
   integer(I4P), intent(inout) :: sfl(s), sx(s), sed(s) !< Fluid, end crossing, edge.
   integer(I4P), intent(inout) :: nxc                !< Crossing slots.
   integer(I4P), intent(in)    :: nxm                !< Crossing slot capacity.
   integer(I4P), intent(inout) :: xsub(2, nxm)       !< Crossing slots.
   real(R8P),    intent(inout) :: flen(4)            !< Fluid length per edge.
   integer(I4P) :: kk, idx, st, en, dir              !< Work.
   real(R8P)    :: pa_(2), pb_(2)                    !< Work.

   if (ed <= 2) then
      st = n0 ; en = n1 ; dir = 1
   else
      st = n1 ; en = n0 ; dir = -1
   endif
   pa_ = cs
   if (n1 >= n0) then
      do kk = st, en, dir
         idx = xord(kk)
         pb_ = xp(:, idx)
         if (ed == 1 .or. ed == 3) pb_(2) = cs(2)
         if (ed == 2 .or. ed == 4) pb_(1) = cs(1)
         call push(ed, pa_, pb_, .not. state, idx, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
         state = .not. state
         pa_ = pb_
      enddo
   endif
   call push(ed, pa_, ce, .not. state, 0_I4P, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   if (check .and. (state .neqv. fend)) ok = .false.
   endsubroutine edge_subs

   pure subroutine push(ed, aa, bb, fluid, xid, ns, s, sa, sb, sfl, sx, sed, nxc, nxm, xsub, flen)
   !< Append a sub-interval (a slot is always free: the cell's crossings + 4 <= s); record the crossing at its end.
   !$acc routine seq
   integer(I4P), intent(in)    :: ed                 !< Edge.
   real(R8P),    intent(in)    :: aa(2), bb(2)       !< End points.
   logical,      intent(in)    :: fluid              !< Fluid.
   integer(I4P), intent(in)    :: xid                !< Crossing at its end (0: corner).
   integer(I4P), intent(inout) :: ns                 !< Sub-intervals.
   integer(I4P), intent(in)    :: s                  !< Slots.
   real(R8P),    intent(inout) :: sa(2, s), sb(2, s) !< Sub-intervals.
   integer(I4P), intent(inout) :: sfl(s), sx(s), sed(s) !< Fluid, end crossing, edge.
   integer(I4P), intent(inout) :: nxc                !< Crossing slots.
   integer(I4P), intent(in)    :: nxm                !< Crossing slot capacity.
   integer(I4P), intent(inout) :: xsub(2, nxm)       !< Crossing slots.
   real(R8P),    intent(inout) :: flen(4)            !< Fluid length per edge.

   ns = ns + 1 ; sa(:, ns) = aa ; sb(:, ns) = bb ; sx(ns) = xid ; sed(ns) = ed
   sfl(ns) = 0_I4P ; if (fluid) sfl(ns) = 1_I4P
   if (fluid) flen(ed) = flen(ed) + norm2_xy(bb(1) - aa(1), bb(2) - aa(2))
   if (xid > 0) then
      nxc = nxc + 1 ; xsub(1, nxc) = xid ; xsub(2, nxc) = ns
   endif
   endsubroutine push

   pure integer(I4P) function fluid_sub_at(xid, nxc, nxm, xsub, nsb, s, sfl)
   !< The fluid sub-interval next to crossing xid (before or after it, counter-clockwise; 0: none).
   !$acc routine seq
   integer(I4P), intent(in) :: xid, nxc, nxm  !< Crossing id, crossing slots, capacity.
   integer(I4P), intent(in) :: xsub(2, nxm)   !< Crossing slots.
   integer(I4P), intent(in) :: nsb, s         !< Boundary sub-intervals, slots.
   integer(I4P), intent(in) :: sfl(s)         !< Fluid.
   integer(I4P) :: kk, a_, b_                 !< Work.

   fluid_sub_at = 0
   do kk = 1, nxc
      if (xsub(1, kk) /= xid) cycle
      a_ = xsub(2, kk) ; b_ = mod(a_, nsb) + 1
      if (sfl(a_) == 1) then
         fluid_sub_at = a_
      elseif (sfl(b_) == 1) then
         fluid_sub_at = b_
      endif
      return
   enddo
   endfunction fluid_sub_at

   pure integer(I4P) function local_start(vtx, nb, lt)
   !< The first wall piece of the cell that starts at vertex vtx (0: none).
   !$acc routine seq
   integer(I4P), intent(in) :: vtx, nb   !< Vertex, wall pieces.
   integer(I4P), intent(in) :: lt(2, nb) !< Tags.
   integer(I4P) :: kk                    !< Counter.

   local_start = 0
   do kk = 1, nb
      if (lt(1, kk) == -vtx) then
         local_start = kk ; return
      endif
   enddo
   endfunction local_start

   subroutine follow(qq, nb, lp, lt, pvis, plast, g, wl, wnn, good, nm, me, sm, hm, area_c, cgm, cwm, ncs)
   !< Follow the chain from wall piece qq to a piece ending at a crossing.
   !$acc routine seq
   integer(I4P), intent(in)    :: qq, nb            !< First piece, wall pieces.
   real(R8P),    intent(in)    :: lp(4, nb)         !< Wall pieces.
   integer(I4P), intent(in)    :: lt(2, nb)         !< Tags.
   integer(I4P), intent(inout) :: pvis(nb)          !< Visited.
   integer(I4P), intent(out)   :: plast             !< Last piece.
   real(R8P),    intent(out)   :: g(3), wl(4), wnn(2) !< Curve Green, wet (len, mx, my, -), normal.
   logical,      intent(out)   :: good              !< The chain closed on a crossing.
   integer(I4P), intent(in)    :: nm                !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)         !< Exponents.
   real(R8P),    intent(in)    :: sm(2), hm(2), area_c !< Cell centre, size, area.
   real(R8P),    intent(inout) :: cgm(nm), cwm(nm)  !< Chain moments.
   integer(I4P), intent(inout) :: ncs               !< Chain segments.
   integer(I4P) :: cur, lq, guard                   !< Work.

   g = 0._R8P ; wl = 0._R8P ; wnn = 0._R8P ; good = .true. ; cgm = 0._R8P ; cwm = 0._R8P
   cur = qq ; guard = 0
   do
      pvis(cur) = 1_I4P
      call add_curve(lp(:, cur), g, wl, wnn, nm, me, sm, hm, area_c, cgm, cwm)
      ncs = ncs + 1_I4P
      if (lt(2, cur) > 0) exit
      lq = local_start(-lt(2, cur), nb, lt)
      guard = guard + 1
      if (lq == 0 .or. guard > nb) then
         good = .false. ; exit
      endif
      cur = lq
   enddo
   plast = cur
   endsubroutine follow

   subroutine follow_loop(qq, nb, lp, lt, pvis, nloopp, loopp, g, wl, wnn, good, nm, me, sm, hm, area_c, cgm, cwm, ncs)
   !< Follow a closed loop of wall pieces starting at qq (a loop slot is always free: each piece enters one loop).
   !$acc routine seq
   integer(I4P), intent(in)    :: qq, nb            !< First piece, wall pieces.
   real(R8P),    intent(in)    :: lp(4, nb)         !< Wall pieces.
   integer(I4P), intent(in)    :: lt(2, nb)         !< Tags.
   integer(I4P), intent(inout) :: pvis(nb)          !< Visited.
   integer(I4P), intent(inout) :: nloopp            !< Loop pieces.
   integer(I4P), intent(inout) :: loopp(nb)         !< Loop pieces.
   real(R8P),    intent(out)   :: g(3), wl(4), wnn(2) !< Accumulators.
   logical,      intent(out)   :: good              !< The loop closed.
   integer(I4P), intent(in)    :: nm                !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)         !< Exponents.
   real(R8P),    intent(in)    :: sm(2), hm(2), area_c !< Cell centre, size, area.
   real(R8P),    intent(inout) :: cgm(nm), cwm(nm)  !< Loop moments.
   integer(I4P), intent(inout) :: ncs               !< Chain segments.
   integer(I4P) :: lq, guard                        !< Work.

   g = 0._R8P ; wl = 0._R8P ; wnn = 0._R8P ; good = .true. ; cgm = 0._R8P ; cwm = 0._R8P
   lq = qq ; guard = 0
   do
      if (pvis(lq) == 1) exit
      pvis(lq) = 1_I4P
      nloopp = nloopp + 1 ; loopp(nloopp) = lq
      call add_curve(lp(:, lq), g, wl, wnn, nm, me, sm, hm, area_c, cgm, cwm)
      ncs = ncs + 1_I4P
      if (lt(2, lq) >= 0) then
         good = .false. ; exit
      endif
      lq = local_start(-lt(2, lq), nb, lt)
      guard = guard + 1
      if (lq == 0 .or. guard > nb) then
         good = .false. ; exit
      endif
   enddo
   endsubroutine follow_loop

   subroutine add_curve(pab, g, wl, wnn, nm, me, sm, hm, area_c, cgm, cwm)
   !< Curve Green, wetted length, moments and normal of the wall piece pab = (a, b).
   !$acc routine seq
   real(R8P),    intent(in)    :: pab(4)            !< Wall piece.
   real(R8P),    intent(inout) :: g(3), wl(4), wnn(2) !< Accumulators.
   integer(I4P), intent(in)    :: nm                !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)         !< Exponents.
   real(R8P),    intent(in)    :: sm(2), hm(2), area_c !< Cell centre, size, area.
   real(R8P),    intent(inout) :: cgm(nm), cwm(nm)  !< Chain moments.
   real(R8P) :: gg(3), d(2), l                      !< Work.

   call green(pab(1:2), pab(3:4), gg)
   g = g + gg
   d = pab(3:4) - pab(1:2) ; l = norm2_xy(d(1), d(2))
   wl(1) = wl(1) + l ; wl(2:3) = wl(2:3) + l * 0.5_R8P * (pab(1:2) + pab(3:4))
   wnn = wnn + [d(2), -d(1)]
   if (nm > 0_I4P) call curve_moments((pab(1:2) - sm) / hm, (pab(3:4) - sm) / hm, l, nm, me, area_c, cgm, cwm)
   endsubroutine add_curve

   subroutine curve_moments(pa, pb, ln, nm, me, area_c, cgm, cwm)
   !< Green moments and wall moments of the wall segment pa -> pb (scaled coordinates, physical length ln) into the
   !< chain accumulators (issue #70): wall moments by 3-point Gauss on the segment (degree <= 4, exact).
   !$acc routine seq
   real(R8P),    intent(in)    :: pa(2), pb(2), ln  !< End points, length.
   integer(I4P), intent(in)    :: nm                !< Moments.
   integer(I4P), intent(in)    :: me(2, nm)         !< Exponents.
   real(R8P),    intent(in)    :: area_c            !< Cell area.
   real(R8P),    intent(inout) :: cgm(nm), cwm(nm)  !< Chain moments.
   real(R8P), parameter  :: GX(3) = [0.5_R8P - sqrt(0.15_R8P), 0.5_R8P, 0.5_R8P + sqrt(0.15_R8P)] !< Gauss nodes.
   real(R8P), parameter  :: GW(3) = [5._R8P / 18._R8P, 8._R8P / 18._R8P, 5._R8P / 18._R8P] !< Weights.
   real(R8P)             :: pt(2) !< Point.
   integer(I4P)          :: gg, kk !< Work.

   call segment_green_moments(pa, pb, me, area_c, cgm)
   do gg = 1, 3
      pt = pa + GX(gg) * (pb - pa)
      do kk = 1, nm
         cwm(kk) = cwm(kk) + GW(gg) * ln * pt(1) ** me(1, kk) * pt(2) ** me(2, kk)
      enddo
   enddo
   endsubroutine curve_moments

   pure subroutine charge(rr, g, wl, wnn, gm, wm, s, nm, acc, wn, accm, wallm)
   !< Charge a chain or loop to the slot rr: fluid Green is minus the curve's.
   !$acc routine seq
   integer(I4P), intent(in)    :: rr                !< Root or hole slot.
   real(R8P),    intent(in)    :: g(3), wl(3), wnn(2) !< Curve Green, wet (length, mx, my), normal.
   integer(I4P), intent(in)    :: s, nm             !< Slots, moments.
   real(R8P),    intent(in)    :: gm(nm), wm(nm)    !< Curve Green moments, wall moments (issue #70).
   real(R8P),    intent(inout) :: acc(6, s), wn(2, s), accm(nm, s), wallm(nm, s) !< Root accumulators.

   acc(1:3, rr) = acc(1:3, rr) - g
   acc(4:6, rr) = acc(4:6, rr) + wl(1:3)
   wn(:, rr) = wn(:, rr) + wnn
   accm(:, rr) = accm(:, rr) - gm ; wallm(:, rr) = wallm(:, rr) + wm
   endsubroutine charge

   subroutine island_host(kk, nb, lp, lt, loopp, ifirst, nh, hfirst, hn, harea, hroot, nsb, s, sfl, par, sa, sb, nxc, &
                          nxm, xsub, host)
   !< Slot that receives island kk: the smallest hole containing it, else the boundary piece (0: none found).
   !$acc routine seq
   integer(I4P), intent(in)  :: kk, nb               !< Island, wall pieces.
   real(R8P),    intent(in)  :: lp(4, nb)            !< Wall pieces.
   integer(I4P), intent(in)  :: lt(2, nb)            !< Tags.
   integer(I4P), intent(in)  :: loopp(nb), ifirst(nb) !< Loop pieces, island loops.
   integer(I4P), intent(in)  :: nh                   !< Holes.
   integer(I4P), intent(in)  :: hfirst(nb), hn(nb), hroot(nb) !< Hole loops, slots.
   real(R8P),    intent(in)  :: harea(nb)            !< Hole areas.
   integer(I4P), intent(in)  :: nsb, s               !< Boundary sub-intervals, slots.
   integer(I4P), intent(in)  :: sfl(s), par(s)       !< Fluid, parents.
   real(R8P),    intent(in)  :: sa(2, s), sb(2, s)   !< Sub-intervals.
   integer(I4P), intent(in)  :: nxc, nxm             !< Crossing slots, capacity.
   integer(I4P), intent(in)  :: xsub(2, nxm)         !< Crossing slots.
   integer(I4P), intent(out) :: host                 !< Slot.
   real(R8P)    :: pt(2), best                       !< Island vertex, smallest hole area.
   integer(I4P) :: hh, kk2, cnt, nroots, rlast       !< Work.

   host = 0
   pt = lp(1:2, loopp(ifirst(kk)))
   best = huge(1._R8P)
   do hh = 1, nh
      if (inside_loop(pt, hn(hh), loopp(hfirst(hh)), nb, lp) .and. harea(hh) < best) then
         best = harea(hh) ; host = hroot(hh)
      endif
   enddo
   if (host /= 0) return
   nroots = 0 ; rlast = 0
   do kk2 = 1, nsb
      if (sfl(kk2) == 0) cycle
      if (find(kk2, s, par) /= kk2) cycle
      nroots = nroots + 1 ; rlast = kk2
   enddo
   if (nroots == 1) then
      host = rlast ; return
   endif
   ! several boundary pieces: parity of a ray along +u against each piece's boundary
   do kk2 = 1, nsb
      if (sfl(kk2) == 0) cycle
      if (find(kk2, s, par) /= kk2) cycle
      cnt = ray_count(pt, kk2, nb, lp, lt, nsb, s, sfl, par, sa, sb, nxc, nxm, xsub)
      if (mod(cnt, 2_I4P) == 1_I4P) then
         host = kk2 ; return
      endif
   enddo
   endsubroutine island_host

   integer(I4P) function ray_count(pt, rt, nb, lp, lt, nsb, s, sfl, par, sa, sb, nxc, nxm, xsub)
   !< Crossings of the ray from pt along +u with the boundary of the piece of root rt (its fluid sub-intervals and the
   !< wall chains charged to it).
   !$acc routine seq
   real(R8P),    intent(in) :: pt(2)              !< Point.
   integer(I4P), intent(in) :: rt, nb             !< Root, wall pieces.
   real(R8P),    intent(in) :: lp(4, nb)          !< Wall pieces.
   integer(I4P), intent(in) :: lt(2, nb)          !< Tags.
   integer(I4P), intent(in) :: nsb, s             !< Boundary sub-intervals, slots.
   integer(I4P), intent(in) :: sfl(s), par(s)     !< Fluid, parents.
   real(R8P),    intent(in) :: sa(2, s), sb(2, s) !< Sub-intervals.
   integer(I4P), intent(in) :: nxc, nxm           !< Crossing slots, capacity.
   integer(I4P), intent(in) :: xsub(2, nxm)       !< Crossing slots.
   integer(I4P) :: kk, qq, rr                     !< Work.

   ray_count = 0
   do kk = 1, nsb
      if (sfl(kk) == 0) cycle
      if (find(kk, s, par) /= rt) cycle
      if (ray_hits(pt, sa(:, kk), sb(:, kk))) ray_count = ray_count + 1
   enddo
   do qq = 1, nb
      if (lt(1, qq) <= 0) cycle
      rr = fluid_sub_at(lt(1, qq), nxc, nxm, xsub, nsb, s, sfl)
      if (rr == 0) cycle
      if (find(rr, s, par) /= rt) cycle
      call ray_chain(pt, qq, nb, lp, lt, ray_count)
   enddo
   endfunction ray_count

   subroutine ray_chain(pt, pfirst, nb, lp, lt, cnt)
   !< Add the ray crossings of the chain starting at wall piece pfirst.
   !$acc routine seq
   real(R8P),    intent(in)    :: pt(2)      !< Point.
   integer(I4P), intent(in)    :: pfirst, nb !< First piece, wall pieces.
   real(R8P),    intent(in)    :: lp(4, nb)  !< Wall pieces.
   integer(I4P), intent(in)    :: lt(2, nb)  !< Tags.
   integer(I4P), intent(inout) :: cnt        !< Count.
   integer(I4P) :: cur, lq, guard            !< Work.

   cur = pfirst ; guard = 0
   do
      if (ray_hits(pt, lp(1:2, cur), lp(3:4, cur))) cnt = cnt + 1
      if (lt(2, cur) > 0) exit
      lq = local_start(-lt(2, cur), nb, lt)
      guard = guard + 1
      if (lq == 0 .or. guard > nb) exit
      cur = lq
   enddo
   endsubroutine ray_chain

   logical function inside_loop(pt, n, list, nb, lp)
   !< Parity of the ray from pt along +u against a closed loop of n wall pieces.
   !$acc routine seq
   real(R8P),    intent(in) :: pt(2)     !< Point.
   integer(I4P), intent(in) :: n, nb     !< Loop pieces, wall pieces.
   integer(I4P), intent(in) :: list(n)   !< Pieces of the loop.
   real(R8P),    intent(in) :: lp(4, nb) !< Wall pieces.
   integer(I4P) :: kk, cnt               !< Work.

   cnt = 0
   do kk = 1, n
      if (ray_hits(pt, lp(1:2, list(kk)), lp(3:4, list(kk)))) cnt = cnt + 1
   enddo
   inside_loop = mod(cnt, 2_I4P) == 1_I4P
   endfunction inside_loop

   pure integer(I4P) function find(kk, s, par) result(rr)
   !< Union-find root.
   !$acc routine seq
   integer(I4P), intent(in) :: kk, s  !< Element, slots.
   integer(I4P), intent(in) :: par(s) !< Parents.

   rr = kk
   do while (par(rr) /= rr)
      rr = par(rr)
   enddo
   endfunction find

   pure subroutine unite(a_, b_, s, par, nm, acc, wn, accm, wallm)
   !< Union (smaller root wins, deterministic). The loser's accumulators move to the winner, moments included (issue
   !< #95: they were dropped, so a root charged with a chain before losing a later union lost its moments).
   !$acc routine seq
   integer(I4P), intent(in)    :: a_, b_, s, nm     !< Elements, slots, moments.
   integer(I4P), intent(inout) :: par(s)            !< Parents.
   real(R8P),    intent(inout) :: acc(6, s), wn(2, s), accm(nm, s), wallm(nm, s) !< Root accumulators.
   integer(I4P) :: ra, rb                           !< Roots.

   ra = find(a_, s, par) ; rb = find(b_, s, par)
   if (ra == rb) return
   if (ra < rb) then
      par(rb) = ra ; acc(:, ra) = acc(:, ra) + acc(:, rb) ; acc(:, rb) = 0._R8P
      wn(:, ra) = wn(:, ra) + wn(:, rb) ; wn(:, rb) = 0._R8P
      accm(:, ra) = accm(:, ra) + accm(:, rb) ; accm(:, rb) = 0._R8P
      wallm(:, ra) = wallm(:, ra) + wallm(:, rb) ; wallm(:, rb) = 0._R8P
   else
      par(ra) = rb ; acc(:, rb) = acc(:, rb) + acc(:, ra) ; acc(:, ra) = 0._R8P
      wn(:, rb) = wn(:, rb) + wn(:, ra) ; wn(:, ra) = 0._R8P
      accm(:, rb) = accm(:, rb) + accm(:, ra) ; accm(:, ra) = 0._R8P
      wallm(:, rb) = wallm(:, rb) + wallm(:, ra) ; wallm(:, ra) = 0._R8P
   endif
   endsubroutine unite

   subroutine quadrature_out(cc, x, y)
   !< Quadrature rules of the pieces from their oriented boundaries (issue #70): area, signed triangles from the cell
   !< centre (7 points each, degree 5), a full fluid cell the 3 x 3 Gauss rule; wall, 3-point Gauss per wall segment.
   type(cut_cells_2d_t), intent(inout) :: cc         !< Cut cells (boundaries set).
   real(R8P),            intent(in)    :: x(:), y(:) !< Cell edges.
   real(R8P), parameter :: GX(3) = [0.5_R8P - sqrt(0.15_R8P), 0.5_R8P, 0.5_R8P + sqrt(0.15_R8P)] !< Gauss nodes.
   real(R8P), parameter :: GW(3) = [5._R8P / 18._R8P, 8._R8P / 18._R8P, 5._R8P / 18._R8P] !< Weights.
   integer(I4P) :: p, k, nq, nw, i, j, g     !< Work.
   real(R8P)    :: c(2), tp(2, 7), tw(7), ln, gp(2, 9), gww(9) !< Work.
   logical      :: full                      !< A fluid cell.

   allocate(cc%quad_start(cc%np + 1), cc%wall_quad_start(cc%np + 1))
   cc%quad_start(1) = 1_I4P ; cc%wall_quad_start(1) = 1_I4P
   do p = 1, cc%np ! counts
      full = .not. any(cc%piece_seg_wall(cc%piece_seg_start(p):cc%piece_seg_start(p + 1) - 1)) .and. &
             cc%piece_seg_start(p + 1) - cc%piece_seg_start(p) == 4_I4P .and. &
             cc%pieces(cc%piece_cell(1, p), cc%piece_cell(2, p)) == 1
      nq = 7_I4P * (cc%piece_seg_start(p + 1) - cc%piece_seg_start(p)) ; if (full) nq = 9_I4P
      nw = 3_I4P * count(cc%piece_seg_wall(cc%piece_seg_start(p):cc%piece_seg_start(p + 1) - 1))
      cc%quad_start(p + 1) = cc%quad_start(p) + nq ; cc%wall_quad_start(p + 1) = cc%wall_quad_start(p) + nw
   enddo
   allocate(cc%quad_point(2, cc%quad_start(cc%np + 1) - 1), cc%quad_weight(cc%quad_start(cc%np + 1) - 1))
   allocate(cc%wall_quad_point(2, cc%wall_quad_start(cc%np + 1) - 1), cc%wall_quad_weight(cc%wall_quad_start(cc%np + 1) - 1))
   do p = 1, cc%np
      i = cc%piece_cell(1, p) ; j = cc%piece_cell(2, p) ; c = 0.5_R8P * [x(i) + x(i + 1), y(j) + y(j + 1)]
      nq = cc%quad_start(p) ; nw = cc%wall_quad_start(p)
      if (cc%quad_start(p + 1) - nq == 9_I4P .and. cc%piece_seg_start(p + 1) - cc%piece_seg_start(p) == 4_I4P) then
         call tensor_rule([x(i), y(j)], [x(i + 1), y(j + 1)], gp, gww)
         cc%quad_point(:, nq:nq + 8) = gp ; cc%quad_weight(nq:nq + 8) = gww
         cycle
      endif
      do k = cc%piece_seg_start(p), cc%piece_seg_start(p + 1) - 1
         call triangle_rule(c, cc%piece_seg(1:2, k), cc%piece_seg(3:4, k), tp, tw, .true.)
         cc%quad_point(:, nq:nq + 6) = tp ; cc%quad_weight(nq:nq + 6) = tw ; nq = nq + 7_I4P
         if (.not. cc%piece_seg_wall(k)) cycle
         ln = norm2(cc%piece_seg(3:4, k) - cc%piece_seg(1:2, k))
         do g = 1, 3
            cc%wall_quad_point(:, nw) = cc%piece_seg(1:2, k) + GX(g) * (cc%piece_seg(3:4, k) - cc%piece_seg(1:2, k))
            cc%wall_quad_weight(nw) = GW(g) * ln ; nw = nw + 1_I4P
         enddo
      enddo
   enddo
   endsubroutine quadrature_out

   pure function norm2_xy(x, y) result(r)
   !< norm2([x, y]) by the scaled algorithm of gfortran's norm2, written out (issue #100 step 2): nvfortran's device norm2
   !< is sqrt(x*x + y*y) and its host norm2 another formula, 1 ulp apart from each other; this one gives the bits of
   !< gfortran's intrinsic (checked on 4 million pairs from 1e-160 to 1e300, zeros included) on the host and the device.
   !$acc routine seq
   real(R8P), intent(in) :: x, y                !< Components.
   real(R8P)             :: r                   !< Norm.
   real(R8P)             :: scale, res, a, val  !< Work.

   scale = 1._R8P ; res = 0._R8P
   if (x /= 0._R8P) then
      a = abs(x)
      if (a > scale) then
         val = scale / a ; res = (val * val) * res + 1._R8P ; scale = a
      else
         val = a / scale ; res = val * val + res
      endif
   endif
   if (y /= 0._R8P) then
      a = abs(y)
      if (a > scale) then
         val = scale / a ; res = (val * val) * res + 1._R8P ; scale = a
      else
         val = a / scale ; res = val * val + res
      endif
   endif
   r = sqrt(res) * scale
   endfunction norm2_xy

   pure subroutine green(a, b, g)
   !< Green integrals along the straight path a -> b: g(1) = 1/2 (u dv - v du) (area), g(2) = 1/2 u^2 dv
   !< (moment in u), g(3) = -1/2 v^2 du (moment in v).
   !$acc routine seq
   real(R8P), intent(in)  :: a(2), b(2) !< End points.
   real(R8P), intent(out) :: g(3)       !< Integrals.

   g(1) = 0.5_R8P * (a(1) * b(2) - b(1) * a(2))
   g(2) = (b(2) - a(2)) * (a(1) * a(1) + a(1) * b(1) + b(1) * b(1)) / 6._R8P
   g(3) = -(b(1) - a(1)) * (a(2) * a(2) + a(2) * b(2) + b(2) * b(2)) / 6._R8P
   endsubroutine green

   pure logical function ray_hits(pt, a, b)
   !< The ray from pt along +u crosses the segment ab (v = pt(2) + eps rule; exact side test).
   !$acc routine seq
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
   !$acc routine seq
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
   !< Inside flags of the grid nodes: parity of the crossings of their horizontal line left of them (`node_parity` over
   !< the segments of the node's row, `node_rows`; the device twin runs the same kernel per node).
   real(R8P),            intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)  :: next(:)    !< Next vertex.
   real(R8P),            intent(in)  :: x(:), y(:) !< Edges.
   logical, allocatable, intent(out) :: node_in(:,:) !< Flags [nx+1, ny+1].
   integer(I4P), allocatable :: rstart(:), rseg(:) !< Segments crossing each row.
   integer(I4P)              :: i, j, nu, nr       !< Work.

   call node_rows(u, v, next, y, rstart, rseg)
   nu = size(u, kind=I4P) ; nr = size(rseg, kind=I4P)
   allocate(node_in(size(x), size(y)))
   do j = 1, size(y)
      do i = 1, size(x)
         node_in(i, j) = node_parity(x(i), y(j), rstart(j), rstart(j + 1) - 1_I4P, nr, rseg, nu, u, v, next)
      enddo
   enddo
   endsubroutine node_flags

   subroutine node_rows(u, v, next, y, rstart, rseg)
   !< The segments whose end points lie on both sides of the row line v = y(j) (one strictly above), in segment order:
   !< rseg(rstart(j) .. rstart(j+1)-1) (rseg has at least one element).
   real(R8P),                 intent(in)  :: u(:), v(:) !< Vertices.
   integer(I4P),              intent(in)  :: next(:)    !< Next vertex.
   real(R8P),                 intent(in)  :: y(:)       !< Row lines.
   integer(I4P), allocatable, intent(out) :: rstart(:)  !< [size(y) + 1].
   integer(I4P), allocatable, intent(out) :: rseg(:)    !< Segments.
   integer(I4P) :: j, s, n                              !< Work.

   allocate(rstart(size(y) + 1))
   n = 0_I4P
   do j = 1, size(y)
      rstart(j) = n + 1_I4P
      do s = 1, size(u)
         if (next(s) == 0) cycle
         if ((v(s) > y(j)) .eqv. (v(next(s)) > y(j))) cycle
         n = n + 1_I4P
      enddo
   enddo
   rstart(size(y) + 1) = n + 1_I4P
   allocate(rseg(max(1_I4P, n))) ; rseg = 0_I4P
   n = 0_I4P
   do j = 1, size(y)
      do s = 1, size(u)
         if (next(s) == 0) cycle
         if ((v(s) > y(j)) .eqv. (v(next(s)) > y(j))) cycle
         n = n + 1_I4P ; rseg(n) = s
      enddo
   enddo
   endsubroutine node_rows

   pure logical function node_parity(xn, yn, k0, k1, nr, rseg, nu, u, v, next)
   !< Inside flag of the node (xn, yn): parity of the crossings of its row segments rseg(k0:k1) left of it, each by the
   !< exact corner test `corner_side` (a crossing is left of the node when the node is right of an upward segment, left
   !< of a downward one).
   !$acc routine seq
   real(R8P),    intent(in) :: xn, yn             !< Node.
   integer(I4P), intent(in) :: k0, k1             !< Row segments.
   integer(I4P), intent(in) :: nr                 !< Size of rseg.
   integer(I4P), intent(in) :: rseg(nr)           !< Segments by row.
   integer(I4P), intent(in) :: nu                 !< Vertices.
   real(R8P),    intent(in) :: u(nu), v(nu)       !< Vertices.
   integer(I4P), intent(in) :: next(nu)           !< Next vertex.
   integer(I4P)             :: k, s, w, cnt, sd   !< Work.

   cnt = 0_I4P
   do k = k0, k1
      s = rseg(k) ; w = next(s)
      sd = corner_side(u(s), v(s), u(w), v(w), xn, yn)
      if ((v(w) > v(s) .and. sd < 0) .or. (v(w) < v(s) .and. sd > 0)) cnt = cnt + 1_I4P
   enddo
   node_parity = mod(cnt, 2_I4P) == 1_I4P
   endfunction node_parity

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
