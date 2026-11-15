!< FOSSIL, the patched interface of a conservative overset assembly in 2-D and 3-D (issue #106, Wang 1995).
module fossil_overset_patch
!< FOSSIL, the patched interface of a conservative overset assembly in 2-D and 3-D (issue #106, Wang 1995).
!<
!< A minor curvilinear grid overlaps a rectilinear major grid. Its outer boundary (the faces of its blocks named
!< `boundary_name`, default 'overlap', or an explicit face list) cuts the major grid exactly:
!< - the major cells inside the minor region are holes;
!< - the major cells it crosses are cut cells (`cut`, with fluid fractions, apertures W_f and pieces);
!< - the boundary, split at the major grid's edges, gives the smallest faces: each one a part of a minor boundary face
!<   shared by one major piece and one minor cell.
!< The two grids exchange fluxes through the smallest faces (Wang 1995), so the composite conserves exactly what each
!< grid conserves: a major piece closes as sum_f W_f S_f n_f + sum of its smallest faces' `sf_area` = 0, and a minor
!< boundary face equals the sum of its smallest faces (`face_covered`) where it lies inside the major box.
!<
!< The minor boundary is a `curve_object` (`boundary`) whose loops have the minor region on their left (the solid of
!< the cut): the faces are oriented by the handedness of their block (the inverse of `wall_cutter_2d`, which leaves the
!< block interior on the right), chained into loops by exact coincidence of their end points, and checked: closed
!< loops, no two segments crossing or touching (`self_intersections`, exact), orientation matching the exact nesting
!< depth (`check_validity`). Several minor blocks are allowed when their loops are disjoint or nested; overlapping minor
!< regions need a minor-minor hierarchy and are refused (OVERSET_PATCH_OVERLAP: the cut cells count regions by parity,
!< so a lens covered twice would come out as fluid). The cut is `boundary%cut_cells(x, y, cut, boundaries=.true.)`, bit
!< for bit; its `piece_seg_src` labels every wall segment with its curve segment, hence its minor face: no geometric
!< matching, no tolerance.
!<
!< Degenerate pieces (a minor boundary along a grid line: the low-side rule of the cut cells puts it into the cells on
!< its solid side, which become cut cells of zero fluid area with a wall and an aperture on the same edge) are flagged,
!< not removed: `degenerate(p)` and `merge_into(p)`, the piece across the coincident edge (the linked piece of largest
!< shared length that is not degenerate itself). merge_into is 0 for an empty piece, one with no fluid edge interval of
!< positive length (a boundary through a cell corner leaves such a point-like piece: no aperture, no smallest face, no
!< flux; 13 on the diamond of the test). The geometry stays exact and every closure holds as is; the
!< fold (move p's smallest faces to merge_into(p) and drop the shared edge) is a choice of the solver, which needs it
!< because a Wang update divides by the piece area. Smallest faces of zero length are dropped; those of ulp length
!< (crossings computed twice at a grid node) are kept, since the closures need their area vectors.
!<
!< Units and signs: lengths and areas in grid units. `sf_area` = length x unit normal, out of the major piece and into
!< the minor cell (a flux F . sf_area leaves the major piece and enters the minor cell `sf_minor`); per piece,
!< sum sf_area = - cut%piece_wet_normal (bit for bit when the piece has one wall chain or loop, the same terms in the
!< same order; to round-off otherwise). `face_vector` = the minor face's length x unit normal out of the minor cell;
!< `face_covered` = - the sum of its smallest faces' sf_area, in order along the face (less than face_vector where the
!< face leaves the major box: the minor grid's own boundary condition applies to the rest).
!<
!< The major grid's own walls are not part of this step (issue #106 step 1). CPU only.
!<
!< `overlap_weights_2d` (issue #106 step 2) gives the intersection weights of Picard et al. (2021) instead: for every
!< receptor cell V_c of a curvilinear block and every cell V_c' of a rectilinear donor grid, |V_c n V_c'|, the weight
!< |V_c n V_c'| / |V_c| and the centroid of V_c n V_c' (a transfer q_c = sum w (q_c' + grad q_c' . (x_cc' - x_c')) is
!< exact for linear fields). It transfers states, not fluxes: not globally conservative (the patched interface is).
!<
!< 3-D (issue #106 step 3): `overset_patch_3d` and `overlap_weights_3d` on the same scheme. The minor boundary faces are
!< quads split along the shorter diagonal (`wall_cutter`'s rule), so the boundary is a closed triangulated surface
!< whose facets are labelled (block, i, j, k, face, triangle 1-2); the major grid is cut by `cut_cells_3d`, whose
!< `walls=` output gives every clipped facet polygon with the major piece it bounds: the smallest faces, with no
!< geometric matching. The receptors of the weights are the hexahedra of 12 triangles (shorter diagonals, ties broken
!< by the global vertex index, so neighbouring receptors split a shared face alike).
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P
use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, TOPO_BOUNDARY
use fossil_curve_object,      only : curve_object, cut_cells_2d_t, CURVE_VALID
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t
use fossil_cut_graph,         only : cut_links_t, cut_links_2d
use fossil_facet_object,      only : facet_object
use fossil_self_intersection, only : intersection_pair_t
use fossil_solid_angle,       only : solid_angle_oac
use fossil_surface_stl_object, only : surface_stl_object

implicit none
private
public :: overset_patch_2d_t, overset_patch_2d
public :: OVERSET_PATCH_OK, OVERSET_PATCH_BAD_INPUT, OVERSET_PATCH_NO_BOUNDARY, OVERSET_PATCH_OPEN, &
          OVERSET_PATCH_OVERLAP, OVERSET_PATCH_NESTING, OVERSET_PATCH_UNTRUSTED
public :: overset_patch_3d_t, overset_patch_3d
public :: overlap_weights_t, overlap_weights_2d, overlap_weights_3d
public :: OVERLAP_WEIGHTS_OK, OVERLAP_WEIGHTS_BAD_INPUT, OVERLAP_WEIGHTS_DEGENERATE, OVERLAP_WEIGHTS_OUTSIDE, &
          OVERLAP_WEIGHTS_OPEN

integer(I4P), parameter :: OVERSET_PATCH_OK          = 0_I4P !< Done.
integer(I4P), parameter :: OVERSET_PATCH_BAD_INPUT   = 1_I4P !< A minor block or face out of range, or bad grid edges.
integer(I4P), parameter :: OVERSET_PATCH_NO_BOUNDARY = 2_I4P !< No boundary face of non-zero length.
integer(I4P), parameter :: OVERSET_PATCH_OPEN        = 3_I4P !< The faces do not chain into closed loops.
integer(I4P), parameter :: OVERSET_PATCH_OVERLAP     = 4_I4P !< Two loops (or one with itself) cross or touch.
integer(I4P), parameter :: OVERSET_PATCH_NESTING     = 5_I4P !< A loop's orientation does not match its nesting (or is
                                                             !< degenerate: fewer than 3 vertices, zero area; 3-D: a
                                                             !< shell of zero volume, or of the wrong sign).
integer(I4P), parameter :: OVERSET_PATCH_UNTRUSTED   = 6_I4P !< 3-D: done, but some cut cells holding smallest faces have no
                                                             !< trusted pieces (pieces = -1): their faces are listed per cell.
integer(I4P), parameter :: OVERLAP_WEIGHTS_OK         = 0_I4P !< Done (receptor: intersections computed).
integer(I4P), parameter :: OVERLAP_WEIGHTS_BAD_INPUT  = 1_I4P !< Bad donor grid edges or a receptor cell out of range.
integer(I4P), parameter :: OVERLAP_WEIGHTS_DEGENERATE = 2_I4P !< A receptor of zero or negative area (volume), or folded.
integer(I4P), parameter :: OVERLAP_WEIGHTS_OUTSIDE    = 3_I4P !< A receptor with no intersection of positive area (it
                                                              !< lies outside the donor box).
integer(I4P), parameter :: OVERLAP_WEIGHTS_OPEN       = 4_I4P !< A donor cell of the receptor's cut did not close.

type :: overset_patch_2d_t
   !< The patched interface of a minor grid in a rectilinear major grid (see the module notes).
   type(curve_object)        :: boundary            !< The minor boundary loops, the minor region on the left.
   integer(I4P), allocatable :: seg_face(:)         !< Face of each curve segment (by its first vertex) [nv].
   type(cut_cells_2d_t)      :: cut                 !< The major grid cut by the minor region (boundaries on).
   logical,      allocatable :: hole(:,:)           !< Major cells inside the minor region (no piece) [nx, ny].
   logical,      allocatable :: degenerate(:)       !< Major pieces of zero (or negative round-off) area [np].
   integer(I4P), allocatable :: merge_into(:)       !< The piece that should take a degenerate piece's fluxes (0: none, or
                                                    !< an empty piece) [np].
   integer(I4P)              :: ns = 0_I4P          !< Smallest faces.
   integer(I4P), allocatable :: piece_sf_start(:)   !< Smallest faces of piece p: start(p) .. start(p + 1) - 1 [np + 1].
   integer(I4P), allocatable :: sf_piece(:)         !< Major piece of each smallest face [ns].
   integer(I4P), allocatable :: sf_face(:)          !< Minor face (index into face_id) [ns].
   integer(I4P), allocatable :: sf_seg(:)           !< Its boundary segment in cut%piece_seg (exact end points) [ns].
   integer(I4P), allocatable :: sf_minor(:,:)       !< Minor cell face (block, i, j, face 1-4) [4, ns].
   real(R8P),    allocatable :: sf_area(:,:)        !< Length x unit normal out of the major piece [2, ns].
   real(R8P),    allocatable :: sf_centroid(:,:)    !< Midpoint [2, ns].
   integer(I4P)              :: nf = 0_I4P          !< Minor boundary faces.
   integer(I4P), allocatable :: face_id(:,:)        !< (block, i, j, face 1-4: -i, +i, -j, +j) [4, nf].
   real(R8P),    allocatable :: face_vector(:,:)    !< Length x unit normal out of the minor cell [2, nf].
   real(R8P),    allocatable :: face_covered(:,:)   !< - sum of its smallest faces' sf_area [2, nf].
   integer(I4P), allocatable :: face_sf_start(:)    !< Smallest faces of face f: face_sf(start(f) .. start(f + 1) - 1),
   integer(I4P), allocatable :: face_sf(:)          !< in order along the face (minor region on the left) [nf + 1], [ns].
endtype overset_patch_2d_t

type :: overlap_weights_t
   !< Intersection weights of receptor cells against the cells of a rectilinear donor grid (Picard et al. 2021; see the
   !< notes of `overlap_weights_2d`). Volumes are areas in 2-D.
   integer(I4P)              :: nr = 0_I4P            !< Receptors.
   integer(I4P)              :: nw = 0_I4P            !< Receptor-donor intersections.
   integer(I4P), allocatable :: start(:)              !< Donors of receptor r: start(r) .. start(r + 1) - 1 [nr + 1].
   integer(I4P), allocatable :: donor(:,:)            !< Donor cell (i, j) of the donor grid [2, nw] ((i, j, k) [3, nw]).
   real(R8P),    allocatable :: volume(:)             !< |V_c n V_c'| [nw].
   real(R8P),    allocatable :: weight(:)             !< volume / receptor_volume [nw].
   real(R8P),    allocatable :: centroid(:,:)         !< Centroid of V_c n V_c' [2, nw] ([3, nw]).
   real(R8P),    allocatable :: receptor_volume(:)    !< Sum of the receptor's volumes (compensated) [nr].
   real(R8P),    allocatable :: covered(:)            !< receptor_volume / the receptor polygon's area (polyhedron's
                                                      !< volume): 1 exactly when
                                                      !< the receptor lies in the closed donor box, below 1 otherwise [nr].
   integer(I4P), allocatable :: flag(:)               !< OVERLAP_WEIGHTS_* per receptor [nr].
endtype overlap_weights_t

type :: receptor_buffer_t
   !< The intersections of one receptor, before the scan (`overlap_weights_2d`).
   integer(I4P)              :: n = 0_I4P     !< Donors.
   integer(I4P), allocatable :: donor(:,:)    !< Donor cells [2 (3), n].
   real(R8P),    allocatable :: volume(:)     !< Intersection volumes [n].
   real(R8P),    allocatable :: centroid(:,:) !< Intersection centroids [2 (3), n].
endtype receptor_buffer_t

type :: overset_patch_3d_t
   !< The patched interface of a minor grid in a rectilinear major grid in 3-D (see `overset_patch_3d`).
   type(surface_stl_object)  :: boundary            !< The minor boundary, minor region solid (outward normals).
   integer(I4P), allocatable :: tri_face(:)         !< Minor face of each boundary facet [nt].
   integer(I4P), allocatable :: tri_id(:)           !< Its triangle in the face's split (1, 2) [nt].
   type(cut_cells_3d_t)      :: cut                 !< The major grid cut by the minor region.
   type(cut_links_t)         :: links               !< Its piece links (the interior fluxes and merges of a solver).
   logical,      allocatable :: hole(:,:,:)         !< Major cells inside the minor region (no piece) [nx, ny, nz].
   logical,      allocatable :: degenerate(:)       !< Major pieces of zero (or negative round-off) volume [np].
   integer(I4P), allocatable :: merge_into(:)       !< The piece that should take a degenerate piece's fluxes (0: none) [np].
   integer(I4P)              :: ns = 0_I4P          !< Smallest faces.
   integer(I4P), allocatable :: piece_sf_start(:)   !< Smallest faces of piece p: start(p) .. start(p + 1) - 1 [np + 1];
                                                    !< those of the cells without trusted pieces follow, from start(np + 1).
   integer(I4P), allocatable :: sf_piece(:)         !< Major piece of each smallest face (0: a cell of the uc_ lists) [ns].
   integer(I4P), allocatable :: sf_cell(:,:)        !< Its major cell (i, j, k) [3, ns].
   integer(I4P), allocatable :: sf_face(:)          !< Minor face (index into face_id) [ns].
   integer(I4P), allocatable :: sf_facet(:)         !< Boundary facet (the triangle it is a part of) [ns].
   integer(I4P), allocatable :: sf_minor(:,:)       !< Minor cell face (block, i, j, k, face 1-6) [5, ns].
   integer(I4P), allocatable :: sf_tri(:)           !< Triangle of the face's split (1, 2) [ns].
   real(R8P),    allocatable :: sf_area(:,:)        !< Area x unit normal out of the major piece [3, ns].
   real(R8P),    allocatable :: sf_centroid(:,:)    !< Centroid [3, ns].
   real(R8P),    allocatable :: sf_wet(:)           !< Area [ns].
   integer(I4P)              :: nuc = 0_I4P         !< Cells without trusted pieces that hold smallest faces.
   integer(I4P), allocatable :: uc_cell(:,:)        !< Those cells (i, j, k) [3, nuc].
   integer(I4P), allocatable :: uc_sf_start(:)      !< Their smallest faces: uc_sf_start(u) .. uc_sf_start(u + 1) - 1 [nuc + 1].
   integer(I4P)              :: nf = 0_I4P          !< Minor boundary faces.
   integer(I4P), allocatable :: face_id(:,:)        !< (block, i, j, k, face 1-6: -i, +i, -j, +j, -k, +k) [5, nf].
   real(R8P),    allocatable :: face_vector(:,:)    !< 1/2 (c3 - c1) x (c4 - c2), out of the minor cell [3, nf].
   real(R8P),    allocatable :: face_covered(:,:)   !< - sum of its smallest faces' sf_area [3, nf].
   integer(I4P), allocatable :: face_sf_start(:)    !< Smallest faces of face f: face_sf(start(f) .. start(f + 1) - 1),
   integer(I4P), allocatable :: face_sf(:)          !< in smallest-face order [nf + 1], [ns].
endtype overset_patch_3d_t

contains
   subroutine overset_patch_2d(blocks, minor, x, y, patch, boundary_name, faces, status)
   !< The patched interface of the minor blocks `blocks(minor)` in the major grid of cell edges x(nx+1), y(ny+1).
   !<
   !< The minor boundary is the faces of the TOPO_BOUNDARY records named `boundary_name` (default 'overlap') of the minor
   !< blocks, in block, record and index order, or the explicit list `faces(4, nf)` = (block, i, j, face 1-4), whose
   !< blocks must be minor ones. On a status other than OVERSET_PATCH_OK nothing is computed past the boundary (the
   !< curve and face list are kept when they were built).
   type(curvilinear_block_2d_t), intent(in)            :: blocks(:)     !< Blocks.
   integer(I4P),                 intent(in)            :: minor(:)      !< Indexes of the minor blocks.
   real(R8P),                    intent(in)            :: x(:), y(:)    !< Major grid edges, strictly increasing.
   type(overset_patch_2d_t),     intent(inout)         :: patch         !< Patched interface (reset first).
   character(*),                 intent(in),  optional :: boundary_name !< Record name (default 'overlap').
   integer(I4P),                 intent(in),  optional :: faces(:,:)    !< Explicit faces (block, i, j, face) [4, nf].
   integer(I4P),                 intent(out), optional :: status        !< OVERSET_PATCH_*.
   type(overset_patch_2d_t)  :: fresh                     !< Default.
   type(cut_links_t)         :: links                     !< Piece links of the cut.
   real(R8P),    allocatable :: fa(:,:), fb(:,:)          !< Oriented face ends [2, nf] (minor region on the left).
   integer(I4P), allocatable :: order(:), succ(:), npred(:), lsize(:) !< Chaining (sorted starts, successor, predecessors).
   real(R8P),    allocatable :: u(:), v(:)                !< Loop vertices.
   logical,      allocatable :: used(:)                   !< Visited.
   real(R8P),    allocatable :: orient(:)                 !< Handedness of the minor blocks.
   integer(I4P), allocatable :: pairs(:,:)                !< Crossing segment pairs.
   real(R8P)                 :: pp(2), qq(2)              !< Wall segment ends.
   character(len=:), allocatable :: bname                 !< Record name.
   integer(I4P) :: st, nf, nseg, f, k, s, nl, nv, p, q, cnt, defect !< Work.

   if (patch%cut%on_device) call patch%cut%exit_device
   patch = fresh
   st = OVERSET_PATCH_OK
   bname = 'overlap' ; if (present(boundary_name)) bname = boundary_name
   if (any(minor < 1_I4P) .or. any(minor > size(blocks, kind=I4P))) st = OVERSET_PATCH_BAD_INPUT
   ! 1. the faces, oriented with the minor region on the left
   if (st == OVERSET_PATCH_OK) then
      if (present(faces)) then
         call explicit_faces(blocks, minor, faces, patch%face_id, st)
      else
         call record_faces(blocks, minor, bname, patch%face_id)
      endif
   endif
   if (st /= OVERSET_PATCH_OK) then
      call put_status(st, status) ; return
   endif
   allocate(orient(size(blocks)))
   orient = 1._R8P
   do k = 1, size(minor, kind=I4P)
      orient(minor(k)) = blocks(minor(k))%orientation()
   enddo
   nf = size(patch%face_id, 2, kind=I4P) ; patch%nf = nf
   allocate(fa(2, nf), fb(2, nf), patch%face_vector(2, nf), patch%face_covered(2, nf))
   nseg = 0_I4P
   do f = 1, nf
      call face_ends(blocks(patch%face_id(1, f)), orient(patch%face_id(1, f)), patch%face_id(2:4, f), fa(:, f), fb(:, f))
      patch%face_vector(:, f) = [fb(2, f) - fa(2, f), -(fb(1, f) - fa(1, f))]
      if (any(fa(:, f) /= fb(:, f))) nseg = nseg + 1_I4P
   enddo
   patch%face_covered = 0._R8P
   if (nseg == 0_I4P) then
      st = OVERSET_PATCH_NO_BOUNDARY ; call put_status(st, status) ; return
   endif
   ! 2. chain the faces of non-zero length into loops by exact coincidence of their ends
   call chain_faces(nf, fa, fb, order, succ, npred, st)
   if (st /= OVERSET_PATCH_OK) then
      call put_status(st, status) ; return
   endif
   allocate(used(nf), u(nseg), v(nseg), patch%seg_face(nseg), lsize(nseg))
   used = .false. ; nl = 0_I4P ; nv = 0_I4P
   do f = 1, nf
      if (used(f) .or. succ(f) == 0_I4P) cycle
      nl = nl + 1_I4P ; lsize(nl) = 0_I4P ; k = f
      do while (.not. used(k))
         used(k) = .true. ; nv = nv + 1_I4P ; lsize(nl) = lsize(nl) + 1_I4P
         u(nv) = fa(1, k) ; v(nv) = fa(2, k) ; patch%seg_face(nv) = k
         k = succ(k)
      enddo
      if (k /= f) st = OVERSET_PATCH_OPEN ! cannot happen: one predecessor each
   enddo
   if (st /= OVERSET_PATCH_OK) then
      call put_status(st, status) ; return
   endif
   call patch%boundary%adopt_loops(u=u, v=v, loop_sizes=lsize(1:nl), status=s)
   if (s /= 0_I4P) then
      st = OVERSET_PATCH_BAD_INPUT ; call put_status(st, status) ; return
   endif
   ! crossings first (an overlap also breaks the nesting test), then the orientation of every loop against its depth
   call patch%boundary%self_intersections(pairs, k)
   if (k > 0_I4P) then
      st = OVERSET_PATCH_OVERLAP ; call put_status(st, status) ; return
   endif
   defect = patch%boundary%check_validity()
   if (defect /= CURVE_VALID) then
      st = OVERSET_PATCH_NESTING ; call put_status(st, status) ; return
   endif
   ! 3. the cut of the major grid, the minor region solid
   call patch%boundary%cut_cells(x, y, patch%cut, status=s, boundaries=.true.)
   if (s /= 0_I4P) then
      st = OVERSET_PATCH_BAD_INPUT ; call put_status(st, status) ; return
   endif
   patch%hole = patch%cut%pieces == 0_I4P .and. patch%cut%closed
   ! 4. the smallest faces, by major piece (count, then fill)
   cnt = 0_I4P
   do q = 1, size(patch%cut%piece_seg_src, kind=I4P)
      if (patch%cut%piece_seg_src(q) <= 0_I4P) cycle
      if (all(patch%cut%piece_seg(1:2, q) == patch%cut%piece_seg(3:4, q))) cycle
      cnt = cnt + 1_I4P
   enddo
   patch%ns = cnt
   allocate(patch%piece_sf_start(patch%cut%np + 1), patch%sf_piece(cnt), patch%sf_face(cnt), patch%sf_minor(4, cnt), &
            patch%sf_area(2, cnt), patch%sf_centroid(2, cnt), patch%sf_seg(cnt))
   cnt = 0_I4P
   do p = 1, patch%cut%np
      patch%piece_sf_start(p) = cnt + 1_I4P
      do q = patch%cut%piece_seg_start(p), patch%cut%piece_seg_start(p + 1) - 1_I4P
         if (patch%cut%piece_seg_src(q) <= 0_I4P) cycle
         pp = patch%cut%piece_seg(1:2, q) ; qq = patch%cut%piece_seg(3:4, q)
         if (all(pp == qq)) cycle
         cnt = cnt + 1_I4P
         f = patch%seg_face(patch%cut%piece_seg_src(q))
         patch%sf_piece(cnt) = p ; patch%sf_face(cnt) = f ; patch%sf_minor(:, cnt) = patch%face_id(:, f)
         patch%sf_seg(cnt) = q
         ! the wall segment runs b -> a of the curve (fluid on the left): minus its wet normal [d2, -d1], d = b - a
         patch%sf_area(:, cnt) = [-(pp(2) - qq(2)), pp(1) - qq(1)]
         patch%sf_centroid(:, cnt) = 0.5_R8P * (pp + qq)
      enddo
   enddo
   patch%piece_sf_start(patch%cut%np + 1) = cnt + 1_I4P
   ! 5. the same by minor face, in order along the face
   call faces_csr(patch, fa, fb)
   ! 6. degenerate pieces and their merge targets
   allocate(patch%degenerate(patch%cut%np), patch%merge_into(patch%cut%np))
   patch%degenerate = patch%cut%piece_area <= 0._R8P
   patch%merge_into = 0_I4P
   if (any(patch%degenerate)) then
      call cut_links_2d(patch%cut, links)
      do p = 1, patch%cut%np
         if (.not. patch%degenerate(p)) cycle
         k = 0_I4P
         do q = links%start(p), links%start(p + 1) - 1_I4P
            if (patch%degenerate(links%piece(q))) cycle
            if (k == 0_I4P) then
               k = q
            elseif (links%area(q) > links%area(k)) then
               k = q
            endif
         enddo
         if (k > 0_I4P) patch%merge_into(p) = links%piece(k)
      enddo
   endif
   call put_status(st, status)
   endsubroutine overset_patch_2d

   subroutine put_status(st, status)
   !< Return the status if asked.
   integer(I4P),           intent(in)  :: st     !< Status.
   integer(I4P), optional, intent(out) :: status !< Returned status.

   if (present(status)) status = st
   endsubroutine put_status

   subroutine record_faces(blocks, minor, bname, face_id)
   !< The faces of the TOPO_BOUNDARY records named bname of the minor blocks, (block, i, j, face), in block, record and
   !< index order.
   type(curvilinear_block_2d_t), intent(in)               :: blocks(:)    !< Blocks.
   integer(I4P),                 intent(in)               :: minor(:)     !< Minor blocks.
   character(*),                 intent(in)               :: bname        !< Record name.
   integer(I4P), allocatable,    intent(out)              :: face_id(:,:) !< Faces [4, nf].
   integer(I4P) :: pass, mb, b, r, nd, d1, m, n, c(2) !< Work.

   n = 0_I4P
   do pass = 1, 2
      if (pass == 2) allocate(face_id(4, n))
      n = 0_I4P
      do mb = 1, size(minor, kind=I4P)
         b = minor(mb)
         do r = 1, blocks(b)%nrec
            if (blocks(b)%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blocks(b)%rec(r)%name) /= bname) cycle
            nd = (blocks(b)%rec(r)%face + 1_I4P) / 2_I4P ; d1 = 3_I4P - nd
            do m = blocks(b)%rec(r)%lo(d1), blocks(b)%rec(r)%hi(d1) - 1_I4P
               n = n + 1_I4P
               if (pass == 1) cycle
               c(d1) = m
               c(nd) = blocks(b)%rec(r)%lo(nd)
               if (mod(blocks(b)%rec(r)%face, 2_I4P) == 0_I4P) c(nd) = c(nd) - 1_I4P ! a high face: the cell below it
               face_id(:, n) = [b, c(1), c(2), blocks(b)%rec(r)%face]
            enddo
         enddo
      enddo
   enddo
   endsubroutine record_faces

   subroutine explicit_faces(blocks, minor, faces, face_id, st)
   !< The explicit face list, checked: minor blocks, cells and faces in range.
   type(curvilinear_block_2d_t), intent(in)    :: blocks(:)    !< Blocks.
   integer(I4P),                 intent(in)    :: minor(:)     !< Minor blocks.
   integer(I4P),                 intent(in)    :: faces(:,:)   !< Faces [4, nf].
   integer(I4P), allocatable,    intent(out)   :: face_id(:,:) !< Faces [4, nf].
   integer(I4P),                 intent(inout) :: st           !< Status.
   integer(I4P) :: f, b !< Work.

   allocate(face_id(4, 0))
   if (size(faces, 1) /= 4) then
      st = OVERSET_PATCH_BAD_INPUT ; return
   endif
   do f = 1, size(faces, 2, kind=I4P)
      b = faces(1, f)
      if (.not. any(minor == b)) then
         st = OVERSET_PATCH_BAD_INPUT ; return
      endif
      if (faces(2, f) < 1_I4P .or. faces(2, f) > blocks(b)%ni - 1_I4P .or. faces(3, f) < 1_I4P .or. &
          faces(3, f) > blocks(b)%nj - 1_I4P .or. faces(4, f) < 1_I4P .or. faces(4, f) > 4_I4P) then
         st = OVERSET_PATCH_BAD_INPUT ; return
      endif
   enddo
   face_id = faces
   endsubroutine explicit_faces

   subroutine face_ends(blk, orient, cf, a, b)
   !< The ends of face cf = (i, j, face) of a block, oriented with the cell on the left: the counter-clockwise traversal
   !< of the cell (i, j) (i+1, j) (i+1, j+1) (i, j+1) on a right-handed block, reversed on a left-handed one.
   type(curvilinear_block_2d_t), intent(in)  :: blk     !< Block.
   real(R8P),                    intent(in)  :: orient  !< Its handedness (`orientation`).
   integer(I4P),                 intent(in)  :: cf(3)   !< Cell (i, j) and face 1-4 (-i, +i, -j, +j).
   real(R8P),                    intent(out) :: a(2), b(2) !< Ends.
   integer(I4P) :: p(2), q(2), i, j   !< Vertices.
   real(R8P)    :: t(2)               !< Work.

   i = cf(1) ; j = cf(2)
   select case(cf(3))
   case(1)
      p = [i, j + 1] ; q = [i, j]
   case(2)
      p = [i + 1, j] ; q = [i + 1, j + 1]
   case(3)
      p = [i, j] ; q = [i + 1, j]
   case default
      p = [i + 1, j + 1] ; q = [i, j + 1]
   endselect
   a = [blk%u(p(1), p(2)), blk%v(p(1), p(2))] ; b = [blk%u(q(1), q(2)), blk%v(q(1), q(2))]
   if (orient < 0._R8P) then
      t = a ; a = b ; b = t
   endif
   endsubroutine face_ends

   subroutine chain_faces(nf, fa, fb, order, succ, npred, st)
   !< The successor of every face of non-zero length: the face whose start is exactly its end (0 for zero-length faces).
   !< OVERSET_PATCH_OPEN unless every such face has exactly one successor and one predecessor.
   integer(I4P),              intent(in)    :: nf             !< Faces.
   real(R8P),                 intent(in)    :: fa(2, nf), fb(2, nf) !< Oriented ends.
   integer(I4P), allocatable, intent(out)   :: order(:)       !< Faces of non-zero length sorted by start.
   integer(I4P), allocatable, intent(out)   :: succ(:)        !< Successor [nf].
   integer(I4P), allocatable, intent(out)   :: npred(:)       !< Predecessors [nf].
   integer(I4P),              intent(inout) :: st             !< Status.
   integer(I4P), allocatable :: tmp(:)  !< Merge sort work.
   integer(I4P) :: n, f, lo, hi, mid, k, w, i, j, i1, j1, e !< Work.

   allocate(succ(nf), npred(nf)) ; succ = 0_I4P ; npred = 0_I4P
   n = 0_I4P
   do f = 1, nf
      if (any(fa(:, f) /= fb(:, f))) n = n + 1_I4P
   enddo
   allocate(order(n), tmp(n)) ; n = 0_I4P
   do f = 1, nf
      if (any(fa(:, f) /= fb(:, f))) then
         n = n + 1_I4P ; order(n) = f
      endif
   enddo
   ! bottom-up merge sort of the starts (lexicographic, then face index: stable)
   w = 1_I4P
   do while (w < n)
      k = 0_I4P
      do lo = 1, n, 2_I4P * w
         mid = min(lo + w, n + 1_I4P) ; hi = min(lo + 2_I4P * w, n + 1_I4P)
         i = lo ; j = mid ; i1 = mid ; j1 = hi
         do e = lo, hi - 1_I4P
            if (i < i1 .and. j < j1) then
               if (.not. less(fa(:, order(j)), fa(:, order(i)))) then
                  tmp(e) = order(i) ; i = i + 1_I4P
               else
                  tmp(e) = order(j) ; j = j + 1_I4P
               endif
            elseif (i < i1) then
               tmp(e) = order(i) ; i = i + 1_I4P
            else
               tmp(e) = order(j) ; j = j + 1_I4P
            endif
         enddo
      enddo
      order = tmp ; w = 2_I4P * w
   enddo
   ! two faces starting at the same point: not a set of disjoint simple loops
   do k = 2, n
      if (all(fa(:, order(k)) == fa(:, order(k - 1)))) then
         st = OVERSET_PATCH_OPEN ; return
      endif
   enddo
   do k = 1, n
      f = order(k)
      lo = 1_I4P ; hi = n ! the first start not less than fb(f)
      do while (lo < hi)
         mid = (lo + hi) / 2_I4P
         if (less(fa(:, order(mid)), fb(:, f))) then
            lo = mid + 1_I4P
         else
            hi = mid
         endif
      enddo
      if (all(fa(:, order(lo)) == fb(:, f))) then
         succ(f) = order(lo) ; npred(order(lo)) = npred(order(lo)) + 1_I4P
      else
         st = OVERSET_PATCH_OPEN ; return
      endif
   enddo
   do k = 1, n
      if (npred(order(k)) /= 1_I4P) then
         st = OVERSET_PATCH_OPEN ; return
      endif
   enddo
   endsubroutine chain_faces

   pure logical function less(p, q)
   !< Lexicographic order of two points.
   real(R8P), intent(in) :: p(2), q(2) !< Points.

   less = p(1) < q(1) .or. (p(1) == q(1) .and. p(2) < q(2))
   endfunction less

   subroutine faces_csr(patch, fa, fb)
   !< The smallest faces by minor face, in order along the face (by the midpoint's parameter; stable), and the covered
   !< vectors summed in that order.
   type(overset_patch_2d_t), intent(inout) :: patch    !< Patch (smallest faces set).
   real(R8P),                intent(in)    :: fa(:,:), fb(:,:) !< Oriented face ends [2, nf].
   real(R8P),    allocatable :: t(:)   !< Parameters along the faces.
   integer(I4P), allocatable :: pos(:) !< Fill positions.
   integer(I4P) :: f, s, k, m, key     !< Work.
   real(R8P)    :: d(2), tk            !< Work.

   allocate(patch%face_sf_start(patch%nf + 1), patch%face_sf(patch%ns), t(patch%ns), pos(patch%nf))
   patch%face_sf_start = 0_I4P
   do s = 1, patch%ns
      f = patch%sf_face(s) ; patch%face_sf_start(f + 1) = patch%face_sf_start(f + 1) + 1_I4P
   enddo
   patch%face_sf_start(1) = 1_I4P
   do f = 1, patch%nf
      patch%face_sf_start(f + 1) = patch%face_sf_start(f + 1) + patch%face_sf_start(f)
   enddo
   pos = patch%face_sf_start(1:patch%nf)
   do s = 1, patch%ns
      f = patch%sf_face(s)
      d = fb(:, f) - fa(:, f)
      patch%face_sf(pos(f)) = s ; t(pos(f)) = dot_product(patch%sf_centroid(:, s) - fa(:, f), d) ; pos(f) = pos(f) + 1_I4P
   enddo
   do f = 1, patch%nf ! insertion sort per face (a face has a few smallest faces)
      do k = patch%face_sf_start(f) + 1_I4P, patch%face_sf_start(f + 1) - 1_I4P
         key = patch%face_sf(k) ; tk = t(k) ; m = k - 1_I4P
         do while (m >= patch%face_sf_start(f))
            if (t(m) <= tk) exit
            patch%face_sf(m + 1) = patch%face_sf(m) ; t(m + 1) = t(m) ; m = m - 1_I4P
         enddo
         patch%face_sf(m + 1) = key ; t(m + 1) = tk
      enddo
      patch%face_covered(:, f) = 0._R8P
      do k = patch%face_sf_start(f), patch%face_sf_start(f + 1) - 1_I4P
         patch%face_covered(:, f) = patch%face_covered(:, f) - patch%sf_area(:, patch%face_sf(k))
      enddo
   enddo
   endsubroutine faces_csr
   subroutine overlap_weights_2d(block, cells, x, y, w, status)
   !< The intersection weights (Picard et al. 2021) of the receptor cells `cells(:, r)` = (i, j) of a curvilinear block
   !< against the cells of the rectilinear donor grid of cell edges x(nx+1), y(ny+1) (issue #106 step 2).
   !<
   !< Each receptor is the closed polygon of its 4 vertices (consecutive coincident vertices merged: a collapsed edge
   !< gives a triangle), counter-clockwise by the block's handedness. It is cut against the donor sub-grid around its
   !< bounding box as a clockwise hole in a solid box (`curve%cut_cells`, the "direct" variant of the step-0 report),
   !< in coordinates shifted to the sub-grid's low corner: the cut cells sum Green's integrals about the origin, and the
   !< shift keeps their relative error at round-off (1e-15) where unshifted it grows with |x| / h (1e-12 at 100). The
   !< fluid of donor cell c' is then V_c n V_c'. Its area and centroid are the cut's cell totals (`alpha` x the shifted
   !< cell area, `centroid` + the shift), never the per-piece arrays.
   !<
   !< The donor list is the exact cut's: every donor cell with a fluid piece, in (j, i) order of the donor grid. Under
   !< the cut's symbolic perturbation (a point on a grid line belongs to the low side) a receptor side lying on an
   !< interior grid line, or through a grid node, can leave a donor whose intersection has zero area: it is listed with
   !< a volume of zero or round-off (a sliver whose Green sum comes out negative counts 0, so every weight lies in
   !< [0, 1]). The sub-grid starts at the last edge <= the receptor's low bound and ends at the first edge >= its high
   !< bound, so a side on the bounding box adds no such donor.
   !<
   !< `receptor_volume` is the compensated sum of the receptor's volumes, so its weights sum to 1 to round-off;
   !< `covered` is 1 when the receptor lies in the closed donor box (decided exactly on its vertices), else
   !< receptor_volume over the polygon's area (shoelace in the shifted coordinates), held below 1.
   !<
   !< Per receptor `flag`: OVERLAP_WEIGHTS_DEGENERATE (fewer than 3 distinct vertices, zero or negative area in the
   !< block's handedness, sides crossing or folding back: `check_validity` with self-intersections, exact) and
   !< OVERLAP_WEIGHTS_OUTSIDE (no intersection of positive volume) have no donors, receptor_volume and covered 0;
   !< OVERLAP_WEIGHTS_OPEN (a donor cell did not close: not expected) has no donors either. `status` is
   !< OVERLAP_WEIGHTS_BAD_INPUT (nothing computed) for donor edges that are not finite and strictly increasing (at least
   !< one cell) or a receptor out of range, OVERLAP_WEIGHTS_DEGENERATE or OVERLAP_WEIGHTS_OPEN if a receptor has that
   !< flag (the others computed), else OVERLAP_WEIGHTS_OK. Serial, CPU only.
   type(curvilinear_block_2d_t), intent(in)            :: block      !< Receptor block.
   integer(I4P),                 intent(in)            :: cells(:,:) !< Receptor cells (i, j) [2, nr].
   real(R8P),                    intent(in)            :: x(:), y(:) !< Donor grid edges, strictly increasing.
   type(overlap_weights_t),      intent(inout)         :: w          !< Weights (reset first).
   integer(I4P),                 intent(out), optional :: status     !< OVERLAP_WEIGHTS_*.
   type(overlap_weights_t)              :: fresh  !< Default.
   type(receptor_buffer_t), allocatable :: buf(:) !< Per receptor intersections.
   real(R8P)    :: orient                         !< Block handedness.
   integer(I4P) :: st, nr, r, k, n                !< Work.

   w = fresh
   st = OVERLAP_WEIGHTS_OK
   if (size(cells, 1) /= 2 .or. size(x) < 2 .or. size(y) < 2) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (.not. (all(abs(x) <= huge(1._R8P)) .and. all(abs(y) <= huge(1._R8P)))) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (any(x(2:) <= x(:size(x) - 1)) .or. any(y(2:) <= y(:size(y) - 1))) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (size(cells, 2) > 0) then
      if (any(cells(1, :) < 1_I4P) .or. any(cells(1, :) > block%ni - 1_I4P) .or. any(cells(2, :) < 1_I4P) .or. &
          any(cells(2, :) > block%nj - 1_I4P)) st = OVERLAP_WEIGHTS_BAD_INPUT
   endif
   if (st /= OVERLAP_WEIGHTS_OK) then
      call put_status(st, status) ; return
   endif
   nr = size(cells, 2, kind=I4P) ; w%nr = nr
   orient = 1._R8P ; if (nr > 0_I4P) orient = block%orientation()
   allocate(buf(nr), w%start(nr + 1), w%receptor_volume(nr), w%covered(nr), w%flag(nr))
   do r = 1, nr
      call receptor_weights_2d(block, cells(1, r), cells(2, r), orient, x, y, buf(r), w%receptor_volume(r), &
                               w%covered(r), w%flag(r))
   enddo
   ! the scan, then the fill
   w%start(1) = 1_I4P
   do r = 1, nr
      w%start(r + 1) = w%start(r) + buf(r)%n
   enddo
   w%nw = w%start(nr + 1) - 1_I4P
   allocate(w%donor(2, w%nw), w%volume(w%nw), w%weight(w%nw), w%centroid(2, w%nw))
   do r = 1, nr
      n = buf(r)%n ; k = w%start(r)
      if (n == 0_I4P) cycle
      w%donor(:, k:k + n - 1) = buf(r)%donor(:, 1:n)
      w%volume(k:k + n - 1) = buf(r)%volume(1:n)
      w%centroid(:, k:k + n - 1) = buf(r)%centroid(:, 1:n)
      if (w%receptor_volume(r) > 0._R8P) then
         w%weight(k:k + n - 1) = buf(r)%volume(1:n) / w%receptor_volume(r)
      else
         w%weight(k:k + n - 1) = 0._R8P
      endif
   enddo
   if (nr > 0_I4P) then
      if (any(w%flag == OVERLAP_WEIGHTS_DEGENERATE)) then
         st = OVERLAP_WEIGHTS_DEGENERATE
      elseif (any(w%flag == OVERLAP_WEIGHTS_OPEN)) then
         st = OVERLAP_WEIGHTS_OPEN
      endif
   endif
   call put_status(st, status)
   endsubroutine overlap_weights_2d

   subroutine receptor_weights_2d(block, i, j, orient, x, y, b, rvol, cov, flag)
   !< The intersections of receptor cell (i, j) with the donor cells (see `overlap_weights_2d`).
   type(curvilinear_block_2d_t), intent(in)    :: block      !< Receptor block.
   integer(I4P),                 intent(in)    :: i, j       !< Receptor cell.
   real(R8P),                    intent(in)    :: orient     !< Block handedness.
   real(R8P),                    intent(in)    :: x(:), y(:) !< Donor grid edges.
   type(receptor_buffer_t),      intent(inout) :: b          !< Intersections.
   real(R8P),                    intent(out)   :: rvol       !< Receptor volume (sum of the intersections).
   real(R8P),                    intent(out)   :: cov        !< Covered fraction.
   integer(I4P),                 intent(out)   :: flag       !< OVERLAP_WEIGHTS_*.
   type(curve_object)     :: c                   !< Solid box with the receptor as a hole.
   type(cut_cells_2d_t)   :: cut                 !< Its cut of the sub-grid.
   real(R8P), allocatable :: xs(:), ys(:)        !< Shifted sub-grid edges.
   real(R8P)    :: p(2, 4), q(2, 4), qs(2, 4)    !< Vertices, distinct ones, shifted.
   real(R8P)    :: lo(2), hi(2), o(2), blo(2), bhi(2), d, area, comp !< Bounds, shift, box, work.
   integer(I4P) :: n, k, nx, ny, i0, i1, j0, j1, ii, jj, st, m !< Work.
   logical      :: inbox                         !< The receptor lies in the closed donor box.

   b%n = 0_I4P ; rvol = 0._R8P ; cov = 0._R8P ; flag = OVERLAP_WEIGHTS_OK
   p(:, 1) = [block%u(i, j), block%v(i, j)]
   p(:, 2) = [block%u(i + 1, j), block%v(i + 1, j)]
   p(:, 3) = [block%u(i + 1, j + 1), block%v(i + 1, j + 1)]
   p(:, 4) = [block%u(i, j + 1), block%v(i, j + 1)]
   if (orient < 0._R8P) p = p(:, [1, 4, 3, 2])
   n = 0_I4P ! consecutive coincident vertices merged (a collapsed edge)
   do k = 1, 4
      if (n > 0_I4P) then
         if (all(p(:, k) == q(:, n))) cycle
      endif
      n = n + 1_I4P ; q(:, n) = p(:, k)
   enddo
   if (n > 1_I4P) then
      if (all(q(:, n) == q(:, 1))) n = n - 1_I4P
   endif
   if (n < 3_I4P) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   lo = minval(q(:, 1:n), dim=2) ; hi = maxval(q(:, 1:n), dim=2)
   inbox = lo(1) >= x(1) .and. hi(1) <= x(nx + 1) .and. lo(2) >= y(1) .and. hi(2) <= y(ny + 1)
   i0 = edge_below(x, lo(1)) ; i1 = edge_above(x, hi(1))
   j0 = edge_below(y, lo(2)) ; j1 = edge_above(y, hi(2))
   ! the receptor, shifted to the sub-grid's low corner, as a clockwise hole in a solid box around it and the sub-grid
   o = [x(i0), y(j0)]
   do k = 1, n
      qs(:, k) = q(:, k) - o
   enddo
   xs = x(i0:i1) - o(1) ; ys = y(j0:j1) - o(2)
   blo = [min(lo(1) - o(1), xs(1)), min(lo(2) - o(2), ys(1))]
   bhi = [max(hi(1) - o(1), xs(size(xs))), max(hi(2) - o(2), ys(size(ys)))]
   d = maxval(bhi - blo) ; blo = blo - d ; bhi = bhi + d
   call c%adopt_loops(u=[blo(1), bhi(1), bhi(1), blo(1), qs(1, n:1:-1)], &
                      v=[blo(2), blo(2), bhi(2), bhi(2), qs(2, n:1:-1)], loop_sizes=[4_I4P, n], status=st)
   if (st /= 0_I4P) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   if (c%check_validity(check_self_intersections=.true.) /= CURVE_VALID) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   if (i1 <= i0 .or. j1 <= j0) then
      flag = OVERLAP_WEIGHTS_OUTSIDE ; return
   endif
   call c%cut_cells(xs, ys, cut, status=st, intervals=.false.)
   if (st /= 0_I4P) then
      flag = OVERLAP_WEIGHTS_OPEN ; return
   endif
   if (.not. all(cut%closed)) then
      flag = OVERLAP_WEIGHTS_OPEN ; return
   endif
   m = count(cut%pieces > 0_I4P, kind=I4P)
   allocate(b%donor(2, m), b%volume(m), b%centroid(2, m))
   comp = 0._R8P
   do jj = 1, cut%ny
      do ii = 1, cut%nx
         if (cut%pieces(ii, jj) <= 0_I4P) cycle
         b%n = b%n + 1_I4P
         b%donor(:, b%n) = [i0 + ii - 1_I4P, j0 + jj - 1_I4P]
         b%volume(b%n) = max(0._R8P, cut%alpha(ii, jj) * ((xs(ii + 1) - xs(ii)) * (ys(jj + 1) - ys(jj))))
         b%centroid(:, b%n) = cut%centroid(:, ii, jj) + o
         call neumaier(rvol, comp, b%volume(b%n))
      enddo
   enddo
   rvol = rvol + comp
   if (.not. (rvol > 0._R8P)) then
      rvol = 0._R8P ; flag = OVERLAP_WEIGHTS_OUTSIDE ; b%n = 0_I4P ; return
   endif
   if (inbox) then
      cov = 1._R8P
   else
      area = 0._R8P
      do k = 2, n - 1
         area = area + ((qs(1, k) - qs(1, 1)) * (qs(2, k + 1) - qs(2, 1)) - (qs(1, k + 1) - qs(1, 1)) * (qs(2, k) - qs(2, 1)))
      enddo
      area = 0.5_R8P * area
      cov = min(rvol / area, nearest(1._R8P, -1._R8P))
   endif
   endsubroutine receptor_weights_2d

   pure integer(I4P) function edge_below(x, t) result(k)
   !< The last edge x(k) <= t (1 if none).
   real(R8P), intent(in) :: x(:) !< Edges, increasing.
   real(R8P), intent(in) :: t    !< Value.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   lo = 1_I4P ; hi = size(x, kind=I4P)
   if (x(1) > t) then
      k = 1_I4P ; return
   endif
   do while (lo < hi) ! invariant: x(lo) <= t
      mid = (lo + hi + 1_I4P) / 2_I4P
      if (x(mid) <= t) then
         lo = mid
      else
         hi = mid - 1_I4P
      endif
   enddo
   k = lo
   endfunction edge_below

   pure integer(I4P) function edge_above(x, t) result(k)
   !< The first edge x(k) >= t (size(x) if none).
   real(R8P), intent(in) :: x(:) !< Edges, increasing.
   real(R8P), intent(in) :: t    !< Value.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   lo = 1_I4P ; hi = size(x, kind=I4P)
   if (x(hi) < t) then
      k = hi ; return
   endif
   do while (lo < hi) ! invariant: x(hi) >= t
      mid = (lo + hi) / 2_I4P
      if (x(mid) >= t) then
         hi = mid
      else
         lo = mid + 1_I4P
      endif
   enddo
   k = hi
   endfunction edge_above

   pure subroutine neumaier(s, comp, a)
   !< s += a with Neumaier's compensation (the sum is s + comp).
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
   ! ------------------------------------------------------------------------------------------------------------- 3-D
   subroutine overset_patch_3d(blocks, minor, x, y, z, patch, boundary_name, faces, status)
   !< The patched interface of the minor blocks `blocks(minor)` in the major grid of cell edges x(nx+1), y(ny+1), z(nz+1)
   !< (issue #106 step 3; the 3-D twin of `overset_patch_2d`).
   !<
   !< The minor boundary is the faces of the TOPO_BOUNDARY records named `boundary_name` (default 'overlap') of the minor
   !< blocks, in block, record and index order (the index along the record's second in-face axis outer), or the explicit
   !< list `faces(5, nf)` = (block, i, j, k, face 1-6), whose blocks must be minor ones. Each face is oriented out of its
   !< cell by the block's handedness (cell (i, j, k), face of axis a: the corners p, p + e1, p + e1 + e2, p + e2 with
   !< (e1, e2) the next two axes in cyclic order, reversed on a low face and on a left-handed block) and split along its
   !< shorter diagonal (c1 c3 on a tie: `wall_cutter`'s rule); exactly collapsed triangles (a pole) are dropped. The
   !< boundary is then a surface with the minor region solid (outward normals: the inverse of `wall_cutter`, which leaves
   !< the block interior fluid), checked: closed (every directed edge, by exact coordinates, has exactly one reverse:
   !< OVERSET_PATCH_OPEN), no two facets crossing or touching (`find_self_intersections` with touching pairs, exact:
   !< OVERSET_PATCH_OVERLAP; overlapping minor regions need a minor-minor hierarchy: the cut counts regions by parity),
   !< every shell's signed volume of the sign of its nesting depth (positive at even depth: OVERSET_PATCH_NESTING; the
   !< depth by the winding number of the other shells at one of its vertices). On a status other than OVERSET_PATCH_OK
   !< or OVERSET_PATCH_UNTRUSTED nothing is computed past the boundary (the face list and surface are kept when built).
   !<
   !< The cut is `cut_cells_3d(boundary, x, y, z, cut, links=, walls=)`, bit for bit what `cut_cells_3d` gives with links
   !< (kept: `links`);
   !< its wall terms are the smallest faces: the part of a boundary facet in one major piece, labelled by the cut itself
   !< (`cut_walls_3d_t`: facet and piece), so the facet gives the minor face (`tri_face`) with no geometric matching.
   !< sf_area = - the term's vector area (out of the major piece, into the minor cell); per piece, sum sf_area =
   !< - cut%piece_wet_normal to round-off, bit for bit when the piece's polygons after its first are single fan
   !< triangles (the cut adds the fans in one running sum). Terms of zero vector area are dropped.
   !<
   !< Untrusted major cut cells (pieces = -1, `cut_cells_3d`'s refusal of a split it cannot certify) have no pieces: their
   !< smallest faces are kept, per cell (`uc_cell`, `uc_sf_start`, sf_piece = 0, after the pieces' faces), and the status
   !< is OVERSET_PATCH_UNTRUSTED with everything else computed: the solver takes such a cell as one piece (its totals,
   !< alpha and centroid, are right) or refines the major grid. None occur on the validation spheres (0 of 25 688 cut
   !< cells at 96^3, issue #106 step 0).
   !<
   !< `degenerate` / `merge_into` as in 2-D (a minor boundary on a grid plane makes zero-volume pieces on its solid side;
   !< merge_into: the linked piece of largest shared area that is not degenerate, 0 if none). `face_vector` is
   !< 1/2 (c3 - c1) x (c4 - c2), the sum of the two triangles of either split, so face_covered equals it to round-off
   !< where the face lies inside the major box. Units: areas in grid units squared. CPU only.
   type(curvilinear_block_t),  intent(in)            :: blocks(:)     !< Blocks.
   integer(I4P),               intent(in)            :: minor(:)      !< Indexes of the minor blocks.
   real(R8P),                  intent(in)            :: x(:), y(:), z(:) !< Major grid edges, strictly increasing.
   type(overset_patch_3d_t),   intent(inout)         :: patch         !< Patched interface (reset first).
   character(*),               intent(in),  optional :: boundary_name !< Record name (default 'overlap').
   integer(I4P),               intent(in),  optional :: faces(:,:)    !< Explicit faces (block, i, j, k, face) [5, nf].
   integer(I4P),               intent(out), optional :: status        !< OVERSET_PATCH_*.
   type(overset_patch_3d_t)               :: fresh      !< Default.
   type(cut_walls_3d_t)                   :: walls      !< Labelled wall terms of the cut.
   type(facet_object),        allocatable :: fct(:)     !< Boundary facets.
   type(intersection_pair_t), allocatable :: pairs(:)   !< Crossing or touching facet pairs.
   real(R8P),                 allocatable :: tri(:,:,:), orient(:) !< Boundary triangles, handedness of the blocks.
   integer(I4P),              allocatable :: tf(:), ti(:) !< Face and triangle id of each triangle.
   character(len=:),          allocatable :: bname      !< Record name.
   real(R8P)    :: c(3, 4), t2(3, 3, 2)                 !< Face corners, its triangles.
   integer(I4P) :: st, nf, nt, f, k, m, n2, tid(2), s   !< Work.
   logical      :: good                                 !< Check.

   if (patch%cut%on_device) call patch%cut%exit_device
   patch = fresh
   st = OVERSET_PATCH_OK
   bname = 'overlap' ; if (present(boundary_name)) bname = boundary_name
   if (any(minor < 1_I4P) .or. any(minor > size(blocks, kind=I4P))) st = OVERSET_PATCH_BAD_INPUT
   ! 1. the faces, oriented out of their cells, split into triangles
   if (st == OVERSET_PATCH_OK) then
      if (present(faces)) then
         call explicit_faces_3d(blocks, minor, faces, patch%face_id, st)
      else
         call record_faces_3d(blocks, minor, bname, patch%face_id)
      endif
   endif
   if (st /= OVERSET_PATCH_OK) then
      call put_status(st, status) ; return
   endif
   allocate(orient(size(blocks)))
   orient = 1._R8P
   do k = 1, size(minor, kind=I4P)
      orient(minor(k)) = blocks(minor(k))%orientation()
   enddo
   nf = size(patch%face_id, 2, kind=I4P) ; patch%nf = nf
   allocate(tri(3, 3, 2 * nf), tf(2 * nf), ti(2 * nf), patch%face_vector(3, nf), patch%face_covered(3, nf))
   patch%face_covered = 0._R8P
   nt = 0_I4P
   do f = 1, nf
      call face_corners_3d(blocks(patch%face_id(1, f)), orient(patch%face_id(1, f)), patch%face_id(2:5, f), c)
      patch%face_vector(:, f) = 0.5_R8P * cross_3(c(:, 3) - c(:, 1), c(:, 4) - c(:, 2))
      call split_quad(c, [0_I4P, 0_I4P, 0_I4P, 0_I4P], t2, tid, n2)
      do m = 1, n2
         nt = nt + 1_I4P ; tri(:, :, nt) = t2(:, :, m) ; tf(nt) = f ; ti(nt) = tid(m)
      enddo
   enddo
   if (nt == 0_I4P) then
      st = OVERSET_PATCH_NO_BOUNDARY ; call put_status(st, status) ; return
   endif
   patch%tri_face = tf(1:nt) ; patch%tri_id = ti(1:nt)
   ! 2. closed: every directed edge has exactly one reverse
   call surface_closed(nt, tri, good)
   if (.not. good) then
      st = OVERSET_PATCH_OPEN ; call put_status(st, status) ; return
   endif
   allocate(fct(nt))
   do m = 1, nt
      do k = 1, 3
         fct(m)%vertex(k) = vector_R8P(tri(1, k, m), tri(2, k, m), tri(3, k, m))
      enddo
   enddo
   call patch%boundary%adopt_facets(facets=fct)
   ! 3. no crossing or touching facets, then the orientation of every shell against its depth
   call patch%boundary%find_self_intersections(pairs, include_touching=.true.)
   if (size(pairs) > 0) then
      st = OVERSET_PATCH_OVERLAP ; call put_status(st, status) ; return
   endif
   call shells_nested(patch%boundary, nt, tri, good)
   if (.not. good) then
      st = OVERSET_PATCH_NESTING ; call put_status(st, status) ; return
   endif
   ! 4. the cut of the major grid, the minor region solid, with its labelled wall terms
   call cut_cells_3d(patch%boundary, x, y, z, patch%cut, status=s, links=patch%links, walls=walls)
   if (s /= 0_I4P) then
      st = OVERSET_PATCH_BAD_INPUT ; call put_status(st, status) ; return
   endif
   patch%hole = patch%cut%pieces == 0_I4P .and. patch%cut%closed
   ! 5. the smallest faces, by piece, then by minor face
   call smallest_faces_3d(patch, walls)
   call faces_csr_3d(patch)
   ! 6. degenerate pieces and their merge targets
   allocate(patch%degenerate(patch%cut%np), patch%merge_into(patch%cut%np))
   patch%degenerate = patch%cut%piece_volume <= 0._R8P
   patch%merge_into = 0_I4P
   if (any(patch%degenerate)) call merge_degenerate(patch%cut%np, patch%degenerate, patch%links, patch%merge_into)
   if (patch%nuc > 0_I4P) st = OVERSET_PATCH_UNTRUSTED
   call put_status(st, status)
   endsubroutine overset_patch_3d

   subroutine merge_degenerate(np, degenerate, links, merge_into)
   !< The merge target of every degenerate piece: the linked piece of largest shared area that is not degenerate.
   integer(I4P),      intent(in)    :: np             !< Pieces.
   logical,           intent(in)    :: degenerate(np) !< Degenerate pieces.
   type(cut_links_t), intent(in)    :: links          !< Links.
   integer(I4P),      intent(inout) :: merge_into(np) !< Targets (0: none).
   integer(I4P)                     :: p, q, k        !< Work.

   do p = 1, np
      if (.not. degenerate(p)) cycle
      k = 0_I4P
      do q = links%start(p), links%start(p + 1) - 1_I4P
         if (degenerate(links%piece(q))) cycle
         if (k == 0_I4P) then
            k = q
         elseif (links%area(q) > links%area(k)) then
            k = q
         endif
      enddo
      if (k > 0_I4P) merge_into(p) = links%piece(k)
   enddo
   endsubroutine merge_degenerate

   subroutine record_faces_3d(blocks, minor, bname, face_id)
   !< The faces of the TOPO_BOUNDARY records named bname of the minor blocks, (block, i, j, k, face), in block, record and
   !< index order (the second in-face axis outer).
   type(curvilinear_block_t), intent(in)  :: blocks(:)    !< Blocks.
   integer(I4P),              intent(in)  :: minor(:)     !< Minor blocks.
   character(*),              intent(in)  :: bname        !< Record name.
   integer(I4P), allocatable, intent(out) :: face_id(:,:) !< Faces [5, nf].
   integer(I4P) :: pass, mb, b, r, nd, d1, d2, m1, m2, n, c(3) !< Work.

   n = 0_I4P
   do pass = 1, 2
      if (pass == 2) allocate(face_id(5, n))
      n = 0_I4P
      do mb = 1, size(minor, kind=I4P)
         b = minor(mb)
         do r = 1, blocks(b)%nrec
            if (blocks(b)%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blocks(b)%rec(r)%name) /= bname) cycle
            nd = (blocks(b)%rec(r)%face + 1_I4P) / 2_I4P ; d1 = mod(nd, 3_I4P) + 1_I4P ; d2 = mod(nd + 1_I4P, 3_I4P) + 1_I4P
            do m2 = blocks(b)%rec(r)%lo(d2), blocks(b)%rec(r)%hi(d2) - 1_I4P
               do m1 = blocks(b)%rec(r)%lo(d1), blocks(b)%rec(r)%hi(d1) - 1_I4P
                  n = n + 1_I4P
                  if (pass == 1) cycle
                  c(d1) = m1 ; c(d2) = m2
                  c(nd) = blocks(b)%rec(r)%lo(nd)
                  if (mod(blocks(b)%rec(r)%face, 2_I4P) == 0_I4P) c(nd) = c(nd) - 1_I4P ! a high face: the cell below it
                  face_id(:, n) = [b, c(1), c(2), c(3), blocks(b)%rec(r)%face]
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine record_faces_3d

   subroutine explicit_faces_3d(blocks, minor, faces, face_id, st)
   !< The explicit face list, checked: minor blocks, cells and faces in range.
   type(curvilinear_block_t), intent(in)    :: blocks(:)    !< Blocks.
   integer(I4P),              intent(in)    :: minor(:)     !< Minor blocks.
   integer(I4P),              intent(in)    :: faces(:,:)   !< Faces [5, nf].
   integer(I4P), allocatable, intent(out)   :: face_id(:,:) !< Faces [5, nf].
   integer(I4P),              intent(inout) :: st           !< Status.
   integer(I4P) :: f, b !< Work.

   allocate(face_id(5, 0))
   if (size(faces, 1) /= 5) then
      st = OVERSET_PATCH_BAD_INPUT ; return
   endif
   do f = 1, size(faces, 2, kind=I4P)
      b = faces(1, f)
      if (.not. any(minor == b)) then
         st = OVERSET_PATCH_BAD_INPUT ; return
      endif
      if (faces(2, f) < 1_I4P .or. faces(2, f) > blocks(b)%ni - 1_I4P .or. faces(3, f) < 1_I4P .or. &
          faces(3, f) > blocks(b)%nj - 1_I4P .or. faces(4, f) < 1_I4P .or. faces(4, f) > blocks(b)%nk - 1_I4P .or. &
          faces(5, f) < 1_I4P .or. faces(5, f) > 6_I4P) then
         st = OVERSET_PATCH_BAD_INPUT ; return
      endif
   enddo
   face_id = faces
   endsubroutine explicit_faces_3d

   pure subroutine face_offsets(fc, orient, off)
   !< The corners of face fc (1-6: -i, +i, -j, +j, -k, +k) of a cell as offsets (0, 1) from its low vertex, in the
   !< order that orients the face out of the cell on a right-handed block: p, p + e1, p + e1 + e2, p + e2 on a high face
   !< ((e1, e2) the next two axes in cyclic order), p, p + e2, p + e1 + e2, p + e1 on a low one; reversed if orient < 0.
   integer(I4P), intent(in)  :: fc        !< Face.
   real(R8P),    intent(in)  :: orient    !< Handedness.
   integer(I4P), intent(out) :: off(3, 4) !< Offsets of the corners.
   integer(I4P)              :: nd, d1, d2, t(3) !< Work.

   nd = (fc + 1_I4P) / 2_I4P ; d1 = mod(nd, 3_I4P) + 1_I4P ; d2 = mod(nd + 1_I4P, 3_I4P) + 1_I4P
   off = 0_I4P
   if (mod(fc, 2_I4P) == 0_I4P) off(nd, :) = 1_I4P
   off(d1, 2) = 1_I4P ; off(d1, 3) = 1_I4P ; off(d2, 3) = 1_I4P ; off(d2, 4) = 1_I4P
   if (mod(fc, 2_I4P) == 1_I4P) then ! low face: p, p + e2, p + e1 + e2, p + e1
      t = off(:, 2) ; off(:, 2) = off(:, 4) ; off(:, 4) = t
   endif
   if (orient < 0._R8P) then
      t = off(:, 2) ; off(:, 2) = off(:, 4) ; off(:, 4) = t
   endif
   endsubroutine face_offsets

   subroutine face_corners_3d(blk, orient, cf, c)
   !< The corners of face cf = (i, j, k, face) of a block, oriented out of the cell (`face_offsets`).
   type(curvilinear_block_t), intent(in)  :: blk     !< Block.
   real(R8P),                 intent(in)  :: orient  !< Its handedness.
   integer(I4P),              intent(in)  :: cf(4)   !< Cell (i, j, k) and face 1-6.
   real(R8P),                 intent(out) :: c(3, 4) !< Corners.
   integer(I4P)                           :: off(3, 4), m, p(3) !< Work.

   call face_offsets(cf(4), orient, off)
   do m = 1, 4
      p = cf(1:3) + off(:, m)
      c(:, m) = [blk%x(p(1), p(2), p(3)), blk%y(p(1), p(2), p(3)), blk%z(p(1), p(2), p(3))]
   enddo
   endsubroutine face_corners_3d

   pure subroutine split_quad(c, g, t2, tid, n2)
   !< The quad c(:, 1:4) (cyclic) split along its shorter diagonal, orientation kept: c1 c3 gives the triangles
   !< (1, 2, 3), (1, 3, 4), c2 c4 gives (1, 2, 4), (2, 3, 4). On a tie of the computed lengths the diagonal through the
   !< least of the vertex indices g (all equal: c1 c3, `wall_cutter`'s rule). Exactly collapsed triangles (zero cross
   !< product: a pole) are dropped; tid: the kept triangles' ids (1, 2).
   real(R8P),    intent(in)  :: c(3, 4)      !< Corners.
   integer(I4P), intent(in)  :: g(4)         !< Global vertex indices (ties).
   real(R8P),    intent(out) :: t2(3, 3, 2)  !< Kept triangles.
   integer(I4P), intent(out) :: tid(2)       !< Their ids.
   integer(I4P), intent(out) :: n2           !< Kept.
   real(R8P)                 :: d13, d24, w(3, 3, 2) !< Work.
   logical                   :: use13        !< Diagonal c1 c3.
   integer(I4P)              :: m            !< Work.

   d13 = norm2(c(:, 3) - c(:, 1)) ; d24 = norm2(c(:, 4) - c(:, 2))
   if (d13 < d24) then
      use13 = .true.
   elseif (d24 < d13) then
      use13 = .false.
   else
      use13 = min(g(1), g(3)) <= min(g(2), g(4))
   endif
   if (use13) then
      w(:, :, 1) = c(:, [1, 2, 3]) ; w(:, :, 2) = c(:, [1, 3, 4])
   else
      w(:, :, 1) = c(:, [1, 2, 4]) ; w(:, :, 2) = c(:, [2, 3, 4])
   endif
   n2 = 0_I4P ; tid = 0_I4P
   do m = 1, 2
      if (norm2(cross_3(w(:, 2, m) - w(:, 1, m), w(:, 3, m) - w(:, 1, m))) == 0._R8P) cycle
      n2 = n2 + 1_I4P ; t2(:, :, n2) = w(:, :, m) ; tid(n2) = m
   enddo
   endsubroutine split_quad

   pure function cross_3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c(1) = a(2) * b(3) - a(3) * b(2) ; c(2) = a(3) * b(1) - a(1) * b(3) ; c(3) = a(1) * b(2) - a(2) * b(1)
   endfunction cross_3

   subroutine surface_closed(nt, tri, ok)
   !< Every directed edge of the triangles (exact coordinates) is unique and has exactly one reverse: a closed, consistently
   !< oriented surface, every edge shared by exactly two triangles.
   integer(I4P), intent(in)  :: nt             !< Triangles.
   real(R8P),    intent(in)  :: tri(3, 3, nt)  !< Triangles.
   logical,      intent(out) :: ok             !< Closed.
   real(R8P),    allocatable :: key(:,:)       !< Edges (from, to) [6, 3 nt].
   integer(I4P), allocatable :: order(:), tmp(:) !< Sorted edges, merge work.
   integer(I4P) :: n, t, m, w, lo, hi, mid, i, j, i1, j1, e, k !< Work.
   real(R8P)    :: rv(6)                        !< Reversed edge.

   n = 3_I4P * nt
   allocate(key(6, n), order(n), tmp(n))
   do t = 1, nt
      do m = 1, 3
         k = 3_I4P * (t - 1_I4P) + m
         key(1:3, k) = tri(:, m, t) ; key(4:6, k) = tri(:, mod(m, 3_I4P) + 1_I4P, t) ; order(k) = k
      enddo
   enddo
   w = 1_I4P ! bottom-up merge sort (stable)
   do while (w < n)
      do lo = 1, n, 2_I4P * w
         mid = min(lo + w, n + 1_I4P) ; hi = min(lo + 2_I4P * w, n + 1_I4P)
         i = lo ; j = mid ; i1 = mid ; j1 = hi
         do e = lo, hi - 1_I4P
            if (i < i1 .and. j < j1) then
               if (.not. less6(key(:, order(j)), key(:, order(i)))) then
                  tmp(e) = order(i) ; i = i + 1_I4P
               else
                  tmp(e) = order(j) ; j = j + 1_I4P
               endif
            elseif (i < i1) then
               tmp(e) = order(i) ; i = i + 1_I4P
            else
               tmp(e) = order(j) ; j = j + 1_I4P
            endif
         enddo
      enddo
      order = tmp ; w = 2_I4P * w
   enddo
   ok = .false.
   do k = 2, n
      if (all(key(:, order(k)) == key(:, order(k - 1)))) return
   enddo
   do k = 1, n
      rv(1:3) = key(4:6, order(k)) ; rv(4:6) = key(1:3, order(k))
      lo = 1_I4P ; hi = n ! the first key not less than rv
      do while (lo < hi)
         mid = (lo + hi) / 2_I4P
         if (less6(key(:, order(mid)), rv)) then
            lo = mid + 1_I4P
         else
            hi = mid
         endif
      enddo
      if (.not. all(key(:, order(lo)) == rv)) return
   enddo
   ok = .true.
   endsubroutine surface_closed

   pure logical function less6(p, q)
   !< Lexicographic order of two 6-tuples.
   real(R8P), intent(in) :: p(6), q(6) !< Keys.
   integer(I4P)          :: k          !< Work.

   less6 = .false.
   do k = 1, 6
      if (p(k) < q(k)) then
         less6 = .true. ; return
      elseif (p(k) > q(k)) then
         return
      endif
   enddo
   endfunction less6

   subroutine shells_nested(srf, nt, tri, ok)
   !< Every shell of the surface has a signed volume of the sign of its nesting depth (positive at even depth: outward
   !< normals around a solid; negative at odd depth: a cavity). The depth of a shell is the number of other shells whose
   !< winding number at its first vertex exceeds 1/2 in magnitude (the shells neither cross nor touch: checked before).
   type(surface_stl_object), intent(in)  :: srf           !< Surface (shells from its analysis).
   integer(I4P),             intent(in)  :: nt            !< Triangles.
   real(R8P),                intent(in)  :: tri(3, 3, nt) !< Triangles (the surface's facets, in order).
   logical,                  intent(out) :: ok            !< Consistent.
   real(R8P),    allocatable :: vol(:), wn(:)             !< Shell volumes, winding numbers.
   integer(I4P), allocatable :: sh(:), first(:), depth(:) !< Shell of each facet, its first facet, depths.
   real(R8P)                 :: r(3), p(3), four_pi       !< Work.
   integer(I4P)              :: ns, t, a, b               !< Work.

   four_pi = 4._R8P * acos(-1._R8P)
   ns = max(1_I4P, srf%get_shells_number())
   allocate(sh(nt), first(ns), vol(ns), depth(ns), wn(ns))
   first = 0_I4P ; vol = 0._R8P ; depth = 0_I4P
   do t = 1, nt
      sh(t) = max(1_I4P, min(ns, srf%facet_shell(t)))
      if (first(sh(t)) == 0_I4P) first(sh(t)) = t
   enddo
   do t = 1, nt ! signed volumes about each shell's first vertex
      r = tri(:, 1, first(sh(t)))
      vol(sh(t)) = vol(sh(t)) + dot_product(tri(:, 1, t) - r, cross_3(tri(:, 2, t) - r, tri(:, 3, t) - r)) / 6._R8P
   enddo
   do a = 1, ns
      if (first(a) == 0_I4P) cycle
      p = tri(:, 1, first(a)) ; wn = 0._R8P
      do t = 1, nt
         b = sh(t)
         if (b == a) cycle
         wn(b) = wn(b) + solid_angle_oac(p(1), p(2), p(3), tri(1, 1, t), tri(2, 1, t), tri(3, 1, t), tri(1, 2, t), &
                                         tri(2, 2, t), tri(3, 2, t), tri(1, 3, t), tri(2, 3, t), tri(3, 3, t))
      enddo
      depth(a) = count(abs(wn) / four_pi > 0.5_R8P, kind=I4P)
   enddo
   ok = .true.
   do a = 1, ns
      if (first(a) == 0_I4P) cycle
      if (mod(depth(a), 2_I4P) == 0_I4P) then
         ok = ok .and. vol(a) > 0._R8P
      else
         ok = ok .and. vol(a) < 0._R8P
      endif
   enddo
   endsubroutine shells_nested

   subroutine smallest_faces_3d(patch, walls)
   !< The smallest faces from the labelled wall terms of the cut: those of non-zero vector area, by piece (CSR, each
   !< piece's in term order: the order in which the cut sums them), then those of the cells without trusted pieces, per
   !< cell in cell order.
   type(overset_patch_3d_t), intent(inout) :: patch !< Patch (cut, face list, triangle labels set).
   type(cut_walls_3d_t),     intent(in)    :: walls !< Wall terms of the cut.
   integer(I4P), allocatable :: tcell(:), pos(:)    !< Cell of each term, fill positions.
   logical,      allocatable :: keep(:)             !< Terms kept.
   integer(I4P) :: np, nx, ny, c, t, p, ns, nu, k, f, last !< Work.

   np = patch%cut%np ; nx = patch%cut%nx ; ny = patch%cut%ny
   allocate(tcell(walls%n), keep(walls%n), patch%piece_sf_start(np + 1), pos(np + 1))
   do c = 1, size(walls%cell_start, kind=I4P) - 1_I4P
      do t = walls%cell_start(c), walls%cell_start(c + 1) - 1_I4P
         tcell(t) = c
      enddo
   enddo
   patch%piece_sf_start = 0_I4P ; ns = 0_I4P ; nu = 0_I4P ; last = 0_I4P
   do t = 1, walls%n
      keep(t) = any(walls%area(:, t) /= 0._R8P)
      if (.not. keep(t)) cycle
      ns = ns + 1_I4P
      if (walls%piece(t) > 0_I4P) then
         patch%piece_sf_start(walls%piece(t) + 1) = patch%piece_sf_start(walls%piece(t) + 1) + 1_I4P
      elseif (tcell(t) /= last) then
         nu = nu + 1_I4P ; last = tcell(t)
      endif
   enddo
   patch%piece_sf_start(1) = 1_I4P
   do p = 1, np
      patch%piece_sf_start(p + 1) = patch%piece_sf_start(p + 1) + patch%piece_sf_start(p)
   enddo
   patch%ns = ns ; patch%nuc = nu
   allocate(patch%sf_piece(ns), patch%sf_cell(3, ns), patch%sf_face(ns), patch%sf_facet(ns), patch%sf_minor(5, ns), &
            patch%sf_tri(ns), patch%sf_area(3, ns), patch%sf_centroid(3, ns), patch%sf_wet(ns), patch%uc_cell(3, nu), &
            patch%uc_sf_start(nu + 1))
   pos = patch%piece_sf_start ; nu = 0_I4P ; last = 0_I4P
   do t = 1, walls%n
      if (.not. keep(t)) cycle
      c = tcell(t)
      if (walls%piece(t) > 0_I4P) then
         k = pos(walls%piece(t)) ; pos(walls%piece(t)) = k + 1_I4P
      else
         k = pos(np + 1) ; pos(np + 1) = k + 1_I4P
         if (c /= last) then
            nu = nu + 1_I4P ; last = c ; patch%uc_sf_start(nu) = k
            patch%uc_cell(:, nu) = [mod(c - 1_I4P, nx) + 1_I4P, mod((c - 1_I4P) / nx, ny) + 1_I4P, (c - 1_I4P) / (nx * ny) + 1_I4P]
         endif
      endif
      f = patch%tri_face(walls%facet(t))
      patch%sf_piece(k) = walls%piece(t) ; patch%sf_facet(k) = walls%facet(t) ; patch%sf_face(k) = f
      patch%sf_cell(:, k) = [mod(c - 1_I4P, nx) + 1_I4P, mod((c - 1_I4P) / nx, ny) + 1_I4P, (c - 1_I4P) / (nx * ny) + 1_I4P]
      patch%sf_minor(:, k) = patch%face_id(:, f) ; patch%sf_tri(k) = patch%tri_id(walls%facet(t))
      patch%sf_area(:, k) = -walls%area(:, t) ; patch%sf_centroid(:, k) = walls%centroid(:, t)
      patch%sf_wet(k) = walls%wet_area(t)
   enddo
   patch%uc_sf_start(nu + 1) = ns + 1_I4P
   endsubroutine smallest_faces_3d

   subroutine faces_csr_3d(patch)
   !< The smallest faces by minor face, in smallest-face order, and the covered vectors summed in that order.
   type(overset_patch_3d_t), intent(inout) :: patch    !< Patch (smallest faces set).
   integer(I4P), allocatable :: pos(:)                 !< Fill positions.
   integer(I4P) :: f, s, k                             !< Work.

   allocate(patch%face_sf_start(patch%nf + 1), patch%face_sf(patch%ns), pos(patch%nf))
   patch%face_sf_start = 0_I4P
   do s = 1, patch%ns
      f = patch%sf_face(s) ; patch%face_sf_start(f + 1) = patch%face_sf_start(f + 1) + 1_I4P
   enddo
   patch%face_sf_start(1) = 1_I4P
   do f = 1, patch%nf
      patch%face_sf_start(f + 1) = patch%face_sf_start(f + 1) + patch%face_sf_start(f)
   enddo
   pos = patch%face_sf_start(1:patch%nf)
   do s = 1, patch%ns
      f = patch%sf_face(s) ; patch%face_sf(pos(f)) = s ; pos(f) = pos(f) + 1_I4P
   enddo
   do f = 1, patch%nf
      patch%face_covered(:, f) = 0._R8P
      do k = patch%face_sf_start(f), patch%face_sf_start(f + 1) - 1_I4P
         patch%face_covered(:, f) = patch%face_covered(:, f) - patch%sf_area(:, patch%face_sf(k))
      enddo
   enddo
   endsubroutine faces_csr_3d

   subroutine overlap_weights_3d(block, cells, x, y, z, w, status)
   !< The intersection weights (Picard et al. 2021) of the receptor cells `cells(:, r)` = (i, j, k) of a curvilinear
   !< block against the cells of the rectilinear donor grid of cell edges x(nx+1), y(ny+1), z(nz+1) (issue #106 step 3,
   !< the 3-D twin of `overlap_weights_2d`).
   !<
   !< Each receptor is the closed polyhedron of its 6 faces, each split along its shorter diagonal into 2 triangles (12
   !< in all), oriented out of the cell by the block's handedness (`face_offsets`). On a tie of the computed diagonal
   !< lengths the diagonal through the least global vertex index (i + ni ((j - 1) + nj (k - 1))) is taken, so the two
   !< receptors sharing a face split it alike and the receptors of a block tile the polyhedral region they cover (no gap,
   !< no overlap); exactly collapsed triangles (a pole) are dropped. It is cut against the donor sub-grid around its
   !< bounding box as a hole in a solid box (`cut_cells_3d`, the direct variant), in coordinates shifted to the
   !< sub-grid's low corner (the cut cells' face sections sum Green's integrals about the origin: the shift keeps their
   !< relative error at round-off), on every axis. Until issue #107 an axis was shifted only where the shift was exact
   !< for every vertex and sub-grid edge (else not shifted): an inexact shift could leave a facet within an ulp of a grid
   !< plane, which the cut cells double-counted (a receptor 1e-16 above z = 0, shifted by 1.3, lost 0.9 % of its
   !< volume; up to 21 % in `fossil_test_overlap_weights_3d`). The cut is right for such facets since #107, so the
   !< fallback is gone: without it the weights move by <= 1e-15 and receptor volumes by <= 1.3e-15 relative, and 17
   !< zero-volume round-off donors leave the lists. It bought no accuracy either: the shift was inexact only near the
   !< origin (coordinates of opposite sign or much smaller than the shift), where the unshifted sums are accurate too.
   !< V_c n V_c' is then the fluid of donor cell c': its volume and centroid are the cut's cell totals (`alpha` x the
   !< shifted cell volume, `centroid` + the shift), never the per-piece arrays (the cut may leave a donor cell untrusted,
   !< pieces = -1, whose totals are right and whose pieces are empty).
   !<
   !< Donors: every donor cell of the exact cut with fluid (a piece, or untrusted with alpha > 0), in (k, j, i) order of
   !< the donor grid; a receptor face on an interior grid plane can leave a donor of zero or round-off volume (listed,
   !< volume clipped at 0, as in 2-D). receptor_volume, covered, flags and status as in `overlap_weights_2d`; the
   !< receptor is DEGENERATE with fewer than 4 triangles, a polyhedron volume <= 0 in the block's handedness, or facets
   !< crossing (`find_self_intersections`, exact: a folded cell). Serial, CPU only.
   type(curvilinear_block_t), intent(in)            :: block      !< Receptor block.
   integer(I4P),              intent(in)            :: cells(:,:) !< Receptor cells (i, j, k) [3, nr].
   real(R8P),                 intent(in)            :: x(:), y(:), z(:) !< Donor grid edges, strictly increasing.
   type(overlap_weights_t),   intent(inout)         :: w          !< Weights (reset first).
   integer(I4P),              intent(out), optional :: status     !< OVERLAP_WEIGHTS_*.
   type(overlap_weights_t)              :: fresh  !< Default.
   type(receptor_buffer_t), allocatable :: buf(:) !< Per receptor intersections.
   real(R8P)    :: orient                         !< Block handedness.
   integer(I4P) :: st, nr, r, k, n                !< Work.

   w = fresh
   st = OVERLAP_WEIGHTS_OK
   if (size(cells, 1) /= 3 .or. size(x) < 2 .or. size(y) < 2 .or. size(z) < 2) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (.not. (all(abs(x) <= huge(1._R8P)) .and. all(abs(y) <= huge(1._R8P)) .and. all(abs(z) <= huge(1._R8P)))) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (any(x(2:) <= x(:size(x) - 1)) .or. any(y(2:) <= y(:size(y) - 1)) .or. any(z(2:) <= z(:size(z) - 1))) then
      st = OVERLAP_WEIGHTS_BAD_INPUT
   elseif (size(cells, 2) > 0) then
      if (any(cells(1, :) < 1_I4P) .or. any(cells(1, :) > block%ni - 1_I4P) .or. any(cells(2, :) < 1_I4P) .or. &
          any(cells(2, :) > block%nj - 1_I4P) .or. any(cells(3, :) < 1_I4P) .or. any(cells(3, :) > block%nk - 1_I4P)) &
         st = OVERLAP_WEIGHTS_BAD_INPUT
   endif
   if (st /= OVERLAP_WEIGHTS_OK) then
      call put_status(st, status) ; return
   endif
   nr = size(cells, 2, kind=I4P) ; w%nr = nr
   orient = 1._R8P ; if (nr > 0_I4P) orient = block%orientation()
   allocate(buf(nr), w%start(nr + 1), w%receptor_volume(nr), w%covered(nr), w%flag(nr))
   do r = 1, nr
      call receptor_weights_3d(block, cells(:, r), orient, x, y, z, buf(r), w%receptor_volume(r), w%covered(r), w%flag(r))
   enddo
   w%start(1) = 1_I4P
   do r = 1, nr
      w%start(r + 1) = w%start(r) + buf(r)%n
   enddo
   w%nw = w%start(nr + 1) - 1_I4P
   allocate(w%donor(3, w%nw), w%volume(w%nw), w%weight(w%nw), w%centroid(3, w%nw))
   do r = 1, nr
      n = buf(r)%n ; k = w%start(r)
      if (n == 0_I4P) cycle
      w%donor(:, k:k + n - 1) = buf(r)%donor(:, 1:n)
      w%volume(k:k + n - 1) = buf(r)%volume(1:n)
      w%centroid(:, k:k + n - 1) = buf(r)%centroid(:, 1:n)
      if (w%receptor_volume(r) > 0._R8P) then
         w%weight(k:k + n - 1) = buf(r)%volume(1:n) / w%receptor_volume(r)
      else
         w%weight(k:k + n - 1) = 0._R8P
      endif
   enddo
   if (nr > 0_I4P) then
      if (any(w%flag == OVERLAP_WEIGHTS_DEGENERATE)) then
         st = OVERLAP_WEIGHTS_DEGENERATE
      elseif (any(w%flag == OVERLAP_WEIGHTS_OPEN)) then
         st = OVERLAP_WEIGHTS_OPEN
      endif
   endif
   call put_status(st, status)
   endsubroutine overlap_weights_3d

   subroutine receptor_weights_3d(block, ijk, orient, x, y, z, b, rvol, cov, flag)
   !< The intersections of receptor cell ijk with the donor cells (see `overlap_weights_3d`).
   type(curvilinear_block_t), intent(in)    :: block      !< Receptor block.
   integer(I4P),              intent(in)    :: ijk(3)     !< Receptor cell.
   real(R8P),                 intent(in)    :: orient     !< Block handedness.
   real(R8P),                 intent(in)    :: x(:), y(:), z(:) !< Donor grid edges.
   type(receptor_buffer_t),   intent(inout) :: b          !< Intersections.
   real(R8P),                 intent(out)   :: rvol       !< Receptor volume (sum of the intersections).
   real(R8P),                 intent(out)   :: cov        !< Covered fraction.
   integer(I4P),              intent(out)   :: flag       !< OVERLAP_WEIGHTS_*.
   integer(I4P), parameter                :: QA(3) = [1_I4P, 2_I4P, 3_I4P], QB(3) = [1_I4P, 3_I4P, 4_I4P] !< Box quad split.
   type(surface_stl_object)               :: srf          !< Solid box with the receptor as a hole.
   type(cut_cells_3d_t)                   :: cut          !< Its cut of the sub-grid.
   type(facet_object),        allocatable :: fct(:)       !< Facets.
   type(intersection_pair_t), allocatable :: pairs(:)     !< Crossing facet pairs.
   real(R8P),                 allocatable :: xs(:), ys(:), zs(:) !< Shifted sub-grid edges.
   real(R8P)    :: v(3, 8), tri(3, 3, 24), c(3, 4), t2(3, 3, 2), lo(3), hi(3), o(3), blo(3), bhi(3), cv(3, 8) !< Work.
   real(R8P)    :: d, vol, comp, h(3)                     !< Work.
   integer(I4P) :: g(8), off(3, 4), gc(4), tid(2), q, fc, m, n2, nt, nx, ny, nz, i0, i1, j0, j1, k0, k1, ii, jj, kk, st !< Work.
   integer(I4P) :: nd, kf                                 !< Work.
   logical      :: inbox                                  !< The receptor lies in the closed donor box.

   b%n = 0_I4P ; rvol = 0._R8P ; cov = 0._R8P ; flag = OVERLAP_WEIGHTS_OK
   do q = 1, 8
      ii = ijk(1) + mod(q - 1_I4P, 2_I4P) ; jj = ijk(2) + mod((q - 1_I4P) / 2_I4P, 2_I4P) ; kk = ijk(3) + (q - 1_I4P) / 4_I4P
      v(:, q) = [block%x(ii, jj, kk), block%y(ii, jj, kk), block%z(ii, jj, kk)]
      g(q) = ii + block%ni * ((jj - 1_I4P) + block%nj * (kk - 1_I4P))
   enddo
   nt = 0_I4P
   do fc = 1, 6
      call face_offsets(fc, orient, off)
      do m = 1, 4
         q = 1_I4P + off(1, m) + 2_I4P * off(2, m) + 4_I4P * off(3, m) ; c(:, m) = v(:, q) ; gc(m) = g(q)
      enddo
      call split_quad(c, gc, t2, tid, n2)
      do m = 1, n2
         nt = nt + 1_I4P ; tri(:, :, nt) = t2(:, :, m)
      enddo
   enddo
   if (nt < 4_I4P) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   lo = minval(v, dim=2) ; hi = maxval(v, dim=2)
   inbox = lo(1) >= x(1) .and. hi(1) <= x(nx + 1) .and. lo(2) >= y(1) .and. hi(2) <= y(ny + 1) .and. &
           lo(3) >= z(1) .and. hi(3) <= z(nz + 1)
   i0 = edge_below(x, lo(1)) ; i1 = edge_above(x, hi(1))
   j0 = edge_below(y, lo(2)) ; j1 = edge_above(y, hi(2))
   k0 = edge_below(z, lo(3)) ; k1 = edge_above(z, hi(3))
   ! the receptor, shifted to the sub-grid's low corner: its volume, then reversed as a hole in a solid box around it
   o(1) = x(i0) ; o(2) = y(j0) ; o(3) = z(k0)
   do m = 1, nt
      do q = 1, 3
         tri(:, q, m) = tri(:, q, m) - o
      enddo
   enddo
   vol = 0._R8P
   do m = 1, nt
      vol = vol + dot_product(tri(:, 1, m), cross_3(tri(:, 2, m), tri(:, 3, m))) / 6._R8P
   enddo
   if (.not. (vol > 0._R8P)) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   xs = x(i0:i1) - o(1) ; ys = y(j0:j1) - o(2) ; zs = z(k0:k1) - o(3)
   blo = [min(lo(1) - o(1), xs(1)), min(lo(2) - o(2), ys(1)), min(lo(3) - o(3), zs(1))]
   bhi = [max(hi(1) - o(1), xs(size(xs))), max(hi(2) - o(2), ys(size(ys))), max(hi(3) - o(3), zs(size(zs)))]
   d = maxval(bhi - blo) ; blo = blo - d ; bhi = bhi + d
   do q = 1, 8
      cv(:, q) = merge(bhi, blo, [btest(q - 1, 0), btest(q - 1, 1), btest(q - 1, 2)])
   enddo
   allocate(fct(nt + 12))
   do m = 1, nt
      fct(m)%vertex(1) = vector_R8P(tri(1, 1, m), tri(2, 1, m), tri(3, 1, m))
      fct(m)%vertex(2) = vector_R8P(tri(1, 3, m), tri(2, 3, m), tri(3, 3, m))
      fct(m)%vertex(3) = vector_R8P(tri(1, 2, m), tri(2, 2, m), tri(3, 2, m))
   enddo
   do fc = 1, 6 ! the box, outward
      call face_offsets(fc, 1._R8P, off)
      do m = 1, 4
         c(:, m) = cv(:, 1_I4P + off(1, m) + 2_I4P * off(2, m) + 4_I4P * off(3, m))
      enddo
      nd = nt + 2_I4P * fc - 1_I4P
      do kf = 1, 3
         q = QA(kf) ; fct(nd)%vertex(kf) = vector_R8P(c(1, q), c(2, q), c(3, q))
         q = QB(kf) ; fct(nd + 1)%vertex(kf) = vector_R8P(c(1, q), c(2, q), c(3, q))
      enddo
   enddo
   call srf%adopt_facets(facets=fct)
   call srf%find_self_intersections(pairs)
   if (size(pairs) > 0) then
      flag = OVERLAP_WEIGHTS_DEGENERATE ; return
   endif
   if (i1 <= i0 .or. j1 <= j0 .or. k1 <= k0) then
      flag = OVERLAP_WEIGHTS_OUTSIDE ; return
   endif
   call cut_cells_3d(srf, xs, ys, zs, cut, status=st)
   if (st /= 0_I4P) then
      flag = OVERLAP_WEIGHTS_OPEN ; return
   endif
   if (.not. all(cut%closed)) then
      flag = OVERLAP_WEIGHTS_OPEN ; return
   endif
   m = count(cut%pieces > 0_I4P .or. (cut%pieces < 0_I4P .and. cut%alpha > 0._R8P), kind=I4P)
   allocate(b%donor(3, m), b%volume(m), b%centroid(3, m))
   comp = 0._R8P
   do kk = 1, cut%nz
      do jj = 1, cut%ny
         do ii = 1, cut%nx
            if (cut%pieces(ii, jj, kk) == 0_I4P) cycle
            if (cut%pieces(ii, jj, kk) < 0_I4P .and. .not. cut%alpha(ii, jj, kk) > 0._R8P) cycle
            h = [xs(ii + 1) - xs(ii), ys(jj + 1) - ys(jj), zs(kk + 1) - zs(kk)]
            b%n = b%n + 1_I4P
            b%donor(:, b%n) = [i0 + ii - 1_I4P, j0 + jj - 1_I4P, k0 + kk - 1_I4P]
            b%volume(b%n) = max(0._R8P, cut%alpha(ii, jj, kk) * (h(1) * h(2) * h(3)))
            b%centroid(:, b%n) = cut%centroid(:, ii, jj, kk) + o
            call neumaier(rvol, comp, b%volume(b%n))
         enddo
      enddo
   enddo
   rvol = rvol + comp
   if (.not. (rvol > 0._R8P)) then
      rvol = 0._R8P ; flag = OVERLAP_WEIGHTS_OUTSIDE ; b%n = 0_I4P ; return
   endif
   if (inbox) then
      cov = 1._R8P
   else
      cov = min(rvol / vol, nearest(1._R8P, -1._R8P))
   endif
   endsubroutine receptor_weights_3d
endmodule fossil_overset_patch
