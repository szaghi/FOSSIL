!< FOSSIL, the patched interface of a conservative overset assembly in 2-D (issue #106, Wang 1995).
module fossil_overset_patch
!< FOSSIL, the patched interface of a conservative overset assembly in 2-D (issue #106, Wang 1995).
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
use penf,                     only : I4P, R8P
use fossil_curvilinear_block, only : curvilinear_block_2d_t, TOPO_BOUNDARY
use fossil_curve_object,      only : curve_object, cut_cells_2d_t, CURVE_VALID
use fossil_cut_graph,         only : cut_links_t, cut_links_2d

implicit none
private
public :: overset_patch_2d_t, overset_patch_2d
public :: OVERSET_PATCH_OK, OVERSET_PATCH_BAD_INPUT, OVERSET_PATCH_NO_BOUNDARY, OVERSET_PATCH_OPEN, &
          OVERSET_PATCH_OVERLAP, OVERSET_PATCH_NESTING

integer(I4P), parameter :: OVERSET_PATCH_OK          = 0_I4P !< Done.
integer(I4P), parameter :: OVERSET_PATCH_BAD_INPUT   = 1_I4P !< A minor block or face out of range, or bad grid edges.
integer(I4P), parameter :: OVERSET_PATCH_NO_BOUNDARY = 2_I4P !< No boundary face of non-zero length.
integer(I4P), parameter :: OVERSET_PATCH_OPEN        = 3_I4P !< The faces do not chain into closed loops.
integer(I4P), parameter :: OVERSET_PATCH_OVERLAP     = 4_I4P !< Two loops (or one with itself) cross or touch.
integer(I4P), parameter :: OVERSET_PATCH_NESTING     = 5_I4P !< A loop's orientation does not match its nesting (or is
                                                             !< degenerate: fewer than 3 vertices, zero area).

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
endmodule fossil_overset_patch
