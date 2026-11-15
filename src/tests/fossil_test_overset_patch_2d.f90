!< FOSSIL test: the patched interface of a conservative overset assembly in 2-D (issue #106 step 1).
program fossil_test_overset_patch_2d
!< FOSSIL test: the patched interface of a conservative overset assembly in 2-D (issue #106 step 1).
!<
!< `overset_patch_2d` on a cylinder O-grid in a Cartesian background at two sizes and on the step-0 edge cases (a
!< boundary along grid lines and through grid nodes, a diamond through grid nodes, a boundary leaving the major box, two
!< disjoint minor blocks), and its refusals (overlapping blocks, nested blocks of equal orientation, an open boundary).
!< Checks per case:
!< - the cut equals `cut_cells(boundaries=.true.)` of the same loops built independently, bit for bit, and the labels
!<   `piece_seg_src` flag exactly the wall segments, each lying on its source segment;
!< - every major piece closes: sum_f W_f S_f n_f + sum of its smallest faces = 0 to 1e-15 h;
!< - per piece, sum sf_area = - piece_wet_normal (bit for bit, or to 1e-15 h where the piece has several wall chains);
!< - a face inside the major box: face_covered = face_vector to 1e-14 |face|, and its smallest faces chain exactly
!<   from its first vertex to its last; a face leaving the box: |face_covered| < |face_vector|;
!< - the Kahan sum of the piece areas = box - minor region to 1e-15 relative;
!< - degenerate pieces (zero area) have a non-degenerate merge target in a neighbouring cell.
use fossil,                   only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, cut_cells_2d_t, &
                                     wall_cutter_2d, overset_patch_2d_t, overset_patch_2d, OVERSET_PATCH_OK, &
                                     OVERSET_PATCH_OPEN, OVERSET_PATCH_OVERLAP, OVERSET_PATCH_NESTING
use fossil_cut_cells_compare, only : cut_cells_2d_mismatches
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P) !< Pi.
type(curvilinear_block_2d_t)  :: blocks(2)  !< Minor blocks.
type(overset_patch_2d_t)      :: patch, p2  !< Patches.
type(curve_object)            :: wc         !< Wall cutter.
real(R8P), allocatable        :: x(:), pu(:), pv(:), qu(:), qv(:) !< Grid, polygons.
integer(I4P), allocatable     :: faces(:,:) !< Explicit faces.
integer(I4P)                  :: k, n, m, st !< Counters.
logical                       :: ok      !< Checks.

ok = .true.
! 1. cylinder O-grid (r in [1, 2], 32 m faces) in the background [-3, 3]^2 of 40 m cells
do m = 1, 2
   call ogrid(m, blocks(1))
   x = grid_1d(-3._R8P, 3._R8P, 40_I4P * m)
   n = 32_I4P * m ; call circle(0._R8P, 0._R8P, 2._R8P, n, pu, pv, ring=.true., blk=blocks(1))
   call check_case('cylinder O-grid, m = '//trim(itoa(m)), blocks(1:1), [1_I4P], x, pu, pv, [n], &
                   36._R8P - shoelace(pu, pv), .true., ok)
enddo
! the orientation is the inverse of the wall cutter's
call wall_cutter_2d(blocks(1:1), wc, wall_name='overlap', status=st)
call overset_patch_2d(blocks(1:1), [1_I4P], x, x, patch, status=st)
ok = ok .and. st == OVERSET_PATCH_OK .and. abs(patch%boundary%get_area() + wc%get_area()) <= 1e-14_R8P * abs(wc%get_area())
print '(A,2F12.7)', '  signed areas (patch, wall_cutter_2d): ', patch%boundary%get_area(), wc%get_area()
! the explicit face list of the same faces gives the same interface
allocate(faces(4, 64))
do k = 1, 64
   faces(:, k) = [1_I4P, k, 4_I4P, 4_I4P]
enddo
call overset_patch_2d(blocks(1:1), [1_I4P], x, x, p2, faces=faces, status=st)
ok = ok .and. st == OVERSET_PATCH_OK .and. p2%ns == patch%ns .and. all(p2%sf_minor == patch%sf_minor) .and. &
     all(p2%sf_area == patch%sf_area) .and. all(p2%face_covered == patch%face_covered)
print '(A,I0,A,L1)', '  explicit faces: status ', st, ', same interface ', p2%ns == patch%ns
x = grid_1d(-3._R8P, 3._R8P, 40_I4P) ! h = 0.15
! 2. a square on grid lines, a vertex at every node it passes: degenerate pieces along its sides
n = 0
deallocate(pu, pv) ; allocate(pu(80), pv(80))
do k = 11, 30
   n = n + 1 ; pu(n) = x(k) ; pv(n) = x(11)
enddo
do k = 11, 30
   n = n + 1 ; pu(n) = x(31) ; pv(n) = x(k)
enddo
do k = 31, 12, -1
   n = n + 1 ; pu(n) = x(k) ; pv(n) = x(31)
enddo
do k = 31, 12, -1
   n = n + 1 ; pu(n) = x(11) ; pv(n) = x(k)
enddo
call ring(pu, pv, blocks(1))
call check_case('square on grid lines, 80 vertices', blocks(1:1), [1_I4P], x, pu, pv, [80_I4P], 36._R8P - shoelace(pu, pv), &
                .true., ok, ndeg=39_I4P)
! 3. the same with 4 vertices (segments along grid lines through nodes)
pu = [x(11), x(31), x(31), x(11)] ; pv = [x(11), x(11), x(31), x(31)]
call ring(pu, pv, blocks(1))
call check_case('square on grid lines, 4 vertices', blocks(1:1), [1_I4P], x, pu, pv, [4_I4P], 36._R8P - shoelace(pu, pv), &
                .true., ok, ndeg=39_I4P)
! 4. a diamond with its vertices on grid nodes, its sides through grid nodes
pu = [x(21), x(31), x(21), x(11)] ; pv = [x(11), x(21), x(31), x(21)]
call ring(pu, pv, blocks(1))
call check_case('diamond through grid nodes', blocks(1:1), [1_I4P], x, pu, pv, [4_I4P], 36._R8P - shoelace(pu, pv), &
                .true., ok)
! 5. a boundary leaving the major box: circle r = 1.5 centred at (2.5, 0.4)
call circle(2.5_R8P, 0.4_R8P, 1.5_R8P, 96_I4P, pu, pv, ring=.true., blk=blocks(1))
call check_case('a boundary leaving the major box', blocks(1:1), [1_I4P], x, pu, pv, [96_I4P], &
                36._R8P - clipped_area(pu, pv, [-3._R8P, 3._R8P, 3._R8P, -3._R8P], [-3._R8P, -3._R8P, 3._R8P, 3._R8P]), &
                .false., ok)
! 6. two disjoint minor blocks
call circle(-1.6_R8P, 0.1_R8P, 1._R8P, 64_I4P, pu, pv, ring=.true., blk=blocks(1))
call circle(1.55_R8P, -0.2_R8P, 1._R8P, 64_I4P, qu, qv, ring=.true., blk=blocks(2))
call check_case('two disjoint minor blocks', blocks, [1_I4P, 2_I4P], x, [pu, qu], [pv, qv], [64_I4P, 64_I4P], &
                36._R8P - shoelace(pu, pv) - shoelace(qu, qv), .true., ok)
! refusals: 7. two overlapping minor blocks
call circle(-0.4_R8P, 0.3_R8P, 1._R8P, 64_I4P, qu, qv, ring=.true., blk=blocks(2))
call overset_patch_2d(blocks, [1_I4P, 2_I4P], x, x, patch, status=st)
print '(A,I0)', '  two overlapping minor blocks: status ', st
ok = ok .and. st == OVERSET_PATCH_OVERLAP .and. .not. allocated(patch%cut%alpha)
! 8. a minor block nested inside another's region with the same orientation (a minor-minor interface)
call circle(-1.6_R8P, 0.1_R8P, 0.4_R8P, 32_I4P, qu, qv, ring=.true., blk=blocks(2))
call overset_patch_2d(blocks, [1_I4P, 2_I4P], x, x, patch, status=st)
print '(A,I0)', '  nested minor blocks of equal orientation: status ', st
ok = ok .and. st == OVERSET_PATCH_NESTING
! 9. an open boundary: the record covers half the outer face
call circle(-1.6_R8P, 0.1_R8P, 1._R8P, 64_I4P, pu, pv, ring=.true., blk=blocks(1))
blocks(1)%rec(1)%hi(1) = 33
call overset_patch_2d(blocks(1:1), [1_I4P], x, x, patch, status=st)
print '(A,I0)', '  an open boundary: status ', st
ok = ok .and. st == OVERSET_PATCH_OPEN
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check_case(tag, blk, minor, x, pu, pv, lsz, expected, inside, ok, ndeg)
   !< overset_patch_2d on the blocks, with the checks of the program notes.
   character(*),                 intent(in)    :: tag          !< Case.
   type(curvilinear_block_2d_t), intent(in)    :: blk(:)       !< Blocks.
   integer(I4P),                 intent(in)    :: minor(:)     !< Minor blocks.
   real(R8P),                    intent(in)    :: x(:)         !< Grid edges (both directions).
   real(R8P),                    intent(in)    :: pu(:), pv(:) !< The expected loops (minor region on the left).
   integer(I4P),                 intent(in)    :: lsz(:)       !< Their sizes.
   real(R8P),                    intent(in)    :: expected     !< Expected fluid area.
   logical,                      intent(in)    :: inside       !< Every face lies inside the major box.
   logical,                      intent(inout) :: ok           !< Checks.
   integer(I4P),       optional, intent(in)    :: ndeg         !< Expected degenerate pieces.
   type(overset_patch_2d_t) :: pt                 !< Patch.
   type(curve_object)       :: c                  !< The loops, built here.
   type(cut_cells_2d_t)     :: ref, cpy           !< Reference cut, copy without boundaries.
   real(R8P)    :: h, r(2), sw(2), dmax, cmax, wmax, fmax, tot, comp, yk, tk, a(2), b(2), e(2), len, rel !< Work.
   integer(I4P) :: st, p, q, s, f, k, i, j, nbit, nrnd, nmis, nchain, nout, nlab, ndg, nbad, sv, t, nlong, nempty !< Work.
   logical      :: inn                            !< Face inside the box.

   call overset_patch_2d(blk, minor, x, x, pt, status=st)
   if (st /= OVERSET_PATCH_OK) then
      print '(3A,I0)', '  ', tag, ': status ', st ; ok = .false. ; return
   endif
   h = x(2) - x(1)
   ! the cut against the same loops built independently
   call c%adopt_loops(u=pu, v=pv, loop_sizes=lsz)
   call c%cut_cells(x, x, ref, boundaries=.true.)
   nmis = 0
   if (any(shape(ref%piece_seg) /= shape(pt%cut%piece_seg))) then
      nmis = nmis + 1
   else
      nmis = nmis + count(ref%piece_seg /= pt%cut%piece_seg) + count(ref%piece_seg_wall .neqv. pt%cut%piece_seg_wall) &
             + count(ref%piece_seg_src /= pt%cut%piece_seg_src) + count(ref%piece_seg_start /= pt%cut%piece_seg_start)
   endif
   cpy = pt%cut
   deallocate(cpy%piece_seg_start, cpy%piece_seg, cpy%piece_seg_wall, cpy%piece_seg_src)
   call c%cut_cells(x, x, ref)
   nmis = nmis + cut_cells_2d_mismatches(cpy, ref)
   ! labels: exactly the wall segments, on their source segment
   nlab = count(pt%cut%piece_seg_wall .neqv. pt%cut%piece_seg_src > 0)
   dmax = 0._R8P
   do q = 1, size(pt%cut%piece_seg_src, kind=I4P)
      s = pt%cut%piece_seg_src(q)
      if (s <= 0) cycle
      a = [pu(s), pv(s)] ; b = [pu(next_of(s, lsz)), pv(next_of(s, lsz))] ; e = b - a
      do k = 0, 1
         r = pt%cut%piece_seg(2 * k + 1:2 * k + 2, q) - a
         dmax = max(dmax, abs(r(1) * e(2) - r(2) * e(1)) / norm2(e) / (x(size(x)) - x(1)))
      enddo
   enddo
   ! per piece: closure, sum of the smallest faces against the wet normal
   cmax = 0._R8P ; wmax = 0._R8P ; nbit = 0 ; nrnd = 0 ; tot = 0._R8P ; comp = 0._R8P
   do p = 1, pt%cut%np
      i = pt%cut%piece_cell(1, p) ; j = pt%cut%piece_cell(2, p)
      sw = 0._R8P
      do k = pt%piece_sf_start(p), pt%piece_sf_start(p + 1) - 1
         sw = sw + pt%sf_area(:, k)
      enddo
      r = [(pt%cut%piece_aperture(2, p) - pt%cut%piece_aperture(1, p)) * (x(j + 1) - x(j)), &
           (pt%cut%piece_aperture(4, p) - pt%cut%piece_aperture(3, p)) * (x(i + 1) - x(i))] + sw
      cmax = max(cmax, norm2(r) / h)
      if (all(sw == -pt%cut%piece_wet_normal(:, p))) then
         nbit = nbit + 1
      else
         nrnd = nrnd + 1 ; wmax = max(wmax, norm2(sw + pt%cut%piece_wet_normal(:, p)) / h)
      endif
      yk = pt%cut%piece_area(p) - comp ; tk = tot + yk ; comp = (tk - tot) - yk ; tot = tk ! Kahan
   enddo
   ! per face: covered against the face, the chain of its smallest faces
   fmax = 0._R8P ; nchain = 0 ; nout = 0 ; nlong = 0
   do f = 1, pt%nf
      len = norm2(pt%face_vector(:, f))
      if (len == 0._R8P) cycle
      s = pt%face_sf_start(f) ; t = pt%face_sf_start(f + 1) - 1
      ! the face's oriented ends: the source segment of any of its smallest faces, or the face vector
      call face_ends_of(pt, f, a, b)
      inn = in_box(a, x) .and. in_box(b, x)
      if (inn) then
         fmax = max(fmax, norm2(pt%face_covered(:, f) - pt%face_vector(:, f)) / len)
         if (s > t) then
            nchain = nchain + 1 ; cycle
         endif
         if (any(pt%cut%piece_seg(3:4, pt%sf_seg(pt%face_sf(s))) /= a)) nchain = nchain + 1
         if (any(pt%cut%piece_seg(1:2, pt%sf_seg(pt%face_sf(t))) /= b)) nchain = nchain + 1
      else
         nout = nout + 1
         if (norm2(pt%face_covered(:, f)) > len * (1._R8P + 1e-14_R8P)) nlong = nlong + 1
      endif
      do k = s, t - 1
         sv = pt%sf_seg(pt%face_sf(k))
         if (any(pt%cut%piece_seg(1:2, sv) /= pt%cut%piece_seg(3:4, pt%sf_seg(pt%face_sf(k + 1))))) nchain = nchain + 1
      enddo
   enddo
   ! degenerate pieces
   ndg = count(pt%degenerate) ; nbad = 0 ; nempty = 0
   do p = 1, pt%cut%np
      if (.not. pt%degenerate(p)) cycle
      q = pt%merge_into(p)
      if (q <= 0) then ! an empty piece only: no aperture, no smallest face
         nempty = nempty + 1
         if (any(pt%cut%piece_aperture(:, p) > 0._R8P) .or. pt%piece_sf_start(p + 1) > pt%piece_sf_start(p)) &
            nbad = nbad + 1
         cycle
      endif
      if (pt%degenerate(q)) nbad = nbad + 1
      if (sum(abs(pt%cut%piece_cell(:, q) - pt%cut%piece_cell(:, p))) /= 1) nbad = nbad + 1
   enddo
   rel = abs(tot - expected) / expected
   print '(3A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)', '  ', tag, ': holes ', count(pt%hole), ', pieces ', pt%cut%np, &
         ', faces ', pt%nf, ', smallest faces ', pt%ns, ', degenerate ', ndg, ' (empty ', nempty, ', bad ', nbad, ')'
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)', '    cut vs independent mismatches ', nmis, ', label errors ', &
         nlab, ', off-line / L ', dmax, '; closure / h ', cmax, '; covered vs face ', fmax, '; area rel ', rel
   print '(A,I0,A,I0,A,ES9.2,A,I0,A,I0,A,I0,A)', '    sum sf = -wet_normal: bit ', nbit, ', round-off ', nrnd, &
         ' (max / h ', wmax, '); chain breaks ', nchain, '; faces leaving the box ', nout, ' (longer than the face ', &
         nlong, ')'
   ok = ok .and. nmis == 0 .and. nlab == 0 .and. dmax <= 1e-15_R8P .and. cmax <= 1e-15_R8P .and. wmax <= 1e-15_R8P .and. &
        fmax <= 1e-14_R8P .and. nchain == 0 .and. nlong == 0 .and. rel <= 1e-15_R8P .and. nbad == 0 .and. &
        all(pt%hole .eqv. pt%cut%pieces == 0)
   if (inside) ok = ok .and. nout == 0
   if (.not. inside) ok = ok .and. nout > 0
   if (present(ndeg)) ok = ok .and. ndg == ndeg
   endsubroutine check_case

   subroutine face_ends_of(pt, f, a, b)
   !< The oriented ends of face f (the minor region on the left), from its block.
   type(overset_patch_2d_t), intent(in)  :: pt   !< Patch.
   integer(I4P),             intent(in)  :: f    !< Face.
   real(R8P),                intent(out) :: a(2), b(2) !< Ends.
   integer(I4P) :: k !< Segment.

   do k = 1, size(pt%seg_face, kind=I4P) ! the curve segment of the face
      if (pt%seg_face(k) == f) exit
   enddo
   a = pt%boundary%get_vertex(k) ; b = pt%boundary%get_vertex(pt%boundary%get_next(k))
   endsubroutine face_ends_of

   pure logical function in_box(p, x)
   !< The point lies strictly inside the box of the grid.
   real(R8P), intent(in) :: p(2), x(:) !< Point, grid.

   in_box = p(1) > x(1) .and. p(1) < x(size(x)) .and. p(2) > x(1) .and. p(2) < x(size(x))
   endfunction in_box

   pure integer(I4P) function next_of(s, lsz)
   !< The next vertex of vertex s in loops of sizes lsz (closed).
   integer(I4P), intent(in) :: s, lsz(:) !< Vertex, loop sizes.
   integer(I4P)             :: l, k0     !< Work.

   k0 = 0
   do l = 1, size(lsz, kind=I4P)
      if (s <= k0 + lsz(l)) exit
      k0 = k0 + lsz(l)
   enddo
   next_of = k0 + mod(s - k0, lsz(l)) + 1_I4P
   endfunction next_of

   subroutine ogrid(m, blk)
   !< The cylinder O-grid: r in [1, 2], 32 m x 4 cells, its outer face (j = nj) named 'overlap'.
   integer(I4P),                 intent(in)    :: m   !< Refinement.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   type(topo_record_t)    :: rec   !< Record.
   real(R8P), allocatable :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)           :: ni, nj, i, j !< Counters.
   real(R8P)              :: r     !< Radius.

   ni = 32 * m + 1 ; nj = 5
   allocate(u(ni, nj), v(ni, nj))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = r * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni, :) = u(1, :) ; v(ni, :) = v(1, :)
   call blk%destroy
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine ogrid

   subroutine circle(cx, cy, rr, n, pu, pv, ring, blk)
   !< The regular n-gon of radius rr about (cx, cy), counter-clockwise from angle 0 (the O-grid's outer ring when rr = 2
   !< about the origin); with ring, the minor block of that boundary (`ring`).
   real(R8P),                    intent(in)    :: cx, cy, rr !< Centre, radius.
   integer(I4P),                 intent(in)    :: n          !< Vertices.
   real(R8P), allocatable,       intent(inout) :: pu(:), pv(:) !< Polygon.
   logical,                      intent(in)    :: ring       !< Build the block too.
   type(curvilinear_block_2d_t), intent(inout) :: blk        !< Block.
   integer(I4P) :: s !< Counter.

   if (allocated(pu)) deallocate(pu, pv)
   allocate(pu(n), pv(n))
   do s = 1, n
      pu(s) = cx + rr * cos(2 * PI * (s - 1) / n) ; pv(s) = cy + rr * sin(2 * PI * (s - 1) / n)
   enddo
   if (ring) then
      if (cx == 0._R8P .and. cy == 0._R8P .and. rr == 2._R8P) return ! the O-grid: same vertices, block already built
      call ring_block(pu, pv, blk)
   endif
   endsubroutine circle

   subroutine ring(pu, pv, blk)
   !< The minor block of the counter-clockwise polygon (pu, pv) (`ring_block`).
   real(R8P),                    intent(in)    :: pu(:), pv(:) !< Polygon.
   type(curvilinear_block_2d_t), intent(inout) :: blk          !< Block.

   call ring_block(pu, pv, blk)
   endsubroutine ring

   subroutine ring_block(pu, pv, blk)
   !< A one-layer block between the polygon (outer, j = 2, named 'overlap') and its copy shrunk by half towards its
   !< vertex mean (j = 1): the polygon is the block's outer boundary.
   real(R8P),                    intent(in)    :: pu(:), pv(:) !< Polygon.
   type(curvilinear_block_2d_t), intent(inout) :: blk          !< Block.
   type(topo_record_t)    :: rec           !< Record.
   real(R8P), allocatable :: u(:,:), v(:,:) !< Vertices.
   real(R8P)              :: c(2)          !< Vertex mean.
   integer(I4P)           :: n             !< Vertices.

   n = size(pu, kind=I4P)
   allocate(u(n + 1, 2), v(n + 1, 2))
   c = [sum(pu), sum(pv)] / n
   u(1:n, 2) = pu ; v(1:n, 2) = pv ; u(n + 1, 2) = pu(1) ; v(n + 1, 2) = pv(1)
   u(:, 1) = c(1) + 0.5_R8P * (u(:, 2) - c(1)) ; v(:, 1) = c(2) + 0.5_R8P * (v(:, 2) - c(2))
   call blk%destroy
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, 2, 1] ; rec%hi = [n + 1, 2, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine ring_block

   pure function grid_1d(a, b, n) result(x)
   !< n + 1 equally spaced edges on [a, b].
   real(R8P),    intent(in) :: a, b     !< Ends.
   integer(I4P), intent(in) :: n        !< Cells.
   real(R8P)                :: x(n + 1) !< Edges.
   integer(I4P)             :: i        !< Counter.

   do i = 1, n + 1
      x(i) = a + (b - a) * (i - 1) / n
   enddo
   x(n + 1) = b
   endfunction grid_1d

   pure real(R8P) function shoelace(px, py) result(a)
   !< Signed area of a polygon.
   real(R8P), intent(in) :: px(:), py(:) !< Vertices.
   integer(I4P)          :: k, n         !< Counters.

   n = size(px, kind=I4P) ; a = 0._R8P
   do k = 1, n
      a = a + px(k) * py(mod(k, n) + 1) - px(mod(k, n) + 1) * py(k)
   enddo
   a = 0.5_R8P * a
   endfunction shoelace

   real(R8P) function clipped_area(px, py, cx, cy) result(a)
   !< Area of polygon p clipped by the convex counter-clockwise polygon c (Sutherland-Hodgman).
   real(R8P), intent(in)  :: px(:), py(:), cx(:), cy(:) !< Polygons.
   real(R8P), allocatable :: wx(:), wy(:), qx(:), qy(:) !< Work.
   real(R8P)              :: s1, s2, t                  !< Work.
   integer(I4P)           :: n, k, k2, nc, l, nq        !< Counters.

   n = size(px, kind=I4P) ; nc = size(cx, kind=I4P)
   allocate(wx(2 * (n + nc) + 8), wy(2 * (n + nc) + 8), qx(2 * (n + nc) + 8), qy(2 * (n + nc) + 8))
   wx(1:n) = px ; wy(1:n) = py
   do l = 1, nc
      nq = 0
      do k = 1, n
         k2 = mod(k, n) + 1
         s1 = (cx(mod(l, nc) + 1) - cx(l)) * (wy(k) - cy(l)) - (cy(mod(l, nc) + 1) - cy(l)) * (wx(k) - cx(l))
         s2 = (cx(mod(l, nc) + 1) - cx(l)) * (wy(k2) - cy(l)) - (cy(mod(l, nc) + 1) - cy(l)) * (wx(k2) - cx(l))
         if (s1 >= 0._R8P) then
            nq = nq + 1 ; qx(nq) = wx(k) ; qy(nq) = wy(k)
         endif
         if ((s1 >= 0._R8P .and. s2 < 0._R8P) .or. (s1 < 0._R8P .and. s2 >= 0._R8P)) then
            t = s1 / (s1 - s2)
            nq = nq + 1 ; qx(nq) = wx(k) + t * (wx(k2) - wx(k)) ; qy(nq) = wy(k) + t * (wy(k2) - wy(k))
         endif
      enddo
      n = nq ; wx(1:n) = qx(1:n) ; wy(1:n) = qy(1:n)
      if (n == 0) exit
   enddo
   a = 0._R8P ; if (n > 2) a = shoelace(wx(1:n), wy(1:n))
   endfunction clipped_area

   pure function itoa(i) result(s)
   !< An integer as a string.
   integer(I4P), intent(in) :: i !< Integer.
   character(16)            :: s !< String.

   write(s, '(I0)') i
   endfunction itoa
endprogram fossil_test_overset_patch_2d
