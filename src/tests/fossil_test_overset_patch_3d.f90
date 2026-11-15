!< FOSSIL test: the patched interface of a conservative overset assembly in 3-D (issue #106 step 3).
program fossil_test_overset_patch_3d
!< FOSSIL test: the patched interface of a conservative overset assembly in 3-D (issue #106 step 3).
!<
!< `overset_patch_3d` on the outer face of a sphere shell (poles collapsed, a left-handed block) in a Cartesian
!< background, on a box block whose faces lie on grid planes (degenerate pieces), on a sphere leaving the major box and
!< on two disjoint shells; its refusals (overlapping shells, nested shells of equal orientation, an open boundary).
!< Checks per case:
!< - the cut equals `cut_cells_3d` of the same surface with links and without `walls`, bit for bit (the labels change
!<   nothing), and every boundary facet points out of the minor region (its signed volume = the polyhedron's);
!< - every major piece closes: sum_f W_f S_f n_f + sum of its smallest faces = 0 to 1e-12 h^2;
!< - per piece, sum sf_area = - piece_wet_normal: bit for bit, or to round-off (counted, max printed);
!< - every smallest face lies in its major cell's box and its label's facet: its face's triangles sum to the face
!<   vector; a face inside the box: face_covered = face_vector to 1e-13 |face|; a face leaving it: |covered| < |face|;
!< - the Kahan sum of the piece volumes = box - minor region to 1e-13 relative;
!< - degenerate pieces (zero volume) have a non-degenerate merge target in a neighbouring cell.
use fossil,                   only : curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, cut_cells_3d_t, cut_cells_3d, &
                                     cut_links_t, overset_patch_3d_t, overset_patch_3d, OVERSET_PATCH_OK, &
                                     OVERSET_PATCH_OPEN, OVERSET_PATCH_OVERLAP, OVERSET_PATCH_NESTING, facet_object
use fossil_cut_cells_compare, only : cut_cells_3d_mismatches
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P) !< Pi.
type(curvilinear_block_t) :: blocks(2)  !< Minor blocks.
type(overset_patch_3d_t)  :: patch, p2  !< Patches.
real(R8P), allocatable    :: x(:)       !< Major grid.
integer(I4P), allocatable :: faces(:,:) !< Explicit faces.
real(R8P)                 :: vs, a      !< Volumes, box side.
integer(I4P)              :: k, st, n   !< Counters.
logical                   :: ok         !< Checks.

ok = .true.
! 1. the outer face (r = 2) of a sphere shell, 16 x 8 faces, in [-2.6, 2.6]^3 of 14^3 cells
x = grid_1d(-2.6_R8P, 2.6_R8P, 14_I4P)
call shell(16_I4P, [0._R8P, 0._R8P, 0._R8P], 1.5_R8P, 2._R8P, blocks(1), vs)
call check_case('sphere shell, 16 x 8 faces', blocks(1:1), [1_I4P], x, 5.2_R8P ** 3 - vs, .true., ok)
! the explicit face list of the same faces gives the same interface
n = 16 * 8 ; allocate(faces(5, n)) ; n = 0
do k = 1, 8
   do st = 1, 16
      n = n + 1 ; faces(:, n) = [1_I4P, st, k, 1_I4P, 6_I4P]
   enddo
enddo
call overset_patch_3d(blocks(1:1), [1_I4P], x, x, x, patch, status=st)
call overset_patch_3d(blocks(1:1), [1_I4P], x, x, x, p2, faces=faces, status=k)
ok = ok .and. st == OVERSET_PATCH_OK .and. k == OVERSET_PATCH_OK .and. p2%ns == patch%ns .and. &
     all(p2%sf_minor == patch%sf_minor) .and. all(p2%sf_area == patch%sf_area) .and. all(p2%face_covered == patch%face_covered)
print '(A,I0,A,L1)', '  explicit faces: status ', k, ', same interface ', p2%ns == patch%ns
! 2. a box block on grid planes (8 cells a side between x(4) and x(12)): degenerate pieces on its high sides
x = grid_1d(-3._R8P, 3._R8P, 15_I4P)
call box_block(x(4), x(12), 4_I4P, blocks(1))
a = x(12) - x(4)
call check_case('box on grid planes', blocks(1:1), [1_I4P], x, 216._R8P - a ** 3, .true., ok, degenerate=.true.)
! 3. a sphere leaving the major box: centre (2.4, 0.3, -0.2), r = 1.2
x = grid_1d(-3._R8P, 3._R8P, 14_I4P)
call shell(16_I4P, [2.4_R8P, 0.3_R8P, -0.2_R8P], 0.9_R8P, 1.2_R8P, blocks(1), vs)
call check_case('a sphere leaving the major box', blocks(1:1), [1_I4P], x, -1._R8P, .false., ok)
! 4. two disjoint shells
call shell(12_I4P, [-1.5_R8P, 0.2_R8P, 0.1_R8P], 0.7_R8P, 1._R8P, blocks(1), vs)
call shell(12_I4P, [1.4_R8P, -0.3_R8P, 0.2_R8P], 0.7_R8P, 1._R8P, blocks(2), a)
call check_case('two disjoint shells', blocks, [1_I4P, 2_I4P], x, 216._R8P - vs - a, .true., ok)
! refusals: 5. two overlapping shells
call shell(12_I4P, [-0.6_R8P, 0.2_R8P, 0.1_R8P], 0.7_R8P, 1._R8P, blocks(2), a)
call overset_patch_3d(blocks, [1_I4P, 2_I4P], x, x, x, patch, status=st)
print '(A,I0)', '  two overlapping shells: status ', st
ok = ok .and. st == OVERSET_PATCH_OVERLAP .and. .not. allocated(patch%cut%alpha)
! 6. a shell nested inside the other's region with the same orientation (a minor-minor interface)
call shell(12_I4P, [-1.5_R8P, 0.2_R8P, 0.1_R8P], 0.3_R8P, 0.45_R8P, blocks(2), a)
call overset_patch_3d(blocks, [1_I4P, 2_I4P], x, x, x, patch, status=st)
print '(A,I0)', '  nested shells of equal orientation: status ', st
ok = ok .and. st == OVERSET_PATCH_NESTING
! 7. an open boundary: the record covers half the outer face
blocks(1)%rec(1)%hi(1) = 7
call overset_patch_3d(blocks(1:1), [1_I4P], x, x, x, patch, status=st)
print '(A,I0)', '  an open boundary: status ', st
ok = ok .and. st == OVERSET_PATCH_OPEN
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check_case(tag, blk, minor, x, expected, inside, ok, degenerate)
   !< overset_patch_3d on the blocks, with the checks of the program notes (expected < 0: no volume check).
   character(*),              intent(in)    :: tag        !< Case.
   type(curvilinear_block_t), intent(in)    :: blk(:)     !< Blocks.
   integer(I4P),              intent(in)    :: minor(:)   !< Minor blocks.
   real(R8P),                 intent(in)    :: x(:)       !< Grid edges (all directions).
   real(R8P),                 intent(in)    :: expected   !< Expected fluid volume (< 0: not checked).
   logical,                   intent(in)    :: inside     !< Every face lies inside the major box.
   logical,                   intent(inout) :: ok         !< Checks.
   logical,         optional, intent(in)    :: degenerate !< Degenerate pieces expected.
   type(overset_patch_3d_t)  :: pt                 !< Patch.
   type(cut_cells_3d_t)      :: ref                !< Reference cut.
   type(cut_links_t)         :: links              !< Its links.
   type(facet_object), pointer :: fp               !< Facet.
   real(R8P),    allocatable :: fsum(:,:)          !< Triangles per face.
   real(R8P)    :: h, r(3), sw(3), cmax, wmax, fmax, tot, comp, yk, tk, len, rel, pf(6), lo(3), hi(3), vb, t(3, 3) !< Work.
   integer(I4P) :: st, p, k, f, i, j, l, nbit, nrnd, nmis, nout, nlong, ndg, nbad, nbox, nin, nfv !< Work.

   call overset_patch_3d(blk, minor, x, x, x, pt, status=st)
   if (st /= OVERSET_PATCH_OK) then
      print '(3A,I0)', '  ', tag, ': status ', st ; ok = .false. ; return
   endif
   h = x(2) - x(1)
   ! the cut against cut_cells_3d of the same surface with links, without walls
   call cut_cells_3d(pt%boundary, x, x, x, ref, links=links)
   nmis = cut_cells_3d_mismatches(pt%cut, ref)
   ! the facets: outward of the minor region (their signed volume), their sum per face = the face vector
   allocate(fsum(3, pt%nf)) ; fsum = 0._R8P ; vb = 0._R8P
   do k = 1, pt%boundary%get_facets_number()
      fp => pt%boundary%facet_at(k)
      do l = 1, 3
         t(:, l) = [fp%vertex(l)%x, fp%vertex(l)%y, fp%vertex(l)%z]
      enddo
      fsum(:, pt%tri_face(k)) = fsum(:, pt%tri_face(k)) + 0.5_R8P * cross(t(:, 2) - t(:, 1), t(:, 3) - t(:, 1))
      vb = vb + dot_product(t(:, 1), cross(t(:, 2), t(:, 3))) / 6._R8P
   enddo
   nfv = 0
   do f = 1, pt%nf
      if (norm2(fsum(:, f) - pt%face_vector(:, f)) > 1e-14_R8P * max(norm2(pt%face_vector(:, f)), h ** 2)) nfv = nfv + 1
   enddo
   ! per piece: closure, sum of the smallest faces against the wet normal
   cmax = 0._R8P ; wmax = 0._R8P ; nbit = 0 ; nrnd = 0 ; tot = 0._R8P ; comp = 0._R8P
   do p = 1, pt%cut%np
      sw = 0._R8P
      do k = pt%piece_sf_start(p), pt%piece_sf_start(p + 1) - 1
         sw = sw + pt%sf_area(:, k)
      enddo
      pf = pt%cut%piece_face(:, p)
      r = h * h * [pf(2) - pf(1), pf(4) - pf(3), pf(6) - pf(5)] + sw
      cmax = max(cmax, norm2(r) / h ** 2)
      if (all(sw == -pt%cut%piece_wet_normal(:, p))) then
         nbit = nbit + 1
      else
         nrnd = nrnd + 1 ; wmax = max(wmax, norm2(sw + pt%cut%piece_wet_normal(:, p)) / h ** 2)
      endif
      yk = pt%cut%piece_volume(p) - comp ; tk = tot + yk ; comp = (tk - tot) - yk ; tot = tk ! Kahan
   enddo
   ! every smallest face inside its major cell's box (the centroid), labels consistent
   nbox = 0
   do k = 1, pt%ns
      i = pt%sf_cell(1, k) ; j = pt%sf_cell(2, k) ; l = pt%sf_cell(3, k)
      lo = [x(i), x(j), x(l)] - 1e-12_R8P * h ; hi = [x(i + 1), x(j + 1), x(l + 1)] + 1e-12_R8P * h
      if (any(pt%sf_centroid(:, k) < lo) .or. any(pt%sf_centroid(:, k) > hi)) nbox = nbox + 1
      if (pt%sf_piece(k) > 0) then
         if (any(pt%cut%piece_cell(:, pt%sf_piece(k)) /= pt%sf_cell(:, k))) nbox = nbox + 1
      endif
      if (any(pt%sf_minor(:, k) /= pt%face_id(:, pt%sf_face(k))) .or. pt%tri_face(pt%sf_facet(k)) /= pt%sf_face(k)) &
         nbox = nbox + 1
   enddo
   ! per face: covered against the face vector
   fmax = 0._R8P ; nout = 0 ; nlong = 0 ; nin = 0
   do f = 1, pt%nf
      len = norm2(pt%face_vector(:, f))
      if (len == 0._R8P) cycle
      if (face_inside(pt, f, x)) then
         nin = nin + 1
         fmax = max(fmax, norm2(pt%face_covered(:, f) - pt%face_vector(:, f)) / len)
      else
         nout = nout + 1
         if (norm2(pt%face_covered(:, f)) > len * (1._R8P + 1e-13_R8P)) nlong = nlong + 1
      endif
   enddo
   ! degenerate pieces
   ndg = count(pt%degenerate) ; nbad = 0
   do p = 1, pt%cut%np
      if (.not. pt%degenerate(p)) cycle
      k = pt%merge_into(p)
      if (k <= 0) then
         nbad = nbad + 1 ; cycle
      endif
      if (pt%degenerate(k)) nbad = nbad + 1
      if (sum(abs(pt%cut%piece_cell(:, k) - pt%cut%piece_cell(:, p))) /= 1) nbad = nbad + 1
   enddo
   rel = 0._R8P ; if (expected > 0._R8P) rel = abs(tot - expected) / expected
   print '(3A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)', '  ', tag, ': holes ', count(pt%hole), ', pieces ', pt%cut%np, &
         ', faces ', pt%nf, ' (', pt%boundary%get_facets_number(), ' facets), smallest faces ', pt%ns, ', degenerate ', &
         ndg, ' (bad ', nbad, ')'
   print '(A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)', '    cut vs cut_cells_3d mismatches ', nmis, &
         ', facet sums vs face ', nfv, ', labels ', nbox, '; closure / h^2 ', cmax, '; covered vs face ', fmax, &
         '; volume rel ', rel, '; boundary vol ', vb
   print '(A,I0,A,I0,A,ES9.2,A,I0,A,I0,A,I0,A,I0)', '    sum sf = -wet_normal: bit ', nbit, ', round-off ', nrnd, &
         ' (max / h^2 ', wmax, '); faces inside ', nin, ', leaving the box ', nout, ' (longer than the face ', nlong, &
         '); untrusted cells ', pt%nuc
   ok = ok .and. nmis == 0 .and. nfv == 0 .and. nbox == 0 .and. cmax <= 1e-12_R8P .and. wmax <= 1e-13_R8P .and. &
        fmax <= 1e-13_R8P .and. nlong == 0 .and. rel <= 1e-13_R8P .and. nbad == 0 .and. pt%nuc == 0 .and. &
        all(pt%hole .eqv. (pt%cut%pieces == 0)) .and. vb > 0._R8P .and. pt%ns > 0
   if (inside) ok = ok .and. nout == 0
   if (.not. inside) ok = ok .and. nout > 0
   if (present(degenerate)) then
      if (degenerate) ok = ok .and. ndg > 0
   endif
   endsubroutine check_case

   logical function face_inside(pt, f, x) result(yes)
   !< Face f of the patch lies in the closed major box.
   type(overset_patch_3d_t), intent(in) :: pt   !< Patch.
   integer(I4P),             intent(in) :: f    !< Face.
   real(R8P),                intent(in) :: x(:) !< Grid edges.
   type(facet_object), pointer          :: fp   !< Facet.
   integer(I4P)                         :: k, l !< Work.
   real(R8P)                            :: v(3) !< Vertex.

   yes = .true.
   do k = 1, pt%boundary%get_facets_number()
      if (pt%tri_face(k) /= f) cycle
      fp => pt%boundary%facet_at(k)
      do l = 1, 3
         v = [fp%vertex(l)%x, fp%vertex(l)%y, fp%vertex(l)%z]
         if (any(v < x(1)) .or. any(v > x(size(x)))) yes = .false.
      enddo
   enddo
   endfunction face_inside

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   subroutine shell(na, c, r0, r1, blk, vol)
   !< A sphere shell block (i around, j pole to pole, k = 1, 2 radially: left-handed) centred at c between r0 and r1, its
   !< outer face (+k) named 'overlap'; vol: the volume of its outer polyhedron (shorter diagonals, poles collapsed).
   integer(I4P),              intent(in)    :: na    !< Cells around (na / 2 pole to pole).
   real(R8P),                 intent(in)    :: c(3)  !< Centre.
   real(R8P),                 intent(in)    :: r0, r1 !< Radii.
   type(curvilinear_block_t), intent(inout) :: blk   !< Block.
   real(R8P),                 intent(out)   :: vol   !< Outer polyhedron volume.
   type(topo_record_t)    :: rec                     !< Record.
   real(R8P), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:) !< Vertices.
   real(R8P)              :: th, ph, q(3, 4), w(3, 3, 2), rr !< Work.
   integer(I4P)           :: np, i, j, k, m              !< Counters.

   np = na / 2
   allocate(xx(na + 1, np + 1, 2), yy(na + 1, np + 1, 2), zz(na + 1, np + 1, 2))
   do k = 1, 2
      rr = merge(r0, r1, k == 1)
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            xx(i, j, k) = c(1) + rr * sin(ph) * cos(th) ; yy(i, j, k) = c(2) + rr * sin(ph) * sin(th)
            zz(i, j, k) = c(3) + rr * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            xx(:, j, k) = c(1) ; yy(:, j, k) = c(2) ; zz(:, j, k) = c(3) + merge(rr, -rr, j == 1)
         endif
      enddo
      xx(na + 1, :, k) = xx(1, :, k) ; yy(na + 1, :, k) = yy(1, :, k) ; zz(na + 1, :, k) = zz(1, :, k)
   enddo
   call blk%destroy
   call blk%initialize(xx, yy, zz)
   rec%kind = TOPO_BOUNDARY ; rec%face = 6 ; rec%lo = [1, 1, 2] ; rec%hi = [na + 1, np + 1, 2] ; rec%name = 'overlap'
   call blk%add_record(rec)
   vol = 0._R8P ! the outer polyhedron, built here (outward from the centre), about the centre
   do j = 1, np
      do i = 1, na
         q(:, 1) = [xx(i, j, 2), yy(i, j, 2), zz(i, j, 2)] - c ; q(:, 2) = [xx(i + 1, j, 2), yy(i + 1, j, 2), zz(i + 1, j, 2)] - c
         q(:, 3) = [xx(i + 1, j + 1, 2), yy(i + 1, j + 1, 2), zz(i + 1, j + 1, 2)] - c
         q(:, 4) = [xx(i, j + 1, 2), yy(i, j + 1, 2), zz(i, j + 1, 2)] - c
         if (norm2(q(:, 3) - q(:, 1)) <= norm2(q(:, 4) - q(:, 2))) then
            w(:, :, 1) = q(:, [1, 2, 3]) ; w(:, :, 2) = q(:, [1, 3, 4])
         else
            w(:, :, 1) = q(:, [1, 2, 4]) ; w(:, :, 2) = q(:, [2, 3, 4])
         endif
         do m = 1, 2
            vol = vol + abs(dot_product(w(:, 1, m), cross(w(:, 2, m), w(:, 3, m)))) / 6._R8P
         enddo
      enddo
   enddo
   endsubroutine shell

   subroutine box_block(a, b, n, blk)
   !< A Cartesian block [a, b]^3 of n cells a side (right-handed), its six faces named 'overlap'.
   real(R8P),                 intent(in)    :: a, b !< Bounds.
   integer(I4P),              intent(in)    :: n    !< Cells a side.
   type(curvilinear_block_t), intent(inout) :: blk  !< Block.
   type(topo_record_t)    :: rec                    !< Record.
   real(R8P), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:) !< Vertices.
   integer(I4P)           :: i, j, k, f             !< Counters.

   allocate(xx(n + 1, n + 1, n + 1), yy(n + 1, n + 1, n + 1), zz(n + 1, n + 1, n + 1))
   do k = 1, n + 1
      do j = 1, n + 1
         do i = 1, n + 1
            xx(i, j, k) = a + (b - a) * (i - 1) / n ; yy(i, j, k) = a + (b - a) * (j - 1) / n
            zz(i, j, k) = a + (b - a) * (k - 1) / n
         enddo
      enddo
   enddo
   xx(n + 1, :, :) = b ; yy(:, n + 1, :) = b ; zz(:, :, n + 1) = b
   call blk%destroy
   call blk%initialize(xx, yy, zz)
   do f = 1, 6
      rec%kind = TOPO_BOUNDARY ; rec%face = f ; rec%lo = [1, 1, 1] ; rec%hi = [n + 1, n + 1, n + 1] ; rec%name = 'overlap'
      rec%lo((f + 1) / 2) = merge(1_I4P, n + 1_I4P, mod(f, 2_I4P) == 1) ; rec%hi((f + 1) / 2) = rec%lo((f + 1) / 2)
      call blk%add_record(rec)
   enddo
   endsubroutine box_block

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
endprogram fossil_test_overset_patch_3d
