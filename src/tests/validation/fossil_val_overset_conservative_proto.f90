!< FOSSIL validation: design probe of the conservative overset interface geometry (issue #106, step 0).
program fossil_val_overset_conservative_proto
!< FOSSIL validation: design probe of the conservative overset interface geometry (issue #106, step 0).
!<
!< No library code: the existing cut cells measured on the overset geometries of `fossil_val_overset`.
!< A (patched interface, Wang 1995): the minor block's outer boundary as a closed curve / surface (solid = the minor
!<   region) cuts a Cartesian major grid: holes, cut cells, apertures W_f, volumes; the closure of every cut cell; the
!<   wall pieces (smallest faces) labelled by the minor boundary face they lie on (post-match), and the sum of a face's
!<   pieces against the face; the total fluid measure against box - minor region.
!< B (intersection weights, Picard et al. 2021): one receptor cell of the minor block as a closed curve / surface cuts
!<   the donor grid restricted to the receptor's box: |V_c n V_c'| per donor cell, their sum against the receptor, the
!<   intersection centroids; cost per receptor and per fringe.
!< Edge cases: a boundary along grid lines and through grid nodes, a receptor on grid lines, two minor blocks
!< (disjoint, overlapping), a boundary leaving the major box.
!< Exit status 0 if every check holds (the overlapping minor boundaries are reported, not checked).

use fossil, only : curvilinear_block_2d_t, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, curve_object, &
                   surface_stl_object, facet_object, wall_cutter, wall_cutter_2d, cut_cells_2d_t, cut_cells_3d_t, &
                   cut_cells_3d
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
integer(I4P)         :: m
logical              :: ok
real(R8P)            :: t_surface = 0._R8P, t_cut = 0._R8P !< B 3-D: time in the receptor surfaces, in the cuts [s].

ok = .true.
print '(A)', '== A 2-D: cylinder O-grid outer boundary (r = 2) in the background [-3, 3]^2'
do m = 1, 8
   if (m == 3 .or. m == 5 .or. m == 6 .or. m == 7) cycle
   call patch_2d_cylinder(m, ok)
enddo
call wall_cutter_orientation_2d(ok)
print '(A)', '== A 2-D: edge cases'
call edge_cases_2d(ok)
print '(A)', '== B 2-D: receptor cells of the O-grid against the background'
do m = 1, 8
   if (m == 3 .or. m == 5 .or. m == 6 .or. m == 7) cycle
   call weights_2d(m, ok)
enddo
call weights_2d_edge(ok)
print '(A)', '== A 3-D: sphere shell outer boundary (r = 2) in the background [-2.6, 2.6]^3'
call wall_cutter_orientation_3d(ok)
call patch_3d(25_I4P, 32_I4P, ok)
call patch_3d(49_I4P, 64_I4P, ok)
call patch_3d(97_I4P, 64_I4P, ok)
call cube_3d(ok)
print '(A)', '== B 3-D: receptor hexes of the shell against the background'
call weights_3d(49_I4P, ok)
call weights_3d(97_I4P, ok)
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   ! ------------------------------------------------------------------------------------------------------ helpers 2-D
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

   pure subroutine clip_halfplane(px, py, n, ax, ay, bx, by)
   !< Keep the part of polygon (px, py)(1:n) on the left of the line a -> b (Sutherland-Hodgman).
   real(R8P),    intent(inout) :: px(:), py(:) !< Vertices (capacity >= 2 n).
   integer(I4P), intent(inout) :: n            !< Vertices.
   real(R8P),    intent(in)    :: ax, ay, bx, by !< Line.
   real(R8P)                   :: qx(size(px)), qy(size(px)), s1, s2, t !< Work.
   integer(I4P)                :: k, k2, nq    !< Counters.

   nq = 0_I4P
   do k = 1, n
      k2 = mod(k, n) + 1_I4P
      s1 = (bx - ax) * (py(k) - ay) - (by - ay) * (px(k) - ax)
      s2 = (bx - ax) * (py(k2) - ay) - (by - ay) * (px(k2) - ax)
      if (s1 >= 0._R8P) then
         nq = nq + 1_I4P ; qx(nq) = px(k) ; qy(nq) = py(k)
      endif
      if ((s1 >= 0._R8P .and. s2 < 0._R8P) .or. (s1 < 0._R8P .and. s2 >= 0._R8P)) then
         t = s1 / (s1 - s2)
         nq = nq + 1_I4P ; qx(nq) = px(k) + t * (px(k2) - px(k)) ; qy(nq) = py(k) + t * (py(k2) - py(k))
      endif
   enddo
   n = nq ; px(1:n) = qx(1:n) ; py(1:n) = qy(1:n)
   endsubroutine clip_halfplane

   real(R8P) function clipped_area(px, py, cx, cy) result(a)
   !< Area of polygon p clipped by the convex CCW polygon c.
   real(R8P), intent(in)  :: px(:), py(:), cx(:), cy(:) !< Polygons.
   real(R8P), allocatable :: wx(:), wy(:) !< Work.
   integer(I4P)           :: n, k, nc     !< Counters.

   n = size(px, kind=I4P) ; nc = size(cx, kind=I4P)
   allocate(wx(2 * (n + nc) + 8), wy(2 * (n + nc) + 8)) ; wx(1:n) = px ; wy(1:n) = py
   do k = 1, nc
      call clip_halfplane(wx, wy, n, cx(k), cy(k), cx(mod(k, nc) + 1), cy(mod(k, nc) + 1))
      if (n == 0) exit
   enddo
   a = 0._R8P ; if (n > 2) a = shoelace(wx(1:n), wy(1:n))
   endfunction clipped_area

   pure function grid_1d(a, b, n) result(x)
   !< n + 1 equally spaced edges on [a, b].
   real(R8P),    intent(in) :: a, b !< Ends.
   integer(I4P), intent(in) :: n    !< Cells.
   real(R8P)                :: x(n + 1) !< Edges.
   integer(I4P)             :: i    !< Counter.

   do i = 1, n + 1
      x(i) = a + (b - a) * (i - 1) / n
   enddo
   x(n + 1) = b
   endfunction grid_1d

   subroutine analyze_2d(tag, pu, pv, lsz, x, y, expected, ok, check, timing)
   !< Cut the grid x, y by the closed loops (pu, pv) of sizes lsz (solid on the left) and measure A's geometry:
   !< counts, closure per cut cell and per piece, the labels of the wall pieces by post-match, the sum of every
   !< boundary face's pieces against the face (faces inside the grid box), and the total fluid area against `expected`.
   character(*), intent(in)    :: tag          !< Case name.
   real(R8P),    intent(in)    :: pu(:), pv(:) !< Loop vertices.
   integer(I4P), intent(in)    :: lsz(:)       !< Loop sizes.
   real(R8P),    intent(in)    :: x(:), y(:)   !< Grid edges.
   real(R8P),    intent(in)    :: expected     !< Expected fluid area.
   logical,      intent(inout) :: ok           !< Checks.
   logical,      intent(in)    :: check        !< Check (else report only).
   real(R8P),    intent(out), optional :: timing !< Cut time [s].
   type(curve_object)          :: c            !< Curve.
   type(cut_cells_2d_t)        :: cut          !< Cut cells.
   integer(I4P), allocatable   :: sa(:), sb(:), cnt(:) !< Segment ends, pieces per segment.
   real(R8P),    allocatable   :: sv(:,:)      !< Sum of the pieces per segment [2, ns].
   integer(I4P)  :: ns, l, k, s, nx, ny, i, j, p, q, nhole, ncut, nfluid, namb, nnone, nmulti, nwall, st, best !< Work.
   integer(I4P)  :: nzero, nfull, nsliver, nzl, ntiny, npart !< Degenerate counts.
   real(R8P)     :: dx, dy, r(2), rmax_piece, rmax_seg, dmax, tot, comp, yk, tk, a(2), b(2), d, t, tol, mid(2) !< Work.
   real(R8P)     :: e(2), rcell_sign(2), apmin, rel, len !< Work.
   integer(I8P)  :: t0, t1, cr !< Clock.

   ns = size(pu, kind=I4P)
   allocate(sa(ns), sb(ns), cnt(ns), sv(2, ns)) ; k = 0_I4P
   do l = 1, size(lsz, kind=I4P)
      do s = 1, lsz(l)
         sa(k + s) = k + s ; sb(k + s) = k + mod(s, lsz(l)) + 1_I4P
      enddo
      k = k + lsz(l)
   enddo
   call c%adopt_loops(u=pu, v=pv, loop_sizes=lsz, status=st)
   call system_clock(t0, cr)
   call c%cut_cells(x, y, cut, status=st, boundaries=.true.)
   call system_clock(t1)
   if (present(timing)) timing = real(t1 - t0, R8P) / cr
   nx = cut%nx ; ny = cut%ny
   nhole = 0 ; ncut = 0 ; nfluid = 0 ; rcell_sign = 0._R8P ; tot = 0._R8P ; comp = 0._R8P ; apmin = 1._R8P
   nzero = 0 ; nfull = 0 ; nsliver = 0
   do j = 1, ny
      do i = 1, nx
         dx = x(i + 1) - x(i) ; dy = y(j + 1) - y(j)
         if (cut%pieces(i, j) == 0) then
            nhole = nhole + 1 ; cycle
         endif
         yk = cut%alpha(i, j) * dx * dy - comp ; tk = tot + yk ; comp = (tk - tot) - yk ; tot = tk ! Kahan
         if (cut%wet_length(i, j) == 0._R8P) then
            nfluid = nfluid + 1 ; cycle
         endif
         ncut = ncut + 1
         if (cut%alpha(i, j) == 0._R8P) nzero = nzero + 1
         if (cut%alpha(i, j) == 1._R8P) nfull = nfull + 1
         if (cut%alpha(i, j) > 0._R8P .and. cut%alpha(i, j) < 1e-12_R8P) nsliver = nsliver + 1
         r = [(cut%aperture(2, i, j) - cut%aperture(1, i, j)) * dy, (cut%aperture(4, i, j) - cut%aperture(3, i, j)) * dx]
         rcell_sign(1) = max(rcell_sign(1), norm2(r + cut%wet_normal(:, i, j)) / max(dx, dy))
         rcell_sign(2) = max(rcell_sign(2), norm2(r - cut%wet_normal(:, i, j)) / max(dx, dy))
         if (any(cut%aperture(:, i, j) > 0._R8P)) &
            apmin = min(apmin, minval(cut%aperture(:, i, j), mask=cut%aperture(:, i, j) > 0._R8P))
      enddo
   enddo
   ! per piece: the oriented boundary closes; labels of the wall segments
   rmax_piece = 0._R8P ; cnt = 0 ; sv = 0._R8P ; namb = 0 ; nnone = 0 ; dmax = 0._R8P ; nwall = 0 ; nmulti = 0
   nzl = 0 ; ntiny = 0
   tol = 1e-12_R8P * max(x(size(x)) - x(1), y(size(y)) - y(1))
   do p = 1, cut%np
      r = 0._R8P
      do q = cut%piece_seg_start(p), cut%piece_seg_start(p + 1) - 1
         a = cut%piece_seg(1:2, q) ; b = cut%piece_seg(3:4, q)
         r = r + [b(2) - a(2), -(b(1) - a(1))]
         if (.not. cut%piece_seg_wall(q)) cycle
         nwall = nwall + 1 ; len = norm2(a - b)
         if (len == 0._R8P) then
            nzl = nzl + 1 ; cycle
         endif
         if (len < tol) ntiny = ntiny + 1
         mid = 0.5_R8P * (a + b) ; best = 0 ; d = huge(1._R8P) ; k = 0
         do s = 1, ns ! post-match: the boundary face whose segment holds the midpoint, same direction (curve: b -> a)
            e = [pu(sb(s)) - pu(sa(s)), pv(sb(s)) - pv(sa(s))]
            t = dot_product(mid - [pu(sa(s)), pv(sa(s))], e) / dot_product(e, e)
            if (t < -1e-9_R8P .or. t > 1._R8P + 1e-9_R8P) cycle
            if (dot_product(a - b, e) <= 0._R8P .and. len >= tol) cycle
            rel = norm2(mid - [pu(sa(s)), pv(sa(s))] - t * e)
            if (rel <= tol) k = k + 1
            if (rel < d) then
               d = rel ; best = s
            endif
         enddo
         if (k == 0) nnone = nnone + 1
         if (k > 1) namb = namb + 1
         if (best > 0) then
            dmax = max(dmax, d) ; cnt(best) = cnt(best) + 1 ; sv(:, best) = sv(:, best) + (a - b)
         endif
      enddo
      rmax_piece = max(rmax_piece, norm2(r) / max(x(2) - x(1), y(2) - y(1)))
   enddo
   rmax_seg = 0._R8P ; npart = 0
   do s = 1, ns
      e = [pu(sb(s)) - pu(sa(s)), pv(sb(s)) - pv(sa(s))]
      if (.not. (inside_box(pu(sa(s)), pv(sa(s)), x, y) .and. inside_box(pu(sb(s)), pv(sb(s)), x, y))) then
         npart = npart + 1 ; cycle
      endif
      rmax_seg = max(rmax_seg, norm2(sv(:, s) - e) / norm2(e))
      nmulti = max(nmulti, cnt(s))
   enddo
   print '(2A)', '  ', tag
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2)', '    cells ', nx * ny, ': holes ', nhole, ', cut ', ncut, &
         ', fluid ', nfluid, '; boundary faces ', ns, ' -> smallest faces ', nwall, '; smallest aperture > 0 ', apmin
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '    degenerate: cut cells with alpha = 0 ', nzero, ', alpha = 1 ', nfull, &
         ', 0 < alpha < 1e-12 ', nsliver, '; smallest faces of length 0 ', nzl, ', below 1e-12 L ', ntiny
   print '(A,2ES10.2,A,ES10.2)', '    closure per cut cell / h (W_f S_f n_f + wet_normal, - wet_normal): ', rcell_sign, &
         '; per piece (oriented boundary) / h: ', rmax_piece
   print '(A,I0,A,I0,A,ES9.2,A,I0,A,I0,A,ES9.2,A,I0,A)', '    labels by post-match (length > 0): ambiguous ', namb, &
         ', none ', nnone, ', max off-line ', dmax, '; faces with no piece ', count(cnt == 0), ', max pieces per face ', &
         nmulti, '; |sum pieces - face| / |face| ', rmax_seg, ' (', npart, ' faces leave the box)'
   print '(A,ES22.14,A,ES22.14,A,ES9.2)', '    fluid area ', tot, ' expected ', expected, ' rel diff ', &
         abs(tot - expected) / expected
   if (check) ok = ok .and. minval(rcell_sign) <= 1e-13_R8P .and. rmax_piece <= 1e-13_R8P .and. rmax_seg <= 1e-12_R8P &
                   .and. abs(tot - expected) <= 1e-13_R8P * expected
   endsubroutine analyze_2d

   pure logical function inside_box(pa, pb, x, y)
   !< The point lies in the box of the grid x, y.
   real(R8P), intent(in) :: pa, pb     !< Point.
   real(R8P), intent(in) :: x(:), y(:) !< Grid.

   inside_box = pa >= x(1) .and. pa <= x(size(x)) .and. pb >= y(1) .and. pb <= y(size(y))
   endfunction inside_box

   ! ----------------------------------------------------------------------------------------------------------- A 2-D
   subroutine patch_2d_cylinder(m, ok)
   !< The O-grid's outer boundary (r = 2, 32 m faces) on the background of 40 m cells per side.
   integer(I4P), intent(in)    :: m  !< Refinement.
   logical,      intent(inout) :: ok !< Checks.
   real(R8P), allocatable :: pu(:), pv(:), x(:) !< Boundary, grid.
   integer(I4P)           :: ns, s, nb          !< Counters.
   real(R8P)              :: tcut               !< Time.
   character(64)          :: tag                !< Name.

   ns = 32 * m ; nb = 40 * m
   allocate(pu(ns), pv(ns))
   do s = 1, ns
      pu(s) = 2._R8P * cos(2 * PI * (s - 1) / ns) ; pv(s) = 2._R8P * sin(2 * PI * (s - 1) / ns)
   enddo
   x = grid_1d(-3._R8P, 3._R8P, nb)
   write(tag, '(A,I0,A,I0,A)') 'm = ', m, ' (', nb, '^2)'
   call analyze_2d(trim(tag), pu, pv, [ns], x, x, 36._R8P - shoelace(pu, pv), ok, .true., tcut)
   print '(A,ES10.3,A)', '    cut cells with boundaries: ', tcut, ' s'
   endsubroutine patch_2d_cylinder

   subroutine wall_cutter_orientation_2d(ok)
   !< What `wall_cutter_2d(wall_name='overlap')` gives for the O-grid: the loop orientation (signed area).
   logical, intent(inout)       :: ok        !< Checks.
   type(curvilinear_block_2d_t) :: blocks(1) !< O-grid.
   type(topo_record_t)          :: rec       !< Record.
   type(curve_object)           :: wc        !< Cutter.
   real(R8P), allocatable       :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)                 :: ni, nj, i, j, st !< Counters.
   real(R8P)                    :: r         !< Radius.

   ni = 33 ; nj = 13
   allocate(u(ni, nj), v(ni, nj))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = r * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   call blocks(1)%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   call wall_cutter_2d(blocks, wc, wall_name='overlap', status=st)
   print '(A,I0,A,I0,A,I0,A,F10.6,A)', '  wall_cutter_2d(overlap): status ', st, ', loops ', wc%get_loops_number(), &
         ', vertices ', wc%get_vertices_number(), ', signed area ', wc%get_area(), &
         ' (negative: solid outside, the minor region is fluid: reverse it for A)'
   ok = ok .and. st == 0
   endsubroutine wall_cutter_orientation_2d

   subroutine edge_cases_2d(ok)
   !< A on the edge cases (module notes).
   logical, intent(inout) :: ok !< Checks.
   real(R8P), allocatable :: x(:), pu(:), pv(:), qu(:), qv(:), bx(:), by(:) !< Grid, loops, box.
   integer(I4P)           :: k, n, s                !< Counters.
   real(R8P)              :: a1, a2, lens, tot      !< Areas.

   x = grid_1d(-3._R8P, 3._R8P, 40_I4P) ! h = 0.15
   ! 1. a square on grid lines through grid nodes: [x(11), x(31)]^2, a vertex at every node it passes
   n = 0 ; allocate(pu(80), pv(80))
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
   call analyze_2d('square on grid lines, a vertex at every node', pu, pv, [n], x, x, 36._R8P - shoelace(pu, pv), ok, &
                   .true.)
   ! 2. the same square with 4 vertices only (long segments along grid lines through nodes)
   call analyze_2d('square on grid lines, 4 vertices', [x(11), x(31), x(31), x(11)], [x(11), x(11), x(31), x(31)], [4], &
                   x, x, 36._R8P - (x(31) - x(11)) ** 2, ok, .true.)
   ! 3. a diamond with its vertices on grid nodes: its edges pass through grid nodes (diagonals)
   call analyze_2d('diamond through grid nodes', [x(21), x(31), x(21), x(11)], [x(11), x(21), x(31), x(21)], [4], x, x, &
                   36._R8P - 0.5_R8P * (x(31) - x(11)) ** 2, ok, .true.)
   ! 4. a polygon (r = 2, 64 faces) whose vertices are snapped to grid lines in x (vertices on grid lines)
   deallocate(pu, pv) ; allocate(pu(64), pv(64))
   do s = 1, 64
      pu(s) = 2._R8P * cos(2 * PI * (s - 1) / 64) ; pv(s) = 2._R8P * sin(2 * PI * (s - 1) / 64)
      pu(s) = x(nint((pu(s) + 3._R8P) / 0.15_R8P) + 1)
   enddo
   call analyze_2d('64-gon, vertices snapped to grid lines in x', pu, pv, [64], x, x, 36._R8P - shoelace(pu, pv), ok, &
                   .true.)
   ! 5. two disjoint minor boundaries
   deallocate(pu, pv) ; allocate(pu(128), pv(128))
   do s = 1, 64
      pu(s) = -1.6_R8P + cos(2 * PI * (s - 1) / 64) ; pv(s) = 0.1_R8P + sin(2 * PI * (s - 1) / 64)
      pu(64 + s) = 1.55_R8P + cos(2 * PI * (s - 1) / 64) ; pv(64 + s) = -0.2_R8P + sin(2 * PI * (s - 1) / 64)
   enddo
   call analyze_2d('two disjoint minor boundaries', pu, pv, [64, 64], x, x, &
                   36._R8P - shoelace(pu(1:64), pv(1:64)) - shoelace(pu(65:128), pv(65:128)), ok, .true.)
   ! 6. two overlapping minor boundaries: the union is the hole; the cut cells' parity makes the lens fluid
   do s = 1, 64
      pu(64 + s) = -0.4_R8P + cos(2 * PI * (s - 1) / 64) ; pv(64 + s) = 0.3_R8P + sin(2 * PI * (s - 1) / 64)
   enddo
   a1 = shoelace(pu(1:64), pv(1:64)) ; a2 = shoelace(pu(65:128), pv(65:128))
   lens = clipped_area(pu(1:64), pv(1:64), pu(65:128), pv(65:128))
   call analyze_2d('two overlapping minor boundaries (expected: box - union)', pu, pv, [64, 64], x, x, &
                   36._R8P - a1 - a2 + lens, ok, .false.)
   print '(A,ES22.14,A,ES22.14)', '    box - union ', 36._R8P - a1 - a2 + lens, '; box - union + lens (parity) ', &
         36._R8P - a1 - a2 + 2 * lens
   ! 7. a minor boundary leaving the major box: circle r = 1.5 centred at (2.5, 0.4)
   deallocate(pu, pv) ; allocate(pu(96), pv(96))
   do s = 1, 96
      pu(s) = 2.5_R8P + 1.5_R8P * cos(2 * PI * (s - 1) / 96) ; pv(s) = 0.4_R8P + 1.5_R8P * sin(2 * PI * (s - 1) / 96)
   enddo
   bx = [-3._R8P, 3._R8P, 3._R8P, -3._R8P] ; by = [-3._R8P, -3._R8P, 3._R8P, 3._R8P]
   tot = 36._R8P - clipped_area(pu, pv, bx, by)
   call analyze_2d('a minor boundary leaving the major box', pu, pv, [96], x, x, tot, ok, .true.)
   endsubroutine edge_cases_2d

   ! ----------------------------------------------------------------------------------------------------------- B 2-D
   subroutine receptor_2d(qu, qv, x, y, mode, area_sum, mom_sum, nover, ok_cut)
   !< Intersection areas and centroids of the quad (qu, qv) (CCW) with the cells of the grid x, y (its box sub-grid is
   !< taken here). mode 1: the quad is the solid, |V_c n V_c'| = |V_c'| - fluid area (complement); mode 2: the quad is
   !< a hole in a solid box around the sub-grid, |V_c n V_c'| = the fluid pieces' areas (direct); mode 3: mode 2 in
   !< coordinates relative to the sub-grid's low corner (the cut cells' Green sums run about the origin).
   real(R8P),    intent(in)  :: qu(4), qv(4)   !< Receptor (CCW).
   real(R8P),    intent(in)  :: x(:), y(:)     !< Donor grid.
   integer(I4P), intent(in)  :: mode           !< 1 complement, 2 direct, 3 direct shifted.
   real(R8P),    intent(out) :: area_sum       !< Sum of the intersection areas.
   real(R8P),    intent(out) :: mom_sum(2)     !< Sum of the intersection first moments.
   integer(I4P), intent(out) :: nover          !< Donor cells with a non-empty intersection.
   logical,      intent(out) :: ok_cut         !< Every cell closed.
   type(curve_object)        :: c              !< Curve.
   type(cut_cells_2d_t)      :: cut            !< Cut.
   integer(I4P)              :: i0, i1, j0, j1, i, j, p, st !< Sub-grid, counters.
   real(R8P)                 :: dx, dy, ai, ci(2), fa, fm(2), bu(4), bv(4), o(2) !< Work.

   i0 = maxloc(x, dim=1, mask=x <= minval(qu)) ; i1 = minloc(x, dim=1, mask=x >= maxval(qu))
   j0 = maxloc(y, dim=1, mask=y <= minval(qv)) ; j1 = minloc(y, dim=1, mask=y >= maxval(qv))
   if (i1 == i0) i1 = i0 + 1
   if (j1 == j0) j1 = j0 + 1
   o = 0._R8P ; if (mode == 3) o = [x(i0), y(j0)]
   if (mode == 1) then
      call c%adopt_loops(u=qu, v=qv, loop_sizes=[4_I4P], status=st)
   else
      dx = x(i1) - x(i0) ; dy = y(j1) - y(j0)
      bu = [x(i0) - dx, x(i1) + dx, x(i1) + dx, x(i0) - dx] ; bv = [y(j0) - dy, y(j0) - dy, y(j1) + dy, y(j1) + dy]
      call c%adopt_loops(u=[bu, qu(4:1:-1)] - o(1), v=[bv, qv(4:1:-1)] - o(2), loop_sizes=[4_I4P, 4_I4P], status=st)
   endif
   call c%cut_cells(x(i0:i1) - o(1), y(j0:j1) - o(2), cut, status=st, intervals=.false.)
   area_sum = 0._R8P ; mom_sum = 0._R8P ; nover = 0 ; ok_cut = all(cut%closed)
   do j = 1, cut%ny
      do i = 1, cut%nx
         dx = x(i0 + i) - x(i0 + i - 1) ; dy = y(j0 + j) - y(j0 + j - 1)
         fa = 0._R8P ; fm = 0._R8P
         do p = cut%first_piece(i, j), cut%first_piece(i, j) + cut%pieces(i, j) - 1
            fa = fa + cut%piece_area(p) ; fm = fm + cut%piece_area(p) * cut%piece_centroid(:, p)
         enddo
         if (mode == 1) then
            ai = dx * dy - fa
            ci = dx * dy * [0.5_R8P * (x(i0 + i) + x(i0 + i - 1)), 0.5_R8P * (y(j0 + j) + y(j0 + j - 1))] - fm
         else
            ai = fa ; ci = fm + fa * o
         endif
         if (ai > 1e-12_R8P * dx * dy) nover = nover + 1 ! round-off slivers are not counted
         area_sum = area_sum + ai ; mom_sum = mom_sum + ci
      enddo
   enddo
   endsubroutine receptor_2d

   subroutine weights_2d(m, ok)
   !< B on the two outer layers of the O-grid (32 m x 12 m cells, r in [1, 2]) against the background (40 m cells).
   integer(I4P), intent(in)    :: m  !< Refinement.
   logical,      intent(inout) :: ok !< Checks.
   real(R8P), allocatable :: x(:)  !< Grid.
   real(R8P)              :: qu(4), qv(4), r0, r1, t0, t1, area, mom(2), asum, msum(2), emax(3), cmax(3) !< Work.
   integer(I4P)           :: ni, nj, i, j, mode, nover, ntot(3), nrec !< Counters.
   integer(I8P)           :: c0, c1, cr !< Clock.
   real(R8P)              :: tm(3)  !< Times.
   logical                :: okc, allc !< Closed.

   ni = 32 * m ; nj = 12 * m
   x = grid_1d(-3._R8P, 3._R8P, 40_I4P * m)
   emax = 0._R8P ; cmax = 0._R8P ; ntot = 0 ; nrec = 0 ; allc = .true.
   do mode = 1, 3
      call system_clock(c0, cr)
      do j = nj - 1, nj
         r0 = 1._R8P + (j - 1._R8P) / nj ; r1 = 1._R8P + real(j, R8P) / nj
         do i = 1, ni
            t0 = 2 * PI * (i - 1) / ni ; t1 = 2 * PI * i / ni
            qu = [r0 * cos(t0), r0 * cos(t1), r1 * cos(t1), r1 * cos(t0)]
            qv = [r0 * sin(t0), r0 * sin(t1), r1 * sin(t1), r1 * sin(t0)] ! (i, j) order: CW; reversed below
            qu = qu(4:1:-1) ; qv = qv(4:1:-1)
            area = shoelace(qu - qu(1), qv - qv(1)) ! about a vertex: the shoelace about the origin loses 1e-12
            call polygon_moments(qu - qu(1), qv - qv(1), mom) ; mom = mom + area * [qu(1), qv(1)]
            call receptor_2d(qu, qv, x, x, mode, asum, msum, nover, okc)
            allc = allc .and. okc
            emax(mode) = max(emax(mode), abs(asum - area) / area)
            cmax(mode) = max(cmax(mode), norm2(msum - mom) / (area * (r1 - r0)))
            ntot(mode) = ntot(mode) + nover
            if (mode == 1) nrec = nrec + 1
         enddo
      enddo
      call system_clock(c1)
      tm(mode) = real(c1 - c0, R8P) / cr
   enddo
   print '(A,I0,A,I0,A,3F6.2,A)', '  m = ', m, ': ', nrec, ' receptors (2 outer layers), ', &
         real(ntot, R8P) / nrec, ' donor cells each (complement, direct, shifted)'
   print '(A,3ES10.2,A,3ES10.2)', '    |sum - |V_c|| / |V_c| (complement, direct, direct shifted) ', emax, &
         '; first moment / (|V_c| dr) ', cmax
   print '(A,3F7.2,A,L1)', '    us per receptor ', 1e6_R8P * tm / nrec, '; all closed ', allc
   ok = ok .and. allc .and. maxval(emax) <= 1e-11_R8P .and. emax(3) <= 1e-14_R8P .and. cmax(3) <= 1e-12_R8P .and. &
        ntot(1) == ntot(2) .and. ntot(1) == ntot(3)
   endsubroutine weights_2d

   pure subroutine polygon_moments(px, py, mom)
   !< First moments int x dA, int y dA of a polygon (Green).
   real(R8P), intent(in)  :: px(:), py(:) !< Vertices.
   real(R8P), intent(out) :: mom(2)       !< Moments.
   integer(I4P)           :: k, k2, n     !< Counters.
   real(R8P)              :: cr           !< Cross.

   n = size(px, kind=I4P) ; mom = 0._R8P
   do k = 1, n
      k2 = mod(k, n) + 1_I4P
      cr = px(k) * py(k2) - px(k2) * py(k)
      mom = mom + cr * [px(k) + px(k2), py(k) + py(k2)]
   enddo
   mom = mom / 6._R8P
   endsubroutine polygon_moments

   subroutine weights_2d_edge(ok)
   !< B edge cases: a receptor on grid lines (vertices on nodes), a receptor with one vertex on a node and edges along
   !< a grid line, a receptor inside one donor cell.
   logical, intent(inout) :: ok !< Checks.
   real(R8P), allocatable :: x(:) !< Grid.
   real(R8P)              :: qu(4), qv(4), asum, msum(2), mom(2), area, e(3) !< Work.
   integer(I4P)           :: mode, nover(3), k !< Counters.
   logical                :: okc(3)  !< Closed.
   character(48)          :: tags(3) !< Names.

   x = grid_1d(-3._R8P, 3._R8P, 40_I4P)
   tags = [character(48) :: 'receptor = 2 x 2 donor cells', 'receptor on a grid line, vertex on a node', &
                            'receptor inside one donor cell']
   do k = 1, 3
      select case(k)
      case(1)
         qu = [x(3), x(5), x(5), x(3)] ; qv = [x(2), x(2), x(4), x(4)]
      case(2)
         qu = [x(3), x(5), x(5) + 0.07_R8P, x(3) - 0.01_R8P] ; qv = [x(2), x(2), x(4) - 0.02_R8P, x(3) + 0.03_R8P]
      case(3)
         qu = x(7) + [0.02_R8P, 0.1_R8P, 0.12_R8P, 0.03_R8P] ; qv = x(9) + [0.01_R8P, 0.02_R8P, 0.13_R8P, 0.1_R8P]
      endselect
      area = shoelace(qu - qu(1), qv - qv(1))
      call polygon_moments(qu - qu(1), qv - qv(1), mom) ; mom = mom + area * [qu(1), qv(1)]
      do mode = 1, 3
         call receptor_2d(qu, qv, x, x, mode, asum, msum, nover(mode), okc(mode))
         e(mode) = abs(asum - area) / area
      enddo
      print '(3A,3(I0,1X),A,3ES10.2,A,3L2)', '  ', tags(k), ': donor cells (complement, direct, shifted) ', nover, &
            '; |sum - |V_c|| / |V_c| ', e, '; closed ', okc
      ok = ok .and. all(okc) .and. maxval(e) <= 1e-13_R8P .and. nover(1) == nover(2) .and. nover(1) == nover(3)
   enddo
   endsubroutine weights_2d_edge

   ! ------------------------------------------------------------------------------------------------------ helpers 3-D
   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   subroutine make_surface(tri, srf)
   !< A surface from triangles [3, 3, n].
   real(R8P),                intent(in)    :: tri(:,:,:) !< Triangles.
   type(surface_stl_object), intent(inout) :: srf        !< Surface.
   type(facet_object), allocatable         :: f(:)       !< Facets.
   integer(I4P)                            :: t, k       !< Counters.

   allocate(f(size(tri, 3)))
   do t = 1, size(tri, 3, kind=I4P)
      do k = 1, 3
         f(t)%vertex(k) = vector_R8P(tri(1, k, t), tri(2, k, t), tri(3, k, t))
      enddo
   enddo
   call srf%adopt_facets(facets=f)
   endsubroutine make_surface

   pure subroutine quad_tris(c, tri, nt)
   !< The two triangles of the quad c(:, 1:4) (cyclic) split along the shorter diagonal (`wall_cutter`'s rule), the
   !< collapsed ones dropped; orientation of the quad kept.
   real(R8P),    intent(in)  :: c(3, 4)      !< Corners.
   real(R8P),    intent(out) :: tri(3, 3, 2) !< Triangles.
   integer(I4P), intent(out) :: nt           !< Triangles kept.
   real(R8P)                 :: t2(3, 3, 2)  !< Work.
   integer(I4P)              :: k            !< Counter.

   if (norm2(c(:, 3) - c(:, 1)) <= norm2(c(:, 4) - c(:, 2))) then
      t2(:, :, 1) = c(:, [1, 2, 3]) ; t2(:, :, 2) = c(:, [1, 3, 4])
   else
      t2(:, :, 1) = c(:, [1, 2, 4]) ; t2(:, :, 2) = c(:, [2, 3, 4])
   endif
   nt = 0
   do k = 1, 2
      if (norm2(cross3(t2(:, 2, k) - t2(:, 1, k), t2(:, 3, k) - t2(:, 1, k))) == 0._R8P) cycle
      nt = nt + 1 ; tri(:, :, nt) = t2(:, :, k)
   enddo
   endsubroutine quad_tris

   pure subroutine clip_tri_box(t, lo, hi, av, n)
   !< Area vector av = 1/2 sum p_k x p_k+1 of triangle t clipped to the box [lo, hi] (Sutherland-Hodgman); n vertices.
   real(R8P),    intent(in)  :: t(3, 3)      !< Triangle.
   real(R8P),    intent(in)  :: lo(3), hi(3) !< Box.
   real(R8P),    intent(out) :: av(3)        !< Area vector.
   integer(I4P), intent(out) :: n            !< Vertices left.
   real(R8P)                 :: p(3, 16), q(3, 16), s1, s2, w !< Work.
   integer(I4P)              :: ax, sd, k, k2, nq, o !< Counters.

   n = 3 ; p(:, 1:3) = t ; av = 0._R8P
   do ax = 1, 3
      do sd = 1, 2
         nq = 0
         do k = 1, n
            k2 = mod(k, n) + 1
            if (sd == 1) then
               s1 = p(ax, k) - lo(ax) ; s2 = p(ax, k2) - lo(ax)
            else
               s1 = hi(ax) - p(ax, k) ; s2 = hi(ax) - p(ax, k2)
            endif
            if (s1 >= 0._R8P) then
               nq = nq + 1 ; q(:, nq) = p(:, k)
            endif
            if ((s1 >= 0._R8P .and. s2 < 0._R8P) .or. (s1 < 0._R8P .and. s2 >= 0._R8P)) then
               w = s1 / (s1 - s2) ; nq = nq + 1 ; q(:, nq) = p(:, k) + w * (p(:, k2) - p(:, k))
               o = ax ; q(o, nq) = merge(lo(ax), hi(ax), sd == 1)
            endif
         enddo
         n = nq ; p(:, 1:n) = q(:, 1:n)
         if (n < 3) then
            n = 0 ; return
         endif
      enddo
   enddo
   do k = 2, n - 1
      av = av + 0.5_R8P * cross3(p(:, k) - p(:, 1), p(:, k + 1) - p(:, 1))
   enddo
   endsubroutine clip_tri_box

   ! ----------------------------------------------------------------------------------------------------------- A 3-D
   subroutine wall_cutter_orientation_3d(ok)
   !< What `wall_cutter(wall_name='overlap')` gives for the shell's outer face: facets and signed volume.
   logical, intent(inout)    :: ok        !< Checks.
   type(curvilinear_block_t) :: blocks(1) !< Shell.
   type(topo_record_t)       :: rec       !< Record.
   type(surface_stl_object)  :: wc        !< Cutter.
   real(R8P), allocatable    :: xi(:,:,:), x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)              :: k, st     !< Counters.
   real(R8P)                 :: gap       !< Bilinear gap.

   allocate(x(65, 33, 3), y(65, 33, 3), z(65, 33, 3))
   do k = 1, 3
      call sphere_vertices(64_I4P, 32_I4P, 1.5_R8P + 0.25_R8P * k, xi)
      x(:, :, k) = xi(1, :, :) ; y(:, :, k) = xi(2, :, :) ; z(:, :, k) = xi(3, :, :)
   enddo
   call blocks(1)%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 6 ; rec%lo = [1, 1, 3] ; rec%hi = [65, 33, 3] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   call wall_cutter(blocks, wc, wall_name='overlap', gap=gap, status=st)
   print '(A,I0,A,I0,A,F10.6,A,ES9.2,A)', '  wall_cutter(overlap): status ', st, ', facets ', wc%get_facets_number(), &
         ', signed volume ', wc%get_volume(), ', bilinear gap ', gap, ' (negative: the minor region is fluid: reverse it)'
   ok = ok .and. st == 0
   endsubroutine wall_cutter_orientation_3d

   subroutine sphere_vertices(na, np, r, x)
   !< The vertices of the sphere of radius r of the shell grid, poles collapsed: x [3, na + 1, np + 1].
   integer(I4P),           intent(in)  :: na, np !< Cells around, pole to pole.
   real(R8P),              intent(in)  :: r      !< Radius.
   real(R8P), allocatable, intent(out) :: x(:,:,:) !< Vertices.
   integer(I4P)                        :: i, j   !< Counters.
   real(R8P)                           :: th, ph !< Angles.

   allocate(x(3, na + 1, np + 1))
   do j = 1, np + 1
      ph = PI * (j - 1) / np
      do i = 1, na + 1
         th = 2 * PI * (i - 1) / na
         x(:, i, j) = r * [sin(ph) * cos(th), sin(ph) * sin(th), cos(ph)]
      enddo
      if (j == 1 .or. j == np + 1) x(:, :, j) = spread([0._R8P, 0._R8P, merge(r, -r, j == 1)], 2, na + 1)
   enddo
   x(:, na + 1, :) = x(:, 1, :)
   endsubroutine sphere_vertices

   subroutine patch_3d(nb, na, ok)
   !< The shell's outer boundary (r = 2, na x na/2 faces, shorter diagonal, outward) on a background of nb - 1 cells per
   !< side over [-2.6, 2.6]^3.
   integer(I4P), intent(in)    :: nb, na !< Background nodes per side, faces around.
   logical,      intent(inout) :: ok     !< Checks.
   type(surface_stl_object)    :: srf    !< Minor boundary.
   type(cut_cells_3d_t)        :: cut    !< Cut.
   real(R8P),    allocatable   :: xs(:,:,:), tri(:,:,:), x(:), fsum(:,:), csum(:,:,:,:) !< Work.
   integer(I4P), allocatable   :: flab(:,:) !< Facet -> (i, j) of the minor face.
   real(R8P)    :: c(3, 4), t2(3, 3, 2), r(3), h, dx, vol, tot, exp_vol, rs(2), av(3), lo(3), hi(3), rlab(2) !< Work.
   real(R8P)    :: comp, yk, tk !< Kahan.
   integer(I4P) :: np, i, j, k, nt, t, n, nh, ncut, nfl, i0(3), i1(3), a, cc, npair, nmax, nv, nun !< Counters.
   integer(I4P), allocatable :: pairs_per_cell(:,:,:) !< Smallest faces per cut cell.
   integer(I8P) :: c0, c1, c2, cr !< Clock.

   np = na / 2
   call sphere_vertices(na, np, 2._R8P, xs)
   allocate(tri(3, 3, 2 * na * np), flab(2, 2 * na * np)) ; n = 0
   do j = 1, np
      do i = 1, na
         c(:, 1) = xs(:, i, j) ; c(:, 2) = xs(:, i + 1, j) ; c(:, 3) = xs(:, i + 1, j + 1) ; c(:, 4) = xs(:, i, j + 1)
         call quad_tris(c, t2, nt)
         do t = 1, nt
            n = n + 1 ; tri(:, :, n) = t2(:, :, t) ; flab(:, n) = [i, j]
            r = cross3(t2(:, 2, t) - t2(:, 1, t), t2(:, 3, t) - t2(:, 1, t))
            if (dot_product(r, sum(t2(:, :, t), dim=2)) < 0._R8P) tri(:, :, n) = t2(:, [1, 3, 2], t) ! outward
         enddo
      enddo
   enddo
   call make_surface(tri(:, :, 1:n), srf)
   exp_vol = srf%get_volume()
   x = grid_1d(-2.6_R8P, 2.6_R8P, nb - 1_I4P) ; h = x(2) - x(1)
   call system_clock(c0, cr)
   call cut_cells_3d(srf, x, x, x, cut)
   call system_clock(c1)
   nh = 0 ; ncut = 0 ; nfl = 0 ; tot = 0._R8P ; comp = 0._R8P ; rs = 0._R8P ; nun = 0
   allocate(csum(3, nb - 1, nb - 1, nb - 1), pairs_per_cell(nb - 1, nb - 1, nb - 1)) ; csum = 0._R8P
   pairs_per_cell = 0
   do k = 1, nb - 1
      do j = 1, nb - 1
         do i = 1, nb - 1
            vol = h ** 3
            if (cut%pieces(i, j, k) < 0) nun = nun + 1
            if (cut%pieces(i, j, k) == 0) then
               nh = nh + 1 ; cycle
            endif
            yk = cut%alpha(i, j, k) * vol - comp ; tk = tot + yk ; comp = (tk - tot) - yk ; tot = tk ! Kahan
            if (cut%wet_area(i, j, k) == 0._R8P) then
               nfl = nfl + 1 ; cycle
            endif
            ncut = ncut + 1
            r = h * h * [cut%face_x(i + 1, j, k) - cut%face_x(i, j, k), cut%face_y(i, j + 1, k) - cut%face_y(i, j, k), &
                         cut%face_z(i, j, k + 1) - cut%face_z(i, j, k)]
            rs(1) = max(rs(1), norm2(r + cut%wet_normal(:, i, j, k)) / h ** 2)
            rs(2) = max(rs(2), norm2(r - cut%wet_normal(:, i, j, k)) / h ** 2)
         enddo
      enddo
   enddo
   ! labels: every facet clipped by the cells of its box (the smallest faces), summed per cell and per facet
   call system_clock(c1)
   allocate(fsum(3, n)) ; fsum = 0._R8P ; npair = 0 ; nmax = 0
   do t = 1, n
      do a = 1, 3
         i0(a) = max(1, min(nb - 1, floor((minval(tri(a, :, t)) - x(1)) / h) + 1))
         i1(a) = max(1, min(nb - 1, floor((maxval(tri(a, :, t)) - x(1)) / h) + 1))
      enddo
      do k = i0(3), i1(3)
         do j = i0(2), i1(2)
            do i = i0(1), i1(1)
               lo = [x(i), x(j), x(k)] ; hi = [x(i + 1), x(j + 1), x(k + 1)]
               call clip_tri_box(tri(:, :, t), lo, hi, av, nv)
               if (nv == 0 .or. norm2(av) == 0._R8P) cycle
               npair = npair + 1 ; pairs_per_cell(i, j, k) = pairs_per_cell(i, j, k) + 1
               csum(:, i, j, k) = csum(:, i, j, k) + av ; fsum(:, t) = fsum(:, t) + av
            enddo
         enddo
      enddo
   enddo
   call system_clock(c2)
   rlab = 0._R8P
   do t = 1, n
      av = 0.5_R8P * cross3(tri(:, 2, t) - tri(:, 1, t), tri(:, 3, t) - tri(:, 1, t))
      rlab(1) = max(rlab(1), norm2(fsum(:, t) - av) / norm2(av))
   enddo
   do k = 1, nb - 1
      do j = 1, nb - 1
         do i = 1, nb - 1
            if (cut%pieces(i, j, k) <= 0) cycle
            rlab(2) = max(rlab(2), norm2(csum(:, i, j, k) - cut%wet_normal(:, i, j, k)) / h ** 2)
         enddo
      enddo
   enddo
   nmax = maxval(pairs_per_cell)
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', '  background ', nb - 1, '^3, minor faces ', na, ' x ', np, &
         ' (', n, ' facets): holes ', nh, ', cut ', ncut, ', fluid ', nfl
   print '(A,2ES10.2,A,ES22.14,A,ES22.14,A,ES9.2)', '    closure per cut cell / h^2 (faces + wet_normal, - wet_normal) ', &
         rs, '; fluid volume ', tot, ' expected ', 5.2_R8P ** 3 - exp_vol, ' rel ', &
         abs(tot - (5.2_R8P ** 3 - exp_vol)) / (5.2_R8P ** 3 - exp_vol)
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2)', '    smallest faces (facet x cell clips) ', npair, ', max per cut cell ', nmax, &
         '; |sum over cells - facet| / |facet| ', rlab(1), '; |sum over facets - wet_normal| / h^2 ', rlab(2)
   print '(A,F8.3,A,F8.3,A)', '    cut_cells_3d ', real(c1 - c0, R8P) / cr, ' s; facet clips ', real(c2 - c1, R8P) / cr, ' s'
   dx = abs(tot - (5.2_R8P ** 3 - exp_vol)) / (5.2_R8P ** 3 - exp_vol)
   ok = ok .and. minval(rs) <= 1e-12_R8P .and. dx <= 1e-13_R8P .and. rlab(1) <= 1e-12_R8P .and. rlab(2) <= 1e-12_R8P
   cc = count(pairs_per_cell > 0 .and. cut%wet_area <= 0._R8P) + count(pairs_per_cell == 0 .and. cut%wet_area > 0._R8P)
   print '(A,I0,A,I0)', '    cells with clips but no wall, or wall but no clips: ', cc, '; untrusted cells (pieces -1) ', nun
   ok = ok .and. cc == 0 .and. nun == 0
   endsubroutine patch_3d

   ! ----------------------------------------------------------------------------------------------------------- B 3-D
   subroutine cube_3d(ok)
   !< A 3-D on a minor boundary along grid planes: the box [x(11), x(31)]^3 (12 triangles) on the 40^3 grid over
   !< [-3, 3]^3: holes, cut cells (zero-volume ones), closure, total volume.
   logical, intent(inout)   :: ok   !< Checks.
   type(surface_stl_object) :: srf  !< Box.
   type(cut_cells_3d_t)     :: cut  !< Cut.
   real(R8P), allocatable   :: x(:) !< Grid.
   real(R8P)                :: cv(3, 8), tt(3, 3, 12), r(3), h, rs, tot, comp, yk, tk, ex !< Work.
   integer(I4P)             :: q, n, i, j, k, nh, ncut, nz, nfl, nun !< Counters.

   x = grid_1d(-3._R8P, 3._R8P, 40_I4P) ; h = x(2) - x(1)
   do q = 1, 8
      cv(:, q) = merge(x(31), x(11), [btest(q - 1, 0), btest(q - 1, 1), btest(q - 1, 2)])
   enddo
   n = 0 ; call box_tris(cv, tt, n) ; call make_surface(tt, srf)
   call cut_cells_3d(srf, x, x, x, cut)
   nh = 0 ; ncut = 0 ; nz = 0 ; nfl = 0 ; nun = 0 ; rs = 0._R8P ; tot = 0._R8P ; comp = 0._R8P
   do k = 1, 40
      do j = 1, 40
         do i = 1, 40
            if (cut%pieces(i, j, k) < 0) nun = nun + 1
            if (cut%pieces(i, j, k) == 0) then
               nh = nh + 1 ; cycle
            endif
            yk = cut%alpha(i, j, k) * h ** 3 - comp ; tk = tot + yk ; comp = (tk - tot) - yk ; tot = tk
            if (cut%wet_area(i, j, k) == 0._R8P) then
               nfl = nfl + 1 ; cycle
            endif
            ncut = ncut + 1 ; if (cut%alpha(i, j, k) == 0._R8P) nz = nz + 1
            r = h * h * [cut%face_x(i + 1, j, k) - cut%face_x(i, j, k), cut%face_y(i, j + 1, k) - cut%face_y(i, j, k), &
                         cut%face_z(i, j, k + 1) - cut%face_z(i, j, k)]
            rs = max(rs, norm2(r - cut%wet_normal(:, i, j, k)) / h ** 2)
         enddo
      enddo
   enddo
   ex = 216._R8P - (x(31) - x(11)) ** 3
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2)', '  box on grid planes (40^3): holes ', nh, ', cut ', ncut, &
         ' (alpha = 0: ', nz, '), fluid ', nfl, ', untrusted ', nun, '; closure / h^2 ', rs, '; volume rel diff ', &
         abs(tot - ex) / ex
   ok = ok .and. rs <= 1e-13_R8P .and. abs(tot - ex) <= 1e-14_R8P * ex
   endsubroutine cube_3d


   subroutine hex_tris(v, tri, nt)
   !< The 12 triangles of the hex with corners v(:, i, j, k) (i, j, k = 1, 2), shorter diagonals, outward if the hex
   !< has positive volume (oriented here), collapsed ones dropped.
   real(R8P),    intent(in)  :: v(3, 2, 2, 2)  !< Corners.
   real(R8P),    intent(out) :: tri(3, 3, 12)  !< Triangles.
   integer(I4P), intent(out) :: nt             !< Kept.
   real(R8P)    :: c(3, 4), t2(3, 3, 2), cen(3), nrm(3), vol !< Work.
   integer(I4P) :: f, k, n2                    !< Counters.

   nt = 0 ; cen = sum(sum(sum(v, dim=4), dim=3), dim=2) / 8._R8P
   do f = 1, 6
      select case(f)
      case(1) ; c(:, 1) = v(:, 1, 1, 1) ; c(:, 2) = v(:, 1, 2, 1) ; c(:, 3) = v(:, 1, 2, 2) ; c(:, 4) = v(:, 1, 1, 2)
      case(2) ; c(:, 1) = v(:, 2, 1, 1) ; c(:, 2) = v(:, 2, 2, 1) ; c(:, 3) = v(:, 2, 2, 2) ; c(:, 4) = v(:, 2, 1, 2)
      case(3) ; c(:, 1) = v(:, 1, 1, 1) ; c(:, 2) = v(:, 2, 1, 1) ; c(:, 3) = v(:, 2, 1, 2) ; c(:, 4) = v(:, 1, 1, 2)
      case(4) ; c(:, 1) = v(:, 1, 2, 1) ; c(:, 2) = v(:, 2, 2, 1) ; c(:, 3) = v(:, 2, 2, 2) ; c(:, 4) = v(:, 1, 2, 2)
      case(5) ; c(:, 1) = v(:, 1, 1, 1) ; c(:, 2) = v(:, 2, 1, 1) ; c(:, 3) = v(:, 2, 2, 1) ; c(:, 4) = v(:, 1, 2, 1)
      case(6) ; c(:, 1) = v(:, 1, 1, 2) ; c(:, 2) = v(:, 2, 1, 2) ; c(:, 3) = v(:, 2, 2, 2) ; c(:, 4) = v(:, 1, 2, 2)
      endselect
      call quad_tris(c, t2, n2)
      do k = 1, n2 ! outward from the hex centre (a convex-enough cell)
         nrm = cross3(t2(:, 2, k) - t2(:, 1, k), t2(:, 3, k) - t2(:, 1, k))
         nt = nt + 1 ; tri(:, :, nt) = t2(:, :, k)
         if (dot_product(nrm, sum(t2(:, :, k), dim=2) / 3._R8P - cen) < 0._R8P) tri(:, :, nt) = t2(:, [1, 3, 2], k)
      enddo
   enddo
   vol = 0._R8P
   do k = 1, nt
      vol = vol + dot_product(tri(:, 1, k), cross3(tri(:, 2, k), tri(:, 3, k))) / 6._R8P
   enddo
   if (vol < 0._R8P) tri(:, :, 1:nt) = tri(:, [1, 3, 2], 1:nt)
   endsubroutine hex_tris

   subroutine receptor_3d(tri, nt, x, mode, vsum, msum, nover, okc, vhex, chex, nunt)
   !< Intersection volumes of the closed triangulated receptor (outward) with the cells of the grid x^3 restricted to
   !< its box. mode 1: the receptor is the solid (complement); mode 2: the receptor is a hole in a solid box (direct).
   real(R8P),    intent(in)  :: tri(:,:,:)   !< Triangles [3, 3, >= nt].
   integer(I4P), intent(in)  :: nt           !< Triangles.
   real(R8P),    intent(in)  :: x(:)         !< Donor grid (same in x, y, z).
   integer(I4P), intent(in)  :: mode         !< 1 complement, 2 direct.
   real(R8P),    intent(out) :: vsum, msum(3) !< Sum of the volumes, first moments.
   integer(I4P), intent(out) :: nover        !< Donor cells with a non-empty intersection.
   logical,      intent(out) :: okc          !< Every cell closed.
   real(R8P),    intent(out) :: vhex, chex(3) !< Receptor polyhedron volume and centroid (surface).
   integer(I4P), intent(out) :: nunt         !< Untrusted cells (pieces = -1) met.
   type(surface_stl_object)  :: srf          !< Surface.
   type(cut_cells_3d_t)      :: cut          !< Cut.
   real(R8P), allocatable    :: tt(:,:,:)    !< Triangles.
   real(R8P)                 :: lo(3), hi(3), bl(3), bh(3), d(3), fv, fm(3), vi, ci(3), h(3), cv(3, 8) !< Work.
   integer(I4P)              :: i0(3), i1(3), a, i, j, k, q, n !< Counters.
   integer(I8P)              :: c0, c1, c2, cr !< Clock.
   type(vector_R8P)          :: cen          !< Centroid.

   do a = 1, 3
      lo(a) = minval(tri(a, :, 1:nt)) ; hi(a) = maxval(tri(a, :, 1:nt))
      i0(a) = maxloc(x, dim=1, mask=x <= lo(a)) ; i1(a) = minloc(x, dim=1, mask=x >= hi(a))
      if (i1(a) == i0(a)) i1(a) = i0(a) + 1
   enddo
   call make_surface(tri(:, :, 1:nt), srf)
   vhex = srf%get_volume() ; cen = srf%get_centroid() ; chex = [cen%x, cen%y, cen%z]
   if (mode == 2) then ! a solid box around the sub-grid, the receptor reversed inside it
      bl = [x(i0(1)), x(i0(2)), x(i0(3))] ; bh = [x(i1(1)), x(i1(2)), x(i1(3))] ; d = bh - bl ; bl = bl - d ; bh = bh + d
      do q = 1, 8
         cv(:, q) = merge(bh, bl, [btest(q - 1, 0), btest(q - 1, 1), btest(q - 1, 2)])
      enddo
      allocate(tt(3, 3, nt + 12))
      tt(:, :, 1:nt) = tri(:, [1, 3, 2], 1:nt) ; n = nt
      call box_tris(cv, tt, n)
   else
      tt = tri(:, :, 1:nt) ; n = nt
   endif
   call system_clock(c0, cr)
   call make_surface(tt(:, :, 1:n), srf)
   call system_clock(c1)
   call cut_cells_3d(srf, x(i0(1):i1(1)), x(i0(2):i1(2)), x(i0(3):i1(3)), cut)
   call system_clock(c2)
   t_surface = t_surface + real(c1 - c0, R8P) / cr ; t_cut = t_cut + real(c2 - c1, R8P) / cr
   okc = all(cut%closed) ; vsum = 0._R8P ; msum = 0._R8P ; nover = 0 ; nunt = 0
   do k = 1, cut%nz
      do j = 1, cut%ny
         do i = 1, cut%nx
            h = [x(i0(1) + i) - x(i0(1) + i - 1), x(i0(2) + j) - x(i0(2) + j - 1), x(i0(3) + k) - x(i0(3) + k - 1)]
            fv = cut%alpha(i, j, k) * product(h) ; fm = fv * cut%centroid(:, i, j, k) ! cell totals (untrusted cells too)
            if (cut%pieces(i, j, k) < 0) nunt = nunt + 1
            if (mode == 1) then
               vi = product(h) - fv
               ci = product(h) * [x(i0(1) + i - 1), x(i0(2) + j - 1), x(i0(3) + k - 1)] + product(h) * 0.5_R8P * h - fm
            else
               vi = fv ; ci = fm
            endif
            if (vi > 1e-12_R8P * product(h)) nover = nover + 1 ! round-off slivers are not counted
            vsum = vsum + vi ; msum = msum + ci
         enddo
      enddo
   enddo
   endsubroutine receptor_3d

   pure subroutine box_tris(cv, tt, n)
   !< Append the 12 outward triangles of the box with corners cv(:, 1 + bits(x, y, z)).
   real(R8P),    intent(in)    :: cv(3, 8)    !< Corners.
   real(R8P),    intent(inout) :: tt(:,:,:)   !< Triangles.
   integer(I4P), intent(inout) :: n           !< Count.
   integer(I4P), parameter     :: fq(4, 6) = reshape([1, 5, 7, 3,  2, 4, 8, 6,  1, 2, 6, 5,  3, 7, 8, 4, &
                                                      1, 3, 4, 2,  5, 6, 8, 7], [4, 6]) !< Outward quads.
   integer(I4P)                :: f            !< Face.

   do f = 1, 6
      n = n + 1 ; tt(:, :, n) = cv(:, fq([1, 2, 3], f))
      n = n + 1 ; tt(:, :, n) = cv(:, fq([1, 3, 4], f))
   enddo
   endsubroutine box_tris

   subroutine weights_3d(nb, ok)
   !< B on the outer layer of the shell (64 x 32 x 16 cells, r in [1, 2]) against the background nb^3 nodes.
   integer(I4P), intent(in)    :: nb !< Background nodes per side.
   logical,      intent(inout) :: ok !< Checks.
   real(R8P), allocatable :: xi(:,:,:), xo(:,:,:), x(:) !< Spheres, grid.
   real(R8P)              :: v(3, 2, 2, 2), tri(3, 3, 12), vs, ms(3), vh, ch(3), e(2), ce(2) !< Work.
   integer(I4P)           :: na, np, i, j, nt, nover(2), mode, nrec, ntot(2), nu, nunt(2) !< Counters.
   integer(I8P)           :: c0, c1, cr !< Clock.
   real(R8P)              :: tm(2) !< Times.
   logical                :: okc, allc !< Closed.

   na = 64 ; np = 32
   call sphere_vertices(na, np, 1._R8P + 15._R8P / 16._R8P, xi) ; call sphere_vertices(na, np, 2._R8P, xo)
   x = grid_1d(-2.6_R8P, 2.6_R8P, nb - 1_I4P)
   e = 0._R8P ; ce = 0._R8P ; ntot = 0 ; nrec = 0 ; allc = .true. ; nunt = 0
   t_surface = 0._R8P ; t_cut = 0._R8P
   do mode = 1, 2
      call system_clock(c0, cr)
      do j = 1, np, 3 ! every third row, the poles included (j = 1: collapsed faces)
         do i = 1, na, 2
            v(:, 1, 1, 1) = xi(:, i, j) ; v(:, 2, 1, 1) = xi(:, i + 1, j)
            v(:, 1, 2, 1) = xi(:, i, j + 1) ; v(:, 2, 2, 1) = xi(:, i + 1, j + 1)
            v(:, 1, 1, 2) = xo(:, i, j) ; v(:, 2, 1, 2) = xo(:, i + 1, j)
            v(:, 1, 2, 2) = xo(:, i, j + 1) ; v(:, 2, 2, 2) = xo(:, i + 1, j + 1)
            call hex_tris(v, tri, nt)
            call receptor_3d(tri, nt, x, mode, vs, ms, nover(mode), okc, vh, ch, nu)
            nunt(mode) = nunt(mode) + nu
            allc = allc .and. okc
            e(mode) = max(e(mode), abs(vs - vh) / vh)
            ce(mode) = max(ce(mode), norm2(ms / vs - ch) / 0.0625_R8P)
            ntot(mode) = ntot(mode) + nover(mode)
            if (mode == 1) nrec = nrec + 1
         enddo
      enddo
      call system_clock(c1)
      tm(mode) = real(c1 - c0, R8P) / cr
   enddo
   print '(A,I0,A,I0,A,2F7.2,A)', '  background ', nb - 1, '^3: ', nrec, ' receptors (outer layer, sampled), ', &
         real(ntot, R8P) / nrec, ' donor cells each (complement, direct)'
   print '(A,2ES10.2,A,2ES10.2,A,2F9.1,A,L1)', '    |sum - |V_c|| / |V_c| (complement, direct) ', e, &
         '; centroid / dr ', ce, '; ms per receptor ', 1e3_R8P * tm / nrec, '; all closed ', allc
   print '(A,2I5)', '    donor cells with pieces = -1 (untrusted; their alpha used) (complement, direct): ', nunt
   print '(A,2F7.3,A)', '    of which: building the receptor surfaces, cut_cells_3d ', 1e3_R8P * [t_surface, t_cut] / (2 * nrec), &
         ' ms per receptor (second surface of each call only)'
   ok = ok .and. allc .and. maxval(e) <= 1e-12_R8P .and. maxval(ce) <= 1e-12_R8P .and. ntot(1) == ntot(2)
   endsubroutine weights_3d
endprogram fossil_val_overset_conservative_proto
