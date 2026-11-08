!< Figures of the reference page "Overset assembly kernels", on the cylinder of the cookbook recipe (an O-grid of
!< 32 x 12 cells to r = 2 over a Cartesian background of spacing 0.15):
!< - mask: a close-up of the overlap, every vertex coloured by its iota (1 field, 0 hole, -k receptor from block k),
!<   the background values written next to the vertices, one donor cell of each block filled;
!< - selection: the explicit rule (offset holes) and the implicit rule (no holes) side by side, from the same candidates;
!< - fringes: a background cut short on the left (the O-grid leaves it): the two INNER levels around the holes, the
!<   two OUTER levels on the O-grid, the orphans with no candidate, and the segment of Chesshire's overlap width d;
!< - connectivity: three receptors of the O-grid's outer ring, each with its 3 x 3 background stencil and its record.
program fig_ref_overset_mask
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                       hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, &
                       SELECT_EXPLICIT, SELECT_IMPLICIT, fringe_levels_2d, overlap_width_2d, connectivity_t, &
                       build_connectivity, FRINGE_INNER_1, FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2, &
                       ORPHAN_NO_CANDIDATE
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_label, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter               :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curvilinear_block_2d_t)       :: blocks(2), short(2)
type(curve_object)                 :: wall
type(block_logical_t), allocatable :: holes(:), holes_off(:)
type(donor_set_t),     allocatable :: donors(:), sdonors(:)
type(overset_mask_t),  allocatable :: masks(:), expl(:), impl(:), smasks(:)
type(connectivity_t)               :: conn
real(R8P),             allocatable :: width(:,:)
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                       :: i, j, uu, it_e, it_i
character(160)                     :: caption

do j = 1, 13
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
do j = 1, 41
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call make_blocks(blocks, x, y)
call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
call hole_flags_2d(blocks, wall, holes_off, offset=0.3_R8P, cuts=[.false., .true.])
call donor_candidates_2d(blocks, donors)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
call classify_2d(blocks, donors, expl, SELECT_EXPLICIT, 3_I4P, holes=holes_off, iterations=it_e)
call classify_2d(blocks, donors, impl, SELECT_IMPLICIT, 3_I4P, iterations=it_i)
print '(A,I0,A,I0,A,I0)', 'mask: holes ', count(masks(2)%iota == 0), ', receptors ', count(masks(1)%iota < 0), &
                          ' + ', count(masks(2)%iota < 0)
print '(A,I0,A,I0,A,I0,A,I0)', 'explicit (offset 0.3): holes ', count(expl(2)%iota == 0), ', receptors ', &
                               count(expl(1)%iota < 0), ' + ', count(expl(2)%iota < 0), ', iterations ', it_e
print '(A,I0,A,I0,A,I0)', 'implicit: receptors ', count(impl(1)%iota < 0), ' + ', count(impl(2)%iota < 0), &
                          ', iterations ', it_i
call draw_mask
call draw_selection
call draw_fringes
call draw_connectivity
contains
   subroutine make_blocks(b, bx, by)
   !< The O-grid (wall on face 3, overlap boundary on face 4) and a background of the given vertices.
   type(curvilinear_block_2d_t), intent(inout) :: b(2)     !< Blocks.
   real(R8P),                    intent(in)    :: bx(:,:), by(:,:) !< Background vertices.
   type(topo_record_t)                         :: rec

   call b(1)%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
   call b(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
   call b(1)%add_record(rec)
   call b(2)%initialize(bx, by)
   call b(1)%build_search ; call b(2)%build_search
   endsubroutine make_blocks

   subroutine draw_mask
   !< Close-up of the overlap: iota per vertex, one donor cell per block.
   type(vector_R8P) :: c
   integer(I4P)     :: k, n, cell(2)
   character(8)     :: t

   call fig_open(uu, 'ref-overset-mask', 'iota: 1 field, 0 hole, -k receptor from block k; '// &
                 'donor cells filled', -90._R8P, 90._R8P, width=640_I4P, height=560_I4P)
   call fig_zoom(uu, pt(1.25_R8P, 0.85_R8P), 0.8_R8P)
   call grid(blocks(2)%u, blocks(2)%v, 'box', 0._R8P, 0._R8P)
   call grid(u, v, 'body', 0._R8P, LIFT)
   ! the O-grid cell that donates to background vertex (27, 24), the background cell that donates to O-grid (4, 13)
   n = 27 + 41 * 23
   if (masks(2)%donor(n) <= 0) error stop 'mask: no donor'
   cell = donors(2)%cell(:, masks(2)%donor(n))
   call fig_quad(uu, 'probe', pt(u(cell(1), cell(2)), v(cell(1), cell(2)), 0.5_R8P * LIFT), &
                 pt(u(cell(1) + 1, cell(2)), v(cell(1) + 1, cell(2)), 0.5_R8P * LIFT), &
                 pt(u(cell(1) + 1, cell(2) + 1), v(cell(1) + 1, cell(2) + 1), 0.5_R8P * LIFT), &
                 pt(u(cell(1), cell(2) + 1), v(cell(1), cell(2) + 1), 0.5_R8P * LIFT))
   print '(A,2I3,A,I0,A,2I3)', 'background vertex', 27, 24, ': iota ', masks(2)%iota(n), ', donor O-grid cell', cell
   n = 4 + 33 * 12
   cell = donors(1)%cell(:, masks(1)%donor(n))
   call fig_quad(uu, 'foot', pt(x(cell(1), cell(2)), y(cell(1), cell(2)), 0.5_R8P * LIFT), &
                 pt(x(cell(1) + 1, cell(2)), y(cell(1) + 1, cell(2)), 0.5_R8P * LIFT), &
                 pt(x(cell(1) + 1, cell(2) + 1), y(cell(1) + 1, cell(2) + 1), 0.5_R8P * LIFT), &
                 pt(x(cell(1), cell(2) + 1), y(cell(1), cell(2) + 1), 0.5_R8P * LIFT))
   print '(A,2I3,A,I0,A,2I3)', 'O-grid vertex', 4, 13, ': iota ', masks(1)%iota(n), ', donor background cell', cell
   call iota_points(blocks(2)%u, blocks(2)%v, masks(2)%iota, 2._R8P * LIFT)
   call iota_points(u, v, masks(1)%iota, 3._R8P * LIFT)
   do j = 1, 41
      do i = 1, 41
         c = pt(x(i, j), y(i, j), 4._R8P * LIFT)
         if (abs(x(i, j) - 1.25_R8P) > 0.8_R8P .or. abs(y(i, j) - 0.85_R8P) > 0.8_R8P) cycle
         k = masks(2)%iota(i + 41 * (j - 1))
         write(t, '(I0)') k
         if (k == 1) then
            call fig_label(uu, 'box', c, trim(t))
         elseif (k == 0) then
            call fig_label(uu, 'hole', c, trim(t))
         else
            call fig_label(uu, 'probe', c, trim(t))
         endif
      enddo
   enddo
   do j = 1, 13
      do i = 1, 32
         k = masks(1)%iota(i + 33 * (j - 1))
         if (k /= -2) cycle
         c = pt(u(i, j), v(i, j), 4._R8P * LIFT)
         if (abs(u(i, j) - 1.25_R8P) > 0.8_R8P .or. abs(v(i, j) - 0.85_R8P) > 0.8_R8P) cycle
         call fig_label(uu, 'foot', c, '-2')
      enddo
   enddo
   call fig_close(uu)
   endsubroutine draw_mask

   subroutine iota_points(pu, pv, iota, dz)
   !< Vertices by class: field (box), hole (hole), receptor from the O-grid (probe), from the background (foot).
   real(R8P),    intent(in) :: pu(:,:), pv(:,:) !< Vertices.
   integer(I4P), intent(in) :: iota(:)          !< Classes.
   real(R8P),    intent(in) :: dz               !< Lift.
   type(vector_R8P)         :: p(size(iota))
   integer(I4P)             :: n1, a, b

   n1 = size(pu, 1)
   do b = 1, size(pu, 2)
      do a = 1, n1
         p(a + n1 * (b - 1)) = pt(pu(a, b), pv(a, b), dz)
      enddo
   enddo
   call fig_points(uu, 'box', pack(p, iota == 1))
   call fig_points(uu, 'hole', pack(p, iota == 0))
   call fig_points(uu, 'probe', pack(p, iota == -1))
   call fig_points(uu, 'foot', pack(p, iota == -2))
   endsubroutine iota_points

   subroutine draw_selection
   !< Explicit (offset holes) and implicit rules side by side.
   real(R8P), parameter :: SH = 6.6_R8P
   type(vector_R8P)     :: p(41 * 41), q(33 * 13)
   integer(I4P)         :: run

   write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'same candidates, explicit rule (left): ', count(expl(1)%iota < 0), &
        ' + ', count(expl(2)%iota < 0), ' receptors; implicit rule (right): ', count(impl(1)%iota < 0), ' + ', &
        count(impl(2)%iota < 0)
   call fig_open(uu, 'ref-overset-selection', trim(caption), -90._R8P, 90._R8P, width=900_I4P, height=470_I4P)
   do run = 1, 2
      call grid(blocks(2)%u + SH * (run - 1), blocks(2)%v, 'box', 0._R8P, 0._R8P)
      call grid(u + SH * (run - 1), v, 'body', 0._R8P, LIFT)
      p = [((pt(x(i, j) + SH * (run - 1), y(i, j), 2 * LIFT), i = 1, 41), j = 1, 41)]
      q = [((pt(u(i, j) + SH * (run - 1), v(i, j), 2 * LIFT), i = 1, 33), j = 1, 13)]
      if (run == 1) then
         call fig_points(uu, 'hole', pack(p, expl(2)%iota == 0))
         call fig_points(uu, 'probe', pack(p, expl(2)%iota < 0))
         call fig_points(uu, 'foot', pack(q, expl(1)%iota < 0))
      else
         call fig_points(uu, 'probe', pack(p, impl(2)%iota < 0))
         call fig_points(uu, 'foot', pack(q, impl(1)%iota < 0))
      endif
   enddo
   call fig_label(uu, 'accent', pt(-3._R8P, 3.25_R8P), 'explicit: offset holes, fringe, overlap boundary')
   call fig_label(uu, 'accent', pt(-3._R8P + SH, 3.25_R8P), 'implicit: smaller donor cell, no holes')
   call fig_close(uu)
   endsubroutine draw_selection

   subroutine draw_fringes
   !< A background cut short on the left: fringe levels, orphans, overlap width.
   real(R8P)                :: sx(33, 41), sy(33, 41), best, dd
   type(vector_R8P)         :: p(33 * 41), q(33 * 13), pa, pb
   integer(I4P)             :: a, b, norph

   sx = x(9:41, :) ; sy = y(9:41, :)                       ! the background from x = -1.8
   call make_blocks(short, sx, sy)
   call donor_candidates_2d(short, sdonors)
   call hole_flags_2d(short, wall, holes, cuts=[.false., .true.])
   call classify_2d(short, sdonors, smasks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   call fringe_levels_2d(short, smasks, levels=2_I4P)
   call overlap_width_2d(short, smasks, width)
   norph = count(smasks(1)%orphan == ORPHAN_NO_CANDIDATE)
   p = [((pt(sx(i, j), sy(i, j), 2 * LIFT), i = 1, 33), j = 1, 41)]
   q = [((pt(u(i, j), v(i, j), 2 * LIFT), i = 1, 33), j = 1, 13)]
   best = huge(1._R8P)                                     ! the pair that sets d, by brute force
   do a = 1, size(q)
      if (smasks(1)%iota(a) /= -2) cycle
      do b = 1, size(p)
         if (smasks(2)%iota(b) /= -1) cycle
         dd = sqrt((q(a)%x - p(b)%x) ** 2 + (q(a)%y - p(b)%y) ** 2)
         if (dd < best) then
            best = dd ; pa = q(a) ; pb = p(b)
         endif
      enddo
   enddo
   print '(A,4I5)', 'fringes INNER 1, INNER 2, OUTER 1, OUTER 2: ', count(smasks(2)%fringe == FRINGE_INNER_1), &
                    count(smasks(2)%fringe == FRINGE_INNER_2), count(smasks(1)%fringe == FRINGE_OUTER_1), &
                    count(smasks(1)%fringe == FRINGE_OUTER_2)
   print '(A,I0,A,F7.4,A,F7.4)', 'orphans (no candidate) ', norph, ', overlap width ', width(1, 2), &
                                 ', brute force ', best
   write(caption, '(A,I0,A,F6.4)') 'fringe levels 1 (teal), 2 (gold); ', norph, &
                                   ' orphans (ember); overlap width d = ', width(1, 2)
   call fig_open(uu, 'ref-overset-fringes', trim(caption), -90._R8P, 90._R8P, width=640_I4P, height=560_I4P)
   call grid(sx, sy, 'box', 0._R8P, 0._R8P)
   call grid(u, v, 'ghost', 0._R8P, LIFT)
   call fig_points(uu, 'cell', pack(p, smasks(2)%iota == 0))
   call fig_points(uu, 'probe', pack(p, smasks(2)%fringe == FRINGE_INNER_1))
   call fig_points(uu, 'foot', pack(p, smasks(2)%fringe == FRINGE_INNER_2))
   call fig_points(uu, 'probe', pack(q, smasks(1)%fringe == FRINGE_OUTER_1 .and. smasks(1)%orphan == 0))
   call fig_points(uu, 'foot', pack(q, smasks(1)%fringe == FRINGE_OUTER_2))
   call fig_points(uu, 'patch', pack(q, smasks(1)%orphan == ORPHAN_NO_CANDIDATE))
   call fig_segments(uu, 'hole', [pa + 3 * LIFT * ez_R8P], [pb + 3 * LIFT * ez_R8P])
   call fig_label(uu, 'hole', 0.5_R8P * (pa + pb) + 4 * LIFT * ez_R8P, 'd')
   call fig_label(uu, 'accent', pt(0.3_R8P, -0.1_R8P, 4 * LIFT), 'INNER')
   call fig_label(uu, 'accent', pt(1.45_R8P, 1.6_R8P, 4 * LIFT), 'OUTER')
   call fig_close(uu)
   endsubroutine draw_fringes

   subroutine draw_connectivity
   !< Three receptors of the O-grid's outer ring with their 3 x 3 background stencils and their records.
   integer(I4P), parameter       :: RV(3) = [2_I4P, 3_I4P, 4_I4P]   ! O-grid vertices (i, 13)
   integer(I4P), allocatable     :: nodes(:,:,:), sst(:), cell(:,:)
   real(R8P),    allocatable     :: w(:,:), xi(:,:)
   type(vector_R8P)              :: rp, sp(9)
   integer(I4P)                  :: k, m, n, c
   character(80)                 :: t

   call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
   allocate(cell(2, 3), xi(2, 3))
   do k = 1, 3
      n = RV(k) + 33 * 12 ; c = masks(1)%donor(n)
      cell(:, k) = donors(1)%cell(:, c) ; xi(:, k) = donors(1)%xi(:, c)
   enddo
   call blocks(2)%stencils(cell, xi, 3_I4P, nodes, w, sst)
   call fig_open(uu, 'ref-overset-connectivity', 'record: block, vertex, donor block, corner cell (filled), '// &
                 'xi, width', -90._R8P, 90._R8P, width=640_I4P, height=520_I4P)
   call fig_zoom(uu, pt(1.9_R8P, 0.75_R8P), 0.62_R8P)
   call grid(x, y, 'box', 0._R8P, 0._R8P)
   call grid(u, v, 'ghost', 0._R8P, LIFT)
   do k = 1, 3
      n = RV(k) + 33 * 12
      rp = pt(u(RV(k), 13), v(RV(k), 13), 10 * LIFT)
      call fig_quad(uu, 'cell', pt(x(cell(1, k), cell(2, k)), y(cell(1, k), cell(2, k)), 0.5_R8P * LIFT), &
                    pt(x(cell(1, k) + 1, cell(2, k)), y(cell(1, k) + 1, cell(2, k)), 0.5_R8P * LIFT), &
                    pt(x(cell(1, k) + 1, cell(2, k) + 1), y(cell(1, k) + 1, cell(2, k) + 1), 0.5_R8P * LIFT), &
                    pt(x(cell(1, k), cell(2, k) + 1), y(cell(1, k), cell(2, k) + 1), 0.5_R8P * LIFT))
      do m = 1, 9
         sp(m) = pt(x(nodes(1, m, k), nodes(2, m, k)), y(nodes(1, m, k), nodes(2, m, k)), 2 * LIFT + k * LIFT)
      enddo
      call fig_segments(uu, 'probe', [(rp, m = 1, 9)], sp)
      call fig_points(uu, 'probe', sp)
      call fig_points(uu, 'foot', [rp])
      do m = 1, conn%n
         if (conn%receptor_block(m) == 1 .and. conn%receptor_vertex(m) == n) exit
      enddo
      write(t, '(I0,1X,I0,1X,I0,2(1X,I0),2F7.3,2(1X,I0))') conn%receptor_block(m), conn%receptor_vertex(m), &
            conn%donor_block(m), conn%corner(1:2, m), conn%xi(1:2, m), conn%width(1:2, m)
      print '(A)', 'record: '//trim(t)
      call fig_label(uu, 'accent', rp + LIFT * ez_R8P, trim(t))
   enddo
   call fig_close(uu)
   endsubroutine draw_connectivity

   subroutine grid(gx, gy, style, dx, dz)
   !< The lines of a block.
   real(R8P),    intent(in) :: gx(:,:), gy(:,:) !< Vertices.
   character(*), intent(in) :: style            !< Style.
   real(R8P),    intent(in) :: dx, dz           !< Shift, lift.
   type(vector_R8P), allocatable :: a(:), b(:)
   integer(I4P)             :: n1, n2, m, p, q

   n1 = size(gx, 1) ; n2 = size(gx, 2)
   allocate(a((n1 - 1) * n2 + n1 * (n2 - 1)), b((n1 - 1) * n2 + n1 * (n2 - 1))) ; m = 0
   do q = 1, n2
      do p = 1, n1 - 1
         m = m + 1 ; a(m) = pt(gx(p, q) + dx, gy(p, q), dz) ; b(m) = pt(gx(p + 1, q) + dx, gy(p + 1, q), dz)
      enddo
   enddo
   do q = 1, n2 - 1
      do p = 1, n1
         m = m + 1 ; a(m) = pt(gx(p, q) + dx, gy(p, q), dz) ; b(m) = pt(gx(p, q + 1) + dx, gy(p, q + 1), dz)
      enddo
   enddo
   call fig_segments(uu, style, a, b)
   endsubroutine grid

   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt
endprogram fig_ref_overset_mask
