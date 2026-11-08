!run anim-ref-cut-cells-plate anim_ref_cut_cells_plate
!fig anim-ref-cut-cells-plate
!< Animation of the reference page "Cut-cell geometry": a plate 0.5 long and 0.03 thick (half a cell) translating and
!< turning through 16 x 10 cells of 1/16; each frame recomputes the exact cut cells (curve%cut_cells). The cut cells are
!< shaded by their fluid fraction; the cells the plate splits into two fluid pieces show the two piece centroids in gold,
!< the cells holding a fluid piece below 0.1 of the cell are outlined in red. The caption gives the counts and the
!< covering sum (solid area of the cells minus the plate area).
program anim_ref_cut_cells_plate
use fossil,     only : curve_object, cut_cells_2d_t
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter :: NX = 16, NY = 10, NF = 24
real(R8P),    parameter :: PI = acos(-1._R8P), H = 1._R8P / 16._R8P, LIFT = 1.e-3_R8P, HL = 0.25_R8P, HT = 0.015_R8P
type(curve_object)      :: plate
type(cut_cells_2d_t)    :: cut
type(vector_R8P)        :: a(4 * NX * NY), b(4 * NX * NY), c(4 * NX * NY), d(4 * NX * NY)
real(R8P)               :: x(NX + 1), y(NY + 1), u(4), v(4), cx, cy, th, t
integer(I4P)            :: i, j, k, f, m, n, un, p
logical                 :: small(NX, NY)
character(100)          :: caption

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
call fig_open(un, 'anim-ref-cut-cells-plate', '', -90._R8P, 90._R8P, width=720_I4P, height=500_I4P)
call fig_anim(un, 0.45_R8P)
m = 0                                                          ! the grid, in every frame
do j = 1, NY + 1
   m = m + 1 ; a(m) = pt(x(1), y(j), 0._R8P) ; b(m) = pt(x(NX + 1), y(j), 0._R8P)
enddo
do i = 1, NX + 1
   m = m + 1 ; a(m) = pt(x(i), y(1), 0._R8P) ; b(m) = pt(x(i), y(NY + 1), 0._R8P)
enddo
call fig_segments(un, 'box', a(1:m), b(1:m))
do f = 0, NF - 1
   t = real(f, R8P) / (NF - 1)
   cx = 0.3_R8P + 0.4_R8P * t ; cy = 0.3_R8P + 0.03_R8P * sin(PI * t) ; th = 0.15_R8P + 0.6_R8P * PI * t
   do k = 1, 4                                                ! the plate corners, counter-clockwise
      u(k) = cx + merge(HL, -HL, k == 2 .or. k == 3) * cos(th) - merge(HT, -HT, k >= 3) * sin(th)
      v(k) = cy + merge(HL, -HL, k == 2 .or. k == 3) * sin(th) + merge(HT, -HT, k >= 3) * cos(th)
   enddo
   call plate%adopt_loops(u=u, v=v, loop_sizes=[4_I4P])
   call plate%cut_cells(x=x, y=y, cut=cut)
   small = .false.
   do p = 1, cut%np                                           ! cells holding a piece below 0.1 of the cell
      if (cut%piece_area(p) > 0._R8P .and. cut%piece_area(p) < 0.1_R8P * H**2) &
         small(cut%piece_cell(1, p), cut%piece_cell(2, p)) = .true.
   enddo
   write(caption, '(A,I2,A,I2,A,I2,A,I2,A,ES9.2)') 'frame ', f, ': cut cells ', &
         count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P), ', two pieces ', count(cut%pieces == 2), ', small ', &
         count(small), ', covering sum ', &
         sum(1._R8P - cut%alpha) * H**2 - 4 * HL * HT
   print '(A)', trim(caption)
   call fig_frame(un, caption=trim(caption))
   do j = 1, NY
      do i = 1, NX
         if (cut%alpha(i, j) >= 1._R8P) cycle
         call fig_quad(un, 'frac', pt(x(i), y(j), 0._R8P), pt(x(i + 1), y(j), 0._R8P), pt(x(i + 1), y(j + 1), 0._R8P), &
                       pt(x(i), y(j + 1), 0._R8P), value=cut%alpha(i, j))
      enddo
   enddo
   call fig_quad(un, 'body', pt(u(1), v(1), LIFT), pt(u(2), v(2), LIFT), pt(u(3), v(3), LIFT), pt(u(4), v(4), LIFT))
   call outline(small, n) ; if (n > 0) call fig_segments(un, 'hole', c(1:n), d(1:n))
   n = 0                                                      ! the two centroids of every split cell
   do p = 1, cut%np
      if (cut%pieces(cut%piece_cell(1, p), cut%piece_cell(2, p)) < 2) cycle
      n = n + 1 ; c(n) = pt(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 3 * LIFT)
   enddo
   if (n > 0) call fig_points(un, 'foot', c(1:n))
enddo
call fig_close(un)
contains
   function pt(px, py, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: px, py, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(px, py, dz)
   endfunction pt

   subroutine outline(mask, nn)
   !< The edges of the cells in mask, into c, d.
   logical,      intent(in)  :: mask(:,:) !< Cells.
   integer(I4P), intent(out) :: nn        !< Segments.
   integer(I4P)              :: ii, jj    !< Cell.

   nn = 0
   do jj = 1, NY
      do ii = 1, NX
         if (.not. mask(ii, jj)) cycle
         c(nn + 1) = pt(x(ii), y(jj), 2 * LIFT)         ; d(nn + 1) = pt(x(ii + 1), y(jj), 2 * LIFT)
         c(nn + 2) = pt(x(ii + 1), y(jj), 2 * LIFT)     ; d(nn + 2) = pt(x(ii + 1), y(jj + 1), 2 * LIFT)
         c(nn + 3) = pt(x(ii + 1), y(jj + 1), 2 * LIFT) ; d(nn + 3) = pt(x(ii), y(jj + 1), 2 * LIFT)
         c(nn + 4) = pt(x(ii), y(jj + 1), 2 * LIFT)     ; d(nn + 4) = pt(x(ii), y(jj), 2 * LIFT)
         nn = nn + 4
      enddo
   enddo
   endsubroutine outline
endprogram anim_ref_cut_cells_plate
