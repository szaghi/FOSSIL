!< Tutorial 10, step 2: the 2-D cut cells of the dragon section (curve%cut_cells): fluid fractions, apertures, pieces.
!< The figures (outside the regions): the block shaded by the fluid fraction with the cells of two pieces outlined, and
!< a zoom on the first of them: its pieces' fluid edge intervals and centroids, one colour per piece.
program tut10_b_cells2d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_box, fig_label, fig_close
use vecfor,     only : vector_R8P
use fossil, only : curve_object, cut_cells_2d_t, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)   :: section
type(cut_cells_2d_t) :: cut
real(R8P)            :: x(49), h
integer(I4P)         :: i, j, k, c(2), p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]                  ! the same 48 x 48 block as step 1
call section%cut_cells(x=x, y=x, cut=cut)
print '(A,I0,A,I0,A,I0)', 'cells: fluid ', count(cut%alpha == 1._R8P), ', solid ', count(cut%alpha == 0._R8P), &
                          ', cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P)
print '(A,F10.7,A,F10.7)', 'solid area from alpha ', sum(1._R8P - cut%alpha) * h**2, ', area of the curve ', &
                           section%get_area()
i = 24                                                ! the -x edges of column 24 lie on the line of step 1
print '(A,F9.6)', 'solid length of the line x = -0.025 from the apertures: ', sum(1._R8P - cut%aperture(1, i, :)) * h
c = maxloc(cut%wet_length)                            ! the cell with the longest wall
print '(A,I0,A,I0,A,F7.4,A,4F7.4)', 'cell (', c(1), ',', c(2), '): alpha ', cut%alpha(c(1), c(2)), &
                                    ', apertures -x +x -y +y', cut%aperture(:, c(1), c(2))
print '(A,F7.4,A,2F8.4)', '  wall length / h ', cut%wet_length(c(1), c(2)) / h, ', normal integral / h ', &
                          cut%wet_normal(:, c(1), c(2)) / h
print '(A,2F8.4)', '  aperture differences (+x - -x, +y - -y)    ', &
                   cut%aperture(2, c(1), c(2)) - cut%aperture(1, c(1), c(2)), &
                   cut%aperture(4, c(1), c(2)) - cut%aperture(3, c(1), c(2))
print '(A,I0,A,I0)', 'pieces in the block ', cut%np, ', cells with two or more ', count(cut%pieces > 1)
do j = 1, 48 ; do i = 1, 48
   if (cut%pieces(i, j) < 2) cycle
   p = cut%first_piece(i, j)                          ! its pieces are p, p + 1, ...
   print '(A,I0,A,I0,A,F7.4)', '  cell (', i, ',', j, ') alpha ', cut%alpha(i, j)
   do k = p, p + cut%pieces(i, j) - 1
      print '(A,I0,A,F7.4,A,4F7.4)', '    piece ', k, ': area / h^2 ', cut%piece_area(k) / h**2, &
                                     ', apertures', cut%piece_aperture(:, k)
   enddo
enddo ; enddo
call figures
contains
   subroutine figures
   !< The block by alpha, and a zoom on the first cell with two pieces.
   type(vector_R8P), allocatable :: a(:), b(:)
   type(vector_R8P)              :: centre
   real(R8P)                     :: lo(2), hi(2), r
   integer(I4P)                  :: u, n, q, e, ci(2)
   character(80)                 :: cap
   character(5), parameter       :: STYLE(2) = ['probe', 'patch']
   real(R8P),    parameter       :: W0 = 0.05_R8P * 0.025_R8P
   real(R8P)                     :: w

   call wall_segments(a, b, n)
   call fig_open(u, 'tut10-cells2d', 'cut_cells: fluid fraction (cream 0, sepia 1), the wall, cells with 2 pieces', &
                 0._R8P, 90._R8P)
   lo = section%get_bmin() ; hi = section%get_bmax()
   centre = vec(0.5_R8P * (lo(1) + hi(1)), 0.5_R8P * (lo(2) + hi(2)), 0._R8P)
   r = max(0.5_R8P * (hi(1) - lo(1)), 0.3_R8P * (hi(2) - lo(2))) + 0.03_R8P
   call fig_zoom(u, centre, r)
   call cells(u, 1, 48, 1, 48)
   call fig_segments(u, 'hole', a(1:n), b(1:n))
   do j = 1, 48 ; do i = 1, 48
      if (cut%pieces(i, j) < 2) cycle
      call fig_box(u, 'accent', vec(x(i), x(j), 2.e-3_R8P), vec(x(i + 1), x(j + 1), 2.e-3_R8P))
      call fig_label(u, 'accent', vec(x(i) - 0.035_R8P, x(j), 2.e-3_R8P), '2 pieces')
   enddo ; enddo
   call fig_close(u)

   ci = findloc(cut%pieces, 2)
   p = cut%first_piece(ci(1), ci(2))
   write(cap, '(A,I0,A,I0,A)') 'cell (', ci(1), ',', ci(2), '): two pieces, their fluid edge intervals and centroids'
   call fig_open(u, 'tut10-pieces2d', trim(cap), 0._R8P, 90._R8P)
   call fig_zoom(u, vec(0.5_R8P * (x(ci(1)) + x(ci(1) + 1)), 0.5_R8P * (x(ci(2)) + x(ci(2) + 1)), 0._R8P), 0.8_R8P * h)
   call cells(u, ci(1) - 3, ci(1) + 3, ci(2) - 4, ci(2) + 4)
   call grid(u, ci(1) - 3, ci(1) + 4, ci(2) - 4, ci(2) + 5)
   call fig_segments(u, 'hole', a(1:n), b(1:n))
   do q = p, p + 1                                     ! each piece: its fluid edge intervals (bands inside the cell)
      do k = cut%piece_iv_start(q), cut%piece_iv_start(q + 1) - 1
         e = cut%piece_iv_edge(k)
         if (e <= 49 * 48) then                        ! vertical edge (i, j): u = x(i), v in [lo, hi]
            i = mod(e - 1, 49) + 1 ; w = merge(W0, -W0, i == cut%piece_cell(1, q))
            call fig_quad(u, STYLE(q - p + 1), vec(x(i), cut%piece_iv_lo(k), 3.e-3_R8P), &
                          vec(x(i), cut%piece_iv_hi(k), 3.e-3_R8P), vec(x(i) + w, cut%piece_iv_hi(k), 3.e-3_R8P), &
                          vec(x(i) + w, cut%piece_iv_lo(k), 3.e-3_R8P))
         else                                          ! horizontal edge (i, j): v = x(j), u in [lo, hi]
            e = e - 49 * 48 ; j = (e - 1) / 48 + 1 ; w = merge(W0, -W0, j == cut%piece_cell(2, q))
            call fig_quad(u, STYLE(q - p + 1), vec(cut%piece_iv_lo(k), x(j), 3.e-3_R8P), &
                          vec(cut%piece_iv_hi(k), x(j), 3.e-3_R8P), vec(cut%piece_iv_hi(k), x(j) + w, 3.e-3_R8P), &
                          vec(cut%piece_iv_lo(k), x(j) + w, 3.e-3_R8P))
         endif
      enddo
      call fig_points(u, STYLE(q - p + 1), [vec(cut%piece_centroid(1, q), cut%piece_centroid(2, q), 4.e-3_R8P)])
      write(cap, '(A,I0,A,F6.3)') 'piece ', q - p + 1, ': ', cut%piece_area(q) / h**2
      call fig_label(u, STYLE(q - p + 1), vec(cut%piece_centroid(1, q) + 0.004_R8P, cut%piece_centroid(2, q), &
                     4.e-3_R8P), trim(cap))
   enddo
   call fig_close(u)
   endsubroutine figures

   subroutine cells(u, i0, i1, j0, j1)
   !< The cells i0..i1, j0..j1 shaded by their fluid fraction.
   integer(I4P), intent(in) :: u              !< Scene.
   integer(I4P), intent(in) :: i0, i1, j0, j1 !< Range.
   integer(I4P)             :: ii, jj

   do jj = max(1, j0), min(48, j1) ; do ii = max(1, i0), min(48, i1)
      call fig_quad(u, 'frac', vec(x(ii), x(jj), 0._R8P), vec(x(ii + 1), x(jj), 0._R8P), &
                    vec(x(ii + 1), x(jj + 1), 0._R8P), vec(x(ii), x(jj + 1), 0._R8P), value=cut%alpha(ii, jj))
   enddo ; enddo
   endsubroutine cells

   subroutine grid(u, i0, i1, j0, j1)
   !< The grid lines x(i0..i1), x(j0..j1) over that range.
   integer(I4P), intent(in) :: u              !< Scene.
   integer(I4P), intent(in) :: i0, i1, j0, j1 !< Range of lines.
   integer(I4P)             :: k

   call fig_segments(u, 'box', [(vec(x(k), x(j0), 1.e-3_R8P), k = i0, i1), (vec(x(i0), x(k), 1.e-3_R8P), k = j0, j1)], &
                     [(vec(x(k), x(j1), 1.e-3_R8P), k = i0, i1), (vec(x(i1), x(k), 1.e-3_R8P), k = j0, j1)])
   endsubroutine grid

   subroutine wall_segments(a, b, n)
   !< The segments of the section, just above the cells.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Ends.
   integer(I4P),                  intent(out) :: n          !< Segments.
   real(R8P)                                  :: uv(2)
   integer(I4P)                               :: k

   allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; n = 0
   do k = 1, section%get_vertices_number()
      if (section%get_next(k) == 0_I4P) cycle
      n = n + 1
      uv = section%get_vertex(k) ; a(n) = vec(uv(1), uv(2), 2.e-3_R8P)
      uv = section%get_vertex(section%get_next(k)) ; b(n) = vec(uv(1), uv(2), 2.e-3_R8P)
   enddo
   endsubroutine wall_segments

   function vec(p, q, z) result(v)
   !< The point (p, q, z).
   real(R8P), intent(in) :: p, q, z !< Coordinates.
   type(vector_R8P)      :: v       !< Point.

   v%x = p ; v%y = q ; v%z = z
   endfunction vec
endprogram tut10_b_cells2d
