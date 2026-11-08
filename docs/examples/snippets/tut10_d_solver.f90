!< Tutorial 10, step 4: what a cut-cell solver needs beyond the fractions, on the 2-D cut cells of the dragon section:
!< the links between pieces (cut_links_2d), the state-redistribution neighbourhoods (srd_neighbourhoods) and the merge
!< targets (merge_targets) of the small pieces; then the links of the 3-D cut cells (cut_cells_3d(links=)).
!< The figure (outside the regions): a zoom on the tail: cut cells by fluid fraction, the links of the pieces in view,
!< the small pieces and their merge targets, the SRD neighbourhood of the smallest piece.
program tut10_d_solver
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_label, fig_close
use vecfor,     only : vector_R8P
use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, cut_links_t, &
                   cut_links_2d, srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, merge_t, &
                   merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)       :: section
type(surface_stl_object) :: surface
type(cut_cells_2d_t)     :: cut
type(cut_cells_3d_t)     :: cut3
type(cut_links_t)        :: links, links3
type(srd_t)              :: srd
type(merge_t)            :: mt
real(R8P)                :: x(49), h
integer(I4P)             :: i, k, p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]
call section%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)                        ! piece p meets links%piece(k) over face links%face(k)
print '(A,I0,A,I0)', 'pieces ', cut%np, ', links (both directions) ', size(links%piece)
p = cut%first_piece(22, 20) + 1                      ! the small second piece of cell (22, 20), step 2
do k = links%start(p), links%start(p + 1) - 1
   print '(A,I0,A,I0,A,2I3,A,I0,A,F7.4)', '  piece ', p, ' -> ', links%piece(k), ' in cell', &
         cut%piece_cell(:, links%piece(k)), ' over face ', links%face(k), ', length / h ', links%area(k) / h
enddo
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, x, srd)
print '(A,I0,A,I0,A,I0)', 'SRD: regular ', count(srd%status == SRD_REGULAR), ', small ', &
                          count(srd%status == SRD_COMPLETE), ', incomplete ', count(srd%status == SRD_INCOMPLETE)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)          ! the smallest piece with area
print '(A,I0,A,2I3,A,ES9.2,A,I0,A)', '  smallest piece ', p, ' in cell', cut%piece_cell(:, p), ', area / h^2 ', &
                                     cut%piece_area(p) / h**2, ', neighbourhood of ', srd%start(p + 1) - srd%start(p), &
                                     ' pieces'
print '(A,F7.4,A,F7.4)', '  neighbourhood area / h^2 ', &
                         sum(cut%piece_area(srd%member(srd%start(p):srd%start(p + 1) - 1))) / h**2, ' >= tau = ', srd%tau
call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, x, mt)
print '(A,I0,A,I0,A,I0)', 'merge: served ', count(mt%status /= MERGE_NONE), ', direct ', count(mt%status == MERGE_OK), &
                          ', through a chain ', count(mt%status == MERGE_CHAIN)
print '(A,2I3,A,I0)', '  smallest piece -> cell', cut%piece_cell(:, mt%final_target(p)), &
                      ', normal-directed targets ', mt%mix_start(p + 1) - mt%mix_start(p)
call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call cut_cells_3d(surface=surface, x=x(1:49:2), y=x(1:49:2), z=x(1:49:2), cut=cut3, links=links3)
print '(A,I0,A,I0,A,I0)', '3-D: pieces ', cut3%np, ', links ', size(links3%piece), ', unmatched regions ', &
                          links3%unmatched
call figure
contains
   subroutine figure
   !< Zoom on the smallest piece: cells by alpha, links, small pieces and targets, its SRD neighbourhood.
   type(vector_R8P), allocatable :: a(:), b(:)
   type(vector_R8P)              :: c
   real(R8P)                     :: uv(2)
   integer(I4P)                  :: u, n, q, ii, jj, m
   character(80)                 :: cap

   write(cap, '(A,ES8.2,A)') 'smallest piece ', cut%piece_area(p) / h**2, &
                             ': SRD cells ember, links gold, small teal to target'
   call fig_open(u, 'tut10-solver', trim(cap), 0._R8P, 90._R8P)
   c = vec(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 0._R8P)
   call fig_zoom(u, c, 2.6_R8P * h)
   do jj = 1, 48 ; do ii = 1, 48
      if (abs(x(ii) - c%x) > 4 * h .or. abs(x(jj) - c%y) > 6 * h) cycle
      call fig_quad(u, 'frac', vec(x(ii), x(jj), 0._R8P), vec(x(ii + 1), x(jj), 0._R8P), &
                    vec(x(ii + 1), x(jj + 1), 0._R8P), vec(x(ii), x(jj + 1), 0._R8P), value=cut%alpha(ii, jj))
   enddo ; enddo
   call fig_segments(u, 'box', [(vec(x(ii), x(1), 1.e-3_R8P), ii = 1, 49), (vec(x(1), x(ii), 1.e-3_R8P), ii = 1, 49)], &
                     [(vec(x(ii), x(49), 1.e-3_R8P), ii = 1, 49), (vec(x(49), x(ii), 1.e-3_R8P), ii = 1, 49)])
   allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; n = 0
   do q = 1, section%get_vertices_number()                                         ! the wall
      if (section%get_next(q) == 0_I4P) cycle
      n = n + 1
      uv = section%get_vertex(q) ; a(n) = vec(uv(1), uv(2), 2.e-3_R8P)
      uv = section%get_vertex(section%get_next(q)) ; b(n) = vec(uv(1), uv(2), 2.e-3_R8P)
   enddo
   call fig_segments(u, 'hole', a(1:n), b(1:n))
   deallocate(a, b) ; allocate(a(4 * 25), b(4 * 25)) ; m = 0                        ! the SRD cells of p
   do k = srd%start(p), srd%start(p + 1) - 1
      ii = cut%piece_cell(1, srd%member(k)) ; jj = cut%piece_cell(2, srd%member(k))
      a(m + 1:m + 4) = [vec(x(ii), x(jj), 3.e-3_R8P), vec(x(ii + 1), x(jj), 3.e-3_R8P), &
                        vec(x(ii + 1), x(jj + 1), 3.e-3_R8P), vec(x(ii), x(jj + 1), 3.e-3_R8P)]
      b(m + 1:m + 4) = [a(m + 2:m + 4), a(m + 1)]
      m = m + 4
   enddo
   call fig_segments(u, 'patch', a(1:m), b(1:m))
   deallocate(a, b) ; allocate(a(size(links%piece)), b(size(links%piece))) ; m = 0   ! the links of cut pieces
   do q = 1, cut%np
      if (cut%alpha(cut%piece_cell(1, q), cut%piece_cell(2, q)) >= 1._R8P) cycle
      do k = links%start(q), links%start(q + 1) - 1
         m = m + 1
         a(m) = vec(cut%piece_centroid(1, q), cut%piece_centroid(2, q), 4.e-3_R8P)
         b(m) = vec(cut%piece_centroid(1, links%piece(k)), cut%piece_centroid(2, links%piece(k)), 4.e-3_R8P)
      enddo
   enddo
   call fig_segments(u, 'accent', a(1:m), b(1:m))
   deallocate(a, b) ; allocate(a(cut%np), b(cut%np)) ; m = 0                        ! small pieces -> final target
   do q = 1, cut%np
      if (mt%status(q) == MERGE_NONE .or. mt%final_target(q) < 1) cycle
      m = m + 1
      a(m) = vec(cut%piece_centroid(1, q), cut%piece_centroid(2, q), 5.e-3_R8P)
      b(m) = vec(cut%piece_centroid(1, mt%final_target(q)), cut%piece_centroid(2, mt%final_target(q)), 5.e-3_R8P)
   enddo
   call fig_segments(u, 'probe', a(1:m), b(1:m))
   call fig_points(u, 'probe', a(1:m))
   call fig_label(u, 'patch', vec(cut%piece_centroid(1, p) + 0.003_R8P, cut%piece_centroid(2, p), 6.e-3_R8P), &
                  'smallest')
   call fig_close(u)
   endsubroutine figure

   function vec(p, q, z) result(v)
   !< The point (p, q, z).
   real(R8P), intent(in) :: p, q, z !< Coordinates.
   type(vector_R8P)      :: v       !< Point.

   v%x = p ; v%y = q ; v%z = z
   endfunction vec
endprogram tut10_d_solver
