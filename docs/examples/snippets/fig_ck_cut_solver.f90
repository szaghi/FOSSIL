!< Figure of the cookbook recipe "Neighbourhoods for a small-cell cure": the disc of ck_cut_solver on its 16 x 16 cells,
!< every cut piece coloured by its area fraction, the smallest piece and the members of its state-redistribution
!< neighbourhood outlined, its merge target marked.
program fig_ck_cut_solver
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, merge_t, &
                       merge_targets
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter    :: lift = 1.e-3_R8P, h = 0.0625_R8P
type(curve_object)      :: disc
type(cut_cells_2d_t)    :: cut
type(cut_links_t)       :: links
type(srd_t)             :: srd
type(merge_t)           :: merge
type(vector_R8P)        :: a(48), b(48)
real(R8P)               :: u(48), v(48), x(17), c0(2)
integer(I4P)            :: i, j, k, p, q, uu

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(h * i, i = 0, 16)]
call disc%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, x, srd)
call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, x, merge)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)

call fig_open(uu, 'ck-cut-solver', 'area fraction of the cut cells; the smallest piece, its SRD neighbourhood and merge target', &
              -90._R8P, 90._R8P, height=520_I4P)
do j = 1, 16 ; do i = 1, 16
   if (cut%alpha(i, j) <= 0._R8P .or. cut%alpha(i, j) >= 1._R8P) cycle
   c0 = [x(i), x(j)]
   call fig_quad(uu, 'frac', pt(c0(1), c0(2)), pt(c0(1) + h, c0(2)), pt(c0(1) + h, c0(2) + h), pt(c0(1), c0(2) + h), &
                 value=cut%alpha(i, j))
enddo ; enddo
a = [(pt(u(i), v(i), lift), i = 1, 48)] ; b = [a(2:48), a(1)]
call fig_segments(uu, 'body', a, b)
do k = srd%start(p), srd%start(p + 1) - 1                               ! the neighbourhood, cell outlines
   q = srd%member(k) ; c0 = [x(cut%piece_cell(1, q)), x(cut%piece_cell(2, q))]
   call fig_segments(uu, 'accent', [pt(c0(1), c0(2), 2 * lift), pt(c0(1) + h, c0(2), 2 * lift), &
                     pt(c0(1) + h, c0(2) + h, 2 * lift), pt(c0(1), c0(2) + h, 2 * lift)], &
                     [pt(c0(1) + h, c0(2), 2 * lift), pt(c0(1) + h, c0(2) + h, 2 * lift), &
                     pt(c0(1), c0(2) + h, 2 * lift), pt(c0(1), c0(2), 2 * lift)])
enddo
call fig_points(uu, 'probe', [pt(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 3 * lift)])
call fig_label(uu, 'accent', pt(cut%piece_centroid(1, p) + 0.02_R8P, cut%piece_centroid(2, p) - 0.03_R8P, 3 * lift), 'small')
q = merge%target(p)
call fig_points(uu, 'foot', [pt(cut%piece_centroid(1, q), cut%piece_centroid(2, q), 3 * lift)])
call fig_label(uu, 'accent', pt(cut%piece_centroid(1, q) + 0.02_R8P, cut%piece_centroid(2, q) + 0.02_R8P, 3 * lift), 'target')
call fig_close(uu)
contains
   function pt(px, py, dz) result(r)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: r      !< Point.

   r = px * ex_R8P + py * ey_R8P
   if (present(dz)) r = r + dz * ez_R8P
   endfunction pt
endprogram fig_ck_cut_solver
