!run anim-cut-solver anim_cut_solver
!fig anim-cut-solver
!< Animation of the cut-solver reference page: an ellipse (semi-axes 0.32 and 0.14) translating and turning through
!< 24 x 12 cells of side 1/12. Each frame recomputes the exact cut cells (curve%cut_cells), their links
!< (cut_links_2d), the state-redistribution neighbourhoods (srd_neighbourhoods) and the merge targets
!< (merge_targets, tau = 1/2). Cut cells are shaded by their fluid fraction; every small piece (teal) points to its
!< merge target (yellow segment and dot); the cells of the smallest piece's SRD neighbourhood are outlined in ember.
program anim_cut_solver
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, merge_t, &
                       merge_targets, MERGE_NONE
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter :: NX = 24, NY = 12, NS = 64, NF = 16
real(R8P),    parameter :: PI = acos(-1._R8P), H = 1._R8P / 12._R8P, LIFT = 1.e-3_R8P
type(curve_object)      :: body
type(cut_cells_2d_t)    :: cut
type(cut_links_t)       :: links
type(srd_t)             :: srd
type(merge_t)           :: merge
type(vector_R8P)        :: a(4 * NX * NY), b(4 * NX * NY)
real(R8P)               :: x(NX + 1), y(NY + 1), u(NS), v(NS), cx, cy, th, t
integer(I4P)            :: i, j, k, f, m, p, q, un, nsmall
character(96)           :: caption

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
call fig_open(un, 'anim-cut-solver', '', -90._R8P, 90._R8P, width=720_I4P, height=420_I4P)
call fig_anim(un, 0.5_R8P)
m = 0                                                          ! the grid, in every frame
do j = 1, NY + 1
   m = m + 1 ; a(m) = pt(x(1), y(j), 0._R8P) ; b(m) = pt(x(NX + 1), y(j), 0._R8P)
enddo
do i = 1, NX + 1
   m = m + 1 ; a(m) = pt(x(i), y(1), 0._R8P) ; b(m) = pt(x(i), y(NY + 1), 0._R8P)
enddo
call fig_segments(un, 'box', a(1:m), b(1:m))

do f = 0, NF
   t = real(f, R8P) / NF
   cx = 0.45_R8P + 1.1_R8P * t ; cy = 0.5_R8P + 0.07_R8P * sin(2 * PI * t) ; th = PI / 3 * t
   do k = 1, NS                                                ! the ellipse, counter-clockwise
      u(k) = cx + 0.32_R8P * cos(2 * PI * (k - 1) / NS) * cos(th) - 0.14_R8P * sin(2 * PI * (k - 1) / NS) * sin(th)
      v(k) = cy + 0.32_R8P * cos(2 * PI * (k - 1) / NS) * sin(th) + 0.14_R8P * sin(2 * PI * (k - 1) / NS) * cos(th)
   enddo
   call body%adopt_loops(u=u, v=v, loop_sizes=[NS])
   call body%cut_cells(x=x, y=y, cut=cut)
   call cut_links_2d(cut, links)
   call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, y, srd)
   call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, y, merge)
   nsmall = count(merge%status /= MERGE_NONE)
   p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)
   write(caption, '(A,I2,A,I3,A,I2,A,ES8.2,A,I0,A)') 'frame ', f, ': pieces ', cut%np, ', small ', nsmall, &
         ', smallest ', cut%piece_area(p) / H**2, ' of its cell, ', srd%start(p + 1) - srd%start(p), ' in its SRD'
   print '(A)', trim(caption)
   call fig_frame(un, caption=trim(caption))
   do j = 1, NY                                                ! cut cells by fluid fraction
      do i = 1, NX
         if (cut%alpha(i, j) <= 0._R8P .or. cut%alpha(i, j) >= 1._R8P) cycle
         call fig_quad(un, 'frac', pt(x(i), y(j), 0._R8P), pt(x(i + 1), y(j), 0._R8P), pt(x(i + 1), y(j + 1), 0._R8P), &
                       pt(x(i), y(j + 1), 0._R8P), value=cut%alpha(i, j))
      enddo
   enddo
   call fig_segments(un, 'body', [(pt(u(k), v(k), LIFT), k = 1, NS)], [(pt(u(mod(k, NS) + 1), v(mod(k, NS) + 1), LIFT), &
                     k = 1, NS)])
   m = 0                                                       ! the smallest piece's SRD neighbourhood
   do k = srd%start(p), srd%start(p + 1) - 1
      q = srd%member(k) ; i = cut%piece_cell(1, q) ; j = cut%piece_cell(2, q)
      a(m + 1) = pt(x(i), y(j), 2 * LIFT)         ; b(m + 1) = pt(x(i + 1), y(j), 2 * LIFT)
      a(m + 2) = pt(x(i + 1), y(j), 2 * LIFT)     ; b(m + 2) = pt(x(i + 1), y(j + 1), 2 * LIFT)
      a(m + 3) = pt(x(i + 1), y(j + 1), 2 * LIFT) ; b(m + 3) = pt(x(i), y(j + 1), 2 * LIFT)
      a(m + 4) = pt(x(i), y(j + 1), 2 * LIFT)     ; b(m + 4) = pt(x(i), y(j), 2 * LIFT)
      m = m + 4
   enddo
   call fig_segments(un, 'patch', a(1:m), b(1:m))
   m = 0                                                       ! small pieces and their merge targets
   do q = 1, cut%np
      if (merge%status(q) == MERGE_NONE .or. merge%target(q) < 1) cycle
      m = m + 1
      a(m) = pt(cut%piece_centroid(1, q), cut%piece_centroid(2, q), 3 * LIFT)
      b(m) = pt(cut%piece_centroid(1, merge%target(q)), cut%piece_centroid(2, merge%target(q)), 3 * LIFT)
   enddo
   call fig_segments(un, 'accent', a(1:m), b(1:m))
   call fig_points(un, 'foot', b(1:m))
   call fig_points(un, 'probe', a(1:m))
enddo
call fig_close(un)
contains
   function pt(px, py, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: px, py, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(px, py, dz)
   endfunction pt
endprogram anim_cut_solver
