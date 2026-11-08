!< Figure of the reference page "2-D geometry (curves)", section Conventions: the orientation rule in the plane of
!< normal_axis = 3 (u = x, v = y). An outer square counter-clockwise, a hole clockwise, an island in the hole
!< counter-clockwise, and an open C-shaped loop beside them. Chevrons show the direction of every segment; the solid
!< (left of every segment) is filled teal, cell by cell, by the fast winding number inside test, which closes the open
!< loop by the chord between its ends.
program fig_ref_2d_geometry_conventions
use fossil,     only : curve_object, SIGN_FAST_WINDING_NUMBER, CURVE_OPEN_LOOP
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter       :: H = 0.1_R8P, S = 0.09_R8P
character(6), parameter       :: ROLE(4) = ['outer ', 'hole  ', 'island', 'open  ']
type(curve_object)            :: c
real(R8P)                     :: u(18), v(18), uv(2), e(2), m(2), nn(2), dl
type(vector_R8P), allocatable :: s1(:), s2(:), c1(:), c2(:)
integer(I4P)                  :: i, j, l, n, k, w, first, last, defect, inside
logical                       :: closed

u(1:4) = [0._R8P, 4._R8P, 4._R8P, 0._R8P]       ; v(1:4) = [0._R8P, 0._R8P, 3._R8P, 3._R8P]       ! outer, CCW
u(5:8) = [0.5_R8P, 0.5_R8P, 3.5_R8P, 3.5_R8P]   ; v(5:8) = [0.5_R8P, 2.5_R8P, 2.5_R8P, 0.5_R8P]   ! hole, CW
u(9:11) = [1.5_R8P, 2.5_R8P, 2._R8P]            ; v(9:11) = [1._R8P, 1._R8P, 2._R8P]              ! island, CCW
u(12:18) = [6.2_R8P, 5.2_R8P, 4.8_R8P, 4.8_R8P, 5.2_R8P, 6.2_R8P, 6.2_R8P]                          ! open C, CCW
v(12:18) = [2.6_R8P, 2.6_R8P, 2.2_R8P, 0.8_R8P, 0.4_R8P, 0.4_R8P, 1.0_R8P]
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 4_I4P, 3_I4P, 7_I4P], closed=[.true., .true., .true., .false.])
defect = c%check_validity(loop=l)
print '(A,I0,A,L1,A,I0)', 'normal axis ', c%get_normal_axis(), ', first defect open loop ', defect == CURVE_OPEN_LOOP, &
                          ' at loop ', l
do l = 1, c%get_loops_number()
   print '(A,I0,3A,I0,A,F8.4,A,L1)', 'loop ', l, ' (', trim(ROLE(l)), '): nesting depth ', c%nesting_depth(l), &
                                     ', signed area ', c%get_loop_area(l), ', closed ', c%is_loop_closed(l)
enddo

call fig_open(w, 'ref-2d-geometry-conventions', 'solid on the left of every segment: outer CCW, hole CW, island CCW; '// &
              'an open loop', -90._R8P, 90._R8P, width=760_I4P, height=420_I4P)
inside = 0
do j = 1, 35
   do i = 1, 70
      uv = [-0.3_R8P + H * (i - 0.5_R8P), -0.3_R8P + H * (j - 0.5_R8P)]
      if (.not. c%is_point_inside(uv, sign_algorithm=SIGN_FAST_WINDING_NUMBER)) cycle
      inside = inside + 1
      call fig_quad(w, 'field', vector_R8P(uv(1) - H / 2, uv(2) - H / 2, 0._R8P), &
                    vector_R8P(uv(1) + H / 2, uv(2) - H / 2, 0._R8P), vector_R8P(uv(1) + H / 2, uv(2) + H / 2, 0._R8P), &
                    vector_R8P(uv(1) - H / 2, uv(2) + H / 2, 0._R8P), -0.75_R8P)
   enddo
enddo
print '(A,I0)', 'cells inside (fast winding number) ', inside
n = c%get_vertices_number()
allocate(s1(n), s2(n), c1(2 * n), c2(2 * n))
n = 0 ; k = 0
do l = 1, c%get_loops_number()
   call c%get_loop_range(l, first, last)
   closed = c%is_loop_closed(l)
   do i = first, merge(last, last - 1, closed)
      n = n + 1 ; s1(n) = vtx(i) ; s2(n) = vtx(merge(first, i + 1, i == last))
      uv = c%get_vertex(i) ; e = c%get_vertex(merge(first, i + 1, i == last)) - uv
      dl = norm2(e) ; e = e / dl ; nn = [-e(2), e(1)] ; m = uv + 0.5_R8P * dl * e
      k = k + 1 ; c1(k) = pv(m + 0.5_R8P * S * e) ; c2(k) = pv(m - 0.5_R8P * S * e + 0.6_R8P * S * nn)
      k = k + 1 ; c1(k) = pv(m + 0.5_R8P * S * e) ; c2(k) = pv(m - 0.5_R8P * S * e - 0.6_R8P * S * nn)
   enddo
enddo
call fig_segments(w, 'accent', s1(1:n), s2(1:n))
call fig_segments(w, 'accent', c1(1:k), c2(1:k))
call fig_label(w, 'accent', pv([0.05_R8P, 3.15_R8P]), 'outer: CCW')
call fig_label(w, 'accent', pv([0.6_R8P, 2.65_R8P]), 'hole: CW')
call fig_label(w, 'accent', pv([1.6_R8P, 0.75_R8P]), 'island: CCW')
call fig_label(w, 'accent', pv([4.9_R8P, 2.75_R8P]), 'open loop')
call fig_segments(w, 'box', [pv([4.6_R8P, -0.1_R8P]), pv([4.6_R8P, -0.1_R8P])], &
                            [pv([5.2_R8P, -0.1_R8P]), pv([4.6_R8P, 0.5_R8P])])
call fig_label(w, 'box', pv([5.25_R8P, -0.15_R8P]), 'u = x')
call fig_label(w, 'box', pv([4.45_R8P, 0.6_R8P]), 'v = y')
call fig_close(w)
contains
   function vtx(i) result(r)
   !< Vertex i of the curve, lifted.
   integer(I4P), intent(in) :: i !< Vertex.
   type(vector_R8P)         :: r !< Point.

   r = pv(c%get_vertex(i))
   endfunction vtx

   function pv(p) result(r)
   !< A point of the plane, lifted over the cells.
   real(R8P), intent(in) :: p(2) !< (u, v).
   type(vector_R8P)      :: r    !< Point.

   r = vector_R8P(p(1), p(2), 0.002_R8P)
   endfunction pv
endprogram fig_ref_2d_geometry_conventions
