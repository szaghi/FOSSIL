!run figref-2d-geometry-queries fig_ref_2d_geometry_queries
!fig ref-2d-geometry-queries
!< Figure of the reference page "2-D geometry (curves)", section What it does: a curve of two bodies, an ellipse with a
!< square hole and a triangle, queried natively in 2-D. The background cells are coloured by the signed distance of
!< their centres from `distance_many` (teal inside, negative; ember outside), the curve is drawn in gold, and six
!< query points (white) are joined to their closest points from `compute_distance`.
program fig_ref_2d_geometry_queries
use fossil,     only : curve_object, CURVE_VALID
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter   :: PI = acos(-1._R8P), H = 0.05_R8P
integer(I4P), parameter   :: NE = 48, MX = 100, MY = 50
type(curve_object)        :: c
real(R8P)                 :: u(NE + 7), v(NE + 7), pts(2, MX * MY), d(MX * MY), q(2, 6), foot(2), dq
type(vector_R8P)          :: a(6), b(6)
integer(I4P)              :: i, j, k, l, n, w, first, last, defect
type(vector_R8P), allocatable :: s1(:), s2(:)

do i = 1, NE                                                   ! ellipse 1.6 x 0.9, counter-clockwise
   u(i) = 1.6_R8P * cos(2 * PI * (i - 1) / NE) ; v(i) = 0.9_R8P * sin(2 * PI * (i - 1) / NE)
enddo
u(NE + 1:NE + 4) = [-0.4_R8P, -0.4_R8P, 0.4_R8P, 0.4_R8P]       ! square hole, clockwise
v(NE + 1:NE + 4) = [-0.3_R8P, 0.3_R8P, 0.3_R8P, -0.3_R8P]
u(NE + 5:NE + 7) = [2.0_R8P, 2.4_R8P, 2.0_R8P]                  ! triangle, counter-clockwise
v(NE + 5:NE + 7) = [-0.6_R8P, 0._R8P, 0.6_R8P]
call c%adopt_loops(u=u, v=v, loop_sizes=[NE, 4_I4P, 3_I4P])
defect = c%check_validity(loop=l, check_self_intersections=.true.)
print '(A,L1,A,I0,A,I0,A,F8.5)', 'valid ', defect == CURVE_VALID, ', loops ', c%get_loops_number(), ', bodies ', &
                                 c%get_bodies_number(), ', area ', c%get_area()
n = 0
do j = 1, MY
   do i = 1, MX
      n = n + 1 ; pts(:, n) = [-2.5_R8P + H * (i - 0.5_R8P), -1.25_R8P + H * (j - 0.5_R8P)]
   enddo
enddo
call c%distance_many(pts, d)
print '(A,I0,A,I0,A,F8.5,A,F8.5)', 'cells ', n, ', inside ', count(d < 0._R8P), ', distance from ', minval(d), ' to ', &
                                   maxval(d)
q(:, 1) = [-2.2_R8P, 0.9_R8P] ; q(:, 2) = [0.1_R8P, 0.05_R8P] ; q(:, 3) = [1.0_R8P, 0.2_R8P]
q(:, 4) = [2.6_R8P, -0.9_R8P] ; q(:, 5) = [1.85_R8P, 0.1_R8P] ; q(:, 6) = [-0.9_R8P, -0.6_R8P]
do k = 1, 6
   call c%compute_distance(q(:, k), distance=dq, closest_point=foot)
   print '(A,2F6.2,A,F8.5,A,2F8.4)', 'p ', q(:, k), ': d ', dq, ', foot ', foot
   a(k) = vector_R8P(q(1, k), q(2, k), 0.003_R8P) ; b(k) = vector_R8P(foot(1), foot(2), 0.003_R8P)
enddo

call fig_open(w, 'ref-2d-geometry-queries', 'signed distance of the cells (teal inside, ember outside), '// &
              'six points and their closest points', -90._R8P, 90._R8P, width=760_I4P, height=440_I4P)
do n = 1, MX * MY
   call fig_quad(w, 'field', vector_R8P(pts(1, n) - H / 2, pts(2, n) - H / 2, 0._R8P), &
                 vector_R8P(pts(1, n) + H / 2, pts(2, n) - H / 2, 0._R8P), &
                 vector_R8P(pts(1, n) + H / 2, pts(2, n) + H / 2, 0._R8P), &
                 vector_R8P(pts(1, n) - H / 2, pts(2, n) + H / 2, 0._R8P), max(-1._R8P, min(1._R8P, d(n) / 0.6_R8P)))
enddo
allocate(s1(c%get_vertices_number()), s2(c%get_vertices_number()))
n = 0
do l = 1, c%get_loops_number()
   call c%get_loop_range(l, first, last)
   do i = first, last
      n = n + 1 ; s1(n) = vtx(i) ; s2(n) = vtx(merge(first, i + 1, i == last))
   enddo
enddo
call fig_segments(w, 'accent', s1(1:n), s2(1:n))
call fig_segments(w, 'probe', a, b)
call fig_points(w, 'probe', a)
call fig_points(w, 'foot', b)
call fig_close(w)
contains
   function vtx(i) result(r)
   !< Vertex i of the curve, lifted.
   integer(I4P), intent(in) :: i !< Vertex.
   type(vector_R8P)         :: r !< Point.
   real(R8P)                :: uv(2)

   uv = c%get_vertex(i) ; r = vector_R8P(uv(1), uv(2), 0.002_R8P)
   endfunction vtx
endprogram fig_ref_2d_geometry_queries
