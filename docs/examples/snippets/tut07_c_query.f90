!< Tutorial 7, third run: signed distance to the section read back from dragon-section.dat: one point in detail, then
!< twelve points on an ellipse by distance_many. The figure, from above, joins each point to its closest point.
program tut07_c_query
use fossil_fig, only : fig_open, fig_segment, fig_points, fig_close
use fossil,     only : curve_object
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P

implicit none
real(R8P), parameter :: PI = 4._R8P * atan(1._R8P)
type(curve_object)   :: section
real(R8P)            :: cp(2), d, ring(2, 12), dr(12), uv(2), cpk(2)
integer(I4P)         :: segment, vertex, k, s, u
type(vector_R8P)     :: a, b

call section%load_from_file(file_name='dragon-section.dat')
call section%compute_distance(p=[0._R8P, 0._R8P], distance=d, closest_point=cp, closest_segment=segment, &
                              closest_vertex=vertex)
print '(A,F10.6,A,2F10.6,A,I0,A,I0)', 'd(0,0) ', d, '  closest ', cp, '  segment ', segment, '  vertex ', vertex
do k = 1, 12                                                    ! an ellipse about the origin
   ring(:, k) = [0.12_R8P * cos(2._R8P * PI * k / 12._R8P), 0.4_R8P * sin(2._R8P * PI * k / 12._R8P)]
enddo
call section%distance_many(points=ring, distances=dr)
print '(A,I0,A,F9.5,A,F9.5)', 'ellipse: inside ', count(dr < 0._R8P), ' of 12, d from ', minval(dr), ' to ', maxval(dr)
call fig_open(u, 'tut07-query', 'section z = 0: points inside (teal), outside (ember), closest points (gold)', &
              0._R8P, 90._R8P)
do s = 1, section%get_vertices_number()
   if (section%get_next(s) == 0_I4P) cycle
   uv = section%get_vertex(s) ; a = uv(1) * ex_R8P + uv(2) * ey_R8P
   uv = section%get_vertex(section%get_next(s)) ; b = uv(1) * ex_R8P + uv(2) * ey_R8P
   call fig_segment(u, 'body', a, b)
enddo
do k = 0, 12
   if (k == 0) then
      uv = 0._R8P ; d = section%distance([0._R8P, 0._R8P])
   else
      uv = ring(:, k) ; d = dr(k)
   endif
   call section%compute_distance(p=uv, distance=d, closest_point=cpk)
   a = uv(1) * ex_R8P + uv(2) * ey_R8P ; b = cpk(1) * ex_R8P + cpk(2) * ey_R8P
   call fig_segment(u, 'accent', a, b)
   call fig_points(u, merge('inside ', 'outside', d < 0._R8P), [a])
   call fig_points(u, 'foot', [b])
enddo
call fig_close(u)
endprogram tut07_c_query
