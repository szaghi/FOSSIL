!run figref-usage-closest fig_ref_usage_closest
!fig ref-usage-closest
!< Figure of the usage reference, distance queries: compute_distance from probe points around the dragon, each with its
!< closest point and the pseudo-normal there (short segment from the closest point); the first probe is (0, 0, 25), the
!< point of ref_usage_distance.
program fig_ref_usage_closest
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
integer(I4P), parameter  :: n = 9
type(surface_stl_object) :: dragon
type(vector_R8P)         :: p(n), foot(n), tip(n), normal
real(R8P)                :: d, a
integer(I4P)             :: u, i

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
p(1) = vector_R8P(0._R8P, 0._R8P, 25._R8P)
do i = 2, n                                                      ! a ring about the dragon, at mid height
   a = 2._R8P * pi * (i - 2) / (n - 1)
   p(i) = vector_R8P(-0.6_R8P + 16._R8P * cos(a), 0.4_R8P + 11._R8P * sin(a), 12.5_R8P)
enddo
do i = 1, n
   call dragon%compute_distance(point=p(i), distance=d, is_square_root=.true., closest_point=foot(i), normal=normal)
   tip(i) = foot(i) + 2._R8P * normal
enddo
call fig_open(u, 'ref-usage-closest', 'compute_distance: probes, closest points, pseudo-normals there', &
              -60._R8P, 55._R8P)
call fig_surface(u, 'body', dragon)
call fig_segments(u, 'probe', p, foot)
call fig_segments(u, 'accent', foot, tip)
call fig_points(u, 'probe', p)
call fig_points(u, 'foot', foot)
call fig_close(u)
endprogram fig_ref_usage_closest
