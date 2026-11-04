!< Figure of tutorial 4: twelve probes around the dragon, their closest points on it and the distance segments.
program fig_tut04
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax, probe(12), foot(12)
real(R8P)                :: d, t
integer(I4P)             :: u, k

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
do k = 1, 12                                                                   ! a wavy ring around the body
   t = 6.283185307179586_R8P * (k - 0.5_R8P) / 12._R8P
   probe(k)%x = 0.5_R8P * (bmin%x + bmax%x) + 0.75_R8P * (bmax%x - bmin%x) * cos(t)
   probe(k)%y = 0.5_R8P * (bmin%y + bmax%y) + 0.75_R8P * (bmax%y - bmin%y) * sin(t)
   probe(k)%z = 0.5_R8P * (bmin%z + bmax%z) + 0.25_R8P * (bmax%z - bmin%z) * sin(3 * t)
   call surface%compute_distance(point=probe(k), distance=d, is_signed=.true., is_square_root=.true., &
                                 closest_point=foot(k))
enddo
call fig_open(u, 'tut04-probes', 'compute_distance: probes, closest points and the distance', -60._R8P, 22._R8P)
call fig_surface(u, 'body', surface)
call fig_segments(u, 'foot', probe, foot)
call fig_points(u, 'probe', probe)
call fig_points(u, 'foot', foot)
call fig_close(u)
endprogram fig_tut04
