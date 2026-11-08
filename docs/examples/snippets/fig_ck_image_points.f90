!< Figure of the cookbook recipe "Image points and probes for a ghost-cell method": the ellipse of ck_image_points over
!< its grid (h = 0.05), seen from above, zoomed on its left end: the ghost nodes (inside, within h of the wall), each
!< joined through its foot point to its image point, and the probe 2h out along the normal from each foot point.
program fig_ck_image_points
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: h = 0.05_R8P
type(curve_object)            :: body
real(R8P)                     :: u(64), v(64), nodes(2, 441), d(441)
real(R8P), allocatable        :: ghosts(:,:), foot(:,:), images(:,:), probes(:,:)
logical,   allocatable        :: visible(:), clear(:)
type(vector_R8P), allocatable :: a(:), b(:)
integer(I4P)                  :: i, j, m, un

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
nodes = reshape([((real(i, R8P) * h, real(j, R8P) * h, i = 0, 20), j = 0, 20)], [2, 441])
call body%distance_many(points=nodes, distances=d)
m = count(d < 0._R8P .and. d > -h)
ghosts = reshape(pack(nodes, spread(d < 0._R8P .and. d > -h, 1, 2)), [2, m])
allocate(foot(2, m), images(2, m), probes(2, m), visible(m), clear(m))
call body%image_points(points=ghosts, images=images, visible=visible, foot_points=foot)
call body%probe_many(points=ghosts, lengths=[(2 * h, i = 1, m)], probe_points=probes, visible=clear)
call fig_open(un, 'ck-image-points', 'ghost nodes (teal), images (ember), probes 2h out (stone)', -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.36_R8P, 0.58_R8P, 0._R8P), 0.22_R8P)
allocate(a(441))
do i = 1, 441
   a(i) = vector_R8P(nodes(1, i), nodes(2, i), 0._R8P)
enddo
call fig_points(un, 'box', a)
a = [(vector_R8P(u(i), v(i), 0._R8P), i = 1, 64)] ; b = [(vector_R8P(u(mod(i, 64) + 1), v(mod(i, 64) + 1), 0._R8P), i = 1, 64)]
call fig_segments(un, 'accent', a, b)
a = [(p(ghosts(:, i)), i = 1, m)] ; b = [(p(images(:, i)), i = 1, m)]
call fig_segments(un, 'ghost', a, b)
call fig_points(un, 'inside', a)
call fig_points(un, 'outside', b)
a = [(p(foot(:, i)), i = 1, m)] ; b = [(p(probes(:, i)), i = 1, m)]
call fig_segments(un, 'ghost', a, b)
call fig_points(un, 'foot', a)
call fig_points(un, 'body', b)
call fig_close(un)
contains
   function p(w) result(r)
   !< A plane point (u, v) as a vector.
   real(R8P), intent(in) :: w(2) !< Point.
   type(vector_R8P)      :: r    !< Vector.

   r = vector_R8P(w(1), w(2), 0._R8P)
   endfunction p
endprogram fig_ck_image_points
