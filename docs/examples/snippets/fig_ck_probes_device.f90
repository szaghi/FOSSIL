!< Figure of the cookbook recipe "Image lines and probes on the GPU" and of the IB-queries section on image lines: the
!< ellipse and the island of ck_probes_device over their grid (h = 0.05), seen from above, zoomed on the gap between
!< them: the ghost nodes, each joined to its foot point and on along the normal through the three samples of its image
!< line; visible samples teal, samples behind the island (the segment from the foot crosses it) ember.
program fig_ck_probes_device
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: h = 0.05_R8P
type(curve_object)            :: body
real(R8P)                     :: u(68), v(68), nodes(2, 441), d(441)
real(R8P), allocatable        :: ghosts(:,:), lengths(:,:), samples(:,:,:), foot(:,:)
logical,   allocatable        :: visible(:,:)
type(vector_R8P), allocatable :: a(:), b(:)
integer(I4P)                  :: i, j, k, m, un, nv, nb

u(1:64) = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
v(1:64) = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
u(65:68) = [0.05_R8P, 0.11_R8P, 0.11_R8P, 0.05_R8P] ; v(65:68) = [0.40_R8P, 0.40_R8P, 0.64_R8P, 0.64_R8P]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 4_I4P])
nodes = reshape([((real(i, R8P) * h, real(j, R8P) * h, i = 0, 20), j = 0, 20)], [2, 441])
call body%distance_many(points=nodes, distances=d)
m = count(d < 0._R8P .and. d > -h)
ghosts = reshape(pack(nodes, spread(d < 0._R8P .and. d > -h, 1, 2)), [2, m])
lengths = spread([(k * sqrt(2._R8P) * h, k = 1, 3)], 2, m)
allocate(samples(2, 3, m), visible(3, m), foot(2, m))
call body%enter_device
call body%probe_many_device(points=ghosts, lengths=lengths, probe_points=samples, visible=visible, foot_points=foot)
call body%exit_device
nv = count(visible) ; nb = count(.not. visible)
call fig_open(un, 'ck-probes-device', 'ghost nodes, image lines of 3 samples: visible (teal), blocked (ember)', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.26_R8P, 0.52_R8P, 0._R8P), 0.26_R8P)
allocate(a(441))
do i = 1, 441
   a(i) = vector_R8P(nodes(1, i), nodes(2, i), 0._R8P)
enddo
call fig_points(un, 'box', a)
a = [(vector_R8P(u(i), v(i), 0._R8P), i = 1, 64), (vector_R8P(u(i), v(i), 0._R8P), i = 65, 68)]
b = [(vector_R8P(u(mod(i, 64) + 1), v(mod(i, 64) + 1), 0._R8P), i = 1, 64), &
     (vector_R8P(u(65 + mod(i - 64, 4)), v(65 + mod(i - 64, 4)), 0._R8P), i = 65, 68)]
call fig_segments(un, 'accent', a, b)
a = [(p(ghosts(:, i)), i = 1, m)] ; b = [(p(samples(:, 3, i)), i = 1, m)]
call fig_segments(un, 'ghost', a, b)
call fig_points(un, 'body', a)
a = [(p(foot(:, i)), i = 1, m)]
call fig_points(un, 'foot', a)
a = [((p(samples(:, k, i)), k = 1, 3), i = 1, m)]
call fig_points(un, 'inside', pack(a, reshape(visible, [3 * m])))
call fig_points(un, 'outside', pack(a, .not. reshape(visible, [3 * m])))
call fig_close(un)
print '(A,I0,A,I0,A,I0,A)', 'ghost nodes ', m, ': ', nv, ' samples visible, ', nb, ' blocked'
contains
   function p(w) result(r)
   !< A plane point (u, v) as a vector.
   real(R8P), intent(in) :: w(2) !< Point.
   type(vector_R8P)      :: r    !< Vector.

   r = vector_R8P(w(1), w(2), 0._R8P)
   endfunction p
endprogram fig_ck_probes_device
