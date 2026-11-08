!run ck_image_points ck_image_points
program ck_image_points
use fossil, only : curve_object
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter    :: h = 0.05_R8P
type(curve_object)      :: body
real(R8P)               :: u(64), v(64), nodes(2, 441), d(441)
real(R8P), allocatable  :: ghosts(:,:), foot(:,:), images(:,:), probes(:,:)
logical,   allocatable  :: visible(:), clear(:)
integer(I4P)            :: i, j, m

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]   ! an ellipse, counter-clockwise
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
nodes = reshape([((real(i, R8P) * h, real(j, R8P) * h, i = 0, 20), j = 0, 20)], [2, 441])   ! grid nodes, h = 0.05
call body%distance_many(points=nodes, distances=d)                               ! negative inside
m = count(d < 0._R8P .and. d > -h)
ghosts = reshape(pack(nodes, spread(d < 0._R8P .and. d > -h, 1, 2)), [2, m])     ! inside, within one spacing
allocate(foot(2, m), images(2, m), probes(2, m), visible(m), clear(m))
call body%image_points(points=ghosts, images=images, visible=visible, foot_points=foot)   ! mirror through the wall
call body%probe_many(points=ghosts, lengths=[(2 * h, i = 1, m)], probe_points=probes, visible=clear)  ! 2h out
print '(A,I0,A,I0,A,I0)', 'ghost points ', m, ', images visible ', count(visible), ', probes visible ', count(clear)
print '(A,2F9.5)', 'ghost 1  ', ghosts(:, 1)
print '(A,2F9.5)', '  foot   ', foot(:, 1)
print '(A,2F9.5)', '  image  ', images(:, 1)
print '(A,2F9.5)', '  probe  ', probes(:, 1)
endprogram ck_image_points
