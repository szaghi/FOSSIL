program ck_probes_device
use fossil, only : curve_object
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter   :: h = 0.05_R8P
type(curve_object)     :: body
real(R8P)              :: u(68), v(68), nodes(2, 441), d(441)
real(R8P), allocatable :: ghosts(:,:), lengths(:,:), samples(:,:,:), samples_dev(:,:,:), images(:,:), images_dev(:,:)
logical,   allocatable :: visible(:,:), visible_dev(:,:), seen(:), seen_dev(:)
integer(I4P)           :: i, j, k, m

u(1:64) = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]   ! an ellipse, counter-clockwise
v(1:64) = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
u(65:68) = [0.05_R8P, 0.11_R8P, 0.11_R8P, 0.05_R8P] ; v(65:68) = [0.40_R8P, 0.40_R8P, 0.64_R8P, 0.64_R8P] ! an island
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 4_I4P])
nodes = reshape([((real(i, R8P) * h, real(j, R8P) * h, i = 0, 20), j = 0, 20)], [2, 441])   ! grid nodes, h = 0.05
call body%distance_many(points=nodes, distances=d)
m = count(d < 0._R8P .and. d > -h)
ghosts = reshape(pack(nodes, spread(d < 0._R8P .and. d > -h, 1, 2)), [2, m])     ! inside, within one spacing
lengths = spread([(k * sqrt(2._R8P) * h, k = 1, 3)], 2, m)                       ! an image line: 3 samples per ghost
allocate(samples(2, 3, m), samples_dev(2, 3, m), visible(3, m), visible_dev(3, m))
allocate(images(2, m), images_dev(2, m), seen(m), seen_dev(m))
call body%probe_many(points=ghosts, lengths=lengths, probe_points=samples, visible=visible)          ! CPU
call body%image_points(points=ghosts, images=images, visible=seen)
call body%enter_device                                                                               ! map the tree
call body%probe_many_device(points=ghosts, lengths=lengths, probe_points=samples_dev, visible=visible_dev)
call body%image_points_device(points=ghosts, images=images_dev, visible=seen_dev)
call body%exit_device
print '(A,I0,A,3(1X,I0),A,I0)', 'ghost points ', m, ', visible samples per position', &
      (count(visible(k, :)), k = 1, 3), ', visible images ', count(seen)
print '(A,I0)', 'device vs CPU mismatches: ', count(samples /= samples_dev) + count(visible .neqv. visible_dev) + &
                                               count(images /= images_dev) + count(seen .neqv. seen_dev)
endprogram ck_probes_device
