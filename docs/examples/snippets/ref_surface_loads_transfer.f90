!< Reference page "Surface loads", moving fields between the wall and the grid: a marching-cubes sphere of radius 0.5 at
!< a pose with a wall temperature per vertex, T = 300 + 20 z in the body frame. `transfer_to_closure` reads it at the feet
!< and the ILW crossings of the closure table (exact for a linear field); the samples of order 2 carry a wall stress
!< back to the vertices with `samples_to_vertices` (a linear one exactly, a curved one to O(h^2)); then the table
!< resident on the device, compared by bits.
program ref_surface_loads_transfer
use fossil, only : surface_stl_object, facet_object, extract_isosurface, ib_closure_t, fd_closure_3d, ib_closure_reset, &
                   load_samples_t, load_samples_3d, transfer_to_closure, samples_to_vertices, rigid_pose_t, &
                   pose_from_axis_angle, to_body, to_world
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object), target  :: sphere
type(ib_closure_t)                :: cl, cr
type(load_samples_t)              :: ls
type(rigid_pose_t)                :: pose
type(vertex_pool_object), pointer :: pool
type(facet_object), allocatable   :: facets(:)
real(R8P),          allocatable   :: tw(:,:), tf(:,:), ti(:,:,:), td(:,:), tid(:,:,:), x(:), tau(:,:), tv(:,:)
type(vector_R8P)                  :: c
real(R8P)                         :: f(33, 33, 33), q(3), w(3), err, el, ec
integer(I4P)                      :: i, j, k, g, sl, ncross

do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
pose = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.5_R8P, [0.3_R8P, 0.2_R8P, 0.1_R8P])
pool => sphere%get_vertex_pool()
allocate(tw(1, pool%vertex_count()))
do i = 1, pool%vertex_count()
   c = pool%coord(i) ; tw(1, i) = 300._R8P + 20._R8P * c%z          ! the wall temperature per vertex (body frame)
enddo
x = [(-0.6_R8P + 1.2_R8P * (i - 0.5_R8P) / 24, i = 1, 24)]          ! 24 cell centres per axis around the body
call fd_closure_3d(sphere, x + 0.3_R8P, x + 0.2_R8P, x + 0.1_R8P, [2_I4P, 2_I4P, 2_I4P], cl, pose=pose, ilw=.true.)
allocate(tf(1, cl%n), ti(1, 6, cl%n))
call transfer_to_closure(sphere, cl, tf, vertex_field=tw, ilw_values=ti, pose=pose)   ! T at the feet, the crossings
err = 0._R8P ; ncross = 0
do g = 1, cl%n
   q = to_body(pose, cl%foot(:, g)) ; err = max(err, abs(tf(1, g) - (300._R8P + 20._R8P * q(3))))
   do sl = 1, 6
      if (cl%ilw_element(sl, g) < 1) cycle
      q = cl%points(:, g) ; q((sl + 1) / 2) = cl%ilw_t(sl, g) ; q = to_body(pose, q) ; ncross = ncross + 1
      err = max(err, abs(ti(1, sl, g) - (300._R8P + 20._R8P * q(3))))
   enddo
enddo
print '(A,I0,A,I0,A,L1)', 'ghosts ', cl%n, ', ILW crossings ', ncross, ';  T at the feet and crossings exact to 1e-12: ', &
      err < 1.e-12_R8P
call load_samples_3d(sphere, samples=ls, order=2_I4P, pose=pose)   ! a wall stress known at the samples (world frame)
allocate(tau(2, ls%n), tv(2, pool%vertex_count()))
tau(1, :) = 0.1_R8P * (1._R8P + ls%points(1, :) - ls%points(3, :))  ! linear
tau(2, :) = 0.1_R8P * cos(4 * ls%points(3, :))                     ! curved
call samples_to_vertices(sphere, ls, tau, tv)                       ! per vertex, area-weighted
el = 0._R8P ; ec = 0._R8P
do i = 1, pool%vertex_count()
   c = pool%coord(i) ; w = to_world(pose, [c%x, c%y, c%z])
   el = max(el, abs(tv(1, i) - 0.1_R8P * (1._R8P + w(1) - w(3)))) ; ec = max(ec, abs(tv(2, i) - 0.1_R8P * cos(4 * w(3))))
enddo
print '(A,I0,A,L1,A,F8.5)', 'vertices ', pool%vertex_count(), ';  linear stress exact to 1e-14: ', el < 1.e-14_R8P, &
      ';  curved stress, largest error ', ec
call sphere%enter_device
call fd_closure_3d(sphere, x + 0.3_R8P, x + 0.2_R8P, x + 0.1_R8P, [2_I4P, 2_I4P, 2_I4P], cr, pose=pose, ilw=.true., &
                   resident=.true.)
allocate(td(1, cl%n), tid(1, 6, cl%n))
call transfer_to_closure(sphere, cr, td, vertex_field=tw, ilw_values=tid, pose=pose, device=.true.)  ! read on the device
print '(A,L1)', 'device = host, bit for bit: ', all(transfer(td, 0_I8P, size(td)) == transfer(tf, 0_I8P, size(tf))) .and. &
      all(transfer(tid, 0_I8P, size(tid)) == transfer(ti, 0_I8P, size(ti)))
call ib_closure_reset(cr)
call sphere%exit_device
endprogram ref_surface_loads_transfer
