program ck_field_transfer
use fossil, only : surface_stl_object, curve_object, ib_closure_t, fd_closure_3d, load_samples_t, load_samples_3d, &
                   transfer_field_3d, transfer_field_2d, transfer_to_closure, transfer_to_samples, samples_to_vertices
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object), target  :: cube
type(curve_object)                :: square
type(ib_closure_t)                :: cl
type(load_samples_t)              :: ls
type(vertex_pool_object), pointer :: pool
real(R8P), allocatable            :: tw(:,:), tf(:,:), tg(:,:), ts(:,:), tv(:,:), tau(:,:), x(:)
type(vector_R8P)                  :: c
real(R8P)                         :: err, t2(1, 1)
integer(I4P)                      :: i

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube [0, 1]^3
pool => cube%get_vertex_pool()
allocate(tw(1, pool%vertex_count()))
do i = 1, pool%vertex_count()
   c = pool%coord(i) ; tw(1, i) = 300._R8P + 10._R8P * c%z              ! a wall temperature per vertex
enddo
x = [(-0.15_R8P + 0.1_R8P * i, i = 0, 13)]                             ! the solver's points, spacing 0.1
call fd_closure_3d(cube, x, x, x, [2_I4P, 2_I4P, 2_I4P], cl)           ! ghost points and their feet
allocate(tf(1, cl%n), tg(1, cl%n))
call transfer_to_closure(cube, cl, tf, vertex_field=tw)                ! T at every foot
call transfer_field_3d(cube, cl%element, cl%foot, tg, vertex_field=tw) ! the same: any points, each with its facet
print '(A,I0,A,F7.2,A,F7.2,A,L1)', 'ghosts ', cl%n, ';  T at the feet from ', minval(tf), ' to ', maxval(tf), &
      ';  the same by transfer_field_3d: ', all(tf == tg)
call load_samples_3d(cube, samples=ls, order=2_I4P)                    ! the way back
allocate(ts(1, ls%n), tv(1, pool%vertex_count()), tau(1, ls%n))
call transfer_to_samples(cube, ls, ts, vertex_field=tw)                ! T at the load samples
call samples_to_vertices(cube, ls, ts, tv)                             ! and back: unchanged
print '(A,L1)', 'T, vertices -> samples -> vertices, unchanged to 1e-12: ', maxval(abs(tv - tw)) < 1.e-12_R8P
tau(1, :) = 0.5_R8P * ls%points(1, :)                                  ! a wall stress computed at the samples
call samples_to_vertices(cube, ls, tau, tv)                            ! per vertex
err = 0._R8P
do i = 1, pool%vertex_count()
   c = pool%coord(i) ; err = max(err, abs(tv(1, i) - 0.5_R8P * c%x))
enddo
print '(A,I0,A,L1)', 'vertices ', pool%vertex_count(), ';  the stress 0.5 x at every vertex to 1e-14: ', err < 1.e-14_R8P
call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call transfer_field_2d(square, [1_I4P], reshape([0.25_R8P, 0._R8P], [2, 1]), t2, &
                       vertex_field=reshape([10._R8P, 20._R8P, 30._R8P, 40._R8P], [1, 4]))   ! 2-D: on segment 1 -> 2
print '(A,F6.2)', '2-D: segment 1 at u = 0.25: ', t2
endprogram ck_field_transfer
