!< Reference page "Immersed-boundary queries": two unit cubes 0.2 apart in one surface. Bodies, probes and images with
!< their visibility across the gap, an ambiguous foot point at mid-gap, wall distance with a cutoff, grid-adapted
!< markers.
program ref_ib_queries
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object)   :: a, b
type(vector_R8P)           :: ghost(3), gap(2), wall(3)
real(R8P)                  :: d(3), probe(3, 3), img(3, 3), second(3, 2), dg(2), dw(3)
real(R8P), allocatable     :: mc(:,:), ma(:), mn(:,:)
real(R8P)                  :: fraction
integer(I4P)               :: bod(3), bw(3), i
logical                    :: vis(3), amb(2)

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! [0, 1]^3
b = a
call b%translate(x=1.2_R8P)                                           ! [1.2, 2.2] x [0, 1]^2
call a%merge_solids(b)
call a%analyze
print '(A,I0,A,I0,A,I0)', 'facets ', a%get_facets_number(), ', shells ', a%get_shells_number(), &
                          ', bodies ', a%get_bodies_number()
! three ghost points inside the first cube: near the face x = 1, deeper, near the top face
ghost = [vector_R8P(0.95_R8P, 0.3_R8P, 0.4_R8P), vector_R8P(0.75_R8P, 0.3_R8P, 0.4_R8P), &
         vector_R8P(0.9_R8P, 0.3_R8P, 0.97_R8P)]
call a%query_many(points=ghost, distances=d, bodies=bod)
call a%probe_many(points=ghost, lengths=[0.3_R8P, 0.3_R8P, 0.3_R8P], probe_points=probe, visible=vis)
do i = 1, 3
   print '(A,3F5.2,A,F6.3,A,I0,A,3F6.3,A,L1)', 'ghost', ghost(i)%x, ghost(i)%y, ghost(i)%z, '  d ', d(i), &
         '  body ', bod(i), '  probe', probe(:, i), '  visible ', vis(i)
enddo
call a%image_points(points=ghost, images=img, visible=vis)
do i = 1, 3
   print '(A,I0,A,3F6.3,A,L1)', 'image ', i, ':', img(:, i), '  visible ', vis(i)
enddo
! a point at mid-gap has two foot points, one on each cube; a point nearer the first cube has one
gap = [vector_R8P(1.1_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.05_R8P, 0.5_R8P, 0.5_R8P)]
call a%query_many(points=gap, distances=dg, ambiguous=amb, second_foot_points=second)
do i = 1, 2
   print '(A,F5.2,A,F6.3,A,L1,A,F5.2)', 'gap point x =', gap(i)%x, '  d ', dg(i), '  ambiguous ', amb(i), &
         '  second foot x =', second(1, i)
enddo
! unsigned wall distance and body, capped at 0.3
wall = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.3_R8P, 0.5_R8P, 0.5_R8P), &
        vector_R8P(1.7_R8P, 0.5_R8P, 1.6_R8P)]
call a%wall_distance_many(points=wall, distances=dw, bodies=bw, cutoff=0.3_R8P)
do i = 1, 3
   print '(A,3F5.2,A,F6.3,A,I0)', 'wall distance at', wall(i)%x, wall(i)%y, wall(i)%z, '  ', dw(i), '  body ', bw(i)
enddo
call b%surface_markers(dx=0.1_R8P, centroids=mc, areas=ma, normals=mn, spacing_fraction=fraction, &
                       preserve_features=.true.)
print '(A,I0,A,F9.6,A,F5.3)', 'markers on one cube, dx 0.1: ', size(ma), ', area ', sum(ma), &
                              ', spacings in [0.8, 1.2] dx: ', fraction
endprogram ref_ib_queries
