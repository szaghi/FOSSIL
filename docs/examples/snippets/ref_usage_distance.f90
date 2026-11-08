!< Usage reference, distance queries: unsigned (squared, then Euclidean), signed with two sign algorithms, the
!< point-in-polyhedron test on the unit cube; compute_distance, distance_many and the two tree kinds on the dragon.
program ref_usage_distance
use fossil, only : surface_stl_object, SIGN_RAY_INTERSECTIONS, AABB_TREE_OCTREE, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface, dragon, dragon_octree
type(vector_R8P)         :: p, closest, normal
type(vector_R8P)         :: points(21)
real(R8P)                :: d2, d, sd, sd_ray, dist, d_bvh(21), d_oct(21)
logical                  :: inside, same
integer(I4P)             :: facet, edge, vertex, i, status

call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! the unit cube [0, 1]^3
d2 = surface%distance(point=2._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P)  ! squared (the default)
print '(A,F10.6)', 'squared distance of (2, 0, 0)         ', d2
p = 2._R8P * ex_R8P + 2._R8P * ey_R8P + 0.5_R8P * ez_R8P
d = surface%distance(point=p, is_square_root=.true.)
print '(A,F10.6)', 'distance of (2, 2, 0.5)               ', d
p = 0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.25_R8P * ez_R8P
sd = surface%distance(point=p, is_signed=.true., is_square_root=.true.)
! Override the sign algorithm for benchmarking or for meshes whose orientation
! cannot be sanitised (e.g. open shells where outward is ill-defined).
sd_ray = surface%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_RAY_INTERSECTIONS)
print '(A,F10.6,A,F10.6)', 'signed distance of (0.5, 0.5, 0.25)   ', sd, ', by ray intersections ', sd_ray
inside = surface%is_point_inside(point=0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.5_R8P * ez_R8P)
if (inside) print '(A)', 'point (0.5, 0.5, 0.5) is inside the surface'
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
p = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 25._R8P * ez_R8P
call dragon%compute_distance(point=p, distance=dist, is_square_root=.true., facet_index=facet, edge_index=edge, &
                             vertex_index=vertex, closest_point=closest, normal=normal)
print '(A,F9.4)',   'dragon, (0, 0, 25): distance ', dist
print '(A,3F9.4)',  '  closest point ', closest%x, closest%y, closest%z
print '(A,3F9.4)',  '  normal        ', normal%x, normal%y, normal%z
print '(A,I0,A,I0,A,I0)', '  facet ', facet, ', local edge ', edge, ', local vertex ', vertex
do i = 1, size(points)                                     ! 21 points on a line through the dragon
   points(i) = (-12._R8P + 1.1_R8P * (i - 1)) * ex_R8P + 0.5_R8P * ey_R8P + 12._R8P * ez_R8P
enddo
call dragon%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true., status=status)
if (status /= STATUS_OK) error stop 'distance_many: points and distances differ in size'
same = .true.
do i = 1, size(points)
   same = same .and. (d_bvh(i) == dragon%distance(point=points(i), is_signed=.true., is_square_root=.true.))
enddo
print '(A,I0,A,L1)', 'distance_many: ', size(points), ' points, equal to distance bit for bit ', same
print '(A,I0)', '  inside ', count(d_bvh < 0._R8P)
call dragon_octree%load_from_file(file_name='dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call dragon_octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,L1)', 'octree and SAH BVH: equal bit for bit ', all(d_bvh == d_oct)
endprogram ref_usage_distance
