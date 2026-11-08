!run ck_tree_kind ck_tree_kind
program ck_tree_kind
use fossil, only : surface_stl_object, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: bvh, octree
type(vector_R8P)         :: p(1000)
real(R8P)                :: d_bvh(1000), d_octree(1000)
integer(I4P)             :: i, j, k, n

call bvh%load_from_file(file_name='dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_SAH_BVH)   ! default
call octree%load_from_file(file_name='dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
n = 0
do k = 1, 10 ; do j = 1, 10 ; do i = 1, 10                      ! 1000 points over the dragon's box
   n = n + 1 ; p(n) = vector_R8P(-11._R8P + 2.2_R8P * i, -4.5_R8P + 1._R8P * j, 4.5_R8P + 1.6_R8P * k)
enddo ; enddo ; enddo
call bvh%distance_many(points=p, distances=d_bvh, is_signed=.true., is_square_root=.true.)
call octree%distance_many(points=p, distances=d_octree, is_signed=.true., is_square_root=.true.)
print '(A,I0,A,ES10.3)', 'points ', n, ', largest |d_bvh - d_octree| ', maxval(abs(d_bvh - d_octree))
print '(A,2F9.4)', 'd range ', minval(d_bvh), maxval(d_bvh)
endprogram ck_tree_kind
