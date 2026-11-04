call bvh%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)       ! default AABB_TREE_SAH_BVH
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call bvh%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true., status=status)
if (status /= STATUS_OK) error stop 'distance_many failed'
call octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,I0)',     'points                ', size(points)
print '(A,ES10.3)', 'max ||bvh| - |oct||  ', maxval(abs(abs(d_bvh) - abs(d_oct)))
print '(A,I0)',     'bvh, oct signs differ ', count(d_bvh * d_oct < 0._R8P)
