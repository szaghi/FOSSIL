call dragon_octree%load_from_file(file_name='dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call dragon_octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,L1)', 'octree and SAH BVH: equal bit for bit ', all(d_bvh == d_oct)
