call body%distance_many_device(points=centre, distances=d, status=status)                     ! not resident
print '(A,L1)', 'not entered:      status invalid ', status /= STATUS_OK
call body%enter_device
call body%distance_many_device(points=centre, distances=d, is_signed=.true., sign_algorithm=SIGN_SOLID_ANGLE, &
                               status=status)                                                 ! no device port
print '(A,L1)', 'solid angle:      status invalid ', status /= STATUS_OK
call body%distance_many_device(points=centre(1:10), distances=d, status=status)                ! sizes differ
print '(A,L1)', 'size mismatch:    status invalid ', status /= STATUS_OK
call body%exit_device
call octree%load_from_file(file_name='cube-01.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call octree%enter_device                                                                      ! refused: BVH only
call octree%distance_many_device(points=centre(1:10), distances=d(1:10), status=status)
print '(A,L1)', 'octree:           status invalid ', status /= STATUS_OK
