call surface%distance_many_device(points=points, distances=d_host, is_signed=.true., &
                                  sign_algorithm=SIGN_RAY_INTERSECTIONS, status=status)  ! not ported: refused
print '(A,L1)', 'ray-intersection sign on the device refused: ', status /= STATUS_OK
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call octree%enter_device                                    ! the device traverses the SAH BVH only: refused
print '(A,L1)', 'octree on the device: ', octree%is_on_device
call surface%exit_device                                    ! release the device copies, before destroy or analyze
print '(A,L1)', 'after exit_device: ', surface%is_on_device
