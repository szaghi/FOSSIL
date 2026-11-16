call body%probe_many(points=ghosts, lengths=lengths, probe_points=samples, visible=visible)          ! CPU
call body%image_points(points=ghosts, images=images, visible=seen)
call body%enter_device                                                                               ! map the tree
call body%probe_many_device(points=ghosts, lengths=lengths, probe_points=samples_dev, visible=visible_dev)
call body%image_points_device(points=ghosts, images=images_dev, visible=seen_dev)
call body%exit_device
