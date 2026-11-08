call slice_surface(surface=a, normal_axis=3_I4P, level=0.75_R8P, curve=section)   ! the section z = 0.75
call lengths(section)                                                         ! each segment has its facet's patch
