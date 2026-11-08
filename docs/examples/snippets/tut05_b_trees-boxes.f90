call draw(bvh, 'tut05-bvh', 'SAH BVH', 4_I4P)
call draw(octree, 'tut05-octree', 'octree', 2_I4P)
contains
   subroutine draw(s, name, title, depth)
   !< The boxes of the tree of s `depth` levels below the root, over the dragon.
   type(surface_stl_object), intent(in) :: s           !< Surface.
   character(*),             intent(in) :: name, title !< Figure name, tree name.
   integer(I4P),             intent(in) :: depth       !< Levels below the root.
   type(vector_R8P), allocatable        :: bmin(:), bmax(:)
   integer(I4P),     allocatable        :: level(:)
   integer(I4P)                         :: c, u, all
   character(96)                        :: caption

   call s%tree_boxes(bmin=bmin, bmax=bmax)                    ! the whole tree
   all = size(bmin)
   call s%tree_boxes(bmin=bmin, bmax=bmax, depth=level, max_depth=depth)
   write(caption, '(A,A,I0,A,I0,A,I0,A)') title, ': ', count(level == depth), ' boxes ', depth, &
                                          ' levels below the root, ', all, ' nodes in all'
   print '(A)', trim(caption)
