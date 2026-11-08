!< Tutorial 5, second run: the same batch on the two kinds of AABB tree. The figures draw the boxes of the BVH four
!< levels below the root and of the octree two levels below (its leaves at the default depth), read from each
!< surface's own tree by tree_boxes (issue #89): the BVH splits the facets, the octree splits the space.
program tut05_b_trees
use fossil_fig,              only : fig_open, fig_surface, fig_box, fig_close
use fossil,                  only : surface_stl_object, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16                      ! points per direction
type(surface_stl_object)      :: bvh, octree
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_bvh(:), d_oct(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p

allocate(points(n**3), d_bvh(n**3), d_oct(n**3))
h = 1.2_R8P / (n - 1)                                        ! the cube [-0.6, 0.6]^3
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call bvh%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)       ! default AABB_TREE_SAH_BVH
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call bvh%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true.)
call octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,I0)',     'points                ', size(points)
print '(A,ES10.3)', 'max ||bvh| - |oct||  ', maxval(abs(abs(d_bvh) - abs(d_oct)))
print '(A,I0)',     'bvh, oct signs differ ', count(d_bvh * d_oct < 0._R8P)
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
   call fig_open(u, name, trim(caption), -60._R8P, 22._R8P)
   call fig_surface(u, 'body', s)
   do c = 1, size(level)
      if (level(c) == depth) call fig_box(u, 'accent', bmin(c), bmax(c))
   enddo
   call fig_close(u)
   endsubroutine draw
endprogram tut05_b_trees
