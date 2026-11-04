!< Tutorial 5: distances of a whole grid of points in one call, SAH BVH and octree.
program tut05_batched
use fossil, only : surface_stl_object, AABB_TREE_OCTREE, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16                      ! points per direction
type(surface_stl_object)      :: bvh, octree
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_bvh(:), d_oct(:), d_wn(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p, status

allocate(points(n**3), d_bvh(n**3), d_oct(n**3), d_wn(n**3))
h = 1.2_R8P / (n - 1)                                        ! the cube [-0.6, 0.6]^3
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call bvh%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)       ! default AABB_TREE_SAH_BVH
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call bvh%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true., status=status)
if (status /= STATUS_OK) error stop 'distance_many failed'
call octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,I0)',     'points                ', size(points)
print '(A,ES10.3)', 'max ||bvh| - |oct||  ', maxval(abs(abs(d_bvh) - abs(d_oct)))
print '(A,I0)',     'bvh, oct signs differ ', count(d_bvh * d_oct < 0._R8P)
call bvh%distance_many(points=points, distances=d_wn, is_signed=.true., is_square_root=.true., &
                       sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I5,A,F10.6)', 'pseudo-normal:  inside ', count(d_bvh < 0._R8P), '  min d ', minval(d_bvh)
print '(A,I5,A,F10.6)', 'winding number: inside ', count(d_wn  < 0._R8P), '  min d ', minval(d_wn)
print '(A,I0,A,ES10.3)', 'signs differ at ', count(d_bvh * d_wn < 0._R8P), ' points, max ||d_pn| - |d_wn|| ', &
                        maxval(abs(abs(d_bvh) - abs(d_wn)))
endprogram tut05_batched
