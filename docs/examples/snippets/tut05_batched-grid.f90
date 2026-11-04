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
