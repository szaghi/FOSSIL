use fossil, only : surface_stl_object, AABB_TREE_OCTREE, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, &
                   SIGN_RAY_INTERSECTIONS, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                    ! points per direction
type(surface_stl_object)      :: surface, octree
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_pn(:), d_fwn(:), d_host(:)
integer(I4P)                  :: i, j, k, p, s, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)  ! SAH BVH: the device tree
call surface%enter_device                                   ! copy the flat BVH arrays to the device, once
call surface%enter_device                                   ! a second call is a no-op
print '(A,L1)', 'on the device: ', surface%is_on_device
