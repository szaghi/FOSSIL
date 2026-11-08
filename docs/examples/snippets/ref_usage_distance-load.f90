use fossil, only : surface_stl_object, SIGN_RAY_INTERSECTIONS, AABB_TREE_OCTREE, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface, dragon, dragon_octree
type(vector_R8P)         :: p, closest, normal
type(vector_R8P)         :: points(21)
real(R8P)                :: d2, d, sd, sd_ray, dist, d_bvh(21), d_oct(21)
logical                  :: inside, same
integer(I4P)             :: facet, edge, vertex, i, status

call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! the unit cube [0, 1]^3
