use fossil, only : surface_stl_object, SIGN_FAST_WINDING_NUMBER, FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 41                         ! grid points per direction
type(surface_stl_object)      :: surface
type(vector_R8P), allocatable :: points(:), ghosts(:)
real(R8P),        allocatable :: d(:), foot(:,:), normal(:,:), image(:,:)
integer(I4P),     allocatable :: feature(:), body(:)
logical,          allocatable :: ambiguous(:), visible(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p, m

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / (n - 1)
allocate(points(n**3), d(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call surface%distance_many(points=points, distances=d, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
ghosts = pack(points, d < 0._R8P .and. d > -h)                  ! inside, within one spacing of the wall
m = size(ghosts)
