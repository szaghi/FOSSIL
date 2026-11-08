use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, FEATURE_FACE, FEATURE_EDGE, &
                       FEATURE_VERTEX, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 12                    ! probes per direction
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object)      :: surface
type(rigid_pose_t)            :: pose
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d(:, :), cp(:, :, :), nm(:, :, :)
integer(I4P),     allocatable :: facet(:, :), feature(:, :), body(:, :)
integer(I4P)                  :: i, j, k, p

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%enter_device
allocate(points(n**3), d(n**3, 2), cp(3, n**3, 2), nm(3, n**3, 2), facet(n**3, 2), feature(n**3, 2), body(n**3, 2))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.55_R8P + 1.1_R8P * (i - 1) / (n - 1)) * ex_R8P + &
               (-0.55_R8P + 1.1_R8P * (j - 1) / (n - 1)) * ey_R8P + (-0.55_R8P + 1.1_R8P * (k - 1) / (n - 1)) * ez_R8P
enddo ; enddo ; enddo
call surface%query_many_device(points=points, distances=d(:, 1), closest_points=cp(:, :, 1), facets=facet(:, 1), &
                               features=feature(:, 1), normals=nm(:, :, 1), bodies=body(:, 1), &
                               sign_algorithm=SIGN_FAST_WINDING_NUMBER)                  ! column 1: device
call surface%query_many(points=points, distances=d(:, 2), closest_points=cp(:, :, 2), facets=facet(:, 2), &
                        features=feature(:, 2), normals=nm(:, :, 2), bodies=body(:, 2), &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)                         ! column 2: host
print '(A,I0,A)', 'query_many_device vs query_many, ', size(points), ' probes, mismatches:'
print '(A,I0,A,I0,A,I0)', '  distances ', count(d(:, 1) /= d(:, 2)), ', closest points ', &
      count(any(cp(:, :, 1) /= cp(:, :, 2), dim=1)), ', facets ', count(facet(:, 1) /= facet(:, 2))
print '(A,I0,A,I0,A,I0)', '  features ', count(feature(:, 1) /= feature(:, 2)), ', normals ', &
      count(any(nm(:, :, 1) /= nm(:, :, 2), dim=1)), ', bodies ', count(body(:, 1) /= body(:, 2))
print '(A,I0,A,I0,A,I0)', 'feet on a face ', count(feature(:, 1) == FEATURE_FACE), ', on an edge ', &
      count(feature(:, 1) == FEATURE_EDGE), ', on a vertex ', count(feature(:, 1) == FEATURE_VERTEX)
