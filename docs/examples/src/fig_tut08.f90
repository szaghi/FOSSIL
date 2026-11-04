!run fig08 fig_tut08
!fig tut08-images
!fig tut08-features
!< Figures of tutorial 8: the ghost points of the 41^3 grid, their wall points and image points, and the feature
!< (face, edge, vertex) each wall point lies on. The dragon is placed as in tutorials 2 and 3 (saved and read back as
!< dragon-placed.stl is); ghosts, query_many and image_points are those of tutorial 8. The first figure shows the
!< ghosts of the grid plane y = 0 and the section of the dragon by that plane, seen along y; the second the features
!< on cube-01.stl, whose regions are clear: a sheet of points over the top face, each joined to its wall point.
program fig_tut08
use fossil,     only : surface_stl_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER, FEATURE_FACE, &
                      FEATURE_EDGE, FEATURE_VERTEX
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
integer(I4P), parameter       :: n = 41
type(surface_stl_object)      :: surface, cube
type(curve_object)            :: section
type(vector_R8P), allocatable :: points(:), ghosts(:), g(:), w(:), im(:)
real(R8P),        allocatable :: d(:), foot(:,:), normal(:,:), image(:,:)
integer(I4P),     allocatable :: feature(:)
logical,          allocatable :: plane(:), visible(:)
type(vector_R8P)              :: bmin, bmax, c
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: h, scale, uv(2), r
integer(I4P)                  :: u, i, j, k, p, m

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)        ! dragon-placed.stl, as tutorials 2-3
call surface%sanitize(fill_holes=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig08-dragon-placed.stl')
call surface%load_from_file(file_name='fig08-dragon-placed.stl', guess_format=.true.)

h = 1.2_R8P / (n - 1)                                                             ! the ghosts of tutorial 8
allocate(points(n**3), d(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call surface%distance_many(points=points, distances=d, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
ghosts = pack(points, d < 0._R8P .and. d > -h)
m = size(ghosts)
allocate(foot(3, m), normal(3, m), image(3, m), feature(m), visible(m))
deallocate(d) ; allocate(d(m))
call surface%query_many(points=ghosts, distances=d, closest_points=foot, features=feature, normals=normal, &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
call surface%image_points(points=ghosts, images=image, visible=visible)

plane = abs(ghosts%y) < 0.5_R8P * h                                               ! the grid plane y = 0
g = pack(ghosts, plane)
w = pack(foot(1, :) * ex_R8P + foot(2, :) * ey_R8P + foot(3, :) * ez_R8P, plane)
im = pack(image(1, :) * ex_R8P + image(2, :) * ey_R8P + image(3, :) * ez_R8P, plane)
call slice_surface(surface, 2_I4P, 0._R8P, section, open_loops=p, status=i) ! the wall in the plane: (u, v) = (z, x)
allocate(a(section%get_vertices_number()), b(section%get_vertices_number()))
do i = 1, section%get_vertices_number()
   uv = section%get_vertex(i) ; a(i) = uv(2) * ex_R8P + uv(1) * ez_R8P
   uv = section%get_vertex(section%get_next(i)) ; b(i) = uv(2) * ex_R8P + uv(1) * ez_R8P
enddo
k = 1                                                                             ! frame the largest loop
do i = 2, section%get_loops_number()
   if (abs(section%get_loop_area(i)) > abs(section%get_loop_area(k))) k = i
enddo
call section%get_loop_range(k, i, j)
c%x = 0.5_R8P * (minval(a(i:j)%x) + maxval(a(i:j)%x)) ; c%y = 0._R8P
c%z = 0.5_R8P * (minval(a(i:j)%z) + maxval(a(i:j)%z))
r = 0.5_R8P * max(maxval(a(i:j)%x) - minval(a(i:j)%x), maxval(a(i:j)%z) - minval(a(i:j)%z)) + 2._R8P * h
print '(A,I0,A,I0,A,I0)', 'ghosts ', m, ', on the plane y = 0: ', size(g), ', open loops ', p
call fig_open(u, 'tut08-images', 'plane y = 0: ghosts (teal), wall points (yellow), images (ember)', -90._R8P, 0._R8P)
call fig_zoom(u, c, r)
call fig_segments(u, 'box', a, b)
call fig_segments(u, 'cell', g, im)
call fig_points(u, 'inside', g)
call fig_points(u, 'foot', w)
call fig_points(u, 'outside', im)
call fig_close(u)

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! the features: a sheet of points over a cube
bmin = cube%get_bmin() ; bmax = cube%get_bmax() ; h = (bmax%x - bmin%x) / 8._R8P
deallocate(points) ; allocate(points(15**2))
p = 0
do j = 1, 15 ; do i = 1, 15
   p = p + 1
   points(p) = (bmin%x + (i - 4.3_R8P) * h) * ex_R8P + (bmin%y + (j - 4.3_R8P) * h) * ey_R8P + (bmax%z + 1.5_R8P * h) * ez_R8P
enddo ; enddo
deallocate(d, foot, feature) ; allocate(d(p), foot(3, p), feature(p))
call cube%query_many(points=points, distances=d, closest_points=foot, features=feature)
w = foot(1, :) * ex_R8P + foot(2, :) * ey_R8P + foot(3, :) * ez_R8P
print '(A,I0,A,I0,A,I0)', 'cube: face ', count(feature == FEATURE_FACE), ', edge ', count(feature == FEATURE_EDGE), &
                          ', vertex ', count(feature == FEATURE_VERTEX)
call fig_open(u, 'tut08-features', 'query_many features: face (yellow), edge (teal), vertex (ember)', -60._R8P, 62._R8P)
call fig_surface(u, 'body', cube)
call fig_segments(u, 'box', points, w)
call fig_points(u, 'foot', pack(points, feature == FEATURE_FACE))
call fig_points(u, 'probe', pack(points, feature == FEATURE_EDGE))
call fig_points(u, 'outside', pack(points, feature == FEATURE_VERTEX))
call fig_close(u)
endprogram fig_tut08
