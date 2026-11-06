!run tut08 tut08_ib_queries
!< Tutorial 8: immersed-boundary queries of the ghost points: wall point, normal, feature, body, image points; walls and
!< gaps thinner than a cell.
program tut08_ib_queries
!region ghosts
use fossil, only : surface_stl_object, SIGN_FAST_WINDING_NUMBER, FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX, &
                  THICKNESS_MEDIAL
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 41                         ! grid points per direction
type(surface_stl_object)      :: surface, wrap
type(vector_R8P), allocatable :: points(:), ghosts(:)
real(R8P),        allocatable :: d(:), foot(:,:), normal(:,:), image(:,:), thickness(:), gap(:), area(:)
integer(I4P),     allocatable :: feature(:), body(:), region(:)
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
!endregion ghosts
!region query
allocate(foot(3, m), normal(3, m), image(3, m), feature(m), body(m), ambiguous(m), visible(m))
deallocate(d) ; allocate(d(m))
call surface%query_many(points=ghosts, distances=d, closest_points=foot, features=feature, normals=normal, bodies=body, &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER, ambiguous=ambiguous)
print '(A,I0,A,I0,A,I0)', 'ghost points ', m, ', shells ', surface%get_shells_number(), ', bodies ', surface%get_bodies_number()
print '(A,I0,A,I0,A,I0)', 'feature: face ', count(feature == FEATURE_FACE), ', edge ', count(feature == FEATURE_EDGE), &
                          ', vertex ', count(feature == FEATURE_VERTEX)
print '(A,I0)', 'ambiguous foot points ', count(ambiguous)
print '(A,3F9.5)', 'ghost 1      ', ghosts(1)%x, ghosts(1)%y, ghosts(1)%z
print '(A,3F9.5)', '  wall point ', foot(:, 1)
print '(A,3F9.5)', '  normal     ', normal(:, 1)
!endregion query
!region images
call surface%image_points(points=ghosts, images=image, visible=visible)
print '(A,3F9.5)', '  image      ', image(:, 1)
print '(A,I0,A,I0)', 'images visible from their wall point ', count(visible), ' of ', m
!endregion images
!region thin
call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap)   ! one closed surface (chapter 11)
call wrap%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL)
call wrap%thin_regions(thickness=thickness, gap=gap, h=h, factor=1._R8P, region=region, region_area=area)
print '(A,ES9.2,A,ES9.2,A,F6.3)', 'wrap: thinnest wall ', minval(thickness), ', narrowest gap ', minval(gap), ', h ', h
print '(A,I0,A,I0,A,I0,A,F7.4)', 'under h: ', count(region > 0), ' of ', size(region), ' facets in ', size(area), &
                                 ' regions, the largest of area ', area(1)
!endregion thin
endprogram tut08_ib_queries
