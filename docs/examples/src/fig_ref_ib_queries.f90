!run figref-ib-queries fig_ref_ib_queries
!fig ref-ib-queries-images
!fig ref-ib-queries-ambiguous
!fig ref-ib-queries-curvature
!fig ref-ib-queries-markers
!< Figures of the reference page "Immersed-boundary queries":
!< - two unit cubes 0.2 apart, a row of ghost points inside the first one at depths 0.01 to 0.39 from its face x = 1,
!<   each mirrored through its wall point: the images within the gap are visible, those inside the second cube not;
!< - the points of the plane z = 0.5 that `query_many` flags as ambiguous (two foot points), joined to both feet: the
!<   medial axis of the two cubes and of the gap between them;
!< - an ellipsoid (marching cubes) coloured by the mean curvature H at the foot of each facet centroid;
!< - grid-adapted markers on the unit cube, dx = 0.1, features kept.
program fig_ref_ib_queries
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segment, fig_points, fig_triangle, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 20, NS = 57
type(surface_stl_object)        :: a, b, ell
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: fp
type(vector_R8P)                :: ghost(NG), cen(1), sheet(NS * NS)
type(vector_R8P), allocatable   :: pts(:)
real(R8P)                       :: img(3, NG), foot(3, NG), dg(NG), ds(NS * NS), second(3, NS * NS), fs(3, NS * NS)
real(R8P)                       :: phi(41, 41, 41), d1(1), hk(2, 1), hmin, hmax, fraction
real(R8P), allocatable          :: h(:), mc(:,:), ma(:), mn(:,:)
logical                         :: vis(NG), amb(NS * NS)
integer(I4P)                    :: i, j, k, u, status
character(100)                  :: caption

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
b = a
call b%translate(x=1.2_R8P)
call a%merge_solids(b)
call a%analyze

! ---- images and their visibility
do i = 1, NG
   ghost(i) = vector_R8P(0.99_R8P - 0.02_R8P * (i - 1), 0.04_R8P + 0.025_R8P * i, 0.5_R8P)
enddo
call a%image_points(points=ghost, images=img, visible=vis, foot_points=foot)
print '(A,I0,A,I0,A)', 'images: ', count(vis), ' visible, ', count(.not. vis), ' inside the other cube'
write(caption, '(A,I0,A,I0,A)') 'images of the ghost points (teal): ', count(vis), ' visible, ', &
                                count(.not. vis), ' behind a wall (ember)'
call fig_open(u, 'ref-ib-queries-images', trim(caption), -90._R8P, 90._R8P, height=420_I4P)
call fig_zoom(u, vector_R8P(1.1_R8P, 0.35_R8P, 0.5_R8P), 0.5_R8P)
call fig_surface(u, 'ghost', a)
call fig_points(u, 'probe', ghost)
do i = 1, NG
   if (vis(i)) then
      call fig_segment(u, 'foot', v(foot(:, i)), v(img(:, i)))
      call fig_points(u, 'foot', [v(img(:, i))])
   else
      call fig_segment(u, 'outside', v(foot(:, i)), v(img(:, i)))
      call fig_points(u, 'outside', [v(img(:, i))])
   endif
   call fig_segment(u, 'probe', ghost(i), v(foot(:, i)))
enddo
call fig_close(u)

! ---- ambiguous foot points on the plane z = 0.5
do j = 1, NS ; do i = 1, NS
   sheet(i + (j - 1) * NS) = vector_R8P(-0.3_R8P + 0.05_R8P * (i - 1), -0.9_R8P + 0.05_R8P * (j - 1), 0.5_R8P)
enddo ; enddo
call a%query_many(points=sheet, distances=ds, closest_points=fs, ambiguous=amb, second_foot_points=second)
print '(A,I0,A,I0)', 'ambiguous: ', count(amb), ' of ', size(amb)
write(caption, '(A,I0,A,I0,A)') 'plane z = 0.5: ', count(amb), ' of ', size(amb), &
                                ' points with two foot points, joined to both'
call fig_open(u, 'ref-ib-queries-ambiguous', trim(caption), -90._R8P, 90._R8P, height=420_I4P)
call fig_surface(u, 'ghost', a)
do i = 1, size(amb)
   if (.not. amb(i)) cycle
   call fig_segment(u, 'foot', sheet(i), v(fs(:, i)))
   call fig_segment(u, 'foot', sheet(i), v(second(:, i)))
enddo
pts = pack(sheet, amb)
call fig_points(u, 'outside', pts)
call fig_close(u)

! ---- foot curvature on an ellipsoid with semi-axes 1, 0.7, 0.45
do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41
   phi(i, j, k) = norm2((-1.2_R8P + 0.06_R8P * ([i, j, k] - 1)) / [1._R8P, 0.7_R8P, 0.45_R8P]) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.2_R8P, -1.2_R8P, -1.2_R8P), &
                        bmax=vector_R8P(1.2_R8P, 1.2_R8P, 1.2_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call ell%adopt_facets(facets=facets)
deallocate(pts) ; allocate(h(ell%get_facets_number()), pts(ell%get_facets_number()))
do i = 1, ell%get_facets_number()
   fp => ell%facet_at(i)
   pts(i) = (fp%vertex(1) + fp%vertex(2) + fp%vertex(3)) / 3._R8P
enddo
block
   real(R8P), allocatable :: dd(:), hh(:,:)
   allocate(dd(size(pts)), hh(2, size(pts)))
   call ell%query_many(points=pts, distances=dd, curvatures=hh)
   h = hh(1, :)
endblock
hmin = minval(h) ; hmax = maxval(h)
print '(A,I0,A,F6.3,A,F6.3)', 'ellipsoid: ', size(h), ' facets, H from ', hmin, ' to ', hmax
write(caption, '(A,F4.2,A,F4.2,A)') 'ellipsoid 1 x 0.7 x 0.45: H at the foot, ', hmin, ' (teal) to ', hmax, &
                                    ' (ember)'
call fig_open(u, 'ref-ib-queries-curvature', trim(caption), -60._R8P, 25._R8P, height=380_I4P)
do i = 1, ell%get_facets_number()
   fp => ell%facet_at(i)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), &
                     value=2._R8P * (h(i) - hmin) / (hmax - hmin) - 1._R8P)
enddo
call fig_close(u)

! ---- markers on the unit cube
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%surface_markers(dx=0.1_R8P, centroids=mc, areas=ma, normals=mn, spacing_fraction=fraction, &
                       preserve_features=.true.)
deallocate(pts) ; allocate(pts(size(ma)))
do i = 1, size(ma)
   pts(i) = v(mc(:, i) + 0.01_R8P * mn(:, i))                   ! lifted off the face to draw
enddo
print '(A,I0,A,F5.3)', 'markers: ', size(ma), ', spacing fraction ', fraction
write(caption, '(A,I0,A,F4.1,A)') 'surface_markers, dx = 0.1: ', size(ma), ' markers, ', 100._R8P * fraction, &
                                  ' % of spacings in [0.8, 1.2] dx'
call fig_open(u, 'ref-ib-queries-markers', trim(caption), -60._R8P, 25._R8P, height=380_I4P)
do k = 0, 1 ; do j = 0, 19 ; do i = 0, 19                       ! the cube's faces as small quads, so that the
   call face_quad(1, k, i, j) ; call face_quad(2, k, i, j) ; call face_quad(3, k, i, j)  ! markers sort right
enddo ; enddo ; enddo
call fig_points(u, 'probe', pts)
call fig_close(u)
contains
   subroutine face_quad(axis, side, i, j)
   !< Quad (i, j) of a 20 x 20 split of the cube face x(axis) = side.
   integer(I4P), intent(in) :: axis, side, i, j !< Face, quad.
   real(R8P)                :: c(3, 4)
   integer(I4P)             :: ua, va

   ua = mod(axis, 3) + 1 ; va = mod(axis + 1, 3) + 1
   c(axis, :) = real(side, R8P)
   c(ua, :) = 0.05_R8P * [i, i + 1, i + 1, i] ; c(va, :) = 0.05_R8P * [j, j, j + 1, j + 1]
   call fig_quad(u, 'body', v(c(:, 1)), v(c(:, 2)), v(c(:, 3)), v(c(:, 4)))
   endsubroutine face_quad

   function v(c) result(p)
   real(R8P), intent(in) :: c(3) !< Coordinates.
   type(vector_R8P)      :: p    !< Point.

   p = vector_R8P(c(1), c(2), c(3))
   endfunction v
endprogram fig_ref_ib_queries
