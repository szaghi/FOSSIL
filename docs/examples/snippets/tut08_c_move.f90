!< Animation of tutorial 8: the ghost points of the grid plane y = 0 (41^3 grid of chapter 8), one per frame in turn
!< around the section of the placed dragon by that plane (seen along y). Each frame draws the ghost point (teal), its
!< wall point (gold), its image point (ember, 2 foot - x) and its probe two grid spacings out along the normal (grey),
!< with the feature of the wall point and the visibility of the image in the caption. Wall points are closest points
!< in 3-D: they may lie off the plane.
program tut08_c_move
use fossil,     only : surface_stl_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER, FEATURE_FACE, &
                       FEATURE_EDGE
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ez_R8P

implicit none
real(R8P),    parameter  :: H = 0.03_R8P                         ! the grid spacing of chapter 8
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(vector_R8P)         :: x(1), foot, image(1), probe(1)
type(vector_R8P), allocatable :: nodes(:), ghosts(:)
real(R8P),        allocatable :: dn(:), angle(:)
integer(I4P),     allocatable :: order(:)
real(R8P)                :: d, im(3, 1), pr(3, 1), uv(2)
integer(I4P)             :: k, s, u, n, feature
logical                  :: seen(1), seen_probe(1)
character(96)            :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 2_I4P, 0._R8P, section)                ! (u, v) = (z, x)
call fig_open(u, 'anim-tut08-move', '', -90._R8P, 0._R8P)
call fig_zoom(u, -0.04_R8P * ez_R8P, 0.22_R8P)
call fig_anim(u, 0.5_R8P)
do s = 1, section%get_vertices_number()
   if (section%get_next(s) == 0_I4P) cycle
   uv = section%get_vertex(s) ; foot = uv(2) * ex_R8P + uv(1) * ez_R8P
   uv = section%get_vertex(section%get_next(s))
   call fig_segment(u, 'body', foot, uv(2) * ex_R8P + uv(1) * ez_R8P)
enddo
n = 41                                                                     ! the plane y = 0 of the grid
allocate(nodes(n * n), dn(n * n))
do k = 1, n ; do s = 1, n
   nodes(s + (k - 1) * n) = (-0.6_R8P + (s - 1) * H) * ex_R8P + (-0.6_R8P + (k - 1) * H) * ez_R8P
enddo ; enddo
call surface%distance_many(points=nodes, distances=dn, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
ghosts = pack(nodes, dn < 0._R8P .and. dn > -H)
allocate(angle(size(ghosts)), order(size(ghosts)))
angle = atan2(ghosts%z - sum(ghosts%z) / size(ghosts), ghosts%x - sum(ghosts%x) / size(ghosts))
do k = 1, size(ghosts)                                                     ! around the section, by angle
   order(k) = count(angle < angle(k)) + 1
enddo
call fig_points(u, 'cell', ghosts)
do k = 1, size(ghosts)
   x(1) = ghosts(findloc(order, k, dim=1))
   call surface%query_point(point=x(1), distance=d, closest_point=foot, feature=feature, &
                            sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   call surface%image_points(points=x, images=im, visible=seen)
   call surface%probe_many(points=x, lengths=[2._R8P * H], probe_points=pr, visible=seen_probe)
   image(1)%x = im(1, 1) ; image(1)%y = im(2, 1) ; image(1)%z = im(3, 1)
   probe(1)%x = pr(1, 1) ; probe(1)%y = pr(2, 1) ; probe(1)%z = pr(3, 1)
   write(caption, '(A,I2,A,F8.4,A,A,A,A)') 'ghost ', k, ': d ', d, ', wall point on ', &
         trim(merge('a face  ', merge('an edge ', 'a vertex', feature == FEATURE_EDGE), feature == FEATURE_FACE)), &
         ', image ', trim(merge('visible', 'hidden ', seen(1)))
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call fig_segment(u, 'accent', x(1), foot)
   call fig_segment(u, 'outside', foot, image(1))
   call fig_segment(u, 'box', foot, probe(1))
   call fig_points(u, 'inside', x)
   call fig_points(u, 'foot', [foot])
   call fig_points(u, 'patch', image)
   call fig_points(u, 'box', probe)
enddo
call fig_close(u)
endprogram tut08_c_move
