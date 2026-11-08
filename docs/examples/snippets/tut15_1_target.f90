!< Tutorial 15, step 1: the projection target of the NACA0012 wing (patches, constraint curves, corners), then points
!< near the trailing edge classified and projected, free and constrained. Figures: the target on the wing seen from
!< behind; the mid-span plane near the trailing edge with the points, their free feet and their constrained feet.
program tut15_1_target
use fossil,     only : surface_stl_object, curve_object, projection_target_t, project_points, classify_points, &
                       slice_surface, PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: N = 6
type(surface_stl_object)  :: wing
type(projection_target_t) :: target
type(curve_object)        :: section
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                 :: px(N), py(N), pz(N), ffx(N), ffy(N), ffz(N), fx(N), fy(N), fz(N)
integer(I4P)              :: kinds(N), ids(N), free(N), i, k, u
character(8)              :: num
character(16), parameter  :: NAME(0:3) = [character(16) :: 'PROJ_FREE', 'PROJ_PATCH', 'PROJ_CURVE', 'PROJ_CORNER']

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
call wing%projection_target(target)                                          ! built once: patches, curves, corners
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'wing: ', wing%get_facets_number(), ' facets; target: patches ', target%np, &
      ', constraint curves ', target%nc, ' (', size(target%seg, 2), ' segments), corners ', target%nk
px = [0.735_R8P, 0.75_R8P, 0.757_R8P, 0.758_R8P, 0.762_R8P, 0.756_R8P]                  ! around the trailing edge, mid-span
py = [0.0075_R8P, 0.0063_R8P, 0.0002_R8P, 0.0040_R8P, -0.0005_R8P, -0.005_R8P]
pz = 0.5_R8P
call classify_points(target, px, py, pz, tolerance=0.01_R8P, kinds=kinds, ids=ids)
free = PROJ_FREE
call project_points(target, px, py, pz, free, ids, ffx, ffy, ffz)   ! anywhere on the wing
call project_points(target, px, py, pz, kinds, ids, fx, fy, fz)     ! on the constraint found
do i = 1, N
   print '(A,I0,A,F6.3,A,F7.4,A,A12,I3,A,F7.4,A,F8.5,A,F7.4,A,F8.5,A)', 'point ', i, ' (', px(i), ',', py(i), '): ', &
         NAME(kinds(i)), ids(i), '  free (', ffx(i), ',', ffy(i), ')  constrained (', fx(i), ',', fy(i), ')'
enddo

call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)  ! for figure 2

! figure 1: the target on the wing, drawn with the span horizontal, (x, y, z) -> (x, -z, y), seen from behind
allocate(a(size(target%seg, 2)), b(size(target%seg, 2)))
do k = 1, size(target%seg, 2)
   a(k) = turn(target%seg(1:3, k)) ; b(k) = turn(target%seg(4:6, k))
enddo
call wing%rotate(axis=vector_R8P(1._R8P, 0._R8P, 0._R8P), angle=acos(-1._R8P) / 2._R8P)
call fig_open(u, 'tut15-target', 'projection_target: 12 constraint curves (gold), 7 corners (teal)', -30._R8P, 30._R8P)
call fig_surface(u, 'body', wing)
call fig_segments(u, 'accent', a, b)
call fig_points(u, 'probe', [(turn(target%corner(:, k)), k = 1, target%nk)])
call fig_label(u, 'accent', turn([0.78_R8P, 0._R8P, 0.5_R8P]), 'trailing-edge lines')
call fig_label(u, 'accent', turn([0.25_R8P, 0.08_R8P, 1.02_R8P]), 'tip outline')
call fig_label(u, 'accent', turn([0.3_R8P, 0.07_R8P, 0._R8P]), 'root outline')
call fig_close(u)

! figure 2: the mid-span plane near the trailing edge, from +z: section, points, free and constrained feet
deallocate(a, b) ; allocate(a(section%get_vertices_number()), b(section%get_vertices_number()))
k = 0
do i = 1, section%get_vertices_number()
   if (section%get_next(i) == 0_I4P) cycle
   k = k + 1 ; a(k) = uv(section%get_vertex(i)) ; b(k) = uv(section%get_vertex(section%get_next(i)))
enddo
call fig_open(u, 'tut15-points', 'points 1-6 near the trailing edge: free feet (grey), constrained feet (yellow)', &
              -90._R8P, 90._R8P)
call fig_zoom(u, vector_R8P(0.749_R8P, 0.001_R8P, 0._R8P), 0.0135_R8P)
call fig_segments(u, 'body', a(1:k), b(1:k))
call fig_points(u, 'accent', [vector_R8P(0.75_R8P, 0.00126_R8P, 0._R8P), vector_R8P(0.75_R8P, -0.00126_R8P, 0._R8P)])
call fig_segments(u, 'box', [(vector_R8P(px(i), py(i), 0._R8P), i = 1, N)], [(vector_R8P(ffx(i), ffy(i), 0._R8P), i = 1, N)])
call fig_segments(u, 'accent', [(vector_R8P(px(i), py(i), 0._R8P), i = 1, N)], [(vector_R8P(fx(i), fy(i), 0._R8P), i = 1, N)])
call fig_points(u, 'probe', [(vector_R8P(px(i), py(i), 0._R8P), i = 1, N)])
do i = 1, N
   write(num, '(I0)') i
   call fig_label(u, 'probe', vector_R8P(px(i) + 0.0006_R8P, py(i) + 0.0004_R8P, 0._R8P), trim(num))
enddo
call fig_label(u, 'body', vector_R8P(0.7395_R8P, 0.0012_R8P, 0._R8P), 'upper side')
call fig_label(u, 'body', vector_R8P(0.7395_R8P, -0.0035_R8P, 0._R8P), 'lower side')
call fig_close(u)
contains
   function turn(p) result(v)
   !< (x, y, z) drawn with the span horizontal: (x, -z, y).
   real(R8P), intent(in) :: p(3) !< Point.
   type(vector_R8P)      :: v    !< Drawn point.

   v = vector_R8P(p(1), -p(3), p(2))
   endfunction turn

   function uv(q) result(v)
   !< A point (x, y) of the mid-span plane, drawn at z = 0.
   real(R8P), intent(in) :: q(2) !< In-plane point.
   type(vector_R8P)      :: v    !< Drawn point.

   v = vector_R8P(q(1), q(2), 0._R8P)
   endfunction uv

endprogram tut15_1_target
