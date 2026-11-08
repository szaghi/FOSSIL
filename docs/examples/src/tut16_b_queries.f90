!run tut16b tut16_b_queries
!fig tut16-queries
!< Tutorial 16, step 2: the rich point queries on the device: query_many_device (distance, closest point, facet,
!< feature, pseudo-normal, body), wall_distance_many_device with a cutoff, and a posed query, each against its host
!< twin. The figure: probes on a circle about the dragon, joined to the closest points query_many_device returns,
!< the feet coloured by the feature they lie on. Built with gfortran the device kernels run on the host.
program tut16_b_queries
!region rich
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
!endregion rich
!region wall
call surface%wall_distance_many_device(points=points, distances=d(:, 1), bodies=body(:, 1), cutoff=0.05_R8P)
call surface%wall_distance_many(points=points, distances=d(:, 2), bodies=body(:, 2), cutoff=0.05_R8P)
print '(A,I0,A,I0,A,I0)', 'wall distance, cutoff 0.05: mismatches ', count(d(:, 1) /= d(:, 2)), &
      ', probes nearer than the cutoff ', count(body(:, 1) /= 0), ' of ', size(points)
!endregion wall
!region posed
pose = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call surface%query_many_device(points=points, distances=d(:, 1), facets=facet(:, 1), pose=pose)  ! chapter 11's step
call surface%query_many(points=points, distances=d(:, 2), facets=facet(:, 2), pose=pose)
print '(A,I0,A,I0)', 'posed query: distance mismatches ', count(d(:, 1) /= d(:, 2)), ', facet mismatches ', &
      count(facet(:, 1) /= facet(:, 2))
call surface%exit_device
!endregion posed
call draw_figure
contains
   subroutine draw_figure
   !< 72 probes (grey) on the circle of radius 0.45 in the plane z = 0, joined to their device closest points; feet on
   !< a face in yellow, on an edge in teal, on a vertex in ember.
   use fossil_fig, only : fig_open, fig_surface, fig_points, fig_segments, fig_close
   integer(I4P), parameter       :: m = 72
   type(vector_R8P)              :: probe(m), foot(m)
   real(R8P)                     :: dd(m), cc(3, m)
   integer(I4P)                  :: ft(m), u, q
   character(100)                :: caption

   do q = 1, m
      probe(q) = 0.45_R8P * (cos(2._R8P * pi * (q - 1) / m) * ex_R8P + sin(2._R8P * pi * (q - 1) / m) * ey_R8P)
   enddo
   call surface%enter_device
   call surface%query_many_device(points=probe, distances=dd, closest_points=cc, features=ft)
   call surface%exit_device
   do q = 1, m
      foot(q) = cc(1, q) * ex_R8P + cc(2, q) * ey_R8P + cc(3, q) * ez_R8P
   enddo
   write(caption, '(A,I0,A,I0,A,I0,A)') 'query_many_device feet: face ', count(ft == FEATURE_FACE), &
         ' (yellow), edge ', count(ft == FEATURE_EDGE), ' (teal), vertex ', count(ft == FEATURE_VERTEX), ' (ember)'
   call fig_open(u, 'tut16-queries', trim(caption), -90._R8P, 60._R8P)
   call fig_surface(u, 'body', surface)
   call fig_segments(u, 'box', probe, foot)
   call fig_points(u, 'box', probe)
   if (count(ft == FEATURE_FACE) > 0) call fig_points(u, 'foot', pack(foot, ft == FEATURE_FACE))
   if (count(ft == FEATURE_EDGE) > 0) call fig_points(u, 'inside', pack(foot, ft == FEATURE_EDGE))
   if (count(ft == FEATURE_VERTEX) > 0) call fig_points(u, 'outside', pack(foot, ft == FEATURE_VERTEX))
   call fig_close(u)
   endsubroutine draw_figure
endprogram tut16_b_queries
