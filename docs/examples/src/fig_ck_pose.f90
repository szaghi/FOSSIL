!run figck-pose fig_ck_pose
!fig ck-pose
!< Figure of the recipe "Query a rotated body without rebuilding it": the unit cube at the recipe's pose (45 degrees
!< about the vertical axis through its centre), the cube at rest as a translucent ghost, the recipe's world-frame
!< query point (teal) and the closest point query_many returns for the cube at rest (yellow; the one with pose= is
!< 0.005 from the query point, under it at this scale: the caption gives both distances). Seen from above with the top
!< faces (z = 1) left out, so the point at z = 0.5 shows; the walls are edge-on, the bottoms (z = 0) seen.
program fig_ck_pose
use fossil,     only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube, moved
type(rigid_pose_t)       :: pose
type(vector_R8P)         :: p(1), rest(1)
real(R8P)                :: d(1), dr(1), cp(3, 1)
character(80)            :: caption
integer(I4P)             :: u

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
pose = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 4)          ! the recipe's pose
pose%t = [0.5_R8P, 0.5_R8P, 0._R8P] - matmul(pose%r, [0.5_R8P, 0.5_R8P, 0._R8P])
p = vector_R8P(1.2_R8P, 0.5_R8P, 0.5_R8P)
call cube%query_many(points=p, distances=d, pose=pose)
call cube%query_many(points=p, distances=dr, closest_points=cp)
rest(1)%x = cp(1, 1) ; rest(1)%y = cp(2, 1) ; rest(1)%z = cp(3, 1)
moved = cube                                                                       ! only to draw it at the pose
call moved%move(pose)
call lidless(moved) ; call lidless(cube)
write(caption, '(A,F7.4,A,F7.4,A)') 'signed distance ', d(1), ' with pose= (stone), ', dr(1), ' at rest (ghost)'
call fig_open(u, 'ck-pose', trim(caption), -90._R8P, 90._R8P)
call fig_surface(u, 'body', moved)
call fig_surface(u, 'ghost', cube)
call fig_segment(u, 'foot', p(1), rest(1))
call fig_points(u, 'foot', rest)
call fig_points(u, 'probe', p)
call fig_close(u)
contains
   subroutine lidless(surface)
   !< Keep only the facets of a surface not facing up (normal z < 0.5).
   type(surface_stl_object), intent(inout) :: surface !< Surface.
   type(facet_object), pointer             :: f
   type(facet_object), allocatable         :: kept(:)
   integer(I4P)                            :: i, n

   n = 0
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i)
      if (f%normal%z <= 0.5_R8P) n = n + 1
   enddo
   allocate(kept(n)) ; n = 0
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i)
      if (f%normal%z > 0.5_R8P) cycle
      n = n + 1 ; kept(n) = f
   enddo
   call surface%adopt_facets(facets=kept)
   endsubroutine lidless
endprogram fig_ck_pose
