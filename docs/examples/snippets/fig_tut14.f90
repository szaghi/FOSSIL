!< Figures of tutorial 14. First the pebble thrown through the dragon within one step: the pebble at the five sampled
!< instants (ghost, every one apart from the dragon), its path, and at the time of first contact time_of_impact
!< returns (solid) with the witness point on the dragon (yellow). Then the 2-D step: the NACA 0012 section at rest and
!< three particles on their paths, drawn at s = 0 (ghost outlines) and at the first contact first_impact_2d returns
!< (the pair in contact in yellow, the others in teal), with the witness point on the section.
program fig_tut14
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                       pose_from_axis_angle, pose_2d, pose_interpolate, first_impact_2d, slice_surface, to_world
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segment, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter          :: PI = acos(-1._R8P)
type(surface_stl_object)      :: dragon, pebble, moved, naca
type(curve_object)            :: wing, disc, curves(4)
type(rigid_pose_t)            :: p0, p1
type(rigid_pose_2d_t)         :: q0(4), q1(4)
type(vector_R8P)              :: bmin, bmax, foot(1)
real(R8P)                     :: scale, toi, pa(3), qa(2)
character(120)                :: caption
integer(I4P)                  :: u, k, status, pair(2)

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! the dragon, placed as in chapter 3
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze
call pebble%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call pebble%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
call pebble%resize(x=0.04_R8P, y=0.04_R8P, z=0.04_R8P)
call pebble%analyze
p0 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi, status, point_self=pa)
print '(A,F9.6)', 'pebble: first contact at s = ', toi

write(caption, '(A,F6.4,A)') 'time_of_impact: contact at s = ', toi, ' (ember); samples every 0.25 apart'
call fig_open(u, 'tut14-throw', trim(caption), -60._R8P, 15._R8P)
call fig_zoom(u, vector_R8P(0._R8P, -0.1_R8P, -0.05_R8P), 0.55_R8P)
call fig_surface(u, 'body', dragon)
call fig_segment(u, 'box', vector_R8P(-0.6_R8P, -0.3_R8P, -0.2_R8P), vector_R8P(0.6_R8P, -0.3_R8P, -0.2_R8P))
do k = 0, 4
   moved = pebble ; call moved%move(pose_interpolate(p0, p1, 0.25_R8P * k))
   call fig_surface(u, 'ghost', moved)
enddo
moved = pebble ; call moved%move(pose_interpolate(p0, p1, toi))
call fig_surface(u, 'patch', moved)
foot(1) = vector_R8P(pa(1), pa(2), pa(3))
call fig_points(u, 'foot', foot)
call fig_close(u)

call naca%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)  ! the 2-D step of the chapter
call slice_surface(naca, 3_I4P, 0.5_R8P, wing)
call make_disc(0.03_R8P, disc)
curves(1) = wing ; curves(2:4) = disc
q0(1) = pose_2d(angle=0._R8P) ; q1(1) = q0(1)
q0(2) = pose_2d(angle=0._R8P, translation=[-0.2_R8P, 0.25_R8P]) ; q1(2) = pose_2d(angle=PI, translation=[0.4_R8P, -0.25_R8P])
q0(3) = pose_2d(angle=0._R8P, translation=[-0.4_R8P, -0.12_R8P]) ; q1(3) = pose_2d(angle=PI, translation=[0.9_R8P, -0.12_R8P])
q0(4) = pose_2d(angle=0._R8P, translation=[-0.45_R8P, 0.2_R8P]) ; q1(4) = pose_2d(angle=PI, translation=[-0.45_R8P, -0.2_R8P])
call first_impact_2d(curves, q0, q1, toi, pair, status)
call wing%time_of_impact(disc, q0(1), q1(1), q0(pair(2)), q1(pair(2)), toi, status, point_self=qa)
print '(A,F9.6,A,2I2)', 'particles: first contact at s = ', toi, ', curves', pair

write(caption, '(A,F6.4,A,I0,A,I0)') 'first_impact_2d: contact at s = ', toi, ', curves ', pair(1), ' and ', pair(2)
call fig_open(u, 'tut14-particles', trim(caption), -90._R8P, 90._R8P)
call draw_curve(u, 'body', wing, q0(1))
do k = 2, 4
   call fig_segment(u, 'box', vector_R8P(q0(k)%t(1), q0(k)%t(2), 0._R8P), vector_R8P(q1(k)%t(1), q1(k)%t(2), 0._R8P))
   call draw_curve(u, 'ghost', disc, q0(k))
   call draw_curve(u, merge('accent', 'probe ', k == pair(2)), disc, pose_interpolate(q0(k), q1(k), toi))
enddo
foot(1) = vector_R8P(qa(1), qa(2), 0.001_R8P)
call fig_points(u, 'foot', foot)
call fig_close(u)
contains
   subroutine make_disc(radius, c)
   !< A regular 32-gon of the given radius about the origin, counter-clockwise (solid inside).
   real(R8P),          intent(in)    :: radius
   type(curve_object), intent(inout) :: c
   real(R8P)                         :: a(32)
   integer(I4P)                      :: m

   a = [(2._R8P * PI * m / 32._R8P, m = 0, 31)]
   call c%adopt_loops(u=radius * cos(a), v=radius * sin(a), loop_sizes=[32_I4P])
   endsubroutine make_disc

   subroutine draw_curve(unit, style, c, pose)
   !< The segments of a curve placed at a pose, in the plane z = 0.
   integer(I4P),          intent(in) :: unit
   character(*),          intent(in) :: style
   type(curve_object),    intent(in) :: c
   type(rigid_pose_2d_t), intent(in) :: pose
   type(vector_R8P), allocatable     :: a(:), b(:)
   real(R8P)                         :: x(2), y(2)
   integer(I4P)                      :: i, n

   allocate(a(c%get_vertices_number()), b(c%get_vertices_number()))
   n = 0
   do i = 1, c%get_vertices_number()
      if (c%get_next(i) == 0) cycle
      x = to_world(pose, c%get_vertex(i)) ; y = to_world(pose, c%get_vertex(c%get_next(i)))
      n = n + 1 ; a(n) = vector_R8P(x(1), x(2), 0._R8P) ; b(n) = vector_R8P(y(1), y(2), 0._R8P)
   enddo
   call fig_segments(unit, trim(style), a(1:n), b(1:n))
   endsubroutine draw_curve
endprogram fig_tut14
