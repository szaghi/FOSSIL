!run fig14 fig_tut14
!fig tut14-gap
!fig tut14-throw
!fig tut14-several
!fig tut14-particles
!fig tut14-sliding
!< Figures of tutorial 14 (self-contained: the dragon is placed here as in the chapter). First the gap d(s) between the
!< dragon and the pebble over the step, sampled at 300 instants (yellow), the instants of contact on the s axis (ember),
!< the five samples of run 1 (teal) and the toi of time_of_impact. Then the pebble at the five samples (ghost) and at
!< the toi (solid) with the witness point. Then the dragon and three pebbles: paths, pebbles at s = 0 (ghost) and at
!< the first contact of first_impact (the pair in ember). Then the 2-D step around the NACA 0012 section. Last, a
!< particle sliding along a wall: the iterates of the plain bound (top) and of the directional bound (bottom).
program fig_tut14
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                       pose_from_axis_angle, pose_2d, pose_interpolate, first_impact, first_impact_2d, slice_surface, &
                       to_world
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ez_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), SX = 1.5_R8P, SY = 1.5_R8P   ! graph scales
integer(I4P), parameter       :: NS = 300
type(surface_stl_object)      :: dragon, pebble, moved, naca, bodies(4)
type(curve_object)            :: wing, disc, wall, curves(4)
type(rigid_pose_t)            :: p0, p1, poses0(4), poses1(4)
type(rigid_pose_2d_t)         :: q0(4), q1(4), w0, w1
type(vector_R8P)              :: bmin, bmax, foot(1), a(NS), b(NS), samples(5)
real(R8P)                     :: scale, toi, pa(3), qa(2), s, d(0:NS), gap, lift
logical                       :: hit(0:NS)
character(120)                :: caption
character(40)                 :: text
integer(I4P)                  :: u, k, m, status, pair(2), steps, panel

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

! the gap over the step: a graph in the plane z = 0, s along x, d along y
do k = 0, NS
   s = real(k, R8P) / NS
   call dragon%distance_to(pebble, distance=d(k), intersecting=hit(k), pose_other=pose_interpolate(p0, p1, s))
enddo
call fig_open(u, 'tut14-gap', 'the gap over the step: the 5 samples (teal) miss both contacts (ember)', &
              -90._R8P, 90._R8P, height=300_I4P)
call fig_segment(u, 'box', vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(SX, 0._R8P, 0._R8P))
call fig_segment(u, 'box', vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(0._R8P, SY * 0.45_R8P, 0._R8P))
do k = 1, NS
   a(k) = vector_R8P(SX * (k - 1) / NS, SY * d(k - 1), 0._R8P) ; b(k) = vector_R8P(SX * k / NS, SY * d(k), 0._R8P)
enddo
call fig_segments(u, 'accent', a, b)
m = 0
do k = 1, NS
   if (.not. (hit(k - 1) .and. hit(k))) cycle
   m = m + 1 ; a(m) = vector_R8P(SX * (k - 1) / NS, 0._R8P, 0._R8P) ; b(m) = vector_R8P(SX * k / NS, 0._R8P, 0._R8P)
enddo
call fig_segments(u, 'hole', a(1:m), b(1:m))
do k = 1, 5
   samples(k) = vector_R8P(SX * 0.25_R8P * (k - 1), SY * d(NS * (k - 1) / 4), 0._R8P)
   write(text, '(A,F5.3)') 'gap ', d(NS * (k - 1) / 4)
   call fig_label(u, 'probe', samples(k) - merge(0.24_R8P, 0._R8P, k == 5) * ex_R8P, trim(text))
enddo
call fig_points(u, 'probe', samples)
foot(1) = vector_R8P(SX * toi, 0._R8P, 0._R8P)
call fig_points(u, 'foot', foot)
write(text, '(A,F6.4)') 'toi ', toi
call fig_label(u, 'foot', vector_R8P(SX * toi - 0.12_R8P, -0.05_R8P, 0._R8P), trim(text))
call fig_label(u, 'box', vector_R8P(-0.02_R8P, -0.05_R8P, 0._R8P), 's = 0')
call fig_label(u, 'box', vector_R8P(SX - 0.08_R8P, -0.05_R8P, 0._R8P), 's = 1')
call fig_label(u, 'box', vector_R8P(0.01_R8P, SY * 0.45_R8P, 0._R8P), 'gap')
call fig_close(u)

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

bodies(1) = dragon ; poses0(1) = pose_identity() ; poses1(1) = poses0(1)       ! several bodies (run 3)
bodies(2:4) = pebble
poses0(2) = p0 ; poses1(2) = p1
poses0(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0._R8P, 0.3_R8P, 0.6_R8P])
poses1(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI, translation=[0._R8P, 0.3_R8P, -0.6_R8P])
poses0(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-0.15_R8P, 0.3_R8P, 0.45_R8P])
poses1(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI, translation=[1.05_R8P, 0.3_R8P, 0.45_R8P])
call first_impact(bodies, poses0, poses1, toi, pair, status, tolerance=1.e-6_R8P)
print '(A,F9.6,A,2I2)', 'several: first contact at s = ', toi, ', bodies', pair
write(caption, '(A,I0,A,I0,A,F6.4,A)') 'first_impact: bodies ', pair(1), ' and ', pair(2), ' touch first, at s = ', &
                                       toi, ' (ember)'
call fig_open(u, 'tut14-several', trim(caption), -60._R8P, 15._R8P)
call fig_surface(u, 'body', dragon)
call fig_label(u, 'body', vector_R8P(0.1_R8P, -0.1_R8P, 0.35_R8P), '1 (at rest)')
do k = 2, 4
   call fig_segment(u, 'box', vector_R8P(poses0(k)%t(1), poses0(k)%t(2), poses0(k)%t(3)), &
                    vector_R8P(poses1(k)%t(1), poses1(k)%t(2), poses1(k)%t(3)))
   moved = pebble ; call moved%move(poses0(k)) ; call fig_surface(u, 'ghost', moved)
   write(text, '(I0,A)') k, ' at s = 0'
   call fig_label(u, 'box', vector_R8P(poses0(k)%t(1), poses0(k)%t(2), poses0(k)%t(3) + 0.04_R8P), trim(text))
   moved = pebble ; call moved%move(pose_interpolate(poses0(k), poses1(k), toi))
   call fig_surface(u, merge('patch', 'probe', any(pair == k)), moved)
enddo
call fig_close(u)

call naca%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)  ! the 2-D step (run 4)
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
call draw_curve(u, 'body', wing, q0(1), 0._R8P)
do k = 2, 4
   call fig_segment(u, 'box', vector_R8P(q0(k)%t(1), q0(k)%t(2), 0._R8P), vector_R8P(q1(k)%t(1), q1(k)%t(2), 0._R8P))
   call draw_curve(u, 'ghost', disc, q0(k), 0._R8P)
   call draw_curve(u, merge('accent', 'probe ', k == pair(2)), disc, pose_interpolate(q0(k), q1(k), toi), 0._R8P)
enddo
foot(1) = vector_R8P(qa(1), qa(2), 0.001_R8P)
call fig_points(u, 'foot', foot)
call fig_close(u)

! a particle sliding along a wall at a gap of 0.02: plain iterates (top), directional (bottom, shifted down by 0.3)
call wall%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[-0.1_R8P, -0.1_R8P, 0._R8P, 0._R8P], loop_sizes=[4_I4P])
gap = 0.02_R8P
w0 = pose_2d(angle=0._R8P, translation=[0.1_R8P, 0.03_R8P + gap])
w1 = pose_2d(angle=0._R8P, translation=[0.9_R8P, 0.03_R8P + gap])
call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), w0, w1, toi, status, &
                         tolerance=1.e-6_R8P, iterations=steps)
write(caption, '(A,I0,A)') 'sliding at a gap of 0.02: plain bound ', steps, ' steps (top), directional 1 (bottom)'
call fig_open(u, 'tut14-sliding', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
do panel = 1, 2
   lift = -0.3_R8P * (panel - 1)
   call draw_curve(u, 'body', wall, pose_2d(angle=0._R8P), lift)
   call draw_curve(u, 'ghost', disc, w0, lift)
   do k = 1, merge(steps, 1_I4P, panel == 1)
      call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), w0, w1, toi, status, &
                               tolerance=1.e-6_R8P, max_iterations=k, directional=panel == 2)
      call draw_curve(u, merge('box   ', 'accent', k < steps .and. panel == 1), disc, pose_interpolate(w0, w1, toi), lift)
   enddo
enddo
call fig_label(u, 'accent', vector_R8P(0.02_R8P, 0.12_R8P, 0._R8P), 'plain bound: steps of (d - tol/2) / mu')
call fig_label(u, 'accent', vector_R8P(0.02_R8P, -0.18_R8P, 0._R8P), 'directional bound: one step')
call fig_close(u)
contains
   subroutine make_disc(radius, c)
   !< A regular 32-gon of the given radius about the origin, counter-clockwise (solid inside).
   real(R8P),          intent(in)    :: radius
   type(curve_object), intent(inout) :: c
   real(R8P)                         :: ang(32)
   integer(I4P)                      :: j

   ang = [(2._R8P * PI * j / 32._R8P, j = 0, 31)]
   call c%adopt_loops(u=radius * cos(ang), v=radius * sin(ang), loop_sizes=[32_I4P])
   endsubroutine make_disc

   subroutine draw_curve(unit, style, c, pose, dy)
   !< The segments of a curve placed at a pose and shifted by dy along y, in the plane z = 0.
   integer(I4P),          intent(in) :: unit
   character(*),          intent(in) :: style
   type(curve_object),    intent(in) :: c
   type(rigid_pose_2d_t), intent(in) :: pose
   real(R8P),             intent(in) :: dy
   type(vector_R8P), allocatable     :: e(:), f(:)
   real(R8P)                         :: x(2), y(2)
   integer(I4P)                      :: i, n

   allocate(e(c%get_vertices_number()), f(c%get_vertices_number()))
   n = 0
   do i = 1, c%get_vertices_number()
      if (c%get_next(i) == 0) cycle
      x = to_world(pose, c%get_vertex(i)) ; y = to_world(pose, c%get_vertex(c%get_next(i)))
      n = n + 1 ; e(n) = vector_R8P(x(1), x(2) + dy, 0._R8P) ; f(n) = vector_R8P(y(1), y(2) + dy, 0._R8P)
   enddo
   call fig_segments(unit, trim(style), e(1:n), f(1:n))
   endsubroutine draw_curve
endprogram fig_tut14
