!run anim14 anim_tut14_advance
!fig anim-tut14-advance
!< Animation of tutorial 14: conservative advancement of the pebble thrown at the dragon, one frame per step. Frame k
!< calls time_of_impact with max_iterations = k, whose toi is the k-th iterate s_k: the pebble at s_k (ember) over the
!< earlier iterates (grey), the closest pair at s_k (yellow), and the next step (d - tol/2) / mu drawn as the stretch
!< of path the pebble's centre may cover before the next distance query (ember line). The dragon is placed as in the
!< chapter (self-contained).
program anim_tut14_advance
use fossil,     only : surface_stl_object, facet_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, &
                       CCD_CONTACT
use fossil,     only : motion_bound_3d
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_triangle, fig_segment, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P

implicit none
real(R8P),    parameter  :: PI = acos(-1._R8P)
type(vector_R8P), parameter :: CENTRE = vector_R8P(-0.37_R8P, -0.3_R8P, -0.17_R8P) ! the view
integer(I4P), parameter  :: KMAX = 5
type(surface_stl_object) :: dragon, pebble, moved
type(rigid_pose_t)       :: p0, p1
type(vector_R8P)         :: bmin, bmax, da, db, c
type(facet_object), pointer :: f
real(R8P), allocatable   :: x(:), y(:), z(:)
real(R8P)                :: mu, scale, s(0:KMAX), d(0:KMAX), pa(3, 0:KMAX), pb(3, 0:KMAX), tol, dt
character(100)           :: caption
integer(I4P)             :: u, k, j, i, status, n

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
da = dragon%get_bmax() - dragon%get_bmin() ; db = pebble%get_bmax() - pebble%get_bmin()
tol = 1.e-6_R8P * max(da%normL2(), db%normL2())                                   ! the default tolerance
allocate(x(3 * pebble%get_facets_number()), y(3 * pebble%get_facets_number()), z(3 * pebble%get_facets_number()))
do i = 1, pebble%get_facets_number()                                              ! the pebble's speed bound
   f => pebble%facet_at(i)
   x(3 * i - 2:3 * i) = f%vertex%x ; y(3 * i - 2:3 * i) = f%vertex%y ; z(3 * i - 2:3 * i) = f%vertex%z
enddo
mu = motion_bound_3d(p0, p1, x, y, z)

s(0) = 0._R8P                                                                     ! the iterates s_k, d(s_k)
call dragon%distance_to(pebble, distance=d(0), point_self=pa(:, 0), point_other=pb(:, 0), pose_other=p0)
n = KMAX
do k = 1, KMAX
   call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, s(k), status, distance=d(k), &
                              point_self=pa(:, k), point_other=pb(:, k), max_iterations=k)
   if (status == CCD_CONTACT) then
      n = k ; exit
   endif
enddo
do k = 0, n
   print '(A,I0,A,F9.6,A,ES10.3)', 'step ', k, ': s = ', s(k), ', gap ', d(k)
enddo

call fig_open(u, 'anim-tut14-advance', 'conservative advancement', -65._R8P, 20._R8P)
call fig_zoom(u, CENTRE, 0.27_R8P)
do i = 1, dragon%get_facets_number()                       ! the part of the dragon in view, its side facing the camera
   f => dragon%facet_at(i)
   c = (f%vertex(1) + f%vertex(2) + f%vertex(3)) / 3._R8P
   if (c%x > CENTRE%x + 0.4_R8P .or. c%y > 0._R8P .or. abs(c%z - CENTRE%z) > 0.3_R8P) cycle
   call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_segment(u, 'box', vector_R8P(-0.6_R8P, -0.3_R8P, -0.2_R8P), vector_R8P(0.6_R8P, -0.3_R8P, -0.2_R8P))
call fig_anim(u, 1.6_R8P)
do k = 0, n
   if (k < n) then
      dt = (d(k) - 0.5_R8P * tol) / mu
      write(caption, '(A,I0,A,F8.6,A,ES9.2,A,ES9.2)') 'step ', k, ': s = ', s(k), ', gap ', d(k), &
                                                     ', next step (d - tol/2)/mu ', dt
      call fig_frame(u, trim(caption))
   else
      write(caption, '(A,I0,A,F8.6,A,ES9.2,A)') 'step ', k, ': s = ', s(k), ', gap ', d(k), ' <= tol: contact'
      call fig_frame(u, trim(caption), seconds=3._R8P)
   endif
   do j = 0, k - 1
      moved = pebble ; call moved%move(pose_interpolate(p0, p1, s(j))) ; call fig_surface(u, 'ghost', moved)
   enddo
   moved = pebble ; call moved%move(pose_interpolate(p0, p1, s(k))) ; call fig_surface(u, 'patch', moved)
   call fig_segment(u, 'accent', vector_R8P(pa(1, k), pa(2, k), pa(3, k)), vector_R8P(pb(1, k), pb(2, k), pb(3, k)))
   if (k < n) call fig_segment(u, 'hole', vector_R8P(-0.6_R8P + 1.2_R8P * s(k), -0.3_R8P, -0.165_R8P), &
                               vector_R8P(-0.6_R8P + 1.2_R8P * (s(k) + dt), -0.3_R8P, -0.165_R8P))
enddo
call fig_close(u)
endprogram anim_tut14_advance
