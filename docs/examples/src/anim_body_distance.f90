!run anim-body-distance anim_body_distance
!fig anim-body-distance
!< Animation of the body-distance reference page and of tutorial 11: a cube approaching the dragon on a screw motion.
!< The placed dragon of chapters 2-3 (rebuilt here with the same calls, decimated to 700 facets for a light drawing)
!< stays at rest; a cube of side 0.25 travels from (1.1, 0.35, 0.3) to the dragon's centre while turning half a turn
!< about (1, 1, 1), its pose interpolated by pose_interpolate (a screw motion). Each frame calls distance_to with
!< pose_other= (the cube's tree as built; move only to draw it) and draws the closest pair joined by a segment, the gap
!< in the caption, until the first frame where the bodies touch or cross.
program anim_body_distance
use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_interpolate, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: NS = 24                                     ! steps
real(R8P),    parameter  :: PI = 4._R8P * atan(1._R8P), SIDE = 0.25_R8P
type(surface_stl_object) :: dragon, cube, moved
type(rigid_pose_t)       :: start, finish, pose
type(vector_R8P)         :: bmin, bmax, qa(1), qb(1)
real(R8P)                :: scale, d, pa(3), pb(3), axis(3)
logical                  :: hit
character(80)            :: caption
integer(I4P)             :: u, s, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)             ! chapter 2: repair
call dragon%sanitize(fill_holes=.true.)
call dragon%save_into_file(file_name='anim-body-distance-dragon-clean.stl')
call dragon%load_from_file(file_name='anim-body-distance-dragon-clean.stl', guess_format=.true.) ! chapter 3: place
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze
call dragon%decimate(target_facets=700_I4P, status=status)                           ! a light drawing
if (status /= DEC_STATUS_OK) print '(A,I0)', 'decimate status ', status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)              ! side SIDE, centred at 0
call cube%translate(delta=vector_R8P(-0.5_R8P, -0.5_R8P, -0.5_R8P))
call cube%resize(x=SIDE, y=SIDE, z=SIDE)
call cube%analyze
axis = [1._R8P, 1._R8P, 1._R8P] / sqrt(3._R8P)
start = pose_from_axis_angle(axis=axis, angle=0._R8P, translation=[1.1_R8P, 0.35_R8P, 0.3_R8P])
finish = pose_from_axis_angle(axis=axis, angle=PI, translation=[0._R8P, 0._R8P, 0._R8P])

call fig_open(u, 'anim-body-distance', '', -35._R8P, 20._R8P)
call fig_anim(u, 0.35_R8P)
do s = 0, NS
   pose = pose_interpolate(start, finish, real(s, R8P) / NS)
   call dragon%distance_to(cube, distance=d, point_self=pa, point_other=pb, intersecting=hit, pose_other=pose)
   if (hit) then
      write(caption, '(A,I2,A)') 'step ', s, ': contact, distance_to = 0'
   else
      write(caption, '(A,I2,A,F7.4)') 'step ', s, ': gap ', d
   endif
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption), seconds=merge(2.0_R8P, merge(1.2_R8P, 0.35_R8P, s == 0), hit))
   moved = cube
   call moved%move(pose)                                                             ! only to draw it
   call fig_surface(u, 'body', dragon)
   call fig_surface(u, 'patch', moved)
   qa(1) = vector_R8P(pa(1), pa(2), pa(3)) ; qb(1) = vector_R8P(pb(1), pb(2), pb(3))
   call fig_segment(u, 'accent', qa(1), qb(1))
   call fig_points(u, 'foot', qa)
   call fig_points(u, 'probe', qb)
   if (hit) exit
enddo
call fig_close(u)
endprogram anim_body_distance
