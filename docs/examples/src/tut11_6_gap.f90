!run tut11-gap tut11_6_gap
!fig tut11-gap
!< Tutorial 11, step 6: the distance between two bodies. The dragon at the pose of step 1 and a unit cube beside it:
!< distance_to with pose_self= on the dragon as built, against the moved copy. The figure: the two bodies, the closest
!< pair joined by a segment and its length.
program tut11_6_gap
!region gap
use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter  :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, moved, box
type(rigid_pose_t)       :: p1
real(R8P)                :: gap, gap1, ps(3), po(3)
integer(I4P)             :: fa, fb
logical                  :: hit

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! a unit cube beside the dragon
call box%translate(x=0.7_R8P, y=-0.5_R8P, z=-0.5_R8P)
call box%analyze
call surface%distance_to(box, distance=gap, point_self=ps, point_other=po, facet_self=fa, facet_other=fb, &
                         intersecting=hit, pose_self=p1)
print '(A,F8.5,A,I0,A,I0,A,L1)', 'gap to the cube ', gap, ' (facets ', fa, ', ', fb, '), touching ', hit
print '(A,3F8.4,A,3F8.4)', 'closest pair', ps, '  -', po
!endregion gap
!region moved
moved = surface
call moved%move(p1)
call moved%distance_to(box, distance=gap1)
print '(A,L1)', 'pose_self = moved copy, bit for bit: ', gap == gap1
!endregion moved
call figure
contains
   subroutine figure
   !< The dragon at pose1, the cube and the closest pair.
   type(vector_R8P) :: a, b
   character(24)    :: text
   integer(I4P)     :: u

   a%x = ps(1) ; a%y = ps(2) ; a%z = ps(3) ; b%x = po(1) ; b%y = po(2) ; b%z = po(3)
   call fig_open(u, 'tut11-gap', 'distance_to: the dragon at pose1, the cube and their closest pair', &
                 -70._R8P, 20._R8P)
   call fig_surface(u, 'body', moved)
   call fig_surface(u, 'ghost', box)
   call fig_segment(u, 'accent', a, b)
   call fig_points(u, 'accent', [a, b])
   write(text, '(A,F7.4)') 'gap ', gap
   call fig_label(u, 'accent', 0.5_R8P * (a + b), trim(text))
   call fig_close(u)
   endsubroutine figure
endprogram tut11_6_gap
