!< Animation of tutorial 3: the placement of the repaired dragon step by step, in CAD units, about the fixed world axes
!< (x and y through the origin, gold): as repaired, its box centre translated to the origin, turned a quarter about z
!< in three steps, resized to a longest side of 1. The dragon is repaired as in tutorial 2 (sanitize with fill_holes);
!< the box and the captions are those of the full surface, the drawing a copy decimated to 900 facets.
program anim_tut03_place
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_box, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = 4._R8P * atan(1._R8P), L = 14._R8P   ! half-length of the drawn axes
type(surface_stl_object) :: surface, light
type(vector_R8P)         :: bmin, bmax, delta
real(R8P)                :: scale
character(96)            :: caption
integer(I4P)             :: u, k, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)  ! tutorial 2
call surface%sanitize(fill_holes=.true.)
light = surface
call light%decimate(target_facets=900_I4P, status=status)
call fig_open(u, 'anim-tut03-place', '', -60._R8P, 30._R8P)
call fig_anim(u, 1.2_R8P)
call fig_segment(u, 'accent', -L * ex_R8P, L * ex_R8P) ; call fig_label(u, 'accent', L * ex_R8P, 'x')
call fig_segment(u, 'accent', -L * ey_R8P, L * ey_R8P) ; call fig_label(u, 'accent', L * ey_R8P, 'y')
call frame('as repaired, in CAD units')
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
delta = -0.5_R8P * (bmin + bmax)
call surface%translate(delta=delta) ; call light%translate(delta=delta)
call surface%analyze
write(caption, '(A,3F7.2,A)') 'translate(delta=', delta%x, delta%y, delta%z, ')'
call frame(trim(caption))
do k = 1, 3
   call surface%rotate(axis=ez_R8P, angle=PI / 6._R8P) ; call light%rotate(axis=ez_R8P, angle=PI / 6._R8P)
   call surface%analyze
   write(caption, '(A,I0,A)') 'rotate(axis=ez_R8P, angle=', 30 * k, ' degrees in all)'
   call frame(trim(caption))
enddo
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale) ; call light%resize(x=scale, y=scale, z=scale)
call surface%analyze
write(caption, '(A,F8.6,A)') 'resize(', scale, '): longest side 1, the small body at the origin'
call frame(trim(caption))
call fig_close(u)
contains
   subroutine frame(caption)
   !< A frame: the body (decimated copy) and the box of the full surface.
   character(*), intent(in) :: caption !< Caption.

   print '(A)', caption
   call fig_frame(u, caption=caption)
   call fig_surface(u, 'body', light)
   call fig_box(u, 'box', surface%get_bmin(), surface%get_bmax())
   endsubroutine frame
endprogram anim_tut03_place
