!< Animation of tutorial 1: the dragon turned about the vertical axis through the centre of its box, 30 degrees per
!< frame, and analyzed again each time. The box analyze computes is axis-aligned, so its sides change with the turn;
!< the volume does not. Statistics from the full surface; the drawing is a copy decimated to 900 facets.
program anim_tut01_turn
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, light, turned, drawn
type(vector_R8P)         :: bmin, bmax, centre
character(96)            :: caption
integer(I4P)             :: u, k, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
centre = 0.5_R8P * (bmin + bmax)
light = surface
call light%decimate(target_facets=900_I4P, status=status)                   ! a light drawing
call fig_open(u, 'anim-tut01-turn', '', -60._R8P, 22._R8P)
call fig_anim(u, 0.6_R8P)
do k = 0, 11
   turned = surface ; drawn = light
   call turned%rotate(axis=ez_R8P, angle=k * PI / 6._R8P, center=centre)
   call turned%analyze                                                       ! box and volume of the turned body
   call drawn%rotate(axis=ez_R8P, angle=k * PI / 6._R8P, center=centre)
   bmin = turned%get_bmin() ; bmax = turned%get_bmax()
   write(caption, '(A,I3,A,F5.2,A,F5.2,A,F5.2,A,F7.2)') 'turn ', 30 * k, ': box ', bmax%x - bmin%x, ' x ', &
                                                       bmax%y - bmin%y, ' x ', bmax%z - bmin%z, ', volume ', &
                                                       turned%get_volume()
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call fig_surface(u, 'body', drawn)
   call fig_box(u, 'box', bmin, bmax)
enddo
call fig_close(u)
endprogram anim_tut01_turn
