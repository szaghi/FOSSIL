!< Animation of the recipe "Decimate a mesh": the dragon as loaded, then decimated again and again to half of its
!< facets by the recipe's call (quadric edge collapses), down to about 200 facets, seen from the same camera.
program anim_decimate
use fossil,     only : surface_stl_object, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object) :: dragon
character(80)            :: caption
integer(I4P)             :: u, k, target, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call fig_open(u, 'anim-decimate', '', -60._R8P, 22._R8P)
call fig_anim(u, 1.0_R8P)
write(caption, '(A,I0,A)') 'dragon.stl as loaded: ', dragon%get_facets_number(), ' facets'
call fig_frame(u, caption=trim(caption))
call fig_surface(u, 'body', dragon)
do k = 1, 5
   target = dragon%get_facets_number() / 2
   call dragon%decimate(target_facets=target, status=status)
   if (status /= DEC_STATUS_OK) error stop 'decimate failed'
   write(caption, '(A,I0,A,I0,A)') 'decimate(target_facets=', target, '): ', dragon%get_facets_number(), ' facets'
   print '(A)', trim(caption)
   if (k == 5) then
      call fig_frame(u, caption=trim(caption), seconds=2.5_R8P)
   else
      call fig_frame(u, caption=trim(caption))
   endif
   call fig_surface(u, 'body', dragon)
enddo
call fig_close(u)
endprogram anim_decimate
