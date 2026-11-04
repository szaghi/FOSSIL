!run figck-decimate fig_ck_decimate
!fig ck-decimate-before
!fig ck-decimate-after
!< Figures of the recipe "Decimate a mesh": the dragon as loaded, then decimated to a quarter of its facets by the
!< recipe's call (quadric edge collapses), seen from the same camera.
program fig_ck_decimate
use fossil,     only : surface_stl_object, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object) :: dragon
character(80)            :: caption
integer(I4P)             :: u, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
write(caption, '(A,I0,A)') 'dragon.stl as loaded: ', dragon%get_facets_number(), ' facets'
call fig_open(u, 'ck-decimate-before', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_close(u)
call dragon%decimate(target_facets=dragon%get_facets_number() / 4, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
write(caption, '(A,I0,A)') 'decimate(target_facets=n/4): ', dragon%get_facets_number(), ' facets'
call fig_open(u, 'ck-decimate-after', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_close(u)
endprogram fig_ck_decimate
