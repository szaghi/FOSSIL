!< Figures of the reference page "Mesh decimation (QEM)": the page decimates the Stanford bunny (69 451 facets, too
!< heavy for a figure); here the dragon as loaded, then decimate(target_facets=n/4), the page's call, seen from the
!< same camera.
program fig_ref_decimate
use fossil,     only : surface_stl_object, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object) :: dragon
character(80)            :: caption
integer(I4P)             :: u, n0, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
n0 = dragon%get_facets_number()
write(caption, '(A,I0,A)') 'dragon.stl as loaded: ', n0, ' facets'
call fig_open(u, 'ref-decimate-dragon-input', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_close(u)
call dragon%decimate(target_facets=n0 / 4_I4P, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
write(caption, '(A,I0,A)') 'decimate(target_facets=n0/4): ', dragon%get_facets_number(), ' facets (DEC_STATUS_OK)'
call fig_open(u, 'ref-decimate-dragon-output', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_close(u)
endprogram fig_ref_decimate
