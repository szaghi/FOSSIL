!run figck-clip fig_ck_clip
!fig ck-clip
!< Figure of the cookbook recipe "Split a surface by a box": the dragon clipped by the box x <= 0 (the head and
!< front, body colour), the remainder (ghost) and the box.
program fig_ck_clip
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object) :: dragon, tail
type(vector_R8P)         :: bmin, bmax
integer(I4P)             :: u
character(90)            :: caption

bmin = vector_R8P(-12._R8P, -6._R8P, 4._R8P) ; bmax = vector_R8P(0._R8P, 6._R8P, 21._R8P)
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%clip(bmin=bmin, bmax=bmax, remainder=tail)
write(caption, '(A,I0,A,I0,A)') 'clip: ', dragon%get_facets_number(), ' facets kept (stone), ', &
                                tail%get_facets_number(), ' in the remainder (ghost)'
call fig_open(u, 'ck-clip', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', tail)
call fig_surface(u, 'body', dragon)
call fig_box(u, 'accent', bmin, bmax)
call fig_close(u)
endprogram fig_ck_clip
