!run fig01 fig_tut01
!fig tut01-dragon
!< Figure of tutorial 1: the dragon as loaded and the bounding box analyze computes.
program fig_tut01
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
call fig_open(u, 'tut01-dragon', 'dragon.stl: 6588 facets and the bounding box analyze computes', -60._R8P, 22._R8P)
call fig_surface(u, 'body', surface)
call fig_box(u, 'box', surface%get_bmin(), surface%get_bmax())
call fig_close(u)
endprogram fig_tut01
