!< Figures of the cookbook recipe "Read an STL file, ASCII or binary": the two files the recipe loads with
!< guess_format, the binary dragon and the ASCII NACA 0012 wing section, each with its facet count.
program fig_ck_read_any
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
implicit none
type(surface_stl_object) :: binary, ascii
integer(I4P)             :: u
character(80)            :: caption

call binary%load_from_file(file_name='dragon.stl', guess_format=.true.)
call ascii%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
write(caption, '(A,I0,A)') 'dragon.stl, binary: ', binary%get_facets_number(), ' facets'
call fig_open(u, 'ck-read-any-binary', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', binary)
call fig_close(u)
write(caption, '(A,I0,A)') 'naca0012-ascii.stl, ASCII: ', ascii%get_facets_number(), ' facets'
call fig_open(u, 'ck-read-any-ascii', trim(caption), -60._R8P, 35._R8P)
call fig_surface(u, 'body', ascii)
call fig_close(u)
endprogram fig_ck_read_any
