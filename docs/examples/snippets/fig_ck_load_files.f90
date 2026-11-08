!< Figure of the cookbook recipe "Load several files as one surface": cube-01.stl (binary, patch 1, translucent) and
!< naca0012-ascii.stl (ASCII, patch 2, ember) read by one load_from_files call.
program fig_ck_load_files
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
implicit none
type(surface_stl_object) :: scene
integer(I4P)             :: u, n1, i
character(90)            :: caption

call scene%load_from_files(file_names=[character(18) :: 'cube-01.stl', 'naca0012-ascii.stl'])
n1 = count([(scene%facet_patch(i) == 1, i = 1, scene%get_facets_number())])
write(caption, '(A,I0,A,I0,A)') 'one surface, two patches: cube-01.stl (', n1, ' facets) and naca0012-ascii.stl (', &
                                scene%get_facets_number() - n1, ')'
call fig_open(u, 'ck-load-files', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'ghost', scene, 1_I4P, n1)
call fig_surface(u, 'patch', scene, n1 + 1_I4P, scene%get_facets_number())
call fig_close(u)
endprogram fig_ck_load_files
