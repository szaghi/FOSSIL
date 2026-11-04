!< Figure of the cookbook recipe "Check that a surface is valid": the NACA 0012 wing section and its boundary edges
!< (facet edges without a neighbour), which make it manifold but not watertight. They are the seams of the two caps:
!< the caps are triangulated on vertices that only partly match the side's rims, so no edge there is shared.
program fig_ck_validity
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_boundary, fig_close
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: wing
type(facet_object), pointer :: f
integer(I4P)                :: u, i, n
character(80)               :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
n = 0
do i = 1, wing%get_facets_number()
   f => wing%facet_at(i) ; n = n + count(f%fcon_edge == 0_I4P)
enddo
write(caption, '(A,I0,A)') 'naca0012-ascii.stl: ', n, ' boundary edges, where the caps meet the side'
call fig_open(u, 'ck-validity', trim(caption), -60._R8P, 35._R8P)
call fig_surface(u, 'body', wing)
call fig_boundary(u, 'hole', wing)
call fig_close(u)
endprogram fig_ck_validity
