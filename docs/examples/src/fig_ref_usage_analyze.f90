!run figref-usage-analyze fig_ref_usage_analyze
!fig ref-usage-analyze
!< Figure of the usage reference, "Printing statistics": what analyze computes on the dragon as loaded: the bounding
!< box (its extents), the centroid, and the boundary edges (edges with one facet) that make is_watertight false.
program fig_ref_usage_analyze
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_box, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: dragon
type(vector_R8P)         :: bmin, bmax, c
integer(I4P)             :: u
character(120)           :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax() ; c = dragon%get_centroid()
write(caption, '(A,I0,A,L1,A)') 'analyze: bounding box, centroid, boundary edges (', dragon%get_facets_number(), &
                                ' facets, is_watertight ', dragon%is_watertight(), ')'
call fig_open(u, 'ref-usage-analyze', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_boundary(u, 'hole', dragon)
call fig_box(u, 'box', bmin, bmax)
call fig_points(u, 'accent', [c])
call fig_label(u, 'accent', c, 'centroid')
call fig_label(u, 'box', bmin, 'bmin')
call fig_label(u, 'box', bmax, 'bmax')
call fig_close(u)
endprogram fig_ref_usage_analyze
