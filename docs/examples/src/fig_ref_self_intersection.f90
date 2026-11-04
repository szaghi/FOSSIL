!run figref-self-intersection fig_ref_self_intersection
!fig ref-self-intersection-dragon
!< Figure of the reference page "Self-intersection": the dragon (ghost) and, in ember, the segment p-q along which
!< each facet pair find_self_intersections reports crosses.
program fig_ref_self_intersection
use fossil,     only : surface_stl_object, intersection_pair_t
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object)               :: dragon
type(intersection_pair_t), allocatable :: pairs(:)
character(80)                          :: caption
integer(I4P)                           :: u, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%find_self_intersections(pairs=pairs, status=status)
write(caption, '(A,I0,A,I0,A)') 'find_self_intersections: ', size(pairs), ' crossing pairs in ', &
                                dragon%get_facets_number(), ' facets'
call fig_open(u, 'ref-self-intersection-dragon', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_segments(u, 'hole', pairs%p, pairs%q)
call fig_close(u)
endprogram fig_ref_self_intersection
