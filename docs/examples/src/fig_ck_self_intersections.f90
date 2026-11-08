!run figck-self-intersections fig_ck_self_intersections
!fig ck-self-intersections-before
!fig ck-self-intersections-after
!< Figures of the cookbook recipe "Find self-intersections and remove them": the two crossing cubes merged into one
!< surface (ghost) with the segments along which find_self_intersections reports the facet pairs crossing, then the
!< surface after resolve_self_intersections, its walls inside the other cube gone.
program fig_ck_self_intersections
use fossil,     only : surface_stl_object, intersection_pair_t
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)               :: a, b, soup
type(intersection_pair_t), allocatable :: pairs(:)
character(80)                          :: caption
integer(I4P)                           :: u

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
b = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.35_R8P, 0.4_R8P)) ; call b%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%analyze
call soup%find_self_intersections(pairs=pairs)
write(caption, '(A,I0,A)') 'merged: ', size(pairs), ' crossing facet pairs, along the ember segments'
call fig_open(u, 'ck-self-intersections-before', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'ghost', soup)
call fig_segments(u, 'hole', pairs%p, pairs%q)
call fig_close(u)
call soup%resolve_self_intersections
write(caption, '(A,I0,A,F6.4)') 'resolved: ', soup%get_facets_number(), ' facets, no crossing, volume ', &
                                soup%get_volume()
call fig_open(u, 'ck-self-intersections-after', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'body', soup)
call fig_close(u)
endprogram fig_ck_self_intersections
