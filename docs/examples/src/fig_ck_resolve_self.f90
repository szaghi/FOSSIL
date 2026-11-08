!run figck-resolve-self fig_ck_resolve_self
!fig ck-resolve-self-arrangement
!fig ck-resolve-self-union
!< Figures of the cookbook recipe "Resolve the self-intersections of one surface": what resolve_self_intersections
!< returns for the two crossing cubes (issue #78: the arrangement, translucent, its open edges in ember) and the union
!< of the two shells the recipe falls back to.
program fig_ck_resolve_self
use fossil,     only : surface_stl_object, facet_object, BOOL_UNION
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object) :: soup, resolved, a, b
integer(I4P)             :: u
character(90)            :: caption

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
b = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)) ; call b%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%analyze
resolved = soup
call resolved%resolve_self_intersections
write(caption, '(A,F6.4,A)') 'resolve_self_intersections: volume ', resolved%get_volume(), &
                             ', open edges in ember (issue #78)'
call fig_open(u, 'ck-resolve-self-arrangement', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'ghost', resolved)
call fig_boundary(u, 'hole', resolved)
call fig_close(u)
call a%boolean(other=b, op=BOOL_UNION)
write(caption, '(A,I0,A,F6.4,A)') 'union of the two shells: ', a%get_facets_number(), ' facets, volume ', &
                                  a%get_volume(), ', watertight'
call fig_open(u, 'ck-resolve-self-union', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'body', a)
call fig_close(u)
endprogram fig_ck_resolve_self
