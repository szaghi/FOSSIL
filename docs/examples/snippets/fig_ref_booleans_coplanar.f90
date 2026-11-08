!< Figure of the reference page "Boolean operations", section "Coplanar faces": two unit cubes offset by 0.3 along x,
!< so that four faces of each lie in the planes of four faces of the other and overlap on half their area, united by
!< a%boolean(b, BOOL_UNION). The result is drawn with its facet edges: on the four merged faces the triangulation is cut
!< along the edges of both inputs (the lines x = 0.3, 0.5, 0.8 and 1, and the diagonals of both cubes), with no crack. The
!< outlines of the two inputs are drawn as boxes.
program fig_ref_booleans_coplanar
use fossil,     only : surface_stl_object, BOOL_UNION, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: a, b
character(96)            :: caption
integer(I4P)             :: u, status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.3_R8P, 0._R8P, 0._R8P)) ; call b%analyze
call a%boolean(other=b, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
write(caption, '(A,I0,A,F6.4,A,L1,A,L1)') 'coplanar union: ', a%get_facets_number(), ' facets, volume ', &
      a%get_volume(), ', watertight ', a%is_watertight(), ', manifold ', a%is_manifold()
call fig_open(u, 'ref-booleans-coplanar', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', a)
call fig_box(u, 'accent', vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(1._R8P, 1._R8P, 1._R8P))
call fig_box(u, 'probe', vector_R8P(0.3_R8P, 0._R8P, 0._R8P), vector_R8P(1.3_R8P, 1._R8P, 1._R8P))
call fig_close(u)
endprogram fig_ref_booleans_coplanar
