!run figck-union fig_ck_union
!fig ck-union-inputs
!fig ck-union-result
!< Figures of the recipe "Unite two bodies": the two unit cubes of the recipe (b translated by 0.5 along the
!< diagonal), then a%boolean(b, BOOL_UNION), seen from the same camera.
program fig_ck_union
use fossil,     only : surface_stl_object, BOOL_UNION, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: a, b
character(80)            :: caption
integer(I4P)             :: u, status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze
call fig_open(u, 'ck-union-inputs', 'the inputs: a (translucent) and b, translated by 0.5 (ember)', -60._R8P, 25._R8P)
call fig_surface(u, 'patch', b)
call fig_surface(u, 'ghost', a)
call fig_close(u)
call a%boolean(other=b, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
write(caption, '(A,I0,A,F6.4)') 'a%boolean(b, BOOL_UNION): ', a%get_facets_number(), ' facets, volume ', a%get_volume()
call fig_open(u, 'ck-union-result', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', a)
call fig_close(u)
endprogram fig_ck_union
