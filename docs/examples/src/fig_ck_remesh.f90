!run figck-remesh fig_ck_remesh
!fig ck-remesh-before
!fig ck-remesh-after
!< Figures of the recipe "Remesh to a uniform edge length": the unit cube as loaded, then remeshed by the recipe's
!< call (target edge 0.1, 10 iterations, sharp edges kept), seen from the same camera.
program fig_ck_remesh
use fossil,     only : surface_stl_object, REM_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P

implicit none
type(surface_stl_object) :: cube
character(80)            :: caption
integer(I4P)             :: u, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
write(caption, '(A,I0,A,F6.4,A,F6.4)') 'cube-01.stl: ', cube%get_facets_number(), ' facets, edges ', &
                                       cube%smallest_edge_len(), ' to ', cube%largest_edge_len()
call fig_open(u, 'ck-remesh-before', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', cube)
call fig_close(u)
call cube%isotropic_remesh(target_length=0.1_R8P, iterations=10_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
write(caption, '(A,I0,A,F6.4,A,F6.4)') 'isotropic_remesh(0.1): ', cube%get_facets_number(), ' facets, edges ', &
                                       cube%smallest_edge_len(), ' to ', cube%largest_edge_len()
call fig_open(u, 'ck-remesh-after', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', cube)
call fig_close(u)
endprogram fig_ck_remesh
