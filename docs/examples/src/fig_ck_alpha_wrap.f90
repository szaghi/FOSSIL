!run figck-alpha-wrap fig_ck_alpha_wrap
!fig ck-alpha-wrap-input
!fig ck-alpha-wrap-wrap
!< Figures of the cookbook recipe "Wrap a triangle soup into a closed surface": the dragon as loaded, its boundary
!< edges (the holes) in ember, and its alpha wrap (alpha 1/20, offset 1/200 of the box diagonal) from the same camera.
program fig_ck_alpha_wrap
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object) :: dragon, wrap
type(vector_R8P)         :: lo, hi
real(R8P)                :: diagonal
integer(I4P)             :: u
character(80)            :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
lo = dragon%get_bmin() ; hi = dragon%get_bmax()
diagonal = norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
call dragon%alpha_wrap(alpha=diagonal / 20._R8P, offset=diagonal / 200._R8P, wrapped=wrap)
write(caption, '(A,I0,A)') 'input: ', dragon%get_facets_number(), ' facets, boundary edges in ember'
call fig_open(u, 'ck-alpha-wrap-input', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_boundary(u, 'hole', dragon)
call fig_close(u)
write(caption, '(A,I0,A)') 'alpha_wrap: ', wrap%get_facets_number(), ' facets, watertight and manifold'
call fig_open(u, 'ck-alpha-wrap-wrap', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', wrap)
call fig_close(u)
endprogram fig_ck_alpha_wrap
