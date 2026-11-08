!< Figure of the cookbook recipe "Extract an isosurface from a sampled field": the torus extract_isosurface draws
!< from the recipe's samples (spacing 0.1) of its signed distance, with the sampling box.
program fig_ck_isosurface
use fossil,     only : surface_stl_object, facet_object, extract_isosurface
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: R = 1._R8P, rt = 0.4_R8P, h = 0.1_R8P
type(surface_stl_object)        :: torus
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: f(31, 31, 11), x(31), z(11)
integer(I4P)                    :: i, j, k, u
character(80)                   :: caption

x = [(-1.5_R8P + h * i, i = 0, 30)] ; z = [(-0.5_R8P + h * i, i = 0, 10)]
do k = 1, 11 ; do j = 1, 31 ; do i = 1, 31
   f(i, j, k) = sqrt((sqrt(x(i)**2 + x(j)**2) - R)**2 + z(k)**2) - rt
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), x(1), z(1)), bmax=vector_R8P(x(31), x(31), z(11)), &
                        iso=0._R8P, surface=facets)
call torus%adopt_facets(facets=facets)
write(caption, '(A,I0,A)') 'extract_isosurface: the zero level of 31 x 31 x 11 samples, ', &
                           torus%get_facets_number(), ' facets'
call fig_open(u, 'ck-isosurface', trim(caption), -60._R8P, 35._R8P)
call fig_box(u, 'box', vector_R8P(x(1), x(1), z(1)), vector_R8P(x(31), x(31), z(11)))
call fig_surface(u, 'body', torus)
call fig_close(u)
endprogram fig_ck_isosurface
