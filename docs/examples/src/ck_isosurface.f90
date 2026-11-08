!run ck_isosurface ck_isosurface
program ck_isosurface
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter            :: R = 1._R8P, rt = 0.4_R8P, h = 0.1_R8P   ! torus radii, sample spacing
type(surface_stl_object)        :: torus
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: f(31, 31, 11), x(31), z(11)
integer(I4P)                    :: i, j, k, status

x = [(-1.5_R8P + h * i, i = 0, 30)] ; z = [(-0.5_R8P + h * i, i = 0, 10)]
do k = 1, 11 ; do j = 1, 31 ; do i = 1, 31                     ! samples of the signed distance to a torus
   f(i, j, k) = sqrt((sqrt(x(i)**2 + x(j)**2) - R)**2 + z(k)**2) - rt
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), x(1), z(1)), bmax=vector_R8P(x(31), x(31), z(11)), &
                        iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'extract_isosurface failed'
call torus%adopt_facets(facets=facets)
print '(A,I0,A,L1)', 'facets ', torus%get_facets_number(), ', watertight ', torus%is_watertight()
print '(A,F8.5,A,F8.5)', 'volume ', torus%get_volume(), ', exact 2 pi^2 R r^2 ', 2 * acos(-1._R8P)**2 * R * rt**2
endprogram ck_isosurface
