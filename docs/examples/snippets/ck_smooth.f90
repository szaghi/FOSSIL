program ck_smooth
use fossil, only : surface_stl_object, facet_object, extract_isosurface, SMOOTH_METHOD_TAUBIN, SMOOTH_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: f(29, 29, 29), x(29)
integer(I4P)                    :: i, j, k, status

x = [(-1.4_R8P + 0.1_R8P * i, i = 0, 28)]
do k = 1, 29 ; do j = 1, 29 ; do i = 1, 29
   f(i, j, k) = norm2([x(i), x(j), x(k)]) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), x(1), x(1)), bmax=vector_R8P(x(29), x(29), x(29)), &
                        iso=0._R8P, surface=facets)               ! a unit sphere ...
do i = 1, size(facets) ; do j = 1, 3                            ! ... with noise: each vertex moved radially by up to 3 %
   facets(i)%vertex(j) = (1._R8P + 0.03_R8P * sin(1.e3_R8P * (facets(i)%vertex(j)%x + 2 * facets(i)%vertex(j)%y + &
                         3 * facets(i)%vertex(j)%z))) * facets(i)%vertex(j)
enddo ; enddo
call sphere%adopt_facets(facets=facets)
call report('noisy:   ')
call sphere%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=10_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'smoothing failed'
call report('smoothed:')
print '(A,F8.5,A,F8.5)', 'unit sphere:                   area ', 4 * acos(-1._R8P), ', volume ', 4 * acos(-1._R8P) / 3
contains
   subroutine report(label)
   character(*), intent(in)    :: label
   type(facet_object), pointer :: g
   real(R8P)                   :: e
   integer(I4P)                :: n, v

   e = 0._R8P
   do n = 1, sphere%get_facets_number()
      g => sphere%facet_at(n)
      do v = 1, 3
         e = e + (norm2([g%vertex(v)%x, g%vertex(v)%y, g%vertex(v)%z]) - 1._R8P)**2
      enddo
   enddo
   print '(2A,F8.5,A,F8.5,A,F8.5)', label, ' rms |r - 1| ', sqrt(e / (3 * sphere%get_facets_number())), &
                                    ', area ', sphere%get_area(), ', volume ', sphere%get_volume()
   endsubroutine report
endprogram ck_smooth
