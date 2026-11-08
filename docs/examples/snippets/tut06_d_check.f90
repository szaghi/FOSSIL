!< Tutorial 6, fourth run: the field against the exact signed distance at every node, in and out of the band. The
!< figure shows the error on the grid plane x = 0 (colour sqrt(error / largest error), cream 0 to sepia).
program tut06_d_check
use fossil,     only : surface_stl_object, distance_field_object, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)    :: surface
type(distance_field_object) :: field
real(R8P), allocatable      :: x(:), err(:,:,:), lo(:), hi(:)
real(R8P)                   :: h, d, err_band, err_far
integer(I4P)                :: i, j, k, m, n, u
character(96)               :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! the grid and field of the first runs
h = 0.04_R8P
do m = 1, 6
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
n = size(x)
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
allocate(err(n, n, n))
err_band = 0._R8P ; err_far = 0._R8P
do k = 1, n ; do j = 1, n ; do i = 1, n
   d = surface%distance(point=x(i) * ex_R8P + x(j) * ey_R8P + x(k) * ez_R8P, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   err(i, j, k) = abs(field%phi(i, j, k) - d)
   if (field%in_band(i, j, k)) then
      err_band = max(err_band, err(i, j, k))
   else
      err_far = max(err_far, err(i, j, k))
   endif
enddo ; enddo ; enddo
print '(A,ES10.3)', 'max error in the band ', err_band
print '(A,ES10.3)', 'max error elsewhere   ', err_far
i = minloc(abs(x), dim=1)                                       ! the plane x = 0, nodes as dual cells
lo = [x(1), (0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1)]
hi = [(0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1), x(n)]
write(caption, '(A,F6.4,A,F6.4,A)') 'error on the plane x = 0: 0 in the band, up to ', maxval(err(i, :, :)), &
                                    ' (largest anywhere ', maxval(err), ')'
print '(A)', trim(caption)
call fig_open(u, 'tut06-error', trim(caption), 0._R8P, 0._R8P, height=500_I4P)
do k = 1, n ; do j = 1, n
   call fig_quad(u, 'frac', p(lo(j), lo(k)), p(hi(j), lo(k)), p(hi(j), hi(k)), p(lo(j), hi(k)), &
                 sqrt(err(i, j, k) / maxval(err)))
enddo ; enddo
call fig_close(u)
contains
   function p(y, z) result(q)
   !< A point of the plane x = x(i).
   real(R8P), intent(in) :: y, z !< In-plane coordinates.
   type(vector_R8P)      :: q    !< Point.

   q%x = x(i) ; q%y = y ; q%z = z
   endfunction p
endprogram tut06_d_check
