!run fig06 fig_tut06
!fig tut06-field
!fig tut06-band
!< Figures of tutorial 6: the signed-distance field of the placed dragon (fig_tut03 writes it) on the grid plane through
!< x = 0, seen along x. Each node is drawn as its dual cell (the midpoints to its neighbours), coloured by phi; the
!< section of the dragon by the same plane (slice_surface) is drawn over it. tut06-field: the whole stretched grid,
!< colour sign(phi) sqrt(|phi| / 0.6). tut06-band: the band nodes only, colour phi / band, zoomed on the body
!< (the dark nodes are outside the band).
program fig_tut06
use fossil,     only : surface_stl_object, distance_field_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER, &
                       STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)    :: surface
type(distance_field_object) :: field
type(curve_object)          :: section
real(R8P), allocatable      :: x(:), lo(:), hi(:)
real(R8P)                   :: h, band, xs
integer(I4P)                :: i, j, k, m, n, u, status
character(80)               :: caption

call surface%load_from_file(file_name='fig-dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! the grid of chapter 6
h = 0.04_R8P
do m = 1, 6
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x, status=status)
call field%compute(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
n = size(x)
i = minloc(abs(x), dim=1)                                       ! the grid plane through x = 0
xs = x(i)
band = field%get_band_width()
lo = [x(1), (0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1)]        ! dual cells
hi = [(0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1), x(n)]
call slice_surface(surface=surface, normal_axis=1_I4P, level=xs, curve=section)

write(caption, '(A,F5.2,A,F4.2,A)') 'phi on the plane x = 0, from ', minval(field%phi(i, :, :)), ' to ', &
                                    maxval(field%phi(i, :, :)), '; the section of the dragon'
call fig_open(u, 'tut06-field', trim(caption), 0._R8P, 0._R8P, height=500_I4P)
do k = 1, n ; do j = 1, n
   call fig_quad(u, 'field', p(lo(j), lo(k)), p(hi(j), lo(k)), p(hi(j), hi(k)), p(lo(j), hi(k)), &
                 sign(sqrt(abs(field%phi(i, j, k)) / 0.6_R8P), field%phi(i, j, k)))
enddo ; enddo
call draw_section
call fig_close(u)

write(caption, '(A,I0,A,F4.2,A)') 'the band: ', count(field%in_band(i, :, :)), ' nodes on this plane, |phi| <= ', &
                                  band, ', exact'
call fig_open(u, 'tut06-band', trim(caption), 0._R8P, 0._R8P, height=500_I4P)
call fig_zoom(u, p(0._R8P, 0._R8P), 0.45_R8P)
do k = 1, n ; do j = 1, n
   if (.not. field%in_band(i, j, k)) cycle
   call fig_quad(u, 'field', p(lo(j), lo(k)), p(hi(j), lo(k)), p(hi(j), hi(k)), p(lo(j), hi(k)), &
                 field%phi(i, j, k) / band)
enddo ; enddo
call draw_section
call fig_close(u)
contains
   function p(y, z) result(q)
   !< A point of the plane x = xs.
   real(R8P), intent(in) :: y, z !< In-plane coordinates.
   type(vector_R8P)      :: q    !< Point.

   q%x = xs ; q%y = y ; q%z = z
   endfunction p

   subroutine draw_section
   !< The section of the dragon, over the cells (moved 1e-3 towards the viewer, along x).
   type(vector_R8P) :: a, b
   real(R8P)        :: uv(2)
   integer(I4P)     :: s

   do s = 1, section%get_vertices_number()
      if (section%get_next(s) == 0_I4P) cycle
      uv = section%get_vertex(s) ; a = p(uv(1), uv(2)) ; a%x = xs + 1.e-3_R8P
      uv = section%get_vertex(section%get_next(s)) ; b = p(uv(1), uv(2)) ; b%x = xs + 1.e-3_R8P
      call fig_segment(u, 'accent', a, b)
   enddo
   endsubroutine draw_section
endprogram fig_tut06
