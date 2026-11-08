!run tut05 tut05_a_grid
!fig anim-tut05-batch
!< Tutorial 5, first run: the signed distance of a 16^3 grid of points over the domain in one distance_many call.
!< The animation walks the batch plane by plane in z: each grid point is drawn as a small square coloured by its
!< distance (teal inside, cream on the wall, ember outside; colour tanh(d / 0.2)) over the dragon, drawn translucent
!< as a copy decimated to 300 facets.
program tut05_a_grid
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_quad, fig_box, fig_close
!region grid
use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16                      ! points per direction
type(surface_stl_object)      :: surface, light
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p, status, u
character(80)                 :: caption

allocate(points(n**3), d(n**3))
h = 1.2_R8P / (n - 1)                                        ! the cube [-0.6, 0.6]^3
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%distance_many(points=points, distances=d, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'distance_many failed'
print '(A,I0)',     'points  ', size(points)
print '(A,I0)',     'inside  ', count(d < 0._R8P)
print '(A,2F10.6)', 'min, max d ', minval(d), maxval(d)
!endregion grid
light = surface
call light%decimate(target_facets=300_I4P, status=status)
call fig_open(u, 'anim-tut05-batch', '', -60._R8P, 40._R8P)
call fig_anim(u, 0.45_R8P)
call fig_surface(u, 'ghost', light)
call fig_box(u, 'box', -0.6_R8P * (ex_R8P + ey_R8P + ez_R8P), 0.6_R8P * (ex_R8P + ey_R8P + ez_R8P))
do k = 1, n
   write(caption, '(A,F6.2,A,I0,A,I0,A)') 'z = ', points(1 + (k - 1) * n**2)%z, ': points ', (k - 1) * n**2 + 1, &
                                          ' to ', k * n**2, ' of the batch'
   call fig_frame(u, caption=trim(caption))
   do p = 1 + (k - 1) * n**2, k * n**2
      call fig_quad(u, 'field', points(p) - 0.4_R8P * h * (ex_R8P + ey_R8P), points(p) + 0.4_R8P * h * (ex_R8P - ey_R8P), &
                    points(p) + 0.4_R8P * h * (ex_R8P + ey_R8P), points(p) - 0.4_R8P * h * (ex_R8P - ey_R8P), &
                    tanh(d(p) / 0.2_R8P))
   enddo
enddo
call fig_close(u)
endprogram tut05_a_grid
