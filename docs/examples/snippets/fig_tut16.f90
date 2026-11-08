!< Figure of tutorial 16: the plane z = const nearest z = 0 of the chapter's 24^3 batch, each point's cell shaded by the
!< signed distance distance_many_device returns (teal inside, ember outside, the colour saturating at half the largest
!< distance), the dragon cut by that plane on top, and the points where it differs from distance_many outlined
!< (the caption counts them). The sign is the winding number, the second of the chapter: on this open dragon the
!< pseudo-normal gives a few far points the wrong sign. The docs are built with gfortran: the !$acc directives are
!< comments and the device kernels run on the host, so the figure shows the device code path, not a GPU run. The
!< placed dragon (chapters 2-3) is rebuilt here with the same calls, since the figure programs run before the tutorial
!< ones.
program fig_tut16
use fossil,     only : surface_stl_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                    ! points per direction
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
real(R8P),    parameter       :: LIFT = 1.e-3_R8P          ! the section drawn just above the cells (seen from +z)
type(surface_stl_object)      :: surface
type(curve_object)            :: section
type(vector_R8P)              :: bmin, bmax, e(4)
type(vector_R8P), allocatable :: points(:), a(:), b(:)
real(R8P),        allocatable :: d_host(:), d_device(:)
real(R8P)                     :: scale, dx, zc, dmax, uv(2)
character(80)                 :: caption
integer(I4P)                  :: u, i, j, k, p, kc, m, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! chapter 2: repair
call surface%sanitize(fill_holes=.true.)
call surface%save_into_file(file_name='fig16-dragon-clean.stl')
call surface%load_from_file(file_name='fig16-dragon-clean.stl', guess_format=.true.) ! chapter 3: place
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig16-dragon-placed.stl')
call surface%load_from_file(file_name='fig16-dragon-placed.stl', guess_format=.true.)

allocate(points(n**3), d_host(n**3), d_device(n**3))                               ! chapter 16: the batch
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + 1.2_R8P * (i - 1) / (n - 1)) * ex_R8P + (-0.6_R8P + 1.2_R8P * (j - 1) / (n - 1)) * ey_R8P + &
               (-0.6_R8P + 1.2_R8P * (k - 1) / (n - 1)) * ez_R8P
enddo ; enddo ; enddo
call surface%enter_device
call surface%distance_many_device(points=points, distances=d_device, is_signed=.true., is_square_root=.true., &
                                  sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'distance_many_device failed'
call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
call surface%exit_device

dx = 1.2_R8P / (n - 1)
kc = n / 2                                                                          ! the plane nearest z = 0
zc = -0.6_R8P + dx * (kc - 1)
p = (kc - 1) * n**2
dmax = maxval(abs(d_device(p + 1:p + n**2)))
m = count(d_device(p + 1:p + n**2) /= d_host(p + 1:p + n**2))
write(caption, '(A,F6.3,A,I0,A,I0,A)') 'distance_many_device, z = ', zc, ': ', m, ' of ', n**2, &
                                       ' differ from distance_many'
call fig_open(u, 'tut16-slice', trim(caption), 0._R8P, 90._R8P)
call fig_zoom(u, zc * ez_R8P, 0.3_R8P)                                             ! the middle of the plane
do j = 1, n ; do i = 1, n                                                           ! a cell about each point
   p = (kc - 1) * n**2 + (j - 1) * n + i
   call corners(points(p))
   call fig_quad(u, 'field', e(1), e(2), e(3), e(4), value=d_device(p) / (0.5_R8P * dmax))
   if (d_device(p) /= d_host(p)) call fig_box(u, 'accent', e(1), e(3))
enddo ; enddo
call slice_surface(surface, 3_I4P, zc, section, status=status)                     ! the dragon at that z
allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; m = 0
do i = 1, section%get_vertices_number()
   if (section%get_next(i) == 0_I4P) cycle
   m = m + 1
   uv = section%get_vertex(i) ; a(m)%x = uv(1) ; a(m)%y = uv(2) ; a(m)%z = zc + LIFT
   uv = section%get_vertex(section%get_next(i)) ; b(m)%x = uv(1) ; b(m)%y = uv(2) ; b(m)%z = zc + LIFT
enddo
call fig_segments(u, 'accent', a(1:m), b(1:m))
call fig_close(u)
contains
   subroutine corners(c)
   !< The corners of the cell of size dx about the point c, in its plane z = const, into e.
   type(vector_R8P), intent(in) :: c !< Point.

   e = c
   e(1)%x = c%x - 0.5_R8P * dx ; e(1)%y = c%y - 0.5_R8P * dx
   e(2)%x = c%x + 0.5_R8P * dx ; e(2)%y = c%y - 0.5_R8P * dx
   e(3)%x = c%x + 0.5_R8P * dx ; e(3)%y = c%y + 0.5_R8P * dx
   e(4)%x = c%x - 0.5_R8P * dx ; e(4)%y = c%y + 0.5_R8P * dx
   endsubroutine corners
endprogram fig_tut16
