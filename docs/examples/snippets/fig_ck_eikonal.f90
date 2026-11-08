!< Figure of the cookbook recipe "Check a distance field by the eikonal diagnostic": q = |1 - |grad phi|| of ck_eikonal
!< (cube-01, 41^3 nodes on [-0.5, 1.5]) on the grid plane z = 0.5, seen from above, each node drawn as its dual cell
!< coloured by min(q / 0.25, 1); the flagged nodes (q > 0.05) are dots, the section of the cube is outlined.
program fig_ck_eikonal
use fossil,     only : distance_field_object, surface_stl_object
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
implicit none
integer(I4P), parameter       :: kp = 21                        ! the grid plane z = 0.5
type(surface_stl_object)      :: cube
type(distance_field_object)   :: field
type(vector_R8P)              :: a(4), b(4)
type(vector_R8P), allocatable :: dots(:)
real(R8P)                     :: x(41), lo(41), hi(41), q(41, 41, 41), zp
integer(I4P)                  :: i, j, u, flagged
character(100)                :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
x = [(-0.5_R8P + 0.05_R8P * i, i = 0, 40)]
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=cube)
call field%gradient_diagnostic(q=q, flagged=flagged)
zp = x(kp)
lo = [x(1), (0.5_R8P * (x(i) + x(i + 1)), i = 1, 40)]
hi = [(0.5_R8P * (x(i) + x(i + 1)), i = 1, 40), x(41)]
write(caption, '(A,I0,A)') 'q = |1 - |grad phi|| on z = 0.5; dots: the ', count(q(:, :, kp) > 0.05_R8P), &
                           ' flagged nodes of the plane'
call fig_open(u, 'ck-eikonal', trim(caption), -90._R8P, 90._R8P, height=520_I4P)
do j = 1, 41 ; do i = 1, 41
   call fig_quad(u, 'frac', p(lo(i), lo(j)), p(hi(i), lo(j)), p(hi(i), hi(j)), p(lo(i), hi(j)), &
                 value=min(q(i, j, kp) / 0.25_R8P, 1._R8P))
enddo ; enddo
allocate(dots(0))
do j = 1, 41 ; do i = 1, 41
   if (q(i, j, kp) > 0.05_R8P) dots = [dots, p(x(i), x(j), 2.e-3_R8P)]
enddo ; enddo
call fig_points(u, 'outside', dots)
a = [p(0._R8P, 0._R8P, 3.e-3_R8P), p(1._R8P, 0._R8P, 3.e-3_R8P), p(1._R8P, 1._R8P, 3.e-3_R8P), p(0._R8P, 1._R8P, 3.e-3_R8P)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a, b)
call fig_close(u)
contains
   function p(xx, yy, dz) result(r)
   !< A point of the plane z = zp, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: xx, yy !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: r      !< Point.

   r = xx * ex_R8P + yy * ey_R8P + zp * ez_R8P
   if (present(dz)) r%z = r%z + dz
   endfunction p
endprogram fig_ck_eikonal
