!< Figure of the cookbook recipe "Compute the signed distance on a whole grid": the field of ck_grid_distance (cube-01,
!< 41^3 nodes on [-0.5, 1.5]) on the grid plane z = 0.5, seen from above. Each node is drawn as its dual cell, coloured
!< by phi / 0.75; the band nodes (exact distances) are framed, the others come from fast sweeping. The section of the
!< cube by the same plane is drawn over it.
program fig_ck_grid_distance
use fossil,     only : distance_field_object, surface_stl_object, STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segment, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter     :: kp = 21                          ! the grid plane z = 0.5
type(surface_stl_object)    :: cube
type(distance_field_object) :: field
type(vector_R8P)            :: a(4), b(4), o
real(R8P)                   :: x(41), lo(41), hi(41), zp
integer(I4P)                :: i, j, u, status
character(80)               :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the recipe's grid and call
x = [(-0.5_R8P + 0.05_R8P * i, i = 0, 40)]
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=cube, status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
zp = x(kp)
lo = [x(1), (0.5_R8P * (x(i) + x(i + 1)), i = 1, 40)]           ! dual cells
hi = [(0.5_R8P * (x(i) + x(i + 1)), i = 1, 40), x(41)]

write(caption, '(A,F6.3,A,F6.3,A)') 'phi on the plane z = 0.5, from ', minval(field%phi(:, :, kp)), ' to ', &
                                    maxval(field%phi(:, :, kp)), '; band nodes framed'
call fig_open(u, 'ck-grid-distance', trim(caption), -90._R8P, 90._R8P, height=520_I4P)
do j = 1, 41 ; do i = 1, 41
   call fig_quad(u, 'field', p(lo(i), lo(j)), p(hi(i), lo(j)), p(hi(i), hi(j)), p(lo(i), hi(j)), &
                 value=field%phi(i, j, kp) / 0.75_R8P)
enddo ; enddo
do j = 1, 41 ; do i = 1, 41                                       ! the band: a small square in each band node
   if (.not. field%in_band(i, j, kp)) cycle
   call fig_segments(u, 'box', [p(x(i) - 0.01_R8P, x(j) - 0.01_R8P, 1.e-3_R8P), p(x(i) + 0.01_R8P, x(j) - 0.01_R8P, &
                     1.e-3_R8P), p(x(i) + 0.01_R8P, x(j) + 0.01_R8P, 1.e-3_R8P), p(x(i) - 0.01_R8P, x(j) + 0.01_R8P, &
                     1.e-3_R8P)], [p(x(i) + 0.01_R8P, x(j) - 0.01_R8P, 1.e-3_R8P), p(x(i) + 0.01_R8P, x(j) + 0.01_R8P, &
                     1.e-3_R8P), p(x(i) - 0.01_R8P, x(j) + 0.01_R8P, 1.e-3_R8P), p(x(i) - 0.01_R8P, x(j) - 0.01_R8P, &
                     1.e-3_R8P)])
enddo ; enddo
a = [p(0._R8P, 0._R8P, 2.e-3_R8P), p(1._R8P, 0._R8P, 2.e-3_R8P), p(1._R8P, 1._R8P, 2.e-3_R8P), p(0._R8P, 1._R8P, 2.e-3_R8P)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a, b)                              ! the section of the cube
write(caption, '(A,F7.4)') 'centre ', field%phi(21, 21, kp)
call fig_label(u, 'accent', p(0.5_R8P, 0.5_R8P, 3.e-3_R8P), trim(caption))
o = p(x(1) + 0.05_R8P, x(1) + 0.05_R8P, 3.e-3_R8P)                ! axes glyph at a corner
call fig_segment(u, 'accent', o, o + 0.2_R8P * ex_R8P) ; call fig_label(u, 'accent', o + 0.2_R8P * ex_R8P, 'x')
call fig_segment(u, 'accent', o, o + 0.2_R8P * ey_R8P) ; call fig_label(u, 'accent', o + 0.2_R8P * ey_R8P, 'y')
call fig_close(u)
contains
   function p(xx, yy, dz) result(q)
   !< A point of the plane z = zp, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: xx, yy !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = xx * ex_R8P + yy * ey_R8P + zp * ez_R8P
   if (present(dz)) q%z = q%z + dz
   endfunction p
endprogram fig_ck_grid_distance
