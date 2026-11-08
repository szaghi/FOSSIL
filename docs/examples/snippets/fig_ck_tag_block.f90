!< Figure of the cookbook recipe "Tag refinement bands on a block": the tags of ck_tag_block (cube-01, 20^3 cells on
!< [-0.55, 1.55], band 3) on the cell layer k = 10 (the one through the cube's middle), seen from above, coloured by
!< tag / 4 (teal inside, cream on the cut cells, ember outside), each cell labelled with its tag; the cube's section
!< outlined.
program fig_ck_tag_block
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object) :: cube
type(vector_R8P)         :: a(4), b(4)
real(R8P)                :: x(21), zp
integer(I4P)             :: tags(20, 20, 20), i, j, u
character(3)             :: lab

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
x = [(-0.55_R8P + 0.105_R8P * i, i = 0, 20)]
call cube%tag_block(x=x, y=x, z=x, band=3, tags=tags)
zp = 0.5_R8P * (x(10) + x(11))
call fig_open(u, 'ck-tag-block', 'tag_block, band 3: the cell layer through the middle of the cube', -90._R8P, 90._R8P, &
              height=520_I4P)
do j = 1, 20 ; do i = 1, 20
   call fig_quad(u, 'field', p(x(i), x(j)), p(x(i + 1), x(j)), p(x(i + 1), x(j + 1)), p(x(i), x(j + 1)), &
                 value=tags(i, j, 10) / 4._R8P)
   write(lab, '(I0)') tags(i, j, 10)
   call fig_label(u, 'body', p(0.5_R8P * (x(i) + x(i + 1)) - 0.03_R8P, 0.5_R8P * (x(j) + x(j + 1)) - 0.02_R8P, 2.e-3_R8P), &
                  trim(lab))
enddo ; enddo
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
endprogram fig_ck_tag_block
