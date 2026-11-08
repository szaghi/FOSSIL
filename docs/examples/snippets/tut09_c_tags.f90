!< Tutorial 9, step 3: signed tag bands on a block of 32^3 cells (tag_block), then the 2-D twin on the section of the
!< dragon at the mid-plane of layer 16 (curve_object%tag_block). The figures: layer 16 of the 3-D tags with that
!< section, and the 2-D tags of the section on the same 32 x 32 cells.
program tut09_c_tags
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
use fossil, only : surface_stl_object, curve_object, slice_surface, STATUS_OK
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nc = 32                     ! cells per direction of the tagged block
integer(I4P), parameter  :: kc = 16                     ! the layer compared with the 2-D twin
type(surface_stl_object) :: surface
type(curve_object)       :: section
real(R8P)                :: x(nc + 1), zc
integer(I4P)             :: tags(nc, nc, nc), tags2(nc, nc), i, status
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: uv(2)
integer(I4P)                  :: u, j, n

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]          ! cell edges of one 32^3 block over the domain
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A)', 'tag        -5     -4     -3     -2     -1      0      1      2      3      4      5'
print '(A,11I7)', '3-D   ', (count(tags == i), i = -5, 5)
zc = 0.5_R8P * (x(kc) + x(kc + 1))                      ! the mid-plane of layer 16
call slice_surface(surface, 3_I4P, zc, section, status=status)
call section%tag_block(x=x, y=x, band=4, tags=tags2, status=status)   ! the same 32 x 32 cells, in the plane
print '(A,11I7)', '2-D   ', (count(tags2 == i), i = -5, 5)
print '(A,I0,A,I0)', 'layer 16: cut cells in 3-D ', count(tags(:, :, kc) == 0), ', on its mid-plane section ', &
                     count(tags2 == 0)
print '(A,I0)', 'cells cut in 2-D but not in the layer: ', count(tags2 == 0 .and. tags(:, :, kc) /= 0)
n = section%get_vertices_number()
allocate(a(n), b(n))
do i = 1, n
   uv = section%get_vertex(i) ; a(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.001_R8P * ez_R8P
   uv = section%get_vertex(section%get_next(i)) ; b(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.001_R8P * ez_R8P
enddo
call fig_open(u, 'tut09-tags', 'tag_block, layer 16 of 32^3: cut (cream), inside (teal), outside rings (ember)', &
              -90._R8P, 90._R8P)
do j = 1, nc ; do i = 1, nc
   call fig_quad(u, 'field', p(x(i), x(j)), p(x(i + 1), x(j)), p(x(i + 1), x(j + 1)), p(x(i), x(j + 1)), &
                 value=shade(tags(i, j, kc)))
enddo ; enddo
call fig_segments(u, 'accent', a, b)
call fig_close(u)
call fig_open(u, 'tut09-tags-2d', 'curve tag_block, the section alone on the same cells: a thinner cut ring', &
              -90._R8P, 90._R8P)
do j = 1, nc ; do i = 1, nc
   call fig_quad(u, 'field', p(x(i), x(j)), p(x(i + 1), x(j)), p(x(i + 1), x(j + 1)), p(x(i), x(j + 1)), &
                 value=shade(tags2(i, j)))
enddo ; enddo
call fig_segments(u, 'accent', a, b)
call fig_close(u)
contains
   pure function shade(tag) result(v)
   !< Colour of a tag: 0 cream, then darker teal (inside) or ember (outside) ring by ring.
   integer(I4P), intent(in) :: tag !< Tag.
   real(R8P)                :: v   !< Value of the field style.

   v = 0._R8P
   if (tag /= 0_I4P) v = sign(0.3_R8P + 0.14_R8P * (abs(tag) - 1), real(tag, R8P))
   endfunction shade

   pure function p(xx, yy) result(w)
   !< The point (xx, yy, 0).
   real(R8P), intent(in) :: xx, yy !< Coordinates.
   type(vector_R8P)      :: w      !< Point.

   w%x = xx ; w%y = yy ; w%z = 0._R8P
   endfunction p
endprogram tut09_c_tags
