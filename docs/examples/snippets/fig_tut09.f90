!< Figures of tutorial 9: the 8^3 blocks box_classify_many finds cut and near the dragon, and one layer of cells of
!< the 32^3 block tag_block tags. The dragon is placed as in tutorials 2 and 3 (saved and read back as
!< dragon-placed.stl is); blocks, radius, cells and band are those of tutorial 9.
program fig_tut09
use fossil,     only : surface_stl_object, curve_object, slice_surface, BOX_CUT, BOX_NEAR
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_segment, fig_segments, fig_label, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
integer(I4P), parameter       :: nb = 8, nc = 32, kc = 16       ! blocks, cells per direction, the layer drawn
type(surface_stl_object)      :: surface
type(curve_object)            :: section
type(vector_R8P)              :: bmin, bmax, o
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: lo(3, nb**3), hi(3, nb**3), x(nc + 1), hb, scale, uv(2), zc
integer(I4P)                  :: codes(nb**3), tags(nc, nc, nc), u, i, j, k, n, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)        ! dragon-placed.stl, as tutorials 2-3
call surface%sanitize(fill_holes=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig09-dragon-placed.stl')
call surface%load_from_file(file_name='fig09-dragon-placed.stl', guess_format=.true.)

hb = 1.2_R8P / nb                                                                 ! the blocks of tutorial 9
n = 0
do k = 1, nb ; do j = 1, nb ; do i = 1, nb
   n = n + 1
   lo(:, n) = -0.6_R8P + [i - 1, j - 1, k - 1] * hb
   hi(:, n) = lo(:, n) + hb
enddo ; enddo ; enddo
call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.05_R8P)
print '(A,I0,A,I0)', 'blocks: cut ', count(codes == BOX_CUT), ', near ', count(codes == BOX_NEAR)
call fig_open(u, 'tut09-blocks', 'box_classify_many: cut blocks (ember), near blocks (teal)', -60._R8P, 22._R8P)
call fig_surface(u, 'body', surface)
do n = 1, nb**3
   if (codes(n) == BOX_CUT) call fig_box(u, 'outside', lo(1, n) * ex_R8P + lo(2, n) * ey_R8P + lo(3, n) * ez_R8P, &
                                                       hi(1, n) * ex_R8P + hi(2, n) * ey_R8P + hi(3, n) * ez_R8P)
   if (codes(n) == BOX_NEAR) call fig_box(u, 'probe', lo(1, n) * ex_R8P + lo(2, n) * ey_R8P + lo(3, n) * ez_R8P, &
                                                      hi(1, n) * ex_R8P + hi(2, n) * ey_R8P + hi(3, n) * ez_R8P)
enddo
call fig_close(u)

x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]                                    ! the tagged block of tutorial 9
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags)
zc = 0.5_R8P * (x(kc) + x(kc + 1))                                                ! the layer kc and its mid-plane
call slice_surface(surface, 3_I4P, zc, section, status=status)
n = section%get_vertices_number()
allocate(a(n), b(n))
do i = 1, n
   uv = section%get_vertex(i) ; a(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.001_R8P * ez_R8P
   uv = section%get_vertex(section%get_next(i)) ; b(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.001_R8P * ez_R8P
enddo
call fig_open(u, 'tut09-tags', 'tag_block, cells of layer 16: cut (cream), inside (teal), outside (ember)', &
              -90._R8P, 90._R8P)
do j = 1, nc
   do i = 1, nc
      call fig_quad(u, 'field', x(i) * ex_R8P + x(j) * ey_R8P, x(i + 1) * ex_R8P + x(j) * ey_R8P, &
                    x(i + 1) * ex_R8P + x(j + 1) * ey_R8P, x(i) * ex_R8P + x(j + 1) * ey_R8P, &
                    value=shade(tags(i, j, kc)))
   enddo
enddo
call fig_segments(u, 'accent', a, b)
o = (x(1) + 0.03_R8P) * ex_R8P + (x(1) + 0.03_R8P) * ey_R8P + 0.001_R8P * ez_R8P   ! axes glyph at a corner
call fig_segment(u, 'box', o, o + 0.1_R8P * ex_R8P) ; call fig_label(u, 'box', o + 0.1_R8P * ex_R8P, 'x')
call fig_segment(u, 'box', o, o + 0.1_R8P * ey_R8P) ; call fig_label(u, 'box', o + 0.1_R8P * ey_R8P, 'y')
call fig_close(u)
contains
   pure function shade(tag) result(v)
   !< Colour of a tag: 0 cream, then darker teal (inside) or ember (outside) ring by ring.
   integer(I4P), intent(in) :: tag !< Tag.
   real(R8P)                :: v   !< Value of the field style.

   v = 0._R8P
   if (tag /= 0_I4P) v = sign(0.3_R8P + 0.14_R8P * (abs(tag) - 1), real(tag, R8P))
   endfunction shade
endprogram fig_tut09
