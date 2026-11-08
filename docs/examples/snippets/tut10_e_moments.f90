!< Tutorial 10, step 5: moments and quadrature rules of the pieces (curve%cut_cells(moments_order=, quadrature=)), then
!< a coarse level built from the fine one (cut_coarsen_2d) against the cut cells computed on the coarse grid.
!< The figure (outside the regions): the first piece of the tail cell (23, 41): its oriented boundary (red), the fan of
!< signed triangles from the cell centre, its area points (teal positive, ember negative weight), its wall points (gold).
program tut10_e_moments
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_close
use vecfor,     only : vector_R8P
use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_coarsen_2d, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)        :: section
type(cut_cells_2d_t)      :: cut, coarse, direct
type(cut_links_t)         :: links, coarse_links
real(R8P)                 :: x(49), h, xc(2), qx
integer(I4P), allocatable :: fine_to_coarse(:)
integer(I4P)              :: i, k, p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]
call section%cut_cells(x=x, y=x, cut=cut, moments_order=2_I4P, boundaries=.true., quadrature=.true.)
p = cut%first_piece(23, 41)                          ! the first piece of the tail cell (23, 41), step 2
xc = 0.5_R8P * [x(23) + x(24), x(41) + x(42)]        ! the cell centre: moments in X = (x - xc) / h, Y = (y - yc) / h
print '(A,6F9.5)', 'moments / h^2 (1, X, Y, X^2, XY, Y^2):', cut%piece_moments(:, p) / h**2
print '(A,2F9.5,A,2F9.5)', '  centroid from the moments ', xc + h * cut%piece_moments(2:3, p) / cut%piece_moments(1, p), &
                           ', piece_centroid ', cut%piece_centroid(:, p)
k = cut%quad_start(p + 1) - cut%quad_start(p)
print '(A,I0,A,I0,A,I0,A)', 'area rule: ', k, ' points (', &
      count(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1) < 0._R8P), ' negative), wall rule: ', &
      cut%wall_quad_start(p + 1) - cut%wall_quad_start(p), ' points'
qx = 0._R8P
do k = cut%quad_start(p), cut%quad_start(p + 1) - 1  ! integrate X^2 with the rule
   qx = qx + cut%quad_weight(k) * ((cut%quad_point(1, k) - xc(1)) / h)**2
enddo
print '(A,ES9.2,A,ES9.2)', '  sum of weights - area ', &
      sum(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1)) - cut%piece_area(p), &
      ', rule on X^2 - moment ', qx - cut%piece_moments(4, p)
call cut_links_2d(cut, links)
call cut_coarsen_2d(cut, links, x, x, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
call section%cut_cells(x=x(1:49:2), y=x(1:49:2), cut=direct)                      ! the coarse grid directly
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'coarse (ratio 2): status ', status, ', pieces ', coarse%np, ', direct ', direct%np, &
                               ', fine piece ', p, ' -> coarse piece ', fine_to_coarse(p)
print '(A,ES9.2)', '  largest |alpha coarsened - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha))
call figure
contains
   subroutine figure
   !< The tail cell: its cell and neighbours by alpha, the wall, the rules of its two pieces.
   type(vector_R8P), allocatable :: a(:), b(:), pts(:)
   real(R8P),        allocatable :: w(:)
   real(R8P)                     :: uv(2)
   integer(I4P)                  :: u, n, q, ii, jj, nn, nneg
   character(80)                 :: cap

   nn = cut%quad_start(p + 1) - cut%quad_start(p)
   nneg = count(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1) < 0._R8P)
   write(cap, '(A,I0,A,I0,A)') 'piece 1 of cell (23,41): ', nn, ' area points, ', nneg, &
                               ' negative (ember), wall points gold'
   call fig_open(u, 'tut10-quadrature', trim(cap), 0._R8P, 90._R8P)
   call fig_zoom(u, vec(xc(1), xc(2), 0._R8P), 0.75_R8P * h)
   do jj = 39, 43 ; do ii = 21, 25
      call fig_quad(u, 'frac', vec(x(ii), x(jj), 0._R8P), vec(x(ii + 1), x(jj), 0._R8P), &
                    vec(x(ii + 1), x(jj + 1), 0._R8P), vec(x(ii), x(jj + 1), 0._R8P), value=cut%alpha(ii, jj))
   enddo ; enddo
   call fig_segments(u, 'box', [(vec(x(ii), x(39), 1.e-3_R8P), ii = 21, 26), (vec(x(21), x(jj), 1.e-3_R8P), jj = 39, 44)], &
                     [(vec(x(ii), x(44), 1.e-3_R8P), ii = 21, 26), (vec(x(26), x(jj), 1.e-3_R8P), jj = 39, 44)])
   allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; n = 0
   do q = 1, section%get_vertices_number()
      if (section%get_next(q) == 0_I4P) cycle
      n = n + 1
      uv = section%get_vertex(q) ; a(n) = vec(uv(1), uv(2), 2.e-3_R8P)
      uv = section%get_vertex(section%get_next(q)) ; b(n) = vec(uv(1), uv(2), 2.e-3_R8P)
   enddo
   call fig_segments(u, 'body', a(1:n), b(1:n))
   n = cut%piece_seg_start(p + 1) - cut%piece_seg_start(p)                          ! the piece's oriented boundary
   a = [(vec(cut%piece_seg(1, q), cut%piece_seg(2, q), 2.5e-3_R8P), q = cut%piece_seg_start(p), cut%piece_seg_start(p + 1) - 1)]
   b = [(vec(cut%piece_seg(3, q), cut%piece_seg(4, q), 2.5e-3_R8P), q = cut%piece_seg_start(p), cut%piece_seg_start(p + 1) - 1)]
   call fig_segments(u, 'hole', a, b)
   call fig_segments(u, 'box', [(vec(xc(1), xc(2), 2.5e-3_R8P), q = 1, n)], a)    ! the fan from the cell centre
   pts = [(vec(cut%quad_point(1, q), cut%quad_point(2, q), 3.e-3_R8P), q = cut%quad_start(p), cut%quad_start(p + 1) - 1)]
   w = cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1)
   call fig_points(u, 'inside', pack(pts, w >= 0._R8P))
   call fig_points(u, 'outside', pack(pts, w < 0._R8P))
   pts = [(vec(cut%wall_quad_point(1, q), cut%wall_quad_point(2, q), 3.e-3_R8P), &
           q = cut%wall_quad_start(p), cut%wall_quad_start(p + 1) - 1)]
   call fig_points(u, 'accent', pts)
   call fig_close(u)
   endsubroutine figure

   function vec(p, q, z) result(v)
   !< The point (p, q, z).
   real(R8P), intent(in) :: p, q, z !< Coordinates.
   type(vector_R8P)      :: v       !< Point.

   v%x = p ; v%y = q ; v%z = z
   endfunction vec
endprogram tut10_e_moments
