!< Figures of the reference page "Finite-difference wall closure":
!< - ghosts: the NACA 0012 nose on cell-centred lines of spacing chord/64, reach 3: the three ghost layers;
!< - masks: a few ghosts on the upper side, each joined to the first fluid point of every direction its mask sets;
!< - low side: the unit square on lines through its sides: the points on u = 0 and v = 0 are solid, those on u = 1 and
!<   v = 1 fluid, as the cut cells (x_i, x_i+1] put the sides; the cells shaded by their fluid fraction;
!< - lines: a circle on node lines and on cell-centre lines of the same spacing: two different ghost sets;
!< - sphere: the ghosts of a marching-cubes sphere on node lines, reach 2, the lower half.
program fig_ref_fd_closure
use fossil,     only : surface_stl_object, curve_object, slice_surface, extract_isosurface, facet_object, cut_cells_2d_t, &
                       ghost_points_t, ghost_points_2d, ghost_points_3d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_label, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),        parameter     :: LIFT = 1.e-3_R8P
type(surface_stl_object)        :: wing, sphere
type(curve_object)              :: section, square, circle
type(ghost_points_t)            :: gp, gq
type(cut_cells_2d_t)            :: cut
type(facet_object), allocatable :: facets(:)
type(vector_R8P),   allocatable :: a(:), b(:), p1(:), p2(:), p3(:)
real(R8P)                       :: h, xc(80), yc(16), e(9), xn(11), xm(10), xs(25), f(33, 33, 33)
integer(I4P)                    :: u, i, j, k, g, d, s, ax, sg, nl(3)

! ---- the NACA 0012 nose: three ghost layers on cell-centred lines
call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)
h = 1._R8P / 64._R8P
xc = [(-0.375_R8P + (i - 0.5_R8P) * h, i = 1, 80)] ; yc = [(-0.125_R8P + (j - 0.5_R8P) * h, j = 1, 16)]
call ghost_points_2d(section, xc, yc, [3_I4P, 3_I4P], gp)
allocate(p1(0), p2(0), p3(0)) ; nl = 0
do g = 1, gp%n
   k = minval(gp%steps(:, g), mask=gp%steps(:, g) > 0) ; nl(k) = nl(k) + 1
   select case (k)
   case (1) ; p1 = [p1, vec(gp%points(1, g), gp%points(2, g), LIFT)]
   case (2) ; p2 = [p2, vec(gp%points(1, g), gp%points(2, g), LIFT)]
   case default ; p3 = [p3, vec(gp%points(1, g), gp%points(2, g), LIFT)]
   endselect
enddo
print '(A,I0,A,3(I0,:,", "))', 'NACA 0012, reach 3: ghosts ', gp%n, ', layers ', nl
call fig_open(u, 'ref-fd-closure-ghosts', 'NACA 0012 nose, h = c/64, reach 3: layers 1 (teal), 2 (gold), 3 (ember)', &
              -90._R8P, 90._R8P)
call fig_zoom(u, vec(-0.17_R8P, 0._R8P), 0.13_R8P)
allocate(a(0))
do j = 1, size(yc) ; do i = 1, size(xc)
   a = [a, vec(xc(i), yc(j))]
enddo ; enddo
call fig_points(u, 'box', a)
call curve_segments(section, a, b)
call fig_segments(u, 'accent', a, b)
call fig_points(u, 'probe', p1)
call fig_points(u, 'foot', p2)
call fig_points(u, 'patch', p3)
call fig_close(u)

! ---- direction masks: ghosts on the upper side near x = -0.1, joined to the first fluid point of each set direction
call fig_open(u, 'ref-fd-closure-masks', 'ghosts (teal) and the first fluid point (ember) of each direction of their mask', &
              -90._R8P, 90._R8P)
call fig_zoom(u, vec(-0.1_R8P, 0.045_R8P), 0.06_R8P)
deallocate(a) ; allocate(a(0))
do j = 1, size(yc) ; do i = 1, size(xc)
   a = [a, vec(xc(i), yc(j))]
enddo ; enddo
call fig_points(u, 'box', a)
call curve_segments(section, a, b)
call fig_segments(u, 'accent', a, b)
deallocate(a, b, p1, p2) ; allocate(a(0), b(0), p1(0), p2(0))
do g = 1, gp%n
   if (gp%points(1, g) < -0.14_R8P .or. gp%points(1, g) > -0.06_R8P .or. gp%points(2, g) < 0._R8P) cycle
   p1 = [p1, vec(gp%points(1, g), gp%points(2, g), 2 * LIFT)]
   do d = 1, 4
      s = gp%steps(d, g) ; if (s == 0) cycle
      ax = (d + 1) / 2 ; sg = merge(-1, 1, mod(d, 2) == 1)
      i = gp%ijk(1, g) ; j = gp%ijk(2, g)
      if (ax == 1) i = i + sg * s
      if (ax == 2) j = j + sg * s
      a = [a, vec(gp%points(1, g), gp%points(2, g), LIFT)] ; b = [b, vec(xc(i), yc(j), LIFT)]
      p2 = [p2, vec(xc(i), yc(j), 2 * LIFT)]
   enddo
enddo
call fig_segments(u, 'ghost', a, b)
call fig_points(u, 'inside', p1)
call fig_points(u, 'outside', p2)
call fig_close(u)

! ---- the low-side rule: the unit square on lines through its sides
call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
e = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call square%cut_cells(x=e, y=e, cut=cut)
call ghost_points_2d(square, e, e, [1_I4P, 1_I4P], gq)
print '(A,I0,A,I0)', 'unit square on its lines: solid points ', count(cut%node_in), ', ghosts (reach 1) ', gq%n
call fig_open(u, 'ref-fd-closure-low-side', 'square on grid lines: solid points (teal: ghosts), fluid points (stone)', &
              -90._R8P, 90._R8P)
do j = 1, 8 ; do i = 1, 8
   call fig_quad(u, 'frac', vec(e(i), e(j)), vec(e(i + 1), e(j)), vec(e(i + 1), e(j + 1)), vec(e(i), e(j + 1)), &
                 value=cut%alpha(i, j))
enddo ; enddo
call grid_segments(e, e, a, b)
call fig_segments(u, 'box', a, b)
a = [vec(0._R8P, 0._R8P, LIFT), vec(1._R8P, 0._R8P, LIFT), vec(1._R8P, 1._R8P, LIFT), vec(0._R8P, 1._R8P, LIFT)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a, b)
call point_sets(e, e, cut%node_in, gq, p1, p2, p3)
call fig_points(u, 'inside', p1)
call fig_points(u, 'body', p2)
call fig_points(u, 'box', p3)
call fig_label(u, 'accent', vec(-0.45_R8P, -0.42_R8P, 3 * LIFT), 'u = 0, v = 0 (left, bottom) at + eps: solid')
call fig_label(u, 'accent', vec(-0.45_R8P, 1.38_R8P, 3 * LIFT), 'u = 1, v = 1 (right, top) at + eps: fluid')
call fig_close(u)

! ---- the solver's own lines: a circle on node lines and on cell-centre lines of the same spacing
call circle%adopt_loops(u=[(0.5_R8P + 0.31_R8P * cos(2 * acos(-1._R8P) * i / 96), i = 0, 95)], &
                        v=[(0.5_R8P + 0.31_R8P * sin(2 * acos(-1._R8P) * i / 96), i = 0, 95)], loop_sizes=[96_I4P])
xn = [(0.1_R8P * i, i = 0, 10)] ; xm = [(0.05_R8P + 0.1_R8P * i, i = 0, 9)]
call ghost_points_2d(circle, xn, xn, [1_I4P, 1_I4P], gp)
call ghost_points_2d(circle, xm, xm, [1_I4P, 1_I4P], gq)
print '(A,I0,A,I0)', 'circle, reach 1: ghosts on node lines ', gp%n, ', on cell-centre lines ', gq%n
call fig_open(u, 'ref-fd-closure-lines', 'one circle, two point sets: node ghosts (teal), cell-centre ghosts (ember)', &
              -90._R8P, 90._R8P)
call grid_segments(xn, xn, a, b)
call fig_segments(u, 'box', a, b)
call curve_segments(circle, a, b)
call fig_segments(u, 'accent', a, b)
deallocate(p1, p2) ; allocate(p1(gp%n), p2(gq%n))
p1 = [(vec(gp%points(1, g), gp%points(2, g), LIFT), g = 1, gp%n)]
p2 = [(vec(gq%points(1, g), gq%points(2, g), LIFT), g = 1, gq%n)]
call fig_points(u, 'inside', p1)
call fig_points(u, 'outside', p2)
call fig_close(u)

! ---- 3-D: a sphere on node lines, reach 2
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vec(-1._R8P, -1._R8P, -1._R8P), bmax=vec(1._R8P, 1._R8P, 1._R8P), iso=0._R8P, &
                        surface=facets)
call sphere%adopt_facets(facets=facets)
xs = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]
call ghost_points_3d(sphere, xs, xs, xs, [2_I4P, 2_I4P, 2_I4P], gp)
deallocate(p1, p2) ; allocate(p1(0), p2(0))
do g = 1, gp%n
   if (gp%points(3, g) > 0._R8P) cycle
   if (any(gp%steps(:, g) == 1)) then
      p1 = [p1, vec(gp%points(1, g), gp%points(2, g), gp%points(3, g))]
   else
      p2 = [p2, vec(gp%points(1, g), gp%points(2, g), gp%points(3, g))]
   endif
enddo
print '(A,I0,A,I0,A,I0)', 'sphere, reach 2: ghosts ', gp%n, ', lower half: layer 1 ', size(p1), ', layer 2 ', size(p2)
call fig_open(u, 'ref-fd-closure-sphere', 'sphere on node lines, reach 2, lower half: layers 1 (teal) and 2 (gold)', &
              -35._R8P, 20._R8P)
call fig_surface(u, 'ghost', sphere)
call fig_points(u, 'probe', p1)
call fig_points(u, 'foot', p2)
call fig_close(u)
contains
   pure function vec(x, y, z) result(v)
   !< A vector, z = 0 by default.
   real(R8P), intent(in)           :: x, y
   real(R8P), intent(in), optional :: z
   type(vector_R8P)                :: v

   v%x = x ; v%y = y ; v%z = 0._R8P
   if (present(z)) v%z = z
   endfunction vec

   subroutine curve_segments(c, a, b)
   !< The segments of a closed curve, lifted.
   type(curve_object),            intent(in)  :: c    !< Curve.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Segment ends.
   real(R8P)                                  :: p(2), q(2) !< Work.
   integer(I4P)                               :: i, n, m    !< Work.

   n = c%get_vertices_number() ; allocate(a(0), b(0))
   do i = 1, n
      m = c%get_next(i) ; if (m == 0) cycle
      p = c%get_vertex(i) ; q = c%get_vertex(m)
      a = [a, vec(p(1), p(2), LIFT)] ; b = [b, vec(q(1), q(2), LIFT)]
   enddo
   endsubroutine curve_segments

   subroutine grid_segments(x, y, a, b)
   !< The grid lines through the points x, y.
   real(R8P),                     intent(in)  :: x(:), y(:) !< Lines.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Segment ends.
   integer(I4P)                               :: i          !< Work.

   allocate(a(0), b(0))
   do i = 1, size(x)
      a = [a, vec(x(i), y(1), 0.5_R8P * LIFT)] ; b = [b, vec(x(i), y(size(y)), 0.5_R8P * LIFT)]
   enddo
   do i = 1, size(y)
      a = [a, vec(x(1), y(i), 0.5_R8P * LIFT)] ; b = [b, vec(x(size(x)), y(i), 0.5_R8P * LIFT)]
   enddo
   endsubroutine grid_segments

   subroutine point_sets(x, y, solid, gp, pg, ps, pf)
   !< The ghost points, the other solid points and the fluid points.
   real(R8P),                     intent(in)  :: x(:), y(:)  !< Lines.
   logical,                       intent(in)  :: solid(:,:)  !< Point flags.
   type(ghost_points_t),          intent(in)  :: gp          !< Ghosts.
   type(vector_R8P), allocatable, intent(out) :: pg(:), ps(:), pf(:) !< Points.
   logical, allocatable                       :: gh(:,:)     !< Ghost flags.
   integer(I4P)                               :: i, j, g     !< Work.

   allocate(gh(size(x), size(y)), pg(0), ps(0), pf(0)) ; gh = .false.
   do g = 1, gp%n
      gh(gp%ijk(1, g), gp%ijk(2, g)) = .true.
   enddo
   do j = 1, size(y) ; do i = 1, size(x)
      if (gh(i, j)) then
         pg = [pg, vec(x(i), y(j), 2 * LIFT)]
      elseif (solid(i, j)) then
         ps = [ps, vec(x(i), y(j), 2 * LIFT)]
      else
         pf = [pf, vec(x(i), y(j), 2 * LIFT)]
      endif
   enddo ; enddo
   endsubroutine point_sets
endprogram fig_ref_fd_closure
