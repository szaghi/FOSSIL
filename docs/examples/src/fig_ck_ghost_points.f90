!run figck-ghost-points fig_ck_ghost_points
!fig ck-ghost-points
!< Figure of the cookbook recipe "Find the ghost points of a finite-difference stencil": the ellipse of ck_ghost_points
!< on its 20 x 20 cell centres (h = 0.05), seen from above: the first ghost layer in teal, the second in gold, the
!< other solid points in stone.
program fig_ck_ghost_points
use fossil,     only : curve_object, ghost_points_t, ghost_points_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: h = 0.05_R8P
type(curve_object)            :: body
type(ghost_points_t)          :: gp
real(R8P)                     :: u(64), v(64), xc(20)
type(vector_R8P), allocatable :: a(:), b(:), p1(:), p2(:), ps(:)
logical                       :: gh(20, 20)
integer(I4P)                  :: i, j, g, un

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xc = [((i - 0.5_R8P) * h, i = 1, 20)]
call ghost_points_2d(body, xc, xc, [2_I4P, 2_I4P], gp)
allocate(p1(0), p2(0), ps(0), a(0)) ; gh = .false.
do g = 1, gp%n
   gh(gp%ijk(1, g), gp%ijk(2, g)) = .true.
   if (any(gp%steps(:, g) == 1)) then
      p1 = [p1, vector_R8P(gp%points(1, g), gp%points(2, g), 1.e-3_R8P)]
   else
      p2 = [p2, vector_R8P(gp%points(1, g), gp%points(2, g), 1.e-3_R8P)]
   endif
enddo
do j = 1, 20 ; do i = 1, 20
   if (gh(i, j)) cycle
   if (body%is_point_inside([xc(i), xc(j)])) then
      ps = [ps, vector_R8P(xc(i), xc(j), 1.e-3_R8P)]
   else
      a = [a, vector_R8P(xc(i), xc(j), 0._R8P)]
   endif
enddo ; enddo
print '(A,I0,A,I0,A,I0)', 'ghosts: layer 1 ', size(p1), ', layer 2 ', size(p2), ', other solid points ', size(ps)
call fig_open(un, 'ck-ghost-points', 'cell centres h = 0.05, reach 2: ghost layers 1 (teal) and 2 (gold), solid (stone)', &
              -90._R8P, 90._R8P)
call fig_points(un, 'box', a)
a = [(vector_R8P(u(i), v(i), 0._R8P), i = 1, 64)] ; b = [(vector_R8P(u(mod(i, 64) + 1), v(mod(i, 64) + 1), 0._R8P), i = 1, 64)]
call fig_segments(un, 'accent', a, b)
call fig_points(un, 'probe', p1)
call fig_points(un, 'foot', p2)
call fig_points(un, 'body', ps)
call fig_close(un)
endprogram fig_ck_ghost_points
