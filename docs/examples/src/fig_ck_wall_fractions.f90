!run figck-wall-fractions fig_ck_wall_fractions
!fig ck-wall-fractions
!< Figure of the cookbook recipe "Find where the wall cuts the links between cell centres": the disc of
!< ck_wall_fractions over its 16 x 16 cell centres, inside centres in teal, outside ones faint; every link the wall
!< cuts drawn from each end to the wall point its fraction gives (teal from inside, ember from outside), the wall
!< points in yellow.
program fig_ck_wall_fractions
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: PI = acos(-1._R8P), H = 0.0625_R8P
integer(I4P), parameter       :: DI(4) = [-1, 1, 0, 0], DJ(4) = [0, 0, -1, 1]
type(curve_object)            :: disc
real(R8P)                     :: u(48), v(48), x(16), frac(4, 16, 16)
logical                       :: inside(16, 16)
type(vector_R8P), allocatable :: a(:), b(:), ci(:), co(:), wi(:), wo(:), pi_(:), po(:), w(:)
type(vector_R8P)              :: c, p
integer(I4P)                  :: i, j, d, un

u = [(0.5_R8P + 0.3_R8P * cos(2 * PI * i / 48), i = 0, 47)]
v = [(0.5_R8P + 0.3_R8P * sin(2 * PI * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(H * (i - 0.5_R8P), i = 1, 16)]
call disc%cell_wall_fractions(x=x, y=x, frac=frac, inside=inside)
allocate(a(48), b(48))
do i = 1, 48
   a(i) = vector_R8P(u(i), v(i), 0._R8P) ; b(i) = vector_R8P(u(mod(i, 48) + 1), v(mod(i, 48) + 1), 0._R8P)
enddo
allocate(ci(0), wi(0), co(0), wo(0), pi_(0), po(0), w(0))
do j = 1, 16 ; do i = 1, 16
   c = vector_R8P(x(i), x(j), 0._R8P)
   if (inside(i, j)) then ; pi_ = [pi_, c] ; else ; po = [po, c] ; endif
   do d = 1, 4
      if (frac(d, i, j) < 0._R8P) cycle
      p = c + vector_R8P(DI(d) * H, DJ(d) * H, 0._R8P) * frac(d, i, j)
      w = [w, p]
      if (inside(i, j)) then ; ci = [ci, c] ; wi = [wi, p] ; else ; co = [co, c] ; wo = [wo, p] ; endif
   enddo
enddo ; enddo
call fig_open(un, 'ck-wall-fractions', 'cell_wall_fractions: the wall point on every link it cuts, from both ends', &
              -90._R8P, 90._R8P, height=440_I4P)
call fig_segments(un, 'body', a, b)
call fig_points(un, 'ghost', po)
call fig_points(un, 'inside', pi_)
call fig_segments(un, 'inside', ci, wi)
call fig_segments(un, 'outside', co, wo)
call fig_points(un, 'foot', w)
call fig_close(un)
endprogram fig_ck_wall_fractions
