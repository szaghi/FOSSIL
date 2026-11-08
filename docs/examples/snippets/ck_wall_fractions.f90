program ck_wall_fractions
use fossil, only : curve_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: disc
real(R8P)            :: u(48), v(48), x(16), frac(4, 16, 16)
logical              :: inside(16, 16)
integer(I4P)         :: i, status, bad

u = [(0.5_R8P + 0.3_R8P * cos(2 * PI * i / 48), i = 0, 47)]   ! a disc of radius 0.3, 48 segments
v = [(0.5_R8P + 0.3_R8P * sin(2 * PI * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(0.0625_R8P * (i - 0.5_R8P), i = 1, 16)]                 ! the cell centres of 16 x 16 cells
call disc%cell_wall_fractions(x=x, y=x, frac=frac, inside=inside, bad_lines=bad, status=status)
if (status /= STATUS_OK) error stop 'cell_wall_fractions failed'
print '(A,I0,A,I0,A,I0)', 'centres inside ', count(inside), ', fractions >= 0 (links cut by the wall) ', count(frac >= 0._R8P), &
                          ', bad lines ', bad
print '(A,L1,A,4F8.4)', 'cell (3, 8): inside ', inside(3, 8), ', fractions -u +u -v +v', frac(:, 3, 8)
print '(A,L1,A,4F8.4)', 'cell (4, 8): inside ', inside(4, 8), ', fractions -u +u -v +v', frac(:, 4, 8)
endprogram ck_wall_fractions
