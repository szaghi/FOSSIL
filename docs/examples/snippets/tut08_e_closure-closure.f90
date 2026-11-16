use fossil,     only : surface_stl_object, ib_closure_t, fd_closure_3d
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter  :: n = 41                                      ! grid points per direction
type(surface_stl_object) :: surface, wrap
type(ib_closure_t)       :: table
real(R8P)                :: h, x(n)
integer(I4P)             :: i, g, u
type(vector_R8P), allocatable :: a(:), b(:), pg(:), pf(:), pp(:)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap)  ! one closed surface (see the last run)
h = 1.2_R8P / (n - 1)
x = [(-0.6_R8P + (i - 1) * h, i = 1, n)]                               ! the grid lines, the same along x, y, z
call fd_closure_3d(wrap, x, x, x, [2_I4P, 2_I4P, 2_I4P], table, lengths=[h, 2 * h], ilw=.true.)
print '(A,I0,A,I0,A,I0)', 'ghost points ', table%n, ', probes visible ', count(table%visible), ' of ', table%n * table%m
print '(A,I0,A,F8.5,A,F8.5)', 'ILW crossings ', count(table%ilw_element > 0), ', distances from ', minval(table%distance), &
      ' to ', maxval(table%distance)
