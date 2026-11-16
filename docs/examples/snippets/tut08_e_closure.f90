!< Tutorial 8, fifth run: the closure table of a finite-difference solver in one call, on the closed wrap of the dragon
!< (the ghosts need an inside: the parities of the exact crossings), the 41^3 grid lines of the chapter, reach 2: the
!< ghost points, their feet and normals, an image line of two probes, the ILW crossings. The figure shows the grid
!< plane y = 0, seen along y: the ghosts of the plane (teal), joined through their feet (gold) to their probes (ember).
program tut08_e_closure
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_close
use vecfor,     only : vector_R8P
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
allocate(a(0), b(0), pg(0), pf(0), pp(0))
do g = 1, table%n
   if (abs(table%points(2, g)) > 0.5_R8P * h) cycle                     ! the plane y = 0
   pg = [pg, vector_R8P(table%points(1, g), 0._R8P, table%points(3, g))]
   pf = [pf, vector_R8P(table%foot(1, g), 0._R8P, table%foot(3, g))]
   pp = [pp, vector_R8P(table%probe(1, 2, g), 0._R8P, table%probe(3, 2, g))]
   a = [a, pg(size(pg))] ; b = [b, pp(size(pp))]
enddo
call fig_open(u, 'tut08-closure', 'plane y = 0: ghosts (teal), feet (gold), image lines of two probes (ember)', &
              -90._R8P, 0._R8P)
call fig_segments(u, 'ghost', a, b)
call fig_points(u, 'probe', pg) ; call fig_points(u, 'foot', pf) ; call fig_points(u, 'patch', pp)
call fig_close(u)
print '(A,I0)', 'ghosts in the plane y = 0: ', size(pg)
endprogram tut08_e_closure
