!< Tutorial 10, step 1: the exact crossings of one grid line with the dragon section of chapter 7 (line_crossings).
!< The figure (outside the regions): the section, the grid, the line with the most crossings, its solid intervals and
!< its crossings, entering in teal and leaving in ember.
program tut10_a_crossings
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_label, fig_close
use fossil,     only : curve_object, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(curve_object)        :: section
real(R8P)                 :: x(49), solid
real(R8P),    allocatable :: t(:)
integer(I4P), allocatable :: segments(:)
logical,      allocatable :: entering(:)
integer(I4P)              :: i, j, best, n, nbest, status

call section%load_from_file(file_name='dragon-section.dat', status=status)       ! chapter 7
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]          ! the cell edges of a 48 x 48 block, h = 0.025
best = 1 ; nbest = 0
do j = 1, 49                                          ! the grid line x = x(j), along y, with the most crossings
   call section%line_crossings(axis=2_I4P, a=x(j), t0=x(1), t1=x(49), t=t, segments=segments, entering=entering, n=n)
   if (n > nbest) then
      best = j ; nbest = n
   endif
enddo
call section%line_crossings(axis=2_I4P, a=x(best), t0=x(1), t1=x(49), t=t, segments=segments, entering=entering, n=n)
print '(A,I0,A,F7.4,A,I0,A)', 'grid line ', best, ' (x = ', x(best), '): ', n, ' crossings'
do i = 1, n
   print '(A,F9.5,A,I4,A,A)', '  y = ', t(i), '  segment ', segments(i), '  ', merge('enters', 'leaves', entering(i))
enddo
solid = sum(t(2:n:2) - t(1:n-1:2))                    ! in +y: enter, leave, enter, leave, ...
print '(A,F9.6)', 'solid length of the line: ', solid
call figure
contains
   subroutine figure
   !< The section, the grid, the line, its solid intervals (yellow) and crossings (teal entering, ember leaving).
   type(vector_R8P), allocatable :: a(:), b(:), p(:)
   type(vector_R8P)              :: c
   real(R8P)                     :: uv(2), lo(2), hi(2)
   integer(I4P)                  :: k, m, u
   character(80)                 :: cap

   write(cap, '(A,F6.3,A)') 'line_crossings on x = ', x(best), ': solid yellow, enters teal, leaves ember'
   call fig_open(u, 'tut10-crossings', trim(cap), 0._R8P, 90._R8P)
   lo = section%get_bmin() ; hi = section%get_bmax()
   c = vec(0.5_R8P * (lo(1) + hi(1)), 0.5_R8P * (lo(2) + hi(2)), 0._R8P)
   call fig_zoom(u, c, max(0.5_R8P * (hi(1) - lo(1)), 0.33_R8P * (hi(2) - lo(2))) + 0.03_R8P)
   allocate(a(98), b(98))                                                    ! the grid
   do k = 1, 49
      a(k) = vec(x(1), x(k), 0._R8P) ; b(k) = vec(x(49), x(k), 0._R8P)
      a(49 + k) = vec(x(k), x(1), 0._R8P) ; b(49 + k) = vec(x(k), x(49), 0._R8P)
   enddo
   call fig_segments(u, 'cell', a, b)
   deallocate(a, b) ; allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; m = 0
   do k = 1, section%get_vertices_number()                                   ! the wall
      if (section%get_next(k) == 0_I4P) cycle
      m = m + 1
      uv = section%get_vertex(k) ; a(m) = vec(uv(1), uv(2), 1.e-3_R8P)
      uv = section%get_vertex(section%get_next(k)) ; b(m) = vec(uv(1), uv(2), 1.e-3_R8P)
   enddo
   call fig_segments(u, 'body', a(1:m), b(1:m))
   call fig_segments(u, 'box', [vec(x(best), x(1), 2.e-3_R8P)], [vec(x(best), x(49), 2.e-3_R8P)])
   call fig_segments(u, 'accent', [(vec(x(best), t(k), 3.e-3_R8P), k = 1, n - 1, 2)], &
                     [(vec(x(best), t(k), 3.e-3_R8P), k = 2, n, 2)])
   p = [(vec(x(best), t(k), 4.e-3_R8P), k = 1, n)]
   call fig_points(u, 'inside', pack(p, entering(1:n)))
   call fig_points(u, 'outside', pack(p, .not. entering(1:n)))
   call fig_label(u, 'accent', vec(x(best) - 0.035_R8P, lo(2) - 0.07_R8P, 4.e-3_R8P), 'the grid line')
   call fig_close(u)
   endsubroutine figure

   function vec(p, q, z) result(v)
   !< The point (p, q, z).
   real(R8P), intent(in) :: p, q, z !< Coordinates.
   type(vector_R8P)      :: v       !< Point.

   v%x = p ; v%y = q ; v%z = z
   endfunction vec
endprogram tut10_a_crossings
