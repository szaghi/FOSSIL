!< Figures of the reference page "Boolean operations": the page's two unit cubes (b translated by 0.5 along the
!< diagonal) combined by a%boolean(b, op) for op = BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, from the same
!< camera, the volume and its error against the analytic one in the caption. The intersection is drawn inside both
!< inputs (translucent), the difference with the subtracted b (translucent).
program fig_ref_booleans
use fossil,     only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
call draw(BOOL_UNION,      'ref-boolean-union',      'BOOL_UNION',      1.875_R8P, .false., .false.)
call draw(BOOL_INTERSECT,  'ref-boolean-intersect',  'BOOL_INTERSECT',  0.125_R8P, .true.,  .true.)
call draw(BOOL_DIFFERENCE, 'ref-boolean-difference', 'BOOL_DIFFERENCE', 0.875_R8P, .false., .true.)
contains
   subroutine draw(op, name, label, exact, ghost_a, ghost_b)
   !< One boolean of the two cubes, drawn over the inputs `ghost_a`, `ghost_b` (translucent).
   integer(I4P), intent(in) :: op      !< Operation.
   character(*), intent(in) :: name    !< Figure name.
   character(*), intent(in) :: label   !< Operation name in the caption.
   real(R8P),    intent(in) :: exact   !< Analytic volume.
   logical,      intent(in) :: ghost_a !< Draw the input a too.
   logical,      intent(in) :: ghost_b !< Draw the input b too.
   type(surface_stl_object) :: a, b, a0, b0
   character(96)            :: caption
   integer(I4P)             :: u, status

   call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze
   call a0%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call b0%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call b0%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
   call a%boolean(other=b, op=op, status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   write(caption, '(A,A,A,I0,A,F6.4,A,ES7.1)') 'a%boolean(b, ', label, '): ', a%get_facets_number(), &
                                                ' facets, volume ', a%get_volume(), ', error ', abs(a%get_volume() - exact)
   call fig_open(u, name, trim(caption), -60._R8P, 25._R8P)
   if (ghost_a) call fig_surface(u, 'ghost', a0)
   if (ghost_b) call fig_surface(u, 'ghost', b0)
   call fig_surface(u, 'body', a)
   call fig_close(u)
   endsubroutine draw
endprogram fig_ref_booleans
