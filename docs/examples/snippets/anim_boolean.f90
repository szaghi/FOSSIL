!< Animation of the recipe "Unite two bodies": the two unit cubes of the recipe (b translated by 0.5 along the
!< diagonal), then a%boolean(b, op) for the four operations, each frame with its computed volume and the exact one.
!< The symmetric difference looks like the union from outside: it is drawn translucent, the facets of its inner void
!< (the cube [0.5, 1]^3 where a and b overlap, found by their centroids) in ember.
program anim_boolean
use fossil,     only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK
use fossil,     only : facet_object
use fossil_fig, only : fig_open, fig_surface, fig_triangle, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
character(*), parameter  :: NAMES(4) = ['BOOL_UNION     ', 'BOOL_INTERSECT ', 'BOOL_DIFFERENCE', 'BOOL_SYMDIFF   ']
integer(I4P), parameter  :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
real(R8P),    parameter  :: EXACT(4) = [1.875_R8P, 0.125_R8P, 0.875_R8P, 1.75_R8P] ! 2 - 1/8, 1/8, 1 - 1/8, 2 - 2/8
type(surface_stl_object)    :: a, b, r
type(facet_object), pointer :: f
type(vector_R8P)            :: c
character(80)               :: caption
integer(I4P)                :: u, k, i, status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze
call fig_open(u, 'anim-boolean', '', -60._R8P, 25._R8P)
call fig_anim(u, 1.6_R8P)
call fig_frame(u, caption='the inputs: a (translucent) and b, translated by 0.5 (ember)')
call fig_surface(u, 'patch', b)
call fig_surface(u, 'ghost', a)
do k = 1, 4
   r = a
   call r%boolean(other=b, op=OPS(k), status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   write(caption, '(A,A,A,F6.4,A,F6.4,A)') 'a%boolean(b, ', trim(NAMES(k)), '): volume ', r%get_volume(), &
                                           ' (exact ', EXACT(k), ')'
   print '(A,A,I0,A,L1)', trim(caption), ', ', r%get_facets_number(), ' facets, manifold ', r%is_manifold()
   if (k < 4) then
      call fig_frame(u, caption=trim(caption))
      call fig_surface(u, 'body', r)
   else
      write(caption, '(A,F6.4,A,F6.4,A)') 'BOOL_SYMDIFF: volume ', r%get_volume(), ' (exact ', EXACT(k), &
                                          '), its inner void in ember'
      call fig_frame(u, caption=trim(caption), seconds=2.6_R8P)
      do i = 1, r%get_facets_number()
         f => r%facet_at(i)
         c = (f%vertex(1) + f%vertex(2) + f%vertex(3)) / 3._R8P
         if (min(c%x, c%y, c%z) > 0.5_R8P - 1.e-9_R8P .and. max(c%x, c%y, c%z) < 1._R8P + 1.e-9_R8P) then
            call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
         else
            call fig_triangle(u, 'ghost', f%vertex(1), f%vertex(2), f%vertex(3))
         endif
      enddo
   endif
enddo
call fig_close(u)
endprogram anim_boolean
