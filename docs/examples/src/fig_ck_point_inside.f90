!run figck-point_inside fig_ck_point_inside
!fig ck-point_inside
!< Figure of the cookbook recipe "Tell whether a point is inside": the cube (ghost) and the two points of the recipe,
!< coloured by is_point_inside (teal inside, ember outside); the caption says whether the four sign algorithms agree.
program fig_ck_point_inside
use fossil,     only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                       SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter  :: algorithms(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                             SIGN_FAST_WINDING_NUMBER]
type(surface_stl_object) :: cube
type(vector_R8P)         :: points(2)
character(1)             :: names(2) = ['p', 'q']
logical                  :: inside(4, 2)
integer(I4P)             :: u, a, k
character(80)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
points(1) = vector_R8P(0.3_R8P, 0.6_R8P, 0.9_R8P) ; points(2) = vector_R8P(1.2_R8P, 0.5_R8P, 0.5_R8P)
do k = 1, 2 ; do a = 1, 4
   inside(a, k) = cube%is_point_inside(point=points(k), sign_algorithm=algorithms(a))
enddo ; enddo
if (all(inside(:, 1) .eqv. inside(1, 1)) .and. all(inside(:, 2) .eqv. inside(1, 2))) then
   caption = 'is_point_inside: the four sign algorithms agree on p and q'
else
   caption = 'is_point_inside: the sign algorithms disagree'
endif
call fig_open(u, 'ck-point_inside', trim(caption), -80._R8P, 18._R8P)
call fig_surface(u, 'ghost', cube)
call fig_box(u, 'box', cube%get_bmin(), cube%get_bmax())
do k = 1, 2
   if (inside(1, k)) then
      call fig_points(u, 'inside', points(k:k))
      call fig_label(u, 'inside', points(k), names(k) // ' inside')
   else
      call fig_points(u, 'outside', points(k:k))
      call fig_label(u, 'outside', points(k), names(k) // ' outside')
   endif
enddo
call fig_close(u)
endprogram fig_ck_point_inside
