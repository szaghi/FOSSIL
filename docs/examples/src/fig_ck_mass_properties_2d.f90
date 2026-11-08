!run figck-mass-properties-2d fig_ck_mass_properties_2d
!fig ck-mass-properties-2d
!< Figure of the cookbook recipe "Second moments of area of a section": the hollow 0.1 x 0.2 m section of the recipe
!< (its wall filled, seen from above), the outer loop and the clockwise hole, its centroid and the two axes through it,
!< each drawn to the radius of gyration sqrt(I / A) on both sides and labelled with its second moment.
program fig_ck_mass_properties_2d
use fossil,     only : curve_object, STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curve_object) :: tube
real(R8P)          :: area, center(2), sm(3), g(2)
real(R8P)          :: uo(4), vo(4), ui(4), vi(4)
integer(I4P)       :: u, k, status
character(48)      :: text

uo = [-0.05_R8P, 0.05_R8P, 0.05_R8P, -0.05_R8P] ; vo = [-0.10_R8P, -0.10_R8P, 0.10_R8P, 0.10_R8P]
ui = [-0.04_R8P, -0.04_R8P, 0.04_R8P, 0.04_R8P] ; vi = [-0.09_R8P, 0.09_R8P, 0.09_R8P, -0.09_R8P]
call tube%adopt_loops(u=[uo, ui], v=[vo, vi], loop_sizes=[4_I4P, 4_I4P])
call tube%mass_properties(area=area, center=center, second_moments=sm, status=status)
if (status /= STATUS_OK) error stop 'not a closed section'
g = sqrt(sm(1:2) / area)                                       ! radii of gyration about v (from int u^2) and u

call fig_open(u, 'ck-mass-properties-2d', 'hollow section: centroid and radii of gyration sqrt(I / A)', &
              -90._R8P, 90._R8P, height=520_I4P)
call fig_quad(u, 'body', p(-0.05_R8P, -0.10_R8P), p(0.05_R8P, -0.10_R8P), p(0.05_R8P, -0.09_R8P), p(-0.05_R8P, -0.09_R8P))
call fig_quad(u, 'body', p(-0.05_R8P,  0.09_R8P), p(0.05_R8P,  0.09_R8P), p(0.05_R8P,  0.10_R8P), p(-0.05_R8P,  0.10_R8P))
call fig_quad(u, 'body', p(-0.05_R8P, -0.09_R8P), p(-0.04_R8P, -0.09_R8P), p(-0.04_R8P, 0.09_R8P), p(-0.05_R8P, 0.09_R8P))
call fig_quad(u, 'body', p( 0.04_R8P, -0.09_R8P), p( 0.05_R8P, -0.09_R8P), p( 0.05_R8P, 0.09_R8P), p( 0.04_R8P, 0.09_R8P))
call fig_segments(u, 'cell', [(p(uo(k), vo(k), 1e-4_R8P), k = 1, 4)], [(p(uo(mod(k, 4) + 1), vo(mod(k, 4) + 1), 1e-4_R8P), &
                  k = 1, 4)])
call fig_segments(u, 'hole', [(p(ui(k), vi(k), 1e-4_R8P), k = 1, 4)], [(p(ui(mod(k, 4) + 1), vi(mod(k, 4) + 1), 1e-4_R8P), &
                  k = 1, 4)])
call fig_segment(u, 'accent', p(center(1) - g(1), center(2), 2e-4_R8P), p(center(1) + g(1), center(2), 2e-4_R8P))
call fig_segment(u, 'accent', p(center(1), center(2) - g(2), 2e-4_R8P), p(center(1), center(2) + g(2), 2e-4_R8P))
write(text, '(A,ES10.3,A)') 'int v^2 dA = ', sm(2), ' m^4'
call fig_label(u, 'accent', p(center(1), center(2) + g(2), 2e-4_R8P), trim(text))
write(text, '(A,ES10.3,A)') 'int u^2 dA = ', sm(1), ' m^4'
call fig_label(u, 'accent', p(center(1) + g(1), center(2), 2e-4_R8P), trim(text))
call fig_points(u, 'patch', [p(center(1), center(2), 3e-4_R8P)])
call fig_label(u, 'patch', p(center(1), center(2) - 0.01_R8P, 3e-4_R8P), 'centroid')
call fig_close(u)
contains
   function p(x, y, dz) result(q)
   !< A point of the plane of the section, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: x, y !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz   !< Lift.
   type(vector_R8P)                :: q    !< Point.

   q = vector_R8P(x, y, 0._R8P)
   if (present(dz)) q%z = dz
   endfunction p
endprogram fig_ck_mass_properties_2d
