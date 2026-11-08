!< Figure of the cookbook recipe "Mass and inertia of a body": the 0.1 x 0.2 x 0.3 m steel block of the recipe
!< (translucent), its centre of mass, its three principal axes through it, each drawn to the block's faces and
!< labelled with its principal moment, and the corner (0, 0, 0) that the recipe's second call takes as `about`.
program fig_ck_mass_properties
use fossil,     only : surface_stl_object, STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object) :: block
real(R8P)                :: mass, center(3), moments(3), axes(3,3), half(3), a(3), b(3)
integer(I4P)             :: u, i, status
character(40)            :: text

call block%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call block%resize(factor=0.1_R8P * ex_R8P + 0.2_R8P * ey_R8P + 0.3_R8P * ez_R8P)
call block%analyze
call block%mass_properties(mass=mass, center=center, principal_moments=moments, principal_axes=axes, &
                           density=7850._R8P, status=status)
if (status /= STATUS_OK) error stop 'not a closed body'
half = [0.05_R8P, 0.1_R8P, 0.15_R8P]

call fig_open(u, 'ck-mass-properties', 'steel block 0.1 x 0.2 x 0.3 m: centre of mass, principal axes and moments', &
              -55._R8P, 25._R8P)
call fig_surface(u, 'ghost', block)
do i = 1, 3                                    ! each axis to the faces it meets: |axis . half| along the block's edges
   a = center - sum(abs(axes(:, i)) * half) * axes(:, i) ; b = center + sum(abs(axes(:, i)) * half) * axes(:, i)
   call fig_segment(u, 'accent', vector_R8P(a(1), a(2), a(3)), vector_R8P(b(1), b(2), b(3)))
   write(text, '(A,I0,A,F8.5,A)') 'I', i, ' = ', moments(i), ' kg m^2'
   if (i == 3) b = a                           ! the third label at the other end, clear of the centre's
   call fig_label(u, 'accent', vector_R8P(b(1), b(2), b(3)), trim(text))
enddo
call fig_points(u, 'patch', [vector_R8P(center(1), center(2), center(3))])
call fig_label(u, 'patch', vector_R8P(center(1), center(2), center(3)), 'centre of mass')
call fig_points(u, 'probe', [vector_R8P(0._R8P, 0._R8P, 0._R8P)])
call fig_label(u, 'probe', vector_R8P(0._R8P, 0._R8P, 0._R8P), 'about = (0, 0, 0)')
call fig_close(u)
endprogram fig_ck_mass_properties
