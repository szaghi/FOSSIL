!< Tutorial 11, step 5: the mass properties of the dragon. The placed dragon of chapters 2-3 is not closed, so its
!< closed surrogate, an alpha wrap, gives the mass, the centre of mass and the inertia tensor; the tensor at a pose
!< (inertia_to_world) equals the one of the moved wrap. The figure: the wrap (translucent), its principal axes at the
!< centre of mass labelled with their moments.
program tut11_5_inertia
use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, inertia_to_world, STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_label, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter  :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, wrap, rest
type(rigid_pose_t)       :: p1
real(R8P)                :: mass, volume, centre(3), inertia(3,3), inertia1(3,3), moments(3), axes(3,3)
integer(I4P)             :: status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%mass_properties(mass=mass, status=status)
print '(A,L1)', 'the placed dragon is a closed body: ', status == STATUS_OK     ! chapter 2: it is not
call surface%alpha_wrap(alpha=0.06_R8P, offset=0.006_R8P, wrapped=wrap, status=status)  ! a closed surrogate
print '(A,I0,A,L1)', 'alpha wrap: ', wrap%get_facets_number(), ' facets, watertight ', wrap%is_watertight()
call wrap%mass_properties(mass=mass, volume=volume, center=centre, inertia=inertia, principal_moments=moments, &
                          principal_axes=axes, density=1000._R8P, status=status)     ! a body of density 1000
if (status /= STATUS_OK) error stop 'the wrap is not closed'
print '(A,ES11.4,A,ES11.4)', 'volume ', volume, ', mass ', mass
print '(A,3F8.4)', 'centre of mass   ', centre
print '(A,3ES11.4)', 'principal moments', moments
rest = wrap                                                                ! kept for the figure
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call wrap%move(p1)
call wrap%mass_properties(inertia=inertia1, density=1000._R8P)
print '(A,ES9.2)', 'moved wrap against inertia_to_world, relative ', &
                   maxval(abs(inertia1 - inertia_to_world(inertia, p1))) / maxval(abs(inertia))
call figure
contains
   subroutine figure
   !< The wrap at rest and its principal axes at the centre of mass.
   type(vector_R8P)         :: c, e
   character(24)            :: text
   integer(I4P)             :: u, k

   c%x = centre(1) ; c%y = centre(2) ; c%z = centre(3)
   call fig_open(u, 'tut11-inertia', 'the alpha wrap of the dragon and its principal axes at the centre of mass', &
                 -35._R8P, 25._R8P)
   call fig_surface(u, 'ghost', rest)
   do k = 1, 3
      e%x = axes(1, k) ; e%y = axes(2, k) ; e%z = axes(3, k)
      e = 0.45_R8P * e
      call fig_segment(u, 'accent', c - e, c + e)
      write(text, '(A,I0,A,ES9.2)') 'I', k, ' = ', moments(k)
      call fig_label(u, 'accent', c + e, trim(text))
   enddo
   call fig_points(u, 'accent', [c])
   call fig_close(u)
   endsubroutine figure
endprogram tut11_5_inertia
