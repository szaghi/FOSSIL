!run figck-inertia-world fig_ck_inertia_world
!fig ck-inertia-world
!< Figure of the cookbook recipe "Inertia of a turned body": the 0.1 x 0.2 x 0.3 m block before (ghost) and after
!< (body) the 30-degree turn about z, seen from above, and the principal frame of the turned block: its three axes
!< from the centre of mass, labelled with their moments.
program fig_ck_inertia_world
use fossil,     only : surface_stl_object, rigid_pose_t, principal_frame
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: block, turned
type(rigid_pose_t)       :: frame
type(vector_R8P)         :: c, a
real(R8P)                :: center(3), moments(3), axes(3,3)
integer(I4P)             :: u, i
character(12)            :: text

call block%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call block%resize(factor=0.1_R8P * ex_R8P + 0.2_R8P * ey_R8P + 0.3_R8P * ez_R8P) ; call block%analyze
turned = block
call turned%rotate(axis=ez_R8P, angle=PI / 6._R8P) ; call turned%analyze
call turned%mass_properties(center=center, principal_moments=moments, principal_axes=axes, density=7850._R8P)
frame = principal_frame(center, axes)
call fig_open(u, 'ck-inertia-world', 'the block turned 30 deg about z and its principal frame, moments in kg m^2', &
              -90._R8P, 90._R8P, height=380_I4P)
call fig_surface(u, 'ghost', block)
call fig_surface(u, 'body', turned)
c = vector_R8P(frame%t(1), frame%t(2), 0.31_R8P)
do i = 1, 3
   a = vector_R8P(frame%r(1,i), frame%r(2,i), frame%r(3,i))
   if (a%x + a%y < 0._R8P) a = a * (-1._R8P)                ! draw each axis to the upper right
   call fig_segment(u, 'accent', c, c + a * 0.12_R8P)
   write(text, '(F7.5)') moments(i)
   if (abs(a%z) > 0.9_R8P) then                              ! the vertical axis points at the viewer
      call fig_label(u, 'accent', c + vector_R8P(0.01_R8P, -0.035_R8P, 0._R8P), trim(text)//' (z)')
   else
      call fig_label(u, 'accent', c + a * 0.135_R8P, trim(text))
   endif
enddo
call fig_points(u, 'foot', [c])
call fig_close(u)
endprogram fig_ck_inertia_world
