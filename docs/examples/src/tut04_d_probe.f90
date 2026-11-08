!run tut04-probe tut04_d_probe
!fig anim-tut04-probe
!< Animation of tutorial 4: a probe travels an ellipse in the plane x = 0 through the placed dragon; each frame draws
!< the probe (teal inside, ember outside, by the fast winding number), its closest point (gold) and the distance
!< segment, with the path (grey) and the signed distance in the caption. The dragon is drawn as a copy decimated to 500 facets.
program tut04_d_probe
use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: NF = 24
real(R8P),    parameter  :: PI = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, light
type(vector_R8P)         :: p, foot
real(R8P)                :: d, t
character(64)            :: caption
integer(I4P)             :: u, k, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
light = surface
call light%decimate(target_facets=500_I4P, status=status)                 ! a light drawing
call fig_open(u, 'anim-tut04-probe', '', -35._R8P, 20._R8P)
call fig_anim(u, 0.5_R8P)
call fig_surface(u, 'ghost', light)
do k = 0, 2 * NF - 1                                                       ! the path
   call fig_segment(u, 'box', path(PI * k / NF), path(PI * (k + 1) / NF))
enddo
do k = 0, NF - 1
   t = 2._R8P * PI * k / NF
   p = path(t)
   call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., closest_point=foot, &
                                 sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   write(caption, '(A,2F7.3,A,F8.4,A)') 'probe (0,', p%y, p%z, '): d = ', d, merge(' inside ', ' outside', d < 0._R8P)
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call fig_segment(u, 'accent', p, foot)
   call fig_points(u, merge('inside ', 'outside', d < 0._R8P), [p])
   call fig_points(u, 'foot', [foot])
enddo
call fig_close(u)
contains
   function path(t) result(q)
   !< The point of the ellipse at angle t (components below 1e-12 set to 0, for the caption).
   real(R8P), intent(in) :: t !< Angle.
   type(vector_R8P)      :: q !< Point.

   q = merge(0._R8P, 0.42_R8P * cos(t), abs(cos(t)) < 1.e-12_R8P) * ey_R8P + &
       merge(0._R8P, 0.25_R8P * sin(t), abs(sin(t)) < 1.e-12_R8P) * ez_R8P
   endfunction path
endprogram tut04_d_probe
