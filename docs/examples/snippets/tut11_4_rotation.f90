!< Tutorial 11, step 4: the swept geometry of a rotation. A section of the placed dragon (chapters 2-3) that does not
!< cross itself, turned by 10 degrees about the origin (the exact rotation path of swept_cells_2d), then the unit cube
!< turned by 15 degrees about its vertical axis on 10^3 cells (swept_cells_3d). The figure: an animation of the section
!< turning in nine steps of 10 degrees, each frame the cells shaded by that step's swept area; every step must take
!< the exact path (issue #86: before its fix, most steps of these sections did not).
program tut11_4_rotation
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                       pose_from_axis_angle, pose_2d, slice_surface, swept_cells_2d_t, swept_cells_3d_t, &
                       swept_cells_2d, swept_cells_3d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(swept_cells_2d_t)   :: sw
real(R8P)                :: x(49)
integer(I4P)             :: i, n, status
integer(I4P), allocatable :: pairs(:,:)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 3_I4P, -0.2_R8P, section, merge_collinear=.true., status=status)   ! z = -0.2
call section%self_intersections(pairs, n)
print '(A,I0,A,I0)', 'section z = -0.2: ', section%get_vertices_number(), ' segments, crossing pairs ', n
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
call swept_cells_2d(section, pose_2d(angle=0._R8P), pose_2d(angle=pi / 18._R8P), x, x, sw)   ! 10 degrees
print '(A,L1,A,I0,A,I0)', 'rotation: exact ', sw%exact, ', cells swept ', count(sw%swept /= 0._R8P), &
                          ', largest evaluations per face ', sw%evaluations
print '(A,ES9.2,A,ES9.2)', '  sum of the swept areas ', sum(sw%swept), ', largest |residual| ', maxval(abs(sw%residual))
call solid
call figure
contains
   subroutine solid
   !< The unit cube turned by 15 degrees about the vertical axis through its centre.
   type(surface_stl_object) :: cube
   type(swept_cells_3d_t)   :: s3
   type(rigid_pose_t)       :: p1
   real(R8P)                :: g(11), c(3)

   call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)    ! [0, 1]^3
   g = [(-0.25_R8P + 0.15_R8P * i, i = 0, 10)]                              ! 10^3 cells of side 0.15
   c = 0.5_R8P                                                               ! the centre stays: t = c - R c
   p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=pi / 12._R8P)
   p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=pi / 12._R8P, translation=c - matmul(p1%r, c))
   call swept_cells_3d(cube, pose_identity(), p1, g, g, g, s3)
   print '(A,L1,A,I0,A,I0)', 'cube, 15 degrees: exact ', s3%exact, ', cells swept ', count(s3%swept /= 0._R8P), &
                             ', cleared ', count(s3%cleared)
   print '(A,ES9.2,A,ES9.2)', '  sum of the swept volumes ', sum(s3%swept), ', largest |residual| ', &
                              maxval(abs(s3%residual))
   endsubroutine solid

   subroutine figure
   !< The section turning in nine steps of 10 degrees, each one exact.
   integer(I4P), parameter       :: NS = 9
   real(R8P),    parameter       :: LIFT = 1.e-3_R8P
   type(curve_object)            :: moved
   type(rigid_pose_2d_t)         :: pa, pb
   type(swept_cells_2d_t)        :: sk
   type(vector_R8P), allocatable :: a(:), b(:)
   type(vector_R8P)              :: centre
   character(80)                 :: caption
   integer(I4P)                  :: u, s, ii, jj

   call fig_open(u, 'tut11-rotation', 'swept_cells_2d, rotation', 0._R8P, 90._R8P)
   centre%x = 0._R8P ; centre%y = 0._R8P ; centre%z = 0._R8P
   call fig_zoom(u, centre, 0.42_R8P)
   call fig_anim(u, 1._R8P)
   do s = 1, NS
      pa = pose_2d(angle=(s - 1) * pi / 18._R8P) ; pb = pose_2d(angle=s * pi / 18._R8P)
      if (s == 1) then
         sk = sw
      else
         call swept_cells_2d(section, pa, pb, x, x, sk)
      endif
      if (.not. sk%exact) error stop 'swept_cells_2d: a step off the exact path'
      print '(A,I0,A,L1,A,ES9.2)', 'step ', s, ': exact ', sk%exact, ', largest |residual| ', maxval(abs(sk%residual))
      write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'step ', s, ': ', 10 * (s - 1), ' to ', 10 * s, ' degrees, ', &
                                                count(sk%swept /= 0._R8P), ' cells swept'
      call fig_frame(u, trim(caption))
      do jj = 1, 48 ; do ii = 1, 48
         if (sk%swept(ii, jj) == 0._R8P) cycle
         call fig_quad(u, 'field', vec(x(ii), x(jj)), vec(x(ii + 1), x(jj)), vec(x(ii + 1), x(jj + 1)), &
                       vec(x(ii), x(jj + 1)), value=sk%swept(ii, jj) / 0.025_R8P**2)
      enddo ; enddo
      moved = section ; call moved%move(pa)
      call curve_segments(moved, LIFT, a, b)
      call fig_segments(u, 'box', a, b)
      moved = section ; call moved%move(pb)
      call curve_segments(moved, LIFT, a, b)
      call fig_segments(u, 'accent', a, b)
   enddo
   call fig_close(u)
   endsubroutine figure

   subroutine curve_segments(curve, level, a, b)
   !< The segments of a curve in the plane z = level.
   type(curve_object),            intent(in)  :: curve      !< Curve.
   real(R8P),                     intent(in)  :: level      !< Plane drawn.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Ends.
   real(R8P)                                  :: uv(2)
   integer(I4P)                               :: s, m

   allocate(a(curve%get_vertices_number()), b(curve%get_vertices_number())) ; m = 0
   do s = 1, curve%get_vertices_number()
      if (curve%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = curve%get_vertex(s) ; a(m) = vec(uv(1), uv(2), level)
      uv = curve%get_vertex(curve%get_next(s)) ; b(m) = vec(uv(1), uv(2), level)
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine curve_segments

   function vec(p, q, z) result(v)
   !< The point (p, q, z), z = 0 by default.
   real(R8P), intent(in)           :: p, q !< In-plane coordinates.
   real(R8P), intent(in), optional :: z    !< Height.
   type(vector_R8P)                :: v    !< Point.

   v%x = p ; v%y = q ; v%z = 0._R8P
   if (present(z)) v%z = z
   endfunction vec
endprogram tut11_4_rotation
