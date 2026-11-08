!< Tutorial 11, step 3: the swept geometry of a step on a fixed grid. The dragon's section of chapter 7 translated by
!< 0.02 (the exact translation path); the residual of its cells, and the check of the input. The figure: each cell
!< shaded by its swept area, the section at pose0 (grey) and at pose1 (yellow).
program tut11_3_swept
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_box, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(curve_object)      :: section
type(rigid_pose_2d_t)   :: q0, q1
type(swept_cells_2d_t)  :: sw
real(R8P)               :: x(49)
integer(I4P)            :: i

call section%load_from_file(file_name='dragon-section.dat')            ! chapter 7: the section z = 0
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]                           ! 48 x 48 cells of side 0.025
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.02_R8P, 0._R8P])            ! a pure translation, 0.8 of a cell
call swept_cells_2d(section, q0, q1, x, x, sw)
print '(A,L1,A,I0,A,I0,A,I0)', 'translation: exact ', sw%exact, ', cells swept ', count(sw%swept /= 0._R8P), &
                               ', cleared ', count(sw%cleared), ', covered ', count(sw%covered)
print '(A,ES9.2)', '  sum of the swept areas ', sum(sw%swept)
print '(A,I0,A,ES9.2,A,2I3)', '  cells with residual > 1e-14: ', count(abs(sw%residual) > 1.e-14_R8P), &
                              ', largest ', maxval(abs(sw%residual)), ' in cell', maxloc(abs(sw%residual))
call swept_cells_2d(section, q0, q1, x, x, sw, check_input=.true.)    ! the section crosses itself (chapter 7)
print '(A,L1)', '  check_input: self-intersecting ', sw%self_intersecting
call figure
contains
   subroutine figure
   !< The translated step: the cells shaded by their swept area, the section at the two poses.
   real(R8P),    parameter       :: LIFT = 1.e-3_R8P
   type(curve_object)            :: moved
   type(vector_R8P), allocatable :: a(:), b(:)
   type(vector_R8P)              :: centre
   character(40)                 :: text
   integer(I4P)                  :: u, k(2)

   call fig_open(u, 'tut11-swept', 'swept_cells_2d: swept area (ember gained, teal lost), pose0 grey, pose1 yellow', &
                 0._R8P, 90._R8P)
   call curve_segments(section, LIFT, a, b)
   centre%x = 0.5_R8P * (minval(a%x) + maxval(a%x) + 0.02_R8P) ; centre%y = 0.5_R8P * (minval(a%y) + maxval(a%y))
   centre%z = 0._R8P
   call fig_zoom(u, centre, max(0.5_R8P * (maxval(a%x) - minval(a%x) + 0.02_R8P), &
                                0.3_R8P * (maxval(a%y) - minval(a%y))) + 0.03_R8P)
   call cells(u, sw)
   call fig_segments(u, 'box', a, b)
   moved = section ; call moved%move(q1)
   call curve_segments(moved, LIFT, a, b)
   call fig_segments(u, 'accent', a, b)
   k = maxloc(abs(sw%residual))                                    ! the cell of the knot
   call fig_box(u, 'probe', vec(x(k(1)), x(k(2)), -LIFT), vec(x(k(1) + 1), x(k(2) + 1), LIFT))
   write(text, '(A,ES8.1)') 'knot: residual ', maxval(abs(sw%residual))
   call fig_label(u, 'probe', vec(x(k(1)), x(k(2) + 1) + 0.01_R8P, LIFT), trim(text))
   call fig_close(u)
   endsubroutine figure

   subroutine cells(u, sk)
   !< The cells with swept area, shaded by it (per cell area).
   integer(I4P),           intent(in) :: u  !< Unit of the scene.
   type(swept_cells_2d_t), intent(in) :: sk !< Swept geometry.
   integer(I4P)                       :: ii, jj

   do jj = 1, 48 ; do ii = 1, 48
      if (sk%swept(ii, jj) == 0._R8P) cycle
      call fig_quad(u, 'field', vec(x(ii), x(jj)), vec(x(ii + 1), x(jj)), vec(x(ii + 1), x(jj + 1)), &
                    vec(x(ii), x(jj + 1)), value=sk%swept(ii, jj) / 0.025_R8P**2)
   enddo ; enddo
   endsubroutine cells

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
endprogram tut11_3_swept
