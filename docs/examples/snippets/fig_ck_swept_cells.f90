!< Figure of the cookbook recipe "Swept volume of a moving body on a grid": the plate of ck_swept_cells at its start pose
!< (blue outline) and end pose (yellow outline), turned 30 degrees and moved in one step, over the 20 x 20 cells
!< coloured by the swept area over the cell area (ember: fluid gained, teal: fluid lost), seen from above and zoomed.
program fig_ck_swept_cells
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curve_object)     :: plate
type(swept_cells_2d_t) :: sw
type(rigid_pose_2d_t)  :: p0, p1
type(vector_R8P)       :: a(21), b(21)
real(R8P)              :: x(21), cu(4), cv(4)
integer(I4P)           :: i, j, un

cu = [-0.15_R8P, 0.15_R8P, 0.15_R8P, -0.15_R8P] ; cv = [-0.05_R8P, -0.05_R8P, 0.05_R8P, 0.05_R8P]
call plate%adopt_loops(u=cu, v=cv, loop_sizes=[4_I4P])
x = [(0.05_R8P * i, i = 0, 20)]
p0 = pose_2d(0._R8P, [0.4_R8P, 0.45_R8P]) ; p1 = pose_2d(acos(-1._R8P) / 6, [0.55_R8P, 0.5_R8P])
call swept_cells_2d(plate, p0, p1, x, x, sw)
call fig_open(un, 'ck-swept-cells', 'swept area: ember gained, teal lost; plate at s = 0 (blue), s = 1 (yellow)', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.5_R8P, 0.5_R8P, 0._R8P), 0.3_R8P)
do j = 1, 20 ; do i = 1, 20
   call fig_quad(un, 'field', vector_R8P(x(i), x(j), 0._R8P), vector_R8P(x(i + 1), x(j), 0._R8P), &
                 vector_R8P(x(i + 1), x(j + 1), 0._R8P), vector_R8P(x(i), x(j + 1), 0._R8P), value=sw%swept(i, j) / 0.05_R8P**2)
enddo ; enddo
a = [(vector_R8P(x(i), 0._R8P, 1.e-3_R8P), i = 1, 21)] ; b = [(vector_R8P(x(i), 1._R8P, 1.e-3_R8P), i = 1, 21)]
call fig_segments(un, 'box', a, b)
a = [(vector_R8P(0._R8P, x(i), 1.e-3_R8P), i = 1, 21)] ; b = [(vector_R8P(1._R8P, x(i), 1.e-3_R8P), i = 1, 21)]
call fig_segments(un, 'box', a, b)
call outline(p0, 'probe') ; call outline(p1, 'accent')
call fig_close(un)
contains
   subroutine outline(p, style)
   !< The plate's four sides at pose p.
   type(rigid_pose_2d_t), intent(in) :: p     !< Pose.
   character(*),          intent(in) :: style !< Style.
   type(vector_R8P)                  :: c(4)  !< Corners.
   integer(I4P)                      :: k     !< Counter.

   do k = 1, 4
      c(k) = vector_R8P(p%r(1, 1) * cu(k) + p%r(1, 2) * cv(k) + p%t(1), p%r(2, 1) * cu(k) + p%r(2, 2) * cv(k) + p%t(2), &
                        2.e-3_R8P)
   enddo
   call fig_segments(un, style, c, [c(2), c(3), c(4), c(1)])
   endsubroutine outline
endprogram fig_ck_swept_cells
