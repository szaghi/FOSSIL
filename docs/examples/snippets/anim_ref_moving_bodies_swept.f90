!< Animation of the reference page "Moving bodies": a square of side 0.5 that rolls along a 24 x 12 grid of spacing 0.1,
!< turning 0.15 rad and advancing 0.14 per step, with a vertical oscillation, over 10 steps. Each step is one
!< `swept_cells_2d` call between consecutive poses (the screw motion, exact rotation path): each cell is shaded by its
!< swept area over the cell area (ember: fluid gained behind the square, teal: fluid lost ahead of it), the square at
!< the start of the step in dark brown and at its end in yellow. The caption gives the sum of the swept areas (0 for a
!< rigid motion) and the largest uniform-flow residual of the step.
program anim_ref_moving_bodies_swept
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, STATUS_OK
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter :: PI = acos(-1._R8P), H = 0.1_R8P
integer(I4P), parameter :: NX = 24, NY = 12, NS = 10
type(curve_object)      :: square
type(rigid_pose_2d_t)   :: p0, p1
type(swept_cells_2d_t)  :: sw
real(R8P)               :: x(NX + 1), y(NY + 1), worst
integer(I4P)            :: i, j, k, u, status
character(100)          :: caption

call square%adopt_loops(u=[-0.25_R8P, 0.25_R8P, 0.25_R8P, -0.25_R8P], v=[-0.25_R8P, -0.25_R8P, 0.25_R8P, 0.25_R8P], &
                        loop_sizes=[4_I4P])
x = [(-0.45_R8P + H * i, i = 0, NX)] ; y = [(-0.6_R8P + H * j, j = 0, NY)]
call fig_open(u, 'anim-ref-moving-bodies-swept', '', -90._R8P, 90._R8P, width=860_I4P, height=480_I4P)
call fig_segments(u, 'box', [(pt(x(i), y(1), 5.e-4_R8P), i = 1, NX + 1), (pt(x(1), y(j), 5.e-4_R8P), j = 1, NY + 1)], &
                  [(pt(x(i), y(NY + 1), 5.e-4_R8P), i = 1, NX + 1), (pt(x(NX + 1), y(j), 5.e-4_R8P), j = 1, NY + 1)])
call fig_anim(u, 0.7_R8P)
worst = 0._R8P
do k = 0, NS - 1
   p0 = pose(k) ; p1 = pose(k + 1)
   call swept_cells_2d(square, p0, p1, x, y, sw, status=status)
   if (status /= STATUS_OK .or. .not. sw%exact) error stop 'swept_cells_2d: no exact result'
   worst = max(worst, maxval(abs(sw%residual)))
   write(caption, '(A,I0,A,ES9.2,A,ES9.2,A,I0,A,I0)') 'step ', k + 1, ': sum of swept ', sum(sw%swept), &
         ', max residual ', maxval(abs(sw%residual)), ', cleared ', count(sw%cleared), ', covered ', count(sw%covered)
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   do j = 1, NY
      do i = 1, NX
         call fig_quad(u, 'field', pt(x(i), y(j)), pt(x(i + 1), y(j)), pt(x(i + 1), y(j + 1)), pt(x(i), y(j + 1)), &
                       value=sw%swept(i, j) / (H * H))
      enddo
   enddo
   call outline(p0, 'cell') ; call outline(p1, 'accent')
enddo
call fig_close(u)
print '(A,ES9.2)', 'largest residual over the 10 steps ', worst
contains
   pure function pose(kk) result(p)
   !< The pose after kk steps: rolling (clockwise) to +x, bobbing up and down.
   integer(I4P), intent(in) :: kk
   type(rigid_pose_2d_t)    :: p

   p = pose_2d(angle=-0.15_R8P * kk, translation=[0.14_R8P * kk, 0.15_R8P * sin(PI * kk / 5._R8P)])
   endfunction pose

   subroutine outline(p, style)
   !< The square at pose p, drawn over the cells.
   type(rigid_pose_2d_t), intent(in) :: p
   character(*),          intent(in) :: style
   real(R8P)                         :: c(2, 4)
   type(vector_R8P)                  :: a(4), b(4)
   integer(I4P)                      :: m

   c = reshape([-0.25_R8P, -0.25_R8P, 0.25_R8P, -0.25_R8P, 0.25_R8P, 0.25_R8P, -0.25_R8P, 0.25_R8P], [2, 4])
   do m = 1, 4
      a(m) = pt(dot_product(p%r(1, :), c(:, m)) + p%t(1), dot_product(p%r(2, :), c(:, m)) + p%t(2), 1.e-3_R8P)
   enddo
   b = [a(2:4), a(1)]
   call fig_segments(u, style, a, b)
   endsubroutine outline

   pure function pt(px, py, pz) result(q)
   !< A point of the plane, lifted by pz (optional).
   real(R8P), intent(in)           :: px, py
   real(R8P), intent(in), optional :: pz
   type(vector_R8P)                :: q

   q%x = px ; q%y = py ; q%z = 0._R8P
   if (present(pz)) q%z = pz
   endfunction pt
endprogram anim_ref_moving_bodies_swept
