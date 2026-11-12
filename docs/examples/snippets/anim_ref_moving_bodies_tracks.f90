!< Animation of the reference page "Moving bodies", section "Piece map": a thin plate (0.9 by 0.03, tilted 60 degrees)
!< crossing a 24 x 12 grid of spacing 0.1. Per step it moves 0.16 to the right, more than a cell plus its own
!< thickness, while turning 0.04 rad. Each step is one `swept_cells_2d(..., map=)` call between consecutive poses. A cell is
!< shaded by its tracks:
!< - teal: one fresh track (fluid appears behind the plate);
!< - ember: one dead track (fluid covered ahead of it);
!< - yellow: several tracks. The plate passed through part of the cell within the step: the fluid ahead of it was pushed
!<   out (a dead track) while the fluid behind it grew (a continuing one), or the plate crossed the whole cell (dead and
!<   fresh). The end-pose overlap would merge the fluid on the two sides of the plate into one piece;
!< - grey: one continuing track, in a cell the wall touches at either end.
!< The plate is drawn at the start of the step in grey and at its end in red. The caption counts the tracks, the cells
!< with several tracks and the largest per-track residual.
program anim_ref_moving_bodies_tracks
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, piece_map_2d_t, &
                       TRACK_FRESH, TRACK_DEAD, TRACK_CONTINUING, TRACK_TRANSIENT, STATUS_OK
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter :: PI = acos(-1._R8P), H = 0.1_R8P, HL = 0.45_R8P, HT = 0.015_R8P
integer(I4P), parameter :: NX = 24, NY = 12, NS = 8
type(curve_object)      :: plate
type(rigid_pose_2d_t)   :: p0, p1
type(swept_cells_2d_t)  :: sw
type(piece_map_2d_t)    :: map
real(R8P)               :: x(NX + 1), y(NY + 1), c(2, 4)
integer(I4P)            :: i, j, k, u, t, status
character(120)          :: caption

c = reshape([-HT, -HL, HT, -HL, HT, HL, -HT, HL], [2, 4])
call plate%adopt_loops(u=c(1, :), v=c(2, :), loop_sizes=[4_I4P])
x = [(-0.3_R8P + H * i, i = 0, NX)] ; y = [(-0.6_R8P + H * j, j = 0, NY)]
call fig_open(u, 'anim-ref-moving-bodies-tracks', '', -90._R8P, 90._R8P, width=860_I4P, height=480_I4P)
call fig_segments(u, 'box', [(pt(x(i), y(1), 5.e-4_R8P), i = 1, NX + 1), (pt(x(1), y(j), 5.e-4_R8P), j = 1, NY + 1)], &
                  [(pt(x(i), y(NY + 1), 5.e-4_R8P), i = 1, NX + 1), (pt(x(NX + 1), y(j), 5.e-4_R8P), j = 1, NY + 1)])
call fig_anim(u, 0.9_R8P)
do k = 0, NS - 1
   p0 = pose(k) ; p1 = pose(k + 1)
   call swept_cells_2d(plate, p0, p1, x, y, sw, map=map, status=status)
   if (status /= STATUS_OK) error stop 'swept_cells_2d: map refused'
   write(caption, '(A,I0,A,I0,A,I0,A,I0,A,ES8.1)') 'step ', k + 1, ': ', map%nt, ' tracks, ', &
         map%multi_track_cells, ' cells with several, ', map%untrusted_cells, ' untrusted, max residual ', &
         maxval(abs(map%residual))
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   do j = 1, NY
      do i = 1, NX
         if (map%tracks(i, j) > 1_I4P) then
            call cell(i, j, 'accent')
         elseif (map%tracks(i, j) == 1_I4P) then
            t = map%first_track(i, j)
            select case (map%status(t))
            case (TRACK_FRESH, TRACK_TRANSIENT)
               call cell(i, j, 'inside')
            case (TRACK_DEAD)
               call cell(i, j, 'outside')
            case (TRACK_CONTINUING)
               if (sw%alpha0(i, j) < 1._R8P .or. sw%alpha1(i, j) < 1._R8P) call cell(i, j, 'ghost')
            endselect
         endif
      enddo
   enddo
   call outline(p0, 'ghost') ; call outline(p1, 'hole')
enddo
call fig_close(u)
contains
   pure function pose(kk) result(p)
   !< The pose after kk steps: tilted 60 degrees, moving right by 0.16 per step and turning 0.04 rad.
   integer(I4P), intent(in) :: kk
   type(rigid_pose_2d_t)    :: p

   p = pose_2d(angle=PI / 6._R8P + 0.04_R8P * kk, translation=[0.16_R8P * kk, 0._R8P])
   endfunction pose

   subroutine cell(ci, cj, style)
   !< Shade cell (ci, cj).
   integer(I4P), intent(in) :: ci, cj
   character(*), intent(in) :: style

   call fig_quad(u, style, pt(x(ci), y(cj)), pt(x(ci + 1), y(cj)), pt(x(ci + 1), y(cj + 1)), pt(x(ci), y(cj + 1)))
   endsubroutine cell

   subroutine outline(p, style)
   !< The plate at pose p, drawn over the cells.
   type(rigid_pose_2d_t), intent(in) :: p
   character(*),          intent(in) :: style
   type(vector_R8P)                  :: a(4), b(4)
   integer(I4P)                      :: m

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
endprogram anim_ref_moving_bodies_tracks
