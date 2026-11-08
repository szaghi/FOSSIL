!run anim-ccd anim_ccd
!fig anim-ccd
!< Animation of the collision reference page: conservative advancement of a square sliding along a floor into a wall,
!< plain bound (top) against the directional bound (bottom). The square (side 0.3) starts 0.18 above the floor of an
!< L-shaped body, translates by 4 along x and turns 45 degrees over the step; it first touches the wall. Frame k calls
!< curve%time_of_impact with max_iterations = k, whose toi is the k-th iterate s_k: the square is drawn at s_k
!< (ember) over the earlier iterates (grey), with the closest pair at s_k (yellow). The plain bound advances by
!< (d - tol/2) / mu, d the floor gap; the directional step sees that the sliding does not close the floor gap.
program anim_ccd
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, pose_interpolate, to_world, CCD_CONTACT
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_segments, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter :: NB = 6, NA = 4, KMAX = 60
real(R8P),    parameter :: PI = acos(-1._R8P)
real(R8P),    parameter :: UB(NB) = [-0.5_R8P, 3.6_R8P, 3.6_R8P, 3.3_R8P, 3.3_R8P, -0.5_R8P], &
                           VB(NB) = [-0.2_R8P, -0.2_R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P],       &
                           UA(NA) = [-0.15_R8P, 0.15_R8P, 0.15_R8P, -0.15_R8P],                   &
                           VA(NA) = [-0.15_R8P, -0.15_R8P, 0.15_R8P, 0.15_R8P],                   &
                           LIFT(2) = [1.7_R8P, 0._R8P]                                          ! panel offsets
type(curve_object)      :: a, b
type(rigid_pose_2d_t)   :: a0, a1, rest
real(R8P)               :: s(0:KMAX, 2), pa(2, 0:KMAX, 2), pb(2, 0:KMAX, 2), toi, d
integer(I4P)            :: n(2), k, m, run, u, status, it
logical                 :: dir
character(100)          :: caption

call b%adopt_loops(u=UB, v=VB, loop_sizes=[NB])                ! floor and wall, at rest
call a%adopt_loops(u=UA, v=VA, loop_sizes=[NA])                ! the square, about its body origin
rest = pose_2d(angle=0._R8P)
a0 = pose_2d(angle=0._R8P, translation=[0._R8P, 0.33_R8P])
a1 = pose_2d(angle=PI / 4, translation=[4._R8P, 0.33_R8P])

do run = 1, 2                                                  ! 1: plain bound, 2: directional
   dir = run == 2
   s(0, run) = 0._R8P ; n(run) = KMAX
   call a%distance_to(b, d, point_self=pa(:, 0, run), point_other=pb(:, 0, run), pose_self=a0)
   do k = 1, KMAX                                              ! the k-th iterate: the step capped at k
      call a%time_of_impact(b, a0, a1, rest, rest, toi, status, point_self=pa(:, k, run), &
                            point_other=pb(:, k, run), iterations=it, max_iterations=k, directional=dir)
      s(k, run) = toi
      if (status == CCD_CONTACT) then
         n(run) = it ; exit
      endif
   enddo
   print '(A,A,I0,A,F7.4)', trim(merge('directional', 'plain      ', dir)), ': contact after ', n(run), &
                            ' steps at s = ', s(n(run), run)
enddo

call fig_open(u, 'anim-ccd', '', -90._R8P, 90._R8P, width=720_I4P, height=420_I4P)
call fig_anim(u, 0.3_R8P)
do run = 1, 2                                                  ! the bodies at rest and the labels, in every frame
   call fig_segments(u, 'body', [(pt(UB(m), VB(m), run), m = 1, NB)], [(pt(UB(mod(m, NB) + 1), VB(mod(m, NB) + 1), &
                     run), m = 1, NB)])
   call fig_label(u, 'accent', pt(-0.4_R8P, 0.75_R8P, run), trim(merge('directional bound', 'plain bound      ', &
                  run == 2)))
enddo
do k = 0, maxval(n)
   write(caption, '(A,I2,A,F6.4,A,F6.4)') 'step ', k, ': plain s = ', s(min(k, n(1)), 1), ', directional s = ', &
                                          s(min(k, n(2)), 2)
   call fig_frame(u, caption=trim(caption), seconds=merge(2._R8P, 0.3_R8P, k == maxval(n)))
   do run = 1, 2
      do m = 0, min(k, n(run)) - 1                             ! the earlier iterates
         call square(m, run, 'box')
      enddo
      call square(min(k, n(run)), run, 'patch')
      call fig_segment(u, 'accent', pt(pa(1, min(k, n(run)), run), pa(2, min(k, n(run)), run), run), &
                       pt(pb(1, min(k, n(run)), run), pb(2, min(k, n(run)), run), run))
      call fig_points(u, 'foot', [pt(pb(1, min(k, n(run)), run), pb(2, min(k, n(run)), run), run)])
   enddo
enddo
call fig_close(u)
contains
   function pt(x, y, panel) result(q)
   !< A point of the plane, shifted to its panel.
   real(R8P),    intent(in) :: x, y  !< Coordinates.
   integer(I4P), intent(in) :: panel !< Panel (1 top, 2 bottom).
   type(vector_R8P)         :: q     !< Point.

   q = vector_R8P(x + LIFT(2), y + merge(LIFT(1), 0._R8P, panel == 1), 0._R8P)
   endfunction pt

   subroutine square(k, panel, style)
   !< The square at iterate k of a run.
   integer(I4P), intent(in) :: k, panel !< Iterate, run.
   character(*), intent(in) :: style    !< Style.
   type(rigid_pose_2d_t)    :: p
   real(R8P)                :: w(2, NA)
   integer(I4P)             :: i

   p = pose_interpolate(a0, a1, s(k, panel))
   do i = 1, NA
      w(:, i) = to_world(p, [UA(i), VA(i)])
   enddo
   call fig_segments(u, style, [(pt(w(1, i), w(2, i), panel), i = 1, NA)], &
                     [(pt(w(1, mod(i, NA) + 1), w(2, mod(i, NA) + 1), panel), i = 1, NA)])
   endsubroutine square
endprogram anim_ccd
