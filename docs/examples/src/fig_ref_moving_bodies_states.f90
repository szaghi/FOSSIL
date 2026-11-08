!run figref-moving-bodies-states fig_ref_moving_bodies_states
!fig ref-moving-bodies-states
!< Figure of the reference page "Moving bodies", section "What it does": `state_changes` over one step of a NACA 0012
!< section (chord 1, pivot at the quarter chord) that pitches nose up by 3 degrees and moves 0.02 forward (-x) and
!< 0.01 up. The cell centres of a grid of spacing 0.02 (outlined) are drawn in grey only where the band filter sends
!< them to the second query (|phi(b0)| <= |b0 - b1|): one query decides the others. Fresh centres (solid to fluid) in
!< teal, dead ones (fluid to solid) in ember; the section at the start in dark brown, at the end in yellow.
program fig_ref_moving_bodies_states
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, to_body, CURVE_VALID
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), H = 0.02_R8P, X0 = -0.45_R8P, Y0 = -0.21_R8P
integer(I4P), parameter       :: NP = 60, NI = 66, NJ = 22
type(curve_object)            :: foil
type(rigid_pose_2d_t)         :: p0, p1
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: su(2 * NP), sv(2 * NP), p(2, NI * NJ), b0(2), b1(2)
logical                       :: fresh(NI * NJ), dead(NI * NJ), band(NI * NJ)
integer(I4P)                  :: i, j, k, u, checked
character(120)                :: caption

call make_foil
call foil%adopt_loops(u=su, v=sv, loop_sizes=[2 * NP])
if (foil%check_validity() /= CURVE_VALID) error stop 'foil: invalid'
do j = 1, NJ
   do i = 1, NI
      p(:, i + NI * (j - 1)) = [X0 + H * (i - 1), Y0 + H * (j - 1)]
   enddo
enddo
!region states
p0 = pose_2d(angle=0._R8P)
p1 = pose_2d(angle=3._R8P * PI / 180._R8P, translation=[-0.02_R8P, 0.01_R8P])
call foil%state_changes(pose0=p0, pose1=p1, points=p, fresh=fresh, dead=dead, checked=checked)
!endregion states
do k = 1, NI * NJ                                   ! the band filter, recomputed to draw it
   b0 = to_body(p0, p(:, k)) ; b1 = to_body(p1, p(:, k))
   band(k) = abs(foil%distance(b0)) <= norm2(b1 - b0)
enddo
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'centres ', NI * NJ, ', second query ', checked, ' (band drawn: ', count(band), &
                                     '), fresh ', count(fresh), ', dead ', count(dead)
write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'state_changes: ', checked, ' of ', NI * NJ, ' centres queried twice (grey); ', &
      count(fresh), ' fresh (teal), ', count(dead), ' dead (ember)'
call fig_open(u, 'ref-moving-bodies-states', trim(caption), -90._R8P, 90._R8P, width=860_I4P, height=330_I4P)
a = [(pt(p(1, k), p(2, k)), k = 1, NI * NJ)]
call fig_points(u, 'box', pack(a, band .and. .not. (fresh .or. dead)))
call fig_segments(u, 'box', [pt(X0, Y0), pt(X0 + H * (NI - 1), Y0), pt(X0 + H * (NI - 1), Y0 + H * (NJ - 1)), &
                              pt(X0, Y0 + H * (NJ - 1))], [pt(X0 + H * (NI - 1), Y0), pt(X0 + H * (NI - 1), Y0 + H * (NJ - 1)), &
                              pt(X0, Y0 + H * (NJ - 1)), pt(X0, Y0)])
call fig_points(u, 'probe', pack(a, fresh))
call fig_points(u, 'outside', pack(a, dead))
allocate(b(2 * NP)) ; deallocate(a) ; allocate(a(2 * NP))
call outline(p0) ; call fig_segments(u, 'cell', a, b)
call outline(p1) ; call fig_segments(u, 'accent', a, b)
call fig_close(u)
contains
   subroutine make_foil
   !< The NACA 0012 section, counter-clockwise from the trailing edge over the upper side, the leading edge at x = -0.25.
   real(R8P)    :: xs(0:NP), ys(0:NP)
   integer(I4P) :: ii, nn

   do ii = 0, NP
      xs(ii) = 0.5_R8P * (1._R8P - cos(PI * ii / NP))
      ys(ii) = 0.6_R8P * (0.2969_R8P * sqrt(xs(ii)) - 0.1260_R8P * xs(ii) - 0.3516_R8P * xs(ii) ** 2 + &
                          0.2843_R8P * xs(ii) ** 3 - 0.1036_R8P * xs(ii) ** 4)
   enddo
   ys(NP) = 0._R8P
   nn = 0
   do ii = NP, 1, -1
      nn = nn + 1 ; su(nn) = xs(ii) - 0.25_R8P ; sv(nn) = ys(ii)
   enddo
   do ii = 0, NP - 1
      nn = nn + 1 ; su(nn) = xs(ii) - 0.25_R8P ; sv(nn) = -ys(ii)
   enddo
   endsubroutine make_foil

   subroutine outline(pose)
   !< The section at a pose, into a, b, lifted over the points.
   type(rigid_pose_2d_t), intent(in) :: pose
   integer(I4P)                      :: ii, i2

   do ii = 1, 2 * NP
      i2 = mod(ii, 2 * NP) + 1
      a(ii) = pt(pose%r(1, 1) * su(ii) + pose%r(1, 2) * sv(ii) + pose%t(1), &
                 pose%r(2, 1) * su(ii) + pose%r(2, 2) * sv(ii) + pose%t(2), 1.e-3_R8P)
      b(ii) = pt(pose%r(1, 1) * su(i2) + pose%r(1, 2) * sv(i2) + pose%t(1), &
                 pose%r(2, 1) * su(i2) + pose%r(2, 2) * sv(i2) + pose%t(2), 1.e-3_R8P)
   enddo
   endsubroutine outline

   pure function pt(x, y, z) result(q)
   !< A point of the plane, lifted by z (optional).
   real(R8P), intent(in)           :: x, y
   real(R8P), intent(in), optional :: z
   type(vector_R8P)                :: q

   q%x = x ; q%y = y ; q%z = 0._R8P
   if (present(z)) q%z = z
   endfunction pt
endprogram fig_ref_moving_bodies_states
