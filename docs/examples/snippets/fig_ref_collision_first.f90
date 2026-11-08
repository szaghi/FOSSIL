!< Figure of the reference page "Continuous collision detection": first_impact_2d over five moving squares. Each square
!< at s = 0 (filled) and s = 1 (outline), its box at s = 0 grown by mu + tol/2 (thin lines, which contain the square
!< over the whole step), the path of its centre, and the first contact: the pair drawn at the time of impact (ember).
program fig_ref_collision_first
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, pose_interpolate, to_world, motion_bound_2d, &
                       first_impact_2d, CCD_CONTACT
use fossil_fig, only : fig_open, fig_segments, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter :: N = 5
real(R8P),    parameter :: PI = acos(-1._R8P), TOL = 1.e-6_R8P
real(R8P),    parameter :: US(4) = [-0.25_R8P, 0.25_R8P, 0.25_R8P, -0.25_R8P], VS(4) = [-0.25_R8P, -0.25_R8P, 0.25_R8P, 0.25_R8P]
type(curve_object)      :: sq(N)
type(rigid_pose_2d_t)   :: p0(N), p1(N)
real(R8P)               :: mu(N), lo(2, N), hi(2, N), toi, y(2)
integer(I4P)            :: pair(2), status, tested, i, j, k, overlaps, u
character(16)           :: text

do i = 1, N
   call sq(i)%adopt_loops(u=US, v=VS, loop_sizes=[4_I4P])
enddo
p0(1) = pose_2d(0._R8P, [0._R8P, 0._R8P])      ; p1(1) = pose_2d(PI / 6, [2._R8P, 0.4_R8P])
p0(2) = pose_2d(0._R8P, [3.2_R8P, 0.3_R8P])    ; p1(2) = pose_2d(-PI / 8, [2.3_R8P, 0.1_R8P])
p0(3) = pose_2d(0._R8P, [0.5_R8P, 2.2_R8P])    ; p1(3) = pose_2d(PI / 4, [2.5_R8P, 2.2_R8P])
p0(4) = pose_2d(0._R8P, [3.5_R8P, 3.0_R8P])    ; p1(4) = pose_2d(0._R8P, [1.5_R8P, 3.0_R8P])
p0(5) = pose_2d(0._R8P, [5.5_R8P, 1.5_R8P])    ; p1(5) = pose_2d(PI / 2, [5.5_R8P, 2.0_R8P])
do i = 1, N ! the broad phase of first_impact_2d, redone for the drawing: boxes at s = 0 grown by mu + tol / 2
   mu(i) = motion_bound_2d(p0(i), p1(i), US, VS)
   do k = 1, 4
      y = to_world(p0(i), [US(k), VS(k)])
      if (k == 1) then
         lo(:, i) = y ; hi(:, i) = y
      endif
      lo(:, i) = min(lo(:, i), y) ; hi(:, i) = max(hi(:, i), y)
   enddo
   lo(:, i) = lo(:, i) - mu(i) - TOL / 2 ; hi(:, i) = hi(:, i) + mu(i) + TOL / 2
enddo
overlaps = 0
do i = 1, N
   do j = i + 1, N
      if (all(lo(:, i) <= hi(:, j)) .and. all(lo(:, j) <= hi(:, i))) overlaps = overlaps + 1
   enddo
enddo
call first_impact_2d(sq, p0, p1, toi, pair, status, tolerance=TOL, tested=tested)
print '(A,5F7.3)', 'mu ', mu
print '(A,I0,A,I0,A,I0)', 'pairs ', N * (N - 1) / 2, ', boxes overlapping ', overlaps, ', time_of_impact runs ', tested
print '(A,L1,A,2I2,A,F8.5)', 'contact ', status == CCD_CONTACT, ', pair', pair, ', toi ', toi

call fig_open(u, 'ref-collision-first', 'first_impact_2d: grown boxes at s = 0 and the first contact (ember)', &
              -90._R8P, 90._R8P, height=440_I4P)
do i = 1, N
   call fig_quad(u, 'ghost', w(i, 0._R8P, 1_I4P), w(i, 0._R8P, 2_I4P), w(i, 0._R8P, 3_I4P), w(i, 0._R8P, 4_I4P))
   call fig_segments(u, 'ghost', [(w(i, 1._R8P, k), k = 1, 4)], [(w(i, 1._R8P, mod(k, 4_I4P) + 1_I4P), k = 1, 4)])
   call fig_segments(u, 'box', [vector_R8P(lo(1, i), lo(2, i), 0._R8P), vector_R8P(hi(1, i), lo(2, i), 0._R8P), &
                                vector_R8P(hi(1, i), hi(2, i), 0._R8P), vector_R8P(lo(1, i), hi(2, i), 0._R8P)], &
                               [vector_R8P(hi(1, i), lo(2, i), 0._R8P), vector_R8P(hi(1, i), hi(2, i), 0._R8P), &
                                vector_R8P(lo(1, i), hi(2, i), 0._R8P), vector_R8P(lo(1, i), lo(2, i), 0._R8P)])
   call fig_segments(u, 'probe', [(c(i, (k - 1) / 20._R8P), k = 1, 20)], [(c(i, k / 20._R8P), k = 1, 20)])
   write(text, '(I0)') i
   y = to_world(p0(i), [0._R8P, 0._R8P])
   call fig_label(u, 'accent', vector_R8P(y(1) - 0.05_R8P, y(2) - 0.05_R8P, 0._R8P), trim(text))
enddo
do k = 1, 2
   i = pair(k)
   call fig_quad(u, 'patch', w(i, toi, 1_I4P), w(i, toi, 2_I4P), w(i, toi, 3_I4P), w(i, toi, 4_I4P))
enddo
call fig_close(u)
contains
   function w(i, s, k) result(q)
   !< Vertex k of square i at s, world frame.
   integer(I4P), intent(in) :: i, k !< Square, vertex.
   real(R8P),    intent(in) :: s    !< Instant.
   type(vector_R8P)         :: q    !< Point.
   real(R8P)                :: z(2)

   z = to_world(pose_interpolate(p0(i), p1(i), s), [US(k), VS(k)])
   q = vector_R8P(z(1), z(2), 0._R8P)
   endfunction w

   function c(i, s) result(q)
   !< Centre of square i at s, world frame.
   integer(I4P), intent(in) :: i !< Square.
   real(R8P),    intent(in) :: s !< Instant.
   type(vector_R8P)         :: q !< Point.
   real(R8P)                :: z(2)

   z = to_world(pose_interpolate(p0(i), p1(i), s), [0._R8P, 0._R8P])
   q = vector_R8P(z(1), z(2), 0._R8P)
   endfunction c
endprogram fig_ref_collision_first
