!< Animation of the reference page "Immersed-boundary queries": a NACA 0012 section (chord 1, 60 points per side,
!< pivot at the quarter chord) pitching by up to 15 degrees while it plunges and crosses a fixed grid of spacing 0.05,
!< in 12 steps. At each step the curve is moved (`move`, the segment tree refitted), the grid nodes are classified by
!< `query_many`, and the ghost nodes (inside, with a fluid neighbour) are mirrored through the wall by `image_points`:
!< ghost nodes in teal, each joined to its image point (yellow); an image whose segment to the wall crosses the curve
!< would be ember. The ghost set changes every step; the images stay on the fluid side, one reflection away.
program anim_ref_ib_queries_ghosts
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, CURVE_VALID
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), H = 0.05_R8P, X0 = -0.85_R8P, Y0 = -0.55_R8P
integer(I4P), parameter       :: NP = 60, NI = 47, NJ = 23, NF = 12
type(curve_object)            :: foil, moved
type(rigid_pose_2d_t)         :: pose
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: su(2 * NP), sv(2 * NP), p(2, NI * NJ), d(NI * NJ), s
real(R8P), allocatable        :: g(:,:), img(:,:)
logical,   allocatable        :: vis(:)
logical                       :: ghost(NI * NJ)
integer(I4P)                  :: i, j, k, m, n, u, ng, total
character(100)                :: caption

call make_foil
call foil%adopt_loops(u=su, v=sv, loop_sizes=[2 * NP])
if (foil%check_validity() /= CURVE_VALID) error stop 'foil: invalid'
do j = 1, NJ
   do i = 1, NI
      p(:, i + NI * (j - 1)) = [X0 + H * (i - 1), Y0 + H * (j - 1)]
   enddo
enddo
call fig_open(u, 'anim-ref-ib-queries-ghosts', '', -90._R8P, 90._R8P, width=860_I4P, height=460_I4P)
allocate(a(2 * NI + 2 * NJ), b(2 * NI + 2 * NJ)) ; m = 0                  ! the grid, common to every frame
do i = 1, NI
   m = m + 1 ; a(m) = pt(X0 + H * (i - 1), Y0) ; b(m) = pt(X0 + H * (i - 1), Y0 + H * (NJ - 1))
enddo
do j = 1, NJ
   m = m + 1 ; a(m) = pt(X0, Y0 + H * (j - 1)) ; b(m) = pt(X0 + H * (NI - 1), Y0 + H * (j - 1))
enddo
call fig_segments(u, 'cell', a(1:m), b(1:m))
call fig_anim(u, 0.45_R8P)
total = 0
do k = 0, NF - 1
   s = real(k, R8P) / NF
   pose = pose_2d(angle=-15._R8P * PI / 180._R8P * sin(2 * PI * s), &
                  translation=[0.6_R8P * s - 0.15_R8P, 0.15_R8P * sin(2 * PI * s + PI / 2)])
   moved = foil
   call moved%move(pose)
   call moved%query_many(points=p, distances=d)
   ghost = .false.
   do j = 1, NJ
      do i = 1, NI
         n = i + NI * (j - 1)
         if (d(n) >= 0._R8P) cycle
         if (i > 1) ghost(n) = ghost(n) .or. d(n - 1) >= 0._R8P
         if (i < NI) ghost(n) = ghost(n) .or. d(n + 1) >= 0._R8P
         if (j > 1) ghost(n) = ghost(n) .or. d(n - NI) >= 0._R8P
         if (j < NJ) ghost(n) = ghost(n) .or. d(n + NI) >= 0._R8P
      enddo
   enddo
   ng = count(ghost)
   g = reshape(pack(p, spread(ghost, 1, 2)), [2, ng])
   if (allocated(img)) deallocate(img, vis)
   allocate(img(2, ng), vis(ng))
   call moved%image_points(points=g, images=img, visible=vis)
   total = total + count(.not. vis)
   write(caption, '(A,I0,A,F5.1,A,I0,A,I0,A)') 'step ', k + 1, ': nose up ', 15._R8P * sin(2 * PI * s), &
         ' deg, ', ng, ' ghost nodes, ', count(.not. vis), ' images not visible'
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   deallocate(a, b) ; allocate(a(2 * NP), b(2 * NP))
   do i = 1, 2 * NP
      a(i) = world(su(i), sv(i)) ; b(i) = world(su(mod(i, 2 * NP) + 1), sv(mod(i, 2 * NP) + 1))
   enddo
   call fig_segments(u, 'accent', a, b)
   deallocate(a, b)
   a = [(pt(g(1, i), g(2, i), 1.e-3_R8P), i = 1, ng)] ; b = [(pt(img(1, i), img(2, i), 1.e-3_R8P), i = 1, ng)]
   call fig_segments(u, 'box', a, b)
   call fig_points(u, 'probe', a)
   call fig_points(u, 'foot', pack(b, vis))
   if (any(.not. vis)) call fig_points(u, 'outside', pack(b, .not. vis))
enddo
call fig_close(u)
print '(A,I0)', 'images not visible over the 12 steps: ', total
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

   function world(x, y) result(q)
   !< A body point of the section at the current pose.
   real(R8P), intent(in) :: x, y
   type(vector_R8P)      :: q

   q = pt(pose%r(1, 1) * x + pose%r(1, 2) * y + pose%t(1), pose%r(2, 1) * x + pose%r(2, 2) * y + pose%t(2), 5.e-4_R8P)
   endfunction world

   pure function pt(x, y, z) result(q)
   !< A point of the plane, lifted by z (optional).
   real(R8P), intent(in)           :: x, y
   real(R8P), intent(in), optional :: z
   type(vector_R8P)                :: q

   q%x = x ; q%y = y ; q%z = 0._R8P
   if (present(z)) q%z = z
   endfunction pt
endprogram anim_ref_ib_queries_ghosts
