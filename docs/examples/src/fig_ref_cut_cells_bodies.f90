!run figref-cut-cells-bodies fig_ref_cut_cells_bodies
!fig ref-cut-cells-bodies
!< Figure of the reference page "Cut-cell geometry", section "Several bodies": a turned ellipse (semi-axes 0.3 and
!< 0.16, a 96-gon) and a 64-gon circle of radius 0.18, each with its own pose, 0.04 apart, cut together on 12 x 7 cells
!< of 0.1 by `cut_cells_2d_bodies`. The cells are shaded by their fluid fraction (cream to sepia); the cells whose
!< fluid piece is bounded by both bodies (two entries in `piece_bodies`) are outlined in red; the wall segments are drawn
!< in the colour of their body (`wall_body`). The run checks the covering sum, then shows a crossing pose refused and
!< two squares touching along an edge accepted.
program fig_ref_cut_cells_bodies
use fossil,     only : curve_object, cut_cells_2d_t, cut_walls_2d_t, piece_bodies_t, rigid_pose_2d_t, pose_2d, &
                       cut_cells_2d_bodies, STATUS_OK
use fossil_fig, only : fig_open, fig_segments, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NX = 12, NY = 7
real(R8P),    parameter       :: H = 0.1_R8P, PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)            :: body(2), sq
type(rigid_pose_2d_t)         :: pose(2)
type(cut_cells_2d_t)          :: cut
type(cut_walls_2d_t)          :: walls
type(piece_bodies_t)          :: pb
real(R8P)                     :: x(NX + 1), y(NY + 1), u(96), v(96), solid
type(vector_R8P), allocatable :: s(:), e(:)
integer(I4P)                  :: i, j, k, m, p, status, nboth, w

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
do k = 1, 96
   u(k) = 0.3_R8P * cos(2 * PI * (k - 1) / 96) ; v(k) = 0.16_R8P * sin(2 * PI * (k - 1) / 96)
enddo
call body(1)%adopt_loops(u=u, v=v, loop_sizes=[96_I4P])
do k = 1, 64
   u(k) = 0.18_R8P * cos(2 * PI * (k - 1) / 64) ; v(k) = 0.18_R8P * sin(2 * PI * (k - 1) / 64)
enddo
call body(2)%adopt_loops(u=u(1:64), v=v(1:64), loop_sizes=[64_I4P])
pose(1) = pose_2d(angle=0.35_R8P, translation=[0.42_R8P, 0.36_R8P])   ! each body in its own pose
pose(2) = pose_2d(angle=0._R8P, translation=[0.93_R8P, 0.33_R8P])
call cut_cells_2d_bodies(body, pose, x, y, cut, status=status, piece_bodies=pb, walls=walls)
if (status /= STATUS_OK) error stop 'cut_cells_2d_bodies failed'
nboth = 0 ; solid = 0._R8P
do p = 1, pb%np
   if (pb%start(p + 1) - pb%start(p) == 2) nboth = nboth + 1
enddo
do j = 1, NY
   do i = 1, NX
      solid = solid + (1._R8P - cut%alpha(i, j)) * H**2
   enddo
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'cells ', NX * NY, ', pieces ', cut%np, ', wall segments ', walls%n, &
                               ', pieces bounded by both bodies ', nboth
print '(A,F10.7,A,F10.7)', 'solid area ', solid, ', A1 + A2 ', body(1)%get_area() + body(2)%get_area()
pose(2) = pose_2d(angle=0._R8P, translation=[0.8_R8P, 0.33_R8P])      ! now the circle crosses the ellipse
call cut_cells_2d_bodies(body, pose, x, y, cut, status=status)
print '(A,I0,A)', 'crossing pose: status ', status, ' (refused)'
call sq%adopt_loops(u=[0._R8P, 0.25_R8P, 0.25_R8P, 0._R8P], v=[0._R8P, 0._R8P, 0.25_R8P, 0.25_R8P], loop_sizes=[4_I4P])
pose(1) = pose_2d(angle=0._R8P, translation=[0.25_R8P, 0.25_R8P])
pose(2) = pose_2d(angle=0._R8P, translation=[0.5_R8P, 0.375_R8P])     ! along the edge x = 0.5
call cut_cells_2d_bodies([sq, sq], pose, x, y, cut, status=status)
print '(A,I0,A,I0)', 'touching squares: status ', status, ', untrusted cells ', count(cut%pieces < 0)

! the figure, on the first cut
pose(1) = pose_2d(angle=0.35_R8P, translation=[0.42_R8P, 0.36_R8P])
pose(2) = pose_2d(angle=0._R8P, translation=[0.93_R8P, 0.33_R8P])
call cut_cells_2d_bodies(body, pose, x, y, cut, status=status, piece_bodies=pb, walls=walls)
call fig_open(w, 'ref-cut-cells-bodies', 'fluid fraction (shading), cells of a piece bounded by both bodies (red), '// &
              'walls by body (gold, teal)', -90._R8P, 90._R8P, width=760_I4P, height=470_I4P)
do j = 1, NY
   do i = 1, NX
      call fig_quad(w, 'frac', q(x(i), y(j), 0._R8P), q(x(i + 1), y(j), 0._R8P), q(x(i + 1), y(j + 1), 0._R8P), &
                    q(x(i), y(j + 1), 0._R8P), value=max(0._R8P, cut%alpha(i, j)))
   enddo
enddo
call fig_segments(w, 'box', [([(q(x(i), y(j), LIFT), i = 1, NX)], j = 1, NY + 1)], &
                  [([(q(x(i + 1), y(j), LIFT), i = 1, NX)], j = 1, NY + 1)])
call fig_segments(w, 'box', [([(q(x(i), y(j), LIFT), j = 1, NY)], i = 1, NX + 1)], &
                  [([(q(x(i), y(j + 1), LIFT), j = 1, NY)], i = 1, NX + 1)])
allocate(s(4 * pb%np), e(4 * pb%np)) ; m = 0
do p = 1, pb%np                                                ! the cells of the pieces bounded by both bodies
   if (pb%start(p + 1) - pb%start(p) < 2) cycle
   i = cut%piece_cell(1, p) ; j = cut%piece_cell(2, p)
   s(m + 1:m + 4) = [q(x(i), y(j), 2 * LIFT), q(x(i + 1), y(j), 2 * LIFT), q(x(i + 1), y(j + 1), 2 * LIFT), &
                     q(x(i), y(j + 1), 2 * LIFT)]
   e(m + 1:m + 4) = [q(x(i + 1), y(j), 2 * LIFT), q(x(i + 1), y(j + 1), 2 * LIFT), q(x(i), y(j + 1), 2 * LIFT), &
                     q(x(i), y(j), 2 * LIFT)]
   m = m + 4
enddo
call fig_segments(w, 'hole', s(1:m), e(1:m))
deallocate(s, e) ; allocate(s(walls%n), e(walls%n))
do k = 1, 2                                                    ! the wall segments, by body
   m = 0
   do i = 1, walls%n
      if (pb%wall_body(i) /= k) cycle
      m = m + 1
      s(m) = q(walls%ends(1, i), walls%ends(2, i), 3 * LIFT) ; e(m) = q(walls%ends(3, i), walls%ends(4, i), 3 * LIFT)
   enddo
   call fig_segments(w, trim(merge('foot ', 'probe', k == 1)), s(1:m), e(1:m))
enddo
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_cut_cells_bodies
