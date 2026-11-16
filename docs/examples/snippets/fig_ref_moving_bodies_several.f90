!< Figure of the reference page "Moving bodies", section "Several bodies": a three-blade rotor (tip radius 0.55, a
!< 96-gon) turning 0.15 rad inside a fixed annular casing (radii 0.6 and 0.75), a tip gap of half a cell, on 16 x 16
!< cells of 0.1, swept by `swept_cells_2d_bodies` on the Gauss path. The cells are shaded by their swept area over the
!< cell area (ember: fluid gained, teal: fluid lost); the rotor at the start (dark) and at the end (yellow), the casing
!< (grey); the cells whose fluid piece is bounded by both bodies at the start or at the end are outlined in red. The run
!< prints the balance, the cells shared by the two bodies, then the refusal of the exact path.
program fig_ref_moving_bodies_several
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d_bodies, cut_cells_2d_t, &
                       cut_cells_2d_bodies, piece_bodies_t, STATUS_OK
use fossil_fig, only : fig_open, fig_segments, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: N = 16
real(R8P),    parameter       :: H = 0.1_R8P, PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)            :: body(2), moved
type(rigid_pose_2d_t)         :: pose0(2), pose1(2)
type(swept_cells_2d_t)        :: sw
type(cut_cells_2d_t)          :: cut
type(piece_bodies_t)          :: pb
real(R8P)                     :: x(N + 1), u(192), v(192), th, r
type(vector_R8P), allocatable :: s(:), e(:)
logical                       :: both(N, N)
integer(I4P)                  :: i, j, k, m, p, status, w, nb0, nb1

x = [(-0.8_R8P + H * i, i = 0, N)]
do k = 1, 96                                                  ! the rotor: r = 0.55 (0.55 + 0.45 cos 3 theta)
   th = 2 * PI * (k - 1) / 96 ; r = 0.55_R8P * (0.55_R8P + 0.45_R8P * cos(3 * th))
   u(k) = r * cos(th) ; v(k) = r * sin(th)
enddo
call body(1)%adopt_loops(u=u(1:96), v=v(1:96), loop_sizes=[96_I4P])
do k = 1, 96                                                  ! the casing: outer loop counter-clockwise, inner clockwise
   u(k) = 0.75_R8P * cos(2 * PI * (k - 1) / 96) ; v(k) = 0.75_R8P * sin(2 * PI * (k - 1) / 96)
   u(96 + k) = 0.6_R8P * cos(-2 * PI * (k - 1) / 96) ; v(96 + k) = 0.6_R8P * sin(-2 * PI * (k - 1) / 96)
enddo
call body(2)%adopt_loops(u=u, v=v, loop_sizes=[96_I4P, 96_I4P])
pose0(1) = pose_2d(angle=0.1_R8P) ; pose1(1) = pose_2d(angle=0.25_R8P)   ! the rotor turns, each body on its own motion
pose0(2) = pose_2d(angle=0._R8P)  ; pose1(2) = pose0(2)                  ! the casing stays
call swept_cells_2d_bodies(body, pose0, pose1, x, x, sw, exact=.false., status=status)
if (status /= STATUS_OK) error stop 'swept_cells_2d_bodies failed'
print '(A,I0,A,I0,A,ES9.2)', 'Gauss nodes ', sw%evaluations, ', cells with a swept area ', count(sw%swept /= 0._R8P), &
                             ', sum of swept ', sum(sw%swept)
print '(A,ES9.2,A,ES9.2)', 'largest residual ', maxval(abs(sw%residual)), ', error estimate ', sw%error_estimate
both = .false.
call shared(pose0, nb0) ; call shared(pose1, nb1)
print '(A,I0,A,I0,A,I0)', 'pieces bounded by both bodies: start ', nb0, ', end ', nb1, ', cells ', count(both)
call swept_cells_2d_bodies(body, pose0, pose1, x, x, sw, status=status)
print '(A,I0,A)', 'exact path asked: status ', status, ' (refused: several bodies have the Gauss path only)'

! the figure
call swept_cells_2d_bodies(body, pose0, pose1, x, x, sw, exact=.false., status=status)
call fig_open(w, 'ref-moving-bodies-several', 'swept area / cell area; rotor start dark, end yellow; '// &
              'shared cells red', -90._R8P, 90._R8P, width=640_I4P, height=640_I4P)
do j = 1, N
   do i = 1, N
      call fig_quad(w, 'field', q(x(i), x(j), 0._R8P), q(x(i + 1), x(j), 0._R8P), q(x(i + 1), x(j + 1), 0._R8P), &
                    q(x(i), x(j + 1), 0._R8P), value=sw%swept(i, j) / H**2)
   enddo
enddo
allocate(s(4 * N * N), e(4 * N * N)) ; m = 0
do j = 1, N
   do i = 1, N
      if (.not. both(i, j)) cycle
      s(m + 1:m + 4) = [q(x(i), x(j), 2 * LIFT), q(x(i + 1), x(j), 2 * LIFT), q(x(i + 1), x(j + 1), 2 * LIFT), &
                        q(x(i), x(j + 1), 2 * LIFT)]
      e(m + 1:m + 4) = [q(x(i + 1), x(j), 2 * LIFT), q(x(i + 1), x(j + 1), 2 * LIFT), q(x(i), x(j + 1), 2 * LIFT), &
                        q(x(i), x(j), 2 * LIFT)]
      m = m + 4
   enddo
enddo
call fig_segments(w, 'hole', s(1:m), e(1:m))
moved = body(1) ; call moved%move(pose0(1)) ; call outline(moved, 'cell')
moved = body(1) ; call moved%move(pose1(1)) ; call outline(moved, 'accent')
call outline(body(2), 'box')
call fig_close(w)
contains
   subroutine shared(pose, nboth)
   !< The pieces bounded by both bodies at a pose set, their cells marked in `both`.
   type(rigid_pose_2d_t), intent(in)  :: pose(2) !< Poses.
   integer(I4P),          intent(out) :: nboth   !< Pieces bounded by both bodies.
   integer(I4P)                       :: p_      !< Piece.

   call cut_cells_2d_bodies(body, pose, x, x, cut, status=status, piece_bodies=pb)
   nboth = 0
   do p_ = 1, pb%np
      if (pb%start(p_ + 1) - pb%start(p_) < 2) cycle
      nboth = nboth + 1 ; both(cut%piece_cell(1, p_), cut%piece_cell(2, p_)) = .true.
   enddo
   endsubroutine shared

   subroutine outline(c, style)
   !< The segments of curve c, just above the cells.
   type(curve_object), intent(in) :: c     !< Curve.
   character(*),       intent(in) :: style !< Style.
   type(vector_R8P), allocatable  :: a(:), b(:)
   real(R8P)                      :: uv(2)
   integer(I4P)                   :: i_, n_, nv

   nv = c%get_vertices_number()
   allocate(a(nv), b(nv)) ; n_ = 0
   do i_ = 1, nv
      if (c%get_next(i_) == 0_I4P) cycle
      n_ = n_ + 1
      uv = c%get_vertex(i_) ; a(n_) = q(uv(1), uv(2), 3 * LIFT)
      uv = c%get_vertex(c%get_next(i_)) ; b(n_) = q(uv(1), uv(2), 3 * LIFT)
   enddo
   call fig_segments(w, style, a(1:n_), b(1:n_))
   endsubroutine outline

   function q(a_, b_, dz) result(r_)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r_         !< Point.

   r_ = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_moving_bodies_several
