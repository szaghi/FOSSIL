!< Animation of the curvilinear-blocks reference page: three points moving through an O-grid around a cylinder (40 x 10
!< cells, radii 1 to 3, the seam at theta = 0 found by detect_topology_2d), located frame after frame by stencil
!< walking: locate_many(walk=.true., start=) starts each point from the cell it was found in at the previous frame.
!< The start cells are teal, the cells found ember, the points with their paths in yellow; the caption gives the cells
!< crossed by the three walks (jumps) and the points handed to the robust path (fallbacks).
program anim_block_walk
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, LOC_INSIDE
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_segments, fig_quad, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter      :: NI = 41, NJ = 11, NP = 3, NF = 16
real(R8P),    parameter      :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curvilinear_block_2d_t) :: grid(1)
type(vector_R8P)             :: a(2 * NI * NJ), b(2 * NI * NJ), path(NP, 0:NF)
real(R8P)                    :: u(NI, NJ), v(NI, NJ), r, px(NP), py(NP)
real(R8P),    allocatable    :: xi(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:)
integer(I4P)                 :: start(2, NP), i, j, k, f, m, un, jumps, fallbacks
character(80)                :: caption

do j = 1, NJ                                                   ! the O-grid
   r = 1._R8P + 2._R8P * (j - 1) / (NJ - 1)
   do i = 1, NI
      u(i, j) = r * cos(2 * PI * (i - 1) / (NI - 1)) ; v(i, j) = r * sin(2 * PI * (i - 1) / (NI - 1))
   enddo
enddo
u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)                            ! the seam, written exactly
call grid(1)%initialize(u, v)
call detect_topology_2d(grid)
call grid(1)%build_search

call fig_open(un, 'anim-block-walk', '', -90._R8P, 90._R8P, width=560_I4P, height=520_I4P)
call fig_anim(un, 0.4_R8P)
m = 0                                                          ! the grid, in every frame
do j = 1, NJ
   do i = 1, NI - 1
      m = m + 1 ; a(m) = pt(u(i, j), v(i, j), 0._R8P) ; b(m) = pt(u(i + 1, j), v(i + 1, j), 0._R8P)
   enddo
enddo
do i = 1, NI - 1                                               ! radial lines: straight
   m = m + 1 ; a(m) = pt(u(i, 1), v(i, 1), 0._R8P) ; b(m) = pt(u(i, NJ), v(i, NJ), 0._R8P)
enddo
call fig_segments(un, 'box', a(1:m), b(1:m))

call position(0_I4P, px, py)                                   ! frame 0: the robust path
call grid(1)%locate_many(px, py, cell, xi, st)
start = cell
do k = 1, NP
   path(k, 0) = pt(px(k), py(k), 3 * LIFT)
enddo
do f = 1, NF
   call position(f, px, py)
   call grid(1)%locate_many(px, py, cell, xi, st, walk=.true., start=start, jumps=jumps, fallbacks=fallbacks)
   write(caption, '(A,I2,A,I3,A,I0,A,I0)') 'frame ', f, ': jumps ', jumps, ', fallbacks ', fallbacks, &
                                            ', inside ', count(st == LOC_INSIDE)
   print '(A,A,6I4)', trim(caption), ', cells', cell
   call fig_frame(un, caption=trim(caption))
   do k = 1, NP
      path(k, f) = pt(px(k), py(k), 3 * LIFT)
      call fill(start(:, k), -0.8_R8P, LIFT)
      call fill(cell(:, k), 0.9_R8P, 2 * LIFT)
      call fig_segments(un, 'accent', path(k, 0:f - 1), path(k, 1:f))
   enddo
   call fig_points(un, 'foot', path(:, f))
   start = cell
enddo
call fig_close(un)
contains
   subroutine position(f, x, y)
   !< The points at frame f: one on a circle across the seam, one on a spiral, one on a straight line.
   integer(I4P), intent(in)  :: f        !< Frame.
   real(R8P),    intent(out) :: x(NP), y(NP) !< Points.
   real(R8P)                 :: t, rr

   t = real(f, R8P) / NF
   rr = 1.5_R8P ; x(1) = rr * cos(-1.6_R8P + 4.2_R8P * t) ; y(1) = rr * sin(-1.6_R8P + 4.2_R8P * t)
   rr = 2.85_R8P - 1.7_R8P * t ; x(2) = rr * cos(1.6_R8P + 3.6_R8P * t) ; y(2) = rr * sin(1.6_R8P + 3.6_R8P * t)
   x(3) = 2.7_R8P - 4.6_R8P * t ; y(3) = 1.1_R8P + 0.5_R8P * t
   endsubroutine position

   subroutine fill(c, value, dz)
   !< Fill cell c = (i, j) with a colour of the field style.
   integer(I4P), intent(in) :: c(2)  !< Cell.
   real(R8P),    intent(in) :: value !< Field value (-1 teal, 1 ember).
   real(R8P),    intent(in) :: dz    !< Lift.

   if (any(c < 1)) return
   call fig_quad(un, 'field', pt(u(c(1), c(2)), v(c(1), c(2)), dz), pt(u(c(1) + 1, c(2)), v(c(1) + 1, c(2)), dz), &
                 pt(u(c(1) + 1, c(2) + 1), v(c(1) + 1, c(2) + 1), dz), pt(u(c(1), c(2) + 1), v(c(1), c(2) + 1), dz), value)
   endsubroutine fill

   function pt(x, y, dz) result(q)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: x, y, dz !< Coordinates, lift.
   type(vector_R8P)      :: q        !< Point.

   q = vector_R8P(x, y, dz)
   endfunction pt
endprogram anim_block_walk
