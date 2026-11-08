!< Figure of the reference page "Cut cells: solver support", section "Edge fractions": cube-01 turned by 0.3 rad about
!< (1, 1, 1) and shifted, on a 4 x 4 x 4 block of cells of 0.4 (cut_cells_3d with edges=.true.). The z edges of the
!< block: fluid parts in teal, solid parts in ember (each cut edge split at its fluid fraction, on the side of its
!< fluid node), and the grid nodes inside the solid (node_in) as gold dots.
program fig_ref_cut_solver_edges
use fossil,     only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: N = 4
real(R8P),    parameter       :: H = 0.4_R8P
type(surface_stl_object)      :: cube
type(cut_cells_3d_t)          :: cut
real(R8P)                     :: g(N + 1), f, zc, dmax
type(vector_R8P), allocatable :: fa(:), fb(:), sa(:), sb(:), nodes(:)
integer(I4P)                  :: i, j, k, mf, ms, ncut, nn, w

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
g = [(-0.12_R8P + H * i, i = 0, N)]
call cut_cells_3d(surface=cube, x=g, y=g, z=g, cut=cut, edges=.true.)
ncut = count(cut%edge_z > 0._R8P .and. cut%edge_z < 1._R8P)
print '(A,I0,A,I0,A,I0,A,I0)', 'z edges ', size(cut%edge_z), ': fluid ', count(cut%edge_z >= 1._R8P), ', cut ', ncut, &
                               ', solid ', count(cut%edge_z <= 0._R8P)
print '(A,I0,A,I0)', 'nodes inside the solid ', count(cut%node_in), ' of ', size(cut%node_in)

allocate(fa(2 * size(cut%edge_z)), fb(2 * size(cut%edge_z)), sa(2 * size(cut%edge_z)), sb(2 * size(cut%edge_z)))
mf = 0 ; ms = 0 ; dmax = 0._R8P
do k = 1, N
   do j = 1, N + 1
      do i = 1, N + 1
         f = cut%edge_z(i, j, k)
         if (f >= 1._R8P) then
            call add(fa, fb, mf, g(k), g(k + 1))
         elseif (f <= 0._R8P) then
            call add(sa, sb, ms, g(k), g(k + 1))
         elseif (.not. cut%node_in(i, j, k)) then             ! fluid at the low end
            zc = g(k) + f * H
            call add(fa, fb, mf, g(k), zc) ; call add(sa, sb, ms, zc, g(k + 1)) ; call on_wall(zc)
         else                                                ! fluid at the high end
            zc = g(k + 1) - f * H
            call add(sa, sb, ms, g(k), zc) ; call add(fa, fb, mf, zc, g(k + 1)) ; call on_wall(zc)
         endif
      enddo
   enddo
enddo
print '(A,ES9.2)', 'cut edges split at their fraction: largest distance of the split point to the surface ', dmax
allocate(nodes(count(cut%node_in))) ; nn = 0
do k = 1, N + 1
   do j = 1, N + 1
      do i = 1, N + 1
         if (.not. cut%node_in(i, j, k)) cycle
         nn = nn + 1 ; nodes(nn) = vector_R8P(g(i), g(j), g(k))
      enddo
   enddo
enddo
call fig_open(w, 'ref-cut-solver-edges', 'z edges: fluid part (teal), solid part (ember); nodes inside the '// &
              'solid (gold)', -62._R8P, 18._R8P, width=640_I4P, height=560_I4P)
call fig_surface(w, 'ghost', cube)
call fig_segments(w, 'probe', fa(1:mf), fb(1:mf))
call fig_segments(w, 'hole', sa(1:ms), sb(1:ms))
call fig_points(w, 'foot', nodes)
call fig_close(w)
contains
   subroutine add(a, b, m, z0, z1)
   !< Append the piece z0..z1 of the current z edge (i, j).
   type(vector_R8P), intent(inout) :: a(:), b(:) !< Segment ends.
   integer(I4P),     intent(inout) :: m          !< Count.
   real(R8P),        intent(in)    :: z0, z1     !< Range.

   m = m + 1 ; a(m) = vector_R8P(g(i), g(j), z0) ; b(m) = vector_R8P(g(i), g(j), z1)
   endsubroutine add

   subroutine on_wall(zz)
   !< Track the distance to the surface of the split point (g(i), g(j), zz).
   real(R8P), intent(in) :: zz !< Split height.

   dmax = max(dmax, cube%distance(point=vector_R8P(g(i), g(j), zz), is_signed=.false., is_square_root=.true.))
   endsubroutine on_wall
endprogram fig_ref_cut_solver_edges
