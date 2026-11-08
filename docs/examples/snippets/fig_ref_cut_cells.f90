!< Figures of the reference page "Cut-cell geometry":
!< - 2-D: a diamond whose vertices are grid nodes and whose sides pass through grid nodes, on an 8 x 8 block; the
!<   cells shaded by their fluid fraction, the crossings of the v lines from `line_crossings` (entering in teal,
!<   leaving in ember). Every crossing lies on a node, and each is counted once per line by the symbolic tie-break;
!< - 3-D: a thin tilted plate (thickness 0.3 of a cell) across a 4 x 4 x 4 block; the cells it splits into two fluid
!<   pieces outlined, and the two piece centroids of each such cell, one on either side of the plate.
program fig_ref_cut_cells
use fossil,     only : curve_object, cut_cells_2d_t, surface_stl_object, cut_cells_3d_t, cut_cells_3d, rigid_pose_t, &
                       pose_from_axis_angle
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_box, fig_label, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter          :: LIFT = 1.e-3_R8P
type(curve_object)            :: diamond
type(cut_cells_2d_t)          :: cut2
type(surface_stl_object)      :: plate
type(cut_cells_3d_t)          :: cut3
type(rigid_pose_t)            :: pose
type(vector_R8P), allocatable :: pin(:), pout(:), a(:), b(:)
real(R8P),        allocatable :: t(:)
integer(I4P),     allocatable :: seg(:)
logical,          allocatable :: ent(:)
real(R8P)                     :: e(9), g(5)
integer(I4P)                  :: u, i, j, k, n, p, two

! ---- 2-D: the diamond (vertices at grid nodes, sides through nodes), its cut cells and the grid-line crossings
e = [(0.125_R8P * i, i = 0, 8)]
call diamond%adopt_loops(u=[0.5_R8P, 0.875_R8P, 0.5_R8P, 0.125_R8P], v=[0.125_R8P, 0.5_R8P, 0.875_R8P, 0.5_R8P], &
                         loop_sizes=[4_I4P])
call diamond%cut_cells(x=e, y=e, cut=cut2)
allocate(pin(0), pout(0))
print '(A,F8.5,A,I0)', 'diamond: solid area ', sum(1._R8P - cut2%alpha) / 64._R8P, ', cut cells ', &
                      count(cut2%alpha > 0._R8P .and. cut2%alpha < 1._R8P)
do i = 1, 9                                                                           ! the horizontal grid lines
   call diamond%line_crossings(axis=1_I4P, a=e(i), t0=0._R8P, t1=1._R8P, t=t, segments=seg, entering=ent, n=n)
   print '(A,F6.3,A,I0,A,*(F6.3,:,","))', '  line v =', e(i), ': ', n, ' crossings at u =', t(1:n)
   do k = 1, n
      if (ent(k)) then ; pin = [pin, vec(t(k), e(i), 2 * LIFT)] ; else ; pout = [pout, vec(t(k), e(i), 2 * LIFT)] ; endif
   enddo
enddo
call fig_open(u, 'ref-cut-cells-crossings', 'diamond on grid nodes: alpha, v-line crossings (teal enter, ember leave)', &
              -90._R8P, 90._R8P)
do j = 1, 8 ; do i = 1, 8
   call fig_quad(u, 'frac', vec(e(i), e(j)), vec(e(i + 1), e(j)), vec(e(i + 1), e(j + 1)), vec(e(i), e(j + 1)), &
                 value=cut2%alpha(i, j))
enddo ; enddo
allocate(a(18), b(18))
do i = 1, 9                                                                           ! the grid lines
   a(i) = vec(e(i), 0._R8P, 0.5_R8P * LIFT) ; b(i) = vec(e(i), 1._R8P, 0.5_R8P * LIFT)
   a(9 + i) = vec(0._R8P, e(i), 0.5_R8P * LIFT) ; b(9 + i) = vec(1._R8P, e(i), 0.5_R8P * LIFT)
enddo
call fig_segments(u, 'box', a, b)
a(1:4) = [vec(0.5_R8P, 0.125_R8P, LIFT), vec(0.875_R8P, 0.5_R8P, LIFT), vec(0.5_R8P, 0.875_R8P, LIFT), &
          vec(0.125_R8P, 0.5_R8P, LIFT)]
b(1:4) = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a(1:4), b(1:4))                                         ! the wall
call fig_points(u, 'inside', pin)
call fig_points(u, 'outside', pout)
call fig_label(u, 'accent', vec(0.53_R8P, 0.09_R8P, 3._R8P * LIFT), 'v = 1/8 taken at v + eps: 2 crossings')
call fig_label(u, 'accent', vec(0.53_R8P, 0.90_R8P, 3._R8P * LIFT), 'v = 7/8 taken at v + eps: none')
call fig_close(u)

! ---- 3-D: a thin plate (0.3 of a cell thick) tilted across a 4 x 4 x 4 block, the cells it splits in two pieces
call plate%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call plate%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)                              ! centred at the origin
call plate%resize(x=1.6_R8P, y=1.6_R8P, z=0.075_R8P)                                  ! 1.6 x 1.6 x 0.075
call plate%analyze
pose = pose_from_axis_angle(axis=[1._R8P, -1._R8P, 0._R8P], angle=0.45_R8P, translation=[0.5_R8P, 0.5_R8P, 0.5_R8P])
call plate%move(pose)
g = [(0.25_R8P * i, i = 0, 4)]
call cut_cells_3d(surface=plate, x=g, y=g, z=g, cut=cut3)
two = count(cut3%pieces == 2)
print '(A,I0,A,I0,A,I0)', 'plate: cells cut ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1._R8P), &
                          ', with 2 pieces ', two, ', flagged ', count(cut3%pieces < 0)
call fig_open(u, 'ref-cut-cells-pieces', 'a plate thinner than the cell: cells split in 2 fluid pieces, their centroids', &
              -45._R8P, 8._R8P)
call fig_surface(u, 'body', plate)
deallocate(pin, pout) ; allocate(pin(0), pout(0))
do k = 1, 4 ; do j = 1, 4 ; do i = 1, 4
   if (cut3%pieces(i, j, k) /= 2) cycle
   call fig_box(u, 'accent', vec(g(i), g(j), g(k)), vec(g(i + 1), g(j + 1), g(k + 1)))
   p = cut3%first_piece(i, j, k)
   pin = [pin, vec(cut3%piece_centroid(1, p), cut3%piece_centroid(2, p), cut3%piece_centroid(3, p))]
   pout = [pout, vec(cut3%piece_centroid(1, p + 1), cut3%piece_centroid(2, p + 1), cut3%piece_centroid(3, p + 1))]
enddo ; enddo ; enddo
call fig_box(u, 'box', vec(0._R8P, 0._R8P, 0._R8P), vec(1._R8P, 1._R8P, 1._R8P))
call fig_points(u, 'probe', pin)
call fig_points(u, 'outside', pout)
call fig_close(u)
contains
   pure function vec(x, y, z) result(v)
   !< A vector, z = 0 by default.
   real(R8P), intent(in)           :: x, y
   real(R8P), intent(in), optional :: z
   type(vector_R8P)                :: v

   v%x = x ; v%y = y ; v%z = 0._R8P
   if (present(z)) v%z = z
   endfunction vec
endprogram fig_ref_cut_cells
