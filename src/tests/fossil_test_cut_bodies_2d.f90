!< FOSSIL test: the 2-D cut cells of several posed bodies (#109 WP3 step 03).
program fossil_test_cut_bodies_2d
!< FOSSIL test: the 2-D cut cells of several posed bodies (#109 WP3 step 03).
!<
!< `cut_cells_2d_bodies` / `merge_posed_bodies` over curves (`fossil_cut_bodies`):
!< 1. two posed polygonal circles and a disc with a hole: the cut equals the cut of the curve built from the vertices
!<    of the copies moved by `curve%move`, bit for bit (without and with the oriented boundaries); the segment map and
!<    the labels of the wall segments agree; the solid area of the cut cells is the sum of the bodies' areas; the
!<    device path (built without OpenACC: the device loops on the host) equals the host path bit for bit, copied back
!<    and resident, with the same labels; boundaries with `device` are refused;
!< 2. the two circles in one cell: the cell's piece lists both bodies;
!< 3. overlaps refused: crossing circles, a circle nested in another (both orders), coincident circles, two squares
!<    overlapping along collinear edges of the same direction (no proper crossing, no vertex strictly inside);
!< 4. touching allowed: two squares along an edge (aligned, shifted) and at a corner: every cell trusted, the solid area
!<    2 to round-off; curves in planes of different normal axes refused.
use fossil,                   only : curve_object, cut_cells_2d_t, rigid_pose_2d_t, pose_2d, piece_bodies_t, &
                                     merge_posed_bodies, cut_cells_2d_bodies, cut_walls_2d_t, STATUS_OK, STATUS_INVALID_INPUT
use fossil_cut_cells_compare, only : cut_cells_2d_mismatches, piece_bodies_mismatches
use penf,                     only : I4P, I8P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P)
type(curve_object)    :: circ, disc, small, square
type(rigid_pose_2d_t) :: p2(2)
real(R8P), allocatable :: x(:), y(:)
integer(I4P)          :: i
logical               :: ok

ok = .true.
call polygon(circ, 0.5_R8P, 48_I4P)
call polygon(small, 0.2_R8P, 24_I4P)
call holed_disc(disc)
call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
p2(1) = pose_2d(0.3_R8P, [-0.62_R8P, 0.02_R8P]) ; p2(2) = pose_2d(-1.1_R8P, [0.61_R8P, -0.05_R8P])
x = [(-1.31_R8P + 0.1_R8P * i, i = 0, 27)] ; y = [(-0.73_R8P + 0.1_R8P * i, i = 0, 15)]
call twin_case([circ, circ], 'two circles')
call twin_case([circ, disc], 'circle and disc with a hole')
call one_cell_case
call overlap_case
call touch_case
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine twin_case(cc, label)
   !< Case 1: equal to the moved copies' curve, labels, area, device twin.
   type(curve_object), intent(in) :: cc(2)  !< Bodies.
   character(*),       intent(in) :: label  !< Case.
   type(curve_object)             :: a, mg  !< Reference, merged.
   type(cut_cells_2d_t)           :: ch, cb, cd !< Reference, bodies, device.
   type(cut_walls_2d_t)           :: wb     !< Wall segments.
   type(piece_bodies_t)           :: pb, pd !< Labels.
   integer(I4P), allocatable      :: sb(:), ss(:) !< Segment map.
   integer(I4P)                   :: m, st, k, r, n1
   real(R8P)                      :: as, a12

   m = 0_I4P
   call moved_union(cc, p2, a)
   call a%cut_cells(x, y, ch, status=st)
   call cut_cells_2d_bodies(cc, p2, x, y, cb, status=st, merged=mg, segment_body=sb, segment_sub_body=ss)
   m = m + abs(st) + cut_cells_2d_mismatches(ch, cb) + vertex_mismatches(a, mg)
   n1 = cc(1)%get_vertices_number()
   if (any(sb(1:n1) /= 1_I4P) .or. any(sb(n1 + 1:) /= 2_I4P)) m = m + 1_I4P
   do k = 1, size(ss)
      if (ss(k) /= 1_I4P) m = m + 1_I4P ! one body per curve (the hole belongs to the disc)
   enddo
   ! with boundaries (the labels turn them on): the same boundaries as the reference's
   call a%cut_cells(x, y, ch, status=st, boundaries=.true.)
   call cut_cells_2d_bodies(cc, p2, x, y, cb, status=st, piece_bodies=pb, walls=wb)
   m = m + abs(st) + seg_mismatches(ch, cb)
   if (pb%np /= cb%np .or. pb%nw /= wb%n) m = m + 1_I4P
   do k = 1, wb%n
      if (pb%wall_body(k) /= sb(wb%segment(k))) m = m + 1_I4P
   enddo
   as = solid_area(cb, x, y) ; a12 = cc(1)%get_area() + cc(2)%get_area()
   if (abs(as - a12) > 1.e-14_R8P) m = m + 1_I4P
   print '(A,A,I0,A,I0,A,I0,A,ES9.2,A,I0)', label, ': cells ', size(cb%alpha), ', pieces ', cb%np, ', wall segments ', &
         wb%n, ', |A_solid - A1 - A2| ', abs(as - a12), ', mismatches ', m
   ok = ok .and. m == 0_I4P
   ! device twin: copied back and resident, labels
   m = 0_I4P
   call cut_cells_2d_bodies(cc, p2, x, y, cb, status=st)
   do r = 1, 2
      call cut_cells_2d_bodies(cc, p2, x, y, cd, status=st, piece_bodies=pd, device=.true., resident=r == 2)
      if (r == 2) call cd%update_host
      m = m + abs(st) + cut_cells_2d_mismatches(cb, cd) + piece_bodies_mismatches(pb, pd)
   enddo
   call cd%exit_device
   call cut_cells_2d_bodies(cc, p2, x, y, cd, status=st, boundaries=.true., device=.true.)
   if (st /= STATUS_INVALID_INPUT .or. cd%nx /= 0_I4P) m = m + 1_I4P
   print '(A,A,I0)', label, ', device path: mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine twin_case

   subroutine one_cell_case
   !< Case 2: both circles in one cell, one piece bounded by both.
   type(cut_cells_2d_t) :: cb
   type(piece_bodies_t) :: pb
   real(R8P)            :: gx(3), gy(3)
   integer(I4P)         :: m, st, p

   m = 0_I4P
   gx = [-1.4_R8P, 1.4_R8P, 1.7_R8P] ; gy = [-0.8_R8P, 0.8_R8P, 1.1_R8P]
   call cut_cells_2d_bodies([circ, circ], p2, gx, gy, cb, status=st, piece_bodies=pb)
   if (st /= STATUS_OK .or. cb%pieces(1, 1) /= 1_I4P) then
      m = m + 1_I4P
   else
      p = cb%first_piece(1, 1)
      if (pb%start(p + 1) - pb%start(p) /= 2_I4P) then
         m = m + 1_I4P
      else
         if (any(pb%body(pb%start(p):pb%start(p) + 1) /= [1_I4P, 2_I4P])) m = m + 1_I4P
      endif
   endif
   print '(A,I0)', 'two circles in one cell: mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine one_cell_case

   subroutine overlap_case
   !< Case 3: overlaps refused.
   type(rigid_pose_2d_t) :: q(2)
   type(curve_object)    :: other
   integer(I4P)          :: m

   m = 0_I4P
   q(1) = pose_2d(0._R8P, [-0.3_R8P, 0._R8P]) ; q(2) = pose_2d(0.2_R8P, [0.3_R8P, 0.01_R8P])
   m = m + refused('crossing circles', [circ, circ], q)
   q(1) = pose_2d(0._R8P) ; q(2) = pose_2d(0.4_R8P, [0.1_R8P, 0.05_R8P])
   m = m + refused('circle in a circle', [circ, small], q)
   m = m + refused('circle in a circle (reversed)', [small, circ], q(2:1:-1))
   q(1) = p2(1) ; q(2) = p2(1)
   m = m + refused('coincident circles', [circ, circ], q)
   q(1) = pose_2d(0._R8P) ; q(2) = pose_2d(0._R8P, [0.5_R8P, 0._R8P])
   m = m + refused('squares on collinear edges', [square, square], q)
   call other%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P], &
                          normal_axis=1_I4P)
   q(2) = pose_2d(0._R8P, [5._R8P, 0._R8P])
   m = m + refused('curves in different planes', [square, other], q)
   ok = ok .and. m == 0_I4P
   endsubroutine overlap_case

   function refused(label, cc, q) result(m)
   !< 1 unless the cut of cc posed by q is refused (and the cut reset).
   character(*),          intent(in) :: label !< Case.
   type(curve_object),    intent(in) :: cc(:) !< Bodies.
   type(rigid_pose_2d_t), intent(in) :: q(:)  !< Poses.
   integer(I4P)                      :: m     !< Mismatches.
   type(cut_cells_2d_t)              :: cb    !< Cut.
   integer(I4P)                      :: st, nc !< Status, contacts.

   call cut_cells_2d_bodies(cc, q, x, y, cb, status=st, contacts=nc)
   m = 0_I4P ; if (st /= STATUS_INVALID_INPUT .or. cb%nx /= 0_I4P .or. allocated(cb%alpha)) m = 1_I4P
   print '(A,A,I0,A,I0)', label, ': refused (contacts ', nc, '), mismatches ', m
   endfunction refused

   subroutine touch_case
   !< Case 4: touching squares allowed, every cell trusted, area 2.
   real(R8P),     parameter :: T(2, 3) = reshape([1._R8P, 0._R8P, 1._R8P, 0.37_R8P, 1._R8P, 1._R8P], [2, 3])
   character(12), parameter :: NAMES(3) = [character(12) :: 'edge', 'edge shifted', 'corner']
   real(R8P),     parameter :: H(2) = [0.13_R8P, 0.125_R8P]
   type(rigid_pose_2d_t)    :: q(2)
   type(cut_cells_2d_t)     :: cb
   type(piece_bodies_t)     :: pb
   real(R8P), allocatable   :: gx(:)
   real(R8P)                :: as
   integer(I4P)             :: c, g, m, st, nc, nun, nboth, p

   do c = 1, 3
      do g = 1, 2
         m = 0_I4P
         q(1) = pose_2d(0._R8P) ; q(2) = pose_2d(0._R8P, T(:, c))
         gx = [(-0.27_R8P + H(g) * i, i = 0, nint(2.6_R8P / H(g)))]
         call cut_cells_2d_bodies([square, square], q, gx, gx, cb, status=st, piece_bodies=pb, contacts=nc)
         if (st /= STATUS_OK) then
            m = m + 1_I4P ; nun = -1_I4P ; as = 0._R8P ; nboth = 0_I4P
         else
            nun = count(cb%pieces < 0_I4P .or. .not. cb%closed)
            as = solid_area(cb, gx, gx)
            if (nun /= 0_I4P .or. abs(as - 2._R8P) > 1.e-14_R8P .or. nc == 0_I4P) m = m + 1_I4P
            nboth = 0_I4P
            do p = 1, pb%np
               if (pb%start(p + 1) - pb%start(p) == 2_I4P) nboth = nboth + 1_I4P
            enddo
            if (nboth == 0_I4P) m = m + 1_I4P
         endif
         print '(A,A,A,F6.3,A,I0,A,I0,A,I0,A,ES9.2,A,I0)', 'touching squares, ', NAMES(c), ', h ', H(g), ': contacts ', nc, &
               ', untrusted ', nun, ', pieces of both ', nboth, ', |A_solid - 2| ', abs(as - 2._R8P), ', mismatches ', m
         ok = ok .and. m == 0_I4P
      enddo
   enddo
   endsubroutine touch_case

   subroutine moved_union(cc, q, a)
   !< The reference: copies moved by `curve%move`, their vertices in one curve (loops in order).
   type(curve_object),    intent(in)    :: cc(:) !< Bodies.
   type(rigid_pose_2d_t), intent(in)    :: q(:)  !< Poses.
   type(curve_object),    intent(inout) :: a     !< Union curve.
   type(curve_object)                   :: w     !< Moved copy.
   real(R8P),    allocatable            :: u(:), v(:) !< Vertices.
   integer(I4P), allocatable            :: sz(:) !< Loop sizes.
   integer(I4P)                         :: c, l, i, f, e !< Work.
   real(R8P)                            :: uv(2) !< Vertex.

   allocate(u(0), v(0), sz(0))
   do c = 1, size(cc)
      w = cc(c) ; call w%move(q(c))
      do l = 1, w%get_loops_number()
         call w%get_loop_range(l, f, e)
         sz = [sz, e - f + 1_I4P]
         do i = f, e
            uv = w%get_vertex(i) ; u = [u, uv(1)] ; v = [v, uv(2)]
         enddo
      enddo
   enddo
   call a%adopt_loops(u=u, v=v, loop_sizes=sz)
   endsubroutine moved_union

   function vertex_mismatches(a, b) result(m)
   !< Vertices of a and b differing in a bit (or a count mismatch).
   type(curve_object), intent(in) :: a, b !< Curves.
   integer(I4P)                   :: m    !< Mismatches.
   integer(I4P)                   :: i    !< Work.

   m = 0_I4P
   if (a%get_vertices_number() /= b%get_vertices_number() .or. a%get_loops_number() /= b%get_loops_number()) then
      m = 1_I4P ; return
   endif
   do i = 1, a%get_vertices_number()
      if (any(transfer(a%get_vertex(i), 0_I8P, 2) /= transfer(b%get_vertex(i), 0_I8P, 2))) m = m + 1_I4P
   enddo
   endfunction vertex_mismatches

   function seg_mismatches(a, b) result(m)
   !< The oriented boundaries of a and b differing in a bit (start, segments, wall flags, sources).
   type(cut_cells_2d_t), intent(in) :: a, b !< Cuts with boundaries.
   integer(I4P)                     :: m    !< Mismatches.

   m = 0_I4P
   if (.not. (allocated(a%piece_seg) .and. allocated(b%piece_seg))) then
      m = 1_I4P ; return
   endif
   if (any(shape(a%piece_seg) /= shape(b%piece_seg)) .or. size(a%piece_seg_start) /= size(b%piece_seg_start)) then
      m = 1_I4P ; return
   endif
   m = count(a%piece_seg_start /= b%piece_seg_start, kind=I4P) + count(a%piece_seg_src /= b%piece_seg_src, kind=I4P) + &
       count(a%piece_seg_wall .neqv. b%piece_seg_wall, kind=I4P) + &
       count(transfer(a%piece_seg, 0_I8P, size(a%piece_seg)) /= transfer(b%piece_seg, 0_I8P, size(b%piece_seg)), kind=I4P)
   endfunction seg_mismatches

   function solid_area(c, gx, gy) result(s)
   !< Solid area of the cut cells: sum of (1 - alpha) |cell| over the cells with solid.
   type(cut_cells_2d_t), intent(in) :: c          !< Cut.
   real(R8P),            intent(in) :: gx(:), gy(:) !< Cell edges.
   real(R8P)                        :: s          !< Solid area.
   integer(I4P)                     :: i, j       !< Work.

   s = 0._R8P
   do j = 1, c%ny
      do i = 1, c%nx
         if (c%alpha(i, j) < 1._R8P) s = s + (1._R8P - c%alpha(i, j)) * (gx(i + 1) - gx(i)) * (gy(j + 1) - gy(j))
      enddo
   enddo
   endfunction solid_area

   subroutine polygon(c, r, n)
   !< Regular n-gon of radius r about the origin, counter-clockwise.
   type(curve_object), intent(inout) :: c !< Curve.
   real(R8P),          intent(in)    :: r !< Radius.
   integer(I4P),       intent(in)    :: n !< Vertices.
   real(R8P)                         :: u(n), v(n) !< Vertices.
   integer(I4P)                      :: k !< Work.

   do k = 1, n
      u(k) = r * cos(2._R8P * PI * (k - 1) / n) ; v(k) = r * sin(2._R8P * PI * (k - 1) / n)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[n])
   endsubroutine polygon

   subroutine holed_disc(c)
   !< A 40-gon of radius 0.45 with a clockwise 16-gon hole of radius 0.15.
   type(curve_object), intent(inout) :: c !< Curve.
   real(R8P)                         :: u(56), v(56) !< Vertices.
   integer(I4P)                      :: k !< Work.

   do k = 1, 40
      u(k) = 0.45_R8P * cos(2._R8P * PI * (k - 1) / 40) ; v(k) = 0.45_R8P * sin(2._R8P * PI * (k - 1) / 40)
   enddo
   do k = 1, 16
      u(40 + k) = 0.15_R8P * cos(-2._R8P * PI * (k - 1) / 16) ; v(40 + k) = 0.15_R8P * sin(-2._R8P * PI * (k - 1) / 16)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[40_I4P, 16_I4P])
   endsubroutine holed_disc
endprogram fossil_test_cut_bodies_2d
