!< FOSSIL test: fluid pieces of the 3-D cut cells (issue #40).
program fossil_test_cut_pieces_3d
!< FOSSIL test: fluid pieces of the 3-D cut cells (issue #40).
!<
!< 1. Prisms (extrusions along z): on every cell strictly inside the height the 3-D pieces equal the
!<    2-D pieces of the section (`slice_surface` + 2-D cut cells), an independent count. Bodies: the
!<    5-degree wedge rotated about z (cells crossed by the thin wedge have 2 pieces), a square tube
!<    with a hole narrower than a cell (the fluid inside the hole meets the cell boundary in two
!<    disks joined inside the cell: one piece, while the cell-boundary fluid alone would give two),
!<    and the L-prism.
!< 2. Thin plate (1e-3) across cells: 2 pieces; cells over the plate edge: 1.
!< 3. Rods along x, y and z through a cell (solid islands traced on two opposite faces): 1 piece.
!< 4. A hollow box (a void inside the solid): a cell holding the void inside the solid wall, 1 piece;
!<    a cell holding the whole body, 2; a small box inside a cell (an inclusion), 1.
!< 6. A square passage along (1, -1, 0) cutting the corner of a cell at its x and y high faces: one
!<    piece. Its wall set touches the cell boundary only around two disks on high faces and has a
!<    positive signed volume, so only the touch rule (not the void rule) counts it.
!< 5. Unit cube on a generic grid: pieces = 1 where alpha > 0, else 0. On a grid through its faces,
!<    edges and corners: no cell flagged, pieces 0 or 1, and equal to the 2-D pieces of the section on
!<    the layer (0, 0.5]. A cube face on the high plane of a cell lies inside that cell (a point on a
!<    plane belongs to the low side), so the open high face (aperture 1) bounds a zero-volume fluid
!<    piece: 1 piece with alpha = 0, as in 2-D. The cube without one facet (open): the cells the
!<    hole's edges overlap, and those whose fraction the open surface corrupts (outside [0, 1]), are
!<    flagged (-1); the others keep their count.

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, rigid_pose_t, pose_from_axis_angle
use fossil_curve_slice,       only : slice_surface
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface, soup_remove
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(cut_cells_3d_t)      :: cut
type(soup_t)              :: so
real(R8P), allocatable    :: x(:), y(:), z(:)
integer(I4P)              :: i, mism(3), ncells(3), n2(3), a, st
logical                   :: ok(6)

call bench_bodies(b)

! 1. prisms against the 2-D pieces of the section
call to_surface(b(2)%soup, s) ! wedge 5, z in [0, 1]
call s%move(pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.37_R8P))
x = [(-0.13_R8P + 0.17_R8P * i, i = 0, 7)] ; y = [(-0.21_R8P + 0.043_R8P * i, i = 0, 21)] ; z = [0.2_R8P, 0.5_R8P, 0.8_R8P]
call prism_check(s, x, y, z, mism(1), ncells(1), n2(1))
so = tube([0._R8P, 3._R8P], [1._R8P, 2._R8P], [0._R8P, 3._R8P]) ; call to_surface(so, s)
x = [-0.5_R8P, 1.2_R8P, 1.8_R8P, 3.4_R8P] ; y = x ; z = [0.5_R8P, 1.5_R8P, 2.5_R8P]
call prism_check(s, x, y, z, mism(2), ncells(2), n2(2)) ! corner cells: outer fluid + a corner of the hole
x = [-0.5_R8P, 0.8_R8P, 2.2_R8P, 3.4_R8P] ; y = x
call cut_cells_3d(s, x, y, z, cut)
ok(1) = all(cut%pieces(2, 2, :) == 1) ! the narrow hole: one piece
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [0.1_R8P, 0.45_R8P, 0.9_R8P]
call prism_check(s, x, y, z, mism(3), ncells(3), n2(3))
ok(1) = ok(1) .and. all(mism == 0) .and. n2(1) > 0 .and. n2(2) > 0
print '(A,3I4,A,3I4,A,3I4,A,L1)', 'prisms (wedge, tube, L): 3-D vs 2-D mismatches', mism, '  cells', ncells, &
      '  cells with 2+ pieces', n2, '  narrow hole 1 piece, ok ', ok(1)

! 2. thin plate
call to_surface(b(5)%soup, s) ! [0, 1]^2 x [0, 1e-3]
call cut_cells_3d(s, [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P, 1.3_R8P], [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P], &
                  [-0.5_R8P, 0.5_R8P], cut)
ok(2) = all(cut%pieces(1:3, :, 1) == 2) .and. all(cut%pieces(4, :, 1) == 1)
print '(A,L1)', 'thin plate: 2 pieces across it, 1 over its edge: ', ok(2)

! 3. rods through a cell along each axis
ok(3) = .true.
do a = 1, 3
   so = box_along(a) ; call to_surface(so, s)
   call cut_cells_3d(s, [1._R8P, 2._R8P], [1._R8P, 2._R8P], [1._R8P, 2._R8P], cut)
   ok(3) = ok(3) .and. cut%pieces(1, 1, 1) == 1
enddo
print '(A,L1)', 'rods along x, y, z through a cell (islands on two faces): 1 piece: ', ok(3)

! 4. void and inclusion
so = hollow_box() ; call to_surface(so, s)
call cut_cells_3d(s, [0.5_R8P, 2.5_R8P], [0.5_R8P, 2.5_R8P], [0.5_R8P, 2.5_R8P], cut)
ok(4) = cut%pieces(1, 1, 1) == 1 .and. abs(s%get_volume() - 26._R8P) < 1.e-12_R8P
call cut_cells_3d(s, [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], cut)
ok(4) = ok(4) .and. cut%pieces(1, 1, 1) == 2
call cut_cells_3d(s, [0.5_R8P, 2.5_R8P], [0.5_R8P, 2.5_R8P], [-1._R8P, 2.5_R8P], cut)
ok(4) = ok(4) .and. cut%pieces(1, 1, 1) == 2
so = box_soup([0.4_R8P, 0.4_R8P, 0.4_R8P], [0.6_R8P, 0.6_R8P, 0.6_R8P], .false.) ; call to_surface(so, s)
call cut_cells_3d(s, [0._R8P, 1._R8P], [0._R8P, 1._R8P], [0._R8P, 1._R8P], cut)
ok(4) = ok(4) .and. cut%pieces(1, 1, 1) == 1
print '(A,L1)', 'void in the solid wall 1, whole hollow body 2, void + fluid below 2, inclusion 1: ', ok(4)

! 5. unit cube, generic and degenerate grids
call to_surface(b(1)%soup, s)
x = [(-0.23_R8P + 0.3_R8P * i, i = 0, 6)]
call cut_cells_3d(s, x, x, x, cut, st)
ok(5) = all(cut%pieces == merge(1, 0, cut%alpha > 0._R8P))
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call cut_cells_3d(s, x, x, x, cut, st)
ok(5) = ok(5) .and. all(cut%pieces >= 0) .and. all(cut%pieces <= 1)
call prism_check(s, x, x, [0._R8P, 0.5_R8P], mism(1), ncells(1), n2(1))
ok(5) = ok(5) .and. mism(1) == 0
so = b(1)%soup ; call soup_remove(so, 1) ; call to_surface(so, s)
x = [(-0.23_R8P + 0.3_R8P * i, i = 0, 6)]
call cut_cells_3d(s, x, x, x, cut, st)
ok(5) = ok(5) .and. count(cut%pieces == -1) > 0 .and. count(cut%pieces == -1) < size(cut%pieces) / 2 .and. &
        all(cut%pieces == -1 .or. cut%pieces == merge(1, 0, cut%alpha > 0._R8P))
print '(A,L1,A,I0,A,I0)', 'unit cube: generic grid pieces = (alpha > 0); aligned grid none flagged, 0 or 1, = 2-D: ', ok(5), &
      '  mismatches ', mism(1), '  open cube: flagged cells ', count(cut%pieces == -1)

! 6. passage through two high faces near their common edge
so = tube([0._R8P, 3._R8P], [1.25_R8P, 1.75_R8P], [-3._R8P, 3._R8P]) ! hole 0.5 wide, clear of every cell edge
call frame_map(so, (2.35_R8P / sqrt(2._R8P) - 0.25_R8P))
call to_surface(so, s)
call cut_cells_3d(s, [0._R8P, 1.2_R8P], [0._R8P, 1.2_R8P], [-0.6_R8P, 0.6_R8P], cut)
ok(6) = cut%pieces(1, 1, 1) == 1 .and. abs(cut%alpha(1, 1, 1) - 0.25_R8P * (0.757106781186548_R8P**2 - &
        0.05_R8P**2) / 1.728_R8P) < 1.e-12_R8P ! fluid: the strip x + y in [2.4 - L1, 2.4 - L2], z in [-0.25, 0.25]
print '(A,I0,A,F8.5,A,L1)', 'passage along (1, -1, 0) through two high faces: pieces ', cut%pieces(1, 1, 1), &
      '  alpha ', cut%alpha(1, 1, 1), '  ok ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine prism_check(sf, xe, ye, ze, nmis, ncell, ntwo)
   !< Compare the 3-D pieces of every cell with the 2-D pieces of the section at mid-height.
   type(surface_stl_object), intent(in)  :: sf         !< Prism (extruded along z).
   real(R8P),                intent(in)  :: xe(:), ye(:), ze(:) !< Edges (ze inside the height).
   integer(I4P),             intent(out) :: nmis, ncell, ntwo !< Mismatches, cells, cells with 2+ pieces.
   type(cut_cells_3d_t)  :: c3   !< 3-D cut cells.
   type(curve_object)    :: sec  !< Section.
   type(cut_cells_2d_t)  :: c2   !< Its cut cells.
   integer(I4P)          :: ii, jj, kk !< Cells.

   call cut_cells_3d(sf, xe, ye, ze, c3)
   call slice_surface(surface=sf, normal_axis=3_I4P, level=0.5_R8P * (ze(1) + ze(2)), curve=sec, on_plane_low=.true.)
   call sec%cut_cells(x=xe, y=ye, cut=c2)
   nmis = 0 ; ncell = 0 ; ntwo = 0
   do kk = 1, size(ze) - 1
      do jj = 1, size(ye) - 1
         do ii = 1, size(xe) - 1
            ncell = ncell + 1
            if (c3%pieces(ii, jj, kk) /= c2%pieces(ii, jj)) nmis = nmis + 1
            if (c3%pieces(ii, jj, kk) >= 2) ntwo = ntwo + 1
         enddo
      enddo
   enddo
   endsubroutine prism_check

   subroutine add_quad(sp, p1, p2, p3, p4, nout)
   !< Append the quad p1 p2 p3 p4 (planar) as two triangles, oriented along nout.
   type(soup_t), intent(inout) :: sp                 !< Soup.
   real(R8P),    intent(in)    :: p1(3), p2(3), p3(3), p4(3), nout(3) !< Corners, outward direction.
   real(R8P)                   :: t(3, 3, 2), nq(3)  !< Triangles, normal.

   nq = [(p2(2) - p1(2)) * (p3(3) - p1(3)) - (p2(3) - p1(3)) * (p3(2) - p1(2)), &
         (p2(3) - p1(3)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(3) - p1(3)), &
         (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))]
   if (dot_product(nq, nout) > 0._R8P) then
      t(:, 1, 1) = p1 ; t(:, 2, 1) = p2 ; t(:, 3, 1) = p3 ; t(:, 1, 2) = p1 ; t(:, 2, 2) = p3 ; t(:, 3, 2) = p4
   else
      t(:, 1, 1) = p1 ; t(:, 2, 1) = p3 ; t(:, 3, 1) = p2 ; t(:, 1, 2) = p1 ; t(:, 2, 2) = p4 ; t(:, 3, 2) = p3
   endif
   if (.not. allocated(sp%v)) allocate(sp%v(3, 3, 0))
   sp%v = reshape([reshape(sp%v, [9 * size(sp%v, 3)]), reshape(t, [18])], [3, 3, size(sp%v, 3) + 2])
   endsubroutine add_quad

   function box_soup(lo, hi, inward) result(sp)
   !< Axis-aligned box [lo, hi], normals outward (or inward).
   real(R8P), intent(in) :: lo(3), hi(3) !< Corners.
   logical,   intent(in) :: inward       !< Inward normals.
   type(soup_t)          :: sp           !< Soup.
   real(R8P)             :: f            !< Orientation.

   f = 1._R8P ; if (inward) f = -1._R8P
   call add_quad(sp, [lo(1), lo(2), lo(3)], [lo(1), hi(2), lo(3)], [lo(1), hi(2), hi(3)], [lo(1), lo(2), hi(3)], &
                 [-f, 0._R8P, 0._R8P])
   call add_quad(sp, [hi(1), lo(2), lo(3)], [hi(1), hi(2), lo(3)], [hi(1), hi(2), hi(3)], [hi(1), lo(2), hi(3)], &
                 [f, 0._R8P, 0._R8P])
   call add_quad(sp, [lo(1), lo(2), lo(3)], [hi(1), lo(2), lo(3)], [hi(1), lo(2), hi(3)], [lo(1), lo(2), hi(3)], &
                 [0._R8P, -f, 0._R8P])
   call add_quad(sp, [lo(1), hi(2), lo(3)], [hi(1), hi(2), lo(3)], [hi(1), hi(2), hi(3)], [lo(1), hi(2), hi(3)], &
                 [0._R8P, f, 0._R8P])
   call add_quad(sp, [lo(1), lo(2), lo(3)], [hi(1), lo(2), lo(3)], [hi(1), hi(2), lo(3)], [lo(1), hi(2), lo(3)], &
                 [0._R8P, 0._R8P, -f])
   call add_quad(sp, [lo(1), lo(2), hi(3)], [hi(1), lo(2), hi(3)], [hi(1), hi(2), hi(3)], [lo(1), hi(2), hi(3)], &
                 [0._R8P, 0._R8P, f])
   endfunction box_soup

   function box_along(ax) result(sp)
   !< Rod [1.4, 1.6]^2 x [-1, 4] along axis ax.
   integer(I4P), intent(in) :: ax     !< Axis.
   type(soup_t)             :: sp     !< Soup.
   real(R8P)                :: lo(3), hi(3) !< Corners.

   lo = 1.4_R8P ; hi = 1.6_R8P ; lo(ax) = -1._R8P ; hi(ax) = 4._R8P
   sp = box_soup(lo, hi, .false.)
   endfunction box_along

   function hollow_box() result(sp)
   !< Box [0, 3]^3 with the void [1, 2]^3 (inner faces facing the void).
   type(soup_t) :: sp    !< Soup.
   type(soup_t) :: inner !< Void surface.

   sp = box_soup([0._R8P, 0._R8P, 0._R8P], [3._R8P, 3._R8P, 3._R8P], .false.)
   inner = box_soup([1._R8P, 1._R8P, 1._R8P], [2._R8P, 2._R8P, 2._R8P], .true.)
   sp%v = reshape([reshape(sp%v, [size(sp%v)]), reshape(inner%v, [size(inner%v)])], [3, 3, size(sp%v, 3) + size(inner%v, 3)])
   endfunction hollow_box

   subroutine frame_map(sp, wc)
   !< Map a tube built along z (cross-section centre (1.5, 1.5)) to the axis along e_u = (1, -1, 0)/sqrt 2
   !< through wc e_w with e_w = (1, 1, 0)/sqrt 2: X -> e_w, Y -> e_z, Z -> e_u (right-handed).
   type(soup_t), intent(inout) :: sp !< Soup.
   real(R8P),    intent(in)    :: wc !< Axis offset along e_w.
   real(R8P)                   :: ew(3), ez(3), eu(3), p(3) !< Frame, point.
   integer(I4P)                :: f, k !< Work.

   ew = [1._R8P, 1._R8P, 0._R8P] / sqrt(2._R8P) ; ez = [0._R8P, 0._R8P, 1._R8P] ; eu = [1._R8P, -1._R8P, 0._R8P] / sqrt(2._R8P)
   do f = 1, size(sp%v, 3)
      do k = 1, 3
         p = sp%v(:, k, f)
         sp%v(:, k, f) = wc * ew + (p(1) - 1.5_R8P) * ew + (p(2) - 1.5_R8P) * ez + p(3) * eu
      enddo
   enddo
   endsubroutine frame_map

   function tube(o, h, zr) result(sp)
   !< Square tube: outer [o(1), o(2)]^2, hole [h(1), h(2)]^2, z in zr.
   real(R8P), intent(in) :: o(2), h(2), zr(2) !< Outer and hole spans, height.
   type(soup_t)          :: sp                !< Soup.
   real(R8P)             :: oc(2, 4), hc(2, 4), c(2), mid(2) !< Corners, centre, edge midpoint.
   integer(I4P)          :: k, k1, m          !< Work.

   oc = reshape([o(1), o(1), o(2), o(1), o(2), o(2), o(1), o(2)], [2, 4])
   hc = reshape([h(1), h(1), h(2), h(1), h(2), h(2), h(1), h(2)], [2, 4])
   c = 0.5_R8P * [o(1) + o(2), o(1) + o(2)]
   do k = 1, 4
      k1 = mod(k, 4) + 1
      mid = 0.5_R8P * (oc(:, k) + oc(:, k1))
      call add_quad(sp, [oc(:, k), zr(1)], [oc(:, k1), zr(1)], [oc(:, k1), zr(2)], [oc(:, k), zr(2)], [mid - c, 0._R8P])
      mid = 0.5_R8P * (hc(:, k) + hc(:, k1))
      call add_quad(sp, [hc(:, k), zr(1)], [hc(:, k1), zr(1)], [hc(:, k1), zr(2)], [hc(:, k), zr(2)], [c - mid, 0._R8P])
      do m = 1, 2
         call add_quad(sp, [oc(:, k), zr(m)], [oc(:, k1), zr(m)], [hc(:, k1), zr(m)], [hc(:, k), zr(m)], &
                       [0._R8P, 0._R8P, merge(-1._R8P, 1._R8P, m == 1)])
      enddo
   enddo
   endfunction tube
endprogram fossil_test_cut_pieces_3d
