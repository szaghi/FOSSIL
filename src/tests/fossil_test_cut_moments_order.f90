!< FOSSIL test: the cut-cell moments do not depend on the order of the loops (issue #95).
program fossil_test_cut_moments_order
!< FOSSIL test: the cut-cell moments do not depend on the order of the loops (issue #95).
!<
!< Four solid squares on the corners of one cell, the loops in all 24 orders: the union of the fluid sub-intervals
!< joins up to four wall chains, and a root that had chains charged to it before losing a later union lost their
!< moments (before #95: zeroth moment 0.58 or 0.52 instead of 0.64 in 16 orders, errors up to 1.2 with the wall
!< moments; the 3-D triangle body 4.0e-2 of the cell). Measured after the fix (gfortran 16.2, -O0 and -O2) next to the
!< bounds:
!< 1. 2-D, `moments_order=2`: per piece, the zeroth moment is the area and the zeroth wall moment the wetted length
!<    (1e-14 of the cell; 2.2e-16); the cell's moments and wall moments equal the closed form (1e-14; 4.4e-16); every
!<    order gives the moments of the first (relative 1e-14; 1.2e-17: the unions add the terms in another order);
!< 2. 3-D: the volume moments come from the x sections, whose loops are sorted (`slice_surface`), so the facet order
!<    does not matter but the geometry does. The four squares extruded along x through a cell sort in a lucky order;
!<    with the square at the corner (1, 1) replaced by a triangle reaching y = -0.9 (it sorts first) the cell's
!<    section lost the triangle's moments. Per piece: zeroth moment = volume (1e-14 of the cell), zeroth wall moment =
!<    wetted area, zeroth face moment = face fraction x face area (1e-14 of the face; 2.4e-15 over the three); first
!<    moments = volume x scaled centroid offset (1e-14 of the cell; 1.8e-16); the four-squares cell in closed form
!<    (1e-14; 5.6e-17);
!< 3. a cell that runs out of sub-interval or loop slots is reported open (`closed` false), every other cell closed,
!<    instead of an out-of-bounds write: 132 square holes in one cell (528 loop pieces > 512: `loopp`), and 488
!<    crossings of the cell's edge plus 30 triangular holes (the hole slots pass 512: `par`, `acc`, ...).

use fossil,              only : surface_stl_object, curve_object, cut_cells_2d_t
use fossil_cut_cells_3d, only : cut_cells_3d_t, cut_cells_3d
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

real(R8P), parameter      :: TOL = 1.e-14_R8P
type(curve_object)        :: c
type(cut_cells_2d_t)      :: c2, c2ref
type(cut_cells_3d_t)      :: c3
type(surface_stl_object)  :: s
real(R8P)                 :: cu(4, 4), cv(4, 4), u(16), v(16), cx(4), cy(4), g(4)
real(R8P), parameter      :: g3(5) = [-1._R8P, 0._R8P, 1._R8P, 2._R8P, 3._R8P]   !< 3-D y, z grid.
real(R8P), parameter      :: g3x(4) = [-0.37_R8P, 0.13_R8P, 0.61_R8P, 1.29_R8P] !< 3-D x grid.
real(R8P)                 :: err(7), m2ref(6), w2ref(6), m3ref(10)
integer(I4P)              :: p(4), a, b, d, e, k, n, norders, nfirst
logical                   :: ok, first

cx = [0._R8P, 1._R8P, 1._R8P, 0._R8P] ; cy = [0._R8P, 0._R8P, 1._R8P, 1._R8P]
do k = 1, 4
   cu(:, k) = cx(k) + 0.3_R8P * [-1._R8P, 1._R8P, 1._R8P, -1._R8P]
   cv(:, k) = cy(k) + 0.3_R8P * [-1._R8P, -1._R8P, 1._R8P, 1._R8P]
enddo
g = [-1._R8P, 0._R8P, 1._R8P, 2._R8P]
! the cell [0, 1]^2 (scaled X, Y in [-1/2, 1/2]) without the four corner squares of side 0.3: moments 1, X, Y, X^2, XY, Y^2
m2ref = [0.64_R8P, 0._R8P, 0._R8P, 1._R8P / 12._R8P - 4._R8P * 0.3_R8P * (0.125_R8P - 0.008_R8P) / 3._R8P, 0._R8P, 0._R8P]
m2ref(6) = m2ref(4)
w2ref = [2.4_R8P, 0._R8P, 0._R8P, 4._R8P * (0.04_R8P * 0.3_R8P + (0.125_R8P - 0.008_R8P) / 3._R8P), 0._R8P, 0._R8P]
w2ref(6) = w2ref(4)
! the 3-D cell [0.13, 0.61] x [0, 1]^2: the 2-D section times h_x = 0.48; moments 1, X, Y, Z, X^2, XY, XZ, Y^2, YZ, Z^2
m3ref = 0._R8P
m3ref(1) = 0.48_R8P * m2ref(1) ; m3ref(5) = 0.48_R8P * m2ref(1) / 12._R8P
m3ref(8) = 0.48_R8P * m2ref(4) ; m3ref(10) = 0.48_R8P * m2ref(6)
err = 0._R8P ; norders = 0 ; first = .true. ; ok = .true.
do a = 1, 4 ; do b = 1, 4 ; do d = 1, 4 ; do e = 1, 4
   p = [a, b, d, e]
   if (any([(count(p == k), k = 1, 4)] /= 1)) cycle
   norders = norders + 1
   ! 1. 2-D
   do k = 1, 4
      u(4 * k - 3:4 * k) = cu(:, p(k)) ; v(4 * k - 3:4 * k) = cv(:, p(k))
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 4_I4P, 4_I4P, 4_I4P])
   call c%cut_cells(g, g, c2, moments_order=2_I4P)
   do n = 1, c2%np
      err(1) = max(err(1), abs(c2%piece_moments(1, n) - c2%piece_area(n)), &
                           abs(c2%piece_wall_moments(1, n) - c2%piece_wet_length(n)))
      if (all(c2%piece_cell(:, n) == [2, 2])) then
         err(2) = max(err(2), maxval(abs(c2%piece_moments(:, n) - m2ref)), maxval(abs(c2%piece_wall_moments(:, n) - w2ref)))
      endif
   enddo
   if (first) then
      c2ref = c2 ; nfirst = c2%np ; first = .false.
   elseif (c2%np /= nfirst) then
      ok = .false.
   else
      err(3) = max(err(3), maxval(abs(c2%piece_moments - c2ref%piece_moments)) / maxval(abs(c2ref%piece_moments)), &
               maxval(abs(c2%piece_wall_moments - c2ref%piece_wall_moments)) / maxval(abs(c2ref%piece_wall_moments)))
   endif
enddo ; enddo ; enddo ; enddo
! 2. 3-D
call prisms(.false., s)
call cut_cells_3d(s, g3x, g3, g3, c3, moments_order=2_I4P)
call check3(.true.)
call prisms(.true., s)
call cut_cells_3d(s, g3x, g3, g3, c3, moments_order=2_I4P)
call check3(.false.)
print '(A,I0)', 'loop orders: ', norders
print '(A,6ES10.2)', 'max errors (2-D zeroth, 2-D closed form, 2-D orders, 3-D zeroth, 3-D closed form, '// &
   '3-D first): ', err(1:6)
! 3. slot overflows
call many_holes(err(7))
print '(A,I0)', 'overflow cells not open or other cells open: ', nint(err(7))
ok = ok .and. norders == 24 .and. all(err(1:6) <= TOL) .and. err(7) == 0._R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check3(squares)
   !< Checks of the 3-D cut c3 (squares: the four-squares body, with the closed form of its cell (2, 2, 2)).
   logical, intent(in)  :: squares !< The four-squares body.
   real(R8P), parameter :: HX(3) = [0.5_R8P, 0.48_R8P, 0.68_R8P] !< x widths of the cells.
   integer(I4P)         :: q, f, ax, b1, b2
   real(R8P)            :: h(3), af, cc(3)

   do q = 1, c3%np
      h = [HX(c3%piece_cell(1, q)), 1._R8P, 1._R8P]
      cc = [g3x(c3%piece_cell(1, q)), g3(c3%piece_cell(2, q)), g3(c3%piece_cell(3, q))] + 0.5_R8P * h
      err(4) = max(err(4), abs(c3%piece_moments(1, q) - c3%piece_volume(q)) / product(h), &
                           abs(c3%piece_wall_moments(1, q) - c3%piece_wet_area(q)) / maxval(h) ** 2)
      do f = 1, 6
         ax = (f + 1) / 2 ; b1 = mod(ax, 3) + 1 ; b2 = mod(ax + 1, 3) + 1 ; af = h(b1) * h(b2)
         err(4) = max(err(4), abs(c3%piece_face_moments(1, f, q) - c3%piece_face(f, q) * af) / af)
      enddo
      err(6) = max(err(6), maxval(abs(c3%piece_moments(2:4, q) - c3%piece_volume(q) * (c3%piece_centroid(:, q) - cc) / h)) &
                   / product(h))
      if (squares .and. all(c3%piece_cell(:, q) == [2, 2, 2])) err(5) = max(err(5), maxval(abs(c3%piece_moments(:, q) - m3ref)))
   enddo
   endsubroutine check3

   subroutine prisms(triangle, surf)
   !< Boxes x in [0, 1] around the corners of the (y, z) cell [0, 1]^2 (half side 0.3); with triangle, the box at the
   !< corner (1, 1) replaced by the prism of the triangle (1.25, 0.7), (1.25, 1.25), (-0.9, 2.5) (it crosses the cell's
   !< top at y = 0.892 and its right side at z = 0.909, and its vertex y = -0.9 sorts its section loop first).
   logical,                  intent(in)    :: triangle !< Replace the box at (1, 1).
   type(surface_stl_object), intent(inout) :: surf     !< Surface.
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: lo(3), hi(3), t(2, 3)
   integer(I4P)                            :: kb, nf, kv, kw

   allocate(f(48)) ; nf = 0
   do kb = 1, 4
      if (triangle .and. kb == 3) then
         t = reshape([1.25_R8P, 0.7_R8P, 1.25_R8P, 1.25_R8P, -0.9_R8P, 2.5_R8P], [2, 3]) ! counter-clockwise in (y, z)
         call tri(f(nf + 1), [1._R8P, t(:, 1)], [1._R8P, t(:, 2)], [1._R8P, t(:, 3)])
         call tri(f(nf + 2), [0._R8P, t(:, 1)], [0._R8P, t(:, 3)], [0._R8P, t(:, 2)])
         nf = nf + 2
         do kv = 1, 3
            kw = mod(kv, 3) + 1
            call tri(f(nf + 1), [0._R8P, t(:, kv)], [0._R8P, t(:, kw)], [1._R8P, t(:, kw)])
            call tri(f(nf + 2), [0._R8P, t(:, kv)], [1._R8P, t(:, kw)], [1._R8P, t(:, kv)])
            nf = nf + 2
         enddo
      else
         lo = [0._R8P, cx(kb) - 0.3_R8P, cy(kb) - 0.3_R8P]
         hi = [1._R8P, cx(kb) + 0.3_R8P, cy(kb) + 0.3_R8P]
         call box(lo, hi, f(nf + 1:nf + 12))
         nf = nf + 12
      endif
   enddo
   f = f(1:nf)
   call surf%adopt_facets(facets=f)
   endsubroutine prisms

   subroutine box(lo, hi, f)
   !< The 12 outward facets of the box [lo, hi].
   real(R8P),          intent(in)    :: lo(3), hi(3) !< Corners.
   type(facet_object), intent(inout) :: f(12)        !< Facets.
   real(R8P)                         :: cn(3, 8)
   integer(I4P), parameter           :: Q(4, 6) = reshape([1, 5, 7, 3,  2, 4, 8, 6,  1, 2, 6, 5,  3, 7, 8, 4, &
                                                           1, 3, 4, 2,  5, 6, 8, 7], [4, 6]) !< Faces, outward.
   integer(I4P)                      :: kc, kf

   do kc = 1, 8 ! corner kc: bit 0 = x, bit 1 = y, bit 2 = z
      cn(:, kc) = merge(hi, lo, [btest(kc - 1, 0), btest(kc - 1, 1), btest(kc - 1, 2)])
   enddo
   do kf = 1, 6
      call tri(f(2 * kf - 1), cn(:, Q(1, kf)), cn(:, Q(2, kf)), cn(:, Q(3, kf)))
      call tri(f(2 * kf), cn(:, Q(1, kf)), cn(:, Q(3, kf)), cn(:, Q(4, kf)))
   enddo
   endsubroutine box

   subroutine tri(f, p1, p2, p3)
   !< Facet p1 p2 p3.
   type(facet_object), intent(inout) :: f              !< Facet.
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3) !< Vertices.

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   endsubroutine tri

   subroutine many_holes(bad)
   !< Check 3: a solid square [0.5, 2.5]^2 on the grid [0, 3]^2 with holes in the cell (2, 2); bad counts the failures.
   real(R8P),    intent(out) :: bad
   real(R8P),    allocatable :: uu(:), vv(:)
   integer(I4P), allocatable :: ls(:)
   real(R8P)                 :: gg(4), a0, b0
   integer(I4P)              :: ih, jh, kz, nv
   integer(I4P), parameter   :: KZIG = 489 !< Zigzag vertices inside x in (1, 2) (488 crossings of y = 1).

   bad = 0._R8P ; gg = [0._R8P, 1._R8P, 2._R8P, 3._R8P]
   ! a. 132 square holes (12 x 11), side 0.02, clockwise: 528 loop pieces in one cell
   allocate(uu(4 + 4 * 132), vv(4 + 4 * 132), ls(133))
   uu(1:4) = [0.5_R8P, 2.5_R8P, 2.5_R8P, 0.5_R8P] ; vv(1:4) = [0.5_R8P, 0.5_R8P, 2.5_R8P, 2.5_R8P] ; ls(1) = 4 ; nv = 4
   do jh = 1, 11
      do ih = 1, 12
         a0 = 1.03_R8P + 0.075_R8P * (ih - 1) ; b0 = 1.03_R8P + 0.085_R8P * (jh - 1)
         uu(nv + 1:nv + 4) = [a0, a0, a0 + 0.02_R8P, a0 + 0.02_R8P] ; vv(nv + 1:nv + 4) = [b0, b0 + 0.02_R8P, b0 + 0.02_R8P, b0]
         nv = nv + 4 ; ls(1 + 12 * (jh - 1) + ih) = 4
      enddo
   enddo
   call c%adopt_loops(u=uu, v=vv, loop_sizes=ls)
   call c%cut_cells(gg, gg, c2, moments_order=2_I4P)
   call tally(bad)
   deallocate(uu, vv, ls)
   ! b. the bottom of the solid zigzags across y = 1 inside the cell (KZIG - 1 crossings: under the 508 sub-intervals),
   !    then 30 triangular holes (90 loop pieces): the hole slots overflow the 512 sub-intervals
   allocate(uu(4 + KZIG + 3 * 30), vv(4 + KZIG + 3 * 30), ls(31))
   uu(1) = 0.5_R8P ; vv(1) = 1.005_R8P ; nv = 1
   do kz = 1, KZIG
      nv = nv + 1 ; uu(nv) = 1._R8P + (kz - 0.5_R8P) / KZIG ; vv(nv) = merge(0.995_R8P, 1.005_R8P, mod(kz, 2) == 1)
   enddo
   ! counter-clockwise: the solid above the zigzag
   uu(nv + 1:nv + 3) = [2.5_R8P, 2.5_R8P, 0.5_R8P] ; vv(nv + 1:nv + 3) = [1.005_R8P, 2.5_R8P, 2.5_R8P]
   nv = nv + 3 ; ls(1) = nv
   do jh = 1, 5
      do ih = 1, 6
         a0 = 1.1_R8P + 0.13_R8P * (ih - 1) ; b0 = 1.3_R8P + 0.11_R8P * (jh - 1)
         uu(nv + 1:nv + 3) = [a0, a0, a0 + 0.05_R8P] ; vv(nv + 1:nv + 3) = [b0, b0 + 0.05_R8P, b0]
         nv = nv + 3 ; ls(1 + 6 * (jh - 1) + ih) = 3
      enddo
   enddo
   call c%adopt_loops(u=uu(1:nv), v=vv(1:nv), loop_sizes=ls)
   call c%cut_cells(gg, gg, c2, moments_order=2_I4P)
   call tally(bad)
   endsubroutine many_holes

   subroutine tally(bad)
   !< The cell (2, 2) must be open, every other cell closed.
   real(R8P), intent(inout) :: bad
   integer(I4P)             :: ic, jc

   do jc = 1, 3
      do ic = 1, 3
         if (c2%closed(ic, jc) .eqv. (ic == 2 .and. jc == 2)) bad = bad + 1._R8P
      enddo
   enddo
   endsubroutine tally
endprogram fossil_test_cut_moments_order
