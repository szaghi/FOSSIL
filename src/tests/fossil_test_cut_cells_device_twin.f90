!< FOSSIL test: the device twin of the 2-D cut cells equals the host path bit for bit (issue #100 step 2).
program fossil_test_cut_cells_device_twin
!< FOSSIL test: the device twin of the 2-D cut cells equals the host path bit for bit (issue #100 step 2).
!<
!< Built without OpenACC (this CI build), `curve%cut_cells_device` runs its loops on the host: the test catches driver,
!< offset and prefix-sum bugs of the twin (`compute_cut_cells_2d_device`), not device codegen (the curve bench does
!< that on the GPU). Every output of the twin (`fossil_cut_cells_compare`: cell and piece arrays, `node_in`, the fluid
!< edge intervals; reals by their bits) is compared against `curve%cut_cells`, with and without the intervals, on
!< stress curves and several blocks (uniform, stretched, coarse, partly outside the curve):
!< - an ellipse, a 7-point star, a self-intersecting Lissajous curve (cells not closed included);
!< - grid ties: a diamond with vertices on grid nodes, a square on grid lines, a triangle with vertices on grid lines;
!< - an open polyline next to a closed loop;
!< - nested loops (solid square, hole, island, hole in the island) and loops inside single cells;
!< - a comb of thin teeth on coarse blocks (many pieces per cell);
!< - one cell holding more than 512 wall pieces in loops (a 700-gon with 100 small holes inside);
!< - no cut cell at all: an empty curve, a curve outside the blocks, the blocks inside the solid.
!< The integer stages of the device path (issue #100 step 6a: `cut_cells_2d_classify_dev`, `cut_cells_2d_scan_dev`,
!< two-level chunk scans) are compared with the host's `cut_cells_2d_classify` / `_scan` integer for integer, on chunks of
!< 1, 3, 7 cells and the default, with and without moments, on synthetic counts for the scan (`stage_check`).
use fossil,                   only : curve_object
use fossil_cut_cells_compare, only : cut_cells_twin_mismatches
use fossil_cut_cells_2d,      only : node_flags, cut_cells_2d_prepare, cut_cells_2d_classify, cut_cells_2d_scan, &
                                     cut_cells_2d_classify_dev, cut_cells_2d_scan_dev
use penf,                     only : I4P, R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P)
type(curve_object)      :: c
real(R8P)               :: e8(9), s8(9), g3(4), gx(6), gy(4), go(12), gf(38), gh(30)
real(R8P), allocatable  :: u(:), v(:)
integer(I4P)            :: k, l, m, n, st, total, nv
logical                 :: ok

do k = 1, 9
   e8(k) = (k - 1) / 8._R8P ; s8(k) = ((k - 1) / 8._R8P)**1.4_R8P
enddo
g3 = [0._R8P, 1._R8P / 3._R8P, 2._R8P / 3._R8P, 1._R8P]
gx = [-0.1_R8P, 0.2_R8P, 0.45_R8P, 0.5_R8P, 0.8_R8P, 1.1_R8P] ; gy = [0.02_R8P, 0.4_R8P, 0.61_R8P, 0.97_R8P]
do k = 1, 12
   go(k) = 0.3_R8P + (k - 1) * 0.1_R8P ! partly outside the curves
enddo
do k = 1, 38
   gf(k) = -0.05_R8P + (k - 1) * 1.08_R8P / 37._R8P
enddo
do k = 1, 30
   gh(k) = -0.02_R8P + ((k - 1) / 29._R8P)**1.2_R8P * 1.05_R8P
enddo
ok = .true. ; total = 0_I4P
! 1. ellipse
nv = 256 ; allocate(u(nv), v(nv))
do k = 1, nv
   u(k) = 0.5_R8P + 0.37_R8P * cos(2 * PI * (k - 1) / nv + 0.1_R8P)
   v(k) = 0.47_R8P + 0.23_R8P * sin(2 * PI * (k - 1) / nv + 0.1_R8P)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[nv])
call check('ellipse')
call stage_check('ellipse', u, v, [nv])
! 2. star
deallocate(u, v) ; nv = 14 ; allocate(u(nv), v(nv))
do k = 1, nv
   u(k) = 0.52_R8P + merge(0.4_R8P, 0.17_R8P, mod(k, 2) == 1) * cos(2 * PI * (k - 1) / nv)
   v(k) = 0.49_R8P + merge(0.4_R8P, 0.17_R8P, mod(k, 2) == 1) * sin(2 * PI * (k - 1) / nv)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[nv])
call check('star')
! 3. self-intersecting Lissajous (3:2)
deallocate(u, v) ; nv = 300 ; allocate(u(nv), v(nv))
do k = 1, nv
   u(k) = 0.5_R8P + 0.4_R8P * sin(3 * 2 * PI * (k - 1) / nv + 0.3_R8P) ; v(k) = 0.5_R8P + 0.35_R8P * sin(2 * 2 * PI * (k - 1) / nv)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[nv])
call check('lissajous')
! 4. grid ties on the 1/8 lattice: vertices on nodes, sides on lines, vertices on lines
call c%adopt_loops(u=[0.5_R8P, 0.75_R8P, 0.5_R8P, 0.25_R8P], v=[0.25_R8P, 0.5_R8P, 0.75_R8P, 0.5_R8P], loop_sizes=[4_I4P])
call check('diamond on grid nodes')
call c%adopt_loops(u=[0.25_R8P, 0.625_R8P, 0.625_R8P, 0.25_R8P], v=[0.25_R8P, 0.25_R8P, 0.875_R8P, 0.875_R8P], &
                   loop_sizes=[4_I4P])
call check('square on grid lines')
call c%adopt_loops(u=[0.375_R8P, 0.9_R8P, 0.375_R8P], v=[0.1_R8P, 0.375_R8P, 0.8_R8P], loop_sizes=[3_I4P])
call check('triangle with vertices on lines')
! 5. an open polyline next to a closed loop
call c%adopt_loops(u=[0.1_R8P, 0.4_R8P, 0.4_R8P, 0.1_R8P, 0.55_R8P, 0.7_R8P, 0.62_R8P, 0.93_R8P], &
                   v=[0.1_R8P, 0.1_R8P, 0.42_R8P, 0.42_R8P, 0.05_R8P, 0.5_R8P, 0.7_R8P, 0.88_R8P], &
                   loop_sizes=[4_I4P, 4_I4P], closed=[.true., .false.])
call check('open + closed')
! 6. nested: solid square, hole, island, hole in the island; an island and a hole inside single cells
deallocate(u, v) ; allocate(u(24), v(24))
u(1:16) = [0.05_R8P, 0.95_R8P, 0.95_R8P, 0.05_R8P, 0.2_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, &
           0.35_R8P, 0.65_R8P, 0.65_R8P, 0.35_R8P, 0.45_R8P, 0.45_R8P, 0.55_R8P, 0.55_R8P]
v(1:16) = [0.05_R8P, 0.05_R8P, 0.95_R8P, 0.95_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P, &
           0.35_R8P, 0.35_R8P, 0.65_R8P, 0.65_R8P, 0.45_R8P, 0.55_R8P, 0.55_R8P, 0.45_R8P]
u(17:24) = [0.27_R8P, 0.29_R8P, 0.29_R8P, 0.27_R8P, 0.07_R8P, 0.07_R8P, 0.09_R8P, 0.09_R8P] ! island in the hole,
v(17:24) = [0.27_R8P, 0.27_R8P, 0.29_R8P, 0.29_R8P, 0.1_R8P, 0.12_R8P, 0.12_R8P, 0.1_R8P]   ! hole in the solid rim
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P])
call check('nested loops')
call stage_check('nested loops', u, v, [4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P])
! 7. a comb of 12 thin teeth (many pieces per cell on coarse blocks)
deallocate(u, v) ; allocate(u(4 * 12 + 2), v(4 * 12 + 2))
n = 0
do k = 1, 12
   n = n + 1 ; u(n) = 0.06_R8P + (k - 1) * 0.075_R8P ; v(n) = 0.1_R8P
   n = n + 1 ; u(n) = u(n - 1) + 0.03_R8P ; v(n) = 0.1_R8P
   n = n + 1 ; u(n) = u(n - 1) ; v(n) = 0.93_R8P
   n = n + 1 ; u(n) = u(n - 1) + 0.045_R8P ; v(n) = 0.93_R8P
enddo
u(n) = 0.95_R8P ; v(n) = 0.93_R8P ; n = n + 1 ; u(n) = 0.95_R8P ; v(n) = 0.03_R8P ; n = n + 1 ; u(n) = 0.06_R8P
v(n) = 0.03_R8P
call c%adopt_loops(u=u(1:n), v=v(1:n), loop_sizes=[n])
call check('comb')
call stage_check('comb', u(1:n), v(1:n), [n])
! 8. one cell of g3 ([1/3, 2/3]^2) with > 512 wall pieces: a 700-gon (solid) with 100 small holes
deallocate(u, v) ; allocate(u(700 + 400), v(700 + 400))
do k = 1, 700
   u(k) = 0.5_R8P + 0.15_R8P * cos(2 * PI * (k - 1) / 700) ; v(k) = 0.5_R8P + 0.15_R8P * sin(2 * PI * (k - 1) / 700)
enddo
n = 700
do l = 1, 10
   do m = 1, 10
      u(n + 1:n + 4) = 0.42_R8P + (l - 1) * 0.016_R8P + [0._R8P, 0._R8P, 0.006_R8P, 0.006_R8P]
      v(n + 1:n + 4) = 0.42_R8P + (m - 1) * 0.016_R8P + [0._R8P, 0.006_R8P, 0.006_R8P, 0._R8P]
      n = n + 4
   enddo
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[700_I4P, [(4_I4P, k = 1, 100)]])
call stage_check('700-gon with 100 holes', u, v, [700_I4P, [(4_I4P, k = 1, 100)]])
k = cut_cells_twin_mismatches(c, g3, g3, st)
total = total + k ; ok = ok .and. k == 0_I4P .and. st == 0_I4P
print '(A,I0)', '700-gon with 100 holes in one cell, 3 x 3: mismatches ', k
k = cut_cells_twin_mismatches(c, e8, e8, st)
total = total + k ; ok = ok .and. k == 0_I4P .and. st == 0_I4P
print '(A,I0)', '700-gon with 100 holes, 8 x 8: mismatches ', k
! 9. no cut cell: an empty curve, a curve outside the blocks, every block inside the solid (empty arrays on the device)
call c%destroy
call check('empty curve')
call c%adopt_loops(u=[2._R8P, 3._R8P, 2.5_R8P], v=[2._R8P, 2._R8P, 3._R8P], loop_sizes=[3_I4P])
call check('curve outside the blocks')
call c%adopt_loops(u=[-1._R8P, 3._R8P, 3._R8P, -1._R8P], v=[-1._R8P, -1._R8P, 3._R8P, 3._R8P], loop_sizes=[4_I4P])
call check('blocks inside the solid')
! 10. bad edges: both refused
k = cut_cells_twin_mismatches(c, [0._R8P, 0._R8P, 1._R8P], e8, st)
ok = ok .and. k == 0_I4P .and. st /= 0_I4P
print '(A,L1)', 'bad edges refused by both: ', st /= 0_I4P
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label)
   !< The twin against the host path on every block.
   character(*), intent(in) :: label !< Curve.
   integer(I4P)             :: nm(6), s(6) !< Mismatches, statuses per block.

   nm(1) = cut_cells_twin_mismatches(c, e8, e8, s(1))
   nm(2) = cut_cells_twin_mismatches(c, s8, e8, s(2))
   nm(3) = cut_cells_twin_mismatches(c, g3, g3, s(3))
   nm(4) = cut_cells_twin_mismatches(c, gx, gy, s(4))
   nm(5) = cut_cells_twin_mismatches(c, go, gf, s(5))
   nm(6) = cut_cells_twin_mismatches(c, gf, gh, s(6))
   total = total + sum(nm) ; ok = ok .and. all(nm == 0_I4P) .and. all(s == 0_I4P)
   print '(A,6(1X,I0))', label//': mismatches on 6 blocks', nm
   endsubroutine check

   subroutine stage_check(label, uu, vv, sizes)
   !< The device-path classification and prefix sums against the host's, integer for integer, on the six blocks of
   !< `check`, closed loops of the given sizes (issue #100 step 6a).
   character(*), intent(in) :: label            !< Curve.
   real(R8P),    intent(in) :: uu(:), vv(:)     !< Vertices.
   integer(I4P), intent(in) :: sizes(:)         !< Loop sizes (closed loops).
   integer(I4P), allocatable :: nxt(:)          !< Next vertex.
   integer(I4P) :: b, f, l, nb                  !< Counters.

   allocate(nxt(size(uu)))
   f = 0_I4P
   do l = 1, size(sizes)
      nxt(f + 1:f + sizes(l)) = [(f + b + 1_I4P, b = 1, sizes(l) - 1), f + 1_I4P]
      f = f + sizes(l)
   enddo
   nb = stage_mismatches(uu, vv, nxt, e8, e8) + stage_mismatches(uu, vv, nxt, s8, e8) + &
        stage_mismatches(uu, vv, nxt, g3, g3) + stage_mismatches(uu, vv, nxt, gx, gy) + &
        stage_mismatches(uu, vv, nxt, go, gf) + stage_mismatches(uu, vv, nxt, gf, gh)
   total = total + nb ; ok = ok .and. nb == 0_I4P
   print '(A,I0)', label//': device-path classification and prefix sums vs host, mismatches ', nb
   endsubroutine stage_check

   function stage_mismatches(uu, vv, nxt, x, y) result(nb)
   !< Mismatches of the device-path integer stages against the host's on the block x, y.
   real(R8P),    intent(in) :: uu(:), vv(:), x(:), y(:) !< Vertices, cell edges.
   integer(I4P), intent(in) :: nxt(:)                   !< Next vertex.
   integer(I4P)             :: nb                       !< Mismatches.
   logical,      allocatable :: node_in(:,:)
   integer(I4P), allocatable :: xstart(:), xord(:), pstart(:), pord(:), ptg(:,:)
   real(R8P),    allocatable :: xp(:,:), pab(:,:)
   integer(I4P), allocatable :: ck1(:), cc1(:), rw1(:), iw1(:), ck2(:), cc2(:), rw2(:), iw2(:), cnt(:,:)
   integer(I4P), allocatable :: fp1(:,:), iv1(:), bd1(:), fp2(:,:), iv2(:), bd2(:)
   integer(I4P) :: nx, ny, ne, ncell, nxs, nps, nc1, nc2, nrt, nit, np1, ni1, nd1, np2, ni2, nd2, q, sg, nm, w
   integer(I4P), parameter :: SEGS(4) = [1_I4P, 3_I4P, 7_I4P, 0_I4P] !< Chunk sizes (0: the default).

   nb = 0_I4P
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   call node_flags(uu, vv, nxt, x, y, node_in)
   call cut_cells_2d_prepare(uu, vv, nxt, x, y, nx, ny, ne, ncell, nxs, nps, xstart, xord, pstart, pord, xp, pab, ptg)
   do nm = 0, 3, 3
      call cut_cells_2d_classify(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ck1, nc1, cc1, rw1, iw1)
      allocate(cnt(3, max(1_I4P, nc1)))
      do q = 1, nc1 ! synthetic counts: the scan is a function of the kinds and the counts only
         cnt(:, q) = [mod(7_I4P * q, 5_I4P), mod(3_I4P * q, 4_I4P), 1_I4P + mod(q, 3_I4P)]
      enddo
      do sg = 1, size(SEGS)
         if (SEGS(sg) > 0_I4P) then
            call cut_cells_2d_classify_dev(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ck2, nc2, cc2, rw2, iw2, nrt, nit, &
                                           seg=SEGS(sg))
         else
            call cut_cells_2d_classify_dev(nx, ny, ne, ncell, xstart, pstart, node_in, nm, ck2, nc2, cc2, rw2, iw2, nrt, nit)
         endif
         if (nc1 /= nc2) then
            nb = nb + 1_I4P ; cycle
         endif
         nb = nb + count(ck1 /= ck2) + count(cc1 /= cc2) + count(rw1 /= rw2) + count(iw1 /= iw2)
         if (nrt /= rw1(nc1 + 1) - 1_I4P .or. nit /= iw1(nc1 + 1) - 1_I4P) nb = nb + 1_I4P
         do w = 0, 3
            allocate(fp1(nx, ny), iv1(ncell), bd1(ncell), fp2(nx, ny), iv2(ncell), bd2(ncell))
            call cut_cells_2d_scan(nx, ny, ncell, ck1, nc1, cnt, mod(w, 2) == 1, w >= 2, fp1, iv1, bd1, np1, ni1, nd1)
            if (SEGS(sg) > 0_I4P) then
               call cut_cells_2d_scan_dev(nx, ny, ncell, ck2, nc2, cnt, mod(w, 2) == 1, w >= 2, fp2, iv2, bd2, np2, ni2, &
                                          nd2, seg=SEGS(sg))
            else
               call cut_cells_2d_scan_dev(nx, ny, ncell, ck2, nc2, cnt, mod(w, 2) == 1, w >= 2, fp2, iv2, bd2, np2, ni2, nd2)
            endif
            nb = nb + count(fp1 /= fp2) + count(iv1 /= iv2) + count(bd1 /= bd2)
            if (np1 /= np2 .or. ni1 /= ni2 .or. nd1 /= nd2) nb = nb + 1_I4P
            deallocate(fp1, iv1, bd1, fp2, iv2, bd2)
         enddo
      enddo
      deallocate(cnt)
   enddo
   endfunction stage_mismatches
endprogram fossil_test_cut_cells_device_twin
