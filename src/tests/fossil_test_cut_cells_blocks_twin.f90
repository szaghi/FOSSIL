!< FOSSIL test: the multi-block 2-D cut cells and the device sections of the 3-D cut cells equal the host paths bit for
!< bit (issue #100 step 3b).
program fossil_test_cut_cells_blocks_twin
!< FOSSIL test: the multi-block 2-D cut cells and the device sections of the 3-D cut cells equal the host paths bit for
!< bit (issue #100 step 3b).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, offset and
!< concatenation bugs of `compute_cut_cells_2d_blocks_device` and of `cut_cells_3d(..., device_sections=.true.)`, not
!< device codegen (the distance bench does that on the GPU, `cut cells 3-D sections device vs serial mismatches`).
!< 1. Twelve stress curves as the blocks of ONE `cut_cells_device_many` call (empty blocks first, in the middle and
!<    last; an ellipse, a star, a self-intersecting Lissajous curve, grid ties, an open polyline, nested loops, a comb, a
!<    curve outside the grid, a square around it) on six grids, against `curve%cut_cells` per block
!<    (`cut_cells_many_mismatches`: with pieces and intervals, pieces only, cells only; copied back and resident).
!< 2. The plane sections of 3-D bodies (a tilted wedge, a cone, an L-prism on an aligned grid, lattice tetrahedra, and
!<    the unit cube `cube-01.stl` on a grid through its faces) as blocks, the same check.
!< 3. On the same bodies, `cut_cells_3d(..., device_sections=.true.)` and `cut_cells_3d_sections(..., device=.true.)`
!<    against the host path: the face arrays bitwise, and the cell fractions and pieces that follow from them.
!< 4. Per-block grids (issue #100 step 4.1, `cut_cells_device_many(curves, x(:, :), y(:, :))`): the stress curves each
!<    on its own shifted grid (shifts by multiples of the cell size put the curves' vertices and sides on the shifted
!<    lines: ties), on stretched and uniform grids of one size mixed, and the translated planes of the swept geometry
!<    (sections of the cone and of cube-01 at the levels a - s d and grids ue - s d of two-point Gauss nodes, as
!<    `plane_average`), against `curve%cut_cells` on each block's grid (`cut_cells_many_grids_mismatches`); a grid count
!<    other than the block count and a bad column give STATUS_INVALID_INPUT.
use fossil,                   only : curve_object, surface_stl_object, slice_surface, cut_cells_3d_t, cut_cells_3d, &
                                     pose_from_axis_angle, facet_object
use fossil_cut_cells_3d,      only : cut_cells_3d_sections
use fossil_cut_cells_compare, only : cut_cells_many_mismatches, cut_cells_3d_face_mismatches, cut_cells_many_grids_mismatches
use fossil_curve_object,      only : cut_cells_2d_blocks_t, cut_cells_device_many
use fossil_surface_stl_object, only : STATUS_INVALID_INPUT
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P),    parameter        :: PI = acos(-1._R8P)
type(curve_object)             :: cv(15)
type(body_t),     allocatable  :: bb(:)
type(surface_stl_object)       :: s
real(R8P)                      :: e8(9), s8(9), g3(4), gx(6), gy(4), go(12), gf(38), gh(30), v(3, 4), vol
real(R8P),        allocatable  :: u(:), w(:), x(:), y(:), z(:), sh(:,:)
integer(I4P)                   :: k, n, nv, total, i
integer(I8P)                   :: seed
logical                        :: ok

do k = 1, 9
   e8(k) = (k - 1) / 8._R8P ; s8(k) = ((k - 1) / 8._R8P)**1.4_R8P
enddo
g3 = [0._R8P, 1._R8P / 3._R8P, 2._R8P / 3._R8P, 1._R8P]
gx = [-0.1_R8P, 0.2_R8P, 0.45_R8P, 0.5_R8P, 0.8_R8P, 1.1_R8P] ; gy = [0.02_R8P, 0.4_R8P, 0.61_R8P, 0.97_R8P]
do k = 1, 12
   go(k) = 0.3_R8P + (k - 1) * 0.1_R8P
enddo
do k = 1, 38
   gf(k) = -0.05_R8P + (k - 1) * 1.08_R8P / 37._R8P
enddo
do k = 1, 30
   gh(k) = -0.02_R8P + ((k - 1) / 29._R8P)**1.2_R8P * 1.05_R8P
enddo
ok = .true. ; total = 0_I4P
! 1. stress curves as blocks; cv(1), cv(8), cv(15) stay empty
nv = 256 ; allocate(u(nv), w(nv))
do k = 1, nv
   u(k) = 0.5_R8P + 0.37_R8P * cos(2 * PI * (k - 1) / nv + 0.1_R8P)
   w(k) = 0.47_R8P + 0.23_R8P * sin(2 * PI * (k - 1) / nv + 0.1_R8P)
enddo
call cv(2)%adopt_loops(u=u, v=w, loop_sizes=[nv])
deallocate(u, w) ; nv = 14 ; allocate(u(nv), w(nv))
do k = 1, nv
   u(k) = 0.52_R8P + merge(0.4_R8P, 0.17_R8P, mod(k, 2) == 1) * cos(2 * PI * (k - 1) / nv)
   w(k) = 0.49_R8P + merge(0.4_R8P, 0.17_R8P, mod(k, 2) == 1) * sin(2 * PI * (k - 1) / nv)
enddo
call cv(3)%adopt_loops(u=u, v=w, loop_sizes=[nv])
deallocate(u, w) ; nv = 300 ; allocate(u(nv), w(nv))
do k = 1, nv
   u(k) = 0.5_R8P + 0.4_R8P * sin(3 * 2 * PI * (k - 1) / nv + 0.3_R8P) ; w(k) = 0.5_R8P + 0.35_R8P * sin(2 * 2 * PI * (k - 1) / nv)
enddo
call cv(4)%adopt_loops(u=u, v=w, loop_sizes=[nv])
call cv(5)%adopt_loops(u=[0.5_R8P, 0.75_R8P, 0.5_R8P, 0.25_R8P], v=[0.25_R8P, 0.5_R8P, 0.75_R8P, 0.5_R8P], loop_sizes=[4_I4P])
call cv(6)%adopt_loops(u=[0.25_R8P, 0.625_R8P, 0.625_R8P, 0.25_R8P], v=[0.25_R8P, 0.25_R8P, 0.875_R8P, 0.875_R8P], &
                       loop_sizes=[4_I4P])
call cv(7)%adopt_loops(u=[0.375_R8P, 0.9_R8P, 0.375_R8P], v=[0.1_R8P, 0.375_R8P, 0.8_R8P], loop_sizes=[3_I4P])
call cv(9)%adopt_loops(u=[0.1_R8P, 0.4_R8P, 0.4_R8P, 0.1_R8P, 0.55_R8P, 0.7_R8P, 0.62_R8P, 0.93_R8P], &
                       v=[0.1_R8P, 0.1_R8P, 0.42_R8P, 0.42_R8P, 0.05_R8P, 0.5_R8P, 0.7_R8P, 0.88_R8P], &
                       loop_sizes=[4_I4P, 4_I4P], closed=[.true., .false.])
deallocate(u, w) ; allocate(u(24), w(24))
u(1:16) = [0.05_R8P, 0.95_R8P, 0.95_R8P, 0.05_R8P, 0.2_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, &
           0.35_R8P, 0.65_R8P, 0.65_R8P, 0.35_R8P, 0.45_R8P, 0.45_R8P, 0.55_R8P, 0.55_R8P]
w(1:16) = [0.05_R8P, 0.05_R8P, 0.95_R8P, 0.95_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P, &
           0.35_R8P, 0.35_R8P, 0.65_R8P, 0.65_R8P, 0.45_R8P, 0.55_R8P, 0.55_R8P, 0.45_R8P]
u(17:24) = [0.27_R8P, 0.29_R8P, 0.29_R8P, 0.27_R8P, 0.07_R8P, 0.07_R8P, 0.09_R8P, 0.09_R8P]
w(17:24) = [0.27_R8P, 0.27_R8P, 0.29_R8P, 0.29_R8P, 0.1_R8P, 0.12_R8P, 0.12_R8P, 0.1_R8P]
call cv(10)%adopt_loops(u=u, v=w, loop_sizes=[4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P, 4_I4P])
deallocate(u, w) ; allocate(u(4 * 12 + 2), w(4 * 12 + 2))
n = 0
do k = 1, 12
   n = n + 1 ; u(n) = 0.06_R8P + (k - 1) * 0.075_R8P ; w(n) = 0.1_R8P
   n = n + 1 ; u(n) = u(n - 1) + 0.03_R8P ; w(n) = 0.1_R8P
   n = n + 1 ; u(n) = u(n - 1) ; w(n) = 0.93_R8P
   n = n + 1 ; u(n) = u(n - 1) + 0.045_R8P ; w(n) = 0.93_R8P
enddo
u(n) = 0.95_R8P ; w(n) = 0.93_R8P ; n = n + 1 ; u(n) = 0.95_R8P ; w(n) = 0.03_R8P ; n = n + 1 ; u(n) = 0.06_R8P
w(n) = 0.03_R8P
call cv(11)%adopt_loops(u=u(1:n), v=w(1:n), loop_sizes=[n])
call cv(12)%adopt_loops(u=[2._R8P, 3._R8P, 2.5_R8P], v=[2._R8P, 2._R8P, 3._R8P], loop_sizes=[3_I4P])
call cv(13)%adopt_loops(u=[-1._R8P, 3._R8P, 3._R8P, -1._R8P], v=[-1._R8P, -1._R8P, 3._R8P, 3._R8P], loop_sizes=[4_I4P])
call cv(14)%adopt_loops(u=[0.3_R8P, 0.7_R8P, 0.5_R8P], v=[0.3_R8P, 0.3_R8P, 0.7_R8P], loop_sizes=[3_I4P])
call blocks('stress curves e8 x e8', cv, e8, e8)
call blocks('stress curves s8 x e8', cv, s8, e8)
call blocks('stress curves g3 x g3', cv, g3, g3)
call blocks('stress curves gx x gy', cv, gx, gy)
call blocks('stress curves go x gf', cv, go, gf)
call blocks('stress curves gf x gh', cv, gf, gh)
call blocks('one empty block', cv(1:1), e8, e8)
! 2-3. bodies: sections as blocks, then the 3-D faces
call bench_bodies(bb)
call to_surface(bb(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call body('tilted wedge')
call to_surface(bb(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call body('cone')
call to_surface(bb(8)%soup, s) ! L-prism, aligned grid
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call body('L-prism, aligned grid')
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)] ; y = x ; z = x
do k = 1, 3
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross3(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(v)
   call body('lattice tetrahedron')
enddo
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)] ; y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)] ; z = x
call body('cube-01, grid through its faces')
! 4. per-block grids
allocate(sh(2, 15))
sh(1, :) = [0._R8P, 0.125_R8P, -0.25_R8P, 0.0625_R8P, 0.375_R8P, -0.125_R8P, 0.5_R8P, 0.03125_R8P, 0.25_R8P, &
            -0.0375_R8P, 0.125_R8P, -2._R8P, 0.75_R8P, 0.1_R8P, 0._R8P]
sh(2, :) = [0._R8P, -0.125_R8P, 0.25_R8P, 0.25_R8P, -0.0625_R8P, 0._R8P, 0.125_R8P, -0.5_R8P, 0.0125_R8P, &
            0.375_R8P, -0.25_R8P, -2._R8P, 0._R8P, 0.2_R8P, 0.0625_R8P]
call shifted('stress curves, e8 x e8 shifted per block', cv, e8, e8, sh)
call shifted('stress curves, gf x gh shifted per block', cv, gf, gh, 0.37_R8P * sh)
call shifted('one empty block, shifted', cv(1:1), e8, e8, sh(:, 2:2))
call mixed('stress curves, s8 and e8 grids mixed', cv, s8, e8)
call bad_grids
call to_surface(bb(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call translated('cone, translated planes', [0.31_R8P, -0.12_R8P, 0.43_R8P])
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)] ; y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)] ; z = x
call translated('cube-01, translated planes through its faces', [0.125_R8P, 0.3_R8P, -0.25_R8P])
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine blocks(label, cc, xx, yy)
   !< The multi-block cut cells of the curves cc against the per-curve host path on the grid xx x yy.
   character(*),       intent(in) :: label   !< Case.
   type(curve_object), intent(in) :: cc(:)   !< Curves.
   real(R8P),          intent(in) :: xx(:), yy(:) !< Cell edges.
   integer(I4P)                   :: m, st   !< Mismatches, status.

   m = cut_cells_many_mismatches(cc, xx, yy, st)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0)', label//': blocks ', size(cc), ', mismatches ', m
   endsubroutine blocks

   subroutine shifted(label, cc, xx, yy, shift)
   !< Block b of the curves cc on the grid (xx - shift(1, b)) x (yy - shift(2, b)), against the per-curve host path.
   character(*),       intent(in) :: label   !< Case.
   type(curve_object), intent(in) :: cc(:)   !< Curves.
   real(R8P),          intent(in) :: xx(:), yy(:) !< Cell edges of the unshifted grid.
   real(R8P),          intent(in) :: shift(:,:) !< Shift of each block's grid [2, size(cc)].
   real(R8P), allocatable         :: gu(:,:), gv(:,:) !< Grids.
   integer(I4P)                   :: m, st, b !< Mismatches, status, block.

   allocate(gu(size(xx), size(cc)), gv(size(yy), size(cc)))
   do b = 1, size(cc)
      gu(:, b) = xx - shift(1, b) ; gv(:, b) = yy - shift(2, b)
   enddo
   m = cut_cells_many_grids_mismatches(cc, gu, gv, st)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0)', label//': blocks ', size(cc), ', mismatches ', m
   endsubroutine shifted

   subroutine mixed(label, cc, xa, xb)
   !< Blocks on different grids of one size: block b on xa x xb if b is odd, else on xb x xa (the stretched one reversed).
   character(*),       intent(in) :: label   !< Case.
   type(curve_object), intent(in) :: cc(:)   !< Curves.
   real(R8P),          intent(in) :: xa(:), xb(:) !< Two grids of one size.
   real(R8P), allocatable         :: gu(:,:), gv(:,:) !< Grids.
   integer(I4P)                   :: m, st, b !< Mismatches, status, block.

   allocate(gu(size(xa), size(cc)), gv(size(xa), size(cc)))
   do b = 1, size(cc)
      if (mod(b, 2) == 1) then
         gu(:, b) = xa ; gv(:, b) = xb
      else
         gu(:, b) = xb ; gv(:, b) = 1._R8P - xa(size(xa):1:-1)
      endif
   enddo
   m = cut_cells_many_grids_mismatches(cc, gu, gv, st)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0)', label//': blocks ', size(cc), ', mismatches ', m
   endsubroutine mixed

   subroutine bad_grids
   !< STATUS_INVALID_INPUT for a grid count other than the block count and for a column that does not increase.
   type(cut_cells_2d_blocks_t) :: cb       !< Result.
   real(R8P)                   :: gu(9, 3), gv(9, 3) !< Grids.
   integer(I4P)                :: st, m, b !< Status, failures, block.

   do b = 1, 3
      gu(:, b) = e8 + 0.1_R8P * b ; gv(:, b) = e8
   enddo
   m = 0_I4P
   call cut_cells_device_many(cv(2:3), gu, gv, cb, st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   gu(5, 3) = gu(4, 3)
   call cut_cells_device_many(cv(2:4), gu, gv, cb, st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call cut_cells_device_many(cv(2:4), gu(:, [1, 2, 1]), gv, cb, st)
   if (st /= 0_I4P) m = m + 1_I4P
   total = total + m ; ok = ok .and. m == 0_I4P
   print '(A,I0)', 'per-block grids, bad input: failures ', m
   endsubroutine bad_grids

   subroutine translated(label, d)
   !< The translated planes of `plane_average` (swept geometry, exact translation by d): for every z plane z_k and two
   !< Gauss nodes s of [0, 1], the section of s at z_k - s d(3) on the grid (x - s d(1)) x (y - s d(2)), all as the
   !< blocks of one call.
   character(*), intent(in) :: label !< Case.
   real(R8P),    intent(in) :: d(3)  !< Displacement.
   real(R8P),    parameter  :: G = 0.5_R8P / sqrt(3._R8P) !< Two-point Gauss offset.
   type(curve_object), allocatable :: secs(:) !< Sections.
   real(R8P), allocatable   :: gu(:,:), gv(:,:) !< Grids.
   real(R8P)                :: sg   !< Node.
   integer(I4P)             :: m, st, k, q, b !< Work.

   allocate(secs(2 * size(z)), gu(size(x), 2 * size(z)), gv(size(y), 2 * size(z)))
   b = 0
   do k = 1, size(z)
      do q = -1, 1, 2
         b = b + 1 ; sg = 0.5_R8P + q * G
         call slice_surface(surface=s, normal_axis=3_I4P, level=z(k) - sg * d(3), curve=secs(b), on_plane_low=.true.)
         gu(:, b) = x - sg * d(1) ; gv(:, b) = y - sg * d(2)
      enddo
   enddo
   m = cut_cells_many_grids_mismatches(secs, gu, gv, st)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0)', label//': blocks ', size(secs), ', mismatches ', m
   endsubroutine translated

   subroutine body(label)
   !< The sections of s along each axis as blocks; the 3-D faces with the device sections against the host's.
   character(*), intent(in) :: label !< Case.
   type(curve_object), allocatable :: secs(:)  !< Sections.
   type(cut_cells_3d_t) :: ch, cd             !< Host, device sections.
   integer(I4P)         :: m, ax, b, st, st2  !< Work.

   m = 0_I4P
   do ax = 1, 3
      select case (ax)
      case (1)
         allocate(secs(size(x)))
         do b = 1, size(x)
            call slice_surface(surface=s, normal_axis=1_I4P, level=x(b), curve=secs(b), on_plane_low=.true.)
         enddo
         m = m + cut_cells_many_mismatches(secs, y, z, st)
      case (2)
         allocate(secs(size(y)))
         do b = 1, size(y)
            call slice_surface(surface=s, normal_axis=2_I4P, level=y(b), curve=secs(b), on_plane_low=.true.)
         enddo
         m = m + cut_cells_many_mismatches(secs, z, x, st)
      case default
         allocate(secs(size(z)))
         do b = 1, size(z)
            call slice_surface(surface=s, normal_axis=3_I4P, level=z(b), curve=secs(b), on_plane_low=.true.)
         enddo
         m = m + cut_cells_many_mismatches(secs, x, y, st)
      endselect
      if (st /= 0_I4P) m = m + 1_I4P
      deallocate(secs)
   enddo
   call cut_cells_3d(s, x, y, z, ch, status=st)
   call cut_cells_3d(s, x, y, z, cd, status=st2, device_sections=.true.)
   m = m + cut_cells_3d_face_mismatches(ch, cd) + abs(st - st2)
   m = m + count(transfer(ch%alpha, 0_I8P, size(ch%alpha)) /= transfer(cd%alpha, 0_I8P, size(cd%alpha)), kind=I4P)
   m = m + count(ch%pieces /= cd%pieces, kind=I4P)
   if (ch%np /= cd%np) m = m + 1_I4P
   call cut_cells_3d_sections(s, x, y, z, cd, device=.true., status=st2)
   m = m + cut_cells_3d_face_mismatches(ch, cd) + abs(st - st2)
   call cut_cells_3d_sections(s, x, y, z, cd, status=st2)
   m = m + cut_cells_3d_face_mismatches(ch, cd) + abs(st - st2)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0)', label//': sections as blocks and 3-D faces, mismatches ', m
   endsubroutine body

   subroutine tetra(vt)
   !< s = the tetrahedron vt (positively oriented).
   real(R8P), intent(in) :: vt(3, 4) !< Vertices.
   type(facet_object), allocatable :: f(:) !< Facets.
   integer(I4P), parameter :: T(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4]) !< Outward faces.
   integer(I4P)          :: q, r     !< Counters.

   allocate(f(4))
   do q = 1, 4
      do r = 1, 3
         f(q)%vertex(r) = vector_R8P(vt(1, T(r, q)), vt(2, T(r, q)), vt(3, T(r, q)))
      enddo
   enddo
   call s%adopt_facets(facets=f)
   endsubroutine tetra

   pure function cross3(a, b) result(c)
   !< a x b.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< Product.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   function rnd() result(r)
   !< Deterministic uniform in (0, 1) (Park-Miller, the generator of `fossil_test_cut_links_3d`).
   real(R8P) :: r !< Value.

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_cells_blocks_twin
