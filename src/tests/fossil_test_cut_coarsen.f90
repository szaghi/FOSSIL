!< FOSSIL test: exact multilevel coarsening of the cut cells (issue #70, step 10).
program fossil_test_cut_coarsen
!< FOSSIL test: exact multilevel coarsening of the cut cells (issue #70, step 10).
!<
!< `cut_coarsen_3d` / `_2d` at ratio 2 against the cut cells computed directly on the coarse grid (every other edge of
!< the fine one). 3-D: the 5-degree wedge tilted off every axis, the L-prism, the thin plate; 2-D: a circle of 64
!< segments and a thin plate.
!< 1. the same number of pieces in every coarse cell, the same untrusted cells, and the same piece volumes (sorted per
!<    cell, 1e-13 of the cell);
!< 2. alpha, the face fractions (apertures in 2-D) and the wetted areas (lengths) agree (1e-13);
!< 3. node flags agree exactly, edge fractions to 1e-14;
!< 4. moments up to order 2 of each coarse cell (sum over its pieces) agree (1e-12 of the cell);
!< 5. the links across every coarse face add up to the same area (1e-13 of the face).

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, &
                                     pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_cut_coarsen,       only : cut_coarsen_3d, cut_coarsen_2d
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
real(R8P), allocatable    :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                 :: err(5)
integer(I4P)              :: i, n, nbad, ncells
logical                   :: ok

err = 0._R8P ; nbad = 0 ; ncells = 0
call bench_bodies(b)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 26)]
call check3(x, y, z)
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.185_R8P * i, i = 0, 14)] ; y = x ; z = [-0.2_R8P, 0.125_R8P, 0.45_R8P, 0.875_R8P, 1.3_R8P]
call check3(x, y, z)
call to_surface(b(5)%soup, s) ! plate
x = [(-0.15_R8P + 0.05_R8P * i, i = 0, 26)] ; z = [-0.3_R8P, -0.2_R8P, -0.1_R8P, 0.00025_R8P, 0.1_R8P, 0.2_R8P, 0.3_R8P]
call check3(x, x, z)
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.035_R8P * (i - 1), i = 1, 29)] ; y = [(-0.021_R8P + 0.0355_R8P * (i - 1), i = 1, 31)]
call check2(x, y)
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P], v=[0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P], loop_sizes=[4_I4P])
x = [(0.05_R8P + 0.05_R8P * (i - 1), i = 1, 19)]
call check2(x, x)
print '(A,I0,A,I0)', 'coarse cells ', ncells, '  mismatches ', nbad
print '(A,5ES10.2)', 'max errors (volumes, fractions, edges, moments, links): ', err
ok = nbad == 0 .and. ncells > 1000 .and. err(1) <= 1.e-13_R8P .and. err(2) <= 1.e-13_R8P .and. err(3) <= 1.e-14_R8P .and. &
     err(4) <= 1.e-12_R8P .and. err(5) <= 1.e-13_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check3(xf, yf, zf)
   !< Checks 1-5 on the fine block xf x yf x zf and its coarse block.
   real(R8P), intent(in)     :: xf(:), yf(:), zf(:)
   type(cut_cells_3d_t)      :: fine, direct, coarse
   type(cut_links_t)         :: lf, ld, lc
   integer(I4P), allocatable :: f2c(:)
   real(R8P), allocatable    :: xc(:), yc(:), zc(:)
   integer(I4P)              :: st, ia, ja, ka, p, np, k0
   real(R8P)                 :: vc, h(3), md(10), mc(10)

   call cut_cells_3d(s, xf, yf, zf, fine, links=lf, edges=.true., moments_order=2_I4P)
   xc = xf(1::2) ; yc = yf(1::2) ; zc = zf(1::2)
   call cut_cells_3d(s, xc, yc, zc, direct, links=ld, edges=.true., moments_order=2_I4P)
   call cut_coarsen_3d(fine, lf, xf, yf, zf, 2_I4P, coarse, lc, f2c, st)
   if (st /= 0) then
      nbad = nbad + 1 ; return
   endif
   if (any(coarse%node_in .neqv. direct%node_in)) nbad = nbad + 1
   err(3) = max(err(3), maxval(abs(coarse%edge_x - direct%edge_x)), maxval(abs(coarse%edge_y - direct%edge_y)), &
                maxval(abs(coarse%edge_z - direct%edge_z)))
   err(2) = max(err(2), maxval(abs(coarse%face_x - direct%face_x)), maxval(abs(coarse%face_y - direct%face_y)), &
                maxval(abs(coarse%face_z - direct%face_z)))
   do ka = 1, coarse%nz
      do ja = 1, coarse%ny
         do ia = 1, coarse%nx
            ncells = ncells + 1
            np = coarse%pieces(ia, ja, ka)
            if (np /= direct%pieces(ia, ja, ka)) then
               nbad = nbad + 1 ; cycle
            endif
            if (np < 0) cycle
            h = [xc(ia + 1) - xc(ia), yc(ja + 1) - yc(ja), zc(ka + 1) - zc(ka)] ; vc = product(h)
            p = coarse%first_piece(ia, ja, ka) ; k0 = direct%first_piece(ia, ja, ka)
            err(1) = max(err(1), maxval(abs(sorted(coarse%piece_volume(p:p + np - 1)) - &
                                             sorted(direct%piece_volume(k0:k0 + np - 1)))) / vc)
            err(2) = max(err(2), abs(coarse%alpha(ia, ja, ka) - direct%alpha(ia, ja, ka)), &
                         abs(coarse%wet_area(ia, ja, ka) - direct%wet_area(ia, ja, ka)) / maxval(h) ** 2)
            md = 0._R8P ; mc = 0._R8P
            do p = direct%first_piece(ia, ja, ka), direct%first_piece(ia, ja, ka) + np - 1
               md = md + direct%piece_moments(:, p)
            enddo
            do p = coarse%first_piece(ia, ja, ka), coarse%first_piece(ia, ja, ka) + np - 1
               mc = mc + coarse%piece_moments(:, p)
            enddo
            err(4) = max(err(4), maxval(abs(md - mc)) / vc)
            err(5) = max(err(5), link_face_err3(coarse, lc, direct, ld, [ia, ja, ka], h))
         enddo
      enddo
   enddo
   endsubroutine check3

   real(R8P) function link_face_err3(cg, lg, dg, ldg, ic, h) result(e)
   !< Largest difference, over the six faces of cell ic, of the summed link areas (relative to the face area).
   type(cut_cells_3d_t), intent(in) :: cg, dg
   type(cut_links_t),    intent(in) :: lg, ldg
   integer(I4P),         intent(in) :: ic(3)
   real(R8P),            intent(in) :: h(3)
   real(R8P)                        :: sc(6), sd(6), af(6)
   integer(I4P)                     :: p, k

   sc = 0._R8P ; sd = 0._R8P
   do p = cg%first_piece(ic(1), ic(2), ic(3)), cg%first_piece(ic(1), ic(2), ic(3)) + cg%pieces(ic(1), ic(2), ic(3)) - 1
      do k = lg%start(p), lg%start(p + 1) - 1
         sc(lg%face(k)) = sc(lg%face(k)) + lg%area(k)
      enddo
   enddo
   do p = dg%first_piece(ic(1), ic(2), ic(3)), dg%first_piece(ic(1), ic(2), ic(3)) + dg%pieces(ic(1), ic(2), ic(3)) - 1
      do k = ldg%start(p), ldg%start(p + 1) - 1
         sd(ldg%face(k)) = sd(ldg%face(k)) + ldg%area(k)
      enddo
   enddo
   af = [h(2) * h(3), h(2) * h(3), h(3) * h(1), h(3) * h(1), h(1) * h(2), h(1) * h(2)]
   e = maxval(abs(sc - sd) / af)
   endfunction link_face_err3

   subroutine check2(xf, yf)
   !< Checks 1-5 on the 2-D fine block xf x yf and its coarse block.
   real(R8P), intent(in)     :: xf(:), yf(:)
   type(cut_cells_2d_t)      :: fine, direct, coarse
   type(cut_links_t)         :: lf, ld, lc
   integer(I4P), allocatable :: f2c(:)
   real(R8P), allocatable    :: xc(:), yc(:)
   integer(I4P)              :: st, ia, ja, p, np, k
   real(R8P)                 :: ac, h(2), md(6), mc(6), sc(4), sd(4)

   call c%cut_cells(x=xf, y=yf, cut=fine, moments_order=2_I4P)
   call cut_links_2d(fine, lf)
   xc = xf(1::2) ; yc = yf(1::2)
   call c%cut_cells(x=xc, y=yc, cut=direct, moments_order=2_I4P)
   call cut_links_2d(direct, ld)
   call cut_coarsen_2d(fine, lf, xf, yf, 2_I4P, coarse, lc, f2c, st)
   if (st /= 0) then
      nbad = nbad + 1 ; return
   endif
   if (any(coarse%node_in .neqv. direct%node_in)) nbad = nbad + 1
   do ja = 1, coarse%ny
      do ia = 1, coarse%nx
         ncells = ncells + 1
         np = coarse%pieces(ia, ja)
         if (np /= direct%pieces(ia, ja)) then
            nbad = nbad + 1 ; cycle
         endif
         if (.not. direct%closed(ia, ja)) cycle
         h = [xc(ia + 1) - xc(ia), yc(ja + 1) - yc(ja)] ; ac = product(h)
         err(1) = max(err(1), maxval(abs(sorted(coarse%piece_area(coarse%first_piece(ia, ja):coarse%first_piece(ia, ja) + np - 1)) &
                  - sorted(direct%piece_area(direct%first_piece(ia, ja):direct%first_piece(ia, ja) + np - 1)))) / ac)
         err(2) = max(err(2), abs(coarse%alpha(ia, ja) - direct%alpha(ia, ja)), &
                      maxval(abs(coarse%aperture(:, ia, ja) - direct%aperture(:, ia, ja))))
         md = 0._R8P ; mc = 0._R8P ; sc = 0._R8P ; sd = 0._R8P
         do p = direct%first_piece(ia, ja), direct%first_piece(ia, ja) + np - 1
            md = md + direct%piece_moments(:, p)
            do k = ld%start(p), ld%start(p + 1) - 1
               sd(ld%face(k)) = sd(ld%face(k)) + ld%area(k)
            enddo
         enddo
         do p = coarse%first_piece(ia, ja), coarse%first_piece(ia, ja) + np - 1
            mc = mc + coarse%piece_moments(:, p)
            do k = lc%start(p), lc%start(p + 1) - 1
               sc(lc%face(k)) = sc(lc%face(k)) + lc%area(k)
            enddo
         enddo
         err(4) = max(err(4), maxval(abs(md - mc)) / ac)
         err(5) = max(err(5), maxval(abs(sc - sd) / [h(2), h(2), h(1), h(1)]))
      enddo
   enddo
   endsubroutine check2

   pure function sorted(a) result(r)
   !< a sorted ascending (insertion).
   real(R8P), intent(in) :: a(:)
   real(R8P)             :: r(size(a)), t
   integer(I4P)          :: i1, j1

   r = a
   do i1 = 2, size(r)
      t = r(i1) ; j1 = i1 - 1
      do while (j1 >= 1)
         if (r(j1) <= t) exit
         r(j1 + 1) = r(j1) ; j1 = j1 - 1
      enddo
      r(j1 + 1) = t
   enddo
   endfunction sorted
endprogram fossil_test_cut_coarsen
