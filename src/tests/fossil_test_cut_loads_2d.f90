!< FOSSIL test: the loads of a curve from the wall segments of its 2-D cut cells (#109 WP2 step 05).
program fossil_test_cut_loads_2d
!< FOSSIL test: the loads of a curve from the wall segments of its 2-D cut cells (#109 WP2 step 05).
!<
!< Curves inside the grid: a 64-gon circle, a disc with a hole (clockwise), a NACA 0012 polygon (200 vertices) and the
!< same turned by 20 degrees, two discs. Per unit depth:
!<
!< 1. Buoyancy through the cut: p = p0 + g . x at the segments' midpoints, tau = 0: F = -g A and, with the gradient,
!<    M = A ((x_c - c) x (-g))_z, A and x_c from the curve's `mass_properties`, each within TOL eps times the sum of the
!<    absolute terms (S_F = sum |p_w| L_w + A |g|, S_M = sum |p_w| L_w |x_w - c| + sum |K_w g| L_w + A |g| |x_c - c|); the
!<    largest ratio printed. Without the gradient the moment error is printed, not checked.
!< 2. Uniform pressure: the vector areas sum to 0 (4 eps sum L_w).
!< 3. The walls of each piece sum to `piece_wet_normal` (8 eps sum L_w of the piece); the bit-equal pieces counted.
!< 4. FD / FV: `integrate_loads` over load samples of order 2 and `cut_loads_2d` over the walls, the same linear p: the
!<    same force and moment within TOL eps (S_FD + S_FV).
!< 5. Uniform shear: F = tau L, the length output = the perimeter (1e-14).
!< 6. Per patch and per body (three patches over two discs): the sums add up; each disc its own buoyancy.
!< 7. `pressure_per_piece` = the per-wall pressure spread by hand, bit for bit; the walls' labels.
!< 8. Refusals: a cut without boundaries, two pressures, sizes, per body without the curve.
use fossil, only : curve_object, cut_cells_2d_t, cut_walls_2d_t, cut_walls_2d, cut_loads_2d, cut_loads_linear_pressure, &
                   load_samples_t, load_samples_2d, integrate_loads, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter    :: PI = acos(-1._R8P), EPS = epsilon(1._R8P), TOL = 1._R8P
real(R8P), parameter    :: G(2) = [0.4_R8P, -0.9_R8P], P0 = 1.5_R8P, C(2) = [0.15_R8P, -0.25_R8P]
real(R8P), parameter    :: TAU(2) = [0.3_R8P, -0.6_R8P]
type(curve_object)      :: circle, holed, naca, turned, discs
type(cut_cells_2d_t)    :: cut
type(cut_walls_2d_t)    :: walls
real(R8P), allocatable  :: u(:), v(:), x(:), y(:), p(:), t(:,:), gr(:,:), ppp(:), pw(:), ppf(:,:), ppm(:), pbf(:,:), &
                           pbm(:), u2(:), v2(:)
real(R8P)               :: f(2), m, f2(2), m2, a, rf, rm, rm0, rfd, ria, sf, sm, ca, sa
integer(I4P)            :: i, w, st, nbit, npc, k
logical                 :: ok

ok = .true. ; rf = 0._R8P ; rm = 0._R8P ; rm0 = 0._R8P ; rfd = 0._R8P ; ria = 0._R8P ; nbit = 0_I4P ; npc = 0_I4P
u = [(0.5_R8P + 0.4_R8P * cos(2 * PI * i / 64), i = 0, 63)] ; v = [(0.5_R8P + 0.4_R8P * sin(2 * PI * i / 64), i = 0, 63)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
u = [u, [(0.5_R8P + 0.15_R8P * cos(-2 * PI * i / 32), i = 0, 31)]]
v = [v, [(0.6_R8P + 0.15_R8P * sin(-2 * PI * i / 32), i = 0, 31)]]
call holed%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 32_I4P])
call naca0012(100_I4P, u, v)
call naca%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ca = cos(20 * PI / 180) ; sa = sin(20 * PI / 180)
u2 = ca * u + sa * v ; v2 = -sa * u + ca * v                 ! nose up about the leading edge
call turned%adopt_loops(u=u2, v=v2, loop_sizes=[200_I4P])
u = [[(0.3_R8P + 0.25_R8P * cos(2 * PI * i / 48), i = 0, 47)], [(1.1_R8P + 0.2_R8P * cos(2 * PI * i / 40), i = 0, 39)]]
v = [[(0.4_R8P + 0.25_R8P * sin(2 * PI * i / 48), i = 0, 47)], [(0.2_R8P + 0.2_R8P * sin(2 * PI * i / 40), i = 0, 39)]]
call discs%adopt_loops(u=u, v=v, loop_sizes=[48_I4P, 40_I4P], &
                       patches=[(1_I4P, i = 1, 20), (2_I4P, i = 21, 48), (3_I4P, i = 1, 40)], &
                       patch_names=['one  ', 'two  ', 'three'])
! 1-4
call body(circle, 'circle 64-gon', 16_I4P)
call body(holed, 'disc with a hole', 20_I4P)
call body(naca, 'NACA 0012', 40_I4P)
call body(turned, 'NACA 0012 turned', 40_I4P)
call body(discs, 'two discs', 24_I4P)
print '(A,F6.2,A,F6.2,A,F6.2,A,ES9.2)', 'largest error / (eps S): force ', rf, ', moment ', rm, ', FD vs FV ', rfd, &
      ';  centroid-only moment error / S_M ', rm0
print '(A,I0,A,I0,A,F6.2,A)', 'piece identity: ', nbit, ' of ', npc, ' pieces bit for bit, largest error ', ria, &
      ' eps sum L'
! 5. uniform shear
call box_lines(naca, 40_I4P)
call naca%cut_cells(x, y, cut, status=st, boundaries=.true.) ; call cut_walls_2d(cut, walls)
allocate(t(2, walls%n))
do w = 1, walls%n
   t(:, w) = TAU
enddo
call cut_loads_2d(walls, f, m, status=st, shear=t, area=a)
call check('uniform shear: F = tau L (4 eps)', st == STATUS_OK .and. all(abs(f - TAU * a) <= 4._R8P * EPS * abs(TAU) * a))
call check('uniform shear: length = perimeter (1e-14)', abs(a - naca%get_perimeter()) <= 1.e-14_R8P * a)
! 6. per patch and per body
call box_lines(discs, 24_I4P)
call discs%cut_cells(x, y, cut, status=st, boundaries=.true.) ; call cut_walls_2d(cut, walls)
p = cut_loads_linear_pressure(walls%centroid, P0, G) ; gr = spread(G, 2, walls%n)
deallocate(t) ; allocate(t(2, walls%n))
do w = 1, walls%n
   t(:, w) = [sin(3 * walls%centroid(2, w)), 0.4_R8P * cos(2 * walls%centroid(1, w))]
enddo
call cut_loads_2d(walls, f, m, status=st, pressure=p, shear=t, gradient=gr, about=C, curve=discs, per_patch_force=ppf, &
                  per_patch_moment=ppm, per_body_force=pbf, per_body_moment=pbm, area=a)
call scales(0._R8P, [0._R8P, 0._R8P], sf, sm)
call check('per patch: three patches, two bodies', st == STATUS_OK .and. size(ppf, 2) == 3 .and. size(ppm) == 3 .and. &
           size(pbf, 2) == 2)
call check('per patch: add up to the totals (4 eps S)', all(abs(sum(ppf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           abs(sum(ppm) - m) <= 4._R8P * EPS * (sm + 4 * a))
call check('per body: add up to the totals (4 eps S)', all(abs(sum(pbf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           abs(sum(pbm) - m) <= 4._R8P * EPS * (sm + 4 * a))
call cut_loads_2d(walls, f, m, status=st, pressure=p, gradient=gr, about=C, curve=discs, per_body_force=pbf, &
                  per_body_moment=pbm)
call body_check(1_I4P, 0_I4P, 48_I4P) ; call body_check(2_I4P, 48_I4P, 40_I4P)
! 7. per piece, labels
allocate(ppp(cut%np))
do k = 1, cut%np
   ppp(k) = 1._R8P + 0.1_R8P * sin(real(k, R8P))
enddo
allocate(pw(walls%n))
do w = 1, walls%n
   pw(w) = ppp(walls%piece(w))
enddo
call cut_loads_2d(walls, f, m, status=st, pressure_per_piece=ppp, shear=t, about=C)
call cut_loads_2d(walls, f2, m2, pressure=pw, shear=t, about=C)
call check('per piece: = the pressure spread by hand, bit for bit', st == STATUS_OK .and. same(f, f2) .and. same([m], [m2]))
call check('labels: every wall has a piece and a curve segment', all(walls%piece >= 1) .and. all(walls%segment >= 1) .and. &
           all(walls%segment <= 88))
! 8. refusals
call discs%cut_cells(x, y, cut, status=st)
call cut_walls_2d(cut, walls, status=st)
call check('refused: a cut without boundaries', st == STATUS_INVALID_INPUT .and. walls%n == 0)
call discs%cut_cells(x, y, cut, status=st, boundaries=.true.) ; call cut_walls_2d(cut, walls)
pw = [(1._R8P, w = 1, walls%n)]
call cut_loads_2d(walls, f, m, status=st)
call check('refused: neither p nor tau', st == STATUS_INVALID_INPUT .and. all(f == 0._R8P) .and. m == 0._R8P)
call cut_loads_2d(walls, f, m, status=st, pressure=pw, pressure_per_piece=ppp)
call check('refused: two pressures', st == STATUS_INVALID_INPUT)
call cut_loads_2d(walls, f, m, status=st, pressure=pw(1:walls%n - 1))
call check('refused: p of the wrong size', st == STATUS_INVALID_INPUT)
call cut_loads_2d(walls, f, m, status=st, shear=spread(TAU, 2, walls%n), gradient=spread(G, 2, walls%n))
call check('refused: gradient without a pressure', st == STATUS_INVALID_INPUT)
call cut_loads_2d(walls, f, m, status=st, pressure=pw, per_body_moment=pbm)
call check('refused: per body without the curve, size 0', st == STATUS_INVALID_INPUT .and. size(pbm) == 0)
call cut_loads_2d(walls, f, m, status=st, pressure=pw, about=[0._R8P, 0._R8P, 0._R8P])
call check('refused: about of size 3', st == STATUS_INVALID_INPUT)

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine body(cur, name, n)
   !< Checks 1-4 on one curve, on n x n cells around it.
   type(curve_object), intent(in) :: cur  !< Curve.
   character(*),       intent(in) :: name !< Label.
   integer(I4P),       intent(in) :: n    !< Cells per axis.
   type(load_samples_t)   :: ls                         !< Samples.
   real(R8P), allocatable :: pl(:), one(:), ac(:,:)     !< Pressures, piece sums.
   real(R8P)    :: aa, x0(2), ff(2), mm, fr(2), mr, ssf, ssm, e1, e2, ffd(2), mfd, sfd, smd, sl, pq !< Work.
   integer(I4P) :: sst, q, pc                           !< Work.

   call cur%mass_properties(area=aa, center=x0, status=sst)
   call box_lines(cur, n)
   call cur%cut_cells(x, y, cut, status=sst, boundaries=.true.)
   call cut_walls_2d(cut, walls, status=sst)
   call check(name//': cut walls', sst == STATUS_OK .and. walls%n > 0)
   pl = cut_loads_linear_pressure(walls%centroid, P0, G) ; gr = spread(G, 2, walls%n)
   fr = -G * aa ; mr = aa * cross2(x0 - C, -G)
   call cut_loads_2d(walls, ff, mm, status=sst, pressure=pl, gradient=gr, about=C)
   call scales(aa, x0, ssf, ssm)
   e1 = maxval(abs(ff - fr)) / (EPS * ssf) ; e2 = abs(mm - mr) / (EPS * ssm)
   rf = max(rf, e1) ; rm = max(rm, e2)
   call check(name//': buoyancy F = -g A', sst == STATUS_OK .and. e1 <= TOL)
   call check(name//': buoyancy M = A (xc - c) x (-g), with the gradient', e2 <= TOL)
   call cut_loads_2d(walls, ff, mm, status=sst, pressure=pl, about=C)
   rm0 = max(rm0, abs(mm - mr) / ssm)
   ! 2. uniform pressure
   allocate(one(walls%n)) ; one = 1._R8P
   call cut_loads_2d(walls, ff, mm, status=sst, pressure=one, area=sl)
   call check(name//': uniform p, sum of the vector areas 0 (4 eps)', all(abs(ff) <= 4._R8P * EPS * sl))
   ! 3. the walls of each piece
   allocate(ac(3, cut%np)) ; ac = 0._R8P
   do q = 1, walls%n
      pc = walls%piece(q)
      ac(1:2, pc) = ac(1:2, pc) + walls%area(:, q) ; ac(3, pc) = ac(3, pc) + walls%length(q)
   enddo
   do pc = 1, cut%np
      if (ac(3, pc) <= 0._R8P) cycle
      npc = npc + 1
      if (same(ac(1:2, pc), cut%piece_wet_normal(:, pc))) nbit = nbit + 1
      ria = max(ria, maxval(abs(ac(1:2, pc) - cut%piece_wet_normal(:, pc))) / (EPS * ac(3, pc)))
   enddo
   call check(name//': walls of a piece = piece_wet_normal (8 eps sum L)', ria <= 8._R8P)
   ! 4. FD / FV
   call load_samples_2d(cur, samples=ls, order=2_I4P)
   call integrate_loads(ls, ffd, mfd, status=sst, pressure=cut_loads_linear_pressure(ls%points, P0, G), about=C)
   call cut_loads_2d(walls, ff, mm, status=sst, pressure=pl, gradient=gr, about=C)
   sfd = aa * norm2(G) ; smd = aa * norm2(G) * norm2(x0 - C)
   do q = 1, ls%n
      pq = abs(P0 + dot_product(G, ls%points(:, q)))
      sfd = sfd + ls%weight(q) * pq ; smd = smd + ls%weight(q) * pq * norm2(ls%points(:, q) - C)
   enddo
   e1 = maxval(abs(ff - ffd)) / (EPS * (ssf + sfd)) ; e2 = abs(mm - mfd) / (EPS * (ssm + smd))
   rfd = max(rfd, e1, e2)
   call check(name//': FD (integrate_loads) = FV (cut_loads), force and moment', e1 <= TOL .and. e2 <= TOL)
   endsubroutine body

   subroutine body_check(b, i0, nv)
   !< Disc b (vertices i0 + 1 ... i0 + nv): its own buoyancy.
   integer(I4P), intent(in) :: b, i0, nv  !< Body, first vertex - 1, vertices.
   type(curve_object) :: one              !< The disc alone.
   real(R8P)          :: aa, x0(2), ssf, ssm !< Area, centroid, scales.

   call one%adopt_loops(u=u(i0 + 1:i0 + nv), v=v(i0 + 1:i0 + nv), loop_sizes=[nv])
   call one%mass_properties(area=aa, center=x0)
   call scales(aa, x0, ssf, ssm)
   call check('per body: disc '//char(48 + b)//' buoyancy', all(abs(pbf(:, b) + G * aa) <= TOL * EPS * ssf) .and. &
              abs(pbm(b) - aa * cross2(x0 - C, -G)) <= TOL * EPS * ssm)
   endsubroutine body_check

   subroutine box_lines(cur, n)
   !< Lines of n cells around the curve's box, uneven margins.
   type(curve_object), intent(in) :: cur !< Curve.
   integer(I4P),       intent(in) :: n   !< Cells per axis.
   real(R8P)    :: l(2), h(2)            !< Corners.
   integer(I4P) :: q                     !< Counter.

   l = cur%get_bmin() - [0.137_R8P, 0.093_R8P] ; h = cur%get_bmax() + [0.119_R8P, 0.157_R8P]
   x = [(l(1) + (h(1) - l(1)) * q / n, q = 0, n)]
   y = [(l(2) + (h(2) - l(2)) * q / n, q = 0, n)]
   endsubroutine box_lines

   subroutine scales(aa, x0, ssf, ssm)
   !< The sums of the absolute terms of the walls (`walls`, p0 + g . x).
   real(R8P),    intent(in)  :: aa, x0(2) !< Area, centroid.
   real(R8P),    intent(out) :: ssf, ssm  !< Scales.
   real(R8P)    :: pq, d(2)               !< Work.
   integer(I4P) :: q                      !< Counter.

   ssf = aa * norm2(G) ; ssm = aa * norm2(G) * norm2(x0 - C)
   do q = 1, walls%n
      pq = abs(P0 + dot_product(G, walls%centroid(:, q))) ; d = walls%ends(3:4, q) - walls%ends(1:2, q)
      ssf = ssf + pq * walls%length(q)
      ssm = ssm + pq * walls%length(q) * norm2(walls%centroid(:, q) - C) + &
            norm2(d) * abs(dot_product(d, G)) / 12._R8P * walls%length(q)
   enddo
   endsubroutine scales

   subroutine naca0012(mh, xx, yy)
   !< NACA 0012 of unit chord, closed trailing edge, 2 mh vertices, cosine spacing, counter-clockwise.
   integer(I4P),           intent(in)  :: mh     !< Points per side.
   real(R8P), allocatable, intent(out) :: xx(:), yy(:) !< Vertices.
   real(R8P)    :: xi                            !< Abscissa.
   integer(I4P) :: q, j                          !< Counters.

   allocate(xx(2 * mh), yy(2 * mh)) ; j = 0
   do q = mh, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = yt(xi)
   enddo
   do q = 1, mh - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = -yt(xi)
   enddo
   endsubroutine naca0012

   pure real(R8P) function yt(xx)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: xx !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.1260_R8P * xx - 0.3516_R8P * xx ** 2 + 0.2843_R8P * xx ** 3 - &
                   0.1036_R8P * xx ** 4)
   endfunction yt

   pure real(R8P) function cross2(a1, b1)
   !< z component of a1 x b1.
   real(R8P), intent(in) :: a1(2), b1(2) !< Vectors.

   cross2 = a1(1) * b1(2) - a1(2) * b1(1)
   endfunction cross2

   logical function same(a1, b1)
   !< Equal bits.
   real(R8P), intent(in) :: a1(:), b1(:) !< Arrays.

   same = size(a1) == size(b1)
   if (same) same = all(transfer(a1, 0_I8P, size(a1)) == transfer(b1, 0_I8P, size(b1)))
   endfunction same

   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (.not. cond) then
      print '(A)', 'FAILED: '//label
      ok = .false.
   else
      print '(A)', 'ok: '//label
   endif
   endsubroutine check
endprogram fossil_test_cut_loads_2d
