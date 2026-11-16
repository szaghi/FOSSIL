!< FOSSIL test: the loads of a body from the wall terms of its 3-D cut cells (#109 WP2 step 05).
program fossil_test_cut_loads_3d
!< FOSSIL test: the loads of a body from the wall terms of its 3-D cut cells (#109 WP2 step 05).
!<
!< Bodies inside the grid (the walls then cover the whole surface): cube-01 on a grid through its faces and on one off
!< them, an icosphere (1280 facets) moved to a pose, two boxes, an irregular tetrahedron.
!<
!< 1. Buoyancy through the cut: p = p0 + g . x at the walls' centroids (`cut_loads_linear_pressure`), tau = 0. The force
!<    is F = -g V (V, x_c from `mass_properties`) and, with the gradient and the walls' gyration (`cut_walls_gyration`),
!<    the moment about c is M = V (x_c - c) x (-g): each component within TOL eps times the sum of the absolute terms
!<    (S_F = sum |p_w| |A_w| + V |g|, S_M = sum |p_w| |A_w| |x_w - c| + sum |K_w g| |A_w| + V |g| |x_c - c|); the largest
!<    ratio is printed. Without the gradient (the centroid value only) the moment error is printed, not checked.
!< 2. Uniform pressure: the walls' vector areas sum to 0 (4 eps sum |A_w|).
!< 3. The walls of each piece sum to `piece_wet_normal` (8 eps sum |A_w| of the piece); the bit-equal pieces are counted.
!< 4. FD / FV: the same linear p integrated by `integrate_loads` over load samples of order 2 and by `cut_loads_3d` over
!<    the walls: the same force and moment within TOL eps (S_FD + S_FV).
!< 5. Uniform shear: F = tau A, the area output = the surface's area (1e-14).
!< 6. Per patch and per body: three patches over two boxes add up to the totals; each box its own buoyancy.
!< 7. Per piece: `pressure_per_piece` = the per-wall pressure spread by hand, bit for bit; an open triangle (untrusted
!<    cells): its terms take no pressure, counted in `unassigned`.
!< 8. Refusals.
use fossil,                  only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t, cut_loads_3d, &
                                    cut_walls_gyration, cut_loads_linear_pressure, load_samples_t, load_samples_3d, &
                                    integrate_loads, rigid_pose_t, pose_from_axis_angle, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_boxes, make_sphere
use penf,                    only : I4P, I8P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: EPS = epsilon(1._R8P), TOL = 1._R8P
real(R8P), parameter            :: G(3) = [0.3_R8P, -0.7_R8P, 1.1_R8P], P0 = 2._R8P, C(3) = [0.1_R8P, -0.2_R8P, 0.3_R8P]
real(R8P), parameter            :: TAU(3) = [0.2_R8P, -0.5_R8P, 0.7_R8P]
type(surface_stl_object)        :: cube, sphere, boxes, tet, plate
type(facet_object), allocatable :: fc(:)
type(cut_cells_3d_t)            :: cut
type(cut_walls_3d_t)            :: walls
type(rigid_pose_t)              :: pose
real(R8P), allocatable          :: x(:), y(:), z(:), p(:), t(:,:), gr(:,:), gy(:,:), ppp(:), pw(:), ppf(:,:), ppm(:,:), &
                                   pbf(:,:), pbm(:,:)
real(R8P)                       :: f(3), m(3), f2(3), m2(3), a, rf, rm, rm0, rfd, ria, sf, sm
integer(I4P)                    :: st, k, w, nbit, npc, nun
logical                         :: ok

ok = .true. ; rf = 0._R8P ; rm = 0._R8P ; rm0 = 0._R8P ; rfd = 0._R8P ; ria = 0._R8P ; nbit = 0_I4P ; npc = 0_I4P
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call make_sphere(3_I4P, fc) ; call sphere%adopt_facets(facets=fc)
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.7_R8P, -0.3_R8P, 1.2_R8P])
call sphere%move(pose)
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, 0.25_R8P, fc) ; call boxes%adopt_facets(facets=fc)
call make_tet([0._R8P, 0._R8P, 0._R8P], [1.3_R8P, 0.1_R8P, 0.2_R8P], [0.2_R8P, 0.9_R8P, -0.1_R8P], [0.4_R8P, 0.3_R8P, 1.1_R8P], fc)
call tet%adopt_facets(facets=fc)
! 1-4. buoyancy, uniform pressure, pieces, FD / FV
x = [(-0.25_R8P + 0.25_R8P * k, k = 0, 6)] ; y = x ; z = x
call body(cube, 'cube-01 on grid planes')
call box_lines(cube, 9_I4P) ; call body(cube, 'cube-01 off the planes')
call box_lines(sphere, 12_I4P) ; call body(sphere, 'icosphere moved')
call box_lines(boxes, 10_I4P) ; call body(boxes, 'two boxes')
call box_lines(tet, 8_I4P) ; call body(tet, 'irregular tetrahedron')
print '(A,F6.2,A,F6.2,A,F6.2,A,ES9.2)', 'largest error / (eps S): force ', rf, ', moment ', rm, ', FD vs FV ', rfd, &
      ';  centroid-only moment error / S_M ', rm0
print '(A,I0,A,I0,A,F6.2,A)', 'piece identity: ', nbit, ' of ', npc, ' pieces bit for bit, largest error ', ria, &
      ' eps sum |A|'
! 5. uniform shear
call box_lines(sphere, 12_I4P)
call cut_cells_3d(sphere, x, y, z, cut, status=st, walls=walls)
allocate(t(3, walls%n))
do w = 1, walls%n
   t(:, w) = TAU
enddo
call cut_loads_3d(walls, f, m, status=st, shear=t, area=a)
call check('uniform shear: F = tau A (4 eps)', st == STATUS_OK .and. all(abs(f - TAU * a) <= 4._R8P * EPS * abs(TAU) * a))
call check('uniform shear: area = the surface area (1e-14)', abs(a - surface_area(sphere)) <= 1.e-14_R8P * a)
! 6. per patch and per body
call boxes%set_patch([1_I4P, 2_I4P, 3_I4P, 4_I4P, 5_I4P], 'two')
call boxes%set_patch([(k, k = 13_I4P, 20_I4P)], 'three')
call box_lines(boxes, 10_I4P)
call cut_cells_3d(boxes, x, y, z, cut, status=st, walls=walls)
p = cut_loads_linear_pressure(walls%centroid, P0, G)
call cut_walls_gyration(walls, boxes, x, y, z, gy, status=st)
gr = spread(G, 2, walls%n)
deallocate(t) ; allocate(t(3, walls%n))
do w = 1, walls%n
   t(:, w) = [sin(3 * walls%centroid(2, w)), cos(2 * walls%centroid(3, w)), 0.3_R8P * walls%centroid(1, w)]
enddo
call cut_loads_3d(walls, f, m, status=st, pressure=p, shear=t, gradient=gr, gyration=gy, about=C, surface=boxes, &
                  per_patch_force=ppf, per_patch_moment=ppm, per_body_force=pbf, per_body_moment=pbm, area=a)
call scales(0._R8P, [0._R8P, 0._R8P, 0._R8P], sf, sm)
call check('per patch: three patches, two bodies', st == STATUS_OK .and. size(ppf, 2) == 3 .and. size(pbf, 2) == 2)
call check('per patch: add up to the totals (4 eps S)', all(abs(sum(ppf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           all(abs(sum(ppm, 2) - m) <= 4._R8P * EPS * (sm + 4 * a)))
call check('per body: add up to the totals (4 eps S)', all(abs(sum(pbf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           all(abs(sum(pbm, 2) - m) <= 4._R8P * EPS * (sm + 4 * a)))
call cut_loads_3d(walls, f, m, status=st, pressure=p, gradient=gr, gyration=gy, about=C, surface=boxes, &
                  per_body_force=pbf, per_body_moment=pbm)
call scales(1._R8P, [1._R8P, 0.5_R8P, 0.5_R8P], sf, sm)
call check('per body: box 1 buoyancy', all(abs(pbf(:, 1) + G) <= TOL * EPS * sf) .and. &
           all(abs(pbm(:, 1) - cross([0.5_R8P, 0.5_R8P, 0.5_R8P] - C, -G)) <= TOL * EPS * sm))
call check('per body: box 2 buoyancy', all(abs(pbf(:, 2) + G) <= TOL * EPS * sf) .and. &
           all(abs(pbm(:, 2) - cross([1.75_R8P, 0.5_R8P, 0.5_R8P] - C, -G)) <= TOL * EPS * sm))
! 7. per piece
allocate(ppp(cut%np))
do k = 1, cut%np
   ppp(k) = 1._R8P + 0.1_R8P * sin(real(k, R8P))
enddo
allocate(pw(walls%n))
do w = 1, walls%n
   pw(w) = 0._R8P ; if (walls%piece(w) >= 1) pw(w) = ppp(walls%piece(w))
enddo
call cut_loads_3d(walls, f, m, status=st, pressure_per_piece=ppp, shear=t, about=C, unassigned=nun)
call cut_loads_3d(walls, f2, m2, pressure=pw, shear=t, about=C)
call check('per piece: = the pressure spread by hand, bit for bit', st == STATUS_OK .and. nun == 0 .and. same(f, f2) .and. &
           same(m, m2))
allocate(fc(1)) ; call put(fc(1), [0.11_R8P, 0.13_R8P, 0.17_R8P], [0.83_R8P, 0.21_R8P, 0.42_R8P], &
                           [0.29_R8P, 0.88_R8P, 0.71_R8P])
call plate%adopt_facets(facets=fc)
x = [(0.05_R8P * k, k = 0, 20)] ; y = x ; z = x
call cut_cells_3d(plate, x, y, z, cut, status=st, walls=walls)
deallocate(ppp) ; allocate(ppp(max(1_I4P, cut%np))) ; ppp = 1._R8P
call cut_loads_3d(walls, f, m, status=st, pressure_per_piece=ppp, unassigned=nun)
call check('open triangle: terms of no piece take no pressure, counted', st == STATUS_OK .and. walls%n > 0 .and. &
           nun == count(walls%piece == 0) .and. nun > 0)
! 8. refusals
call cut_loads_3d(walls, f, m, status=st)
call check('refused: neither p nor tau', st == STATUS_INVALID_INPUT .and. all(f == 0._R8P) .and. all(m == 0._R8P))
pw = [(1._R8P, w = 1, walls%n)]
call cut_loads_3d(walls, f, m, status=st, pressure=pw, pressure_per_piece=ppp)
call check('refused: two pressures', st == STATUS_INVALID_INPUT)
call cut_loads_3d(walls, f, m, status=st, pressure=pw(1:walls%n - 1))
call check('refused: p of the wrong size', st == STATUS_INVALID_INPUT)
call cut_loads_3d(walls, f, m, status=st, pressure=pw, gradient=spread(G, 2, walls%n))
call check('refused: gradient without gyration', st == STATUS_INVALID_INPUT)
call cut_loads_3d(walls, f, m, status=st, pressure=pw, per_body_force=pbf)
call check('refused: per body without the surface, size 0', st == STATUS_INVALID_INPUT .and. size(pbf, 2) == 0)
call cut_loads_3d(walls, f(1:2), m, status=st, pressure=pw)
call check('refused: force of size 2', st == STATUS_INVALID_INPUT)
call cut_walls_gyration(walls, plate, x(1:10), y, z, gy, status=st)
call check('refused: gyration on other cell edges', st == STATUS_INVALID_INPUT .and. size(gy, 2) == 0)

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine body(surf, name)
   !< Checks 1-4 on one body inside the lines x, y, z.
   type(surface_stl_object), intent(in) :: surf !< Body.
   character(*),             intent(in) :: name !< Label.
   type(load_samples_t)     :: ls                         !< Samples.
   real(R8P), allocatable   :: pl(:), one(:), ac(:,:)     !< Pressures, piece sums.
   real(R8P)                :: v, x0(3), ff(3), mm(3), fr(3), mr(3), ssf, ssm, e1, e2, ffd(3), mfd(3), sfd, smd, sa
   integer(I4P)             :: sst, q, pc                 !< Work.
   logical                  :: bits                       !< Bit-equal pieces.

   call surf%mass_properties(volume=v, center=x0, status=sst)
   call cut_cells_3d(surf, x, y, z, cut, status=sst, walls=walls)
   call check(name//': cut with walls', sst == STATUS_OK .and. walls%n > 0)
   pl = cut_loads_linear_pressure(walls%centroid, P0, G)
   call cut_walls_gyration(walls, surf, x, y, z, gy, status=sst)
   call check(name//': gyration', sst == STATUS_OK .and. size(gy, 2) == walls%n)
   gr = spread(G, 2, walls%n)
   fr = -G * v ; mr = cross(v * (x0 - C), -G)
   call cut_loads_3d(walls, ff, mm, status=sst, pressure=pl, gradient=gr, gyration=gy, about=C)
   call scales(v, x0, ssf, ssm)
   e1 = maxval(abs(ff - fr)) / (EPS * ssf) ; e2 = maxval(abs(mm - mr)) / (EPS * ssm)
   rf = max(rf, e1) ; rm = max(rm, e2)
   call check(name//': buoyancy F = -g V', sst == STATUS_OK .and. e1 <= TOL)
   call check(name//': buoyancy M = V (xc - c) x (-g), with the gyration', e2 <= TOL)
   call cut_loads_3d(walls, ff, mm, status=sst, pressure=pl, about=C)
   rm0 = max(rm0, maxval(abs(mm - mr)) / ssm)
   ! 2. uniform pressure
   allocate(one(walls%n)) ; one = 1._R8P
   call cut_loads_3d(walls, ff, mm, status=sst, pressure=one, area=sa)
   call check(name//': uniform p, sum of the vector areas 0 (4 eps)', all(abs(ff) <= 4._R8P * EPS * sa))
   ! 3. the walls of each piece
   allocate(ac(4, cut%np)) ; ac = 0._R8P
   do q = 1, walls%n
      pc = walls%piece(q)
      if (pc < 1) cycle
      ac(1:3, pc) = ac(1:3, pc) + walls%area(:, q) ; ac(4, pc) = ac(4, pc) + walls%wet_area(q)
   enddo
   do pc = 1, cut%np
      if (ac(4, pc) <= 0._R8P) cycle
      npc = npc + 1
      bits = same(ac(1:3, pc), cut%piece_wet_normal(:, pc))
      if (bits) nbit = nbit + 1
      ria = max(ria, maxval(abs(ac(1:3, pc) - cut%piece_wet_normal(:, pc))) / (EPS * ac(4, pc)))
   enddo
   call check(name//': walls of a piece = piece_wet_normal (8 eps sum |A|)', ria <= 8._R8P)
   ! 4. FD / FV
   call load_samples_3d(surf, samples=ls, order=2_I4P)
   call integrate_loads(ls, ffd, mfd, status=sst, pressure=cut_loads_linear_pressure(ls%points, P0, G), about=C)
   call cut_loads_3d(walls, ff, mm, status=sst, pressure=pl, gradient=gr, gyration=gy, about=C)
   sfd = v * norm2(G) ; smd = v * norm2(G) * norm2(x0 - C)
   do q = 1, ls%n
      sfd = sfd + ls%weight(q) * abs(P0 + dot_product(G, ls%points(:, q)))
      smd = smd + ls%weight(q) * abs(P0 + dot_product(G, ls%points(:, q))) * norm2(ls%points(:, q) - C)
   enddo
   e1 = maxval(abs(ff - ffd)) / (EPS * (ssf + sfd)) ; e2 = maxval(abs(mm - mfd)) / (EPS * (ssm + smd))
   rfd = max(rfd, e1, e2)
   call check(name//': FD (integrate_loads) = FV (cut_loads), force and moment', e1 <= TOL .and. e2 <= TOL)
   endsubroutine body

   subroutine box_lines(surf, n)
   !< Lines of n cells around the body's box, with uneven margins (no plane through a vertex by design).
   type(surface_stl_object), intent(in) :: surf !< Body.
   integer(I4P),             intent(in) :: n    !< Cells per axis.
   type(vector_R8P) :: lo, hi                   !< Box.
   real(R8P)        :: l(3), h(3)               !< Corners.
   integer(I4P)     :: q                        !< Counter.

   lo = surf%get_bmin() ; hi = surf%get_bmax()
   l = [lo%x, lo%y, lo%z] - [0.137_R8P, 0.093_R8P, 0.171_R8P] ; h = [hi%x, hi%y, hi%z] + [0.119_R8P, 0.157_R8P, 0.083_R8P]
   x = [(l(1) + (h(1) - l(1)) * q / n, q = 0, n)]
   y = [(l(2) + (h(2) - l(2)) * q / n, q = 0, n)]
   z = [(l(3) + (h(3) - l(3)) * q / n, q = 0, n)]
   endsubroutine box_lines

   subroutine scales(v, x0, ssf, ssm)
   !< The sums of the absolute terms of the walls (`walls`, pressure p0 + g . x, gyration gy).
   real(R8P),    intent(in)  :: v, x0(3) !< Volume, centroid.
   real(R8P),    intent(out) :: ssf, ssm !< Scales.
   real(R8P)    :: pp, kg(3)             !< Work.
   integer(I4P) :: q                     !< Counter.

   ssf = v * norm2(G) ; ssm = v * norm2(G) * norm2(x0 - C)
   do q = 1, walls%n
      pp = abs(P0 + dot_product(G, walls%centroid(:, q)))
      ssf = ssf + pp * walls%wet_area(q) ; ssm = ssm + pp * walls%wet_area(q) * norm2(walls%centroid(:, q) - C)
      if (size(gy, 2) == walls%n) then
         kg = [gy(1, q) * G(1) + gy(4, q) * G(2) + gy(6, q) * G(3), gy(4, q) * G(1) + gy(2, q) * G(2) + gy(5, q) * G(3), &
               gy(6, q) * G(1) + gy(5, q) * G(2) + gy(3, q) * G(3)]
         ssm = ssm + norm2(kg) * walls%wet_area(q)
      endif
   enddo
   endsubroutine scales

   real(R8P) function surface_area(surf) result(aa)
   !< Sum of the facet areas.
   type(surface_stl_object), intent(in) :: surf  !< Surface.
   real(R8P)    :: vv(3, 3)                      !< Vertices.
   integer(I4P) :: ff                            !< Counter.

   aa = 0._R8P
   do ff = 1, surf%get_facets_number()
      vv = reshape(surf%tri_vertices(:, ff), [3, 3])
      aa = aa + 0.5_R8P * norm2(cross(vv(:, 2) - vv(:, 1), vv(:, 3) - vv(:, 1)))
   enddo
   endfunction surface_area

   subroutine make_tet(a, b, cc, d, ff)
   !< The tetrahedron a b cc d, facets outward.
   real(R8P),                       intent(in)  :: a(3), b(3), cc(3), d(3) !< Vertices.
   type(facet_object), allocatable, intent(out) :: ff(:)                   !< Facets.

   allocate(ff(4))
   call put(ff(1), a, cc, b) ; call put(ff(2), a, b, d) ; call put(ff(3), b, cc, d) ; call put(ff(4), cc, a, d)
   endsubroutine make_tet

   subroutine put(fo, a, b, cc)
   !< One facet.
   type(facet_object), intent(inout) :: fo                !< Facet.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo%vertex(1) = vector_R8P(a(1), a(2), a(3)) ; fo%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo%vertex(3) = vector_R8P(cc(1), cc(2), cc(3)) ; call fo%compute_metrix
   endsubroutine put

   logical function same(u, v)
   !< Equal bits.
   real(R8P), intent(in) :: u(:), v(:) !< Arrays.

   same = size(u) == size(v)
   if (same) same = all(transfer(u, 0_I8P, size(u)) == transfer(v, 0_I8P, size(v)))
   endfunction same

   pure function cross(u, v) result(r)
   !< Cross product.
   real(R8P), intent(in) :: u(3), v(3) !< Vectors.
   real(R8P)             :: r(3)       !< u cross v.

   r = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
   endfunction cross

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
endprogram fossil_test_cut_loads_3d
