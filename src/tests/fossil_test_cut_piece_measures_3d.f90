!< FOSSIL test: measures of the fluid pieces of the 3-D cut cells (issue #40).
program fossil_test_cut_piece_measures_3d
!< FOSSIL test: measures of the fluid pieces of the 3-D cut cells (issue #40).
!<
!< 1. Closure: on every trusted cell the sums over its pieces of the volume, first moments, wetted
!<    area, wetted moments and normal integral, and of each of the six face fractions, equal the cell
!<    totals (round-off). Bodies: the 5-degree wedge tilted off every axis (thin cells with two
!<    pieces cut obliquely), the cone, the tube, the L-prism, the thin plate, the hollow box, a plate
!<    with an inclusion; generic and grid-aligned grids. The index arrays (first_piece, piece_cell,
!<    np) are consistent.
!< 2. Prisms against the 2-D pieces of the section (independent): every 3-D piece of a cell inside
!<    the height has the volume (2-D area x dz), the centroid (2-D centroid, mid-height) and the -z,
!<    +z face fractions (2-D area / cell area) of a 2-D piece of the same cell (matched by centroid).
!<    The centroid is compared as a first moment (error x volume / (cell volume h)): a sliver of
!<    fraction 3e-4 has a centroid good to 2e-11 h only, the moment to 5e-15.
!< 3. Analytic: the thin plate (below 0.5 hx hy, above (0.5 - 1e-3) hx hy, face fractions); the
!<    hollow box in one cell (outer fluid 98, void 1) and over the void (fluid below 4, void 1); an
!<    inclusion above a plate (attached by the ray through the +x face) and left of a wall (attached
!<    by the ray hitting the wall): its volume leaves the piece around it; an inclusion inside the
!<    void (the void piece loses its volume).

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, pose_from_axis_angle
use fossil_curve_slice,       only : slice_surface
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(cut_cells_3d_t)      :: cut
type(soup_t)              :: so
real(R8P), allocatable    :: x(:), y(:), z(:)
real(R8P)                 :: err(4), emax(4), v(2)
integer(I4P)              :: i, nmulti, ntot, nmis, nmatch, p
logical                   :: ok(3)

call bench_bodies(b)
emax = 0._R8P ; nmulti = 0 ; ntot = 0

! 1. closure
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call run(s, x, y, z)
call to_surface(b(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call run(s, x, y, z)
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [-0.2_R8P, 0.45_R8P, 1.3_R8P]
call run(s, x, y, z)
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call run(s, x, y, z) ! aligned
so = tube([0._R8P, 3._R8P], [1._R8P, 2._R8P], [0._R8P, 3._R8P]) ; call to_surface(so, s)
x = [-0.5_R8P, 1.2_R8P, 1.8_R8P, 3.4_R8P] ; y = x ; z = [-0.5_R8P, 1.5_R8P, 3.5_R8P]
call run(s, x, y, z)
call s%move(pose_from_axis_angle([1._R8P, 0.2_R8P, 0.1_R8P], 0.3_R8P))
x = [(-0.9_R8P + 0.31_R8P * i, i = 0, 16)] ; y = x ; z = x
call run(s, x, y, z)
call to_surface(b(5)%soup, s) ! plate
call run(s, [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P, 1.3_R8P], [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P], [-0.5_R8P, 0.5_R8P])
so = hollow_box([1._R8P, 2._R8P], .false.) ; call to_surface(so, s)
call run(s, [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], [-1._R8P, 4._R8P])
x = [(-0.7_R8P + 0.47_R8P * i, i = 0, 10)]
call run(s, x, x, x)
so = plate_with_inclusion(.false.) ; call to_surface(so, s)
call run(s, [0._R8P, 1._R8P], [0._R8P, 1._R8P], [0._R8P, 1._R8P])
ok(1) = all(emax < 1.e-13_R8P) .and. nmulti > 50
print '(A,4ES10.2,A,I0,A,I0,A,L1)', 'closure: max rel err (volume, moments, wall, faces)', emax, '  cells ', ntot, &
      '  with 2+ pieces ', nmulti, '  ok ', ok(1)

! 2. prisms against the 2-D pieces of the section
call to_surface(b(2)%soup, s) ! wedge 5 rotated about z
call s%move(pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.37_R8P))
x = [(-0.13_R8P + 0.17_R8P * i, i = 0, 7)] ; y = [(-0.21_R8P + 0.043_R8P * i, i = 0, 21)] ; z = [0.2_R8P, 0.5_R8P, 0.8_R8P]
call prism_check(s, x, y, z, nmis, nmatch, err(1))
ok(2) = nmis == 0 .and. err(1) < 1.e-12_R8P ; p = nmatch
so = tube([0._R8P, 3._R8P], [1._R8P, 2._R8P], [0._R8P, 3._R8P]) ; call to_surface(so, s)
x = [-0.5_R8P, 1.2_R8P, 1.8_R8P, 3.4_R8P] ; y = x ; z = [0.5_R8P, 1.5_R8P, 2.5_R8P]
call prism_check(s, x, y, z, nmis, nmatch, err(2))
ok(2) = ok(2) .and. nmis == 0 .and. err(2) < 1.e-12_R8P ; p = p + nmatch
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [0.1_R8P, 0.45_R8P, 0.9_R8P]
call prism_check(s, x, y, z, nmis, nmatch, err(3))
ok(2) = ok(2) .and. nmis == 0 .and. err(3) < 1.e-12_R8P .and. p > 20 ; p = p + nmatch
print '(A,I0,A,3ES10.2,A,L1)', 'prisms (wedge, tube, L) vs 2-D pieces: matched pieces ', p, '  max err', err(1:3), &
      '  ok ', ok(2)

! 3. analytic
ok(3) = .true.
call to_surface(b(5)%soup, s) ! plate [0, 1]^2 x [0, 1e-3]
call cut_cells_3d(s, [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P, 1.3_R8P], [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P], &
                  [-0.5_R8P, 0.5_R8P], cut)
call two_vols(cut, 2, 2, 1, v)
ok(3) = ok(3) .and. cut%pieces(2, 2, 1) == 2 .and. abs(v(1) - 0.09_R8P * 0.5_R8P) < 1.e-15_R8P .and. &
        abs(v(2) - 0.09_R8P * 0.499_R8P) < 1.e-15_R8P
p = cut%first_piece(2, 2, 1)
if (cut%piece_volume(p) < cut%piece_volume(p + 1)) p = p + 1 ! the piece below
ok(3) = ok(3) .and. all(abs(cut%piece_face(:, p) - [0.5_R8P, 0.5_R8P, 0.5_R8P, 0.5_R8P, 1._R8P, 0._R8P]) < 1.e-14_R8P) .and. &
        abs(cut%piece_centroid(3, p) + 0.25_R8P) < 1.e-14_R8P .and. abs(cut%piece_wet_area(p) - 0.09_R8P) < 1.e-15_R8P .and. &
        abs(cut%piece_wet_normal(3, p) + 0.09_R8P) < 1.e-15_R8P
print '(A,2F12.8,A,L1)', 'thin plate: volumes ', v, '  ok ', ok(3)
so = hollow_box([1._R8P, 2._R8P], .false.) ; call to_surface(so, s)
call cut_cells_3d(s, [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], cut)
call two_vols(cut, 1, 1, 1, v)
ok(3) = ok(3) .and. abs(v(1) - 98._R8P) < 1.e-12_R8P .and. abs(v(2) - 1._R8P) < 1.e-13_R8P
print '(A,2F12.8)', 'hollow box in one cell: outer fluid, void ', v
call cut_cells_3d(s, [0.5_R8P, 2.5_R8P], [0.5_R8P, 2.5_R8P], [-1._R8P, 2.5_R8P], cut)
call two_vols(cut, 1, 1, 1, v)
ok(3) = ok(3) .and. abs(v(1) - 4._R8P) < 1.e-13_R8P .and. abs(v(2) - 1._R8P) < 1.e-13_R8P
print '(A,2F12.8)', 'over the void: fluid below, void ', v
so = hollow_box([1._R8P, 2._R8P], .true.) ; call to_surface(so, s)
call cut_cells_3d(s, [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], [-1._R8P, 4._R8P], cut)
call two_vols(cut, 1, 1, 1, v)
ok(3) = ok(3) .and. cut%pieces(1, 1, 1) == 2 .and. abs(v(1) - 98._R8P) < 1.e-12_R8P .and. &
        abs(v(2) - (1._R8P - 0.008_R8P)) < 1.e-13_R8P
print '(A,2F12.8)', 'inclusion inside the void: outer fluid, void ', v
so = plate_with_inclusion(.false.) ; call to_surface(so, s)
call cut_cells_3d(s, [0._R8P, 1._R8P], [0._R8P, 1._R8P], [0._R8P, 1._R8P], cut)
call two_vols(cut, 1, 1, 1, v)
ok(3) = ok(3) .and. cut%pieces(1, 1, 1) == 2 .and. abs(v(1) - (0.79_R8P - 0.008_R8P)) < 1.e-14_R8P .and. &
        abs(v(2) - 0.2_R8P) < 1.e-14_R8P
print '(A,2F12.8)', 'inclusion above a plate (ray through the +x face): above, below ', v
so = plate_with_inclusion(.true.) ; call to_surface(so, s)
call cut_cells_3d(s, [0._R8P, 1._R8P], [0._R8P, 1._R8P], [0._R8P, 1._R8P], cut)
call two_vols(cut, 1, 1, 1, v)
ok(3) = ok(3) .and. cut%pieces(1, 1, 1) == 2 .and. abs(v(1) - (0.8_R8P - 0.008_R8P)) < 1.e-14_R8P .and. &
        abs(v(2) - 0.19_R8P) < 1.e-14_R8P
print '(A,2F12.8,A,L1)', 'inclusion left of a wall (ray hits the wall): left, right ', v, '  ok ', ok(3)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine run(sf, xe, ye, ze)
   !< Cut and check the closure (updates emax, ntot, nmulti).
   type(surface_stl_object), intent(in) :: sf                !< Surface.
   real(R8P),                intent(in) :: xe(:), ye(:), ze(:) !< Edges.
   type(cut_cells_3d_t) :: c3   !< Cut cells.
   real(R8P)            :: e(4) !< Errors.
   integer(I4P)         :: nm, nt_ !< Counts.

   call cut_cells_3d(sf, xe, ye, ze, c3)
   call closure(c3, xe, ye, ze, e, nm, nt_)
   emax = max(emax, e) ; nmulti = nmulti + nm ; ntot = ntot + nt_
   endsubroutine run

   subroutine closure(c3, xe, ye, ze, e, nm, nt_)
   !< Max relative closure errors over the trusted cells: volume (/ cell volume), first moments
   !< (/ cell volume h), wall terms (area, moments, normal; / h^2, h^3), faces (fractions).
   type(cut_cells_3d_t), intent(in)  :: c3                  !< Cut cells.
   real(R8P),            intent(in)  :: xe(:), ye(:), ze(:) !< Edges.
   real(R8P),            intent(out) :: e(4)                !< Errors.
   integer(I4P),         intent(out) :: nm, nt_             !< Cells with 2+ pieces, trusted cells.
   integer(I4P) :: ii, jj, kk, q, q0, nexp !< Work.
   real(R8P)    :: h(3), vc, hm, vs, ms(3), ws, wm(3), wn(3), fs(6), ft(6) !< Work.

   e = 0._R8P ; nm = 0 ; nt_ = 0 ; nexp = 0
   do kk = 1, c3%nz
      do jj = 1, c3%ny
         do ii = 1, c3%nx
            if (c3%first_piece(ii, jj, kk) /= nexp + 1) e = huge(1._R8P)
            if (c3%pieces(ii, jj, kk) < 1) cycle
            nexp = nexp + c3%pieces(ii, jj, kk)
            nt_ = nt_ + 1 ; if (c3%pieces(ii, jj, kk) >= 2) nm = nm + 1
            h = [xe(ii + 1) - xe(ii), ye(jj + 1) - ye(jj), ze(kk + 1) - ze(kk)] ; vc = product(h) ; hm = maxval(h)
            vs = 0._R8P ; ms = 0._R8P ; ws = 0._R8P ; wm = 0._R8P ; wn = 0._R8P ; fs = 0._R8P
            q0 = c3%first_piece(ii, jj, kk)
            do q = q0, q0 + c3%pieces(ii, jj, kk) - 1
               if (any(c3%piece_cell(:, q) /= [ii, jj, kk])) e = huge(1._R8P)
               vs = vs + c3%piece_volume(q) ; ms = ms + c3%piece_volume(q) * c3%piece_centroid(:, q)
               ws = ws + c3%piece_wet_area(q) ; wm = wm + c3%piece_wet_area(q) * c3%piece_wet_centroid(:, q)
               wn = wn + c3%piece_wet_normal(:, q) ; fs = fs + c3%piece_face(:, q)
            enddo
            ft = [c3%face_x(ii, jj, kk), c3%face_x(ii + 1, jj, kk), c3%face_y(ii, jj, kk), c3%face_y(ii, jj + 1, kk), &
                  c3%face_z(ii, jj, kk), c3%face_z(ii, jj, kk + 1)]
            e(1) = max(e(1), abs(vs - c3%alpha(ii, jj, kk) * vc) / vc)
            e(2) = max(e(2), maxval(abs(ms - c3%alpha(ii, jj, kk) * vc * c3%centroid(:, ii, jj, kk))) / (vc * hm))
            e(3) = max(e(3), abs(ws - c3%wet_area(ii, jj, kk)) / hm**2, &
                       maxval(abs(wm - c3%wet_area(ii, jj, kk) * c3%wet_centroid(:, ii, jj, kk))) / hm**3, &
                       maxval(abs(wn - c3%wet_normal(:, ii, jj, kk))) / hm**2)
            e(4) = max(e(4), maxval(abs(fs - ft)))
         enddo
      enddo
   enddo
   if (nexp /= c3%np) e = huge(1._R8P)
   endsubroutine closure

   subroutine prism_check(sf, xe, ye, ze, nmis, nmat, er)
   !< Every 3-D piece of the cells against the 2-D pieces of the section at mid-height.
   type(surface_stl_object), intent(in)  :: sf         !< Prism (extruded along z).
   real(R8P),                intent(in)  :: xe(:), ye(:), ze(:) !< Edges (ze inside the height).
   integer(I4P),             intent(out) :: nmis, nmat !< Unmatched pieces, matched pieces.
   real(R8P),                intent(out) :: er         !< Max relative error (volume, first moments, z faces).
   type(cut_cells_3d_t) :: c3  !< 3-D cut cells.
   type(curve_object)   :: sec !< Section.
   type(cut_cells_2d_t) :: c2  !< Its cut cells.
   integer(I4P)         :: ii, jj, kk, q, r, rb !< Work.
   real(R8P)            :: h(3), d, db, a2 !< Work.

   call cut_cells_3d(sf, xe, ye, ze, c3)
   call slice_surface(surface=sf, normal_axis=3_I4P, level=0.5_R8P * (ze(1) + ze(2)), curve=sec, on_plane_low=.true.)
   call sec%cut_cells(x=xe, y=ye, cut=c2)
   nmis = 0 ; nmat = 0 ; er = 0._R8P
   do kk = 1, size(ze) - 1
      do jj = 1, size(ye) - 1
         do ii = 1, size(xe) - 1
            if (c3%pieces(ii, jj, kk) /= c2%pieces(ii, jj)) then
               nmis = nmis + 1 ; cycle
            endif
            h = [xe(ii + 1) - xe(ii), ye(jj + 1) - ye(jj), ze(kk + 1) - ze(kk)]
            do q = c3%first_piece(ii, jj, kk), c3%first_piece(ii, jj, kk) + c3%pieces(ii, jj, kk) - 1
               db = huge(1._R8P) ; rb = 0
               do r = c2%first_piece(ii, jj), c2%first_piece(ii, jj) + c2%pieces(ii, jj) - 1
                  d = norm2(c3%piece_centroid(1:2, q) - c2%piece_centroid(:, r))
                  if (d < db) then
                     db = d ; rb = r
                  endif
               enddo
               if (rb == 0) then
                  nmis = nmis + 1 ; cycle
               endif
               a2 = c2%piece_area(rb)
               nmat = nmat + 1
               er = max(er, abs(c3%piece_volume(q) - a2 * h(3)) / product(h), db * c3%piece_volume(q) / (product(h) * maxval(h)), &
                        abs(c3%piece_centroid(3, q) - 0.5_R8P * (ze(kk) + ze(kk + 1))) / maxval(h), &
                        abs(c3%piece_face(5, q) - a2 / (h(1) * h(2))), abs(c3%piece_face(6, q) - a2 / (h(1) * h(2))))
            enddo
         enddo
      enddo
   enddo
   endsubroutine prism_check

   subroutine two_vols(c3, ii, jj, kk, vv)
   !< The volumes of the two pieces of a cell, the larger first.
   type(cut_cells_3d_t), intent(in)  :: c3         !< Cut cells.
   integer,              intent(in)  :: ii, jj, kk !< Cell.
   real(R8P),            intent(out) :: vv(2)      !< Volumes.
   integer(I4P)                      :: q          !< First piece.

   vv = -1._R8P
   if (c3%pieces(ii, jj, kk) /= 2) return
   q = c3%first_piece(ii, jj, kk)
   vv = [maxval(c3%piece_volume(q:q + 1)), minval(c3%piece_volume(q:q + 1))]
   endsubroutine two_vols

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

   subroutine append(sp, add)
   !< Append the triangles of add to sp.
   type(soup_t), intent(inout) :: sp  !< Soup.
   type(soup_t), intent(in)    :: add !< Added.

   sp%v = reshape([reshape(sp%v, [size(sp%v)]), reshape(add%v, [size(add%v)])], [3, 3, size(sp%v, 3) + size(add%v, 3)])
   endsubroutine append

   function hollow_box(vr, inclusion) result(sp)
   !< Box [0, 3]^3 with the void [vr(1), vr(2)]^3 (inner faces facing the void); optionally the solid
   !< box [1.4, 1.6]^3 inside the void.
   real(R8P), intent(in) :: vr(2)     !< Void span.
   logical,   intent(in) :: inclusion !< Add the inclusion.
   type(soup_t)          :: sp        !< Soup.

   sp = box_soup([0._R8P, 0._R8P, 0._R8P], [3._R8P, 3._R8P, 3._R8P], .false.)
   call append(sp, box_soup([vr(1), vr(1), vr(1)], [vr(2), vr(2), vr(2)], .true.))
   if (inclusion) call append(sp, box_soup([1.4_R8P, 1.4_R8P, 1.4_R8P], [1.6_R8P, 1.6_R8P, 1.6_R8P], .false.))
   endfunction hollow_box

   function plate_with_inclusion(wall) result(sp)
   !< The box [0.4, 0.6] x [0.4, 0.6] x [0.5, 0.7] with a plate across the cell [0, 1]^3: the slab
   !< z in [0.2, 0.21] below it, or (wall) the slab x in [0.8, 0.81] right of it.
   logical, intent(in) :: wall !< Vertical plate.
   type(soup_t)        :: sp   !< Soup.

   if (wall) then
      sp = box_soup([0.8_R8P, -1._R8P, -1._R8P], [0.81_R8P, 2._R8P, 2._R8P], .false.)
   else
      sp = box_soup([-1._R8P, -1._R8P, 0.2_R8P], [2._R8P, 2._R8P, 0.21_R8P], .false.)
   endif
   call append(sp, box_soup([0.4_R8P, 0.4_R8P, 0.5_R8P], [0.6_R8P, 0.6_R8P, 0.7_R8P], .false.))
   endfunction plate_with_inclusion

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
endprogram fossil_test_cut_piece_measures_3d
