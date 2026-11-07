!< FOSSIL test: geometric moments of the cut-cell pieces (issue #70, step 8).
program fossil_test_cut_moments
!< FOSSIL test: geometric moments of the cut-cell pieces (issue #70, step 8).
!<
!< Moments up to order 4 (`moments_order=4`), int X^a Y^b (Z^c) over each piece, in the scaled coordinates of its cell.
!< 3-D: the 5-degree wedge tilted off every axis, the L-prism on a grid-aligned grid, the thin plate.
!< 1. summed over the pieces, re-centred to the global monomials x^a y^b z^c, they equal the block's minus the body's,
!<    the latter integrated independently (tetrahedra from the origin, collapsed 4x4x4 Gauss) (relative 1e-12);
!< 2. the zeroth moment is the piece volume, the first ones its centroid offset (1e-13 of the cell volume); the zeroth
!<    wall moment is the wetted area, the first ones its centroid (1e-13 of the cell);
!< 3. the face moments of the pieces of a cell add up to the face moments of the plane section, and their zeroth is
!<    the piece face fraction times the face area (1e-13 of the face);
!< 4. 2-D (a circle of 64 segments): the same sum against the block minus the polygon (relative 1e-13; 3.6e-15
!<    measured), the zeroth wall moment against the wetted length; the zeroth and first moments against the areas and
!<    centroids (1e-12: those centroids come from Green in absolute coordinates, 1.7e-13 of the cell area here).

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_cut_moments,       only : mono_count, mono_exponents
use fossil_facet_object,      only : facet_object
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
type(cut_cells_2d_t)      :: c2
real(R8P), allocatable    :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                 :: err(5)
integer(I4P)              :: i, n
logical                   :: ok

err = 0._R8P
call bench_bodies(b)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.142_R8P * i, i = 0, 11)] ; y = [(-0.67_R8P + 0.126_R8P * i, i = 0, 12)]
z = [(-0.59_R8P + 0.154_R8P * i, i = 0, 13)] ! the wedge reaches z = 1.27
call check3(x, y, z)
call to_surface(b(8)%soup, s) ! L-prism, aligned
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call check3(x, y, z)
call to_surface(b(5)%soup, s) ! plate
x = [(-0.15_R8P + 0.2_R8P * i, i = 0, 7)] ; z = [-0.3_R8P, -0.1_R8P, 0.1_R8P, 0.3_R8P]
call check3(x, x, z)
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call check2(x, y)
print '(A,5ES10.2)', 'max errors (sums 3-D, piece measures, faces, 2-D sums and wall, 2-D low orders): ', err
ok = err(1) <= 1.e-12_R8P .and. err(2) <= 1.e-13_R8P .and. err(3) <= 1.e-13_R8P .and. err(4) <= 1.e-13_R8P .and. &
     err(5) <= 1.e-12_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check3(xx, yy, zz)
   !< Checks 1-3 on the block xx x yy x zz.
   real(R8P), intent(in)     :: xx(:), yy(:), zz(:)
   type(cut_cells_3d_t)      :: cut
   integer(I4P), allocatable :: e(:,:), e2(:,:)
   real(R8P),    allocatable :: glob(:), ref(:), fs(:)
   integer(I4P)              :: p, k, l, nm, ic(3), f, ax, b1, b2, q
   real(R8P)                 :: cc(3), h(3), vc, af, scale
   type(facet_object), pointer :: fp

   call cut_cells_3d(s, xx, yy, zz, cut, moments_order=4_I4P)
   nm = mono_count(3_I4P, 4_I4P) ; e = mono_exponents(3_I4P, 4_I4P) ; e2 = mono_exponents(2_I4P, 4_I4P)
   allocate(glob(nm), ref(nm), fs(size(e2, 2))) ; glob = 0._R8P ; ref = 0._R8P
   do p = 1, cut%np
      ic = cut%piece_cell(:, p)
      h = [xx(ic(1) + 1) - xx(ic(1)), yy(ic(2) + 1) - yy(ic(2)), zz(ic(3) + 1) - zz(ic(3))] ; vc = product(h)
      cc = [xx(ic(1)), yy(ic(2)), zz(ic(3))] + 0.5_R8P * h
      do k = 1, nm
         do l = 1, nm
            if (any(e(:, l) > e(:, k))) cycle
            glob(k) = glob(k) + product([(binom(e(q, k), e(q, l)) * cc(q) ** (e(q, k) - e(q, l)) * h(q) ** e(q, l), q = 1, 3)]) &
                      * cut%piece_moments(l, p)
         enddo
      enddo
      ! 2. low orders
      err(2) = max(err(2), abs(cut%piece_moments(1, p) - cut%piece_volume(p)) / vc)
      do q = 1, 3
         err(2) = max(err(2), abs(cut%piece_moments(1 + q, p) - cut%piece_volume(p) * (cut%piece_centroid(q, p) - cc(q)) / h(q)) &
                  / vc)
      enddo
      scale = maxval(h) ** 2
      err(2) = max(err(2), abs(cut%piece_wall_moments(1, p) - cut%piece_wet_area(p)) / scale)
      do q = 1, 3
         err(2) = max(err(2), abs(cut%piece_wall_moments(1 + q, p) - cut%piece_wet_area(p) * &
                  (cut%piece_wet_centroid(q, p) - cc(q)) / h(q)) / scale)
      enddo
      do f = 1, 6 ! 3. zeroth face moment
         ax = (f + 1) / 2 ; b1 = mod(ax, 3) + 1 ; b2 = mod(ax + 1, 3) + 1 ; af = h(b1) * h(b2)
         err(3) = max(err(3), abs(cut%piece_face_moments(1, f, p) - cut%piece_face(f, p) * af) / af)
      enddo
   enddo
   do k = 1, nm ! block minus the body
      ref(k) = (xx(size(xx)) ** (e(1, k) + 1) - xx(1) ** (e(1, k) + 1)) / (e(1, k) + 1) * &
               (yy(size(yy)) ** (e(2, k) + 1) - yy(1) ** (e(2, k) + 1)) / (e(2, k) + 1) * &
               (zz(size(zz)) ** (e(3, k) + 1) - zz(1) ** (e(3, k) + 1)) / (e(3, k) + 1)
      do f = 1, s%get_facets_number()
         fp => s%facet_at(f)
         ref(k) = ref(k) - tet_int([fp%vertex(1)%x, fp%vertex(1)%y, fp%vertex(1)%z], &
                                   [fp%vertex(2)%x, fp%vertex(2)%y, fp%vertex(2)%z], &
                                   [fp%vertex(3)%x, fp%vertex(3)%y, fp%vertex(3)%z], e(:, k))
      enddo
      err(1) = max(err(1), abs(glob(k) - ref(k)) / max(abs(ref(k)), 1.e-3_R8P * maxval(abs(ref))))
   enddo
   do k = 1, cut%nz ! 3. multi-piece cells: the face moments of the pieces add up to the section's (+x faces)
      do l = 1, cut%ny
         do p = 1, cut%nx
            if (cut%pieces(p, l, k) < 2) cycle
            fs = 0._R8P
            do q = cut%first_piece(p, l, k), cut%first_piece(p, l, k) + cut%pieces(p, l, k) - 1
               fs = fs + cut%piece_face_moments(:, 2, q)
            enddo
            af = (yy(l + 1) - yy(l)) * (zz(k + 1) - zz(k))
            err(3) = max(err(3), maxval(abs(fs - cut%face_x_moments(:, p + 1, l, k))) / af)
         enddo
      enddo
   enddo
   endsubroutine check3

   subroutine check2(xx, yy)
   !< Check 4 on the block xx x yy.
   real(R8P), intent(in)     :: xx(:), yy(:)
   integer(I4P), allocatable :: e(:,:)
   real(R8P),    allocatable :: glob(:), ref(:)
   integer(I4P)              :: p, k, l, nm, ic(2), q
   real(R8P)                 :: cc(2), h(2), ac

   call c%cut_cells(x=xx, y=yy, cut=c2, moments_order=4_I4P)
   nm = mono_count(2_I4P, 4_I4P) ; e = mono_exponents(2_I4P, 4_I4P)
   allocate(glob(nm), ref(nm)) ; glob = 0._R8P ; ref = 0._R8P
   do p = 1, c2%np
      ic = c2%piece_cell(:, p) ; h = [xx(ic(1) + 1) - xx(ic(1)), yy(ic(2) + 1) - yy(ic(2))] ; ac = product(h)
      cc = [xx(ic(1)), yy(ic(2))] + 0.5_R8P * h
      do k = 1, nm
         do l = 1, nm
            if (any(e(:, l) > e(:, k))) cycle
            glob(k) = glob(k) + product([(binom(e(q, k), e(q, l)) * cc(q) ** (e(q, k) - e(q, l)) * h(q) ** e(q, l), q = 1, 2)]) &
                      * c2%piece_moments(l, p)
         enddo
      enddo
      err(5) = max(err(5), abs(c2%piece_moments(1, p) - c2%piece_area(p)) / ac)
      if (c2%piece_area(p) > 0._R8P) then
         do q = 1, 2
            err(5) = max(err(5), abs(c2%piece_moments(1 + q, p) - c2%piece_area(p) * (c2%piece_centroid(q, p) - cc(q)) / h(q)) / ac)
         enddo
      endif
      err(4) = max(err(4), abs(c2%piece_wall_moments(1, p) - c2%piece_wet_length(p)) / maxval(h))
   enddo
   do k = 1, nm
      ref(k) = (xx(size(xx)) ** (e(1, k) + 1) - xx(1) ** (e(1, k) + 1)) / (e(1, k) + 1) * &
               (yy(size(yy)) ** (e(2, k) + 1) - yy(1) ** (e(2, k) + 1)) / (e(2, k) + 1)
      do p = 1, size(u)
         ref(k) = ref(k) - tri_int([0.5_R8P, 0.5_R8P], [u(p), v(p)], [u(mod(p, size(u)) + 1), v(mod(p, size(u)) + 1)], e(:, k))
      enddo
      err(4) = max(err(4), abs(glob(k) - ref(k)) / max(abs(ref(k)), 1.e-3_R8P * maxval(abs(ref))))
   enddo
   endsubroutine check2

   pure real(R8P) function binom(n_, k_)
   !< Binomial coefficient.
   integer(I4P), intent(in) :: n_, k_
   integer(I4P)             :: j

   binom = 1._R8P
   do j = 1, k_
      binom = binom * (n_ - k_ + j) / j
   enddo
   endfunction binom

   pure real(R8P) function tri_int(a, bb, cc, ex)
   !< int x^ex(1) y^ex(2) over the triangle a bb cc (collapsed 4x4 Gauss, exact for degree <= 5).
   real(R8P),    intent(in) :: a(2), bb(2), cc(2)
   integer(I4P), intent(in) :: ex(2)
   real(R8P)                :: gx(4), gw(4), ss, tt, pt(2), jac
   integer(I4P)             :: i1, i2

   call gauss4(gx, gw)
   jac = abs((bb(1) - a(1)) * (cc(2) - a(2)) - (bb(2) - a(2)) * (cc(1) - a(1)))
   tri_int = 0._R8P
   do i1 = 1, 4
      do i2 = 1, 4
         ss = 0.5_R8P * (gx(i1) + 1._R8P) ; tt = 0.5_R8P * (gx(i2) + 1._R8P)
         pt = a + ss * (bb - a) + ss * tt * (cc - bb)
         tri_int = tri_int + 0.25_R8P * gw(i1) * gw(i2) * ss * jac * pt(1) ** ex(1) * pt(2) ** ex(2)
      enddo
   enddo
   endfunction tri_int

   pure real(R8P) function tet_int(p1, p2, p3, ex)
   !< Signed int x^a y^b z^c over the tetrahedron (0, p1, p2, p3) (collapsed 4x4x4 Gauss, exact for degree <= 4).
   real(R8P),    intent(in) :: p1(3), p2(3), p3(3)
   integer(I4P), intent(in) :: ex(3)
   real(R8P)                :: gx(4), gw(4), r, ss, tt, pt(3), det
   integer(I4P)             :: i1, i2, i3

   call gauss4(gx, gw)
   det = p1(1) * (p2(2) * p3(3) - p2(3) * p3(2)) - p1(2) * (p2(1) * p3(3) - p2(3) * p3(1)) + &
         p1(3) * (p2(1) * p3(2) - p2(2) * p3(1))
   tet_int = 0._R8P
   do i1 = 1, 4
      do i2 = 1, 4
         do i3 = 1, 4
            r = 0.5_R8P * (gx(i1) + 1._R8P) ; ss = 0.5_R8P * (gx(i2) + 1._R8P) ; tt = 0.5_R8P * (gx(i3) + 1._R8P)
            pt = r * (p1 + ss * (p2 - p1) + ss * tt * (p3 - p2)) ! Duffy: x = r (p1 + s (p2 - p1) + s t (p3 - p2))
            tet_int = tet_int + 0.125_R8P * gw(i1) * gw(i2) * gw(i3) * r * r * ss * det * &
                      pt(1) ** ex(1) * pt(2) ** ex(2) * pt(3) ** ex(3)
         enddo
      enddo
   enddo
   endfunction tet_int

   pure subroutine gauss4(gx, gw)
   !< 4-point Gauss-Legendre on [-1, 1] (degree 7: with the Duffy factors r^2 s, exact for the order-4 monomials).
   real(R8P), intent(out) :: gx(4), gw(4)

   gx = [-0.8611363115940526_R8P, -0.3399810435848563_R8P, 0.3399810435848563_R8P, 0.8611363115940526_R8P]
   gw = [0.3478548451374538_R8P, 0.6521451548625461_R8P, 0.6521451548625461_R8P, 0.3478548451374538_R8P]
   endsubroutine gauss4
endprogram fossil_test_cut_moments
