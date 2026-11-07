!< FOSSIL test: quadrature rules on the cut-cell pieces (issue #70, step 9).
program fossil_test_cut_quadrature
!< FOSSIL test: quadrature rules on the cut-cell pieces (issue #70, step 9).
!<
!< Volume rules: signed cones from the cell centre over the oriented fluid boundary of each piece (the wall fans, and on
!< every face the triangles from the face centre: the plane section's boundary for a one-piece cell, the face arcs of the
!< pieces otherwise), 36 points per cone; fluid cells the 3^d Gauss rule. Wall rules: the 7-point rule on the wall fans
!< (3-D), 3-point Gauss on the wall segments (2-D). Bodies: the L-prism on a generic and on a grid-aligned grid,
!< the thin plate, 4 lattice tetrahedra (3-D); a circle of 64 segments (2-D).
!< 1. every rule with quad_status 0 integrates the monomials of degree <= 4 to the piece's moments (1e-13 of the cell),
!<    every wall rule to its wall moments (1e-13 of the cell size squared, or length in 2-D); at most NSTATUS pieces with
!<    quad_status 1: measured 1 of 5028, a two-piece cell of a lattice tetrahedron whose face arcs (from the clipped
!<    traces, as its face areas) do not close; a one-piece cell takes the plane sections and is always exact;
!< 2. the wall weights are positive; the volume weights of some pieces are negative (cones of a piece that is not
!<    star-shaped about the cell centre): counted and reported, not an error.

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_cut_moments,       only : mono_exponents
use fossil_facet_object,      only : facet_object
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable       :: b(:)
type(surface_stl_object)        :: s
type(curve_object)              :: c
type(facet_object), allocatable :: f(:)
real(R8P), allocatable          :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                       :: err(2), vv(3, 4), vol
integer(I4P)                    :: i, k, n, nbad, nneg, npieces, nstat
integer(I4P), parameter         :: NSTATUS = 1 !< Pieces with an inexact rule (measured, see the header).
integer(I8P)                    :: seed
logical                         :: ok

err = 0._R8P ; nbad = 0 ; nneg = 0 ; npieces = 0 ; nstat = 0
call bench_bodies(b)
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [-0.2_R8P, 0.45_R8P, 1.3_R8P]
call check3(x, y, z)
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call check3(x, y, z) ! aligned
call to_surface(b(5)%soup, s) ! plate
x = [(-0.15_R8P + 0.1_R8P * i, i = 0, 13)] ; z = [-0.3_R8P, -0.1_R8P, 0.1_R8P, 0.3_R8P]
call check3(x, x, z)
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)]
do k = 1, 4
   do
      vv = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(vv(:, 2) - vv(:, 1), cross(vv(:, 3) - vv(:, 1), vv(:, 4) - vv(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) vv(:, [3, 4]) = vv(:, [4, 3])
   call tetra(s, vv)
   call check3(x, x, x)
enddo
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call check2(x, y)
print '(A,I0,A,I0,A,I0,A,I0)', 'pieces ', npieces, '  bad ', nbad, '  pieces with negative volume weights ', nneg, &
      '  inexact (quad_status 1) ', nstat
print '(A,2ES10.2)', 'max errors (volume rules, wall rules): ', err
ok = nbad == 0 .and. nstat <= NSTATUS .and. npieces > 500 .and. all(err <= 1.e-13_R8P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check3(xx, yy, zz)
   !< Checks 1-2 on the block xx x yy x zz.
   real(R8P), intent(in)     :: xx(:), yy(:), zz(:)
   type(cut_cells_3d_t)      :: cut
   integer(I4P), allocatable :: e(:,:)
   integer(I4P)              :: p, kk, q, ic(3)
   real(R8P)                 :: h(3), cc(3), sv, sw, pt(3)

   call cut_cells_3d(s, xx, yy, zz, cut, moments_order=4_I4P, quadrature=.true.)
   e = mono_exponents(3_I4P, 4_I4P)
   npieces = npieces + cut%np ; nstat = nstat + count(cut%quad_status /= 0)
   if (any(cut%wall_quad_weight < 0._R8P)) nbad = nbad + 1
   do p = 1, cut%np
      ic = cut%piece_cell(:, p)
      h = [xx(ic(1) + 1) - xx(ic(1)), yy(ic(2) + 1) - yy(ic(2)), zz(ic(3) + 1) - zz(ic(3))]
      cc = [xx(ic(1)), yy(ic(2)), zz(ic(3))] + 0.5_R8P * h
      if (any(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1) < 0._R8P)) nneg = nneg + 1
      do kk = 1, size(e, 2)
         sv = 0._R8P ; sw = 0._R8P
         do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
            pt = (cut%quad_point(:, q) - cc) / h ; sv = sv + cut%quad_weight(q) * product(pt ** e(:, kk))
         enddo
         do q = cut%wall_quad_start(p), cut%wall_quad_start(p + 1) - 1
            pt = (cut%wall_quad_point(:, q) - cc) / h ; sw = sw + cut%wall_quad_weight(q) * product(pt ** e(:, kk))
         enddo
         if (cut%quad_status(p) == 0) err(1) = max(err(1), abs(sv - cut%piece_moments(kk, p)) / product(h))
         err(2) = max(err(2), abs(sw - cut%piece_wall_moments(kk, p)) / maxval(h) ** 2)
      enddo
   enddo
   endsubroutine check3

   subroutine check2(xx, yy)
   !< Checks 1-2 on the 2-D block xx x yy.
   real(R8P), intent(in)     :: xx(:), yy(:)
   type(cut_cells_2d_t)      :: cut
   integer(I4P), allocatable :: e(:,:)
   integer(I4P)              :: p, kk, q
   real(R8P)                 :: h(2), cc(2), sv, sw, pt(2)

   call c%cut_cells(x=xx, y=yy, cut=cut, moments_order=4_I4P, quadrature=.true.)
   e = mono_exponents(2_I4P, 4_I4P)
   npieces = npieces + cut%np
   if (any(cut%wall_quad_weight < 0._R8P)) nbad = nbad + 1
   do p = 1, cut%np
      h = [xx(cut%piece_cell(1, p) + 1) - xx(cut%piece_cell(1, p)), yy(cut%piece_cell(2, p) + 1) - yy(cut%piece_cell(2, p))]
      cc = [xx(cut%piece_cell(1, p)), yy(cut%piece_cell(2, p))] + 0.5_R8P * h
      if (any(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1) < 0._R8P)) nneg = nneg + 1
      do kk = 1, size(e, 2)
         sv = 0._R8P ; sw = 0._R8P
         do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
            pt = (cut%quad_point(:, q) - cc) / h ; sv = sv + cut%quad_weight(q) * product(pt ** e(:, kk))
         enddo
         do q = cut%wall_quad_start(p), cut%wall_quad_start(p + 1) - 1
            pt = (cut%wall_quad_point(:, q) - cc) / h ; sw = sw + cut%wall_quad_weight(q) * product(pt ** e(:, kk))
         enddo
         err(1) = max(err(1), abs(sv - cut%piece_moments(kk, p)) / product(h))
         err(2) = max(err(2), abs(sw - cut%piece_wall_moments(kk, p)) / maxval(h))
      enddo
   enddo
   endsubroutine check2

   subroutine tetra(sf, vt)
   !< Tetrahedron vt (positive orientation) with outward facets.
   type(surface_stl_object), intent(inout) :: sf
   real(R8P),                intent(in)    :: vt(3, 4)
   integer(I4P), parameter                 :: T(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4])
   integer(I4P)                            :: q, j

   if (allocated(f)) deallocate(f)
   allocate(f(4))
   do q = 1, 4
      do j = 1, 3
         f(q)%vertex(j) = vector_R8P(vt(1, T(j, q)), vt(2, T(j, q)), vt(3, T(j, q)))
      enddo
   enddo
   call sf%adopt_facets(facets=f)
   endsubroutine tetra

   pure function cross(a, cb) result(w)
   real(R8P), intent(in) :: a(3), cb(3)
   real(R8P)             :: w(3)

   w = [a(2) * cb(3) - a(3) * cb(2), a(3) * cb(1) - a(1) * cb(3), a(1) * cb(2) - a(2) * cb(1)]
   endfunction cross

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_quadrature
