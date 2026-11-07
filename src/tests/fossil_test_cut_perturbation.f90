!< FOSSIL test: one global symbolic perturbation for the exact line queries (issue #69, step 6).
!<
!< Every grid node P is taken at P + (eps1, eps2, eps3), eps1 >> eps2 >> eps3, and its inside flag is the parity of
!< the crossings before it along its x, y and z lines (`line_crossings` + `line_crossing_before`). With one global
!< perturbation the three parities are the point-in-solid test of the same perturbed point, so they agree at every
!< node, on the surface included. Bodies with facets on grid planes, edges on grid lines and vertices on grid nodes:
!< 1. cube-01 on a grid of step 0.25 over [-0.5, 1.5];
!< 2. an octahedron with vertices on the nodes (0.5 +- 0.5, 0.5, 0.5), ...;
!< 3. 20 fixed-seed tetrahedra with vertices on a dyadic lattice (step 1/8), positive volume.
!< 2-D twin: node flags of a curve along u and along v (`line_crossings` + `line_crossing_before_2d`) for a square, a
!< diamond and 20 lattice triangles. Negative control (before step 6): the per-axis cyclic tie-break disagrees.
program fossil_test_cut_perturbation

use fossil,              only : surface_stl_object, curve_object
use fossil_facet_object, only : facet_object
use fossil_line_query,   only : line_crossing_before, line_crossing_before_2d
use penf,                only : I4P, I8P, R8P
use vecfor,              only : vector_R8P

implicit none

type(surface_stl_object)        :: s
type(curve_object)              :: c
type(facet_object), allocatable :: f(:)
real(R8P)                       :: g(17), v(3, 4), w(2, 3), vol
integer(I4P)                    :: k, nbad, nnode, nbad2, nnode2, i
integer(I8P)                    :: seed
logical                         :: ok

g = [(-0.5_R8P + 0.125_R8P * (i - 1), i = 1, 17)]
nbad = 0 ; nnode = 0
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call check3(s, g)
call octahedron(s)
call check3(s, g)
seed = 69_I8P
do k = 1, 20
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(s, v)
   call check3(s, g)
enddo
print '(A,I0,A,I0)', '3-D nodes ', nnode, '  disagreements ', nbad

nbad2 = 0 ; nnode2 = 0
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call check2(c, g)
call c%adopt_loops(u=[0.5_R8P, 1._R8P, 0.5_R8P, 0._R8P], v=[0._R8P, 0.5_R8P, 1._R8P, 0.5_R8P], loop_sizes=[4_I4P])
call check2(c, g)
do k = 1, 20
   do
      w = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 6)], [2, 3])
      vol = (w(1, 2) - w(1, 1)) * (w(2, 3) - w(2, 1)) - (w(2, 2) - w(2, 1)) * (w(1, 3) - w(1, 1))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) w(:, [2, 3]) = w(:, [3, 2])
   call c%adopt_loops(u=w(1, :), v=w(2, :), loop_sizes=[3_I4P])
   call check2(c, g)
enddo
print '(A,I0,A,I0)', '2-D nodes ', nnode2, '  disagreements ', nbad2
ok = nbad == 0 .and. nbad2 == 0 .and. nnode > 0 .and. nnode2 > 0
print '(A,L1)', 'Are all tests passed? ', ok
contains
   subroutine check3(sf, gg)
   !< Parities of every node of the grid gg^3 along x, y, z: count the nodes where they disagree.
   type(surface_stl_object), intent(in) :: sf
   real(R8P),                intent(in) :: gg(:)
   integer(I4P)                         :: i1, i2, i3
   logical                              :: px, py, pz

   do i3 = 1, size(gg) ; do i2 = 1, size(gg) ; do i1 = 1, size(gg)
      px = parity3(sf, 1_I4P, gg(i2), gg(i3), gg(i1))
      py = parity3(sf, 2_I4P, gg(i3), gg(i1), gg(i2))
      pz = parity3(sf, 3_I4P, gg(i1), gg(i2), gg(i3))
      nnode = nnode + 1
      if ((px .neqv. py) .or. (px .neqv. pz)) nbad = nbad + 1
   enddo ; enddo ; enddo
   endsubroutine check3

   logical function parity3(sf, ax, a, b, cc)
   !< Parity of the crossings before the point at cc on the line along ax at (a, b).
   type(surface_stl_object), intent(in) :: sf
   integer(I4P),             intent(in) :: ax
   real(R8P),                intent(in) :: a, b, cc
   real(R8P),    allocatable            :: t(:)
   integer(I4P), allocatable            :: fc(:)
   logical,      allocatable            :: en(:)
   integer(I4P)                         :: n, q, m

   call sf%line_crossings(ax, a, b, -huge(1._R8P), huge(1._R8P), t, fc, en, n)
   m = 0
   do q = 1, n
      if (line_crossing_before(reshape(sf%tri_vertices(:, fc(q)), [3, 3]), ax, a, b, cc, en(q))) m = m + 1
   enddo
   parity3 = mod(m, 2) == 1
   endfunction parity3

   subroutine check2(cv, gg)
   !< Parities of every node of the grid gg^2 along u and v.
   type(curve_object), intent(in) :: cv
   real(R8P),          intent(in) :: gg(:)
   integer(I4P)                   :: i1, i2

   do i2 = 1, size(gg) ; do i1 = 1, size(gg)
      nnode2 = nnode2 + 1
      if (parity2(cv, 1_I4P, gg(i2), gg(i1)) .neqv. parity2(cv, 2_I4P, gg(i1), gg(i2))) nbad2 = nbad2 + 1
   enddo ; enddo
   endsubroutine check2

   logical function parity2(cv, ax, a, cc)
   !< Parity of the crossings before the point at cc on the line along ax at a.
   type(curve_object), intent(in) :: cv
   integer(I4P),       intent(in) :: ax
   real(R8P),          intent(in) :: a, cc
   real(R8P),    allocatable      :: t(:)
   integer(I4P), allocatable      :: sg(:)
   logical,      allocatable      :: en(:)
   integer(I4P)                   :: n, q, m
   real(R8P)                      :: p(2), r(2)

   call cv%line_crossings(ax, a, -huge(1._R8P), huge(1._R8P), t, sg, en, n)
   m = 0
   do q = 1, n
      p = cv%get_vertex(sg(q)) ; r = cv%get_vertex(cv%get_next(sg(q)))
      if (line_crossing_before_2d(p(1), p(2), r(1), r(2), ax, a, cc, en(q))) m = m + 1
   enddo
   parity2 = mod(m, 2) == 1
   endfunction parity2

   subroutine octahedron(sf)
   !< Octahedron with vertices on grid nodes, outward facets.
   type(surface_stl_object), intent(inout) :: sf
   real(R8P)                               :: o(3, 6)
   integer(I4P), parameter                 :: T(3, 8) = reshape([1, 3, 5, 3, 2, 5, 2, 4, 5, 4, 1, 5, &
                                                                 3, 1, 6, 2, 3, 6, 4, 2, 6, 1, 4, 6], [3, 8])
   integer(I4P)                            :: q, j

   o = reshape([1._R8P, 0.5_R8P, 0.5_R8P, 0._R8P, 0.5_R8P, 0.5_R8P, 0.5_R8P, 1._R8P, 0.5_R8P, &
                0.5_R8P, 0._R8P, 0.5_R8P, 0.5_R8P, 0.5_R8P, 1._R8P, 0.5_R8P, 0.5_R8P, 0._R8P], [3, 6])
   if (allocated(f)) deallocate(f)
   allocate(f(8))
   do q = 1, 8
      do j = 1, 3
         f(q)%vertex(j) = vector_R8P(o(1, T(j, q)), o(2, T(j, q)), o(3, T(j, q)))
      enddo
   enddo
   call sf%adopt_facets(facets=f)
   if (sf%get_volume() < 0._R8P) error stop 'octahedron wound inwards'
   endsubroutine octahedron

   subroutine tetra(sf, vv)
   !< Tetrahedron vv (positive orientation) with outward facets.
   type(surface_stl_object), intent(inout) :: sf
   real(R8P),                intent(in)    :: vv(3, 4)
   integer(I4P), parameter                 :: T(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4])
   integer(I4P)                            :: q, j

   if (allocated(f)) deallocate(f)
   allocate(f(4))
   do q = 1, 4
      do j = 1, 3
         f(q)%vertex(j) = vector_R8P(vv(1, T(j, q)), vv(2, T(j, q)), vv(3, T(j, q)))
      enddo
   enddo
   call sf%adopt_facets(facets=f)
   endsubroutine tetra

   pure function cross(x, y) result(z)
   real(R8P), intent(in) :: x(3), y(3)
   real(R8P)             :: z(3)

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_perturbation
