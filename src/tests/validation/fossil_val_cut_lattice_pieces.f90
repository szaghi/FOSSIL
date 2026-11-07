!< FOSSIL validation: the 3-D cut-cell pieces of lattice bodies (issue #70, step 1).
program fossil_val_cut_lattice_pieces
!< FOSSIL validation: the 3-D cut-cell pieces of lattice bodies (issue #70, step 1).
!<
!< VALIDATION version (issue #45 tiers): 42 bodies on three grids, 1.57 million cells (75 s in the debug coverage
!< build). The regression subset is src/tests/fossil_test_cut_lattice_pieces.f90.
!< Bodies whose vertices sit on grid nodes, edges on grid lines and facets on grid planes: cube-01, an octahedron
!< with vertices on the nodes, and 40 fixed-seed tetrahedra with vertices on the dyadic lattice of step 1/8, on
!< grids of step 1/8 and 1/4 (nodes on the lattice) and 1/16. Before #70 the piercings of a cell edge were the
!< clipped polygon vertices on it: a mesh edge or vertex on the edge gave one per incident facet, the corner states
!< did not alternate, and up to 38 cells of a tetrahedron were untrusted (pieces = -1): 947 of these 1.57 million
!< cells. Now the corner states are the exact node parities and the crossings of each edge the exact ones under the
!< global perturbation (#69), the on-plane flags come from the coordinates, and a facet lying in or ending on a face
!< plane is read on the perturbed face; a multi-piece cell whose split does not close is not trusted.
!< Checks, every cell of every body and grid:
!< 1. alpha in [0, 1] up to 1e-14 (round-off of the face sums: 1 + 4e-16 seen), and at most NFLAG untrusted cells
!<    (pieces = -1): 19 measured, all on tetrahedra at h <= 1/8 where a facet touches the cell along a segment of a
!<    cell edge only (its clipped polygon is empty, so its trace on the face is missing and the split cannot close);
!< 2. closure per piece (1e-13; 3.5e-14 measured, as on base): the fluid boundary of a piece has a zero normal integral,
!<    (f+ - f-) A - wet normal = 0 per axis (relative to the cell face area);
!< 3. the pieces of a cell add up to the cell volume, alpha |cell| (1e-13);
!< 4. the fluid volume of the block is the block volume minus the body volume (1e-14).

use fossil,              only : surface_stl_object
use fossil_cut_cells_3d, only : cut_cells_3d_t, cut_cells_3d
use fossil_facet_object, only : facet_object
use penf,                only : I4P, I8P, R8P
use vecfor,              only : vector_R8P

implicit none

type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
real(R8P)                       :: v(3, 4), vol, emax(4)
integer(I4P)                    :: k, i, nbad, nflag, ncells, nmulti
integer(I4P), parameter         :: NFLAG_MAX = 19 !< Untrusted cells measured (see the header).
integer(I8P)                    :: seed
logical                         :: ok

nbad = 0 ; nflag = 0 ; ncells = 0 ; nmulti = 0 ; emax = 0._R8P
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call check_grids(s)
call octahedron(s)
call check_grids(s)
seed = 70_I8P
do k = 1, 40
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(s, v)
   call check_grids(s)
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'cells ', ncells, '  alpha outside [0, 1]: ', nbad, '  untrusted: ', nflag, &
      '  multi-piece cells ', nmulti
print '(A,4ES10.2)', 'max closure / cell-sum / block-volume errors (relative): ', emax(1:3)
ok = nbad == 0 .and. nflag <= NFLAG_MAX .and. emax(1) <= 1.e-13_R8P .and. emax(2) <= 1.e-13_R8P .and. emax(3) <= 1.e-14_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check_grids(sf)
   !< Check the body on the three grids.
   type(surface_stl_object), intent(in) :: sf
   integer(I4P)                         :: n, ig
   real(R8P)                            :: h
   real(R8P), allocatable               :: g(:)
   real(R8P), parameter                 :: STEPS(3) = [0.125_R8P, 0.25_R8P, 0.0625_R8P]

   do ig = 1, 3
      h = STEPS(ig)
      n = nint(2._R8P / h) + 1
      g = [(-0.5_R8P + h * (i - 1), i = 1, n)]
      call check(sf, g, h)
   enddo
   endsubroutine check_grids

   subroutine check(sf, g, h)
   !< Checks 1-4 on the grid g x g x g (uniform step h).
   type(surface_stl_object), intent(in) :: sf
   real(R8P),                intent(in) :: g(:), h
   type(cut_cells_3d_t)                 :: cut
   integer(I4P)                         :: a, b, c, p, n
   real(R8P)                            :: fl, sp, cl(3), hv, ha
   real(R8P), parameter                 :: TOL = 1.e-14_R8P !< alpha round-off (1 + 4e-16 seen: the face sums).

   call cut_cells_3d(sf, g, g, g, cut)
   n = size(g) - 1 ; hv = h ** 3 ; ha = h ** 2 ; fl = 0._R8P
   do c = 1, n
      do b = 1, n
         do a = 1, n
            ncells = ncells + 1
            if (cut%alpha(a, b, c) < -TOL .or. cut%alpha(a, b, c) > 1._R8P + TOL) then
               nbad = nbad + 1 ; cycle
            endif
            fl = fl + cut%alpha(a, b, c) * hv
            if (cut%pieces(a, b, c) < 0) then
               nflag = nflag + 1 ; cycle
            endif
            if (cut%pieces(a, b, c) > 1) nmulti = nmulti + 1
            sp = 0._R8P
            do p = cut%first_piece(a, b, c), cut%first_piece(a, b, c) + cut%pieces(a, b, c) - 1
               sp = sp + cut%piece_volume(p)
               cl = [cut%piece_face(2, p) - cut%piece_face(1, p), cut%piece_face(4, p) - cut%piece_face(3, p), &
                     cut%piece_face(6, p) - cut%piece_face(5, p)] * ha - cut%piece_wet_normal(:, p)
               emax(1) = max(emax(1), maxval(abs(cl)) / ha)
            enddo
            emax(2) = max(emax(2), abs(sp - cut%alpha(a, b, c) * hv) / hv)
         enddo
      enddo
   enddo
   emax(3) = max(emax(3), abs(fl - ((g(n + 1) - g(1)) ** 3 - sf%get_volume())) / (g(n + 1) - g(1)) ** 3)
   endsubroutine check

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
endprogram fossil_val_cut_lattice_pieces
