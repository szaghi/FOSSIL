!< FOSSIL test: exact inside flags of cell centres near the wall (issue #69, step 5).
!<
!< A tetrahedron with tilted faces. Along 200 fixed-seed x lines through it, the centres are each rounded crossing
!< parameter t and its neighbours up to 3 ulp on either side (plus two far points); `cell_wall_fractions` gives their
!< inside flags. Reference: the point-in-tetrahedron test in quadruple precision (real128; the four orientations of
!< points with double coordinates, far from exact ties at these offsets). Every flag must match. Before issue #69 the
!< rounded parameter decided whether a crossing was before a centre: centres within round-off of a face landed on
!< either side (negative control on 04e5bbf: mismatches > 0).
!< 2-D twin: a triangle with tilted sides as a `curve_object`, lines along u, `curve%cell_wall_fractions`, reference in
!< quadruple precision.
program fossil_test_cut_centre_exact

use, intrinsic :: iso_fortran_env, only : real128
use fossil,              only : surface_stl_object, curve_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, I8P, R8P
use vecfor,              only : vector_R8P

implicit none

integer,   parameter            :: QP = real128
real(R8P), parameter            :: V(3, 4) = reshape([0._R8P, 0._R8P, 0._R8P, 1._R8P, 0.13_R8P, 0.07_R8P, &
                                                      0.21_R8P, 1._R8P, 0.11_R8P, 0.17_R8P, 0.19_R8P, 1._R8P], [3, 4])
integer(I4P), parameter         :: FACES(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4])
real(R8P), parameter            :: W(2, 3) = reshape([0._R8P, 0._R8P, 1._R8P, 0.17_R8P, 0.23_R8P, 1._R8P], [2, 3])
type(surface_stl_object)        :: tet
type(curve_object)              :: tri
integer(I4P)                    :: nbad2, npts2
real(R8P),    allocatable       :: frac2(:,:,:)
logical,      allocatable       :: ins2(:,:)
type(facet_object), allocatable :: f(:)
real(R8P),    allocatable       :: t(:), x(:), frac(:,:,:,:)
integer(I4P), allocatable       :: fc(:)
logical,      allocatable       :: en(:), ins(:,:,:)
real(R8P)                       :: a, b, y(1), z(1)
integer(I4P)                    :: l, k, j, n, m, nbad, npts
integer(I8P)                    :: seed

allocate(f(4))
do k = 1, 4
   do j = 1, 3
      f(k)%vertex(j) = vector_R8P(V(1, FACES(j, k)), V(2, FACES(j, k)), V(3, FACES(j, k)))
   enddo
enddo
call tet%adopt_facets(facets=f)
if (tet%get_volume() < 0._R8P) error stop 'tetrahedron wound inwards'
seed = 20261106_I8P ; nbad = 0 ; npts = 0
do l = 1, 200
   a = 0.05_R8P + 0.4_R8P * rnd() ; b = 0.05_R8P + 0.4_R8P * rnd()
   call tet%line_crossings(1_I4P, a, b, -10._R8P, 10._R8P, t, fc, en, n)
   if (n == 0) cycle
   allocate(x(0))
   x = [-5._R8P]
   do k = 1, n
      x = [x, around(t(k))]
   enddo
   x = [x, 5._R8P]
   call sort_unique(x)
   m = size(x)
   y = a ; z = b
   allocate(frac(6, m, 1, 1), ins(m, 1, 1))
   call tet%cell_wall_fractions(x, y, z, frac, ins)
   do k = 1, m
      npts = npts + 1
      if (ins(k, 1, 1) .neqv. inside_ref([x(k), a, b])) nbad = nbad + 1
   enddo
   deallocate(x, frac, ins)
enddo
print '(A,I0,A,I0)', '3-D points ', npts, '  mismatches ', nbad

call tri%adopt_loops(u=W(1, :), v=W(2, :), loop_sizes=[3_I4P])
nbad2 = 0 ; npts2 = 0
do l = 1, 200
   a = 0.05_R8P + 0.6_R8P * rnd()
   call tri%line_crossings(1_I4P, a, -10._R8P, 10._R8P, t, fc, en, n)
   if (n == 0) cycle
   x = [-5._R8P]
   do k = 1, n
      x = [x, around(t(k))]
   enddo
   x = [x, 5._R8P]
   call sort_unique(x)
   m = size(x)
   allocate(frac2(4, m, 1), ins2(m, 1))
   call tri%cell_wall_fractions(x, [a], frac2, ins2)
   do k = 1, m
      npts2 = npts2 + 1
      if (ins2(k, 1) .neqv. inside_ref_2d(x(k), a)) nbad2 = nbad2 + 1
   enddo
   deallocate(x, frac2, ins2)
enddo
print '(A,I0,A,I0)', '2-D points ', npts2, '  mismatches ', nbad2
print '(A,L1)', 'Are all tests passed? ', nbad == 0 .and. npts > 1000 .and. nbad2 == 0 .and. npts2 > 1000
contains
   function around(tc) result(r)
   !< tc and its neighbours up to 3 ulp on either side.
   real(R8P), intent(in) :: tc
   real(R8P)             :: r(7)
   integer(I4P)          :: i

   r(4) = tc
   do i = 1, 3
      r(4 + i) = nearest(r(3 + i), 1._R8P) ; r(4 - i) = nearest(r(5 - i), -1._R8P)
   enddo
   endfunction around

   subroutine sort_unique(w)
   !< Sort ascending and drop repeats.
   real(R8P), allocatable, intent(inout) :: w(:)
   real(R8P)                             :: s
   integer(I4P)                          :: i, k, nn

   nn = size(w)
   do i = 2, nn
      s = w(i) ; k = i - 1
      do while (k >= 1)
         if (w(k) <= s) exit
         w(k + 1) = w(k) ; k = k - 1
      enddo
      w(k + 1) = s
   enddo
   k = 1
   do i = 2, nn
      if (w(i) /= w(k)) then
         k = k + 1 ; w(k) = w(i)
      endif
   enddo
   w = w(1:k)
   endsubroutine sort_unique

   logical function inside_ref(p)
   !< p strictly inside the tetrahedron, in quadruple precision: on the inner side of its four outward faces.
   real(R8P), intent(in) :: p(3)
   real(QP)              :: q(3), a1(3), a2(3), a3(3), d
   integer(I4P)          :: k

   inside_ref = .true.
   q = real(p, QP)
   do k = 1, 4
      a1 = real(V(:, FACES(1, k)), QP) - q ; a2 = real(V(:, FACES(2, k)), QP) - q ; a3 = real(V(:, FACES(3, k)), QP) - q
      d = a1(1) * (a2(2) * a3(3) - a2(3) * a3(2)) - a1(2) * (a2(1) * a3(3) - a2(3) * a3(1)) + &
          a1(3) * (a2(1) * a3(2) - a2(2) * a3(1))
      if (d <= 0._QP) inside_ref = .false. ! outward faces: det[a - p, b - p, c - p] > 0 for p inside
   enddo
   endfunction inside_ref

   logical function inside_ref_2d(pu, pv)
   !< (pu, pv) strictly inside the counter-clockwise triangle W, in quadruple precision.
   real(R8P), intent(in) :: pu, pv
   real(QP)              :: d
   integer(I4P)          :: k, k1

   inside_ref_2d = .true.
   do k = 1, 3
      k1 = mod(k, 3) + 1
      d = (real(W(1, k1), QP) - real(W(1, k), QP)) * (real(pv, QP) - real(W(2, k), QP)) - &
          (real(W(2, k1), QP) - real(W(2, k), QP)) * (real(pu, QP) - real(W(1, k), QP))
      if (d <= 0._QP) inside_ref_2d = .false.
   enddo
   endfunction inside_ref_2d

   function rnd() result(r)
   !< Park-Miller uniform in (0, 1).
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_centre_exact
