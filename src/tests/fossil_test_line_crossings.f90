!< FOSSIL test: exact crossings of axis-aligned lines, wall fractions (issue #40 step 1).
program fossil_test_line_crossings
!< FOSSIL test: exact crossings of axis-aligned lines, wall fractions (issue #40 step 1).
!<
!< 3-D (`line_crossings`, `cell_wall_fractions`):
!< - every sign-benchmark body on a 13^3 block of centres covering its box: the parity inside flag
!<   equals the solid-angle test at every centre farther than 1e-9 L from the surface; flags
!<   alternate on every line;
!< - the unit cube on centres k/8, k = -2..10 (lines along faces, through edges, vertices and the
!<   diagonals of the triangulation): every line has an even number of crossings with alternating
!<   flags, the parity is right at every centre off the surface, and the fractions toward x = 0 and
!<   x = 1 are exact;
!< - the 1e-3 plate: two crossings 1e-3 apart on lines through it;
!< - the same block on the octree (brute force) gives the BVH crossings.
!< 2-D: the square on centres k/8 (lines along its sides and through its corners), the square
!< annulus (four crossings per line through the hole), a NACA 0012: parity = exact winding number
!< off the curve, alternating flags, exact fractions on the square. Refusals.

use fossil,                   only : surface_stl_object, curve_object, SIGN_SOLID_ANGLE, AABB_TREE_OCTREE, STATUS_INVALID_INPUT
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s, so
type(curve_object)        :: c
real(R8P),    allocatable :: t(:)
integer(I4P), allocatable :: f(:)
logical,      allocatable :: e(:)
real(R8P)                 :: g(13), u(200), v(200)
integer(I4P)              :: i, k, n, nbad, st
logical                   :: ok(6)
real(R8P), parameter      :: PI = acos(-1._R8P)

call bench_bodies(b)
! 1. every body: parity = solid angle off the surface
ok(1) = .true.
do i = 1, size(b)
   call to_surface(b(i)%soup, s)
   ok(1) = ok(1) .and. parity3(s, b(i)%diag, trim(b(i)%name))
enddo
! 2. grid-aligned unit cube
call to_surface(b(1)%soup, s)
do k = 1, 13
   g(k) = (k - 3) / 8._R8P
enddo
ok(2) = cube_aligned(s, g)
! 3. plate: two crossings 1e-3 apart
call to_surface(b(5)%soup, s)
call s%line_crossings(axis=3_I4P, a=0.3_R8P, b=0.6_R8P, t0=-1._R8P, t1=1._R8P, t=t, facets=f, entering=e, n=n)
ok(3) = n == 2 .and. e(1) .and. .not. e(2) .and. abs(t(1)) < 1.e-15_R8P .and. abs(t(2) - 1.e-3_R8P) < 1.e-15_R8P
print '(A,L1)', 'plate: two crossings 1e-3 apart, entering then leaving: ', ok(3)
! 4. octree = BVH
call to_surface(b(8)%soup, s)
call to_surface(b(8)%soup, so)
call so%analyze(aabb_tree_kind=AABB_TREE_OCTREE)
ok(4) = same_crossings(s, so)
print '(A,L1)', 'L-prism: octree (brute force) crossings = BVH crossings: ', ok(4)
! 5. 2-D
ok(5) = square_aligned()
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
call c%line_crossings(axis=1_I4P, a=3._R8P, t0=-1._R8P, t1=7._R8P, t=t, segments=f, entering=e, n=n)
ok(5) = ok(5) .and. n == 4 .and. all(t == [0._R8P, 1._R8P, 5._R8P, 6._R8P]) .and. all(e .eqv. [.true., .false., .true., .false.])
ok(5) = ok(5) .and. parity2(c, -1._R8P, 7._R8P, 'annulus')
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(5) = ok(5) .and. parity2(c, -0.2_R8P, 1.2_R8P, 'NACA 0012')
! 6. refusals
call s%line_crossings(axis=4_I4P, a=0._R8P, b=0._R8P, t0=0._R8P, t1=1._R8P, t=t, facets=f, entering=e, n=n, status=st)
ok(6) = st == STATUS_INVALID_INPUT
call c%line_crossings(axis=1_I4P, a=0._R8P, t0=1._R8P, t1=0._R8P, t=t, segments=f, entering=e, n=n, status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function parity3(sf, diag, label) result(good)
   !< Parity inside = solid angle at the centres off the surface; flags alternate.
   type(surface_stl_object), intent(in) :: sf    !< Surface.
   real(R8P),                intent(in) :: diag  !< Diagonal.
   character(*),             intent(in) :: label !< Label.
   logical                              :: good  !< Result.
   real(R8P)        :: x(13), y(13), z(13), fr(6, 13, 13, 13), d, lo(3), hi(3) !< Work.
   logical          :: ins(13, 13, 13)  !< Parity flags.
   integer(I4P)     :: i, j, k, nb, nwrong, ncheck !< Work.
   type(vector_R8P) :: p, bmin, bmax    !< Work.

   bmin = sf%get_bmin() ; bmax = sf%get_bmax()
   lo = [bmin%x, bmin%y, bmin%z] ; hi = [bmax%x, bmax%y, bmax%z]
   do i = 1, 13
      x(i) = lo(1) - 0.1_R8P * (hi(1) - lo(1)) + 1.2_R8P * (hi(1) - lo(1)) * (i - 0.5_R8P) / 13._R8P
      y(i) = lo(2) - 0.1_R8P * (hi(2) - lo(2)) + 1.2_R8P * (hi(2) - lo(2)) * (i - 0.47_R8P) / 13._R8P
      z(i) = lo(3) - 0.1_R8P * (hi(3) - lo(3)) + 1.2_R8P * (hi(3) - lo(3)) * (i - 0.53_R8P) / 13._R8P
   enddo
   call sf%cell_wall_fractions(x=x, y=y, z=z, frac=fr, inside=ins, bad_lines=nb)
   nwrong = 0 ; ncheck = 0
   do k = 1, 13
      do j = 1, 13
         do i = 1, 13
            p = vector_R8P(x(i), y(j), z(k))
            call sf%compute_distance(point=p, distance=d, is_square_root=.true.)
            if (d < 1.e-9_R8P * diag) cycle
            ncheck = ncheck + 1
            if (ins(i, j, k) .neqv. sf%is_point_inside(p, sign_algorithm=SIGN_SOLID_ANGLE)) nwrong = nwrong + 1
         enddo
      enddo
   enddo
   good = nb == 0 .and. nwrong == 0 .and. count(ins) > 0
   print '(A,A,I0,A,I0,A,I0,A,I0)', label, ': centres ', ncheck, '  inside ', count(ins), '  parity errors ', nwrong, &
         '  bad lines ', nb
   endfunction parity3

   function cube_aligned(sf, g) result(good)
   !< Unit cube, centres on k/8: even crossings, alternating flags, parity off the surface, fractions.
   type(surface_stl_object), intent(in) :: sf   !< Unit cube [0, 1]^3.
   real(R8P),                intent(in) :: g(:) !< Centres (k / 8).
   logical                              :: good !< Result.
   real(R8P)    :: fr(6, size(g), size(g), size(g)), ex !< Fractions, expected.
   logical      :: ins(size(g), size(g), size(g)), on, truth !< Work.
   integer(I4P) :: i, j, k, ax, nb, nodd, nwrong, nfrac !< Work.

   nodd = 0
   do ax = 1, 3
      do j = 1, size(g)
         do k = 1, size(g)
            call sf%line_crossings(axis=int(ax, I4P), a=g(j), b=g(k), t0=-1._R8P, t1=2._R8P, t=t, facets=f, entering=e, n=n)
            if (mod(n, 2) /= 0) nodd = nodd + 1
         enddo
      enddo
   enddo
   call sf%cell_wall_fractions(x=g, y=g, z=g, frac=fr, inside=ins, bad_lines=nb)
   nwrong = 0 ; nfrac = 0
   do k = 1, size(g)
      do j = 1, size(g)
         do i = 1, size(g)
            on = any([g(i), g(j), g(k)] == 0._R8P .or. [g(i), g(j), g(k)] == 1._R8P)
            truth = all([g(i), g(j), g(k)] > 0._R8P .and. [g(i), g(j), g(k)] < 1._R8P)
            if (.not. on .and. (ins(i, j, k) .neqv. truth)) nwrong = nwrong + 1
            ! along +x from a centre strictly inside the yz square: the wall x = 0 or x = 1
            if (g(j) > 0._R8P .and. g(j) < 1._R8P .and. g(k) > 0._R8P .and. g(k) < 1._R8P .and. i < size(g)) then
               ex = -1._R8P
               if (g(i) <= 0._R8P .and. g(i + 1) >= 0._R8P) ex = (0._R8P - g(i)) / (g(i + 1) - g(i))
               if (g(i) <= 1._R8P .and. g(i + 1) >= 1._R8P .and. g(i) > 0._R8P) ex = (1._R8P - g(i)) / (g(i + 1) - g(i))
               if (abs(fr(2, i, j, k) - ex) > 1.e-14_R8P) nfrac = nfrac + 1
            endif
         enddo
      enddo
   enddo
   good = nodd == 0 .and. nb == 0 .and. nwrong == 0 .and. nfrac == 0
   print '(A,I0,A,I0,A,I0,A,I0)', 'grid-aligned cube: lines with odd crossings ', nodd, '  bad lines ', nb, &
         '  parity errors ', nwrong, '  fraction errors ', nfrac
   endfunction cube_aligned

   function same_crossings(s1, s2) result(same)
   !< Crossings of the two surfaces equal on a set of lines along each axis.
   type(surface_stl_object), intent(in) :: s1, s2 !< Surfaces.
   logical                              :: same   !< Result.
   real(R8P),    allocatable :: t2(:)  !< Work.
   integer(I4P), allocatable :: f2(:)  !< Work.
   logical,      allocatable :: e2(:)  !< Work.
   integer(I4P)              :: ax, j, k, n2 !< Work.

   same = .true.
   do ax = 1, 3
      do j = 0, 8
         do k = 0, 8
            call s1%line_crossings(axis=int(ax, I4P), a=j / 4._R8P, b=k / 8._R8P, t0=-1._R8P, t1=3._R8P, t=t, facets=f, &
                                   entering=e, n=n)
            call s2%line_crossings(axis=int(ax, I4P), a=j / 4._R8P, b=k / 8._R8P, t0=-1._R8P, t1=3._R8P, t=t2, facets=f2, &
                                   entering=e2, n=n2)
            same = same .and. n == n2
            if (n == n2) same = same .and. all(t == t2) .and. all(f == f2) .and. all(e .eqv. e2)
         enddo
      enddo
   enddo
   endfunction same_crossings

   function square_aligned() result(good)
   !< Unit square, centres on k/8: parity off the curve, alternating flags, exact fractions.
   logical      :: good !< Result.
   real(R8P)    :: q(13), fr(4, 13, 13), ex !< Work.
   logical      :: ins(13, 13), on, truth !< Work.
   integer(I4P) :: i, j, kk, nb, nwrong, nfrac !< Work.

   do kk = 1, 13
      q(kk) = (kk - 3) / 8._R8P
   enddo
   call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
   call c%cell_wall_fractions(x=q, y=q, frac=fr, inside=ins, bad_lines=nb)
   nwrong = 0 ; nfrac = 0
   do j = 1, 13
      do i = 1, 13
         on = any([q(i), q(j)] == 0._R8P .or. [q(i), q(j)] == 1._R8P)
         truth = all([q(i), q(j)] > 0._R8P .and. [q(i), q(j)] < 1._R8P)
         if (.not. on .and. (ins(i, j) .neqv. truth)) nwrong = nwrong + 1
         if (q(j) > 0._R8P .and. q(j) < 1._R8P .and. i < 13) then
            ex = -1._R8P
            if (q(i) <= 0._R8P .and. q(i + 1) >= 0._R8P) ex = (0._R8P - q(i)) / (q(i + 1) - q(i))
            if (q(i) <= 1._R8P .and. q(i + 1) >= 1._R8P .and. q(i) > 0._R8P) ex = (1._R8P - q(i)) / (q(i + 1) - q(i))
            if (abs(fr(2, i, j) - ex) > 1.e-14_R8P) nfrac = nfrac + 1
         endif
      enddo
   enddo
   good = nb == 0 .and. nwrong == 0 .and. nfrac == 0
   print '(A,I0,A,I0,A,I0)', 'grid-aligned square: bad lines ', nb, '  parity errors ', nwrong, '  fraction errors ', nfrac
   endfunction square_aligned

   function parity2(cur, lo, hi, label) result(good)
   !< Parity inside = exact winding number on a 41 x 41 block; flags alternate.
   type(curve_object), intent(in) :: cur    !< Curve.
   real(R8P),          intent(in) :: lo, hi !< Block range (both axes).
   character(*),       intent(in) :: label  !< Label.
   logical                        :: good   !< Result.
   real(R8P)    :: q(41), fr(4, 41, 41), d !< Work.
   logical      :: ins(41, 41) !< Work.
   integer(I4P) :: i, j, nb, nwrong !< Work.

   do i = 1, 41
      q(i) = lo + (hi - lo) * (i - 0.5_R8P) / 41._R8P
   enddo
   call cur%cell_wall_fractions(x=q, y=q, frac=fr, inside=ins, bad_lines=nb)
   nwrong = 0
   do j = 1, 41
      do i = 1, 41
         d = abs(cur%distance([q(i), q(j)]))
         if (d < 1.e-12_R8P) cycle
         if (ins(i, j) .neqv. (cur%winding_number_exact([q(i), q(j)]) /= 0)) nwrong = nwrong + 1
      enddo
   enddo
   good = nb == 0 .and. nwrong == 0 .and. count(ins) > 0
   print '(A,A,I0,A,I0,A,I0)', label, ': inside ', count(ins), '  parity errors ', nwrong, '  bad lines ', nb
   endfunction parity2

   pure function thickness(x) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: x !< Chord position.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.126_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_line_crossings
