!< FOSSIL test: signed tag bands on a block of cells (issue #38 step 4).
program fossil_test_tag_block
!< FOSSIL test: signed tag bands on a block of cells (issue #38 step 4).
!<
!< Reference: cut cells from `box_classify_many` (one box per cell), the band by brute force (the
!< minimum over all cut cells of the Chebyshev index distance, capped at band + 1), the sign by
!< `is_point_inside` at the centre of *every* uncut cell (not one per component, so the component
!< argument is checked too). `tag_block` must equal it cell by cell.
!< - 3-D: cube-02 on a uniform 16^3 block straddling a corner, on a stretched 12^3 block, a block
!<   inside (all -(band + 1)), a block outside (all band + 1); bunny on a 12^3 block; band 4.
!< - 2-D: NACA 0012 on a stretched 48 x 20 block, a 200-gon on a uniform 32 x 32 block, a block
!<   inside the 200-gon, a block outside; band 3.
!< - Refusals: non-increasing edges, wrong shape, negative band.
!< - The default sign is the fast winding number (issue #39): every block gives the same tags with
!<   `sign_algorithm=SIGN_FAST_WINDING_NUMBER`.

use fossil, only : surface_stl_object, curve_object, BOX_CUT, STATUS_INVALID_INPUT, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube, bunny
type(curve_object)       :: c
real(R8P)                :: e16(17), s12(13), u(200), v(200), g(33), sx(49), sy(21), bmin(3), bmax(3)
integer(I4P)             :: k, st, st2, st3, t3(2, 2, 2), t2(2, 2)
logical                  :: ok(5)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bbox(cube, bmin, bmax)
! uniform block over the corner region of the cube, stretched block across the whole cube
do k = 1, 17
   e16(k) = bmin(1) + (bmax(1) - bmin(1)) * (-0.3_R8P + 0.8_R8P * (k - 1) / 16._R8P)
enddo
do k = 1, 13
   s12(k) = bmin(1) + (bmax(1) - bmin(1)) * (-0.2_R8P + 1.4_R8P * ((k - 1) / 12._R8P)**1.6_R8P)
enddo
ok(1) = check3(cube, e16, e16, e16, 4, 'cube-02 uniform 16^3') .and. check3(cube, s12, s12, s12, 4, 'cube-02 stretched 12^3')
ok(1) = ok(1) .and. check3(cube, lin(bmin(1), bmax(1), 0.3_R8P, 0.6_R8P, 7), lin(bmin(2), bmax(2), 0.3_R8P, 0.6_R8P, 7), &
                           lin(bmin(3), bmax(3), 0.3_R8P, 0.6_R8P, 7), 4, 'cube-02 block inside', -5)
ok(1) = ok(1) .and. check3(cube, lin(bmin(1), bmax(1), 1.2_R8P, 1.5_R8P, 7), lin(bmin(2), bmax(2), 0._R8P, 1._R8P, 7), &
                           lin(bmin(3), bmax(3), 0._R8P, 1._R8P, 7), 4, 'cube-02 block outside', 5)
call bbox(bunny, bmin, bmax)
ok(2) = check3(bunny, lin(bmin(1), bmax(1), -0.1_R8P, 1.1_R8P, 13), lin(bmin(2), bmax(2), -0.1_R8P, 1.1_R8P, 13), &
               lin(bmin(3), bmax(3), -0.1_R8P, 1.1_R8P, 13), 4, 'bunny 12^3')

! 2-D
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
do k = 1, 49
   sx(k) = -0.2_R8P + 1.4_R8P * ((k - 1) / 48._R8P)**1.3_R8P
enddo
do k = 1, 21
   sy(k) = 0.25_R8P * sinh(2._R8P * (k - 11) / 10._R8P) / sinh(2._R8P)
enddo
ok(3) = check2(c, sx, sy, 3, 'NACA 0012 stretched 48 x 20')
do k = 1, 200
   u(k) = cos(2 * PI * (k - 1) / 200) ; v(k) = sin(2 * PI * (k - 1) / 200)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
do k = 1, 33
   g(k) = -1.3_R8P + 2.6_R8P * (k - 1) / 32._R8P
enddo
ok(4) = check2(c, g, g, 3, '200-gon uniform 32 x 32') .and. &
        check2(c, lin(-1._R8P, 1._R8P, 0.35_R8P, 0.65_R8P, 9), lin(-1._R8P, 1._R8P, 0.35_R8P, 0.65_R8P, 9), 3, &
               '200-gon block inside', -4) .and. &
        check2(c, lin(-1._R8P, 1._R8P, 1.2_R8P, 1.6_R8P, 9), lin(-1._R8P, 1._R8P, 0._R8P, 1._R8P, 9), 3, &
               '200-gon block outside', 4)

call cube%tag_block(x=[0._R8P, 1._R8P, 1._R8P], y=[0._R8P, 1._R8P, 2._R8P], z=[0._R8P, 1._R8P, 2._R8P], band=2, tags=t3, &
                    status=st)
call cube%tag_block(x=[0._R8P, 1._R8P], y=[0._R8P, 1._R8P, 2._R8P], z=[0._R8P, 1._R8P, 2._R8P], band=2, tags=t3, status=st2)
call c%tag_block(x=[0._R8P, 1._R8P, 2._R8P], y=[0._R8P, 1._R8P, 2._R8P], band=-1, tags=t2, status=st3)
ok(5) = st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT .and. st3 == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check3(s, x, y, z, band, label, uniform) result(good)
   !< tag_block vs the brute-force reference; optionally every tag must equal `uniform`.
   type(surface_stl_object), intent(in)           :: s          !< Surface.
   real(R8P),                intent(in)           :: x(:), y(:), z(:) !< Cell edges.
   integer,                  intent(in)           :: band       !< Band.
   character(*),             intent(in)           :: label      !< Label.
   integer,                  intent(in), optional :: uniform    !< Expected uniform tag.
   logical                                        :: good       !< Result.
   integer(I4P), allocatable :: tags(:,:,:), ref(:,:,:), codes(:), cuts(:,:), tf(:,:,:) !< Work.
   real(R8P),    allocatable :: lo(:,:), hi(:,:) !< Cell boxes.
   integer(I4P)              :: nx, ny, nz, i, j, k, m, d, ncut, st !< Work.

   nx = size(x) - 1 ; ny = size(y) - 1 ; nz = size(z) - 1
   allocate(tags(nx, ny, nz), ref(nx, ny, nz), codes(nx * ny * nz), lo(3, nx * ny * nz), hi(3, nx * ny * nz))
   call s%tag_block(x=x, y=y, z=z, band=band, tags=tags, status=st)
   allocate(tf(nx, ny, nz)) ! the default sign is the fast winding number (issue #39)
   call s%tag_block(x=x, y=y, z=z, band=band, tags=tf, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   if (any(tf /= tags)) st = -1
   m = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            m = m + 1 ; lo(:, m) = [x(i), y(j), z(k)] ; hi(:, m) = [x(i+1), y(j+1), z(k+1)]
         enddo
      enddo
   enddo
   call s%box_classify_many(lo=lo, hi=hi, codes=codes)
   ncut = count(codes == BOX_CUT)
   allocate(cuts(3, ncut))
   m = 0 ; ncut = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            m = m + 1
            if (codes(m) == BOX_CUT) then
               ncut = ncut + 1 ; cuts(:, ncut) = [i, j, k]
            endif
         enddo
      enddo
   enddo
   m = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            m = m + 1
            if (codes(m) == BOX_CUT) then
               ref(i, j, k) = 0
               cycle
            endif
            d = band + 1
            if (ncut > 0) d = min(d, minval(max(abs(cuts(1, :) - i), abs(cuts(2, :) - j), abs(cuts(3, :) - k))))
            if (s%is_point_inside(vector_R8P(0.5_R8P * (x(i) + x(i+1)), 0.5_R8P * (y(j) + y(j+1)), &
                                             0.5_R8P * (z(k) + z(k+1))))) d = -d
            ref(i, j, k) = d
         enddo
      enddo
   enddo
   good = st == 0 .and. all(tags == ref)
   if (present(uniform)) good = good .and. all(tags == uniform)
   print '(A,I0,A,I0,A,I0,A,I0,A,L1)', label//': cut ', ncut, '  inside ', count(tags < 0), '  outside ', count(tags > 0), &
         '  differences ', count(tags /= ref), '  ok ', good
   endfunction check3

   function check2(cur, x, y, band, label, uniform) result(good)
   !< Curve tag_block vs the brute-force reference.
   type(curve_object), intent(in)           :: cur       !< Curve.
   real(R8P),          intent(in)           :: x(:), y(:) !< Cell edges.
   integer,            intent(in)           :: band      !< Band.
   character(*),       intent(in)           :: label     !< Label.
   integer,            intent(in), optional :: uniform   !< Expected uniform tag.
   logical                                  :: good      !< Result.
   integer(I4P), allocatable :: tags(:,:), ref(:,:), codes(:), cuts(:,:), tf(:,:) !< Work.
   real(R8P),    allocatable :: lo(:,:), hi(:,:) !< Cell boxes.
   integer(I4P)              :: nx, ny, i, j, m, d, ncut, st !< Work.

   nx = size(x) - 1 ; ny = size(y) - 1
   allocate(tags(nx, ny), ref(nx, ny), codes(nx * ny), lo(2, nx * ny), hi(2, nx * ny))
   call cur%tag_block(x=x, y=y, band=band, tags=tags, status=st)
   allocate(tf(nx, ny)) ! the default sign is the fast winding number (issue #39)
   call cur%tag_block(x=x, y=y, band=band, tags=tf, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   if (any(tf /= tags)) st = -1
   m = 0
   do j = 1, ny
      do i = 1, nx
         m = m + 1 ; lo(:, m) = [x(i), y(j)] ; hi(:, m) = [x(i+1), y(j+1)]
      enddo
   enddo
   call cur%box_classify_many(lo=lo, hi=hi, codes=codes)
   ncut = count(codes == BOX_CUT)
   allocate(cuts(2, ncut))
   m = 0 ; ncut = 0
   do j = 1, ny
      do i = 1, nx
         m = m + 1
         if (codes(m) == BOX_CUT) then
            ncut = ncut + 1 ; cuts(:, ncut) = [i, j]
         endif
      enddo
   enddo
   m = 0
   do j = 1, ny
      do i = 1, nx
         m = m + 1
         if (codes(m) == BOX_CUT) then
            ref(i, j) = 0
            cycle
         endif
         d = band + 1
         if (ncut > 0) d = min(d, minval(max(abs(cuts(1, :) - i), abs(cuts(2, :) - j))))
         if (cur%is_point_inside([0.5_R8P * (x(i) + x(i+1)), 0.5_R8P * (y(j) + y(j+1))])) d = -d
         ref(i, j) = d
      enddo
   enddo
   good = st == 0 .and. all(tags == ref)
   if (present(uniform)) good = good .and. all(tags == uniform)
   print '(A,I0,A,I0,A,I0,A,I0,A,L1)', label//': cut ', ncut, '  inside ', count(tags < 0), '  outside ', count(tags > 0), &
         '  differences ', count(tags /= ref), '  ok ', good
   endfunction check2

   pure function lin(a, b, f0, f1, n) result(e)
   !< n edges from a + f0 (b - a) to a + f1 (b - a).
   real(R8P), intent(in) :: a, b, f0, f1 !< Range, fractions.
   integer,   intent(in) :: n            !< Edges.
   real(R8P)             :: e(n)         !< Edges.
   integer               :: k            !< Counter.

   do k = 1, n
      e(k) = a + (b - a) * (f0 + (f1 - f0) * (k - 1) / real(n - 1, R8P))
   enddo
   endfunction lin

   subroutine bbox(s, bmin, bmax)
   !< Bounding box of a surface as arrays.
   type(surface_stl_object), intent(in)  :: s                !< Surface.
   real(R8P),                intent(out) :: bmin(3), bmax(3) !< Box.
   type(vector_R8P)                      :: p                !< Work.

   p = s%get_bmin() ; bmin = [p%x, p%y, p%z]
   p = s%get_bmax() ; bmax = [p%x, p%y, p%z]
   endsubroutine bbox

   pure function thickness(x) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: x !< Chord position.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.126_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_tag_block
