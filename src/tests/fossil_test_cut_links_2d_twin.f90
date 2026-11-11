!< FOSSIL test: the 2-D piece links on the device equal the host's bit for bit (issue #100 step 5).
program fossil_test_cut_links_2d_twin
!< FOSSIL test: the 2-D piece links on the device equal the host's bit for bit (issue #100 step 5).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, offset and
!< prefix-sum bugs of the device stages, not device codegen (the curve bench does that on the GPU, `cut links 2-D device
!< vs serial mismatches`).
!< `cut_links_2d_device` against `cut_links_2d` (`cut_links_mismatches`: np, unmatched, start, piece, face, area
!< bitwise) on a host cut and on a resident one (`cut_cells_device(..., resident=.true.)`), the links copied back and
!< resident (then `update_host`; the next call releases them); a copy of resident links is not resident. Curves: an
!< ellipse, a NACA 0012, a disc with a hole, 50 circles, a thin plate (two pieces per cut cell, links between pieces of
!< one edge) on generic and grid-aligned grids; a cut without intervals (no links).
use fossil,                   only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_links_2d_device
use fossil_cut_cells_compare, only : cut_links_mismatches
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter   :: PI = acos(-1._R8P)
type(curve_object)     :: c
type(cut_cells_2d_t)   :: cut
type(cut_links_t)      :: lk, lb
real(R8P), allocatable :: u(:), v(:), x(:), y(:)
integer(I4P)           :: i, k, m, total, nlinks
integer(I4P)           :: nl(50)
logical                :: ok

ok = .true. ; total = 0_I4P ; nlinks = 0_I4P
u = [(0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 256), k = 1, 256)]
v = [(0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 256), k = 1, 256)]
call c%adopt_loops(u=u, v=v, loop_sizes=[256_I4P])
call grid(48_I4P) ; call case('ellipse')
call naca0012(128_I4P, u, v)
u = 0.15_R8P + 0.7_R8P * u ; v = 0.5_R8P + 0.7_R8P * v
call c%adopt_loops(u=u, v=v, loop_sizes=[256_I4P])
call grid(64_I4P) ; call case('NACA 0012')
deallocate(u, v) ; allocate(u(300), v(300))
do k = 1, 200
   u(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 200) ; v(k) = 0.5_R8P + 0.4_R8P * sin(2._R8P * PI * (k - 1) / 200)
enddo
do k = 1, 100
   u(200 + k) = 0.5_R8P + 0.12_R8P * cos(-2._R8P * PI * (k - 1) / 100)
   v(200 + k) = 0.5_R8P + 0.12_R8P * sin(-2._R8P * PI * (k - 1) / 100)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P, 100_I4P])
call grid(40_I4P) ; call case('disc with a hole')
deallocate(u, v) ; allocate(u(800), v(800))
do m = 1, 50
   do k = 1, 16
      u(16 * (m - 1) + k) = 0.1_R8P + 0.8_R8P * mod(m - 1, 10) / 9._R8P + 0.03_R8P * cos(2._R8P * PI * (k - 1) / 16)
      v(16 * (m - 1) + k) = 0.1_R8P + 0.8_R8P * ((m - 1) / 10) / 4._R8P + 0.03_R8P * sin(2._R8P * PI * (k - 1) / 16)
   enddo
enddo
nl = 16_I4P
call c%adopt_loops(u=u, v=v, loop_sizes=nl)
call grid(64_I4P) ; call case('50 circles')
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P], v=[0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P], &
                   loop_sizes=[4_I4P])
x = [(0.05_R8P + 0.1_R8P * (i - 1), i = 1, 10)] ; y = x
call case('thin plate, grid-aligned')
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call case('thin plate, generic grid')
call c%cut_cells(x=x, y=y, cut=cut, intervals=.false.) ! no intervals: the host's links on that cut (none)
call cut_links_2d(cut, lk)
call cut_links_2d_device(cut, lb)
total = total + cut_links_mismatches(lk, lb)
if (nlinks == 0_I4P) ok = .false.
ok = ok .and. total == 0_I4P
print '(A,I0,A,I0)', 'total mismatches: ', total, ', links: ', nlinks
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine grid(n)
   !< n x n cells over [0.003, 0.997]^2 grown by 2 %.
   integer(I4P), intent(in) :: n !< Cells per axis.
   real(R8P)                :: a, b !< Box.
   integer(I4P)             :: q    !< Counter.

   a = 0.0031_R8P - 0.02_R8P * 0.9938_R8P ; b = 0.9969_R8P + 0.02_R8P * 0.9938_R8P
   x = [(a + (b - a) * q / real(n, R8P), q = 0, n)] ; y = x
   endsubroutine grid

   subroutine case(label)
   !< The links of the device path on c against the host's on the grid x, y.
   character(*), intent(in) :: label      !< Case.
   type(cut_cells_2d_t)     :: ch, cr     !< Host cut, resident cut.
   type(cut_links_t)        :: lh, ld, lb !< Host, device links, a copy of resident ones.
   integer(I4P)             :: mm, st, r, q !< Work.

   mm = 0_I4P
   call c%cut_cells(x=x, y=y, cut=ch, status=st)
   call cut_links_2d(ch, lh)
   call c%cut_cells_device(x=x, y=y, cut=cr, status=st, resident=.true.)
   do q = 1, 2 ! host cut, resident cut
      do r = 1, 2 ! links copied back, resident
         if (q == 1) then
            call cut_links_2d_device(ch, ld, resident=r == 2)
         else
            call cut_links_2d_device(cr, ld, resident=r == 2)
         endif
         if (r == 2) then
            if (.not. ld%on_device) mm = mm + 1_I4P
            lb = ld ; call lb%exit_device ! a copy: the original stays mapped
            if (lb%on_device .or. .not. ld%on_device) mm = mm + 1_I4P
            call ld%update_host
         elseif (ld%on_device) then
            mm = mm + 1_I4P
         endif
         mm = mm + cut_links_mismatches(lh, ld)
      enddo
   enddo
   call ld%exit_device ; call cr%exit_device
   if (ld%on_device) mm = mm + 1_I4P
   nlinks = nlinks + size(lh%piece, kind=I4P)
   total = total + mm
   print '(A,I0,A,I0,A,I0,A,I0)', label//': cells ', size(ch%alpha), ', pieces ', ch%np, ', links ', size(lh%piece), &
         ', device path mismatches ', mm
   endsubroutine case

   subroutine naca0012(mh, xx, yy)
   !< NACA 0012 with a closed trailing edge, 2 mh vertices, cosine spacing, counter-clockwise (the curve bench's).
   integer(I4P),           intent(in)  :: mh      !< Points per side.
   real(R8P), allocatable, intent(out) :: xx(:), yy(:) !< Vertices.
   real(R8P)                           :: xi      !< Abscissa.
   integer(I4P)                        :: kk, j   !< Counters.

   allocate(xx(2 * mh), yy(2 * mh))
   j = 0
   do kk = mh, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / mh))
      j = j + 1 ; xx(j) = xi ; yy(j) = thick(xi)
   enddo
   do kk = 1, mh - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / mh))
      j = j + 1 ; xx(j) = xi ; yy(j) = -thick(xi)
   enddo
   endsubroutine naca0012

   pure function thick(xc) result(t)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge.
   real(R8P), intent(in) :: xc !< Chord position.
   real(R8P)             :: t  !< Half thickness.

   t = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc ** 2 + 0.2843_R8P * xc ** 3 - 0.1036_R8P * xc ** 4)
   endfunction thick
endprogram fossil_test_cut_links_2d_twin
