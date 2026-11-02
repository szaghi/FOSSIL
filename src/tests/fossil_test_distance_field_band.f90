!< FOSSIL test: narrow band of `distance_field_object` (issue #32, steps 1-2).
program fossil_test_distance_field_band
!< FOSSIL test: narrow band of `distance_field_object` (issue #32, steps 1-2).
!<
!< On cube-02 (unit cube), for a uniform grid, a 1.1-stretched grid with a halo, and an
!< explicit band:
!<  - band membership equals the brute-force set of interior cells with |d| <= band, cell by
!<    cell (so the counts match too);
!<  - band values are bit-exact against a per-cell `surface%distance` with the same options;
!<  - halo and non-band cells hold FAR_VALUE.
!< Invalid grids and a non-positive band return STATUS_INVALID_INPUT.
!< `use_device=.true.` (host loop under gfortran, GPU under nvfortran) gives the host band bit
!< for bit and leaves the surface's residency as it found it; an unsupported device sign
!< algorithm is refused.

use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_SOLID_ANGLE, SIGN_FAST_WINDING_NUMBER, STATUS_OK, &
                  STATUS_INVALID_INPUT
use fossil_distance_field, only : distance_field_object, FAR_VALUE
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)    :: cube   !< Unit cube.
type(distance_field_object) :: field  !< Field under test.
real(R8P), allocatable      :: xu(:)  !< Uniform axis.
real(R8P), allocatable      :: xs(:)  !< Stretched axis.
logical                     :: ok     !< Running verdict.
integer(I4P)                :: st     !< Status.
integer(I4P)                :: i      !< Counter.

ok = .true.
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)

! uniform grid: 24 points spanning [-0.5, 1.5]
allocate(xu(24))
do i = 1, 24
   xu(i) = -0.5_R8P + 2._R8P * real(i - 1, R8P) / 23._R8P
enddo
call check_case(label='uniform, default band', x=xu, ng=0_I4P)
call check_case(label='uniform, band 0.13', x=xu, ng=0_I4P, band=0.13_R8P)

! stretched grid: uniform 0.05 spacing on [-0.1, 1.1], then 1.1 growth for 10 points per side
xs = stretched_axis()
call check_case(label='stretched, default band, ng=2', x=xs, ng=2_I4P)

! invalid inputs
call field%initialize(x=[0._R8P, 1._R8P, 1._R8P], y=xu, z=xu, status=st)
ok = ok .and. st == STATUS_INVALID_INPUT .and. .not. allocated(field%phi)
call field%initialize(x=xu(1:4), y=xu, z=xu, ng=2_I4P, status=st)
ok = ok .and. st == STATUS_INVALID_INPUT
call field%initialize(x=xu, y=xu, z=xu, status=st)
call field%compute_band(surface=cube, band=0._R8P, status=st)
ok = ok .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'invalid inputs refused: ', ok

call device_case

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine device_case
   !< use_device band == host band bit for bit (PN and FWN); refusal of SIGN_SOLID_ANGLE.
   type(distance_field_object) :: host, dev
   integer(I4P)                :: st_h, st_d, st_r
   logical                     :: case_ok

   call host%initialize(x=xs, y=xs, z=xs, ng=2_I4P)
   call dev%initialize(x=xs, y=xs, z=xs, ng=2_I4P)
   call host%compute_band(surface=cube, sign_algorithm=SIGN_PSEUDO_NORMAL, status=st_h)
   call dev%compute_band(surface=cube, sign_algorithm=SIGN_PSEUDO_NORMAL, use_device=.true., status=st_d)
   case_ok = st_h == STATUS_OK .and. st_d == STATUS_OK .and. all(host%phi == dev%phi) .and. &
             all(host%in_band .eqv. dev%in_band) .and. .not. cube%is_on_device
   call host%compute_band(surface=cube, sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=st_h)
   call cube%enter_device   ! already resident: compute_band must leave it resident
   call dev%compute_band(surface=cube, sign_algorithm=SIGN_FAST_WINDING_NUMBER, use_device=.true., status=st_d)
   case_ok = case_ok .and. st_h == STATUS_OK .and. st_d == STATUS_OK .and. all(host%phi == dev%phi) .and. &
             cube%is_on_device
   call dev%compute_band(surface=cube, sign_algorithm=SIGN_SOLID_ANGLE, use_device=.true., status=st_r)
   call cube%exit_device
   case_ok = case_ok .and. st_r == STATUS_INVALID_INPUT .and. dev%get_band_cells() == 0
   print '(A,L1)', 'device band == host band (PN, FWN), residency kept, SIGN_SOLID_ANGLE refused: ', case_ok
   ok = ok .and. case_ok
   endsubroutine device_case

   function stretched_axis() result(c)
   !< Uniform core on [-0.1, 1.1] (spacing 0.05) plus 10 points per side growing by 1.1.
   real(R8P), allocatable :: c(:)       !< Axis.
   real(R8P)              :: core(25)   !< Core points.
   real(R8P)              :: left(10)   !< Left extension.
   real(R8P)              :: right(10)  !< Right extension.
   real(R8P)              :: h          !< Spacing.
   integer(I4P)           :: m          !< Counter.

   do m = 1, 25
      core(m) = -0.1_R8P + 0.05_R8P * real(m - 1, R8P)
   enddo
   h = 0.05_R8P
   right(1) = core(25) + h * 1.1_R8P
   left(10) = core(1)  - h * 1.1_R8P
   do m = 2, 10
      h = h * 1.1_R8P
      right(m)      = right(m - 1) + h * 1.1_R8P
      left(11 - m)  = left(12 - m) - h * 1.1_R8P
   enddo
   c = [left, core, right]
   endfunction stretched_axis

   subroutine check_case(label, x, ng, band)
   !< Build the band on an x = y = z grid and compare with brute force.
   character(*), intent(in)           :: label  !< Case label.
   real(R8P),    intent(in)           :: x(:)   !< Axis (used for x, y, z).
   integer(I4P), intent(in)           :: ng     !< Halo width.
   real(R8P),    intent(in), optional :: band   !< Explicit band.
   real(R8P)                          :: b      !< Band used.
   real(R8P)                          :: d      !< Brute-force distance.
   type(vector_R8P)                   :: p      !< Point.
   integer(I4P)                       :: n, a, c, e, n_brute, n_mismatch, n_value, n_far
   logical                            :: interior, case_ok

   n = size(x, kind=I4P)
   call field%initialize(x=x, y=x, z=x, ng=ng, status=st)
   case_ok = st == STATUS_OK
   call field%compute_band(surface=cube, band=band, is_signed=.true., sign_algorithm=SIGN_PSEUDO_NORMAL, status=st)
   case_ok = case_ok .and. st == STATUS_OK
   b = field%get_band_width()
   n_brute = 0 ; n_mismatch = 0 ; n_value = 0 ; n_far = 0
   do e = 1, n
      do c = 1, n
         do a = 1, n
            interior = min(a, c, e) > ng .and. max(a, c, e) <= n - ng
            p%x = x(a) ; p%y = x(c) ; p%z = x(e)
            if (interior) then
               d = cube%distance(point=p, is_signed=.true., sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
               if (abs(d) <= b) n_brute = n_brute + 1
               if ((abs(d) <= b) .neqv. field%in_band(a, c, e)) n_mismatch = n_mismatch + 1
               if (field%in_band(a, c, e)) then
                  if (field%phi(a, c, e) /= d) n_value = n_value + 1
               elseif (field%phi(a, c, e) /= FAR_VALUE) then
                  n_far = n_far + 1
               endif
            else
               if (field%in_band(a, c, e) .or. field%phi(a, c, e) /= FAR_VALUE) n_far = n_far + 1
            endif
         enddo
      enddo
   enddo
   case_ok = case_ok .and. n_brute > 0 .and. n_brute == field%get_band_cells() .and. &
             n_mismatch == 0 .and. n_value == 0 .and. n_far == 0
   print '(A,A,ES10.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,L1)', label, ': band ', b, '  cells ', field%get_band_cells(), &
         '  brute ', n_brute, '  membership mismatches ', n_mismatch, '  value mismatches ', n_value, &
         '  bad far/halo ', n_far, '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine check_case
endprogram fossil_test_distance_field_band
