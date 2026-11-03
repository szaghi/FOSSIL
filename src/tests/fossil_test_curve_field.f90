!< FOSSIL test: 2-D signed-distance field from a curve (issue #36 step 8).
program fossil_test_curve_field
!< FOSSIL test: 2-D signed-distance field from a curve (issue #36 step 8).
!<
!< - A 4096-gon of radius 0.5 on uniform grids of 32^2, 64^2, 128^2 cells over [-1, 1]^2: the band
!<   cells equal `curve%distance` bit for bit; every cell has the right sign; the far-field
!<   error against the exact distance falls with h and is below 0.02 at 128^2.
!< - The same on a stretched grid (tanh clustering at the circle): band bit-exact, signs right,
!<   far-field error below 3 times the largest spacing.
!< - One block entirely inside and one entirely outside the curve (a quadtree block that does
!<   not see it): empty band, one point query, the whole block -FAR_VALUE / +FAR_VALUE.
!< - An open circle (a gap of 1/16 of the perimeter) with SIGN_FAST_WINDING_NUMBER: every cell
!<   farther than 0.2 from the gap gets the circle's inside/outside right.
!< - Refusals: a curve on a 3-D field, a surface on a 2-D field.

use fossil, only : curve_object, distance_field_object, surface_stl_object, FAR_VALUE, SIGN_FAST_WINDING_NUMBER, &
                   STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter        :: PI = acos(-1._R8P), R = 0.5_R8P
type(curve_object)          :: c, co
type(distance_field_object) :: f
type(surface_stl_object)    :: cube
real(R8P)                   :: u(4096), v(4096), err(3), x(129), hmax
integer(I4P)                :: k, m, st, fb, i, j, wrong
logical                     :: ok(5)

do k = 1, 4096
   u(k) = R * cos(2._R8P * PI * (k - 1) / 4096) ; v(k) = R * sin(2._R8P * PI * (k - 1) / 4096)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[4096_I4P])

! uniform grids
ok(1) = .true.
do m = 1, 3
   call uniform(32 * 2**(m - 1), x)
   call f%initialize(x=x(1:32 * 2**(m - 1)), y=x(1:32 * 2**(m - 1)))
   call f%compute(curve=c, status=st)
   ok(1) = ok(1) .and. st == STATUS_OK .and. f%get_is_2d() .and. check(f, c, x(1:32 * 2**(m - 1)), err(m))
enddo
ok(1) = ok(1) .and. err(2) < err(1) .and. err(3) < err(2) .and. err(3) < 0.02_R8P
print '(A,L1,3ES11.3)', 'uniform grids 32, 64, 128, far-field max error: ', ok(1), err

! stretched grid
do k = 1, 129
   x(k) = tanh(2._R8P * (-1._R8P + 2._R8P * (k - 1) / 128._R8P)) / tanh(2._R8P)
enddo
hmax = maxval(x(2:129) - x(1:128))
call f%initialize(x=x, y=x)
call f%compute(curve=c, status=st)
ok(2) = st == STATUS_OK .and. check(f, c, x, err(1)) .and. err(1) < 3._R8P * hmax
print '(A,L1,2ES11.3)', 'stretched grid, far-field error, max spacing: ', ok(2), err(1), hmax

! blocks that do not see the curve
call uniform(16, x)
call f%initialize(x=0.2_R8P * x(1:16), y=0.2_R8P * x(1:16))
call f%compute(curve=c, fallback_cells=fb, status=st)
ok(3) = st == STATUS_OK .and. f%get_band_cells() == 0 .and. fb == 1 .and. all(f%phi == -FAR_VALUE)
call f%initialize(x=3._R8P + 0.2_R8P * x(1:16), y=0.2_R8P * x(1:16))
call f%compute(curve=c, fallback_cells=fb, status=st)
ok(3) = ok(3) .and. st == STATUS_OK .and. f%get_band_cells() == 0 .and. fb == 1 .and. all(f%phi == FAR_VALUE)
print '(A,L1)', 'blocks inside / outside, no band: ', ok(3)

! open circle, fast winding number sign
call co%adopt_loops(u=u(1:3840), v=v(1:3840), loop_sizes=[3840_I4P], closed=[.false.])
call uniform(64, x)
call f%initialize(x=x(1:64), y=x(1:64))
call f%compute(curve=co, sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=st)
wrong = 0
do j = 1, 64
   do i = 1, 64
      if (hypot(x(i) - R * cos(-PI / 16._R8P), x(j) - R * sin(-PI / 16._R8P)) < 0.2_R8P) cycle
      if ((f%phi(i, j, 1) < 0._R8P) .neqv. (hypot(x(i), x(j)) < R)) wrong = wrong + 1
   enddo
enddo
ok(4) = st == STATUS_OK .and. wrong == 0
print '(A,I0)', 'open circle (FWN), wrong signs away from the gap: ', wrong

! refusals
call f%initialize(x=x(1:8), y=x(1:8), z=x(1:8))
call f%compute_band(curve=c, status=st)
ok(5) = st == STATUS_INVALID_INPUT
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call f%initialize(x=x(1:8), y=x(1:8))
call f%compute_band(surface=cube, status=st)
ok(5) = ok(5) .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals (curve on 3-D, surface on 2-D): ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine uniform(n, xc)
   !< n cell centres of a uniform grid over [-1, 1].
   integer(I4P), intent(in)  :: n     !< Cells.
   real(R8P),    intent(out) :: xc(:) !< Centres (first n set).
   integer(I4P)              :: q     !< Counter.

   do q = 1, n
      xc(q) = -1._R8P + (2._R8P * q - 1._R8P) / n
   enddo
   endsubroutine uniform

   function check(fld, cur, xc, far_err) result(good)
   !< Band cells bit-exact with cur%distance, every sign right, far-field max |phi - exact|.
   type(distance_field_object), intent(in)  :: fld     !< Field.
   type(curve_object),          intent(in)  :: cur     !< Curve.
   real(R8P),                   intent(in)  :: xc(:)   !< Coordinates (both axes).
   real(R8P),                   intent(out) :: far_err !< Far-field error.
   logical                                  :: good    !< Result.
   real(R8P)                                :: dex     !< Exact distance.
   integer(I4P)                             :: p, q, nb, ns !< Counters.

   far_err = 0._R8P ; nb = 0 ; ns = 0
   do q = 1, size(xc)
      do p = 1, size(xc)
         dex = cur%distance([xc(p), xc(q)])
         if (fld%in_band(p, q, 1)) then
            if (fld%phi(p, q, 1) /= dex) nb = nb + 1
         else
            far_err = max(far_err, abs(fld%phi(p, q, 1) - dex))
         endif
         if ((fld%phi(p, q, 1) < 0._R8P) .neqv. (hypot(xc(p), xc(q)) < R)) ns = ns + 1
      enddo
   enddo
   good = nb == 0 .and. ns == 0 .and. fld%get_band_cells() > 0
   if (.not. good) print '(A,I0,A,I0)', '   band mismatches ', nb, '  wrong signs ', ns
   endfunction check
endprogram fossil_test_curve_field
