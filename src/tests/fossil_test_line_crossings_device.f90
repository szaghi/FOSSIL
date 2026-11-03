!< FOSSIL test: batched and device grid-line crossings (issue #40 step 2).
program fossil_test_line_crossings_device
!< FOSSIL test: batched and device grid-line crossings (issue #40 step 2).
!<
!< - The 9 sign-benchmark bodies, 17 x 17 lines along each axis over their box (through edges and
!<   vertices where the grid hits them): `line_crossings_many` equals `line_crossings` line by line,
!<   and `line_crossings_many_device` (after `enter_device`; on the host in tests-gnu) equals
!<   `line_crossings_many` bit for bit, buffers included.
!< - The same for curves: the square annulus, a NACA 0012, a 64-gon, along u and v.
!< - A buffer of 1 on a line with 2 crossings: overflow reported, one crossing kept.
!< - Use before `enter_device` is refused.

use fossil,                   only : surface_stl_object, curve_object, STATUS_INVALID_INPUT
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
real(R8P)                 :: u(200), v(200), t1(1, 1)
integer(I4P)              :: i, k, nbad3, nbad2, st, f1(1, 1), c1(1)
logical                   :: ok(3), e1(1, 1), o1(1)

call bench_bodies(b)
nbad3 = 0
do i = 1, size(b)
   call to_surface(b(i)%soup, s)
   nbad3 = nbad3 + check3(s)
enddo
nbad2 = 0
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
nbad2 = nbad2 + check2(c)
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
nbad2 = nbad2 + check2(c)
do k = 1, 64
   u(k) = cos(2 * PI * (k - 1) / 64) ; v(k) = sin(2 * PI * (k - 1) / 64)
enddo
call c%adopt_loops(u=u(1:64), v=v(1:64), loop_sizes=[64_I4P])
nbad2 = nbad2 + check2(c)
ok(1) = nbad3 == 0 .and. nbad2 == 0
print '(A,I0,A,I0)', 'many = single and device = many: 3-D mismatches ', nbad3, '  2-D mismatches ', nbad2

! overflow: the cube, a line with 2 crossings, buffer 1
call to_surface(b(1)%soup, s)
call s%line_crossings_many(axis=1_I4P, a=[0.3_R8P], b=[0.6_R8P], t0=-1._R8P, t1=2._R8P, t=t1, facets=f1, entering=e1, &
                           counts=c1, overflow=o1)
ok(2) = o1(1) .and. c1(1) == 1 .and. (t1(1, 1) == 0._R8P .or. t1(1, 1) == 1._R8P)
print '(A,L1)', 'buffer of 1 on 2 crossings: overflow reported, one crossing kept: ', ok(2)

call s%line_crossings_many_device(axis=1_I4P, a=[0.3_R8P], b=[0.6_R8P], t0=-1._R8P, t1=2._R8P, t=t1, facets=f1, &
                                  entering=e1, counts=c1, overflow=o1, status=st)
ok(3) = st == STATUS_INVALID_INPUT
call c%line_crossings_many_device(axis=1_I4P, a=[0._R8P], t0=-1._R8P, t1=2._R8P, t=t1, segments=f1, entering=e1, &
                                  counts=c1, overflow=o1, status=st)
ok(3) = ok(3) .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'device use before enter_device refused: ', ok(3)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check3(sf) result(nbad)
   !< many vs single and device vs many on 17 x 17 lines per axis.
   type(surface_stl_object), intent(inout) :: sf   !< Surface.
   integer                                 :: nbad !< Mismatches.
   integer, parameter        :: NL = 17 * 17, NMAX = 16 !< Lines, buffer.
   real(R8P)                 :: a(NL), bb(NL), lo(3), hi(3), tc(NMAX, NL), td(NMAX, NL) !< Work.
   integer(I4P)              :: fc(NMAX, NL), fd(NMAX, NL), cc(NL), cd(NL), ax, i1, i2, j, k, m, n !< Work.
   logical                   :: ec(NMAX, NL), ed(NMAX, NL), oc(NL), od(NL) !< Work.
   real(R8P),    allocatable :: ts(:)  !< Single-line crossings.
   integer(I4P), allocatable :: fs(:)  !< Work.
   logical,      allocatable :: es(:)  !< Work.
   type(vector_R8P)          :: p      !< Work.

   p = sf%get_bmin() ; lo = [p%x, p%y, p%z] ; p = sf%get_bmax() ; hi = [p%x, p%y, p%z]
   nbad = 0
   call sf%enter_device
   do ax = 1, 3
      i1 = mod(ax, 3) + 1 ; i2 = mod(ax + 1, 3) + 1
      m = 0
      do k = 0, 16
         do j = 0, 16
            m = m + 1
            a(m) = lo(i1) + (hi(i1) - lo(i1)) * j / 16._R8P ; bb(m) = lo(i2) + (hi(i2) - lo(i2)) * k / 16._R8P
         enddo
      enddo
      call sf%line_crossings_many(axis=int(ax, I4P), a=a, b=bb, t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=tc, facets=fc, &
                                  entering=ec, counts=cc, overflow=oc)
      call sf%line_crossings_many_device(axis=int(ax, I4P), a=a, b=bb, t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=td, &
                                         facets=fd, entering=ed, counts=cd, overflow=od)
      nbad = nbad + count(tc /= td) + count(fc /= fd) + count(ec .neqv. ed) + count(cc /= cd) + count(oc .neqv. od)
      do m = 1, NL
         call sf%line_crossings(axis=int(ax, I4P), a=a(m), b=bb(m), t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=ts, &
                                facets=fs, entering=es, n=n)
         if (oc(m) .or. n /= cc(m)) then
            nbad = nbad + 1
         elseif (n > 0) then
            if (any(ts /= tc(1:n, m)) .or. any(fs /= fc(1:n, m)) .or. any(es .neqv. ec(1:n, m))) nbad = nbad + 1
         endif
      enddo
   enddo
   call sf%exit_device
   endfunction check3

   function check2(cur) result(nbad)
   !< many vs single and device vs many on 33 lines per axis.
   type(curve_object), intent(inout) :: cur  !< Curve.
   integer                           :: nbad !< Mismatches.
   integer, parameter        :: NL = 33, NMAX = 16 !< Lines, buffer.
   real(R8P)                 :: a(NL), lo(2), hi(2), tc(NMAX, NL), td(NMAX, NL) !< Work.
   integer(I4P)              :: fc(NMAX, NL), fd(NMAX, NL), cc(NL), cd(NL), ax, m, n !< Work.
   logical                   :: ec(NMAX, NL), ed(NMAX, NL), oc(NL), od(NL) !< Work.
   real(R8P),    allocatable :: ts(:) !< Work.
   integer(I4P), allocatable :: fs(:) !< Work.
   logical,      allocatable :: es(:) !< Work.

   lo = cur%get_bmin() ; hi = cur%get_bmax()
   nbad = 0
   call cur%enter_device
   do ax = 1, 2
      do m = 1, NL
         a(m) = lo(3 - ax) + (hi(3 - ax) - lo(3 - ax)) * (m - 1) / real(NL - 1, R8P)
      enddo
      call cur%line_crossings_many(axis=int(ax, I4P), a=a, t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=tc, segments=fc, &
                                   entering=ec, counts=cc, overflow=oc)
      call cur%line_crossings_many_device(axis=int(ax, I4P), a=a, t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=td, &
                                          segments=fd, entering=ed, counts=cd, overflow=od)
      nbad = nbad + count(tc /= td) + count(fc /= fd) + count(ec .neqv. ed) + count(cc /= cd) + count(oc .neqv. od)
      do m = 1, NL
         call cur%line_crossings(axis=int(ax, I4P), a=a(m), t0=lo(ax) - 1._R8P, t1=hi(ax) + 1._R8P, t=ts, segments=fs, &
                                 entering=es, n=n)
         if (oc(m) .or. n /= cc(m)) then
            nbad = nbad + 1
         elseif (n > 0) then
            if (any(ts /= tc(1:n, m)) .or. any(fs /= fc(1:n, m)) .or. any(es .neqv. ec(1:n, m))) nbad = nbad + 1
         endif
      enddo
   enddo
   call cur%exit_device
   endfunction check2

   pure function thickness(x) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: x !< Chord position.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.126_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_line_crossings_device
