!< VALIDATION version (issue #45): the full, slow test, run by scripts/run_validation.sh at -O2.
!< The fast regression variant run per push is src/tests/fossil_test_distance_field_fsm.f90.
!<
!< FOSSIL test: far field of `distance_field_object` by fast sweeping (issue #32, step 4).
program fossil_val_distance_field_fsm
!< FOSSIL test: far field of `distance_field_object` by fast sweeping (issue #32, step 4).
!<
!< Against the analytic signed distance of the unit cube [0, 1]^3 (cube-02):
!<  - uniform grids on [-1, 2]^3 with 25, 49, 97 points per axis: the far-field max error
!<    falls like h |log h|, the FSM bound near point-like features (the cube's corners and
!<    edges; Zhao 2005): observed order against h |log h| in [0.8, 1.2] for both pairs. The
!<    order against h alone is printed too (about 0.6-0.7 for that reason), and the mean
!<    error's order against h;
!<  - the band cells keep their exact values bit for bit through the sweeps;
!<  - the sweeps converge within 5 cycles (40 directional sweeps);
!<  - the field is symmetric on the symmetric grid (to 1e-12);
!<  - a 1.15-stretched grid: far-field max error below the largest spacing;
!<  - two blocks with a one-cell halo, swept with caller-driven halo exchange, reproduce the
!<    single-block field (to 1e-10);
!<  - OpenMP: the sweep with all threads equals the sweep with one thread bit for bit (97^3;
!<    trivially true without OpenMP), timings printed.

use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, STATUS_OK
use fossil_distance_field, only : distance_field_object, FAR_VALUE
use penf, only : I4P, I8P, R8P

implicit none

type(surface_stl_object)    :: cube                !< Unit cube.
real(R8P)                   :: err(3), h(3), order !< Errors, spacings, observed order.
real(R8P)                   :: mean(3), order_c   !< Mean errors, coarse-pair order.
integer(I4P), parameter     :: NS(3) = [25, 49, 97]
logical                     :: ok
integer(I4P)                :: r

ok = .true.
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)

do r = 1, 3
   call uniform_case(n=NS(r), err=err(r), mean=mean(r), h=h(r))
enddo
order_c = log(err(1) / err(2)) / log(hlogh(h(1)) / hlogh(h(2)))
order   = log(err(2) / err(3)) / log(hlogh(h(2)) / hlogh(h(3)))
print '(A,F6.3,A,F6.3)', 'max error, observed order vs h |log h|: coarse pair ', order_c, '  fine pair ', order
print '(A,F6.3,A,F6.3)', 'max error, observed order vs h:         coarse pair ', log(err(1) / err(2)) / log(h(1) / h(2)), &
                         '  fine pair ', log(err(2) / err(3)) / log(h(2) / h(3))
print '(A,F6.3,A,F6.3)', 'mean error, observed order vs h:        coarse pair ', log(mean(1) / mean(2)) / log(h(1) / h(2)), &
                         '  fine pair ', log(mean(2) / mean(3)) / log(h(2) / h(3))
ok = ok .and. err(3) < err(2) .and. err(2) < err(1) .and. &
     order_c >= 0.8_R8P .and. order_c <= 1.2_R8P .and. order >= 0.8_R8P .and. order <= 1.2_R8P

call stretched_case
call two_block_case
call threads_case

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   pure function hlogh(hh) result(v)
   !< h |log h|.
   real(R8P), intent(in) :: hh
   real(R8P)             :: v

   v = hh * abs(log(hh))
   endfunction hlogh

   pure elemental function box_sdf(px, py, pz) result(d)
   !< Signed distance to the unit cube [0, 1]^3 (negative inside).
   real(R8P), intent(in) :: px, py, pz !< Point.
   real(R8P)             :: d          !< Distance.
   real(R8P)             :: qx, qy, qz !< |p - centre| - half side.

   qx = abs(px - 0.5_R8P) - 0.5_R8P
   qy = abs(py - 0.5_R8P) - 0.5_R8P
   qz = abs(pz - 0.5_R8P) - 0.5_R8P
   d = sqrt(max(qx, 0._R8P)**2 + max(qy, 0._R8P)**2 + max(qz, 0._R8P)**2) + min(max(qx, qy, qz), 0._R8P)
   endfunction box_sdf

   subroutine far_error(field, x, y, z, err, mean)
   !< Max and mean |phi - exact| over the non-band cells.
   type(distance_field_object), intent(in)            :: field
   real(R8P),                   intent(in)            :: x(:), y(:), z(:)
   real(R8P),                   intent(out)           :: err
   real(R8P),                   intent(out), optional :: mean
   integer(I4P)                                       :: i, j, k, n
   real(R8P)                                          :: total, e

   err = 0._R8P ; total = 0._R8P ; n = 0
   do k = 1, size(z)
      do j = 1, size(y)
         do i = 1, size(x)
            if (field%in_band(i, j, k)) cycle
            e = abs(field%phi(i, j, k) - box_sdf(x(i), y(j), z(k)))
            err = max(err, e) ; total = total + e ; n = n + 1
         enddo
      enddo
   enddo
   if (present(mean)) mean = total / real(max(n, 1), R8P)
   endsubroutine far_error

   subroutine uniform_case(n, err, mean, h)
   !< Uniform grid on [-1, 2]^3 with n points per axis.
   integer(I4P), intent(in)    :: n
   real(R8P),    intent(out)   :: err, mean, h
   type(distance_field_object) :: field
   real(R8P),    allocatable   :: x(:), band_phi(:,:,:)
   real(R8P)                   :: sym
   integer(I4P)                :: i, sweeps, fb, st
   logical                     :: case_ok

   allocate(x(n))
   do i = 1, n
      x(i) = -1._R8P + 3._R8P * real(i - 1, R8P) / real(n - 1, R8P)
   enddo
   h = x(2) - x(1)
   call field%initialize(x=x, y=x, z=x)
   call field%compute_band(surface=cube, sign_algorithm=SIGN_PSEUDO_NORMAL, status=st)
   call field%compute_signs(surface=cube, fallback_cells=fb)
   band_phi = field%phi
   call field%sweep(sweeps=sweeps)
   call far_error(field, x, x, x, err, mean)
   sym = max(maxval(abs(field%phi - field%phi(n:1:-1, :, :))), maxval(abs(field%phi - field%phi(:, n:1:-1, :))), &
             maxval(abs(field%phi - field%phi(:, :, n:1:-1))))
   case_ok = st == STATUS_OK .and. fb == 0 .and. sweeps <= 40 .and. all(field%phi < FAR_VALUE) .and. &
             all(merge(field%phi == band_phi, .true., field%in_band)) .and. sym <= 1.e-12_R8P
   print '(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A,ES10.3,A,L1)', 'uniform ', n, '^3: h ', h, '  far max err ', err, &
         '  mean ', mean, '  sweeps ', sweeps, '  asymmetry ', sym, '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine uniform_case

   function stretched_axis(h0, growth, n_ext) result(c)
   !< Uniform core on [-0.2, 1.2] with spacing h0, plus n_ext points per side growing by `growth`.
   real(R8P),    intent(in) :: h0, growth
   integer(I4P), intent(in) :: n_ext
   real(R8P), allocatable   :: c(:), core(:), left(:), right(:)
   real(R8P)                :: hh
   integer(I4P)             :: m, nc

   nc = nint(1.4_R8P / h0) + 1
   allocate(core(nc), left(n_ext), right(n_ext))
   do m = 1, nc
      core(m) = -0.2_R8P + h0 * real(m - 1, R8P)
   enddo
   hh = h0
   do m = 1, n_ext
      hh = hh * growth
      if (m == 1) then
         right(1) = core(nc) + hh ; left(n_ext) = core(1) - hh
      else
         right(m) = right(m - 1) + hh ; left(n_ext + 1 - m) = left(n_ext + 2 - m) - hh
      endif
   enddo
   c = [left, core, right]
   endfunction stretched_axis

   subroutine stretched_case
   !< 1.15-stretched grid: far-field error below the largest spacing.
   type(distance_field_object) :: field
   real(R8P),    allocatable   :: x(:)
   real(R8P)                   :: e, hmax
   integer(I4P)                :: sweeps, st
   logical                     :: case_ok

   x = stretched_axis(h0=0.05_R8P, growth=1.15_R8P, n_ext=14)
   hmax = maxval(x(2:) - x(:size(x) - 1))
   call field%initialize(x=x, y=x, z=x)
   call field%compute(surface=cube, sweeps=sweeps, status=st)
   call far_error(field, x, x, x, e)
   case_ok = st == STATUS_OK .and. e < hmax .and. all(field%phi < FAR_VALUE)
   print '(A,I0,A,F8.4,A,ES10.3,A,ES10.3,A,I0,A,L1)', 'stretched ', size(x), '^3 (to ', x(size(x)), '): h max ', hmax, &
         '  far max err ', e, '  sweeps ', sweeps, '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine stretched_case

   subroutine two_block_case
   !< Two blocks split along x with ng = 1, swept with halo exchange, vs one block.
   type(distance_field_object) :: one, a, b
   real(R8P),    allocatable   :: g(:), xa(:), xb(:), yz(:)
   real(R8P)                   :: dg, ca, cb, diff
   integer(I4P), parameter     :: N = 41, M = 17
   integer(I4P)                :: i, it
   logical                     :: case_ok

   allocate(g(N))
   do i = 1, N
      g(i) = -1._R8P + 3._R8P * real(i - 1, R8P) / real(N - 1, R8P)
   enddo
   dg = g(2) - g(1)
   call one%initialize(x=g, y=g, z=g)
   call one%compute(surface=cube, band=0.2_R8P)

   yz = [g(1) - dg, g, g(N) + dg]
   xa = [g(1) - dg, g(1:M), g(M + 1)]          ! interior g(1:M), halo high = g(M+1)
   xb = [g(M), g(M + 1:N), g(N) + dg]          ! halo low = g(M), interior g(M+1:N)
   call a%initialize(x=xa, y=yz, z=yz, ng=1_I4P)
   call b%initialize(x=xb, y=yz, z=yz, ng=1_I4P)
   call a%compute(surface=cube, band=0.2_R8P)
   call b%compute(surface=cube, band=0.2_R8P)
   do it = 1, 50
      a%phi(M + 2, :, :) = b%phi(2, :, :)
      b%phi(1, :, :)     = a%phi(M + 1, :, :)
      call a%sweep(max_sweeps=8_I4P, max_change=ca)
      call b%sweep(max_sweeps=8_I4P, max_change=cb)
      if (max(ca, cb) == 0._R8P) exit
   enddo
   diff = max(maxval(abs(a%phi(2:M + 1, 2:N + 1, 2:N + 1) - one%phi(1:M, :, :))), &
              maxval(abs(b%phi(2:N - M + 1, 2:N + 1, 2:N + 1) - one%phi(M + 1:N, :, :))))
   case_ok = diff <= 1.e-10_R8P
   print '(A,I0,A,ES10.3,A,L1)', 'two blocks, halo exchange: iterations ', it, '  max diff vs one block ', diff, &
         '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine two_block_case

   subroutine threads_case
   !< Sweep from the same signed band with 1 thread and with all threads: bit-exact.
   type(distance_field_object) :: f1, fn
   real(R8P),    allocatable   :: x(:)
   integer(I4P), parameter     :: N = 97
   integer(I4P)                :: i, nthreads
   integer(I8P)                :: c0, c1, rate
   real(R8P)                   :: t1, tn
   logical                     :: case_ok
   !$ integer, external        :: omp_get_max_threads

   allocate(x(N))
   do i = 1, N
      x(i) = -1._R8P + 3._R8P * real(i - 1, R8P) / real(N - 1, R8P)
   enddo
   call f1%initialize(x=x, y=x, z=x)
   call f1%compute_band(surface=cube)
   call f1%compute_signs(surface=cube)
   fn = f1
   nthreads = 1
   !$ nthreads = omp_get_max_threads()
   !$ call omp_set_num_threads(1)
   call system_clock(c0, rate)
   call f1%sweep
   call system_clock(c1)
   t1 = real(c1 - c0, R8P) / real(rate, R8P)
   !$ call omp_set_num_threads(nthreads)
   call system_clock(c0, rate)
   call fn%sweep
   call system_clock(c1)
   tn = real(c1 - c0, R8P) / real(rate, R8P)
   case_ok = all(f1%phi == fn%phi)
   print '(A,I0,A,F8.4,A,I0,A,F8.4,A,L1)', 'threads: 1 thread ', N, '^3 sweep ', t1, ' s,  ', nthreads, ' threads ', tn, &
         ' s,  bit-exact -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine threads_case
endprogram fossil_val_distance_field_fsm
