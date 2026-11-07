!< FOSSIL validation: time of first contact at size (issue #73, step 7).
!<
!< 1. The bunny (decimated to 5 000 facets) against the dragon (2 000), both scaled to a unit diagonal, over 8 fixed-seed
!<    steps that turn and translate the dragon across the bunny: time_of_impact with the plain and the directional
!<    bound against a brute force (200 samples of d, then the first time d <= tol / 2 bisected): no sample of [0, toi) is
!<    below tol / 2, d(toi) <= tol, toi <= the bisected time; iterations, pair evaluations and times reported.
!< 2. Sliding passes: a cube sliding past another, face parallel to face, at a gap tol 2^k, k = -2..8: contact at and
!<    below tol; above it the plain bound advances by about (gap - tol / 2) / mu per step and hits the cap
!<    (CCD_NOT_CONVERGED, a certified lower bound), the directional bound (n . dt = 0 for the parallel faces) separates
!<    in a few steps (both reported, the statuses asserted).
!< 3. Several bodies: 6 spheres (10^3 marching cubes) with fixed-seed velocities over 10 sub-steps (passing through
!<    each other); per sub-step first_impact equals the loop over all pairs, bit for bit (toi, pair, status); pruning reported.
!< Exit status 0 if every check holds.
program fossil_val_ccd

use fossil,                only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, &
                                  first_impact, CCD_CONTACT, CCD_SEPARATE, CCD_NOT_CONVERGED
use fossil_boolean_bodies, only : BODY_SPHERE, mc_body
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

real(R8P), parameter     :: TOL = 1e-6_R8P
type(surface_stl_object) :: a, b, s(6)
type(rigid_pose_t)       :: a0, a1, b0, p0(6), p1(6)
real(R8P)                :: t0, t1, tp, td, tbrute, toi, tb, d, v(3, 6)
integer(I4P)             :: k, i, st, sb, ip, id, np, pair(2), pb(2), nbad, ntest, nall
integer(I8P)             :: seed
logical                  :: all_ok

all_ok = .true.
b0 = pose_identity()

! 1. bunny vs dragon
call a%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.) ; call a%decimate(target_facets=5000_I4P)
call b%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.) ; call b%decimate(target_facets=2000_I4P)
call scale_to_unit(a) ; call scale_to_unit(b)
seed = 20261107_I8P ; tp = 0._R8P ; td = 0._R8P ; tbrute = 0._R8P
print '(A)', 'step  status  toi           plain steps  directional steps  pairs   bisected - toi'
do k = 1, 8
   a0 = pose_from_axis_angle([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P], 6.28_R8P * rnd(), &
                             [-1.2_R8P, 0.4_R8P * (rnd() - 0.5_R8P), 0.4_R8P * (rnd() - 0.5_R8P)])
   a1 = pose_from_axis_angle([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P], 6.28_R8P * rnd(), &
                             [1.2_R8P, 0.4_R8P * (rnd() - 0.5_R8P), 0.4_R8P * (rnd() - 0.5_R8P)])
   call cpu_time(t0)
   call b%time_of_impact(a, a0, a1, b0, b0, toi, st, tolerance=TOL, distance=d, iterations=ip)
   call cpu_time(t1) ; tp = tp + (t1 - t0)
   call b%time_of_impact(a, a0, a1, b0, b0, tb, sb, tolerance=TOL, iterations=id, directional=.true., pairs=np)
   call cpu_time(t0) ; td = td + (t0 - t1)
   call brute(toi, st, d, ip, id, np, k)
   call cpu_time(t1) ; tbrute = tbrute + (t1 - t0)
   if (sb /= st .or. id > ip) then
      all_ok = .false. ; print '(A)', '   directional status differs or more steps'
   endif
enddo
print '(A,3F9.3)', 'times: plain, directional, brute force [s] ', tp, td, tbrute

! 2. grazing passes
call cube(a, [0._R8P, 0._R8P, 0._R8P])
print '(A)', 'gap / tol   plain status  steps   directional status  steps'
do k = -2, 8
   call cube(b, [2._R8P, 1._R8P + TOL * 2._R8P ** k, 0._R8P])
   a0 = pose_identity() ; a1 = a0 ; a1%t = [4._R8P, 0._R8P, 0._R8P]
   call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=TOL, iterations=ip)
   call a%time_of_impact(b, a0, a1, b0, b0, tb, sb, tolerance=TOL, iterations=id, directional=.true.)
   print '(F9.2,I15,I7,I21,I7)', 2._R8P ** k, st, ip, sb, id
   if (k <= 0 .and. (st /= CCD_CONTACT .or. sb /= CCD_CONTACT)) all_ok = .false.
   if (k > 0 .and. (st /= CCD_NOT_CONVERGED .or. sb /= CCD_SEPARATE)) all_ok = .false.
enddo

! 3. several bodies
do i = 1, 6
   call mc_body(BODY_SPHERE, 10_I4P, s(i))
   p0(i) = pose_identity() ; p0(i)%t = [1.3_R8P * mod(i - 1, 3), 1.3_R8P * ((i - 1) / 3), 0._R8P]
   v(:, i) = 0.25_R8P * ([1.3_R8P, 0.65_R8P, 0._R8P] - p0(i)%t) + 0.2_R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, 0._R8P]
enddo
nbad = 0 ; nall = 0
do k = 1, 10
   do i = 1, 6
      p1(i) = p0(i) ; p1(i)%t = p0(i)%t + v(:, i)
   enddo
   call first_impact(s, p0, p1, toi, pair, st, tolerance=TOL, tested=ntest)
   call brute_bodies(tb, pb, sb)
   if (toi /= tb .or. any(pair /= pb) .or. st /= sb) nbad = nbad + 1
   nall = nall + ntest
   print '(A,I3,A,I2,A,F9.6,A,2I2,A,I0)', 'sub-step ', k, ': status ', st, ', toi ', toi, ', pair ', pair, &
         ', pairs tested ', ntest
   p0 = p1 ! the bodies pass through each other: every sub-step is a new configuration
enddo
print '(A,I0,A,I0,A)', 'several bodies: mismatches with the loop over all pairs ', nbad, ', pairs tested ', nall, ' of 150'
if (nbad /= 0) all_ok = .false.

print '(A,L1)', 'Are all tests passed? ', all_ok
if (.not. all_ok) error stop 1

contains
   subroutine brute(toi_, st_, d_, ip_, id_, np_, k_)
   !< 200 samples of d over the step, the first time d <= tol / 2 bisected; check toi against them.
   real(R8P),    intent(in) :: toi_, d_       !< Result, its distance.
   integer(I4P), intent(in) :: st_, ip_, id_, np_, k_ !< Status, steps, pairs, step number.
   real(R8P)                :: t, dt, lo, hi, mid, dm, first
   integer(I4P)             :: jj, n, nlow
   logical                  :: ok

   nlow = 0 ; first = 2._R8P
   do jj = 0, 200
      t = real(jj, R8P) / 200
      call b%distance_to(a, dt, pose_self=pose_interpolate(a0, a1, t))
      if (t < toi_ .and. dt < 0.5_R8P * TOL * (1._R8P - 1e-9_R8P)) nlow = nlow + 1
      if (dt <= 0.5_R8P * TOL .and. first > 1._R8P) then
         lo = max(0._R8P, t - 1._R8P / 200) ; hi = t
         do n = 1, 60
            mid = 0.5_R8P * (lo + hi)
            call b%distance_to(a, dm, pose_self=pose_interpolate(a0, a1, mid))
            if (dm <= 0.5_R8P * TOL) then
               hi = mid
            else
               lo = mid
            endif
         enddo
         first = hi
      endif
   enddo
   ok = nlow == 0
   if (st_ == CCD_CONTACT) ok = ok .and. d_ <= TOL .and. toi_ <= first
   if (st_ == CCD_SEPARATE) ok = ok .and. first > 1._R8P
   print '(I4,I8,ES14.6,I13,I19,I7,ES17.6,A,L1)', k_, st_, toi_, ip_, id_, np_, first - toi_, '  ok ', ok
   all_ok = all_ok .and. ok
   endsubroutine brute

   subroutine brute_bodies(t_, p_, s_)
   !< The loop over all pairs, ties to the smallest pair.
   real(R8P),    intent(out) :: t_    !< Earliest toi.
   integer(I4P), intent(out) :: p_(2) !< Pair.
   integer(I4P), intent(out) :: s_    !< Status.
   integer(I4P)              :: ii, jj, stt
   real(R8P)                 :: t

   t_ = 1._R8P ; p_ = 0 ; s_ = CCD_SEPARATE
   do ii = 1, 6
      do jj = ii + 1, 6
         call s(ii)%time_of_impact(s(jj), p0(ii), p1(ii), p0(jj), p1(jj), t, stt, tolerance=TOL)
         if (stt == CCD_SEPARATE) cycle
         if (t < t_ .or. (t == t_ .and. s_ == CCD_SEPARATE)) then
            t_ = t ; p_ = [ii, jj] ; s_ = stt
         endif
      enddo
   enddo
   endsubroutine brute_bodies

   subroutine cube(sf, o)
   !< The unit cube [0, 1]^3 moved by o.
   type(surface_stl_object), intent(inout) :: sf   !< Cube.
   real(R8P),                intent(in)    :: o(3) !< Offset.
   type(surface_stl_object)                :: fresh

   sf = fresh
   call sf%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call sf%translate(x=o(1), y=o(2), z=o(3)) ; call sf%analyze
   endsubroutine cube

   subroutine scale_to_unit(sf)
   !< Centre at the origin, unit bounding-box diagonal.
   type(surface_stl_object), intent(inout) :: sf !< Body.
   type(vector_R8P)                        :: lo, hi
   real(R8P)                               :: f

   lo = sf%get_bmin() ; hi = sf%get_bmax()
   call sf%translate(delta=-0.5_R8P * (lo + hi))
   f = 1._R8P / norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call sf%resize(factor=vector_R8P(f, f, f)) ; call sf%analyze
   endsubroutine scale_to_unit

   function rnd() result(r)
   !< Park-Miller, fixed seed.
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_val_ccd
