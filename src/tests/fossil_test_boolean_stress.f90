!< FOSSIL test: boolean stress on snapped box pairs with coplanar faces (issue #33, step 6).
program fossil_test_boolean_stress
!< FOSSIL test: boolean stress on snapped box pairs with coplanar faces (issue #33, step 6).
!<
!< 100 configurations from a fixed-seed generator: A = unit cube, B = a box whose sizes (0.5, 1,
!< 1.5) and offsets (-0.5 to 1, step 0.25) are snapped to a 0.25 grid, so coplanar faces, shared
!< edges and touching contacts are frequent; every other configuration rotates both bodies by
!< the same random proper rotation. The four ops run on each (400 booleans) against the analytic
!< box volumes.
!<
!< Pass criterion (the acceptance of issue #33): no silent wrong volume. A result with status
!< OK must have the analytic volume (1e-10 relative to max(V, 1)); a non-OK status is allowed and
!< counted. Watertightness and manifoldness are reported, and required for the union,
!< intersection and difference of configurations without contact (every overlap length either
!< positive or the boxes apart by a positive gap on some axis): with contact the result can be
!< non-manifold by definition (boxes sharing only an edge), and so can a symmetric difference.

use fossil, only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none

integer(I4P), parameter :: N_CASES = 100
integer(I4P), parameter :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
integer(I8P)            :: seed
integer(I4P)            :: c, k, st, n_ok_status, n_bad_status, n_wrong, n_topo_req, n_topo_bad, n_wt, n_mf
integer(I4P)            :: n_overlap, n_contact, n_apart, n_coplanar !< Configuration mix.
real(R8P)               :: s(3), t(3), lo(3), hi(3), ov(3), vi, va, vb, expected(4), v, axis(3), angle
logical                 :: rotated, contact, wt, mf

seed = 20261102_I8P
n_ok_status = 0 ; n_bad_status = 0 ; n_wrong = 0 ; n_topo_req = 0 ; n_topo_bad = 0 ; n_wt = 0 ; n_mf = 0
n_overlap = 0 ; n_contact = 0 ; n_apart = 0 ; n_coplanar = 0
do c = 1, N_CASES
   do k = 1, 3
      s(k) = 0.5_R8P * real(1 + int(3 * rnd()), R8P)          ! 0.5, 1, 1.5
      t(k) = -0.5_R8P + 0.25_R8P * real(int(7 * rnd()), R8P) ! -0.5 .. 1
   enddo
   rotated = mod(c, 2) == 0
   axis  = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
   angle = 6.283185307179586_R8P * rnd()
   lo = max(0._R8P, t) ; hi = min(1._R8P, t + s)
   ov = hi - lo
   contact = any(ov == 0._R8P) .and. .not. any(ov < 0._R8P)
   vi = product(max(ov, 0._R8P))
   if (any(ov < 0._R8P)) then
      n_apart = n_apart + 1
   elseif (contact) then
      n_contact = n_contact + 1
   else
      n_overlap = n_overlap + 1
   endif
   ! a face plane of B on a face plane of A
   if (any(t == 0._R8P) .or. any(t == 1._R8P) .or. any(t + s == 0._R8P) .or. any(t + s == 1._R8P)) &
      n_coplanar = n_coplanar + 1
   va = 1._R8P ; vb = product(s)
   expected = [va + vb - vi, vi, va - vi, va + vb - 2._R8P * vi]
   do k = 1, 4
      call one(k)
   enddo
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'configurations: overlapping ', n_overlap, '  touching ', n_contact, '  apart ', n_apart, &
      '  with coplanar face planes ', n_coplanar
print '(A,I0,A,I0,A,I0)', 'booleans: ', 4 * N_CASES, '  status OK: ', n_ok_status, '  non-OK status: ', n_bad_status
print '(A,I0)',           'silent wrong volumes (status OK, volume off): ', n_wrong
print '(A,I0,A,I0)',      'watertight: ', n_wt, '  manifold: ', n_mf
print '(A,I0,A,I0)',      'topology required (no contact, U/I/D, non-empty): ', n_topo_req, '  failing: ', n_topo_bad
print '(A,L1)', 'Are all tests passed? ', n_wrong == 0 .and. n_topo_bad == 0
if (.not. (n_wrong == 0 .and. n_topo_bad == 0)) error stop 1

contains
   function rnd() result(r)
   !< Fixed-seed linear congruential generator in [0, 1) (portable, no intrinsic seeding).
   real(R8P) :: r

   seed = mod(seed * 6364136223846793005_I8P + 1442695040888963407_I8P, huge(seed))
   if (seed < 0_I8P) seed = -seed
   r = real(mod(seed, 1000003_I8P), R8P) / 1000003._R8P
   endfunction rnd

   subroutine one(k)
   integer(I4P), intent(in) :: k
   type(surface_stl_object) :: a, b
   real(R8P)                :: err

   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%resize(x=s(1), y=s(2), z=s(3))
   call b%translate(x=t(1), y=t(2), z=t(3))
   if (rotated) then
      call a%rotate(axis=vector_R8P(axis(1), axis(2), axis(3)), angle=angle)
      call b%rotate(axis=vector_R8P(axis(1), axis(2), axis(3)), angle=angle)
   endif
   call a%analyze ; call b%analyze
   call a%boolean(other=b, op=OPS(k), status=st)
   if (st /= BOOL_STATUS_OK) then
      n_bad_status = n_bad_status + 1
      print '(A,I0,A,I0,A,I0,A,L1)', 'non-OK status: case ', c, ' op ', OPS(k), ' status ', st, ' rotated ', rotated
   else
      n_ok_status = n_ok_status + 1
      v = a%get_volume()
      err = abs(v - expected(k)) / max(expected(k), 1._R8P)
      if (err > 1.e-10_R8P) then
         n_wrong = n_wrong + 1
         print '(A,I0,A,I0,A,F10.6,A,F10.6,A,3F6.2,A,3F6.2,A,L1)', 'WRONG volume: case ', c, ' op ', OPS(k), ' got ', v, &
               ' expected ', expected(k), ' s ', s, ' t ', t, ' rotated ', rotated
         if (rotated) print '(A,3ES24.16,A,ES24.16)', '   rotation axis ', axis, ' angle ', angle
      endif
      wt = a%is_watertight() ; mf = a%is_manifold()
      if (wt) n_wt = n_wt + 1
      if (mf) n_mf = n_mf + 1
      if (.not. contact .and. k /= 4 .and. expected(k) > 0._R8P) then
         n_topo_req = n_topo_req + 1
         if (.not. (wt .and. mf)) then
            n_topo_bad = n_topo_bad + 1
            print '(A,I0,A,I0,A,2L1,A,3F6.2,A,3F6.2,A,L1)', 'topology: case ', c, ' op ', OPS(k), ' wt/mf ', wt, mf, &
                  ' s ', s, ' t ', t, ' rotated ', rotated
            if (rotated) print '(A,3ES24.16,A,ES24.16)', '   rotation axis ', axis, ' angle ', angle
         endif
      endif
   endif
   call a%destroy ; call b%destroy
   endsubroutine one
endprogram fossil_test_boolean_stress
