!< FOSSIL, test the SIGN_FAST_WINDING_NUMBER sign algorithm (issue #29, step 5).
!<
!<   1. Watertight cube-01, cube-02 and a rotated cube-02: the FWN sign equals the pseudo-normal
!<      sign at every 13^3 grid point farther than 1e-9 from the surface.
!<   2. Holed cube-02 (16 facets removed): the FWN sign is correct at all 729 interior/exterior
!<      grid points. For reference the test prints the errors of the other algorithms on the
!<      same points (measured: pseudo-normal 10, ray parity 15, solid angle 0).
!<   3. Inside-out cube put through sanitize: FWN sign correct at the centre.
!<   4. distance_many with FWN equals the serial distance loop bit for bit, and
!<      is_point_inside(FWN) agrees with the sign of distance(FWN), on dragon (it replaced
!<      dragon-fine, issue #45).

program fossil_test_fwn_sign

use fossil, only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                   SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

logical :: ok(6)

ok(1) = check_watertight('src/tests/cube-01.stl', .false., 'cube-01        ')
ok(2) = check_watertight('src/tests/cube-02.stl', .false., 'cube-02        ')
ok(3) = check_watertight('src/tests/cube-02.stl', .true.,  'cube-02 rotated')
ok(4) = check_holed()
ok(5) = check_inside_out()
ok(6) = check_batch()

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check_watertight(file_name, rotated, label) result(passed)
   character(*), intent(in) :: file_name, label
   logical,      intent(in) :: rotated
   logical                  :: passed
   type(surface_stl_object) :: s
   type(vector_R8P)         :: bmin, h, p
   real(R8P)                :: d_pn, d_fwn
   integer(I4P)             :: i, j, k, n_checked, n_mismatch

   call s%load_from_file(file_name=file_name, guess_format=.true.)
   call s%sanitize
   if (rotated) then
      call s%rotate(axis=ex_R8P + 0.7_R8P * ey_R8P + 0.3_R8P * ez_R8P, angle=0.61_R8P)
      call s%analyze
   endif
   bmin = s%get_bmin() ; h = s%get_bmax() - bmin
   bmin = bmin - 0.25_R8P * h ; h = 1.5_R8P * h / 12._R8P
   n_checked = 0 ; n_mismatch = 0
   do k = 0, 12 ; do j = 0, 12 ; do i = 0, 12
      p = vector_R8P(bmin%x + i * h%x, bmin%y + j * h%y, bmin%z + k * h%z)
      d_pn = s%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_PSEUDO_NORMAL)
      if (abs(d_pn) < 1.e-9_R8P) cycle
      d_fwn = s%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
      n_checked = n_checked + 1
      if (sign(1._R8P, d_pn) /= sign(1._R8P, d_fwn)) n_mismatch = n_mismatch + 1
   enddo ; enddo ; enddo
   passed = (n_mismatch == 0) .and. (n_checked > 0)
   print '(A,A,I0,A,I0,A,L1)', label, ': points=', n_checked, ' FWN vs PN sign mismatches=', n_mismatch, ' ok=', passed
   endfunction check_watertight

   function check_holed() result(passed)
   logical                         :: passed
   type(surface_stl_object)        :: s, holed
   type(facet_object), allocatable :: f(:)
   type(facet_object), pointer     :: fp
   type(vector_R8P)                :: p
   integer(I4P), parameter         :: N_REMOVED = 16_I4P
   integer(I4P)                    :: algos(4), errors(4), i, j, k, a, n_checked
   real(R8P)                       :: d
   logical                         :: truth

   algos = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_FAST_WINDING_NUMBER]
   call s%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   call s%sanitize
   allocate(f(s%get_facets_number() - N_REMOVED))
   do i = 1, size(f)
      fp => s%facet_at(i + N_REMOVED)
      f(i) = fp
   enddo
   call holed%adopt_facets(facets=f)
   errors = 0 ; n_checked = 0
   do k = 1, 9 ; do j = 1, 9 ; do i = 1, 9
      p = vector_R8P(-0.25_R8P + 0.15_R8P * i, -0.25_R8P + 0.15_R8P * j, -0.25_R8P + 0.15_R8P * k)
      truth = all([p%x, p%y, p%z] > 0._R8P) .and. all([p%x, p%y, p%z] < 1._R8P)
      n_checked = n_checked + 1
      do a = 1, 4
         d = holed%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=algos(a))
         if ((d < 0._R8P) .neqv. truth) errors(a) = errors(a) + 1
      enddo
   enddo ; enddo ; enddo
   passed = (errors(4) == 0)
   print '(A,I0,A,I0,A,4(I0,1X),A,L1)', 'cube-02 with ', N_REMOVED, ' facets removed, ', n_checked, &
         ' points; wrong signs PN RI SA FWN: ', errors, ' ok=', passed
   endfunction check_holed

   function check_inside_out() result(passed)
   logical                  :: passed
   type(surface_stl_object) :: s
   real(R8P)                :: d

   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%reverse_normals
   call s%analyze
   call s%sanitize
   d = s%distance(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), is_signed=.true., is_square_root=.true., &
                  sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   passed = abs(d + 0.5_R8P) < 1.e-12_R8P
   print '(A,F8.4,A,L1)', 'inside-out cube after sanitize: d(centre) by FWN = ', d, ' ok=', passed
   endfunction check_inside_out

   function check_batch() result(passed)
   logical                       :: passed
   type(surface_stl_object)      :: s
   type(vector_R8P), allocatable :: pts(:)
   real(R8P),        allocatable :: d_many(:)
   real(R8P)                     :: d_one
   integer(I4P)                  :: i, n_diff, n_inside_diff, status

   call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
   call s%sanitize
   allocate(pts(1000), d_many(1000))
   do i = 1, 1000
      pts(i) = s%get_bmin() + (s%get_bmax() - s%get_bmin()) * (real(mod(i * 37, 1000), R8P) / 1000._R8P)
      pts(i)%y = pts(i)%y + 0.3_R8P * real(mod(i * 11, 7), R8P)
   enddo
   call s%distance_many(points=pts, distances=d_many, is_signed=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER, &
                        is_square_root=.true., status=status)
   n_diff = 0 ; n_inside_diff = 0
   do i = 1, 1000
      d_one = s%distance(point=pts(i), is_signed=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER, is_square_root=.true.)
      if (d_one /= d_many(i)) n_diff = n_diff + 1
      if (s%is_point_inside(point=pts(i), sign_algorithm=SIGN_FAST_WINDING_NUMBER) .neqv. (d_one < 0._R8P)) &
         n_inside_diff = n_inside_diff + 1
   enddo
   passed = (n_diff == 0) .and. (n_inside_diff == 0) .and. (status == 0)
   print '(A,I0,A,I0,A,L1)', 'dragon 1000 points: distance_many vs serial diffs=', n_diff, &
         '  is_point_inside vs sign diffs=', n_inside_diff, ' ok=', passed
   endfunction check_batch
endprogram fossil_test_fwn_sign
