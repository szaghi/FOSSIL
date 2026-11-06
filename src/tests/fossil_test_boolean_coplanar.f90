!< FOSSIL test: mesh booleans on coplanar overlaps (issue #33).
program fossil_test_boolean_coplanar
!< FOSSIL test: mesh booleans on coplanar overlaps (issue #33).
!<
!< A = unit cube, B = a box at an offset; the four ops against analytic volumes, axis-aligned
!< and under a proper rotation of both bodies (axis (1, 2, 3), angle 0.7):
!<  - offset x 0.5 (B's y and z faces coplanar with A's);
!<  - offset x and y 0.5;
!<  - face touching (offset x 1);
!<  - identical cubes;
!<  - a 1 x 0.5 x 0.5 box on A's z = 0 plane, partially overlapping it;
!<  - generic offset (0.5, 0.5, 0.5), no coplanar faces (control).
!< A check passes when the status is OK, the volume is within 1e-10 (relative to max(V, 1)) of
!< the analytic one, and a non-empty result is watertight and manifold. Exception: a symmetric
!< difference A^B = (A\B) u (B\A) whose two pieces touch along the curves where the surfaces cross
!< (offset xy, box on face plane, generic) is non-manifold by definition (four facets at each
!< such edge), so there only its volume is checked. So is the rotated offset x since issue #65: its
!< coplanar faces are off by round-off once rotated, and the exact arrangement crosses them (no snap).
!<
!< Started as fossil_xfail_boolean_coplanar (4 of 48 checks passing); a regular test since the
!< coplanar cuts, the exact pair classification and the T-junction closing made all 48 pass.
use fossil, only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
integer(I4P), parameter :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
character(1), parameter :: OPN(4) = ['U', 'I', 'D', 'S']
character(256) :: mesh
integer(I4P) :: n_ok, n_all, r
n_ok = 0 ; n_all = 0
mesh = 'src/tests/cube.stl'
if (command_argument_count() >= 1) call get_command_argument(1, mesh)
do r = 0, 1
   ! A = unit cube; B = box of size s at offset t; expected volumes [U, I, D, S]
   call run('offset x 0.5', r == 0, &
            [0.5_R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], &
            [1.5_R8P, 0.5_R8P, 0.5_R8P, 1._R8P], r == 1)
   call run('offset xy 0.5', .false., &
            [0.5_R8P, 0.5_R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], &
            [1.75_R8P, 0.25_R8P, 0.75_R8P, 1.5_R8P], r == 1)
   call run('face touching', .true., &
            [1._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], &
            [2._R8P, 0._R8P, 1._R8P, 2._R8P], r == 1)
   call run('identical', .true., &
            [0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], &
            [1._R8P, 1._R8P, 0._R8P, 0._R8P], r == 1)
   call run('box on face plane', .false., &
            [0.5_R8P, 0.25_R8P, 0._R8P], [1._R8P, 0.5_R8P, 0.5_R8P], &
            [1.125_R8P, 0.125_R8P, 0.875_R8P, 1._R8P], r == 1)
   call run('generic xyz 0.5', .false., &
            [0.5_R8P, 0.5_R8P, 0.5_R8P], [1._R8P, 1._R8P, 1._R8P], &
            [1.875_R8P, 0.125_R8P, 0.875_R8P, 1.75_R8P], r == 1)
enddo
print '(A,I0,A,I0)', 'checks passed: ', n_ok, ' of ', n_all
print '(A,L1)', 'Are all tests passed? ', n_ok == n_all
if (n_ok /= n_all) error stop 1
contains
   subroutine run(label, symdiff_manifold, t, s, expected, rotated)
   character(*), intent(in) :: label
   logical,      intent(in) :: symdiff_manifold !< A\B and B\A disjoint or sharing whole faces.
   real(R8P),    intent(in) :: t(3), s(3), expected(4)
   logical,      intent(in) :: rotated
   type(surface_stl_object) :: a, b
   integer(I4P) :: k, st
   real(R8P)    :: v, err
   logical      :: ok, wt, mf
   character(64) :: line
   line = ''
   do k = 1, 4
      call a%load_from_file(file_name=trim(mesh), guess_format=.true.)
      call b%load_from_file(file_name=trim(mesh), guess_format=.true.)
      call b%resize(x=s(1), y=s(2), z=s(3))
      call b%translate(x=t(1), y=t(2), z=t(3))
      if (rotated) then
         call a%rotate(axis=vector_R8P(1._R8P, 2._R8P, 3._R8P), angle=0.7_R8P)
         call b%rotate(axis=vector_R8P(1._R8P, 2._R8P, 3._R8P), angle=0.7_R8P)
      endif
      call a%analyze ; call b%analyze
      call a%boolean(other=b, op=OPS(k), status=st)
      v = a%get_volume()
      err = abs(v - expected(k)) / max(expected(k), 1._R8P)
      wt = a%is_watertight() ; mf = a%is_manifold()
      ok = st == BOOL_STATUS_OK .and. err <= 1.e-10_R8P .and. &
           (expected(k) == 0._R8P .or. (wt .and. mf) .or. (OPS(k) == BOOL_SYMDIFF .and. .not. symdiff_manifold))
      if (expected(k) == 0._R8P) ok = st == BOOL_STATUS_OK .and. abs(v) <= 1.e-10_R8P
      n_all = n_all + 1 ; if (ok) n_ok = n_ok + 1
      write(line(1+15*(k-1):), '(A1,A1,F9.5,A1,L1,L1)') OPN(k), merge('+', 'x', ok), v, ' ', wt, mf
      call a%destroy ; call b%destroy
   enddo
   print '(A,A,A,I2,A,A)', merge('rot ', '    ', rotated), label, repeat(' ', 20 - len(label)), st, '  ', trim(line)
   endsubroutine run
endprogram fossil_test_boolean_coplanar
