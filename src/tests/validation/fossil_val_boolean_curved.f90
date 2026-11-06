!< FOSSIL validation: boolean identities on curved and dense bodies (issue #63, roadmap #58 WP5, step 2).
!<
!< 1. Marching-cubes bodies (sphere, torus, capped cylinder) and the unit cube, every pair of kinds (10), at 24^3 and
!<    40^3 samples, each under 4 fixed-seed proper rigid motions of B (random axis and angle, translation components within +-0.2):
!<    80 pairs, 320 booleans.
!< 2. The Stanford bunny and the dragon, not watertight, enter through alpha_wrap (alpha 1/50, offset 1/500 of the
!<    diagonal), scaled to a unit diagonal and centred; 4 motions of the dragon against the bunny.
!< The harness (fossil_boolean_harness) checks every identity of the four ops; prints one line per pair, the counts of
!< non-OK statuses and failing pairs, and the time. Pass criterion (#33's): no pair with status OK fails an identity.
program fossil_val_boolean_curved

use fossil,                 only : surface_stl_object, BOOL_STATUS_OK
use fossil_boolean_bodies,  only : BODY_SPHERE, BODY_TORUS, BODY_CYLINDER, mc_body, move_body
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

integer(I4P), parameter  :: RES(2) = [24_I4P, 40_I4P]
character(8), parameter  :: NAMES(0:3) = ['cube    ', 'sphere  ', 'torus   ', 'cylinder']
type(surface_stl_object) :: a, b, bunny, dragon
type(boolean_report_t)   :: rep
real(R8P)                :: axis(3), angle, t(3)
integer(I8P)             :: seed, c0, c1, rate
integer(I4P)             :: r, i, j, m, npass, nall, nbad_status
character(48)            :: label

call system_clock(c0, rate)
seed = 20261106_I8P ; npass = 0 ; nall = 0 ; nbad_status = 0
do r = 1, size(RES)
   do i = 0, 3
      do j = i, 3
         do m = 1, 4
            call body(i, RES(r), a) ; call body(j, RES(r), b)
            call motion
            call move_body(b, axis, angle, t)
            call check_boolean(a, b, rep, 200_I4P, seed, .true., .false.)
            write(label, '(I0,A,3A,I0)') RES(r), '^3 ', trim(NAMES(i)), ' x ', trim(NAMES(j)) // ' m', m
            call tally(label)
         enddo
      enddo
   enddo
enddo
call wrapped('src/tests/bunny.stl', bunny) ; call wrapped('src/tests/dragon.stl', dragon)
do m = 1, 4
   a = bunny ; b = dragon
   call motion
   call move_body(b, axis, angle, t)
   call check_boolean(a, b, rep, 200_I4P, seed, .true., .false.)
   write(label, '(A,I0)') 'wrapped bunny x dragon m', m
   call tally(label)
enddo
call system_clock(c1)
print '(A,I0,A,I0,A,I0,A,F8.1,A)', 'pairs ', nall, ', passing ', npass, ', non-OK statuses ', nbad_status, ', ', &
                                   real(c1 - c0, R8P) / rate, ' s'
print '(A,L1)', 'Are all tests passed? ', npass == nall
if (npass /= nall) error stop 1

contains
   subroutine tally(lab)
   !< Print the report, count.
   character(*), intent(in) :: lab

   call print_report('   '//lab(1:34), rep)
   nall = nall + 1 ; if (boolean_pass(rep)) npass = npass + 1
   nbad_status = nbad_status + count(rep%status /= BOOL_STATUS_OK)
   endsubroutine tally

   subroutine motion()
   !< A fixed-seed proper rigid motion.
   axis = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] ; angle = 6.283185307179586_R8P * rnd()
   t = 0.4_R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
   endsubroutine motion

   subroutine body(kind, n, s)
   !< A marching-cubes body, or the unit cube centred at the origin (kind 0).
   integer(I4P),             intent(in)    :: kind, n
   type(surface_stl_object), intent(inout) :: s

   if (kind == 0) then
      call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
      call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
      call s%analyze
   else
      call mc_body(kind, n, s)
   endif
   endsubroutine body

   subroutine wrapped(file_name, s)
   !< The alpha wrap of a mesh, scaled to a unit diagonal and centred at the origin.
   character(*),             intent(in)    :: file_name
   type(surface_stl_object), intent(inout) :: s
   type(surface_stl_object)                :: raw
   type(vector_R8P)                        :: lo, hi
   real(R8P)                               :: d
   integer(I4P)                            :: st

   call raw%load_from_file(file_name=file_name, guess_format=.true.)
   lo = raw%get_bmin() ; hi = raw%get_bmax()
   d = norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call raw%alpha_wrap(alpha=d / 50._R8P, offset=d / 500._R8P, wrapped=s, status=st)
   if (.not. s%is_watertight()) error stop 'alpha wrap not closed'
   lo = s%get_bmin() ; hi = s%get_bmax()
   call s%translate(x=-0.5_R8P * (lo%x + hi%x), y=-0.5_R8P * (lo%y + hi%y), z=-0.5_R8P * (lo%z + hi%z))
   call s%resize(factor=vector_R8P(1._R8P / d, 1._R8P / d, 1._R8P / d))
   call s%analyze
   print '(3A,I0,A)', '   ', file_name, ' wrapped: ', s%get_facets_number(), ' facets'
   endsubroutine wrapped

   function rnd() result(x)
   !< Park-Miller in (0, 1).
   real(R8P) :: x

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   x = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_val_boolean_curved
