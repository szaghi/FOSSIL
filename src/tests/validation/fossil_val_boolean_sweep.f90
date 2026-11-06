!< FOSSIL validation: dense degenerate sweep of the exact boolean (issue #65, step 6).
!<
!< Bodies A, each centred at the origin, at two resolutions: the unit cube (cube.stl, 12 facets; cube-02.stl), the
!< sphere (marching cubes, 10^3 and 24^3 samples), the torus (12^3 and 24^3). B is A translated along (0.31, 0.53, 0.79) by
!< v times the size, or rotated about (0.3, -0.5, 0.8) by v rad, for v = 1e-4 .. 1e-14 in steps of 10^0.5 (21
!< values); and, for v <= 1e-10, the same with the opt-in weld (weld = 1e-9 of the size: B's vertices take A's). The
!< harness (fossil_boolean_harness) checks every identity of the four ops, validity included (self-intersections only
!< as rounding folds). Pass criterion: every case passes with status OK; a non-OK status is counted and reported, and
!< fails the program too (none is expected). One line per failing case, then the counts per body and motion.
!< The direction is generic on purpose: along (1, 2, 3), parallel to the marching-cubes facets of normal (-1, -1, 1), a
!< translated edge of B can stay exactly in a facet plane of A (rounding keeps it there), and B then touches A along
!< it without crossing: the difference pinches there, non-manifold by definition, like any contact (measured on the
!< 10^3 sphere at 3.2e-5, 3.2e-6, 3.2e-7; the exact topology is right).
program fossil_val_boolean_sweep

use fossil,                 only : surface_stl_object, BOOL_STATUS_OK
use fossil_boolean_bodies,  only : BODY_SPHERE, BODY_TORUS, mc_body
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

character(6), parameter  :: BNAME(3) = ['cube  ', 'sphere', 'torus ']
type(surface_stl_object) :: a, b
type(vector_R8P)         :: dir, axis
real(R8P)                :: v, t0, t1
integer(I4P)             :: body, res, motion, i, w, nall, npass, nnonok, ncase, ncase_pass
character(80)            :: label
logical                  :: all_ok

dir = vector_R8P(0.31_R8P, 0.53_R8P, 0.79_R8P) ; dir = dir / dir%normL2()
axis = vector_R8P(0.3_R8P, -0.5_R8P, 0.8_R8P) ; axis = axis / axis%normL2()
nall = 0 ; npass = 0 ; nnonok = 0
call cpu_time(t0)
do body = 1, 3
   do res = 1, 2
      call make(body, res, a)
      do motion = 1, 2
         ncase = 0 ; ncase_pass = 0
         do i = 0, 20
            v = 10._R8P ** (-4._R8P - 0.5_R8P * real(i, R8P))
            do w = 0, 1
               if (w == 1 .and. v > 1.0001e-10_R8P) cycle
               call make(body, res, b)
               if (motion == 1) then
                  call b%translate(delta=v * dir)
               else
                  call b%rotate(axis=axis, angle=v)
               endif
               call b%analyze
               write(label, '(A,A,I0,A,A,ES8.1,A)') trim(BNAME(body)), ' res ', res, merge(' translated ', ' rotated    ', &
                     motion == 1), ' v ', v, merge(' weld', '     ', w == 1)
               call run(trim(label), w == 1)
            enddo
         enddo
         print '(A,A,I0,A,A,I0,A,I0)', BNAME(body), ' res ', res, merge(' translated', ' rotated   ', motion == 1), &
               ': passing ', ncase_pass, ' of ', ncase
      enddo
   enddo
enddo
call cpu_time(t1)
all_ok = npass == nall .and. nnonok == 0
print '(A,I0,A,I0,A,I0,A,F8.1,A)', 'cases ', nall, ', passing ', npass, ', non-OK statuses ', nnonok, ', ', t1 - t0, ' s'
print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine make(kind, r, s)
   !< Body kind at resolution r, centred at the origin, size ~1.
   integer(I4P),             intent(in)    :: kind, r !< Body, resolution.
   type(surface_stl_object), intent(inout) :: s       !< Surface.
   type(surface_stl_object)                :: fresh   !< Reset.

   s = fresh
   select case (kind)
   case (1)
      call s%load_from_file(file_name=trim(merge('src/tests/cube.stl   ', 'src/tests/cube-02.stl', r == 1)), &
                            guess_format=.true.)
      call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
      call s%analyze
   case (2)
      call mc_body(BODY_SPHERE, merge(10_I4P, 24_I4P, r == 1), s)
   case default
      call mc_body(BODY_TORUS, merge(12_I4P, 24_I4P, r == 1), s)
   endselect
   endsubroutine make

   subroutine run(lab, welded)
   !< One case: harness, count, print if it fails.
   character(*), intent(in) :: lab    !< Case name.
   logical,      intent(in) :: welded !< With the opt-in weld.
   type(boolean_report_t)   :: rep    !< Report.
   logical                  :: ok     !< Pass.

   if (welded) then
      call check_boolean(a, b, rep, 100_I4P, 20261106_I8P, .true., .false., weld=1.e-9_R8P)
   else
      call check_boolean(a, b, rep, 100_I4P, 20261106_I8P, .true., .false.)
   endif
   ok = boolean_pass(rep) .and. all(rep%status == BOOL_STATUS_OK)
   if (any(rep%status /= BOOL_STATUS_OK)) nnonok = nnonok + 1
   nall = nall + 1 ; ncase = ncase + 1
   if (ok) then
      npass = npass + 1 ; ncase_pass = ncase_pass + 1
   else
      call print_report('   '//lab, rep)
      print '(A,4I3,A,4L2,A,4L2,A,4I5,A,4ES9.1)', '      status', rep%status, ' wt', rep%watertight, ' mf', &
            rep%manifold, ' self', rep%self_pairs, ' depth', rep%fold_depth
   endif
   endsubroutine run
endprogram fossil_val_boolean_sweep
