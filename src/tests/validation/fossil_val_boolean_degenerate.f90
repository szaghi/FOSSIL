!< FOSSIL validation: boolean identities on near-degenerate configurations (issue #63, roadmap #58 WP5, step 3).
!<
!< Two bodies A: the unit cube and a marching-cubes sphere (24^3 samples), both centred at the origin. Families:
!< 1. A against A (coincident everywhere);
!< 2. A against A translated by 1e-4, 1e-6, 1e-8, 1e-10, 1e-12 and 1e-14 of the size along a generic direction;
!< 3. A against A rotated about a generic axis through the origin by 1e-4, 1e-6, 1e-8, 1e-10 and 1e-12 rad;
!< 4. nested: A against A scaled by 0.5 (strictly inside), and the cube against a half cube sharing part of a face;
!< 5. contacts of the cube with a translated cube: at one vertex, along one edge, on a face patch;
!< 6. coincident curved faces: the sphere against the same sphere re-facetted (25^3 samples).
!< The harness (fossil_boolean_harness) checks every identity of the four ops. One line per case, then the outcome map:
!< OK statuses, passes, and for a failing case the identities it misses. Every case must pass except the curved body
!< against itself moved by 1e-10 .. 1e-6 (translated) or 1e-10 .. 1e-4 rad (rotated), the documented limit of
!< booleans.md: features between the 1e-10 snap and the 1e-9 tolerances of the retriangulation (measured, not gated).
program fossil_val_boolean_degenerate

use fossil,                 only : surface_stl_object, BOOL_STATUS_OK
use fossil_boolean_bodies,  only : BODY_SPHERE, mc_body
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report, TOL_IDENTITY
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

real(R8P),    parameter  :: OFFSETS(6) = [1e-4_R8P, 1e-6_R8P, 1e-8_R8P, 1e-10_R8P, 1e-12_R8P, 1e-14_R8P]
real(R8P),    parameter  :: ANGLES(5) = [1e-4_R8P, 1e-6_R8P, 1e-8_R8P, 1e-10_R8P, 1e-12_R8P]
type(surface_stl_object) :: a, b
type(vector_R8P)         :: dir, axis
integer(I4P)             :: body, i, nall, npass, nmust, nmust_pass
character(64)            :: label
character(6)             :: bname

nall = 0 ; npass = 0 ; nmust = 0 ; nmust_pass = 0
dir = vector_R8P(1._R8P, 2._R8P, 3._R8P) ; dir = dir / dir%normL2()
axis = vector_R8P(0.3_R8P, -0.5_R8P, 0.8_R8P) ; axis = axis / axis%normL2()
do body = 1, 2
   bname = merge('cube  ', 'sphere', body == 1)
   ! 1. coincident
   call make(body, a) ; b = a
   call run(trim(bname)//' x itself', a, b, .false., .true., .false.)
   ! 2. translated copies
   do i = 1, size(OFFSETS)
      call make(body, b) ; call b%translate(delta=OFFSETS(i) * dir) ; call b%analyze
      write(label, '(2A,ES7.0)') trim(bname), ' translated ', OFFSETS(i)
      call run(label, a, b, OFFSETS(i) >= 1e-8_R8P, .false., body == 1 .or. OFFSETS(i) >= 1e-6_R8P .or. OFFSETS(i) < 1e-10_R8P)
   enddo
   ! 3. rotated copies
   do i = 1, size(ANGLES)
      call make(body, b) ; call b%rotate(axis=axis, angle=ANGLES(i)) ; call b%analyze
      write(label, '(2A,ES7.0)') trim(bname), ' rotated ', ANGLES(i)
      call run(label, a, b, ANGLES(i) >= 1e-8_R8P, .false., body == 1 .or. ANGLES(i) >= 1e-4_R8P .or. ANGLES(i) < 1e-10_R8P)
   enddo
   ! 4. nested, strictly inside
   call make(body, b) ; call b%resize(factor=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), respect_centroid=.true.)
   call b%analyze
   call run(trim(bname)//' nested', a, b, .true., .false., .true.)
enddo
! 4. nested with a shared face patch: half cube on the face x = 0.5 of the cube
call make(1_I4P, a)
call make(1_I4P, b) ; call b%resize(factor=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), respect_centroid=.true.)
call b%translate(x=0.25_R8P) ; call b%analyze
call run('cube nested, shared face', a, b, .false., .true., .true.)
! 5. contacts
call make(1_I4P, b) ; call b%translate(x=1._R8P, y=1._R8P, z=1._R8P) ; call b%analyze
call run('cube contact vertex', a, b, .true., .true., .true.)
call make(1_I4P, b) ; call b%translate(x=1._R8P, y=1._R8P) ; call b%analyze
call run('cube contact edge', a, b, .true., .true., .true.)
call make(1_I4P, b) ; call b%translate(x=1._R8P, y=0.3_R8P, z=0.2_R8P) ; call b%analyze
call run('cube contact face patch', a, b, .false., .true., .true.)
! 6. re-facetted sphere
call make(2_I4P, a) ; call mc_body(BODY_SPHERE, 25_I4P, b)
call run('sphere x sphere re-facetted', a, b, .true., .false., .true.)

print '(A,I0,A,I0,A,I0,A,I0)', 'cases ', nall, ', passing ', npass, '; required ', nmust, ', passing ', nmust_pass
print '(A,L1)', 'Are all tests passed? ', nmust_pass == nmust

contains
   subroutine make(kind, s)
   !< The unit cube centred at the origin (kind 1) or the 24^3 sphere (kind 2).
   integer(I4P),             intent(in)    :: kind !< Body.
   type(surface_stl_object), intent(inout) :: s    !< Surface.
   type(surface_stl_object)                :: fresh !< Reset.

   s = fresh
   if (kind == 1_I4P) then
      call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
      call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
      call s%analyze
   else
      call mc_body(BODY_SPHERE, 24_I4P, s)
   endif
   endsubroutine make

   subroutine run(lab, a, b, area_identity, contact, required)
   !< Check one case, print it, count it; on failure, name the identities it misses.
   character(*),             intent(in) :: lab           !< Case name.
   type(surface_stl_object), intent(in) :: a, b          !< Inputs.
   logical,                  intent(in) :: area_identity !< Check the area identity.
   logical,                  intent(in) :: contact       !< The inputs touch.
   logical,                  intent(in) :: required      !< Must pass (outside the snap's range).
   type(boolean_report_t)               :: rep           !< Report.
   character(64)                        :: lab64         !< Padded label.
   integer(I4P)                         :: k             !< Counter.
   logical                              :: ok            !< Pass.

   call check_boolean(a, b, rep, 200_I4P, 20261106_I8P, area_identity, contact)
   lab64 = lab
   call print_report('   '//lab64(1:34), rep)
   ok = boolean_pass(rep)
   nall = nall + 1 ; if (ok) npass = npass + 1
   if (required) then
      nmust = nmust + 1 ; if (ok) nmust_pass = nmust_pass + 1
   endif
   if (ok) return
   if (any(rep%vol_res > TOL_IDENTITY)) print '(A,3ES9.1)', '      volume residuals ', rep%vol_res
   if (any(rep%mom_res > TOL_IDENTITY)) print '(A,3ES9.1)', '      moment residuals ', rep%mom_res
   if (rep%area_res > TOL_IDENTITY) print '(A,ES9.1)', '      area residual ', rep%area_res
   do k = 1, 4
      if (rep%status(k) /= BOOL_STATUS_OK) cycle
      if (rep%member_bad(k) > 0 .or. rep%dev_lower(k) > TOL_IDENTITY) &
         print '(A,I0,A,I0,A,ES9.1)', '      op ', k, ': membership mismatches ', rep%member_bad(k), ', deviation ', &
                                     rep%dev_lower(k)
      if (rep%validity_checked(k)) then
         if (.not. (rep%watertight(k) .and. rep%manifold(k) .and. rep%volume(k) > 0._R8P .and. rep%self_pairs(k) == 0)) &
            print '(A,I0,A,2L2,ES10.2,A,I0)', '      op ', k, ': watertight, manifold, V ', rep%watertight(k), &
                                            rep%manifold(k), rep%volume(k), ', self pairs ', rep%self_pairs(k)
      endif
   enddo
   endsubroutine run
endprogram fossil_val_boolean_degenerate
