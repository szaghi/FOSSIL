!< FOSSIL test: boolean CSG chains (issue #63, roadmap #58 WP5, step 4).
!<
!< Short chains of fossil_boolean_csg: the unit cube drilled by 4 cylinders and the union of 4 overlapping spheres
!< (marching cubes, 10^3). Every step is checked against volume and second moments from independent booleans; every
!< intermediate body must be watertight and manifold. The long chains (9 to 49 operations) are
!< src/tests/validation/fossil_val_boolean_csg.f90.
program fossil_test_boolean_csg

use fossil_boolean_csg, only : chain_report_t, drill_chain, lattice_chain
use penf,               only : I4P, R8P

implicit none

type(chain_report_t) :: rep
logical              :: ok

ok = .true.
call drill_chain(2_I4P, 10_I4P, rep, .false.)
call summary('drilled block, 4 cylinders')
call lattice_chain(2_I4P, 10_I4P, rep, .false.)
call summary('sphere lattice, 4 spheres ')
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine summary(label)
   !< One line per chain; fold into the verdict.
   character(*), intent(in) :: label
   logical                  :: pass

   pass = rep%bad_status == 0 .and. rep%max_vol <= 1e-9_R8P .and. rep%max_mom <= 1e-9_R8P .and. &
          rep%final_vol <= 1e-9_R8P .and. rep%final_mom <= 1e-9_R8P .and. rep%valid
   print '(2A,I0,A,ES8.1,A,ES8.1,A,L1,A,L1)', label, ': steps ', rep%steps, ', vol ', max(rep%max_vol, rep%final_vol), &
         ', mom ', max(rep%max_mom, rep%final_mom), ', valid ', rep%valid, ', pass ', pass
   ok = ok .and. pass
   endsubroutine summary
endprogram fossil_test_boolean_csg
