!< FOSSIL validation: boolean CSG chains (issue #63, roadmap #58 WP5, step 4).
!<
!< The chains of fossil_boolean_csg at three lengths: the unit cube drilled by 9, 25 and 49 cylinders, and the union of
!< 9, 25 and 49 overlapping spheres (marching cubes, 12^3). Each step is checked against volume and second moments
!< computed from independent booleans; drift shows as a final residual growing with the length. Pass: every boolean
!< status OK, every step and final residual <= 1e-9, every intermediate body watertight and manifold.
program fossil_val_boolean_csg

use fossil_boolean_csg, only : chain_report_t, drill_chain, lattice_chain
use penf,               only : I4P, R8P

implicit none

integer(I4P), parameter :: KS(3) = [3_I4P, 5_I4P, 7_I4P]
type(chain_report_t)    :: rep
integer(I4P)            :: i
logical                 :: ok

ok = .true.
do i = 1, size(KS)
   print '(A,I0,A)', 'drilled block, ', KS(i) ** 2, ' cylinders'
   call drill_chain(KS(i), 12_I4P, rep, .true.)
   call summary
enddo
do i = 1, size(KS)
   print '(A,I0,A)', 'sphere lattice, ', KS(i) ** 2, ' spheres'
   call lattice_chain(KS(i), 12_I4P, rep, .true.)
   call summary
enddo
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine summary()
   !< One line per chain; fold into the verdict.
   logical :: pass

   pass = rep%bad_status == 0 .and. rep%max_vol <= 1e-9_R8P .and. rep%max_mom <= 1e-9_R8P .and. &
          rep%final_vol <= 1e-9_R8P .and. rep%final_mom <= 1e-9_R8P .and. rep%valid
   print '(A,I0,A,I0,A,ES8.1,A,ES8.1,A,ES8.1,A,ES8.1,A,L1,A,I0,A,F8.2,A,L1)', '   steps ', rep%steps, ', non-OK ', &
         rep%bad_status, ', max step vol ', rep%max_vol, ' mom ', rep%max_mom, ', final vol ', rep%final_vol, ' mom ', &
         rep%final_mom, ', valid ', rep%valid, ', facets ', rep%facets, ', ', rep%seconds, ' s, pass ', pass
   ok = ok .and. pass
   endsubroutine summary
endprogram fossil_val_boolean_csg
