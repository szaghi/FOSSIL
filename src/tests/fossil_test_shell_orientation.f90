!< FOSSIL, assert sanitize orients every closed shell by its nesting parity.
!<
!< sanitize_normals used to decide "outward" by one volume check over the whole surface,
!< so a shell oriented against the rest stayed wrong: two disjoint cubes, one inverted,
!< sum to volume 0 and nothing was flipped. Each closed shell is now oriented by nesting
!< parity: an outer boundary encloses positive volume, a cavity wall negative. Cases
!< (unit cube U = [0,1]^3, big cube B = [-1,2]^3), each merged, analyzed, sanitized:
!<   A. U and U translated by 5 along x and inverted: both outward, volume 2.
!<   B. hollow body, correct (B outward, U inward): unchanged, volume 26.
!<   C. hollow body, fully inverted (B inward, U outward): fixed, volume 26.
!<   D. hollow body, cavity wall wrong (B and U outward): U flipped, volume 26.
!< Signed distances: inside the solid negative, inside the cavity and outside positive.

program fossil_test_shell_orientation

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

logical :: ok(4) !< Per-case results.

ok(1) = check_disjoint()
ok(2) = check_hollow(invert_outer=.false., invert_inner=.true.,  label='B hollow, correct          ')
ok(3) = check_hollow(invert_outer=.true.,  invert_inner=.false., label='C hollow, fully inverted   ')
ok(4) = check_hollow(invert_outer=.false., invert_inner=.false., label='D hollow, cavity wall wrong')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check_disjoint() result(passed)
   !< Case A: two disjoint unit cubes, the second inverted.
   logical                  :: passed !< Check result.
   type(surface_stl_object) :: a, b   !< Cubes.
   real(R8P)                :: d1, d2 !< Signed distances at the two centres.

   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%translate(delta=5._R8P * ex_R8P)
   call b%reverse_normals
   call a%merge_solids(other=b)
   call a%analyze
   call a%sanitize
   d1 = a%distance(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), is_signed=.true., is_square_root=.true.)
   d2 = a%distance(point=vector_R8P(5.5_R8P, 0.5_R8P, 0.5_R8P), is_signed=.true., is_square_root=.true.)
   passed = (abs(a%get_volume() - 2._R8P) < 1.e-12_R8P) .and. (d1 < 0._R8P) .and. (d2 < 0._R8P)
   print '(A,F8.4,A,F8.4,A,F8.4,A,L1)', 'A disjoint, one inverted   : volume=', a%get_volume(), ' d(c1)=', d1, &
                                        ' d(c2)=', d2, ' ok=', passed
   endfunction check_disjoint

   function check_hollow(invert_outer, invert_inner, label) result(passed)
   !< Cases B-D: big cube [-1,2]^3 with the unit cube as a cavity.
   logical,      intent(in) :: invert_outer, invert_inner !< Initial orientation of each shell.
   character(*), intent(in) :: label                      !< Case label.
   logical                  :: passed                     !< Check result.
   type(surface_stl_object) :: outer, inner               !< Shells.
   real(R8P)                :: d_solid, d_cavity, d_out   !< Signed distances.

   call outer%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call outer%resize(factor=3._R8P * (ex_R8P + ey_R8P + ez_R8P))
   call outer%translate(delta=-1._R8P * (ex_R8P + ey_R8P + ez_R8P))
   call inner%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   if (invert_outer) call outer%reverse_normals
   if (invert_inner) call inner%reverse_normals
   call outer%merge_solids(other=inner)
   call outer%analyze
   call outer%sanitize
   d_solid  = outer%distance(point=vector_R8P(-0.5_R8P, -0.5_R8P, -0.5_R8P), is_signed=.true., is_square_root=.true.)
   d_cavity = outer%distance(point=vector_R8P( 0.5_R8P,  0.5_R8P,  0.5_R8P), is_signed=.true., is_square_root=.true.)
   d_out    = outer%distance(point=vector_R8P( 3.0_R8P,  0.5_R8P,  0.5_R8P), is_signed=.true., is_square_root=.true.)
   passed = (abs(outer%get_volume() - 26._R8P) < 1.e-10_R8P) .and. &
            (d_solid < 0._R8P) .and. (d_cavity > 0._R8P) .and. (d_out > 0._R8P)
   print '(A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,L1)', label, ': volume=', outer%get_volume(), ' d(solid)=', d_solid, &
                                                 ' d(cavity)=', d_cavity, ' d(outside)=', d_out, ' ok=', passed
   endfunction check_hollow
endprogram fossil_test_shell_orientation
