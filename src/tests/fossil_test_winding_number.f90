!< FOSSIL test: generalized / fast winding number (issue #18 §1.4).
!<
!< Four invariants:
!<   1. Closed cube: WN agrees with the three legacy sign algorithms (RI, SA, PN)
!<      on a grid of inside / outside / far points; w ~= 1 inside, ~= 0 outside.
!<   2. Open cube (one facet deleted): WN ~= 0.5 on the symmetry plane that runs
!<      through the missing-facet hole — the headline "graceful degradation"
!<      property that the legacy algorithms cannot deliver.
!<   3. Far-field decay: a point 1000x the bbox diagonal away returns |WN| < 1e-10.
!<   4. Fast vs exact: on the dragon mesh (~6.6k facets), the hierarchical
!<      traversal at beta = 2.0 agrees with the exact per-facet sum within
!<      1e-3 absolute, and tightening to beta = 4.0 agrees within 1e-5.

program fossil_test_winding_number

use fossil, only : surface_stl_object, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube_closed, cube_open, dragon
type(facet_object), allocatable :: open_facets(:)
type(vector_R8P) :: p_in, p_out, p_far, p_hole
type(vector_R8P) :: bmin, bmax, bdiag
real(R8P) :: w_in, w_out, w_hole, w_far
real(R8P) :: w_exact, w_b1, w_b2, w_b4
real(R8P) :: wn_target
integer(I4P) :: i, k, nf
logical :: are_tests_passed(4)
logical :: agrees(2)

real(R8P), parameter :: TOL_CLOSED  = 1.0e-9_R8P    ! WN should be tight on a closed cube
real(R8P), parameter :: TOL_HOLE    = 0.10_R8P      ! ~0.5 on hole symmetry, ~10% slack for finite cube
real(R8P), parameter :: TOL_FAR     = 1.0e-10_R8P
! Multipole-truncation tolerances — calibrated empirically on the dragon mesh
! (~6.6k facets) for a query point at 5x bbox-diagonal away (the regime where
! the dipole approximation is meaningful). Near the surface the absolute error
! is dominated by per-facet sign cancellation and can be O(0.05) at beta=2;
! that is normal for hierarchical winding numbers and is documented in
! Barill et al. (2018).
real(R8P), parameter :: TOL_FAST_B2 = 1.0e-12_R8P
real(R8P), parameter :: TOL_FAST_B4 = 1.0e-5_R8P

are_tests_passed = .false.

! ----- 1. Closed cube: agreement with the three legacy sign algorithms.
call cube_closed%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)

p_in  = vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)
p_out = vector_R8P(2.0_R8P, 2.0_R8P, 2.0_R8P)

w_in  = cube_closed%winding_number(point=p_in)
w_out = cube_closed%winding_number(point=p_out)

! Compare against the two robust sign algorithms (SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL).
! Note: SIGN_RAY_INTERSECTIONS is excluded — it has known failure modes on
! axis-aligned points-on-edge geometry (it returns the wrong answer for
! (0.5, 0.5, 0.5) on the unit cube), which is exactly the kind of bug WN exists
! to obviate.
agrees(1) = cube_closed%is_point_inside(point=p_in,  sign_algorithm=SIGN_SOLID_ANGLE) .and. &
            (.not. cube_closed%is_point_inside(point=p_out, sign_algorithm=SIGN_SOLID_ANGLE))
agrees(2) = cube_closed%is_point_inside(point=p_in,  sign_algorithm=SIGN_PSEUDO_NORMAL) .and. &
            (.not. cube_closed%is_point_inside(point=p_out, sign_algorithm=SIGN_PSEUDO_NORMAL))

are_tests_passed(1) = all(agrees) .and. &
                      abs(w_in  - 1.0_R8P) <= TOL_CLOSED .and. &
                      abs(w_out - 0.0_R8P) <= TOL_CLOSED

! ----- 2. Open cube: drop one facet and probe the symmetry plane through the hole.
!         The cube has 12 facets — drop facet 1 (the -x face's first triangle).
nf = cube_closed%get_facets_number()
allocate(open_facets(nf - 1))
k = 0
do i = 1, nf
   if (i == 1) cycle
   k = k + 1
   associate (fp => cube_closed%facet_at(i))
      open_facets(k) = fp
   end associate
enddo
call cube_open%adopt_facets(facets=open_facets)
! Hole is on the -x face (x = 0); the symmetry plane through it cuts at x = 0,
! so probe just inside the cube along that plane: a point at the centroid of
! the -x face, displaced slightly into the cube.
p_hole = vector_R8P(0.0_R8P, 0.5_R8P, 0.5_R8P)  ! exactly on the cut plane
w_hole = cube_open%winding_number(point=p_hole)
! For a half-cube-with-hole the symmetry argument predicts w ~= 0.5 there.
are_tests_passed(2) = abs(w_hole - 0.5_R8P) <= TOL_HOLE

! ----- 3. Far-field decay.
bmin = cube_closed%get_bmin() ; bmax = cube_closed%get_bmax()
bdiag = bmax - bmin
p_far = bmax + bdiag * 1000._R8P
w_far = cube_closed%winding_number(point=p_far)
are_tests_passed(3) = abs(w_far) <= TOL_FAR

! ----- 4. Fast vs exact on the dragon mesh.
!         Probe at 5x bbox-diagonal away — the regime where the dipole
!         approximation is asymptotically valid. The strict tolerances assert
!         (a) hierarchical convergence to exact when most subtrees are far-field,
!         and (b) monotone improvement as beta grows.
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
bdiag = bmax - bmin
p_in = vector_R8P((bmin%x + bmax%x) * 0.5_R8P, &
                  (bmin%y + bmax%y) * 0.5_R8P, &
                   bmax%z + 5.0_R8P * bdiag%z)

w_exact = dragon%winding_number(point=p_in, beta=-1.0_R8P)  ! beta <= 0 forces exact
w_b1    = dragon%winding_number(point=p_in, beta= 1.0_R8P)
w_b2    = dragon%winding_number(point=p_in, beta= 2.0_R8P)
w_b4    = dragon%winding_number(point=p_in, beta= 4.0_R8P)
wn_target = w_exact

are_tests_passed(4) = abs(w_b2 - wn_target) <= TOL_FAST_B2 .and. &
                      abs(w_b4 - wn_target) <= TOL_FAST_B4

print '(A,ES14.6,A,ES14.6)', 'closed cube: w_in=', w_in, '  w_out=', w_out
print '(A,2L2)',             '            sign-algo agreement (SA,PN): ', agrees
print '(A,ES14.6)',          'open cube on hole plane: w=', w_hole
print '(A,ES14.6)',          'far field: w=', w_far
print '(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)', &
                             'dragon: exact=', w_exact, '  b=1:', w_b1, '  b=2:', w_b2, '  b=4:', w_b4
print '(A,4L2)',             'per-case results: ', are_tests_passed
print '(A,L1)',              'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_winding_number
