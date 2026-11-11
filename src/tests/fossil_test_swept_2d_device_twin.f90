!< FOSSIL test: `swept_cells_2d_device` equals `swept_cells_2d` bit for bit (issue #100 step 4).
program fossil_test_swept_2d_device_twin
!< FOSSIL test: `swept_cells_2d_device` equals `swept_cells_2d` bit for bit (issue #100 step 4).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver bugs of the
!< device twin (the two end poses as the blocks of one multi-block call, the alpha taken out of it), not device codegen
!< (the curve bench does that on the GPU, `swept cells 2-D device vs serial mismatches`). Every output compared by its
!< bits (`swept_twin_mismatches`):
!< - exact translation: an ellipse on a stretched grid, a square moving by exactly one cell (its sides on grid lines at
!<   both poses: ties), a heptagon translating at a rotated pose;
!< - exact rotation: the ellipse about an off-centre point, the square about a grid node;
!< - Gauss-Legendre (`exact=.false.`, nq = 3) on a rotation, and explicit node poses;
!< - an empty curve, and the refusals (decreasing edges, nq = 0, a poses size other than nq): the same statuses;
!< - negative control: a pose one ulp off gives mismatches.
use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, pose_interpolate, gauss_legendre_01, swept_cells_2d_t, &
                   swept_cells_2d
use fossil_swept_compare, only : swept_twin_mismatches, swept_cells_2d_mismatches
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P)
type(curve_object)    :: ce, sq, hp, em
type(rigid_pose_2d_t) :: p0, p1, pr, pq(3)
type(swept_cells_2d_t) :: a, b
real(R8P)             :: x(25), y(21), xu(9), u(32), v(32), xd(9)
real(R8P), allocatable :: sn(:), wn(:)
integer(I4P)          :: k, m, total
logical               :: ok

do k = 1, 25
   x(k) = -0.05_R8P + 1.1_R8P * ((k - 1) / 24._R8P)**1.2_R8P
enddo
do k = 1, 21
   y(k) = -0.03_R8P + 1.06_R8P * (k - 1) / 20._R8P
enddo
do k = 1, 9
   xu(k) = (k - 1) / 8._R8P
enddo
do k = 1, 32
   u(k) = 0.5_R8P + 0.31_R8P * cos(2 * PI * (k - 1) / 32) ; v(k) = 0.48_R8P + 0.2_R8P * sin(2 * PI * (k - 1) / 32)
enddo
call ce%adopt_loops(u=u, v=v, loop_sizes=[32_I4P])
call sq%adopt_loops(u=[0.25_R8P, 0.625_R8P, 0.625_R8P, 0.25_R8P], v=[0.25_R8P, 0.25_R8P, 0.75_R8P, 0.75_R8P], &
                    loop_sizes=[4_I4P])
call hp%adopt_loops(u=[0.3_R8P, 0.55_R8P, 0.7_R8P, 0.62_R8P, 0.5_R8P, 0.38_R8P, 0.25_R8P], &
                    v=[0.3_R8P, 0.25_R8P, 0.4_R8P, 0.62_R8P, 0.45_R8P, 0.7_R8P, 0.5_R8P], loop_sizes=[7_I4P])
ok = .true. ; total = 0_I4P
p0 = pose_2d(0._R8P, [0._R8P, 0._R8P])
! exact translation
call check('ellipse, translation', swept_twin_mismatches(ce, p0, pose_2d(0._R8P, [0.071_R8P, -0.043_R8P]), x, y))
call check('square, one-cell translation (ties)', swept_twin_mismatches(sq, p0, pose_2d(0._R8P, [0.125_R8P, 0._R8P]), &
                                                                        xu, xu))
pr = pose_2d(0.4_R8P, [0.1_R8P, -0.2_R8P])
call check('heptagon, translation at a rotated pose', swept_twin_mismatches(hp, pr, pose_2d(0.4_R8P, [0.13_R8P, -0.17_R8P]), &
                                                                            x, y))
! exact rotation
pr = pose_2d(0.15_R8P, [0._R8P, 0._R8P]) ; pr%t = [0.6_R8P, 0.4_R8P] - matmul(pr%r, [0.6_R8P, 0.4_R8P])
call check('ellipse, rotation about an off-centre point', swept_twin_mismatches(ce, p0, pr, x, y))
pr = pose_2d(0.2_R8P, [0._R8P, 0._R8P]) ; pr%t = [0.5_R8P, 0.5_R8P] - matmul(pr%r, [0.5_R8P, 0.5_R8P])
call check('square, rotation about a grid node', swept_twin_mismatches(sq, p0, pr, xu, xu))
! Gauss-Legendre, explicit poses
call check('ellipse, Gauss nq = 3', swept_twin_mismatches(ce, p0, pr, x, y, nq=3_I4P, exact=.false.))
call gauss_legendre_01(3_I4P, sn, wn)
do k = 1, 3
   pq(k) = pose_interpolate(p0, pr, sn(k))
enddo
call check('ellipse, explicit node poses', swept_twin_mismatches(ce, p0, pr, x, y, nq=3_I4P, poses=pq))
! empty curve, refusals
call check('empty curve, translation', swept_twin_mismatches(em, p0, pose_2d(0._R8P, [0.1_R8P, 0._R8P]), xu, xu))
xd = xu ; xd(4) = xd(3)
call check('refusal: decreasing edges', swept_twin_mismatches(ce, p0, pr, xd, xu))
call check('refusal: nq = 0', swept_twin_mismatches(ce, p0, pr, x, y, nq=0_I4P, exact=.false.))
call check('refusal: poses size', swept_twin_mismatches(ce, p0, pr, x, y, nq=2_I4P, poses=pq))
! negative control
call swept_cells_2d(ce, p0, pose_2d(0._R8P, [0.071_R8P, -0.043_R8P]), x, y, a)
call swept_cells_2d(ce, p0, pose_2d(0._R8P, [nearest(0.071_R8P, 1._R8P), -0.043_R8P]), x, y, b)
m = swept_cells_2d_mismatches(a, b)
print '(A,I0)', 'negative control (one ulp): mismatches ', m
ok = ok .and. m > 0_I4P
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label, n)
   !< One case.
   character(*), intent(in) :: label !< Case.
   integer(I4P), intent(in) :: n     !< Mismatches.

   print '(A,I0)', label//': mismatches ', n
   total = total + n ; ok = ok .and. n == 0_I4P
   endsubroutine check
endprogram fossil_test_swept_2d_device_twin
