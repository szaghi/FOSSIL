!< FOSSIL, test the shared solid-angle core `solid_angle_oac` (issue #29, step 1).
!<
!< Checks:
!<   1. Agreement with the former vector_R8P-based Van Oosterom-Strackee formula (kept here as
!<      the reference) on 10^4 pseudo-random triangles and points: relative error <= 1e-13.
!<   2. `facet%solid_angle` delegates: identical to `solid_angle_oac` bit for bit.
!<   3. Sign: a counter-clockwise triangle seen from the side opposite its normal is positive.
!<   4. Degenerate triangle and point on the triangle plane (outside it): 0.
!<   5. Point just above the centre of a large triangle: close to 2*pi.
!<   6. Unit cube: sum over its 12 facets = 4*pi at an interior point, 0 at an exterior one.

program fossil_test_solid_angle

use fossil,             only : surface_stl_object, facet_object
use fossil_solid_angle, only : solid_angle_oac
use penf,               only : I4P, I8P, R8P
use vecfor,             only : vector_R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P)
logical                  :: ok(6)
type(vector_R8P)         :: p, a, b, c
type(facet_object)       :: f
type(surface_stl_object) :: cube
type(facet_object), pointer :: fp
real(R8P)                :: w, w_ref, err_max, s_in, s_out
integer(I4P)             :: k, i
integer(I8P)             :: seed

! 1-2. Random triangles vs the reference formula; facet delegation.
seed = 12345_I8P
err_max = 0._R8P
ok(2) = .true.
do k = 1, 10000
   p = rnd_vec() ; a = rnd_vec() ; b = rnd_vec() ; c = rnd_vec()
   w     = solid_angle_oac(p%x, p%y, p%z, a%x, a%y, a%z, b%x, b%y, b%z, c%x, c%y, c%z)
   w_ref = reference(p, a, b, c)
   err_max = max(err_max, abs(w - w_ref) / max(abs(w_ref), 1.e-3_R8P))
   f%vertex(1) = a ; f%vertex(2) = b ; f%vertex(3) = c
   if (f%solid_angle(point=p) /= w) ok(2) = .false.
enddo
ok(1) = err_max <= 1.e-13_R8P

! 3. Triangle in z = 0, counter-clockwise from +z (normal +z); seen from below it is positive.
w = solid_angle_oac(0._R8P, 0._R8P, -1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P)
ok(3) = w > 0._R8P .and. solid_angle_oac(0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, &
                                         1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P) < 0._R8P

! 4. Degenerate (collinear) triangle; point in the plane, outside the triangle. The triple product
!    is zero only up to rounding (1.4e-16 for the collinear case), so compare with a tolerance.
ok(4) = abs(solid_angle_oac(0.3_R8P, 0.2_R8P, 0.7_R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P, 1._R8P, &
                            2._R8P, 2._R8P, 2._R8P)) < 1.e-14_R8P .and. &
        abs(solid_angle_oac(5._R8P, 5._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, &
                            0._R8P, 1._R8P, 0._R8P)) < 1.e-14_R8P

! 5. Point just below the centre of a huge counter-clockwise triangle: about 2*pi.
w = solid_angle_oac(0._R8P, 0._R8P, -1.e-6_R8P, -1.e3_R8P, -1.e3_R8P, 0._R8P, 1.e3_R8P, -1.e3_R8P, 0._R8P, &
                    0._R8P, 1.e3_R8P, 0._R8P)
ok(5) = abs(w - 2._R8P * PI) < 1.e-6_R8P

! 6. Closed unit cube.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%sanitize
s_in = 0._R8P ; s_out = 0._R8P
do i = 1, cube%get_facets_number()
   fp => cube%facet_at(i)
   s_in  = s_in  + fp%solid_angle(point=vector_R8P(0.3_R8P, 0.6_R8P, 0.4_R8P))
   s_out = s_out + fp%solid_angle(point=vector_R8P(2.0_R8P, 0.6_R8P, 0.4_R8P))
enddo
ok(6) = abs(s_in - 4._R8P * PI) < 1.e-12_R8P .and. abs(s_out) < 1.e-12_R8P

print '(A,ES10.3,A,6L2)', 'max relative error vs reference: ', err_max, '  checks:', ok
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function reference(p, a, b, c) result(omega)
   !< The former facet%solid_angle formula, written with vector_R8P operations.
   type(vector_R8P), intent(in) :: p, a, b, c
   real(R8P)                    :: omega
   type(vector_R8P)             :: r1, r2, r3
   real(R8P)                    :: n1, n2, n3, num, den

   r1 = a - p ; n1 = r1%normL2()
   r2 = b - p ; n2 = r2%normL2()
   r3 = c - p ; n3 = r3%normL2()
   num = r1%dotproduct(rhs=r2%crossproduct(rhs=r3))
   den = n1 * n2 * n3 + (r1%dotproduct(rhs=r2)) * n3 + (r1%dotproduct(rhs=r3)) * n2 + (r2%dotproduct(rhs=r3)) * n1
   omega = 2._R8P * atan2(num, den)
   endfunction reference

   function rnd_vec() result(v)
   !< Deterministic pseudo-random vector in [-1, 1]^3 (64-bit LCG).
   type(vector_R8P) :: v
   v%x = rnd() ; v%y = rnd() ; v%z = rnd()
   endfunction rnd_vec

   function rnd() result(r)
   real(R8P) :: r
   seed = 6364136223846793005_I8P * seed + 1442695040888963407_I8P
   r = 2._R8P * real(ishft(seed, -11), R8P) / 9007199254740992._R8P - 1._R8P
   endfunction rnd
endprogram fossil_test_solid_angle
