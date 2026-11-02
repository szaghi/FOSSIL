!< FOSSIL, signed solid angle of a triangle seen from a point (shared, device-callable core).
module fossil_solid_angle
!< Signed solid angle subtended by a triangle at a query point, Van Oosterom-Strackee:
!<
!<     Omega = 2 atan2( R1 . (R2 x R3),
!<                      |R1||R2||R3| + (R1.R2)|R3| + (R1.R3)|R2| + (R2.R3)|R1| )
!<
!< with Ri = vi - p. Positive when the triangle is seen counter-clockwise from p (p on the
!< side its normal points away from), so a closed outward surface sums to 4 pi at interior
!< points and 0 outside.
!<
!< Plain-real arguments, no derived types: the routine is `!$acc routine seq` and is used by
!< the fast-winding-number traversal on host and device alike (issue #29, roadmap #34), so the
!< two paths run the same arithmetic. It lives in its own module because adding a device
!< routine to `fossil_facet_object` triggers the nvfortran NVVM duplicate-type bug (see
!< `fossil_triangle_closest`). `facet%solid_angle` delegates here.
use penf, only : R8P
implicit none
private
public :: solid_angle_oac
contains
   pure function solid_angle_oac(px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz) result(omega)
   !< Signed solid angle of triangle (a, b, c) seen from point p.
   !$acc routine seq
   real(R8P), intent(in) :: px, py, pz       !< Query point.
   real(R8P), intent(in) :: ax, ay, az       !< Vertex 1.
   real(R8P), intent(in) :: bx, by, bz       !< Vertex 2.
   real(R8P), intent(in) :: cx, cy, cz       !< Vertex 3.
   real(R8P)             :: omega            !< Solid angle.
   real(R8P)             :: r1x, r1y, r1z    !< Vertex 1 minus point.
   real(R8P)             :: r2x, r2y, r2z    !< Vertex 2 minus point.
   real(R8P)             :: r3x, r3y, r3z    !< Vertex 3 minus point.
   real(R8P)             :: n1, n2, n3       !< Norms of r1, r2, r3.
   real(R8P)             :: numerator        !< r1 . (r2 x r3).
   real(R8P)             :: denominator      !< Van Oosterom-Strackee denominator.

   r1x = ax - px ; r1y = ay - py ; r1z = az - pz
   r2x = bx - px ; r2y = by - py ; r2z = bz - pz
   r3x = cx - px ; r3y = cy - py ; r3z = cz - pz
   n1 = sqrt(r1x * r1x + r1y * r1y + r1z * r1z)
   n2 = sqrt(r2x * r2x + r2y * r2y + r2z * r2z)
   n3 = sqrt(r3x * r3x + r3y * r3y + r3z * r3z)
   numerator = r1x * (r2y * r3z - r2z * r3y) + &
               r1y * (r2z * r3x - r2x * r3z) + &
               r1z * (r2x * r3y - r2y * r3x)
   denominator = n1 * n2 * n3 + (r1x * r2x + r1y * r2y + r1z * r2z) * n3 + &
                                (r1x * r3x + r1y * r3y + r1z * r3z) * n2 + &
                                (r2x * r3x + r2y * r3y + r2z * r3z) * n1
   omega = 2._R8P * atan2(numerator, denominator)
   endfunction solid_angle_oac
endmodule fossil_solid_angle
