!< FOSSIL, utils library.

module fossil_utils
!< FOSSIL, utils library.

use penf, only : I4P, R4P, R8P, ZeroR4P
use vecfor, only : vector_R8P

implicit none
private
public :: EPS
public :: FRLEN
public :: PI
public :: is_inside_bb
public :: triangle_overlaps_aabb

real(R8P), protected    :: EPS=real(ZeroR4P, R8P)     !< Small EPSILON to avoid rund off errors.
integer(I4P), parameter :: FRLEN=80                   !< Maximum length of fossil STL label string.
real(R8P),    parameter :: PI = 4._R8P * atan(1._R8P) !< Pi greek.

contains
   pure function is_inside_bb(bmin, bmax, point)
   !< Return the true if point is inside AABB.
   type(vector_R8P), intent(in) :: bmin, bmax   !< Bounding box extents.
   type(vector_R8P), intent(in) :: point        !< Point reference.
   logical                      :: is_inside_bb !< Check result.

   is_inside_bb = ((point%x >= bmin%x.and.point%x <= bmax%x).and.&
                   (point%y >= bmin%y.and.point%y <= bmax%y).and.&
                   (point%z >= bmin%z.and.point%z <= bmax%z))
   endfunction is_inside_bb

   pure function triangle_overlaps_aabb(bmin, bmax, v1, v2, v3) result(overlap)
   !< Conservative triangle-AABB overlap test (Akenine-Möller 2001, 13-axis SAT).
   !<
   !< Returns .true. iff the closed triangle (v1,v2,v3) and the closed AABB share
   !< at least one point. Never rejects a truly overlapping facet; touching contact
   !< is reported as overlap to match the closed-set semantics of `is_inside_bb`.
   !<
   !< Reference: T. Akenine-Möller, "Fast 3D Triangle-Box Overlap Testing",
   !< Journal of Graphics Tools, 6(1), 2001.
   type(vector_R8P), intent(in) :: bmin, bmax    !< AABB extents.
   type(vector_R8P), intent(in) :: v1, v2, v3    !< Triangle vertices.
   logical                      :: overlap       !< Overlap result.
   type(vector_R8P)             :: c, h          !< Box centre and half-extents.
   type(vector_R8P)             :: t1, t2, t3    !< Triangle recentred on box.
   type(vector_R8P)             :: e1, e2, e3    !< Triangle edge vectors.
   type(vector_R8P)             :: n             !< Triangle normal.
   real(R8P)                    :: r             !< Projected box radius.

   c = 0.5_R8P * (bmin + bmax)
   h = 0.5_R8P * (bmax - bmin)
   t1 = v1 - c
   t2 = v2 - c
   t3 = v3 - c

   overlap = .false.

   ! Category 1: 3 AABB face normals -- cheapest, rejects the majority of facets.
   if (min(t1%x, t2%x, t3%x) >  h%x .or. max(t1%x, t2%x, t3%x) < -h%x) return
   if (min(t1%y, t2%y, t3%y) >  h%y .or. max(t1%y, t2%y, t3%y) < -h%y) return
   if (min(t1%z, t2%z, t3%z) >  h%z .or. max(t1%z, t2%z, t3%z) < -h%z) return

   ! Category 2: triangle plane vs box.
   e1 = t2 - t1
   e2 = t3 - t2
   n  = e1.cross.e2
   r  = h%x*abs(n%x) + h%y*abs(n%y) + h%z*abs(n%z)
   if (abs(n.dot.t1) > r) return

   ! Category 3: 9 edge-edge cross-product axes (AABB edge_i x triangle edge_j).
   e3 = t1 - t3
   if (sat_separated(0._R8P, -e1%z,  e1%y, t1, t2, t3, h)) return
   if (sat_separated( e1%z, 0._R8P, -e1%x, t1, t2, t3, h)) return
   if (sat_separated(-e1%y,  e1%x, 0._R8P, t1, t2, t3, h)) return
   if (sat_separated(0._R8P, -e2%z,  e2%y, t1, t2, t3, h)) return
   if (sat_separated( e2%z, 0._R8P, -e2%x, t1, t2, t3, h)) return
   if (sat_separated(-e2%y,  e2%x, 0._R8P, t1, t2, t3, h)) return
   if (sat_separated(0._R8P, -e3%z,  e3%y, t1, t2, t3, h)) return
   if (sat_separated( e3%z, 0._R8P, -e3%x, t1, t2, t3, h)) return
   if (sat_separated(-e3%y,  e3%x, 0._R8P, t1, t2, t3, h)) return

   overlap = .true.
   endfunction triangle_overlaps_aabb

   pure function sat_separated(ax, ay, az, t1, t2, t3, h) result(separated)
   !< One SAT axis test. Returns .true. iff axis (ax,ay,az) separates triangle from box.
   real(R8P),        intent(in) :: ax, ay, az    !< Test axis components.
   type(vector_R8P), intent(in) :: t1, t2, t3    !< Triangle vertices (box-centred).
   type(vector_R8P), intent(in) :: h             !< Box half-extents.
   logical                      :: separated     !< Result: .true. iff axis is a separator.
   real(R8P)                    :: p1, p2, p3    !< Triangle vertex projections on axis.
   real(R8P)                    :: pmin, pmax, r !< Projection bounds and box radius.

   p1 = ax*t1%x + ay*t1%y + az*t1%z
   p2 = ax*t2%x + ay*t2%y + az*t2%z
   p3 = ax*t3%x + ay*t3%y + az*t3%z
   pmin = min(p1, p2, p3)
   pmax = max(p1, p2, p3)
   r = h%x*abs(ax) + h%y*abs(ay) + h%z*abs(az)
   separated = (pmin > r .or. pmax < -r)
   endfunction sat_separated
endmodule fossil_utils
