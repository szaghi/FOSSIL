!< FOSSIL, local thickness and gaps: pure helpers (issue #62, roadmap #58 WP4).
module fossil_thickness
!< FOSSIL, local thickness and gaps: pure helpers (issue #62, roadmap #58 WP4).
!<
!< Plain reals, no surface type (the surface and curve modules wire the type-bound procedures):
!< - `cone_cosines` / `cone_directions`: the fixed ray set of the ray method, `m` unit directions in a cone of full
!<   aperture `cone_angle` [degrees] around a unit axis. Direction 1 is the axis; direction k has the polar cosine
!<   1 - (1 - cos(cone_angle/2)) (k - 1)/(m - 1) (equal-area steps up to the rim) and the azimuth (k - 1) times the
!<   golden angle. No random numbers: the same facet gives the same rays on every run.
!< - `lower_median`: the median of a small sample, the lower one for an even count.
!< - `ball_radius`: the radius of the ball tangent at p, centred along the unit direction a from p, that passes
!<   through q: r = |p - q|^2 / (2 a.(q - p)) (the shrinking-ball update).
use penf, only : I4P, R8P

implicit none
private
public :: THICKNESS_RAY, THICKNESS_MEDIAL, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS
public :: cone_cosines, cone_directions, lower_median, ball_radius

integer(I4P), parameter :: THICKNESS_RAY          = 1_I4P   !< Ray method (shape diameter).
integer(I4P), parameter :: THICKNESS_MEDIAL       = 2_I4P   !< Medial method (shrinking ball).
real(R8P),    parameter :: THICKNESS_DEFAULT_CONE = 60._R8P !< Default full cone aperture [degrees].
integer(I4P), parameter :: THICKNESS_DEFAULT_RAYS = 31_I4P  !< Default rays per facet.
real(R8P),    parameter :: PI = 3.141592653589793_R8P       !< Pi.
real(R8P),    parameter :: GOLDEN_ANGLE = PI * (3._R8P - 2.23606797749979_R8P) !< Golden angle [rad].

contains
   pure function cone_cosines(cone_angle, m) result(c)
   !< Polar cosines of the m directions of the ray set (direction 1 the axis, direction m on the rim).
   real(R8P),    intent(in) :: cone_angle !< Full aperture [degrees].
   integer(I4P), intent(in) :: m          !< Number of directions.
   real(R8P)                :: c(m)       !< Polar cosines.
   real(R8P)                :: ch         !< Cosine of the half aperture.
   integer(I4P)             :: k          !< Counter.

   ch = cos(0.5_R8P * cone_angle * PI / 180._R8P)
   c(1) = 1._R8P
   do k = 2, m
      c(k) = 1._R8P - (1._R8P - ch) * real(k - 1, R8P) / real(m - 1, R8P)
   enddo
   endfunction cone_cosines

   pure subroutine cone_directions(axis, cone_angle, m, dirs)
   !< The m unit directions of the ray set around the unit `axis`.
   real(R8P),    intent(in)  :: axis(3)    !< Unit axis.
   real(R8P),    intent(in)  :: cone_angle !< Full aperture [degrees].
   integer(I4P), intent(in)  :: m          !< Number of directions.
   real(R8P),    intent(out) :: dirs(3, m) !< Unit directions.
   real(R8P)                 :: c(m), s, e1(3), e2(3), phi !< Work.
   integer(I4P)              :: k          !< Counter.

   ! orthonormal basis (axis, e1, e2): e1 from the coordinate axis least parallel to `axis`
   if (abs(axis(1)) <= abs(axis(2)) .and. abs(axis(1)) <= abs(axis(3))) then
      e1 = [0._R8P, -axis(3), axis(2)]
   elseif (abs(axis(2)) <= abs(axis(3))) then
      e1 = [-axis(3), 0._R8P, axis(1)]
   else
      e1 = [-axis(2), axis(1), 0._R8P]
   endif
   e1 = e1 / sqrt(e1(1) ** 2 + e1(2) ** 2 + e1(3) ** 2)
   e2 = [axis(2) * e1(3) - axis(3) * e1(2), axis(3) * e1(1) - axis(1) * e1(3), axis(1) * e1(2) - axis(2) * e1(1)]
   c = cone_cosines(cone_angle, m)
   do k = 1, m
      s = sqrt(max(0._R8P, 1._R8P - c(k) * c(k)))
      phi = real(k - 1, R8P) * GOLDEN_ANGLE
      dirs(:, k) = c(k) * axis + s * (cos(phi) * e1 + sin(phi) * e2)
   enddo
   endsubroutine cone_directions

   pure function lower_median(x) result(med)
   !< Median of x (the lower one for an even count; 0 for an empty sample), by insertion sort of a copy.
   real(R8P), intent(in) :: x(:) !< Sample.
   real(R8P)             :: med  !< Median.
   real(R8P)             :: w(size(x)), key !< Work.
   integer(I4P)          :: n, i, j         !< Counters.

   n = size(x, kind=I4P)
   med = 0._R8P
   if (n == 0) return
   w = x
   do i = 2, n
      key = w(i) ; j = i - 1
      do while (j >= 1)
         if (w(j) <= key) exit
         w(j + 1) = w(j) ; j = j - 1
      enddo
      w(j + 1) = key
   enddo
   med = w((n + 1) / 2)
   endfunction lower_median

   pure function ball_radius(p, a, q) result(r)
   !< Radius of the ball tangent at p, centred at p + r a (a unit), through q: |p - q|^2 / (2 a.(q - p)); huge when q is
   !< not on the side of a (no finite ball reaches it).
   real(R8P), intent(in) :: p(3), a(3), q(3) !< Tangency point, direction to the centre, point on the sphere.
   real(R8P)             :: r                !< Radius.
   real(R8P)             :: w(3), den        !< Work.

   w = q - p
   den = 2._R8P * (a(1) * w(1) + a(2) * w(2) + a(3) * w(3))
   r = huge(1._R8P)
   if (den > 0._R8P) r = (w(1) * w(1) + w(2) * w(2) + w(3) * w(3)) / den
   endfunction ball_radius
endmodule fossil_thickness
