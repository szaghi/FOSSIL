!< FOSSIL test-side module: closed curved bodies for the boolean validation (issue #63).
module fossil_boolean_bodies
!< FOSSIL test-side module: closed curved bodies for the boolean validation (issue #63).
!<
!< `mc_body(kind, n, s)`: the zero level set of a signed field on n^3 samples of its box, by marching cubes, then
!< `sanitize(fill_holes=.true.)` (drops the slivers marching cubes leaves, closes the holes that leaves); the result
!< must be watertight with V > 0 (else error stop). Bodies, centred at the origin:
!< - BODY_SPHERE: radius 0.5;
!< - BODY_TORUS: radii 0.4 (tube centre) and 0.15 (tube), axis z;
!< - BODY_CYLINDER: radius 0.3, height 1, axis z, flat caps (the min of the two distances).
!< `move_body(s, axis, angle, t)`: a proper rigid motion (rotation about the origin, then translation).
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: BODY_SPHERE, BODY_TORUS, BODY_CYLINDER, mc_body, move_body

integer(I4P), parameter :: BODY_SPHERE = 1_I4P, BODY_TORUS = 2_I4P, BODY_CYLINDER = 3_I4P !< Body kinds.

contains
   subroutine mc_body(kind, n, s)
   !< A closed body of the given kind from n^3 samples.
   integer(I4P),             intent(in)    :: kind !< BODY_*.
   integer(I4P),             intent(in)    :: n    !< Samples per direction.
   type(surface_stl_object), intent(inout) :: s    !< Body.
   type(facet_object), allocatable         :: f(:) !< Facets.
   real(R8P)                               :: phi(n, n, n), x(3), h, ext, q !< Field.
   integer(I4P)                            :: i, j, k, status !< Counters.

   ext = 0.65_R8P ; h = 2._R8P * ext / (n - 1)
   do k = 1, n ; do j = 1, n ; do i = 1, n
      x = -ext + h * [i - 1, j - 1, k - 1]
      select case (kind)
      case (BODY_TORUS)
         q = sqrt(x(1) ** 2 + x(2) ** 2) - 0.4_R8P
         phi(i, j, k) = sqrt(q * q + x(3) ** 2) - 0.15_R8P
      case (BODY_CYLINDER)
         phi(i, j, k) = max(sqrt(x(1) ** 2 + x(2) ** 2) - 0.3_R8P, abs(x(3)) - 0.5_R8P)
      case default
         phi(i, j, k) = norm2(x) - 0.5_R8P
      endselect
   enddo ; enddo ; enddo
   call extract_isosurface(values=phi, bmin=vector_R8P(-ext, -ext, -ext), bmax=vector_R8P(ext, ext, ext), iso=0._R8P, &
                           surface=f, status=status)
   if (status /= MC_STATUS_OK) error stop 'mc_body: marching cubes failed'
   call s%adopt_facets(facets=f)
   call s%sanitize(fill_holes=.true.)
   if (.not. s%is_watertight() .or. s%get_volume() <= 0._R8P) error stop 'mc_body: not a closed body'
   endsubroutine mc_body

   subroutine move_body(s, axis, angle, t)
   !< Rotate s about the origin by angle around axis, then translate by t; rebuild.
   type(surface_stl_object), intent(inout) :: s       !< Body.
   real(R8P),                intent(in)    :: axis(3) !< Rotation axis.
   real(R8P),                intent(in)    :: angle   !< Rotation angle [rad].
   real(R8P),                intent(in)    :: t(3)    !< Translation.

   call s%rotate(axis=vector_R8P(axis(1), axis(2), axis(3)), angle=angle)
   call s%translate(x=t(1), y=t(2), z=t(3))
   call s%analyze
   endsubroutine move_body
endmodule fossil_boolean_bodies
