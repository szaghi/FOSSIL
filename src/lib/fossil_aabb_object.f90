!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition.

module fossil_aabb_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition.

use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none
private
public :: aabb_object

type :: aabb_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) class.
   type(vector_R8P)                :: bmin            !< Minimum point of AABB.
   type(vector_R8P)                :: bmax            !< Maximum point of AABB.
   integer(I4P)                    :: facets_number=0 !< Facets number.
   type(facet_object), allocatable :: facet(:)        !< Facets.
   contains
      ! public methods
      procedure, pass(self) :: closest_point    !< Return closest point on AABB from point reference.
      procedure, pass(self) :: compute_octant   !< Compute AABB octants.
      procedure, pass(self) :: destroy          !< Destroy AABB.
      procedure, pass(self) :: distance         !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: do_ray_intersect !< Return true if AABB is intersected by ray.
      procedure, pass(self) :: initialize       !< Initialize AABB.
      procedure, pass(self) :: vertex           !< Return AABB vertices.
      ! operators
      generic :: assignment(=) => aabb_assign_aabb      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_assign_aabb !< Operator `=`.
endtype aabb_object

contains
   ! public methods
   pure function closest_point(self, point) result(closest)
   !< Return closest point on (or in) AABB from point reference.
   class(aabb_object), intent(in) :: self    !< AABB.
   type(vector_R8P),   intent(in) :: point   !< Point reference.
   type(vector_R8P)               :: closest !< Closest point on (on in) aabb to point.

   closest = point
   closest%x = max(closest%x, self%bmin%x) ; closest%x = min(closest%x, self%bmax%x)
   closest%y = max(closest%y, self%bmin%y) ; closest%y = min(closest%y, self%bmax%y)
   closest%z = max(closest%z, self%bmin%z) ; closest%z = min(closest%z, self%bmax%z)
   endfunction closest_point

   pure subroutine compute_octant(self, octant)
   !< Return AABB octants.
   class(aabb_object), intent(in)  :: self      !< AABB.
   type(aabb_object),  intent(out) :: octant(8) !< AABB octants.
   type(vector_R8P)                :: vertex(8) !< AABB vertices.
   integer(I4P)                    :: o         !< Counter.

   vertex = self%vertex()
   octant(1)%bmin = self%bmin      ; octant(1)%bmax = 0.5_R8P * (self%bmin + self%bmax)
   octant(8)%bmin = octant(1)%bmax ; octant(8)%bmax = self%bmax
   do o=2, 7 ! loop over remaining octants
      octant(o)%bmin = 0.5_R8P * (self%bmin + vertex(o)) ; octant(o)%bmax = 0.5_R8P * (vertex(o) + self%bmax)
   enddo
   endsubroutine compute_octant

   elemental subroutine destroy(self)
   !< Destroy AABB.
   class(aabb_object), intent(inout) :: self  !< AABB.
   type(aabb_object)                 :: fresh !< Fresh instance of AABB box.

   self = fresh
   endsubroutine destroy

   pure function distance(self, point)
   !< Return the (square) distance from point to AABB.
   class(aabb_object), intent(in) :: self     !< AABB.
   type(vector_R8P),   intent(in) :: point    !< Point reference.
   real(R8P)                      :: distance !< Distance from point to AABB.

   distance = 0._R8P
   if (point%x < self%bmin%x) distance = distance + (self%bmin%x - point%x) * (self%bmin%x - point%x)
   if (point%y < self%bmin%y) distance = distance + (self%bmin%y - point%y) * (self%bmin%y - point%y)
   if (point%z < self%bmin%z) distance = distance + (self%bmin%z - point%z) * (self%bmin%z - point%z)
   if (point%x > self%bmax%x) distance = distance + (point%x - self%bmax%x) * (point%x - self%bmax%x)
   if (point%y > self%bmax%y) distance = distance + (point%y - self%bmax%y) * (point%y - self%bmax%y)
   if (point%z > self%bmax%z) distance = distance + (point%z - self%bmax%z) * (point%z - self%bmax%z)
   endfunction distance

   pure function do_ray_intersect(self, ray_origin, ray_direction) result(do_intersect)
   !< Return true if AABB is intersected by ray from origin and oriented as ray direction vector.
   class(aabb_object), intent(in) :: self          !< AABB box.
   type(vector_R8P),   intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),   intent(in) :: ray_direction !< Ray direction.
   logical                        :: do_intersect  !< Test result.
   logical                        :: must_return   !< Flag to check when to return from procedure.
   real(R8P)                      :: tmin, tmax    !< Minimum maximum ray intersections with box slabs.

   do_intersect = .false.
   must_return = .false.
   tmin = 0._R8P
   tmax = MaxR8P
   call check_slab(aabb_min=self%bmin%x, aabb_max=self%bmax%x, &
                   o=ray_origin%x, d=ray_direction%x, must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   call check_slab(aabb_min=self%bmin%y, aabb_max=self%bmax%y, &
                   o=ray_origin%y, d=ray_direction%y, must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   call check_slab(aabb_min=self%bmin%z, aabb_max=self%bmax%z, &
                   o=ray_origin%z, d=ray_direction%z, must_return=must_return, tmin=tmin, tmax=tmax)
   if (must_return) return
   ! ray intersects all 3 slabs
   do_intersect = .true.
   contains
      pure subroutine check_slab(aabb_min, aabb_max, o, d, must_return, tmin, tmax)
      !< Perform ray intersection check in a direction-split fashion over slabs.
      real(R8P), intent(in)    :: aabb_min     !< Box minimum bound in the current direction.
      real(R8P), intent(in)    :: aabb_max     !< Box maximum bound in the current direction.
      real(R8P), intent(in)    :: o            !< Ray origin in the current direction.
      real(R8P), intent(in)    :: d            !< Ray slope in the current direction.
      logical,   intent(inout) :: must_return  !< Flag to check when to return from procedure.
      real(R8P), intent(inout) :: tmin, tmax   !< Minimum maximum ray intersections with box slabs.
      real(R8P)                :: ood, t1, t2  !< Intersection coefficients.
      real(R8P)                :: tmp          !< Temporary buffer.
      real(R8P), parameter     :: eps=1e-6_R8P !< Small epsilon to control round off errors.

      if ((d) < eps) then
         ! ray is parallel to slab, no hit if origin not within slab
         if ((o < aabb_min).or.(o > aabb_max)) then
            must_return = .true.
            return
         endif
      else
         ! compute intersection t value of ray with near and far plane of slab
         ood = 1._R8P / d
         t1 = (aabb_min - o) * ood
         t2 = (aabb_max - o) * ood
         ! make t1 be intersection with near plane, t2 with far plane
         if (t1 > t2) then
            tmp = t1
            t1 = t2
            t2 = tmp
         endif
         ! compute the intersection of slab intersection intervals
         if (t1 > tmin) tmin = t1
         if (t2 > tmax) tmax = t2
         ! exit with no collision as soon as slab intersection becomes empty
         if (tmin > tmax) then
            must_return = .true.
            return
         endif
      endif
      endsubroutine check_slab
   endfunction do_ray_intersect

   pure subroutine initialize(self, facet, bmin, bmax)
   !< Initialize AABB.
   class(aabb_object), intent(inout)        :: self     !< AABB.
   type(facet_object), intent(in), optional :: facet(:) !< Facets list.
   type(vector_R8P),   intent(in), optional :: bmin     !< Minimum point of AABB.
   type(vector_R8P),   intent(in), optional :: bmax     !< Maximum point of AABB.

   call self%destroy
   if (present(facet)) then
      self%bmin%x = minval(facet(:)%bb(1)%x)
      self%bmin%y = minval(facet(:)%bb(1)%y)
      self%bmin%z = minval(facet(:)%bb(1)%z)
      self%bmax%x = maxval(facet(:)%bb(2)%x)
      self%bmax%y = maxval(facet(:)%bb(2)%y)
      self%bmax%z = maxval(facet(:)%bb(2)%z)
   elseif (present(bmin).and.present(bmax)) then
      self%bmin = bmin
      self%bmax = bmax
   endif
   endsubroutine initialize

   pure function vertex(self)
   !< Return AABB vertices.
   class(aabb_object), intent(in) :: self      !< AABB.
   type(vector_R8P)               :: vertex(8) !< AABB vertices.

   vertex(1) = self%bmin
   vertex(2) = self%bmax%x * ex_R8P + self%bmin%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(3) = self%bmin%x * ex_R8P + self%bmax%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(4) = self%bmax%x * ex_R8P + self%bmax%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(5) = self%bmin%x * ex_R8P + self%bmin%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(6) = self%bmax%x * ex_R8P + self%bmin%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(7) = self%bmin%x * ex_R8P + self%bmax%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(8) = self%bmax
   endfunction vertex

   ! operators
   ! =
   pure subroutine aabb_assign_aabb(lhs, rhs)
   !< Operator `=`.
   class(aabb_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_object),  intent(in)    :: rhs !< Right hand side.

   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
   lhs%facets_number = rhs%facets_number
   if (allocated(lhs%facet)) deallocate(lhs%facet)
   if (allocated(rhs%facet)) allocate(lhs%facet(1:lhs%facets_number), source=rhs%facet)
   endsubroutine aabb_assign_aabb
endmodule fossil_aabb_object
