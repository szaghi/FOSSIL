!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition.

module fossil_aabb_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition.

use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : FR8P, I4P, R8P, MaxR8P, str
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
      procedure, pass(self) :: add_facets                  !< Add facets to AABB.
      procedure, pass(self) :: closest_point               !< Return closest point on AABB from point reference.
      procedure, pass(self) :: compute_octants             !< Compute AABB octants.
      procedure, pass(self) :: destroy                     !< Destroy AABB.
      procedure, pass(self) :: distance                    !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: distance_from_facets        !< Return the (square) distance from point to AABB's facets.
      procedure, pass(self) :: do_ray_intersect            !< Return true if AABB is intersected by ray.
      procedure, pass(self) :: has_facets                  !< Return true if AABB has facets.
      procedure, pass(self) :: initialize                  !< Initialize AABB.
      procedure, pass(self) :: is_inside                   !< Return the true if point is inside ABB.
      procedure, pass(self) :: ray_intersections_number    !< Return ray intersections number.
      procedure, pass(self) :: save_geometry_tecplot_ascii !< Save AABB geometry into Tecplot ascii file.
      procedure, pass(self) :: save_facets_into_file_stl   !< Save facets into file STL.
      procedure, pass(self) :: update_extents              !< Update AABB bounding box extents.
      procedure, pass(self) :: vertex                      !< Return AABB vertices.
      ! operators
      generic :: assignment(=) => aabb_assign_aabb      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_assign_aabb !< Operator `=`.
endtype aabb_object

contains
   ! public methods
   subroutine add_facets(self, facet)
   !< Add facets to AABB.
   !<
   !< @note Previously stored facets are lost.
   !<
   !< @note Facets added to AABB are removed to facets list that is also returned.
   class(aabb_object),              intent(inout) :: self              !< AABB.
   type(facet_object), allocatable, intent(inout) :: facet(:)          !< Facets list.
   integer(I4P)                                   :: scratch_unit_add  !< Scratch unit file for added facets.
   integer(I4P)                                   :: scratch_unit_rem  !< Scratch unit file for remaining facets.
   integer(I4P)                                   :: rem_facets_number !< Remaining facets number.
   integer(I4P)                                   :: f                 !< Counter.

   self%facets_number = 0
   rem_facets_number = 0
   if (allocated(self%facet)) deallocate(self%facet)
   open(newunit=scratch_unit_add, status='scratch', access='stream', form='unformatted')
   open(newunit=scratch_unit_rem, status='scratch', access='stream', form='unformatted')
   do f=1, size(facet, dim=1)
      if (self%is_inside(point=facet(f)%vertex_1).and.&
          self%is_inside(point=facet(f)%vertex_2).and.&
          self%is_inside(point=facet(f)%vertex_3)) then
         self%facets_number = self%facets_number + 1
         call facet(f)%save_into_file_binary(file_unit=scratch_unit_add)
      else
         rem_facets_number = rem_facets_number + 1
         call facet(f)%save_into_file_binary(file_unit=scratch_unit_rem)
      endif
   enddo
   if (self%facets_number > 0) then
      allocate(self%facet(1:self%facets_number))
      rewind(unit=scratch_unit_add)
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_binary(file_unit=scratch_unit_add)
         call self%facet(f)%compute_metrix
      enddo
   endif
   close(unit=scratch_unit_add)
   deallocate(facet)
   if (rem_facets_number > 0) then
      allocate(facet(1:rem_facets_number))
      rewind(unit=scratch_unit_rem)
      do f=1, rem_facets_number
         call facet(f)%load_from_file_binary(file_unit=scratch_unit_rem)
         call facet(f)%compute_metrix
      enddo
   endif
   close(unit=scratch_unit_rem)
   endsubroutine add_facets

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

   pure subroutine compute_octants(self, octant)
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
   endsubroutine compute_octants

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
   real(R8P)                      :: dx, dy, dz !< Distance components.

   dx = max(self%bmin%x - point%x, 0._R8P, point%x - self%bmax%x)
   dy = max(self%bmin%y - point%y, 0._R8P, point%y - self%bmax%y)
   dz = max(self%bmin%z - point%z, 0._R8P, point%z - self%bmax%z)
   distance = dx * dx + dy * dy + dz * dz
   endfunction distance

   pure function distance_from_facets(self, point) result(distance)
   !< Return the (square) distance from point to AABB's facets.
   class(aabb_object), intent(in) :: self      !< AABB.
   type(vector_R8P),   intent(in) :: point     !< Point reference.
   real(R8P)                      :: distance  !< Distance from point to AABB's facets.
   real(R8P)                      :: distance_ !< Distance from point to AABB's facets, local variable.
   integer(I4P)                   :: f         !< Counter.

   distance = MaxR8P
   if (self%facets_number > 0) then
      do f=1, self%facets_number
         distance_ = self%facet(f)%distance(point=point)
         if (abs(distance_) <= abs(distance)) distance = distance_
      enddo
   endif
   endfunction distance_from_facets

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

   pure function has_facets(self)
   !< Return true if AABB has facets.
   class(aabb_object), intent(in) :: self       !< AABB box.
   logical                        :: has_facets !< Check result.

   has_facets = self%facets_number > 0
   endfunction has_facets

   pure subroutine initialize(self, facet, bmin, bmax)
   !< Initialize AABB.
   class(aabb_object), intent(inout)        :: self     !< AABB.
   type(facet_object), intent(in), optional :: facet(:) !< Facets list.
   type(vector_R8P),   intent(in), optional :: bmin     !< Minimum point of AABB.
   type(vector_R8P),   intent(in), optional :: bmax     !< Maximum point of AABB.

   call self%destroy
   if (present(facet)) then
      call compute_bb_from_facets(facet=facet, bmin=self%bmin, bmax=self%bmax)
   elseif (present(bmin).and.present(bmax)) then
      self%bmin = bmin
      self%bmax = bmax
   endif
   endsubroutine initialize

   pure function is_inside(self, point)
   !< Return the true if point is inside ABB.
   class(aabb_object), intent(in) :: self      !< AABB.
   type(vector_R8P),   intent(in) :: point     !< Point reference.
   logical                        :: is_inside !< Check result.

   is_inside = ((point%x >= self%bmin%x.and.point%x <= self%bmax%x).and.&
                (point%y >= self%bmin%y.and.point%y <= self%bmax%y).and.&
                (point%z >= self%bmin%z.and.point%z <= self%bmax%z))
   endfunction is_inside

   pure function ray_intersections_number(self, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number.
   class(aabb_object), intent(in) :: self                 !< AABB.
   type(vector_R8P),   intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),   intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                   :: intersections_number !< Intersection number.
   integer(I4P)                   :: f                    !< Counter.

   intersections_number = 0
   if (self%facets_number > 0) then
      do f=1, self%facets_number
         if (self%facet(f)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) &
            intersections_number = intersections_number + 1
      enddo
   endif
   endfunction ray_intersections_number

   subroutine  save_geometry_tecplot_ascii(self, file_unit, aabb_name)
   !< Save AABB geometry into Tecplot ascii file.
   class(aabb_object), intent(in)           :: self       !< AABB.
   integer(I4P),       intent(in)           :: file_unit  !< File unit.
   character(*),       intent(in), optional :: aabb_name  !< Name of AABB.
   character(len=:), allocatable            :: aabb_name_ !< Name of AABB, local variable.
   type(vector_R8P)                         :: vertex(8)  !< AABB vertices.
   integer(I4P)                             :: v          !< Counter.

   aabb_name_ = 'AABB' ; if (present(aabb_name)) aabb_name_ = aabb_name
   write(file_unit, '(A)') 'ZONE T="'//aabb_name//'", I=2, J=2, K=2'
   vertex = self%vertex()
   do v=1, 8
      write(file_unit, '(3('//FR8P//',1X))') vertex(v)%x, vertex(v)%y, vertex(v)%z
   enddo
   endsubroutine  save_geometry_tecplot_ascii

   subroutine save_facets_into_file_stl(self, file_name, is_ascii)
   !< Save facets into file STL.
   class(aabb_object), intent(in) :: self      !< AABB.
   character(*),       intent(in) :: file_name !< File name.
   logical,            intent(in) :: is_ascii  !< Sentinel to check if file is ASCII.
   integer(I4P)                   :: file_unit !< File unit.
   integer(I4P)                   :: f         !< Counter.

   if (self%facets_number > 0) then
      call open_file
      if (is_ascii) then
         do f=1, self%facets_number
            call self%facet(f)%save_into_file_ascii(file_unit=file_unit)
         enddo
      else
         do f=1, self%facets_number
            call self%facet(f)%save_into_file_binary(file_unit=file_unit)
         enddo
      endif
      call close_file
   endif
   contains
      subroutine open_file()
      !< Open STL file.

      if (is_ascii) then
         open(newunit=file_unit, file=trim(adjustl(file_name)),                  form='formatted')
         write(file_unit, '(A)') 'solid '//trim(adjustl(file_name))
      else
         open(newunit=file_unit, file=trim(adjustl(file_name)), access='stream', form='unformatted')
         write(file_unit) repeat('a', 80)
         write(file_unit) self%facets_number
      endif
      endsubroutine open_file

      subroutine close_file()
      !< Close STL file.

      if (is_ascii) write(file_unit, '(A)') 'endsolid '//trim(adjustl(file_name))
      close(unit=file_unit)
      endsubroutine close_file
   endsubroutine save_facets_into_file_stl

   pure subroutine update_extents(self)
   !< Update AABB bounding box extents.
   class(aabb_object), intent(inout) :: self !< AABB.

   if (self%facets_number > 0) call compute_bb_from_facets(facet=self%facet, bmin=self%bmin, bmax=self%bmax)
   endsubroutine update_extents

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

   ! non TBP
   pure subroutine compute_bb_from_facets(facet, bmin, bmax)
   !< Compute AABB extents (minimum and maximum bounding box) from facets list.
   !<
   !< @note Facets' metrix must be already computed.
   type(facet_object), intent(in)    :: facet(:) !< Facets list.
   type(vector_R8P),   intent(inout) :: bmin     !< Minimum point of AABB.
   type(vector_R8P),   intent(inout) :: bmax     !< Maximum point of AABB.
   real(R8P)                         :: eps(3) !< Small epsilon.

   eps(1) = (maxval(facet(:)%bb(2)%x) - minval(facet(:)%bb(1)%x)) / 100._R8P
   eps(2) = (maxval(facet(:)%bb(2)%y) - minval(facet(:)%bb(1)%y)) / 100._R8P
   eps(3) = (maxval(facet(:)%bb(2)%z) - minval(facet(:)%bb(1)%z)) / 100._R8P
   bmin%x = minval(facet(:)%bb(1)%x) - eps(1)
   bmin%y = minval(facet(:)%bb(1)%y) - eps(2)
   bmin%z = minval(facet(:)%bb(1)%z) - eps(3)
   bmax%x = maxval(facet(:)%bb(2)%x) + eps(1)
   bmax%y = maxval(facet(:)%bb(2)%y) + eps(2)
   bmax%z = maxval(facet(:)%bb(2)%z) + eps(3)
   endsubroutine compute_bb_from_facets
endmodule fossil_aabb_object
