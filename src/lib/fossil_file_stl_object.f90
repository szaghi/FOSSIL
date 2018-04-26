!< FOSSIL,  STL file class definition.

module fossil_file_stl_object
!< FOSSIL,  STL file class definition.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none
private
public :: file_stl_object

real(R8P), parameter :: PI = 4._R8P * atan(1._R8P)

type :: file_stl_object
   !< FOSSIL STL file class.
   character(len=:), allocatable   :: file_name       !< File name
   integer(I4P)                    :: file_unit=0     !< File unit.
   character(FRLEN)                :: header          !< File header.
   integer(I4P)                    :: facets_number=0 !< Facets number.
   type(facet_object), allocatable :: facet(:)        !< Facets.
   type(aabb_tree_object)          :: aabb            !< AABB tree.
   logical                         :: is_ascii=.true. !< Sentinel to check if file is ASCII.
   logical                         :: is_open=.false. !< Sentinel to check if file is open.
   contains
      ! public methods
      procedure, pass(self) :: close_file                    !< Close file.
      procedure, pass(self) :: compute_metrix                !< Compute facets metrix.
      procedure, pass(self) :: create_aabb_tree              !< Create the AABB tree.
      procedure, pass(self) :: destroy                       !< Destroy file.
      procedure, pass(self) :: distance                      !< Compute the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: initialize                    !< Initialize file.
      procedure, pass(self) :: is_point_inside_polyhedron_ri !< Determinate is a point is inside or not STL facets by ray intersect.
      procedure, pass(self) :: is_point_inside_polyhedron_sa !< Determinate is a point is inside or not STL facets by solid angle.
      procedure, pass(self) :: load_from_file                !< Load from file.
      procedure, pass(self) :: open_file                     !< Open file, once initialized.
      procedure, pass(self) :: sanitize_normals              !< Sanitize normals, make normals consistent with vertices.
      procedure, pass(self) :: save_into_file                !< Save into file.
      ! operators
      generic :: assignment(=) => file_stl_assign_file_stl       !< Overload `=`.
      procedure, pass(lhs),  private :: file_stl_assign_file_stl !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: allocate_facets              !< Allocate facets.
      procedure, pass(self), private :: load_facets_number_from_file !< Load facets number from file.
      procedure, pass(self), private :: load_header_from_file        !< Load header from file.
      procedure, pass(self), private :: save_header_into_file        !< Save header into file.
endtype file_stl_object

contains
   ! public methods
   subroutine close_file(self)
   !< Close file.
   class(file_stl_object), intent(inout) :: self       !< File STL.
   logical                               :: file_exist !< Sentinel to check if file exist.

   if (self%is_open) then
      close(unit=self%file_unit)
      self%file_unit = 0
      self%is_open = .false.
   endif
   endsubroutine close_file

   pure subroutine compute_metrix(self)
   !< Compute facets metrix.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) then
      call self%facet%compute_metrix
   endif
   endsubroutine compute_metrix

   subroutine create_aabb_tree(self, refinement_levels)
   !< Create AABB tree.
   !<
   !< @note Facets metrix must be already computed.
   class(file_stl_object), intent(inout)        :: self               !< File STL.
   integer(I4P),           intent(in), optional :: refinement_levels  !< Total number of refinement levels used.
   integer(I4P)                                 :: refinement_levels_ !< Total number of refinement levels used, local variable.

   refinement_levels_ = 2 ; if (present(refinement_levels)) refinement_levels_ = refinement_levels
   call self%aabb%initialize(refinement_levels=refinement_levels_, facet=self%facet)
   endsubroutine create_aabb_tree

   elemental subroutine destroy(self)
   !< Destroy file.
   class(file_stl_object), intent(inout) :: self  !< File STL.
   type(file_stl_object)                 :: fresh !< Fresh instance of file STL.

   self = fresh
   endsubroutine destroy

   pure function distance(self, point, is_signed, sign_algorithm, is_square_root)
   !< Compute the (minimum) distance from a point to the triangulated surface.
   !<
   !< @note STL's metrix must be already computed.
   class(file_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),       intent(in)           :: point           !< Point coordinates.
   logical,                intent(in), optional :: is_signed       !< Sentinel to trigger signed distance.
   character(*),           intent(in), optional :: sign_algorithm  !< Algorithm used for "point in polyhedron" test.
   logical,                intent(in), optional :: is_square_root  !< Sentinel to trigger square-root distance.
   real(R8P)                                    :: distance        !< Minimum distance from point to the triangulated surface.
   real(R8P)                                    :: distance_       !< Minimum distance, temporary buffer.
   character(len=:), allocatable                :: sign_algorithm_ !< Algorithm used for "point in polyhedron" test, local variable.
   integer(I4P)                                 :: f               !< Counter.

   if (self%facets_number > 0) then
      if (self%aabb%is_initialized) then
         ! exploit AABB refinement levels
         distance = self%aabb%distance(point=point)
      else
         ! brute search over all facets
         distance = MaxR8P
         do f=1, self%facets_number
            distance_ = self%facet(f)%distance(point=point)
            if (abs(distance_) <= abs(distance)) distance = distance_
         enddo
      endif
   endif

   if (present(is_square_root)) then
      if (is_square_root) distance = sqrt(distance)
   endif

   if (present(is_signed)) then
      if (is_signed) then
        sign_algorithm_ = 'ray_intersections' ; if (present(sign_algorithm)) sign_algorithm_ = sign_algorithm
        select case(sign_algorithm_)
        case('solid_angle')
           if (self%is_point_inside_polyhedron_sa(point=point)) distance = -distance
        case('ray_intersections')
           if (self%is_point_inside_polyhedron_ri(point=point)) distance = -distance
        case default
          ! raise error: "unknown point in polyhedron algorithm"
        endselect
      endif
   endif
   endfunction distance

   pure function is_point_inside_polyhedron_ri(self, point) result(is_inside)
   !< Determinate is a point is inside or not to a polyhedron described by STL facets by means ray intersections count.
   !<
   !< @note STL's metrix must be already computed.
   class(file_stl_object), intent(in) :: self           !< File STL.
   type(vector_R8P),       intent(in) :: point          !< Point coordinates.
   logical                            :: is_inside      !< Check result.
   logical                            :: is_inside_by_x !< Test result by x-aligned ray intersections.
   logical                            :: is_inside_by_y !< Test result by y-aligned ray intersections.
   logical                            :: is_inside_by_z !< Test result by z-aligned ray intersections.

   is_inside_by_x = is_inside_by_ray_intersect(ray_origin=point, ray_direction=ex_R8P)
   is_inside_by_y = is_inside_by_ray_intersect(ray_origin=point, ray_direction=ey_R8P)
   if (is_inside_by_x.and.is_inside_by_y) then
     is_inside = .true.
   else
      is_inside_by_z = is_inside_by_ray_intersect(ray_origin=point, ray_direction=ez_R8P)
      is_inside = ((is_inside_by_x.and.is_inside_by_y).or.&
                   (is_inside_by_x.and.is_inside_by_z).or.&
                   (is_inside_by_y.and.is_inside_by_z))
   endif
   contains
      pure function is_inside_by_ray_intersect(ray_origin, ray_direction) result(is_inside_by)
      !< Generic line intersect test.
      type(vector_R8P), intent(in) :: ray_origin           !< Ray origin.
      type(vector_R8P), intent(in) :: ray_direction        !< Ray direction.
      integer(I4P)                 :: intersections_number !< Winding number of STL polyhedra with respect point.
      integer(I4P)                 :: f                    !< Counter.
      logical                      :: is_inside_by         !< Test result.

      intersections_number = 0
      do f=1, self%facets_number
         if (self%facet(f)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) &
            intersections_number = intersections_number + 1
      enddo
      if (mod(intersections_number, 2) == 0) then
        is_inside_by = .false.
      else
        is_inside_by = .true.
      endif
      endfunction is_inside_by_ray_intersect
   endfunction is_point_inside_polyhedron_ri

   pure function is_point_inside_polyhedron_sa(self, point) result(is_inside)
   !< Determinate is a point is inside or not to a polyhedron described by STL facets by means of the solid angle criteria.
   !<
   !< @note STL's metrix must be already computed.
   class(file_stl_object), intent(in) :: self        !< File STL.
   type(vector_R8P),       intent(in) :: point       !< Point coordinates.
   logical                            :: is_inside   !< Check result.
   real(R8P)                          :: solid_angle !< Solid angle of STL polyhedra projected on point unit sphere.
   integer(I4P)                       :: f           !< Counter.

   solid_angle = 0._R8P
   do f=1, self%facets_number
      solid_angle = solid_angle + self%facet(f)%solid_angle(point=point)
   enddo
   if (solid_angle < -2._R8P * PI .or. 2._R8P * PI < solid_angle) then
     is_inside = .true.
   else
     is_inside = .false.
   endif
   endfunction is_point_inside_polyhedron_sa

   elemental subroutine initialize(self, skip_destroy, file_name, is_ascii)
   !< Initialize file.
   class(file_stl_object), intent(inout)        :: self          !< File STL.
   logical,                intent(in), optional :: skip_destroy  !< Flag to skip destroy file.
   character(*),           intent(in), optional :: file_name     !< File name.
   logical,                intent(in), optional :: is_ascii      !< Sentinel to check if file is ASCII.
   logical                                      :: skip_destroy_ !< Flag to skip destroy file, local variable.

   skip_destroy_ = .false. ; if (present(skip_destroy)) skip_destroy_ = skip_destroy
   if (.not.skip_destroy_) call self%destroy
   if (present(file_name)) self%file_name = trim(adjustl(file_name))
   if (present(is_ascii)) self%is_ascii = is_ascii
   endsubroutine initialize

   subroutine load_from_file(self, file_name, is_ascii, guess_format)
   !< Load from file.
   class(file_stl_object), intent(inout)        :: self         !< File STL.
   character(*),           intent(in), optional :: file_name    !< File name.
   logical,                intent(in), optional :: is_ascii     !< Sentinel to check if file is ASCII.
   logical,                intent(in), optional :: guess_format !< Sentinel to try to guess format directly from file.
   integer(I4P)                                 :: f            !< Counter.

   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='read', guess_format=guess_format)
   call self%load_facets_number_from_file
   call self%allocate_facets
   call self%load_header_from_file
   if (self%is_ascii) then
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_ascii(file_unit=self%file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_binary(file_unit=self%file_unit)
      enddo
   endif
   call self%close_file
   endsubroutine load_from_file

   subroutine open_file(self, file_action, guess_format)
   !< Open file, once initialized.
   class(file_stl_object), intent(inout)        :: self          !< File STL.
   character(*),           intent(in)           :: file_action   !< File action, "read" or "write".
   logical,                intent(in), optional :: guess_format  !< Sentinel to try to guess format directly from file.
   logical                                      :: guess_format_ !< Sentinel to try to guess format directly from file, local var.
   logical                                      :: file_exist    !< Sentinel to check if file exist.
   character(5)                                 :: ascii_header  !< Ascii header sentinel.

   if (allocated(self%file_name)) then
      select case(trim(adjustl(file_action)))
      case('read')
         guess_format_ = .false. ; if (present(guess_format)) guess_format_ = guess_format
         inquire(file=self%file_name, exist=file_exist)
         if (file_exist) then
            if (guess_format_) then
               open(newunit=self%file_unit, file=self%file_name, form='formatted')
               read(self%file_unit, '(A)') ascii_header
               close(self%file_unit)
               if (ascii_header=='solid') then
                  self%is_ascii = .true.
               else
                  self%is_ascii = .false.
               endif
            endif
            if (self%is_ascii) then
               open(newunit=self%file_unit, file=self%file_name,                  form='formatted')
            else
               open(newunit=self%file_unit, file=self%file_name, access='stream', form='unformatted')
            endif
            self%is_open = .true.
         else
            write(stderr, '(A)') 'error: file "'//self%file_name//'" does not exist, impossible to open file!'
         endif
      case('write')
         if (self%is_ascii) then
            open(newunit=self%file_unit, file=self%file_name,                  form='formatted')
         else
            open(newunit=self%file_unit, file=self%file_name, access='stream', form='unformatted')
         endif
         self%is_open = .true.
      case default
         write(stderr, '(A)') 'error: file action "'//trim(adjustl(file_action))//'" unknown!'
      endselect
   else
      write(stderr, '(A)') 'error: file name has not be initialized, impossible to open file!'
   endif
   endsubroutine open_file

   elemental subroutine sanitize_normals(self)
   !< Sanitize normals, make normals consistent with vertices.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%sanitize_normal
   endsubroutine sanitize_normals

   subroutine save_into_file(self, file_name, is_ascii)
   !< Save into file.
   class(file_stl_object), intent(inout)        :: self      !< File STL.
   character(*),           intent(in), optional :: file_name !< File name.
   logical,                intent(in), optional :: is_ascii  !< Sentinel to check if file is ASCII.
   integer(I4P)                                 :: f         !< Counter.

   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='write')
   call self%save_header_into_file
   if (self%is_ascii) then
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_ascii(file_unit=self%file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_binary(file_unit=self%file_unit)
      enddo
   endif
   call self%close_file
   endsubroutine save_into_file

   ! operators
   ! =
   pure subroutine file_stl_assign_file_stl(lhs, rhs)
   !< Operator `=`.
   class(file_stl_object), intent(inout) :: lhs !< Left hand side.
   type(file_stl_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(lhs%file_name)) deallocate(lhs%file_name)
   if (allocated(rhs%file_name)) lhs%file_name = rhs%file_name
   lhs%file_unit = rhs%file_unit
   lhs%header = rhs%header
   lhs%facets_number = rhs%facets_number
   if (allocated(lhs%facet)) deallocate(lhs%facet)
   if (allocated(rhs%facet)) allocate(lhs%facet(1:lhs%facets_number), source=rhs%facet)
   lhs%aabb = rhs%aabb
   lhs%is_ascii = rhs%is_ascii
   lhs%is_open = rhs%is_open
   endsubroutine file_stl_assign_file_stl

   ! private methods
   elemental subroutine allocate_facets(self)
   !< Allocate facets.
   !<
   !< @note Facets previously allocated are lost.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) then
      if (allocated(self%facet))  deallocate(self%facet)
      allocate(self%facet(1:self%facets_number))
   endif
   endsubroutine allocate_facets

   subroutine load_facets_number_from_file(self)
   !< Load facets number from file.
   !<
   !< @note File is rewinded.
   class(file_stl_object), intent(inout) :: self         !< File STL.
   character(FRLEN)                      :: facet_record !< Facet record string buffer.

   if (self%is_open) then
      self%facets_number = 0
      rewind(self%file_unit)
      if (self%is_ascii) then
         do
            read(self%file_unit, '(A)', end=10, err=10) facet_record
            if (index(string=facet_record, substring='facet normal') > 0) self%facets_number = self%facets_number + 1
         enddo
      else
         read(self%file_unit, end=10, err=10) facet_record
         read(self%file_unit, end=10, err=10) self%facets_number
      endif
      10 rewind(self%file_unit)
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load facets number from file!'
   endif
   endsubroutine load_facets_number_from_file

   subroutine load_header_from_file(self)
   !< Load header from file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         read(self%file_unit, '(A)') self%header
         self%header = trim(self%header(index(self%header, 'solid')+1:))
      else
         read(self%file_unit) self%header
         read(self%file_unit) self%facets_number
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine load_header_from_file

   subroutine save_header_into_file(self)
   !< Save header into file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         write(self%file_unit, '(A)') 'solid '//trim(self%header)
      else
         write(self%file_unit) self%header
         write(self%file_unit) self%facets_number
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine save_header_into_file

   subroutine save_trailer_into_file(self)
   !< Save trailer into file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      if (self%is_ascii) write(self%file_unit, '(A)') 'endsolid '//trim(self%header)
   else
      write(stderr, '(A)') 'error: file is not open, impossible to write trailer into file!'
   endif
   endsubroutine save_trailer_into_file
endmodule fossil_file_stl_object
