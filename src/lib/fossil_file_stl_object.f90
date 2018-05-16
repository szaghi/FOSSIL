!< FOSSIL,  STL file class definition.

module fossil_file_stl_object
!< FOSSIL,  STL file class definition.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_utils, only : EPS, PI, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: file_stl_object

type :: file_stl_object
   !< FOSSIL STL file class.
   character(len=:), allocatable   :: file_name            !< File name
   integer(I4P)                    :: file_unit=0          !< File unit.
   character(FRLEN)                :: header               !< File header.
   integer(I4P)                    :: facets_number=0      !< Facets number.
   type(facet_object), allocatable :: facet(:)             !< Facets.
   integer(I4P), allocatable       :: facet_1_de(:)        !< Facets with one disconnected edge.
   integer(I4P)                    :: facets_1_de_number=0 !< Facets number with one disconnected edge.
   integer(I4P), allocatable       :: facet_2_de(:)        !< Facets with two disconnected edges.
   integer(I4P)                    :: facets_2_de_number=0 !< Facets number with two disconnected edge.
   integer(I4P), allocatable       :: facet_3_de(:)        !< Facets with three disconnected edges.
   integer(I4P)                    :: facets_3_de_number=0 !< Facets number with three disconnected edge.
   type(aabb_tree_object)          :: aabb                 !< AABB tree.
   type(vector_R8P)                :: bmin                 !< Minimum point of STL.
   type(vector_R8P)                :: bmax                 !< Maximum point of STL.
   real(R8P)                       :: volume=0._R8P        !< Volume bounded by STL surface.
   logical                         :: is_ascii=.true.      !< Sentinel to check if file is ASCII.
   logical                         :: is_open=.false.      !< Sentinel to check if file is open.
   contains
      ! public methods
      procedure, pass(self) :: analize                         !< Analize STL.
      procedure, pass(self) :: build_connectivity              !< Build facets connectivity.
      procedure, pass(self) :: close_file                      !< Close file.
      procedure, pass(self) :: compute_metrix                  !< Compute facets metrix.
      procedure, pass(self) :: compute_normals                 !< Compute facets normals by means of vertices data.
      procedure, pass(self) :: compute_volume                  !< Compute volume bounded by STL surface.
      procedure, pass(self) :: connect_nearby_vertices         !< Connect nearby vertices of disconnected edges.
      procedure, pass(self) :: create_aabb_tree                !< Create the AABB tree.
      procedure, pass(self) :: destroy                         !< Destroy file.
      procedure, pass(self) :: distance                        !< Compute the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: initialize                      !< Initialize file.
      procedure, pass(self) :: is_point_inside_polyhedron_ri   !< Determinate is point is inside or not STL facets by ray intersect.
      procedure, pass(self) :: is_point_inside_polyhedron_sa   !< Determinate is point is inside or not STL facets by solid angle.
      procedure, pass(self) :: load_from_file                  !< Load from file.
      generic               :: mirror => mirror_by_normal, &
                                         mirror_by_matrix      !< Mirror facets.
      procedure, pass(self) :: open_file                       !< Open file, once initialized.
      procedure, pass(self) :: reverse_normals                 !< Reverse facets normals.
      procedure, pass(self) :: resize                          !< Resize (scale) facets by x or y or z or vectorial factors.
      generic               :: rotate => rotate_by_axis_angle, &
                                         rotate_by_matrix      !< Rotate facets.
      procedure, pass(self) :: sanitize                        !< Sanitize STL.
      procedure, pass(self) :: sanitize_normals                !< Sanitize facets normals, make them consistent.
      procedure, pass(self) :: save_into_file                  !< Save into file.
      procedure, pass(self) :: smallest_edge_len               !< Return the smallest edge length.
      procedure, pass(self) :: statistics                      !< Return STL statistics.
      procedure, pass(self) :: translate                       !< Translate facet given vectorial delta.
      ! operators
      generic :: assignment(=) => file_stl_assign_file_stl       !< Overload `=`.
      procedure, pass(lhs),  private :: file_stl_assign_file_stl !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: allocate_facets              !< Allocate facets.
      procedure, pass(self), private :: compute_facets_disconnected  !< Compute facets with disconnected edges.
      procedure, pass(self), private :: load_facets_number_from_file !< Load facets number from file.
      procedure, pass(self), private :: load_header_from_file        !< Load header from file.
      procedure, pass(self), private :: mirror_by_normal             !< Mirror facets given normal of mirroring plane.
      procedure, pass(self), private :: mirror_by_matrix             !< Mirror facets given matrix.
      procedure, pass(self), private :: rotate_by_axis_angle         !< Rotate facets given axis and angle.
      procedure, pass(self), private :: rotate_by_matrix             !< Rotate facets given matrix.
      procedure, pass(self), private :: save_header_into_file        !< Save header into file.
      procedure, pass(self), private :: save_trailer_into_file       !< Save trailer into file.
endtype file_stl_object

contains
   ! public methods
   pure subroutine analize(self)
   !< Analize STL.
   !<
   !< Buil connectivity, compute metrix, compute volume.
   class(file_stl_object), intent(inout) :: self   !< File STL.

   if (self%facets_number>0) then
      call self%build_connectivity
      call self%compute_metrix
      call self%compute_facets_disconnected
      call self%compute_volume
   endif
   endsubroutine analize

   pure subroutine build_connectivity(self)
   !< Build facets connectivity.
   class(file_stl_object), intent(inout) :: self              !< File STL.
   real(R8P)                             :: smallest_edge_len !< Smallest edge length.
   integer(I4P)                          :: f1, f2            !< Counter.

   if (self%facets_number>0) then
      call self%facet%destroy_connectivity
      smallest_edge_len = self%smallest_edge_len() * 0.9_R8P
      do f1=1, self%facets_number - 1
         do f2=f1 + 1, self%facets_number
            call self%facet(f1)%compute_vertices_nearby(other=self%facet(f2), &
                                                        tolerance_to_be_identical=EPS, tolerance_to_be_nearby=smallest_edge_len)
         enddo
      enddo
      do f1=1, self%facets_number
         call self%facet(f1)%update_connectivity
      enddo
   endif
   endsubroutine build_connectivity

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
      ! computing bounding box extents
      self%bmin%x = minval(self%facet(:)%bb(1)%x)
      self%bmin%y = minval(self%facet(:)%bb(1)%y)
      self%bmin%z = minval(self%facet(:)%bb(1)%z)
      self%bmax%x = maxval(self%facet(:)%bb(2)%x)
      self%bmax%y = maxval(self%facet(:)%bb(2)%y)
      self%bmax%z = maxval(self%facet(:)%bb(2)%z)
   endif
   endsubroutine compute_metrix

   elemental subroutine compute_normals(self)
   !< Compute facets normals by means of vertices data.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%compute_normal
   endsubroutine compute_normals

   elemental subroutine compute_volume(self)
   !< Compute volume bounded by STL surface.
   class(file_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                          :: f    !< Counter.

   if (self%facets_number>0) then
      self%volume = 0._R8P
      do f=1, self%facets_number
         self%volume = self%volume + self%facet(f)%tetrahedron_volume(apex=self%facet(1)%vertex(1))
      enddo
   endif
   endsubroutine compute_volume

   pure subroutine connect_nearby_vertices(self)
   !< Connect nearby vertices of disconnected edges.
   class(file_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                          :: f    !< Counter.

   if (self%facets_number>0) then
      if     (self%facets_1_de_number>0) then
         do f=1, self%facets_1_de_number
            call self%facet(self%facet_1_de(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
      if (self%facets_2_de_number>0) then
         do f=1, self%facets_2_de_number
            call self%facet(self%facet_2_de(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
      if (self%facets_3_de_number>0) then
         do f=1, self%facets_3_de_number
            call self%facet(self%facet_3_de(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
   endif
   endsubroutine connect_nearby_vertices

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
         ! brute-force search over all facets
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

   is_inside_by_x = is_inside_by_ray_intersect(ray_origin=point, ray_direction=      ex_R8P + EPS * ey_R8P + EPS * ez_R8P)
   is_inside_by_y = is_inside_by_ray_intersect(ray_origin=point, ray_direction=EPS * ex_R8P +       ey_R8P + EPS * ez_R8P)
   if (is_inside_by_x.and.is_inside_by_y) then
     is_inside = .true.
   else
      is_inside_by_z = is_inside_by_ray_intersect(ray_origin=point, ray_direction=EPS * ex_R8P + EPS * ey_R8P + ez_R8P)
      is_inside = ((is_inside_by_x.and.is_inside_by_y).or.&
                   (is_inside_by_x.and.is_inside_by_z).or.&
                   (is_inside_by_y.and.is_inside_by_z))
   endif
   contains
      pure function is_inside_by_ray_intersect(ray_origin, ray_direction) result(is_inside_by)
      !< Generic line intersect test.
      type(vector_R8P), intent(in) :: ray_origin           !< Ray origin.
      type(vector_R8P), intent(in) :: ray_direction        !< Ray direction.
      integer(I4P)                 :: intersections_number !< Ray intersections number of STL polyhedra with respect point.
      integer(I4P)                 :: f                    !< Counter.
      logical                      :: is_inside_by         !< Test result.

      intersections_number = 0

      if (self%aabb%is_initialized) then
         ! exploit AABB refinement levels
         intersections_number = self%aabb%ray_intersections_number(ray_origin=ray_origin, ray_direction=ray_direction)
      else
         ! brute-force search over all facets
         do f=1, self%facets_number
            if (self%facet(f)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) &
               intersections_number = intersections_number + 1
         enddo
      endif

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

   subroutine load_from_file(self, file_name, is_ascii, guess_format, disable_analysis)
   !< Load from file.
   class(file_stl_object), intent(inout)        :: self              !< File STL.
   character(*),           intent(in), optional :: file_name         !< File name.
   logical,                intent(in), optional :: is_ascii          !< Sentinel to check if file is ASCII.
   logical,                intent(in), optional :: guess_format      !< Sentinel to try to guess format directly from file.
   logical,                intent(in), optional :: disable_analysis  !< Sentinel to disable STL analysis.
   logical                                      :: disable_analysis_ !< Sentinel to disable STL analysis, local variable.
   integer(I4P)                                 :: f                 !< Counter.

   disable_analysis_ = .false. ; if (present(disable_analysis)) disable_analysis_ = disable_analysis
   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='read', guess_format=guess_format)
   call self%load_facets_number_from_file
   call self%allocate_facets
   call self%load_header_from_file
   if (self%is_ascii) then
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_ascii(file_unit=self%file_unit)
         self%facet(f)%id = f
      enddo
   else
      do f=1, self%facets_number
         self%facet(f)%id = f
         call self%facet(f)%load_from_file_binary(file_unit=self%file_unit)
      enddo
   endif
   call self%close_file
   if (.not.disable_analysis_) call self%analize
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

   elemental subroutine resize(self, x, y, z, factor, recompute_metrix)
   !< Resize (scale) facets by x or y or z or vectorial factors.
   !<
   !< @note The name `scale` has not been used, it been a Fortran built-in.
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),              intent(in), optional :: x                !< Factor along x axis.
   real(R8P),              intent(in), optional :: y                !< Factor along y axis.
   real(R8P),              intent(in), optional :: z                !< Factor along z axis.
   type(vector_R8P),       intent(in), optional :: factor           !< Vectorial factor.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   type(vector_R8P)                             :: factor_          !< Vectorial factor, local variable.

   if (self%facets_number>0) then
      factor_ = 1._R8P
      if (present(factor)) then
         factor_ = factor
      else
         if (present(x)) factor_%x = x
         if (present(y)) factor_%y = y
         if (present(z)) factor_%z = z
      endif
      call self%facet%resize(factor=factor_, recompute_metrix=recompute_metrix)
   endif
   endsubroutine resize

   elemental subroutine reverse_normals(self)
   !< Reverse facets normals.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%reverse_normal
   endsubroutine reverse_normals

   pure subroutine sanitize(self, do_analysis)
   !< Sanitize STL.
   class(file_stl_object), intent(inout)        :: self !< File STL.
   logical,                intent(in), optional :: do_analysis !< Sentil for performing a first analysis.

   if (self%facets_number>0) then
      if (present(do_analysis)) then
         if (do_analysis) call self%analize
      endif
      if (self%facets_1_de_number>0.or.self%facets_2_de_number>0.or.self%facets_3_de_number>0) call self%connect_nearby_vertices
      call self%analize
      call self%sanitize_normals
   endif
   endsubroutine sanitize

   pure subroutine sanitize_normals(self)
   !< Sanitize facets normals, make them consistent.
   !<
   !< @note Facets connectivity and normals must be already computed.
   class(file_stl_object), intent(inout) :: self             !< File STL.
   logical, allocatable                  :: facet_checked(:) !< List of facets checked.
   integer(I4P)                          :: f, ff            !< Counter.

   if (self%facets_number>0) then
      allocate(facet_checked(1:self%facets_number))
      facet_checked = .false.
      f = 1
      facet_checked(f) = .true.
      do
         ff = 0
         if (self%facet(f)%fcon_edge_12>0) then
            if (.not.facet_checked(self%facet(f)%fcon_edge_12)) then
               call self%facet(f)%make_normal_consistent(edge_dir='edge_12', other=self%facet(self%facet(f)%fcon_edge_12))
               facet_checked(self%facet(f)%fcon_edge_12) = .true.
               ff = self%facet(f)%fcon_edge_12
            endif
         endif
         if (self%facet(f)%fcon_edge_23>0) then
            if (.not.facet_checked(self%facet(f)%fcon_edge_23)) then
               call self%facet(f)%make_normal_consistent(edge_dir='edge_23', other=self%facet(self%facet(f)%fcon_edge_23))
               facet_checked(self%facet(f)%fcon_edge_23) = .true.
               ff = self%facet(f)%fcon_edge_23
            endif
         endif
         if (self%facet(f)%fcon_edge_31>0) then
            if (.not.facet_checked(self%facet(f)%fcon_edge_31)) then
               call self%facet(f)%make_normal_consistent(edge_dir='edge_31', other=self%facet(self%facet(f)%fcon_edge_31))
               facet_checked(self%facet(f)%fcon_edge_31) = .true.
               ff = self%facet(f)%fcon_edge_31
            endif
         endif
         if (ff==0) then
            exit
         else
            f = ff
         endif
      enddo
   endif
   call self%compute_volume
   if (self%volume < 0) call self%reverse_normals
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
   call self%save_trailer_into_file
   call self%close_file
   endsubroutine save_into_file

   pure function smallest_edge_len(self) result(smallest)
   !< Return the smallest edge length.
   class(file_stl_object), intent(in) :: self     !< File STL.
   real(R8P)                          :: smallest !< Smallest edge length.
   integer(I4P)                       :: f        !< Counter.

   smallest = MaxR8P
   if (self%facets_number>0) then
      do f=1, self%facets_number
         smallest = min(smallest, self%facet(f)%smallest_edge_len())
      enddo
   endif
   endfunction smallest_edge_len

   pure function statistics(self, prefix) result(stats)
   !< Return STL statistics.
   class(file_stl_object), intent(in)           :: self             !< File STL.
   character(*),           intent(in), optional :: prefix           !< Lines prefix.
   character(len=:), allocatable                :: stats            !< STL statistics.
   character(len=:), allocatable                :: prefix_          !< Lines prefix, local variable.
   character(1), parameter                      :: NL=new_line('a') !< Line terminator.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   stats = prefix_//self%header//NL
   if (allocated(self%file_name)) stats=stats//prefix_//'file name:   '//self%file_name//NL
   if (self%is_ascii) then
      stats=stats//prefix_//'file format: ascii'//NL
   else
      stats=stats//prefix_//'file format: binary'//NL
   endif
   if (self%facets_number > 0) then
      stats=stats//prefix_//'X extents: ['//trim(str(self%bmin%x))//', '//trim(str(self%bmax%x))//']'//NL
      stats=stats//prefix_//'Y extents: ['//trim(str(self%bmin%y))//', '//trim(str(self%bmax%y))//']'//NL
      stats=stats//prefix_//'Z extents: ['//trim(str(self%bmin%z))//', '//trim(str(self%bmax%z))//']'//NL
      stats=stats//prefix_//'volume: '//trim(str(self%volume))//NL
      stats=stats//prefix_//'number of facets: '//trim(str(self%facets_number))//NL
      stats=stats//prefix_//'number of facets with 1 edges disconnected: '//trim(str(self%facets_1_de_number))//NL
      stats=stats//prefix_//'number of facets with 2 edges disconnected: '//trim(str(self%facets_2_de_number))//NL
      stats=stats//prefix_//'number of facets with 3 edges disconnected: '//trim(str(self%facets_3_de_number))!//NL
   endif
   endfunction statistics

   elemental subroutine translate(self, x, y, z, delta, recompute_metrix)
   !< Translate facets x or y or z or vectorial delta increments.
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),              intent(in), optional :: x                !< Increment along x axis.
   real(R8P),              intent(in), optional :: y                !< Increment along y axis.
   real(R8P),              intent(in), optional :: z                !< Increment along z axis.
   type(vector_R8P),       intent(in), optional :: delta            !< Vectorial increment.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   type(vector_R8P)                             :: delta_           !< Vectorial increment, local variable.

   if (self%facets_number>0) then
      delta_ = 0._R8P
      if (present(delta)) then
         delta_ = delta
      else
         if (present(x)) delta_%x = x
         if (present(y)) delta_%y = y
         if (present(z)) delta_%z = z
      endif
      call self%facet%translate(delta=delta_, recompute_metrix=recompute_metrix)
   endif
   endsubroutine translate

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
   lhs%facets_1_de_number = rhs%facets_1_de_number
   if (allocated(lhs%facet_1_de)) deallocate(lhs%facet_1_de)
   if (allocated(rhs%facet_1_de)) allocate(lhs%facet_1_de, source=rhs%facet_1_de)
   lhs%facets_2_de_number = rhs%facets_2_de_number
   if (allocated(lhs%facet_2_de)) deallocate(lhs%facet_2_de)
   if (allocated(rhs%facet_2_de)) allocate(lhs%facet_2_de, source=rhs%facet_2_de)
   lhs%facets_3_de_number = rhs%facets_3_de_number
   if (allocated(lhs%facet_3_de)) deallocate(lhs%facet_3_de)
   if (allocated(rhs%facet_3_de)) allocate(lhs%facet_3_de, source=rhs%facet_3_de)
   lhs%aabb = rhs%aabb
   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
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

   pure subroutine compute_facets_disconnected(self)
   !< Compute facets with disconnected edges.
   class(file_stl_object), intent(inout) :: self  !< File STL.
   logical                               :: de(3) !< Flag to check edges disconnection.
   integer(I4P)                          :: f     !< Counter.

   if (self%facets_number>0) then
      self%facets_1_de_number = 0
      self%facets_2_de_number = 0
      self%facets_3_de_number = 0
      do f=1, self%facets_number
         de = .false.
         if (self%facet(f)%fcon_edge_12==0_I4P) de(1) = .true.
         if (self%facet(f)%fcon_edge_23==0_I4P) de(2) = .true.
         if (self%facet(f)%fcon_edge_31==0_I4P) de(3) = .true.
         select case(count(de))
         case(1_I4P)
            self%facets_1_de_number = self%facets_1_de_number + 1
            call add_to_de_list(id=self%facet(f)%id, list=self%facet_1_de)
         case(2_I4P)
            self%facets_2_de_number = self%facets_2_de_number + 1
            call add_to_de_list(id=self%facet(f)%id, list=self%facet_2_de)
         case(3_I4P)
            self%facets_3_de_number = self%facets_3_de_number + 1
            call add_to_de_list(id=self%facet(f)%id, list=self%facet_3_de)
         endselect
      enddo
   endif
   contains
      pure subroutine add_to_de_list(id, list)
      !< Add given facet ID to a list of disconnected facets.
      integer(I4P),              intent(in)    :: id      !< Facet global ID.
      integer(I4P), allocatable, intent(inout) :: list(:) !< Disconnected facets list.
      integer(I4P), allocatable                :: tmp(:)  !< Disconnected facets list, temporary variable.
      integer(I4P)                             :: n       !< Size of input list.

      if (allocated(list)) then
         n = size(list, dim=1)
         allocate(tmp(1:n+1))
         tmp(1:n) = list
         tmp(n+1) = id
         call move_alloc(from=tmp, to=list)
      else
         allocate(list(1))
         list(1) = id
      endif
      endsubroutine add_to_de_list
   endsubroutine compute_facets_disconnected

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
         self%header = trim(adjustl(self%header(index(self%header, 'solid')+6:)))
      else
         read(self%file_unit) self%header
         read(self%file_unit) self%facets_number
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine load_header_from_file

   elemental subroutine mirror_by_normal(self, normal, recompute_metrix)
   !< Mirror facets given normal of mirroring plane.
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),       intent(in)           :: normal           !< Normal of mirroring plane.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                    :: matrix(3,3)      !< Mirroring matrix.
   integer(I4P)                                 :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = mirror_matrix_R8P(normal=normal)
      do f=1, self%facets_number
         call self%facet(f)%mirror(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine mirror_by_normal

   pure subroutine mirror_by_matrix(self, matrix, recompute_metrix)
   !< Mirror facet given matrix (of mirroring).
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),              intent(in)           :: matrix(3,3)      !< Mirroring matrix.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                 :: f                !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         call self%facet(f)%mirror(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine mirror_by_matrix

   elemental subroutine rotate_by_axis_angle(self, axis, angle, recompute_metrix)
   !< Rotate facets given axis and angle.
   !<
   !< Angle must be in radiants.
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),       intent(in)           :: axis             !< Axis of rotation.
   real(R8P),              intent(in)           :: angle            !< Angle of rotation.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                    :: matrix(3,3)      !< Rotation matrix.
   integer(I4P)                                 :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = rotation_matrix_R8P(axis=axis, angle=angle)
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_axis_angle

   pure subroutine rotate_by_matrix(self, matrix, recompute_metrix)
   !< Rotate facet given matrix (of ratation).
   class(file_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),              intent(in)           :: matrix(3,3)      !< Rotation matrix.
   logical,                intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                 :: f                !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_matrix

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
