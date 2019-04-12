!< FOSSIL, STL surface class definition.

module fossil_surface_stl_object
!< FOSSIL, STL surface class definition.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use fossil_utils, only : EPS, PI, is_inside_bb
use penf, only : I4P, R8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: surface_stl_object

type :: surface_stl_object
   !< FOSSIL STL surface class.
   integer(I4P)                    :: facets_number=0 !< Facets number.
   type(facet_object), allocatable :: facet(:)        !< Facets.
   type(list_id_object)            :: facet_1_de      !< Facets with one disconnected edges.
   type(list_id_object)            :: facet_2_de      !< Facets with two disconnected edges.
   type(list_id_object)            :: facet_3_de      !< Facets with three disconnected edges.
   type(aabb_tree_object)          :: aabb            !< AABB tree.
   type(vector_R8P)                :: bmin            !< Minimum point of STL.
   type(vector_R8P)                :: bmax            !< Maximum point of STL.
   real(R8P)                       :: volume=0._R8P   !< Volume bounded by STL surface.
   type(vector_R8P)                :: centroid        !< Centroid of STL surface.
   contains
      ! public methods
      procedure, pass(self) :: allocate_facets                 !< Allocate facets.
      procedure, pass(self) :: analize                         !< Analize STL.
      procedure, pass(self) :: build_connectivity              !< Build facets connectivity.
      procedure, pass(self) :: clip                            !< Clip triangulated surface given an AABB.
      procedure, pass(self) :: compute_centroid                !< Compute centroid of STL surface.
      procedure, pass(self) :: compute_distance                !< Compute the (minimum) distance returning also the closest point.
      procedure, pass(self) :: compute_metrix                  !< Compute facets metrix.
      procedure, pass(self) :: compute_normals                 !< Compute facets normals by means of vertices data.
      procedure, pass(self) :: compute_volume                  !< Compute volume bounded by STL surface.
      procedure, pass(self) :: connect_nearby_vertices         !< Connect nearby vertices of disconnected edges.
      procedure, pass(self) :: destroy                         !< Destroy file.
      procedure, pass(self) :: distance                        !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: initialize                      !< Initialize file.
      procedure, pass(self) :: is_point_inside_polyhedron_ri   !< Determinate is point is inside or not STL facets by ray intersect.
      procedure, pass(self) :: is_point_inside_polyhedron_sa   !< Determinate is point is inside or not STL facets by solid angle.
      procedure, pass(self) :: merge_solids                    !< Merge facets with ones of other STL file.
      generic               :: mirror => mirror_by_normal, &
                                         mirror_by_matrix      !< Mirror facets.
      procedure, pass(self) :: reverse_normals                 !< Reverse facets normals.
      procedure, pass(self) :: resize                          !< Resize (scale) facets by x or y or z or vectorial factors.
      generic               :: rotate => rotate_by_axis_angle, &
                                         rotate_by_matrix      !< Rotate facets.
      procedure, pass(self) :: sanitize                        !< Sanitize STL.
      procedure, pass(self) :: sanitize_normals                !< Sanitize facets normals, make them consistent.
      procedure, pass(self) :: smallest_edge_len               !< Return the smallest edge length.
      procedure, pass(self) :: statistics                      !< Return STL statistics.
      procedure, pass(self) :: translate                       !< Translate facet given vectorial delta.
      ! operators
      generic :: assignment(=) => surface_stl_assign_surface_stl       !< Overload `=`.
      procedure, pass(lhs),  private :: surface_stl_assign_surface_stl !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: compute_facets_disconnected !< Compute facets with disconnected edges.
      procedure, pass(self), private :: mirror_by_normal            !< Mirror facets given normal of mirroring plane.
      procedure, pass(self), private :: mirror_by_matrix            !< Mirror facets given matrix.
      procedure, pass(self), private :: rotate_by_axis_angle        !< Rotate facets given axis and angle.
      procedure, pass(self), private :: rotate_by_matrix            !< Rotate facets given matrix.
      procedure, pass(self), private :: set_facets_id               !< (Re)set facets ID.
endtype surface_stl_object

contains
   ! public methods
   elemental subroutine allocate_facets(self, facets_number)
   !< Allocate facets.
   !<
   !< @note Facets previously allocated are lost.
   class(surface_stl_object), intent(inout) :: self          !< File STL.
   integer(I4P),              intent(in)    :: facets_number !< Facets number.

   if (allocated(self%facet)) then
      call self%facet%destroy
      deallocate(self%facet)
   endif
   self%facets_number = facets_number
   if (self%facets_number>0) then
      allocate(self%facet(1:self%facets_number))
   endif
   endsubroutine allocate_facets

   elemental subroutine analize(self, aabb_refinement_levels)
   !< Analize STL.
   !<
   !< Buil connectivity, compute metrix, compute volume.
   class(surface_stl_object), intent(inout)        :: self                   !< File STL.
   integer(I4P),              intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.

   self%facets_number = 0
   if (allocated(self%facet)) self%facets_number = size(self%facet, dim=1)
   if (self%facets_number>0) then
      call self%set_facets_id
      call self%compute_metrix
      call self%aabb%initialize(refinement_levels=aabb_refinement_levels, facet=self%facet)
      call self%build_connectivity
      call self%compute_facets_disconnected
      call self%compute_volume
      call self%compute_centroid
   endif
   endsubroutine analize

   pure subroutine build_connectivity(self)
   !< Build facets connectivity.
   class(surface_stl_object), intent(inout) :: self              !< File STL.
   real(R8P)                                :: smallest_edge_len !< Smallest edge length.
   integer(I4P)                             :: f1, f2            !< Counter.
   type(aabb_tree_object)                   :: aabb              !< Temporary AABB tree.

   if (self%facets_number>0) then
      call self%facet%destroy_connectivity
      smallest_edge_len = self%smallest_edge_len() * 0.9_R8P
      if (self%aabb%is_initialized) then
         ! exploit AABB structure
         call aabb%initialize(facet=self%facet, refinement_levels=self%aabb%refinement_levels, do_facets_distribute=.false.)
         call aabb%distribute_facets(facet=self%facet, is_exclusive=.false., do_update_extents=.false.)
         call aabb%compute_vertices_nearby(facet=self%facet,              &
                                           tolerance_to_be_identical=EPS, &
                                           tolerance_to_be_nearby=smallest_edge_len)
      else
         ! brute-force search over all facets
         do f1=1, self%facets_number - 1
            do f2=f1 + 1, self%facets_number
               call self%facet(f1)%compute_vertices_nearby(other=self%facet(f2),          &
                                                           tolerance_to_be_identical=EPS, &
                                                           tolerance_to_be_nearby=smallest_edge_len)
            enddo
         enddo
      endif
      do f1=1, self%facets_number
         call self%facet(f1)%update_connectivity
      enddo
   endif
   endsubroutine build_connectivity

   subroutine clip(self, bmin, bmax, remainder)
   !< Clip triangulated surface given an AABB.
   class(surface_stl_object), intent(inout)         :: self              !< File STL.
   type(vector_R8P),          intent(in)            :: bmin, bmax        !< Bounding box extents.
   type(surface_stl_object),  intent(out), optional :: remainder         !< Remainder part of the triangulated surface.
   type(facet_object), allocatable                  :: facet(:)          !< Clipped facets.
   integer(I4P)                                     :: facets_in_number  !< Number of facets inside bounding box.
   integer(I4P)                                     :: facets_out_number !< Number of facets outside bounding box.
   integer(I4P)                                     :: f, fi, fo         !< Counter.

   if (self%facets_number>0) then
      facets_in_number = 0
      facets_out_number = 0
      do f=1, self%facets_number
         if (is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(1)).and.&
             is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(2)).and.&
             is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(3))) then
            facets_in_number = facets_in_number + 1
         else
            facets_out_number = facets_out_number + 1
         endif
      enddo
      if (facets_in_number>0) then
         allocate(facet(1:facets_in_number))
         if (present(remainder)) then
            remainder%facets_number = facets_out_number
            allocate(remainder%facet(1:facets_out_number))
         endif
         fi = 0
         fo = 0
         do f=1, self%facets_number
            if (is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(1)).and.&
                is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(2)).and.&
                is_inside_bb(bmin=bmin, bmax=bmax, point=self%facet(f)%vertex(3))) then
               fi = fi + 1
               facet(fi) = self%facet(f)
               facet(fi)%id = fi
            else
               fo = fo + 1
               if (present(remainder)) then
                  remainder%facet(fo) = self%facet(f)
                  remainder%facet(fo)%id = fo
               endif
            endif
         enddo
         call move_alloc(from=facet, to=self%facet)
         self%facets_number = facets_in_number
         call self%analize(aabb_refinement_levels=self%aabb%refinement_levels)
         if (present(remainder)) call remainder%analize(aabb_refinement_levels=self%aabb%refinement_levels)
      endif
   endif
   endsubroutine clip

   pure subroutine compute_centroid(self)
   !< Compute centroid of STL surface.
   !<
   !< @note Metrix and volume must be already computed.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      self%centroid = 0._R8P
      do f=1, self%facets_number
         self%centroid = self%centroid - self%facet(f)%centroid_part()
      enddo
      self%centroid = self%centroid / (48 * self%volume)
   endif
   endsubroutine compute_centroid

   pure subroutine compute_distance(self, point, distance, is_signed, sign_algorithm, is_square_root, &
                                    facet_index, edge_index, vertex_index)
   !< Compute the (minimum) distance returning also the closest point.
   class(surface_stl_object), intent(in)            :: self            !< File STL.
   type(vector_R8P),          intent(in)            :: point           !< Point coordinates.
   real(R8P),                 intent(out)           :: distance        !< Minimum distance.
   logical,                   intent(in),  optional :: is_signed       !< Sentinel to trigger signed distance.
   character(*),              intent(in),  optional :: sign_algorithm  !< Algorithm used for "point in polyhedron" test.
   logical,                   intent(in),  optional :: is_square_root  !< Sentinel to trigger square-root distance.
   integer(I4P),              intent(out), optional :: facet_index     !< Index of facet containing the closest point.
   integer(I4P),              intent(out), optional :: edge_index      !< Index of edge on facet containing the closest point.
   integer(I4P),              intent(out), optional :: vertex_index    !< Index of vertex on facet containing the closest point.
   real(R8P)                                        :: distance_       !< Minimum distance, temporary buffer.
   character(len=:), allocatable                    :: sign_algorithm_ !< Algorithm used for "point in polyhedron" test, local var.
   integer(I4P)                                     :: facet_index_    !< Index of facet containing the closest point, local var.
   integer(I4P)                                     :: f               !< Counter.

   if (self%facets_number > 0) then
      if (self%aabb%is_initialized) then
         ! exploit AABB refinement levels
         distance = self%aabb%distance(facet=self%facet, point=point)
      else
         ! brute-force search over all facets
         distance = MaxR8P
         do f=1, self%facets_number
            call self%facet(f)%compute_distance(point=point, distance=distance_)
            if (abs(distance_) <= abs(distance)) then
               facet_index_ = facet_index
               distance = distance_
            endif
         enddo
      endif
   endif
   if (present(facet_index)) facet_index = facet_index_

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
   endsubroutine compute_distance

   pure subroutine compute_metrix(self)
   !< Compute facets metrix.
   class(surface_stl_object), intent(inout) :: self !< File STL.

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
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%compute_normal
   endsubroutine compute_normals

   elemental subroutine compute_volume(self)
   !< Compute volume bounded by STL surface.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      self%volume = 0._R8P
      do f=1, self%facets_number
         self%volume = self%volume + self%facet(f)%tetrahedron_volume(apex=self%facet(1)%vertex(1))
      enddo
   endif
   endsubroutine compute_volume

   pure subroutine connect_nearby_vertices(self)
   !< Connect nearby vertices of disconnected edges.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      if (self%facet_1_de%ids_number>0) then
         do f=1, self%facet_1_de%ids_number
            call self%facet(self%facet_1_de%id(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
      if (self%facet_2_de%ids_number>0) then
         do f=1, self%facet_2_de%ids_number
            call self%facet(self%facet_2_de%id(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
      if (self%facet_3_de%ids_number>0) then
         do f=1, self%facet_3_de%ids_number
            call self%facet(self%facet_3_de%id(f))%connect_nearby_vertices(facet=self%facet)
         enddo
      endif
   endif
   endsubroutine connect_nearby_vertices

   elemental subroutine destroy(self)
   !< Destroy file.
   class(surface_stl_object), intent(inout) :: self  !< File STL.
   type(surface_stl_object)                 :: fresh !< Fresh instance of file STL.

   self = fresh
   endsubroutine destroy

   pure function distance(self, point, is_signed, sign_algorithm, is_square_root)
   !< Return the (minimum) distance from a point to the triangulated surface.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical,                   intent(in), optional :: is_signed       !< Sentinel to trigger signed distance.
   character(*),              intent(in), optional :: sign_algorithm  !< Algorithm used for "point in polyhedron" test.
   logical,                   intent(in), optional :: is_square_root  !< Sentinel to trigger square-root distance.
   real(R8P)                                       :: distance        !< Minimum distance from point to the triangulated surface.

   call self%compute_distance(point=point, distance=distance, &
                              is_signed=is_signed, sign_algorithm=sign_algorithm, is_square_root=is_square_root)
   endfunction distance

   pure function is_point_inside_polyhedron_ri(self, point) result(is_inside)
   !< Determinate is a point is inside or not to a polyhedron described by STL facets by means ray intersections count.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in) :: self           !< File STL.
   type(vector_R8P),          intent(in) :: point          !< Point coordinates.
   logical                               :: is_inside      !< Check result.
   logical                               :: is_inside_by_x !< Test result by x-aligned ray intersections.
   logical                               :: is_inside_by_y !< Test result by y-aligned ray intersections.
   logical                               :: is_inside_by_z !< Test result by z-aligned ray intersections.

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
         intersections_number = self%aabb%ray_intersections_number(facet=self%facet, &
                                                                   ray_origin=ray_origin, ray_direction=ray_direction)
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
   class(surface_stl_object), intent(in) :: self        !< File STL.
   type(vector_R8P),          intent(in) :: point       !< Point coordinates.
   logical                               :: is_inside   !< Check result.
   real(R8P)                             :: solid_angle !< Solid angle of STL polyhedra projected on point unit sphere.
   integer(I4P)                          :: f           !< Counter.

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

   elemental subroutine initialize(self, aabb_refinement_levels)
   !< Initialize file.
   class(surface_stl_object), intent(inout)        :: self                   !< File STL.
   integer(I4P),              intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.

   call self%destroy
   if (present(aabb_refinement_levels)) self%aabb%refinement_levels = aabb_refinement_levels
   endsubroutine initialize

   pure subroutine merge_solids(self, other)
   !< Merge facets with ones of other STL file.
   class(surface_stl_object), intent(inout) :: self     !< File STL.
   type(surface_stl_object),  intent(in)    :: other    !< Other file STL.
   type(facet_object), allocatable          :: facet(:) !< Facets temporary list.
   integer(I4P)                             :: f        !< Counter.

   if (other%facets_number > 0) then
      if (self%facets_number > 0) then
         allocate(facet(1:self%facets_number + other%facets_number))
         do f=1, self%facets_number
            facet(f)  =  self%facet(f)
         enddo
         do f=1, other%facets_number
            facet(self%facets_number+f) = other%facet(f)
         enddo
         call move_alloc(from=facet, to=self%facet)
         self%facets_number = self%facets_number + other%facets_number
      else
         allocate(self%facet(1:other%facets_number))
         do f=1, other%facets_number
            self%facet(f) = other%facet(f)
         enddo
         self%facets_number = other%facets_number
      endif
      call self%analize(aabb_refinement_levels=self%aabb%refinement_levels)
   endif
   endsubroutine merge_solids

   elemental subroutine resize(self, x, y, z, factor, respect_centroid, recompute_metrix)
   !< Resize (scale) facets by x or y or z or vectorial factors.
   !<
   !< @note The name `scale` has not been used, it been a Fortran built-in.
   !<
   !< @note If centroid must be used for center of resize it must be already computed.
   class(surface_stl_object), intent(inout)        :: self              !< File STL.
   real(R8P),                 intent(in), optional :: x                 !< Factor along x axis.
   real(R8P),                 intent(in), optional :: y                 !< Factor along y axis.
   real(R8P),                 intent(in), optional :: z                 !< Factor along z axis.
   type(vector_R8P),          intent(in), optional :: factor            !< Vectorial factor.
   logical,                   intent(in), optional :: respect_centroid  !< Sentinel to activate centroid as resize center.
   logical,                   intent(in), optional :: recompute_metrix  !< Sentinel to activate metrix recomputation.
   type(vector_R8P)                                :: factor_           !< Vectorial factor, local variable.
   logical                                         :: respect_centroid_ !< Sentinel to activate centroid as resize center, local v.

   respect_centroid_ = .false. ; if (present(respect_centroid)) respect_centroid_ = respect_centroid
   if (self%facets_number>0) then
      factor_ = 1._R8P
      if (present(factor)) then
         factor_ = factor
      else
         if (present(x)) factor_%x = x
         if (present(y)) factor_%y = y
         if (present(z)) factor_%z = z
      endif
      if (respect_centroid_) then
         call self%facet%resize(factor=factor_, center=self%centroid)
      else
         call self%facet%resize(factor=factor_, center=0 * ex_R8P)
      endif
      if (present(recompute_metrix)) then
         if (recompute_metrix) call self%compute_metrix
      endif
   endif
   endsubroutine resize

   elemental subroutine reverse_normals(self)
   !< Reverse facets normals.
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%reverse_normal
   endsubroutine reverse_normals

   pure subroutine sanitize(self, do_analysis)
   !< Sanitize STL.
   class(surface_stl_object), intent(inout)        :: self        !< File STL.
   logical,                   intent(in), optional :: do_analysis !< Sentil for performing a first analysis.

   if (self%facets_number>0) then
      if (present(do_analysis)) then
         if (do_analysis) call self%analize(aabb_refinement_levels=self%aabb%refinement_levels)
      endif
      if (self%facet_1_de%ids_number>0.or.&
          self%facet_2_de%ids_number>0.or.&
          self%facet_3_de%ids_number>0) call self%connect_nearby_vertices
      call self%analize(aabb_refinement_levels=self%aabb%refinement_levels)
      call self%sanitize_normals
   endif
   endsubroutine sanitize

   pure subroutine sanitize_normals(self)
   !< Sanitize facets normals, make them consistent.
   !<
   !< @note Facets connectivity and normals must be already computed.
   class(surface_stl_object), intent(inout) :: self             !< File STL.
   logical, allocatable                     :: facet_checked(:) !< List of facets checked.
   integer(I4P)                             :: f, ff            !< Counter.

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

   pure function smallest_edge_len(self) result(smallest)
   !< Return the smallest edge length.
   class(surface_stl_object), intent(in) :: self     !< File STL.
   real(R8P)                             :: smallest !< Smallest edge length.
   integer(I4P)                          :: f        !< Counter.

   smallest = MaxR8P
   if (self%facets_number>0) then
      do f=1, self%facets_number
         smallest = min(smallest, self%facet(f)%smallest_edge_len())
      enddo
   endif
   endfunction smallest_edge_len

   pure function statistics(self, prefix) result(stats)
   !< Return STL statistics.
   class(surface_stl_object), intent(in)           :: self             !< File STL.
   character(*),              intent(in), optional :: prefix           !< Lines prefix.
   character(len=:), allocatable                   :: stats            !< STL statistics.
   character(len=:), allocatable                   :: prefix_          !< Lines prefix, local variable.
   character(1), parameter                         :: NL=new_line('a') !< Line terminator.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   stats = ''
   if (self%facets_number > 0) then
      stats=stats//prefix_//'X extents: ['//trim(str(self%bmin%x))//', '//trim(str(self%bmax%x))//']'//NL
      stats=stats//prefix_//'Y extents: ['//trim(str(self%bmin%y))//', '//trim(str(self%bmax%y))//']'//NL
      stats=stats//prefix_//'Z extents: ['//trim(str(self%bmin%z))//', '//trim(str(self%bmax%z))//']'//NL
      stats=stats//prefix_//'volume: '//trim(str(self%volume))//NL
      stats=stats//prefix_//'centroid: ['//trim(str(self%centroid%x))//', '//&
                                           trim(str(self%centroid%y))//', '//&
                                           trim(str(self%centroid%z))//']'//NL
      stats=stats//prefix_//'number of facets: '//trim(str(self%facets_number))//NL
      stats=stats//prefix_//'number of facets with 1 edges disconnected: '//trim(str(self%facet_1_de%ids_number))//NL
      stats=stats//prefix_//'number of facets with 2 edges disconnected: '//trim(str(self%facet_2_de%ids_number))//NL
      stats=stats//prefix_//'number of facets with 3 edges disconnected: '//trim(str(self%facet_3_de%ids_number))//NL
      stats=stats//prefix_//'number of AABB refinement levels: '//trim(str(self%aabb%refinement_levels))!//NL
   endif
   endfunction statistics

   elemental subroutine translate(self, x, y, z, delta, recompute_metrix)
   !< Translate facets x or y or z or vectorial delta increments.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in), optional :: x                !< Increment along x axis.
   real(R8P),                 intent(in), optional :: y                !< Increment along y axis.
   real(R8P),                 intent(in), optional :: z                !< Increment along z axis.
   type(vector_R8P),          intent(in), optional :: delta            !< Vectorial increment.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   type(vector_R8P)                                :: delta_           !< Vectorial increment, local variable.

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
      if (present(recompute_metrix)) then
         if (recompute_metrix) call self%compute_metrix
      endif
   endif
   endsubroutine translate

   ! operators
   ! =
   pure subroutine surface_stl_assign_surface_stl(lhs, rhs)
   !< Operator `=`.
   class(surface_stl_object), intent(inout) :: lhs !< Left hand side.
   type(surface_stl_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                             :: f   !< Counter.

   lhs%facets_number = rhs%facets_number
   if (allocated(lhs%facet)) then
      call lhs%facet%destroy
      deallocate(lhs%facet)
   endif
   if (allocated(rhs%facet)) then
      allocate(lhs%facet(1:lhs%facets_number))
      do f=1, lhs%facets_number
         lhs%facet(f) = rhs%facet(f)
      enddo
   endif
   lhs%facet_1_de = rhs%facet_1_de
   lhs%facet_2_de = rhs%facet_2_de
   lhs%facet_3_de = rhs%facet_3_de
   lhs%aabb = rhs%aabb
   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
   lhs%volume = rhs%volume
   lhs%centroid = rhs%centroid
   endsubroutine surface_stl_assign_surface_stl

   ! private methods
   pure subroutine compute_facets_disconnected(self)
   !< Compute facets with disconnected edges.
   class(surface_stl_object), intent(inout) :: self  !< File STL.
   logical                                  :: de(3) !< Flag to check edges disconnection.
   integer(I4P)                             :: f     !< Counter.

   call self%facet_1_de%destroy
   call self%facet_2_de%destroy
   call self%facet_3_de%destroy
   if (self%facets_number>0) then
      do f=1, self%facets_number
         de = .false.
         if (self%facet(f)%fcon_edge_12==0_I4P) de(1) = .true.
         if (self%facet(f)%fcon_edge_23==0_I4P) de(2) = .true.
         if (self%facet(f)%fcon_edge_31==0_I4P) de(3) = .true.
         select case(count(de))
         case(1_I4P)
            call self%facet_1_de%put(id=self%facet(f)%id)
         case(2_I4P)
            call self%facet_2_de%put(id=self%facet(f)%id)
         case(3_I4P)
            call self%facet_3_de%put(id=self%facet(f)%id)
         endselect
      enddo
   endif
   endsubroutine compute_facets_disconnected

   elemental subroutine mirror_by_normal(self, normal, recompute_metrix)
   !< Mirror facets given normal of mirroring plane.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),          intent(in)           :: normal           !< Normal of mirroring plane.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                       :: matrix(3,3)      !< Mirroring matrix.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = mirror_matrix_R8P(normal=normal)
      do f=1, self%facets_number
         call self%facet(f)%mirror(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine mirror_by_normal

   pure subroutine mirror_by_matrix(self, matrix, recompute_metrix)
   !< Mirror facet given matrix (of mirroring).
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in)           :: matrix(3,3)      !< Mirroring matrix.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                    :: f                !< Counter.

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
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),          intent(in)           :: axis             !< Axis of rotation.
   real(R8P),                 intent(in)           :: angle            !< Angle of rotation.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                       :: matrix(3,3)      !< Rotation matrix.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = rotation_matrix_R8P(axis=axis, angle=angle)
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_axis_angle

   pure subroutine rotate_by_matrix(self, matrix, recompute_metrix)
   !< Rotate facet given matrix (of ratation).
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in)           :: matrix(3,3)      !< Rotation matrix.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_matrix

   elemental subroutine set_facets_id(self)
   !< (Re)set facets ID.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         self%facet(f)%id = f
      enddo
   endif
   endsubroutine set_facets_id
endmodule fossil_surface_stl_object
