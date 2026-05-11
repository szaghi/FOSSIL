!< FOSSIL, STL surface class definition.

module fossil_surface_stl_object
!< FOSSIL, STL surface class definition.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use fossil_utils, only : EPS, FRLEN, PI, is_inside_bb
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: surface_stl_object
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE
public :: sign_algorithm_from_string

! Point-in-polyhedron algorithm selector for `is_point_inside`, `compute_distance`,
! and `distance` (when `is_signed=.true.`):
!   SIGN_RAY_INTERSECTIONS — count intersections of an axis-aligned ray with the
!                            polyhedron; odd = inside, even = outside.
!   SIGN_SOLID_ANGLE       — sum the projected solid angles of all facets;
!                            ~±4π = inside, ~0 = outside.
integer(I4P), parameter :: SIGN_RAY_INTERSECTIONS = 1_I4P
integer(I4P), parameter :: SIGN_SOLID_ANGLE       = 2_I4P

type :: surface_stl_object
   !< FOSSIL STL surface class.
   !<
   !< All components are `private`. External code interacts via accessors:
   !<  - scalar/value getters: `get_facets_number`, `get_bmin`, `get_bmax`, `get_volume`,
   !<    `get_centroid`, `get_header`, `set_header`
   !<  - facet access: `facet_at(i)` returns a `pointer` to a single facet (or `null()`
   !<    if `i` is out of range), `facets_ref()` returns a `pointer` to the whole array.
   !<    Both are intended as read-only views — the library does not enforce this through
   !<    the pointer (Fortran has no `const`-equivalent on components), so by convention
   !<    callers do not write through these accessors. Mutation goes via the surface's
   !<    own TBPs (`translate`, `rotate`, `mirror`, etc.).
   !<  - ownership transfer: `adopt_facets(arr)` moves an allocatable array into the
   !<    surface via `move_alloc`, then auto-runs `analize`. Used internally by
   !<    `load_from_file`.
   !<  - `aabb` is technically public so that callers can invoke its TBPs
   !<    (e.g. `surface%aabb%set_use_index(...)`); its own components are private.
   integer(I4P),                    private :: facets_number=0 !< Facets number (== size(facet)).
   type(facet_object), allocatable, private :: facet(:)        !< Facets.
   type(list_id_object),            private :: facet_1_de      !< Facets with one disconnected edges.
   type(list_id_object),            private :: facet_2_de      !< Facets with two disconnected edges.
   type(list_id_object),            private :: facet_3_de      !< Facets with three disconnected edges.
   type(aabb_tree_object)                   :: aabb            !< AABB tree handle (its own state is private).
   type(vector_R8P),                private :: bmin            !< Bounding-box min.
   type(vector_R8P),                private :: bmax            !< Bounding-box max.
   real(R8P),                       private :: volume=0._R8P   !< Volume bounded by STL surface.
   type(vector_R8P),                private :: centroid        !< Centroid of STL surface.
   character(FRLEN),                private :: header=''       !< STL file header (preserved across load/save).
   contains
      ! read-only accessors (pure, inlined at -O2, zero data copy for scalars)
      procedure, pass(self) :: get_facets_number !< Return facets_number.
      procedure, pass(self) :: get_bmin          !< Return bmin (bounding-box minimum).
      procedure, pass(self) :: get_bmax          !< Return bmax (bounding-box maximum).
      procedure, pass(self) :: get_volume        !< Return volume.
      procedure, pass(self) :: get_centroid      !< Return centroid.
      procedure, pass(self) :: get_header        !< Return STL header.
      ! mutator (the only externally-permitted direct writes; aabb has its own mutators)
      procedure, pass(self) :: set_header        !< Set STL header.
      ! facet access (pointer-returning; treat as read-only views)
      procedure, pass(self) :: facet_at          !< Return pointer to facet(i), null() if i is out of range.
      procedure, pass(self) :: facets_ref        !< Return pointer to the whole facet array.
      ! ownership transfer
      procedure, pass(self) :: adopt_facets      !< Take ownership of an allocatable facet array via move_alloc.
      ! file I/O (moved from file_stl_object)
      procedure, pass(self) :: load_from_file       !< Load STL from a file path (ASCII or binary, auto-detected).
      procedure, pass(self) :: save_into_file       !< Save STL to a file path.
      procedure, pass(self) :: save_aabb_into_file  !< Save the AABB tree leaves as separate STL files.
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
      procedure, pass(self) :: is_point_inside                 !< Determinate if point is inside or not STL.
      procedure, pass(self) :: is_point_inside_polyhedron_ri   !< Determinate if point is inside or not STL facets by ray intersect.
      procedure, pass(self) :: is_point_inside_polyhedron_sa   !< Determinate if point is inside or not STL facets by solid angle.
      procedure, pass(self) :: largest_edge_len                !< Return the largest edge length.
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
      ! finaliser — releases facet(:) and nested allocatables (aabb%node, facet_*_de%id)
      ! when the instance goes out of scope or is wrapped in an array container.
      final :: surface_stl_finalize
      ! private methods
      procedure, pass(self), private :: compute_facets_disconnected !< Compute facets with disconnected edges.
      procedure, pass(self), private :: mirror_by_normal            !< Mirror facets given normal of mirroring plane.
      procedure, pass(self), private :: mirror_by_matrix            !< Mirror facets given matrix.
      procedure, pass(self), private :: rotate_by_axis_angle        !< Rotate facets given axis and angle.
      procedure, pass(self), private :: rotate_by_matrix            !< Rotate facets given matrix.
      procedure, pass(self), private :: set_facets_id               !< (Re)set facets ID.
endtype surface_stl_object

contains
   ! accessors (pure, scalar/vector return — inlined by gfortran/ifort at -O2, no data copy)
   pure function get_facets_number(self) result(n)
   !< Return facets_number.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Facets number.

   n = self%facets_number
   endfunction get_facets_number

   pure function get_bmin(self) result(b)
   !< Return bmin (bounding-box minimum).
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: b    !< Bounding-box minimum.

   b = self%bmin
   endfunction get_bmin

   pure function get_bmax(self) result(b)
   !< Return bmax (bounding-box maximum).
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: b    !< Bounding-box maximum.

   b = self%bmax
   endfunction get_bmax

   pure function get_volume(self) result(v)
   !< Return volume bounded by STL surface.
   class(surface_stl_object), intent(in) :: self !< File STL.
   real(R8P)                             :: v    !< Volume.

   v = self%volume
   endfunction get_volume

   pure function get_centroid(self) result(c)
   !< Return centroid of STL surface.
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: c    !< Centroid.

   c = self%centroid
   endfunction get_centroid

   pure function get_header(self) result(h)
   !< Return STL header (the 80-char "solid <name>" string).
   class(surface_stl_object), intent(in) :: self !< File STL.
   character(FRLEN)                      :: h    !< Header.

   h = self%header
   endfunction get_header

   pure subroutine set_header(self, header)
   !< Set STL header (truncated/padded to FRLEN).
   class(surface_stl_object), intent(inout) :: self   !< File STL.
   character(*),              intent(in)    :: header !< New header text.

   self%header = header
   endsubroutine set_header

   ! facet access (pointer-returning, zero-copy; treat the returned pointer as a read-only view)

   function facet_at(self, i) result(p)
   !< Return a pointer to facet `i`, or `null()` if `i` is out of range.
   !<
   !< The returned pointer aliases the surface's internal storage — do not write through
   !< it; use the surface's TBPs (`translate`, `rotate`, etc.) for mutation. Bounds-test
   !< at the call site: `p => surface%facet_at(i); if (.not. associated(p)) ...`.
   class(surface_stl_object), intent(in), target :: self !< File STL.
   integer(I4P),              intent(in)         :: i    !< Facet index (1..facets_number).
   type(facet_object),        pointer            :: p    !< Pointer to facet, or null() if out of range.

   if (allocated(self%facet) .and. i >= 1 .and. i <= self%facets_number) then
      p => self%facet(i)
   else
      p => null()
   endif
   endfunction facet_at

   function facets_ref(self) result(p)
   !< Return a pointer to the whole facet array (length = facets_number).
   !<
   !< Returns `null()` if no facets are allocated. Same read-only-view convention as
   !< `facet_at`: do not mutate through this pointer.
   class(surface_stl_object), intent(in), target :: self !< File STL.
   type(facet_object),        pointer            :: p(:) !< Pointer to facet array, or null() if unallocated.

   if (allocated(self%facet)) then
      p => self%facet
   else
      p => null()
   endif
   endfunction facets_ref

   ! ownership transfer

   subroutine adopt_facets(self, facets, aabb_refinement_levels)
   !< Take ownership of an allocatable facet array via `move_alloc`, then `analize`.
   !<
   !< The caller's `facets(:)` becomes unallocated on return — this is a zero-copy
   !< handoff. Used internally by `load_from_file` and available to external code that
   !< builds facets procedurally.
   class(surface_stl_object),       intent(inout)        :: self                   !< File STL.
   type(facet_object), allocatable, intent(inout)        :: facets(:)              !< Facets to adopt.
   integer(I4P),                    intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.

   if (allocated(self%facet)) deallocate(self%facet)
   if (allocated(facets)) then
      call move_alloc(from=facets, to=self%facet)
      self%facets_number = size(self%facet, dim=1)
   else
      self%facets_number = 0
   endif
   call self%analize(aabb_refinement_levels=aabb_refinement_levels)
   endsubroutine adopt_facets

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

   ! elemental subroutine analize(self, aabb_refinement_levels)
   subroutine analize(self, aabb_refinement_levels)
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
      call self%aabb%initialize(refinement_levels=aabb_refinement_levels, facet=self%facet,largest_edge_len=self%largest_edge_len())
      call self%build_connectivity
      call self%compute_facets_disconnected
      call self%compute_volume
      call self%compute_centroid
   endif
   endsubroutine analize

   ! pure subroutine build_connectivity(self)
   subroutine build_connectivity(self)
   !< Build facets connectivity.
   class(surface_stl_object), intent(inout) :: self              !< File STL.
   real(R8P)                                :: smallest_edge_len !< Smallest edge length.
   integer(I4P)                             :: f1, f2            !< Counter.
   type(aabb_tree_object)                   :: aabb              !< Temporary AABB tree.

   if (self%facets_number>0) then
      call self%facet%destroy_connectivity
      smallest_edge_len = self%smallest_edge_len() * 0.9_R8P
      if (self%aabb%get_is_initialized()) then
         ! exploit AABB structure
         call aabb%initialize(facet=self%facet, refinement_levels=self%aabb%get_refinement_levels(), do_facets_distribute=.false.)
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
         call self%analize(aabb_refinement_levels=self%aabb%get_refinement_levels())
         if (present(remainder)) call remainder%analize(aabb_refinement_levels=self%aabb%get_refinement_levels())
      endif
   endif
   endsubroutine clip

   pure subroutine compute_centroid(self)
   !< Compute centroid of STL surface.
   !<
   !< @note Metrix and volume must be already computed.
   !<
   !< Degenerate / non-watertight surfaces can have `volume == 0`; the centroid
   !< formula divides by `48 * volume`, so we guard against producing NaN/Inf.
   !< In the degenerate case we fall back to the origin — a defined value the
   !< caller can detect via `volume == 0` if needed.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      self%centroid = 0._R8P
      if (abs(self%volume) > tiny(0._R8P)) then
         do f=1, self%facets_number
            self%centroid = self%centroid - self%facet(f)%centroid_part()
         enddo
         self%centroid = self%centroid / (48 * self%volume)
      endif
   endif
   endsubroutine compute_centroid

   subroutine compute_distance(self, point, distance, is_signed, sign_algorithm, is_square_root, &
                               facet_index, edge_index, vertex_index)
   !< Compute the (minimum) distance returning also the closest point.
   class(surface_stl_object), intent(in)            :: self            !< File STL.
   type(vector_R8P),          intent(in)            :: point           !< Point coordinates.
   real(R8P),                 intent(out)           :: distance        !< Minimum distance.
   logical,                   intent(in),  optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in),  optional :: sign_algorithm  !< Algorithm code (SIGN_RAY_INTERSECTIONS or SIGN_SOLID_ANGLE).
   logical,                   intent(in),  optional :: is_square_root  !< Sentinel to trigger square-root distance.
   integer(I4P),              intent(out), optional :: facet_index     !< Index of facet containing the closest point.
   integer(I4P),              intent(out), optional :: edge_index      !< Index of edge on facet containing the closest point.
   integer(I4P),              intent(out), optional :: vertex_index    !< Index of vertex on facet containing the closest point.
   real(R8P)                                        :: distance_       !< Minimum distance, temporary buffer.
   integer(I4P)                                     :: facet_index_    !< Index of facet containing the closest point, local var.
   integer(I4P)                                     :: f               !< Counter.

   ! initialize intent(out) outputs so they are defined on every path:
   ! - empty surface (facets_number == 0) — distance stays MaxR8P, facet_index_ stays 0
   ! - AABB path — facet_index_ is not (yet) populated by distance_tree; keep 0 as sentinel
   distance     = MaxR8P
   facet_index_ = 0
   if (self%facets_number > 0) then
      if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
         ! exploit AABB refinement levels
         ! distance = self%aabb%distance(facet=self%facet, point=point)
         distance = self%aabb%distance_tree(facet=self%facet, point=point)
      else
         ! brute-force search over all facets
         do f=1, self%facets_number
            call self%facet(f)%compute_distance(point=point, distance=distance_)
            if (abs(distance_) <= abs(distance)) then
               facet_index_ = self%facet(f)%id
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
        if (self%is_point_inside(point=point, sign_algorithm=sign_algorithm)) distance = -distance
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

   function distance(self, point, is_signed, sign_algorithm, is_square_root)
   !< Return the (minimum) distance from a point to the triangulated surface.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical,                   intent(in), optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_RAY_INTERSECTIONS or SIGN_SOLID_ANGLE).
   logical,                   intent(in), optional :: is_square_root  !< Sentinel to trigger square-root distance.
   real(R8P)                                       :: distance        !< Minimum distance from point to the triangulated surface.

   call self%compute_distance(point=point, distance=distance, &
                              is_signed=is_signed, sign_algorithm=sign_algorithm, is_square_root=is_square_root)
   endfunction distance

   function is_point_inside(self, point, sign_algorithm) result(is_inside)
   !< Compute sign.
   !<
   !< `sign_algorithm` selects the point-in-polyhedron test; pass `SIGN_RAY_INTERSECTIONS`
   !< (default) or `SIGN_SOLID_ANGLE`. Out-of-range values raise `error stop`.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical                                         :: is_inside       !< Return true if point is inside surface.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_RAY_INTERSECTIONS or SIGN_SOLID_ANGLE).
   integer(I4P)                                    :: algo            !< Effective algorithm code.

   algo = SIGN_RAY_INTERSECTIONS ; if (present(sign_algorithm)) algo = sign_algorithm
   select case(algo)
   case(SIGN_SOLID_ANGLE)
      is_inside = self%is_point_inside_polyhedron_sa(point=point)
   case(SIGN_RAY_INTERSECTIONS)
      is_inside = self%is_point_inside_polyhedron_ri(point=point)
   case default
      error stop 'fossil_surface_stl_object%is_point_inside: unknown sign_algorithm code &
                 &(valid: SIGN_RAY_INTERSECTIONS=1, SIGN_SOLID_ANGLE=2)'
   endselect
   endfunction is_point_inside

   function is_point_inside_polyhedron_ri(self, point) result(is_inside)
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
      function is_inside_by_ray_intersect(ray_origin, ray_direction) result(is_inside_by)
      !< Generic line intersect test.
      type(vector_R8P), intent(in) :: ray_origin           !< Ray origin.
      type(vector_R8P), intent(in) :: ray_direction        !< Ray direction.
      integer(I4P)                 :: intersections_number !< Ray intersections number of STL polyhedra with respect point.
      integer(I4P)                 :: f                    !< Counter.
      logical                      :: is_inside_by         !< Test result.

      intersections_number = 0

      if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
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
   if (present(aabb_refinement_levels)) call self%aabb%set_refinement_levels(aabb_refinement_levels)
   endsubroutine initialize

   pure function largest_edge_len(self) result(largest)
   !< Return the largest edge length.
   class(surface_stl_object), intent(in) :: self    !< File STL.
   real(R8P)                             :: largest !< Largest edge length.
   integer(I4P)                          :: f       !< Counter.

   largest = 0._R8P
   if (self%facets_number>0) then
      do f=1, self%facets_number
         largest = max(largest, self%facet(f)%largest_edge_len())
      enddo
   endif
   endfunction largest_edge_len

   ! pure subroutine merge_solids(self, other)
   subroutine merge_solids(self, other)
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
      call self%analize(aabb_refinement_levels=self%aabb%get_refinement_levels())
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

   ! pure subroutine sanitize(self, do_analysis)
   subroutine sanitize(self, do_analysis)
   !< Sanitize STL.
   class(surface_stl_object), intent(inout)        :: self        !< File STL.
   logical,                   intent(in), optional :: do_analysis !< Sentil for performing a first analysis.

   if (self%facets_number>0) then
      if (present(do_analysis)) then
         if (do_analysis) call self%analize(aabb_refinement_levels=self%aabb%get_refinement_levels())
      endif
      if (self%facet_1_de%ids_number>0.or.&
          self%facet_2_de%ids_number>0.or.&
          self%facet_3_de%ids_number>0) call self%connect_nearby_vertices
      call self%analize(aabb_refinement_levels=self%aabb%get_refinement_levels())
      call self%sanitize_normals
   endif
   endsubroutine sanitize

   pure subroutine sanitize_normals(self)
   !< Sanitize facets normals, make them consistent.
   !<
   !< @note Facets connectivity and normals must be already computed.
   class(surface_stl_object), intent(inout) :: self             !< File STL.
   logical, allocatable                     :: facet_checked(:) !< List of facets checked.
   integer(I4P)                             :: f, ff, e, neighbour !< Counters.

   if (self%facets_number>0) then
      allocate(facet_checked(1:self%facets_number))
      facet_checked = .false.
      f = 1
      facet_checked(f) = .true.
      do
         ff = 0
         do e=1, 3
            neighbour = self%facet(f)%fcon_edge(e)
            if (neighbour <= 0) cycle
            if (facet_checked(neighbour)) cycle
            call self%facet(f)%make_normal_consistent(edge=e, other=self%facet(neighbour))
            facet_checked(neighbour) = .true.
            ff = neighbour
         enddo
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
      stats=stats//prefix_//'number of AABB refinement levels: '//trim(str(self%aabb%get_refinement_levels()))!//NL
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

   ! finaliser
   subroutine surface_stl_finalize(self)
   !< Release facet(:) and reset state.
   !<
   !< Nested components (aabb%node, facet_*_de%id, per-facet allocatables) have their
   !< own finalisers / `=` semantics; releasing facet(:) here cascades correctly.
   type(surface_stl_object), intent(inout) :: self !< Surface.

   if (allocated(self%facet)) deallocate(self%facet)
   self%facets_number = 0
   endsubroutine surface_stl_finalize

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
         de = (self%facet(f)%fcon_edge == 0_I4P)
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

   ! non-TBP helpers

   pure function sign_algorithm_from_string(name) result(code)
   !< Translate a human-readable algorithm name into a `SIGN_*` code.
   !<
   !< Intended for CLI / config-file boundaries where users type a string and the library
   !< wants the integer dispatch code. Unknown names raise `error stop` with the valid set.
   !< Comparison is case-sensitive — keep names in canonical form (`'ray_intersections'`,
   !< `'solid_angle'`) at call sites.
   character(*), intent(in) :: name !< Algorithm name string.
   integer(I4P)             :: code !< Corresponding SIGN_* code.

   select case(trim(adjustl(name)))
   case('ray_intersections')
      code = SIGN_RAY_INTERSECTIONS
   case('solid_angle')
      code = SIGN_SOLID_ANGLE
   case default
      error stop 'fossil_surface_stl_object%sign_algorithm_from_string: unknown name "'//trim(adjustl(name))// &
                 '" (valid: "ray_intersections", "solid_angle")'
   endselect
   endfunction sign_algorithm_from_string

   ! file I/O (migrated from the deleted fossil_file_stl_object module)
   !
   ! These TBPs and their internal helpers replace what used to live on `file_stl_object`.
   ! Callers now do `call surface%load_from_file('foo.stl')` directly — there is no
   ! intermediate file handle to manage. Internally the routines use a local `file_unit`
   ! variable; format-mode (ASCII vs binary) is a local flag, not stored state.

   subroutine load_from_file(self, file_name, is_ascii, guess_format, clip_min, clip_max, aabb_refinement_levels)
   !< Load an STL file into the surface.
   !<
   !< Builds a local facet array, then transfers ownership via `adopt_facets` (which
   !< runs `analize`). Auto-detects ASCII vs binary when `guess_format=.true.` (size
   !< identity-check; see the binary-header trap discussion in audit #14 S3).
   !< When `clip_min`/`clip_max` are present, only facets entirely inside the AABB
   !< are loaded.
   class(surface_stl_object),       intent(inout)        :: self                   !< Surface STL.
   character(*),                    intent(in)           :: file_name              !< STL file path.
   logical,                         intent(in), optional :: is_ascii               !< Force ASCII (default .true. if guess_format=.false.).
   logical,                         intent(in), optional :: guess_format           !< Auto-detect format from file size.
   type(vector_R8P),                intent(in), optional :: clip_min, clip_max     !< AABB clip extents (facets inside only).
   integer(I4P),                    intent(in), optional :: aabb_refinement_levels !< AABB refinement levels passed to analize.
   type(facet_object), allocatable                       :: facets(:)              !< Local buffer for ownership transfer.
   integer(I4P)                                          :: file_unit              !< File unit.
   logical                                               :: is_ascii_              !< Effective ASCII flag.
   integer(I4P)                                          :: facets_number          !< Facet count from header.
   type(facet_object)                                    :: facet_clip             !< Buffer for clipped loading.
   integer(I4P)                                          :: f, ff                  !< Counters.

   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_read(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_, guess_format=guess_format)
   call stl_load_facets_number(file_unit=file_unit, is_ascii=is_ascii_, facets_number=facets_number)
   call stl_load_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header)
   if (present(clip_min).and.present(clip_max)) then
      ! two-pass: count, then re-read and keep facets inside the AABB
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call facet_clip%load_from_file_ascii(file_unit=file_unit)
         else
            call facet_clip%load_from_file_binary(file_unit=file_unit)
         endif
         if (is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(1)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(2)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(3))) ff = ff + 1
      enddo
      call stl_load_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header)
      allocate(facets(1:ff))
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call facet_clip%load_from_file_ascii(file_unit=file_unit)
         else
            call facet_clip%load_from_file_binary(file_unit=file_unit)
         endif
         if (is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(1)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(2)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(3))) then
            ff = ff + 1
            facets(ff) = facet_clip
            facets(ff)%id = ff
         endif
      enddo
   else
      allocate(facets(1:facets_number))
      do f=1, facets_number
         if (is_ascii_) then
            call facets(f)%load_from_file_ascii(file_unit=file_unit)
         else
            call facets(f)%load_from_file_binary(file_unit=file_unit)
         endif
         facets(f)%id = f
      enddo
   endif
   close(file_unit)
   call self%adopt_facets(facets=facets, aabb_refinement_levels=aabb_refinement_levels)
   endsubroutine load_from_file

   subroutine save_into_file(self, file_name, is_ascii)
   !< Save the surface to an STL file.
   class(surface_stl_object), intent(in)           :: self      !< Surface STL.
   character(*),              intent(in)           :: file_name !< STL file path.
   logical,                   intent(in), optional :: is_ascii  !< Write as ASCII (default .true.).
   integer(I4P)                                    :: file_unit !< File unit.
   logical                                         :: is_ascii_ !< Effective ASCII flag.
   integer(I4P)                                    :: f         !< Counter.

   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_write(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_)
   call stl_save_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header, facets_number=self%facets_number)
   if (is_ascii_) then
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_ascii(file_unit=file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_binary(file_unit=file_unit)
      enddo
   endif
   if (is_ascii_) write(file_unit, '(A)') 'endsolid '//trim(self%header)
   close(file_unit)
   endsubroutine save_into_file

   subroutine save_aabb_into_file(self, base_file_name, is_ascii)
   !< Save each AABB tree leaf as a separate STL file.
   !<
   !< File name pattern: `<base_file_name>aabb-l_<level>-b_<box>.stl`.
   class(surface_stl_object), intent(in)           :: self           !< Surface STL.
   character(*),              intent(in)           :: base_file_name !< Base file-name prefix.
   logical,                   intent(in), optional :: is_ascii       !< Write as ASCII (default .true.).
   type(facet_object), allocatable                 :: aabb_facet(:)  !< AABB leaf facets.
   type(surface_stl_object)                        :: leaf           !< Scratch surface for save.
   integer(I4P)                                    :: b, l           !< Counters.
   character(len=:), allocatable                   :: file_name      !< Composed file name.

   if (.not. self%aabb%get_is_initialized()) return
   do while (self%aabb%loop_node(facet=self%facet, aabb_facet=aabb_facet, b=b, l=l))
      file_name = trim(adjustl(base_file_name))//'aabb-l_'//trim(str(l, .true.))//'-b_'//trim(str(b, .true.))//'.stl'
      call leaf%adopt_facets(facets=aabb_facet)
      call leaf%set_header(self%header)
      call leaf%save_into_file(file_name=file_name, is_ascii=is_ascii)
      call leaf%destroy
   enddo
   endsubroutine save_aabb_into_file

   ! internal I/O helpers (module-private; not exposed as TBPs)

   subroutine stl_open_for_read(file_name, file_unit, is_ascii, guess_format)
   !< Open an STL file for reading; auto-detect format via size identity if requested.
   character(*),  intent(in)              :: file_name    !< File path.
   integer(I4P),  intent(out)             :: file_unit    !< Newunit-assigned unit.
   logical,       intent(inout)           :: is_ascii     !< In/out: format flag (rewritten by guess).
   logical,       intent(in),   optional  :: guess_format !< Auto-detect via file-size identity.
   logical                                :: guess_format_, file_exist
   integer(I4P)                           :: file_size, facets_count, ios
   integer(I4P), parameter                :: BINARY_HEADER_BYTES = 80_I4P
   integer(I4P), parameter                :: BINARY_FACET_BYTES  = 50_I4P
   integer(I4P), parameter                :: BINARY_COUNT_BYTES  =  4_I4P

   guess_format_ = .false. ; if (present(guess_format)) guess_format_ = guess_format
   inquire(file=file_name, exist=file_exist, size=file_size)
   if (.not. file_exist) then
      write(stderr, '(A)') 'error: file "'//file_name//'" does not exist, impossible to open file!'
      error stop 'fossil_surface_stl_object%load_from_file: file not found'
   endif
   if (guess_format_) then
      is_ascii = .true.
      if (file_size > BINARY_HEADER_BYTES + BINARY_COUNT_BYTES) then
         open(newunit=file_unit, file=file_name, access='stream', form='unformatted', action='read', iostat=ios)
         if (ios == 0) then
            read(file_unit, pos=BINARY_HEADER_BYTES + 1, iostat=ios) facets_count
            close(file_unit)
            if (ios == 0 .and. facets_count > 0 .and. &
                file_size == BINARY_HEADER_BYTES + BINARY_COUNT_BYTES + facets_count * BINARY_FACET_BYTES) then
               is_ascii = .false.
            endif
         endif
      endif
   endif
   if (is_ascii) then
      open(newunit=file_unit, file=file_name,                  form='formatted',   action='read')
   else
      open(newunit=file_unit, file=file_name, access='stream', form='unformatted', action='read')
   endif
   endsubroutine stl_open_for_read

   subroutine stl_open_for_write(file_name, file_unit, is_ascii)
   !< Open an STL file for writing.
   character(*), intent(in)  :: file_name !< File path.
   integer(I4P), intent(out) :: file_unit !< Newunit-assigned unit.
   logical,      intent(in)  :: is_ascii  !< Format flag.

   if (is_ascii) then
      open(newunit=file_unit, file=file_name,                  form='formatted')
   else
      open(newunit=file_unit, file=file_name, access='stream', form='unformatted')
   endif
   endsubroutine stl_open_for_write

   subroutine stl_load_facets_number(file_unit, is_ascii, facets_number)
   !< Count facets in the file (ASCII: count "facet normal" lines; binary: read uint32).
   !< File is rewound on exit.
   integer(I4P), intent(in)  :: file_unit     !< Open file unit.
   logical,      intent(in)  :: is_ascii      !< Format flag.
   integer(I4P), intent(out) :: facets_number !< Count.
   character(FRLEN)          :: facet_record  !< Line buffer for ASCII scan.

   facets_number = 0
   rewind(file_unit)
   if (is_ascii) then
      do
         read(file_unit, '(A)', end=10, err=10) facet_record
         if (index(string=facet_record, substring='facet normal') > 0) facets_number = facets_number + 1
      enddo
   else
      read(file_unit, end=10, err=10) facet_record
      read(file_unit, end=10, err=10) facets_number
   endif
   10 rewind(file_unit)
   endsubroutine stl_load_facets_number

   subroutine stl_load_header(file_unit, is_ascii, header)
   !< Read the STL header into `header`. File is rewound before reading.
   integer(I4P),     intent(in)    :: file_unit !< Open file unit.
   logical,          intent(in)    :: is_ascii  !< Format flag.
   character(FRLEN), intent(out)   :: header    !< Header text.
   integer(I4P)                    :: facets_number_dummy

   rewind(file_unit)
   if (is_ascii) then
      read(file_unit, '(A)') header
      header = trim(adjustl(header(index(header, 'solid')+6:)))
   else
      read(file_unit) header
      read(file_unit) facets_number_dummy   ! count read elsewhere; consume here to advance the stream
   endif
   endsubroutine stl_load_header

   subroutine stl_save_header(file_unit, is_ascii, header, facets_number)
   !< Write the STL header. File is rewound before writing.
   integer(I4P),     intent(in) :: file_unit     !< Open file unit.
   logical,          intent(in) :: is_ascii      !< Format flag.
   character(FRLEN), intent(in) :: header        !< Header text.
   integer(I4P),     intent(in) :: facets_number !< Facet count (binary only).

   rewind(file_unit)
   if (is_ascii) then
      write(file_unit, '(A)') 'solid '//trim(header)
   else
      write(file_unit) header
      write(file_unit) facets_number
   endif
   endsubroutine stl_save_header
endmodule fossil_surface_stl_object
