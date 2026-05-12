!< FOSSIL, vertex pool object.
!<
!< Stage 1 of the indexed-vertex refactor (issue #5): a persisted, queryable
!< unique-vertex pool built from a facet array via EPS coincidence. In this
!< stage the pool is a **derived artifact** -- facets still own their inline
!< `vertex(3)` coordinates. Consumers can query the pool to check structural
!< coincidence without re-running EPS comparisons. Future stages will move
!< facets to carry only `vertex_id(3)` and drop the inline coordinates, at
!< which point the pool becomes the authoritative storage.
!<
!< API:
!<   call pool%initialize_from_facets(facet, status)
!<   n = pool%vertex_count()
!<   v = pool%coord(vid)
!<   vid = pool%facet_vid(facet_id, local_v)
!<   call pool%destroy
!<
!< Build: O(N^2) worst case in stage 1 (mirrors the union-find pass already in
!< compute_connectivity). Stage 5 will replace with a spatial hash.

module fossil_vertex_pool_object

use fossil_facet_object, only : facet_object
use fossil_utils,        only : EPS
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: vertex_pool_object

type :: vertex_pool_object
   !< Unique-vertex coordinate pool with (facet,local_v) -> pool_id mapping.
   !<
   !< Invariants after `initialize_from_facets`:
   !<   - `coord(1:n_vertices)` contains pairwise-distinct vertex coordinates
   !<     (pairwise distance > EPS on at least one component).
   !<   - `facet_to_pool(1:3, 1:n_facets)` maps every (local_v, facet_id) onto
   !<     a valid pool index in [1, n_vertices].
   !<   - For every (f, v): `coord(facet_to_pool(v, f))` and `facet(f)%vertex(v)`
   !<     are EPS-coincident.
   private
   type(vector_R8P), allocatable :: coord_(:)         !< Unique vertex coordinates.
   integer(I4P),     allocatable :: facet_to_pool(:,:) !< (local_v=1..3, facet_id=1..n_facets) -> pool id.
   integer(I4P)                  :: n_vertices = 0    !< Pool size.
   integer(I4P)                  :: n_facets   = 0    !< Source facet count.
   logical                       :: is_initialized = .false.
   contains
      procedure, pass(self) :: initialize_from_facets
      procedure, pass(self) :: vertex_count
      procedure, pass(self) :: facets_count
      procedure, pass(self) :: coord
      procedure, pass(self) :: facet_vid
      procedure, pass(self) :: get_is_initialized
      procedure, pass(self) :: destroy
      final :: finalize
endtype vertex_pool_object

contains

   subroutine initialize_from_facets(self, facet, status)
   !< Build the pool from a facet array via EPS coincidence (union-find).
   !<
   !< Algorithm:
   !<   1. Treat each (f, v) as a candidate vertex; initialize one union-find
   !<      class per (f, v) over global ids gid = (f-1)*3 + v.
   !<   2. Pairwise compare all (f1,v1) against (f2>=f1, v2) at strict EPS; union.
   !<   3. Collapse each gid to its root; assign a dense pool id per root.
   !<
   !< Complexity: O((3N)^2) pairwise tests. Acceptable for stage 1 (matches the
   !< complexity of the union-find pass already in compute_connectivity).
   class(vertex_pool_object),   intent(inout)         :: self      !< Vertex pool.
   type(facet_object),          intent(in)            :: facet(:)  !< Source facets.
   integer(I4P),                intent(out), optional :: status    !< 0=ok, 1=alloc fail.
   integer(I4P), allocatable                          :: parent(:), rank_(:) !< Union-find arrays.
   integer(I4P), allocatable                          :: root_to_pid(:)      !< Root gid -> dense pool id (0 if unassigned).
   integer(I4P)                                       :: nf, nv_total
   integer(I4P)                                       :: f1, f2, v1, v2
   integer(I4P)                                       :: gid1, gid2, root, pid
   integer(I4P)                                       :: istat

   call self%destroy
   if (present(status)) status = 0_I4P

   nf = size(facet)
   self%n_facets = nf
   if (nf <= 0) then
      self%is_initialized = .true.
      return
   endif

   nv_total = 3 * nf
   allocate(parent(nv_total), rank_(nv_total), stat=istat)
   if (istat /= 0) then
      if (present(status)) status = 1_I4P
      return
   endif
   do gid1 = 1, nv_total
      parent(gid1) = gid1
      rank_(gid1)  = 0_I4P
   enddo

   ! Pairwise EPS union. Inner loop walks f2 > f1 only, but for each (f1, v1) we
   ! also union with (f1, v2 > v1) -- a degenerate triangle can have coincident
   ! local vertices and we want to detect that (it would otherwise produce two
   ! pool entries for the same point).
   do f1 = 1, nf
      do v1 = 1, 3
         gid1 = (f1 - 1) * 3 + v1
         ! Within-facet coincidences (degenerate-triangle case).
         do v2 = v1 + 1, 3
            if (vertices_eps_match(facet(f1)%vertex(v1), facet(f1)%vertex(v2))) then
               call uf_union(gid1, (f1 - 1) * 3 + v2, parent, rank_)
            endif
         enddo
         ! Cross-facet coincidences.
         do f2 = f1 + 1, nf
            do v2 = 1, 3
               if (vertices_eps_match(facet(f1)%vertex(v1), facet(f2)%vertex(v2))) then
                  call uf_union(gid1, (f2 - 1) * 3 + v2, parent, rank_)
               endif
            enddo
         enddo
      enddo
   enddo

   ! Collapse to roots, assign dense pool ids in order of first appearance.
   allocate(root_to_pid(nv_total), self%facet_to_pool(3, nf), stat=istat)
   if (istat /= 0) then
      if (present(status)) status = 1_I4P
      deallocate(parent, rank_)
      return
   endif
   root_to_pid = 0_I4P
   pid = 0_I4P
   do f1 = 1, nf
      do v1 = 1, 3
         gid1 = (f1 - 1) * 3 + v1
         root = uf_find(gid1, parent)
         if (root_to_pid(root) == 0_I4P) then
            pid = pid + 1_I4P
            root_to_pid(root) = pid
         endif
         self%facet_to_pool(v1, f1) = root_to_pid(root)
      enddo
   enddo
   self%n_vertices = pid

   ! Materialize unique coordinates in the same dense order.
   allocate(self%coord_(self%n_vertices), stat=istat)
   if (istat /= 0) then
      if (present(status)) status = 1_I4P
      deallocate(parent, rank_, root_to_pid)
      return
   endif
   do f1 = 1, nf
      do v1 = 1, 3
         pid = self%facet_to_pool(v1, f1)
         self%coord_(pid) = facet(f1)%vertex(v1)
      enddo
   enddo

   deallocate(parent, rank_, root_to_pid)
   self%is_initialized = .true.
   endsubroutine initialize_from_facets

   pure function vertex_count(self) result(n)
   !< Return number of unique vertices in the pool.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P)                          :: n
   n = self%n_vertices
   endfunction vertex_count

   pure function facets_count(self) result(n)
   !< Return the source facet count the pool was built from.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P)                          :: n
   n = self%n_facets
   endfunction facets_count

   pure function coord(self, vid) result(v)
   !< Return the coordinate of pool vertex `vid`. Returns origin if out of range.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P),              intent(in) :: vid
   type(vector_R8P)                      :: v
   if (vid >= 1 .and. vid <= self%n_vertices) then
      v = self%coord_(vid)
   else
      v = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   endif
   endfunction coord

   pure function facet_vid(self, facet_id, local_v) result(vid)
   !< Return the pool vertex id for (facet_id, local_v). Returns 0 if out of range.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P),              intent(in) :: facet_id
   integer(I4P),              intent(in) :: local_v
   integer(I4P)                          :: vid
   if (facet_id >= 1 .and. facet_id <= self%n_facets .and. &
       local_v  >= 1 .and. local_v  <= 3) then
      vid = self%facet_to_pool(local_v, facet_id)
   else
      vid = 0_I4P
   endif
   endfunction facet_vid

   pure function get_is_initialized(self) result(yes)
   !< Has the pool been built?
   class(vertex_pool_object), intent(in) :: self
   logical                               :: yes
   yes = self%is_initialized
   endfunction get_is_initialized

   pure subroutine destroy(self)
   !< Release all storage and reset to uninitialized.
   class(vertex_pool_object), intent(inout) :: self
   if (allocated(self%coord_))        deallocate(self%coord_)
   if (allocated(self%facet_to_pool)) deallocate(self%facet_to_pool)
   self%n_vertices     = 0_I4P
   self%n_facets       = 0_I4P
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine finalize(self)
   !< Finaliser. `class` not allowed on `final`; takes `type`.
   type(vertex_pool_object), intent(inout) :: self
   if (allocated(self%coord_))        deallocate(self%coord_)
   if (allocated(self%facet_to_pool)) deallocate(self%facet_to_pool)
   endsubroutine finalize

   ! ---- internal helpers ----

   pure function vertices_eps_match(a, b) result(yes)
   !< Strict EPS coincidence -- same rule as `vertices_match` inside
   !< compute_connectivity. Kept private to the pool so the coincidence
   !< definition lives in one place.
   type(vector_R8P), intent(in) :: a, b
   logical                      :: yes
   yes = (abs(a%x - b%x) <= EPS) .and. &
         (abs(a%y - b%y) <= EPS) .and. &
         (abs(a%z - b%z) <= EPS)
   endfunction vertices_eps_match

   pure function uf_find(x, parent) result(root)
   !< Union-find root with path halving. Pure: parent is intent(in), so the
   !< usual path-compression write is replaced by re-traversal -- O(log* N)
   !< amortized in practice for our densities, acceptable in stage 1.
   integer(I4P), intent(in) :: x
   integer(I4P), intent(in) :: parent(:)
   integer(I4P)             :: root
   root = x
   do while (parent(root) /= root)
      root = parent(root)
   enddo
   endfunction uf_find

   subroutine uf_union(a, b, parent, rank_)
   !< Union-by-rank on `parent(:)`. Path compression on lookup.
   integer(I4P), intent(in)    :: a, b
   integer(I4P), intent(inout) :: parent(:), rank_(:)
   integer(I4P)                :: ra, rb

   ra = uf_find(a, parent)
   rb = uf_find(b, parent)
   if (ra == rb) return
   if (rank_(ra) < rank_(rb)) then
      parent(ra) = rb
   else if (rank_(ra) > rank_(rb)) then
      parent(rb) = ra
   else
      parent(rb) = ra
      rank_(ra)  = rank_(ra) + 1_I4P
   endif
   endsubroutine uf_union

endmodule fossil_vertex_pool_object
