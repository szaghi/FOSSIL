!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.

module fossil_aabb_node_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.
!<
!< This is just a *container* for AABB tree's nodes.

use fossil_aabb_object, only : aabb_object
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use fossil_ray_query, only : ray_hit_t
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : vector_R8P

implicit none
private
public :: aabb_node_object

type :: aabb_node_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree-node class.
   !<
   !< Stores both the AABB payload and the child links used by the SAH BVH path.
   !< The octree path leaves the child links at zero and instead uses implicit
   !< parent/child indexing (`first_child_node = TREE_RATIO * parent + 1`), so the
   !< two trees share the same node type without conflict. A BVH leaf has
   !< `left_child = right_child = 0`; an internal node has both non-zero.
   private
   type(aabb_object), allocatable :: aabb           !< AABB data.
   integer(I4P)                   :: left_child=0   !< BVH left child node index into the tree's node array; 0 = leaf.
   integer(I4P)                   :: right_child=0  !< BVH right child node index into the tree's node array; 0 = leaf.
   contains
      ! public methods
      procedure, pass(self) :: add_facets                  !< Add facets to AABB.
      procedure, pass(self) :: get_left_child              !< Return left_child (BVH); 0 if leaf or octree.
      procedure, pass(self) :: get_right_child             !< Return right_child (BVH); 0 if leaf or octree.
      procedure, pass(self) :: set_children                !< Set (left_child, right_child) — BVH builder.
      procedure, pass(self) :: set_facet_ids               !< Bulk-load the leaf's facet_id list (BVH builder).
      procedure, pass(self) :: bmin                        !< Return AABB bmin.
      procedure, pass(self) :: bmax                        !< Return AABB bmax.
      procedure, pass(self) :: closest_point               !< Return closest point on AABB from point reference.
      procedure, pass(self) :: compute_octants             !< Compute AABB octants.
      procedure, pass(self) :: compute_vertices_nearby     !< Compute vertices nearby.
      procedure, pass(self) :: destroy                     !< Destroy AABB.
      procedure, pass(self) :: distance                    !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: distance_from_facets        !< Return the (square) distance from point to AABB's facets.
      procedure, pass(self) :: update_best_from_facets     !< Update (best d^2, best facet id, best region) over node's facets.
      procedure, pass(self) :: do_ray_intersect            !< Return true if AABB is intersected by ray.
      procedure, pass(self) :: ray_slab_interval           !< Return AABB-ray slab-intersection interval (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_all_facets    !< Append every t >= 0 hit on this node's facets (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_first_facets  !< Update closest hit using this node's facets (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_any_facets    !< Early-exit any-hit on this node's facets (issue #18 §2.5).
      procedure, pass(self) :: facet_id                    !< Return the facets IDs list.
      procedure, pass(self) :: get_aabb_facets             !< Get AABB facets list.
      procedure, pass(self) :: has_facets                  !< Return true if AABB has facets.
      procedure, pass(self) :: initialize                  !< Initialize AABB.
      procedure, pass(self) :: is_allocated                !< Return true is node is allocated.
      procedure, pass(self) :: ray_intersections_number    !< Return ray intersections number.
      procedure, pass(self) :: save_geometry_tecplot_ascii !< Save AABB geometry into Tecplot ascii file.
      procedure, pass(self) :: save_facets_into_file_stl   !< Save facets into file STL.
      procedure, pass(self) :: translate                   !< Translate AABB by delta.
      procedure, pass(self) :: union                       !< Make AABB the union of other AABBs.
      procedure, pass(self) :: update_extents              !< Update AABB bounding box extents.
      ! operators
      generic :: assignment(=) => aabb_node_assign_aabb_node      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_node_assign_aabb_node !< Operator `=`.
endtype aabb_node_object

contains
   ! public methods
   pure subroutine add_facets(self, facet_id, facet, is_exclusive)
   !< Add facets to AABB.
   !<
   !< @note Facets added to AABB are removed to facets list that is also returned.
   class(aabb_node_object), intent(inout)        :: self         !< AABB.
   type(list_id_object),    intent(inout)        :: facet_id     !< List of facets IDs.
   type(facet_object),      intent(in)           :: facet(:)     !< Facets list.
   logical,                 intent(in), optional :: is_exclusive !< Sentinel to enable/disable exclusive addition.

   if (allocated(self%aabb)) call self%aabb%add_facets(facet_id=facet_id, facet=facet, is_exclusive=is_exclusive)
   endsubroutine add_facets

   pure function bmin(self)
   !< Return AABB bmin.
   class(aabb_node_object), intent(in) :: self !< AABB box.
   type(vector_R8P)                    :: bmin !< AABB bmin.

   if (allocated(self%aabb)) bmin = self%aabb%bmin
   endfunction bmin

   pure function get_left_child(self) result(idx)
   !< Return the BVH left-child node index. 0 means leaf (or octree node, which
   !< ignores these fields entirely and uses implicit indexing).
   class(aabb_node_object), intent(in) :: self !< AABB node.
   integer(I4P)                        :: idx  !< Left-child index into the tree's node array.

   idx = self%left_child
   endfunction get_left_child

   pure function get_right_child(self) result(idx)
   !< Return the BVH right-child node index. 0 means leaf (or octree node).
   class(aabb_node_object), intent(in) :: self !< AABB node.
   integer(I4P)                        :: idx  !< Right-child index into the tree's node array.

   idx = self%right_child
   endfunction get_right_child

   pure subroutine set_children(self, left_child, right_child)
   !< Set both child links in one call. Used by the BVH builder; the octree
   !< path never invokes this.
   class(aabb_node_object), intent(inout) :: self        !< AABB node.
   integer(I4P),            intent(in)    :: left_child  !< Left-child node index (0 = leaf).
   integer(I4P),            intent(in)    :: right_child !< Right-child node index (0 = leaf).

   self%left_child  = left_child
   self%right_child = right_child
   endsubroutine set_children

   pure subroutine set_facet_ids(self, ids)
   !< Bulk-load this node's facet_id list with the given ids array, replacing any
   !< previous contents. Used by the BVH builder to populate a leaf without going
   !< through `add_facets`'s centroid-inside-bbox filter (the partition step in the
   !< builder has already done that work). No-op if the underlying aabb has not
   !< been initialized (e.g. an empty node).
   class(aabb_node_object), intent(inout) :: self   !< AABB node.
   integer(I4P),            intent(in)    :: ids(:) !< Facet ids to store.

   if (allocated(self%aabb)) call self%aabb%facet_id%initialize(id=ids)
   endsubroutine set_facet_ids

   pure function bmax(self)
   !< Return AABB bmax.
   class(aabb_node_object), intent(in) :: self !< AABB box.
   type(vector_R8P)                    :: bmax !< AABB bmax.

   if (allocated(self%aabb)) bmax = self%aabb%bmax
   endfunction bmax

   pure function closest_point(self, point) result(closest)
   !< Return closest point on (or in) AABB from point reference.
   class(aabb_node_object), intent(in) :: self    !< AABB box.
   type(vector_R8P),        intent(in) :: point   !< Point reference.
   type(vector_R8P)                    :: closest !< Closest point on (on in) aabb to point.

   closest = MaxR8P
   if (allocated(self%aabb)) closest = self%aabb%closest_point(point=point)
   endfunction closest_point

   pure subroutine compute_octants(self, octant)
   !< Return AABB octants.
   class(aabb_node_object), intent(in)  :: self      !< AABB.
   type(aabb_object),       intent(out) :: octant(8) !< AABB octants.
   type(vector_R8P)                     :: vertex(8) !< AABB vertices.
   integer(I4P)                         :: o         !< Counter.

   call self%aabb%compute_octants(octant=octant)
   endsubroutine compute_octants

   pure subroutine compute_vertices_nearby(self, facet, tolerance_to_be_nearby)
   !< Compute vertices nearby (loose tolerance only; see facet variant).
   class(aabb_node_object), intent(in)    :: self                   !< AABB.
   type(facet_object),      intent(inout) :: facet(:)               !< Facets list.
   real(R8P),               intent(in)    :: tolerance_to_be_nearby !< Tolerance to identify nearby vertices.

   if (allocated(self%aabb)) call self%aabb%compute_vertices_nearby(facet=facet, tolerance_to_be_nearby=tolerance_to_be_nearby)
   endsubroutine compute_vertices_nearby

   elemental subroutine destroy(self)
   !< Destroy AABB.
   class(aabb_node_object), intent(inout) :: self  !< AABB.

   if (allocated(self%aabb)) then
      call self%aabb%destroy
      deallocate(self%aabb)
   endif
   self%left_child  = 0_I4P
   self%right_child = 0_I4P
   endsubroutine destroy

   pure function distance(self, point)
   !< Return the (square) distance from point to AABB.
   class(aabb_node_object), intent(in) :: self     !< AABB.
   type(vector_R8P),        intent(in) :: point    !< Point reference.
   real(R8P)                           :: distance !< Distance from point to AABB.

   distance = MaxR8P
   if (allocated(self%aabb)) distance = self%aabb%distance(point=point)
   endfunction distance

   pure function distance_from_facets(self, facet, point) result(distance)
   !< Return the (square) distance from point to AABB's facets.
   class(aabb_node_object), intent(in) :: self      !< AABB.
   type(facet_object),      intent(in) :: facet(:)  !< Facets list.
   type(vector_R8P),        intent(in) :: point     !< Point reference.
   real(R8P)                           :: distance  !< Distance from point to AABB's facets.

   distance = MaxR8P
   if (allocated(self%aabb)) distance = self%aabb%distance_from_facets(facet=facet, point=point)
   endfunction distance_from_facets

   pure subroutine update_best_from_facets(self, facet, point, best, best_facet, best_region)
   !< Pass-through to the underlying AABB's update_best_from_facets — no-op when the node is unallocated.
   class(aabb_node_object), intent(in)    :: self        !< AABB node.
   type(facet_object),      intent(in)    :: facet(:)    !< Facets list.
   type(vector_R8P),        intent(in)    :: point       !< Point reference.
   real(R8P),               intent(inout) :: best        !< Running best squared distance.
   integer(I4P),            intent(inout) :: best_facet  !< Running best facet id.
   integer(I4P),            intent(inout) :: best_region !< Running best Voronoi region tag.

   if (allocated(self%aabb)) call self%aabb%update_best_from_facets(facet=facet, point=point, &
                                                                    best=best, best_facet=best_facet, best_region=best_region)
   endsubroutine update_best_from_facets

   pure function do_ray_intersect(self, ray_origin, ray_direction) result(do_intersect)
   !< Return true if AABB is intersected by ray from origin and oriented as ray direction vector.
   class(aabb_node_object), intent(in) :: self          !< AABB.
   type(vector_R8P),        intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction !< Ray direction.
   logical                             :: do_intersect  !< Test result.

   do_intersect = .false.
   if (allocated(self%aabb)) do_intersect = self%aabb%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)
   endfunction do_ray_intersect

   subroutine intersect_ray_all_facets(self, facet, ray_origin, ray_direction, tmp, n_hits, n_overflow)
   !< Append every facet hit on this node into `tmp` (issue #18 §2.5).
   !<
   !< `n_hits` is the running write index, advanced once per accepted hit. If
   !< `n_hits` would exceed `size(tmp)`, the hit is counted in `n_overflow`
   !< instead — the caller checks this flag, grows `tmp`, and retries (or
   !< pre-sizes `tmp` generously).
   class(aabb_node_object),      intent(in)    :: self           !< AABB node.
   type(facet_object),           intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),             intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),             intent(in)    :: ray_direction  !< Ray direction.
   type(ray_hit_t), allocatable, intent(inout) :: tmp(:)         !< Hit accumulator.
   integer(I4P),                 intent(inout) :: n_hits         !< Write index (advanced).
   integer(I4P),                 intent(inout) :: n_overflow     !< Hits dropped because tmp was full.
   integer(I4P)                                :: nf, f, fid     !< Facet-list counters.
   real(R8P)                                   :: t, u, v        !< Per-facet output.
   logical                                     :: facet_hit      !< Per-facet flag.

   if (.not. allocated(self%aabb)) return
   nf = self%aabb%facet_id%ids_number
   do f = 1_I4P, nf
      fid = self%aabb%facet_id%id(f)
      call facet(fid)%intersect_ray(ray_origin=ray_origin, ray_direction=ray_direction, &
                                    t=t, u=u, v=v, hit=facet_hit)
      if (.not. facet_hit) cycle
      if (t < 0._R8P) cycle
      if (n_hits + 1_I4P > size(tmp, kind=I4P)) then
         n_overflow = n_overflow + 1_I4P
         cycle
      endif
      n_hits = n_hits + 1_I4P
      tmp(n_hits)%facet_id = fid
      tmp(n_hits)%t        = t
      tmp(n_hits)%point    = ray_origin + t * ray_direction
   enddo
   endsubroutine intersect_ray_all_facets

   subroutine intersect_ray_first_facets(self, facet, ray_origin, ray_direction, hit, has_hit, t_best)
   !< Update closest hit using this node's facets (issue #18 §2.5).
   !<
   !< Mutates `hit / has_hit / t_best` if any facet on this node hits at
   !< `0 <= t < t_best`. Caller is responsible for the AABB-level prune (don't
   !< call this on a node whose AABB the ray missed).
   class(aabb_node_object), intent(in)    :: self           !< AABB node.
   type(facet_object),      intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)    :: ray_direction  !< Ray direction.
   type(ray_hit_t),         intent(inout) :: hit            !< Running closest hit.
   logical,                 intent(inout) :: has_hit        !< Set true on first hit.
   real(R8P),               intent(inout) :: t_best         !< Running best t.
   integer(I4P)                           :: nf, f, fid     !< Counters.
   real(R8P)                              :: t, u, v        !< Per-facet output.
   logical                                :: facet_hit      !< Per-facet flag.

   if (.not. allocated(self%aabb)) return
   nf = self%aabb%facet_id%ids_number
   do f = 1_I4P, nf
      fid = self%aabb%facet_id%id(f)
      call facet(fid)%intersect_ray(ray_origin=ray_origin, ray_direction=ray_direction, &
                                    t=t, u=u, v=v, hit=facet_hit)
      if (.not. facet_hit) cycle
      if (t < 0._R8P) cycle
      if (t >= t_best) cycle
      t_best = t
      has_hit = .true.
      hit%facet_id = fid
      hit%t        = t
      hit%point    = ray_origin + t * ray_direction
   enddo
   endsubroutine intersect_ray_first_facets

   subroutine intersect_ray_any_facets(self, facet, ray_origin, ray_direction, max_t, found)
   !< Set `found = .true.` and return as soon as any facet on this node hits at
   !< `0 <= t <= max_t` (issue #18 §2.5). No record kept; the caller only cares
   !< about existence (shadow ray, occlusion test).
   class(aabb_node_object), intent(in)    :: self           !< AABB node.
   type(facet_object),      intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)    :: ray_direction  !< Ray direction.
   real(R8P),               intent(in)    :: max_t          !< Upper t bound (use MaxR8P for "any hit").
   logical,                 intent(inout) :: found          !< Set true on first hit; left untouched otherwise.
   integer(I4P)                           :: nf, f, fid     !< Counters.
   real(R8P)                              :: t, u, v        !< Per-facet output.
   logical                                :: facet_hit      !< Per-facet flag.

   if (.not. allocated(self%aabb)) return
   nf = self%aabb%facet_id%ids_number
   do f = 1_I4P, nf
      fid = self%aabb%facet_id%id(f)
      call facet(fid)%intersect_ray(ray_origin=ray_origin, ray_direction=ray_direction, &
                                    t=t, u=u, v=v, hit=facet_hit)
      if (.not. facet_hit) cycle
      if (t < 0._R8P) cycle
      if (t > max_t)   cycle
      found = .true.
      return
   enddo
   endsubroutine intersect_ray_any_facets

   pure subroutine ray_slab_interval(self, ray_origin, ray_direction, t_near, t_far, do_intersect)
   !< Forward to `aabb%ray_slab_interval`. Returns `do_intersect = .false.` and a
   !< degenerate interval when this node carries no AABB (unallocated container).
   class(aabb_node_object), intent(in)  :: self          !< AABB node.
   type(vector_R8P),        intent(in)  :: ray_origin    !< Ray origin.
   type(vector_R8P),        intent(in)  :: ray_direction !< Ray direction.
   real(R8P),               intent(out) :: t_near        !< Slab interval near (entry) parameter.
   real(R8P),               intent(out) :: t_far         !< Slab interval far (exit) parameter.
   logical,                 intent(out) :: do_intersect  !< True if half-ray hits the box.

   t_near = 0._R8P
   t_far  = 0._R8P
   do_intersect = .false.
   if (allocated(self%aabb)) call self%aabb%ray_slab_interval(ray_origin=ray_origin, ray_direction=ray_direction, &
                                                              t_near=t_near, t_far=t_far, do_intersect=do_intersect)
   endsubroutine ray_slab_interval

   pure function facet_id(self)
   !< Return facets IDs list.
   class(aabb_node_object), intent(in) :: self     !< AABB box.
   type(list_id_object)                :: facet_id !< List of facets IDs contained into AABB.

   if (self%is_allocated()) facet_id = self%aabb%facet_id
   endfunction facet_id

   pure subroutine get_aabb_facets(self, facet, aabb_facet)
   !< Get AABB facets list.
   class(aabb_node_object), intent(in)               :: self          !< AABB.
   type(facet_object),      intent(in)               :: facet(:)      !< Whole facets list.
   type(facet_object),      intent(out), allocatable :: aabb_facet(:) !< AABB facets list.

   if (allocated(self%aabb)) call self%aabb%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
   endsubroutine get_aabb_facets

   pure function has_facets(self)
   !< Return true if AABB has facets.
   class(aabb_node_object), intent(in) :: self       !< AABB box.
   logical                             :: has_facets !< Check result.

   has_facets = allocated(self%aabb)
   if (has_facets) has_facets = self%aabb%has_facets()
   endfunction has_facets

   pure subroutine initialize(self, facet, bmin, bmax)
   !< Initialize AABB.
   class(aabb_node_object), intent(inout)        :: self     !< AABB.
   type(facet_object),      intent(in), optional :: facet(:) !< Facets list.
   type(vector_R8P),        intent(in), optional :: bmin     !< Minimum point of AABB.
   type(vector_R8P),        intent(in), optional :: bmax     !< Maximum point of AABB.

   call self%destroy
   if (present(facet).or.(present(bmin).and.present(bmin))) then
      allocate(self%aabb)
      call self%aabb%initialize(facet=facet, bmin=bmin, bmax=bmax)
   endif
   endsubroutine initialize

   pure function is_allocated(self)
   !< Return true if node is allocated.
   class(aabb_node_object), intent(in) :: self         !< AABB box.
   logical                             :: is_allocated !< Check result.

   is_allocated = allocated(self%aabb)
   endfunction is_allocated

   pure function ray_intersections_number(self, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number.
   class(aabb_node_object), intent(in) :: self                 !< AABB.
   type(facet_object),      intent(in) :: facet(:)             !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                        :: intersections_number !< Intersection number.

   intersections_number = 0
   if (allocated(self%aabb)) &
      intersections_number = self%aabb%ray_intersections_number(facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
   endfunction ray_intersections_number

   subroutine  save_geometry_tecplot_ascii(self, file_unit, aabb_name)
   !< Save AABB geometry into Tecplot ascii file.
   class(aabb_node_object), intent(in)           :: self       !< AABB.
   integer(I4P),            intent(in)           :: file_unit  !< File unit.
   character(*),            intent(in), optional :: aabb_name  !< Name of AABB.

   if (allocated(self%aabb)) call self%aabb%save_geometry_tecplot_ascii(file_unit=file_unit, aabb_name=aabb_name)
   endsubroutine  save_geometry_tecplot_ascii

   subroutine save_facets_into_file_stl(self, facet, file_name, is_ascii)
   !< Save facets into file STL.
   class(aabb_node_object), intent(in) :: self      !< AABB.
   type(facet_object),      intent(in) :: facet(:)  !< Facets list.
   character(*),            intent(in) :: file_name !< File name.
   logical,                 intent(in) :: is_ascii  !< Sentinel for file format.

   if (allocated(self%aabb)) call self%aabb%save_facets_into_file_stl(facet=facet, file_name=file_name, is_ascii=is_ascii)
   endsubroutine save_facets_into_file_stl

   elemental subroutine translate(self, delta)
   !< Translate AABB by delta.
   class(aabb_node_object), intent(inout) :: self  !< AABB.
   type(vector_R8P),        intent(in)    :: delta !< Delta of translation.

   if (allocated(self%aabb)) call self%aabb%translate(delta=delta)
   endsubroutine translate

   pure subroutine union(self, node, id)
   !< Make AABB the union of other AABBs.
   class(aabb_node_object), intent(inout) :: self    !< AABB.
   type(aabb_node_object),  intent(in)    :: node(:) !< Nodes list.
   integer(I4P),            intent(in)    :: id(:)   !< Nodes ID list.
   integer(I4P)                           :: i       !< Counter.

   call self%destroy
   allocate(self%aabb)
   do i=1, size(id, dim=1)
      if (allocated(node(id(i))%aabb)) call self%aabb%union(other=node(id(i))%aabb)
   enddo
   endsubroutine union

   pure subroutine update_extents(self, facet)
   !< Update AABB bounding box extents.
   class(aabb_node_object), intent(inout) :: self     !< AABB.
   type(facet_object),      intent(in)    :: facet(:) !< Facets list.

   if (allocated(self%aabb)) call self%aabb%update_extents(facet=facet)
   endsubroutine update_extents

   ! operators
   ! =
   pure subroutine aabb_node_assign_aabb_node(lhs, rhs)
   !< Operator `=`.
   class(aabb_node_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_node_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(lhs%aabb)) then
      call lhs%aabb%destroy
      deallocate(lhs%aabb)
   endif
   if (allocated(rhs%aabb)) then
      allocate(lhs%aabb)
      lhs%aabb = rhs%aabb
   endif
   lhs%left_child  = rhs%left_child
   lhs%right_child = rhs%right_child
   endsubroutine aabb_node_assign_aabb_node
endmodule fossil_aabb_node_object
