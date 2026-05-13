!< FOSSIL, mesh-mesh boolean operations (issue #18 §1.1).

module fossil_boolean
!< FOSSIL, mesh-mesh boolean operations (union / intersect / difference / symdiff).
!<
!< This module implements the **selection + stitch** stages of the Mesh
!< Arrangements pipeline (Zhou et al. SIGGRAPH 2016). The earlier stages
!< (arrangement initialization, intersection collection, retriangulation)
!< are in `fossil_arrangement`; the underlying winding-number primitive is
!< in `fossil_winding_number`.
!<
!< Per-op semantics — what the WN-tag selection truth table preserves:
!<
!<     op           keep A-facets where    keep B-facets where    flip B?
!<     ------------ ---------------------- ---------------------- -------
!<     UNION        outside B              outside A              no
!<     INTERSECT    inside  B              inside  A              no
!<     DIFFERENCE   outside B              inside  A              YES
!<     SYMDIFF      outside B              outside A              YES
!<                  + inside  B            + inside  A
!<
!< The flip on B (in DIFFERENCE / SYMDIFF) is needed because B's outward
!< normal points away from B's interior, but those B-facets are now serving
!< as the *inner* surface of the result (the cavity walls where B carved
!< into A). Flipping the normal points it into A's exterior, which is the
!< new outside of the result.
!<
!< Stage 1 of this PR ships only `BOOL_DIFFERENCE` end-to-end and tested.
!< The other ops have parameter constants exposed for API stability but
!< return STATUS_BOOL_NOT_IMPLEMENTED at present.

use fossil_arrangement,    only : arrangement_t, arrangement_initialize, &
                                  arrangement_collect_intersections, &
                                  arrangement_retriangulate, &
                                  ARR_STATUS_OK, ARR_STATUS_CDT_FAILED
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,   only : facet_object
use fossil_winding_number, only : winding_number
use penf,                  only : I4P, R8P, MaxR8P
use vecfor,                only : vector_R8P

implicit none
private
public :: BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF
public :: BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT
public :: boolean_compute

! Op codes — the API surface for `surface%boolean(other, op, status)`.
integer(I4P), parameter :: BOOL_UNION      = 1_I4P  !< A ∪ B  : keep everything that is inside A or inside B.
integer(I4P), parameter :: BOOL_INTERSECT  = 2_I4P  !< A ∩ B  : keep only what is inside both.
integer(I4P), parameter :: BOOL_DIFFERENCE = 3_I4P  !< A \ B  : keep what is inside A but not inside B (this PR's MVP).
integer(I4P), parameter :: BOOL_SYMDIFF    = 4_I4P  !< A △ B  : keep what is inside exactly one of A, B.

! Status codes.
integer(I4P), parameter :: BOOL_STATUS_OK              = 0_I4P
integer(I4P), parameter :: BOOL_STATUS_CDT_FAILED      = 1_I4P  !< Retriangulation could not recover a cut segment.
integer(I4P), parameter :: BOOL_STATUS_NOT_IMPLEMENTED = 2_I4P  !< Op other than DIFFERENCE in this PR.
integer(I4P), parameter :: BOOL_STATUS_EMPTY_INPUT     = 3_I4P  !< One of the input surfaces has zero facets.

! Three-state classification thresholds for the WN test:
!   wn > 0.5 + WN_BOUNDARY_TOL  → strictly inside
!   wn < 0.5 - WN_BOUNDARY_TOL  → strictly outside
!   |wn - 0.5| <= WN_BOUNDARY_TOL  → on the boundary (centroid lies on the surface)
!
! WN_BOUNDARY_TOL is the slack to declare "near-0.5" as boundary. The hierarchical
! WN with default beta=2 has typical absolute error around 0.05 for points near the
! surface (Barill et al. 2018), so 0.1 is a comfortable cutoff that handles the
! per-facet-centroid-on-boundary case without false positives in the interior.
real(R8P), parameter :: WN_INSIDE_THRESHOLD = 0.5_R8P
real(R8P), parameter :: WN_BOUNDARY_TOL     = 0.1_R8P

contains

   subroutine boolean_compute(facet_a, tree_a, facet_b, tree_b, op, kept_facet, status)
   !< Compute the result of `op(A, B)` and emit the result as a flat facet
   !< array in `kept_facet`.
   !<
   !< Pipeline:
   !<   1. arrangement_initialize  — concat facets, allocate cut lists.
   !<   2. arrangement_collect_intersections — populate cuts via tri-tri.
   !<   3. arrangement_retriangulate — split cut facets via CDT.
   !<   4. tag_and_select  — per-sub-triangle WN query + op-specific keep / flip.
   !<
   !< Caller is responsible for deciding what to do with `kept_facet` (typically
   !< adopt it into a `surface_stl_object`).
   type(facet_object),              intent(in)            :: facet_a(:)
   type(aabb_tree_object),          intent(in), target    :: tree_a
   type(facet_object),              intent(in)            :: facet_b(:)
   type(aabb_tree_object),          intent(in), target    :: tree_b
   integer(I4P),                    intent(in)            :: op
   type(facet_object), allocatable, intent(out)           :: kept_facet(:)
   integer(I4P),                    intent(out), optional :: status
   type(arrangement_t)                                    :: arr
   type(facet_object), allocatable                        :: sub_facet(:)
   integer(I4P),       allocatable                        :: sub_owner(:), sub_source(:)
   integer(I4P)                                           :: arr_status

   if (present(status)) status = BOOL_STATUS_OK

   if (size(facet_a, dim=1) == 0 .or. size(facet_b, dim=1) == 0) then
      if (present(status)) status = BOOL_STATUS_EMPTY_INPUT
      allocate(kept_facet(0))
      return
   endif

   ! Stages 1-3.
   call arrangement_initialize(arr=arr, facet_a=facet_a, tree_a=tree_a, &
                                        facet_b=facet_b, tree_b=tree_b)
   call arrangement_collect_intersections(arr=arr)
   call arrangement_retriangulate(arr=arr, sub_facet=sub_facet, sub_owner=sub_owner, &
                                  sub_source=sub_source, status=arr_status)
   if (arr_status /= ARR_STATUS_OK) then
      if (present(status)) status = BOOL_STATUS_CDT_FAILED
      allocate(kept_facet(0))
      return
   endif

   ! Stage 4 — tag + select.
   call tag_and_select(sub_facet=sub_facet, sub_owner=sub_owner, &
                       facet_a=facet_a, tree_a=tree_a, &
                       facet_b=facet_b, tree_b=tree_b, &
                       op=op, kept_facet=kept_facet, status=status)
   endsubroutine boolean_compute

   subroutine tag_and_select(sub_facet, sub_owner, facet_a, tree_a, facet_b, tree_b, &
                             op, kept_facet, status)
   !< For each sub-triangle, query both WN fields at its centroid, apply the
   !< op-specific selection rule, emit the keepers into `kept_facet`.
   !<
   !< Per-sub-triangle WN query is the dominant cost: one `winding_number`
   !< call per sub-triangle per surface. The hierarchical Barnes-Hut path in
   !< `fossil_winding_number` keeps each call sub-linear in facet count, so
   !< the boolean's overall cost is approximately
   !<     O(N log N)  where N = total sub-triangles.
   type(facet_object),              intent(in)            :: sub_facet(:)
   integer(I4P),                    intent(in)            :: sub_owner(:)
   type(facet_object),              intent(in)            :: facet_a(:)
   type(aabb_tree_object),          intent(in), target    :: tree_a
   type(facet_object),              intent(in)            :: facet_b(:)
   type(aabb_tree_object),          intent(in), target    :: tree_b
   integer(I4P),                    intent(in)            :: op
   type(facet_object), allocatable, intent(out)           :: kept_facet(:)
   integer(I4P),                    intent(out), optional :: status
   type(facet_object), allocatable                        :: buf(:), tmp(:)
   integer(I4P)                                           :: cap, used, k, n_sub
   real(R8P)                                              :: wn_a, wn_b
   logical                                                :: in_a, in_b, keep, flip
   type(vector_R8P)                                       :: c

   n_sub = size(sub_facet, dim=1)
   cap = max(64_I4P, n_sub)
   allocate(buf(cap))
   used = 0_I4P

   do k = 1, n_sub
      ! Centroid in 3D — the facet's `centroid` is set by compute_metrix,
      ! which arrangement_retriangulate calls on every emitted sub-facet.
      c    = sub_facet(k)%centroid
      wn_a = winding_number(facet=facet_a, tree=tree_a, point=c)
      wn_b = winding_number(facet=facet_b, tree=tree_b, point=c)

      ! Detect "shared boundary" sub-triangles: those whose centroid lies on
      ! a coplanar facet of the OTHER surface.
      !
      ! Detection rule: distance to the closest other-surface facet is below
      ! tolerance AND the normals are (anti-)parallel. The WN-based detection
      ! that would naively work in symbolic geometry (WN ≈ 0.5 on the boundary)
      ! does not work in our hierarchical implementation: WN at a point exactly
      ! on the boundary surface evaluates to ~1 (treated as interior) rather
      ! than 0.5, because Bærentzen's solid-angle convention plus hierarchical
      ! summation produce a 4π total at boundary points rather than 2π. So we
      ! use the geometric test directly.
      block
         logical   :: shared
         real(R8P) :: orient_dot
         if (sub_owner(k) == 1_I4P) then
            shared = is_on_other_surface(sub=sub_facet(k), other_facet=facet_b, &
                                          other_tree=tree_b, orient_dot=orient_dot)
         else
            shared = is_on_other_surface(sub=sub_facet(k), other_facet=facet_a, &
                                          other_tree=tree_a, orient_dot=orient_dot)
         endif
         in_a = (wn_a > WN_INSIDE_THRESHOLD)
         in_b = (wn_b > WN_INSIDE_THRESHOLD)

         call apply_selection(op=op, owner=sub_owner(k), in_a=in_a, in_b=in_b, &
                              shared=shared, orient_dot=orient_dot, &
                              keep=keep, flip=flip, status=status)
      endblock

      if (present(status)) then
         if (status /= BOOL_STATUS_OK) then
            allocate(kept_facet(0))
            return
         endif
      endif
      if (.not. keep) cycle
      if (used == cap) then
         cap = 2 * cap
         allocate(tmp(cap))
         tmp(1:used) = buf(1:used)
         call move_alloc(from=tmp, to=buf)
      endif
      used = used + 1
      buf(used) = sub_facet(k)
      if (flip) then
         call buf(used)%reverse_normal
         call buf(used)%compute_metrix
      endif
   enddo

   allocate(kept_facet(used))
   kept_facet(1:used) = buf(1:used)
   endsubroutine tag_and_select

   pure subroutine apply_selection(op, owner, in_a, in_b, shared, orient_dot, &
                                   keep, flip, status)
   !< Per-op selection truth table.
   !<
   !< For **non-shared** sub-triangles (the common case — sub-triangle's
   !< centroid is strictly inside or outside the other surface):
   !<   the standard mesh-arrangement table applies, indexed by (owner, in_other).
   !<
   !< For **shared boundary** sub-triangles (centroid on BOTH surfaces — i.e.
   !< the two surfaces have coplanar facets that overlap in this region):
   !<   the answer depends on whether the two surfaces' outward normals point
   !<   the same way at this point.
   !<     - same orientation (orient_dot > 0): both surfaces locally bound the
   !<       same volume. For DIFFERENCE the shared face does NOT bound the
   !<       result (just below: inside both → outside A\B; just above: outside
   !<       both → outside A\B), so drop both copies. By convention we drop
   !<       both A and B copies (rather than keeping one) to avoid downstream
   !<       deduplication.
   !<     - opposite orientation (orient_dot < 0): one surface's interior
   !<       coincides with the other's exterior. For DIFFERENCE the shared
   !<       face IS a boundary of the result; keep A's copy (B's would be
   !<       inverted-orientation, which `reverse_normal` could fix but for
   !<       the simpler symmetric case we let A own it).
   integer(I4P), intent(in)            :: op, owner
   logical,      intent(in)            :: in_a, in_b
   logical,      intent(in)            :: shared      !< True if centroid is on the OTHER surface (shared boundary).
   real(R8P),    intent(in)            :: orient_dot  !< Dot product of sub-tri normal with closest other-surface normal.
   logical,      intent(out)           :: keep, flip
   integer(I4P), intent(out), optional :: status

   keep = .false. ; flip = .false.
   if (present(status)) status = BOOL_STATUS_OK

   select case (op)
   case (BOOL_DIFFERENCE)
      if (shared) then
         ! Shared boundary: opposite-orientation faces ARE boundaries of A\B
         ! (kept by A's copy; B's copy would be opposite, drop). Same-orientation
         ! faces are NOT boundaries of A\B (both surfaces enclose the same side
         ! locally, so the shared face is interior of both A and B); drop both.
         if (owner == 1_I4P) then
            keep = (orient_dot < 0._R8P)
            flip = .false.
         else  ! owner == 2 (B): always drop on shared boundary
            keep = .false.
         endif
      else
         ! Non-shared: standard truth table.
         ! A \ B = A's outer surface (the parts not inside B) ∪ B's surface
         ! flipped (the parts that lie inside A become the cavity walls).
         if (owner == 1_I4P) then
            keep = .not. in_b   ! A-facets outside B
            flip = .false.
         else  ! owner == 2 (B)
            keep = in_a         ! B-facets inside A
            flip = .true.       ! flip to point outward (into A's exterior)
         endif
      endif

   case (BOOL_UNION, BOOL_INTERSECT, BOOL_SYMDIFF)
      if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED

   case default
      if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED
   endselect
   endsubroutine apply_selection

   function is_on_other_surface(sub, other_facet, other_tree, orient_dot) result(yes)
   !< True iff `sub`'s centroid lies on the OTHER surface (within tolerance)
   !< AND the closest other-surface facet's normal is (anti-)parallel to `sub`'s
   !< normal. Returns the dot product of normals via `orient_dot` for
   !< `apply_selection` to disambiguate same-vs-opposite orientation.
   !<
   !< Distance and normal-parallelism tolerances:
   !<   - Distance: `SHARED_DIST_TOL` ≈ 1e-9 in the surface's coordinate units.
   !<     Any sub-triangle centroid this close to another surface's facet is
   !<     treated as coincident. Note that `distance_tree_with_region` returns
   !<     SQUARED distance, so we compare against `SHARED_DIST_TOL**2`.
   !<   - Normal parallelism: `|orient_dot| > SHARED_NORM_TOL` ≈ 0.99 (i.e.,
   !<     angles within ~8° of parallel/antiparallel). Looser than 1.0 to
   !<     accommodate per-facet normal jitter on retriangulated cuts.
   type(facet_object),     intent(in)         :: sub
   type(facet_object),     intent(in)         :: other_facet(:)
   type(aabb_tree_object), intent(in), target :: other_tree
   real(R8P),              intent(out)        :: orient_dot
   logical                                    :: yes
   real(R8P)                                  :: distance_sq
   integer(I4P)                               :: closest_id, region_id
   real(R8P), parameter :: SHARED_DIST_TOL = 1.0e-9_R8P
   real(R8P), parameter :: SHARED_NORM_TOL = 0.99_R8P

   yes = .false. ; orient_dot = 0._R8P
   distance_sq = MaxR8P
   closest_id  = 0_I4P
   region_id   = 0_I4P
   call other_tree%distance_tree_with_region(facet=other_facet, point=sub%centroid, &
                                              distance=distance_sq, closest_facet=closest_id, &
                                              closest_region=region_id)
   if (closest_id <= 0_I4P .or. closest_id > size(other_facet, dim=1)) return
   if (distance_sq > SHARED_DIST_TOL * SHARED_DIST_TOL) return
   orient_dot = sub%normal%dotproduct(rhs=other_facet(closest_id)%normal)
   yes = (abs(orient_dot) > SHARED_NORM_TOL)
   endfunction is_on_other_surface

endmodule fossil_boolean
