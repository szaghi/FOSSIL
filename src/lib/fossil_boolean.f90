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
use penf,                  only : I4P, R8P
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

! WN threshold for in/out classification: WN ~ 1 inside, ~ 0 outside.
real(R8P), parameter :: WN_INSIDE_THRESHOLD = 0.5_R8P

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
      in_a = (wn_a > WN_INSIDE_THRESHOLD)
      in_b = (wn_b > WN_INSIDE_THRESHOLD)
      call apply_selection(op=op, owner=sub_owner(k), in_a=in_a, in_b=in_b, &
                           keep=keep, flip=flip, status=status)
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

   pure subroutine apply_selection(op, owner, in_a, in_b, keep, flip, status)
   !< Per-op selection truth table — given a sub-triangle's owner (1=A, 2=B)
   !< and inside-flags, decide whether to keep it and whether to flip its
   !< normal. See module docstring for the full table.
   !<
   !< Sub-triangles whose owner-side WN is "inside" (e.g. an A-owned sub-tri
   !< with in_a=true) are sub-triangles of A's surface that lie *inside* A.
   !< Those are the sub-triangles that are interior to the source body and
   !< therefore not on the result's boundary — by construction A's surface
   !< doesn't have facets inside A. But the **WN test** can return the wrong
   !< answer for a sub-triangle whose centroid lies on A's boundary (WN ~= 0.5);
   !< those edge cases are why production implementations use exact predicates
   !< on the boundary. For the §1.1 MVP the WN_INSIDE_THRESHOLD = 0.5 cutoff
   !< handles the common case; pathological cases would manifest as missing
   !< or duplicated boundary facets in the output.
   integer(I4P), intent(in)            :: op, owner
   logical,      intent(in)            :: in_a, in_b
   logical,      intent(out)           :: keep, flip
   integer(I4P), intent(out), optional :: status

   keep = .false. ; flip = .false.
   if (present(status)) status = BOOL_STATUS_OK

   select case (op)
   case (BOOL_DIFFERENCE)
      ! A \ B = A's outer surface (the parts not inside B) ∪ B's surface
      ! flipped (the parts that lie inside A become the cavity walls).
      if (owner == 1_I4P) then
         keep = .not. in_b   ! A-facets outside B
         flip = .false.
      else  ! owner == 2 (B)
         keep = in_a         ! B-facets inside A
         flip = .true.       ! flip to point outward (into A's exterior)
      endif

   case (BOOL_UNION, BOOL_INTERSECT, BOOL_SYMDIFF)
      if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED

   case default
      if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED
   endselect
   endsubroutine apply_selection

endmodule fossil_boolean
