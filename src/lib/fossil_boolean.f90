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
!< All four ops are implemented. Coplanar overlaps, shared edges and
!< contacts are handled since issue #33 (shared-boundary sub-triangles,
!< see `apply_selection`). `BOOL_STATUS_NOT_IMPLEMENTED` is returned only
!< for an op code that is none of the four.
!<
!< Since issue #65 the arrangement is exact (implicit points, no snap) and the selection goes by patches: the
!< sub-triangles of a surface connected across edges on no cut are all inside or all outside the other surface, decided
!< once, exactly at a cut edge (`local_inside`) or, for a patch without cuts, by the winding number at a large
!< sub-triangle. Shared (coplanar) sub-triangles are decided exactly (`coplanar_partner`). Faces coplanar by
!< construction but off by round-off are no longer merged: they cross exactly, and the result keeps the slivers
!< (watertight, the right volume; its rounding to doubles can fold a sliver by an ulp).

use fossil_arrangement,    only : arrangement_t, arrangement_initialize, &
                                  arrangement_collect_intersections, &
                                  arrangement_retriangulate, coplanar_partner, cut_sources, &
                                  ARR_STATUS_OK, ARR_STATUS_CDT_FAILED
use fossil_aabb_tree_object, only : aabb_tree_object, wn_node_t
use fossil_facet_object,   only : facet_object
use fossil_implicit_points, only : ip_orient3d
use fossil_predicates,     only : orient3d
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
integer(I4P), parameter :: BOOL_DIFFERENCE = 3_I4P  !< A \ B  : keep what is inside A but not inside B.
integer(I4P), parameter :: BOOL_SYMDIFF    = 4_I4P  !< A △ B  : keep what is inside exactly one of A, B.

! Status codes.
integer(I4P), parameter :: BOOL_STATUS_OK              = 0_I4P
integer(I4P), parameter :: BOOL_STATUS_CDT_FAILED      = 1_I4P  !< Retriangulation could not recover a cut segment.
integer(I4P), parameter :: BOOL_STATUS_NOT_IMPLEMENTED = 2_I4P  !< Op code is none of BOOL_UNION/INTERSECT/DIFFERENCE/SYMDIFF.
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
   integer(I4P),       allocatable                        :: sub_owner(:), sub_source(:), sub_pid(:, :)
   logical,            allocatable                        :: sub_cut(:, :)
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
                                  sub_source=sub_source, status=arr_status, sub_pid=sub_pid, &
                                  sub_cut=sub_cut)
   if (arr_status /= ARR_STATUS_OK) then
      if (present(status)) status = BOOL_STATUS_CDT_FAILED
      allocate(kept_facet(0))
      return
   endif

   ! Stage 4 — tag + select.
   call tag_and_select(arr=arr, sub_facet=sub_facet, sub_owner=sub_owner, sub_source=sub_source, sub_pid=sub_pid, &
                       sub_cut=sub_cut, &
                       facet_a=facet_a, tree_a=tree_a, &
                       facet_b=facet_b, tree_b=tree_b, &
                       op=op, kept_facet=kept_facet, status=status)
   endsubroutine boolean_compute

   subroutine tag_and_select(arr, sub_facet, sub_owner, sub_source, sub_pid, sub_cut, facet_a, tree_a, facet_b, &
                             tree_b, op, kept_facet, status)
   !< Classify the sub-triangles and keep those the op selects (issue #65: by patches).
   !<
   !< - Shared: a sub-triangle inside a coplanar partner of its parent (exact, `coplanar_partner`).
   !< - Patches: the sub-triangles of one surface connected across edges that lie on no cut (registry indices,
   !<   `sub_cut` from the triangulations). The other surface crosses a surface only along cuts, so a patch is all
   !<   inside or all outside it (or touching it without crossing): one decision per patch, by the other surface's
   !<   winding number at the centroid of a large sub-triangle of the patch (the most decisive of its `N_PROBES`
   !<   largest non-shared ones). A sliver 1e-12 thick takes its patch's decision; its own centroid, at round-off
   !<   distance from the other surface, would not settle it.
   type(arrangement_t),             intent(in)            :: arr
   type(facet_object),              intent(in)            :: sub_facet(:)
   integer(I4P),                    intent(in)            :: sub_owner(:)
   integer(I4P),                    intent(in)            :: sub_source(:)
   integer(I4P),                    intent(in)            :: sub_pid(:, :)
   logical,                         intent(in)            :: sub_cut(:, :)
   type(facet_object),              intent(in)            :: facet_a(:)
   type(aabb_tree_object),          intent(in), target    :: tree_a
   type(facet_object),              intent(in)            :: facet_b(:)
   type(aabb_tree_object),          intent(in), target    :: tree_b
   integer(I4P),                    intent(in)            :: op
   type(facet_object), allocatable, intent(out)           :: kept_facet(:)
   integer(I4P),                    intent(out), optional :: status
   integer(I4P), parameter                                :: N_PROBES = 3_I4P
   type(facet_object), allocatable                        :: buf(:)
   integer(I4P)                                           :: used, k, n_sub, partner, p, j, np
   integer(I4P), allocatable                              :: patch(:), part(:), probe(:, :)
   real(R8P),    allocatable                              :: probe_area(:, :), best(:)
   real(R8P)                                              :: wn, orient_dot, ar
   logical,      allocatable                              :: inside(:), decided(:)
   type(wn_node_t), allocatable                           :: cache_a(:), cache_b(:) !< Winding-number caches, built once.
   logical                                                :: keep, flip, shared, in_a, in_b

   if (present(status)) status = BOOL_STATUS_OK
   n_sub = size(sub_facet, dim=1)
   call tree_a%build_winding_payload(facet=facet_a, wn=cache_a)
   call tree_b%build_winding_payload(facet=facet_b, wn=cache_b)
   if (.not. allocated(cache_a)) allocate(cache_a(0:-1))
   if (.not. allocated(cache_b)) allocate(cache_b(0:-1))
   allocate(part(n_sub))
   do k = 1, n_sub
      part(k) = coplanar_partner(arr=arr, gi=sub_source(k), pid=sub_pid(:, k))
   enddo
   call build_patches(sub_owner, sub_pid, sub_cut, patch, np)
   ! probes: the N_PROBES largest non-shared sub-triangles of each patch
   allocate(probe(N_PROBES, np), probe_area(N_PROBES, np))
   probe = 0_I4P ; probe_area = -1._R8P
   do k = 1, n_sub
      if (part(k) > 0_I4P) cycle
      p = patch(k) ; ar = sub_facet(k)%area()
      do j = 1, N_PROBES
         if (ar > probe_area(j, p)) then
            probe(j + 1:N_PROBES, p) = probe(j:N_PROBES - 1, p)
            probe_area(j + 1:N_PROBES, p) = probe_area(j:N_PROBES - 1, p)
            probe(j, p) = k ; probe_area(j, p) = ar
            exit
         endif
      enddo
   enddo
   allocate(inside(np), decided(np), best(np))
   inside = .false. ; decided = .false. ; best = -1._R8P
   ! exact decision at a cut edge of the patch, when it has one
   do k = 1, n_sub
      if (part(k) > 0_I4P) cycle
      p = patch(k)
      if (decided(p)) cycle
      do j = 1, 3
         if (.not. sub_cut(j, k)) cycle
         call local_inside(arr, sub_source(k), sub_pid(:, k), j, in_a, decided(p))
         if (decided(p)) then
            inside(p) = in_a ; exit
         endif
      enddo
   enddo
   do p = 1, np
      if (decided(p)) cycle
      do j = 1, N_PROBES
         k = probe(j, p)
         if (k == 0_I4P) exit
         if (sub_owner(k) == 1_I4P) then
            wn = winding_number(facet=facet_b, tree=tree_b, wn=cache_b, point=sub_facet(k)%centroid)
         else
            wn = winding_number(facet=facet_a, tree=tree_a, wn=cache_a, point=sub_facet(k)%centroid)
         endif
         if (abs(wn - WN_INSIDE_THRESHOLD) > best(p)) then
            best(p) = abs(wn - WN_INSIDE_THRESHOLD) ; inside(p) = wn > WN_INSIDE_THRESHOLD ; decided(p) = .true.
         endif
         if (best(p) > 0.5_R8P - WN_BOUNDARY_TOL) exit ! decisive: no need for the next probe
      enddo
   enddo

   allocate(buf(max(1_I4P, n_sub)))
   used = 0_I4P
   do k = 1, n_sub
      partner = part(k)
      shared = partner > 0_I4P
      orient_dot = 0._R8P
      if (shared) orient_dot = arr%facet(sub_source(k))%normal%dotproduct(rhs=arr%facet(partner)%normal)
      in_a = .false. ; in_b = .false.
      if (.not. shared) then
         if (sub_owner(k) == 1_I4P) then
            in_b = inside(patch(k))
         else
            in_a = inside(patch(k))
         endif
      endif
      call apply_selection(op=op, owner=sub_owner(k), in_a=in_a, in_b=in_b, shared=shared, orient_dot=orient_dot, &
                           keep=keep, flip=flip, status=status)
      if (present(status)) then
         if (status /= BOOL_STATUS_OK) then
            allocate(kept_facet(0))
            return
         endif
      endif
      if (.not. keep) cycle
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

   subroutine local_inside(arr, gi, pid, e, inside, ok)
   !< Exact side of the sub-triangle (registry vertices `pid`) of facet gi beside its cut edge e (opposite vertex e)
   !< with respect to the other surface (issue #65). The edge lies on the other surface: inside one facet g of it, or
   !< on the common edge of two, g1 and g2 (`cut_sources`). With r the sub-triangle's opposite vertex and in(r, g) the
   !< exact orient3d sign (+1 on g's inner side, outward normals): one facet, inside iff in(r, g) > 0; two facets,
   !< convex edge (g2's far vertex on g1's inner side) inside iff both signs >= 0 and one > 0, reflex edge inside iff
   !< one sign > 0, flat edge as one facet. `ok` is false when the sources do not settle it.
   type(arrangement_t), intent(in)  :: arr
   integer(I4P),        intent(in)  :: gi, pid(3), e
   logical,             intent(out) :: inside, ok
   integer(I4P)                     :: g(4), ng, r, s1, s2, cvx, k
   real(R8P)                        :: a(3), b(3), c(3)

   inside = .false. ; ok = .false.
   call cut_sources(arr, gi, pid(mod(e, 3) + 1), pid(mod(e + 1, 3) + 1), g, ng)
   if (ng == 0_I4P) return
   r = pid(e)
   s1 = side(g(1), r)
   if (ng == 1_I4P) then
      if (s1 == 0_I4P) return
      inside = s1 > 0_I4P ; ok = .true. ; return
   endif
   s2 = side(g(2), r)
   ! convexity: the vertex of g2 not on g1's edge, against g1
   cvx = 0_I4P
   do k = 1, 3
      cvx = vertex_side(g(1), arr%facet(g(2))%vertex(k)%x, arr%facet(g(2))%vertex(k)%y, arr%facet(g(2))%vertex(k)%z)
      if (cvx /= 0_I4P) exit
   enddo
   if (cvx > 0_I4P) then
      inside = s1 >= 0_I4P .and. s2 >= 0_I4P .and. (s1 > 0_I4P .or. s2 > 0_I4P)
   elseif (cvx < 0_I4P) then
      inside = s1 > 0_I4P .or. s2 > 0_I4P
   else
      if (s1 == 0_I4P) return
      inside = s1 > 0_I4P
   endif
   ok = s1 /= 0_I4P .or. s2 /= 0_I4P
   contains
      function side(gg, rr) result(sg)
      integer(I4P), intent(in) :: gg, rr
      integer(I4P)             :: sg

      a = [arr%facet(gg)%vertex(1)%x, arr%facet(gg)%vertex(1)%y, arr%facet(gg)%vertex(1)%z]
      b = [arr%facet(gg)%vertex(2)%x, arr%facet(gg)%vertex(2)%y, arr%facet(gg)%vertex(2)%z]
      c = [arr%facet(gg)%vertex(3)%x, arr%facet(gg)%vertex(3)%y, arr%facet(gg)%vertex(3)%z]
      sg = ip_orient3d(arr%pts(rr), a, b, c)
      endfunction side

      function vertex_side(gg, x, y, z) result(sg)
      integer(I4P), intent(in) :: gg
      real(R8P),    intent(in) :: x, y, z
      integer(I4P)             :: sg
      real(R8P)                :: d

      d = orient3d(arr%facet(gg)%vertex(1), arr%facet(gg)%vertex(2), arr%facet(gg)%vertex(3), vector_R8P(x, y, z))
      sg = 0_I4P
      if (d > 0._R8P) sg = 1_I4P
      if (d < 0._R8P) sg = -1_I4P
      endfunction vertex_side
   endsubroutine local_inside

   subroutine build_patches(owner, pid, cut, patch, np)
   !< Patches: union-find of the sub-triangles of each surface over the edges (registry index pairs) that no
   !< triangulation marked as a cut; patch(k) in 1..np.
   integer(I4P),              intent(in)  :: owner(:), pid(:, :)
   logical,                   intent(in)  :: cut(:, :)
   integer(I4P), allocatable, intent(out) :: patch(:)
   integer(I4P),              intent(out) :: np
   integer(I4P), allocatable              :: key(:, :), ord(:), parent(:), root_id(:)
   integer(I4P)                           :: n, k, e, h, i, j, a, b
   logical,      allocatable              :: hcut(:)
   logical                                :: any_cut

   n = size(owner)
   allocate(key(3, 3 * n), hcut(3 * n), ord(3 * n), parent(n))
   h = 0_I4P
   do k = 1, n
      do e = 1, 3
         a = pid(mod(e, 3) + 1, k) ; b = pid(mod(e + 1, 3) + 1, k)
         h = h + 1_I4P
         key(:, h) = [owner(k), min(a, b), max(a, b)] ; hcut(h) = cut(e, k) ; ord(h) = h
      enddo
   enddo
   call sort_keys(key, ord)
   parent = [(k, k = 1, n)]
   i = 1_I4P
   do while (i <= h)
      j = i
      do while (j < h)
         if (any(key(:, ord(j + 1)) /= key(:, ord(i)))) exit
         j = j + 1_I4P
      enddo
      any_cut = any(hcut(ord(i:j)))
      if (.not. any_cut) then
         do k = i + 1, j
            call unite((ord(i) - 1_I4P) / 3_I4P + 1_I4P, (ord(k) - 1_I4P) / 3_I4P + 1_I4P)
         enddo
      endif
      i = j + 1_I4P
   enddo
   allocate(patch(n), root_id(n))
   root_id = 0_I4P ; np = 0_I4P
   do k = 1, n
      a = find(k)
      if (root_id(a) == 0_I4P) then
         np = np + 1_I4P ; root_id(a) = np
      endif
      patch(k) = root_id(a)
   enddo
   contains
      recursive function find(x) result(r)
      integer(I4P), intent(in) :: x
      integer(I4P)             :: r

      r = x
      do while (parent(r) /= r)
         parent(r) = parent(parent(r)) ; r = parent(r)
      enddo
      endfunction find

      subroutine unite(x, y)
      integer(I4P), intent(in) :: x, y
      integer(I4P)             :: rx, ry

      rx = find(x) ; ry = find(y)
      if (rx /= ry) parent(max(rx, ry)) = min(rx, ry)
      endsubroutine unite
   endsubroutine build_patches

   subroutine sort_keys(key, ord)
   !< Bottom-up merge sort of `ord` by the integer triples key(:, ord(:)) (lexicographic).
   integer(I4P), intent(in)    :: key(:, :)
   integer(I4P), intent(inout) :: ord(:)
   integer(I4P), allocatable   :: work(:)
   integer(I4P)                :: n, width, lo, mid, hi, i, j, k

   n = size(ord)
   allocate(work(n))
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (less(ord(j), ord(i))) then
               work(k) = ord(j) ; j = j + 1_I4P
            else
               work(k) = ord(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i <= mid)
            work(k) = ord(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j <= hi)
            work(k) = ord(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         lo = lo + 2_I4P * width
      enddo
      ord = work
      width = 2_I4P * width
   enddo
   contains
      pure function less(x, y) result(l)
      integer(I4P), intent(in) :: x, y
      logical                  :: l
      integer(I4P)             :: c

      l = .false.
      do c = 1, 3
         if (key(c, x) /= key(c, y)) then
            l = key(c, x) < key(c, y) ; return
         endif
      enddo
      endfunction less
   endsubroutine sort_keys

   pure subroutine apply_selection(op, owner, in_a, in_b, shared, orient_dot, &
                                   keep, flip, status)
   !< Per-op selection truth table.
   !<
   !< For **non-shared** sub-triangles (the common case — sub-triangle's
   !< centroid is strictly inside or outside the other surface):
   !<   the standard mesh-arrangement table applies, indexed by (owner, in_other).
   !<   Each row below describes whether to keep the sub-triangle and whether
   !<   to flip its normal (B's facets need flipping when they become the
   !<   *inner* surface of the result, e.g. cavity walls in A\B).
   !<
   !<     op             owner  in_other   keep   flip
   !<     -------------- -----  --------   -----  ----
   !<     UNION          A      out B      yes    no
   !<     UNION          A      in  B      no     -
   !<     UNION          B      out A      yes    no
   !<     UNION          B      in  A      no     -
   !<     INTERSECT      A      out B      no     -
   !<     INTERSECT      A      in  B      yes    no
   !<     INTERSECT      B      out A      no     -
   !<     INTERSECT      B      in  A      yes    no
   !<     DIFFERENCE     A      out B      yes    no
   !<     DIFFERENCE     A      in  B      no     -
   !<     DIFFERENCE     B      out A      no     -
   !<     DIFFERENCE     B      in  A      yes    YES   (cavity wall, flip outward)
   !<     SYMDIFF        A      out B      yes    no    (outer surface)
   !<     SYMDIFF        A      in  B      yes    YES   (cavity wall)
   !<     SYMDIFF        B      out A      yes    no    (outer surface)
   !<     SYMDIFF        B      in  A      yes    YES   (cavity wall)
   !<
   !< For **shared boundary** sub-triangles (centroid on BOTH surfaces — the
   !< two surfaces have coplanar facets that overlap in this region):
   !<   - same orientation (orient_dot > 0): both surfaces locally bound the
   !<     same volume. The shared face is or is not a boundary of the result
   !<     depending on the op. By convention we let A own any kept copy and
   !<     drop B's, to avoid downstream coincident-facet deduplication.
   !<   - opposite orientation (orient_dot < 0): one surface's interior
   !<     coincides with the other's exterior. Rare on ordinary watertight
   !<     meshes; treated symmetrically with same-orientation per op below.
   !<
   !<     op             shared   orient   owner=A keep   owner=B keep
   !<     -------------- ------   ------   -------------  -------------
   !<     UNION          yes      same     yes            no
   !<     UNION          yes      opp      no             no    (cancel)
   !<     INTERSECT      yes      same     yes            no
   !<     INTERSECT      yes      opp      no             no
   !<     DIFFERENCE     yes      same     no             no
   !<     DIFFERENCE     yes      opp      yes (no flip)  no
   !<     SYMDIFF        yes      any      no             no    (shared cancels)
   integer(I4P), intent(in)            :: op, owner
   logical,      intent(in)            :: in_a, in_b
   logical,      intent(in)            :: shared      !< True if centroid is on the OTHER surface (shared boundary).
   real(R8P),    intent(in)            :: orient_dot  !< Dot product of sub-tri normal with closest other-surface normal.
   logical,      intent(out)           :: keep, flip
   integer(I4P), intent(out), optional :: status
   logical                             :: in_other    !< Inside the OTHER surface from owner's perspective.

   keep = .false. ; flip = .false.
   if (present(status)) status = BOOL_STATUS_OK

   ! "in_other" is the inside-flag of whichever surface is NOT the owner.
   if (owner == 1_I4P) then
      in_other = in_b
   else
      in_other = in_a
   endif

   if (shared) then
      ! All four ops drop B's copy on shared boundaries (A owns the boundary
      ! by convention). DIFFERENCE alone keeps A's copy when orientations
      ! oppose. SYMDIFF drops both unconditionally (boundary cancels).
      if (owner == 2_I4P) then
         keep = .false. ; return
      endif
      ! owner == 1 (A) on shared boundary:
      select case (op)
      case (BOOL_UNION)
         keep = (orient_dot > 0._R8P)
      case (BOOL_INTERSECT)
         keep = (orient_dot > 0._R8P)
      case (BOOL_DIFFERENCE)
         keep = (orient_dot < 0._R8P)
      case (BOOL_SYMDIFF)
         keep = .false.
      case default
         if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED
      endselect
      return
   endif

   ! Non-shared sub-triangle — apply the per-op truth table.
   select case (op)
   case (BOOL_UNION)
      keep = .not. in_other        ! keep the parts of each surface that lie outside the other
      flip = .false.
   case (BOOL_INTERSECT)
      keep = in_other              ! keep the parts of each surface that lie inside the other
      flip = .false.
   case (BOOL_DIFFERENCE)
      ! A's outer surface (parts not inside B) ∪ B's surface flipped (parts
      ! inside A — the cavity walls).
      if (owner == 1_I4P) then
         keep = .not. in_b
         flip = .false.
      else
         keep = in_a
         flip = .true.
      endif
   case (BOOL_SYMDIFF)
      ! Always keep; flip when the sub-triangle is inside the other body
      ! (it serves as the inner cavity wall).
      keep = .true.
      flip = in_other
   case default
      if (present(status)) status = BOOL_STATUS_NOT_IMPLEMENTED
   endselect
   endsubroutine apply_selection

endmodule fossil_boolean
