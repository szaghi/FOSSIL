!< FOSSIL, mesh decimation via Quadric Error Metrics (issue #18 §1.3).

module fossil_decimate
!< FOSSIL, mesh decimation via Quadric Error Metrics (Garland & Heckbert 1997).
!<
!< Reduces a triangle mesh's facet count by iteratively collapsing edges,
!< choosing each collapse to minimize the geometric error (squared distance
!< from the collapsed vertex to all incident face planes). The algorithm:
!<
!<   1. Per-vertex 4×4 quadric `Q_v = Σ K_p` over incident face planes,
!<      where `K_p = p p^T` and `p = [a, b, c, d]` are the plane coefficients.
!<   2. Per-edge optimal collapse target `v* = argmin v^T (Q_u + Q_v) v`,
!<      found by solving the 3×3 linear system `[Q_uv]_{1:3,1:3} v* = -[Q_uv]_{1:3,4}`,
!<      with midpoint fallback when the matrix is singular.
!<   3. Per-edge cost `cost = v*^T Q_uv v*`.
!<   4. Min-heap keyed by cost. Pop cheapest, apply the patch rules (the target of an edge touching a patch rim is
!<      constrained, its cost re-evaluated lazily: the free optimum is a lower bound), run safety checks (link
!<      condition, normal-flip), collapse if accepted. Update neighbor quadrics and re-cost.
!<   5. Stop when facet count reaches `target_facets`, the heap is empty, or (optional `max_error`) the cheapest
!<      collapse costs more than max_error**2.
!<
!< Safety checks implemented:
!<   - **Link condition** (Dey, Edelsbrunner, Guha, Nekhayev 1999; issue #81):
!<     edge (a, b) collapses only if Lk(a) ∩ Lk(b) = Lk(ab), i.e. the common
!<     neighbours of a and b are the apexes of the edge's facets and no
!<     tetrahedron (a, c, d), (b, c, d) closes over it; on an open surface
!<     the links are taken with a dummy vertex joined to the boundary (two
!<     boundary vertices merge only along a boundary edge). The result stays a
!<     2-manifold of the input's topology (a tetrahedron is not collapsed
!<     into two coincident facets).
!<   - **Normal-flip rejection**: a candidate collapse that would reverse the
!<     normal of any surviving incident facet, or collapse it to zero area
!<     (issue #81: two pool vertices at one place, merged into a pinch by the
!<     re-analysis; a flat cube ended with non-manifold edges and two
!<     shells), is rejected. Without this,
!<     decimation produces inverted facets that downstream rendering and
!<     volume computation interpret backwards.
!<   - **Degenerate-facet rejection**: rejected if a surviving facet would
!<     get two equal vertices.
!<
!< Not implemented (deferred):
!<   - Boundary-preservation weights on mesh boundaries (open surfaces) by
!<     default — the original Garland-Heckbert paper adds constraint planes
!<     along them. FOSSIL inputs are typically closed solids (no boundary); the
!<     planes are added along patch rims, and along open edges only with
!<     `max_error` (issue #81).
!<   - Texture / attribute preservation — out of scope for STL.
!<
!< Pre-conditions:
!<   - Input facets must have `vertex_id(:)` populated (canonical pool ids).
!<     `surface_stl_object%adopt_facets`/`load_from_file` ensure this.
!<   - `fcon_edge(:)` populated (run `analyze` first). The decimate driver
!<     uses both for incidence queries.

use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: decimate
public :: DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS

integer(I4P), parameter :: DEC_STATUS_OK          = 0_I4P
integer(I4P), parameter :: DEC_STATUS_BAD_INPUT   = 1_I4P  !< Empty facet array, vertex_ids unset or target <= 0.
integer(I4P), parameter :: DEC_STATUS_NO_PROGRESS = 2_I4P  !< Target not reached: heap exhausted, or max_error.

! Quadric is a 4×4 symmetric matrix; store as 10-entry upper-triangle vector.
! Index convention: Q(i,j) with i<=j packs into the linear array as
!   (1,1)=1  (1,2)=2  (1,3)=3  (1,4)=4
!            (2,2)=5  (2,3)=6  (2,4)=7
!                     (3,3)=8  (3,4)=9
!                              (4,4)=10
integer(I4P), parameter :: QSIZE = 10_I4P

! Heap sentinel for "this edge is not in the heap".
integer(I4P), parameter :: NOT_IN_HEAP = 0_I4P

! Collapse target rules of an edge (patches, issues #66 and #81).
integer(I4P), parameter :: MODE_FREE   = 0_I4P !< Both interior: the quadric optimum.
integer(I4P), parameter :: MODE_AT_A   = 1_I4P !< va on a rim, vb interior: at va.
integer(I4P), parameter :: MODE_AT_B   = 2_I4P !< vb on a rim, va interior: at vb.
integer(I4P), parameter :: MODE_RIM    = 3_I4P !< A rim edge, no corner: at the cheaper endpoint.
integer(I4P), parameter :: MODE_LOCKED = 4_I4P !< Two rim vertices otherwise: never.

contains

   subroutine decimate(facet, target_facets, max_error, rim_collapses, status)
   !< Reduce `facet` to at most `target_facets` triangles via QEM edge collapse.
   !<
   !< Patches (issues #66, #81): facets keep their patch. A patch-boundary (rim) vertex stays where it is when an
   !< interior vertex merges onto it. Two rim vertices merge only along their own rim edge (the edge between two facets
   !< of different patches), onto the cheaper of its two endpoints, or onto the corner if one is a corner (a rim vertex
   !< with other than two rim edges: three or more patches meet there; two corners never merge), and never so as to
   !< empty a patch: a rim coarsens as a curve whose vertices are a subset of the input rim's, its corners kept.
   !< `rim_collapses=.false.` locks the rims (two rim vertices never merge). Each rim edge adds to its endpoints'
   !< quadrics the plane through it perpendicular to the mean normal of its two facets (Garland & Heckbert's boundary
   !< constraint), so a rim within a plane is not free to straighten.
   !<
   !< `max_error` (a distance): the cost of a collapse is a sum of squared distances to planes of the input (unit
   !< normals), so sqrt(cost) bounds the distance of the new vertex from each of those planes; when the cheapest
   !< admissible collapse has sqrt(cost) > max_error, the loop stops with DEC_STATUS_NO_PROGRESS. With max_error the open
   !< edges (one facet) get constraint planes too (through the edge, perpendicular to its facet), else a vertex of an
   !< open edge slides along its plane for free. It bounds planes, not the surface: no certificate of the deviation.
   !<
   !< On success, `facet` is reallocated to the surviving (active) facets only,
   !< with their `compute_metrix` already invoked. On failure (`status /= 0`)
   !< the input is left in a possibly partially-decimated state — caller should
   !< rebuild from the original input if recovery is needed.
   !<
   !< Allocation strategy: we work on copies of the input geometry (vertex
   !< coords + vertex incidence lists + edge list + heap), leaving the input
   !< facet array as the source of truth. At the end we materialize a new
   !< facet array from the surviving alive-facets and `move_alloc` it back.
   type(facet_object), allocatable, intent(inout)        :: facet(:)        !< Input facets; replaced by decimated set on success.
   integer(I4P),                    intent(in)           :: target_facets   !< Stop when active facet count reaches this.
   real(R8P),                       intent(in),  optional :: max_error      !< Stop when the cheapest collapse exceeds it.
   logical,                         intent(in),  optional :: rim_collapses  !< Rim vertices merge along rims (default).
   integer(I4P),                    intent(out), optional :: status         !< DEC_STATUS_*.
   integer(I4P)                                          :: nf_in           !< Input facet count.
   integer(I4P)                                          :: nv              !< Vertex count (max pool id).
   integer(I4P)                                          :: ne              !< Edge count.
   real(R8P),     allocatable                            :: vcoord(:, :)    !< Vertex coordinates (3, nv).
   real(R8P),     allocatable                            :: Q(:, :)         !< Per-vertex quadric (QSIZE, nv).
   logical,       allocatable                            :: vertex_alive(:) !< (nv).
   integer(I4P),  allocatable                            :: f_v(:, :)       !< Per-facet vertex pool ids (3, nf_in).
   logical,       allocatable                            :: facet_alive(:)  !< (nf_in).
   ! Per-vertex incidence lists, packed as (head, count) into a flat values array.
   integer(I4P),  allocatable                            :: v2f_head(:)     !< (nv) starting index in v2f_val.
   integer(I4P),  allocatable                            :: v2f_count(:)    !< (nv) live count of facets touching this vertex.
   integer(I4P),  allocatable                            :: v2f_cap(:)      !< (nv) allocated capacity.
   integer(I4P),  allocatable                            :: v2f_val(:)      !< Flat per-vertex facet-id buffer.
   ! Edge list.
   integer(I4P),  allocatable                            :: e_v(:, :)       !< (2, ne) endpoints (v_min, v_max).
   real(R8P),     allocatable                            :: e_target(:, :)  !< (3, ne) optimal collapse target.
   real(R8P),     allocatable                            :: e_cost(:)       !< (ne) cost.
   logical,       allocatable                            :: e_alive(:)      !< (ne).
   ! Per-vertex edge incidence (analogous to v2f_*).
   integer(I4P),  allocatable                            :: v2e_head(:)     !< (nv).
   integer(I4P),  allocatable                            :: v2e_count(:)    !< (nv).
   integer(I4P),  allocatable                            :: v2e_cap(:)      !< (nv).
   integer(I4P),  allocatable                            :: v2e_val(:)      !< Flat buffer.
   ! Heap state.
   integer(I4P),  allocatable                            :: heap(:)         !< Heap slots → edge ids.
   integer(I4P),  allocatable                            :: e_heappos(:)    !< (ne) position of edge in heap, or NOT_IN_HEAP.
   integer(I4P)                                          :: heap_size
   integer(I4P)                                          :: nf_alive, st
   integer(I4P)                                          :: e_top, va, vb
   logical,       allocatable                            :: v_bnd(:)        !< On a patch boundary (issue #66).
   real(R8P)                                             :: tgt(3)          !< Collapse target.
   real(R8P)                                             :: cost            !< Collapse cost of the popped edge.
   integer(I4P),  allocatable                            :: f_p(:)          !< Patch of each facet (nf_in).
   integer(I4P),  allocatable                            :: p_live(:)       !< Live facets per patch.
   integer(I4P)                                          :: f, k, mode
   logical                                               :: rims            !< Rim collapses on.

   if (present(status)) status = DEC_STATUS_OK
   if (.not. allocated(facet)) then
      if (present(status)) status = DEC_STATUS_BAD_INPUT
      return
   endif
   nf_in = size(facet, dim=1)
   if (nf_in == 0_I4P .or. target_facets <= 0_I4P) then
      if (present(status)) status = DEC_STATUS_BAD_INPUT
      return
   endif
   if (nf_in <= target_facets) return  ! nothing to do

   ! --- Stage 1: gather per-facet vertex IDs and derive vertex count. ---
   call collect_vertex_ids(facet=facet, f_v=f_v, nv=nv, st=st)
   if (st /= DEC_STATUS_OK) then
      if (present(status)) status = st
      return
   endif

   allocate(vcoord(3, nv), source=0._R8P)
   allocate(vertex_alive(nv), source=.true.)
   call collect_vertex_coords(facet=facet, f_v=f_v, vcoord=vcoord)

   allocate(facet_alive(nf_in), source=.true.)
   nf_alive = nf_in
   rims = .true. ; if (present(rim_collapses)) rims = rim_collapses
   allocate(f_p(nf_in))
   do f = 1, nf_in
      f_p(f) = max(1_I4P, facet(f)%patch)
   enddo
   allocate(p_live(maxval(f_p)), source=0_I4P)
   do f = 1, nf_in
      p_live(f_p(f)) = p_live(f_p(f)) + 1
   enddo
   ! patches (issue #66): a vertex whose facets are in two patches is on a patch boundary (a rim); an interior vertex
   ! merges onto a rim one at its place, two rim vertices only along their rim edge (issue #81, see edge_mode)
   allocate(v_bnd(nv), source=.false.)
   block
      integer(I4P), allocatable :: vp(:)
      allocate(vp(nv), source=0_I4P)
      do f = 1, nf_in
         do k = 1, 3
            if (vp(f_v(k, f)) == 0_I4P) then
               vp(f_v(k, f)) = facet(f)%patch
            elseif (vp(f_v(k, f)) /= facet(f)%patch) then
               v_bnd(f_v(k, f)) = .true.
            endif
         enddo
      enddo
   endblock

   ! --- Stage 2: build per-vertex incidence (which facets touch each vertex). ---
   call build_v2f(f_v=f_v, nv=nv, &
                  v2f_head=v2f_head, v2f_count=v2f_count, v2f_cap=v2f_cap, v2f_val=v2f_val)

   ! --- Stage 3: build per-vertex quadrics from incident face planes. ---
   allocate(Q(QSIZE, nv), source=0._R8P)
   call build_quadrics(facet=facet, f_v=f_v, vcoord=vcoord, Q=Q)
   ! constraint planes along the patch rims, and along the open edges with a cost ceiling (issue #81)
   if (any(v_bnd) .or. present(max_error)) &
      call add_constraint_planes(facet=facet, f_v=f_v, f_p=f_p, vcoord=vcoord, v2f_head=v2f_head, v2f_count=v2f_count, &
                                 v2f_val=v2f_val, rims=any(v_bnd), open=present(max_error), Q=Q)

   ! --- Stage 4: build edge list (canonical, unique). ---
   call build_edges(facet=facet, f_v=f_v, e_v=e_v, ne=ne)
   allocate(e_target(3, ne), source=0._R8P)
   allocate(e_cost(ne), source=0._R8P)
   allocate(e_alive(ne), source=.true.)
   call build_v2e(e_v=e_v, ne=ne, nv=nv, &
                  v2e_head=v2e_head, v2e_count=v2e_count, v2e_cap=v2e_cap, v2e_val=v2e_val)

   ! --- Stage 5: cost every edge. ---
   call cost_all_edges(e_v=e_v, ne=ne, vcoord=vcoord, Q=Q, e_target=e_target, e_cost=e_cost)

   ! --- Stage 6: build heap. ---
   allocate(heap(ne), source=0_I4P)
   allocate(e_heappos(ne), source=NOT_IN_HEAP)
   heap_size = 0_I4P
   call heap_build(ne=ne, e_alive=e_alive, e_cost=e_cost, heap=heap, e_heappos=e_heappos, &
                   heap_size=heap_size)

   ! --- Stage 7: main collapse loop. ---
   do
      if (nf_alive <= target_facets) exit
      if (heap_size == 0_I4P) then
         if (present(status)) status = DEC_STATUS_NO_PROGRESS
         exit
      endif
      e_top = heap_pop(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size)
      if (.not. e_alive(e_top)) cycle  ! stale entry
      va = e_v(1, e_top) ; vb = e_v(2, e_top)
      if (.not. (vertex_alive(va) .and. vertex_alive(vb))) then
         e_alive(e_top) = .false.
         cycle
      endif
      ! the stored cost is the free optimum (a lower bound); the patch rules constrain the target: re-cost lazily
      mode = edge_mode(va=va, vb=vb, v_bnd=v_bnd, rims=rims, f_v=f_v, f_p=f_p, facet_alive=facet_alive, &
                       v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val)
      if (mode == MODE_LOCKED) then
         e_alive(e_top) = .false.
         cycle
      endif
      call cost_one_edge(va=va, vb=vb, vcoord=vcoord, Q=Q, target=tgt, cost=cost, mode=mode)
      if (cost > e_cost(e_top) .or. cost < e_cost(e_top)) then ! changed (never true for a NaN)
         e_cost(e_top) = cost ; e_target(:, e_top) = tgt
         call heap_insert(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, e=e_top)
         cycle
      endif
      if (present(max_error)) then
         if (cost > max_error**2) then
            if (present(status)) status = DEC_STATUS_NO_PROGRESS
            exit
         endif
      endif
      ! Try the collapse; if rejected by safety checks, mark this edge dead
      ! and continue (its cost would not change, so no point re-queuing).
      if (.not. try_collapse(va=va, vb=vb, target=tgt, &
                              f_v=f_v, vcoord=vcoord, vertex_alive=vertex_alive, &
                              facet_alive=facet_alive, &
                              v2f_head=v2f_head, v2f_count=v2f_count, v2f_cap=v2f_cap, v2f_val=v2f_val, &
                              v2e_head=v2e_head, v2e_count=v2e_count, v2e_cap=v2e_cap, v2e_val=v2e_val, &
                              e_v=e_v, e_alive=e_alive, &
                              Q=Q, &
                              nf_alive=nf_alive, f_p=f_p, p_live=p_live, &
                              e_target=e_target, e_cost=e_cost, &
                              heap=heap, e_heappos=e_heappos, heap_size=heap_size)) then
         e_alive(e_top) = .false.
         cycle
      endif
      v_bnd(va) = v_bnd(va) .or. v_bnd(vb) ! va survives at the boundary vertex's place
      e_alive(e_top) = .false.
   enddo

   ! --- Stage 8: rebuild facet array from survivors. ---
   call materialize_output(facet=facet, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   endsubroutine decimate

   ! ===========================================================================
   ! Stage helpers
   ! ===========================================================================

   subroutine collect_vertex_ids(facet, f_v, nv, st)
   !< Read the canonical pool vertex IDs from each input facet, copy into a
   !< flat (3, nf) array, and compute the maximum pool id (= vertex count).
   !< Verifies all vertex_ids are populated (>0).
   type(facet_object),         intent(in)            :: facet(:)
   integer(I4P),  allocatable, intent(out)           :: f_v(:, :)
   integer(I4P),               intent(out)           :: nv
   integer(I4P),               intent(out)           :: st
   integer(I4P)                                      :: f, k

   st = DEC_STATUS_OK
   nv = 0_I4P
   allocate(f_v(3, size(facet, dim=1)))
   do f = 1, size(facet, dim=1)
      do k = 1, 3
         if (facet(f)%vertex_id(k) <= 0_I4P) then
            st = DEC_STATUS_BAD_INPUT
            return
         endif
         f_v(k, f) = facet(f)%vertex_id(k)
         if (f_v(k, f) > nv) nv = f_v(k, f)
      enddo
   enddo
   endsubroutine collect_vertex_ids

   subroutine collect_vertex_coords(facet, f_v, vcoord)
   !< Resolve the 3D coordinates of each vertex by reading them off the first
   !< facet that references that vertex. (All facets sharing a pool id have
   !< the same coords by construction.)
   type(facet_object), intent(in)    :: facet(:)
   integer(I4P),       intent(in)    :: f_v(:, :)
   real(R8P),          intent(inout) :: vcoord(:, :)
   logical, allocatable              :: seen(:)
   integer(I4P)                      :: f, k, v

   allocate(seen(size(vcoord, dim=2)), source=.false.)
   do f = 1, size(facet, dim=1)
      do k = 1, 3
         v = f_v(k, f)
         if (.not. seen(v)) then
            vcoord(1, v) = facet(f)%vertex(k)%x
            vcoord(2, v) = facet(f)%vertex(k)%y
            vcoord(3, v) = facet(f)%vertex(k)%z
            seen(v) = .true.
         endif
      enddo
   enddo
   endsubroutine collect_vertex_coords

   subroutine build_v2f(f_v, nv, v2f_head, v2f_count, v2f_cap, v2f_val)
   !< Build the per-vertex facet-incidence list as a flat packed array.
   !< Initial allocation gives each vertex 6 slots (typical valence on a
   !< regular triangulation is 6); per-vertex grow doubles capacity on demand
   !< via `v2f_grow_for_vertex`.
   integer(I4P),               intent(in)  :: f_v(:, :)
   integer(I4P),               intent(in)  :: nv
   integer(I4P), allocatable,  intent(out) :: v2f_head(:), v2f_count(:), v2f_cap(:), v2f_val(:)
   integer(I4P)                            :: f, k, v, total_cap
   integer(I4P), parameter                 :: INIT_CAP = 6_I4P

   allocate(v2f_head(nv), v2f_count(nv), v2f_cap(nv))
   v2f_count = 0_I4P
   v2f_cap   = INIT_CAP
   total_cap = nv * INIT_CAP
   allocate(v2f_val(total_cap), source=0_I4P)
   ! Lay out heads.
   v2f_head(1) = 1_I4P
   do v = 2, nv
      v2f_head(v) = v2f_head(v - 1) + INIT_CAP
   enddo
   ! Populate.
   do f = 1, size(f_v, dim=2)
      do k = 1, 3
         v = f_v(k, f)
         call v2x_append(head=v2f_head, count=v2f_count, cap=v2f_cap, val=v2f_val, &
                         vertex=v, item=f)
      enddo
   enddo
   endsubroutine build_v2f

   subroutine v2x_append(head, count, cap, val, vertex, item)
   !< Append `item` to vertex `vertex`'s incidence list. Grows the per-vertex
   !< slice by re-laying-out the global flat array if the local capacity is
   !< exhausted (rare; only when a vertex's valence exceeds the initial cap).
   integer(I4P), allocatable, intent(inout) :: head(:)
   integer(I4P),              intent(inout) :: count(:)
   integer(I4P),              intent(inout) :: cap(:)
   integer(I4P), allocatable, intent(inout) :: val(:)
   integer(I4P),              intent(in)    :: vertex, item
   integer(I4P), allocatable                :: new_val(:)
   integer(I4P)                             :: nv, v, new_cap, write_pos

   if (count(vertex) < cap(vertex)) then
      val(head(vertex) + count(vertex)) = item
      count(vertex) = count(vertex) + 1
      return
   endif
   ! Need to grow: double THIS vertex's slice and re-layout the entire
   ! flat array. Cheaper would be to keep per-vertex linked lists, but
   ! this only fires for high-valence vertices which are rare.
   nv = size(head)
   new_cap = 2 * cap(vertex)
   allocate(new_val(sum(cap) - cap(vertex) + new_cap), source=0_I4P)
   write_pos = 1_I4P
   do v = 1, nv
      if (count(v) > 0) new_val(write_pos:write_pos + count(v) - 1) = val(head(v):head(v) + count(v) - 1)
      head(v) = write_pos
      if (v == vertex) then
         cap(v) = new_cap
      endif
      write_pos = write_pos + cap(v)
   enddo
   call move_alloc(from=new_val, to=val)
   val(head(vertex) + count(vertex)) = item
   count(vertex) = count(vertex) + 1
   endsubroutine v2x_append

   subroutine build_quadrics(facet, f_v, vcoord, Q)
   !< For each input facet, compute its plane coefficients [a, b, c, d] and
   !< accumulate `K_p = p p^T` into the quadrics of all three of its vertices.
   type(facet_object), intent(in)    :: facet(:)
   integer(I4P),       intent(in)    :: f_v(:, :)
   real(R8P),          intent(in)    :: vcoord(:, :)
   real(R8P),          intent(inout) :: Q(:, :)
   integer(I4P)                      :: f, k, v
   real(R8P)                         :: a, b, c, dpl
   real(R8P)                         :: Kp(QSIZE)

   do f = 1, size(facet, dim=1)
      ! Plane equation: n . x = d (d = n . v1). Standard form: a x + b y + c z + d_pl = 0
      ! with [a, b, c, d_pl] = [n.x, n.y, n.z, -d].
      a   = facet(f)%normal%x
      b   = facet(f)%normal%y
      c   = facet(f)%normal%z
      dpl = -facet(f)%d
      Kp(1)  = a * a   ! (1,1)
      Kp(2)  = a * b   ! (1,2)
      Kp(3)  = a * c   ! (1,3)
      Kp(4)  = a * dpl ! (1,4)
      Kp(5)  = b * b   ! (2,2)
      Kp(6)  = b * c   ! (2,3)
      Kp(7)  = b * dpl ! (2,4)
      Kp(8)  = c * c   ! (3,3)
      Kp(9)  = c * dpl ! (3,4)
      Kp(10) = dpl * dpl ! (4,4)
      do k = 1, 3
         v = f_v(k, f)
         Q(:, v) = Q(:, v) + Kp
      enddo
   enddo
   endsubroutine build_quadrics

   subroutine build_edges(facet, f_v, e_v, ne)
   !< Build the unique edge list. Each edge is canonicalized as (min, max)
   !< vertex IDs. Walks all (3 × nf) half-edges, sorts by packed key, scans
   !< for unique keys.
   !<
   !< Same sort-and-pair approach as `build_connectivity` in surface_stl, but
   !< we don't need the facet-pair links here (those come from `fcon_edge`
   !< on the facets themselves and we use them only via per-vertex incidence).
   type(facet_object), intent(in)            :: facet(:)
   integer(I4P),       intent(in)            :: f_v(:, :)
   integer(I4P),       allocatable, intent(out) :: e_v(:, :)
   integer(I4P),                    intent(out) :: ne
   integer(I4P)                              :: f, k, h, va, vb, total_he
   integer(I4P), allocatable                 :: he_lo(:), he_hi(:), order(:)
   integer(I4P)                              :: prev_lo, prev_hi
   integer(I4P), allocatable                 :: e_lo(:), e_hi(:)

   total_he = 3 * size(f_v, dim=2)
   allocate(he_lo(total_he), he_hi(total_he))
   h = 0_I4P
   do f = 1, size(f_v, dim=2)
      do k = 1, 3
         va = f_v(k, f) ; vb = f_v(mod(k, 3) + 1, f)
         h = h + 1
         he_lo(h) = min(va, vb)
         he_hi(h) = max(va, vb)
      enddo
   enddo
   ! Sort by packed (lo, hi). Use a simple key-array sort: build (lo * (max_v+1) + hi) keys.
   call sort_he(he_lo=he_lo, he_hi=he_hi, order=order)
   ! Scan for unique runs.
   allocate(e_lo(total_he), e_hi(total_he))
   ne = 0_I4P
   prev_lo = -1_I4P ; prev_hi = -1_I4P
   do h = 1, total_he
      if (he_lo(order(h)) /= prev_lo .or. he_hi(order(h)) /= prev_hi) then
         ne = ne + 1
         e_lo(ne) = he_lo(order(h))
         e_hi(ne) = he_hi(order(h))
         prev_lo = e_lo(ne) ; prev_hi = e_hi(ne)
      endif
   enddo
   allocate(e_v(2, ne))
   e_v(1, 1:ne) = e_lo(1:ne)
   e_v(2, 1:ne) = e_hi(1:ne)
   endsubroutine build_edges

   subroutine sort_he(he_lo, he_hi, order)
   !< Sort half-edges by (he_lo, he_hi) into the permutation `order`.
   !< Simple O(n log n) merge sort would be ideal, but a small-allocation
   !< quicksort over packed I8P keys is shorter; for typical n_he < 1M this
   !< is still well within tolerance.
   integer(I4P),              intent(in)  :: he_lo(:), he_hi(:)
   integer(I4P), allocatable, intent(out) :: order(:)
   integer(I4P)                           :: n, i

   n = size(he_lo)
   allocate(order(n))
   do i = 1, n
      order(i) = i
   enddo
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=1_I4P, hi=n)
   endsubroutine sort_he

   recursive subroutine qsort_he(he_lo, he_hi, order, lo, hi)
   !< Recursive quicksort with median-of-three pivot. Sorts `order` in place
   !< by the (he_lo, he_hi) lexicographic key.
   integer(I4P), intent(in)    :: he_lo(:), he_hi(:)
   integer(I4P), intent(inout) :: order(:)
   integer(I4P), intent(in)    :: lo, hi
   integer(I4P)                :: i, j, mid, pivot_lo, pivot_hi, tmp

   if (hi - lo < 1) return
   mid = (lo + hi) / 2
   pivot_lo = he_lo(order(mid))
   pivot_hi = he_hi(order(mid))
   i = lo ; j = hi
   do
      do while (he_lo(order(i)) < pivot_lo .or. &
                (he_lo(order(i)) == pivot_lo .and. he_hi(order(i)) < pivot_hi))
         i = i + 1
      enddo
      do while (he_lo(order(j)) > pivot_lo .or. &
                (he_lo(order(j)) == pivot_lo .and. he_hi(order(j)) > pivot_hi))
         j = j - 1
      enddo
      if (i >= j) exit
      tmp = order(i) ; order(i) = order(j) ; order(j) = tmp
      i = i + 1 ; j = j - 1
   enddo
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=lo, hi=j)
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=j + 1, hi=hi)
   endsubroutine qsort_he

   subroutine build_v2e(e_v, ne, nv, v2e_head, v2e_count, v2e_cap, v2e_val)
   !< Per-vertex edge incidence, same packing as v2f.
   integer(I4P),              intent(in)  :: e_v(:, :)
   integer(I4P),              intent(in)  :: ne, nv
   integer(I4P), allocatable, intent(out) :: v2e_head(:), v2e_count(:), v2e_cap(:), v2e_val(:)
   integer(I4P)                           :: e, v
   integer(I4P), parameter                :: INIT_CAP = 6_I4P

   allocate(v2e_head(nv), v2e_count(nv), v2e_cap(nv))
   v2e_count = 0_I4P
   v2e_cap   = INIT_CAP
   allocate(v2e_val(nv * INIT_CAP), source=0_I4P)
   v2e_head(1) = 1_I4P
   do v = 2, nv
      v2e_head(v) = v2e_head(v - 1) + INIT_CAP
   enddo
   do e = 1, ne
      call v2x_append(head=v2e_head, count=v2e_count, cap=v2e_cap, val=v2e_val, &
                      vertex=e_v(1, e), item=e)
      call v2x_append(head=v2e_head, count=v2e_count, cap=v2e_cap, val=v2e_val, &
                      vertex=e_v(2, e), item=e)
   enddo
   endsubroutine build_v2e

   subroutine cost_all_edges(e_v, ne, vcoord, Q, e_target, e_cost)
   !< Initial costing of every edge. Per-edge costing is also called after
   !< each successful collapse for incident edges (`recost_edge`).
   integer(I4P), intent(in)    :: e_v(:, :)
   integer(I4P), intent(in)    :: ne
   real(R8P),    intent(in)    :: vcoord(:, :)
   real(R8P),    intent(in)    :: Q(:, :)
   real(R8P),    intent(inout) :: e_target(:, :)
   real(R8P),    intent(inout) :: e_cost(:)
   integer(I4P)                :: e

   do e = 1, ne
      call cost_one_edge(va=e_v(1, e), vb=e_v(2, e), vcoord=vcoord, Q=Q, &
                         target=e_target(:, e), cost=e_cost(e))
   enddo
   endsubroutine cost_all_edges

   subroutine cost_one_edge(va, vb, vcoord, Q, target, cost, mode)
   !< Compute the optimal collapse target and cost for edge (va, vb).
   !<
   !< Algorithm (Garland-Heckbert §4):
   !<   Q_uv = Q[va] + Q[vb]
   !<   Solve A v* = -b where A = [Q_uv]_{1:3, 1:3}, b = [Q_uv]_{1:3, 4}.
   !<   If A is singular, fall back to the midpoint of the edge.
   !<   cost = v*^T Q_uv v* (compute with v*_homog = [v*, 1]).
   integer(I4P), intent(in)  :: va, vb
   real(R8P),    intent(in)  :: vcoord(:, :)
   real(R8P),    intent(in)  :: Q(:, :)
   real(R8P),    intent(out) :: target(3)
   real(R8P),    intent(out) :: cost
   integer(I4P), intent(in), optional :: mode !< MODE_*, default MODE_FREE.
   real(R8P)                 :: Quv(QSIZE)
   real(R8P)                 :: A(3, 3), b(3), det, cost_b

   Quv = Q(:, va) + Q(:, vb)
   if (present(mode)) then
      select case (mode)
      case (MODE_AT_A)
         target = vcoord(:, va) ; cost = quadric_eval(Q=Quv, x=target) ; return
      case (MODE_AT_B)
         target = vcoord(:, vb) ; cost = quadric_eval(Q=Quv, x=target) ; return
      case (MODE_RIM)
         target = vcoord(:, va) ; cost = quadric_eval(Q=Quv, x=target)
         cost_b = quadric_eval(Q=Quv, x=vcoord(:, vb))
         if (cost_b < cost) then
            target = vcoord(:, vb) ; cost = cost_b
         endif
         return
      endselect
   endif
   ! Unpack the upper-3x3 block.
   A(1, 1) = Quv(1)  ; A(1, 2) = Quv(2)  ; A(1, 3) = Quv(3)
   A(2, 1) = Quv(2)  ; A(2, 2) = Quv(5)  ; A(2, 3) = Quv(6)
   A(3, 1) = Quv(3)  ; A(3, 2) = Quv(6)  ; A(3, 3) = Quv(8)
   b(1) = -Quv(4) ; b(2) = -Quv(7) ; b(3) = -Quv(9)
   call solve_3x3(A=A, b=b, x=target, det=det)
   if (abs(det) < 1.0e-14_R8P) then
      ! Singular: fall back to midpoint.
      target(1) = 0.5_R8P * (vcoord(1, va) + vcoord(1, vb))
      target(2) = 0.5_R8P * (vcoord(2, va) + vcoord(2, vb))
      target(3) = 0.5_R8P * (vcoord(3, va) + vcoord(3, vb))
   endif
   cost = quadric_eval(Q=Quv, x=target)
   endsubroutine cost_one_edge

   pure subroutine solve_3x3(A, b, x, det)
   !< Cramer's rule for 3×3 system. Returns the determinant; caller decides
   !< singularity threshold.
   real(R8P), intent(in)  :: A(3, 3), b(3)
   real(R8P), intent(out) :: x(3)
   real(R8P), intent(out) :: det
   real(R8P)              :: m(3, 3), inv_det

   det = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
       - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
       + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
   if (abs(det) < tiny(1._R8P)) then
      x = 0._R8P ; return
   endif
   inv_det = 1._R8P / det
   m = A
   m(:, 1) = b ; x(1) = (m(1,1) * (m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
                       - m(1,2) * (m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
                       + m(1,3) * (m(2,1)*m(3,2) - m(2,2)*m(3,1))) * inv_det
   m = A ; m(:, 2) = b
   x(2) = (m(1,1) * (m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
         - m(1,2) * (m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
         + m(1,3) * (m(2,1)*m(3,2) - m(2,2)*m(3,1))) * inv_det
   m = A ; m(:, 3) = b
   x(3) = (m(1,1) * (m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
         - m(1,2) * (m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
         + m(1,3) * (m(2,1)*m(3,2) - m(2,2)*m(3,1))) * inv_det
   endsubroutine solve_3x3

   pure function quadric_eval(Q, x) result(c)
   !< Evaluate `[x, 1]^T Q [x, 1]` for a 4×4 quadric stored as the 10-entry
   !< upper triangle.
   real(R8P), intent(in) :: Q(QSIZE), x(3)
   real(R8P)             :: c

   c =       Q(1) * x(1) * x(1) &
       + 2 * Q(2) * x(1) * x(2) &
       + 2 * Q(3) * x(1) * x(3) &
       + 2 * Q(4) * x(1)        &
       +     Q(5) * x(2) * x(2) &
       + 2 * Q(6) * x(2) * x(3) &
       + 2 * Q(7) * x(2)        &
       +     Q(8) * x(3) * x(3) &
       + 2 * Q(9) * x(3)        &
       +     Q(10)
   endfunction quadric_eval

   ! ===========================================================================
   ! Heap (min-heap with decrease-key)
   ! ===========================================================================

   subroutine heap_build(ne, e_alive, e_cost, heap, e_heappos, heap_size)
   !< Heapify all alive edges. Floyd's algorithm: insert all edges into the
   !< array in order, then sift down from the last interior node.
   integer(I4P), intent(in)    :: ne
   logical,      intent(in)    :: e_alive(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   integer(I4P), intent(inout) :: heap_size
   integer(I4P)                :: e, i

   heap_size = 0_I4P
   do e = 1, ne
      if (.not. e_alive(e)) cycle
      heap_size = heap_size + 1
      heap(heap_size) = e
      e_heappos(e) = heap_size
   enddo
   do i = heap_size / 2, 1, -1
      call heap_sift_down(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, i=i)
   enddo
   endsubroutine heap_build

   function heap_pop(heap, e_heappos, e_cost, heap_size) result(top)
   !< Remove and return the minimum-cost edge id.
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(inout) :: heap_size
   integer(I4P)                :: top

   top = heap(1)
   e_heappos(top) = NOT_IN_HEAP
   if (heap_size == 1_I4P) then
      heap_size = 0_I4P
      return
   endif
   heap(1) = heap(heap_size)
   e_heappos(heap(1)) = 1_I4P
   heap_size = heap_size - 1
   call heap_sift_down(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, i=1_I4P)
   endfunction heap_pop

   subroutine heap_decrease_key(heap, e_heappos, e_cost, heap_size, e)
   !< The cost of edge `e` has decreased; bubble it up toward the root.
   !< If `e` is not currently in the heap (e.g. it was popped already),
   !< this is a no-op — caller should re-insert via `heap_insert` instead.
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(in)    :: heap_size
   integer(I4P), intent(in)    :: e

   if (e_heappos(e) == NOT_IN_HEAP) return
   call heap_sift_up(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, &
                     i=e_heappos(e))
   endsubroutine heap_decrease_key

   subroutine heap_increase_key(heap, e_heappos, e_cost, heap_size, e)
   !< Cost of edge `e` has increased; sift it down.
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(in)    :: heap_size
   integer(I4P), intent(in)    :: e

   if (e_heappos(e) == NOT_IN_HEAP) return
   call heap_sift_down(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, &
                       i=e_heappos(e))
   endsubroutine heap_increase_key

   subroutine heap_insert(heap, e_heappos, e_cost, heap_size, e)
   !< Insert edge `e` (assumed not currently in heap).
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(inout) :: heap_size
   integer(I4P), intent(in)    :: e

   heap_size = heap_size + 1
   heap(heap_size) = e
   e_heappos(e) = heap_size
   call heap_sift_up(heap=heap, e_heappos=e_heappos, e_cost=e_cost, heap_size=heap_size, &
                     i=heap_size)
   endsubroutine heap_insert

   subroutine heap_sift_up(heap, e_heappos, e_cost, heap_size, i)
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(in)    :: heap_size
   integer(I4P), intent(in)    :: i
   integer(I4P)                :: cur, parent, tmp

   cur = i
   do while (cur > 1_I4P)
      parent = cur / 2
      if (e_cost(heap(cur)) >= e_cost(heap(parent))) exit
      tmp = heap(cur) ; heap(cur) = heap(parent) ; heap(parent) = tmp
      e_heappos(heap(cur))    = cur
      e_heappos(heap(parent)) = parent
      cur = parent
   enddo
   endsubroutine heap_sift_up

   subroutine heap_sift_down(heap, e_heappos, e_cost, heap_size, i)
   integer(I4P), intent(inout) :: heap(:), e_heappos(:)
   real(R8P),    intent(in)    :: e_cost(:)
   integer(I4P), intent(in)    :: heap_size
   integer(I4P), intent(in)    :: i
   integer(I4P)                :: cur, left, right, smallest, tmp

   cur = i
   do
      left  = 2 * cur
      right = 2 * cur + 1
      smallest = cur
      ! Nested guards (not compound `.and.`) — Fortran's `.and.` is not
      ! short-circuit, so a compound test would read heap(left/right)
      ! even when the index is out of range. -fcheck=bounds catches it.
      if (left  <= heap_size) then
         if (e_cost(heap(left )) < e_cost(heap(smallest))) smallest = left
      endif
      if (right <= heap_size) then
         if (e_cost(heap(right)) < e_cost(heap(smallest))) smallest = right
      endif
      if (smallest == cur) exit
      tmp = heap(cur) ; heap(cur) = heap(smallest) ; heap(smallest) = tmp
      e_heappos(heap(cur))      = cur
      e_heappos(heap(smallest)) = smallest
      cur = smallest
   enddo
   endsubroutine heap_sift_down

   ! ===========================================================================
   ! The collapse itself
   ! ===========================================================================

   function try_collapse(va, vb, target, &
                          f_v, vcoord, vertex_alive, facet_alive, &
                          v2f_head, v2f_count, v2f_cap, v2f_val, &
                          v2e_head, v2e_count, v2e_cap, v2e_val, &
                          e_v, e_alive, &
                          Q, &
                          nf_alive, f_p, p_live, &
                          e_target, e_cost, &
                          heap, e_heappos, heap_size) result(accepted)
   !< Attempt to collapse edge (va, vb) to position `target`. Returns .true.
   !< if the collapse passes the safety checks and was applied; .false. if
   !< rejected.
   !<
   !< Convention: we keep vertex `va` (the lower-id endpoint) and remove `vb`.
   !< The two facets incident to the edge (va, vb) are deleted; all other
   !< facets touching `vb` get their `vb` reference rewritten to `va`.
   integer(I4P),              intent(in)    :: va, vb
   real(R8P),                 intent(in)    :: target(3)
   integer(I4P),              intent(inout) :: f_v(:, :)
   real(R8P),                 intent(inout) :: vcoord(:, :)
   logical,                   intent(inout) :: vertex_alive(:)
   logical,                   intent(inout) :: facet_alive(:)
   integer(I4P), allocatable, intent(inout) :: v2f_head(:), v2f_count(:), v2f_cap(:), v2f_val(:)
   integer(I4P), allocatable, intent(inout) :: v2e_head(:), v2e_count(:), v2e_cap(:), v2e_val(:)
   integer(I4P),              intent(inout) :: e_v(:, :)
   logical,                   intent(inout) :: e_alive(:)
   real(R8P),                 intent(inout) :: Q(:, :)
   integer(I4P),              intent(inout) :: nf_alive
   integer(I4P),              intent(in)    :: f_p(:)    !< Patch of each facet.
   integer(I4P),              intent(inout) :: p_live(:) !< Live facets per patch.
   real(R8P),                 intent(inout) :: e_target(:, :)
   real(R8P),                 intent(inout) :: e_cost(:)
   integer(I4P),              intent(inout) :: heap(:), e_heappos(:)
   integer(I4P),              intent(inout) :: heap_size
   logical                                  :: accepted
   integer(I4P), allocatable                :: facets_b(:)  !< Facets touching vb (will be rewritten or deleted).
   integer(I4P), allocatable                :: edges_b(:)   !< Edges touching vb (will be rewritten).
   integer(I4P)                             :: i, j, f, e, k, w
   integer(I4P)                             :: shared_facets(2), n_shared
   real(R8P)                                :: old_xyz_a(3)
   logical                                  :: would_flip, would_dup

   accepted = .false.

   ! Collect vb's incidence lists into local copies (we'll mutate them below).
   call gather_incidence(head=v2f_head, count=v2f_count, val=v2f_val, vertex=vb, out=facets_b)
   call gather_incidence(head=v2e_head, count=v2e_count, val=v2e_val, vertex=vb, out=edges_b)

   ! Identify the (up to 2) facets shared between va and vb — those are the
   ! ones the collapse deletes. The remaining facets in facets_b need their
   ! vb slot rewritten to va.
   call find_shared_facets(va=va, vb=vb, f_v=f_v, facet_alive=facet_alive, &
                           shared_facets=shared_facets, n_shared=n_shared)

   ! --- Safety check 0: the link condition (issue #81): the collapse keeps the surface a 2-manifold (with boundary). ---
   if (.not. link_condition_holds(va=va, vb=vb, v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, &
                                  f_v=f_v, facet_alive=facet_alive, shared_facets=shared_facets, n_shared=n_shared)) return
   ! --- a patch keeps at least one facet (issue #81: a rim collapse deletes facets of two patches). ---
   if (n_shared == 2_I4P) then
      if (f_p(shared_facets(1)) == f_p(shared_facets(2))) then
         if (p_live(f_p(shared_facets(1))) <= 2_I4P) return
      endif
   endif
   do i = 1, n_shared
      if (p_live(f_p(shared_facets(i))) <= 1_I4P) return
   enddo

   ! --- Safety check 1: would the collapse flip any non-deleted facet's normal? ---
   old_xyz_a(1) = vcoord(1, va) ; old_xyz_a(2) = vcoord(2, va) ; old_xyz_a(3) = vcoord(3, va)
   would_flip = check_normal_flip(va=va, vb=vb, target=target, &
                                  facets_b=facets_b, n_b=size(facets_b), &
                                  shared_facets=shared_facets, n_shared=n_shared, &
                                  f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   if (would_flip) return
   ! Also check facets touching va (only) — same risk if va moves.
   would_flip = check_normal_flip_around_va(va=va, vb=vb, target=target, &
                                            v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, &
                                            shared_facets=shared_facets, n_shared=n_shared, &
                                            f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   if (would_flip) return

   ! --- Safety check 2: would any facet become degenerate (two vertices equal)? ---
   would_dup = .false.
   do i = 1, size(facets_b)
      f = facets_b(i)
      if (.not. facet_alive(f)) cycle
      if (any(shared_facets(1:n_shared) == f)) cycle
      ! After rewrite vb→va, count how many slots now equal va.
      k = 0_I4P
      do j = 1, 3
         if (f_v(j, f) == va .or. f_v(j, f) == vb) k = k + 1
      enddo
      if (k >= 2) then
         would_dup = .true. ; exit
      endif
   enddo
   if (would_dup) return


   ! --- Apply the collapse. ---
   ! 1. Move va to target.
   vcoord(1, va) = target(1) ; vcoord(2, va) = target(2) ; vcoord(3, va) = target(3)
   ! 2. Update Q[va] := Q[va] + Q[vb].
   Q(:, va) = Q(:, va) + Q(:, vb)
   ! 3. For every facet touching vb that's NOT a shared facet: rewrite the vb slot to va.
   do i = 1, size(facets_b)
      f = facets_b(i)
      if (.not. facet_alive(f)) cycle
      if (any(shared_facets(1:n_shared) == f)) cycle
      do j = 1, 3
         if (f_v(j, f) == vb) f_v(j, f) = va
      enddo
      ! Append this facet to va's incidence list.
      call v2x_append(head=v2f_head, count=v2f_count, cap=v2f_cap, val=v2f_val, &
                      vertex=va, item=f)
   enddo
   ! 4. Delete shared facets.
   do i = 1, n_shared
      facet_alive(shared_facets(i)) = .false.
      nf_alive = nf_alive - 1
      p_live(f_p(shared_facets(i))) = p_live(f_p(shared_facets(i))) - 1
   enddo
   ! 5. Mark vb as not-alive.
   vertex_alive(vb) = .false.
   ! 6. Rewrite vb's edges → va, possibly merging with existing va-edges.
   do i = 1, size(edges_b)
      e = edges_b(i)
      if (.not. e_alive(e)) cycle
      ! What's the OTHER endpoint of this edge?
      if (e_v(1, e) == vb) then
         w = e_v(2, e)
      else
         w = e_v(1, e)
      endif
      if (w == va) then
         ! This was the collapsed edge itself — kill it.
         e_alive(e) = .false.
         cycle
      endif
      ! Does an edge (va, w) already exist?
      if (find_edge(va=va, vw=w, v2e_head=v2e_head, v2e_count=v2e_count, v2e_val=v2e_val, &
                    e_v=e_v, e_alive=e_alive) /= 0_I4P) then
         ! Merge: kill this vb-side edge.
         e_alive(e) = .false.
      else
         ! Rewrite the vb slot to va; canonicalize endpoints; append to va's edge list.
         if (e_v(1, e) == vb) then
            e_v(1, e) = va
         else
            e_v(2, e) = va
         endif
         if (e_v(1, e) > e_v(2, e)) then
            k = e_v(1, e) ; e_v(1, e) = e_v(2, e) ; e_v(2, e) = k
         endif
         call v2x_append(head=v2e_head, count=v2e_count, cap=v2e_cap, val=v2e_val, &
                         vertex=va, item=e)
      endif
   enddo
   ! 7. Re-cost all edges incident to va and bubble in the heap.
   call recost_va_edges(va=va, &
                        v2e_head=v2e_head, v2e_count=v2e_count, v2e_val=v2e_val, &
                        e_v=e_v, e_alive=e_alive, vcoord=vcoord, Q=Q, &
                        e_target=e_target, e_cost=e_cost, &
                        heap=heap, e_heappos=e_heappos, heap_size=heap_size)
   accepted = .true.
   endfunction try_collapse

   subroutine gather_incidence(head, count, val, vertex, out)
   !< Copy vertex `vertex`'s current incidence list into a fresh allocatable
   !< array `out`. Used so the caller can mutate the underlying buffer
   !< (append entries on `va`'s list) while iterating over `vb`'s list.
   integer(I4P), allocatable, intent(in)  :: head(:), val(:)
   integer(I4P),              intent(in)  :: count(:)
   integer(I4P),              intent(in)  :: vertex
   integer(I4P), allocatable, intent(out) :: out(:)

   allocate(out(count(vertex)))
   if (count(vertex) > 0) out = val(head(vertex):head(vertex) + count(vertex) - 1)
   endsubroutine gather_incidence

   subroutine find_shared_facets(va, vb, f_v, facet_alive, shared_facets, n_shared)
   !< Find the (≤ 2) live facets that contain both va and vb.
   integer(I4P), intent(in)  :: va, vb
   integer(I4P), intent(in)  :: f_v(:, :)
   logical,      intent(in)  :: facet_alive(:)
   integer(I4P), intent(out) :: shared_facets(2)
   integer(I4P), intent(out) :: n_shared
   integer(I4P)              :: f, has_a, has_b, j

   n_shared = 0_I4P
   shared_facets = 0_I4P
   do f = 1, size(f_v, dim=2)
      if (.not. facet_alive(f)) cycle
      has_a = 0_I4P ; has_b = 0_I4P
      do j = 1, 3
         if (f_v(j, f) == va) has_a = 1_I4P
         if (f_v(j, f) == vb) has_b = 1_I4P
      enddo
      if (has_a == 1_I4P .and. has_b == 1_I4P) then
         n_shared = n_shared + 1
         if (n_shared > 2_I4P) return  ! non-manifold; caller handles
         shared_facets(n_shared) = f
      endif
   enddo
   endsubroutine find_shared_facets

   function check_normal_flip(va, vb, target, facets_b, n_b, shared_facets, n_shared, &
                              f_v, vcoord, facet_alive) result(would_flip)
   !< For every live facet incident to vb (other than the shared ones which
   !< get deleted), simulate the vb→va rewrite + va move and check whether
   !< the new normal flips relative to the old.
   integer(I4P), intent(in) :: va, vb
   real(R8P),    intent(in) :: target(3)
   integer(I4P), intent(in) :: facets_b(:)
   integer(I4P), intent(in) :: n_b
   integer(I4P), intent(in) :: shared_facets(2)
   integer(I4P), intent(in) :: n_shared
   integer(I4P), intent(in) :: f_v(:, :)
   real(R8P),    intent(in) :: vcoord(:, :)
   logical,      intent(in) :: facet_alive(:)
   logical                  :: would_flip
   integer(I4P)             :: i, f, j, vrt(3)
   real(R8P)                :: p1(3), p2(3), p3(3), n_old(3), n_new(3), p_new(3)

   would_flip = .false.
   do i = 1, n_b
      f = facets_b(i)
      if (.not. facet_alive(f)) cycle
      if (any(shared_facets(1:n_shared) == f)) cycle
      ! Old normal:
      do j = 1, 3
         vrt(j) = f_v(j, f)
      enddo
      p1 = vcoord(:, vrt(1)) ; p2 = vcoord(:, vrt(2)) ; p3 = vcoord(:, vrt(3))
      n_old = cross(p2 - p1, p3 - p1)
      ! New positions: vb is replaced by target (since va moves AND vb→va,
      ! the "old vb" slot becomes the new va = target).
      do j = 1, 3
         if (vrt(j) == vb) then
            p_new = target
         else if (vrt(j) == va) then
            p_new = target
         else
            p_new = vcoord(:, vrt(j))
         endif
         if (j == 1) p1 = p_new
         if (j == 2) p2 = p_new
         if (j == 3) p3 = p_new
      enddo
      n_new = cross(p2 - p1, p3 - p1)
      ! a facet collapsed to zero area counts as flipped: it would put two pool vertices at one place, which the
      ! re-analysis merges into a pinch (issue #81: cube-01 to 36 facets, 2 non-manifold edges, 2 shells)
      if (dot_product(n_old, n_new) < 0._R8P .or. norm2(n_new) <= 1.e-12_R8P * norm2(n_old) .or. &
          all(n_new == 0._R8P)) then
         would_flip = .true. ; return
      endif
   enddo
   endfunction check_normal_flip

   function check_normal_flip_around_va(va, vb, target, &
                                         v2f_head, v2f_count, v2f_val, &
                                         shared_facets, n_shared, &
                                         f_v, vcoord, facet_alive) result(would_flip)
   !< Same as check_normal_flip but for facets touching only va (not vb).
   !< These get unaffected by the vb rewrite but DO see va move to target.
   integer(I4P),              intent(in) :: va, vb
   real(R8P),                 intent(in) :: target(3)
   integer(I4P), allocatable, intent(in) :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in) :: v2f_count(:)
   integer(I4P),              intent(in) :: shared_facets(2)
   integer(I4P),              intent(in) :: n_shared
   integer(I4P),              intent(in) :: f_v(:, :)
   real(R8P),                 intent(in) :: vcoord(:, :)
   logical,                   intent(in) :: facet_alive(:)
   logical                               :: would_flip
   integer(I4P)                          :: i, f, j, has_b
   real(R8P)                             :: p1(3), p2(3), p3(3), n_old(3), n_new(3), p_new(3)

   would_flip = .false.
   do i = 1, v2f_count(va)
      f = v2f_val(v2f_head(va) + i - 1)
      if (.not. facet_alive(f)) cycle
      if (any(shared_facets(1:n_shared) == f)) cycle
      ! Skip facets that ALSO touch vb — those were checked already.
      has_b = 0_I4P
      do j = 1, 3
         if (f_v(j, f) == vb) has_b = 1_I4P
      enddo
      if (has_b == 1_I4P) cycle
      ! Old normal.
      p1 = vcoord(:, f_v(1, f))
      p2 = vcoord(:, f_v(2, f))
      p3 = vcoord(:, f_v(3, f))
      n_old = cross(p2 - p1, p3 - p1)
      ! New: va → target.
      do j = 1, 3
         if (f_v(j, f) == va) then
            p_new = target
         else
            p_new = vcoord(:, f_v(j, f))
         endif
         if (j == 1) p1 = p_new
         if (j == 2) p2 = p_new
         if (j == 3) p3 = p_new
      enddo
      n_new = cross(p2 - p1, p3 - p1)
      ! a facet collapsed to zero area counts as flipped: it would put two pool vertices at one place, which the
      ! re-analysis merges into a pinch (issue #81: cube-01 to 36 facets, 2 non-manifold edges, 2 shells)
      if (dot_product(n_old, n_new) < 0._R8P .or. norm2(n_new) <= 1.e-12_R8P * norm2(n_old) .or. &
          all(n_new == 0._R8P)) then
         would_flip = .true. ; return
      endif
   enddo
   endfunction check_normal_flip_around_va

   pure function cross(a, b) result(c)
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: c(3)
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
   endfunction cross

   function link_condition_holds(va, vb, v2f_head, v2f_count, v2f_val, f_v, facet_alive, &
                                 shared_facets, n_shared) result(holds)
   !< Link condition (Dey, Edelsbrunner, Guha, Nekhayev 1999) for the collapse of edge (va, vb): the collapse keeps a
   !< 2-manifold a 2-manifold of the same topology iff Lk(va) ∩ Lk(vb) = Lk(va vb). A surface with boundary is coned
   !< to a dummy vertex omega joined to every boundary vertex (a boundary edge has one incident facet), so that:
   !<   - the edge has one or two incident facets, with distinct apexes (else no collapse);
   !<   - every common neighbour of va and vb is an apex of the edge (vertex part; omega: va and vb are both on the
   !<     boundary only if (va, vb) is a boundary edge);
   !<   - no edge is in both links: not both facets (va, c, d) and (vb, c, d), c and d the apexes (a tetrahedron), and
   !<     no vertex u with (va, u) and (vb, u) both boundary edges (the edge (u, omega)).
   integer(I4P),              intent(in) :: va, vb
   integer(I4P), allocatable, intent(in) :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in) :: v2f_count(:)
   integer(I4P),              intent(in) :: f_v(:, :)
   logical,                   intent(in) :: facet_alive(:)
   integer(I4P),              intent(in) :: shared_facets(2)
   integer(I4P),              intent(in) :: n_shared
   logical                               :: holds
   integer(I4P), allocatable             :: nb_a(:), cnt_a(:), nb_b(:), cnt_b(:) !< 1-rings and facets per neighbour.
   integer(I4P)                          :: n_a, n_b, apex(2), i, j, s, k
   logical                               :: bnd_a, bnd_b

   holds = .false.
   if (n_shared < 1_I4P .or. n_shared > 2_I4P) return
   apex = 0_I4P
   do s = 1, n_shared
      do k = 1, 3
         if (f_v(k, shared_facets(s)) /= va .and. f_v(k, shared_facets(s)) /= vb) apex(s) = f_v(k, shared_facets(s))
      enddo
   enddo
   if (n_shared == 2_I4P .and. apex(1) == apex(2)) return
   call vertex_ring(v=va, v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, f_v=f_v, facet_alive=facet_alive, &
                    nb=nb_a, cnt=cnt_a, n=n_a)
   call vertex_ring(v=vb, v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, f_v=f_v, facet_alive=facet_alive, &
                    nb=nb_b, cnt=cnt_b, n=n_b)
   bnd_a = any(cnt_a(1:n_a) == 1_I4P)
   bnd_b = any(cnt_b(1:n_b) == 1_I4P)
   if (n_shared == 2_I4P .and. bnd_a .and. bnd_b) return
   do i = 1, n_a
      if (nb_a(i) == vb) cycle
      do j = 1, n_b
         if (nb_b(j) /= nb_a(i)) cycle
         if (all(apex(1:n_shared) /= nb_a(i))) return
         if (cnt_a(i) == 1_I4P .and. cnt_b(j) == 1_I4P) return
      enddo
   enddo
   if (n_shared == 2_I4P) then
      if (has_facet(v=va, c=apex(1), d=apex(2), v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, f_v=f_v, &
                    facet_alive=facet_alive) .and. &
          has_facet(v=vb, c=apex(1), d=apex(2), v2f_head=v2f_head, v2f_count=v2f_count, v2f_val=v2f_val, f_v=f_v, &
                    facet_alive=facet_alive)) return
   endif
   holds = .true.
   endfunction link_condition_holds

   subroutine vertex_ring(v, v2f_head, v2f_count, v2f_val, f_v, facet_alive, nb, cnt, n)
   !< The 1-ring of vertex v: its n neighbours nb(1:n) and the number of live facets cnt(1:n) on each edge (v, nb).
   integer(I4P),              intent(in)  :: v
   integer(I4P), allocatable, intent(in)  :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in)  :: v2f_count(:)
   integer(I4P),              intent(in)  :: f_v(:, :)
   logical,                   intent(in)  :: facet_alive(:)
   integer(I4P), allocatable, intent(out) :: nb(:), cnt(:)
   integer(I4P),              intent(out) :: n
   integer(I4P)                           :: i, f, k, w, j

   allocate(nb(2 * v2f_count(v)), cnt(2 * v2f_count(v)))
   n = 0_I4P
   do i = 1, v2f_count(v)
      f = v2f_val(v2f_head(v) + i - 1)
      if (.not. facet_alive(f)) cycle
      if (all(f_v(:, f) /= v)) cycle
      do k = 1, 3
         w = f_v(k, f)
         if (w == v) cycle
         do j = 1, n
            if (nb(j) == w) exit
         enddo
         if (j > n) then
            n = n + 1 ; nb(n) = w ; cnt(n) = 0_I4P
         endif
         cnt(j) = cnt(j) + 1
      enddo
   enddo
   endsubroutine vertex_ring

   function has_facet(v, c, d, v2f_head, v2f_count, v2f_val, f_v, facet_alive) result(found)
   !< True if a live facet has the vertices v, c and d.
   integer(I4P),              intent(in) :: v, c, d
   integer(I4P), allocatable, intent(in) :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in) :: v2f_count(:)
   integer(I4P),              intent(in) :: f_v(:, :)
   logical,                   intent(in) :: facet_alive(:)
   logical                               :: found
   integer(I4P)                          :: i, f

   found = .false.
   do i = 1, v2f_count(v)
      f = v2f_val(v2f_head(v) + i - 1)
      if (.not. facet_alive(f)) cycle
      if (any(f_v(:, f) == v) .and. any(f_v(:, f) == c) .and. any(f_v(:, f) == d)) then
         found = .true. ; return
      endif
   enddo
   endfunction has_facet

   function find_edge(va, vw, v2e_head, v2e_count, v2e_val, e_v, e_alive) result(e_found)
   !< Look up an existing alive edge (va, vw) in va's incidence list.
   !< Returns 0 if none.
   integer(I4P),              intent(in) :: va, vw
   integer(I4P), allocatable, intent(in) :: v2e_head(:), v2e_val(:)
   integer(I4P),              intent(in) :: v2e_count(:)
   integer(I4P),              intent(in) :: e_v(:, :)
   logical,                   intent(in) :: e_alive(:)
   integer(I4P)                          :: e_found
   integer(I4P)                          :: i, e

   e_found = 0_I4P
   do i = 1, v2e_count(va)
      e = v2e_val(v2e_head(va) + i - 1)
      if (.not. e_alive(e)) cycle
      if ((e_v(1, e) == va .and. e_v(2, e) == vw) .or. &
          (e_v(2, e) == va .and. e_v(1, e) == vw)) then
         e_found = e ; return
      endif
   enddo
   endfunction find_edge

   subroutine recost_va_edges(va, v2e_head, v2e_count, v2e_val, e_v, e_alive, vcoord, Q, &
                              e_target, e_cost, heap, e_heappos, heap_size)
   !< Re-cost every alive edge incident to vertex va and update its heap entry.
   integer(I4P),              intent(in)    :: va
   integer(I4P), allocatable, intent(in)    :: v2e_head(:), v2e_val(:)
   integer(I4P),              intent(in)    :: v2e_count(:)
   integer(I4P),              intent(in)    :: e_v(:, :)
   logical,                   intent(in)    :: e_alive(:)
   real(R8P),                 intent(in)    :: vcoord(:, :)
   real(R8P),                 intent(in)    :: Q(:, :)
   real(R8P),                 intent(inout) :: e_target(:, :), e_cost(:)
   integer(I4P),              intent(inout) :: heap(:), e_heappos(:)
   integer(I4P),              intent(inout) :: heap_size
   integer(I4P)                             :: i, e
   real(R8P)                                :: old_cost, new_cost, new_target(3)

   do i = 1, v2e_count(va)
      e = v2e_val(v2e_head(va) + i - 1)
      if (.not. e_alive(e)) cycle
      old_cost = e_cost(e)
      call cost_one_edge(va=e_v(1, e), vb=e_v(2, e), vcoord=vcoord, Q=Q, &
                         target=new_target, cost=new_cost)
      e_target(:, e) = new_target
      e_cost(e)      = new_cost
      if (e_heappos(e) == NOT_IN_HEAP) then
         call heap_insert(heap=heap, e_heappos=e_heappos, e_cost=e_cost, &
                          heap_size=heap_size, e=e)
      else if (new_cost < old_cost) then
         call heap_decrease_key(heap=heap, e_heappos=e_heappos, e_cost=e_cost, &
                                heap_size=heap_size, e=e)
      else
         call heap_increase_key(heap=heap, e_heappos=e_heappos, e_cost=e_cost, &
                                heap_size=heap_size, e=e)
      endif
   enddo
   endsubroutine recost_va_edges

   ! ===========================================================================
   ! Patch rims (issues #66, #81)
   ! ===========================================================================

   function edge_mode(va, vb, v_bnd, rims, f_v, f_p, facet_alive, v2f_head, v2f_count, v2f_val) result(mode)
   !< The target rule of edge (va, vb): free if both are interior, at the rim vertex if one is on a rim; two rim
   !< vertices merge only if rims are on and (va, vb) is a rim edge (its facets in two patches): at the cheaper endpoint
   !< (MODE_RIM) if neither is a corner (a rim vertex with other than two rim edges), at the corner if one is (corners
   !< stay); else MODE_LOCKED.
   integer(I4P),              intent(in) :: va, vb
   logical,                   intent(in) :: v_bnd(:)
   logical,                   intent(in) :: rims
   integer(I4P),              intent(in) :: f_v(:, :), f_p(:)
   logical,                   intent(in) :: facet_alive(:)
   integer(I4P), allocatable, intent(in) :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in) :: v2f_count(:)
   integer(I4P)                          :: mode
   integer(I4P)                          :: da, db !< rim_degree of va and vb.

   if (.not. v_bnd(va)) then
      mode = MODE_FREE ; if (v_bnd(vb)) mode = MODE_AT_B
      return
   endif
   if (.not. v_bnd(vb)) then
      mode = MODE_AT_A ; return
   endif
   mode = MODE_LOCKED
   if (.not. rims) return
   da = rim_degree(v=va, w=vb, f_v=f_v, f_p=f_p, facet_alive=facet_alive, v2f_head=v2f_head, v2f_count=v2f_count, &
                   v2f_val=v2f_val)
   if (da >= 0_I4P) return ! not a rim edge
   db = rim_degree(v=vb, w=va, f_v=f_v, f_p=f_p, facet_alive=facet_alive, v2f_head=v2f_head, v2f_count=v2f_count, &
                   v2f_val=v2f_val)
   if (da == -2_I4P .and. db == -2_I4P) then
      mode = MODE_RIM
   elseif (da == -2_I4P) then
      mode = MODE_AT_B ! vb is a corner: it stays
   elseif (db == -2_I4P) then
      mode = MODE_AT_A
   endif
   endfunction edge_mode

   function rim_degree(v, w, f_v, f_p, facet_alive, v2f_head, v2f_count, v2f_val) result(deg)
   !< The number of rim edges at v (edges whose live facets are in two or more patches), negated if (v, w) is one of
   !< them: -2 means v is on a rim, not a corner, and (v, w) is a rim edge.
   integer(I4P),              intent(in) :: v, w
   integer(I4P),              intent(in) :: f_v(:, :), f_p(:)
   logical,                   intent(in) :: facet_alive(:)
   integer(I4P), allocatable, intent(in) :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in) :: v2f_count(:)
   integer(I4P)                          :: deg
   integer(I4P), allocatable             :: nb(:), pat(:)
   logical,      allocatable             :: rim(:)
   integer(I4P)                          :: i, f, k, u, j, n

   allocate(nb(2 * v2f_count(v)), pat(2 * v2f_count(v)), rim(2 * v2f_count(v)))
   n = 0_I4P
   do i = 1, v2f_count(v)
      f = v2f_val(v2f_head(v) + i - 1)
      if (.not. facet_alive(f)) cycle
      if (all(f_v(:, f) /= v)) cycle
      do k = 1, 3
         u = f_v(k, f)
         if (u == v) cycle
         do j = 1, n
            if (nb(j) == u) exit
         enddo
         if (j > n) then
            n = n + 1 ; nb(n) = u ; pat(n) = f_p(f) ; rim(n) = .false.
         elseif (pat(j) /= f_p(f)) then
            rim(j) = .true.
         endif
      enddo
   enddo
   deg = count(rim(1:n))
   do j = 1, n
      if (nb(j) == w .and. rim(j)) deg = -deg
   enddo
   endfunction rim_degree

   subroutine add_constraint_planes(facet, f_v, f_p, vcoord, v2f_head, v2f_count, v2f_val, rims, open, Q)
   !< The boundary constraints of Garland & Heckbert 1997, weight 1 like the facet planes, added to the quadrics of the
   !< edge's ends: for every rim edge (p, r) (two facets f, g of different patches; if `rims`) the plane through it with
   !< normal (r - p) x (n_f + n_g); for every open edge (one facet f; if `open`) the plane with normal (r - p) x n_f.
   type(facet_object),        intent(in)    :: facet(:)
   integer(I4P),              intent(in)    :: f_v(:, :), f_p(:)
   real(R8P),                 intent(in)    :: vcoord(:, :)
   integer(I4P), allocatable, intent(in)    :: v2f_head(:), v2f_val(:)
   integer(I4P),              intent(in)    :: v2f_count(:)
   logical,                   intent(in)    :: rims, open
   real(R8P),                 intent(inout) :: Q(:, :)
   integer(I4P)                             :: f, g, k, i, p, r
   real(R8P)                                :: nf(3)
   logical                                  :: shared

   do f = 1, size(f_v, dim=2)
      nf = [facet(f)%normal%x, facet(f)%normal%y, facet(f)%normal%z]
      do k = 1, 3
         p = f_v(k, f) ; r = f_v(mod(k, 3) + 1, f)
         shared = .false.
         do i = 1, v2f_count(p)
            g = v2f_val(v2f_head(p) + i - 1)
            if (g == f .or. all(f_v(:, g) /= r)) cycle
            shared = .true.
            if (.not. rims .or. g < f .or. f_p(g) == f_p(f)) cycle ! a rim edge once, from its lower facet
            call add_edge_plane(p=p, r=r, n=nf + [facet(g)%normal%x, facet(g)%normal%y, facet(g)%normal%z], &
                                n_alt=nf, vcoord=vcoord, Q=Q)
         enddo
         if (open .and. .not. shared) call add_edge_plane(p=p, r=r, n=nf, n_alt=nf, vcoord=vcoord, Q=Q)
      enddo
   enddo
   endsubroutine add_constraint_planes

   subroutine add_edge_plane(p, r, n, n_alt, vcoord, Q)
   !< Add the plane through the edge (p, r) with normal (r - p) x n (n_alt if n vanishes), unit, to Q(:, p), Q(:, r).
   integer(I4P), intent(in)    :: p, r
   real(R8P),    intent(in)    :: n(3), n_alt(3)
   real(R8P),    intent(in)    :: vcoord(:, :)
   real(R8P),    intent(inout) :: Q(:, :)
   real(R8P)                   :: m(3), dpl, Kp(QSIZE)

   if (norm2(n) > 1.e-12_R8P) then
      m = cross(vcoord(:, r) - vcoord(:, p), n)
   else
      m = cross(vcoord(:, r) - vcoord(:, p), n_alt)
   endif
   if (norm2(m) <= 0._R8P) return
   m = m / norm2(m)
   dpl = -dot_product(m, vcoord(:, p))
   Kp = [m(1) * m(1), m(1) * m(2), m(1) * m(3), m(1) * dpl, m(2) * m(2), m(2) * m(3), m(2) * dpl, &
         m(3) * m(3), m(3) * dpl, dpl * dpl]
   Q(:, p) = Q(:, p) + Kp
   Q(:, r) = Q(:, r) + Kp
   endsubroutine add_edge_plane

   subroutine materialize_output(facet, f_v, vcoord, facet_alive)
   !< Rebuild the public facet array from surviving (alive) facets, reading
   !< the (possibly moved) vertex positions from `vcoord`.
   type(facet_object), allocatable, intent(inout) :: facet(:)
   integer(I4P),                    intent(in)    :: f_v(:, :)
   real(R8P),                       intent(in)    :: vcoord(:, :)
   logical,                         intent(in)    :: facet_alive(:)
   type(facet_object), allocatable                :: out(:)
   integer(I4P)                                   :: f, n_out, k, j

   n_out = count(facet_alive)
   allocate(out(n_out))
   k = 0_I4P
   do f = 1, size(facet, dim=1)
      if (.not. facet_alive(f)) cycle
      k = k + 1
      do j = 1, 3
         out(k)%vertex(j)%x = vcoord(1, f_v(j, f))
         out(k)%vertex(j)%y = vcoord(2, f_v(j, f))
         out(k)%vertex(j)%z = vcoord(3, f_v(j, f))
      enddo
      ! Reset connectivity / pool ids — caller's adopt_facets re-runs analyze.
      out(k)%vertex_id = 0_I4P
      out(k)%fcon_edge = 0_I4P
      out(k)%patch = facet(f)%patch ! surviving facets keep their patch (issue #66)
      call out(k)%compute_metrix
   enddo
   call move_alloc(from=out, to=facet)
   endsubroutine materialize_output

endmodule fossil_decimate
