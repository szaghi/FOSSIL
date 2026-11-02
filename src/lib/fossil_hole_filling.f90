!< FOSSIL, hole filling: boundary loops and their triangulation (issue #31).
module fossil_hole_filling
!< Close small holes of a triangulated surface.
!<
!< A hole is a closed loop of boundary edges: facet edges with no neighbour (`fcon_edge == 0`
!< after `analyze`). Loops are built on pool vertex ids (welded vertices), so coincident
!< coordinates count as one vertex.
!<
!< Orientation: a boundary edge of facet f runs `vertex(e) -> vertex(e+1)` in f's winding. A
!< patch triangle sharing that edge must run it the other way to keep the surface consistently
!< oriented, so loops are stored in patch order, each edge reversed with respect to its owner
!< facet. A triangle (v_i, v_m, v_k) with i < m < k along the loop is then oriented like the
!< surrounding mesh.
!<
!< Pinch vertices (a vertex on more than one hole, or met twice by one boundary) are handled
!< without heuristics: boundary edges are chained greedily into closed walks, and every closed
!< walk is split into simple loops at its repeated vertices.
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: boundary_loop_t
public :: extract_boundary_loops

type :: boundary_loop_t
   !< One closed boundary loop, in patch orientation.
   integer(I4P),     allocatable :: vid(:)   !< Pool vertex ids; edge k joins vid(k) -> vid(k+1) (cyclic).
   type(vector_R8P), allocatable :: x(:)     !< Coordinates of vid(:).
   integer(I4P),     allocatable :: owner(:) !< Mesh facet owning edge k.
endtype boundary_loop_t

contains
   subroutine extract_boundary_loops(facet, loops, n_open)
   !< Extract every closed boundary loop of the facet set.
   !<
   !< Requires pool vertex ids (`vertex_id`, set by `analyze`) and connectivity (`fcon_edge`).
   !< Boundary edges that cannot be chained into a closed walk (inconsistent orientation or a
   !< non-manifold boundary) are counted in `n_open` and left out.
   type(facet_object),                 intent(in)  :: facet(:)  !< Facets (analyzed).
   type(boundary_loop_t), allocatable, intent(out) :: loops(:)  !< Closed simple loops.
   integer(I4P),                       intent(out) :: n_open    !< Boundary edges left in open chains.
   integer(I4P),     allocatable :: e_start(:), e_end(:), e_owner(:) !< Directed boundary edges (patch order).
   type(vector_R8P), allocatable :: e_xs(:)                        !< Start coordinates of each edge.
   integer(I4P),     allocatable :: head(:), next(:)               !< Per-start-vertex linked lists of edges.
   logical,          allocatable :: used(:)                        !< Edge already walked.
   integer(I4P),     allocatable :: walk(:)                        !< Edge indices of the current walk.
   type(boundary_loop_t), allocatable :: buf(:), tmp(:)            !< Loop accumulator.
   integer(I4P) :: n_edges, f, e, k, v_max, e0, cur, nxt, n_walk, n_loops

   n_open = 0_I4P
   ! 1. Directed boundary edges, reversed to patch orientation.
   n_edges = 0_I4P
   do f = 1, size(facet)
      n_edges = n_edges + count(facet(f)%fcon_edge(:) == 0_I4P)
   enddo
   allocate(loops(0))
   if (n_edges == 0_I4P) return
   allocate(e_start(n_edges), e_end(n_edges), e_owner(n_edges), e_xs(n_edges))
   k = 0_I4P
   do f = 1, size(facet)
      do e = 1, 3
         if (facet(f)%fcon_edge(e) /= 0_I4P) cycle
         k = k + 1_I4P
         e_start(k) = facet(f)%vertex_id(mod(e, 3) + 1)
         e_end(k)   = facet(f)%vertex_id(e)
         e_owner(k) = f
         e_xs(k)    = facet(f)%vertex(mod(e, 3) + 1)
      enddo
   enddo
   if (any(e_start <= 0_I4P) .or. any(e_end <= 0_I4P)) then
      n_open = n_edges   ! no vertex pool: nothing can be chained
      return
   endif
   ! 2. Linked lists of outgoing edges per start vertex.
   v_max = max(maxval(e_start), maxval(e_end))
   allocate(head(v_max), next(n_edges), used(n_edges), walk(n_edges))
   head = 0_I4P ; used = .false.
   do k = n_edges, 1, -1
      next(k) = head(e_start(k))
      head(e_start(k)) = k
   enddo
   ! 3. Greedy closed walks, each split into simple loops.
   allocate(buf(8)) ; n_loops = 0_I4P
   do e0 = 1, n_edges
      if (used(e0)) cycle
      n_walk = 1_I4P ; walk(1) = e0 ; used(e0) = .true. ; cur = e0
      do while (e_end(cur) /= e_start(e0))
         nxt = head(e_end(cur))
         do while (nxt > 0_I4P)
            if (.not. used(nxt)) exit
            nxt = next(nxt)
         enddo
         if (nxt == 0_I4P) exit          ! stuck: open chain
         n_walk = n_walk + 1_I4P ; walk(n_walk) = nxt ; used(nxt) = .true. ; cur = nxt
      enddo
      if (e_end(cur) /= e_start(e0)) then
         n_open = n_open + n_walk
         cycle
      endif
      call split_walk(walk(1:n_walk))
   enddo
   call move_alloc(from=buf, to=tmp)
   deallocate(loops)               ! the empty array allocated for the early returns
   allocate(loops(n_loops))
   do k = 1, n_loops
      call move_alloc(from=tmp(k)%vid,   to=loops(k)%vid)
      call move_alloc(from=tmp(k)%x,     to=loops(k)%x)
      call move_alloc(from=tmp(k)%owner, to=loops(k)%owner)
   enddo
   contains
      subroutine split_walk(w)
      !< Split a closed walk (edge indices) into simple loops at repeated vertices.
      integer(I4P), intent(in) :: w(:)
      integer(I4P)             :: stack(size(w)), sp, i, j, v

      sp = 0_I4P
      do i = 1, size(w)
         v = e_start(w(i))
         ! If v is already on the stack, the edges since its occurrence form a simple loop.
         do j = sp, 1, -1
            if (e_start(stack(j)) == v) exit
         enddo
         if (j >= 1) then
            call emit(stack(j:sp))
            sp = j - 1_I4P
         endif
         sp = sp + 1_I4P ; stack(sp) = w(i)
      enddo
      if (sp > 0_I4P) call emit(stack(1:sp))
      endsubroutine split_walk

      subroutine emit(edges)
      !< Append one simple loop built from consecutive edges.
      integer(I4P), intent(in) :: edges(:)
      integer(I4P)             :: m

      if (size(edges) < 3) then
         n_open = n_open + size(edges, kind=I4P)  ! degenerate (two-edge) cycle: not a fillable hole
         return
      endif
      if (n_loops == size(buf)) then
         allocate(tmp(2 * size(buf)))
         do m = 1, n_loops
            call move_alloc(from=buf(m)%vid,   to=tmp(m)%vid)
            call move_alloc(from=buf(m)%x,     to=tmp(m)%x)
            call move_alloc(from=buf(m)%owner, to=tmp(m)%owner)
         enddo
         call move_alloc(from=tmp, to=buf)
      endif
      n_loops = n_loops + 1_I4P
      allocate(buf(n_loops)%vid(size(edges)), buf(n_loops)%x(size(edges)), buf(n_loops)%owner(size(edges)))
      do m = 1, size(edges)
         buf(n_loops)%vid(m)   = e_start(edges(m))
         buf(n_loops)%x(m)     = e_xs(edges(m))
         buf(n_loops)%owner(m) = e_owner(edges(m))
      enddo
      endsubroutine emit
   endsubroutine extract_boundary_loops
endmodule fossil_hole_filling
