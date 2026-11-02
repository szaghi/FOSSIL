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
public :: triangulate_loop

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

   subroutine triangulate_loop(loop, facet, patch, status)
   !< Triangulate one boundary loop by Liepa's minimum-weight dynamic programming (Liepa 2003).
   !<
   !< Over sub-polygons (i, k) of the loop, the best apex m minimizes, lexicographically:
   !<  1. the largest dihedral angle (angle between unit normals) of triangle (i, m, k) with its
   !<     neighbours across edges (i, m) and (m, k): the sub-problems' best triangles, or the
   !<     owning mesh facet when the edge is a loop edge; the top triangle also checks the closing
   !<     loop edge (n -> 1) against its owner. Degenerate triangles count as angle pi;
   !<  2. the total area.
   !< O(n^3) time, O(n^2) memory. `patch` gets n - 2 facets with metrix computed, oriented like the
   !< mesh (loops are stored in patch orientation). `status` = 0, or 1 if the loop has < 3 vertices.
   type(boundary_loop_t),           intent(in)  :: loop      !< Loop to fill.
   type(facet_object),              intent(in)  :: facet(:)  !< Mesh facets (owners of the loop edges).
   type(facet_object), allocatable, intent(out) :: patch(:)  !< Patch triangles.
   integer(I4P),                    intent(out) :: status    !< 0 = ok.
   real(R8P),    parameter   :: PI_ = acos(-1._R8P)
   real(R8P),    parameter   :: TIE = 1.e-12_R8P           !< Angle tie tolerance before comparing areas.
   real(R8P),    allocatable :: w_ang(:,:), w_area(:,:)    !< DP weights per sub-polygon (i, k).
   integer(I4P), allocatable :: apex(:,:)                  !< Best apex per sub-polygon.
   integer(I4P)              :: n, i, k, m, gap, sp, np
   integer(I4P), allocatable :: stack(:,:)
   real(R8P)                 :: ang, area, nt(3), tarea, best_ang, best_area

   n = size(loop%vid, kind=I4P)
   status = 1_I4P
   allocate(patch(0))
   if (n < 3_I4P) return
   status = 0_I4P
   allocate(w_ang(n, n), w_area(n, n), apex(n, n))
   w_ang = 0._R8P ; w_area = 0._R8P ; apex = 0_I4P
   do gap = 2, n - 1
      do i = 1, n - gap
         k = i + gap
         best_ang = huge(1._R8P) ; best_area = huge(1._R8P)
         do m = i + 1, k - 1
            call tri_normal(loop%x(i), loop%x(m), loop%x(k), nt, tarea)
            if (tarea <= 0._R8P) then
               ang = PI_
            else
               ang = max(w_ang(i, m), w_ang(m, k))
               ang = max(ang, across(i, m, nt), across(m, k, nt))
               if (i == 1 .and. k == n) ang = max(ang, angle_to(facet(loop%owner(n))%normal, nt))
            endif
            area = w_area(i, m) + w_area(m, k) + tarea
            if (ang < best_ang - TIE .or. (abs(ang - best_ang) <= TIE .and. area < best_area)) then
               best_ang = ang ; best_area = area ; apex(i, k) = m
            endif
         enddo
         w_ang(i, k) = best_ang ; w_area(i, k) = best_area
      enddo
   enddo
   ! Recover the triangles from the apex table.
   deallocate(patch) ; allocate(patch(n - 2)) ; np = 0
   allocate(stack(2, n)) ; sp = 1 ; stack(:, 1) = [1_I4P, n]
   do while (sp > 0)
      i = stack(1, sp) ; k = stack(2, sp) ; sp = sp - 1
      if (k - i < 2) cycle
      m = apex(i, k)
      np = np + 1
      patch(np)%vertex(1) = loop%x(i) ; patch(np)%vertex(2) = loop%x(m) ; patch(np)%vertex(3) = loop%x(k)
      patch(np)%vertex_id = [loop%vid(i), loop%vid(m), loop%vid(k)]
      call patch(np)%compute_metrix
      sp = sp + 1 ; stack(:, sp) = [i, m]
      sp = sp + 1 ; stack(:, sp) = [m, k]
   enddo
   contains
      function across(a, b, nt) result(ang)
      !< Dihedral angle of a triangle with normal nt across its edge (a, b): against the mesh owner
      !< if (a, b) is a loop edge, else against the sub-problem's best triangle on (a, b).
      integer(I4P), intent(in) :: a, b
      real(R8P),    intent(in) :: nt(3)
      real(R8P)                :: ang, nb(3), ab
      if (b == a + 1) then
         ang = angle_to(facet(loop%owner(a))%normal, nt)
      else
         call tri_normal(loop%x(a), loop%x(apex(a, b)), loop%x(b), nb, ab)
         if (ab <= 0._R8P) then
            ang = PI_
         else
            ang = acos(max(-1._R8P, min(1._R8P, dot_product(nb, nt))))
         endif
      endif
      endfunction across

      function angle_to(v, nt) result(ang)
      type(vector_R8P), intent(in) :: v
      real(R8P),        intent(in) :: nt(3)
      real(R8P)                    :: ang
      ang = acos(max(-1._R8P, min(1._R8P, v%x * nt(1) + v%y * nt(2) + v%z * nt(3))))
      endfunction angle_to
   endsubroutine triangulate_loop

   pure subroutine tri_normal(a, b, c, nrm, area)
   !< Unit normal and area of triangle (a, b, c); area 0 and normal 0 if degenerate.
   type(vector_R8P), intent(in)  :: a, b, c
   real(R8P),        intent(out) :: nrm(3), area
   real(R8P)                     :: u(3), v(3), w(3), len2

   u = [b%x - a%x, b%y - a%y, b%z - a%z]
   v = [c%x - a%x, c%y - a%y, c%z - a%z]
   w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
   len2 = sum(w * w)
   if (len2 <= tiny(1._R8P)) then
      nrm = 0._R8P ; area = 0._R8P
   else
      nrm = w / sqrt(len2) ; area = 0.5_R8P * sqrt(len2)
   endif
   endsubroutine tri_normal
endmodule fossil_hole_filling
