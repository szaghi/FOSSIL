!< FOSSIL, hole filling: boundary loops and their triangulation (issue #31).
module fossil_hole_filling
!< Close small holes of a triangulated surface.
!<
!< A hole is a closed loop of facet edges without a unique neighbour (`fcon_edge == 0` after
!< `analyze`). Such an edge is a boundary edge (one incident facet) or lies on a non-manifold edge
!< (3+ incident facets, none linked); loops are chained over both kinds, as before issue #88, but
!< the counts tell them apart. Loops are built on pool vertex ids (welded vertices), so coincident
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
use fossil_csr_matrix,   only : csr_matrix_t
use fossil_linear_solver, only : cg_solve, CG_STATUS_OK
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: boundary_loop_t
public :: extract_boundary_loops
public :: triangulate_loop
public :: refine_patch
public :: fair_patch

type :: boundary_loop_t
   !< One closed boundary loop, in patch orientation.
   integer(I4P),     allocatable :: vid(:)   !< Pool vertex ids; edge k joins vid(k) -> vid(k+1) (cyclic).
   type(vector_R8P), allocatable :: x(:)     !< Coordinates of vid(:).
   integer(I4P),     allocatable :: owner(:) !< Mesh facet owning edge k.
endtype boundary_loop_t

contains
   subroutine extract_boundary_loops(facet, loops, n_open, n_nonmanifold)
   !< Extract every closed boundary loop of the facet set.
   !<
   !< Requires pool vertex ids (`vertex_id`, set by `analyze`) and connectivity (`fcon_edge`).
   !< A facet edge without a unique neighbour is a boundary edge if no other facet has it, else it
   !< lies on a non-manifold edge (issue #88). Both kinds are chained into loops. Boundary edges
   !< left out of every closed loop (open chains, inconsistent orientation, two-edge cycles) are
   !< counted in `n_open`; the non-manifold edges met are counted in `n_nonmanifold`, once per edge
   !< (before #88 their facet edges left unchained were counted in `n_open` as boundary edges).
   type(facet_object),                 intent(in)            :: facet(:)      !< Facets (analyzed).
   type(boundary_loop_t), allocatable, intent(out)           :: loops(:)      !< Closed simple loops.
   integer(I4P),                       intent(out)           :: n_open        !< Boundary edges left out of loops.
   integer(I4P),                       intent(out), optional :: n_nonmanifold !< Non-manifold edges met (3+ facets).
   integer(I4P),     allocatable :: e_start(:), e_end(:), e_owner(:) !< Directed boundary edges (patch order).
   type(vector_R8P), allocatable :: e_xs(:)                        !< Start coordinates of each edge.
   integer(I4P),     allocatable :: head(:), next(:)               !< Per-start-vertex linked lists of edges.
   logical,          allocatable :: used(:)                        !< Edge already walked.
   integer(I4P),     allocatable :: walk(:)                        !< Edge indices of the current walk.
   type(boundary_loop_t), allocatable :: buf(:), tmp(:)            !< Loop accumulator.
   integer(I4P) :: n_edges, f, e, k, v_max, e0, cur, nxt, n_walk, n_loops
   integer(I4P) :: q, mult, n_nm                                   !< Same-edge search.
   logical      :: first                                           !< Lowest slot of its edge.
   logical,      allocatable :: nonman(:)                          !< Edge k lies on a non-manifold edge.

   n_open = 0_I4P ; n_nm = 0_I4P
   if (present(n_nonmanifold)) n_nonmanifold = 0_I4P
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
   ! 2b. A facet edge whose edge another unlinked slot shares (either direction) lies on a
   ! non-manifold edge (two linked facets never leave a 0).
   allocate(nonman(n_edges)) ; nonman = .false.
   do k = 1, n_edges
      mult = 0_I4P ; first = .true.
      q = head(e_start(k))
      do while (q > 0_I4P)
         if (q /= k .and. e_end(q) == e_end(k)) then
            mult = mult + 1_I4P ; if (q < k) first = .false.
         endif
         q = next(q)
      enddo
      q = head(e_end(k))
      do while (q > 0_I4P)
         if (q /= k .and. e_end(q) == e_start(k)) then
            mult = mult + 1_I4P ; if (q < k) first = .false.
         endif
         q = next(q)
      enddo
      if (mult > 0_I4P) then
         nonman(k) = .true.
         if (first) n_nm = n_nm + 1_I4P
      endif
   enddo
   if (present(n_nonmanifold)) n_nonmanifold = n_nm
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
         n_open = n_open + count(.not. nonman(walk(1:n_walk)), kind=I4P)
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
         n_open = n_open + count(.not. nonman(edges), kind=I4P)  ! degenerate (two-edge) cycle: not a fillable hole
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

   subroutine refine_patch(loop, patch)
   !< Refine a hole patch so its triangles match the density of the surrounding mesh (Liepa 2003,
   !< issue #31 step 5). The loop vertices never move and no triangle is added outside the patch.
   !<
   !<  - Scale: each loop vertex gets the mean length of its two loop edges; a new vertex gets the
   !<    mean scale of the triangle it splits.
   !<  - Split: a triangle with centroid c and centroid scale s_c is split into three at c when
   !<    sqrt(2) |c - v| > max(s_c, s_v) for each of its vertices v.
   !<  - Relax: an interior edge (any patch edge but a loop edge) shared by triangles (a, b, c) and
   !<    (b, a, d) is flipped to (c, d) when the angles at c and d sum to more than pi (the edge is not
   !<    locally Delaunay), provided both new triangles keep the patch's orientation.
   !< Split and relax alternate until no triangle splits (at most 50 rounds, 20000 vertices).
   !< On return `patch` holds the refined triangles; loop vertices keep their pool ids, new vertices
   !< have id 0 and identical coordinates in every triangle that shares them.
   type(boundary_loop_t),           intent(in)    :: loop     !< Loop the patch fills.
   type(facet_object), allocatable, intent(inout) :: patch(:) !< Patch triangles (from triangulate_loop).
   type(vector_R8P), allocatable :: x(:)        !< Patch vertices (1:nb are the loop's).
   real(R8P),        allocatable :: sc(:)       !< Vertex scales.
   integer(I4P),     allocatable :: t(:,:)      !< Triangles (3, nt), patch orientation.
   integer(I4P)                  :: nb, nv, nt, i, k, round, j, m
   logical                       :: split
   real(R8P)                     :: c(3), s_c
   integer(I4P), parameter       :: MAX_ROUNDS = 50_I4P, MAX_VERTICES = 20000_I4P

   nb = size(loop%vid, kind=I4P)
   if (nb < 3_I4P .or. size(patch) == 0) return
   allocate(x(nb), sc(nb), t(3, size(patch)))
   x = loop%x
   do i = 1, nb
      sc(i) = 0.5_R8P * (dist(x(i), x(prev(i))) + dist(x(i), x(next(i))))
   enddo
   nv = nb ; nt = size(patch, kind=I4P)
   do k = 1, nt
      do j = 1, 3
         t(j, k) = loop_index(patch(k)%vertex_id(j))
         if (t(j, k) == 0_I4P) return ! not a loop patch: leave it untouched
      enddo
   enddo

   do round = 1, MAX_ROUNDS
      split = .false.
      m = nt
      do k = 1, m
         c = ([x(t(1, k))%x, x(t(1, k))%y, x(t(1, k))%z] + [x(t(2, k))%x, x(t(2, k))%y, x(t(2, k))%z] + &
              [x(t(3, k))%x, x(t(3, k))%y, x(t(3, k))%z]) / 3._R8P
         s_c = (sc(t(1, k)) + sc(t(2, k)) + sc(t(3, k))) / 3._R8P
         if (all([(sqrt(2._R8P) * dist3(c, x(t(j, k))) > max(s_c, sc(t(j, k))), j = 1, 3)])) then
            if (nv >= MAX_VERTICES) exit
            call push_vertex(c, s_c)
            call push_triangle(t(2, k), t(3, k), nv)
            call push_triangle(t(3, k), t(1, k), nv)
            t(3, k) = nv
            split = .true.
         endif
      enddo
      call relax
      if (.not. split) exit
   enddo

   deallocate(patch) ; allocate(patch(nt))
   do k = 1, nt
      do j = 1, 3
         patch(k)%vertex(j) = x(t(j, k))
         patch(k)%vertex_id(j) = 0_I4P
         if (t(j, k) <= nb) patch(k)%vertex_id(j) = loop%vid(t(j, k))
      enddo
      call patch(k)%compute_metrix
   enddo
   contains
      pure function prev(i) result(j)
      integer(I4P), intent(in) :: i
      integer(I4P)             :: j

      j = i - 1_I4P ; if (j < 1_I4P) j = nb
      endfunction prev

      pure function next(i) result(j)
      integer(I4P), intent(in) :: i
      integer(I4P)             :: j

      j = i + 1_I4P ; if (j > nb) j = 1_I4P
      endfunction next

      pure function loop_index(vid) result(i)
      integer(I4P), intent(in) :: vid
      integer(I4P)             :: i

      do i = 1, nb
         if (loop%vid(i) == vid) return
      enddo
      i = 0_I4P
      endfunction loop_index

      pure function is_loop_edge(a, b) result(yes)
      integer(I4P), intent(in) :: a, b
      logical                  :: yes

      yes = a <= nb .and. b <= nb .and. (next(a) == b .or. next(b) == a)
      endfunction is_loop_edge

      subroutine push_vertex(p, s)
      real(R8P), intent(in)         :: p(3), s
      type(vector_R8P), allocatable :: xt(:)
      real(R8P),        allocatable :: st(:)

      if (nv == size(x)) then
         allocate(xt(2 * nv), st(2 * nv))
         xt(1:nv) = x(1:nv) ; st(1:nv) = sc(1:nv)
         call move_alloc(from=xt, to=x) ; call move_alloc(from=st, to=sc)
      endif
      nv = nv + 1_I4P
      x(nv)%x = p(1) ; x(nv)%y = p(2) ; x(nv)%z = p(3)
      sc(nv) = s
      endsubroutine push_vertex

      subroutine push_triangle(a, b, cc)
      integer(I4P), intent(in)  :: a, b, cc
      integer(I4P), allocatable :: tt(:,:)

      if (nt == size(t, dim=2)) then
         allocate(tt(3, 2 * nt))
         tt(:, 1:nt) = t(:, 1:nt)
         call move_alloc(from=tt, to=t)
      endif
      nt = nt + 1_I4P
      t(:, nt) = [a, b, cc]
      endsubroutine push_triangle

      subroutine relax
      !< Flip non-locally-Delaunay interior edges, sweep after sweep, until a sweep flips none (at
      !< most 100 sweeps). After a flip the scan moves on to the next triangle.
      integer(I4P) :: sweep, k1, k2, e1, e2, a, b, cc, dd
      logical      :: flipped, done_k1
      real(R8P)    :: n_old(3), n1(3), n2(3), ar, a1, a2

      do sweep = 1, 100
         flipped = .false.
         do k1 = 1, nt
            done_k1 = .false.
            do e1 = 1, 3
               a = t(e1, k1) ; b = t(mod(e1, 3) + 1, k1) ; cc = t(mod(e1 + 1, 3) + 1, k1)
               if (is_loop_edge(a, b)) cycle
               ! the triangle across (a, b) runs it as (b, a)
               do k2 = 1, nt
                  if (k2 == k1) cycle
                  do e2 = 1, 3
                     if (t(e2, k2) == b .and. t(mod(e2, 3) + 1, k2) == a) then
                        dd = t(mod(e2 + 1, 3) + 1, k2)
                        if (angle_at(cc, a, b) + angle_at(dd, a, b) > acos(-1._R8P) + 1.e-12_R8P .and. &
                            .not. has_edge(cc, dd)) then
                           call tri_normal(x(a), x(b), x(cc), n_old, ar)
                           call tri_normal(x(a), x(dd), x(cc), n1, a1)
                           call tri_normal(x(dd), x(b), x(cc), n2, a2)
                           if (a1 > 0._R8P .and. a2 > 0._R8P .and. sum(n1 * n_old) > 0._R8P .and. &
                               sum(n2 * n_old) > 0._R8P) then
                              t(:, k1) = [a, dd, cc]
                              t(:, k2) = [dd, b, cc]
                              flipped = .true. ; done_k1 = .true.
                           endif
                        endif
                        exit
                     endif
                  enddo
                  if (done_k1 .or. e2 <= 3) exit ! the neighbour across (a, b) was found
               enddo
               if (done_k1) exit
            enddo
         enddo
         if (.not. flipped) exit
      enddo
      endsubroutine relax

      pure function has_edge(a, b) result(yes)
      integer(I4P), intent(in) :: a, b
      logical                  :: yes
      integer(I4P)             :: kk, ee

      yes = .false.
      do kk = 1, nt
         do ee = 1, 3
            if ((t(ee, kk) == a .and. t(mod(ee, 3) + 1, kk) == b) .or. (t(ee, kk) == b .and. t(mod(ee, 3) + 1, kk) == a)) then
               yes = .true. ; return
            endif
         enddo
      enddo
      endfunction has_edge

      pure function angle_at(apex, a, b) result(ang)
      !< Angle at `apex` of the triangle (apex, a, b).
      integer(I4P), intent(in) :: apex, a, b
      real(R8P)                :: ang
      real(R8P)                :: u(3), v(3), nu, nv_

      u = [x(a)%x - x(apex)%x, x(a)%y - x(apex)%y, x(a)%z - x(apex)%z]
      v = [x(b)%x - x(apex)%x, x(b)%y - x(apex)%y, x(b)%z - x(apex)%z]
      nu = sqrt(sum(u * u)) ; nv_ = sqrt(sum(v * v))
      if (nu * nv_ <= 0._R8P) then
         ang = 0._R8P
      else
         ang = acos(max(-1._R8P, min(1._R8P, sum(u * v) / (nu * nv_))))
      endif
      endfunction angle_at

      pure function dist(p1, p2) result(dd)
      type(vector_R8P), intent(in) :: p1, p2
      real(R8P)                    :: dd

      dd = sqrt((p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2)
      endfunction dist

      pure function dist3(p1, p2) result(dd)
      real(R8P),        intent(in) :: p1(3)
      type(vector_R8P), intent(in) :: p2
      real(R8P)                    :: dd

      dd = sqrt((p1(1) - p2%x)**2 + (p1(2) - p2%y)**2 + (p1(3) - p2%z)**2)
      endfunction dist3
   endsubroutine refine_patch

   subroutine fair_patch(loop, patch, facet, status)
   !< Fair a refined hole patch (issue #31 step 5): move its interior vertices to the thin-plate
   !< (bi-Laplacian) surface spanned by the loop, so the patch continues the surrounding surface
   !< with matching tangents instead of lying flat.
   !<
   !< Unknowns: the patch vertices that are not loop vertices (refine_patch adds them; a plain
   !< patch has none and is returned unchanged). Fixed: the loop vertices and their 1-ring in the
   !< surface (the vertices of the mesh facets around each loop vertex), which carries the
   !< tangent information. Feature-aware ring: at each loop vertex only the mesh facets whose
   !< normal is within 60 degrees of the patch's own normal there (area-weighted normal of the
   !< patch triangles at that vertex, before fairing) contribute. A hole along a sharp edge (a CAD
   !< crease, the edge of a box) then keeps the crease instead of being bent toward the face
   !< across it; the owners of the loop edges are no reference, since along a crease they are the
   !< facets across it. With the uniform (graph) Laplacian (L x)_v = sum_{u ~ v} (x_u - x_v),
   !< the equations (L (L x))_i = 0 for the interior vertices i form a symmetric positive-definite
   !< system (L^2 restricted to the interior: each interior vertex reaches the fixed boundary), solved
   !< by `cg_solve` once per coordinate. Uniform weights suit a refined patch, whose edges are close
   !< to uniform by construction, and need no cotangent of a possibly degenerate initial triangle.
   !< Loop vertices are copied, so the patch still welds to the mesh bit for bit.
   type(boundary_loop_t),           intent(in)    :: loop     !< Loop the patch fills.
   type(facet_object), allocatable, intent(inout) :: patch(:) !< Refined patch triangles.
   type(facet_object),              intent(in)    :: facet(:) !< Mesh facets (with pool vertex ids).
   integer(I4P),                    intent(out)   :: status   !< 0 = ok (or nothing to do); CG status otherwise.
   type(vector_R8P), allocatable :: x(:)          !< Graph nodes: 1:nb loop, then interior, then ring.
   integer(I4P),     allocatable :: node_vid(:)   !< Pool id of loop and ring nodes (0 for interior).
   integer(I4P),     allocatable :: tri(:,:)      !< Patch triangles as node indices.
   integer(I4P),     allocatable :: adj(:,:), deg(:) !< Adjacency lists (capacity grows) and degrees.
   integer(I4P),     allocatable :: unk(:)        !< Unknown index of a node (0 = fixed).
   type(csr_matrix_t)            :: a
   real(R8P),        allocatable :: b(:,:), sol(:)
   integer(I4P)                  :: nb, nn, ni, nr, k, j, f, i, kk, jj, c, st, n_unk
   real(R8P)                     :: lik, lkj
   real(R8P),        allocatable :: pn(:,:)       !< Patch normal at each loop vertex.

   status = 0_I4P
   nb = size(loop%vid, kind=I4P)
   if (nb < 3_I4P .or. size(patch) == 0) return
   ! nodes: loop vertices, then the patch's interior vertices (id 0, matched by coordinates)
   allocate(x(nb), node_vid(nb), tri(3, size(patch)))
   x = loop%x ; node_vid = loop%vid ; nn = nb
   do k = 1, size(patch)
      do j = 1, 3
         if (patch(k)%vertex_id(j) > 0_I4P) then
            tri(j, k) = find_vid(patch(k)%vertex_id(j))
         else
            tri(j, k) = find_or_add_coord(patch(k)%vertex(j))
         endif
         if (tri(j, k) == 0_I4P) return ! a patch vertex outside the loop: leave the patch untouched
      enddo
   enddo
   ni = nn - nb
   if (ni == 0_I4P) return
   ! the patch's normal at each loop vertex (area-weighted, unnormalized cross products)
   allocate(pn(3, nb)) ; pn = 0._R8P
   do k = 1, size(patch)
      do j = 1, 3
         if (tri(j, k) <= nb) pn(:, tri(j, k)) = pn(:, tri(j, k)) + cross3(x(tri(1, k)), x(tri(2, k)), x(tri(3, k)))
      enddo
   enddo
   do i = 1, nb
      lik = sqrt(sum(pn(:, i)**2))
      if (lik > 0._R8P) pn(:, i) = pn(:, i) / lik
   enddo
   ! ring: the other vertices of the mesh facets around each loop vertex
   allocate(deg(nn), adj(16, nn)) ; deg = 0_I4P
   do k = 1, size(patch)
      do j = 1, 3
         call link(tri(j, k), tri(mod(j, 3) + 1, k))
      enddo
   enddo
   nr = 0_I4P
   do f = 1, size(facet)
      do j = 1, 3
         i = find_vid(facet(f)%vertex_id(j))
         if (i == 0_I4P .or. i > nb) cycle
         if (.not. smooth_with_loop(f, i)) cycle ! across a sharp edge: not part of the ring
         do jj = 1, 3
            if (jj == j) cycle
            c = find_vid(facet(f)%vertex_id(jj))
            if (c == 0_I4P) then
               c = add_ring(facet(f)%vertex(jj), facet(f)%vertex_id(jj))
               nr = nr + 1_I4P
            endif
            call link(i, c)
         enddo
      enddo
   enddo
   ! unknowns: interior nodes nb+1 .. nb+ni
   allocate(unk(nn)) ; unk = 0_I4P
   n_unk = 0_I4P
   do i = nb + 1, nb + ni
      n_unk = n_unk + 1_I4P ; unk(i) = n_unk
   enddo
   ! A = (L^2) restricted to the unknowns; b = - (L^2) applied to the fixed values
   call a%initialize(n_unk, n_unk)
   allocate(b(n_unk, 3), sol(n_unk)) ; b = 0._R8P
   do i = nb + 1, nb + ni
      do kk = 0, deg(i)
         k = i ; if (kk > 0) k = adj(kk, i)
         lik = 1._R8P ; if (k == i) lik = -real(deg(i), R8P)
         do jj = 0, deg(k)
            j = k ; if (jj > 0) j = adj(jj, k)
            lkj = 1._R8P ; if (j == k) lkj = -real(deg(k), R8P)
            if (unk(j) > 0_I4P) then
               call a%append(unk(i), unk(j), lik * lkj)
            else
               b(unk(i), :) = b(unk(i), :) - lik * lkj * [x(j)%x, x(j)%y, x(j)%z]
            endif
         enddo
      enddo
   enddo
   call a%finalize
   do c = 1, 3
      do i = nb + 1, nb + ni ! initial guess: current position
         sol(unk(i)) = merge(x(i)%x, merge(x(i)%y, x(i)%z, c == 2), c == 1)
      enddo
      call cg_solve(a, b(:, c), sol, tol=1.e-12_R8P, status=st)
      if (st /= CG_STATUS_OK) then
         status = st
         return
      endif
      do i = nb + 1, nb + ni
         select case (c)
         case (1) ; x(i)%x = sol(unk(i))
         case (2) ; x(i)%y = sol(unk(i))
         case (3) ; x(i)%z = sol(unk(i))
         endselect
      enddo
   enddo
   do k = 1, size(patch)
      do j = 1, 3
         if (tri(j, k) > nb) patch(k)%vertex(j) = x(tri(j, k))
      enddo
      call patch(k)%compute_metrix
   enddo
   contains
      function smooth_with_loop(ff, li) result(yes)
      !< Facet ff's normal is within 60 degrees of the patch's normal at loop vertex li.
      integer(I4P), intent(in) :: ff, li
      logical                  :: yes
      real(R8P), parameter     :: COS_FEATURE = 0.5_R8P ! cos(60 deg)

      associate(nf => facet(ff)%normal)
         yes = pn(1, li) * nf%x + pn(2, li) * nf%y + pn(3, li) * nf%z >= COS_FEATURE
      end associate
      endfunction smooth_with_loop

      pure function cross3(p1, p2, p3) result(w)
      !< (p2 - p1) x (p3 - p1).
      type(vector_R8P), intent(in) :: p1, p2, p3
      real(R8P)                    :: w(3)
      real(R8P)                    :: u(3), v(3)

      u = [p2%x - p1%x, p2%y - p1%y, p2%z - p1%z]
      v = [p3%x - p1%x, p3%y - p1%y, p3%z - p1%z]
      w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
      endfunction cross3

      function find_vid(vid) result(idx)
      integer(I4P), intent(in) :: vid
      integer(I4P)             :: idx

      idx = 0_I4P
      if (vid <= 0_I4P) return
      do idx = 1, size(node_vid)
         if (node_vid(idx) == vid) return
      enddo
      idx = 0_I4P
      endfunction find_vid

      function find_or_add_coord(p) result(idx)
      type(vector_R8P), intent(in) :: p
      integer(I4P)                 :: idx

      do idx = nb + 1, nn
         if (x(idx)%x == p%x .and. x(idx)%y == p%y .and. x(idx)%z == p%z) return
      enddo
      x = [x, p] ; node_vid = [node_vid, 0_I4P] ; nn = nn + 1_I4P
      idx = nn
      endfunction find_or_add_coord

      function add_ring(p, vid) result(idx)
      type(vector_R8P), intent(in) :: p
      integer(I4P),     intent(in) :: vid
      integer(I4P)                 :: idx
      integer(I4P), allocatable    :: tmp(:,:)

      x = [x, p] ; node_vid = [node_vid, vid] ; nn = nn + 1_I4P
      deg = [deg, 0_I4P]
      allocate(tmp(size(adj, dim=1), nn)) ; tmp = 0_I4P
      tmp(:, 1:nn - 1) = adj
      call move_alloc(from=tmp, to=adj)
      idx = nn
      endfunction add_ring

      subroutine link(p, q)
      !< Undirected edge p - q in the adjacency lists (no duplicates).
      integer(I4P), intent(in)  :: p, q
      integer(I4P), allocatable :: tmp(:,:)

      if (p == q) return
      if (any(adj(1:deg(p), p) == q)) return
      if (deg(p) == size(adj, dim=1) .or. deg(q) == size(adj, dim=1)) then
         allocate(tmp(2 * size(adj, dim=1), size(adj, dim=2))) ; tmp = 0_I4P
         tmp(1:size(adj, dim=1), :) = adj
         call move_alloc(from=tmp, to=adj)
      endif
      deg(p) = deg(p) + 1_I4P ; adj(deg(p), p) = q
      deg(q) = deg(q) + 1_I4P ; adj(deg(q), q) = p
      endsubroutine link
   endsubroutine fair_patch

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
