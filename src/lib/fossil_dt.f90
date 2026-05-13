!< FOSSIL, 2D Delaunay triangulation via Bowyer-Watson (issue #18 §1.1 stage 2).

module fossil_dt
!< FOSSIL, 2D Delaunay triangulation via incremental Bowyer-Watson insertion.
!<
!< This is **Phase 1** of the eventual constrained Delaunay triangulator (CDT)
!< that §1.1 booleans needs. The full Sloan 1993 CDT layers constraint recovery
!< and Delaunay restoration on top of an unconstrained Delaunay foundation;
!< this module builds and tests that foundation in isolation. Phase 2
!< (constraint recovery via edge-walking + diagonal swapping) lands in a
!< follow-up against the data structures defined here.
!<
!< Phase-1 scope (this PR):
!<   - Bowyer-Watson incremental insertion of a 2D point cloud
!<   - Standard super-triangle bracketing + post-strip
!<   - Triangle / vertex / neighbor pointers maintained through every insertion
!<
!< Independently useful: an unconstrained Delaunay triangulation is also the
!< primitive that §1.7 isotropic remeshing's tangential-relaxation step needs
!< when reprojecting smoothed vertices onto the local 1-ring.
!<
!< Algorithm (Bowyer-Watson incremental):
!<
!<   1. Compute a bounding super-triangle large enough to contain every input
!<      point strictly inside it.
!<   2. Initialise the triangulation with the single super-triangle.
!<   3. For each input point P (in input order):
!<      a. Locate the triangle currently containing P by orient2d-walk.
!<      b. BFS outward from that triangle, collecting the **cavity** — the
!<         set of triangles whose circumcircle contains P (incircle test).
!<      c. Delete the cavity triangles; the cavity boundary is a star-shaped
!<         polygon visible from P.
!<      d. Re-triangulate the cavity by connecting P to every boundary edge.
!<      e. Update neighbor pointers across every new triangle.
!<   4. Strip super-triangle: deactivate every triangle incident to a
!<      super-vertex. The remaining active triangles index only into the
!<      original input point set.
!<
!< Edge convention (Sloan):
!<   - Edge `i` of triangle `t` is the edge **opposite** vertex `t%v(i)`.
!<     i.e. edge 1 connects v(2) and v(3); edge 2 connects v(3) and v(1);
!<     edge 3 connects v(1) and v(2).
!<   - `t%nbr(i)` is the triangle index sharing edge `i`, or 0 if boundary.
!<
!< Numeric robustness: predicates are called via `fossil_predicates`, which
!< currently provides `orient2d` and `incircle` as plain double-precision
!< stubs (Phase-1 acceptable per §1.1 stage 1's deferral). Upgrading those
!< predicates to adaptive Shewchuk later improves CDT robustness for free
!< through this module's existing call sites.

use fossil_predicates, only : orient2d, incircle
use penf,              only : I4P, R8P
use vecfor,            only : vector_R8P

implicit none
private
public :: triangulation_t
public :: dt_build

integer(I4P), parameter :: SUPER_OFFSET = 3_I4P  !< Number of super-triangle vertices appended at end of point list.

type :: triangle_t
   !< One triangle in the triangulation. Vertex indices reference the parent
   !< `triangulation_t%px / py` arrays; neighbor indices reference the parent's
   !< `tri` array, with 0 meaning "boundary edge, no neighbor".
   integer(I4P) :: v(3)   = 0_I4P  !< Vertex indices into px/py (CCW).
   integer(I4P) :: nbr(3) = 0_I4P  !< Neighbor triangle indices; 0 = boundary.
   logical      :: active = .false. !< False for triangles deleted during a cavity step.
endtype triangle_t

type :: triangulation_t
   !< Result of an incremental Delaunay build over an input point cloud.
   !<
   !< `px`, `py` hold every vertex coordinate including the three super-
   !< triangle vertices appended at the end (indices `n_input + 1 .. n_input + 3`).
   !< `tri` is the triangle array, with `n_tri_alloc` entries of which only
   !< those with `tri(i)%active == .true.` are valid; deleted triangles are
   !< left in place to keep indices stable across the build.
   integer(I4P)                       :: n_input = 0_I4P  !< Number of input points (excludes super-triangle).
   integer(I4P)                       :: n_tri_alloc = 0_I4P !< Allocated length of `tri`.
   integer(I4P)                       :: n_tri_used = 0_I4P !< Highest used index in `tri` (may exceed live triangle count).
   real(R8P), allocatable             :: px(:), py(:)     !< Point coordinates (input + super-triangle).
   type(triangle_t), allocatable      :: tri(:)           !< Triangle array; check %active before use.
   contains
      procedure, pass(self) :: num_triangles  !< Count of active triangles excluding super-triangle ones.
      procedure, pass(self) :: triangle_vertices !< Return v(3) for the i-th active triangle (1-indexed over active ones).
endtype triangulation_t

contains

   subroutine dt_build(tri, points)
   !< Build the Delaunay triangulation of an input 2D point cloud.
   !<
   !< `points(2, N)` holds (x, y) for each of N input points in any order; on
   !< return `tri` contains the Delaunay triangulation of those N points (the
   !< super-triangle has been stripped — every active triangle's vertices are
   !< in `[1, N]`).
   !<
   !< @note Input points must be **distinct** within EPS. Coincident or
   !<       near-coincident points produce degenerate triangles whose
   !<       circumcircle is ill-defined; current behaviour on duplicates is
   !<       to attempt insertion and let the FP predicate decide, which can
   !<       leave the triangulation in a non-Delaunay state. Caller is
   !<       responsible for deduplication; future work could add a snap-round
   !<       pass here.
   type(triangulation_t), intent(out) :: tri      !< Output triangulation.
   real(R8P),             intent(in)  :: points(:,:) !< Input points, shape (2, N).
   integer(I4P)                       :: n, i, locate, host, p_idx
   integer(I4P), allocatable          :: cavity(:)
   integer(I4P)                       :: n_cavity

   n = size(points, dim=2)
   tri%n_input = n

   ! Allocate point arrays with room for the 3 super-triangle vertices appended at the end.
   allocate(tri%px(n + SUPER_OFFSET), tri%py(n + SUPER_OFFSET))
   tri%px(1:n) = points(1, 1:n)
   tri%py(1:n) = points(2, 1:n)

   ! Worst-case triangle count for N points + 3 super-vertices is ~2*(N+3)-5
   ! Delaunay triangles, but during construction we hold both deleted and new
   ! ones; over-allocate ~6N to be safe.
   tri%n_tri_alloc = max(64_I4P, 6_I4P * (n + SUPER_OFFSET))
   allocate(tri%tri(tri%n_tri_alloc))
   tri%n_tri_used = 0_I4P

   call install_super_triangle(tri=tri)

   ! Cavity workspace; over-allocated to hold any plausible cavity size.
   allocate(cavity(tri%n_tri_alloc))

   ! Insert each input point, in order.
   do i = 1, n
      p_idx = i
      ! Locate the active triangle containing this point.
      host = locate_triangle(tri=tri, x=tri%px(p_idx), y=tri%py(p_idx))
      if (host == 0_I4P) cycle  ! defensive: should not happen if super-triangle covers input
      ! Collect the cavity (triangles violating Delaunay w.r.t. P) by BFS.
      call collect_cavity(tri=tri, p_idx=p_idx, host=host, cavity=cavity, n_cavity=n_cavity)
      ! Re-triangulate the cavity by fanning P to every boundary edge.
      call retriangulate_cavity(tri=tri, p_idx=p_idx, cavity=cavity, n_cavity=n_cavity)
   enddo

   call strip_super_triangle(tri=tri)
   endsubroutine dt_build

   subroutine install_super_triangle(tri)
   !< Append three super-triangle vertices to the point array and seed the
   !< triangulation with the single super-triangle.
   !<
   !< Sizing: the super-triangle must strictly contain every input point.
   !< Standard trick: take the bounding box of the input, blow it up by 10x
   !< the diagonal, and emit an equilateral-ish triangle around it. The
   !< constants here are conservative enough to handle any IEEE-double
   !< magnitude input the predicates can resolve.
   type(triangulation_t), intent(inout) :: tri
   real(R8P)                            :: xmin, xmax, ymin, ymax, dx, dy, dmax, midx, midy
   integer(I4P)                         :: n, sa, sb, sc

   n = tri%n_input
   if (n == 0_I4P) then
      xmin = -1._R8P ; xmax = 1._R8P ; ymin = -1._R8P ; ymax = 1._R8P
   else
      xmin = minval(tri%px(1:n)) ; xmax = maxval(tri%px(1:n))
      ymin = minval(tri%py(1:n)) ; ymax = maxval(tri%py(1:n))
   endif
   dx = xmax - xmin ; dy = ymax - ymin
   dmax = max(dx, dy, 1._R8P) * 100._R8P  ! generous margin
   midx = 0.5_R8P * (xmin + xmax)
   midy = 0.5_R8P * (ymin + ymax)

   sa = n + 1 ; sb = n + 2 ; sc = n + 3
   tri%px(sa) = midx - 2._R8P * dmax ; tri%py(sa) = midy - dmax
   tri%px(sb) = midx + 2._R8P * dmax ; tri%py(sb) = midy - dmax
   tri%px(sc) = midx                 ; tri%py(sc) = midy + 2._R8P * dmax

   tri%n_tri_used = 1_I4P
   tri%tri(1)%v   = [sa, sb, sc]
   tri%tri(1)%nbr = [0_I4P, 0_I4P, 0_I4P]
   tri%tri(1)%active = .true.
   endsubroutine install_super_triangle

   function locate_triangle(tri, x, y) result(host)
   !< Linear scan over active triangles to find one containing (x, y).
   !<
   !< Bowyer-Watson canonically uses an oriented walk (start from the
   !< previous insertion, walk through neighbour faces toward the target),
   !< which is asymptotically O(sqrt(N)) per insertion. The linear scan is
   !< O(N) per insertion, giving O(N^2) total — fine for the small cavities
   !< (≤ a few dozen points) that the §1.1 booleans triangulator actually
   !< feeds in, and far simpler to verify correct.
   !<
   !< @note Replace with the oriented walk if profiling on real CDT inputs
   !<       shows the locate step dominating.
   type(triangulation_t), intent(in) :: tri
   real(R8P),             intent(in) :: x, y
   integer(I4P)                      :: host
   integer(I4P)                      :: t

   host = 0_I4P
   do t = 1, tri%n_tri_used
      if (.not. tri%tri(t)%active) cycle
      if (point_in_triangle(tri=tri, t=t, x=x, y=y)) then
         host = t ; return
      endif
   enddo
   endfunction locate_triangle

   pure function point_in_triangle(tri, t, x, y) result(yes)
   !< True if (x, y) lies inside triangle `t` (or on its boundary).
   !< Uses three orient2d sign tests; consistent CCW orientation of triangles
   !< (maintained by Bowyer-Watson insertion) means all three signs must be
   !< non-negative for an interior point.
   type(triangulation_t), intent(in) :: tri
   integer(I4P),          intent(in) :: t
   real(R8P),             intent(in) :: x, y
   logical                           :: yes
   type(vector_R8P)                  :: a, b, c, p
   real(R8P)                         :: s1, s2, s3

   a = vector_R8P(tri%px(tri%tri(t)%v(1)), tri%py(tri%tri(t)%v(1)), 0._R8P)
   b = vector_R8P(tri%px(tri%tri(t)%v(2)), tri%py(tri%tri(t)%v(2)), 0._R8P)
   c = vector_R8P(tri%px(tri%tri(t)%v(3)), tri%py(tri%tri(t)%v(3)), 0._R8P)
   p = vector_R8P(x, y, 0._R8P)
   s1 = orient2d(a=a, b=b, c=p)
   s2 = orient2d(a=b, b=c, c=p)
   s3 = orient2d(a=c, b=a, c=p)
   yes = (s1 >= 0._R8P .and. s2 >= 0._R8P .and. s3 >= 0._R8P) .or. &
         (s1 <= 0._R8P .and. s2 <= 0._R8P .and. s3 <= 0._R8P)
   endfunction point_in_triangle

   subroutine collect_cavity(tri, p_idx, host, cavity, n_cavity)
   !< BFS outward from `host`, collecting every triangle whose circumcircle
   !< contains the point at index `p_idx`. The Delaunay invariant is that no
   !< point lies strictly inside any triangle's circumcircle; a triangle
   !< failing this test must be deleted and replaced.
   !<
   !< Marks collected triangles as inactive **before** adding their neighbors
   !< to the BFS queue, so each triangle is visited exactly once.
   type(triangulation_t), intent(inout) :: tri
   integer(I4P),          intent(in)    :: p_idx, host
   integer(I4P),          intent(out)   :: cavity(:)
   integer(I4P),          intent(out)   :: n_cavity
   integer(I4P)                         :: head, tail, t, e, nbr

   ! Use cavity(:) as both the BFS queue AND the output list — one pop at
   ! `head` advances the queue; pushes append at `tail`. After BFS the
   ! cavity contents are exactly cavity(1:n_cavity).
   n_cavity      = 1_I4P
   cavity(1)     = host
   tri%tri(host)%active = .false.

   head = 1_I4P ; tail = 1_I4P
   do while (head <= tail)
      t = cavity(head) ; head = head + 1
      do e = 1, 3
         nbr = tri%tri(t)%nbr(e)
         if (nbr == 0_I4P) cycle  ! boundary edge
         if (.not. tri%tri(nbr)%active) cycle  ! already in cavity
         if (in_circumcircle(tri=tri, t=nbr, p_idx=p_idx)) then
            tri%tri(nbr)%active = .false.
            tail = tail + 1
            cavity(tail) = nbr
         endif
      enddo
   enddo
   n_cavity = tail
   endsubroutine collect_cavity

   pure function in_circumcircle(tri, t, p_idx) result(yes)
   !< True if the point at index `p_idx` lies strictly inside triangle `t`'s
   !< circumcircle. The `incircle` predicate returns positive for inside when
   !< (a, b, c) is in CCW order (which Bowyer-Watson maintains).
   type(triangulation_t), intent(in) :: tri
   integer(I4P),          intent(in) :: t, p_idx
   logical                           :: yes
   type(vector_R8P)                  :: a, b, c, p

   a = vector_R8P(tri%px(tri%tri(t)%v(1)), tri%py(tri%tri(t)%v(1)), 0._R8P)
   b = vector_R8P(tri%px(tri%tri(t)%v(2)), tri%py(tri%tri(t)%v(2)), 0._R8P)
   c = vector_R8P(tri%px(tri%tri(t)%v(3)), tri%py(tri%tri(t)%v(3)), 0._R8P)
   p = vector_R8P(tri%px(p_idx), tri%py(p_idx), 0._R8P)
   yes = incircle(a=a, b=b, c=c, d=p) > 0._R8P
   endfunction in_circumcircle

   subroutine retriangulate_cavity(tri, p_idx, cavity, n_cavity)
   !< Re-triangulate the cavity by emitting one new triangle (P, u, v) for
   !< each cavity-boundary edge (u, v).
   !<
   !< Cavity-boundary edges are those edges of cavity triangles whose
   !< neighbor is **outside** the cavity (i.e. the neighbor is either active
   !< or 0 = boundary). For each such edge (u, v):
   !<   - emit triangle (P, u, v) — note vertex order matters for CCW
   !<     consistency: the cavity triangle had this edge as (u, v) opposite
   !<     a vertex w that's inside the cavity, so going u → v keeps the
   !<     traversal CCW around P.
   !<   - the new triangle's neighbor on its u-v edge is the cavity-external
   !<     triangle (the original neighbor); its other two edges' neighbors
   !<     are the two adjacent newly-emitted triangles, which we patch up
   !<     in a second pass.
   type(triangulation_t), intent(inout) :: tri
   integer(I4P),          intent(in)    :: p_idx
   integer(I4P),          intent(in)    :: cavity(:)
   integer(I4P),          intent(in)    :: n_cavity
   integer(I4P)                         :: c, t, e, nbr, e_op, u, v
   integer(I4P)                         :: first_new, last_new, new_t
   ! For neighbor stitching: per new triangle, its (u, v) so we can match
   ! adjacent boundary edges. Indexed [1..n_new]; n_new = number of
   ! boundary edges, equals length of the cavity's outer polygon.
   integer(I4P), allocatable            :: new_u(:), new_v(:), new_idx(:)
   integer(I4P)                         :: n_new, i, j

   ! Worst-case n_new is 3 * n_cavity (every cavity triangle contributes
   ! all three edges as boundary — happens when n_cavity == 1).
   allocate(new_u(3 * n_cavity), new_v(3 * n_cavity), new_idx(3 * n_cavity))
   n_new = 0_I4P
   first_new = tri%n_tri_used + 1
   last_new  = first_new - 1

   do c = 1, n_cavity
      t = cavity(c)
      do e = 1, 3
         nbr = tri%tri(t)%nbr(e)
         if (nbr /= 0_I4P) then
            if (.not. tri%tri(nbr)%active) cycle  ! both inside cavity, not a boundary edge
         endif
         ! Edge `e` of triangle `t`: opposite vertex v(e), so the edge is
         ! (v(e_next), v(e_next2)). With CCW convention (1,2,3,1,2,...):
         !   edge 1 -> v(2), v(3)
         !   edge 2 -> v(3), v(1)
         !   edge 3 -> v(1), v(2)
         u = tri%tri(t)%v(mod(e, 3) + 1)
         v = tri%tri(t)%v(mod(e + 1, 3) + 1)
         ! Emit new triangle (P, u, v); its edge opposite P is (u, v),
         ! adjacent to the external neighbor `nbr`.
         tri%n_tri_used = tri%n_tri_used + 1
         new_t = tri%n_tri_used
         tri%tri(new_t)%v   = [p_idx, u, v]
         tri%tri(new_t)%nbr = [nbr, 0_I4P, 0_I4P]  ! nbr opposite P; sides patched below
         tri%tri(new_t)%active = .true.
         last_new = new_t
         ! Patch the external neighbor's back-pointer to this new triangle.
         if (nbr /= 0_I4P) call patch_neighbor(tri=tri, side=nbr, old=t, new=new_t)
         ! Record (u, v) for the side-stitching pass below.
         n_new = n_new + 1
         new_u(n_new)   = u
         new_v(n_new)   = v
         new_idx(n_new) = new_t
      enddo
   enddo

   ! Side-stitching: for each new triangle T_i with vertices (P, u_i, v_i),
   ! identify the other two new triangles that share its non-(u_i, v_i) edges.
   !
   ! Edge convention recap (Sloan): edge `k` of a triangle is opposite vertex k.
   !   T_i has v(1)=P, v(2)=u_i, v(3)=v_i.
   !   edge 1 = opposite P  = (u_i, v_i)         -> external neighbor (set above)
   !   edge 2 = opposite u_i = (v_i, P)          -> shared with the next CCW T_j
   !   edge 3 = opposite v_i = (P, u_i)          -> shared with the prev CCW T_k
   !
   ! Going CCW around P, the next triangle T_j has u_j == v_i (its "P→u" edge
   ! connects P back to v_i). So T_i's edge 2 ↔ T_j's edge 3, with
   ! match condition `new_u(j) == new_v(i)`.
   ! The previous triangle T_k symmetrically has v_k == u_i, matching with
   ! condition `new_v(k) == new_u(i)` for T_i's edge 3 ↔ T_k's edge 2.
   do i = 1, n_new
      ! Edge 2 of T_i (the (v_i, P) edge): find T_j with u_j == v_i.
      do j = 1, n_new
         if (j == i) cycle
         if (new_u(j) == new_v(i)) then
            tri%tri(new_idx(i))%nbr(2) = new_idx(j)
            exit
         endif
      enddo
      ! Edge 3 of T_i (the (P, u_i) edge): find T_k with v_k == u_i.
      do j = 1, n_new
         if (j == i) cycle
         if (new_v(j) == new_u(i)) then
            tri%tri(new_idx(i))%nbr(3) = new_idx(j)
            exit
         endif
      enddo
   enddo
   endsubroutine retriangulate_cavity

   pure subroutine patch_neighbor(tri, side, old, new)
   !< When a cavity triangle `old` is replaced by `new` triangles, the
   !< external triangle `side` (which used to point to `old` via one of its
   !< three nbr slots) must be redirected to point to `new`.
   type(triangulation_t), intent(inout) :: tri
   integer(I4P),          intent(in)    :: side, old, new
   integer(I4P)                         :: e

   do e = 1, 3
      if (tri%tri(side)%nbr(e) == old) then
         tri%tri(side)%nbr(e) = new
         return
      endif
   enddo
   endsubroutine patch_neighbor

   subroutine strip_super_triangle(tri)
   !< Deactivate every triangle that has at least one super-triangle vertex.
   !< Also clear neighbor pointers from surviving triangles back to the
   !< stripped ones (so external code iterating neighbours doesn't follow
   !< a dangling reference).
   type(triangulation_t), intent(inout) :: tri
   integer(I4P)                         :: t, e, n, sa, sb, sc

   n = tri%n_input
   sa = n + 1 ; sb = n + 2 ; sc = n + 3

   do t = 1, tri%n_tri_used
      if (.not. tri%tri(t)%active) cycle
      if (any(tri%tri(t)%v == sa) .or. any(tri%tri(t)%v == sb) .or. any(tri%tri(t)%v == sc)) then
         tri%tri(t)%active = .false.
      endif
   enddo

   ! Clean up dangling neighbor pointers in surviving triangles.
   do t = 1, tri%n_tri_used
      if (.not. tri%tri(t)%active) cycle
      do e = 1, 3
         if (tri%tri(t)%nbr(e) /= 0_I4P) then
            if (.not. tri%tri(tri%tri(t)%nbr(e))%active) tri%tri(t)%nbr(e) = 0_I4P
         endif
      enddo
   enddo
   endsubroutine strip_super_triangle

   pure function num_triangles(self) result(n)
   !< Count of active triangles (post-strip, so super-triangle ones excluded).
   class(triangulation_t), intent(in) :: self
   integer(I4P)                       :: n
   integer(I4P)                       :: t

   n = 0_I4P
   do t = 1, self%n_tri_used
      if (self%tri(t)%active) n = n + 1
   enddo
   endfunction num_triangles

   pure subroutine triangle_vertices(self, k, v)
   !< Return the vertex indices of the `k`-th active triangle (1-indexed
   !< over active triangles only). Vertex indices are into `self%px / py`
   !< and lie in `[1, self%n_input]` after strip_super_triangle.
   class(triangulation_t), intent(in)  :: self
   integer(I4P),           intent(in)  :: k       !< 1-based index over active triangles.
   integer(I4P),           intent(out) :: v(3)    !< Vertex indices.
   integer(I4P)                        :: t, count

   v = 0_I4P
   count = 0_I4P
   do t = 1, self%n_tri_used
      if (.not. self%tri(t)%active) cycle
      count = count + 1
      if (count == k) then
         v = self%tri(t)%v
         return
      endif
   enddo
   endsubroutine triangle_vertices

endmodule fossil_dt
