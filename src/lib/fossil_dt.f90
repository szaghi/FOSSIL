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
public :: dt_build, cdt_build
public :: DT_STATUS_OK, DT_STATUS_CROSSING_CONSTRAINTS, DT_STATUS_BAD_SEGMENT, DT_STATUS_RECOVERY_FAILED

integer(I4P), parameter :: SUPER_OFFSET = 3_I4P  !< Number of super-triangle vertices appended at end of point list.

! Status codes for cdt_build.
integer(I4P), parameter :: DT_STATUS_OK                   = 0_I4P  !< Successful build.
integer(I4P), parameter :: DT_STATUS_CROSSING_CONSTRAINTS = 1_I4P  !< Two constraint segments cross each other.
integer(I4P), parameter :: DT_STATUS_BAD_SEGMENT          = 2_I4P  !< Constraint segment endpoint out of range or coincident.
integer(I4P), parameter :: DT_STATUS_RECOVERY_FAILED      = 3_I4P  !< Recovery loop did not converge (bug or pathological input).

type :: triangle_t
   !< One triangle in the triangulation. Vertex indices reference the parent
   !< `triangulation_t%px / py` arrays; neighbor indices reference the parent's
   !< `tri` array, with 0 meaning "boundary edge, no neighbor". The
   !< `constrained` flag, parallel to `nbr`, marks edges that must be preserved
   !< across Delaunay-restoration sweeps (set by `cdt_build` after constraint
   !< recovery; left .false. by `dt_build`).
   integer(I4P) :: v(3)         = 0_I4P  !< Vertex indices into px/py (CCW).
   integer(I4P) :: nbr(3)       = 0_I4P  !< Neighbor triangle indices; 0 = boundary.
   logical      :: constrained(3) = .false. !< True for edges that recovery has marked as fixed.
   logical      :: active       = .false. !< False for triangles deleted during a cavity step.
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

   call dt_insert_points(tri=tri, points=points)
   call strip_super_triangle(tri=tri)
   endsubroutine dt_build

   subroutine cdt_build(tri, points, segments, status)
   !< Build the **constrained** Delaunay triangulation of an input 2D point
   !< cloud + a list of constraint segments (Sloan 1993).
   !<
   !< `points(2, N)` is the input point cloud; `segments(2, M)` holds pairs of
   !< point-array indices `(va, vb)` for each constraint segment. Both
   !< endpoints must lie in `[1, N]`. On return every constraint segment
   !< appears as an edge in the output triangulation.
   !<
   !< Algorithm (Sloan 1993, three phases on top of dt_build's foundation):
   !<   1. **Unconstrained Delaunay**: insert all points, leaving the super-
   !<      triangle in place so constraint endpoints near the input bbox don't
   !<      walk off a stripped boundary.
   !<   2. **Constraint recovery**: for each segment (va, vb): if it's already
   !<      a triangulation edge, mark it constrained; otherwise locate the
   !<      crossing-edges chain and swap diagonals on convex quads until vb
   !<      becomes reachable from va along an edge.
   !<   3. **Delaunay restoration**: sweep all non-constrained edges; flip any
   !<      that violate the local empty-circumcircle property; iterate until
   !<      stable. Constrained edges are exempt — they may end up non-Delaunay,
   !<      which is the price of honoring them.
   !<   4. **Strip super-triangle** at the very end (after recovery, since the
   !<      walk needs the super-bounded triangulation to be intact).
   !<
   !< Crossing-constraints policy: this MVP rejects with `DT_STATUS_CROSSING_CONSTRAINTS`
   !< when the recovery walk for a new segment hits an already-marked
   !< constrained edge. A future enhancement could insert a Steiner point at
   !< the intersection and split both segments; the issue text for §1.1 notes
   !< the snap-rounding requirement that would underpin such a path.
   type(triangulation_t), intent(out)           :: tri          !< Output CDT.
   real(R8P),             intent(in)            :: points(:,:)  !< Input points, shape (2, N).
   integer(I4P),          intent(in)            :: segments(:,:) !< Constraint endpoints, shape (2, M).
   integer(I4P),          intent(out), optional :: status       !< DT_STATUS_* code.
   integer(I4P)                                 :: m, k, va, vb, st_local, n_iter
   integer(I4P), parameter                      :: MAX_RESTORATION_ITER = 1000_I4P

   if (present(status)) status = DT_STATUS_OK

   ! Phase 1: unconstrained Delaunay (super-triangle still installed).
   call dt_insert_points(tri=tri, points=points)

   ! Phase 2: recover each constraint.
   m = size(segments, dim=2)
   do k = 1, m
      va = segments(1, k) ; vb = segments(2, k)
      if (va < 1 .or. va > tri%n_input .or. vb < 1 .or. vb > tri%n_input .or. va == vb) then
         if (present(status)) status = DT_STATUS_BAD_SEGMENT
         call strip_super_triangle(tri=tri)
         return
      endif
      call recover_segment(tri=tri, va=va, vb=vb, status=st_local)
      if (st_local /= DT_STATUS_OK) then
         if (present(status)) status = st_local
         call strip_super_triangle(tri=tri)
         return
      endif
   enddo

   ! Phase 3: restore Delaunay in non-constrained regions.
   ! Bounded loop guards against pathological non-convergence; in practice
   ! restoration converges in O(N) sweeps per Lawson's algorithm.
   do n_iter = 1, MAX_RESTORATION_ITER
      if (.not. restoration_sweep(tri=tri)) exit
   enddo

   ! Phase 4: strip the super-triangle.
   call strip_super_triangle(tri=tri)
   endsubroutine cdt_build

   subroutine dt_insert_points(tri, points)
   !< Internal: build a Delaunay triangulation **without stripping the super-
   !< triangle** at the end. Used by both `dt_build` (which then strips) and
   !< `cdt_build` (which strips after constraint recovery).
   type(triangulation_t), intent(out) :: tri
   real(R8P),             intent(in)  :: points(:,:)
   integer(I4P)                       :: n, i, host, p_idx
   integer(I4P), allocatable          :: cavity(:)
   integer(I4P)                       :: n_cavity

   n = size(points, dim=2)
   tri%n_input = n

   allocate(tri%px(n + SUPER_OFFSET), tri%py(n + SUPER_OFFSET))
   tri%px(1:n) = points(1, 1:n)
   tri%py(1:n) = points(2, 1:n)

   tri%n_tri_alloc = max(64_I4P, 6_I4P * (n + SUPER_OFFSET))
   allocate(tri%tri(tri%n_tri_alloc))
   tri%n_tri_used = 0_I4P

   call install_super_triangle(tri=tri)

   allocate(cavity(tri%n_tri_alloc))

   do i = 1, n
      p_idx = i
      host = locate_triangle(tri=tri, x=tri%px(p_idx), y=tri%py(p_idx))
      if (host == 0_I4P) cycle
      call collect_cavity(tri=tri, p_idx=p_idx, host=host, cavity=cavity, n_cavity=n_cavity)
      call retriangulate_cavity(tri=tri, p_idx=p_idx, cavity=cavity, n_cavity=n_cavity)
   enddo
   endsubroutine dt_insert_points

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

   subroutine recover_segment(tri, va, vb, status)
   !< Recover a single constraint segment (va, vb) so it appears as an edge
   !< in `tri`. Marks both directions of the resulting edge as constrained.
   !<
   !< Algorithm (Sloan 1993, MVP variant):
   !<   1. If (va, vb) is already an edge, mark and return.
   !<   2. Outer loop: rebuild the crossing list from scratch (the cheapest
   !<      way to keep references valid across edge flips, which can re-label
   !<      vertex sets and neighbor pointers); scan the list for the first
   !<      convex-quad crossing edge and flip it.
   !<   3. Repeat until the segment becomes an edge or no convex flip is
   !<      available.
   !<   4. Crossing-constraints policy: if the list contains an already-
   !<      constrained edge, return DT_STATUS_CROSSING_CONSTRAINTS.
   !<   5. If the loop terminates without installing the segment, return
   !<      DT_STATUS_RECOVERY_FAILED.
   !<
   !< Known MVP limitation: this convex-flip-greedy variant can spin
   !< indefinitely on inputs where every available convex flip produces a
   !< new diagonal that also crosses (va, vb), because the chain length
   !< doesn't strictly decrease. The full Sloan procedure handles this by
   !< processing the chain in *segment-order* (closest crossing to va first)
   !< and choosing flips that demonstrably advance toward vb; that variant
   !< is left for a future enhancement. For the §1.1 MVP this is acceptable:
   !< the boolean driver can detect DT_STATUS_RECOVERY_FAILED and either
   !< snap-round its cut endpoints or fall back to a coarser approach. Hand-
   !< crafted constraints (the cubes-style boolean test cases) succeed
   !< robustly; arbitrary random-cloud inputs may not.
   type(triangulation_t), intent(inout)   :: tri
   integer(I4P),          intent(in)      :: va, vb
   integer(I4P),          intent(out)     :: status
   integer(I4P)                           :: t, e, t2, e2, n_outer, k
   integer(I4P), allocatable              :: chain_t(:), chain_e(:)
   integer(I4P)                           :: n_chain
   logical                                :: progressed
   integer(I4P), parameter                :: MAX_OUTER = 100000_I4P

   status = DT_STATUS_OK

   ! Quick check: is (va, vb) already an edge?
   call find_edge(tri=tri, va=va, vb=vb, t=t, e=e)
   if (t /= 0_I4P) then
      call mark_edge_constrained(tri=tri, t=t, e=e)
      return
   endif

   allocate(chain_t(tri%n_tri_used), chain_e(tri%n_tri_used))

   do n_outer = 1, MAX_OUTER
      ! Has the segment become an edge? Check before rebuilding the chain
      ! so that the success case doesn't depend on `build_crossing_list`'s
      ! handling of boundary configurations.
      call find_edge(tri=tri, va=va, vb=vb, t=t, e=e)
      if (t /= 0_I4P) then
         call mark_edge_constrained(tri=tri, t=t, e=e)
         return
      endif
      ! Rebuild crossing list each pass — this is the simplest way to keep
      ! triangle / edge references valid across flips, which mutate vertex
      ! lists and neighbor pointers in-place. Asymptotically O(crossings^2)
      ! per recovery, but for small chains (typical §1.1 use case) it is
      ! negligible.
      call build_crossing_list(tri=tri, va=va, vb=vb, &
                               chain_t=chain_t, chain_e=chain_e, n_chain=n_chain, status=status)
      if (status /= DT_STATUS_OK) return
      if (n_chain == 0_I4P) then
         ! No crossings but segment isn't an edge — geometric inconsistency.
         status = DT_STATUS_RECOVERY_FAILED
         return
      endif
      ! Find the first crossing whose quad is convex (flip-valid).
      progressed = .false.
      do k = 1, n_chain
         t = chain_t(k) ; e = chain_e(k)
         if (tri%tri(t)%constrained(e)) then
            status = DT_STATUS_CROSSING_CONSTRAINTS
            return
         endif
         t2 = tri%tri(t)%nbr(e)
         if (t2 == 0_I4P) then
            status = DT_STATUS_RECOVERY_FAILED
            return
         endif
         e2 = neighbor_edge(tri=tri, t=t, e=e)
         if (is_quad_convex(tri=tri, t=t, e=e, t2=t2, e2=e2)) then
            call flip_edge(tri=tri, t=t, e=e)
            progressed = .true.
            exit  ! restart outer loop with fresh chain
         endif
      enddo
      if (.not. progressed) then
         status = DT_STATUS_RECOVERY_FAILED
         return
      endif
   enddo
   status = DT_STATUS_RECOVERY_FAILED
   endsubroutine recover_segment

   subroutine build_crossing_list(tri, va, vb, chain_t, chain_e, n_chain, status)
   !< Walk from va toward vb across the triangulation, recording every edge
   !< the segment crosses. Output is a list of (triangle, local-edge) pairs.
   !<
   !< Strategy: starting from a triangle incident to va whose interior
   !< direction contains vb, find the edge AB exits through (the one whose
   !< endpoints are on opposite sides of AB). Cross to the next triangle.
   !< Repeat until the next triangle is incident to vb.
   type(triangulation_t), intent(in)  :: tri
   integer(I4P),          intent(in)  :: va, vb
   integer(I4P),          intent(out) :: chain_t(:), chain_e(:)
   integer(I4P),          intent(out) :: n_chain, status
   integer(I4P)                       :: t, e, t_next, k, e_next
   real(R8P)                          :: ax, ay, bx, by
   integer(I4P)                       :: u, v, found_e, prev_t
   integer(I4P), parameter            :: MAX_CHAIN_LEN = 100000_I4P
   real(R8P)                          :: s_u, s_v

   status = DT_STATUS_OK
   n_chain = 0_I4P
   ax = tri%px(va) ; ay = tri%py(va)
   bx = tri%px(vb) ; by = tri%py(vb)

   ! Find the starting triangle: incident to va, with its edge opposite va
   ! crossed by AB.
   call find_crossing_edge(tri=tri, va=va, vb=vb, t=t, e=e, status=status)
   if (status /= DT_STATUS_OK) return
   if (t == 0_I4P) return  ! caller should have handled "already an edge" case

   prev_t = 0_I4P
   do k = 1, MAX_CHAIN_LEN
      n_chain = n_chain + 1
      if (n_chain > size(chain_t)) then
         status = DT_STATUS_RECOVERY_FAILED  ! chain too long for buffer
         return
      endif
      chain_t(n_chain) = t
      chain_e(n_chain) = e
      ! Step to the neighbor across edge e.
      t_next = tri%tri(t)%nbr(e)
      if (t_next == 0_I4P) then
         status = DT_STATUS_RECOVERY_FAILED  ! AB exits the triangulation
         return
      endif
      ! If t_next is incident to vb, we are done — no more edges to cross
      ! beyond this one.
      if (any(tri%tri(t_next)%v == vb)) return
      ! Otherwise find the next exit edge of t_next: the edge whose endpoints
      ! are on opposite sides of AB AND which is not the entry edge.
      found_e = 0_I4P
      do e_next = 1, 3
         u = tri%tri(t_next)%v(mod(e_next, 3) + 1)
         v = tri%tri(t_next)%v(mod(e_next + 1, 3) + 1)
         ! Skip the entry edge (which is the one shared with t).
         if (tri%tri(t_next)%nbr(e_next) == t) cycle
         s_u = side2d(ax, ay, bx, by, tri%px(u), tri%py(u))
         s_v = side2d(ax, ay, bx, by, tri%px(v), tri%py(v))
         if (s_u * s_v < 0._R8P) then
            ! Crossing — but also confirm AB actually reaches this edge:
            ! both endpoints u,v are on opposite sides of AB AND vb is on
            ! the *far* side of (u,v) from va (otherwise vb is in t_next
            ! and we should stop).
            ! Test by checking sign of vb wrt the edge UV and comparing to
            ! the sign of va.
            if (segment_crosses(ax=ax, ay=ay, bx=bx, by=by, &
                                ux=tri%px(u), uy=tri%py(u), vx=tri%px(v), vy=tri%py(v))) then
               found_e = e_next ; exit
            endif
         endif
      enddo
      if (found_e == 0_I4P) then
         ! No further crossing — vb lies in the interior of t_next. The
         ! chain is complete; return what we have so far.
         return
      endif
      prev_t = t
      t = t_next
      e = found_e
   enddo
   status = DT_STATUS_RECOVERY_FAILED
   ! Suppress unused-variable warning.
   k = prev_t  ! never used downstream
   endsubroutine build_crossing_list

   pure function segment_crosses(ax, ay, bx, by, ux, uy, vx, vy) result(yes)
   !< True if segment AB crosses segment UV (proper crossing, endpoints excluded
   !< from being on each other's lines).
   real(R8P), intent(in) :: ax, ay, bx, by, ux, uy, vx, vy
   logical               :: yes
   real(R8P)             :: s_u, s_v, s_a, s_b

   s_u = side2d(ax, ay, bx, by, ux, uy)
   s_v = side2d(ax, ay, bx, by, vx, vy)
   s_a = side2d(ux, uy, vx, vy, ax, ay)
   s_b = side2d(ux, uy, vx, vy, bx, by)
   yes = (s_u * s_v < 0._R8P) .and. (s_a * s_b < 0._R8P)
   endfunction segment_crosses

   pure subroutine find_edge(tri, va, vb, t, e)
   !< Search for an edge with endpoints (va, vb). Returns the first triangle
   !< containing this edge and the edge index (1..3). Returns (0, 0) if not
   !< found.
   type(triangulation_t), intent(in)  :: tri
   integer(I4P),          intent(in)  :: va, vb
   integer(I4P),          intent(out) :: t, e
   integer(I4P)                       :: it, ie, u, v

   t = 0_I4P ; e = 0_I4P
   do it = 1, tri%n_tri_used
      if (.not. tri%tri(it)%active) cycle
      do ie = 1, 3
         u = tri%tri(it)%v(mod(ie, 3) + 1)
         v = tri%tri(it)%v(mod(ie + 1, 3) + 1)
         if ((u == va .and. v == vb) .or. (u == vb .and. v == va)) then
            t = it ; e = ie ; return
         endif
      enddo
   enddo
   endsubroutine find_edge

   subroutine mark_edge_constrained(tri, t, e)
   !< Mark edge `e` of triangle `t` as constrained, and the corresponding
   !< edge of its neighbor (so the flag is symmetric across the shared edge).
   type(triangulation_t), intent(inout) :: tri
   integer(I4P),          intent(in)    :: t, e
   integer(I4P)                         :: t2, e2

   tri%tri(t)%constrained(e) = .true.
   t2 = tri%tri(t)%nbr(e)
   if (t2 /= 0_I4P) then
      e2 = neighbor_edge(tri=tri, t=t, e=e)
      if (e2 /= 0_I4P) tri%tri(t2)%constrained(e2) = .true.
   endif
   endsubroutine mark_edge_constrained

   pure function neighbor_edge(tri, t, e) result(e2)
   !< Given that triangle `t` has neighbor `tri%nbr(t,e)` across edge `e`,
   !< find which of the neighbor's three edges (1..3) is the shared one.
   !< Returns 0 if no neighbor (boundary edge).
   type(triangulation_t), intent(in) :: tri
   integer(I4P),          intent(in) :: t, e
   integer(I4P)                      :: e2
   integer(I4P)                      :: t2, k

   e2 = 0_I4P
   t2 = tri%tri(t)%nbr(e)
   if (t2 == 0_I4P) return
   do k = 1, 3
      if (tri%tri(t2)%nbr(k) == t) then
         e2 = k ; return
      endif
   enddo
   endfunction neighbor_edge

   subroutine find_crossing_edge(tri, va, vb, t, e, status)
   !< Find the **starting triangle** for a segment-crossing walk from va
   !< toward vb: the unique triangle incident to va whose interior contains
   !< the ray va→vb. Returns that triangle and the index of its edge
   !< opposite va (which is the first edge AB crosses).
   !<
   !< Selection rule: among the triangles incident to va, the correct one
   !< is the triangle (va, u, v) such that vb lies in the angular sector
   !< bounded by rays va→u (CCW boundary) and va→v (CW boundary). Equivalent
   !< 2D test using orient2d:
   !<   - vb is on the CCW side of (va→u) ⇔ orient2d(va, u, vb) > 0
   !<   - vb is on the CW  side of (va→v) ⇔ orient2d(va, v, vb) < 0
   !< Both conditions hold simultaneously for exactly one incident triangle
   !< (assuming va is not on the boundary and vb is in general position).
   !<
   !< Boundary handling: if va is on the convex hull, vb may lie *outside*
   !< the angular fan of incident triangles — in that case we return t=0,
   !< letting the caller declare the constraint unrecoverable.
   type(triangulation_t), intent(in)  :: tri
   integer(I4P),          intent(in)  :: va, vb
   integer(I4P),          intent(out) :: t, e, status
   integer(I4P)                       :: it, ie, vt(3), va_pos, u, v
   real(R8P)                          :: ax, ay, bx, by, ux, uy, vx, vy
   real(R8P)                          :: s_left, s_right

   status = DT_STATUS_OK
   t = 0_I4P ; e = 0_I4P
   ax = tri%px(va) ; ay = tri%py(va)
   bx = tri%px(vb) ; by = tri%py(vb)

   ! Iterate active triangles incident to va.
   do it = 1, tri%n_tri_used
      if (.not. tri%tri(it)%active) cycle
      vt = tri%tri(it)%v
      ! Find local index of va in this triangle.
      va_pos = 0_I4P
      do ie = 1, 3
         if (vt(ie) == va) then ; va_pos = ie ; exit ; endif
      enddo
      if (va_pos == 0_I4P) cycle
      ! In CCW order with va at va_pos, the next vertex (mod 3) is the "u"
      ! that bounds the CCW side of va, and the previous vertex bounds the
      ! CW side. The edge opposite va is the one between these two.
      u = vt(mod(va_pos, 3) + 1)        ! next CCW from va
      v = vt(mod(va_pos + 1, 3) + 1)    ! the one after that (= prev CW from va)
      ux = tri%px(u) ; uy = tri%py(u)
      vx = tri%px(v) ; vy = tri%py(v)
      ! vb lies in the wedge va→u to va→v (CCW) iff:
      !   orient2d(va, u, vb) >= 0 (vb on CCW side of va→u, including the ray)
      !   orient2d(va, v, vb) <= 0 (vb on CW  side of va→v, including the ray)
      ! Use side2d (positive = left of directed line).
      s_left  = side2d(ax, ay, ux, uy, bx, by)  ! sign of vb wrt va→u
      s_right = side2d(ax, ay, vx, vy, bx, by)  ! sign of vb wrt va→v
      if (s_left >= 0._R8P .and. s_right <= 0._R8P) then
         ! vb is in this triangle's angular sector at va. The edge opposite
         ! va is the one to cross; its index is the position OPPOSITE va_pos.
         t = it
         e = va_pos  ! edge opposite vertex va = edge index va_pos in our convention
         return
      endif
   enddo
   endsubroutine find_crossing_edge

   pure function side2d(ax, ay, bx, by, px, py) result(s)
   !< Sign-of-side: positive if (px, py) is left of directed line (ax,ay)->(bx,by).
   real(R8P), intent(in) :: ax, ay, bx, by, px, py
   real(R8P)             :: s

   s = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
   endfunction side2d

   pure function is_quad_convex(tri, t, e, t2, e2) result(yes)
   !< Test whether the quadrilateral formed by triangles t and t2 (sharing
   !< edge e of t / e2 of t2) is strictly convex. A diagonal flip is valid
   !< iff the quad is convex.
   !<
   !< Quad vertices in CCW order: a (apex of t opposite e), b, d (apex of t2
   !< opposite e2), c — where (b, c) is the shared edge as seen from t.
   !< Strictly convex iff a and d lie on opposite sides of (b, c) AND
   !< b and c lie on opposite sides of (a, d). The first condition is
   !< guaranteed by the triangulation invariant; we test the second.
   type(triangulation_t), intent(in) :: tri
   integer(I4P),          intent(in) :: t, e, t2, e2
   logical                           :: yes
   integer(I4P)                      :: a, b, c, d
   real(R8P)                         :: ax, ay, bx, by, cx, cy, dx_, dy_
   real(R8P)                         :: s_b, s_c

   a  = tri%tri(t )%v(e )
   d  = tri%tri(t2)%v(e2)
   b  = tri%tri(t )%v(mod(e , 3) + 1)
   c  = tri%tri(t )%v(mod(e + 1, 3) + 1)
   ax = tri%px(a) ; ay = tri%py(a)
   bx = tri%px(b) ; by = tri%py(b)
   cx = tri%px(c) ; cy = tri%py(c)
   dx_ = tri%px(d) ; dy_ = tri%py(d)
   s_b = side2d(ax, ay, dx_, dy_, bx, by)
   s_c = side2d(ax, ay, dx_, dy_, cx, cy)
   yes = (s_b * s_c < 0._R8P)
   endfunction is_quad_convex

   subroutine flip_edge(tri, t, e)
   !< Flip the diagonal of the quad formed by triangle t (with edge e shared)
   !< and its neighbor t2 across that edge. After the flip, the shared edge
   !< connects the apices a and d (instead of b and c).
   !<
   !< Before:                    After:
   !<     a                          a
   !<    /|\                        / \
   !<   / | \                      /   \
   !<  b--+--c   (e is bc)        b-----c
   !<   \ | /                      \   /
   !<    \|/                        \ /
   !<     d                          d
   !<
   !< T = (a, b, c) with edge e = (b,c) opposite a.
   !< T2 = (b, d, c)? — depends on neighbor's vertex order. We extract the
   !< apex d explicitly via (t2, e2).
   !<
   !< New triangles:
   !<   T' = (a, b, d) -- shares edge (b, d) with caller's neighbor across the
   !<                     b-side of T2; shares (d, a) with new T2'; shares
   !<                     (a, b) with caller's neighbor across the b-side of T.
   !<   T2' = (a, d, c) -- analogous.
   !<
   !< Neighbor pointers: there are 4 "outside" neighbors that need re-pointing
   !< plus the new shared edge between T' and T2'.
   type(triangulation_t), intent(inout) :: tri
   integer(I4P),          intent(in)    :: t, e
   integer(I4P)                         :: t2, e2
   integer(I4P)                         :: a, b, c, d
   integer(I4P)                         :: nb_t_b, nb_t_c, nb_t2_b, nb_t2_c
   integer(I4P)                         :: e_t_b, e_t_c, e_t2_b, e_t2_c
   logical                              :: cb_t_b, cb_t_c, cb_t2_b, cb_t2_c

   t2 = tri%tri(t)%nbr(e)
   e2 = neighbor_edge(tri=tri, t=t, e=e)

   ! Identify the four quad vertices.
   a = tri%tri(t )%v(e)
   b = tri%tri(t )%v(mod(e , 3) + 1)
   c = tri%tri(t )%v(mod(e + 1, 3) + 1)
   d = tri%tri(t2)%v(e2)

   ! Identify the four outside neighbors (and which edge they're connected via)
   ! BEFORE we mutate t and t2.
   !   In T, edge opposite b is the one whose endpoints are (c, a). Index =
   !     position k where v(k) == b → the edge opposite that vertex.
   e_t_b = which_vertex(tri=tri, t=t, vid=b)
   e_t_c = which_vertex(tri=tri, t=t, vid=c)
   nb_t_b = tri%tri(t)%nbr(e_t_b) ; cb_t_b = tri%tri(t)%constrained(e_t_b)
   nb_t_c = tri%tri(t)%nbr(e_t_c) ; cb_t_c = tri%tri(t)%constrained(e_t_c)
   e_t2_b = which_vertex(tri=tri, t=t2, vid=b)
   e_t2_c = which_vertex(tri=tri, t=t2, vid=c)
   nb_t2_b = tri%tri(t2)%nbr(e_t2_b) ; cb_t2_b = tri%tri(t2)%constrained(e_t2_b)
   nb_t2_c = tri%tri(t2)%nbr(e_t2_c) ; cb_t2_c = tri%tri(t2)%constrained(e_t2_c)

   ! Reconfigure T as T' = (a, b, d) in CCW order (assuming original CCW).
   ! Edge opposite a = (b, d), opposite b = (d, a), opposite d = (a, b).
   tri%tri(t)%v = [a, b, d]
   ! Edge 1 of new T = opposite a = (b, d) → originally an edge of T2 between
   ! b and d. That was edge e_t2_c of T2 (the edge opposite c in T2). Its
   ! external neighbor was nb_t2_c; constrained flag cb_t2_c.
   tri%tri(t)%nbr(1) = nb_t2_c ; tri%tri(t)%constrained(1) = cb_t2_c
   ! Edge 2 = opposite b = (d, a) → the new flipped diagonal. Neighbor is T2'.
   tri%tri(t)%nbr(2) = t2 ; tri%tri(t)%constrained(2) = .false.
   ! Edge 3 = opposite d = (a, b) → originally edge of T between a and b.
   ! That was edge e_t_c of T (opposite c).
   tri%tri(t)%nbr(3) = nb_t_c ; tri%tri(t)%constrained(3) = cb_t_c

   ! Reconfigure T2 as T2' = (a, d, c) in CCW order.
   tri%tri(t2)%v = [a, d, c]
   ! Edge 1 = opposite a = (d, c) → originally edge of T2 between d and c.
   ! That was edge e_t2_b of T2 (opposite b).
   tri%tri(t2)%nbr(1) = nb_t2_b ; tri%tri(t2)%constrained(1) = cb_t2_b
   ! Edge 2 = opposite d = (c, a) → originally edge of T between c and a.
   ! That was edge e_t_b of T (opposite b).
   tri%tri(t2)%nbr(2) = nb_t_b ; tri%tri(t2)%constrained(2) = cb_t_b
   ! Edge 3 = opposite c = (a, d) → the new flipped diagonal, shared with T.
   tri%tri(t2)%nbr(3) = t ; tri%tri(t2)%constrained(3) = .false.

   ! Patch external neighbors' back-pointers. Each external neighbor used to
   ! point to t or t2 via some edge; we redirect it to the new triangle that
   ! now owns that boundary edge.
   if (nb_t2_c /= 0_I4P) call patch_neighbor(tri=tri, side=nb_t2_c, old=t2, new=t)
   if (nb_t_c  /= 0_I4P) call patch_neighbor(tri=tri, side=nb_t_c , old=t , new=t)
   if (nb_t2_b /= 0_I4P) call patch_neighbor(tri=tri, side=nb_t2_b, old=t2, new=t2)
   if (nb_t_b  /= 0_I4P) call patch_neighbor(tri=tri, side=nb_t_b , old=t , new=t2)
   endsubroutine flip_edge

   pure function which_vertex(tri, t, vid) result(k)
   !< Return the local index (1..3) of vertex `vid` in triangle `t`.
   !< Returns 0 if not present (defensive).
   type(triangulation_t), intent(in) :: tri
   integer(I4P),          intent(in) :: t, vid
   integer(I4P)                      :: k

   do k = 1, 3
      if (tri%tri(t)%v(k) == vid) return
   enddo
   k = 0_I4P
   endfunction which_vertex

   function restoration_sweep(tri) result(any_flipped)
   !< Single Lawson restoration sweep: visit each non-constrained interior
   !< edge; if its quad is flippable AND the alternate diagonal is preferred
   !< by the empty-circumcircle test, flip. Returns .true. if any flip
   !< happened (caller iterates until no flips occur).
   type(triangulation_t), intent(inout) :: tri
   logical                              :: any_flipped
   integer(I4P)                         :: t, e, t2, e2
   integer(I4P)                         :: a, b, c, d
   type(vector_R8P)                     :: va, vb, vc, vd
   real(R8P)                            :: incirc

   any_flipped = .false.

   do t = 1, tri%n_tri_used
      if (.not. tri%tri(t)%active) cycle
      do e = 1, 3
         if (tri%tri(t)%constrained(e)) cycle
         t2 = tri%tri(t)%nbr(e)
         if (t2 == 0_I4P) cycle  ! boundary
         ! Process each edge once: only handle when t < t2 to avoid double-visiting.
         if (t2 < t) cycle
         e2 = neighbor_edge(tri=tri, t=t, e=e)
         if (.not. is_quad_convex(tri=tri, t=t, e=e, t2=t2, e2=e2)) cycle
         ! Check if flipping improves Delaunay. Quad apex d should NOT lie
         ! inside circumcircle of (a, b, c); if it does, flip.
         a = tri%tri(t )%v(e )
         d = tri%tri(t2)%v(e2)
         b = tri%tri(t )%v(mod(e , 3) + 1)
         c = tri%tri(t )%v(mod(e + 1, 3) + 1)
         va = vector_R8P(tri%px(a), tri%py(a), 0._R8P)
         vb = vector_R8P(tri%px(b), tri%py(b), 0._R8P)
         vc = vector_R8P(tri%px(c), tri%py(c), 0._R8P)
         vd = vector_R8P(tri%px(d), tri%py(d), 0._R8P)
         incirc = incircle(a=va, b=vb, c=vc, d=vd)
         if (incirc > 0._R8P) then
            call flip_edge(tri=tri, t=t, e=e)
            any_flipped = .true.
         endif
      enddo
   enddo
   endfunction restoration_sweep

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
