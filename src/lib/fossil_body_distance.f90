!< FOSSIL, distance between two bodies (issue #67, roadmap #58 WP7).
module fossil_body_distance
!< FOSSIL, distance between two bodies (issue #67, roadmap #58 WP7).
!<
!< `tri_tri_distance`: the squared distance between two triangles with its witness points. It is 0 when the exact
!< classification says they meet (`pair_cut` of `fossil_tritri_exact`: crossing, touching or coplanar overlap, no
!< tolerance); otherwise the minimum of the 6 vertex-triangle and the 9 edge-edge distances (the closest pair of two
!< disjoint triangles is one of these), each as the norm of a residual, evaluated in a fixed order: the first of equal
!< values wins.
!<
!< `body_distance`: the minimum over the facet pairs of two surfaces, by a best-first descent of their two AABB trees: a
!< heap of node pairs ordered by the squared distance between their boxes (0 when they overlap), a pair dropped when its
!< bound exceeds the best distance found (with a relative margin far above round-off, so a pair whose computed distance
!< ties the best is never dropped), leaf pairs evaluated in a fixed order. Of facet pairs at the same distance the
!< lexicographically smallest (facet_a, facet_b) wins, so the result (distance, facets, points) is bit-identical to the
!< brute-force loop `body_distance_brute`. Works for both tree kinds (it uses the generic node API).
!<
!< Poses (step 2): with `pose_a` / `pose_b` the facets are evaluated where `move` would put them (`to_world` of each
!< vertex, the same bits) and each tree node is bounded by the box of its 8 moved corners, a conservative bound: the
!< trees are not refitted, and the result equals the one after `move`, bit for bit.
!<
!< Per call (issue #94) each facet is posed and boxed once (`posed_facets_t`: the same `to_world`, `minval`, `maxval` as
!< the per-pair evaluation of `body_distance_brute`). In a leaf pair, a facet's whole row is skipped when its box is
!< farther than the pruning bound from the union of the other leaf's facet boxes: rounded subtraction is monotone, so
!< that gap is <= the gap to each facet box of the row, and every pair of the row would be pruned (d2 only decreases
!< along the row). The exact `pair_cut` is skipped when the two triangle boxes are strictly disjoint on some axis and
!< both triangles are regular (`facet_regular`: exactly non-degenerate, coordinates in the predicates' range): such
!< triangles share no point, so `pair_cut` (exact on regular triangles) says SEPARATED, or COPLANAR with an exact
!< `coplanar_meet` that is false. A degenerate triangle keeps `pair_cut` (its sections are not exact). Same bits.
!<
!< `body_safe_step` (issue #73, the directional bound of C2A): over the facet pairs of two bodies moving along
!< `pose_interpolate` (axis k, angle theta, translation dt per body), the largest step that keeps every pair at least
!< tol / 2 apart, by the separating-plane bound of each pair (a triangle is convex): with the witness direction n from a
!< to b at distance d_ij, the gap along n bounds d_ij from below at every s and closes at most at c_ij = n . (dt_a -
!< dt_b) + theta_a rho_a + theta_b rho_b (rho: the largest distance of the triangle's vertices to the rotation axis
!< through the body origin); the pair's step is (d_ij - tol / 2) / c_ij (unbounded if c_ij <= 0). Node pairs whose
!< plain bound (box distance - tol / 2) / (mu_a + mu_b), mu = |dt| + theta rho_max(node box), is already >= the best
!< step are pruned: a pair's directional step is never below its plain one. Depth-first with an explicit stack; the
!< result (a minimum) does not depend on the visiting order.
use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,     only : facet_object
use fossil_implicit_points,  only : implicit_point_t, ip_approx
use fossil_list_id_object,   only : list_id_object
use fossil_predicates,       only : orient2d_xy
use fossil_rigid_pose,       only : rigid_pose_t, to_world, rotate
use fossil_triangle_closest, only : triangle_closest_st
use fossil_tritri_exact,     only : pair_cut, facet_drop, PAIR_SEPARATED, PAIR_COPLANAR
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: tri_tri_distance, body_distance, body_distance_brute, closest_on_triangle, closest_segment_segment, body_safe_step

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P   !< Children buffer (octree fan-out).
real(R8P),    parameter :: PRUNE_REL = 1.e-10_R8P !< Relative margin of the pruning test (bounds are not compared
                                                  !< bit for bit with computed distances).
integer(I4P), parameter :: POSED_NONE      = 0_I4P !< Facet not posed yet.
integer(I4P), parameter :: POSED_NONFINITE = 1_I4P !< Posed, a coordinate is not finite (no skip uses it).
integer(I4P), parameter :: POSED_FINITE    = 2_I4P !< Posed, finite, regularity not tested yet.
integer(I4P), parameter :: POSED_REGULAR   = 3_I4P !< Posed, finite, exactly non-degenerate, in the predicates' range.
integer(I4P), parameter :: POSED_IRREGULAR = 4_I4P !< Posed, finite, but degenerate or out of the predicates' range.
real(R8P),    parameter :: REGULAR_MIN = 1.e-60_R8P !< Smallest nonzero |coordinate| of a regular facet (no underflow in
                                                    !< the exact predicates).
real(R8P),    parameter :: REGULAR_MAX = 1.e60_R8P  !< Largest |coordinate| of a regular facet (no overflow).

type :: posed_facets_t
   !< Facets posed once per `body_distance` call (issue #94): vertices as columns and box, filled on first use.
   real(R8P),    allocatable :: v(:, :, :) !< Posed vertices (3, 3, n), as `try_pair` builds them.
   real(R8P),    allocatable :: box(:, :)  !< Box (6, n): min x, y, z, max x, y, z of the posed vertices.
   integer(I4P), allocatable :: state(:)   !< POSED_* code.
endtype posed_facets_t

contains
   subroutine tri_tri_distance(a, b, d2, pa, pb, meet)
   !< Squared distance between triangles a and b (vertices as columns), witness points pa on a and pb on b.
   real(R8P), intent(in)  :: a(3, 3), b(3, 3) !< Triangles.
   real(R8P), intent(out) :: d2               !< Squared distance.
   real(R8P), intent(out) :: pa(3), pb(3)     !< Witness points.
   logical,   intent(out) :: meet             !< The triangles meet (exact).

   call tri_tri_distance_core(a, b, .false., d2, pa, pb, meet)
   endsubroutine tri_tri_distance

   subroutine tri_tri_distance_core(a, b, apart, d2, pa, pb, meet)
   !< `tri_tri_distance`; with `apart` the caller has proved that the triangles share no point (strictly disjoint boxes,
   !< both non-degenerate): the exact classification, which would find them separated, is skipped.
   real(R8P), intent(in)  :: a(3, 3), b(3, 3) !< Triangles.
   logical,   intent(in)  :: apart            !< The triangles share no point (proved by the caller).
   real(R8P), intent(out) :: d2               !< Squared distance.
   real(R8P), intent(out) :: pa(3), pb(3)     !< Witness points.
   logical,   intent(out) :: meet             !< The triangles meet (exact).
   type(vector_R8P)       :: ta(3), tb(3)
   type(implicit_point_t) :: p, q
   real(R8P)              :: c(3), c2(3), dd
   integer(I4P)           :: i, j, kind

   meet = .false.
   if (.not. apart) then
      do i = 1, 3
         ta(i) = vector_R8P(a(1, i), a(2, i), a(3, i)) ; tb(i) = vector_R8P(b(1, i), b(2, i), b(3, i))
      enddo
      call pair_cut(ta, tb, kind, p, q)
      meet = kind /= PAIR_SEPARATED
      if (kind == PAIR_COPLANAR) meet = coplanar_meet(a, b, facet_drop(ta)) ! same plane: do they overlap?
   endif
   if (meet) then
      d2 = 0._R8P
      if (p%kind /= 0_I4P) then
         pa = ip_approx(p)
      else ! coplanar: a vertex of one inside the other, or the closest pair below (0 anyway)
         pa = a(:, 1)
      endif
      pb = pa
      return
   endif
   d2 = huge(1._R8P)
   do i = 1, 3 ! vertices of a against b
      c = closest_on_triangle(a(:, i), b)
      dd = sum((a(:, i) - c)**2)
      if (dd < d2) then
         d2 = dd ; pa = a(:, i) ; pb = c
      endif
   enddo
   do i = 1, 3 ! vertices of b against a
      c = closest_on_triangle(b(:, i), a)
      dd = sum((b(:, i) - c)**2)
      if (dd < d2) then
         d2 = dd ; pa = c ; pb = b(:, i)
      endif
   enddo
   do i = 1, 3 ! edges against edges
      do j = 1, 3
         call closest_segment_segment(a(:, i), a(:, mod(i, 3) + 1), b(:, j), b(:, mod(j, 3) + 1), c, c2)
         dd = sum((c - c2)**2)
         if (dd < d2) then
            d2 = dd ; pa = c ; pb = c2
         endif
      enddo
   enddo
   endsubroutine tri_tri_distance_core

   function coplanar_meet(a, b, drop) result(meet)
   !< Two exactly coplanar triangles share a point (exact orient2d in the projection that drops `drop`): a vertex of one
   !< in the closed other, or two edges that meet.
   real(R8P),    intent(in) :: a(3, 3), b(3, 3)
   integer(I4P), intent(in) :: drop
   logical                  :: meet
   integer(I4P)             :: i, j, u, v

   meet = .false.
   if (drop == 0_I4P) return
   u = mod(drop, 3_I4P) + 1_I4P ; v = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   do i = 1, 3
      if (inside(a(:, i), b) .or. inside(b(:, i), a)) then
         meet = .true. ; return
      endif
   enddo
   do i = 1, 3
      do j = 1, 3
         if (segments_meet(a(:, i), a(:, mod(i, 3) + 1), b(:, j), b(:, mod(j, 3) + 1))) then
            meet = .true. ; return
         endif
      enddo
   enddo
   contains
      function o2(p, q, r) result(sg)
      real(R8P), intent(in) :: p(3), q(3), r(3)
      integer(I4P)          :: sg
      real(R8P)             :: d

      d = orient2d_xy(p(u), p(v), q(u), q(v), r(u), r(v))
      sg = 0_I4P
      if (d > 0._R8P) sg = 1_I4P
      if (d < 0._R8P) sg = -1_I4P
      endfunction o2

      function inside(p, t) result(yes)
      real(R8P), intent(in) :: p(3), t(3, 3)
      logical               :: yes
      integer(I4P)          :: s1, s2, s3

      s1 = o2(t(:, 1), t(:, 2), p) ; s2 = o2(t(:, 2), t(:, 3), p) ; s3 = o2(t(:, 3), t(:, 1), p)
      yes = (s1 >= 0 .and. s2 >= 0 .and. s3 >= 0) .or. (s1 <= 0 .and. s2 <= 0 .and. s3 <= 0)
      endfunction inside

      function segments_meet(p1, q1, p2, q2) result(yes)
      real(R8P), intent(in) :: p1(3), q1(3), p2(3), q2(3)
      logical               :: yes
      integer(I4P)          :: d1, d2, d3, d4

      d1 = o2(p2, q2, p1) ; d2 = o2(p2, q2, q1) ; d3 = o2(p1, q1, p2) ; d4 = o2(p1, q1, q2)
      yes = d1 * d2 < 0 .and. d3 * d4 < 0
      if (yes) return
      yes = (d1 == 0 .and. on(p2, q2, p1)) .or. (d2 == 0 .and. on(p2, q2, q1)) .or. &
            (d3 == 0 .and. on(p1, q1, p2)) .or. (d4 == 0 .and. on(p1, q1, q2))
      endfunction segments_meet

      pure function on(p, q, r) result(yes)
      !< r (collinear with p, q) on the closed segment (p, q).
      real(R8P), intent(in) :: p(3), q(3), r(3)
      logical               :: yes

      yes = r(u) >= min(p(u), q(u)) .and. r(u) <= max(p(u), q(u)) .and. r(v) >= min(p(v), q(v)) .and. &
            r(v) <= max(p(v), q(v))
      endfunction on
   endfunction coplanar_meet

   pure function closest_on_triangle(p, t) result(c)
   !< Closest point of triangle t (vertices as columns) to p: the point kernel of the distance queries
   !< (`triangle_closest_st`: Eberly's solve, or for needle triangles the residual form `robust_st`, issue #39), so the
   !< distance is accurate to ~eps L on needles too (issue #69; the Voronoi-region form of Ericson 5.1.5 used before
   !< lost up to 1e-10 on needles of aspect 1e8).
   real(R8P), intent(in) :: p(3), t(3, 3)
   real(R8P)             :: c(3)
   type(vector_R8P)      :: v1, e12, e13, q
   real(R8P)             :: a, b, cc, det, sq, tq, d2

   v1%x = t(1, 1) ; v1%y = t(2, 1) ; v1%z = t(3, 1)
   e12%x = t(1, 2) - t(1, 1) ; e12%y = t(2, 2) - t(2, 1) ; e12%z = t(3, 2) - t(3, 1)
   e13%x = t(1, 3) - t(1, 1) ; e13%y = t(2, 3) - t(2, 1) ; e13%z = t(3, 3) - t(3, 1)
   q%x = p(1) ; q%y = p(2) ; q%z = p(3)
   a = e12%x * e12%x + e12%y * e12%y + e12%z * e12%z
   b = e12%x * e13%x + e12%y * e13%y + e12%z * e13%z
   cc = e13%x * e13%x + e13%y * e13%y + e13%z * e13%z
   det = a * cc - b * b
   call triangle_closest_st(v1=v1, e12=e12, e13=e13, a=a, b=b, c=cc, det=det, point=q, sq=sq, tq=tq, distance=d2)
   c = t(:, 1) + sq * (t(:, 2) - t(:, 1)) + tq * (t(:, 3) - t(:, 1))
   endfunction closest_on_triangle

   pure subroutine closest_segment_segment(p1, q1, p2, q2, c1, c2)
   !< Closest points c1 on segment (p1, q1) and c2 on (p2, q2), stable for near-parallel segments (issue #69; Ericson
   !< 5.1.9 used before solves through a e - b^2 = |d1 x d2|^2 and lost ~eps L / angle). The longer segment is the
   !< reference (say 1, direction d1). The distance from x2(t) to segment 1 is |P (x2(t) - p1)|, P the projector on the
   !< plane normal to d1, while the foot s(t) = d1 . (x2(t) - p1) / |d1|^2 (linear in t) is in [0, 1]; elsewhere it
   !< is the distance to p1 or q1. So the minimum is the least of: a point-to-segment distance in that plane, over
   !< the t interval where s(t) is in [0, 1] (projections cost ~eps L, and a tiny projected segment makes the error
   !< of t harmless), and the distances from p1 and q1 to segment 2. Ties keep the first candidate in that order.
   real(R8P), intent(in)  :: p1(3), q1(3), p2(3), q2(3)
   real(R8P), intent(out) :: c1(3), c2(3)
   real(R8P)              :: a1(3), b1(3), a2(3), b2(3), d1(3), d2(3), aa, s0, s1, tlo, thi, t, s, w0(3), e(3), ee
   real(R8P)              :: best, dd, x(3)
   logical                :: swap

   swap = sum((q2 - p2)**2) > sum((q1 - p1)**2)
   if (swap) then
      a1 = p2 ; b1 = q2 ; a2 = p1 ; b2 = q1
   else
      a1 = p1 ; b1 = q1 ; a2 = p2 ; b2 = q2
   endif
   d1 = b1 - a1 ; d2 = b2 - a2 ; aa = dot_product(d1, d1)
   if (aa <= 0._R8P) then ! both degenerate
      c1 = p1 ; c2 = p2 ; return
   endif
   best = huge(1._R8P)
   ! interior candidate: t where the foot on segment 1 is inside, s(t) = s0 + t (s1 - s0)
   s0 = dot_product(d1, a2 - a1) / aa ; s1 = dot_product(d1, b2 - a1) / aa
   tlo = 0._R8P ; thi = 1._R8P
   if (s1 /= s0) then
      tlo = max(0._R8P, min((0._R8P - s0) / (s1 - s0), (1._R8P - s0) / (s1 - s0)))
      thi = min(1._R8P, max((0._R8P - s0) / (s1 - s0), (1._R8P - s0) / (s1 - s0)))
   elseif (s0 < 0._R8P .or. s0 > 1._R8P) then
      tlo = 1._R8P ; thi = 0._R8P ! never inside
   endif
   if (tlo <= thi) then
      w0 = proj(a2 - a1) ; e = proj(d2) ; ee = dot_product(e, e)
      t = tlo
      if (ee > 0._R8P) t = min(thi, max(tlo, -dot_product(w0, e) / ee))
      x = a2 + t * d2 ; s = min(1._R8P, max(0._R8P, dot_product(d1, x - a1) / aa))
      dd = sum(((a1 + s * d1) - x)**2)
      if (dd < best) then
         best = dd ; c1 = a1 + s * d1 ; c2 = x
      endif
   endif
   ! the end points of segment 1 against segment 2
   call point_segment(a1, a2, d2, x) ; dd = sum((a1 - x)**2)
   if (dd < best) then
      best = dd ; c1 = a1 ; c2 = x
   endif
   call point_segment(b1, a2, d2, x) ; dd = sum((b1 - x)**2)
   if (dd < best) then
      best = dd ; c1 = b1 ; c2 = x
   endif
   if (swap) then
      x = c1 ; c1 = c2 ; c2 = x
   endif
   contains
      pure function proj(v) result(w)
      !< Component of v normal to d1.
      real(R8P), intent(in) :: v(3)
      real(R8P)             :: w(3)

      w = v - (dot_product(v, d1) / aa) * d1
      endfunction proj

      pure subroutine point_segment(q, o, d, f)
      !< Foot f of q on the segment o + [0, 1] d.
      real(R8P), intent(in)  :: q(3), o(3), d(3)
      real(R8P), intent(out) :: f(3)
      real(R8P)              :: dd2, u

      dd2 = dot_product(d, d) ; u = 0._R8P
      if (dd2 > 0._R8P) u = min(1._R8P, max(0._R8P, dot_product(q - o, d) / dd2))
      f = o + u * d
      endsubroutine point_segment
   endsubroutine closest_segment_segment

   subroutine body_distance(fa, ta, fb, tb, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
   !< Minimum squared distance between the facets fa (tree ta) and fb (tree tb), the facets realizing it and the
   !< witness points; `meet` if they touch or intersect. Best-first descent of the two trees (see the module).
   type(facet_object),     intent(in)         :: fa(:), fb(:) !< Facets.
   type(aabb_tree_object), intent(in), target :: ta, tb       !< Their trees.
   real(R8P),              intent(out)        :: d2           !< Squared distance.
   integer(I4P),           intent(out)        :: ia, ib       !< Closest facets.
   real(R8P),              intent(out)        :: pa(3), pb(3) !< Witness points.
   logical,                intent(out)        :: meet         !< Contact.
   type(rigid_pose_t),     intent(in), optional :: pose_a, pose_b !< Poses of the two bodies (default: as they are).
   real(R8P),    allocatable                  :: hb(:)        !< Heap: bounds.
   integer(I4P), allocatable                  :: hn(:, :)     !< Heap: node pairs.
   integer(I4P)                               :: hs, na, nb, ca(MAX_CHILDREN), cb(MAX_CHILDREN), nca, ncb, i, j
   real(R8P)                                  :: bound
   type(aabb_node_object), pointer            :: nda, ndb
   type(posed_facets_t)                       :: xa, xb       !< Posed facets of the two bodies (filled on use).

   d2 = huge(1._R8P) ; ia = 0_I4P ; ib = 0_I4P ; pa = 0._R8P ; pb = 0._R8P ; meet = .false.
   if (size(fa) == 0 .or. size(fb) == 0) return
   if (ta%get_nodes_number() <= 0 .or. tb%get_nodes_number() <= 0) then
      call body_distance_brute(fa, fb, d2, ia, ib, pa, pb, meet, pose_a, pose_b) ; return
   endif
   call posed_allocate(xa, size(fa, kind=I4P)) ; call posed_allocate(xb, size(fb, kind=I4P))
   allocate(hb(64), hn(2, 64)) ; hs = 0_I4P
   call push(node_bound(0_I4P, 0_I4P), 0_I4P, 0_I4P)
   do while (hs > 0_I4P)
      call pop(bound, na, nb)
      if (bound > d2 * (1._R8P + PRUNE_REL)) exit ! every pair left is farther
      nda => ta%node_at(i=na) ; ndb => tb%node_at(i=nb)
      call ta%enumerate_children(n=na, out_idx=ca, nchild=nca)
      call tb%enumerate_children(n=nb, out_idx=cb, nchild=ncb)
      if (nca == 0 .and. ncb == 0) then
         call leaf_pair(nda, ndb, fa, fb, xa, xb, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
      elseif (ncb == 0 .or. (nca > 0 .and. volume(nda) >= volume(ndb))) then
         do i = 1, nca
            if (.not. live(ta, ca(i))) cycle
            call push(node_bound(ca(i), nb), ca(i), nb)
         enddo
      else
         do j = 1, ncb
            if (.not. live(tb, cb(j))) cycle
            call push(node_bound(na, cb(j)), na, cb(j))
         enddo
      endif
   enddo
   contains
      function live(t, n) result(ok)
      type(aabb_tree_object), intent(in), target :: t
      integer(I4P),           intent(in)         :: n
      logical                                    :: ok
      type(aabb_node_object), pointer            :: nd

      nd => t%node_at(i=n)
      ok = associated(nd)
      if (ok) ok = nd%is_allocated()
      endfunction live

      function node_bound(m, n) result(bb)
      integer(I4P), intent(in)        :: m, n
      real(R8P)                       :: bb
      type(aabb_node_object), pointer :: x, y

      type(vector_R8P)                :: l1, h1, l2, h2

      x => ta%node_at(i=m) ; y => tb%node_at(i=n)
      l1 = x%bmin() ; h1 = x%bmax() ; l2 = y%bmin() ; h2 = y%bmax()
      if (present(pose_a)) call moved_box(pose_a, l1, h1)
      if (present(pose_b)) call moved_box(pose_b, l2, h2)
      bb = box_distance2(l1, h1, l2, h2)
      endfunction node_bound

      function volume(x) result(v)
      type(aabb_node_object), intent(in) :: x
      real(R8P)                          :: v
      type(vector_R8P)                   :: lo, hi

      lo = x%bmin() ; hi = x%bmax()
      v = (hi%x - lo%x) * (hi%y - lo%y) * (hi%z - lo%z)
      endfunction volume

      subroutine push(b, m, n)
      real(R8P),    intent(in) :: b
      integer(I4P), intent(in) :: m, n
      real(R8P),    allocatable :: tb_(:)
      integer(I4P), allocatable :: tn(:, :)
      integer(I4P)              :: k, up
      real(R8P)                 :: sb
      integer(I4P)              :: sn(2)

      if (b > d2 * (1._R8P + PRUNE_REL)) return
      if (hs == size(hb)) then
         allocate(tb_(2 * hs), tn(2, 2 * hs)) ; tb_(1:hs) = hb ; tn(:, 1:hs) = hn
         call move_alloc(tb_, hb) ; call move_alloc(tn, hn)
      endif
      hs = hs + 1_I4P ; hb(hs) = b ; hn(:, hs) = [m, n]
      k = hs
      do while (k > 1)
         up = k / 2
         if (.not. heap_less(k, up)) exit
         sb = hb(up) ; hb(up) = hb(k) ; hb(k) = sb
         sn = hn(:, up) ; hn(:, up) = hn(:, k) ; hn(:, k) = sn
         k = up
      enddo
      endsubroutine push

      subroutine pop(b, m, n)
      real(R8P),    intent(out) :: b
      integer(I4P), intent(out) :: m, n
      integer(I4P)              :: k, c, sn(2)
      real(R8P)                 :: sb

      b = hb(1) ; m = hn(1, 1) ; n = hn(2, 1)
      hb(1) = hb(hs) ; hn(:, 1) = hn(:, hs) ; hs = hs - 1_I4P
      k = 1
      do
         c = 2 * k
         if (c > hs) exit
         if (c + 1 <= hs) then
            if (heap_less(c + 1, c)) c = c + 1
         endif
         if (.not. heap_less(c, k)) exit
         sb = hb(c) ; hb(c) = hb(k) ; hb(k) = sb
         sn = hn(:, c) ; hn(:, c) = hn(:, k) ; hn(:, k) = sn
         k = c
      enddo
      endsubroutine pop

      pure function heap_less(x, y) result(l)
      !< Order of heap entries: bound, then node pair (deterministic).
      integer(I4P), intent(in) :: x, y
      logical                  :: l

      if (hb(x) /= hb(y)) then
         l = hb(x) < hb(y)
      elseif (hn(1, x) /= hn(1, y)) then
         l = hn(1, x) < hn(1, y)
      else
         l = hn(2, x) < hn(2, y)
      endif
      endfunction heap_less
   endsubroutine body_distance

   subroutine body_safe_step(fa, ta, fb, tb, pose_a, pose_b, ka, tha, dta, kb, thb, dtb, tolerance, cap, step, pairs)
   !< The directional safe step of two moving bodies at the poses pose_a, pose_b (see the module notes), at most cap.
   type(facet_object),     intent(in)            :: fa(:), fb(:)   !< Facets (body frames).
   type(aabb_tree_object), intent(in), target    :: ta, tb         !< Their trees (body frames).
   type(rigid_pose_t),     intent(in)            :: pose_a, pose_b !< Current poses.
   real(R8P),              intent(in)            :: ka(3), kb(3)   !< Unit rotation axes over the step.
   real(R8P),              intent(in)            :: tha, thb       !< Rotation angles over the step (>= 0).
   real(R8P),              intent(in)            :: dta(3), dtb(3) !< Translations over the step.
   real(R8P),              intent(in)            :: tolerance      !< Contact tolerance.
   real(R8P),              intent(in)            :: cap            !< Largest step wanted (1 - s).
   real(R8P),              intent(out)           :: step           !< Safe step (0 if a pair is within tol / 2).
   integer(I4P),           intent(out), optional :: pairs          !< Facet pairs evaluated exactly.
   integer(I4P), allocatable                     :: stk(:, :), tmp(:, :) !< Node-pair stack.
   integer(I4P)                                  :: ns, na, nb, ca(MAX_CHILDREN), cb(MAX_CHILDREN), nca, ncb, i, j, k, l
   integer(I4P)                                  :: np             !< Pairs evaluated.
   type(aabb_node_object), pointer               :: nda, ndb       !< Nodes.
   type(list_id_object)                          :: la, lb         !< Facet lists.
   real(R8P)                                     :: half, bound    !< tol / 2, node bound.

   step = cap ; np = 0_I4P ; half = 0.5_R8P * tolerance
   if (size(fa) > 0 .and. size(fb) > 0) then
      if (ta%get_nodes_number() <= 0 .or. tb%get_nodes_number() <= 0) then ! no trees: every pair
         do i = 1, size(fa, kind=I4P)
            do j = 1, size(fb, kind=I4P)
               call pair_safe_step(fa, fb, i, j, pose_a, pose_b, ka, tha, dta, kb, thb, dtb, half, step, np)
            enddo
         enddo
      else
         allocate(stk(2, 64)) ; ns = 1_I4P ; stk(:, 1) = [0_I4P, 0_I4P]
         do while (ns > 0_I4P .and. step > 0._R8P)
            na = stk(1, ns) ; nb = stk(2, ns) ; ns = ns - 1_I4P
            nda => ta%node_at(i=na) ; ndb => tb%node_at(i=nb)
            bound = node_plain_step(nda, ndb, pose_a, pose_b, ka, tha, dta, kb, thb, dtb, half)
            if (bound >= step) cycle
            call ta%enumerate_children(n=na, out_idx=ca, nchild=nca)
            call tb%enumerate_children(n=nb, out_idx=cb, nchild=ncb)
            if (nca == 0 .and. ncb == 0) then
               if (.not. nda%has_facets() .or. .not. ndb%has_facets()) cycle
               la = nda%facet_id() ; lb = ndb%facet_id()
               do k = 1, la%ids_number
                  do l = 1, lb%ids_number
                     call pair_safe_step(fa, fb, la%id(k), lb%id(l), pose_a, pose_b, ka, tha, dta, kb, thb, dtb, half, &
                                         step, np)
                  enddo
               enddo
            else
               if (ns + max(nca, ncb) > size(stk, 2)) then
                  allocate(tmp(2, 2 * size(stk, 2) + max(nca, ncb))) ; tmp(:, 1:ns) = stk(:, 1:ns)
                  call move_alloc(tmp, stk)
               endif
               if (ncb == 0 .or. (nca > 0 .and. node_box_volume(nda) >= node_box_volume(ndb))) then
                  do i = 1, nca
                     if (.not. node_live(ta, ca(i))) cycle
                     ns = ns + 1_I4P ; stk(:, ns) = [ca(i), nb]
                  enddo
               else
                  do j = 1, ncb
                     if (.not. node_live(tb, cb(j))) cycle
                     ns = ns + 1_I4P ; stk(:, ns) = [na, cb(j)]
                  enddo
               endif
            endif
         enddo
      endif
   endif
   if (present(pairs)) pairs = np
   endsubroutine body_safe_step

   subroutine pair_safe_step(fa, fb, i, j, pose_a, pose_b, ka, tha, dta, kb, thb, dtb, half, step, np)
   !< Lower step to facet pair (i, j)'s directional step if smaller (plain pre-filter first).
   type(facet_object), intent(in)    :: fa(:), fb(:)   !< Facets.
   integer(I4P),       intent(in)    :: i, j           !< Pair.
   type(rigid_pose_t), intent(in)    :: pose_a, pose_b !< Current poses.
   real(R8P),          intent(in)    :: ka(3), kb(3), tha, thb, dta(3), dtb(3) !< Motion over the step.
   real(R8P),          intent(in)    :: half           !< tol / 2.
   real(R8P),          intent(inout) :: step           !< Best step.
   integer(I4P),       intent(inout) :: np             !< Pairs evaluated.
   real(R8P)                         :: a(3, 3), b(3, 3), ra, rb, dd, qa(3), qb(3), d, n(3), c, g2
   logical                           :: mt
   integer(I4P)                      :: k

   if (i < 1 .or. i > size(fa) .or. j < 1 .or. j > size(fb)) return
   ra = 0._R8P ; rb = 0._R8P
   do k = 1, 3
      a(:, k) = [fa(i)%vertex(k)%x, fa(i)%vertex(k)%y, fa(i)%vertex(k)%z]
      b(:, k) = [fb(j)%vertex(k)%x, fb(j)%vertex(k)%y, fb(j)%vertex(k)%z]
      ra = max(ra, axis_distance(ka, rotate(pose_a, a(:, k)))) ; rb = max(rb, axis_distance(kb, rotate(pose_b, b(:, k))))
      a(:, k) = to_world(pose_a, a(:, k)) ; b(:, k) = to_world(pose_b, b(:, k))
   enddo
   ! plain pre-filter: box distance over the pair's plain speed
   g2 = sqrt(max(0._R8P, minval(b(1, :)) - maxval(a(1, :)), minval(a(1, :)) - maxval(b(1, :)))**2 + &
             max(0._R8P, minval(b(2, :)) - maxval(a(2, :)), minval(a(2, :)) - maxval(b(2, :)))**2 + &
             max(0._R8P, minval(b(3, :)) - maxval(a(3, :)), minval(a(3, :)) - maxval(b(3, :)))**2)
   c = norm2(dta) + norm2(dtb) + tha * ra + thb * rb
   if (c > 0._R8P) then
      if ((g2 - half) / c >= step) return
   else
      return ! the pair does not move relative to each other: its distance stays d_ij > tol / 2 (d_k > tol)
   endif
   np = np + 1_I4P
   call tri_tri_distance(a, b, dd, qa, qb, mt)
   d = sqrt(dd)
   if (mt .or. d <= half) then
      step = 0._R8P ; return
   endif
   n = (qb - qa) / d
   c = dot_product(n, dta - dtb) + tha * ra + thb * rb
   if (c > 0._R8P) step = min(step, (d - half) / c)
   endsubroutine pair_safe_step

   function node_plain_step(x, y, pose_a, pose_b, ka, tha, dta, kb, thb, dtb, half) result(st)
   !< (box distance - tol / 2) / (mu_x + mu_y) of two nodes at the current poses (huge if they do not move).
   type(aabb_node_object), intent(in) :: x, y           !< Nodes.
   type(rigid_pose_t),     intent(in) :: pose_a, pose_b !< Current poses.
   real(R8P),              intent(in) :: ka(3), kb(3), tha, thb, dta(3), dtb(3) !< Motion over the step.
   real(R8P),              intent(in) :: half           !< tol / 2.
   real(R8P)                          :: st             !< Plain step bound.
   type(vector_R8P)                   :: l1, h1, l2, h2 !< Boxes.
   real(R8P)                          :: mu             !< Speed bound.

   l1 = x%bmin() ; h1 = x%bmax() ; l2 = y%bmin() ; h2 = y%bmax()
   mu = norm2(dta) + norm2(dtb) + tha * box_rho(ka, pose_a, l1, h1) + thb * box_rho(kb, pose_b, l2, h2)
   call moved_box(pose_a, l1, h1) ; call moved_box(pose_b, l2, h2)
   if (mu > 0._R8P) then
      st = (sqrt(box_distance2(l1, h1, l2, h2)) - half) / mu
   else
      st = huge(1._R8P)
   endif
   endfunction node_plain_step

   pure function box_rho(k, pose, lo, hi) result(r)
   !< Largest distance to the axis k (through the body origin) of the 8 corners of a body-frame box, rotated by the pose.
   real(R8P),          intent(in) :: k(3)   !< Unit axis.
   type(rigid_pose_t), intent(in) :: pose   !< Pose.
   type(vector_R8P),   intent(in) :: lo, hi !< Box.
   real(R8P)                      :: r      !< Distance.
   integer(I4P)                   :: c      !< Corner.

   r = 0._R8P
   do c = 0, 7
      r = max(r, axis_distance(k, rotate(pose, [merge(hi%x, lo%x, btest(c, 0)), merge(hi%y, lo%y, btest(c, 1)), &
                                                merge(hi%z, lo%z, btest(c, 2))])))
   enddo
   endfunction box_rho

   pure function axis_distance(k, y) result(r)
   !< |k x y|: distance of y to the axis k through the origin.
   real(R8P), intent(in) :: k(3), y(3) !< Unit axis, point.
   real(R8P)             :: r          !< Distance.

   r = norm2([k(2) * y(3) - k(3) * y(2), k(3) * y(1) - k(1) * y(3), k(1) * y(2) - k(2) * y(1)])
   endfunction axis_distance

   function node_live(t, n) result(ok)
   !< Node n of t exists and is allocated.
   type(aabb_tree_object), intent(in), target :: t  !< Tree.
   integer(I4P),           intent(in)         :: n  !< Node.
   logical                                    :: ok !< Live.
   type(aabb_node_object), pointer            :: nd !< Node.

   nd => t%node_at(i=n)
   ok = associated(nd)
   if (ok) ok = nd%is_allocated()
   endfunction node_live

   pure function node_box_volume(x) result(v)
   !< Volume of a node's box.
   type(aabb_node_object), intent(in) :: x  !< Node.
   real(R8P)                          :: v  !< Volume.
   type(vector_R8P)                   :: lo, hi

   lo = x%bmin() ; hi = x%bmax()
   v = (hi%x - lo%x) * (hi%y - lo%y) * (hi%z - lo%z)
   endfunction node_box_volume

   subroutine body_distance_brute(fa, fb, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
   !< The reference: every facet pair, the same pair function and tie rule as `body_distance`.
   type(facet_object), intent(in)  :: fa(:), fb(:)
   real(R8P),          intent(out) :: d2
   integer(I4P),       intent(out) :: ia, ib
   real(R8P),          intent(out) :: pa(3), pb(3)
   logical,            intent(out) :: meet
   type(rigid_pose_t), intent(in), optional :: pose_a, pose_b
   integer(I4P)                    :: i, j

   d2 = huge(1._R8P) ; ia = 0_I4P ; ib = 0_I4P ; pa = 0._R8P ; pb = 0._R8P ; meet = .false.
   do i = 1, size(fa)
      do j = 1, size(fb)
         call try_pair(fa, fb, i, j, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
      enddo
   enddo
   endsubroutine body_distance_brute

   subroutine leaf_pair(x, y, fa, fb, xa, xb, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
   !< Every facet pair of two leaves in order (row k of x, then column l of y), as `try_pair` would evaluate them, on the
   !< facets posed once per call; a row is skipped when its facet box is farther than the pruning bound from the union of
   !< y's facet boxes (every pair of the row would be pruned, see the module notes).
   type(aabb_node_object), intent(in)    :: x, y           !< Leaves of the trees of fa and fb.
   type(facet_object),     intent(in)    :: fa(:), fb(:)   !< Facets.
   type(posed_facets_t),   intent(inout) :: xa, xb         !< Posed facets of fa and fb.
   real(R8P),              intent(inout) :: d2             !< Best squared distance.
   integer(I4P),           intent(inout) :: ia, ib         !< Best pair.
   real(R8P),              intent(inout) :: pa(3), pb(3)   !< Its witness points.
   logical,                intent(inout) :: meet           !< Its contact.
   type(rigid_pose_t),     intent(in), optional :: pose_a, pose_b !< Poses.
   type(list_id_object)                  :: la, lb         !< Facet lists of the leaves.
   real(R8P)                             :: ub(6)          !< Union box of y's facets.
   logical                               :: finite_union   !< Every facet of y in range is finite.
   integer(I4P)                          :: k, l, i, j

   if (.not. x%has_facets() .or. .not. y%has_facets()) return
   la = x%facet_id() ; lb = y%facet_id()
   do k = 1, la%ids_number
      call pose_facet(xa, fa, la%id(k), pose_a)
   enddo
   ub(1:3) = huge(1._R8P) ; ub(4:6) = -huge(1._R8P) ; finite_union = .true.
   do l = 1, lb%ids_number
      j = lb%id(l)
      if (j < 1 .or. j > size(fb)) cycle
      call pose_facet(xb, fb, j, pose_b)
      if (xb%state(j) == POSED_NONFINITE) then
         finite_union = .false.
      else
         ub(1:3) = min(ub(1:3), xb%box(1:3, j)) ; ub(4:6) = max(ub(4:6), xb%box(4:6, j))
      endif
   enddo
   do k = 1, la%ids_number
      i = la%id(k)
      if (i < 1 .or. i > size(fa)) cycle ! try_pair returns at once for every pair of the row
      if (finite_union .and. xa%state(i) /= POSED_NONFINITE) then
         if (box_gap2(xa%box(:, i), ub) > d2 * (1._R8P + PRUNE_REL)) cycle ! every pair of the row is pruned
      endif
      do l = 1, lb%ids_number
         call try_pair_posed(xa, xb, i, lb%id(l), d2, ia, ib, pa, pb, meet)
      enddo
   enddo
   endsubroutine leaf_pair

   subroutine try_pair_posed(xa, xb, i, j, d2, ia, ib, pa, pb, meet)
   !< `try_pair` on facets posed once per call (same pruning test, same pair function, same tie rule); the exact
   !< classification is skipped for strictly disjoint boxes of two finite, non-degenerate triangles.
   type(posed_facets_t), intent(inout) :: xa, xb         !< Posed facets (i and j already posed when in range).
   integer(I4P),         intent(in)    :: i, j           !< Pair.
   real(R8P),            intent(inout) :: d2             !< Best squared distance.
   integer(I4P),         intent(inout) :: ia, ib         !< Best pair.
   real(R8P),            intent(inout) :: pa(3), pb(3)   !< Its witness points.
   logical,              intent(inout) :: meet           !< Its contact.
   real(R8P)                           :: dd, qa(3), qb(3)
   logical                             :: mt, apart

   if (i < 1 .or. i > size(xa%state) .or. j < 1 .or. j > size(xb%state)) return
   if (box_gap2(xa%box(:, i), xb%box(:, j)) > d2 * (1._R8P + PRUNE_REL)) return
   apart = xb%box(1, j) > xa%box(4, i) .or. xa%box(1, i) > xb%box(4, j) .or. &
           xb%box(2, j) > xa%box(5, i) .or. xa%box(2, i) > xb%box(5, j) .or. &
           xb%box(3, j) > xa%box(6, i) .or. xa%box(3, i) > xb%box(6, j)
   if (apart) apart = facet_regular(xa, i)
   if (apart) apart = facet_regular(xb, j)
   call tri_tri_distance_core(xa%v(:, :, i), xb%v(:, :, j), apart, dd, qa, qb, mt)
   if (dd < d2 .or. (dd == d2 .and. (i < ia .or. (i == ia .and. j < ib)))) then
      d2 = dd ; ia = i ; ib = j ; pa = qa ; pb = qb ; meet = mt
   endif
   endsubroutine try_pair_posed

   subroutine posed_allocate(x, n)
   !< Room for n posed facets, none posed yet.
   type(posed_facets_t), intent(inout) :: x !< Posed facets.
   integer(I4P),         intent(in)    :: n !< Facets.

   allocate(x%v(3, 3, n), x%box(6, n), x%state(n))
   x%state = POSED_NONE
   endsubroutine posed_allocate

   subroutine pose_facet(x, f, i, pose)
   !< Pose facet i of f (if in range and not posed yet): vertices and box as `try_pair` computes them.
   type(posed_facets_t), intent(inout)        :: x    !< Posed facets.
   type(facet_object),   intent(in)           :: f(:) !< Facets.
   integer(I4P),         intent(in)           :: i    !< Facet.
   type(rigid_pose_t),   intent(in), optional :: pose !< Pose.
   real(R8P)                                  :: a(3, 3)
   integer(I4P)                               :: k

   if (i < 1 .or. i > size(f)) return
   if (x%state(i) /= POSED_NONE) return
   do k = 1, 3
      a(:, k) = [f(i)%vertex(k)%x, f(i)%vertex(k)%y, f(i)%vertex(k)%z]
      if (present(pose)) a(:, k) = to_world(pose, a(:, k))
   enddo
   x%v(:, :, i) = a
   x%box(:, i) = [minval(a(1, :)), minval(a(2, :)), minval(a(3, :)), maxval(a(1, :)), maxval(a(2, :)), maxval(a(3, :))]
   if (all(abs(a) <= huge(1._R8P))) then
      x%state(i) = POSED_FINITE
   else
      x%state(i) = POSED_NONFINITE
   endif
   endsubroutine pose_facet

   function facet_regular(x, i) result(ok)
   !< Posed facet i is regular: every coordinate 0 or of magnitude in [REGULAR_MIN, REGULAR_MAX] (the exact predicates
   !< neither overflow nor underflow) and the triangle exactly non-degenerate (a nonzero exact orient2d in some coordinate
   !< projection, the test of `facet_drop`). Decided once per call. `pair_cut` is exact on two regular triangles.
   type(posed_facets_t), intent(inout) :: x  !< Posed facets.
   integer(I4P),         intent(in)    :: i  !< Facet (posed).
   logical                             :: ok !< Regular.
   integer(I4P)                        :: d, u, w

   if (x%state(i) == POSED_FINITE) then
      x%state(i) = POSED_IRREGULAR
      if (all(x%v(:, :, i) == 0._R8P .or. (abs(x%v(:, :, i)) >= REGULAR_MIN .and. abs(x%v(:, :, i)) <= REGULAR_MAX))) then
         do d = 1, 3
            u = mod(d, 3_I4P) + 1_I4P ; w = mod(d + 1_I4P, 3_I4P) + 1_I4P
            if (orient2d_xy(x%v(u, 1, i), x%v(w, 1, i), x%v(u, 2, i), x%v(w, 2, i), x%v(u, 3, i), x%v(w, 3, i)) /= &
                0._R8P) then
               x%state(i) = POSED_REGULAR ; exit
            endif
         enddo
      endif
   endif
   ok = x%state(i) == POSED_REGULAR
   endfunction facet_regular

   pure function box_gap2(p, q) result(d2)
   !< `box_distance2` of the boxes p and q stored as [min x, y, z, max x, y, z]: the same operations in the same order.
   real(R8P), intent(in) :: p(6), q(6) !< Boxes.
   real(R8P)             :: d2         !< Squared distance (0 if they overlap).

   d2 = max(0._R8P, q(1) - p(4), p(1) - q(4))**2 + max(0._R8P, q(2) - p(5), p(2) - q(5))**2 + &
        max(0._R8P, q(3) - p(6), p(3) - q(6))**2
   endfunction box_gap2

   subroutine try_pair(fa, fb, i, j, d2, ia, ib, pa, pb, meet, pose_a, pose_b)
   !< Evaluate facets i of fa and j of fb; keep them if closer, or as close with a smaller (i, j).
   type(facet_object), intent(in)    :: fa(:), fb(:)
   integer(I4P),       intent(in)    :: i, j
   real(R8P),          intent(inout) :: d2
   integer(I4P),       intent(inout) :: ia, ib
   real(R8P),          intent(inout) :: pa(3), pb(3)
   logical,            intent(inout) :: meet
   type(rigid_pose_t), intent(in), optional :: pose_a, pose_b
   real(R8P)                         :: a(3, 3), b(3, 3), dd, qa(3), qb(3)
   type(vector_R8P)                  :: l1, h1, l2, h2
   logical                           :: mt
   integer(I4P)                      :: k

   if (i < 1 .or. i > size(fa) .or. j < 1 .or. j > size(fb)) return
   do k = 1, 3
      a(:, k) = [fa(i)%vertex(k)%x, fa(i)%vertex(k)%y, fa(i)%vertex(k)%z]
      b(:, k) = [fb(j)%vertex(k)%x, fb(j)%vertex(k)%y, fb(j)%vertex(k)%z]
      if (present(pose_a)) a(:, k) = to_world(pose_a, a(:, k))
      if (present(pose_b)) b(:, k) = to_world(pose_b, b(:, k))
   enddo
   l1 = vector_R8P(minval(a(1, :)), minval(a(2, :)), minval(a(3, :))) ; h1 = vector_R8P(maxval(a(1, :)), maxval(a(2, :)), &
                                                                                        maxval(a(3, :)))
   l2 = vector_R8P(minval(b(1, :)), minval(b(2, :)), minval(b(3, :))) ; h2 = vector_R8P(maxval(b(1, :)), maxval(b(2, :)), &
                                                                                        maxval(b(3, :)))
   if (box_distance2(l1, h1, l2, h2) > d2 * (1._R8P + PRUNE_REL)) return
   call tri_tri_distance(a, b, dd, qa, qb, mt)
   if (dd < d2 .or. (dd == d2 .and. (i < ia .or. (i == ia .and. j < ib)))) then
      d2 = dd ; ia = i ; ib = j ; pa = qa ; pb = qb ; meet = mt
   endif
   endsubroutine try_pair

   subroutine moved_box(pose, lo, hi)
   !< The box [lo, hi] replaced by the box of its 8 corners moved by the pose (a bound of the moved content).
   type(rigid_pose_t), intent(in)    :: pose
   type(vector_R8P),   intent(inout) :: lo, hi
   real(R8P)                         :: c(3), w(3), mn(3), mx(3)
   integer(I4P)                      :: k

   mn = huge(1._R8P) ; mx = -huge(1._R8P)
   do k = 0, 7
      c = [merge(hi%x, lo%x, btest(k, 0)), merge(hi%y, lo%y, btest(k, 1)), merge(hi%z, lo%z, btest(k, 2))]
      w = to_world(pose, c)
      mn = min(mn, w) ; mx = max(mx, w)
   enddo
   ! a corner rounded inward by an ulp must not shrink the box below the moved facets (also rounded): widen by an ulp
   mn = mn - spacing(abs(mn)) * 4._R8P ; mx = mx + spacing(abs(mx)) * 4._R8P
   lo = vector_R8P(mn(1), mn(2), mn(3)) ; hi = vector_R8P(mx(1), mx(2), mx(3))
   endsubroutine moved_box

   pure function box_distance2(lo1, hi1, lo2, hi2) result(d2)
   !< Squared distance between two axis-aligned boxes (0 if they overlap).
   type(vector_R8P), intent(in) :: lo1, hi1, lo2, hi2
   real(R8P)                    :: d2

   d2 = gap(lo1%x, hi1%x, lo2%x, hi2%x)**2 + gap(lo1%y, hi1%y, lo2%y, hi2%y)**2 + gap(lo1%z, hi1%z, lo2%z, hi2%z)**2
   contains
      pure function gap(a0, a1, b0, b1) result(g)
      real(R8P), intent(in) :: a0, a1, b0, b1
      real(R8P)             :: g

      g = max(0._R8P, b0 - a1, a0 - b1)
      endfunction gap
   endfunction box_distance2
endmodule fossil_body_distance
