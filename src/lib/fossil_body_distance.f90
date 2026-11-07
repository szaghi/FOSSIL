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
use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,     only : facet_object
use fossil_implicit_points,  only : implicit_point_t, ip_approx
use fossil_list_id_object,   only : list_id_object
use fossil_predicates,       only : orient2d_xy
use fossil_rigid_pose,       only : rigid_pose_t, to_world
use fossil_triangle_closest, only : triangle_closest_st
use fossil_tritri_exact,     only : pair_cut, facet_drop, PAIR_SEPARATED, PAIR_COPLANAR
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: tri_tri_distance, body_distance, body_distance_brute, closest_on_triangle, closest_segment_segment

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P   !< Children buffer (octree fan-out).
real(R8P),    parameter :: PRUNE_REL = 1.e-10_R8P !< Relative margin of the pruning test (bounds are not compared
                                                  !< bit for bit with computed distances).

contains
   subroutine tri_tri_distance(a, b, d2, pa, pb, meet)
   !< Squared distance between triangles a and b (vertices as columns), witness points pa on a and pb on b.
   real(R8P), intent(in)  :: a(3, 3), b(3, 3) !< Triangles.
   real(R8P), intent(out) :: d2               !< Squared distance.
   real(R8P), intent(out) :: pa(3), pb(3)     !< Witness points.
   logical,   intent(out) :: meet             !< The triangles meet (exact).
   type(vector_R8P)       :: ta(3), tb(3)
   type(implicit_point_t) :: p, q
   real(R8P)              :: c(3), c2(3), dd
   integer(I4P)           :: i, j, kind

   do i = 1, 3
      ta(i) = vector_R8P(a(1, i), a(2, i), a(3, i)) ; tb(i) = vector_R8P(b(1, i), b(2, i), b(3, i))
   enddo
   call pair_cut(ta, tb, kind, p, q)
   meet = kind /= PAIR_SEPARATED
   if (kind == PAIR_COPLANAR) meet = coplanar_meet(a, b, facet_drop(ta)) ! same plane: do they overlap?
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
   endsubroutine tri_tri_distance

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

   d2 = huge(1._R8P) ; ia = 0_I4P ; ib = 0_I4P ; pa = 0._R8P ; pb = 0._R8P ; meet = .false.
   if (size(fa) == 0 .or. size(fb) == 0) return
   if (ta%get_nodes_number() <= 0 .or. tb%get_nodes_number() <= 0) then
      call body_distance_brute(fa, fb, d2, ia, ib, pa, pb, meet, pose_a, pose_b) ; return
   endif
   allocate(hb(64), hn(2, 64)) ; hs = 0_I4P
   call push(node_bound(0_I4P, 0_I4P), 0_I4P, 0_I4P)
   do while (hs > 0_I4P)
      call pop(bound, na, nb)
      if (bound > d2 * (1._R8P + PRUNE_REL)) exit ! every pair left is farther
      nda => ta%node_at(i=na) ; ndb => tb%node_at(i=nb)
      call ta%enumerate_children(n=na, out_idx=ca, nchild=nca)
      call tb%enumerate_children(n=nb, out_idx=cb, nchild=ncb)
      if (nca == 0 .and. ncb == 0) then
         call leaf_pair(nda, ndb)
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

      subroutine leaf_pair(x, y)
      type(aabb_node_object), intent(in) :: x, y
      type(list_id_object)               :: la, lb
      integer(I4P)                       :: k, l

      if (.not. x%has_facets() .or. .not. y%has_facets()) return
      la = x%facet_id() ; lb = y%facet_id()
      do k = 1, la%ids_number
         do l = 1, lb%ids_number
            call try_pair(fa, fb, la%id(k), lb%id(l), d2, ia, ib, pa, pb, meet, pose_a, pose_b)
         enddo
      enddo
      endsubroutine leaf_pair

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
