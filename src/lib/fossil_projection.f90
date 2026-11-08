!< FOSSIL, constrained projection onto bodies for body-fitted grids (issue #74).
module fossil_projection
!< FOSSIL, constrained projection onto bodies for body-fitted grids (issue #74).
!<
!< Design: docs/guide/advanced/body-fitted.md. A `projection_target_t` is built once from a mesh (vertices, facets and
!< their patches; 2-D: the segments of a curve and their patches) and every query reads it:
!< - the elements (triangles; 2-D segments) sorted by patch, one bounding-volume tree per patch;
!< - 3-D: the constraint curves, i.e. the feature edges of #61 (sharp, boundary, non-manifold) and the patch boundaries
!<   (edges between facets of different patches), chained into polylines between corners, with one segment tree;
!< - the corners: the vertices where constraint edges meet in a number other than two, or turn by more than
!<   corner_angle (3-D); the feature corners and the points where two patches meet (2-D).
!< Trees: median split on the element centroids along the longest box axis (deterministic: ties by element index),
!< leaves of at most LEAF elements; node(:) = [left, right, first, count] over a permutation of the elements.
!< This module knows no surface or curve type: the surface and the curve build their target through it.
use penf,                 only : I4P, I8P, R8P
use fossil_feature_edges, only : find_feature_edges

implicit none
private
public :: projection_target_t, build_target_3d, build_target_2d
public :: PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
public :: bvh_build

integer(I4P), parameter :: PROJ_FREE   = 0_I4P !< The whole body.
integer(I4P), parameter :: PROJ_PATCH  = 1_I4P !< One patch (id).
integer(I4P), parameter :: PROJ_CURVE  = 2_I4P !< One constraint curve (id > 0), or any (id 0); 3-D only.
integer(I4P), parameter :: PROJ_CORNER = 3_I4P !< One corner (id > 0), or the nearest (id 0); fixed.
integer(I4P), parameter :: LEAF        = 4_I4P !< Elements per leaf.

type :: projection_target_t
   !< What the constrained projection reads.
   integer(I4P)              :: dims = 0_I4P        !< 3 (surface) or 2 (curve); 0: not built.
   integer(I4P)              :: ne = 0_I4P          !< Elements (triangles or segments).
   real(R8P),    allocatable :: el(:,:)             !< Elements sorted by patch [9, ne] (2-D: [u0, v0, u1, v1, ...]).
   integer(I4P), allocatable :: el_id(:)            !< Original facet / segment index [ne].
   integer(I4P), allocatable :: el_patch(:)         !< Patch of each element [ne].
   integer(I4P)              :: np = 0_I4P          !< Patches (ids 1..np).
   integer(I4P), allocatable :: patch_root(:)       !< Root node of each patch's tree [np] (0: empty patch).
   real(R8P),    allocatable :: node_box(:,:)       !< Element-tree node boxes [6, nn].
   integer(I4P), allocatable :: node(:,:)           !< Element-tree nodes [4, nn]: left, right, first, count.
   integer(I4P), allocatable :: perm(:)             !< Element order of the trees [ne].
   integer(I4P)              :: nc = 0_I4P          !< Constraint curves (3-D).
   integer(I4P), allocatable :: curve_first(:)      !< Segments of curve c: curve_first(c):curve_first(c+1)-1 [nc+1].
   real(R8P),    allocatable :: seg(:,:)            !< Constraint segments sorted by curve [6, ns].
   integer(I4P), allocatable :: seg_curve(:)        !< Curve of each segment [ns].
   integer(I4P)              :: seg_root = 0_I4P    !< Root of the segment tree (0: none).
   real(R8P),    allocatable :: seg_box(:,:)        !< Segment-tree node boxes [6, nsn].
   integer(I4P), allocatable :: seg_node(:,:)       !< Segment-tree nodes [4, nsn].
   integer(I4P), allocatable :: seg_perm(:)         !< Segment order of the tree [ns].
   integer(I4P)              :: nk = 0_I4P          !< Corners.
   real(R8P),    allocatable :: corner(:,:)         !< Corner coordinates [3, nk] (2-D: [u, v, 0]).
endtype projection_target_t

contains
   subroutine build_target_3d(v, fv, fpatch, angle, corner_angle, target)
   !< The target of a triangle mesh: vertices v [3, nv], facets fv [3, nf] (vertex indices), patch of each facet.
   real(R8P),                 intent(in)  :: v(:,:)       !< Vertices [3, nv].
   integer(I4P),              intent(in)  :: fv(:,:)      !< Facets [3, nf].
   integer(I4P),              intent(in)  :: fpatch(:)    !< Patch of each facet [nf] (>= 1).
   real(R8P),                 intent(in)  :: angle        !< Sharp-edge threshold [degrees].
   real(R8P),                 intent(in)  :: corner_angle !< Corner turning threshold [degrees].
   type(projection_target_t), intent(out) :: target       !< Target.
   integer(I4P), allocatable              :: fe(:,:), fk(:), fc(:), ps(:), pv(:), pe(:,:), ed(:,:), order(:)
   integer(I8P), allocatable              :: key(:)
   real(R8P),    allocatable              :: boxes(:,:)
   integer(I4P)                           :: nf, i, k, n, ndeg

   nf = size(fv, 2, kind=I4P)
   target%dims = 3_I4P ; target%ne = nf
   ! elements sorted by patch (stable)
   target%np = maxval([0_I4P, fpatch])
   allocate(key(nf), order(nf))
   do i = 1, nf
      key(i) = int(fpatch(i), I8P) * (int(nf, I8P) + 1_I8P) + i
   enddo
   call sort_keys(key, order)
   allocate(target%el(9, nf), target%el_id(nf), target%el_patch(nf))
   do k = 1, nf
      i = order(k)
      target%el(:, k) = [v(:, fv(1, i)), v(:, fv(2, i)), v(:, fv(3, i))]
      target%el_id(k) = i ; target%el_patch(k) = fpatch(i)
   enddo
   allocate(boxes(6, nf))
   do k = 1, nf
      boxes(1:3, k) = min(target%el(1:3, k), target%el(4:6, k), target%el(7:9, k))
      boxes(4:6, k) = max(target%el(1:3, k), target%el(4:6, k), target%el(7:9, k))
   enddo
   call build_patch_trees(target, boxes)
   ! constraint edges: features of #61 and patch boundaries
   call find_feature_edges(fv, v, angle, corner_angle, .true., .true., fe, fk, fc, ps, pv, ndeg)
   call patch_boundary_edges(fv, fpatch, pe)
   n = size(fe, 2, kind=I4P) + size(pe, 2, kind=I4P)
   allocate(ed(2, n)) ; ed(:, 1:size(fe, 2)) = fe ; ed(:, size(fe, 2) + 1:) = pe
   call unique_edges(ed, size(v, 2, kind=I4P))
   call chain_curves(ed, v, corner_angle, target)
   endsubroutine build_target_3d

   subroutine build_target_2d(u, v, next, prev, spatch, corner_angle, target)
   !< The target of a curve: vertices (u, v) [nv], segment i = i -> next(i) (0: none), its patch spatch(i) (>= 1).
   real(R8P),                 intent(in)  :: u(:), v(:)    !< Vertices [nv].
   integer(I4P),              intent(in)  :: next(:), prev(:) !< Loop links [nv] (0: none).
   integer(I4P),              intent(in)  :: spatch(:)     !< Patch of each segment [nv].
   real(R8P),                 intent(in)  :: corner_angle  !< Corner turning threshold [degrees].
   type(projection_target_t), intent(out) :: target        !< Target.
   integer(I4P), allocatable              :: order(:)
   integer(I8P), allocatable              :: key(:)
   real(R8P),    allocatable              :: boxes(:,:), ck(:,:)
   real(R8P)                              :: a(2), b(2), cosc
   integer(I4P)                           :: nv, ns, i, k, j

   nv = size(u, kind=I4P)
   target%dims = 2_I4P
   ns = count(next(1:nv) > 0)
   target%ne = ns
   target%np = 0_I4P
   do i = 1, nv
      if (next(i) > 0) target%np = max(target%np, spatch(i))
   enddo
   allocate(key(ns), order(ns)) ; k = 0
   do i = 1, nv
      if (next(i) == 0) cycle
      k = k + 1 ; key(k) = int(spatch(i), I8P) * (int(nv, I8P) + 1_I8P) + i
   enddo
   call sort_keys(key, order)
   allocate(target%el(9, ns), target%el_id(ns), target%el_patch(ns))
   target%el = 0._R8P
   do k = 1, ns
      i = int(mod(key(order(k)), int(nv, I8P) + 1_I8P), I4P)
      target%el(1:4, k) = [u(i), v(i), u(next(i)), v(next(i))]
      target%el_id(k) = i ; target%el_patch(k) = spatch(i)
   enddo
   allocate(boxes(6, ns))
   do k = 1, ns
      boxes(:, k) = [min(target%el(1, k), target%el(3, k)), min(target%el(2, k), target%el(4, k)), 0._R8P, &
                     max(target%el(1, k), target%el(3, k)), max(target%el(2, k), target%el(4, k)), 0._R8P]
   enddo
   call build_patch_trees(target, boxes)
   ! corners: open ends, turning above corner_angle, patch junctions
   cosc = cos(corner_angle * acos(-1._R8P) / 180._R8P)
   allocate(ck(3, nv)) ; j = 0
   do i = 1, nv
      if (next(i) == 0 .and. prev(i) == 0) cycle
      if (next(i) == 0 .or. prev(i) == 0) then
         j = j + 1 ; ck(:, j) = [u(i), v(i), 0._R8P] ; cycle
      endif
      a = [u(i) - u(prev(i)), v(i) - v(prev(i))] ; b = [u(next(i)) - u(i), v(next(i)) - v(i)]
      if (spatch(prev(i)) /= spatch(i) .or. dot_product(a, b) < cosc * norm2(a) * norm2(b)) then
         j = j + 1 ; ck(:, j) = [u(i), v(i), 0._R8P]
      endif
   enddo
   target%nk = j ; target%corner = ck(:, 1:j)
   target%nc = 0_I4P ; allocate(target%curve_first(1), target%seg(6, 0), target%seg_curve(0)) ; target%curve_first = 1
   allocate(target%seg_box(6, 0), target%seg_node(4, 0), target%seg_perm(0))
   endsubroutine build_target_2d

   subroutine build_patch_trees(target, boxes)
   !< One tree per patch over the elements (sorted by patch) with the given boxes [6, ne].
   type(projection_target_t), intent(inout) :: target     !< Target (el_patch set).
   real(R8P),                 intent(in)    :: boxes(:,:) !< Element boxes [6, ne].
   real(R8P),    allocatable                :: nb(:,:)
   integer(I4P), allocatable                :: nd(:,:), pm(:)
   integer(I4P)                             :: p, first, last, nn, k, m

   allocate(target%patch_root(target%np), target%perm(target%ne)) ; target%patch_root = 0
   allocate(target%node_box(6, max(1_I4P, 2 * target%ne)), target%node(4, max(1_I4P, 2 * target%ne)))
   nn = 0 ; first = 1
   do p = 1, target%np
      last = first - 1
      do while (last < target%ne)
         if (target%el_patch(last + 1) /= p) exit
         last = last + 1
      enddo
      if (last >= first) then
         call bvh_build(boxes(:, first:last), nb, nd, pm)
         m = size(nd, 2, kind=I4P)
         do k = 1, m ! to the global numbering
            if (nd(1, k) > 0) nd(1:2, k) = nd(1:2, k) + nn
            nd(3, k) = nd(3, k) + first - 1
         enddo
         target%perm(first:last) = pm + first - 1
         target%patch_root(p) = nn + 1
         target%node_box(:, nn + 1:nn + m) = nb ; target%node(:, nn + 1:nn + m) = nd
         nn = nn + m
      endif
      first = last + 1
   enddo
   endsubroutine build_patch_trees

   subroutine bvh_build(boxes, nbox, nodes, perm)
   !< A tree over items with boxes [6, n]: median split on the centroids along the longest axis of the node box (ties by
   !< item index), leaves of at most LEAF items. nodes(:, k) = [left, right, first, count] (left = 0: leaf); the items
   !< of a node are perm(first:first+count-1). Node 1 is the root.
   real(R8P),                 intent(in)  :: boxes(:,:) !< Item boxes [6, n].
   real(R8P),    allocatable, intent(out) :: nbox(:,:)  !< Node boxes [6, nn].
   integer(I4P), allocatable, intent(out) :: nodes(:,:) !< Nodes [4, nn].
   integer(I4P), allocatable, intent(out) :: perm(:)    !< Item order [n].
   real(R8P),    allocatable              :: b(:,:)
   integer(I4P), allocatable              :: l(:,:), todo(:)
   integer(I4P)                           :: n, nn, top, nd, first, cnt, axis, mid, i

   n = size(boxes, 2, kind=I4P)
   allocate(perm(n)) ; perm = [(i, i = 1, n)]
   if (n == 0) then
      allocate(nbox(6, 0), nodes(4, 0)) ; return
   endif
   allocate(b(6, 2 * n), l(4, 2 * n), todo(2 * n))
   nn = 1 ; l(:, 1) = [0_I4P, 0_I4P, 1_I4P, n] ; top = 1 ; todo(1) = 1
   do while (top > 0)
      nd = todo(top) ; top = top - 1
      first = l(3, nd) ; cnt = l(4, nd)
      b(1:3, nd) = minval(boxes(1:3, perm(first:first + cnt - 1)), dim=2)
      b(4:6, nd) = maxval(boxes(4:6, perm(first:first + cnt - 1)), dim=2)
      if (cnt <= LEAF) cycle
      axis = maxloc(b(4:6, nd) - b(1:3, nd), dim=1)
      mid = cnt / 2
      call select_items(perm(first:first + cnt - 1), boxes, axis, mid)
      l(1, nd) = nn + 1 ; l(2, nd) = nn + 2
      l(:, nn + 1) = [0_I4P, 0_I4P, first, mid] ; l(:, nn + 2) = [0_I4P, 0_I4P, first + mid, cnt - mid]
      top = top + 1 ; todo(top) = nn + 1 ; top = top + 1 ; todo(top) = nn + 2
      nn = nn + 2
   enddo
   nbox = b(:, 1:nn) ; nodes = l(:, 1:nn)
   endsubroutine bvh_build

   subroutine select_items(idx, boxes, axis, k)
   !< Quickselect of idx by (centroid on axis, item index) so the first k come first.
   integer(I4P), intent(inout) :: idx(:)     !< Items.
   real(R8P),    intent(in)    :: boxes(:,:) !< Boxes.
   integer(I4P), intent(in)    :: axis       !< Axis.
   integer(I4P), intent(in)    :: k          !< Split.
   integer(I4P)                :: lo, hi, i, j, t, pv

   lo = 1 ; hi = size(idx, kind=I4P)
   do while (lo < hi)
      pv = idx((lo + hi) / 2)
      i = lo ; j = hi
      do while (i <= j)
         do while (item_before(boxes, axis, idx(i), pv))
            i = i + 1
         enddo
         do while (item_before(boxes, axis, pv, idx(j)))
            j = j - 1
         enddo
         if (i <= j) then
            t = idx(i) ; idx(i) = idx(j) ; idx(j) = t
            i = i + 1 ; j = j - 1
         endif
      enddo
      if (k <= j) then
         hi = j
      elseif (k >= i) then
         lo = i
      else
         exit
      endif
   enddo
   endsubroutine select_items

   pure logical function item_before(boxes, axis, x, y)
   !< Item x before item y: centroid on axis, then index.
   real(R8P),    intent(in) :: boxes(:,:) !< Boxes.
   integer(I4P), intent(in) :: axis, x, y !< Axis, items.
   real(R8P)                :: kx, ky

   kx = boxes(axis, x) + boxes(axis + 3, x) ; ky = boxes(axis, y) + boxes(axis + 3, y)
   item_before = kx < ky .or. (kx == ky .and. x < y)
   endfunction item_before

   subroutine patch_boundary_edges(fv, fpatch, pe)
   !< Edges shared by exactly two facets of different patches [2, n] (low id first).
   integer(I4P),              intent(in)  :: fv(:,:)   !< Facets [3, nf].
   integer(I4P),              intent(in)  :: fpatch(:) !< Patches [nf].
   integer(I4P), allocatable, intent(out) :: pe(:,:)   !< Edges.
   integer(I8P), allocatable              :: key(:)
   integer(I4P), allocatable              :: order(:), ef(:), tmp(:,:)
   integer(I4P)                           :: nf, i, k, a, b, m, n, j
   integer(I8P)                           :: nv

   nf = size(fv, 2, kind=I4P) ; nv = int(maxval([1_I4P, fv]), I8P) + 1_I8P
   allocate(key(3 * nf), order(3 * nf), ef(3 * nf))
   do i = 1, nf
      do k = 1, 3
         a = fv(k, i) ; b = fv(mod(k, 3) + 1, i)
         key(3 * (i - 1) + k) = int(min(a, b), I8P) * nv + max(a, b) ; ef(3 * (i - 1) + k) = i
      enddo
   enddo
   call sort_keys(key, order)
   allocate(tmp(2, 3 * nf)) ; n = 0 ; i = 1
   do while (i <= 3 * nf)
      j = i
      do while (j < 3 * nf)
         if (key(order(j + 1)) /= key(order(i))) exit
         j = j + 1
      enddo
      m = j - i + 1
      if (m == 2) then
         if (fpatch(ef(order(i))) /= fpatch(ef(order(j)))) then
            n = n + 1
            tmp(:, n) = [int(key(order(i)) / nv, I4P), int(mod(key(order(i)), nv), I4P)]
         endif
      endif
      i = j + 1
   enddo
   pe = tmp(:, 1:n)
   endsubroutine patch_boundary_edges

   subroutine unique_edges(ed, nv)
   !< Sort the edges (low id first) and drop duplicates.
   integer(I4P), allocatable, intent(inout) :: ed(:,:) !< Edges [2, n].
   integer(I4P),              intent(in)    :: nv      !< Vertices.
   integer(I8P), allocatable                :: key(:)
   integer(I4P), allocatable                :: order(:), out(:,:)
   integer(I4P)                             :: n, i, m
   integer(I8P)                             :: w

   n = size(ed, 2, kind=I4P) ; w = int(nv, I8P) + 1_I8P
   allocate(key(n), order(n), out(2, n))
   do i = 1, n
      key(i) = int(min(ed(1, i), ed(2, i)), I8P) * w + max(ed(1, i), ed(2, i))
   enddo
   call sort_keys(key, order)
   m = 0
   do i = 1, n
      if (i > 1) then
         if (key(order(i)) == key(order(i - 1))) cycle
      endif
      m = m + 1 ; out(:, m) = [int(key(order(i)) / w, I4P), int(mod(key(order(i)), w), I4P)]
   enddo
   ed = out(:, 1:m)
   endsubroutine unique_edges

   subroutine chain_curves(ed, v, corner_angle, target)
   !< Chain the constraint edges into curves between corners (valence /= 2, or turning above corner_angle); loops
   !< without a corner are closed curves, started at their lowest vertex. Sets the curves, segments, segment tree and
   !< corners of the target.
   integer(I4P),              intent(in)    :: ed(:,:)      !< Unique edges [2, ne].
   real(R8P),                 intent(in)    :: v(:,:)       !< Vertices [3, nv].
   real(R8P),                 intent(in)    :: corner_angle !< Corner threshold [degrees].
   type(projection_target_t), intent(inout) :: target       !< Target.
   integer(I4P), allocatable                :: val(:), adj_start(:), adj(:), fill(:), cf(:), sc(:)
   logical,      allocatable                :: corner(:), used(:)
   real(R8P),    allocatable                :: sg(:,:), boxes(:,:)
   integer(I4P)                             :: ne, nv, i, k, e, a, b, c, ns, nc, cur, prv
   real(R8P)                                :: d1(3), d2(3), cosc

   ne = size(ed, 2, kind=I4P) ; nv = size(v, 2, kind=I4P)
   allocate(val(nv), adj_start(nv + 1), fill(nv), adj(2 * ne), corner(nv), used(ne)) ; val = 0
   do e = 1, ne
      val(ed(1, e)) = val(ed(1, e)) + 1 ; val(ed(2, e)) = val(ed(2, e)) + 1
   enddo
   adj_start(1) = 1
   do i = 1, nv
      adj_start(i + 1) = adj_start(i) + val(i)
   enddo
   fill = adj_start(1:nv)
   do e = 1, ne
      adj(fill(ed(1, e))) = e ; fill(ed(1, e)) = fill(ed(1, e)) + 1
      adj(fill(ed(2, e))) = e ; fill(ed(2, e)) = fill(ed(2, e)) + 1
   enddo
   cosc = cos(corner_angle * acos(-1._R8P) / 180._R8P)
   corner = .false.
   do i = 1, nv
      if (val(i) == 0) cycle
      if (val(i) /= 2) then
         corner(i) = .true. ; cycle
      endif
      a = other_end(ed(:, adj(adj_start(i))), i) ; b = other_end(ed(:, adj(adj_start(i) + 1)), i)
      d1 = v(:, i) - v(:, a) ; d2 = v(:, b) - v(:, i)
      corner(i) = dot_product(d1, d2) < cosc * norm2(d1) * norm2(d2)
   enddo
   allocate(sg(6, ne), sc(ne), cf(ne + 1)) ; used = .false. ; ns = 0 ; nc = 0
   ! open curves from the corners (vertices in increasing order, edges in adjacency order)
   do i = 1, nv
      if (.not. corner(i)) cycle
      do k = adj_start(i), adj_start(i + 1) - 1
         e = adj(k)
         if (used(e)) cycle
         nc = nc + 1 ; cf(nc) = ns + 1 ; prv = i
         do
            used(e) = .true. ; cur = other_end(ed(:, e), prv)
            ns = ns + 1 ; sg(:, ns) = [v(:, prv), v(:, cur)] ; sc(ns) = nc
            if (corner(cur)) exit
            e = next_edge(cur, e, val, adj_start, adj, used)
            if (e == 0) exit
            prv = cur
         enddo
      enddo
   enddo
   ! closed curves without a corner
   do e = 1, ne
      if (used(e)) cycle
      nc = nc + 1 ; cf(nc) = ns + 1
      prv = min(ed(1, e), ed(2, e)) ; c = e
      do
         used(c) = .true. ; cur = other_end(ed(:, c), prv)
         ns = ns + 1 ; sg(:, ns) = [v(:, prv), v(:, cur)] ; sc(ns) = nc
         c = next_edge(cur, c, val, adj_start, adj, used)
         if (c == 0) exit
         if (used(c)) exit
         prv = cur
      enddo
   enddo
   cf(nc + 1) = ns + 1
   target%nc = nc ; target%curve_first = cf(1:nc + 1) ; target%seg = sg(:, 1:ns) ; target%seg_curve = sc(1:ns)
   allocate(boxes(6, ns))
   do k = 1, ns
      boxes(1:3, k) = min(sg(1:3, k), sg(4:6, k)) ; boxes(4:6, k) = max(sg(1:3, k), sg(4:6, k))
   enddo
   call bvh_build(boxes, target%seg_box, target%seg_node, target%seg_perm)
   target%seg_root = merge(1_I4P, 0_I4P, ns > 0)
   target%nk = count(corner)
   allocate(target%corner(3, target%nk)) ; k = 0
   do i = 1, nv
      if (.not. corner(i)) cycle
      k = k + 1 ; target%corner(:, k) = v(:, i)
   enddo
   endsubroutine chain_curves

   pure integer(I4P) function other_end(edge, x)
   !< The end of edge that is not x.
   integer(I4P), intent(in) :: edge(2), x !< Edge, vertex.

   other_end = merge(edge(2), edge(1), edge(1) == x)
   endfunction other_end

   pure integer(I4P) function next_edge(x, from, val, adj_start, adj, used)
   !< The other unused edge at a valence-2 vertex x (0 if none).
   integer(I4P), intent(in) :: x, from              !< Vertex, edge arrived by.
   integer(I4P), intent(in) :: val(:), adj_start(:), adj(:) !< Adjacency.
   logical,      intent(in) :: used(:)              !< Edges chained.
   integer(I4P)             :: kk

   next_edge = 0
   if (val(x) /= 2) return
   do kk = adj_start(x), adj_start(x + 1) - 1
      if (adj(kk) /= from .and. .not. used(adj(kk))) next_edge = adj(kk)
   enddo
   endfunction next_edge

   subroutine sort_keys(key, order)
   !< Stable merge sort: order such that key(order) is non-decreasing.
   integer(I8P), intent(in)  :: key(:)   !< Keys.
   integer(I4P), intent(out) :: order(:) !< Order.
   integer(I4P), allocatable :: tmp(:)
   integer(I4P)              :: n, w, lo, mid, hi, i, j, k

   n = size(key, kind=I4P) ; order = [(i, i = 1, n)]
   allocate(tmp(n)) ; w = 1
   do while (w < n)
      lo = 1
      do while (lo <= n)
         mid = min(lo + w - 1, n) ; hi = min(lo + 2 * w - 1, n)
         i = lo ; j = mid + 1 ; k = lo
         do while (i <= mid .and. j <= hi)
            if (key(order(j)) < key(order(i))) then
               tmp(k) = order(j) ; j = j + 1
            else
               tmp(k) = order(i) ; i = i + 1
            endif
            k = k + 1
         enddo
         tmp(k:k + mid - i) = order(i:mid) ; k = k + mid - i + 1
         tmp(k:k + hi - j) = order(j:hi)
         order(lo:hi) = tmp(lo:hi)
         lo = lo + 2 * w
      enddo
      w = 2 * w
   enddo
   endsubroutine sort_keys
endmodule fossil_projection
