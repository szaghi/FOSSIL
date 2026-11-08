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
!<
!< Projection (`project_points`, `project_points_2d`): each point goes to the closest point of its target: the whole
!< body (PROJ_FREE), one patch (PROJ_PATCH, id), one constraint curve or any (PROJ_CURVE, id or 0), one corner or the
!< nearest (PROJ_CORNER, id or 0). Feet on triangles by `closest_on_triangle` (#69), on segments by the clamped
!< projection; ties in the distance go to the lowest facet / segment index (deterministic, independent of the tree).
!< Corners are returned bit for bit. `classify_points`: a corner within the tolerance, else a constraint curve within
!< it, else the patch of the closest facet.
!<
!< Block faces (`project_block_face`, 3-D faces 1-6 = -i, +i, -j, +j, -k, +k; `project_block_side`, 2-D sides 1-4): every
!< vertex of the face is classified with the tolerance (so face corners near body corners go to them, edge vertices
!< near constraint curves to the curves, the others to the patch of their closest facet) and moved to its foot, in
!< place; `moved` gives the distance each vertex moved.
use penf,                 only : I4P, I8P, R8P
use fossil_feature_edges, only : find_feature_edges
use fossil_body_distance, only : closest_on_triangle
use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t

implicit none
private
public :: projection_target_t, build_target_3d, build_target_2d
public :: PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
public :: bvh_build
public :: project_points, project_points_2d, classify_points, classify_points_2d
public :: PROJ_STATUS_OK, PROJ_STATUS_BAD_INPUT
public :: project_block_face, project_block_side

integer(I4P), parameter :: PROJ_FREE   = 0_I4P !< The whole body.
integer(I4P), parameter :: PROJ_PATCH  = 1_I4P !< One patch (id).
integer(I4P), parameter :: PROJ_CURVE  = 2_I4P !< One constraint curve (id > 0), or any (id 0); 3-D only.
integer(I4P), parameter :: PROJ_CORNER = 3_I4P !< One corner (id > 0), or the nearest (id 0); fixed.
integer(I4P), parameter :: LEAF        = 4_I4P !< Elements per leaf.
integer(I4P), parameter :: PROJ_STATUS_OK        = 0_I4P !< Every point projected.
integer(I4P), parameter :: PROJ_STATUS_BAD_INPUT = 1_I4P !< A kind or id out of range, or a target not built: those
                                                         !< points keep their position, distance -1.

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

   subroutine project_points(target, x, y, z, kinds, ids, fx, fy, fz, distances, status)
   !< Feet of the points on their targets (3-D; see the module notes).
   type(projection_target_t), intent(in)            :: target           !< Target (dims 3).
   real(R8P),                 intent(in)            :: x(:), y(:), z(:) !< Points [n].
   integer(I4P),              intent(in)            :: kinds(:), ids(:) !< PROJ_* and ids [n].
   real(R8P),                 intent(out)           :: fx(:), fy(:), fz(:) !< Feet [n].
   real(R8P),                 intent(out), optional :: distances(:)     !< Distances [n] (-1: not projected).
   integer(I4P),              intent(out), optional :: status           !< PROJ_STATUS_*.
   real(R8P)                                        :: p(3), f(3), d2
   integer(I4P)                                     :: i, el, st
   logical                                          :: ok

   st = PROJ_STATUS_OK
   do i = 1, size(x, kind=I4P)
      p = [x(i), y(i), z(i)]
      call project_one(target, p, kinds(i), ids(i), f, d2, el, ok)
      if (.not. ok) then
         st = PROJ_STATUS_BAD_INPUT ; f = p ; d2 = -1._R8P
      endif
      fx(i) = f(1) ; fy(i) = f(2) ; fz(i) = f(3)
      if (present(distances)) distances(i) = merge(sqrt(max(0._R8P, d2)), -1._R8P, d2 >= 0._R8P)
   enddo
   if (present(status)) status = st
   endsubroutine project_points

   subroutine project_points_2d(target, u, v, kinds, ids, fu, fv, distances, status)
   !< Feet of the points on their targets (2-D: PROJ_FREE, PROJ_PATCH, PROJ_CORNER).
   type(projection_target_t), intent(in)            :: target           !< Target (dims 2).
   real(R8P),                 intent(in)            :: u(:), v(:)       !< Points [n].
   integer(I4P),              intent(in)            :: kinds(:), ids(:) !< PROJ_* and ids [n].
   real(R8P),                 intent(out)           :: fu(:), fv(:)     !< Feet [n].
   real(R8P),                 intent(out), optional :: distances(:)     !< Distances [n] (-1: not projected).
   integer(I4P),              intent(out), optional :: status           !< PROJ_STATUS_*.
   real(R8P)                                        :: p(3), f(3), d2
   integer(I4P)                                     :: i, el, st
   logical                                          :: ok

   st = PROJ_STATUS_OK
   do i = 1, size(u, kind=I4P)
      p = [u(i), v(i), 0._R8P]
      ok = kinds(i) /= PROJ_CURVE
      if (ok) call project_one(target, p, kinds(i), ids(i), f, d2, el, ok)
      if (.not. ok) then
         st = PROJ_STATUS_BAD_INPUT ; f = p ; d2 = -1._R8P
      endif
      fu(i) = f(1) ; fv(i) = f(2)
      if (present(distances)) distances(i) = merge(sqrt(max(0._R8P, d2)), -1._R8P, d2 >= 0._R8P)
   enddo
   if (present(status)) status = st
   endsubroutine project_points_2d

   subroutine classify_points(target, x, y, z, tolerance, kinds, ids)
   !< Constraints of the points: PROJ_CORNER (a corner within tolerance), else PROJ_CURVE (a constraint curve within
   !< it), else PROJ_PATCH with the patch of the closest facet.
   type(projection_target_t), intent(in)  :: target           !< Target (dims 3).
   real(R8P),                 intent(in)  :: x(:), y(:), z(:) !< Points [n].
   real(R8P),                 intent(in)  :: tolerance        !< Tolerance.
   integer(I4P),              intent(out) :: kinds(:), ids(:) !< PROJ_* and ids [n].
   integer(I4P)                           :: i

   do i = 1, size(x, kind=I4P)
      call classify_one(target, [x(i), y(i), z(i)], tolerance, kinds(i), ids(i))
   enddo
   endsubroutine classify_points

   subroutine classify_points_2d(target, u, v, tolerance, kinds, ids)
   !< 2-D twin of classify_points: PROJ_CORNER or PROJ_PATCH.
   type(projection_target_t), intent(in)  :: target           !< Target (dims 2).
   real(R8P),                 intent(in)  :: u(:), v(:)       !< Points [n].
   real(R8P),                 intent(in)  :: tolerance        !< Tolerance.
   integer(I4P),              intent(out) :: kinds(:), ids(:) !< PROJ_* and ids [n].
   integer(I4P)                           :: i

   do i = 1, size(u, kind=I4P)
      call classify_one(target, [u(i), v(i), 0._R8P], tolerance, kinds(i), ids(i))
   enddo
   endsubroutine classify_points_2d

   subroutine project_block_face(block, face, target, tolerance, moved, kinds, status)
   !< Classify and project the vertices of one face of a 3-D block onto the body, in place (see the module notes).
   class(curvilinear_block_t), intent(inout)         :: block      !< Block.
   integer(I4P),               intent(in)            :: face       !< Face 1-6 (-i, +i, -j, +j, -k, +k).
   type(projection_target_t),  intent(in)            :: target     !< Target (dims 3).
   real(R8P),                  intent(in)            :: tolerance  !< Classification tolerance.
   real(R8P), allocatable,     intent(out), optional :: moved(:,:) !< Distance moved per face vertex.
   integer(I4P), allocatable,  intent(out), optional :: kinds(:,:) !< Constraint per face vertex (PROJ_*).
   integer(I4P),               intent(out), optional :: status     !< PROJ_STATUS_*.
   integer(I4P)                                      :: lo(3), hi(3), a, n1, n2, i, j, kk, ii, st
   integer(I4P)                                      :: idx(3)
   real(R8P)                                         :: p(3), f(3), d2
   integer(I4P)                                      :: el
   logical                                           :: ok

   st = PROJ_STATUS_OK
   if (face < 1 .or. face > 6 .or. target%dims /= 3 .or. .not. allocated(block%x)) then
      if (present(status)) status = PROJ_STATUS_BAD_INPUT
      return
   endif
   lo = 1 ; hi = [block%ni, block%nj, block%nk]
   a = (face + 1) / 2 ! the fixed direction
   if (mod(face, 2) == 1) then
      hi(a) = 1
   else
      lo(a) = hi(a)
   endif
   select case (a)
   case (1) ; n1 = block%nj ; n2 = block%nk
   case (2) ; n1 = block%ni ; n2 = block%nk
   case default ; n1 = block%ni ; n2 = block%nj
   endselect
   if (present(moved)) allocate(moved(n1, n2))
   if (present(kinds)) allocate(kinds(n1, n2))
   do j = 1, n2
      do i = 1, n1
         select case (a)
         case (1) ; idx = [lo(1), i, j]
         case (2) ; idx = [i, lo(2), j]
         case default ; idx = [i, j, lo(3)]
         endselect
         p = [block%x(idx(1), idx(2), idx(3)), block%y(idx(1), idx(2), idx(3)), block%z(idx(1), idx(2), idx(3))]
         call classify_one(target, p, tolerance, kk, ii)
         call project_one(target, p, kk, ii, f, d2, el, ok)
         if (.not. ok) then
            st = PROJ_STATUS_BAD_INPUT ; f = p
         endif
         block%x(idx(1), idx(2), idx(3)) = f(1) ; block%y(idx(1), idx(2), idx(3)) = f(2)
         block%z(idx(1), idx(2), idx(3)) = f(3)
         if (present(moved)) moved(i, j) = norm2(f - p)
         if (present(kinds)) kinds(i, j) = kk
      enddo
   enddo
   if (present(status)) status = st
   endsubroutine project_block_face

   subroutine project_block_side(block, side, target, tolerance, moved, kinds, status)
   !< Classify and project the vertices of one side of a 2-D block onto the curve, in place (sides 1-4 = -i, +i, -j, +j).
   class(curvilinear_block_2d_t), intent(inout)         :: block     !< Block.
   integer(I4P),                  intent(in)            :: side      !< Side 1-4.
   type(projection_target_t),     intent(in)            :: target    !< Target (dims 2).
   real(R8P),                     intent(in)            :: tolerance !< Classification tolerance.
   real(R8P), allocatable,        intent(out), optional :: moved(:)  !< Distance moved per side vertex.
   integer(I4P), allocatable,     intent(out), optional :: kinds(:)  !< Constraint per side vertex (PROJ_*).
   integer(I4P),                  intent(out), optional :: status    !< PROJ_STATUS_*.
   integer(I4P)                                         :: n, i, iu, iv, kk, ii, el, st
   real(R8P)                                            :: p(3), f(3), d2
   logical                                              :: ok

   st = PROJ_STATUS_OK
   if (side < 1 .or. side > 4 .or. target%dims /= 2 .or. .not. allocated(block%u)) then
      if (present(status)) status = PROJ_STATUS_BAD_INPUT
      return
   endif
   n = merge(block%nj, block%ni, side <= 2)
   if (present(moved)) allocate(moved(n))
   if (present(kinds)) allocate(kinds(n))
   do i = 1, n
      select case (side)
      case (1) ; iu = 1 ; iv = i
      case (2) ; iu = block%ni ; iv = i
      case (3) ; iu = i ; iv = 1
      case default ; iu = i ; iv = block%nj
      endselect
      p = [block%u(iu, iv), block%v(iu, iv), 0._R8P]
      call classify_one(target, p, tolerance, kk, ii)
      call project_one(target, p, kk, ii, f, d2, el, ok)
      if (.not. ok) then
         st = PROJ_STATUS_BAD_INPUT ; f = p
      endif
      block%u(iu, iv) = f(1) ; block%v(iu, iv) = f(2)
      if (present(moved)) moved(i) = norm2(f - p)
      if (present(kinds)) kinds(i) = kk
   enddo
   if (present(status)) status = st
   endsubroutine project_block_side

   subroutine classify_one(target, p, tolerance, kind, id)
   !< Constraint of one point.
   type(projection_target_t), intent(in)  :: target    !< Target.
   real(R8P),                 intent(in)  :: p(3)      !< Point.
   real(R8P),                 intent(in)  :: tolerance !< Tolerance.
   integer(I4P),              intent(out) :: kind, id  !< Constraint.
   real(R8P)                              :: f(3), d2, t2
   integer(I4P)                           :: el, k
   logical                                :: ok

   t2 = tolerance * tolerance
   call nearest_corner(target, p, k, d2)
   if (k > 0 .and. d2 <= t2) then
      kind = PROJ_CORNER ; id = k ; return
   endif
   if (target%dims == 3 .and. target%seg_root > 0) then
      call closest_on_curves(target, p, 0_I4P, f, d2, el)
      if (el > 0 .and. d2 <= t2) then
         kind = PROJ_CURVE ; id = target%seg_curve(el) ; return
      endif
   endif
   call project_one(target, p, PROJ_FREE, 0_I4P, f, d2, el, ok)
   kind = PROJ_PATCH ; id = 0
   if (ok .and. el > 0) id = target%el_patch(el)
   endsubroutine classify_one

   subroutine project_one(target, p, kind, id, f, d2, el, ok)
   !< Foot of one point; el = the element (or segment, or corner) reached.
   type(projection_target_t), intent(in)  :: target !< Target.
   real(R8P),                 intent(in)  :: p(3)   !< Point.
   integer(I4P),              intent(in)  :: kind, id !< Constraint.
   real(R8P),                 intent(out) :: f(3)   !< Foot.
   real(R8P),                 intent(out) :: d2     !< Squared distance.
   integer(I4P),              intent(out) :: el     !< Element / segment / corner.
   logical,                   intent(out) :: ok     !< Projected.
   integer(I4P)                           :: q

   f = p ; d2 = huge(1._R8P) ; el = 0 ; ok = target%dims > 0
   if (.not. ok) return
   select case (kind)
   case (PROJ_FREE)
      do q = 1, target%np
         if (target%patch_root(q) > 0) call closest_in_tree(target, target%patch_root(q), p, f, d2, el)
      enddo
   case (PROJ_PATCH)
      ok = id >= 1 .and. id <= target%np
      if (ok) ok = target%patch_root(id) > 0
      if (ok) call closest_in_tree(target, target%patch_root(id), p, f, d2, el)
   case (PROJ_CURVE)
      ok = target%dims == 3 .and. id >= 0 .and. id <= target%nc .and. target%nc > 0
      if (ok) call closest_on_curves(target, p, id, f, d2, el)
   case (PROJ_CORNER)
      ok = id >= 0 .and. id <= target%nk .and. target%nk > 0
      if (ok) then
         if (id > 0) then
            el = id ; d2 = sum((target%corner(:, id) - p)**2)
         else
            call nearest_corner(target, p, el, d2)
         endif
         f = target%corner(:, el)
      endif
   case default
      ok = .false.
   endselect
   if (ok) ok = el > 0
   endsubroutine project_one

   subroutine closest_in_tree(target, root, p, f, d2, el)
   !< Closest element of the tree under root (depth first, box pruning; ties to the lowest original index).
   type(projection_target_t), intent(in)    :: target !< Target.
   integer(I4P),              intent(in)    :: root   !< Root node.
   real(R8P),                 intent(in)    :: p(3)   !< Point.
   real(R8P),                 intent(inout) :: f(3)   !< Best foot.
   real(R8P),                 intent(inout) :: d2     !< Best squared distance.
   integer(I4P),              intent(inout) :: el     !< Best element.
   integer(I4P)                             :: stk(128), ns, nd, k, e
   real(R8P)                                :: c(3), dd

   ns = 1 ; stk(1) = root
   do while (ns > 0)
      nd = stk(ns) ; ns = ns - 1
      if (box_d2(target%node_box(:, nd), p) > d2) cycle
      if (target%node(1, nd) == 0) then
         do k = target%node(3, nd), target%node(3, nd) + target%node(4, nd) - 1
            e = target%perm(k)
            if (target%dims == 3) then
               c = closest_on_triangle(p, reshape(target%el(:, e), [3, 3]))
            else
               c = closest_on_segment(p, [target%el(1:2, e), 0._R8P], [target%el(3:4, e), 0._R8P])
            endif
            dd = sum((c - p)**2)
            if (el == 0) then
               f = c ; d2 = dd ; el = e
            elseif (dd < d2 .or. (dd == d2 .and. target%el_id(e) < target%el_id(el))) then
               f = c ; d2 = dd ; el = e
            endif
         enddo
      else
         stk(ns + 1) = target%node(2, nd) ; stk(ns + 2) = target%node(1, nd) ; ns = ns + 2
      endif
   enddo
   endsubroutine closest_in_tree

   subroutine closest_on_curves(target, p, id, f, d2, sg)
   !< Closest constraint segment: of curve id (all its segments), or of any curve (id 0, the segment tree); ties to the
   !< lowest segment index.
   type(projection_target_t), intent(in)  :: target !< Target.
   real(R8P),                 intent(in)  :: p(3)   !< Point.
   integer(I4P),              intent(in)  :: id     !< Curve (0: any).
   real(R8P),                 intent(out) :: f(3)   !< Foot.
   real(R8P),                 intent(out) :: d2     !< Squared distance.
   integer(I4P),              intent(out) :: sg     !< Segment.
   integer(I4P)                           :: stk(128), ns, nd, k, e
   real(R8P)                              :: c(3), dd

   f = p ; d2 = huge(1._R8P) ; sg = 0
   if (id > 0) then
      do e = target%curve_first(id), target%curve_first(id + 1) - 1
         c = closest_on_segment(p, target%seg(1:3, e), target%seg(4:6, e)) ; dd = sum((c - p)**2)
         if (sg == 0 .or. dd < d2) then
            f = c ; d2 = dd ; sg = e
         endif
      enddo
      return
   endif
   if (target%seg_root == 0) return
   ns = 1 ; stk(1) = target%seg_root
   do while (ns > 0)
      nd = stk(ns) ; ns = ns - 1
      if (box_d2(target%seg_box(:, nd), p) > d2) cycle
      if (target%seg_node(1, nd) == 0) then
         do k = target%seg_node(3, nd), target%seg_node(3, nd) + target%seg_node(4, nd) - 1
            e = target%seg_perm(k)
            c = closest_on_segment(p, target%seg(1:3, e), target%seg(4:6, e)) ; dd = sum((c - p)**2)
            if (sg == 0) then
               f = c ; d2 = dd ; sg = e
            elseif (dd < d2 .or. (dd == d2 .and. e < sg)) then
               f = c ; d2 = dd ; sg = e
            endif
         enddo
      else
         stk(ns + 1) = target%seg_node(2, nd) ; stk(ns + 2) = target%seg_node(1, nd) ; ns = ns + 2
      endif
   enddo
   endsubroutine closest_on_curves

   subroutine nearest_corner(target, p, k, d2)
   !< Nearest corner (ties to the lowest index; 0 if none).
   type(projection_target_t), intent(in)  :: target !< Target.
   real(R8P),                 intent(in)  :: p(3)   !< Point.
   integer(I4P),              intent(out) :: k      !< Corner.
   real(R8P),                 intent(out) :: d2     !< Squared distance.
   integer(I4P)                           :: i
   real(R8P)                              :: dd

   k = 0 ; d2 = huge(1._R8P)
   do i = 1, target%nk
      dd = sum((target%corner(:, i) - p)**2)
      if (dd < d2) then
         k = i ; d2 = dd
      endif
   enddo
   endsubroutine nearest_corner

   pure function closest_on_segment(p, a, b) result(c)
   !< Closest point of the segment [a, b] to p (clamped projection).
   real(R8P), intent(in) :: p(3), a(3), b(3) !< Point, ends.
   real(R8P)             :: c(3)             !< Foot.
   real(R8P)             :: e(3), l2, t

   e = b - a ; l2 = dot_product(e, e)
   if (l2 > 0._R8P) then
      t = max(0._R8P, min(1._R8P, dot_product(p - a, e) / l2))
   else
      t = 0._R8P
   endif
   c = a + t * e
   if (t == 1._R8P) c = b
   endfunction closest_on_segment

   pure function box_d2(box, p) result(d2)
   !< Squared distance from p to a box [6].
   real(R8P), intent(in) :: box(6), p(3) !< Box, point.
   real(R8P)             :: d2           !< Squared distance.

   d2 = sum(max(0._R8P, box(1:3) - p, p - box(4:6))**2)
   endfunction box_d2
endmodule fossil_projection
