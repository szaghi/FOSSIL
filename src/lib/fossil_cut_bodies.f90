!< FOSSIL, cut cells of several posed bodies (#109 WP3 step 03).
module fossil_cut_bodies
!< FOSSIL, cut cells of several posed bodies (#109 WP3 step 03).
!<
!< Bodies (decision B1 of WP3). An array of surfaces (curves in 2-D), each with its own rigid pose; one surface may hold
!< several bodies of its own (`facet_body`, `loop_body`), which keep their ids. The body of a result is the pair (surface
!< index, sub-body): `body` is the surface (curve) index, `sub_body` the `facet_body` (`loop_body`) of that surface.
!<
!< Merge (B7). `merge_posed_bodies` places every body as `move` places it and merges the placed copies as
!< `merge_solids` merges (3-D: a copy of each surface moved by `move`, the copies merged in order into the first one, so
!< the merged surface is the one a caller gets with `a%move(p1)`, `b%move(p2)`, `a%merge_solids(b)`, bit for bit: same
!< facets in the same order, patches merged by name; 2-D: the vertices placed by `to_world`, as `curve%move` places
!< them, the loops concatenated in order, patches merged by name, then `adopt_loops`). It returns the map facet
!< (segment) -> (body, sub_body). The merged object is host-only (not on the device): `enter_device` it if wanted.
!<
!< Cut (B7). The cut cells of several bodies are the cut of the merged surface (`cut_cells_3d`, `cut_cells_3d_device`;
!< `curve%cut_cells`, `curve%cut_cells_device`), bit for bit, and their wall terms are attributed through the map:
!< every wall term (3-D: `cut_walls_3d_t`, one per (cell, facet); 2-D: `cut_walls_2d_t`, one per wall segment of a
!< piece's boundary) gets the (body, sub_body) of its facet (segment), and every piece the sorted list of the distinct
!< bodies its wall terms belong to (`piece_bodies_t`, CSR). A piece bounded by two bodies (two spheres in one cell)
!< lists both. A piece whose wall reaches it only along an edge or a vertex (no wall term) lists no body for it.
!< `cut_cells_3d_bodies` / `cut_cells_2d_bodies` do merge, check, cut and attribution in one call; the parts are
!< public (`merge_posed_bodies`, `label_piece_bodies`) for a caller that keeps the merged surface or its own walls.
!< Why a wrapper and not a recipe: the overlap check (below) must not be skippable by forgetting a call, and the
!< attribution needs the wall terms, which the wrapper requests itself when `piece_bodies` is asked.
!<
!< Overlap (B2). The cut of a union of overlapping bodies would need the boolean union (the merged surface is then not
!< the boundary of a solid), so the cut refuses an overlap (STATUS_INVALID_INPUT). Touching bodies (a contact of their
!< surfaces, no common volume) are allowed: the cut handles them (cubes touching along a face, an edge, a corner, on and
!< off the grid planes: every cell trusted, the solid volume the sum of the bodies'; `fossil_test_cut_bodies`). The
!< check, between bodies of different surfaces (the sub-bodies of one surface are that surface's business, as for one
!< surface alone), uses exact predicates only:
!< 1. broad phase: the closed boxes of the placed vertices; two surfaces whose boxes do not meet cannot overlap (and the
!<    check costs nothing more when no boxes meet);
!< 2. crossings, exact: an edge of one triangle piercing the interior of a triangle of another surface transversally
!<    (orient3d of the edge's ends strictly opposite, the three orient3d of the edge with the triangle's edges strictly
!<    of one sign: the edge passes from outside the other body to inside it), or two coplanar triangles of the same
!<    orientation whose interiors overlap (no edge line of either weakly separates them, exact orient2d: both solids on
!<    the same side). The candidates are every inter-surface pair of triangles whose closed boxes meet, by a sweep along
!<    x (`sweep_overlap`). Also an overlap: two sheets crossing along a shared edge (SELFINT_SHEETS) and the coplanar /
!<    folded pairs of the same orientation among the inter-surface pairs of `find_self_intersections(include_touching=)`
!<    of the merged surface, whose count is `contacts`. The sweep is ours because `find_self_intersections` misses pairs
!<    in degenerate configurations: two unit cubes offset by exactly 0.5 along y give no pair at all (0.25 or 0.7 give
!<    240 / 128). 2-D: a proper crossing of two segments, or two collinear segments overlapping over a positive length in
!<    the same direction, among the inter-curve pairs of `curve%self_intersections` (exact orient2d);
!< 3. containment: with no such crossing, every shell (loop) of one body lies inside or outside the other. One vertex
!<    per shell that lies on no facet of the other body (exact: orient3d zero and inside the closed triangle; 2-D
!<    orient2d zero and inside the closed segment) is tested by the exact parity of the crossings of the +z line
!<    (2-D: +u) below it with the other body (`tri_line_hit_oac` and `line_crossing_before`, the cut's own kernels):
!<    inside means an overlap. A shell all of whose vertices lie on the other body is counted an overlap (coincident
!<    bodies; conservative).
!< Every positive answer is a true overlap, but for the conservative one of step 3 (a shell on the other body). A
!< crossing made only of edges and vertices lying exactly on the other body (no transversal piercing of a facet
!< interior, no vertex strictly inside) is decided by step 3 alone; such configurations need coordinates aligned to the
!< last bit and were not met by the tests. Cost, when boxes meet: the sweep (O(N log N) plus the candidate pairs), the
!< pairs of `find_self_intersections`, and O(N) per shell for step 3.
!< Open loops (2-D) are not tested for containment (they bound no area).
!<
!< Device. `device=.true.` cuts on the device (`cut_cells_3d_device`, `curve%cut_cells_device`), bit for bit the host
!< cut; the merge, the check and the attribution run on the host. 3-D: the wall terms come from the device cut
!< (`walls=`); with `resident=.true.` they are copied back for the attribution (`update_host`; they stay mapped when the
!< caller passed `walls`). 2-D: the device cut has no oriented boundaries, so the attribution runs a host cut with
!< boundaries for the labels (the same pieces, bit for bit). Moments and quadrature (2-D also boundaries) are host-only,
!< as for one body: asked with `device` they are refused (STATUS_INVALID_INPUT).
!<
!< CPU only for the merge and the check; free routines, re-exported by `use fossil`.
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_cut_cells_2d,       only : cut_cells_2d_reset
use fossil_cut_cells_3d,       only : cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, cut_cells_3d_reset, cut_walls_3d_t, &
                                      cut_walls_reset
use fossil_cut_graph,          only : cut_links_t, cut_links_reset
use fossil_cut_loads,          only : cut_walls_2d_t, cut_walls_2d
use fossil_facet_object,       only : facet_object
use fossil_line_query,         only : line_crossing_before, line_crossing_before_2d, seg_line_hit_oac, tri_line_hit_oac
use fossil_predicates,         only : orient2d_xy, orient3d_xyz
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t, to_world
use fossil_self_intersection,  only : intersection_pair_t, SELFINT_COPLANAR, SELFINT_EDGE_FOLD, SELFINT_SHEETS
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P

implicit none
private
public :: piece_bodies_t, merge_posed_bodies, label_piece_bodies, cut_cells_3d_bodies, cut_cells_2d_bodies

type :: piece_bodies_t
   !< The bodies of the wall terms and of the pieces of a cut of several bodies (see the module header).
   !<
   !< Wall terms in the order of the cut's (3-D: `cut_walls_3d_t`, cell order then facet order; 2-D: `cut_walls_2d_t`,
   !< piece order then boundary order). Piece p lists the distinct bodies of its wall terms, sorted by (body, sub_body),
   !< as entries start(p) .. start(p + 1) - 1 of `body` and `sub_body`; a piece without wall term lists none.
   integer(I4P)              :: np = 0_I4P        !< Pieces.
   integer(I4P)              :: nw = 0_I4P        !< Wall terms.
   integer(I4P), allocatable :: wall_body(:)      !< Body (surface or curve index) of each wall term [nw].
   integer(I4P), allocatable :: wall_sub_body(:)  !< Its sub-body (`facet_body` / `loop_body` of that surface) [nw].
   integer(I4P), allocatable :: start(:)          !< Bodies of piece p: start(p) .. start(p + 1) - 1 [np + 1].
   integer(I4P), allocatable :: body(:)           !< Body of each entry.
   integer(I4P), allocatable :: sub_body(:)       !< Sub-body of each entry.
endtype piece_bodies_t

interface merge_posed_bodies
   !< One surface (curve) with every body placed by its pose, and the map facet (segment) -> (body, sub_body).
   module procedure merge_posed_surfaces, merge_posed_curves
endinterface merge_posed_bodies

interface label_piece_bodies
   !< The bodies of the wall terms and pieces of a cut of the merged bodies (`piece_bodies_t`).
   module procedure label_piece_bodies_3d, label_piece_bodies_2d
endinterface label_piece_bodies

contains
   subroutine merge_posed_surfaces(surfaces, poses, merged, facet_body, facet_sub_body, status, overlap, contacts)
   !< The surfaces placed by their poses and merged (see the module header): facet f of the result is facet
   !< f - offset(s) of surface s, `facet_body(f)` = s, `facet_sub_body(f)` = `surfaces(s)%facet_body(f - offset(s))`.
   !< With `overlap` (or `contacts`) the bodies are checked (the module header): `overlap` true if two surfaces' volumes
   !< intersect; `contacts` the inter-surface facet pairs in contact. The merge itself never refuses an overlap (the FD
   !< closure reports overlaps, B2). `status = STATUS_INVALID_INPUT` (merged empty) if there is no surface or the poses
   !< do not match them. A merged surface mapped on the device is released first.
   type(surface_stl_object),  intent(in)            :: surfaces(:)       !< Bodies (body frame).
   type(rigid_pose_t),        intent(in)            :: poses(:)          !< Pose of each surface.
   type(surface_stl_object),  intent(inout)         :: merged            !< Merged posed surface.
   integer(I4P), allocatable, intent(out)           :: facet_body(:)     !< Surface of each facet.
   integer(I4P), allocatable, intent(out)           :: facet_sub_body(:) !< Its facet_body in that surface.
   integer(I4P),              intent(out), optional :: status            !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                   intent(out), optional :: overlap           !< Two surfaces' volumes intersect.
   integer(I4P),              intent(out), optional :: contacts          !< Inter-surface facet pairs in contact.
   type(surface_stl_object)  :: work, empty                              !< Placed copy, empty surface.
   integer(I4P), allocatable :: offset(:)                                !< First facet - 1 of each surface.
   integer(I4P)              :: ns, nf, s, f, st                         !< Work.
   logical                   :: ov                                       !< Overlap.

   if (present(status)) status = STATUS_OK
   if (present(overlap)) overlap = .false.
   if (present(contacts)) contacts = 0_I4P
   if (merged%is_on_device) call merged%exit_device
   ns = size(surfaces, kind=I4P)
   if (ns < 1_I4P .or. size(poses, kind=I4P) /= ns) then
      merged = empty
      allocate(facet_body(0), facet_sub_body(0))
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(offset(ns + 1))
   offset(1) = 0_I4P
   do s = 1, ns
      offset(s + 1) = offset(s) + surfaces(s)%get_facets_number()
   enddo
   nf = offset(ns + 1)
   allocate(facet_body(nf), facet_sub_body(nf))
   do s = 1, ns
      do f = 1, surfaces(s)%get_facets_number()
         facet_body(offset(s) + f) = s ; facet_sub_body(offset(s) + f) = surfaces(s)%facet_body(f)
      enddo
      work = surfaces(s)
      work%is_on_device = .false. ! a copy: its host arrays are not the mapped ones (move would release the original's)
      call work%move(poses(s))    ! the placement of `move`, bit for bit
      if (s == 1_I4P) then
         merged = work
      else
         call merged%merge_solids(work, status=st)
         if (st /= STATUS_OK) then
            merged = empty
            if (present(status)) status = STATUS_INVALID_INPUT
            return
         endif
      endif
   enddo
   merged%is_on_device = .false.
   if (present(overlap) .or. present(contacts)) then
      call overlap_3d(surfaces, merged, offset, ov, st)
      if (present(overlap)) overlap = ov
      if (present(contacts)) contacts = st
   endif
   endsubroutine merge_posed_surfaces

   subroutine merge_posed_curves(curves, poses, merged, segment_body, segment_sub_body, status, overlap, contacts)
   !< The curves placed by their poses and merged (see the module header): the loops of curve 1, then of curve 2, ...;
   !< segment i of the result (its first vertex) is segment i - offset(c) of curve c, `segment_body(i)` = c,
   !< `segment_sub_body(i)` its `loop_body`. With `overlap` (or `contacts`) the curves are checked (module header).
   !< `status = STATUS_INVALID_INPUT` (merged empty) if there is no curve, the poses do not match them or the curves lie
   !< in planes of different normal axes.
   type(curve_object),        intent(in)            :: curves(:)           !< Bodies (body frame).
   type(rigid_pose_2d_t),     intent(in)            :: poses(:)            !< Pose of each curve.
   type(curve_object),        intent(inout)         :: merged              !< Merged posed curve.
   integer(I4P), allocatable, intent(out)           :: segment_body(:)     !< Curve of each segment.
   integer(I4P), allocatable, intent(out)           :: segment_sub_body(:) !< Its loop_body in that curve.
   integer(I4P),              intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                   intent(out), optional :: overlap             !< Two curves' areas intersect.
   integer(I4P),              intent(out), optional :: contacts            !< Inter-curve segment pairs in contact.
   real(R8P),        allocatable :: u(:), v(:)                             !< Placed vertices.
   integer(I4P),     allocatable :: sizes(:), patches(:), offset(:), pmap(:) !< Loop sizes, patches, offsets, patch map.
   logical,          allocatable :: closed(:)                              !< Loops closed.
   character(64),    allocatable :: names(:)                               !< Patch names.
   character(len=:), allocatable :: nm                                     !< Name.
   real(R8P)                     :: w(2)                                   !< Placed vertex.
   integer(I4P)                  :: nc, nv, nl, c, l, i, first, last, kl, p, q, st, axis !< Work.
   logical                       :: ov                                     !< Overlap.

   if (present(status)) status = STATUS_OK
   if (present(overlap)) overlap = .false.
   if (present(contacts)) contacts = 0_I4P
   if (merged%get_is_on_device()) call merged%exit_device
   nc = size(curves, kind=I4P)
   st = STATUS_OK
   if (nc < 1_I4P .or. size(poses, kind=I4P) /= nc) then
      st = STATUS_INVALID_INPUT
   else
      axis = curves(1)%get_normal_axis()
      do c = 2, nc
         if (curves(c)%get_normal_axis() /= axis) st = STATUS_INVALID_INPUT
      enddo
   endif
   if (st /= STATUS_OK) then
      call merged%destroy
      allocate(segment_body(0), segment_sub_body(0))
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(offset(nc + 1))
   offset(1) = 0_I4P ; nl = 0_I4P
   do c = 1, nc
      offset(c + 1) = offset(c) + curves(c)%get_vertices_number()
      nl = nl + curves(c)%get_loops_number()
   enddo
   nv = offset(nc + 1)
   allocate(u(nv), v(nv), segment_body(nv), segment_sub_body(nv), patches(nv), sizes(nl), closed(nl))
   ! patch names: curve 1's table, then the new names of the others (equal names merge)
   allocate(names(0))
   kl = 0_I4P
   do c = 1, nc
      allocate(pmap(curves(c)%get_patches_number()))
      do p = 1, size(pmap)
         nm = curves(c)%patch_name(p)
         pmap(p) = 0_I4P
         if (c > 1_I4P) then
            do q = 1, size(names)
               if (trim(names(q)) == nm) then
                  pmap(p) = q ; exit
               endif
            enddo
         endif
         if (pmap(p) == 0_I4P) then
            names = [names, repeat(' ', 64)] ; names(size(names)) = nm ; pmap(p) = size(names, kind=I4P)
         endif
      enddo
      do l = 1, curves(c)%get_loops_number()
         call curves(c)%get_loop_range(l, first, last)
         kl = kl + 1_I4P
         sizes(kl) = max(0_I4P, last - first + 1_I4P) ; closed(kl) = curves(c)%is_loop_closed(l)
         do i = first, last
            w = to_world(poses(c), curves(c)%get_vertex(i)) ! the placement of `curve%move`, bit for bit
            u(offset(c) + i) = w(1) ; v(offset(c) + i) = w(2)
            segment_body(offset(c) + i) = c ; segment_sub_body(offset(c) + i) = curves(c)%loop_body(l)
            p = curves(c)%segment_patch(i)
            patches(offset(c) + i) = pmap(max(1_I4P, min(p, size(pmap, kind=I4P))))
         enddo
      enddo
      deallocate(pmap)
   enddo
   call merged%adopt_loops(u, v, sizes, normal_axis=axis, closed=closed, status=st, patches=patches, patch_names=names)
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(overlap) .or. present(contacts)) then
      call overlap_2d(merged, offset, segment_body, u, v, ov, st)
      if (present(overlap)) overlap = ov
      if (present(contacts)) contacts = st
   endif
   endsubroutine merge_posed_curves

   subroutine overlap_3d(surfaces, merged, offset, ov, contacts)
   !< The overlap check of the module header on the merged surface (surface s = facets offset(s) + 1 .. offset(s + 1)).
   type(surface_stl_object), intent(in)  :: surfaces(:) !< Bodies (for their shells).
   type(surface_stl_object), intent(in)  :: merged      !< Merged posed surface.
   integer(I4P),             intent(in)  :: offset(:)   !< Facet offsets [ns + 1].
   logical,                  intent(out) :: ov          !< Overlap.
   integer(I4P),             intent(out) :: contacts    !< Inter-surface pairs in contact.
   type(intersection_pair_t), allocatable :: pairs(:)   !< Contact pairs.
   type(facet_object), pointer            :: fp         !< Facet.
   real(R8P),    allocatable :: tri(:,:), bmin(:,:), bmax(:,:) !< Placed triangles [9, nf], surface boxes [3, ns].
   integer(I4P), allocatable :: owner(:)                !< Surface of each facet.
   logical,      allocatable :: meet(:,:)               !< Boxes of two surfaces meet.
   integer(I4P)              :: ns, nf, s, t, f, m, k, st !< Work.

   ov = .false. ; contacts = 0_I4P
   ns = size(surfaces, kind=I4P) ; nf = offset(ns + 1)
   if (ns < 2_I4P) return
   allocate(tri(9, nf), owner(nf), bmin(3, ns), bmax(3, ns), meet(ns, ns))
   do s = 1, ns
      bmin(:, s) = huge(1._R8P) ; bmax(:, s) = -huge(1._R8P)
      do f = offset(s) + 1, offset(s + 1)
         owner(f) = s
         fp => merged%facet_at(f)
         do m = 1, 3
            tri(3 * m - 2:3 * m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
            bmin(:, s) = min(bmin(:, s), tri(3 * m - 2:3 * m, f)) ; bmax(:, s) = max(bmax(:, s), tri(3 * m - 2:3 * m, f))
         enddo
      enddo
   enddo
   meet = .false.
   do s = 1, ns
      do t = 1, ns
         if (s == t .or. offset(s + 1) == offset(s) .or. offset(t + 1) == offset(t)) cycle
         meet(s, t) = all(bmin(:, s) <= bmax(:, t)) .and. all(bmin(:, t) <= bmax(:, s))
      enddo
   enddo
   if (.not. any(meet)) return
   ! contacts: the inter-surface pairs, exact
   call merged%find_self_intersections(pairs, status=st, include_touching=.true.)
   do k = 1, size(pairs)
      if (owner(pairs(k)%a) == owner(pairs(k)%b)) cycle
      contacts = contacts + 1_I4P
      if (ov) cycle
      if (pairs(k)%kind == SELFINT_SHEETS) then
         ov = .true.
      elseif (pierces(tri(:, pairs(k)%a), tri(:, pairs(k)%b)) .or. pierces(tri(:, pairs(k)%b), tri(:, pairs(k)%a))) then
         ov = .true.
      elseif (pairs(k)%kind == SELFINT_COPLANAR .or. pairs(k)%kind == SELFINT_EDGE_FOLD) then
         ov = dot_product(tri_normal(tri(:, pairs(k)%a)), tri_normal(tri(:, pairs(k)%b))) > 0._R8P
      endif
   enddo
   if (ov) return
   ! crossings: an exact narrow phase of our own over the inter-surface pairs of overlapping boxes (sweep along x), as
   ! `find_self_intersections` misses pairs in degenerate configurations (two unit cubes offset by exactly 0.5 along y
   ! give no pair at all)
   ov = sweep_overlap(tri, owner)
   if (ov) return
   ! containment: one vertex per shell, off the other surface, by the exact parity
   do s = 1, ns
      do t = 1, ns
         if (.not. meet(s, t)) cycle
         if (shell_inside(surfaces(s), offset(s), tri, offset(t), offset(t + 1), bmin(:, t), bmax(:, t))) then
            ov = .true. ; return
         endif
      enddo
   enddo
   endsubroutine overlap_3d

   function sweep_overlap(tri, owner) result(ov)
   !< True if two triangles of different owners cross (`pierces`) or overlap coplanar with the same orientation
   !< (`coplanar_same_side`): the candidates are the pairs of closed boxes that meet, by a sweep along x.
   real(R8P),    intent(in) :: tri(:,:)  !< Triangles [9, nf].
   integer(I4P), intent(in) :: owner(:)  !< Owner of each triangle.
   logical                  :: ov        !< Overlap.
   real(R8P),    allocatable :: lo(:,:), hi(:,:) !< Triangle boxes [3, nf].
   integer(I4P), allocatable :: order(:), active(:) !< Triangles by lo x, active list.
   integer(I4P)              :: nf, i, a, b, k, na, m !< Work.

   ov = .false.
   nf = size(tri, 2, kind=I4P)
   allocate(lo(3, nf), hi(3, nf), order(nf), active(nf))
   do i = 1, nf
      lo(:, i) = min(tri(1:3, i), tri(4:6, i), tri(7:9, i)) ; hi(:, i) = max(tri(1:3, i), tri(4:6, i), tri(7:9, i))
      order(i) = i
   enddo
   call sort_by_key(lo(1, :), order)
   na = 0_I4P
   do k = 1, nf
      a = order(k)
      m = 0_I4P
      do i = 1, na ! drop the triangles that end before this one starts, test the others
         b = active(i)
         if (hi(1, b) < lo(1, a)) cycle
         m = m + 1_I4P ; active(m) = b
         if (owner(a) == owner(b)) cycle
         if (any(hi(2:3, b) < lo(2:3, a)) .or. any(hi(2:3, a) < lo(2:3, b))) cycle
         if (pierces(tri(:, a), tri(:, b)) .or. pierces(tri(:, b), tri(:, a)) .or. coplanar_same_side(tri(:, a), tri(:, b))) then
            ov = .true. ; return
         endif
      enddo
      na = m + 1_I4P ; active(na) = a
   enddo
   endfunction sweep_overlap

   subroutine sort_by_key(key, idx)
   !< Stable merge sort of idx by key(idx) (bottom-up).
   real(R8P),    intent(in)    :: key(:) !< Keys.
   integer(I4P), intent(inout) :: idx(:) !< Indexes.
   integer(I4P), allocatable   :: tmp(:) !< Buffer.
   integer(I4P)                :: n, w, l, mid, r, i, j, k !< Work.

   n = size(idx, kind=I4P)
   allocate(tmp(n))
   w = 1_I4P
   do while (w < n)
      l = 1_I4P
      do while (l <= n)
         mid = min(l + w, n + 1_I4P) ; r = min(l + 2_I4P * w, n + 1_I4P)
         i = l ; j = mid ; k = l
         do while (i < mid .and. j < r)
            if (key(idx(j)) < key(idx(i))) then
               tmp(k) = idx(j) ; j = j + 1_I4P
            else
               tmp(k) = idx(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i < mid)
            tmp(k) = idx(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j < r)
            tmp(k) = idx(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         l = r
      enddo
      idx = tmp
      w = 2_I4P * w
   enddo
   endsubroutine sort_by_key

   pure function coplanar_same_side(ta, tb) result(yes)
   !< True if the triangles are coplanar (exact), have the same orientation and their interiors overlap (exact: no edge
   !< line of either weakly separates them in the projection dropping the largest normal component of ta).
   real(R8P), intent(in) :: ta(9), tb(9) !< Triangles.
   logical               :: yes          !< Coplanar overlap, same orientation.
   real(R8P)             :: na(3), nb(3), sa, sb !< Normals, projected orientations.
   integer(I4P)          :: k, i1, i2   !< Projection.

   yes = .false.
   if (orient3d_xyz(ta(1:3), ta(4:6), ta(7:9), tb(1:3)) /= 0._R8P) return
   if (orient3d_xyz(ta(1:3), ta(4:6), ta(7:9), tb(4:6)) /= 0._R8P) return
   if (orient3d_xyz(ta(1:3), ta(4:6), ta(7:9), tb(7:9)) /= 0._R8P) return
   na = tri_normal(ta) ; nb = tri_normal(tb)
   k = maxloc(abs(na), dim=1)
   if (na(k) == 0._R8P) return
   i1 = mod(k, 3_I4P) + 1_I4P ; i2 = mod(k + 1_I4P, 3_I4P) + 1_I4P
   sa = orient2d_xy(ta(i1), ta(i2), ta(3 + i1), ta(3 + i2), ta(6 + i1), ta(6 + i2))
   sb = orient2d_xy(tb(i1), tb(i2), tb(3 + i1), tb(3 + i2), tb(6 + i1), tb(6 + i2))
   if (sa == 0._R8P .or. sb == 0._R8P) return
   if ((sa > 0._R8P) .neqv. (sb > 0._R8P)) return ! opposite orientation: the solids on the two sides (touching)
   if (separates(ta, tb, sa) .or. separates(tb, ta, sb)) return
   yes = .true.
   contains
      pure function separates(t, o, st) result(sep)
      !< True if the line of an edge of t has every vertex of o on its closed outer side.
      real(R8P), intent(in) :: t(9), o(9), st !< Triangle, other, orientation of t.
      logical               :: sep            !< Weakly separated.
      real(R8P)             :: d              !< Orientation.
      integer(I4P)          :: e, i, j, v     !< Work.

      sep = .false.
      do e = 1, 3
         i = 3 * (e - 1) ; j = 3 * mod(e, 3_I4P)
         sep = .true.
         do v = 0, 2
            d = orient2d_xy(t(i + i1), t(i + i2), t(j + i1), t(j + i2), o(3 * v + i1), o(3 * v + i2))
            if ((st > 0._R8P .and. d > 0._R8P) .or. (st < 0._R8P .and. d < 0._R8P)) then
               sep = .false. ; exit
            endif
         enddo
         if (sep) return
      enddo
      endfunction separates
   endfunction coplanar_same_side

   function shell_inside(surface, off, tri, t0, t1, tmin, tmax) result(inside)
   !< True if a shell of `surface` (its facets at off + 1 .. in `tri`) lies inside the body of facets t0 + 1 .. t1, tested
   !< at one vertex per shell that lies on none of those facets; a shell with every vertex on them counts as inside.
   type(surface_stl_object), intent(in) :: surface     !< Body whose shells are tested (body frame: shells only).
   integer(I4P),             intent(in) :: off         !< Its first facet - 1 in `tri`.
   real(R8P),                intent(in) :: tri(:,:)    !< Placed triangles [9, nf].
   integer(I4P),             intent(in) :: t0, t1      !< The other body: facets t0 + 1 .. t1.
   real(R8P),                intent(in) :: tmin(3), tmax(3) !< Its box.
   logical                              :: inside      !< Some shell inside (or all on) the other body.
   logical, allocatable                 :: decided(:)  !< Shell decided.
   integer(I4P)                         :: f, m, sh, nsh !< Work.
   real(R8P)                            :: p(3)        !< Vertex.

   inside = .false.
   nsh = surface%get_shells_number()
   if (nsh < 1_I4P) return
   allocate(decided(nsh)) ; decided = .false.
   do f = 1, surface%get_facets_number()
      sh = surface%facet_shell(f)
      if (sh < 1_I4P .or. sh > nsh) cycle
      if (decided(sh)) cycle
      do m = 1, 3
         p = tri(3 * m - 2:3 * m, off + f)
         if (any(p < tmin) .or. any(p > tmax)) then ! outside the box: off the body and outside it
            decided(sh) = .true. ; exit
         endif
         if (on_triangles(p, tri, t0, t1)) cycle
         decided(sh) = .true.
         if (parity_3d(p, tri, t0, t1)) then
            inside = .true. ; return
         endif
         exit
      enddo
   enddo
   inside = .not. all(decided)
   endfunction shell_inside

   pure function on_triangles(p, tri, t0, t1) result(on)
   !< True if the point p lies on one of the closed triangles t0 + 1 .. t1 (exact).
   real(R8P),    intent(in) :: p(3)     !< Point.
   real(R8P),    intent(in) :: tri(:,:) !< Triangles [9, nf].
   integer(I4P), intent(in) :: t0, t1   !< Range.
   logical                  :: on       !< On a triangle.
   real(R8P)                :: n(3), o1, o2, o3 !< Normal, orientations.
   integer(I4P)             :: f, k, i1, i2 !< Work.

   on = .false.
   do f = t0 + 1, t1
      if (p(1) < min(tri(1, f), tri(4, f), tri(7, f)) .or. p(1) > max(tri(1, f), tri(4, f), tri(7, f))) cycle
      if (p(2) < min(tri(2, f), tri(5, f), tri(8, f)) .or. p(2) > max(tri(2, f), tri(5, f), tri(8, f))) cycle
      if (p(3) < min(tri(3, f), tri(6, f), tri(9, f)) .or. p(3) > max(tri(3, f), tri(6, f), tri(9, f))) cycle
      if (orient3d_xyz(tri(1:3, f), tri(4:6, f), tri(7:9, f), p) /= 0._R8P) cycle
      n = tri_normal(tri(:, f))
      k = maxloc(abs(n), dim=1)
      if (n(k) == 0._R8P) cycle ! degenerate facet
      i1 = mod(k, 3_I4P) + 1_I4P ; i2 = mod(k + 1_I4P, 3_I4P) + 1_I4P ! projection dropping axis k
      o1 = orient2d_xy(tri(i1, f), tri(i2, f), tri(3 + i1, f), tri(3 + i2, f), p(i1), p(i2))
      o2 = orient2d_xy(tri(3 + i1, f), tri(3 + i2, f), tri(6 + i1, f), tri(6 + i2, f), p(i1), p(i2))
      o3 = orient2d_xy(tri(6 + i1, f), tri(6 + i2, f), tri(i1, f), tri(i2, f), p(i1), p(i2))
      if ((o1 >= 0._R8P .and. o2 >= 0._R8P .and. o3 >= 0._R8P) .or. (o1 <= 0._R8P .and. o2 <= 0._R8P .and. o3 <= 0._R8P)) then
         on = .true. ; return
      endif
   enddo
   endfunction on_triangles

   pure function parity_3d(p, tri, t0, t1) result(inside)
   !< True if the point p (on none of the triangles) is inside the closed body of triangles t0 + 1 .. t1: an odd number
   !< of their exact crossings with the +z line through p lies below p (the cut's line kernels, global perturbation).
   real(R8P),    intent(in) :: p(3)     !< Point.
   real(R8P),    intent(in) :: tri(:,:) !< Triangles [9, nf].
   integer(I4P), intent(in) :: t0, t1   !< Range.
   logical                  :: inside   !< Inside.
   real(R8P)                :: t, tv(3, 3) !< Crossing, triangle.
   logical                  :: hit, enter  !< Kernel outputs.
   integer(I4P)             :: f, n        !< Work.

   n = 0_I4P
   do f = t0 + 1, t1
      if (p(1) < min(tri(1, f), tri(4, f), tri(7, f)) .or. p(1) > max(tri(1, f), tri(4, f), tri(7, f))) cycle
      if (p(2) < min(tri(2, f), tri(5, f), tri(8, f)) .or. p(2) > max(tri(2, f), tri(5, f), tri(8, f))) cycle
      call tri_line_hit_oac(tri(:, f), 3_I4P, p(1), p(2), hit, t, enter)
      if (.not. hit) cycle
      tv = reshape(tri(:, f), [3, 3])
      if (line_crossing_before(tv, 3_I4P, p(1), p(2), p(3), enter)) n = n + 1_I4P
   enddo
   inside = mod(n, 2_I4P) == 1_I4P
   endfunction parity_3d

   pure function pierces(ta, tb) result(yes)
   !< True if an edge of triangle ta passes through the interior of triangle tb transversally (exact): its ends strictly
   !< on the two sides of tb's plane, and the line of the edge strictly inside tb (the three orient3d of the edge with
   !< tb's edges of one strict sign).
   real(R8P), intent(in) :: ta(9), tb(9) !< Triangles.
   logical               :: yes          !< Transversal piercing.
   real(R8P)             :: s1, s2, o1, o2, o3 !< Orientations.
   integer(I4P)          :: e, i, j      !< Work.

   yes = .false.
   do e = 1, 3
      i = 3 * (e - 1) + 1 ; j = 3 * mod(e, 3_I4P) + 1
      s1 = orient3d_xyz(tb(1:3), tb(4:6), tb(7:9), ta(i:i + 2))
      s2 = orient3d_xyz(tb(1:3), tb(4:6), tb(7:9), ta(j:j + 2))
      if (.not. ((s1 > 0._R8P .and. s2 < 0._R8P) .or. (s1 < 0._R8P .and. s2 > 0._R8P))) cycle
      o1 = orient3d_xyz(ta(i:i + 2), ta(j:j + 2), tb(1:3), tb(4:6))
      o2 = orient3d_xyz(ta(i:i + 2), ta(j:j + 2), tb(4:6), tb(7:9))
      o3 = orient3d_xyz(ta(i:i + 2), ta(j:j + 2), tb(7:9), tb(1:3))
      if ((o1 > 0._R8P .and. o2 > 0._R8P .and. o3 > 0._R8P) .or. (o1 < 0._R8P .and. o2 < 0._R8P .and. o3 < 0._R8P)) then
         yes = .true. ; return
      endif
   enddo
   endfunction pierces

   pure function tri_normal(t) result(n)
   !< (B - A) x (C - A) of the triangle t = [A, B, C].
   real(R8P), intent(in) :: t(9) !< Triangle.
   real(R8P)             :: n(3) !< Normal (not unit).
   real(R8P)             :: e1(3), e2(3) !< Edges.

   e1 = t(4:6) - t(1:3) ; e2 = t(7:9) - t(1:3)
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   endfunction tri_normal

   subroutine overlap_2d(merged, offset, owner, u, v, ov, contacts)
   !< The overlap check of the module header on the merged curve (curve c = vertices offset(c) + 1 .. offset(c + 1)).
   type(curve_object), intent(in)  :: merged    !< Merged posed curve.
   integer(I4P),       intent(in)  :: offset(:) !< Vertex offsets [nc + 1].
   integer(I4P),       intent(in)  :: owner(:)  !< Curve of each vertex.
   real(R8P),          intent(in)  :: u(:), v(:) !< Placed vertices.
   logical,            intent(out) :: ov        !< Overlap.
   integer(I4P),       intent(out) :: contacts  !< Inter-curve segment pairs in contact.
   integer(I4P), allocatable :: nxt(:), pairs(:,:) !< Next vertex, contact pairs.
   real(R8P),    allocatable :: bmin(:,:), bmax(:,:) !< Curve boxes [2, nc].
   logical,      allocatable :: meet(:,:)        !< Boxes of two curves meet.
   integer(I4P)              :: nc, nv, c, d, i, k, n, a, b, first, last, l !< Work.
   real(R8P)                 :: o1, o2, o3, o4  !< Orientations.
   logical                   :: decided, cross  !< Work.

   ov = .false. ; contacts = 0_I4P
   nc = size(offset, kind=I4P) - 1_I4P ; nv = offset(nc + 1)
   if (nc < 2_I4P) return
   allocate(nxt(nv), bmin(2, nc), bmax(2, nc), meet(nc, nc))
   do i = 1, nv
      nxt(i) = merged%get_next(i)
   enddo
   do c = 1, nc
      bmin(:, c) = huge(1._R8P) ; bmax(:, c) = -huge(1._R8P)
      do i = offset(c) + 1, offset(c + 1)
         bmin(:, c) = min(bmin(:, c), [u(i), v(i)]) ; bmax(:, c) = max(bmax(:, c), [u(i), v(i)])
      enddo
   enddo
   meet = .false.
   do c = 1, nc
      do d = 1, nc
         if (c == d .or. offset(c + 1) == offset(c) .or. offset(d + 1) == offset(d)) cycle
         meet(c, d) = all(bmin(:, c) <= bmax(:, d)) .and. all(bmin(:, d) <= bmax(:, c))
      enddo
   enddo
   if (.not. any(meet)) return
   ! contacts: the inter-curve pairs, exact
   call merged%self_intersections(pairs, n)
   do k = 1, n
      a = pairs(1, k) ; b = pairs(2, k)
      if (owner(a) == owner(b)) cycle
      contacts = contacts + 1_I4P
      if (ov .or. nxt(a) == 0_I4P .or. nxt(b) == 0_I4P) cycle
      o1 = orient2d_xy(u(a), v(a), u(nxt(a)), v(nxt(a)), u(b), v(b))
      o2 = orient2d_xy(u(a), v(a), u(nxt(a)), v(nxt(a)), u(nxt(b)), v(nxt(b)))
      o3 = orient2d_xy(u(b), v(b), u(nxt(b)), v(nxt(b)), u(a), v(a))
      o4 = orient2d_xy(u(b), v(b), u(nxt(b)), v(nxt(b)), u(nxt(a)), v(nxt(a)))
      cross = ((o1 > 0._R8P .and. o2 < 0._R8P) .or. (o1 < 0._R8P .and. o2 > 0._R8P)) .and. &
              ((o3 > 0._R8P .and. o4 < 0._R8P) .or. (o3 < 0._R8P .and. o4 > 0._R8P))
      if (cross) then
         ov = .true.
      elseif (o1 == 0._R8P .and. o2 == 0._R8P) then ! collinear: overlapping over a positive length in the same direction
         ov = collinear_overlap(u(a), v(a), u(nxt(a)), v(nxt(a)), u(b), v(b), u(nxt(b)), v(nxt(b)))
      endif
   enddo
   if (ov) return
   ! containment: one vertex per closed loop, off the other curve, by the exact parity
   do l = 1, merged%get_loops_number()
      if (.not. merged%is_loop_closed(l)) cycle
      call merged%get_loop_range(l, first, last)
      if (last < first) cycle
      c = owner(first)
      do d = 1, nc
         if (.not. meet(c, d)) cycle
         decided = .false.
         do i = first, last
            if (u(i) < bmin(1, d) .or. u(i) > bmax(1, d) .or. v(i) < bmin(2, d) .or. v(i) > bmax(2, d)) then
               decided = .true. ; exit
            endif
            if (on_segments(u(i), v(i), u, v, nxt, offset(d), offset(d + 1))) cycle
            decided = .true.
            if (parity_2d(u(i), v(i), u, v, nxt, offset(d), offset(d + 1))) then
               ov = .true. ; return
            endif
            exit
         enddo
         if (.not. decided) then ! every vertex on the other curve: coincident (conservative)
            ov = .true. ; return
         endif
      enddo
   enddo
   endsubroutine overlap_2d

   pure function collinear_overlap(au, av, bu, bv, cu, cv, du, dv) result(yes)
   !< True if the collinear segments a -> b and c -> d overlap over a positive length with the same direction (exact:
   !< compared on the coordinate along which a -> b is not constant).
   real(R8P), intent(in) :: au, av, bu, bv, cu, cv, du, dv !< Segments.
   logical               :: yes                            !< Same-direction overlap.
   real(R8P)             :: a0, a1, c0, c1                 !< Coordinates along the line.

   if (au /= bu) then
      a0 = au ; a1 = bu ; c0 = cu ; c1 = du
   else
      a0 = av ; a1 = bv ; c0 = cv ; c1 = dv
   endif
   yes = max(min(a0, a1), min(c0, c1)) < min(max(a0, a1), max(c0, c1)) .and. ((a1 > a0) .eqv. (c1 > c0))
   endfunction collinear_overlap

   pure function on_segments(pu, pv, u, v, nxt, i0, i1) result(on)
   !< True if (pu, pv) lies on one of the closed segments i -> nxt(i), i = i0 + 1 .. i1 (exact).
   real(R8P),    intent(in) :: pu, pv     !< Point.
   real(R8P),    intent(in) :: u(:), v(:) !< Vertices.
   integer(I4P), intent(in) :: nxt(:)     !< Next vertex (0: none).
   integer(I4P), intent(in) :: i0, i1     !< Range.
   logical                  :: on         !< On a segment.
   integer(I4P)             :: i, j       !< Work.

   on = .false.
   do i = i0 + 1, i1
      j = nxt(i)
      if (j == 0_I4P) cycle
      if (pu < min(u(i), u(j)) .or. pu > max(u(i), u(j)) .or. pv < min(v(i), v(j)) .or. pv > max(v(i), v(j))) cycle
      if (orient2d_xy(u(i), v(i), u(j), v(j), pu, pv) == 0._R8P) then
         on = .true. ; return
      endif
   enddo
   endfunction on_segments

   pure function parity_2d(pu, pv, u, v, nxt, i0, i1) result(inside)
   !< True if (pu, pv) (on none of the segments) is inside the closed loops of segments i0 + 1 .. i1: an odd number of
   !< their exact crossings with the +u line through it lies before it (the cut's 2-D line kernels).
   real(R8P),    intent(in) :: pu, pv     !< Point.
   real(R8P),    intent(in) :: u(:), v(:) !< Vertices.
   integer(I4P), intent(in) :: nxt(:)     !< Next vertex (0: none).
   integer(I4P), intent(in) :: i0, i1     !< Range.
   logical                  :: inside     !< Inside.
   real(R8P)                :: t          !< Crossing.
   logical                  :: hit, enter !< Kernel outputs.
   integer(I4P)             :: i, j, n    !< Work.

   n = 0_I4P
   do i = i0 + 1, i1
      j = nxt(i)
      if (j == 0_I4P) cycle
      call seg_line_hit_oac(u(i), v(i), u(j), v(j), 1_I4P, pv, hit, t, enter)
      if (.not. hit) cycle
      if (line_crossing_before_2d(u(i), v(i), u(j), v(j), 1_I4P, pv, pu, enter)) n = n + 1_I4P
   enddo
   inside = mod(n, 2_I4P) == 1_I4P
   endfunction parity_2d

   subroutine label_piece_bodies_3d(walls, facet_body, facet_sub_body, np, piece_bodies, status)
   !< The bodies of the wall terms of a 3-D cut of the merged surface and of its np pieces (`piece_bodies_t`), from the
   !< map of `merge_posed_bodies`. The terms must be on the host. `status = STATUS_INVALID_INPUT` (that term labelled 0)
   !< if a term's facet is outside the map.
   type(cut_walls_3d_t), intent(in)            :: walls             !< Wall terms (host copies).
   integer(I4P),         intent(in)            :: facet_body(:)     !< Body of each facet.
   integer(I4P),         intent(in)            :: facet_sub_body(:) !< Its sub-body.
   integer(I4P),         intent(in)            :: np                !< Pieces of the cut.
   type(piece_bodies_t), intent(out)           :: piece_bodies      !< Bodies.
   integer(I4P),         intent(out), optional :: status            !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                                :: st                !< Status.

   if (walls%n > 0_I4P) then
      call label_core(walls%n, walls%piece, walls%facet, facet_body, facet_sub_body, np, piece_bodies, st)
   else
      block
         integer(I4P), allocatable :: none(:) !< No terms.
         allocate(none(0))
         call label_core(0_I4P, none, none, facet_body, facet_sub_body, np, piece_bodies, st)
      end block
   endif
   if (present(status)) status = st
   endsubroutine label_piece_bodies_3d

   subroutine label_piece_bodies_2d(walls, segment_body, segment_sub_body, np, piece_bodies, status)
   !< The bodies of the wall segments of a 2-D cut of the merged curve (`cut_walls_2d`) and of its np pieces, from the
   !< map of `merge_posed_bodies`. `status = STATUS_INVALID_INPUT` (that term labelled 0) if a segment is outside it.
   type(cut_walls_2d_t), intent(in)            :: walls               !< Wall segments.
   integer(I4P),         intent(in)            :: segment_body(:)     !< Body of each segment.
   integer(I4P),         intent(in)            :: segment_sub_body(:) !< Its sub-body.
   integer(I4P),         intent(in)            :: np                  !< Pieces of the cut.
   type(piece_bodies_t), intent(out)           :: piece_bodies        !< Bodies.
   integer(I4P),         intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                                :: st                  !< Status.

   if (walls%n > 0_I4P) then
      call label_core(walls%n, walls%piece, walls%segment, segment_body, segment_sub_body, np, piece_bodies, st)
   else
      block
         integer(I4P), allocatable :: none(:) !< No terms.
         allocate(none(0))
         call label_core(0_I4P, none, none, segment_body, segment_sub_body, np, piece_bodies, st)
      end block
   endif
   if (present(status)) status = st
   endsubroutine label_piece_bodies_2d

   subroutine label_core(nw, piece, src, src_body, src_sub, np, pb, st)
   !< Labels of nw wall terms (piece, source element) and the sorted distinct bodies of each of np pieces.
   integer(I4P),         intent(in)  :: nw          !< Wall terms.
   integer(I4P),         intent(in)  :: piece(:)    !< Piece of each term (0: none).
   integer(I4P),         intent(in)  :: src(:)      !< Source facet / segment of each term.
   integer(I4P),         intent(in)  :: src_body(:) !< Body of each source.
   integer(I4P),         intent(in)  :: src_sub(:)  !< Sub-body of each source.
   integer(I4P),         intent(in)  :: np          !< Pieces.
   type(piece_bodies_t), intent(out) :: pb          !< Labels.
   integer(I4P),         intent(out) :: st          !< Status.
   integer(I4P), allocatable :: first(:), order(:), kb(:), ks(:) !< Terms by piece, entries.
   integer(I4P) :: w, p, k, m, n, b, s, q, j !< Work.
   logical      :: found                     !< Entry known.

   st = STATUS_OK
   pb%np = max(0_I4P, np) ; pb%nw = max(0_I4P, nw)
   allocate(pb%wall_body(pb%nw), pb%wall_sub_body(pb%nw), pb%start(pb%np + 1))
   do w = 1, pb%nw
      pb%wall_body(w) = 0_I4P ; pb%wall_sub_body(w) = 0_I4P
      if (src(w) >= 1_I4P .and. src(w) <= size(src_body, kind=I4P) .and. src(w) <= size(src_sub, kind=I4P)) then
         pb%wall_body(w) = src_body(src(w)) ; pb%wall_sub_body(w) = src_sub(src(w))
      else
         st = STATUS_INVALID_INPUT
      endif
   enddo
   ! terms by piece (counting sort, stable)
   allocate(first(pb%np + 2)) ; first = 0_I4P
   do w = 1, pb%nw
      p = piece(w)
      if (p >= 1_I4P .and. p <= pb%np .and. pb%wall_body(w) > 0_I4P) first(p + 2) = first(p + 2) + 1_I4P
   enddo
   first(1) = 1_I4P ; first(2) = 1_I4P
   do p = 2, pb%np + 1
      first(p + 1) = first(p + 1) + first(p)
   enddo
   allocate(order(first(pb%np + 2) - 1_I4P))
   do w = 1, pb%nw
      p = piece(w)
      if (p >= 1_I4P .and. p <= pb%np .and. pb%wall_body(w) > 0_I4P) then
         order(first(p + 1)) = w ; first(p + 1) = first(p + 1) + 1_I4P
      endif
   enddo
   ! now first(p) .. first(p + 1) - 1 are the terms of piece p; the distinct (body, sub_body) of each, sorted
   allocate(kb(size(order)), ks(size(order)))
   n = 0_I4P
   do p = 1, pb%np
      pb%start(p) = n + 1_I4P
      m = n
      do k = first(p), first(p + 1) - 1_I4P
         b = pb%wall_body(order(k)) ; s = pb%wall_sub_body(order(k))
         found = .false.
         do q = n + 1, m
            if (kb(q) == b .and. ks(q) == s) then
               found = .true. ; exit
            endif
         enddo
         if (found) cycle
         ! insertion in (body, sub_body) order
         q = m
         do while (q > n)
            if (kb(q) < b .or. (kb(q) == b .and. ks(q) < s)) exit
            q = q - 1_I4P
         enddo
         do j = m, q + 1, -1
            kb(j + 1) = kb(j) ; ks(j + 1) = ks(j)
         enddo
         kb(q + 1) = b ; ks(q + 1) = s ; m = m + 1_I4P
      enddo
      n = m
   enddo
   pb%start(pb%np + 1) = n + 1_I4P
   allocate(pb%body(n), pb%sub_body(n))
   pb%body = kb(1:n) ; pb%sub_body = ks(1:n)
   endsubroutine label_core

   subroutine cut_cells_3d_bodies(surfaces, poses, x, y, z, cut, status, links, walls, edges, moments_order, quadrature, &
                                  device, resident, piece_bodies, merged, facet_body, facet_sub_body, contacts)
   !< The cut cells of several posed bodies on the block x(nx+1), y(ny+1), z(nz+1) (see the module header): merge
   !< (`merge_posed_bodies`), check (an overlap is refused), cut of the merged surface (`cut_cells_3d`, or
   !< `cut_cells_3d_device` with `device`), bit for bit the cut of the merged surface, and with `piece_bodies` the bodies
   !< of the wall terms and of the pieces (the wall terms are then computed, in the caller's `walls` if present).
   !< `links`, `walls`, `edges`, `moments_order`, `quadrature`, `resident` as in the cut (moments and quadrature are
   !< refused with `device`). `merged` returns the merged surface (the cut's surface: `cut_walls_gyration` and
   !< `cut_loads_3d` need it), `facet_body` / `facet_sub_body` its map, `contacts` the inter-surface contacts.
   !< `status = STATUS_INVALID_INPUT` (stderr message, the cut reset) on an overlap, bad poses, bad edges or a refused
   !< option.
   type(surface_stl_object),  intent(in)              :: surfaces(:)       !< Bodies (body frame).
   type(rigid_pose_t),        intent(in)              :: poses(:)          !< Pose of each surface.
   real(R8P),                 intent(in)              :: x(:), y(:), z(:)  !< Cell edges (strictly increasing).
   type(cut_cells_3d_t),      intent(inout)           :: cut               !< Cut cells (released, reset).
   integer(I4P),              intent(out),   optional :: status            !< STATUS_OK or STATUS_INVALID_INPUT.
   type(cut_links_t),         intent(inout), optional :: links             !< Links between the pieces.
   type(cut_walls_3d_t),      intent(inout), optional :: walls             !< The wall terms, labelled by piece and facet.
   logical,                   intent(in),    optional :: edges             !< Edge fractions too.
   integer(I4P),              intent(in),    optional :: moments_order     !< Moments of the pieces (host only).
   logical,                   intent(in),    optional :: quadrature        !< Quadrature rules (host only).
   logical,                   intent(in),    optional :: device            !< Cut on the device (default .false.).
   logical,                   intent(in),    optional :: resident          !< Device outputs left resident.
   type(piece_bodies_t),      intent(out),   optional :: piece_bodies      !< Bodies of the wall terms and pieces.
   type(surface_stl_object),  intent(inout), optional :: merged            !< Merged posed surface.
   integer(I4P), allocatable, intent(out),   optional :: facet_body(:)     !< Body of each merged facet.
   integer(I4P), allocatable, intent(out),   optional :: facet_sub_body(:) !< Its sub-body.
   integer(I4P),              intent(out),   optional :: contacts          !< Inter-surface facet pairs in contact.
   type(surface_stl_object)  :: lm                                        !< Merged surface (local).
   type(cut_walls_3d_t)      :: lw                                        !< Wall terms (local).
   integer(I4P), allocatable :: fb(:), fs(:)                              !< Facet map.
   integer(I4P)              :: st, nc                                    !< Status, contacts.
   logical                   :: dev, ov                                   !< Options.

   call cut_cells_3d_reset(cut)
   if (present(links)) call cut_links_reset(links)
   if (present(walls)) call cut_walls_reset(walls)
   if (present(status)) status = STATUS_OK
   if (present(contacts)) contacts = 0_I4P
   dev = .false. ; if (present(device)) dev = device
   if (dev) then
      st = STATUS_OK
      if (present(moments_order)) then
         if (moments_order > 0_I4P) st = STATUS_INVALID_INPUT
      endif
      if (present(quadrature)) then
         if (quadrature) st = STATUS_INVALID_INPUT
      endif
      if (st /= STATUS_OK) then
         write(stderr, '(A)') 'fossil cut_cells_3d_bodies: refused -- moments and quadrature are host-only (device)'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   if (present(merged)) then
      call merge_cut(merged)
   else
      call merge_cut(lm)
   endif
   if (present(facet_body)) call move_alloc(from=fb, to=facet_body)
   if (present(facet_sub_body)) call move_alloc(from=fs, to=facet_sub_body)
   contains
      subroutine merge_cut(m)
      !< Merge into m, check, cut, label.
      type(surface_stl_object), intent(inout) :: m !< Merged surface.

      call merge_posed_surfaces(surfaces, poses, m, fb, fs, status=st, overlap=ov, contacts=nc)
      if (present(contacts)) contacts = nc
      if (st /= STATUS_OK) then
         write(stderr, '(A)') 'fossil cut_cells_3d_bodies: refused -- no surface, or poses not matching the surfaces'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      if (ov) then
         write(stderr, '(A)') 'fossil cut_cells_3d_bodies: refused -- two bodies overlap (their volumes intersect)'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      if (present(walls)) then
         call cut_one(m, walls)
      elseif (present(piece_bodies)) then
         call cut_one(m, lw)
         if (lw%on_device) call lw%exit_device
      else
         call cut_one(m)
      endif
      endsubroutine merge_cut

      subroutine cut_one(m, w)
      !< The cut of m (host or device), then the labels.
      type(surface_stl_object), intent(in)              :: m !< Merged surface.
      type(cut_walls_3d_t),     intent(inout), optional :: w !< Wall terms.

      if (dev) then
         call cut_cells_3d_device(m, x, y, z, cut, status=st, edges=edges, resident=resident, links=links, walls=w)
      else
         call cut_cells_3d(m, x, y, z, cut, status=st, links=links, edges=edges, moments_order=moments_order, &
                           quadrature=quadrature, walls=w)
      endif
      if (present(status)) status = st
      if (st /= STATUS_OK .or. .not. present(piece_bodies) .or. .not. present(w)) return
      if (w%on_device) call w%update_host
      call label_piece_bodies_3d(w, fb, fs, cut%np, piece_bodies)
      endsubroutine cut_one
   endsubroutine cut_cells_3d_bodies

   subroutine cut_cells_2d_bodies(curves, poses, x, y, cut, status, intervals, moments_order, boundaries, quadrature, &
                                  device, resident, piece_bodies, walls, merged, segment_body, segment_sub_body, contacts)
   !< The 2-D twin of `cut_cells_3d_bodies` on the block x(nx+1), y(ny+1) (u, v): merge, check, cut of the merged curve
   !< (`curve%cut_cells`, or `curve%cut_cells_device` with `device`), bit for bit, and with `piece_bodies` (or `walls`,
   !< the wall segments of `cut_walls_2d`) the bodies of the wall segments and of the pieces. On the host the labels
   !< need the oriented boundaries: `piece_bodies` or `walls` turn `boundaries` on. On the device the cut has none: the
   !< labels come from a host cut with boundaries (the same pieces). `intervals`, `moments_order`, `boundaries`,
   !< `quadrature`, `resident` as in the cut (moments, boundaries and quadrature refused with `device`). `merged`,
   !< `segment_body`, `segment_sub_body`, `contacts` as in 3-D. `status = STATUS_INVALID_INPUT` (stderr message, the cut
   !< reset) on an overlap, bad poses or planes, bad edges or a refused option.
   type(curve_object),        intent(in)              :: curves(:)           !< Bodies (body frame).
   type(rigid_pose_2d_t),     intent(in)              :: poses(:)            !< Pose of each curve.
   real(R8P),                 intent(in)              :: x(:), y(:)          !< Cell edges (u, v), strictly increasing.
   type(cut_cells_2d_t),      intent(inout)           :: cut                 !< Cut cells (released, reset).
   integer(I4P),              intent(out),   optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                   intent(in),    optional :: intervals           !< Fluid edge intervals (default .true.).
   integer(I4P),              intent(in),    optional :: moments_order       !< Moments of the pieces (host only).
   logical,                   intent(in),    optional :: boundaries          !< Oriented boundaries (host only).
   logical,                   intent(in),    optional :: quadrature          !< Quadrature rules (host only).
   logical,                   intent(in),    optional :: device              !< Cut on the device (default .false.).
   logical,                   intent(in),    optional :: resident            !< Device outputs left resident.
   type(piece_bodies_t),      intent(out),   optional :: piece_bodies        !< Bodies of the wall segments and pieces.
   type(cut_walls_2d_t),      intent(out),   optional :: walls               !< Wall segments (`cut_walls_2d`).
   type(curve_object),        intent(inout), optional :: merged              !< Merged posed curve.
   integer(I4P), allocatable, intent(out),   optional :: segment_body(:)     !< Body of each merged segment.
   integer(I4P), allocatable, intent(out),   optional :: segment_sub_body(:) !< Its sub-body.
   integer(I4P),              intent(out),   optional :: contacts            !< Inter-curve segment pairs in contact.
   type(curve_object)        :: lm                                          !< Merged curve (local).
   integer(I4P), allocatable :: sb(:), ss(:)                                !< Segment map.
   integer(I4P)              :: st, nc                                      !< Status, contacts.
   logical                   :: dev, ov, wb, lab                            !< Options.

   call cut_cells_2d_reset(cut)
   if (present(status)) status = STATUS_OK
   if (present(contacts)) contacts = 0_I4P
   dev = .false. ; if (present(device)) dev = device
   wb = .false. ; if (present(boundaries)) wb = boundaries
   lab = present(piece_bodies) .or. present(walls)
   if (dev) then
      st = STATUS_OK
      if (present(moments_order)) then
         if (moments_order > 0_I4P) st = STATUS_INVALID_INPUT
      endif
      if (present(quadrature)) then
         if (quadrature) st = STATUS_INVALID_INPUT
      endif
      if (wb) st = STATUS_INVALID_INPUT
      if (st /= STATUS_OK) then
         write(stderr, '(A)') 'fossil cut_cells_2d_bodies: refused -- moments, boundaries and quadrature are host-only '// &
                              '(device)'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   if (present(merged)) then
      call merge_cut(merged)
   else
      call merge_cut(lm)
   endif
   if (present(segment_body)) call move_alloc(from=sb, to=segment_body)
   if (present(segment_sub_body)) call move_alloc(from=ss, to=segment_sub_body)
   contains
      subroutine merge_cut(m)
      !< Merge into m, check, cut, label.
      type(curve_object), intent(inout) :: m !< Merged curve.
      type(cut_cells_2d_t)              :: ch !< Host cut with boundaries (device labels).
      type(cut_walls_2d_t)              :: lw !< Wall segments (local).

      call merge_posed_curves(curves, poses, m, sb, ss, status=st, overlap=ov, contacts=nc)
      if (present(contacts)) contacts = nc
      if (st /= STATUS_OK) then
         write(stderr, '(A)') 'fossil cut_cells_2d_bodies: refused -- no curve, poses not matching the curves, or '// &
                              'curves in different planes'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      if (ov) then
         write(stderr, '(A)') 'fossil cut_cells_2d_bodies: refused -- two bodies overlap (their areas intersect)'
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      if (dev) then
         call m%cut_cells_device(x, y, cut, status=st, intervals=intervals, resident=resident)
      else
         call m%cut_cells(x, y, cut, status=st, intervals=intervals, moments_order=moments_order, boundaries=wb .or. lab, &
                          quadrature=quadrature)
      endif
      if (present(status)) status = st
      if (st /= STATUS_OK .or. .not. lab) return
      if (dev) then
         call m%cut_cells(x, y, ch, status=st, intervals=.false., boundaries=.true.)
         call cut_walls_2d(ch, lw)
      else
         call cut_walls_2d(cut, lw)
      endif
      if (present(piece_bodies)) then
         if (dev) then
            call label_piece_bodies_2d(lw, sb, ss, ch%np, piece_bodies)
         else
            call label_piece_bodies_2d(lw, sb, ss, cut%np, piece_bodies)
         endif
      endif
      if (present(walls)) walls = lw
      endsubroutine merge_cut
   endsubroutine cut_cells_2d_bodies
endmodule fossil_cut_bodies
