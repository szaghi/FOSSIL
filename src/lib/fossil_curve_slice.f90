!< FOSSIL, section of a surface by an axis-aligned plane: STL -> curve (issue #36 step 5).
module fossil_curve_slice
!< FOSSIL, section of a surface by an axis-aligned plane: STL -> curve (issue #36 step 5).
!<
!< `slice_surface` cuts a `surface_stl_object` with the plane x(normal_axis) = level and returns the
!< section as a `curve_object` in that plane, so a 2-D run and a 3-D run can share one body.
!<
!< Exactness, without tolerances:
!< - a vertex is classified by the exact comparison of its coordinate with `level`; a vertex ON
!<   the plane counts as above it (symbolic perturbation), so every facet is crossed by the plane
!<   along exactly two edges or not at all, and a facet lying in the plane is not crossed;
!< - the crossing point of an edge is a canonical construction: from the edge endpoints in
!<   lexicographic order, with the endpoint itself when the parameter is exactly 0 or 1, so the
!<   two facets sharing an edge produce the same bits and the loops close by exact equality;
!< - the direction of each segment comes from the facet winding (from the crossing of the edge
!<   going above -> below to the crossing of the edge going below -> above, along the facet's
!<   vertex order), so with outward-oriented facets the solid is on the left in the (u, v) plane
!<   of `curve_object`: no normal vector is used.
!< Segments are chained into loops by exact endpoint equality; chains that do not close (a surface
!< with holes) become open loops. Zero-length segments (a facet touching the plane at a vertex)
!< are dropped. Optionally, vertices between two collinear segments (exact `orient2d_xy`) are
!< removed.
!<
!< The plane at a face that bounds the solid from below (e.g. z = 0 for the unit cube) gives an
!< empty section: its vertices count as above, so nothing is below the plane. At a face bounding
!< it from above (z = 1) the section is that face's outline. With `on_plane_low` (issue #40: the
!< plane taken at level + eps, as the cut cells do) a vertex on the plane counts as below, and the
!< roles swap: z = 0 gives the face outline, z = 1 an empty section.
use fossil_curve_object,       only : curve_object
use fossil_facet_object,       only : facet_object
use fossil_predicates,         only : orient2d_xy
use fossil_rigid_pose,         only : rigid_pose_t, to_world
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
private
public :: slice_surface
public :: slice_planes_t, slice_planes_prepare, slice_surface_plane, slice_surface_planes

type :: slice_planes_t
   !< Plane sweep over a set of parallel planes (issue #100 step 3c): the facets that can meet each plane.
   !<
   !< Built once per axis by `slice_planes_prepare` from the facet spans [min, max] along the axis; plane b then lists,
   !< in increasing facet index, every facet whose closed span contains its level (a facet with a NaN coordinate along
   !< the axis is listed on every plane). `slice_surface` keeps a facet only if a vertex is strictly on one side and one
   !< on the other or on the plane, so every facet it keeps has the level in its span: `slice_surface_plane` visits a
   !< superset of those facets in the order of `slice_surface` and builds the same segments in the same order.
   integer(I4P)              :: normal_axis = 0_I4P !< Plane normal axis (0: not prepared or invalid).
   logical                   :: on_plane_low = .false. !< A vertex on the plane counts as below.
   real(R8P),    allocatable :: levels(:)           !< Plane positions.
   integer(I4P), allocatable :: first(:)            !< Facets of plane b: facet(first(b):first(b+1)-1).
   integer(I4P), allocatable :: facet(:)            !< Facet ids per plane, increasing.
endtype slice_planes_t

contains
   subroutine slice_surface(surface, normal_axis, level, curve, merge_collinear, open_loops, on_plane_low, status, pose)
   !< Section of `surface` by the plane x(normal_axis) = level, as a curve in that plane.
   !<
   !< `open_loops` returns the number of chains that did not close (0 for a watertight surface).
   !< `status = STATUS_INVALID_INPUT` (curve empty) if normal_axis is not 1-3 or level is not
   !< finite. An empty section is not an error: the curve has no loops.
   !<
   !< With `pose` (issue #94) the section is that of the surface moved by `pose`, without moving it: every
   !< vertex is placed by `to_world(pose, vertex)` as it is read, the operation `surface%move(pose)` applies,
   !< so the curve equals the section of a moved copy bit for bit (the slice reads nothing else that `move`
   !< changes: the facet vertices only, besides the patches).
   type(surface_stl_object), intent(in)            :: surface         !< Surface.
   integer(I4P),             intent(in)            :: normal_axis     !< Plane normal axis (1, 2, 3).
   real(R8P),                intent(in)            :: level           !< Plane position.
   type(curve_object),       intent(inout)         :: curve           !< Section.
   logical,                  intent(in),  optional :: merge_collinear !< Remove collinear vertices (default .false.).
   integer(I4P),             intent(out), optional :: open_loops      !< Chains that did not close.
   logical,                  intent(in),  optional :: on_plane_low    !< A vertex on the plane counts as below (default .false.).
   integer(I4P),             intent(out), optional :: status          !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_t),       intent(in),  optional :: pose            !< Pose the surface is sliced at (default: as is).
   type(facet_object), pointer                     :: f               !< Facet cursor.
   real(R8P),    allocatable                       :: s(:,:)          !< Segments [su, sv, eu, ev; n].
   integer(I4P), allocatable                       :: sp(:)           !< Patch of each segment (#66).
   real(R8P)                                       :: x(3,3), pa(2), pb(2) !< Facet vertices, crossings.
   integer(I4P)                                    :: n, i, nopen     !< Work.
   logical                                         :: hit, merge_, low_ !< Work.

   call curve%destroy
   if (present(status)) status = STATUS_OK
   if (present(open_loops)) open_loops = 0_I4P
   if (normal_axis < 1_I4P .or. normal_axis > 3_I4P .or. .not. abs(level) <= huge(1._R8P)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   merge_ = .false. ; if (present(merge_collinear)) merge_ = merge_collinear
   low_ = .false. ; if (present(on_plane_low)) low_ = on_plane_low
   allocate(s(4, max(1, surface%get_facets_number())), sp(max(1, surface%get_facets_number())))
   n = 0_I4P
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i)
      x(:, 1) = [f%vertex(1)%x, f%vertex(1)%y, f%vertex(1)%z]
      x(:, 2) = [f%vertex(2)%x, f%vertex(2)%y, f%vertex(2)%z]
      x(:, 3) = [f%vertex(3)%x, f%vertex(3)%y, f%vertex(3)%z]
      if (present(pose)) then ! as `move` places them (issue #94)
         x(:, 1) = to_world(pose, x(:, 1)) ; x(:, 2) = to_world(pose, x(:, 2)) ; x(:, 3) = to_world(pose, x(:, 3))
      endif
      call facet_segment(x, normal_axis, level, low_, pa, pb, hit)
      if (.not. hit) cycle
      n = n + 1_I4P
      s(:, n) = [pa, pb]
      sp(n) = f%patch ! each segment takes the patch of the facet it is cut from (issue #66)
   enddo
   call finish_section(surface, normal_axis, s(:, 1:n), sp(1:n), merge_, curve, nopen)
   if (present(open_loops)) open_loops = nopen
   endsubroutine slice_surface

   subroutine slice_planes_prepare(surface, normal_axis, levels, planes, on_plane_low, status)
   !< Plane sweep index of `surface` for the planes x(normal_axis) = levels(b) (issue #100 step 3c), any order of the
   !< levels. Pass it to `slice_surface_plane`, with the same surface unchanged. `status = STATUS_INVALID_INPUT` if
   !< normal_axis is not 1-3 (every plane then gives an empty curve and that status, as `slice_surface`).
   type(surface_stl_object), intent(in)            :: surface      !< Surface.
   integer(I4P),             intent(in)            :: normal_axis  !< Plane normal axis (1, 2, 3).
   real(R8P),                intent(in)            :: levels(:)    !< Plane positions.
   type(slice_planes_t),     intent(out)           :: planes       !< Index.
   logical,                  intent(in),  optional :: on_plane_low !< A vertex on the plane counts as below (default .false.).
   integer(I4P),             intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   type(facet_object), pointer :: f                 !< Facet cursor.
   integer(I4P), allocatable   :: ord(:), lo(:), hi(:), cnt(:), pos(:) !< Sorted planes, facet ranges, counts, cursors.
   real(R8P),    allocatable   :: sl(:)             !< Sorted finite levels.
   real(R8P)                   :: c(3), mn, mx      !< Facet coordinates along the axis, span.
   integer(I4P)                :: nf, nl, ns, i, k, b, a, z, mid !< Work.

   if (present(status)) status = STATUS_OK
   nl = size(levels, kind=I4P) ; nf = surface%get_facets_number()
   planes%on_plane_low = .false. ; if (present(on_plane_low)) planes%on_plane_low = on_plane_low
   allocate(planes%levels(nl), planes%first(nl + 1_I4P))
   planes%levels = levels
   planes%first = 1_I4P
   allocate(planes%facet(0))
   if (normal_axis < 1_I4P .or. normal_axis > 3_I4P) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   planes%normal_axis = normal_axis
   ! the finite levels in increasing order (a non-finite level lists no facet: `slice_surface_plane` rejects it)
   allocate(ord(nl))
   ns = 0_I4P
   do b = 1, nl
      if (abs(levels(b)) <= huge(1._R8P)) then
         ns = ns + 1_I4P ; ord(ns) = b
      endif
   enddo
   call sort_levels(levels, ord(1:ns))
   allocate(sl(ns))
   do k = 1, ns
      sl(k) = levels(ord(k))
   enddo
   ! per facet, the range lo..hi of sorted planes with mn <= level <= mx
   allocate(lo(max(1, nf)), hi(max(1, nf)), cnt(ns + 1_I4P))
   cnt = 0_I4P
   do i = 1, nf
      f => surface%facet_at(i)
      select case (normal_axis)
      case (1_I4P)
         c = [f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x]
      case (2_I4P)
         c = [f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y]
      case default
         c = [f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z]
      endselect
      if (any(c /= c)) then ! a NaN coordinate: every comparison false, listed everywhere
         lo(i) = 1_I4P ; hi(i) = ns
      else
         mn = c(1) ; mx = c(1)
         if (c(2) < mn) mn = c(2)
         if (c(3) < mn) mn = c(3)
         if (c(2) > mx) mx = c(2)
         if (c(3) > mx) mx = c(3)
         a = 1_I4P ; z = ns + 1_I4P ! first sorted plane with level >= mn
         do while (a < z)
            mid = (a + z) / 2_I4P
            if (sl(mid) < mn) then
               a = mid + 1_I4P
            else
               z = mid
            endif
         enddo
         lo(i) = a
         a = 0_I4P ; z = ns ! last sorted plane with level <= mx
         do while (a < z)
            mid = (a + z + 1_I4P) / 2_I4P
            if (sl(mid) > mx) then
               z = mid - 1_I4P
            else
               a = mid
            endif
         enddo
         hi(i) = a
      endif
      if (lo(i) <= hi(i)) then
         cnt(lo(i)) = cnt(lo(i)) + 1_I4P ; cnt(hi(i) + 1_I4P) = cnt(hi(i) + 1_I4P) - 1_I4P
      endif
   enddo
   ! facets per plane (difference array), then the CSR in the caller's plane order
   allocate(pos(nl + 1_I4P))
   pos = 0_I4P
   k = 0_I4P
   do a = 1, ns
      k = k + cnt(a)
      pos(ord(a)) = k
   enddo
   planes%first(1) = 1_I4P
   do b = 1, nl
      planes%first(b + 1_I4P) = planes%first(b) + pos(b)
   enddo
   deallocate(planes%facet) ; allocate(planes%facet(planes%first(nl + 1_I4P) - 1_I4P))
   pos(1:nl) = planes%first(1:nl)
   do i = 1, nf ! increasing facet index: every list comes out sorted
      do a = lo(i), hi(i)
         b = ord(a)
         planes%facet(pos(b)) = i ; pos(b) = pos(b) + 1_I4P
      enddo
   enddo
   endsubroutine slice_planes_prepare

   subroutine slice_surface_plane(surface, planes, b, curve, merge_collinear, open_loops, status)
   !< `slice_surface(surface, planes%normal_axis, planes%levels(b), curve, merge_collinear, open_loops,
   !< planes%on_plane_low, status)` through the plane sweep index (issue #100 step 3c): the same curve bit for bit, the
   !< facets visited being only those listed for the plane. `surface` must be the one the index was prepared from.
   type(surface_stl_object), intent(in)            :: surface         !< Surface.
   type(slice_planes_t),     intent(in)            :: planes          !< Index (`slice_planes_prepare`).
   integer(I4P),             intent(in)            :: b               !< Plane.
   type(curve_object),       intent(inout)         :: curve           !< Section.
   logical,                  intent(in),  optional :: merge_collinear !< Remove collinear vertices (default .false.).
   integer(I4P),             intent(out), optional :: open_loops      !< Chains that did not close.
   integer(I4P),             intent(out), optional :: status          !< STATUS_OK or STATUS_INVALID_INPUT.
   type(facet_object), pointer :: f                 !< Facet cursor.
   real(R8P),    allocatable   :: s(:,:)            !< Segments [su, sv, eu, ev; n].
   integer(I4P), allocatable   :: sp(:)             !< Patch of each segment (#66).
   real(R8P)                   :: x(3,3), pa(2), pb(2), level !< Facet vertices, crossings, plane position.
   integer(I4P)                :: n, i, k, nopen, axis !< Work.
   logical                     :: hit, merge_       !< Work.

   call curve%destroy
   if (present(status)) status = STATUS_OK
   if (present(open_loops)) open_loops = 0_I4P
   axis = planes%normal_axis
   if (axis < 1_I4P .or. axis > 3_I4P .or. b < 1_I4P .or. b > size(planes%levels, kind=I4P)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   level = planes%levels(b)
   if (.not. abs(level) <= huge(1._R8P)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   merge_ = .false. ; if (present(merge_collinear)) merge_ = merge_collinear
   allocate(s(4, max(1, planes%first(b + 1_I4P) - planes%first(b))), sp(max(1, planes%first(b + 1_I4P) - planes%first(b))))
   n = 0_I4P
   do k = planes%first(b), planes%first(b + 1_I4P) - 1_I4P
      i = planes%facet(k)
      f => surface%facet_at(i)
      x(:, 1) = [f%vertex(1)%x, f%vertex(1)%y, f%vertex(1)%z]
      x(:, 2) = [f%vertex(2)%x, f%vertex(2)%y, f%vertex(2)%z]
      x(:, 3) = [f%vertex(3)%x, f%vertex(3)%y, f%vertex(3)%z]
      call facet_segment(x, axis, level, planes%on_plane_low, pa, pb, hit)
      if (.not. hit) cycle
      n = n + 1_I4P
      s(:, n) = [pa, pb]
      sp(n) = f%patch
   enddo
   call finish_section(surface, axis, s(:, 1:n), sp(1:n), merge_, curve, nopen)
   if (present(open_loops)) open_loops = nopen
   endsubroutine slice_surface_plane

   subroutine slice_surface_planes(surface, normal_axis, levels, curves, merge_collinear, on_plane_low, status)
   !< The sections of `surface` by the planes x(normal_axis) = levels(b), any order (issue #100 step 3c): curves(b)
   !< equals `slice_surface(surface, normal_axis, levels(b), curves(b), merge_collinear, on_plane_low=on_plane_low)` bit
   !< for bit, by one plane sweep index. `status` is STATUS_INVALID_INPUT if any plane is (bad axis, non-finite level).
   type(surface_stl_object),        intent(in)            :: surface         !< Surface.
   integer(I4P),                    intent(in)            :: normal_axis     !< Plane normal axis (1, 2, 3).
   real(R8P),                       intent(in)            :: levels(:)       !< Plane positions.
   type(curve_object), allocatable, intent(inout)         :: curves(:)       !< Sections [size(levels)].
   logical,                         intent(in),  optional :: merge_collinear !< Remove collinear vertices.
   logical,                         intent(in),  optional :: on_plane_low    !< A vertex on the plane counts as below.
   integer(I4P),                    intent(out), optional :: status          !< STATUS_OK or STATUS_INVALID_INPUT.
   type(slice_planes_t) :: planes !< Index.
   integer(I4P)         :: b, st  !< Work.

   if (present(status)) status = STATUS_OK
   if (allocated(curves)) deallocate(curves)
   allocate(curves(size(levels)))
   call slice_planes_prepare(surface, normal_axis, levels, planes, on_plane_low=on_plane_low)
   do b = 1_I4P, size(levels, kind=I4P)
      call slice_surface_plane(surface, planes, b, curves(b), merge_collinear=merge_collinear, status=st)
      if (present(status) .and. st /= STATUS_OK) status = st
   enddo
   endsubroutine slice_surface_planes

   ! private
   pure subroutine facet_segment(x, normal_axis, level, low, pa, pb, hit)
   !< The directed segment the plane cuts from the facet x(:, 1:3) (the classification and canonical crossings of
   !< `slice_surface`): hit false if the plane misses the facet or touches it at a vertex only.
   real(R8P),    intent(in)  :: x(3,3)       !< Facet vertices.
   integer(I4P), intent(in)  :: normal_axis  !< Plane normal axis.
   real(R8P),    intent(in)  :: level        !< Plane position.
   logical,      intent(in)  :: low          !< A vertex on the plane counts as below.
   real(R8P),    intent(out) :: pa(2), pb(2) !< Segment start, end.
   logical,      intent(out) :: hit          !< A segment.
   logical                   :: above(3)     !< Vertex sides.
   integer(I4P)              :: e, en, ia, ib !< Work.

   hit = .false. ; pa = 0._R8P ; pb = 0._R8P
   if (low) then
      above = x(normal_axis, :) > level  ! plane at level + eps (the cut-cell rule, issue #40)
   else
      above = x(normal_axis, :) >= level ! plane at level - eps
   endif
   if (all(above) .or. .not. any(above)) return
   ia = 0_I4P ; ib = 0_I4P
   do e = 1, 3
      en = mod(e, 3_I4P) + 1_I4P
      if (above(e) .and. .not. above(en)) ia = e ! edge e -> en goes above -> below
      if (.not. above(e) .and. above(en)) ib = e ! edge e -> en goes below -> above
   enddo
   pa = crossing(x(:, ia), x(:, mod(ia, 3_I4P) + 1_I4P), normal_axis, level)
   pb = crossing(x(:, ib), x(:, mod(ib, 3_I4P) + 1_I4P), normal_axis, level)
   hit = .not. all(pa == pb) ! false: the facet touches the plane at a vertex
   endsubroutine facet_segment

   subroutine finish_section(surface, normal_axis, s, sp, merge_, curve, nopen)
   !< Chain the segments of a section into loops (optionally dropping collinear vertices) and adopt them into curve,
   !< with the surface's patch names.
   type(surface_stl_object), intent(in)    :: surface     !< Surface (patch names).
   integer(I4P),             intent(in)    :: normal_axis !< Plane normal axis.
   real(R8P),                intent(in)    :: s(:,:)      !< Segments [su, sv, eu, ev; n].
   integer(I4P),             intent(in)    :: sp(:)       !< Patch of each segment.
   logical,                  intent(in)    :: merge_      !< Remove collinear vertices.
   type(curve_object),       intent(inout) :: curve       !< Section.
   integer(I4P),             intent(out)   :: nopen       !< Chains that did not close.
   integer(I4P), allocatable :: lp(:)        !< Patch of each loop vertex (#66).
   real(R8P),    allocatable :: lu(:), lv(:) !< Loop vertices.
   integer(I4P), allocatable :: sizes(:)     !< Loop sizes.
   logical,      allocatable :: closed(:)    !< Loop closed flags.
   character(64), allocatable :: names(:)    !< Patch names.
   integer(I4P)              :: i            !< Work.

   call chain(s, lu, lv, sizes, closed, nopen, sp, lp)
   if (merge_) call drop_collinear(lu, lv, sizes, closed, lp)
   allocate(names(surface%get_patches_number()))
   do i = 1, size(names)
      names(i) = surface%patch_name(i)
   enddo
   call curve%adopt_loops(u=lu, v=lv, loop_sizes=sizes, normal_axis=normal_axis, closed=closed, patches=lp, &
                          patch_names=names)
   endsubroutine finish_section

   pure subroutine sort_levels(levels, ord)
   !< Sort the plane indices ord by level, ties by index (bottom-up merge sort, deterministic).
   real(R8P),    intent(in)    :: levels(:) !< Plane positions (finite at ord).
   integer(I4P), intent(inout) :: ord(:)    !< Plane indices.
   integer(I4P), allocatable   :: tmp(:)    !< Buffer.
   integer(I4P)                :: width, lo, mid, hi, i, j, k, n !< Work.

   n = size(ord, kind=I4P)
   allocate(tmp(n))
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (levels(ord(j)) < levels(ord(i))) then
               tmp(k) = ord(j) ; j = j + 1_I4P
            else
               tmp(k) = ord(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i <= mid)
            tmp(k) = ord(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j <= hi)
            tmp(k) = ord(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         lo = lo + 2_I4P * width
      enddo
      ord = tmp
      width = 2_I4P * width
   enddo
   endsubroutine sort_levels

   pure function crossing(p, q, axis, level) result(c)
   !< In-plane coordinates of the crossing of the edge pq with the plane (canonical).
   !<
   !< The endpoints are ordered lexicographically first, so (p, q) and (q, p) give the same bits;
   !< a parameter of exactly 0 or 1 returns the endpoint itself.
   real(R8P),    intent(in) :: p(3), q(3) !< Edge endpoints (one below, one at or above).
   integer(I4P), intent(in) :: axis       !< Plane normal axis.
   real(R8P),    intent(in) :: level      !< Plane position.
   real(R8P)                :: c(2)       !< Crossing (u, v).
   real(R8P)                :: a(3), b(3), t, x(3) !< Ordered endpoints, parameter, point.

   if (lex_less(p, q)) then
      a = p ; b = q
   else
      a = q ; b = p
   endif
   t = (level - a(axis)) / (b(axis) - a(axis))
   if (t == 0._R8P) then
      x = a
   elseif (t == 1._R8P) then
      x = b
   else
      x = a + t * (b - a)
   endif
   c = in_plane(x, axis)
   endfunction crossing

   pure function in_plane(x, axis) result(uv)
   !< (u, v) of x in the plane normal to axis: (y, z), (z, x), (x, y) for axis 1, 2, 3.
   real(R8P),    intent(in) :: x(3)  !< Point.
   integer(I4P), intent(in) :: axis  !< Normal axis.
   real(R8P)                :: uv(2) !< In-plane coordinates.

   uv = [x(mod(axis, 3_I4P) + 1_I4P), x(mod(axis + 1_I4P, 3_I4P) + 1_I4P)]
   endfunction in_plane

   pure function lex_less(a, b) result(less)
   !< Lexicographic a < b.
   real(R8P), intent(in) :: a(:), b(:) !< Points.
   logical               :: less       !< Result.
   integer(I4P)          :: k          !< Counter.

   less = .false.
   do k = 1, size(a)
      if (a(k) < b(k)) then
         less = .true.
         return
      elseif (a(k) > b(k)) then
         return
      endif
   enddo
   endfunction lex_less

   subroutine chain(s, lu, lv, sizes, closed, nopen, sp, lp)
   !< Chain directed segments into loops by exact endpoint equality.
   !<
   !< Segment starts are sorted lexicographically; the successor of a segment is a not yet used
   !< segment starting at its end (binary search, first match). Chains are started first at
   !< segments that no segment ends at (open polylines), then at any unused segment (closed loops);
   !< a chain that comes back to its first point is closed. Deterministic.
   real(R8P),                 intent(in)  :: s(:,:)       !< Segments [su, sv, eu, ev; n].
   real(R8P),    allocatable, intent(out) :: lu(:), lv(:) !< Loop vertices.
   integer(I4P), allocatable, intent(out) :: sizes(:)     !< Loop sizes.
   logical,      allocatable, intent(out) :: closed(:)    !< Loop closed flags.
   integer(I4P),              intent(out) :: nopen        !< Open chains.
   integer(I4P),              intent(in)  :: sp(:)        !< Patch of each segment.
   integer(I4P), allocatable, intent(out) :: lp(:)        !< Patch of the segment starting at each loop vertex.
   integer(I4P), allocatable              :: ord(:)       !< Segments sorted by start.
   logical,      allocatable              :: used(:), has_pred(:) !< Flags.
   integer(I4P)                           :: n, i, k, cur, nv, nl, pass, start !< Work.

   n = size(s, dim=2, kind=I4P)
   allocate(lu(2 * n), lv(2 * n), sizes(n), closed(n), used(n), has_pred(n), lp(2 * n)) ! an open chain adds its end point
   allocate(ord(n))
   ord = [(i, i = 1, n)]
   call sort_starts(s, ord)
   has_pred = .false.
   used = .false. ! nothing used yet: find_start returns the first match
   do i = 1, n
      k = find_start(s, ord, s(3:4, i), used)
      if (k > 0_I4P) has_pred(k) = .true.
   enddo
   nv = 0_I4P ; nl = 0_I4P ; nopen = 0_I4P
   do pass = 1, 2
      do i = 1, n
         start = ord(i)
         if (used(start)) cycle
         if (pass == 1 .and. has_pred(start)) cycle
         nl = nl + 1_I4P ; sizes(nl) = 0_I4P
         cur = start
         do
            used(cur) = .true.
            nv = nv + 1_I4P ; lu(nv) = s(1, cur) ; lv(nv) = s(2, cur) ; lp(nv) = sp(cur) ; sizes(nl) = sizes(nl) + 1_I4P
            k = find_start(s, ord, s(3:4, cur), used)
            if (k == 0_I4P) exit
            cur = k
         enddo
         closed(nl) = all(s(3:4, cur) == s(1:2, start))
         if (.not. closed(nl)) then ! keep the final end point of an open chain
            nv = nv + 1_I4P ; lu(nv) = s(3, cur) ; lv(nv) = s(4, cur) ; lp(nv) = sp(cur) ; sizes(nl) = sizes(nl) + 1_I4P
            nopen = nopen + 1_I4P
         endif
      enddo
   enddo
   lu = lu(1:nv) ; lv = lv(1:nv) ; lp = lp(1:nv) ; sizes = sizes(1:nl) ; closed = closed(1:nl)
   endsubroutine chain

   pure function find_start(s, ord, p, used) result(k)
   !< First unused segment (in sorted order) starting exactly at p, or 0.
   real(R8P),    intent(in) :: s(:,:)  !< Segments.
   integer(I4P), intent(in) :: ord(:)  !< Segments sorted by start.
   real(R8P),    intent(in) :: p(2)    !< Point.
   logical,      intent(in) :: used(:) !< Used flags.
   integer(I4P)             :: k       !< Segment, or 0.
   integer(I4P)             :: lo, hi, mid !< Binary search.

   k = 0_I4P
   lo = 1_I4P ; hi = size(ord, kind=I4P)
   do while (lo < hi) ! first position whose start is not less than p
      mid = (lo + hi) / 2_I4P
      if (lex_less(s(1:2, ord(mid)), p)) then
         lo = mid + 1_I4P
      else
         hi = mid
      endif
   enddo
   do mid = lo, size(ord, kind=I4P)
      if (any(s(1:2, ord(mid)) /= p)) exit
      if (.not. used(ord(mid))) then
         k = ord(mid)
         return
      endif
   enddo
   endfunction find_start

   pure subroutine sort_starts(s, ord)
   !< Sort ord by segment start, lexicographically, ties by index (merge sort, deterministic).
   real(R8P),    intent(in)    :: s(:,:) !< Segments.
   integer(I4P), intent(inout) :: ord(:) !< Indices.
   integer(I4P), allocatable   :: tmp(:) !< Buffer.
   integer(I4P)                :: width, lo, mid, hi, i, j, k, n !< Work.

   n = size(ord, kind=I4P)
   allocate(tmp(n))
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (before(ord(j), ord(i))) then
               tmp(k) = ord(j) ; j = j + 1_I4P
            else
               tmp(k) = ord(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i <= mid)
            tmp(k) = ord(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j <= hi)
            tmp(k) = ord(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         lo = lo + 2_I4P * width
      enddo
      ord = tmp
      width = 2_I4P * width
   enddo
   contains
      pure function before(a, b) result(yes)
      !< Segment a sorts strictly before b (start, then index).
      integer(I4P), intent(in) :: a, b !< Segments.
      logical                  :: yes  !< Result.

      if (lex_less(s(1:2, a), s(1:2, b))) then
         yes = .true.
      elseif (lex_less(s(1:2, b), s(1:2, a))) then
         yes = .false.
      else
         yes = a < b
      endif
      endfunction before
   endsubroutine sort_starts

   pure subroutine drop_collinear(lu, lv, sizes, closed, lp)
   !< Remove vertices between two collinear segments running the same way (exact orient2d_xy);
   !< the end vertices of open loops are kept, as is any loop that would drop below 3 vertices.
   real(R8P),    allocatable, intent(inout) :: lu(:), lv(:) !< Loop vertices.
   integer(I4P), allocatable, intent(inout) :: sizes(:)     !< Loop sizes.
   logical,                   intent(in)    :: closed(:)    !< Loop closed flags.
   integer(I4P), allocatable, intent(inout) :: lp(:)        !< Segment patches: a vertex where the patch changes stays.
   real(R8P),    allocatable                :: nu(:), nv(:) !< Kept vertices.
   integer(I4P), allocatable                :: np_(:)       !< Kept patches.
   logical,      allocatable                :: keep(:)      !< Keep flags of a loop.
   integer(I4P)                             :: l, first, m, i, ip, in, k, kept !< Work.

   allocate(nu(size(lu)), nv(size(lv)), np_(size(lp)))
   first = 1_I4P ; k = 0_I4P
   do l = 1, size(sizes)
      m = sizes(l)
      allocate(keep(m))
      keep = .true.
      if (m > 3_I4P) then
         do i = 1, m
            if (.not. closed(l) .and. (i == 1_I4P .or. i == m)) cycle
            ip = first + mod(i - 2_I4P + m, m) ; in = first + mod(i, m)
            if (lp(ip) /= lp(first + i - 1)) cycle
            if (orient2d_xy(lu(ip), lv(ip), lu(first + i - 1), lv(first + i - 1), lu(in), lv(in)) == 0._R8P .and. &
                (lu(first + i - 1) - lu(ip)) * (lu(in) - lu(first + i - 1)) + &
                (lv(first + i - 1) - lv(ip)) * (lv(in) - lv(first + i - 1)) > 0._R8P) keep(i) = .false.
         enddo
         if (count(keep) < 3) keep = .true.
      endif
      kept = 0_I4P
      do i = 1, m
         if (.not. keep(i)) cycle
         k = k + 1_I4P ; kept = kept + 1_I4P
         nu(k) = lu(first + i - 1) ; nv(k) = lv(first + i - 1) ; np_(k) = lp(first + i - 1)
      enddo
      first = first + m
      sizes(l) = kept
      deallocate(keep)
   enddo
   lu = nu(1:k) ; lv = nv(1:k) ; lp = np_(1:k)
   endsubroutine drop_collinear
endmodule fossil_curve_slice
