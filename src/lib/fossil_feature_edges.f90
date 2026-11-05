!< FOSSIL, feature edges: sharp, boundary and non-manifold edges, corners and polylines (issue #61, roadmap #58 WP3).
module fossil_feature_edges
!< FOSSIL, feature edges: sharp, boundary and non-manifold edges, corners and polylines (issue #61, roadmap #58 WP3).
!<
!< Pure routines on plain arrays: facets as vertex ids f_v(3, nf), vertex coordinates vcoord(3, nv).
!< - Sharp edge: two incident facets whose unit normals make an angle above `angle` (`is_sharp_dihedral`; the facets must
!<   be consistently oriented, as `sanitize` leaves them). A facet of zero area makes no normal: its edges are not
!<   classified sharp, and are counted in `degenerate`.
!< - Boundary edge: one incident facet. Non-manifold edge: three or more.
!< - Corner: a feature vertex with a number of feature edges other than 2, or with 2 whose turning angle exceeds
!<   `corner_angle`.
!< - Polylines: the feature edges chained through the other (valence-2, non-corner) feature vertices; a polyline with no
!<   corner on it is closed (its first vertex repeated at its end).
!< `isotropic_remesh` locks its features through `is_sharp_dihedral` and `FEATURE_EDGE_DEFAULT_ANGLE` too.
use penf, only : I4P, R8P

implicit none
private
public :: FEATURE_EDGE_SHARP, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD, FEATURE_EDGE_DEFAULT_ANGLE
public :: is_sharp_dihedral, feature_cos, find_feature_edges

integer(I4P), parameter :: FEATURE_EDGE_SHARP       = 1_I4P  !< Two facets, dihedral above the threshold.
integer(I4P), parameter :: FEATURE_EDGE_BOUNDARY    = 2_I4P  !< One facet.
integer(I4P), parameter :: FEATURE_EDGE_NONMANIFOLD = 3_I4P  !< Three or more facets.
real(R8P),    parameter :: FEATURE_EDGE_DEFAULT_ANGLE = 30._R8P !< Default threshold [degrees].

contains
   pure function feature_cos(angle) result(c)
   !< Cosine of an angle in degrees, the one threshold both `find_feature_edges` and `isotropic_remesh` use.
   real(R8P), intent(in) :: angle !< Angle [degrees].
   real(R8P)             :: c     !< Cosine.

   c = cos(angle * 3.141592653589793_R8P / 180._R8P)
   endfunction feature_cos

   pure function is_sharp_dihedral(n1, n2, cos_threshold) result(sharp)
   !< True when the unit normals n1, n2 make an angle above acos(cos_threshold): n1 . n2 < cos_threshold.
   real(R8P), intent(in) :: n1(3), n2(3)   !< Unit normals.
   real(R8P), intent(in) :: cos_threshold  !< Cosine of the threshold angle.
   logical               :: sharp          !< Sharp.

   sharp = n1(1) * n2(1) + n1(2) * n2(2) + n1(3) * n2(3) < cos_threshold
   endfunction is_sharp_dihedral

   subroutine find_feature_edges(f_v, vcoord, angle, corner_angle, include_boundary, include_nonmanifold, edges, kinds, &
                                 corners, polyline_start, polyline_vertices, degenerate)
   !< Feature edges of the mesh (f_v, vcoord), their kinds, the corners and the polylines (CSR: polyline p is
   !< polyline_vertices(polyline_start(p):polyline_start(p+1)-1)). Edges are listed in increasing (low id, high id) order.
   integer(I4P),              intent(in)  :: f_v(:,:)             !< Facets [3, nf].
   real(R8P),                 intent(in)  :: vcoord(:,:)          !< Vertices [3, nv].
   real(R8P),                 intent(in)  :: angle                !< Sharp threshold [degrees].
   real(R8P),                 intent(in)  :: corner_angle         !< Corner turning threshold [degrees].
   logical,                   intent(in)  :: include_boundary     !< Boundary edges are features.
   logical,                   intent(in)  :: include_nonmanifold  !< Non-manifold edges are features.
   integer(I4P), allocatable, intent(out) :: edges(:,:)           !< Feature edges [2, ne].
   integer(I4P), allocatable, intent(out) :: kinds(:)             !< FEATURE_EDGE_* [ne].
   integer(I4P), allocatable, intent(out) :: corners(:)           !< Corner vertices.
   integer(I4P), allocatable, intent(out) :: polyline_start(:)    !< Polyline offsets [np + 1].
   integer(I4P), allocatable, intent(out) :: polyline_vertices(:) !< Polyline vertices.
   integer(I4P),              intent(out) :: degenerate           !< Two-facet edges skipped for a zero-area facet.
   integer(I4P), allocatable              :: lo(:), hi(:), fac(:), order(:) !< Half-edges.
   integer(I4P), allocatable              :: e_lo(:), e_hi(:), e_n(:), e_f1(:), e_f2(:), e_kind(:) !< Unique edges.
   integer(I4P), allocatable              :: valence(:), adj_start(:), adj(:), fill(:) !< Vertex -> feature edges.
   logical,      allocatable              :: is_corner(:), used(:) !< Corners, edges already chained.
   integer(I4P), allocatable              :: pv(:), ps(:)        !< Polyline buffers.
   real(R8P)                              :: n1(3), n2(3), cos_thr, cos_corner, a(3), b(3) !< Work.
   integer(I4P)                           :: nf, nv, nh, ne, nfe, h, k, f, e, v, w, np, nps, ecur, vcur !< Counters.

   nf = size(f_v, 2) ; nv = size(vcoord, 2) ; nh = 3 * nf
   cos_thr = feature_cos(angle)
   cos_corner = feature_cos(corner_angle)
   degenerate = 0_I4P
   ! Unique edges with their incidence count and first two facets (half-edges sorted by (lo, hi, facet)).
   allocate(lo(nh), hi(nh), fac(nh), order(nh))
   h = 0
   do f = 1, nf
      do k = 1, 3
         h = h + 1
         lo(h) = min(f_v(k, f), f_v(mod(k, 3) + 1, f)) ; hi(h) = max(f_v(k, f), f_v(mod(k, 3) + 1, f)) ; fac(h) = f
         order(h) = h
      enddo
   enddo
   call sort_half_edges(lo, hi, fac, order)
   allocate(e_lo(nh), e_hi(nh), e_n(nh), e_f1(nh), e_f2(nh))
   ne = 0
   do k = 1, nh
      h = order(k)
      if (ne > 0) then
         if (lo(h) == e_lo(ne) .and. hi(h) == e_hi(ne)) then
            e_n(ne) = e_n(ne) + 1
            if (e_n(ne) == 2) e_f2(ne) = fac(h)
            cycle
         endif
      endif
      ne = ne + 1
      e_lo(ne) = lo(h) ; e_hi(ne) = hi(h) ; e_n(ne) = 1 ; e_f1(ne) = fac(h) ; e_f2(ne) = 0
   enddo
   ! Classification.
   allocate(e_kind(ne))
   e_kind = 0
   do e = 1, ne
      if (e_n(e) == 1) then
         if (include_boundary) e_kind(e) = FEATURE_EDGE_BOUNDARY
      elseif (e_n(e) >= 3) then
         if (include_nonmanifold) e_kind(e) = FEATURE_EDGE_NONMANIFOLD
      else
         n1 = unit_normal(e_f1(e)) ; n2 = unit_normal(e_f2(e))
         if (all(n1 == 0._R8P) .or. all(n2 == 0._R8P)) then
            degenerate = degenerate + 1
         elseif (is_sharp_dihedral(n1, n2, cos_thr)) then
            e_kind(e) = FEATURE_EDGE_SHARP
         endif
      endif
   enddo
   nfe = count(e_kind > 0)
   allocate(edges(2, nfe), kinds(nfe))
   k = 0
   do e = 1, ne
      if (e_kind(e) == 0) cycle
      k = k + 1 ; edges(:, k) = [e_lo(e), e_hi(e)] ; kinds(k) = e_kind(e)
   enddo
   ! Vertex -> feature edges (CSR), valence, corners.
   allocate(valence(nv), adj_start(nv + 1), adj(2 * nfe), fill(nv), is_corner(nv))
   valence = 0
   do e = 1, nfe
      valence(edges(1, e)) = valence(edges(1, e)) + 1 ; valence(edges(2, e)) = valence(edges(2, e)) + 1
   enddo
   adj_start(1) = 1
   do v = 1, nv
      adj_start(v + 1) = adj_start(v) + valence(v)
   enddo
   fill = adj_start(1:nv)
   do e = 1, nfe
      do k = 1, 2
         v = edges(k, e) ; adj(fill(v)) = e ; fill(v) = fill(v) + 1
      enddo
   enddo
   is_corner = valence > 0 .and. valence /= 2
   do v = 1, nv
      if (valence(v) /= 2) cycle
      a = vcoord(:, v) - vcoord(:, other_end(adj(adj_start(v)), v))
      b = vcoord(:, other_end(adj(adj_start(v) + 1), v)) - vcoord(:, v)
      if (norm2(a) > 0._R8P .and. norm2(b) > 0._R8P) then
         if (dot_product(a, b) / (norm2(a) * norm2(b)) < cos_corner) is_corner(v) = .true.
      endif
   enddo
   corners = pack([(v, v = 1, nv)], is_corner)
   ! Polylines: from every corner along each unused edge, then the closed loops left.
   allocate(used(nfe), pv(2 * nfe + 1), ps(nfe + 1))
   used = .false. ; np = 0 ; nps = 0
   do v = 1, nv
      if (.not. is_corner(v)) cycle
      do k = adj_start(v), adj_start(v + 1) - 1
         if (used(adj(k))) cycle
         call walk(v, adj(k))
      enddo
   enddo
   do e = 1, nfe
      if (.not. used(e)) call walk(edges(1, e), e)
   enddo
   ps(np + 1) = nps + 1
   polyline_start = ps(1:np + 1)
   polyline_vertices = pv(1:nps)
   contains
      function unit_normal(ff) result(n)
      !< Unit normal of facet ff from the vertex order (zero for a zero-area facet).
      integer(I4P), intent(in) :: ff
      real(R8P)                :: n(3), u(3), w_(3), l

      u = vcoord(:, f_v(2, ff)) - vcoord(:, f_v(1, ff)) ; w_ = vcoord(:, f_v(3, ff)) - vcoord(:, f_v(1, ff))
      n = [u(2) * w_(3) - u(3) * w_(2), u(3) * w_(1) - u(1) * w_(3), u(1) * w_(2) - u(2) * w_(1)]
      l = norm2(n)
      if (l > 0._R8P) then
         n = n / l
      else
         n = 0._R8P
      endif
      endfunction unit_normal

      function other_end(ee, vv) result(ww)
      !< The end of feature edge ee other than vv.
      integer(I4P), intent(in) :: ee, vv
      integer(I4P)             :: ww

      ww = merge(edges(2, ee), edges(1, ee), edges(1, ee) == vv)
      endfunction other_end

      subroutine walk(v0, e0)
      !< Chain from vertex v0 along edge e0 until a corner, a dead end or back to v0.
      integer(I4P), intent(in) :: v0, e0
      integer(I4P)             :: kk

      np = np + 1 ; ps(np) = nps + 1
      nps = nps + 1 ; pv(nps) = v0
      vcur = v0 ; ecur = e0
      do
         used(ecur) = .true.
         w = other_end(ecur, vcur)
         nps = nps + 1 ; pv(nps) = w
         if (is_corner(w) .or. w == v0) exit
         ecur = 0
         do kk = adj_start(w), adj_start(w + 1) - 1
            if (.not. used(adj(kk))) then
               ecur = adj(kk) ; exit
            endif
         enddo
         if (ecur == 0) exit
         vcur = w
      enddo
      endsubroutine walk
   endsubroutine find_feature_edges

   subroutine sort_half_edges(lo, hi, fac, order)
   !< Sort the permutation `order` by (lo, hi, fac), ascending: a bottom-up merge sort (stable, O(n log n), deterministic).
   integer(I4P), intent(in)    :: lo(:), hi(:), fac(:) !< Keys.
   integer(I4P), intent(inout) :: order(:)             !< Permutation.
   integer(I4P), allocatable   :: tmp(:)               !< Merge buffer.
   integer(I4P)                :: n, width, i, l, m, r, p, q, k !< Counters.

   n = size(order)
   allocate(tmp(n))
   width = 1
   do while (width < n)
      i = 1
      do while (i <= n)
         l = i ; m = min(i + width, n + 1) ; r = min(i + 2 * width, n + 1)
         p = l ; q = m ; k = l
         do while (p < m .and. q < r)
            if (less(order(q), order(p))) then
               tmp(k) = order(q) ; q = q + 1
            else
               tmp(k) = order(p) ; p = p + 1
            endif
            k = k + 1
         enddo
         do while (p < m)
            tmp(k) = order(p) ; p = p + 1 ; k = k + 1
         enddo
         do while (q < r)
            tmp(k) = order(q) ; q = q + 1 ; k = k + 1
         enddo
         i = i + 2 * width
      enddo
      order = tmp
      width = 2 * width
   enddo
   contains
      pure function less(x, y) result(lt)
      !< (lo, hi, fac) of half-edge x before that of y.
      integer(I4P), intent(in) :: x, y
      logical                  :: lt

      if (lo(x) /= lo(y)) then
         lt = lo(x) < lo(y)
      elseif (hi(x) /= hi(y)) then
         lt = hi(x) < hi(y)
      else
         lt = fac(x) < fac(y)
      endif
      endfunction less
   endsubroutine sort_half_edges
endmodule fossil_feature_edges
