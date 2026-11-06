!< FOSSIL, constrained triangulation of a facet on implicit points (issue #65).
module fossil_cdt_implicit
!< FOSSIL, constrained triangulation of a facet on implicit points (issue #65).
!<
!< `cdt_facet` triangulates a triangle (the cut facet: points 1, 2, 3) with points inside it or on its edges (the cut
!< points, distinct) and constraint segments between them (the cuts). Every validity decision is the exact sign of
!< `ip_orient2d` in the coordinate plane that drops axis `drop` (where the facet's projection is not degenerate):
!< - insertion: each point splits the triangle containing it in three, or the two triangles beside the edge it lies
!<   on in two (no super-triangle, no tolerance: the facet's own edges are the boundary);
!< - constraints: each is first split at the points lying on it (exact collinearity, exact betweenness by
!<   `ip_compare`), then recovered by Sloan's flips of the edges crossing it (exact crossing and convexity tests);
!< - quality: Lawson flips of the unconstrained edges by `incircle` on the rounded coordinates, a heuristic that
!<   only chooses among valid triangulations (every flip is of an exactly convex quadrilateral), bounded in sweeps.
!< The output triangles have the orientation of (1, 2, 3).
use fossil_implicit_points, only : implicit_point_t, ip_orient2d, ip_compare, ip_approx
use fossil_predicates,      only : incircle
use penf,                   only : I4P, R8P
use vecfor,                 only : vector_R8P

implicit none
private
public :: cdt_facet
public :: CDTI_OK, CDTI_DEGENERATE, CDTI_OUTSIDE, CDTI_DUPLICATE, CDTI_CROSSING, CDTI_RECOVERY

integer(I4P), parameter :: CDTI_OK         = 0_I4P !< Success.
integer(I4P), parameter :: CDTI_DEGENERATE = 1_I4P !< Points 1, 2, 3 collinear in the projection.
integer(I4P), parameter :: CDTI_OUTSIDE    = 2_I4P !< A point outside the triangle (1, 2, 3).
integer(I4P), parameter :: CDTI_DUPLICATE  = 3_I4P !< Two coincident points.
integer(I4P), parameter :: CDTI_CROSSING   = 4_I4P !< Two constraints cross.
integer(I4P), parameter :: CDTI_RECOVERY   = 5_I4P !< A constraint could not be recovered.

integer(I4P), parameter :: MAX_SWEEPS = 30_I4P     !< Bound on the Lawson quality sweeps.

type :: mesh_t
   !< Triangles: vertices `v` (counter-clockwise in the oriented projection), neighbours `nb` (edge k is opposite
   !< v(k), 0 on the facet boundary), constrained flags `c`.
   integer(I4P)              :: n = 0_I4P     !< Triangles in use.
   integer(I4P), allocatable :: v(:, :)       !< Vertices (3, cap).
   integer(I4P), allocatable :: nb(:, :)      !< Neighbours (3, cap).
   logical,      allocatable :: c(:, :)       !< Constrained edges (3, cap).
   integer(I4P)              :: last = 1_I4P  !< Start of the next point location.
endtype mesh_t

contains
   subroutine cdt_facet(pts, drop, segs, tris, ntri, status, tcon)
   !< Constrained triangulation of the triangle pts(1:3) with the points pts(4:) and the constraints `segs`.
   type(implicit_point_t),    intent(in)  :: pts(:)     !< Points: 1..3 the triangle, then distinct points in it.
   integer(I4P),              intent(in)  :: drop       !< Dropped axis of the projection.
   integer(I4P),              intent(in)  :: segs(:, :) !< Constraints (2, m), local indices.
   integer(I4P), allocatable, intent(out) :: tris(:, :) !< Triangles (3, ntri), with the orientation of (1, 2, 3).
   integer(I4P),              intent(out) :: ntri       !< Number of triangles.
   integer(I4P),              intent(out) :: status     !< CDTI_*.
   logical, allocatable,      intent(out), optional :: tcon(:, :) !< Constrained edges (3, ntri), edge k opposite vertex k.
   type(mesh_t)                           :: m
   integer(I4P)                           :: s0, n, ip, t, e, kind
   integer(I4P), allocatable              :: sub(:, :)
   integer(I4P)                           :: nsub, k
   real(R8P), allocatable                 :: xy(:, :)

   status = CDTI_OK ; ntri = 0_I4P
   allocate(tris(3, 0))
   n = size(pts)
   s0 = ip_orient2d(pts(1), pts(2), pts(3), drop)
   if (s0 == 0_I4P) then
      status = CDTI_DEGENERATE ; return
   endif
   allocate(m%v(3, 2 * n), m%nb(3, 2 * n), m%c(3, 2 * n))
   m%n = 1_I4P ; m%v(:, 1) = [1_I4P, 2_I4P, 3_I4P] ; m%nb(:, 1) = 0_I4P ; m%c(:, 1) = .false.
   do ip = 4, n
      call locate(m, pts, drop, s0, ip, t, kind, e)
      select case (kind)
      case (0)
         call split3(m, t, ip)
      case (1)
         call split_edge(m, t, e, ip)
      case (2)
         status = CDTI_DUPLICATE ; return
      case default
         status = CDTI_OUTSIDE ; return
      endselect
   enddo
   allocate(xy(2, n))
   do ip = 1, n
      xy(:, ip) = projected(ip_approx(pts(ip)), drop)
   enddo
   call split_constraints(pts, drop, s0, xy, segs, sub, nsub)
   do k = 1, nsub
      call recover(m, pts, drop, s0, sub(1, k), sub(2, k), status)
      if (status /= CDTI_OK) return
   enddo
   call lawson(m, pts, drop, s0, xy)
   deallocate(tris)
   allocate(tris(3, m%n))
   tris = m%v(:, 1:m%n)
   ntri = m%n
   if (present(tcon)) then
      allocate(tcon(3, m%n)) ; tcon = m%c(:, 1:m%n)
   endif
   endsubroutine cdt_facet

   ! orientation ------------------------------------------------------------------------------------------------------
   function o2(pts, drop, s0, i, j, k) result(s)
   !< Exact orientation of points i, j, k, positive as (1, 2, 3).
   type(implicit_point_t), intent(in) :: pts(:)
   integer(I4P),           intent(in) :: drop, s0, i, j, k
   integer(I4P)                       :: s

   s = s0 * ip_orient2d(pts(i), pts(j), pts(k), drop)
   endfunction o2

   pure function projected(x, drop) result(r)
   !< Kept coordinates (i, j) = (drop + 1, drop + 2) mod 3 of x.
   real(R8P),    intent(in) :: x(3)
   integer(I4P), intent(in) :: drop
   real(R8P)                :: r(2)

   r = [x(mod(drop, 3_I4P) + 1_I4P), x(mod(drop + 1_I4P, 3_I4P) + 1_I4P)]
   endfunction projected

   ! insertion --------------------------------------------------------------------------------------------------------
   subroutine locate(m, pts, drop, s0, ip, t, kind, e)
   !< The triangle containing point ip: kind 0 inside, 1 on its edge e, 2 on a vertex, 3 outside the facet. A walk
   !< from the last triangle (start edge rotated each step, so it cannot cycle forever), then a scan if it is long.
   type(mesh_t),           intent(inout) :: m
   type(implicit_point_t), intent(in)    :: pts(:)
   integer(I4P),           intent(in)    :: drop, s0, ip
   integer(I4P),           intent(out)   :: t, kind, e
   integer(I4P)                          :: steps, kk, ee, s(3), a, b
   logical                               :: moved

   t = min(max(m%last, 1_I4P), m%n)
   kind = 3_I4P ; e = 0_I4P
   do steps = 1, 4 * m%n + 8
      moved = .false.
      do kk = 0, 2
         ee = mod(steps + kk, 3_I4P) + 1_I4P
         a = m%v(mod(ee, 3_I4P) + 1_I4P, t) ; b = m%v(mod(ee + 1_I4P, 3_I4P) + 1_I4P, t)
         if (o2(pts, drop, s0, a, b, ip) < 0_I4P) then
            if (m%nb(ee, t) == 0_I4P) return ! outside the facet
            t = m%nb(ee, t) ; moved = .true. ; exit
         endif
      enddo
      if (.not. moved) then
         call classify(t) ; return
      endif
   enddo
   do t = 1, m%n ! fallback scan
      do ee = 1, 3
         a = m%v(mod(ee, 3_I4P) + 1_I4P, t) ; b = m%v(mod(ee + 1_I4P, 3_I4P) + 1_I4P, t)
         s(ee) = o2(pts, drop, s0, a, b, ip)
      enddo
      if (all(s >= 0_I4P)) then
         call classify(t) ; return
      endif
   enddo
   t = 0_I4P
   contains
      subroutine classify(tt)
      integer(I4P), intent(in) :: tt
      integer(I4P)             :: k, z

      z = 0_I4P
      do k = 1, 3
         a = m%v(mod(k, 3_I4P) + 1_I4P, tt) ; b = m%v(mod(k + 1_I4P, 3_I4P) + 1_I4P, tt)
         s(k) = o2(pts, drop, s0, a, b, ip)
         if (s(k) == 0_I4P) then
            z = z + 1_I4P ; e = k
         endif
      enddo
      kind = min(z, 2_I4P)
      m%last = tt
      endsubroutine classify
   endsubroutine locate

   subroutine split3(m, t, p)
   !< Split triangle t = (a, b, c) at the interior point p into (p, b, c), (a, p, c), (a, b, p).
   type(mesh_t), intent(inout) :: m
   integer(I4P), intent(in)    :: t, p
   integer(I4P)                :: a, b, c, n1, n2, n3, t2, t3
   logical                     :: c1, c2, c3

   a = m%v(1, t) ; b = m%v(2, t) ; c = m%v(3, t)
   n1 = m%nb(1, t) ; n2 = m%nb(2, t) ; n3 = m%nb(3, t)
   c1 = m%c(1, t) ; c2 = m%c(2, t) ; c3 = m%c(3, t)
   t2 = m%n + 1_I4P ; t3 = m%n + 2_I4P ; m%n = m%n + 2_I4P
   m%v(:, t)  = [p, b, c] ; m%nb(:, t)  = [n1, t2, t3] ; m%c(:, t)  = [c1, .false., .false.]
   m%v(:, t2) = [a, p, c] ; m%nb(:, t2) = [t, n2, t3]  ; m%c(:, t2) = [.false., c2, .false.]
   m%v(:, t3) = [a, b, p] ; m%nb(:, t3) = [t, t2, n3]  ; m%c(:, t3) = [.false., .false., c3]
   if (n2 /= 0_I4P) call repoint(m, n2, t, t2)
   if (n3 /= 0_I4P) call repoint(m, n3, t, t3)
   endsubroutine split3

   subroutine split_edge(m, t, e, p)
   !< Split the edge e of t at p (on it, strictly inside): t = (a, b, c) with e = (b, c) becomes (a, b, p), (a, p, c);
   !< the triangle across, (d, c, b), becomes (d, c, p), (d, p, b).
   type(mesh_t), intent(inout) :: m
   integer(I4P), intent(in)    :: t, e, p
   integer(I4P)                :: a, b, c, d, nca, nab, u, f, ubd, udc, t2, u2
   logical                     :: cbc, cca, cab, cbd, cdc

   a = m%v(e, t) ; b = m%v(mod(e, 3_I4P) + 1_I4P, t) ; c = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
   u   = m%nb(e, t)
   nca = m%nb(mod(e, 3_I4P) + 1_I4P, t)        ; cca = m%c(mod(e, 3_I4P) + 1_I4P, t)
   nab = m%nb(mod(e + 1_I4P, 3_I4P) + 1_I4P, t) ; cab = m%c(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
   cbc = m%c(e, t)
   t2 = m%n + 1_I4P ; m%n = m%n + 1_I4P
   if (u == 0_I4P) then
      m%v(:, t)  = [a, b, p] ; m%nb(:, t)  = [0_I4P, t2, nab] ; m%c(:, t)  = [cbc, .false., cab]
      m%v(:, t2) = [a, p, c] ; m%nb(:, t2) = [0_I4P, nca, t]  ; m%c(:, t2) = [cbc, cca, .false.]
      if (nca /= 0_I4P) call repoint(m, nca, t, t2)
      return
   endif
   f = opposite(m, u, b, c)
   d = m%v(f, u)
   ubd = m%nb(mod(f, 3_I4P) + 1_I4P, u)        ; cbd = m%c(mod(f, 3_I4P) + 1_I4P, u)
   udc = m%nb(mod(f + 1_I4P, 3_I4P) + 1_I4P, u) ; cdc = m%c(mod(f + 1_I4P, 3_I4P) + 1_I4P, u)
   u2 = m%n + 1_I4P ; m%n = m%n + 1_I4P
   m%v(:, t)  = [a, b, p] ; m%nb(:, t)  = [u2, t2, nab] ; m%c(:, t)  = [cbc, .false., cab]
   m%v(:, t2) = [a, p, c] ; m%nb(:, t2) = [u, nca, t]   ; m%c(:, t2) = [cbc, cca, .false.]
   m%v(:, u)  = [d, c, p] ; m%nb(:, u)  = [t2, u2, udc] ; m%c(:, u)  = [cbc, .false., cdc]
   m%v(:, u2) = [d, p, b] ; m%nb(:, u2) = [t, ubd, u]   ; m%c(:, u2) = [cbc, cbd, .false.]
   if (nca /= 0_I4P) call repoint(m, nca, t, t2)
   if (ubd /= 0_I4P) call repoint(m, ubd, u, u2)
   endsubroutine split_edge

   pure function opposite(m, t, b, c) result(f)
   !< Local index of the vertex of t that is neither b nor c.
   type(mesh_t), intent(in) :: m
   integer(I4P), intent(in) :: t, b, c
   integer(I4P)             :: f

   do f = 1, 3
      if (m%v(f, t) /= b .and. m%v(f, t) /= c) return
   enddo
   f = 0_I4P
   endfunction opposite

   pure subroutine repoint(m, side, old, new)
   !< The neighbour pointer of triangle `side` to `old` now points to `new`.
   type(mesh_t), intent(inout) :: m
   integer(I4P), intent(in)    :: side, old, new
   integer(I4P)                :: k

   do k = 1, 3
      if (m%nb(k, side) == old) then
         m%nb(k, side) = new ; return
      endif
   enddo
   endsubroutine repoint

   ! constraints ------------------------------------------------------------------------------------------------------
   subroutine split_constraints(pts, drop, s0, xy, segs, sub, nsub)
   !< Every constraint replaced by the chain of its pieces between the points lying on it (exactly collinear and
   !< strictly between, ordered by exact coordinate comparison), duplicates dropped. Points whose rounded coordinates
   !< are far outside the segment's box are skipped before the exact tests (the margin is far above the rounding).
   type(implicit_point_t),    intent(in)  :: pts(:)
   integer(I4P),              intent(in)  :: drop, s0
   real(R8P),                 intent(in)  :: xy(:, :)
   integer(I4P),              intent(in)  :: segs(:, :)
   integer(I4P), allocatable, intent(out) :: sub(:, :)
   integer(I4P),              intent(out) :: nsub
   integer(I4P), allocatable              :: on(:)
   integer(I4P)                           :: s, a, b, c, n_on, ax, sab, i, j, it, n
   real(R8P)                              :: margin, lo(2), hi(2)

   n = size(pts)
   allocate(on(n), sub(2, max(1_I4P, size(segs, dim=2) * 4_I4P)))
   nsub = 0_I4P
   margin = 1.e-9_R8P * max(1.e-300_R8P, maxval(abs(xy)))
   do s = 1, size(segs, dim=2)
      a = segs(1, s) ; b = segs(2, s)
      if (a == b) cycle
      ax = mod(drop, 3_I4P) + 1_I4P
      sab = ip_compare(pts(a), pts(b), ax)
      if (sab == 0_I4P) then
         ax = mod(drop + 1_I4P, 3_I4P) + 1_I4P
         sab = ip_compare(pts(a), pts(b), ax)
      endif
      lo = min(xy(:, a), xy(:, b)) - margin ; hi = max(xy(:, a), xy(:, b)) + margin
      n_on = 0_I4P
      do c = 1, n
         if (c == a .or. c == b) cycle
         if (any(xy(:, c) < lo) .or. any(xy(:, c) > hi)) cycle
         if (o2(pts, drop, s0, a, b, c) /= 0_I4P) cycle
         if (ip_compare(pts(a), pts(c), ax) /= sab .or. ip_compare(pts(c), pts(b), ax) /= sab) cycle
         n_on = n_on + 1_I4P ; on(n_on) = c
      enddo
      do i = 2, n_on ! insertion sort from a to b
         j = i
         do while (j > 1)
            if (ip_compare(pts(on(j - 1)), pts(on(j)), ax) == sab) exit
            it = on(j) ; on(j) = on(j - 1) ; on(j - 1) = it
            j = j - 1
         enddo
      enddo
      if (n_on == 0_I4P) then
         call emit(a, b)
      else
         call emit(a, on(1))
         do i = 1, n_on - 1
            call emit(on(i), on(i + 1))
         enddo
         call emit(on(n_on), b)
      endif
   enddo
   contains
      subroutine emit(p1, p2)
      integer(I4P), intent(in)  :: p1, p2
      integer(I4P), allocatable :: tmp(:, :)
      integer(I4P)              :: k

      do k = 1, nsub
         if ((sub(1, k) == p1 .and. sub(2, k) == p2) .or. (sub(1, k) == p2 .and. sub(2, k) == p1)) return
      enddo
      if (nsub == size(sub, dim=2)) then
         allocate(tmp(2, 2 * nsub)) ; tmp(:, 1:nsub) = sub(:, 1:nsub) ; call move_alloc(tmp, sub)
      endif
      nsub = nsub + 1_I4P ; sub(:, nsub) = [p1, p2]
      endsubroutine emit
   endsubroutine split_constraints

   subroutine recover(m, pts, drop, s0, a, b, status)
   !< Make (a, b) an edge and mark it constrained (Sloan 1993): the edges crossing it are flipped when their
   !< quadrilateral is convex, else put back in the queue; a flipped edge still crossing goes back too.
   type(mesh_t),           intent(inout) :: m
   type(implicit_point_t), intent(in)    :: pts(:)
   integer(I4P),           intent(in)    :: drop, s0, a, b
   integer(I4P),           intent(out)   :: status
   integer(I4P), allocatable             :: q(:, :)
   integer(I4P)                          :: head, tail, t, e, u, w, x, y, iter, cap

   status = CDTI_OK
   call find_edge(m, a, b, t, e)
   if (t /= 0_I4P) then
      call mark(m, t, e) ; return
   endif
   cap = 3_I4P * m%n + 8_I4P
   allocate(q(2, cap))
   head = 1_I4P ; tail = 0_I4P
   do t = 1, m%n
      do e = 1, 3
         if (m%nb(e, t) == 0_I4P .or. m%nb(e, t) < t) cycle
         u = m%v(mod(e, 3_I4P) + 1_I4P, t) ; w = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
         if (crosses(u, w)) call push(u, w)
      enddo
   enddo
   iter = 0_I4P
   do while (head <= tail)
      iter = iter + 1_I4P
      if (iter > 100000_I4P) then
         status = CDTI_RECOVERY ; return
      endif
      u = q(1, head) ; w = q(2, head) ; head = head + 1_I4P
      call find_edge(m, u, w, t, e)
      if (t == 0_I4P) then
         status = CDTI_RECOVERY ; return
      endif
      if (m%c(e, t)) then
         status = CDTI_CROSSING ; return
      endif
      if (.not. convex(m, pts, drop, s0, t, e)) then
         call push(u, w) ; cycle
      endif
      x = m%v(e, t)
      call flip(m, t, e)
      y = m%v(3, t) ! after the flip t = (x, b, d): the new diagonal is (x, d)
      if (crosses(x, y)) call push(x, y)
   enddo
   call find_edge(m, a, b, t, e)
   if (t == 0_I4P) then
      status = CDTI_RECOVERY ; return
   endif
   call mark(m, t, e)
   contains
      function crosses(p1, p2) result(yes)
      integer(I4P), intent(in) :: p1, p2
      logical                  :: yes

      yes = .false.
      if (p1 == a .or. p1 == b .or. p2 == a .or. p2 == b) return
      if (o2(pts, drop, s0, a, b, p1) * o2(pts, drop, s0, a, b, p2) >= 0_I4P) return
      yes = o2(pts, drop, s0, p1, p2, a) * o2(pts, drop, s0, p1, p2, b) < 0_I4P
      endfunction crosses

      subroutine push(p1, p2)
      integer(I4P), intent(in)  :: p1, p2
      integer(I4P), allocatable :: tmp(:, :)

      if (tail == size(q, dim=2)) then
         allocate(tmp(2, 2 * size(q, dim=2))) ; tmp(:, 1:tail) = q(:, 1:tail) ; call move_alloc(tmp, q)
      endif
      tail = tail + 1_I4P ; q(:, tail) = [p1, p2]
      endsubroutine push
   endsubroutine recover

   pure subroutine find_edge(m, a, b, t, e)
   !< A triangle t with edge e = (a, b), or t = 0.
   type(mesh_t), intent(in)  :: m
   integer(I4P), intent(in)  :: a, b
   integer(I4P), intent(out) :: t, e
   integer(I4P)              :: u, w

   do t = 1, m%n
      do e = 1, 3
         u = m%v(mod(e, 3_I4P) + 1_I4P, t) ; w = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
         if ((u == a .and. w == b) .or. (u == b .and. w == a)) return
      enddo
   enddo
   t = 0_I4P ; e = 0_I4P
   endsubroutine find_edge

   pure subroutine mark(m, t, e)
   !< Edge e of t constrained, on both sides.
   type(mesh_t), intent(inout) :: m
   integer(I4P), intent(in)    :: t, e
   integer(I4P)                :: t2, k

   m%c(e, t) = .true.
   t2 = m%nb(e, t)
   if (t2 == 0_I4P) return
   do k = 1, 3
      if (m%nb(k, t2) == t) m%c(k, t2) = .true.
   enddo
   endsubroutine mark

   function convex(m, pts, drop, s0, t, e) result(yes)
   !< The quadrilateral of t and its neighbour across edge e is strictly convex (exact).
   type(mesh_t),           intent(in) :: m
   type(implicit_point_t), intent(in) :: pts(:)
   integer(I4P),           intent(in) :: drop, s0, t, e
   logical                            :: yes
   integer(I4P)                       :: t2, a, b, c, d

   yes = .false.
   t2 = m%nb(e, t)
   if (t2 == 0_I4P) return
   a = m%v(e, t) ; b = m%v(mod(e, 3_I4P) + 1_I4P, t) ; c = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
   d = m%v(opposite(m, t2, b, c), t2)
   yes = o2(pts, drop, s0, a, d, b) * o2(pts, drop, s0, a, d, c) < 0_I4P
   endfunction convex

   subroutine flip(m, t, e)
   !< Flip edge e = (b, c) of t = (a, b, c) with the neighbour (d, c, b): t becomes (a, b, d), the neighbour (a, d, c).
   type(mesh_t), intent(inout) :: m
   integer(I4P), intent(in)    :: t, e
   integer(I4P)                :: t2, f, a, b, c, d, nca, nab, nbd, ndc
   logical                     :: cca, cab, cbd, cdc

   t2 = m%nb(e, t)
   a = m%v(e, t) ; b = m%v(mod(e, 3_I4P) + 1_I4P, t) ; c = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
   nca = m%nb(mod(e, 3_I4P) + 1_I4P, t)        ; cca = m%c(mod(e, 3_I4P) + 1_I4P, t)
   nab = m%nb(mod(e + 1_I4P, 3_I4P) + 1_I4P, t) ; cab = m%c(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
   f = opposite(m, t2, b, c)
   d = m%v(f, t2)
   ! t2 = (d, c, b): edge opposite c is (b, d), opposite b is (d, c)
   nbd = m%nb(mod(f, 3_I4P) + 1_I4P, t2)        ; cbd = m%c(mod(f, 3_I4P) + 1_I4P, t2)
   ndc = m%nb(mod(f + 1_I4P, 3_I4P) + 1_I4P, t2) ; cdc = m%c(mod(f + 1_I4P, 3_I4P) + 1_I4P, t2)
   m%v(:, t)  = [a, b, d] ; m%nb(:, t)  = [nbd, t2, nab] ; m%c(:, t)  = [cbd, .false., cab]
   m%v(:, t2) = [a, d, c] ; m%nb(:, t2) = [ndc, nca, t]  ; m%c(:, t2) = [cdc, cca, .false.]
   if (nbd /= 0_I4P) call repoint(m, nbd, t2, t)
   if (nca /= 0_I4P) call repoint(m, nca, t, t2)
   endsubroutine flip

   ! quality ----------------------------------------------------------------------------------------------------------
   subroutine lawson(m, pts, drop, s0, xy)
   !< Lawson flips of the unconstrained interior edges whose quadrilateral is exactly convex and whose opposite apex
   !< is inside the circumcircle (rounded coordinates: a quality choice, never a validity one).
   type(mesh_t),           intent(inout) :: m
   type(implicit_point_t), intent(in)    :: pts(:)
   integer(I4P),           intent(in)    :: drop, s0
   real(R8P),              intent(in)    :: xy(:, :)
   integer(I4P)                          :: sweep, t, e, t2, a, b, c, d
   logical                               :: flipped
   real(R8P)                             :: ic

   do sweep = 1, MAX_SWEEPS
      flipped = .false.
      do t = 1, m%n
         do e = 1, 3
            if (m%c(e, t)) cycle
            t2 = m%nb(e, t)
            if (t2 == 0_I4P .or. t2 < t) cycle
            a = m%v(e, t) ; b = m%v(mod(e, 3_I4P) + 1_I4P, t) ; c = m%v(mod(e + 1_I4P, 3_I4P) + 1_I4P, t)
            d = m%v(opposite(m, t2, b, c), t2)
            ic = incircle(a=v2(a), b=v2(b), c=v2(c), d=v2(d))
            if (real(s0, R8P) * ic <= 0._R8P) cycle
            if (.not. convex(m, pts, drop, s0, t, e)) cycle
            call flip(m, t, e)
            flipped = .true.
         enddo
      enddo
      if (.not. flipped) exit
   enddo
   contains
      pure function v2(k) result(r)
      integer(I4P), intent(in) :: k
      type(vector_R8P)         :: r

      r%x = xy(1, k) ; r%y = xy(2, k) ; r%z = 0._R8P
      endfunction v2
   endsubroutine lawson
endmodule fossil_cdt_implicit
