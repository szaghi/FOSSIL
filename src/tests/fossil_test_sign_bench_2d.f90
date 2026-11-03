!< FOSSIL test/bench: sign algorithms at 2-D corners and cusps (issue #39 step 3).
program fossil_test_sign_bench_2d
!< FOSSIL test/bench: sign algorithms at 2-D corners and cusps (issue #39 step 3).
!<
!< Curves: a square, wedges of 1, 5, 30, 90 and 150 degrees, an L-shape, a cusp (two circular arcs
!< tangent at the origin, corner angle ~0.01 degrees after discretization), a C-shape with a 1e-3
!< slot, a square annulus (outer loop and hole), and the sections at mid-height of the L-prism and
!< of the 5-degree wedge of `fossil_sign_bench_bodies` (`slice_surface`).
!<
!< Features are generic: every vertex (a corner of interior angle theta, measured counter-clockwise
!< from the outgoing to the incoming segment through the solid) and the mid-point of every segment.
!< Query points x + delta d, d at theta k / 6 (inside) and theta + (360 - theta) k / 6 (outside),
!< k = 1..5, or -n / +n on segments; delta = 1e-1 ... 1e-12 L, at most a quarter of the distance
!< to the nearest segment not containing x. The side is known by construction, and the exact
!< winding number (`winding_number_exact`, exact `orient2d`) must agree with it on every clean point.
!<
!< Algorithms: pseudo-normal (`is_point_inside`), fast winding number (SIGN_FAST_WINDING_NUMBER), and
!< the exact winding number. Dirty variants of the first loop (ground truth the clean curve): a gap
!< (the closing segment removed: open loop), a duplicated vertex (a zero-length segment), a
!< micro-segment (a vertex at 1e-9 of the first segment).
!<
!< Asserts: the exact winding number equals the construction on every clean point; every algorithm
!< right at each feature's largest delta; errors within the upper bounds of the first run.

use fossil,                   only : curve_object, surface_stl_object, slice_surface, SIGN_FAST_WINDING_NUMBER
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P

implicit none

integer, parameter :: NA = 3, NK = 5, ND = 12, NV = 3, NC = 12 !< Algorithms, kinds, decades, dirty variants, curves.
character(14), parameter :: ANAME(NA) = ['pseudo-normal ', 'fast winding  ', 'exact winding ']
character(16), parameter :: KNAME(NK) = ['segment         ', 'flat vertex     ', 'convex corner   ', 'acute corner    ', &
                                         'reentrant corner']
character(16), parameter :: VNAME(NV) = ['gap (open loop) ', 'duplicate vertex', 'micro-segment   ']
! upper bounds from the first run (issue #39 step 3), never relaxed: clean and dirty errors per algorithm
integer, parameter :: BOUND_CLEAN(NA) = [0, 0, 0]
integer, parameter :: BOUND_DIRTY(NA) = [584, 5, 7893]
real(R8P), parameter :: PI = acos(-1._R8P)

type :: crv_t
   !< A benchmark curve as loops of vertices.
   character(24)          :: name = '' !< Name.
   real(R8P), allocatable :: u(:), v(:) !< Vertices.
   integer(I4P), allocatable :: sizes(:) !< Loop sizes.
endtype crv_t

type(crv_t)               :: cv(NC)
type(curve_object)        :: c, cd
real(R8P),    allocatable :: p(:,:)
logical,      allocatable :: ins(:), got(:,:)
integer,      allocatable :: knd(:), dec(:), fid(:)
integer      :: err(NA, NK, ND), tot(NK, ND), errd(NA, NV), totd(NV), i, m, a, v, k, nbig, nref, firstfail
integer(I8P) :: c0, c1, rate
real(R8P)    :: tcost(NA), diag
logical      :: ok(3)
character(ND * 9) :: line
character(8) :: tag

call build_curves(cv)
err = 0 ; tot = 0 ; errd = 0 ; totd = 0 ; tcost = 0._R8P ; nbig = 0 ; nref = 0
do i = 1, NC
   call c%adopt_loops(u=cv(i)%u, v=cv(i)%v, loop_sizes=cv(i)%sizes)
   call make_points(cv(i), p, ins, knd, dec, fid, diag)
   call classify(c, p, got, tcost)
   do m = 1, size(ins)
      tot(knd(m), dec(m)) = tot(knd(m), dec(m)) + 1
      if (got(3, m) .neqv. ins(m)) nref = nref + 1
      do a = 1, NA
         if (got(a, m) .neqv. ins(m)) then
            err(a, knd(m), dec(m)) = err(a, knd(m), dec(m)) + 1
            if (dec(m) == minval(dec, mask=fid == fid(m))) nbig = nbig + 1
         endif
      enddo
   enddo
   do v = 1, NV
      call dirty(cv(i), v, cd)
      call classify(cd, p, got)
      totd(v) = totd(v) + size(ins)
      do a = 1, NA
         errd(a, v) = errd(a, v) + count(got(a, :) .neqv. ins)
      enddo
   enddo
   print '(A,A,I0,A,I0,A,ES9.2)', cv(i)%name, ' vertices ', size(cv(i)%u), '  points ', size(ins), '  L ', diag
enddo

print '(A)', 'clean curves: errors / points per feature kind and decade delta/L = 1e-1 ... 1e-12'
do a = 1, NA
   print '(A)', trim(ANAME(a))//':'
   do k = 1, NK
      if (sum(tot(k, :)) == 0) cycle
      line = ''
      do m = 1, ND
         write(line(9 * (m - 1) + 1:9 * m), '(I4,A,I4)') err(a, k, m), '/', tot(k, m)
      enddo
      firstfail = 0
      do m = 1, ND
         if (err(a, k, m) > 0) then
            firstfail = m ; exit
         endif
      enddo
      write(tag, '(A,I0)') '1e-', firstfail ; if (firstfail == 0) tag = 'none'
      print '(3X,A16,A,A,A)', KNAME(k), line, '  first failing ', trim(tag)
   enddo
enddo
print '(A)', 'cost per query [us], clean curves:'
do a = 1, NA
   print '(3X,A14,F10.3)', ANAME(a), 1.e6_R8P * tcost(a) / max(sum(tot), 1)
enddo
print '(A)', 'dirty variants of the first loop: errors / points'
do v = 1, NV
   print '(3X,A16,3(A,A,I0,A,I0))', VNAME(v), (' | ', trim(ANAME(a))//' ', errd(a, v), '/', totd(v), a = 1, NA)
enddo

ok(1) = nref == 0
ok(2) = nbig == 0
ok(3) = .true.
do a = 1, NA
   if (BOUND_CLEAN(a) >= 0) ok(3) = ok(3) .and. sum(err(a, :, :)) <= BOUND_CLEAN(a)
   if (BOUND_DIRTY(a) >= 0) ok(3) = ok(3) .and. sum(errd(a, :)) <= BOUND_DIRTY(a)
enddo
print '(A,I0,A,L1)', 'exact winding number = construction on every clean point (', sum(tot), ' points): ', ok(1)
print '(A,L1)', 'every algorithm right at the largest delta of every feature (clean): ', ok(2)
print '(A,3I7,A,3I7,A,L1)', 'errors clean ', [(sum(err(a, :, :)), a = 1, NA)], '  dirty ', [(sum(errd(a, :)), a = 1, NA)], &
      '  within bounds: ', ok(3)
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine build_curves(cv)
   !< The benchmark curves.
   type(crv_t), intent(out) :: cv(NC) !< Curves.
   type(body_t), allocatable :: b(:)  !< 3-D bodies (for the sections).
   type(surface_stl_object)  :: s     !< Surface.
   type(curve_object)        :: sec   !< Section.
   real(R8P)                 :: x(33), uv(2) !< Work.
   integer                   :: k, j, l, i1, i2 !< Counters.

   call set(cv(1), 'square', [0._R8P, 1._R8P, 1._R8P, 0._R8P], [0._R8P, 0._R8P, 1._R8P, 1._R8P], [4])
   call wedge(cv(2), 1._R8P) ; call wedge(cv(3), 5._R8P) ; call wedge(cv(4), 30._R8P)
   call wedge(cv(5), 90._R8P) ; call wedge(cv(6), 150._R8P)
   call set(cv(7), 'L-shape', [1._R8P, 1._R8P, 0._R8P, 0._R8P, 2._R8P, 2._R8P], &
           [1._R8P, 2._R8P, 2._R8P, 0._R8P, 0._R8P, 1._R8P], [6])
   ! cusp: lower arc (circle of radius 1 tangent to the u axis at the origin), upper arc (radius 0.5),
   ! u in [0, 0.4] with quadratic spacing, so the first segments are short and the corner very sharp
   do k = 1, 33
      x(k) = 0.4_R8P * ((k - 1) / 32._R8P)**2
   enddo
   allocate(cv(8)%u(65), cv(8)%v(65)) ! the origin once: x(1) = 0 starts the lower arc only
   do k = 1, 33
      cv(8)%u(k) = x(k) ; cv(8)%v(k) = 1._R8P - sqrt(1._R8P - x(k)**2)
   enddo
   do k = 1, 32
      cv(8)%u(33 + k) = x(34 - k) ; cv(8)%v(33 + k) = 0.5_R8P - sqrt(0.25_R8P - x(34 - k)**2)
   enddo
   cv(8)%name = 'cusp' ; cv(8)%sizes = [65_I4P]
   call set(cv(9), 'C-shape slot 1e-3', [0._R8P, 0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P, 1._R8P, 1._R8P, 0._R8P], &
            [0._R8P, 0._R8P, 0.5_R8P, 0.5_R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P], [8])
   call set(cv(10), 'square annulus', [0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
            [0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], [4, 4])
   call bench_bodies(b)
   do j = 1, 2
      l = merge(8, 2, j == 1) ! L-prism, wedge 5
      call to_surface(b(l)%soup, s)
      call slice_surface(surface=s, normal_axis=3_I4P, level=0.5_R8P, curve=sec)
      call sec%get_loop_range(1_I4P, i1, i2)
      allocate(cv(10 + j)%u(i2 - i1 + 1), cv(10 + j)%v(i2 - i1 + 1))
      do k = i1, i2
         uv = sec%get_vertex(k) ; cv(10 + j)%u(k - i1 + 1) = uv(1) ; cv(10 + j)%v(k - i1 + 1) = uv(2)
      enddo
      cv(10 + j)%sizes = [i2 - i1 + 1]
      cv(10 + j)%name = 'section of '//trim(b(l)%name)
   enddo
   endsubroutine build_curves

   subroutine set(cc, name, u, v, sizes)
   !< Fill a curve.
   type(crv_t),  intent(inout) :: cc      !< Curve.
   character(*), intent(in)    :: name    !< Name.
   real(R8P),    intent(in)    :: u(:), v(:) !< Vertices.
   integer,      intent(in)    :: sizes(:) !< Loop sizes.

   cc%name = name ; cc%u = u ; cc%v = v ; cc%sizes = int(sizes, I4P)
   endsubroutine set

   subroutine wedge(cc, theta)
   !< Triangle with apex angle theta at the origin, sides of length 1, symmetric about +u.
   type(crv_t), intent(inout) :: cc    !< Curve.
   real(R8P),   intent(in)    :: theta !< Apex angle [deg].
   real(R8P)                  :: h     !< Half angle.
   character(24)              :: nm    !< Name.

   h = 0.5_R8P * theta * PI / 180._R8P
   write(nm, '(A,I0)') 'wedge ', nint(theta)
   call set(cc, trim(nm), [0._R8P, cos(h), cos(h)], [0._R8P, -sin(h), sin(h)], [3])
   endsubroutine wedge

   subroutine make_points(cc, p, inside, kind, decade, feature, diag)
   !< Generic features (every vertex, every segment mid-point) and their query points.
   type(crv_t),          intent(in)  :: cc         !< Curve.
   real(R8P), allocatable, intent(out) :: p(:,:)   !< Points [2, m].
   logical,   allocatable, intent(out) :: inside(:) !< Side by construction.
   integer,   allocatable, intent(out) :: kind(:), decade(:), feature(:) !< Kind, decade, feature id.
   real(R8P),            intent(out) :: diag       !< Bounding-box diagonal L.
   real(R8P), allocatable :: pb(:,:), fx(:,:), fdir(:,:,:)
   logical,   allocatable :: ib(:), fin(:,:)
   integer,   allocatable :: kb(:), db(:), fb(:), fk(:), fnd(:)
   integer   :: l, j, j0, n, jp, jn, nf, f, k, m, q
   real(R8P) :: u1(2), u2(2), th, phi, lfs, delta, nn(2), x(2), d

   diag = norm2([maxval(cc%u) - minval(cc%u), maxval(cc%v) - minval(cc%v)])
   nf = 2 * size(cc%u)
   allocate(fx(2, nf), fdir(2, 10, nf), fin(10, nf), fk(nf), fnd(nf))
   f = 0 ; j0 = 0
   do l = 1, size(cc%sizes)
      n = cc%sizes(l)
      do j = 1, n
         jn = j0 + mod(j, n) + 1 ; jp = j0 + mod(j + n - 2, n) + 1
         x = [cc%u(j0 + j), cc%v(j0 + j)]
         ! corner at vertex j
         u1 = [cc%u(jn), cc%v(jn)] - x ; u1 = u1 / norm2(u1)
         u2 = [cc%u(jp), cc%v(jp)] - x ; u2 = u2 / norm2(u2)
         th = atan2(u1(1) * u2(2) - u1(2) * u2(1), dot_product(u1, u2))
         if (th <= 0._R8P) th = th + 2._R8P * PI
         f = f + 1 ; fx(:, f) = x ; fnd(f) = 10
         if (abs(th - PI) < 1.e-9_R8P) then
            fk(f) = 2
         elseif (th > PI) then
            fk(f) = 5
         elseif (th <= 30._R8P * PI / 180._R8P) then
            fk(f) = 4
         else
            fk(f) = 3
         endif
         do k = 1, 5
            phi = th * k / 6._R8P
            fdir(:, k, f) = [cos(phi) * u1(1) - sin(phi) * u1(2), sin(phi) * u1(1) + cos(phi) * u1(2)] ; fin(k, f) = .true.
            phi = th + (2._R8P * PI - th) * k / 6._R8P
            fdir(:, 5 + k, f) = [cos(phi) * u1(1) - sin(phi) * u1(2), sin(phi) * u1(1) + cos(phi) * u1(2)]
            fin(5 + k, f) = .false.
         enddo
         ! mid-point of segment j -> jn: outward normal (du, dv) -> (dv, -du)
         f = f + 1 ; fx(:, f) = 0.5_R8P * (x + [cc%u(jn), cc%v(jn)]) ; fk(f) = 1 ; fnd(f) = 2
         nn = [u1(2), -u1(1)]
         fdir(:, 1, f) = -nn ; fin(1, f) = .true. ; fdir(:, 2, f) = nn ; fin(2, f) = .false.
      enddo
      j0 = j0 + n
   enddo
   allocate(pb(2, 10 * ND * f), ib(10 * ND * f), kb(10 * ND * f), db(10 * ND * f), fb(10 * ND * f))
   m = 0
   do q = 1, f
      lfs = huge(1._R8P)
      j0 = 0
      do l = 1, size(cc%sizes)
         n = cc%sizes(l)
         do j = 1, n
            jn = j0 + mod(j, n) + 1
            d = seg_dist(fx(:, q), [cc%u(j0 + j), cc%v(j0 + j)], [cc%u(jn), cc%v(jn)])
            if (d > 1.e-13_R8P * diag) lfs = min(lfs, d)
         enddo
         j0 = j0 + n
      enddo
      lfs = 0.25_R8P * lfs
      do k = 1, ND
         delta = 10._R8P**(-k) * diag
         if (delta > lfs) cycle
         do j = 1, fnd(q)
            m = m + 1 ; pb(:, m) = fx(:, q) + delta * fdir(:, j, q) ; ib(m) = fin(j, q) ; kb(m) = fk(q) ; db(m) = k ; fb(m) = q
         enddo
      enddo
   enddo
   p = pb(:, 1:m) ; inside = ib(1:m) ; kind = kb(1:m) ; decade = db(1:m) ; feature = fb(1:m)
   endsubroutine make_points

   pure function seg_dist(x, a, b) result(d)
   !< Distance from x to the segment ab.
   real(R8P), intent(in) :: x(2), a(2), b(2) !< Point, segment.
   real(R8P)             :: d                !< Distance.
   real(R8P)             :: t                !< Parameter.

   t = 0._R8P
   if (sum((b - a)**2) > 0._R8P) t = max(0._R8P, min(1._R8P, dot_product(x - a, b - a) / sum((b - a)**2)))
   d = norm2(x - a - t * (b - a))
   endfunction seg_dist

   subroutine dirty(cc, v, out)
   !< Dirty variant v of the first loop of cc.
   type(crv_t),        intent(in)    :: cc  !< Clean curve.
   integer,            intent(in)    :: v   !< Variant.
   type(curve_object), intent(inout) :: out !< Dirty curve.
   real(R8P), allocatable :: u(:), w(:) !< Vertices.
   integer(I4P), allocatable :: sz(:)  !< Sizes.
   logical,   allocatable :: cl(:)    !< Closed flags.

   u = cc%u ; w = cc%v ; sz = cc%sizes ; allocate(cl(size(sz))) ; cl = .true.
   select case (v)
   case (1) ; cl(1) = .false.
   case (2) ; u = [u(1:2), u(2:)] ; w = [w(1:2), w(2:)] ; sz(1) = sz(1) + 1
   case (3)
      u = [u(1), u(1) + 1.e-9_R8P * (u(2) - u(1)), u(2:)] ; w = [w(1), w(1) + 1.e-9_R8P * (w(2) - w(1)), w(2:)]
      sz(1) = sz(1) + 1
   endselect
   call out%adopt_loops(u=u, v=w, loop_sizes=sz, closed=cl)
   endsubroutine dirty

   subroutine classify(cur, pts, inside, cost)
   !< Inside flags of every point by every algorithm; optionally accumulate the cost.
   type(curve_object), intent(in)              :: cur         !< Curve.
   real(R8P),          intent(in)              :: pts(:,:)    !< Points [2, n].
   logical, allocatable, intent(out)           :: inside(:,:) !< Flags [NA, n].
   real(R8P),          intent(inout), optional :: cost(NA)    !< Accumulated seconds.
   integer :: mm !< Counter.

   allocate(inside(NA, size(pts, 2)))
   call system_clock(c0, rate)
   do mm = 1, size(pts, 2)
      inside(1, mm) = cur%is_point_inside(pts(:, mm))
   enddo
   call system_clock(c1) ; if (present(cost)) cost(1) = cost(1) + real(c1 - c0, R8P) / rate
   call system_clock(c0)
   do mm = 1, size(pts, 2)
      inside(2, mm) = cur%is_point_inside(pts(:, mm), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   enddo
   call system_clock(c1) ; if (present(cost)) cost(2) = cost(2) + real(c1 - c0, R8P) / rate
   call system_clock(c0)
   do mm = 1, size(pts, 2)
      inside(3, mm) = cur%winding_number_exact(pts(:, mm)) /= 0
   enddo
   call system_clock(c1) ; if (present(cost)) cost(3) = cost(3) + real(c1 - c0, R8P) / rate
   endsubroutine classify
endprogram fossil_test_sign_bench_2d
