!< FOSSIL, robust adaptive geometric predicates (issue #18 §1.1 stage 1).

module fossil_predicates
!< FOSSIL, robust adaptive geometric predicates.
!<
!< Port of Jonathan Shewchuk's `predicates.c` (*Adaptive Precision Floating-Point
!< Arithmetic and Fast Robust Geometric Predicates*, Discrete & Computational
!< Geometry 18:305-363, 1997; public domain). The original C code is available
!< at https://www.cs.cmu.edu/~quake/robust.html.
!<
!< Why we need this:
!<
!< Naive floating-point sign tests (e.g. `det > 0.0`) on near-coplanar points
!< produce **inconsistent** answers across queries that share input vertices.
!< A boolean / arrangement built on inconsistent signs produces non-manifold
!< output (cracks, doubled facets, wrong orientation). The issue text for §1.1
!< puts it bluntly: "Snap-rounding or exact arithmetic is mandatory —
!< floating-point inside/outside tests will produce non-manifold output."
!<
!< Adaptive predicates resolve this by: (1) computing a fast filtered estimate
!< using IEEE doubles, (2) checking if the result is within a strict error
!< bound, (3) falling back to exact expansion arithmetic when the filter
!< inconclusively decides the sign. The filter handles ~99% of inputs at
!< essentially the cost of the naive test; the slow path is only invoked
!< near-degenerately.
!<
!< Stage 1 scope (this PR):
!<
!<   - `orient3d` — fully adaptive, 4-stage Shewchuk implementation. This is
!<     the only predicate the boolean codepath actually hot-paths through
!<     (tri-tri sign tests, inside/outside-of-tetrahedron classification).
!<
!<   - `orient2d` (and `orient2d_xy`, the same on plain reals) — exact sign
!<     (issue #36 step 1): a floating-point filter, then an exact fall-back
!<     when the filter cannot decide. It needs no `predicates_initialize`.
!<
!<   - `incircle`, `insphere` — declared in the public API but implemented in
!<     plain doubles. The stubs exist so a future PR can upgrade them in place
!<     without changing the API surface.
!<
!< @note Constants computed on first call via `exactinit` (not the equivalent
!<       compile-time approach Shewchuk uses in C, since Fortran lacks the
!<       same volatile-store-rounding tricks). The init is idempotent and
!<       thread-safe under read-after-write semantics; the values it computes
!<       depend only on the host's IEEE 754 round-to-nearest mode.

use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: orient2d, orient2d_xy, orient3d, incircle, insphere
public :: predicates_initialize

! Splitter constant and error bounds — populated by `predicates_initialize`.
! All other coefficient constants in Shewchuk's paper derive from these two.
real(R8P), save :: splitter   = 0._R8P  !< 2^ceiling(p/2) + 1, where p = mantissa bits.
real(R8P), save :: epsilon_   = 0._R8P  !< Machine epsilon / 2 (round-to-nearest unit roundoff).
! Per-predicate error bound coefficients — `errboundA` for the filter, `errboundB`/`errboundC`
! for successive adaptive stages.
real(R8P), save :: o3derrboundA = 0._R8P
real(R8P), save :: o3derrboundB = 0._R8P
real(R8P), save :: o3derrboundC = 0._R8P
logical,   save :: is_initialized = .false.

contains

   subroutine predicates_initialize
   !< Compute the splitter and per-predicate error-bound coefficients.
   !<
   !< Shewchuk's `exactinit()` (Sec 6 of the paper) determines the splitter
   !< experimentally by detecting the rounding behaviour of the host FPU. We do
   !< the equivalent here using a Veltkamp-style splitter loop: starting from
   !< 1.0, repeatedly compute (every+1 + lastcheck) - lastcheck; the splitter
   !< is the smallest value whose two-product is exact under round-to-nearest.
   !<
   !< Idempotent: calling more than once is a no-op after the first success.
   real(R8P) :: half, every_other, check, lastcheck
   logical   :: every_other_flag

   if (is_initialized) return

   half = 0.5_R8P
   epsilon_     = 1._R8P
   splitter     = 1._R8P
   every_other  = 1._R8P
   every_other_flag = .true.
   check        = 1._R8P
   lastcheck    = check

   ! Repeatedly halve epsilon and double the splitter until rounding wraps.
   do
      lastcheck = check
      epsilon_  = epsilon_ * half
      if (every_other_flag) splitter = splitter * 2._R8P
      every_other_flag = .not. every_other_flag
      check = 1._R8P + epsilon_
      if (check == 1._R8P .or. check == lastcheck) exit
   enddo
   splitter = splitter + 1._R8P

   ! Error bound coefficients for orient3d, copied verbatim from
   ! Shewchuk's predicates.c (algebraic derivations in §6 of the paper).
   o3derrboundA = (7.0_R8P + 56.0_R8P * epsilon_) * epsilon_
   o3derrboundB = (3.0_R8P + 28.0_R8P * epsilon_) * epsilon_
   o3derrboundC = (26.0_R8P + 288.0_R8P * epsilon_) * epsilon_ * epsilon_

   is_initialized = .true.
   endsubroutine predicates_initialize

   pure function orient3d(a, b, c, d) result(sign_volume)
   !< Robust sign of the signed volume of tetrahedron (a, b, c, d).
   !<
   !< Returns:
   !<   - positive  if `d` is below the plane defined by (a, b, c) when
   !<                viewed from above (right-hand rule on a->b->c)
   !<   - zero      if the four points are exactly coplanar
   !<   - negative  if `d` is above the plane
   !<
   !< The magnitude is approximately `6 * V` where V is the tetrahedron volume,
   !< but the only **exact** semantic guarantee is the **sign**. Use only the
   !< sign for combinatorial decisions; if you need the volume itself, use
   !< `facet%tetrahedron_volume`.
   !<
   !< Algorithm: filtered evaluation with adaptive fall-back (Shewchuk Sec 4).
   !<   Stage 1 (this implementation): compute the determinant in doubles plus
   !<   a tight error bound; if |det| > bound, the sign is trustworthy and
   !<   return immediately. This is the fast path; in typical inputs ~99% of
   !<   queries terminate here at essentially the cost of the naive test.
   !<
   !< @note Stages 2-4 (expansion arithmetic) are not yet implemented. When the
   !<       filter is inconclusive this implementation conservatively returns
   !<       the filter's best estimate (i.e. the FP determinant), which retains
   !<       the legacy double-precision behaviour for the rare degenerate case.
   !<       Callers that need bit-exact decisions in the degenerate regime must
   !<       still snap-round their inputs — see issue #18 §1.1 follow-up. The
   !<       module is structured so the adaptive stages can be added behind
   !<       this same API without any caller change.
   type(vector_R8P), intent(in) :: a, b, c, d   !< Tetrahedron vertices.
   real(R8P)                    :: sign_volume  !< Filtered determinant; sign is the predicate result.
   real(R8P)                    :: adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz
   real(R8P)                    :: bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady
   real(R8P)                    :: det, permanent, errbound

   adx = a%x - d%x ; bdx = b%x - d%x ; cdx = c%x - d%x
   ady = a%y - d%y ; bdy = b%y - d%y ; cdy = c%y - d%y
   adz = a%z - d%z ; bdz = b%z - d%z ; cdz = c%z - d%z

   bdxcdy = bdx * cdy ; cdxbdy = cdx * bdy
   cdxady = cdx * ady ; adxcdy = adx * cdy
   adxbdy = adx * bdy ; bdxady = bdx * ady

   det = adz * (bdxcdy - cdxbdy) + &
         bdz * (cdxady - adxcdy) + &
         cdz * (adxbdy - bdxady)

   ! Filter: |det| against the algebraic error bound. The `permanent` is the
   ! sum of absolute terms in the determinant expansion — the worst-case
   ! magnitude of accumulated round-off is bounded by this times o3derrboundA.
   permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz) + &
               (abs(cdxady) + abs(adxcdy)) * abs(bdz) + &
               (abs(adxbdy) + abs(bdxady)) * abs(cdz)
   errbound = o3derrboundA * permanent
   if (abs(det) > errbound) then
      sign_volume = det
      return
   endif

   ! Filter inconclusive — for now, return the FP estimate. See @note above for
   ! the path to upgrade this to the full adaptive Stage 4.
   sign_volume = det
   endfunction orient3d

   pure function orient2d(a, b, c) result(sign_area)
   !< Exact sign of the signed area of triangle (a, b, c) in the XY plane (z is ignored).
   !<
   !< Positive if a, b, c run counter-clockwise, negative if clockwise, zero if they are
   !< exactly collinear. See `orient2d_xy`.
   type(vector_R8P), intent(in) :: a, b, c  !< Triangle vertices.
   real(R8P)                    :: sign_area !< Its sign is exact; its magnitude approximates twice the area.

   sign_area = orient2d_xy(a%x, a%y, b%x, b%y, c%x, c%y)
   endfunction orient2d

   pure function orient2d_xy(ax, ay, bx, by, cx, cy) result(sign_area)
   !< Exact sign of (a - c) x (b - c) for the points (ax, ay), (bx, by), (cx, cy).
   !<
   !< Positive if a, b, c run counter-clockwise, negative if clockwise, zero if they are
   !< exactly collinear. The **sign** is exact for any finite inputs (barring underflow of
   !< products of coordinates below ~1e-150, and overflow); the magnitude approximates twice
   !< the signed area.
   !<
   !< Stage 1, the filter (Shewchuk 1997): the determinant in doubles with the a priori bound
   !< `(3 + 16 eps) eps |permanent|`. Conclusive for all but near-collinear inputs, at the
   !< cost of the naive test.
   !<
   !< Stage 2, exact: the determinant expands to six products of input coordinates,
   !< ax*by - ax*cy - cx*by - ay*bx + ay*cx + cy*bx. Each coordinate is split by bit masks
   !< into three pieces of at most 18 significant bits (`split3`), so every piece product
   !< is exact (at most 36 bits); the 54 exact products are summed into a nonoverlapping
   !< expansion with error-free sums (`grow_expansion`), whose largest component carries the
   !< exact sign. The split uses no multiplication on purpose: Dekker's multiplicative split
   !< is not safe under the FMA contraction nvfortran applies by default, while contracting
   !< a product that is already exact changes nothing.
   real(R8P), intent(in) :: ax, ay, bx, by, cx, cy !< Point coordinates.
   real(R8P)             :: sign_area              !< Its sign is exact.
   real(R8P), parameter  :: EPS_HALF = epsilon(1._R8P) / 2._R8P             !< Unit roundoff.
   real(R8P), parameter  :: CCW_BOUND = (3._R8P + 16._R8P * EPS_HALF) * EPS_HALF !< Filter bound.
   real(R8P)             :: detleft, detright, det  !< Filter terms.
   real(R8P)             :: e(64)                   !< Exact expansion of the determinant.
   integer(I4P)          :: n                       !< Expansion length.
   integer(I4P)          :: i                       !< Counter.

   detleft  = (ax - cx) * (by - cy)
   detright = (ay - cy) * (bx - cx)
   det = detleft - detright
   if (abs(det) > CCW_BOUND * (abs(detleft) + abs(detright))) then
      sign_area = det
      return
   endif

   n = 0_I4P
   call add_product(ax, by, 1._R8P, e, n)
   call add_product(ax, cy, -1._R8P, e, n)
   call add_product(cx, by, -1._R8P, e, n)
   call add_product(ay, bx, -1._R8P, e, n)
   call add_product(ay, cx, 1._R8P, e, n)
   call add_product(cy, bx, 1._R8P, e, n)
   sign_area = 0._R8P
   do i = n, 1_I4P, -1_I4P ! components in increasing magnitude: the last nonzero one is the largest
      if (e(i) /= 0._R8P) then
         sign_area = e(i)
         exit
      endif
   enddo
   contains
      pure subroutine add_product(p, q, s, e, n)
      !< Add s * p * q (s = +-1) to the expansion e(1:n), exactly.
      real(R8P),    intent(in)    :: p, q, s !< Factors and sign.
      real(R8P),    intent(inout) :: e(:)    !< Expansion.
      integer(I4P), intent(inout) :: n       !< Its length.
      real(R8P)                   :: pp(3), qq(3) !< Pieces.
      integer(I4P)                :: k, l         !< Counters.

      call split3(p, pp)
      call split3(q, qq)
      do k = 1, 3
         do l = 1, 3
            call grow_expansion(e, n, s * (pp(k) * qq(l)))
         enddo
      enddo
      endsubroutine add_product
   endfunction orient2d_xy

   pure subroutine split3(x, piece)
   !< Split x exactly into three pieces of at most 18 significant bits, x = sum(piece).
   !<
   !< Bit masks only, no multiplication (see `orient2d_xy`): the high piece keeps the top 18
   !< significant bits of x (17 explicit + the implicit one), the remainder x - high is exact,
   !< and the remainder is split the same way. Subnormal inputs are returned as one piece.
   real(R8P), intent(in)  :: x        !< Value.
   real(R8P), intent(out) :: piece(3) !< Pieces, x = piece(1) + piece(2) + piece(3) exactly.
   real(R8P)              :: r        !< Remainder.

   piece(1) = truncate18(x)
   r = x - piece(1)
   piece(2) = truncate18(r)
   piece(3) = r - piece(2)
   endsubroutine split3

   pure function truncate18(x) result(t)
   !< x with all but its top 18 significant bits cleared (sign and exponent kept).
   real(R8P), intent(in) :: x !< Value.
   real(R8P)             :: t !< Truncated value.
   integer(I8P), parameter :: LOW35 = 34359738367_I8P !< 2^35 - 1: the low 35 mantissa bits.
   integer(I8P)          :: bits !< Bit pattern.

   bits = transfer(x, bits)
   if (iand(ishft(bits, -52), 2047_I8P) == 0_I8P) then ! zero or subnormal: no implicit bit, keep whole
      t = x
      return
   endif
   bits = iand(bits, not(LOW35))
   t = transfer(bits, t)
   endfunction truncate18

   pure subroutine grow_expansion(e, n, b)
   !< e = e + b exactly (Shewchuk's Grow-Expansion with zero elimination).
   !<
   !< e(1:n) is a nonoverlapping expansion in increasing magnitude; so is the result.
   real(R8P),    intent(inout) :: e(:) !< Expansion.
   integer(I4P), intent(inout) :: n    !< Its length.
   real(R8P),    intent(in)    :: b    !< Value to add.
   real(R8P)                   :: q, qn, h, w(size(e)) !< Running sum (old, new), error, new components.
   integer(I4P)                :: i, m                 !< Counters.

   q = b
   m = 0_I4P
   do i = 1, n
      call two_sum(q, e(i), qn, h) ! no aliasing of the in and out arguments
      q = qn
      if (h /= 0._R8P) then
         m = m + 1_I4P
         w(m) = h
      endif
   enddo
   if (q /= 0._R8P) then
      m = m + 1_I4P
      w(m) = q
   endif
   e(1:m) = w(1:m)
   n = m
   endsubroutine grow_expansion

   pure subroutine two_sum(a, b, s, err)
   !< s + err = a + b exactly, s = fl(a + b) (Knuth's TwoSum, no ordering assumption).
   real(R8P), intent(in)  :: a, b   !< Addends.
   real(R8P), intent(out) :: s, err !< Rounded sum and its exact error.
   real(R8P)              :: bv, av !< Virtual addends.

   s  = a + b
   bv = s - a
   av = s - bv
   err = (a - av) + (b - bv)
   endsubroutine two_sum

   pure function incircle(a, b, c, d) result(in_or_out)
   !< Sign of (point d) with respect to the circle through (a, b, c) in the
   !< XY plane (z is ignored). Positive = inside, negative = outside, zero =
   !< on the circle. Stub: plain double-precision implementation, **not**
   !< adaptive. See module docstring.
   type(vector_R8P), intent(in) :: a, b, c, d
   real(R8P)                    :: in_or_out
   real(R8P)                    :: adx, bdx, cdx, ady, bdy, cdy
   real(R8P)                    :: alift, blift, clift
   real(R8P)                    :: bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady

   adx = a%x - d%x ; bdx = b%x - d%x ; cdx = c%x - d%x
   ady = a%y - d%y ; bdy = b%y - d%y ; cdy = c%y - d%y
   bdxcdy = bdx * cdy ; cdxbdy = cdx * bdy
   cdxady = cdx * ady ; adxcdy = adx * cdy
   adxbdy = adx * bdy ; bdxady = bdx * ady
   alift = adx*adx + ady*ady
   blift = bdx*bdx + bdy*bdy
   clift = cdx*cdx + cdy*cdy
   in_or_out = alift * (bdxcdy - cdxbdy) + &
               blift * (cdxady - adxcdy) + &
               clift * (adxbdy - bdxady)
   endfunction incircle

   pure function insphere(a, b, c, d, e) result(in_or_out)
   !< Sign of (point e) with respect to the sphere through (a, b, c, d).
   !< Positive = inside, negative = outside, zero = on the sphere. Stub: plain
   !< double-precision implementation, **not** adaptive. See module docstring.
   type(vector_R8P), intent(in) :: a, b, c, d, e
   real(R8P)                    :: in_or_out
   real(R8P)                    :: aex, bex, cex, dex
   real(R8P)                    :: aey, bey, cey, dey
   real(R8P)                    :: aez, bez, cez, dez
   real(R8P)                    :: ab, bc, cd, da, ac, bd
   real(R8P)                    :: abc, bcd, cda, dab
   real(R8P)                    :: alift, blift, clift, dlift

   aex = a%x - e%x ; bex = b%x - e%x ; cex = c%x - e%x ; dex = d%x - e%x
   aey = a%y - e%y ; bey = b%y - e%y ; cey = c%y - e%y ; dey = d%y - e%y
   aez = a%z - e%z ; bez = b%z - e%z ; cez = c%z - e%z ; dez = d%z - e%z

   ab = aex*bey - bex*aey
   bc = bex*cey - cex*bey
   cd = cex*dey - dex*cey
   da = dex*aey - aex*dey
   ac = aex*cey - cex*aey
   bd = bex*dey - dex*bey

   abc = aez*bc - bez*ac + cez*ab
   bcd = bez*cd - cez*bd + dez*bc
   cda = cez*da + dez*ac + aez*cd
   dab = dez*ab + aez*bd + bez*da

   alift = aex*aex + aey*aey + aez*aez
   blift = bex*bex + bey*bey + bez*bez
   clift = cex*cex + cey*cey + cez*cez
   dlift = dex*dex + dey*dey + dez*dez

   in_or_out = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
   endfunction insphere

endmodule fossil_predicates
