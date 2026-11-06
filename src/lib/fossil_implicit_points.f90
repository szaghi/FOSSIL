!< FOSSIL, implicit points and exact predicates on them (issue #65).
module fossil_implicit_points
!< FOSSIL, implicit points and exact predicates on them (issue #65).
!<
!< An intersection point of the boolean arrangement is kept as its construction, not rounded when it is made
!< (Attene 2020, "Indirect predicates for geometric constructions"; Cherchi et al. 2020, "Fast and robust mesh
!< arrangements using floating-point arithmetic"):
!< - IPK_EXPLICIT: an input vertex x;
!< - IPK_LPI: the line through p, q meets the plane through r, s, t;
!< - IPK_SSI: segment (a, b) meets segment (c, d), both in one plane, solved in the coordinate plane that drops axis
!<   `drop`: the point of the line (a, b) whose projection lies on the projected line (c, d).
!< Each point has homogeneous coordinates (X, W), polynomials in the input coordinates:
!< - explicit: X = x, W = 1;
!< - LPI: n = (s - r) x (t - r), W = n . (q - p), X = W p + (n . (r - p)) (q - p);
!< - SSI: with u = b - a, v = d - c, w = c - a in the kept axes (i, j): W = u_i v_j - u_j v_i,
!<   X = W a + (w_i v_j - w_j v_i) u.
!< The defining inputs are put in a canonical order first (endpoints and segments lexicographically, the plane's
!< three points lexicographically), so one geometric construction is one point, and `ip_approx` rounds it by one
!< formula: equal identities give equal bits.
!<
!< Predicates (`ip_orient3d`, `ip_orient2d`, `ip_compare`) return the exact sign of a polynomial in the input
!< coordinates. Stage 1 evaluates it in doubles with a running bound of the rounding error of every operation (a
!< dynamic filter, conservative); when the bound does not settle the sign, stage 2 evaluates the same polynomial in
!< expansion arithmetic (nonoverlapping sums of doubles, Shewchuk 1997), exactly. Products are made exact by
!< splitting each factor by bit masks into three 18-bit pieces, so FMA contraction cannot change a result (as in
!< `fossil_predicates`). The sign is exact barring overflow and underflow of the intermediate products.
use penf, only : I4P, I8P, R8P

implicit none
private
public :: implicit_point_t, IPK_EXPLICIT, IPK_LPI, IPK_SSI
public :: ip_explicit, ip_lpi, ip_ssi, ip_same, ip_approx
public :: ip_orient3d, ip_orient2d, ip_compare

integer(I4P), parameter :: IPK_EXPLICIT = 1_I4P !< Input vertex.
integer(I4P), parameter :: IPK_LPI      = 2_I4P !< Line-plane intersection.
integer(I4P), parameter :: IPK_SSI      = 3_I4P !< Coplanar segment-segment intersection.

type :: implicit_point_t
   !< An implicit point: its kind, its defining input points (canonical order) and, for IPK_SSI, the dropped axis.
   integer(I4P) :: kind = 0_I4P            !< IP_*.
   real(R8P)    :: v(3, 5) = 0._R8P        !< Defining points (columns): x; p, q, r, s, t; a, b, c, d.
   integer(I4P) :: drop = 0_I4P            !< Dropped axis of an IPK_SSI.
endtype implicit_point_t

type :: rerr_t
   !< A double with a bound of its absolute error (stage 1).
   real(R8P) :: v = 0._R8P !< Value.
   real(R8P) :: e = 0._R8P !< Error bound.
endtype rerr_t

type :: expn_t
   !< A nonoverlapping expansion, components in increasing magnitude (stage 2).
   real(R8P), allocatable :: c(:)       !< Components.
   integer(I4P)           :: n = 0_I4P  !< Number of components.
endtype expn_t

real(R8P), parameter :: U_ROUND = epsilon(1._R8P) / 2._R8P !< Unit roundoff.

interface operator(+)
   module procedure r_add, x_add
endinterface
interface operator(-)
   module procedure r_sub, x_sub
endinterface
interface operator(*)
   module procedure r_mul, x_mul
endinterface

contains
   ! constructors -------------------------------------------------------------------------------------------------------
   pure function ip_explicit(x) result(p)
   !< An input vertex.
   real(R8P), intent(in)  :: x(3) !< Coordinates.
   type(implicit_point_t) :: p    !< Point.

   p%kind = IPK_EXPLICIT ; p%v(:, 1) = x
   endfunction ip_explicit

   pure function ip_lpi(p0, q0, r0, s0, t0) result(p)
   !< The line through p0, q0 meets the plane through r0, s0, t0 (canonical order of both).
   real(R8P), intent(in)  :: p0(3), q0(3), r0(3), s0(3), t0(3) !< Line, plane.
   type(implicit_point_t) :: p                                  !< Point.
   real(R8P)              :: w(3, 3)

   p%kind = IPK_LPI
   if (lex_less(q0, p0)) then
      p%v(:, 1) = q0 ; p%v(:, 2) = p0
   else
      p%v(:, 1) = p0 ; p%v(:, 2) = q0
   endif
   w(:, 1) = r0 ; w(:, 2) = s0 ; w(:, 3) = t0
   call sort3(w)
   p%v(:, 3:5) = w
   endfunction ip_lpi

   pure function ip_ssi(a0, b0, c0, d0, drop) result(p)
   !< Segment (a0, b0) meets segment (c0, d0) in the coordinate plane dropping axis `drop` (canonical order).
   real(R8P),    intent(in) :: a0(3), b0(3), c0(3), d0(3) !< Segments.
   integer(I4P), intent(in) :: drop                        !< Dropped axis (1..3).
   type(implicit_point_t)   :: p                           !< Point.
   real(R8P)                :: s1(3, 2), s2(3, 2), t(3, 2)

   p%kind = IPK_SSI ; p%drop = drop
   s1(:, 1) = a0 ; s1(:, 2) = b0 ; s2(:, 1) = c0 ; s2(:, 2) = d0
   if (lex_less(s1(:, 2), s1(:, 1))) s1 = s1(:, [2, 1])
   if (lex_less(s2(:, 2), s2(:, 1))) s2 = s2(:, [2, 1])
   if (lex_less(s2(:, 1), s1(:, 1)) .or. (all(s2(:, 1) == s1(:, 1)) .and. lex_less(s2(:, 2), s1(:, 2)))) then
      t = s1 ; s1 = s2 ; s2 = t
   endif
   p%v(:, 1:2) = s1 ; p%v(:, 3:4) = s2
   endfunction ip_ssi

   pure function ip_same(a, b) result(same)
   !< The same construction (identity): same kind, defining points and axis, bit for bit.
   type(implicit_point_t), intent(in) :: a, b !< Points.
   logical                            :: same !< Same identity.

   same = a%kind == b%kind .and. a%drop == b%drop .and. all(a%v == b%v)
   endfunction ip_same

   pure function ip_approx(p) result(x)
   !< The point rounded to doubles, by one formula of its identity (equal identities give equal bits).
   type(implicit_point_t), intent(in) :: p    !< Point.
   real(R8P)                          :: x(3) !< Coordinates.
   type(expn_t)                       :: hx(3), hw
   integer(I4P)                       :: k

   if (p%kind == IPK_EXPLICIT) then
      x = p%v(:, 1) ; return
   endif
   call homog_x(p, hx, hw)
   do k = 1, 3
      x(k) = x_estimate(hx(k)) / x_estimate(hw)
   enddo
   endfunction ip_approx

   ! predicates ---------------------------------------------------------------------------------------------------------
   function ip_orient3d(p, a, b, c) result(sg)
   !< Exact sign of orient3d(a, b, c, p) (Shewchuk's convention: positive when p is below the plane a, b, c seen
   !< counter-clockwise from above), for an implicit p and explicit a, b, c.
   type(implicit_point_t), intent(in) :: p          !< Point.
   real(R8P),              intent(in) :: a(3), b(3), c(3) !< Plane.
   integer(I4P)                       :: sg         !< -1, 0, 1.
   type(rerr_t)                       :: rx(3), rw, rd
   type(expn_t)                       :: hx(3), hw, d

   ! det[a - P; b - P; c - P] W^3 = det[W a - X; W b - X; W c - X]; sign(W^3) = sign(W)
   call homog_r(p, rx, rw)
   rd = det3_r(rw, rx, a, b, c)
   sg = r_sign(rd) * r_sign(rw)
   if (r_sign(rd) /= 2_I4P .and. r_sign(rw) /= 2_I4P) return
   call homog_x(p, hx, hw)
   d = det3_x(hw, hx, a, b, c)
   sg = x_sign(d) * x_sign(hw)
   endfunction ip_orient3d

   function ip_orient2d(p, q, r, drop) result(sg)
   !< Exact sign of the 2D orientation of p, q, r in the coordinate plane dropping axis `drop` (positive:
   !< counter-clockwise in the kept axes (i, j) = (drop + 1, drop + 2) mod 3).
   type(implicit_point_t), intent(in) :: p, q, r !< Points.
   integer(I4P),           intent(in) :: drop    !< Dropped axis.
   integer(I4P)                       :: sg      !< -1, 0, 1.
   type(rerr_t)                       :: px(3), pw, qx(3), qw, rx(3), rw, d
   type(expn_t)                       :: ex(3), ew, fx(3), fw, gx(3), gw, e
   integer(I4P)                       :: i, j, s1, s2, s3, s4

   i = mod(drop, 3_I4P) + 1_I4P ; j = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   call homog_r(p, px, pw) ; call homog_r(q, qx, qw) ; call homog_r(r, rx, rw)
   ! det [[Xp_i, Xp_j, Wp], [Xq_i, Xq_j, Wq], [Xr_i, Xr_j, Wr]] = orient2d Wp Wq Wr
   d = px(i) * (qx(j) * rw - rx(j) * qw) - qx(i) * (px(j) * rw - rx(j) * pw) + rx(i) * (px(j) * qw - qx(j) * pw)
   s1 = r_sign(d) ; s2 = r_sign(pw) ; s3 = r_sign(qw) ; s4 = r_sign(rw)
   if (all([s1, s2, s3, s4] /= 2_I4P)) then
      sg = s1 * s2 * s3 * s4 ; return
   endif
   call homog_x(p, ex, ew) ; call homog_x(q, fx, fw) ; call homog_x(r, gx, gw)
   e = ex(i) * (fx(j) * gw - gx(j) * fw) - fx(i) * (ex(j) * gw - gx(j) * ew) + gx(i) * (ex(j) * fw - fx(j) * ew)
   sg = x_sign(e) * x_sign(ew) * x_sign(fw) * x_sign(gw)
   endfunction ip_orient2d

   function ip_compare(p, q, axis) result(sg)
   !< Exact sign of p_axis - q_axis.
   type(implicit_point_t), intent(in) :: p, q !< Points.
   integer(I4P),           intent(in) :: axis !< Axis (1..3).
   integer(I4P)                       :: sg   !< -1, 0, 1.
   type(rerr_t)                       :: px(3), pw, qx(3), qw, d
   type(expn_t)                       :: ex(3), ew, fx(3), fw, e
   integer(I4P)                       :: s1, s2, s3

   call homog_r(p, px, pw) ; call homog_r(q, qx, qw)
   d = px(axis) * qw - qx(axis) * pw
   s1 = r_sign(d) ; s2 = r_sign(pw) ; s3 = r_sign(qw)
   if (all([s1, s2, s3] /= 2_I4P)) then
      sg = s1 * s2 * s3 ; return
   endif
   call homog_x(p, ex, ew) ; call homog_x(q, fx, fw)
   e = ex(axis) * fw - fx(axis) * ew
   sg = x_sign(e) * x_sign(ew) * x_sign(fw)
   endfunction ip_compare

   ! homogeneous coordinates -----------------------------------------------------------------------------------------
   pure subroutine homog_r(p, x, w)
   !< (X, W) of p with error bounds.
   type(implicit_point_t), intent(in)  :: p    !< Point.
   type(rerr_t),           intent(out) :: x(3) !< X.
   type(rerr_t),           intent(out) :: w    !< W.
   type(rerr_t)                        :: n(3), d(3), e(3), t
   integer(I4P)                        :: k, i, j

   select case (p%kind)
   case (IPK_EXPLICIT)
      do k = 1, 3 ; x(k) = rerr_t(p%v(k, 1), 0._R8P) ; enddo
      w = rerr_t(1._R8P, 0._R8P)
   case (IPK_LPI)
      do k = 1, 3
         d(k) = r_of(p%v(k, 4)) - r_of(p%v(k, 3)) ! s - r
         e(k) = r_of(p%v(k, 5)) - r_of(p%v(k, 3)) ! t - r
      enddo
      n(1) = d(2) * e(3) - d(3) * e(2) ; n(2) = d(3) * e(1) - d(1) * e(3) ; n(3) = d(1) * e(2) - d(2) * e(1)
      w = n(1) * (r_of(p%v(1, 2)) - r_of(p%v(1, 1))) + n(2) * (r_of(p%v(2, 2)) - r_of(p%v(2, 1))) + &
          n(3) * (r_of(p%v(3, 2)) - r_of(p%v(3, 1)))
      t = n(1) * (r_of(p%v(1, 3)) - r_of(p%v(1, 1))) + n(2) * (r_of(p%v(2, 3)) - r_of(p%v(2, 1))) + &
          n(3) * (r_of(p%v(3, 3)) - r_of(p%v(3, 1)))
      do k = 1, 3
         x(k) = w * r_of(p%v(k, 1)) + t * (r_of(p%v(k, 2)) - r_of(p%v(k, 1)))
      enddo
   case (IPK_SSI)
      i = mod(p%drop, 3_I4P) + 1_I4P ; j = mod(p%drop + 1_I4P, 3_I4P) + 1_I4P
      w = (r_of(p%v(i, 2)) - r_of(p%v(i, 1))) * (r_of(p%v(j, 4)) - r_of(p%v(j, 3))) - &
          (r_of(p%v(j, 2)) - r_of(p%v(j, 1))) * (r_of(p%v(i, 4)) - r_of(p%v(i, 3)))
      t = (r_of(p%v(i, 3)) - r_of(p%v(i, 1))) * (r_of(p%v(j, 4)) - r_of(p%v(j, 3))) - &
          (r_of(p%v(j, 3)) - r_of(p%v(j, 1))) * (r_of(p%v(i, 4)) - r_of(p%v(i, 3)))
      do k = 1, 3
         x(k) = w * r_of(p%v(k, 1)) + t * (r_of(p%v(k, 2)) - r_of(p%v(k, 1)))
      enddo
   endselect
   endsubroutine homog_r

   pure subroutine homog_x(p, x, w)
   !< (X, W) of p, exactly.
   type(implicit_point_t), intent(in)  :: p    !< Point.
   type(expn_t),           intent(out) :: x(3) !< X.
   type(expn_t),           intent(out) :: w    !< W.
   type(expn_t)                        :: n(3), d(3), e(3), t
   integer(I4P)                        :: k, i, j

   select case (p%kind)
   case (IPK_EXPLICIT)
      do k = 1, 3 ; x(k) = x_of(p%v(k, 1)) ; enddo
      w = x_of(1._R8P)
   case (IPK_LPI)
      do k = 1, 3
         d(k) = x_of(p%v(k, 4)) - x_of(p%v(k, 3))
         e(k) = x_of(p%v(k, 5)) - x_of(p%v(k, 3))
      enddo
      n(1) = d(2) * e(3) - d(3) * e(2) ; n(2) = d(3) * e(1) - d(1) * e(3) ; n(3) = d(1) * e(2) - d(2) * e(1)
      w = n(1) * (x_of(p%v(1, 2)) - x_of(p%v(1, 1))) + n(2) * (x_of(p%v(2, 2)) - x_of(p%v(2, 1))) + &
          n(3) * (x_of(p%v(3, 2)) - x_of(p%v(3, 1)))
      t = n(1) * (x_of(p%v(1, 3)) - x_of(p%v(1, 1))) + n(2) * (x_of(p%v(2, 3)) - x_of(p%v(2, 1))) + &
          n(3) * (x_of(p%v(3, 3)) - x_of(p%v(3, 1)))
      do k = 1, 3
         x(k) = w * x_of(p%v(k, 1)) + t * (x_of(p%v(k, 2)) - x_of(p%v(k, 1)))
      enddo
   case (IPK_SSI)
      i = mod(p%drop, 3_I4P) + 1_I4P ; j = mod(p%drop + 1_I4P, 3_I4P) + 1_I4P
      w = (x_of(p%v(i, 2)) - x_of(p%v(i, 1))) * (x_of(p%v(j, 4)) - x_of(p%v(j, 3))) - &
          (x_of(p%v(j, 2)) - x_of(p%v(j, 1))) * (x_of(p%v(i, 4)) - x_of(p%v(i, 3)))
      t = (x_of(p%v(i, 3)) - x_of(p%v(i, 1))) * (x_of(p%v(j, 4)) - x_of(p%v(j, 3))) - &
          (x_of(p%v(j, 3)) - x_of(p%v(j, 1))) * (x_of(p%v(i, 4)) - x_of(p%v(i, 3)))
      do k = 1, 3
         x(k) = w * x_of(p%v(k, 1)) + t * (x_of(p%v(k, 2)) - x_of(p%v(k, 1)))
      enddo
   endselect
   endsubroutine homog_x

   pure function det3_r(w, x, a, b, c) result(d)
   !< det[W a - X; W b - X; W c - X], with error bound.
   type(rerr_t), intent(in) :: w, x(3)
   real(R8P),    intent(in) :: a(3), b(3), c(3)
   type(rerr_t)             :: d
   type(rerr_t)             :: ra(3), rb(3), rc(3)
   integer(I4P)             :: k

   do k = 1, 3
      ra(k) = w * r_of(a(k)) - x(k) ; rb(k) = w * r_of(b(k)) - x(k) ; rc(k) = w * r_of(c(k)) - x(k)
   enddo
   d = ra(1) * (rb(2) * rc(3) - rb(3) * rc(2)) - ra(2) * (rb(1) * rc(3) - rb(3) * rc(1)) + &
       ra(3) * (rb(1) * rc(2) - rb(2) * rc(1))
   endfunction det3_r

   pure function det3_x(w, x, a, b, c) result(d)
   !< det[W a - X; W b - X; W c - X], exactly.
   type(expn_t), intent(in) :: w, x(3)
   real(R8P),    intent(in) :: a(3), b(3), c(3)
   type(expn_t)             :: d
   type(expn_t)             :: ra(3), rb(3), rc(3)
   integer(I4P)             :: k

   do k = 1, 3
      ra(k) = w * x_of(a(k)) - x(k) ; rb(k) = w * x_of(b(k)) - x(k) ; rc(k) = w * x_of(c(k)) - x(k)
   enddo
   d = ra(1) * (rb(2) * rc(3) - rb(3) * rc(2)) - ra(2) * (rb(1) * rc(3) - rb(3) * rc(1)) + &
       ra(3) * (rb(1) * rc(2) - rb(2) * rc(1))
   endfunction det3_x

   ! stage 1: doubles with error bounds ---------------------------------------------------------------------------------
   pure elemental function r_of(x) result(r)
   real(R8P), intent(in) :: x
   type(rerr_t)          :: r

   r%v = x ; r%e = 0._R8P
   endfunction r_of

   pure elemental function r_add(a, b) result(r)
   type(rerr_t), intent(in) :: a, b
   type(rerr_t)             :: r

   r%v = a%v + b%v
   r%e = (a%e + b%e + U_ROUND * abs(r%v)) * (1._R8P + 4._R8P * U_ROUND)
   endfunction r_add

   pure elemental function r_sub(a, b) result(r)
   type(rerr_t), intent(in) :: a, b
   type(rerr_t)             :: r

   r%v = a%v - b%v
   r%e = (a%e + b%e + U_ROUND * abs(r%v)) * (1._R8P + 4._R8P * U_ROUND)
   endfunction r_sub

   pure elemental function r_mul(a, b) result(r)
   type(rerr_t), intent(in) :: a, b
   type(rerr_t)             :: r

   r%v = a%v * b%v
   r%e = (abs(a%v) * b%e + abs(b%v) * a%e + a%e * b%e + U_ROUND * abs(r%v)) * (1._R8P + 8._R8P * U_ROUND)
   endfunction r_mul

   pure function r_sign(a) result(s)
   !< Sign of a if its error bound settles it, else 2.
   type(rerr_t), intent(in) :: a
   integer(I4P)             :: s

   if (a%v > a%e) then
      s = 1_I4P
   elseif (-a%v > a%e) then
      s = -1_I4P
   elseif (a%v == 0._R8P .and. a%e == 0._R8P) then
      s = 0_I4P
   else
      s = 2_I4P
   endif
   endfunction r_sign

   ! stage 2: expansions ------------------------------------------------------------------------------------------------
   pure function x_of(x) result(e)
   real(R8P), intent(in) :: x
   type(expn_t)          :: e

   allocate(e%c(1)) ; e%n = 0_I4P
   if (x /= 0._R8P) then
      e%c(1) = x ; e%n = 1_I4P
   endif
   endfunction x_of

   pure function x_add(a, b) result(r)
   type(expn_t), intent(in) :: a, b
   type(expn_t)             :: r
   integer(I4P)             :: k

   allocate(r%c(a%n + b%n + 1)) ; r%n = a%n
   if (a%n > 0) r%c(1:a%n) = a%c(1:a%n)
   do k = 1, b%n
      call grow(r%c, r%n, b%c(k))
   enddo
   endfunction x_add

   pure function x_sub(a, b) result(r)
   type(expn_t), intent(in) :: a, b
   type(expn_t)             :: r
   integer(I4P)             :: k

   allocate(r%c(a%n + b%n + 1)) ; r%n = a%n
   if (a%n > 0) r%c(1:a%n) = a%c(1:a%n)
   do k = 1, b%n
      call grow(r%c, r%n, -b%c(k))
   enddo
   endfunction x_sub

   pure function x_mul(a, b) result(r)
   !< Exact product: every pair of components split into 18-bit pieces, 9 exact piece products each.
   type(expn_t), intent(in) :: a, b
   type(expn_t)             :: r
   real(R8P)                :: pa(3), pb(3)
   integer(I4P)             :: k, l, m, o

   allocate(r%c(9 * a%n * b%n + 1)) ; r%n = 0_I4P
   do k = 1, a%n
      call split3(a%c(k), pa)
      do l = 1, b%n
         call split3(b%c(l), pb)
         do m = 1, 3
            if (pa(m) == 0._R8P) cycle
            do o = 1, 3
               if (pb(o) == 0._R8P) cycle
               call grow(r%c, r%n, pa(m) * pb(o))
            enddo
         enddo
      enddo
   enddo
   endfunction x_mul

   pure function x_sign(a) result(s)
   !< Exact sign: the largest (last) nonzero component.
   type(expn_t), intent(in) :: a
   integer(I4P)             :: s
   integer(I4P)             :: k

   s = 0_I4P
   do k = a%n, 1, -1
      if (a%c(k) > 0._R8P) then
         s = 1_I4P ; return
      elseif (a%c(k) < 0._R8P) then
         s = -1_I4P ; return
      endif
   enddo
   endfunction x_sign

   pure function x_estimate(a) result(v)
   !< The expansion's value rounded (sum from the smallest component, deterministic).
   type(expn_t), intent(in) :: a
   real(R8P)                :: v
   integer(I4P)             :: k

   v = 0._R8P
   do k = 1, a%n
      v = v + a%c(k)
   enddo
   endfunction x_estimate

   pure subroutine grow(e, n, b)
   !< e = e + b exactly (Shewchuk's Grow-Expansion with zero elimination); e has room for n + 1 components.
   real(R8P),    intent(inout) :: e(:)
   integer(I4P), intent(inout) :: n
   real(R8P),    intent(in)    :: b
   real(R8P)                   :: q, s, h, bv, av
   integer(I4P)                :: i, m

   q = b ; m = 0_I4P
   do i = 1, n
      s = q + e(i) ; bv = s - q ; av = s - bv ; h = (q - av) + (e(i) - bv)
      q = s
      if (h /= 0._R8P) then
         m = m + 1_I4P ; e(m) = h
      endif
   enddo
   if (q /= 0._R8P) then
      m = m + 1_I4P ; e(m) = q
   endif
   n = m
   endsubroutine grow

   pure subroutine split3(x, piece)
   !< x = piece(1) + piece(2) + piece(3) exactly, each of at most 18 significant bits (bit masks, no multiplication).
   real(R8P), intent(in)  :: x
   real(R8P), intent(out) :: piece(3)
   real(R8P)              :: r

   piece(1) = truncate18(x) ; r = x - piece(1)
   piece(2) = truncate18(r) ; piece(3) = r - piece(2)
   endsubroutine split3

   pure function truncate18(x) result(t)
   real(R8P), intent(in) :: x
   real(R8P)             :: t
   integer(I8P), parameter :: LOW35 = 34359738367_I8P
   integer(I8P)          :: bits

   bits = transfer(x, bits)
   if (iand(ishft(bits, -52), 2047_I8P) == 0_I8P) then
      t = x ; return
   endif
   bits = iand(bits, not(LOW35))
   t = transfer(bits, t)
   endfunction truncate18

   ! ordering helpers ---------------------------------------------------------------------------------------------------
   pure function lex_less(a, b) result(less)
   real(R8P), intent(in) :: a(3), b(3)
   logical               :: less

   if (a(1) /= b(1)) then
      less = a(1) < b(1)
   elseif (a(2) /= b(2)) then
      less = a(2) < b(2)
   else
      less = a(3) < b(3)
   endif
   endfunction lex_less

   pure subroutine sort3(w)
   real(R8P), intent(inout) :: w(3, 3)
   real(R8P)                :: t(3)

   if (lex_less(w(:, 2), w(:, 1))) then ; t = w(:, 1) ; w(:, 1) = w(:, 2) ; w(:, 2) = t ; endif
   if (lex_less(w(:, 3), w(:, 2))) then ; t = w(:, 2) ; w(:, 2) = w(:, 3) ; w(:, 3) = t ; endif
   if (lex_less(w(:, 2), w(:, 1))) then ; t = w(:, 1) ; w(:, 1) = w(:, 2) ; w(:, 2) = t ; endif
   endsubroutine sort3
endmodule fossil_implicit_points
