!< FOSSIL, 2-D hyperbolic grid extrusion from a curve (issue #74).
module fossil_hyperbolic
!< FOSSIL, 2-D hyperbolic grid extrusion from a curve (issue #74).
!<
!< Design: docs/guide/advanced/body-fitted.md. `hyperbolic_grid_2d` marches layers from a curve (closed: an O-grid;
!< open: a C- or H-type sheet), to the right of the curve's direction (outward for a counter-clockwise closed curve).
!< The conditions of the hyperbolic system (gk: Steger & Chaussee 1980), orthogonality r_xi . r_eta = 0 and cell size
!< |r_xi x r_eta| = dA, give at each point the step r_eta = h n, normal to the current layer (central differences;
!< one-sided at the ends of an open curve), with h the layer's height. The step directions are smoothed by the implicit
!< dissipation (I - eps_k delta_xixi) D = h n (a cyclic tridiagonal solve for a closed curve, Neumann ends for an
!< open one), eps_k = epsilon (k - 1) / (layers - 1) growing with the layer. Each step has the mitre length
!< h / cos(beta / 2), beta the layer's turning angle at the point (capped at 4 h): the exact offset of a polygon, so a
!< layer stays at height h from both sides of a corner (a plain h along the bisector makes the corner lag and its cells
!< fold); from a regular polygon the layers are its exact offsets, vertices on the rays. The cells are checked with
!< WP2's quality (scaled certified min J): `invalid` counts the cells not certified valid (folds).
use penf,                     only : I4P, R8P
use fossil_curvilinear_block, only : curvilinear_block_2d_t

implicit none
private
public :: hyperbolic_grid_2d, geometric_spacing

contains
   pure function geometric_spacing(h1, ratio, layers) result(h)
   !< Layer heights h1, h1 ratio, h1 ratio^2, ...
   real(R8P),    intent(in) :: h1, ratio !< First height, growth ratio.
   integer(I4P), intent(in) :: layers    !< Layers.
   real(R8P)                :: h(layers) !< Heights.
   integer(I4P)             :: k

   do k = 1, layers
      h(k) = h1 * ratio**(k - 1)
   enddo
   endfunction geometric_spacing

   subroutine hyperbolic_grid_2d(u, v, spacing, block, closed, epsilon, invalid, status)
   !< March len(spacing) layers from the curve (u, v) [n] (closed: the last point is not repeated; the block repeats the
   !< first point at i = n + 1). status: 0, 1 (bad input), 2 (cells not certified valid: a fold).
   real(R8P),                    intent(in)            :: u(:), v(:) !< Curve points [n].
   real(R8P),                    intent(in)            :: spacing(:) !< Layer heights [layers] (> 0).
   type(curvilinear_block_2d_t), intent(out)           :: block      !< Block [n (+1), layers + 1].
   logical,                      intent(in),  optional :: closed     !< Closed curve (default .true.).
   real(R8P),                    intent(in),  optional :: epsilon    !< Largest smoothing (default 0.5).
   integer(I4P),                 intent(out), optional :: invalid    !< Cells not certified valid.
   integer(I4P),                 intent(out), optional :: status     !< 0, 1, 2.
   real(R8P), allocatable                              :: x(:,:), y(:,:), dx(:), dy(:), mj(:,:)
   real(R8P)                                           :: eps0, ek, t(2), l, nl
   integer(I4P)                                        :: n, nl_, k, i, ni, bad
   logical                                             :: cl

   cl = .true. ; if (present(closed)) cl = closed
   eps0 = 0.5_R8P ; if (present(epsilon)) eps0 = epsilon
   n = size(u, kind=I4P) ; nl_ = size(spacing, kind=I4P)
   if (present(invalid)) invalid = 0
   if (n < 3 .or. size(v) /= n .or. nl_ < 1 .or. any(.not. (spacing > 0._R8P)) .or. eps0 < 0._R8P) then
      if (present(status)) status = 1
      return
   endif
   allocate(x(n, nl_ + 1), y(n, nl_ + 1), dx(n), dy(n))
   x(:, 1) = u ; y(:, 1) = v
   do k = 1, nl_
      do i = 1, n ! normals of the current layer (right of the direction)
         t = tangent(x(:, k), y(:, k), i, cl)
         l = norm2(t)
         if (l > 0._R8P) then
            dx(i) = t(2) / l ; dy(i) = -t(1) / l
         else
            dx(i) = 0._R8P ; dy(i) = 0._R8P
         endif
      enddo
      ek = 0._R8P
      if (nl_ > 1) ek = eps0 * real(k - 1, R8P) / real(nl_ - 1, R8P)
      if (ek > 0._R8P) then
         call smooth(dx, ek, cl) ; call smooth(dy, ek, cl)
      endif
      do i = 1, n
         nl = hypot(dx(i), dy(i)) / mitre(x(:, k), y(:, k), i, cl)
         if (nl > 0._R8P) then
            x(i, k + 1) = x(i, k) + spacing(k) * dx(i) / nl ; y(i, k + 1) = y(i, k) + spacing(k) * dy(i) / nl
         else
            x(i, k + 1) = x(i, k) ; y(i, k + 1) = y(i, k)
         endif
      enddo
   enddo
   ni = merge(n + 1_I4P, n, cl)
   if (cl) then
      call block%initialize(reshape([(x(:, k), x(1, k), k = 1, nl_ + 1)], [ni, nl_ + 1]), &
                            reshape([(y(:, k), y(1, k), k = 1, nl_ + 1)], [ni, nl_ + 1]))
   else
      call block%initialize(x, y)
   endif
   call block%quality(min_jacobian=mj)
   bad = count(.not. (mj > 0._R8P))
   if (present(invalid)) invalid = bad
   if (present(status)) status = merge(2_I4P, 0_I4P, bad > 0)
   endsubroutine hyperbolic_grid_2d

   pure function mitre(x, y, i, cl) result(f)
   !< 1 / cos(beta / 2), beta the turning angle at point i (1 at the ends of an open curve), at most 4.
   real(R8P),    intent(in) :: x(:), y(:) !< Layer.
   integer(I4P), intent(in) :: i          !< Point.
   logical,      intent(in) :: cl         !< Closed.
   real(R8P)                :: f          !< Factor.
   real(R8P)                :: a(2), b(2), c
   integer(I4P)             :: n, ia, ib

   n = size(x, kind=I4P) ; f = 1._R8P
   if (.not. cl .and. (i == 1 .or. i == n)) return
   ia = merge(n, i - 1_I4P, i == 1) ; ib = merge(1_I4P, i + 1_I4P, i == n)
   a = [x(i) - x(ia), y(i) - y(ia)] ; b = [x(ib) - x(i), y(ib) - y(i)]
   if (norm2(a) == 0._R8P .or. norm2(b) == 0._R8P) return
   c = dot_product(a, b) / (norm2(a) * norm2(b))
   f = min(4._R8P, sqrt(2._R8P / max(1._R8P / 8._R8P, 1._R8P + c)))
   endfunction mitre

   pure function tangent(x, y, i, cl) result(t)
   !< Central difference x_{i+1} - x_{i-1} (periodic if cl; one-sided at the ends otherwise).
   real(R8P),    intent(in) :: x(:), y(:) !< Layer.
   integer(I4P), intent(in) :: i          !< Point.
   logical,      intent(in) :: cl         !< Closed.
   real(R8P)                :: t(2)       !< Tangent.
   integer(I4P)             :: n, a, b

   n = size(x, kind=I4P)
   if (cl) then
      a = merge(n, i - 1_I4P, i == 1) ; b = merge(1_I4P, i + 1_I4P, i == n)
   else
      a = max(1_I4P, i - 1_I4P) ; b = min(n, i + 1_I4P)
   endif
   t = [x(b) - x(a), y(b) - y(a)]
   endfunction tangent

   subroutine smooth(d, e, cl)
   !< Solve (1 + 2e) D_i - e (D_{i-1} + D_{i+1}) = d_i in place: cyclic (Sherman-Morrison) if cl, else Neumann ends
   !< (D_0 = D_1, D_{n+1} = D_n).
   real(R8P), intent(inout) :: d(:) !< Right-hand side, then the solution.
   real(R8P), intent(in)    :: e    !< Smoothing.
   logical,   intent(in)    :: cl   !< Closed.
   real(R8P), allocatable   :: a(:), b(:), c(:), z(:), w(:)
   real(R8P)                :: gam, fac
   integer(I4P)             :: n

   n = size(d, kind=I4P)
   allocate(a(n), b(n), c(n)) ; a = -e ; c = -e ; b = 1._R8P + 2._R8P * e
   if (.not. cl) then
      b(1) = 1._R8P + e ; b(n) = 1._R8P + e
      call thomas(a, b, c, d)
      return
   endif
   ! cyclic: A = T + u v^T with corners a(1) = c(n) = -e; u = [gam, 0.., c(n)], v = [1, 0.., a(1) / gam]
   gam = -b(1)
   b(1) = b(1) - gam ; b(n) = b(n) - c(n) * a(1) / gam
   allocate(z(n), w(n)) ; z = d ; w = 0._R8P ; w(1) = gam ; w(n) = c(n)
   call thomas(a, b, c, z) ; call thomas(a, b, c, w)
   fac = (z(1) + a(1) * z(n) / gam) / (1._R8P + w(1) + a(1) * w(n) / gam)
   d = z - fac * w
   endsubroutine smooth

   subroutine thomas(a, b, c, d)
   !< Tridiagonal solve (sub a, diagonal b, super c; a(1), c(n) ignored), d overwritten with the solution.
   real(R8P), intent(in)    :: a(:), b(:), c(:) !< Bands.
   real(R8P), intent(inout) :: d(:)             !< Right-hand side / solution.
   real(R8P), allocatable   :: cp(:)
   real(R8P)                :: m
   integer(I4P)             :: i, n

   n = size(d, kind=I4P) ; allocate(cp(n))
   cp(1) = c(1) / b(1) ; d(1) = d(1) / b(1)
   do i = 2, n
      m = b(i) - a(i) * cp(i - 1)
      cp(i) = c(i) / m ; d(i) = (d(i) - a(i) * d(i - 1)) / m
   enddo
   do i = n - 1, 1, -1
      d(i) = d(i) - cp(i) * d(i + 1)
   enddo
   endsubroutine thomas
endmodule fossil_hyperbolic
