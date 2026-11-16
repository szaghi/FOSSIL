!< Figure of the reference page "Finite-difference wall closure", What accuracy to expect: the 2-D half of the validation
!< fossil_val_fd_closure. A circle of radius R = 0.5 on uniform point lines of spacing h = R/8 .. R/128, reach 2; at each
!< ghost of fd_closure_2d the image-line closure of the smooth field u = cos(1.1 x - 0.7 y + 0.3) exp(0.5 x) + sin(0.9 y)
!< with three probes (l = sqrt(2) h), u at the probes by tensor Lagrange on 5 x 5 grid values from the probe's cell, the
!< wall datum of the true circle at the foot (Dirichlet u, Neumann d_n u), the cubic along the normal evaluated at the
!< ghost's signed distance. The largest error over the ghosts in three regimes: exact geometry (foot and normal of the
!< circle), facets refined with the grid (N = 48 .. 768 segments, h_s / h = 1.05), facets frozen (N = 48). The figure plots
!< the Dirichlet errors against h, with the floor estimate of the frozen facets.
program fig_ref_fd_closure_accuracy
use fossil,     only : curve_object, ib_closure_t, fd_closure_2d
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter :: PI = acos(-1._R8P), R = 0.5_R8P
integer(I4P), parameter :: NL = 5, M = 3
type(curve_object)      :: wall, frozen
type(ib_closure_t)      :: t
real(R8P), allocatable  :: xs(:), ys(:)
real(R8P)               :: h(NL), ed(NL, 3), en(NL, 3), lam, gmax, sag, est
type(vector_R8P)        :: pts(NL)
character(16)           :: text
integer(I4P)            :: k, reg, u, j

ed = 0._R8P ; en = 0._R8P ; lam = 0._R8P ; gmax = 0._R8P
call circle(48_I4P, frozen)
do k = 1, NL
   h(k) = R / (8._R8P * 2._R8P**(k - 1))
   xs = [(-1.1_R8P + (j - 1) * h(k) + 0.0371_R8P * h(k), j = 1, nint(2.2_R8P / h(k)) + 1)]
   ys = [(-1.1_R8P + (j - 1) * h(k) + 0.0613_R8P * h(k), j = 1, nint(2.2_R8P / h(k)) + 1)]
   call circle(48_I4P * 2_I4P**(k - 1), wall)
   call fd_closure_2d(wall, xs, ys, [2_I4P, 2_I4P], t, lengths=[(j * sqrt(2._R8P) * h(k), j = 1, M)])
   call closure_errors(1_I4P)                                     ! exact geometry, on the same ghosts
   call closure_errors(2_I4P)                                     ! the facets of the grid's size
   call fd_closure_2d(frozen, xs, ys, [2_I4P, 2_I4P], t, lengths=[(j * sqrt(2._R8P) * h(k), j = 1, M)])
   call closure_errors(3_I4P)                                     ! the facets frozen
enddo
sag = R * (1._R8P - cos(PI / 48._R8P)) ; est = lam * gmax * sag
print '(A)', 'circle R = 0.5, reach 2, image line of 3 probes (l = sqrt(2) h), max |u_ghost - u(ghost)|'
print '(A)', '                       Dirichlet                           Neumann'
print '(A)', '          h      exact  h_s = c h     frozen        exact  h_s = c h     frozen'
do k = 1, NL
   print '(2X,ES9.2,3ES11.2,2X,3ES11.2)', h(k), ed(k, :), en(k, :)
enddo
print '(A11,3F11.2,2X,3F11.2)', 'order', (log(ed(NL - 1, reg) / ed(NL, reg)) / log(2._R8P), reg = 1, 3), &
                                             (log(en(NL - 1, reg) / en(NL, reg)) / log(2._R8P), reg = 1, 3)
print '(A,ES9.2,A,ES9.2,A,F5.2,A,F5.2,A,ES9.2,A)', 'frozen facets: Dirichlet floor ', ed(NL, 3), &
      ', estimate Lambda |d_n u| s = ', est, ' (Lambda ', lam, ', |d_n u| ', gmax, ', sag s ', sag, ')'

call fig_open(u, 'ref-fd-closure-accuracy', 'Dirichlet ghost error against h, three probes', -90._R8P, 90._R8P)
call fig_segment(u, 'box', p(h(1), 1.e-10_R8P), p(h(NL), 1.e-10_R8P))
call fig_segment(u, 'box', p(h(1), 1.e-10_R8P), p(h(1), 1.e-1_R8P))
do j = -9, -1, 2                                                  ! y ticks: 1e-9 .. 1e-1
   write(text, '(A,I0)') '1e', j
   call fig_label(u, 'box', vector_R8P(-1.1_R8P, real(j, R8P), 0._R8P), trim(text))
   call fig_segment(u, 'box', vector_R8P(-0.2_R8P, real(j, R8P), 0._R8P), vector_R8P(0._R8P, real(j, R8P), 0._R8P))
enddo
do k = 1, NL                                                      ! x ticks: h = R/8 .. R/128
   write(text, '(A,I0)') 'R/', 8 * 2**(k - 1)
   call fig_label(u, 'box', p(h(k), 10._R8P**(-10.6_R8P)) - vector_R8P(0.3_R8P, 0._R8P, 0._R8P), trim(text))
enddo
call fig_label(u, 'box', p(h(NL), 10._R8P**(-10.6_R8P)) + vector_R8P(0.9_R8P, 0._R8P, 0._R8P), 'h')
call fig_segment(u, 'box', p(h(1), est), p(h(NL), est))
call fig_label(u, 'box', p(h(2), est) + vector_R8P(0._R8P, 0.35_R8P, 0._R8P), 'floor estimate Lambda |d_n u| s')
do reg = 1, 3
   do k = 1, NL
      pts(k) = p(h(k), ed(k, reg))
   enddo
   select case (reg)
   case (1)
      call fig_segments(u, 'probe', pts(1:NL - 1), pts(2:NL)) ; call fig_points(u, 'probe', pts)
      call fig_label(u, 'probe', pts(NL - 1) + vector_R8P(0.6_R8P, 0.3_R8P, 0._R8P), 'exact geometry, order 4')
   case (2)
      call fig_segments(u, 'accent', pts(1:NL - 1), pts(2:NL)) ; call fig_points(u, 'accent', pts)
      call fig_label(u, 'accent', pts(NL - 1) + vector_R8P(0.6_R8P, 0.3_R8P, 0._R8P), 'h_s = c h, order 2')
   case (3)
      call fig_segments(u, 'outside', pts(1:NL - 1), pts(2:NL)) ; call fig_points(u, 'outside', pts)
      call fig_label(u, 'outside', pts(NL - 1) + vector_R8P(0._R8P, -0.55_R8P, 0._R8P), 'h_s frozen: the floor')
   endselect
enddo
call fig_close(u)
contains
   subroutine closure_errors(regime)
   !< The largest Dirichlet and Neumann ghost errors of the table t in a regime (1 exact, 2 and 3 the table's geometry).
   integer(I4P), intent(in) :: regime !< Regime.
   real(R8P)                :: foot(2), nrm(2), eta, probe(2, M), ui(M), wall, dwall, pw(2), q, l0
   integer(I4P)             :: g, kk, cell(2, M)

   do g = 1, t%n
      if (regime == 1_I4P) then
         nrm = t%points(:, g) / norm2(t%points(:, g)) ; foot = R * nrm ; eta = norm2(t%points(:, g)) - R
         do kk = 1, M
            probe(:, kk) = foot + t%lengths(kk) * nrm
            cell(:, kk) = [locate(xs, probe(1, kk)), locate(ys, probe(2, kk))]
         enddo
      else
         foot = t%foot(:, g) ; nrm = t%normal(:, g) ; eta = t%distance(g)
         probe = t%probe(:, :, g) ; cell = t%probe_cell(:, :, g)
      endif
      pw = R * foot / norm2(foot)                                 ! the true wall point
      do kk = 1, M                                                ! u at the probes: the solver's interpolation
         ui(kk) = interp(probe(:, kk), cell(:, kk))
      enddo
      wall = field(pw)                                            ! Dirichlet datum of the wall at the foot
      call lagrange([0._R8P, t%lengths], [wall, ui], eta, q, l0)  ! along the normal, at the ghost's distance
      ed(k, regime) = max(ed(k, regime), abs(q - field(t%points(:, g))))
      dwall = dot_product(grad(pw), pw / R)                       ! Neumann datum: d_n u along the closure's normal
      q = neumann(t%lengths, ui, dwall, eta)
      en(k, regime) = max(en(k, regime), abs(q - field(t%points(:, g))))
      if (regime == 3_I4P .and. k == NL) then
         lam = max(lam, abs(l0)) ; gmax = max(gmax, abs(dwall))
      endif
   enddo
   endsubroutine closure_errors

   subroutine lagrange(tn, fv, s, q, l0)
   !< The interpolating polynomial through (tn(j), fv(j)) at s, and the weight of the first datum.
   real(R8P), intent(in)  :: tn(:), fv(:), s !< Nodes, values, abscissa.
   real(R8P), intent(out) :: q, l0           !< Value, weight of fv(1).
   real(R8P)              :: w
   integer(I4P)           :: i, jj

   q = 0._R8P
   do jj = 1, size(tn)
      w = 1._R8P
      do i = 1, size(tn)
         if (i /= jj) w = w * (s - tn(i)) / (tn(jj) - tn(i))
      enddo
      q = q + w * fv(jj)
      if (jj == 1) l0 = w
   enddo
   endsubroutine lagrange

   function neumann(tn, fv, gn, s) result(q)
   !< The cubic q with q'(0) = gn and q(tn(i)) = fv(i), i = 1..3, at s (Cramer-free: q = q0 + gn t + a t^2 + b t^3).
   real(R8P), intent(in) :: tn(M), fv(M), gn, s !< Probe lengths, values, normal derivative, abscissa.
   real(R8P)             :: q, a(M, M), b(M), c(M), f
   integer(I4P)          :: i, jj

   do i = 1, M                                                    ! unknowns q0, a, b; rows q(tn(i)) - gn tn(i)
      a(i, :) = [1._R8P, tn(i)**2, tn(i)**3] ; b(i) = fv(i) - gn * tn(i)
   enddo
   do i = 1, M - 1                                                ! Gaussian elimination (the nodes are distinct)
      do jj = i + 1, M
         f = a(jj, i) / a(i, i) ; a(jj, :) = a(jj, :) - f * a(i, :) ; b(jj) = b(jj) - f * b(i)
      enddo
   enddo
   do i = M, 1, -1
      c(i) = (b(i) - dot_product(a(i, i + 1:), c(i + 1:))) / a(i, i)
   enddo
   q = c(1) + gn * s + c(2) * s**2 + c(3) * s**3
   endfunction neumann

   function interp(q, cell) result(val)
   !< Tensor Lagrange interpolation of the field's grid values on 5 x 5 points, centred on the nearer end of the cell.
   real(R8P),    intent(in) :: q(2)    !< Point.
   integer(I4P), intent(in) :: cell(2) !< Enclosing cell (low-side rule).
   real(R8P)                :: val, wx(5), wy(5)
   integer(I4P)             :: sx, sy, i, jj

   call weights(xs, q(1), cell(1), sx, wx) ; call weights(ys, q(2), cell(2), sy, wy)
   val = 0._R8P
   do jj = 1, 5 ; do i = 1, 5
      val = val + wx(i) * wy(jj) * field([xs(sx + i - 1), ys(sy + jj - 1)])
   enddo ; enddo
   endfunction interp

   subroutine weights(lines, q, i, s, w)
   !< The 5 Lagrange weights at q from the cell i, centred on the nearer end of the cell.
   real(R8P),    intent(in)  :: lines(:), q !< Lines, coordinate.
   integer(I4P), intent(in)  :: i           !< Cell.
   integer(I4P), intent(out) :: s           !< First point.
   real(R8P),    intent(out) :: w(5)        !< Weights.
   integer(I4P)              :: jj, kk

   if (i < 1) error stop 'probe outside the lines'
   s = i - 2 ; if (q - lines(i) > lines(i + 1) - q) s = s + 1
   if (s < 1 .or. s + 4 > size(lines)) error stop 'stencil outside the lines'
   do jj = 1, 5
      w(jj) = 1._R8P
      do kk = 1, 5
         if (kk /= jj) w(jj) = w(jj) * (q - lines(s + kk - 1)) / (lines(s + jj - 1) - lines(s + kk - 1))
      enddo
   enddo
   endsubroutine weights

   pure function locate(lines, q) result(i)
   !< The enclosing cell by the low-side rule, lines(i) < q <= lines(i+1), 0 outside.
   real(R8P), intent(in) :: lines(:), q !< Lines, coordinate.
   integer(I4P)          :: i

   i = 0
   do while (i < size(lines) - 1)
      if (lines(i + 1) < q .and. q <= lines(i + 2)) then
         i = i + 1 ; return
      endif
      i = i + 1
   enddo
   i = 0
   endfunction locate

   pure function field(x) result(f)
   !< The smooth field.
   real(R8P), intent(in) :: x(2) !< Point.
   real(R8P)             :: f

   f = cos(1.1_R8P * x(1) - 0.7_R8P * x(2) + 0.3_R8P) * exp(0.5_R8P * x(1)) + sin(0.9_R8P * x(2))
   endfunction field

   pure function grad(x) result(g)
   !< Its gradient.
   real(R8P), intent(in) :: x(2) !< Point.
   real(R8P)             :: g(2), a, e

   a = 1.1_R8P * x(1) - 0.7_R8P * x(2) + 0.3_R8P ; e = exp(0.5_R8P * x(1))
   g = [(-1.1_R8P * sin(a) + 0.5_R8P * cos(a)) * e, 0.7_R8P * sin(a) * e + 0.9_R8P * cos(0.9_R8P * x(2))]
   endfunction grad

   function p(hk, e) result(v)
   !< A plot point: x = 8 log10(h(1) / hk) (h decreasing to the right), y = log10(e).
   real(R8P), intent(in) :: hk, e !< Spacing, error.
   type(vector_R8P)      :: v

   v = vector_R8P(8._R8P * log10(h(1) / hk), log10(e), 0._R8P)
   endfunction p

   subroutine circle(n, c)
   !< The circle of radius R as n segments, vertices on it (turned by 0.1 rad).
   integer(I4P),       intent(in)    :: n !< Segments.
   type(curve_object), intent(inout) :: c !< Curve.
   integer(I4P)                      :: i

   call c%adopt_loops(u=[(R * cos(0.1_R8P + 2 * PI * i / n), i = 0, n - 1)], &
                      v=[(R * sin(0.1_R8P + 2 * PI * i / n), i = 0, n - 1)], loop_sizes=[n])
   endsubroutine circle
endprogram fig_ref_fd_closure_accuracy
