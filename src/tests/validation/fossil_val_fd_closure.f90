!< FOSSIL validation: what accuracy the closure table of a finite-difference wall closure supports with planar facets
!< (#109 WP1 step 05). A geometry-closure study, no PDE solve.
program fossil_val_fd_closure
!< FOSSIL validation: what accuracy the closure table of a finite-difference wall closure supports with planar facets
!< (#109 WP1 step 05). A geometry-closure study, no PDE solve.
!<
!< Bodies: a circle of radius R = 1/2 (a curve of N segments, vertices on the circle) and a sphere of radius 1/2 (an
!< icosphere, vertices on the sphere), centred at the origin, on uniform point lines of spacing h (h halved each level:
!< 2-D 8 levels R/8 .. R/1024, 3-D 5 levels R/4 .. R/64), shifted off the centre, stencil half-width 2 per axis.
!<
!< For every ghost of `fd_closure_2d` / `fd_closure_3d` the ghost value of the smooth field
!<   2-D  u = cos(1.1 x - 0.7 y + 0.3) exp(0.5 x) + sin(0.9 y)
!<   3-D  u = cos(1.1 x - 0.7 y + 0.3) exp(0.6 z) + sin(0.9 y + 1.2 z)
!< is built by an image-line closure: the m probes of lengths l, 2 l, .., m l (l = sqrt(2) h, the table's `lengths`), u at
!< each probe interpolated from the grid values by tensor Lagrange on 4 points per axis (m = 1, 2) or 5 (m = 3) placed
!< from the probe's enclosing cell (`probe_cell`), the wall datum at the foot, the polynomial of degree m along the
!< normal through the m + 1 data, evaluated at the ghost's signed distance. Error = |u_ghost - u(ghost)|, the maximum
!< over the ghosts. Dirichlet: the datum is u at the wall; Neumann: the normal derivative there.
!<
!< The wall data are those of the true wall: at the faceted foot x_f the solver imposes g_D = u(pi(x_f)) (Dirichlet) or
!< g_N = grad u(pi(x_f)) . n_true (Neumann, along the closure's normal), pi(x_f) = R x_f / |x_f| the true wall point.
!< The difference between the faceted and the true wall (the sag, O(h_s^2 / R), and the normal tilt, O(h_s / R)) is the
!< geometry error of planar facets. Three regimes:
!<   exact geometry  the ghosts of the h_s = c h table, foot R p / |p|, normal p / |p|, distance |p| - R, probes on that
!<                   line (cells by the same low-side rule): the closure's own error, design order m + 1;
!<   h_s = c h       facets refined with the grid (2-D N = 48 2^k, 3-D icosphere of 2 + k subdivisions; h_s / h about
!<                   1.05 in 2-D, 1.3 in 3-D): the table's foot, normal, distance, probes and cells; geometry error O(h_s^2) = O(h^2);
!<   h_s frozen      the coarsest mesh of the h_s = c h regime at every level: the error stops at the geometry floor.
!<
!< Printed checks: the closure's curvature at the feet (h_s = c h: k1, k2 against 1/R) and, for the body at a pose
!< spinning about the pose's origin, the wall velocity against the rigid field V + omega x (x - c) at the foot (to
!< round-off) and at the true wall point (the difference is at most |omega| times the sag).
!<
!< Pass criteria (`Are all tests passed? T` iff all hold):
!<   1. exact geometry: the observed order of the two finest levels within 0.3 of m + 1, Dirichlet and Neumann, m = 1..3;
!<   2. h_s = c h: the observed Dirichlet order of the two finest levels at least 2 - 0.3 (the planar-facet bound O(h^2)
!<      with the tolerance of criterion 1), m = 1..3;
!<   3. h_s frozen: the Dirichlet error stops decreasing (observed order of the two finest levels below 0.5) and its
!<      floor is at most the estimate Lambda_m |d_n u|_max s, s the largest sag of the mesh (O(h_s^2 / R)), Lambda_m the
!<      largest weight of the wall datum in the extrapolation to a ghost;
!<   4. curvature at the feet: 2-D 1/R to round-off, 100 eps (R / h_s)^2 (exact on a uniform polygon of a circle, the
!<      round-off of a three-point circle); 3-D within 1 % at the finest level and decreasing with the level;
!<   5. wall velocity: equal to V + omega x (x_f - c) within 1e-14, and within |omega| s (1 + 1e-6) + 1e-14 of the rigid
!<      field at the true wall point;
!<   6. every call STATUS_OK and every probe inside the block with its interpolation stencil.
use fossil,              only : surface_stl_object, curve_object, ib_closure_t, fd_closure_3d, fd_closure_2d, &
                                rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, pose_2d, STATUS_OK
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P), R = 0.5_R8P
integer(I4P), parameter :: NL = 8_I4P, NLV(2) = [8_I4P, 5_I4P], NM = 3_I4P, LEVEL0_3D = 2_I4P, N0_2D = 48_I4P
integer(I4P), parameter :: EXACT = 1_I4P, CH = 2_I4P, FROZEN = 3_I4P, DIR = 1_I4P, NEU = 2_I4P
character(14), parameter :: REG_NAME(3) = [character(14) :: 'exact geometry', 'h_s = c h', 'h_s frozen']
real(R8P)               :: err(NM, 2, 3, NL, 2) = 0._R8P !< Max ghost error (m, BC, regime, level, 2-D/3-D).
real(R8P)               :: lam(NM, 3, NL, 2) = 0._R8P    !< Largest Dirichlet weight of the wall datum.
real(R8P)               :: gmx(3, NL, 2) = 0._R8P        !< Largest |d_n u| at the true wall points.
integer(I4P)            :: ngh(3, NL, 2) = 0_I4P         !< Ghosts.
real(R8P)               :: hh(NL, 2), hsc(NL, 2), sagc(NL, 2), hsf(2), sagf(2)
real(R8P)               :: curv(NL, 2) = 0._R8P, vel1(NL, 2) = 0._R8P, vel2(NL, 2) = 0._R8P, vtol(NL, 2) = 0._R8P
real(R8P), allocatable  :: xs(:), ys(:), zs(:), ug(:,:,:)
integer(I4P)            :: nbad = 0_I4P, nstat = 0_I4P
logical                 :: ok = .true.

call study_2d
call study_3d
call report(1_I4P, '2-D circle R = 0.5 (curve of N segments)')
call report(2_I4P, '3-D sphere R = 0.5 (icosphere)')
call check('every call STATUS_OK, every probe and stencil inside the block', nstat == 0_I4P .and. nbad == 0_I4P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine study_2d()
   !< The 2-D study: 5 levels, three regimes, curvature and velocity checks.
   type(curve_object)    :: cv, fr
   type(ib_closure_t)    :: t
   type(rigid_pose_2d_t) :: pose
   real(R8P)             :: h, c(2), rv(2), vf(2), vt(2), w
   integer(I4P)          :: k, nv, st, g

   call circle(N0_2D, fr) ; hsf(1) = 2._R8P * R * sin(PI / N0_2D) ; sagf(1) = R * (1._R8P - cos(PI / N0_2D))
   do k = 1_I4P, NLV(1)
      h = R / (8._R8P * 2._R8P**(k - 1)) ; hh(k, 1) = h
      call make_lines(2_I4P, h)
      nv = N0_2D * 2_I4P**(k - 1)
      call circle(nv, cv) ; call cv%build_curvature(status=st) ; call stat(st)
      hsc(k, 1) = 2._R8P * R * sin(PI / nv) ; sagc(k, 1) = R * (1._R8P - cos(PI / nv))
      call fd_closure_2d(cv, xs, ys, [2_I4P, 2_I4P], t, status=st, lengths=sqrt(2._R8P) * h * [1._R8P, 2._R8P, 3._R8P])
      call stat(st)
      call table_errors(t, 2_I4P, EXACT, k, 1_I4P)
      call table_errors(t, 2_I4P, CH, k, 1_I4P)
      curv(k, 1) = maxval(abs(t%kappa(1, :) * R - 1._R8P))
      call fd_closure_2d(fr, xs, ys, [2_I4P, 2_I4P], t, status=st, lengths=sqrt(2._R8P) * h * [1._R8P, 2._R8P, 3._R8P])
      call stat(st)
      call table_errors(t, 2_I4P, FROZEN, k, 1_I4P)
      ! the body at a pose, spinning about the pose's origin and moving
      pose = pose_2d(0.4_R8P, [0.03_R8P, -0.02_R8P]) ; w = 2._R8P
      call fd_closure_2d(cv, xs, ys, [2_I4P, 2_I4P], t, status=st, pose=pose, linear_velocity=[0.1_R8P, -0.2_R8P], &
                         angular_velocity=w)
      call stat(st)
      c = pose%t
      do g = 1_I4P, t%n
         rv = t%foot(:, g) - c
         vf = [0.1_R8P - w * rv(2), -0.2_R8P + w * rv(1)]
         rv = R * rv / norm2(rv)                                        ! the true wall point
         vt = [0.1_R8P - w * rv(2), -0.2_R8P + w * rv(1)]
         vel1(k, 1) = max(vel1(k, 1), maxval(abs(t%velocity(:, g) - vf)))
         vel2(k, 1) = max(vel2(k, 1), norm2(t%velocity(:, g) - vt))
      enddo
      vtol(k, 1) = abs(w) * sagc(k, 1) * (1._R8P + 1.e-6_R8P) + 1.e-14_R8P
   enddo
   endsubroutine study_2d

   subroutine study_3d()
   !< The 3-D study: 5 levels, three regimes, curvature and velocity checks.
   type(surface_stl_object)        :: sc, fr
   type(facet_object), allocatable :: fc(:)
   type(ib_closure_t)              :: t
   type(rigid_pose_t)              :: pose
   real(R8P)                       :: h, c(3), rv(3), vf(3), vt(3), w(3), v0(3)
   integer(I4P)                    :: k, st, g

   call icosphere(LEVEL0_3D, R, fc, hsf(2), sagf(2)) ; call fr%adopt_facets(facets=fc)
   do k = 1_I4P, NLV(2)
      h = R / (4._R8P * 2._R8P**(k - 1)) ; hh(k, 2) = h
      call make_lines(3_I4P, h)
      call icosphere(LEVEL0_3D + k - 1_I4P, R, fc, hsc(k, 2), sagc(k, 2))
      call sc%adopt_facets(facets=fc) ; call sc%build_curvature(status=st) ; call stat(st)
      call fd_closure_3d(sc, xs, ys, zs, [2_I4P, 2_I4P, 2_I4P], t, status=st, &
                         lengths=sqrt(2._R8P) * h * [1._R8P, 2._R8P, 3._R8P])
      call stat(st)
      call table_errors(t, 3_I4P, EXACT, k, 2_I4P)
      call table_errors(t, 3_I4P, CH, k, 2_I4P)
      curv(k, 2) = maxval(abs(t%kappa * R - 1._R8P))
      call fd_closure_3d(fr, xs, ys, zs, [2_I4P, 2_I4P, 2_I4P], t, status=st, &
                         lengths=sqrt(2._R8P) * h * [1._R8P, 2._R8P, 3._R8P])
      call stat(st)
      call table_errors(t, 3_I4P, FROZEN, k, 2_I4P)
      ! the body at a pose, spinning about the pose's origin and moving
      pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.03_R8P, -0.02_R8P, 0.01_R8P])
      w = [0.5_R8P, -1._R8P, 2._R8P] ; v0 = [0.1_R8P, -0.2_R8P, 0.3_R8P]
      call fd_closure_3d(sc, xs, ys, zs, [2_I4P, 2_I4P, 2_I4P], t, status=st, pose=pose, linear_velocity=v0, &
                         angular_velocity=w)
      call stat(st)
      c = pose%t
      do g = 1_I4P, t%n
         rv = t%foot(:, g) - c
         vf = v0 + cross(w, rv)
         rv = R * rv / norm2(rv)                                        ! the true wall point
         vt = v0 + cross(w, rv)
         vel1(k, 2) = max(vel1(k, 2), maxval(abs(t%velocity(:, g) - vf)))
         vel2(k, 2) = max(vel2(k, 2), norm2(t%velocity(:, g) - vt))
      enddo
      vtol(k, 2) = norm2(w) * sagc(k, 2) * (1._R8P + 1.e-6_R8P) + 1.e-14_R8P
   enddo
   endsubroutine study_3d

   subroutine table_errors(t, dims, reg, k, id)
   !< The ghost errors of one table in regime `reg` at level k (id = 1: 2-D, 2: 3-D).
   type(ib_closure_t), intent(in) :: t         !< Closure table.
   integer(I4P),       intent(in) :: dims, reg, k, id !< Dimensions, regime, level, study.
   real(R8P)    :: p(3), foot(3), nrm(3), eta, prb(3, NM), pw(3), gd, gn, ue, ui4(NM), ui5(NM), ui(NM), fv(0:NM)
   real(R8P)    :: tn(0:NM), l0, q
   integer(I4P) :: cel(3, NM), g, kk, m
   logical      :: bad

   ngh(reg, k, id) = t%n
   do g = 1_I4P, t%n
      p = 0._R8P ; p(1:dims) = t%points(:, g)
      if (reg == EXACT) then
         nrm = p / norm2(p) ; foot = R * nrm ; eta = norm2(p) - R
         do kk = 1_I4P, NM
            prb(:, kk) = foot + t%lengths(kk) * nrm
            cel(1, kk) = locate(xs, prb(1, kk)) ; cel(2, kk) = locate(ys, prb(2, kk)) ; cel(3, kk) = 1_I4P
            if (dims == 3_I4P) cel(3, kk) = locate(zs, prb(3, kk))
         enddo
         gd = field(dims, foot) ; gn = dot_product(grad(dims, foot), nrm)
      else
         foot = 0._R8P ; nrm = 0._R8P ; prb = 0._R8P ; cel = 1_I4P
         foot(1:dims) = t%foot(:, g) ; nrm(1:dims) = t%normal(:, g) ; eta = t%distance(g)
         prb(1:dims, :) = t%probe(:, :, g) ; cel(1:dims, :) = t%probe_cell(:, :, g)
         pw = R * foot / norm2(foot)                                 ! the true wall point
         gd = field(dims, pw) ; gn = dot_product(grad(dims, pw), pw / R)
      endif
      gmx(reg, k, id) = max(gmx(reg, k, id), abs(dot_product(grad(dims, R * foot / norm2(foot)), foot / norm2(foot))))
      ue = field(dims, p)
      bad = .false.
      do kk = 1_I4P, NM
         ui4(kk) = interp(dims, prb(:, kk), cel(:, kk), 4_I4P, bad)
         ui5(kk) = interp(dims, prb(:, kk), cel(:, kk), 5_I4P, bad)
      enddo
      if (bad) then
         nbad = nbad + 1_I4P ; cycle
      endif
      tn(0) = 0._R8P ; tn(1:NM) = t%lengths(1:NM)
      do m = 1_I4P, NM
         ui = ui4 ; if (m == 3_I4P) ui = ui5
         fv(0) = gd ; fv(1:m) = ui(1:m)
         call lagrange(tn(0:m), fv(0:m), eta, q, l0)
         err(m, DIR, reg, k, id) = max(err(m, DIR, reg, k, id), abs(q - ue))
         lam(m, reg, k, id) = max(lam(m, reg, k, id), abs(l0))
         q = neumann(tn(1:m), ui(1:m), gn, eta)
         err(m, NEU, reg, k, id) = max(err(m, NEU, reg, k, id), abs(q - ue))
      enddo
   enddo
   endsubroutine table_errors

   subroutine lagrange(tn, fv, s, q, l0)
   !< The interpolating polynomial through (tn(j), fv(j)) at s, and the weight of the first datum.
   real(R8P), intent(in)  :: tn(0:), fv(0:) !< Nodes, values.
   real(R8P), intent(in)  :: s              !< Abscissa.
   real(R8P), intent(out) :: q, l0          !< Value, weight of fv(0).
   real(R8P)              :: w
   integer(I4P)           :: j, i

   q = 0._R8P
   do j = 0, ubound(tn, 1)
      w = 1._R8P
      do i = 0, ubound(tn, 1)
         if (i /= j) w = w * (s - tn(i)) / (tn(j) - tn(i))
      enddo
      q = q + w * fv(j)
      if (j == 0) l0 = w
   enddo
   endsubroutine lagrange

   function neumann(tn, fv, gn, s) result(q)
   !< The polynomial q of degree m = size(tn) with q'(0) = gn and q(tn(k)) = fv(k), at s.
   real(R8P), intent(in) :: tn(:), fv(:), gn, s !< Probe lengths, values, normal derivative, abscissa.
   real(R8P)             :: q
   real(R8P)             :: a(size(tn) + 1, size(tn) + 1), b(size(tn) + 1), c(size(tn) + 1), sc, piv
   integer(I4P)          :: m, i, j, p

   m = size(tn) ; sc = tn(1)
   a = 0._R8P ; a(1, 2) = 1._R8P / sc ; b(1) = gn
   do i = 1, m
      do j = 0, m
         a(i + 1, j + 1) = (tn(i) / sc)**j
      enddo
      b(i + 1) = fv(i)
   enddo
   do i = 1, m + 1                                                    ! Gaussian elimination, partial pivoting
      p = i - 1 + maxloc(abs(a(i:, i)), dim=1)
      if (p /= i) then
         c = a(i, :) ; a(i, :) = a(p, :) ; a(p, :) = c
         piv = b(i) ; b(i) = b(p) ; b(p) = piv
      endif
      do j = i + 1, m + 1
         piv = a(j, i) / a(i, i)
         a(j, i:) = a(j, i:) - piv * a(i, i:) ; b(j) = b(j) - piv * b(i)
      enddo
   enddo
   do i = m + 1, 1, -1
      c(i) = (b(i) - dot_product(a(i, i + 1:), c(i + 1:))) / a(i, i)
   enddo
   q = 0._R8P
   do j = m, 0, -1
      q = q * (s / sc) + c(j + 1)
   enddo
   endfunction neumann

   function interp(dims, q, cel, np, bad) result(val)
   !< Tensor Lagrange interpolation of the grid values on np points per axis, placed from the enclosing cell.
   integer(I4P), intent(in)    :: dims, cel(3), np !< Dimensions, enclosing cell, points per axis.
   real(R8P),    intent(in)    :: q(3)             !< Point.
   logical,      intent(inout) :: bad              !< Set when the stencil leaves the block.
   real(R8P)                   :: val
   real(R8P)                   :: w(8, 3)
   integer(I4P)                :: s(3), i, j, kk

   s = 1_I4P ; w = 0._R8P ; w(1, 3) = 1._R8P ; val = 0._R8P
   call weights(xs, q(1), cel(1), np, s(1), w(:, 1), bad)
   call weights(ys, q(2), cel(2), np, s(2), w(:, 2), bad)
   if (dims == 3_I4P) call weights(zs, q(3), cel(3), np, s(3), w(:, 3), bad)
   if (bad) return
   if (dims == 2_I4P) then
      do j = 1, np ; do i = 1, np
         val = val + w(i, 1) * w(j, 2) * ug(s(1) + i - 1, s(2) + j - 1, 1)
      enddo ; enddo
   else
      do kk = 1, np ; do j = 1, np ; do i = 1, np
         val = val + w(i, 1) * w(j, 2) * w(kk, 3) * ug(s(1) + i - 1, s(2) + j - 1, s(3) + kk - 1)
      enddo ; enddo ; enddo
   endif
   endfunction interp

   subroutine weights(t, q, i, np, s, w, bad)
   !< The np Lagrange weights at q on the lines t from the cell i (t(i) < q <= t(i+1)): np even, the points i - np/2 + 1
   !< .. i + np/2; np odd, centred on the nearer end of the cell.
   real(R8P),    intent(in)    :: t(:), q  !< Lines, coordinate.
   integer(I4P), intent(in)    :: i, np    !< Cell, points.
   integer(I4P), intent(out)   :: s        !< First point.
   real(R8P),    intent(inout) :: w(:)     !< Weights.
   logical,      intent(inout) :: bad      !< Set when the stencil leaves the lines.
   integer(I4P)                :: j, kk

   s = 1_I4P
   if (i < 1_I4P .or. i >= size(t)) then
      bad = .true. ; return
   endif
   if (mod(np, 2_I4P) == 0_I4P) then
      s = i - (np / 2_I4P - 1_I4P)
   else
      s = i - (np - 1_I4P) / 2_I4P ; if (q - t(i) > t(i + 1) - q) s = s + 1_I4P
   endif
   if (s < 1_I4P .or. s + np - 1_I4P > size(t)) then
      bad = .true. ; return
   endif
   do j = 1, np
      w(j) = 1._R8P
      do kk = 1, np
         if (kk /= j) w(j) = w(j) * (q - t(s + kk - 1)) / (t(s + j - 1) - t(s + kk - 1))
      enddo
   enddo
   endsubroutine weights

   pure function locate(t, q) result(i)
   !< The enclosing cell by the low-side rule, t(i) < q <= t(i+1), 0 outside.
   real(R8P), intent(in) :: t(:), q !< Lines, coordinate.
   integer(I4P)          :: i

   i = 0_I4P
   if (.not. (q > t(1) .and. q <= t(size(t)))) return
   i = max(1_I4P, min(int(size(t), I4P) - 1_I4P, int((q - t(1)) / (t(2) - t(1)), I4P) + 1_I4P))
   do while (i > 1_I4P .and. q <= t(i))
      i = i - 1_I4P
   enddo
   do while (i < size(t) - 1 .and. q > t(i + 1))
      i = i + 1_I4P
   enddo
   endfunction locate

   subroutine make_lines(dims, h)
   !< Uniform point lines of spacing h covering the body, its probes and their stencils, shifted off the centre; the
   !< grid values of the field.
   integer(I4P), intent(in) :: dims !< Dimensions.
   real(R8P),    intent(in) :: h    !< Spacing.
   real(R8P)                :: b
   integer(I4P)             :: n, i, j, kk

   b = R + 0.05_R8P + (3._R8P * sqrt(2._R8P) + 6._R8P) * h
   n = ceiling(2._R8P * b / h) + 1_I4P
   xs = [(-b + (i - 1) * h + 0.0371_R8P * h, i = 1, n)]
   ys = [(-b + (i - 1) * h + 0.0613_R8P * h, i = 1, n)]
   if (dims == 3_I4P) then
      zs = [(-b + (i - 1) * h + 0.0829_R8P * h, i = 1, n)]
   else
      zs = [0._R8P]
   endif
   if (allocated(ug)) deallocate(ug)
   allocate(ug(size(xs), size(ys), size(zs)))
   do kk = 1, size(zs) ; do j = 1, size(ys) ; do i = 1, size(xs)
      ug(i, j, kk) = field(dims, [xs(i), ys(j), zs(kk)])
   enddo ; enddo ; enddo
   endsubroutine make_lines

   pure function field(dims, p) result(u)
   !< The manufactured field.
   integer(I4P), intent(in) :: dims !< Dimensions.
   real(R8P),    intent(in) :: p(3) !< Point (z ignored in 2-D).
   real(R8P)                :: u

   if (dims == 2_I4P) then
      u = cos(1.1_R8P * p(1) - 0.7_R8P * p(2) + 0.3_R8P) * exp(0.5_R8P * p(1)) + sin(0.9_R8P * p(2))
   else
      u = cos(1.1_R8P * p(1) - 0.7_R8P * p(2) + 0.3_R8P) * exp(0.6_R8P * p(3)) + sin(0.9_R8P * p(2) + 1.2_R8P * p(3))
   endif
   endfunction field

   pure function grad(dims, p) result(gr)
   !< The gradient of the manufactured field (z component 0 in 2-D).
   integer(I4P), intent(in) :: dims !< Dimensions.
   real(R8P),    intent(in) :: p(3) !< Point.
   real(R8P)                :: gr(3), a, e, b

   a = 1.1_R8P * p(1) - 0.7_R8P * p(2) + 0.3_R8P
   if (dims == 2_I4P) then
      e = exp(0.5_R8P * p(1))
      gr = [(-1.1_R8P * sin(a) + 0.5_R8P * cos(a)) * e, 0.7_R8P * sin(a) * e + 0.9_R8P * cos(0.9_R8P * p(2)), 0._R8P]
   else
      e = exp(0.6_R8P * p(3)) ; b = 0.9_R8P * p(2) + 1.2_R8P * p(3)
      gr = [-1.1_R8P * sin(a) * e, 0.7_R8P * sin(a) * e + 0.9_R8P * cos(b), 0.6_R8P * cos(a) * e + 1.2_R8P * cos(b)]
   endif
   endfunction grad

   subroutine report(id, title)
   !< Print the tables of one study and check its criteria.
   integer(I4P), intent(in) :: id    !< 1: 2-D, 2: 3-D.
   character(*), intent(in) :: title !< Title.
   real(R8P)                :: o, est, fl
   integer(I4P)             :: reg, k, m, bc, n
   character(160)           :: line
   logical                  :: pass

   n = NLV(id)
   print '(A)', ''
   print '(A)', title//', reach 2, probes l, 2 l, 3 l (l = sqrt(2) h), interpolation on 4 (m = 1, 2) / 5 (m = 3) points'
   print '(A)', 'max |u_ghost - u(ghost)| over the ghosts; D = Dirichlet, N = Neumann, m probes (design order m + 1)'
   do reg = 1_I4P, 3_I4P
      print '(A)', '  '//trim(REG_NAME(reg))
      print '(A)', '    h          h_s        ghosts   D m=1      D m=2      D m=3      N m=1      N m=2      N m=3'
      do k = 1_I4P, n
         if (reg == EXACT) then
            write(line, '(4X,ES9.3,2X,A9,I8,6(2X,ES9.3))') hh(k, id), '        -', ngh(reg, k, id), &
                                                           ((err(m, bc, reg, k, id), m = 1, NM), bc = 1, 2)
         elseif (reg == CH) then
            write(line, '(4X,ES9.3,2X,ES9.3,I8,6(2X,ES9.3))') hh(k, id), hsc(k, id), ngh(reg, k, id), &
                                                              ((err(m, bc, reg, k, id), m = 1, NM), bc = 1, 2)
         else
            write(line, '(4X,ES9.3,2X,ES9.3,I8,6(2X,ES9.3))') hh(k, id), hsf(id), ngh(reg, k, id), &
                                                              ((err(m, bc, reg, k, id), m = 1, NM), bc = 1, 2)
         endif
         print '(A)', trim(line)
      enddo
      write(line, '(4X,A,6(2X,F9.2))') 'order, two finest levels  ', &
                                       ((order(err(m, bc, reg, n - 1, id), err(m, bc, reg, n, id)), m = 1, NM), bc = 1, 2)
      print '(A)', trim(line)
   enddo
   ! criterion 1
   pass = .true.
   do bc = 1_I4P, 2_I4P
      do m = 1_I4P, NM
         o = order(err(m, bc, EXACT, n - 1, id), err(m, bc, EXACT, n, id))
         if (abs(o - (m + 1)) > 0.3_R8P) pass = .false.
      enddo
   enddo
   call check(title//': exact geometry, orders within 0.3 of m + 1 (D and N)', pass)
   ! criterion 2
   pass = .true.
   do m = 1_I4P, NM
      if (order(err(m, DIR, CH, n - 1, id), err(m, DIR, CH, n, id)) < 1.7_R8P) pass = .false.
   enddo
   call check(title//': h_s = c h, Dirichlet order >= 2 - 0.3', pass)
   ! criterion 3
   pass = .true.
   print '(A,ES9.3,A,ES9.3,A,ES9.3)', '  h_s frozen: h_s ', hsf(id), ', largest sag s ', sagf(id), &
                                      ', h_s^2 / (8 R) ', hsf(id)**2 / (8._R8P * R)
   do m = 1_I4P, NM
      o = order(err(m, DIR, FROZEN, n - 1, id), err(m, DIR, FROZEN, n, id))
      fl = err(m, DIR, FROZEN, n, id)
      est = lam(m, FROZEN, n, id) * gmx(FROZEN, n, id) * sagf(id)
      print '(A,I0,A,ES9.3,A,F5.2,A,ES9.3,A,F5.2,A,ES9.3,A,F5.2)', '    D m=', m, ': floor ', fl, ' (order ', o, &
            '); estimate Lambda |d_n u| s = ', est, ' (Lambda ', lam(m, FROZEN, n, id), ', |d_n u| ', &
            gmx(FROZEN, n, id), '); floor / estimate ', fl / est
      if (o >= 0.5_R8P .or. fl > est) pass = .false.
   enddo
   call check(title//': h_s frozen, the Dirichlet error stops at a floor below the O(h_s^2 / R) estimate', pass)
   ! criterion 4
   print '(A)', '  curvature at the feet (h_s = c h), max |k R - 1|:'
   print '(4X,8ES11.3)', curv(1:n, id)
   if (id == 1_I4P) then                                              ! the round-off of a three-point circle
      pass = all(curv(1:n, id) <= 100._R8P * epsilon(1._R8P) * (R / hsc(1:n, id))**2)
   else
      pass = curv(n, id) < 1.e-2_R8P .and. curv(n, id) < curv(n - 1, id)
   endif
   call check(title//': curvature at the feet', pass)
   ! criterion 5
   print '(A)', '  wall velocity, body at a pose: max |v - (V + omega x (x_f - c))|, max |v - v(true wall point)|, '// &
                '|omega| s:'
   do k = 1_I4P, n
      print '(4X,3ES11.3)', vel1(k, id), vel2(k, id), vtol(k, id)
   enddo
   call check(title//': wall velocity, the rigid field to round-off and to |omega| s', &
              all(vel1(1:n, id) <= 1.e-14_R8P) .and. all(vel2(1:n, id) <= vtol(1:n, id)))
   endsubroutine report

   pure function order(ea, eb) result(o)
   !< Observed order between two levels of halved spacing.
   real(R8P), intent(in) :: ea, eb !< Errors at h and h / 2.
   real(R8P)             :: o

   o = log(ea / eb) / log(2._R8P)
   endfunction order

   subroutine stat(st)
   !< Count a failed status.
   integer(I4P), intent(in) :: st !< Status.

   if (st /= STATUS_OK) nstat = nstat + 1_I4P
   endsubroutine stat

   subroutine check(label, cond)
   !< Report one criterion.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (cond) then
      print '(A)', 'ok: '//label
   else
      print '(A)', 'FAILED: '//label
      ok = .false.
   endif
   endsubroutine check

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   subroutine circle(n, cv)
   !< A circle of radius R as a closed curve of n segments, vertices on the circle (turned by 0.1 rad).
   integer(I4P),       intent(in)    :: n  !< Segments.
   type(curve_object), intent(inout) :: cv !< Curve.
   real(R8P)                         :: u(n), v(n)
   integer(I4P)                      :: i, st

   u = [(R * cos(0.1_R8P + 2._R8P * PI * i / n), i = 0, n - 1)]
   v = [(R * sin(0.1_R8P + 2._R8P * PI * i / n), i = 0, n - 1)]
   call cv%adopt_loops(u=u, v=v, loop_sizes=[n], status=st) ; call stat(st)
   endsubroutine circle

   subroutine icosphere(levels, radius, fo, edge, sag)
   !< Icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere; its longest edge and its
   !< largest sag (radius minus the least distance of a facet plane from the centre).
   integer(I4P),                    intent(in)  :: levels    !< Subdivisions.
   real(R8P),                       intent(in)  :: radius    !< Radius.
   type(facet_object), allocatable, intent(out) :: fo(:)     !< Facets.
   real(R8P),                       intent(out) :: edge, sag !< Longest edge, largest sag.
   real(R8P), allocatable  :: tri(:,:,:), nxt(:,:,:)
   real(R8P)               :: vv(3, 12), gr, nn(3), dmin
   integer(I4P), parameter :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                    2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                    4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                    5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   integer(I4P)            :: kk, l, m, j

   gr = (1._R8P + sqrt(5._R8P)) / 2._R8P
   vv = reshape([-1._R8P, gr, 0._R8P, 1._R8P, gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P, -gr, 0._R8P, &
                 0._R8P, -1._R8P, gr, 0._R8P, 1._R8P, gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P, -gr, &
                 gr, 0._R8P, -1._R8P, gr, 0._R8P, 1._R8P, -gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P], [3, 12])
   do kk = 1, 12
      vv(:, kk) = vv(:, kk) / norm2(vv(:, kk))
   enddo
   allocate(tri(3, 3, 20))
   do kk = 1, 20
      tri(:, :, kk) = vv(:, IDX(:, kk))
   enddo
   do l = 1, levels
      m = size(tri, dim=3)
      allocate(nxt(3, 3, 4 * m))
      do kk = 1, m
         nxt(:, :, 4 * kk - 3) = reshape([tri(:, 1, kk), mid(tri(:, 1, kk), tri(:, 2, kk)), &
                                          mid(tri(:, 3, kk), tri(:, 1, kk))], [3, 3])
         nxt(:, :, 4 * kk - 2) = reshape([mid(tri(:, 1, kk), tri(:, 2, kk)), tri(:, 2, kk), &
                                          mid(tri(:, 2, kk), tri(:, 3, kk))], [3, 3])
         nxt(:, :, 4 * kk - 1) = reshape([mid(tri(:, 3, kk), tri(:, 1, kk)), mid(tri(:, 2, kk), tri(:, 3, kk)), &
                                          tri(:, 3, kk)], [3, 3])
         nxt(:, :, 4 * kk)     = reshape([mid(tri(:, 1, kk), tri(:, 2, kk)), mid(tri(:, 2, kk), tri(:, 3, kk)), &
                                          mid(tri(:, 3, kk), tri(:, 1, kk))], [3, 3])
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   m = size(tri, dim=3)
   allocate(fo(m))
   edge = 0._R8P ; dmin = huge(1._R8P)
   do kk = 1, m
      do j = 1, 3
         fo(kk)%vertex(j) = vector_R8P(radius * tri(1, j, kk), radius * tri(2, j, kk), radius * tri(3, j, kk))
         edge = max(edge, radius * norm2(tri(:, mod(j, 3) + 1, kk) - tri(:, j, kk)))
      enddo
      call fo(kk)%compute_metrix
      nn = cross(tri(:, 2, kk) - tri(:, 1, kk), tri(:, 3, kk) - tri(:, 1, kk))
      dmin = min(dmin, abs(dot_product(nn, tri(:, 1, kk))) / norm2(nn))
   enddo
   sag = radius * (1._R8P - dmin)
   endsubroutine icosphere

   pure function mid(a, b) result(m)
   !< Unit midpoint.
   real(R8P), intent(in) :: a(3), b(3) !< Points.
   real(R8P)             :: m(3)       !< Midpoint on the unit sphere.

   m = (a + b) / norm2(a + b)
   endfunction mid
endprogram fossil_val_fd_closure
