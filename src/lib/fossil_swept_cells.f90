!< FOSSIL, swept geometry of a rigidly moving body on a fixed grid (issue #41 step 5).
module fossil_swept_cells
!< FOSSIL, swept geometry of a rigidly moving body on a fixed grid (issue #41 step 5).
!<
!< A conservative cut-cell update over a time step with a moving wall needs, per cell, the old and
!< new fluid volumes and, per grid face, the fluid aperture averaged over the step (Falcovitz et al.
!< 1997). On a fixed grid nothing sweeps through a face; the volume the wall sweeps through a cell is
!< its fluid volume change (Reynolds transport), so it is taken exactly from the cut cells at the two
!< poses (decision A of the #41 step-5 plan). The only time integral is the face aperture
!< (decision B); the time-averaged projected wall area of a cell follows from the closed-surface
!< identity, B = -sum_f A_f n_f, so the wall term uses the same quantities as the face fluxes.
!<
!< Time is normalized: s in [0, 1] from pose0 to pose1. The motion is the screw interpolation of
!< `pose_interpolate` (decision C), or the poses the caller gives at the Gauss nodes.
!<
!< Integration of the face apertures:
!< - 2-D pure translation (the two poses share their rotation): exact. On a copy at pose0 the motion
!<   is a translation by s D; a face is an axis-aligned segment moving by -s D in that frame, so its
!<   fluid length is piecewise linear in s. The breakpoints are the instants when a wall vertex
!<   crosses the face's line, and when one of the face's two end nodes crosses a wall segment; the
!<   two-point Gauss rule on every piece is then exact.
!< - 3-D pure translation: exact. On a copy at pose0 a grid plane x_ax = a moves to a - s d_ax and its
!<   grid lines by -s d. Every vertex of the fluid polygon of a face (a mesh edge on the plane, a grid
!<   line on a triangle, a grid corner) moves linearly in s while the combinatorics hold, so the area
!<   is quadratic in s. The combinatorics change only when a mesh vertex crosses the moving plane, a
!<   mesh edge crosses a moving grid line of the plane, or a grid node (moving by -s d) crosses a
!<   triangle. These instants, collected per plane, split [0, 1]; on every piece two-point Gauss is
!<   exact, and one evaluation (`slice_surface` of the shifted plane, 2-D cut cells on the shifted
!<   grid) gives every face of the plane.
!< - otherwise (rotation): Gauss-Legendre at nq nodes, the apertures from a copy moved to each node
!<   pose (2-D: exact line crossings of every face; 3-D: `cut_cells_3d` face fractions); with the screw
!<   motion, `error_estimate` is the largest face difference between nq and 2 nq nodes.
!<
!< Check (`residual`): the volume balance of the rigid velocity field of the screw motion, a field
!< that moves with the body (Falcovitz's design requirement: uniform flow moving with the body is
!< exact). Per cell, residual = swept + int_0^1 sum_f (flux of the field out of the fluid part of
!< face f) ds; the field is linear in x, so a face flux is its fluid length (area) times the field at
!< the fluid centroid. It is round-off for the exact translation paths and the quadrature error
!< otherwise. Not computed (0) when the caller gives the node poses. The exact paths assume a wall that
!< does not cross itself: on a self-intersecting curve (surface) the swept volumes come from the cut
!< cells (Green's theorem, regions counted by winding number) and the face apertures from line
!< crossings (parity), which count a knot differently, so the residual of the cells around the knot is
!< not round-off (a dragon section crossing itself once: 5.3e-6, the knot removed: 7.8e-18). Pass
!< `check_input=.true.` to have the input checked first (`sw%self_intersecting`, a warning on stderr;
!< issue #44), or call `self_intersections` / `check_validity(check_self_intersections=.true.)`.
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_curve_slice,        only : slice_surface
use fossil_cut_cells_3d,       only : cut_cells_3d_t, cut_cells_3d
use fossil_facet_object,       only : facet_object
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t, pose_interpolate, pose_angle_2d, pose_log_3d
use fossil_self_intersection,  only : intersection_pair_t
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P
use, intrinsic :: iso_fortran_env, only : stderr => error_unit

implicit none
private
public :: swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, swept_cells_3d, gauss_legendre_01

type :: swept_cells_2d_t
   !< Swept geometry of a curve between two poses on a 2-D block (fluid = outside the solid).
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P !< Cells.
   real(R8P),    allocatable :: alpha0(:,:)   !< Fluid area fraction at pose0 [nx, ny] (-1 if not closed).
   real(R8P),    allocatable :: alpha1(:,:)   !< Fluid area fraction at pose1 [nx, ny] (-1 if not closed).
   real(R8P),    allocatable :: swept(:,:)    !< Fluid area change (alpha1 - alpha0) |cell| [nx, ny].
   real(R8P),    allocatable :: face_u(:,:)   !< Time-averaged fluid length of the faces u = x_i [nx+1, ny].
   real(R8P),    allocatable :: face_v(:,:)   !< Time-averaged fluid length of the faces v = y_j [nx, ny+1].
   real(R8P),    allocatable :: wall(:,:,:)   !< Time-averaged projected wall, -sum_f A_f n_f [2, nx, ny].
   real(R8P),    allocatable :: residual(:,:) !< Volume balance of the rigid field [nx, ny].
   logical,      allocatable :: cleared(:,:)  !< No fluid at pose0, fluid at pose1 [nx, ny].
   logical,      allocatable :: covered(:,:)  !< Fluid at pose0, none at pose1 [nx, ny].
   logical                   :: exact = .false. !< Exact (2-D translation) integration used.
   real(R8P)                 :: error_estimate = 0._R8P !< Max face |A(nq) - A(2 nq)| (Gauss, screw motion).
   integer(I4P)              :: evaluations = 0_I4P !< Face-aperture evaluations per face (largest).
   logical                   :: self_intersecting = .false. !< The curve crosses itself (only with `check_input`).
endtype swept_cells_2d_t

type :: swept_cells_3d_t
   !< Swept geometry of a surface between two poses on a 3-D block (fluid = outside the solid).
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P, nz = 0_I4P !< Cells.
   real(R8P),    allocatable :: alpha0(:,:,:)   !< Fluid volume fraction at pose0 (-1 if not closed).
   real(R8P),    allocatable :: alpha1(:,:,:)   !< Fluid volume fraction at pose1 (-1 if not closed).
   real(R8P),    allocatable :: swept(:,:,:)    !< Fluid volume change (alpha1 - alpha0) |cell|.
   real(R8P),    allocatable :: face_x(:,:,:)   !< Time-averaged fluid area of the faces x = x_i [nx+1, ny, nz].
   real(R8P),    allocatable :: face_y(:,:,:)   !< Time-averaged fluid area of the faces y = y_j [nx, ny+1, nz].
   real(R8P),    allocatable :: face_z(:,:,:)   !< Time-averaged fluid area of the faces z = z_k [nx, ny, nz+1].
   real(R8P),    allocatable :: wall(:,:,:,:)   !< Time-averaged projected wall, -sum_f A_f n_f [3, nx, ny, nz].
   real(R8P),    allocatable :: residual(:,:,:) !< Volume balance of the rigid field.
   logical,      allocatable :: cleared(:,:,:)  !< No fluid at pose0, fluid at pose1.
   logical,      allocatable :: covered(:,:,:)  !< Fluid at pose0, none at pose1.
   logical                   :: exact = .false. !< Exact (translation) integration used.
   real(R8P)                 :: error_estimate = 0._R8P !< Max face |A(nq) - A(2 nq)| (Gauss, screw motion).
   integer(I4P)              :: evaluations = 0_I4P !< Face-aperture evaluations per grid plane (largest).
   logical                   :: self_intersecting = .false. !< The surface crosses itself (only with `check_input`).
endtype swept_cells_3d_t

real(R8P), parameter :: FLOOR_H = 1e-6_R8P !< Pieces shorter than this converge to 64 eps of the measure times FLOOR_H
                                         !< (issue #69: slices within 1e-10 of a face carry ~1e-7 relative noise).

type :: slist_t
   !< A list of breakpoints.
   real(R8P), allocatable :: s(:) !< Breakpoints.
endtype slist_t

contains
   subroutine swept_cells_2d(curve, pose0, pose1, x, y, sw, nq, poses, exact, status, check_input)
   !< Swept geometry of `curve` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1) (issue #41 step 5; see the module notes). `nq` (default 4) Gauss-Legendre nodes for the
   !< rotating case; `poses(nq)` the poses at the nodes (`gauss_legendre_01`) instead of the screw
   !< motion; `exact` (default .true.) uses the exact integration for a pure translation and, since issue #69,
   !< for the screw motion with a rotation (certified breakpoints, Gauss doubled to round-off on every piece;
   !< `sw%exact` false if a piece did not converge at 64 nodes).
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   !< `check_input` (default .false., issue #44): look for self-intersections of the curve first
   !< (`self_intersections`, one BVH pass); if any, `sw%self_intersecting` and a warning on stderr: the
   !< residual of the cells around a knot is not round-off (see the module notes). The result is
   !< computed anyway.
   type(curve_object),     intent(in)            :: curve    !< Curve (body frame).
   type(rigid_pose_2d_t),  intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_2d_t),  intent(in)            :: pose1    !< Pose at the end.
   real(R8P),              intent(in)            :: x(:), y(:) !< Cell edges (strictly increasing).
   type(swept_cells_2d_t), intent(out)           :: sw       !< Swept geometry.
   integer(I4P),           intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_2d_t),  intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                intent(in),  optional :: exact    !< Exact translation path (default .true.).
   integer(I4P),           intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                intent(in),  optional :: check_input !< Check the curve for self-intersections.
   integer(I4P), allocatable :: pairs(:,:)             !< Crossing pairs (check_input).
   type(curve_object)       :: c0, c1                !< Copies at the two poses.
   type(cut_cells_2d_t)     :: k0, k1                !< Their cut cells.
   real(R8P), allocatable   :: fu2(:,:), fv2(:,:), fl2(:,:), res(:,:) !< 2 nq face lengths, work, residual.
   integer(I4P)             :: nx, ny, n, i, j       !< Work.
   real(R8P)                :: d(2), dth             !< Displacement, rotation angle.
   logical                  :: translation, screw    !< Paths.
   logical                  :: rotation_exact        !< Exact rotation path requested (default).

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   n = 4_I4P ; if (present(nq)) n = nq
   if (nx < 1 .or. ny < 1 .or. n < 1) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(poses)) then
      if (size(poses) /= n) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   sw%nx = nx ; sw%ny = ny
   allocate(sw%alpha0(nx, ny), sw%alpha1(nx, ny), sw%swept(nx, ny), sw%face_u(nx + 1, ny), sw%face_v(nx, ny + 1), &
            sw%wall(2, nx, ny), sw%residual(nx, ny), sw%cleared(nx, ny), sw%covered(nx, ny))
   if (present(check_input)) then
      if (check_input) then
         call curve%self_intersections(pairs, n)
         sw%self_intersecting = n > 0_I4P
         if (sw%self_intersecting) write(stderr, '(A,I0,A)') 'WARNING: swept_cells_2d: the curve crosses itself (', n, &
                                                             ' segment pairs): residuals around the knots are not round-off'
         n = 4_I4P ; if (present(nq)) n = nq
      endif
   endif
   ! the volumes: exact cut cells at the two poses
   c0 = curve ; call c0%move(pose0) ; call c0%cut_cells(x, y, k0)
   c1 = curve ; call c1%move(pose1) ; call c1%cut_cells(x, y, k1)
   sw%alpha0 = k0%alpha ; sw%alpha1 = k1%alpha
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         sw%swept(i, j) = (sw%alpha1(i, j) - sw%alpha0(i, j)) * (x(i + 1) - x(i)) * (y(j + 1) - y(j))
      enddo
   enddo
   sw%cleared = sw%alpha0 == 0._R8P .and. sw%alpha1 > 0._R8P
   sw%covered = sw%alpha0 > 0._R8P .and. sw%alpha1 == 0._R8P
   ! the face apertures
   d = pose1%t - pose0%t
   dth = pose_angle_2d(pose1) - pose_angle_2d(pose0) ; dth = atan2(sin(dth), cos(dth))
   translation = all(pose1%r == pose0%r) .and. .not. present(poses)
   rotation_exact = .true.
   if (present(exact)) then
      translation = translation .and. exact ; rotation_exact = exact
   endif
   screw = .not. present(poses)
   allocate(res(nx, ny)) ; res = 0._R8P
   if (translation) then
      sw%exact = .true.
      call faces_translation_2d(c0, d, x, y, sw%face_u, sw%face_v, sw%evaluations)
      do j = 1_I4P, ny ! the field is uniform (d): only the fluid lengths enter
         do i = 1_I4P, nx
            res(i, j) = d(1) * (sw%face_u(i + 1, j) - sw%face_u(i, j)) + d(2) * (sw%face_v(i, j + 1) - sw%face_v(i, j))
         enddo
      enddo
   elseif (screw .and. rotation_exact .and. dth /= 0._R8P) then
      ! exact rotation (issue #69): certified breakpoints per face, Gauss doubled to round-off on every piece
      call faces_rotation_2d(curve, pose0, pose1, x, y, sw%face_u, sw%face_v, res, sw%evaluations, sw%exact)
   else
      call faces_gauss_2d(curve, pose0, pose1, x, y, n, sw%face_u, sw%face_v, res, poses)
      sw%evaluations = n
      if (screw) then
         allocate(fu2(nx + 1, ny), fv2(nx, ny + 1), fl2(nx, ny))
         call faces_gauss_2d(curve, pose0, pose1, x, y, 2_I4P * n, fu2, fv2, fl2)
         sw%error_estimate = max(maxval(abs(fu2 - sw%face_u)), maxval(abs(fv2 - sw%face_v)))
      else
         res = 0._R8P
      endif
   endif
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         sw%wall(:, i, j) = -[sw%face_u(i + 1, j) - sw%face_u(i, j), sw%face_v(i, j + 1) - sw%face_v(i, j)]
         sw%residual(i, j) = 0._R8P
         if (screw) sw%residual(i, j) = sw%swept(i, j) + res(i, j)
      enddo
   enddo
   endsubroutine swept_cells_2d

   subroutine faces_rotation_2d(curve, pose0, pose1, x, y, fu, fv, flux, nev, exact)
   !< Exact time-averaged fluid lengths and rigid-field outflows of every face under the screw motion of
   !< `pose_interpolate` with a rotation (issue #69; see the moving-bodies reference page, "Exact rotation"). A face
   !< outside the box swept by the body's bounding circle is fluid throughout (closed form). Every other face gets its
   !< breakpoints, the zeros of g(s) = (p0 + p1 s) cos phi + (q0 + q1 s) sin phi + k0 + k1 s isolated with the bound
   !< on |g''| (`isolate_zeros`): a wall vertex crossing the face's line inside its span, an end node crossing a wall
   !< segment inside it. On every piece between breakpoints the fluid length and the outflow are analytic: Gauss-
   !< Legendre doubled from 2 nodes until two values agree to round-off (at most 64; else exact = .false.).
   type(curve_object),    intent(in)  :: curve        !< Curve (body frame).
   type(rigid_pose_2d_t), intent(in)  :: pose0, pose1 !< End poses.
   real(R8P),             intent(in)  :: x(:), y(:)   !< Cell edges.
   real(R8P),             intent(out) :: fu(:,:), fv(:,:) !< Time-averaged fluid lengths.
   real(R8P),             intent(out) :: flux(:,:)    !< Time integral of the rigid-field outflow per cell.
   integer(I4P),          intent(out) :: nev          !< Largest number of evaluations per face.
   logical,               intent(out) :: exact        !< Every piece converged.
   real(R8P),    allocatable :: vu(:), vv(:), ou(:,:), ov(:,:), bp(:)
   integer(I4P), allocatable :: nxt(:)
   real(R8P)                 :: th0, dth, t0(2), dt(2), cen(2), rho, cbox(2, 2), uv(2), ss, cc(2), cmax
   integer(I4P)              :: nx, ny, nv, i, j, k, nb, ne
   logical                   :: ok

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   nv = curve%get_vertices_number()
   allocate(vu(nv), vv(nv), nxt(nv), ou(nx + 1, ny), ov(nx, ny + 1))
   do k = 1_I4P, nv
      uv = curve%get_vertex(k) ; vu(k) = uv(1) ; vv(k) = uv(2) ; nxt(k) = curve%get_next(k)
   enddo
   th0 = pose_angle_2d(pose0) ; dth = pose_angle_2d(pose1) - th0 ; dth = atan2(sin(dth), cos(dth))
   t0 = pose0%t ; dt = pose1%t - pose0%t
   ! the box swept by the body's bounding circle: centre path C(s) = R(th(s)) cen + t(s), |C''| <= dth^2 |cen|
   cen = 0.5_R8P * ([minval(vu), minval(vv)] + [maxval(vu), maxval(vv)])
   rho = sqrt(maxval((vu - cen(1))**2 + (vv - cen(2))**2))
   cbox(:, 1) = huge(1._R8P) ; cbox(:, 2) = -huge(1._R8P)
   do k = 0_I4P, 64_I4P
      ss = real(k, R8P) / 64._R8P
      cc = rot2(th0 + ss * dth, cen) + t0 + ss * dt
      cbox(:, 1) = min(cbox(:, 1), cc) ; cbox(:, 2) = max(cbox(:, 2), cc)
   enddo
   cmax = rho + dth**2 * norm2(cen) / (8._R8P * 64._R8P**2) + 1e-12_R8P * (rho + norm2(cen) + norm2(t0) + norm2(dt))
   cbox(:, 1) = cbox(:, 1) - cmax ; cbox(:, 2) = cbox(:, 2) + cmax
   nev = 0_I4P ; exact = .true.
   ! faces u = x_i, v in [y_j, y_j+1]
   do j = 1_I4P, ny
      do i = 1_I4P, nx + 1_I4P
         if (x(i) < cbox(1, 1) .or. x(i) > cbox(1, 2) .or. y(j + 1) < cbox(2, 1) .or. y(j) > cbox(2, 2)) then
            fu(i, j) = y(j + 1) - y(j) ! fluid throughout
            ou(i, j) = closed_outflow(2_I4P, y(j), y(j + 1))
            cycle
         endif
         call face_events(1_I4P, x(i), y(j), y(j + 1), bp, nb)
         call integrate_face(1_I4P, x(i), y(j), y(j + 1), bp, nb, fu(i, j), ou(i, j), ne, ok)
         nev = max(nev, ne) ; exact = exact .and. ok
      enddo
   enddo
   ! faces v = y_j, u in [x_i, x_i+1]
   do j = 1_I4P, ny + 1_I4P
      do i = 1_I4P, nx
         if (y(j) < cbox(2, 1) .or. y(j) > cbox(2, 2) .or. x(i + 1) < cbox(1, 1) .or. x(i) > cbox(1, 2)) then
            fv(i, j) = x(i + 1) - x(i)
            ov(i, j) = closed_outflow(1_I4P, x(i), x(i + 1))
            cycle
         endif
         call face_events(2_I4P, y(j), x(i), x(i + 1), bp, nb)
         call integrate_face(2_I4P, y(j), x(i), x(i + 1), bp, nb, fv(i, j), ov(i, j), ne, ok)
         nev = max(nev, ne) ; exact = exact .and. ok
      enddo
   enddo
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         flux(i, j) = ou(i + 1, j) - ou(i, j) + ov(i, j + 1) - ov(i, j)
      enddo
   enddo
   contains
      function closed_outflow(lax, lo, hi) result(o)
      !< Time integral of the rigid-field outflow through a face fluid throughout: the face along `lax` from lo to hi.
      !< Face u = x_i (lax = 2): int u_x dl = -dth (m - t_v(s) l) + d_u l; face v = y_j (lax = 1): dth (m - t_u l) +
      !< d_v l; with l = hi - lo, m = (hi^2 - lo^2) / 2 and the mean of t(s) over the step.
      integer(I4P), intent(in) :: lax
      real(R8P),    intent(in) :: lo, hi
      real(R8P)                :: o, l, m

      l = hi - lo ; m = 0.5_R8P * (hi - lo) * (hi + lo)
      if (lax == 2_I4P) then
         o = -dth * (m - (t0(2) + 0.5_R8P * dt(2)) * l) + dt(1) * l
      else
         o = dth * (m - (t0(1) + 0.5_R8P * dt(1)) * l) + dt(2) * l
      endif
      endfunction closed_outflow

      subroutine face_events(nax, a, lo, hi, b, m)
      !< Breakpoints in (0, 1) of the face on the line normal to `nax` (1: u = a, 2: v = a), span [lo, hi].
      integer(I4P),           intent(in)  :: nax
      real(R8P),              intent(in)  :: a, lo, hi
      real(R8P), allocatable, intent(out) :: b(:)
      integer(I4P),           intent(out) :: m
      real(R8P), allocatable              :: tmp(:), z(:)
      real(R8P)                           :: co(6), w0(2), w1(2), d(2), xs(2), r(2), lam, node(2)
      integer(I4P)                        :: q, nz, kk, e
      logical                             :: cert

      allocate(tmp(0)) ; m = 0_I4P
      ! wall vertices crossing the line: X_nax(s) = a, the other coordinate inside [lo, hi]
      do q = 1_I4P, nv
         if (nax == 1_I4P) then
            co = [vu(q), 0._R8P, -vv(q), 0._R8P, t0(1) - a, dt(1)]
         else
            co = [vv(q), 0._R8P, vu(q), 0._R8P, t0(2) - a, dt(2)]
         endif
         call isolate_zeros(co, th0, dth, abs(vu(q)) + abs(vv(q)) + norm2(t0) + abs(a) + norm2(dt), z, nz, cert)
         if (.not. cert) exact = .false.
         do kk = 1_I4P, nz
            xs = rot2(th0 + z(kk) * dth, [vu(q), vv(q)]) + t0 + z(kk) * dt
            if (xs(3 - nax) >= lo .and. xs(3 - nax) <= hi) tmp = [tmp, z(kk)]
         enddo
      enddo
      ! the end nodes crossing wall segments: orient2d(X_q, X_q', N) = (R d) x (N - t(s)) - d x q
      do e = 1_I4P, 2_I4P
         if (nax == 1_I4P) then
            node = [a, merge(lo, hi, e == 1_I4P)]
         else
            node = [merge(lo, hi, e == 1_I4P), a]
         endif
         w0 = node - t0 ; w1 = -dt
         do q = 1_I4P, nv
            if (nxt(q) == 0_I4P) cycle
            d = [vu(nxt(q)) - vu(q), vv(nxt(q)) - vv(q)]
            co = [d(1) * w0(2) - d(2) * w0(1), d(1) * w1(2) - d(2) * w1(1), &
                  -(d(1) * w0(1) + d(2) * w0(2)), -(d(1) * w1(1) + d(2) * w1(2)), &
                  -(d(1) * vv(q) - d(2) * vu(q)), 0._R8P]
            call isolate_zeros(co, th0, dth, norm2(d) * (norm2(w0) + norm2(w1) + abs(vu(q)) + abs(vv(q))), z, nz, cert)
            if (.not. cert) exact = .false.
            do kk = 1_I4P, nz
               xs = rot2(th0 + z(kk) * dth, [vu(q), vv(q)]) + t0 + z(kk) * dt
               r = rot2(th0 + z(kk) * dth, d)
               lam = dot_product(node - xs, r) / dot_product(r, r)
               if (lam >= 0._R8P .and. lam <= 1._R8P) tmp = [tmp, z(kk)]
            enddo
         enddo
      enddo
      m = size(tmp, kind=I4P)
      call sort_unique(tmp, m)
      call merge_close(tmp, m)
      b = tmp(1:m)
      endsubroutine face_events

      subroutine integrate_face(nax, a, lo, hi, b, m, avg, out, ne, conv)
      !< Time averages over [0, 1] of the fluid length and of the rigid-field outflow of one face, piece by piece.
      integer(I4P), intent(in)  :: nax
      real(R8P),    intent(in)  :: a, lo, hi
      real(R8P),    intent(in)  :: b(:)
      integer(I4P), intent(in)  :: m
      real(R8P),    intent(out) :: avg, out
      integer(I4P), intent(out) :: ne
      logical,      intent(out) :: conv
      real(R8P)                 :: sa, sb, l1, o1, l2, o2, scl, sco
      integer(I4P)              :: q, ng

      avg = 0._R8P ; out = 0._R8P ; ne = 0_I4P ; conv = .true.
      scl = hi - lo ; sco = scl * (abs(dth) * (maxval(abs([x, y])) + norm2(t0) + norm2(dt)) + norm2(dt))
      do q = 0_I4P, m
         sa = 0._R8P ; if (q > 0_I4P) sa = b(q)
         sb = 1._R8P ; if (q < m) sb = b(q + 1)
         if (sb <= sa) cycle
         ng = 2_I4P
         call piece(nax, a, lo, hi, sa, sb, ng, l1, o1) ; ne = ne + ng
         do
            ng = 2_I4P * ng
            call piece(nax, a, lo, hi, sa, sb, ng, l2, o2) ; ne = ne + ng
            if (abs(l2 - l1) <= 64._R8P * epsilon(1._R8P) * scl * max(sb - sa, FLOOR_H) .and. &
                abs(o2 - o1) <= 64._R8P * epsilon(1._R8P) * sco * max(sb - sa, FLOOR_H)) exit
            if (ng >= 64_I4P) then
               conv = .false. ; exit
            endif
            l1 = l2 ; o1 = o2
         enddo
         avg = avg + l2 ; out = out + o2
      enddo
      endsubroutine integrate_face

      subroutine piece(nax, a, lo, hi, sa, sb, ng, l, o)
      !< ng-point Gauss-Legendre over [sa, sb] of the fluid length and the outflow of one face (moved copies).
      integer(I4P), intent(in)  :: nax, ng
      real(R8P),    intent(in)  :: a, lo, hi, sa, sb
      real(R8P),    intent(out) :: l, o
      type(curve_object)        :: cq
      real(R8P), allocatable    :: sq(:), wq(:)
      real(R8P)                 :: ll, mm, sn, tq(2)
      integer(I4P)              :: q

      call gauss_legendre_01(ng, sq, wq)
      l = 0._R8P ; o = 0._R8P
      do q = 1_I4P, ng
         sn = sa + (sb - sa) * sq(q)
         cq = curve ; call cq%move(pose_interpolate(pose0, pose1, sn))
         tq = t0 + sn * dt
         if (nax == 1_I4P) then ! face u = a along v
            call face_measure(cq, 2_I4P, a, lo, hi, ll, mm)
            o = o + (sb - sa) * wq(q) * (-dth * (mm - tq(2) * ll) + dt(1) * ll)
         else                   ! face v = a along u
            call face_measure(cq, 1_I4P, a, lo, hi, ll, mm)
            o = o + (sb - sa) * wq(q) * (dth * (mm - tq(1) * ll) + dt(2) * ll)
         endif
         l = l + (sb - sa) * wq(q) * ll
      enddo
      endsubroutine piece
   endsubroutine faces_rotation_2d

   pure subroutine merge_close(a, m)
   !< Drop sorted breakpoints within 1e-10 of the previous kept one, of 0 or of 1 (see `isolate_zeros`).
   real(R8P),    intent(inout) :: a(:) !< Breakpoints, sorted, a(1:m) used.
   integer(I4P), intent(inout) :: m    !< Their number.
   real(R8P)                   :: last
   integer(I4P)                :: i, k

   k = 0_I4P ; last = 0._R8P
   do i = 1_I4P, m
      if (a(i) - last < 1e-10_R8P .or. a(i) > 1._R8P - 1e-10_R8P) cycle
      k = k + 1_I4P ; a(k) = a(i) ; last = a(i)
   enddo
   m = k
   endsubroutine merge_close

   pure function rot2(th, p) result(r)
   !< p rotated by th.
   real(R8P), intent(in) :: th, p(2)
   real(R8P)             :: r(2)

   r = [cos(th) * p(1) - sin(th) * p(2), sin(th) * p(1) + cos(th) * p(2)]
   endfunction rot2

   subroutine isolate_zeros(co, phi0, dphi, scale, z, nz, certified)
   !< Zeros in (0, 1) of g(s) = (p0 + p1 s) cos phi + (q0 + q1 s) sin phi + k0 + k1 s, phi = phi0 + s dphi, with
   !< co = [p0, p1, q0, q1, k0, k1] (issue #69). |g''| <= M = 2 |dphi| |(p1, q1)| + dphi^2 max|(p, q)| (the norm of
   !< an affine function is largest at an end): on [sa, sb] (h = sb - sa) with end values of one sign, both larger
   !< than M h^2/8 plus the evaluation error, g has no zero; otherwise the interval is halved. A sign change is
   !< refined by bisection to one ulp. An interval narrower than GAP = 1e-10 without a certificate (a zero at an end of
   !< the step, a tangency, g within round-off of zero) gives one breakpoint at its middle: a kink within GAP of a
   !< breakpoint changes the Gauss value of its piece by ~slope GAP^2, far below round-off. Zeros within GAP of 0 or
   !< 1, or of each other, are merged by the callers (`merge_close`). An interval on which |g| is at round-off
   !< throughout (end values and the curvature bound below `err`: a zero zone, e.g. a mesh face starting on a grid
   !< plane) gets breakpoints at its two ends only. A g with every coefficient at
   !< round-off of `scale` (identically zero: the element never leaves the face's line or plane) has no event. At most 4096
   !< intervals are examined: beyond, `certified` is false and the caller must not claim exactness.
   real(R8P),              intent(in)  :: co(6), phi0, dphi
   real(R8P),              intent(in)  :: scale !< Magnitude of the terms that make up the coefficients.
   real(R8P), allocatable, intent(out) :: z(:)
   integer(I4P),           intent(out) :: nz
   logical,                intent(out) :: certified
   integer(I4P), parameter             :: NSTACK = 4096, BUDGET = 4096
   real(R8P),    parameter             :: GAP = 1e-10_R8P !< Uncertified intervals stop here (see above).
   real(R8P)                           :: big_m, err, st(2, NSTACK), sa, sb, ga, gb, sm, gm, lo, hi, glo
   integer(I4P)                        :: top, used

   allocate(z(0)) ; nz = 0_I4P ; certified = .true. ; used = 0_I4P
   big_m = 2._R8P * abs(dphi) * norm2(co([2, 4])) + &
           dphi**2 * max(norm2([co(1), co(3)]), norm2([co(1) + co(2), co(3) + co(4)]))
   err = 32._R8P * epsilon(1._R8P) * max(scale, abs(co(1)) + abs(co(2)) + abs(co(3)) + abs(co(4)) + abs(co(5)) + &
                                          abs(co(6)))
   if (maxval(abs(co)) <= 64._R8P * epsilon(1._R8P) * scale) return ! identically zero (to round-off): no event
   top = 1_I4P ; st(:, 1) = [0._R8P, 1._R8P]
   do while (top > 0_I4P)
      used = used + 1_I4P
      if (used > BUDGET) then
         certified = .false. ; exit
      endif
      sa = st(1, top) ; sb = st(2, top) ; top = top - 1_I4P
      ga = g(sa) ; gb = g(sb)
      if (ga * gb < 0._R8P .and. abs(ga) > err .and. abs(gb) > err .and. big_m * (sb - sa)**2 / 8._R8P < err) then
         ! one sign change, g nearly linear here: bisection to one ulp
         lo = sa ; hi = sb ; glo = ga
         do while (nearest(lo, 1._R8P) < hi)
            sm = 0.5_R8P * (lo + hi)
            if (sm <= lo .or. sm >= hi) exit
            gm = g(sm)
            if (gm == 0._R8P) then
               lo = sm ; hi = sm ; exit
            endif
            if ((gm > 0._R8P) .eqv. (glo > 0._R8P)) then
               lo = sm ; glo = gm
            else
               hi = sm
            endif
         enddo
         call add(0.5_R8P * (lo + hi))
         cycle
      endif
      if (ga * gb > 0._R8P .and. min(abs(ga), abs(gb)) > big_m * (sb - sa)**2 / 8._R8P + err) cycle ! no zero
      if (abs(ga) <= err .and. abs(gb) <= err .and. big_m * (sb - sa)**2 / 8._R8P <= err) then
         ! a zero zone: |g| at round-off over the whole interval, the event undetermined in double precision
         call add(sa) ; call add(sb) ; cycle
      endif
      if (sb - sa < GAP) then
         call add(0.5_R8P * (sa + sb)) ; cycle
      endif
      if (top + 2_I4P > NSTACK) then ! cannot happen for |g''| bounded; keep a breakpoint rather than drop one
         call add(0.5_R8P * (sa + sb)) ; cycle
      endif
      sm = 0.5_R8P * (sa + sb)
      top = top + 1_I4P ; st(:, top) = [sm, sb]
      top = top + 1_I4P ; st(:, top) = [sa, sm]
   enddo
   contains
      function g(s) result(v)
      real(R8P), intent(in) :: s
      real(R8P)             :: v, ph

      ph = phi0 + s * dphi
      v = (co(1) + co(2) * s) * cos(ph) + (co(3) + co(4) * s) * sin(ph) + co(5) + co(6) * s
      endfunction g

      subroutine add(s)
      real(R8P), intent(in) :: s

      if (s > GAP .and. s < 1._R8P - GAP) then
         z = [z, s] ; nz = nz + 1_I4P
      endif
      endsubroutine add
   endsubroutine isolate_zeros

   subroutine faces_translation_2d(c0, d, x, y, fu, fv, nev)
   !< Exact time-averaged fluid lengths of every face under the translation s d of the curve c0
   !< (at pose0): piecewise linear in s between breakpoints, two-point Gauss on every piece.
   type(curve_object), intent(in)  :: c0        !< Curve at pose0.
   real(R8P),          intent(in)  :: d(2)      !< Displacement over the step.
   real(R8P),          intent(in)  :: x(:), y(:) !< Cell edges.
   real(R8P),          intent(out) :: fu(:,:), fv(:,:) !< Fluid lengths [nx+1, ny], [nx, ny+1].
   integer(I4P),       intent(out) :: nev       !< Largest number of evaluations per face.
   type(slist_t), allocatable :: node(:,:)      !< Node-path breakpoints [nx+1, ny+1].
   real(R8P),     allocatable :: vu(:), vv(:)   !< Vertices.
   integer(I4P),  allocatable :: nxt(:)         !< Segment ends.
   real(R8P),     allocatable :: bp(:)          !< Breakpoints of one face.
   real(R8P)                  :: uv(2)          !< Vertex.
   integer(I4P)               :: nx, ny, nv, i, j, k, nb !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   nv = c0%get_vertices_number()
   allocate(vu(nv), vv(nv), nxt(nv), node(nx + 1, ny + 1))
   do k = 1_I4P, nv
      uv = c0%get_vertex(k) ; vu(k) = uv(1) ; vv(k) = uv(2) ; nxt(k) = c0%get_next(k)
   enddo
   ! node paths N - s d against every segment
   do j = 1_I4P, ny + 1_I4P
      do i = 1_I4P, nx + 1_I4P
         call node_breakpoints([x(i), y(j)], d, vu, vv, nxt, node(i, j)%s)
      enddo
   enddo
   nev = 0_I4P
   ! faces u = x_i, v in [y_j, y_j+1]: vertices crossing u = x_i - s d1
   do j = 1_I4P, ny
      do i = 1_I4P, nx + 1_I4P
         call collect(1_I4P, x(i), y(j), y(j + 1), node(i, j)%s, node(i, j + 1)%s, bp, nb)
         call integrate(2_I4P, x(i), y(j), y(j + 1), [-d(1), -d(2)], bp, nb, fu(i, j))
         nev = max(nev, 2_I4P * (nb + 1_I4P))
      enddo
   enddo
   ! faces v = y_j, u in [x_i, x_i+1]: vertices crossing v = y_j - s d2
   do j = 1_I4P, ny + 1_I4P
      do i = 1_I4P, nx
         call collect(2_I4P, y(j), x(i), x(i + 1), node(i, j)%s, node(i + 1, j)%s, bp, nb)
         call integrate(1_I4P, y(j), x(i), x(i + 1), [-d(2), -d(1)], bp, nb, fv(i, j))
         nev = max(nev, 2_I4P * (nb + 1_I4P))
      enddo
   enddo
   contains
      subroutine collect(ax, a, lo, hi, n1, n2, b, m)
      !< Breakpoints of the face on the line (u = a for ax = 1, v = a for ax = 2), [lo, hi] along it.
      integer(I4P),           intent(in)  :: ax       !< Normal axis of the face.
      real(R8P),              intent(in)  :: a, lo, hi !< Line, span.
      real(R8P),              intent(in)  :: n1(:), n2(:) !< End-node breakpoints.
      real(R8P), allocatable, intent(out) :: b(:)     !< Sorted breakpoints in (0, 1).
      integer(I4P),           intent(out) :: m        !< Their number.
      real(R8P)                           :: dn, sk, w !< Motion along the normal, instant, position.
      real(R8P), allocatable              :: tmp(:)   !< Work.
      integer(I4P)                        :: q        !< Vertex.

      allocate(tmp(nv + size(n1) + size(n2)))
      m = 0_I4P
      dn = d(ax)
      if (dn /= 0._R8P) then
         do q = 1_I4P, nv ! vertex q on the moving line: (u_q or v_q) = a - s dn
            if (ax == 1_I4P) then
               sk = (a - vu(q)) / dn ; w = vv(q) + sk * d(2)
            else
               sk = (a - vv(q)) / dn ; w = vu(q) + sk * d(1)
            endif
            if (sk > 0._R8P .and. sk < 1._R8P .and. w >= lo .and. w <= hi) then
               m = m + 1_I4P ; tmp(m) = sk
            endif
         enddo
      endif
      do q = 1_I4P, size(n1, kind=I4P)
         m = m + 1_I4P ; tmp(m) = n1(q)
      enddo
      do q = 1_I4P, size(n2, kind=I4P)
         m = m + 1_I4P ; tmp(m) = n2(q)
      enddo
      call sort_unique(tmp, m)
      b = tmp(1:m)
      endsubroutine collect

      subroutine integrate(lax, a, lo, hi, shift, b, m, avg)
      !< Time average over [0, 1] of the fluid length of the face (line along `lax` at `a`, [lo, hi]),
      !< the face moving by s shift (shift(1) normal to it, shift(2) along it) in the frame of c0.
      integer(I4P), intent(in)  :: lax      !< Axis along the face (1: u, 2: v).
      real(R8P),    intent(in)  :: a, lo, hi !< Line, span.
      real(R8P),    intent(in)  :: shift(2) !< Motion of the face per unit s (normal, along).
      real(R8P),    intent(in)  :: b(:)     !< Breakpoints.
      integer(I4P), intent(in)  :: m        !< Their number.
      real(R8P),    intent(out) :: avg      !< Time average.
      real(R8P), parameter      :: G = 0.5_R8P / sqrt(3._R8P) !< Two-point Gauss offset on a unit piece.
      real(R8P)                 :: sa, sb, sm, h, l1, l2, mm !< Work.
      integer(I4P)              :: q        !< Piece.

      avg = 0._R8P
      do q = 0_I4P, m
         sa = 0._R8P ; if (q > 0_I4P) sa = b(q)
         sb = 1._R8P ; if (q < m) sb = b(q + 1)
         h = sb - sa
         if (h <= 0._R8P) cycle
         sm = 0.5_R8P * (sa + sb)
         call face_measure(c0, lax, a + (sm - G * h) * shift(1), lo + (sm - G * h) * shift(2), &
                           hi + (sm - G * h) * shift(2), l1, mm)
         call face_measure(c0, lax, a + (sm + G * h) * shift(1), lo + (sm + G * h) * shift(2), &
                           hi + (sm + G * h) * shift(2), l2, mm)
         avg = avg + 0.5_R8P * h * (l1 + l2)
      enddo
      endsubroutine integrate
   endsubroutine faces_translation_2d

   pure subroutine node_breakpoints(p, d, vu, vv, nxt, b)
   !< Instants s in (0, 1) when the point p - s d lies on a segment of the curve.
   real(R8P),              intent(in)  :: p(2), d(2)  !< Node, displacement.
   real(R8P),              intent(in)  :: vu(:), vv(:) !< Vertices.
   integer(I4P),           intent(in)  :: nxt(:)      !< Segment ends.
   real(R8P), allocatable, intent(out) :: b(:)        !< Instants.
   real(R8P), allocatable              :: tmp(:)      !< Work.
   real(R8P)                           :: ex, ey, ax, ay, det, s, lam, plo(2), phi(2) !< Work.
   integer(I4P)                        :: k, m        !< Work.

   allocate(tmp(size(vu)))
   m = 0_I4P
   plo = min(p, p - d) ; phi = max(p, p - d)
   do k = 1_I4P, size(vu, kind=I4P)
      if (nxt(k) <= 0_I4P) cycle
      if (max(vu(k), vu(nxt(k))) < plo(1) .or. min(vu(k), vu(nxt(k))) > phi(1) .or. &
          max(vv(k), vv(nxt(k))) < plo(2) .or. min(vv(k), vv(nxt(k))) > phi(2)) cycle
      ex = vu(nxt(k)) - vu(k) ; ey = vv(nxt(k)) - vv(k)
      ax = vu(k) - p(1) ; ay = vv(k) - p(2)
      ! p - s d = a + lam e  <=>  s (-d) - lam e = a - p
      det = d(1) * ey - ex * d(2)
      if (det == 0._R8P) cycle
      s = (ex * ay - ey * ax) / det
      lam = (d(2) * ax - d(1) * ay) / det
      if (s > 0._R8P .and. s < 1._R8P .and. lam >= 0._R8P .and. lam <= 1._R8P) then
         m = m + 1_I4P ; tmp(m) = s
      endif
   enddo
   b = tmp(1:m)
   endsubroutine node_breakpoints

   pure subroutine sort_unique(a, m)
   !< Sort a(1:m) ascending (insertion) and drop exact duplicates.
   real(R8P),    intent(inout) :: a(:) !< Values.
   integer(I4P), intent(inout) :: m    !< Count.
   real(R8P)                   :: v    !< Work.
   integer(I4P)                :: i, j, k !< Work.

   do i = 2_I4P, m
      v = a(i) ; j = i - 1_I4P
      do while (j >= 1_I4P)
         if (a(j) <= v) exit
         a(j + 1) = a(j) ; j = j - 1_I4P
      enddo
      a(j + 1) = v
   enddo
   k = min(m, 1_I4P)
   do i = 2_I4P, m
      if (a(i) /= a(k)) then
         k = k + 1_I4P ; a(k) = a(i)
      endif
   enddo
   m = k
   endsubroutine sort_unique

   subroutine face_measure(c, lax, a, lo, hi, len, mom)
   !< Fluid length and first moment (int s ds over the fluid part) of the axis-aligned segment along
   !< `lax` (1: u at v = a; 2: v at u = a), parameter in [lo, hi]: exact crossings, the state from the
   !< parity from below the curve box (fluid outside).
   type(curve_object), intent(in)  :: c        !< Curve.
   integer(I4P),       intent(in)  :: lax      !< Axis along the segment.
   real(R8P),          intent(in)  :: a, lo, hi !< Line, span.
   real(R8P),          intent(out) :: len, mom !< Fluid length, first moment.
   real(R8P),    allocatable :: t(:)   !< Crossings.
   integer(I4P), allocatable :: g(:)   !< Segments.
   logical,      allocatable :: e(:)   !< Entering flags.
   real(R8P)                 :: bmin(2), start, p0, p1 !< Work.
   integer(I4P)              :: n, q   !< Work.
   logical                   :: fluid  !< State.

   bmin = c%get_bmin()
   start = min(lo, bmin(lax)) - 1._R8P
   call c%line_crossings(lax, a, start, hi, t, g, e, n)
   len = 0._R8P ; mom = 0._R8P
   fluid = .true. ; p0 = start
   do q = 1_I4P, n + 1_I4P
      p1 = hi ; if (q <= n) p1 = t(q)
      if (fluid) call add(max(p0, lo), min(p1, hi))
      if (q <= n) then
         fluid = .not. e(q) ; p0 = t(q)
      endif
   enddo
   contains
      subroutine add(sa, sb)
      !< Add the fluid interval [sa, sb] if not empty.
      real(R8P), intent(in) :: sa, sb !< Interval.

      if (sb > sa) then
         len = len + (sb - sa) ; mom = mom + 0.5_R8P * (sb - sa) * (sb + sa)
      endif
      endsubroutine add
   endsubroutine face_measure

   subroutine faces_gauss_2d(curve, pose0, pose1, x, y, n, fu, fv, flux, poses)
   !< Time-averaged fluid lengths by n-point Gauss-Legendre, and the time integral of the outflow of the
   !< rigid field of the screw motion per cell (`flux`, zero if `poses` given).
   type(curve_object),    intent(in)           :: curve      !< Curve (body frame).
   type(rigid_pose_2d_t), intent(in)           :: pose0, pose1 !< End poses.
   real(R8P),             intent(in)           :: x(:), y(:) !< Cell edges.
   integer(I4P),          intent(in)           :: n          !< Nodes.
   real(R8P),             intent(out)          :: fu(:,:), fv(:,:) !< Fluid lengths.
   real(R8P),             intent(out)          :: flux(:,:)  !< Outflow of the rigid field [nx, ny].
   type(rigid_pose_2d_t), intent(in), optional :: poses(:)   !< Node poses.
   type(curve_object)       :: cq                  !< Copy at a node.
   type(rigid_pose_2d_t)    :: pq                  !< Node pose.
   real(R8P), allocatable   :: sq(:), wq(:), ou(:,:), ov(:,:) !< Nodes, weights, face outflows.
   real(R8P)                :: l, m, d(2), dth, tq(2) !< Work.
   integer(I4P)             :: nx, ny, q, i, j     !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   call gauss_legendre_01(n, sq, wq)
   allocate(ou(nx + 1, ny), ov(nx, ny + 1))
   fu = 0._R8P ; fv = 0._R8P ; flux = 0._R8P
   d = pose1%t - pose0%t
   dth = pose_angle_2d(pose1) - pose_angle_2d(pose0) ; dth = atan2(sin(dth), cos(dth))
   do q = 1_I4P, n
      if (present(poses)) then
         pq = poses(q)
      else
         pq = pose_interpolate(pose0, pose1, sq(q))
      endif
      tq = pose0%t + sq(q) * d
      cq = curve ; call cq%move(pq)
      ! velocity of the screw motion: dth J (x - t(s)) + d, J = [[0, -1], [1, 0]]
      do j = 1_I4P, ny
         do i = 1_I4P, nx + 1_I4P
            call face_measure(cq, 2_I4P, x(i), y(j), y(j + 1), l, m)
            fu(i, j) = fu(i, j) + wq(q) * l
            ou(i, j) = -dth * (m - tq(2) * l) + d(1) * l ! int u_x dl over the fluid part
         enddo
      enddo
      do j = 1_I4P, ny + 1_I4P
         do i = 1_I4P, nx
            call face_measure(cq, 1_I4P, y(j), x(i), x(i + 1), l, m)
            fv(i, j) = fv(i, j) + wq(q) * l
            ov(i, j) = dth * (m - tq(1) * l) + d(2) * l ! int u_y dl over the fluid part
         enddo
      enddo
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            flux(i, j) = flux(i, j) + wq(q) * (ou(i + 1, j) - ou(i, j) + ov(i, j + 1) - ov(i, j))
         enddo
      enddo
   enddo
   endsubroutine faces_gauss_2d

   subroutine swept_cells_3d(surface, pose0, pose1, x, y, z, sw, nq, poses, exact, status, check_input)
   !< Swept geometry of `surface` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1), z(nz+1) (issue #41 step 5): volumes from `cut_cells_3d` at the two poses; face apertures
   !< exact for a pure translation and, since issue #69, for the screw motion with a rotation (`exact`, default
   !< .true.; see the module notes and the moving-bodies page), else by `nq`-point Gauss-Legendre (default 4) on
   !< copies moved to the screw poses or to `poses(nq)`.
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   type(surface_stl_object), intent(in)            :: surface  !< Surface (body frame).
   type(rigid_pose_t),       intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_t),       intent(in)            :: pose1    !< Pose at the end.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(swept_cells_3d_t),   intent(out)           :: sw       !< Swept geometry.
   integer(I4P),             intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_t),       intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                  intent(in),  optional :: exact    !< Exact translation path (default .true.).
   logical,                  intent(in),  optional :: check_input !< Check the surface for self-intersections
                                                                  !< (`find_self_intersections`; see swept_cells_2d).
   integer(I4P),             intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   type(surface_stl_object) :: s0                 !< Moved copy.
   type(intersection_pair_t), allocatable :: pairs(:) !< Crossing facet pairs (check_input).
   real(R8P)                :: d(3)               !< Displacement.
   logical                  :: translation        !< Exact path.
   logical                  :: rotation_exact     !< Exact rotation path (issue #69).
   type(cut_cells_3d_t)     :: k0, k1             !< Cut cells at the end poses.
   real(R8P), allocatable   :: fx2(:,:,:), fy2(:,:,:), fz2(:,:,:), flux(:,:,:), fl2(:,:,:) !< Work.
   integer(I4P)             :: nx, ny, nz, n, i, j, k, st !< Work.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   n = 4_I4P ; if (present(nq)) n = nq
   if (nx < 1 .or. ny < 1 .or. nz < 1 .or. n < 1) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny)) .or. any(z(2:) <= z(:nz))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(poses)) then
      if (size(poses) /= n) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   sw%nx = nx ; sw%ny = ny ; sw%nz = nz
   allocate(sw%swept(nx, ny, nz), sw%face_x(nx + 1, ny, nz), sw%face_y(nx, ny + 1, nz), sw%face_z(nx, ny, nz + 1), &
            sw%wall(3, nx, ny, nz), sw%residual(nx, ny, nz), flux(nx, ny, nz))
   if (present(check_input)) then
      if (check_input) then
         call surface%find_self_intersections(pairs)
         sw%self_intersecting = size(pairs) > 0
         if (sw%self_intersecting) write(stderr, '(A,I0,A)') 'WARNING: swept_cells_3d: the surface crosses itself (', &
                                                             size(pairs), ' facet pairs): residuals around them are not round-off'
      endif
   endif
   s0 = surface ; call s0%move(pose1) ; call cut_cells_3d(s0, x, y, z, k1, st)
   s0 = surface ; call s0%move(pose0) ; call cut_cells_3d(s0, x, y, z, k0, st)
   sw%alpha0 = k0%alpha ; sw%alpha1 = k1%alpha
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            sw%swept(i, j, k) = (sw%alpha1(i, j, k) - sw%alpha0(i, j, k)) * (x(i + 1) - x(i)) * (y(j + 1) - y(j)) * &
                                (z(k + 1) - z(k))
         enddo
      enddo
   enddo
   sw%cleared = sw%alpha0 == 0._R8P .and. sw%alpha1 > 0._R8P
   sw%covered = sw%alpha0 > 0._R8P .and. sw%alpha1 == 0._R8P
   translation = all(pose1%r == pose0%r) .and. .not. present(poses)
   rotation_exact = .not. present(poses) .and. .not. translation
   if (present(exact)) then
      translation = translation .and. exact ; rotation_exact = rotation_exact .and. exact
   endif
   if (rotation_exact) then
      ! exact rotation (issue #69): certified breakpoints per grid plane, Gauss doubled to round-off on every piece
      call faces_rotation_3d(surface, pose0, pose1, x, y, z, sw%face_x, sw%face_y, sw%face_z, flux, sw%evaluations, &
                             sw%exact)
   elseif (translation) then
      sw%exact = .true.
      d = pose1%t - pose0%t
      call faces_translation_3d(s0, d, x, y, z, sw%face_x, sw%face_y, sw%face_z, sw%evaluations)
      do k = 1_I4P, nz ! the field is uniform (d): only the fluid areas enter
         do j = 1_I4P, ny
            do i = 1_I4P, nx
               flux(i, j, k) = d(1) * (sw%face_x(i + 1, j, k) - sw%face_x(i, j, k)) + &
                               d(2) * (sw%face_y(i, j + 1, k) - sw%face_y(i, j, k)) + &
                               d(3) * (sw%face_z(i, j, k + 1) - sw%face_z(i, j, k))
            enddo
         enddo
      enddo
   else
      call faces_gauss_3d(surface, pose0, pose1, x, y, z, n, sw%face_x, sw%face_y, sw%face_z, flux, poses)
      sw%evaluations = n
   endif
   if (.not. present(poses) .and. .not. translation .and. .not. rotation_exact) then
      allocate(fx2(nx + 1, ny, nz), fy2(nx, ny + 1, nz), fz2(nx, ny, nz + 1), fl2(nx, ny, nz))
      call faces_gauss_3d(surface, pose0, pose1, x, y, z, 2_I4P * n, fx2, fy2, fz2, fl2)
      sw%error_estimate = max(maxval(abs(fx2 - sw%face_x)), maxval(abs(fy2 - sw%face_y)), maxval(abs(fz2 - sw%face_z)))
   endif
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            sw%wall(:, i, j, k) = -[sw%face_x(i + 1, j, k) - sw%face_x(i, j, k), sw%face_y(i, j + 1, k) - sw%face_y(i, j, k), &
                                    sw%face_z(i, j, k + 1) - sw%face_z(i, j, k)]
            sw%residual(i, j, k) = 0._R8P
            if (.not. present(poses)) sw%residual(i, j, k) = sw%swept(i, j, k) + flux(i, j, k)
         enddo
      enddo
   enddo
   endsubroutine swept_cells_3d

   subroutine faces_rotation_3d(surface, pose0, pose1, x, y, z, fx, fy, fz, flux, nev, exact)
   !< Exact time-averaged fluid face areas and rigid-field outflows under the screw motion with a rotation (issue #69;
   !< moving-bodies page, "Exact rotation"), one grid plane at a time. The body point w = R0 q moves to
   !< X(s) = Rot(s angle) w + t0 + s dt about the fixed unit axis k (Rodrigues: Rot(phi) w = A cos phi + B sin phi + C,
   !< C = k (k . w), A = w - C, B = k x w). A plane's breakpoints are the zeros (`isolate_zeros`) of: a mesh vertex
   !< on the plane, a mesh edge meeting one of its grid lines (det[D, e_l, N - t(s)] - det[d, Rot^T e_l, w]), one of
   !< its grid nodes on a triangle ((Rot n0) . (N - t(s)) - n0 . w_a), each checked to lie on the plane's grid and
   !< prefiltered by the box the element sweeps. A plane outside the box the body sweeps is fluid throughout (closed
   !< form). On every piece, Gauss-Legendre doubled from 2 nodes until the plane's face areas and outflows change by
   !< round-off (at most 64 nodes, else exact = .false.); one evaluation slices a moved copy (`slice_surface` with
   !< `on_plane_low` and the 2-D cut cells, as `cut_cells_3d`).
   type(surface_stl_object), intent(in)  :: surface        !< Surface (body frame).
   type(rigid_pose_t),       intent(in)  :: pose0, pose1   !< End poses.
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   real(R8P),                intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Time-averaged fluid areas.
   real(R8P),                intent(out) :: flux(:,:,:)    !< Time integral of the rigid-field outflow per cell.
   integer(I4P),             intent(out) :: nev            !< Largest number of evaluations per plane.
   logical,                  intent(out) :: exact          !< Every piece converged.
   type(facet_object), pointer :: f
   real(R8P), allocatable      :: w(:,:,:), vlo(:,:,:), vhi(:,:,:), ox(:,:,:), oy(:,:,:), oz(:,:,:), bp(:)
   real(R8P), allocatable      :: area(:,:), outf(:,:)
   real(R8P)                   :: kax(3), ang, t0(3), dt(3), om(3), blo(3), bhi(3), xs(3), mrg, g(3, 3)
   integer(I4P)                :: nf, nx, ny, nz, i, j, k, q, c, ne, nb
   logical                     :: ok

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   call pose_log_3d(pose0, pose1, kax, ang)
   t0 = pose0%t ; dt = pose1%t - pose0%t ; om = ang * kax
   nf = surface%get_facets_number()
   allocate(w(3, 3, nf), vlo(3, 3, nf), vhi(3, 3, nf))
   blo = huge(1._R8P) ; bhi = -huge(1._R8P)
   do i = 1_I4P, nf
      f => surface%facet_at(i)
      do c = 1, 3
         w(:, c, i) = matmul(pose0%r, [f%vertex(c)%x, f%vertex(c)%y, f%vertex(c)%z])
         ! box of the vertex path: 17 samples, |X''| <= ang^2 |w|
         vlo(:, c, i) = huge(1._R8P) ; vhi(:, c, i) = -huge(1._R8P)
         do q = 0, 16
            xs = at(w(:, c, i), real(q, R8P) / 16._R8P)
            vlo(:, c, i) = min(vlo(:, c, i), xs) ; vhi(:, c, i) = max(vhi(:, c, i), xs)
         enddo
         mrg = ang**2 * norm2(w(:, c, i)) / (8._R8P * 16._R8P**2) + 1e-12_R8P * (norm2(w(:, c, i)) + norm2(t0) + &
               norm2(dt) + 1._R8P)
         vlo(:, c, i) = vlo(:, c, i) - mrg ; vhi(:, c, i) = vhi(:, c, i) + mrg
         blo = min(blo, vlo(:, c, i)) ; bhi = max(bhi, vhi(:, c, i))
      enddo
   enddo
   allocate(ox(nx + 1, ny, nz), oy(nx, ny + 1, nz), oz(nx, ny, nz + 1))
   nev = 0_I4P ; exact = .true.
   g(:, 1) = [minval(x), minval(y), minval(z)] ; g(:, 2) = [maxval(x), maxval(y), maxval(z)]
   do i = 1_I4P, nx + 1_I4P ! planes x = x_i: (u, v) = (y, z)
      call one_plane(1_I4P, x(i), y, z, area, outf)
      fx(i, :, :) = area ; ox(i, :, :) = outf
   enddo
   do j = 1_I4P, ny + 1_I4P ! planes y = y_j: (u, v) = (z, x)
      call one_plane(2_I4P, y(j), z, x, area, outf)
      do k = 1_I4P, nz
         do i = 1_I4P, nx
            fy(i, j, k) = area(k, i) ; oy(i, j, k) = outf(k, i)
         enddo
      enddo
   enddo
   do k = 1_I4P, nz + 1_I4P ! planes z = z_k: (u, v) = (x, y)
      call one_plane(3_I4P, z(k), x, y, area, outf)
      fz(:, :, k) = area ; oz(:, :, k) = outf
   enddo
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            flux(i, j, k) = ox(i + 1, j, k) - ox(i, j, k) + oy(i, j + 1, k) - oy(i, j, k) + oz(i, j, k + 1) - oz(i, j, k)
         enddo
      enddo
   enddo
   contains
      pure function at(wv, sv) result(xv)
      !< Position at s of the body point with R0 q = wv.
      real(R8P), intent(in) :: wv(3), sv
      real(R8P)             :: xv(3), cv(3)

      cv = kax * dot_product(kax, wv)
      xv = (wv - cv) * cos(sv * ang) + cross(kax, wv) * sin(sv * ang) + cv + t0 + sv * dt
      endfunction at

      pure subroutine rodrigues(wv, a3, b3, c3)
      !< Rot(phi) wv = a3 cos phi + b3 sin phi + c3.
      real(R8P), intent(in)  :: wv(3)
      real(R8P), intent(out) :: a3(3), b3(3), c3(3)

      c3 = kax * dot_product(kax, wv) ; a3 = wv - c3 ; b3 = cross(kax, wv)
      endsubroutine rodrigues

      subroutine one_plane(na, lev, gu, gv, ar, ou)
      !< Time-averaged fluid areas ar(nu, nv) and outflows ou(nu, nv) of the faces of the plane normal to na at lev,
      !< whose in-plane axes (cyclic) carry the cell edges gu, gv.
      integer(I4P),           intent(in)  :: na
      real(R8P),              intent(in)  :: lev, gu(:), gv(:)
      real(R8P), allocatable, intent(out) :: ar(:,:), ou(:,:)
      real(R8P), allocatable              :: a1(:,:), o1(:,:), a2(:,:), o2(:,:)
      real(R8P)                           :: sa, sb, sca, sco
      integer(I4P)                        :: iu, iv, nu, nv, ng, m

      nu = size(gu, kind=I4P) - 1_I4P ; nv = size(gv, kind=I4P) - 1_I4P
      allocate(ar(nu, nv), ou(nu, nv))
      iu = mod(na, 3_I4P) + 1_I4P ; iv = mod(na + 1_I4P, 3_I4P) + 1_I4P
      if (lev < blo(na) .or. lev > bhi(na) .or. gu(nu + 1) < blo(iu) .or. gu(1) > bhi(iu) .or. &
          gv(nv + 1) < blo(iv) .or. gv(1) > bhi(iv)) then ! fluid throughout: closed form
         do m = 1_I4P, nv
            do ne = 1_I4P, nu
               ar(ne, m) = (gu(ne + 1) - gu(ne)) * (gv(m + 1) - gv(m))
               xs(na) = lev ; xs(iu) = 0.5_R8P * (gu(ne) + gu(ne + 1)) ; xs(iv) = 0.5_R8P * (gv(m) + gv(m + 1))
               ou(ne, m) = ar(ne, m) * (cross_k(om, xs - (t0 + 0.5_R8P * dt), na) + dt(na))
            enddo
         enddo
         return
      endif
      call plane_events(na, lev, gu, gv, iu, iv, bp, nb)
      ar = 0._R8P ; ou = 0._R8P ; ne = 0_I4P
      sca = maxval(gu(2:) - gu(:nu)) * maxval(gv(2:) - gv(:nv))
      sco = sca * (norm2(om) * (maxval(abs([x, y, z])) + norm2(t0) + norm2(dt)) + norm2(dt))
      do m = 0_I4P, nb
         sa = 0._R8P ; if (m > 0_I4P) sa = bp(m)
         sb = 1._R8P ; if (m < nb) sb = bp(m + 1)
         if (sb <= sa) cycle
         ng = 2_I4P
         call plane_piece(na, lev, gu, gv, iu, iv, sa, sb, ng, a1, o1) ; ne = ne + ng
         do
            ng = 2_I4P * ng
            call plane_piece(na, lev, gu, gv, iu, iv, sa, sb, ng, a2, o2) ; ne = ne + ng
            if (maxval(abs(a2 - a1)) <= 64._R8P * epsilon(1._R8P) * sca * max(sb - sa, FLOOR_H) .and. &
                maxval(abs(o2 - o1)) <= 64._R8P * epsilon(1._R8P) * sco * max(sb - sa, FLOOR_H)) exit
            if (ng >= 64_I4P) then
               exact = .false. ; exit
            endif
            a1 = a2 ; o1 = o2
         enddo
         ar = ar + a2 ; ou = ou + o2
      enddo
      nev = max(nev, ne)
      endsubroutine one_plane

      subroutine plane_piece(na, lev, gu, gv, iu, iv, sa, sb, ng, ar, ou)
      !< ng-point Gauss-Legendre over [sa, sb] of the fluid areas and outflows of one plane (moved copies).
      integer(I4P),           intent(in)  :: na, iu, iv, ng
      real(R8P),              intent(in)  :: lev, gu(:), gv(:), sa, sb
      real(R8P), allocatable, intent(out) :: ar(:,:), ou(:,:)
      type(surface_stl_object)            :: sq_
      type(curve_object)                  :: sec
      type(cut_cells_2d_t)                :: c2
      real(R8P), allocatable              :: sq(:), wq(:)
      real(R8P)                           :: sn, tq(3), a, cw(3)
      integer(I4P)                        :: q, ju, jv, st

      allocate(ar(size(gu) - 1, size(gv) - 1), ou(size(gu) - 1, size(gv) - 1)) ; ar = 0._R8P ; ou = 0._R8P
      call gauss_legendre_01(ng, sq, wq)
      do q = 1_I4P, ng
         sn = sa + (sb - sa) * sq(q)
         sq_ = surface ; call sq_%move(pose_interpolate(pose0, pose1, sn))
         tq = t0 + sn * dt
         call slice_surface(surface=sq_, normal_axis=na, level=lev, curve=sec, on_plane_low=.true.)
         call sec%cut_cells(x=gu, y=gv, cut=c2, status=st)
         do jv = 1_I4P, size(gv, kind=I4P) - 1_I4P
            do ju = 1_I4P, size(gu, kind=I4P) - 1_I4P
               a = c2%alpha(ju, jv) * (gu(ju + 1) - gu(ju)) * (gv(jv + 1) - gv(jv))
               cw(na) = lev ; cw(iu) = c2%centroid(1, ju, jv) ; cw(iv) = c2%centroid(2, ju, jv)
               ar(ju, jv) = ar(ju, jv) + (sb - sa) * wq(q) * a
               ou(ju, jv) = ou(ju, jv) + (sb - sa) * wq(q) * a * (cross_k(om, cw - tq, na) + dt(na))
            enddo
         enddo
      enddo
      endsubroutine plane_piece

      subroutine plane_events(na, lev, gu, gv, iu, iv, b, m)
      !< Breakpoints in (0, 1) of the plane normal to na at lev with grid gu x gv.
      integer(I4P),           intent(in)  :: na, iu, iv
      real(R8P),              intent(in)  :: lev, gu(:), gv(:)
      real(R8P), allocatable, intent(out) :: b(:)
      integer(I4P),           intent(out) :: m
      real(R8P), allocatable              :: tmp(:), zz(:)
      real(R8P)                           :: co(6), a3(3), b3(3), c3(3), ae(3), be(3), ce(3), d3(3), f0(3), f1(3)
      real(R8P)                           :: nn(3), xq(3), dd(3), lam, pp(3), tol, el(3), ebox(3, 2), tb(3, 2), bary(3)
      real(R8P)                           :: xa(3), xb(3), xc(3), wxd(3)
      integer(I4P)                        :: ii, cc, nz, kk, li, ldir, lfix, mm, e1, e2
      logical                             :: cert

      allocate(tmp(0))
      tol = 1e-12_R8P * (maxval(abs([x, y, z])) + 1._R8P)
      do ii = 1_I4P, nf
         tb(:, 1) = minval(vlo(:, :, ii), dim=2) ; tb(:, 2) = maxval(vhi(:, :, ii), dim=2)
         if (lev < tb(na, 1) .or. lev > tb(na, 2)) cycle
         ! (a) vertices on the plane
         do cc = 1, 3
            if (lev < vlo(na, cc, ii) .or. lev > vhi(na, cc, ii)) cycle
            call rodrigues(w(:, cc, ii), a3, b3, c3)
            co = [a3(na), 0._R8P, b3(na), 0._R8P, c3(na) + t0(na) - lev, dt(na)]
            call isolate_zeros(co, 0._R8P, ang, norm2(w(:, cc, ii)) + norm2(t0) + abs(lev) + norm2(dt), zz, nz, cert)
            if (.not. cert) exact = .false.
            do kk = 1_I4P, nz
               xq = at(w(:, cc, ii), zz(kk))
               if (xq(iu) >= gu(1) - tol .and. xq(iu) <= gu(size(gu)) + tol .and. &
                   xq(iv) >= gv(1) - tol .and. xq(iv) <= gv(size(gv)) + tol) tmp = [tmp, zz(kk)]
            enddo
         enddo
         ! (b) edges meeting grid lines of the plane: lines along iv at u = gu(li), along iu at v = gv(li)
         do cc = 1, 3
            e1 = cc ; e2 = mod(cc, 3) + 1
            ebox(:, 1) = min(vlo(:, e1, ii), vlo(:, e2, ii)) ; ebox(:, 2) = max(vhi(:, e1, ii), vhi(:, e2, ii))
            if (lev < ebox(na, 1) .or. lev > ebox(na, 2)) cycle
            d3 = w(:, e2, ii) - w(:, e1, ii) ; wxd = cross(w(:, e1, ii), d3)
            call rodrigues(d3, a3, b3, c3)
            do ldir = 1, 2
               if (ldir == 1) then
                  lfix = iu ; el = 0._R8P ; el(iv) = 1._R8P
               else
                  lfix = iv ; el = 0._R8P ; el(iu) = 1._R8P
               endif
               call rodrigues(el, ae, be, ce)
               do li = 1, merge(size(gu), size(gv), ldir == 1)
                  pp = 0._R8P ; pp(na) = lev ; pp(lfix) = merge(gu(li), gv(li), ldir == 1)
                  if (pp(lfix) < ebox(lfix, 1) .or. pp(lfix) > ebox(lfix, 2)) cycle
                  f0 = cross(el, pp - t0) ; f1 = cross(el, -dt)
                  co = [dot_product(a3, f0) - dot_product(ae, wxd), dot_product(a3, f1), &
                        dot_product(b3, f0) + dot_product(be, wxd), dot_product(b3, f1), &
                        dot_product(c3, f0) - dot_product(ce, wxd), dot_product(c3, f1)]
                  call isolate_zeros(co, 0._R8P, ang, norm2(d3) * (norm2(pp) + norm2(t0) + norm2(dt) + &
                                     norm2(w(:, e1, ii))), zz, nz, cert)
                  if (.not. cert) exact = .false.
                  do kk = 1_I4P, nz
                     xq = at(w(:, e1, ii), zz(kk)) ; dd = at(w(:, e2, ii), zz(kk)) - xq
                     if (dd(na) == 0._R8P) cycle
                     lam = (lev - xq(na)) / dd(na)
                     if (lam < -1e-12_R8P .or. lam > 1._R8P + 1e-12_R8P) cycle
                     pp = xq + lam * dd
                     if (abs(pp(lfix) - merge(gu(li), gv(li), ldir == 1)) > tol * 1e3_R8P) cycle
                     mm = merge(iv, iu, ldir == 1)
                     if (ldir == 1) then
                        if (pp(mm) < gv(1) - tol .or. pp(mm) > gv(size(gv)) + tol) cycle
                     else
                        if (pp(mm) < gu(1) - tol .or. pp(mm) > gu(size(gu)) + tol) cycle
                     endif
                     tmp = [tmp, zz(kk)]
                  enddo
               enddo
            enddo
         enddo
         ! (c) grid nodes of the plane on the triangle
         nn = cross(w(:, 2, ii) - w(:, 1, ii), w(:, 3, ii) - w(:, 1, ii))
         call rodrigues(nn, a3, b3, c3)
         do li = 1, size(gu)
            if (gu(li) < tb(iu, 1) .or. gu(li) > tb(iu, 2)) cycle
            do mm = 1, size(gv)
               if (gv(mm) < tb(iv, 1) .or. gv(mm) > tb(iv, 2)) cycle
               pp(na) = lev ; pp(iu) = gu(li) ; pp(iv) = gv(mm)
               f0 = pp - t0 ; f1 = -dt
               co = [dot_product(a3, f0), dot_product(a3, f1), dot_product(b3, f0), dot_product(b3, f1), &
                     dot_product(c3, f0) - dot_product(nn, w(:, 1, ii)), dot_product(c3, f1)]
               call isolate_zeros(co, 0._R8P, ang, norm2(nn) * (norm2(pp) + norm2(t0) + norm2(dt) + norm2(w(:, 1, ii))), &
                                  zz, nz, cert)
               if (.not. cert) exact = .false.
               do kk = 1_I4P, nz
                  xa = at(w(:, 1, ii), zz(kk)) ; xb = at(w(:, 2, ii), zz(kk)) ; xc = at(w(:, 3, ii), zz(kk))
                  el = cross(xb - xa, xc - xa)
                  if (dot_product(el, el) <= 0._R8P) cycle
                  bary = [dot_product(cross(xc - xb, pp - xb), el), dot_product(cross(xa - xc, pp - xc), el), &
                          dot_product(cross(xb - xa, pp - xa), el)] / dot_product(el, el)
                  if (all(bary >= -1e-12_R8P)) tmp = [tmp, zz(kk)]
               enddo
            enddo
         enddo
      enddo
      m = size(tmp, kind=I4P)
      call sort_unique(tmp, m)
      call merge_close(tmp, m)
      b = tmp(1:m)
      endsubroutine plane_events
   endsubroutine faces_rotation_3d

   pure function cross_k(a, b, k) result(c)
   !< Component k of a x b.
   real(R8P),    intent(in) :: a(3), b(3)
   integer(I4P), intent(in) :: k
   real(R8P)                :: c, v(3)

   v = cross(a, b) ; c = v(k)
   endfunction cross_k

   subroutine faces_translation_3d(s0, d, x, y, z, fx, fy, fz, nev)
   !< Exact time-averaged fluid areas of every face under the translation s d of the surface s0 (at
   !< pose0): per grid plane, quadratic in s between the breakpoints of `plane_breakpoints`, two-point
   !< Gauss on every piece.
   type(surface_stl_object), intent(in)  :: s0        !< Surface at pose0.
   real(R8P),                intent(in)  :: d(3)      !< Displacement over the step.
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   real(R8P),                intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Fluid areas.
   integer(I4P),             intent(out) :: nev       !< Largest number of evaluations per plane.
   type(facet_object), pointer :: fp                  !< Facet.
   real(R8P), allocatable      :: tri(:,:,:)          !< Triangles [3, 3, facets].
   real(R8P), allocatable      :: ue(:), ve(:), bp(:), area(:,:) !< In-plane edges, breakpoints, areas.
   real(R8P)                   :: a, dl(3)            !< Plane, displacement in the plane frame.
   integer(I4P)                :: nf, f, m, ax, a1, a2, ip, np, nb !< Work.

   nf = s0%get_facets_number()
   allocate(tri(3, 3, nf))
   do f = 1_I4P, nf
      fp => s0%facet_at(f)
      do m = 1_I4P, 3_I4P
         tri(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
   enddo
   nev = 0_I4P
   do ax = 1_I4P, 3_I4P
      a1 = mod(ax, 3_I4P) + 1_I4P ; a2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P ! in-plane axes, cyclic (slice_surface)
      select case (ax)
      case (1_I4P) ; np = size(x, kind=I4P) ; ue = y ; ve = z
      case (2_I4P) ; np = size(y, kind=I4P) ; ue = z ; ve = x
      case default ; np = size(z, kind=I4P) ; ue = x ; ve = y
      endselect
      dl = [d(ax), d(a1), d(a2)]
      do ip = 1_I4P, np
         select case (ax)
         case (1_I4P) ; a = x(ip)
         case (2_I4P) ; a = y(ip)
         case default ; a = z(ip)
         endselect
         call plane_breakpoints(tri, [ax, a1, a2], a, dl, ue, ve, bp, nb)
         call plane_average(s0, ax, a, dl, ue, ve, bp, nb, area)
         select case (ax)
         case (1_I4P) ; fx(ip, :, :) = area
         case (2_I4P) ; fy(:, ip, :) = transpose(area)
         case default ; fz(:, :, ip) = area
         endselect
         nev = max(nev, 2_I4P * (nb + 1_I4P))
      enddo
   enddo
   endsubroutine faces_translation_3d

   subroutine plane_breakpoints(tri, perm, a, dl, ue, ve, b, m)
   !< Instants s in (0, 1) when the combinatorics of the fluid polygons of the grid plane x_ax = a,
   !< moving by -s d in the frame of the surface (in-plane grid ue - s du, ve - s dv), change: a mesh
   !< vertex on the plane inside the block's in-plane extent, a mesh edge on a grid line of the plane,
   !< a grid node on a triangle. Coordinates in the plane frame (ax, a1, a2) = perm.
   real(R8P),              intent(in)  :: tri(:,:,:) !< Triangles [3, 3, facets], global axes.
   integer(I4P),           intent(in)  :: perm(3)    !< Plane frame axes (normal, u, v).
   real(R8P),              intent(in)  :: a          !< Plane.
   real(R8P),              intent(in)  :: dl(3)      !< Displacement in the plane frame.
   real(R8P),              intent(in)  :: ue(:), ve(:) !< In-plane grid edges.
   real(R8P), allocatable, intent(out) :: b(:)       !< Sorted unique breakpoints.
   integer(I4P),           intent(out) :: m          !< Their number.
   real(R8P), allocatable              :: tmp(:)     !< Work.
   real(R8P)                           :: w(3, 3), e1(3), e2(3), g(3), r(3), det, s, be, ga, slo, shi !< Work.
   integer(I4P)                        :: nu, nv, f, q, i, j, i1, i2, j1, j2 !< Work.

   nu = size(ue, kind=I4P) ; nv = size(ve, kind=I4P)
   allocate(tmp(64)) ; m = 0_I4P
   slo = a - max(dl(1), 0._R8P) ; shi = a - min(dl(1), 0._R8P) ! slab swept by the plane
   do f = 1_I4P, size(tri, 3, kind=I4P)
      do q = 1_I4P, 3_I4P
         w(:, q) = tri(perm, q, f)
      enddo
      if (maxval(w(1, :)) < slo .or. minval(w(1, :)) > shi) cycle
      ! mesh vertices crossing the plane
      if (dl(1) /= 0._R8P) then
         do q = 1_I4P, 3_I4P
            s = (a - w(1, q)) / dl(1)
            if (s > 0._R8P .and. s < 1._R8P .and. inside_extent(w(2, q) + s * dl(2), w(3, q) + s * dl(3))) call add(s)
         enddo
      endif
      ! mesh edges crossing the grid lines u = ue(i) - s du and v = ve(j) - s dv
      do q = 1_I4P, 3_I4P
         call edge_lines(w(:, q), w(:, mod(q, 3_I4P) + 1_I4P) - w(:, q), 2_I4P, 3_I4P, ue, ve)
         call edge_lines(w(:, q), w(:, mod(q, 3_I4P) + 1_I4P) - w(:, q), 3_I4P, 2_I4P, ve, ue)
      enddo
      ! grid nodes G - s d crossing the triangle: w1 + be e1 + ga e2 + s dl = G
      e1 = w(:, 2) - w(:, 1) ; e2 = w(:, 3) - w(:, 1)
      det = dot_product(e1, cross(e2, dl))
      if (det == 0._R8P) cycle
      call index_range(ue, minval(w(2, :)) + min(0._R8P, dl(2)), maxval(w(2, :)) + max(0._R8P, dl(2)), i1, i2)
      call index_range(ve, minval(w(3, :)) + min(0._R8P, dl(3)), maxval(w(3, :)) + max(0._R8P, dl(3)), j1, j2)
      do j = j1, j2
         do i = i1, i2
            g = [a, ue(i), ve(j)] ; r = g - w(:, 1)
            s = dot_product(e1, cross(e2, r)) / det
            if (s <= 0._R8P .or. s >= 1._R8P) cycle
            be = dot_product(r, cross(e2, dl)) / det ; ga = dot_product(e1, cross(r, dl)) / det
            if (be >= 0._R8P .and. ga >= 0._R8P .and. be + ga <= 1._R8P) call add(s)
         enddo
      enddo
   enddo
   call sort_unique(tmp, m)
   b = tmp(1:m)
   contains
      logical function inside_extent(pu, pv)
      !< The in-plane point (relative to the moving grid, i.e. + s d) inside the block's extent.
      real(R8P), intent(in) :: pu, pv !< Point.

      inside_extent = pu >= ue(1) .and. pu <= ue(nu) .and. pv >= ve(1) .and. pv <= ve(nv)
      endfunction inside_extent

      subroutine edge_lines(p, e, k, o, lk, lo)
      !< Crossings of the edge p + lam e with the moving grid lines w(k) = lk(i) - s dl(k) of the plane,
      !< the other in-plane coordinate (axis o) inside [lo(1), lo(end)] relative to the moving grid.
      real(R8P),    intent(in) :: p(3), e(3)    !< Edge.
      integer(I4P), intent(in) :: k, o          !< Line normal and line direction (plane frame).
      real(R8P),    intent(in) :: lk(:), lo(:)  !< Grid edges along k and along o.
      real(R8P)                :: dd, lam, sk, wo !< Work.
      integer(I4P)             :: ii, il, ih    !< Work.

      dd = e(1) * dl(k) - dl(1) * e(k)
      if (dd == 0._R8P) return
      call index_range(lk, min(p(k), p(k) + e(k)) + min(0._R8P, dl(k)), max(p(k), p(k) + e(k)) + max(0._R8P, dl(k)), il, ih)
      do ii = il, ih
         lam = ((a - p(1)) * dl(k) - dl(1) * (lk(ii) - p(k))) / dd
         sk = (e(1) * (lk(ii) - p(k)) - e(k) * (a - p(1))) / dd
         if (lam < 0._R8P .or. lam > 1._R8P .or. sk <= 0._R8P .or. sk >= 1._R8P) cycle
         wo = p(o) + lam * e(o) + sk * dl(o)
         if (wo >= lo(1) .and. wo <= lo(size(lo))) call add(sk)
      enddo
      endsubroutine edge_lines

      subroutine add(sv)
      !< Append a breakpoint.
      real(R8P), intent(in) :: sv !< Instant.

      if (m == size(tmp, kind=I4P)) tmp = [tmp, tmp]
      m = m + 1_I4P ; tmp(m) = sv
      endsubroutine add
   endsubroutine plane_breakpoints

   pure subroutine index_range(e, lo, hi, i1, i2)
   !< Indexes i1..i2 of the sorted edges e with lo <= e(i) <= hi (empty if i1 > i2), by bisection.
   real(R8P),    intent(in)  :: e(:)   !< Sorted edges.
   real(R8P),    intent(in)  :: lo, hi !< Range.
   integer(I4P), intent(out) :: i1, i2 !< Indexes.
   integer(I4P)              :: l, h, c !< Work.

   l = 1_I4P ; h = size(e, kind=I4P) + 1_I4P ! first index with e >= lo
   do while (l < h)
      c = (l + h) / 2_I4P
      if (e(c) < lo) then
         l = c + 1_I4P
      else
         h = c
      endif
   enddo
   i1 = l
   l = 0_I4P ; h = size(e, kind=I4P) ! last index with e <= hi
   do while (l < h)
      c = (l + h + 1_I4P) / 2_I4P
      if (e(c) > hi) then
         h = c - 1_I4P
      else
         l = c
      endif
   enddo
   i2 = l
   endsubroutine index_range

   subroutine plane_average(s0, ax, a, dl, ue, ve, b, m, area)
   !< Time-averaged fluid areas of the faces of the plane x_ax = a moving by -s d (plane frame dl):
   !< two-point Gauss on every piece between the breakpoints b(1:m), one section and one 2-D cut-cell
   !< evaluation per node.
   type(surface_stl_object), intent(in)  :: s0         !< Surface at pose0.
   integer(I4P),             intent(in)  :: ax         !< Plane normal axis.
   real(R8P),                intent(in)  :: a          !< Plane.
   real(R8P),                intent(in)  :: dl(3)      !< Displacement in the plane frame.
   real(R8P),                intent(in)  :: ue(:), ve(:) !< In-plane grid edges.
   real(R8P),                intent(in)  :: b(:)       !< Breakpoints.
   integer(I4P),             intent(in)  :: m          !< Their number.
   real(R8P), allocatable,   intent(out) :: area(:,:)  !< Time-averaged fluid areas [nu, nv].
   real(R8P), parameter     :: G = 0.5_R8P / sqrt(3._R8P) !< Two-point Gauss offset on a unit piece.
   type(curve_object)       :: sec                     !< Section.
   type(cut_cells_2d_t)     :: c2                      !< Its cut cells.
   real(R8P)                :: sa, sb, h, sg           !< Work.
   integer(I4P)             :: q, g2, j, k, st         !< Work.

   allocate(area(size(ue) - 1, size(ve) - 1)) ; area = 0._R8P
   do q = 0_I4P, m
      sa = 0._R8P ; if (q > 0_I4P) sa = b(q)
      sb = 1._R8P ; if (q < m) sb = b(q + 1)
      h = sb - sa
      if (h <= 0._R8P) cycle
      do g2 = -1_I4P, 1_I4P, 2_I4P
         sg = 0.5_R8P * (sa + sb) + g2 * G * h
         call slice_surface(surface=s0, normal_axis=ax, level=a - sg * dl(1), curve=sec, on_plane_low=.true.)
         call sec%cut_cells(x=ue - sg * dl(2), y=ve - sg * dl(3), cut=c2, status=st)
         area = area + 0.5_R8P * h * c2%alpha
      enddo
   enddo
   do k = 1_I4P, size(ve, kind=I4P) - 1_I4P
      do j = 1_I4P, size(ue, kind=I4P) - 1_I4P
         area(j, k) = area(j, k) * (ue(j + 1) - ue(j)) * (ve(k + 1) - ve(k))
      enddo
   enddo
   endsubroutine plane_average

   subroutine faces_gauss_3d(surface, pose0, pose1, x, y, z, n, fx, fy, fz, flux, poses)
   !< Time-averaged fluid face areas by n-point Gauss-Legendre (`cut_cells_3d` at each node pose), and
   !< the time integral of the outflow of the rigid field of the screw motion per cell.
   type(surface_stl_object), intent(in)           :: surface      !< Surface (body frame).
   type(rigid_pose_t),       intent(in)           :: pose0, pose1 !< End poses.
   real(R8P),                intent(in)           :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),             intent(in)           :: n            !< Nodes.
   real(R8P),                intent(out)          :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Fluid areas.
   real(R8P),                intent(out)          :: flux(:,:,:)  !< Outflow of the rigid field.
   type(rigid_pose_t),       intent(in), optional :: poses(:)     !< Node poses.
   type(surface_stl_object) :: sq_                 !< Copy at a node.
   type(cut_cells_3d_t)     :: kc                  !< Its cut cells.
   type(rigid_pose_t)       :: pq                  !< Node pose.
   real(R8P), allocatable   :: sq(:), wq(:), ox(:,:,:), oy(:,:,:), oz(:,:,:) !< Nodes, weights, face outflows.
   real(R8P)                :: d(3), om(3), axis(3), angle, tq(3), a, v(3) !< Work.
   integer(I4P)             :: nx, ny, nz, q, i, j, k, st !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   call gauss_legendre_01(n, sq, wq)
   allocate(ox(nx + 1, ny, nz), oy(nx, ny + 1, nz), oz(nx, ny, nz + 1))
   fx = 0._R8P ; fy = 0._R8P ; fz = 0._R8P ; flux = 0._R8P
   d = pose1%t - pose0%t
   call pose_log_3d(pose0, pose1, axis, angle) ; om = angle * axis
   do q = 1_I4P, n
      if (present(poses)) then
         pq = poses(q)
      else
         pq = pose_interpolate(pose0, pose1, sq(q))
      endif
      tq = pose0%t + sq(q) * d
      sq_ = surface ; call sq_%move(pq) ; call cut_cells_3d(sq_, x, y, z, kc, st)
      ! velocity of the screw motion: om x (x - t(s)) + d; a face flux is area * u(centroid) . n
      do k = 1_I4P, nz
         do j = 1_I4P, ny
            do i = 1_I4P, nx + 1_I4P
               a = kc%face_x(i, j, k) * (y(j + 1) - y(j)) * (z(k + 1) - z(k))
               fx(i, j, k) = fx(i, j, k) + wq(q) * a
               v = cross(om, kc%face_x_centroid(:, i, j, k) - tq) + d
               ox(i, j, k) = a * v(1)
            enddo
         enddo
      enddo
      do k = 1_I4P, nz
         do j = 1_I4P, ny + 1_I4P
            do i = 1_I4P, nx
               a = kc%face_y(i, j, k) * (z(k + 1) - z(k)) * (x(i + 1) - x(i))
               fy(i, j, k) = fy(i, j, k) + wq(q) * a
               v = cross(om, kc%face_y_centroid(:, i, j, k) - tq) + d
               oy(i, j, k) = a * v(2)
            enddo
         enddo
      enddo
      do k = 1_I4P, nz + 1_I4P
         do j = 1_I4P, ny
            do i = 1_I4P, nx
               a = kc%face_z(i, j, k) * (x(i + 1) - x(i)) * (y(j + 1) - y(j))
               fz(i, j, k) = fz(i, j, k) + wq(q) * a
               v = cross(om, kc%face_z_centroid(:, i, j, k) - tq) + d
               oz(i, j, k) = a * v(3)
            enddo
         enddo
      enddo
      do k = 1_I4P, nz
         do j = 1_I4P, ny
            do i = 1_I4P, nx
               flux(i, j, k) = flux(i, j, k) + wq(q) * (ox(i + 1, j, k) - ox(i, j, k) + oy(i, j + 1, k) - oy(i, j, k) + &
                                                        oz(i, j, k + 1) - oz(i, j, k))
            enddo
         enddo
      enddo
   enddo
   endsubroutine faces_gauss_3d

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   pure subroutine gauss_legendre_01(n, s, w)
   !< Gauss-Legendre nodes and weights on [0, 1] (Newton on P_n from the Chebyshev guesses).
   integer(I4P),           intent(in)  :: n    !< Nodes (>= 1).
   real(R8P), allocatable, intent(out) :: s(:) !< Nodes, ascending.
   real(R8P), allocatable, intent(out) :: w(:) !< Weights (sum 1).
   real(R8P), parameter                :: PI = acos(-1._R8P) !< Pi.
   real(R8P)                           :: z, z1, p0, p1, p2, dp !< Work.
   integer(I4P)                        :: i, k, it !< Work.

   allocate(s(n), w(n))
   do i = 1_I4P, n
      z = cos(PI * (i - 0.25_R8P) / (n + 0.5_R8P))
      do it = 1, 100
         p0 = 1._R8P ; p1 = 0._R8P
         do k = 1_I4P, n
            p2 = p1 ; p1 = p0
            p0 = ((2 * k - 1) * z * p1 - (k - 1) * p2) / k
         enddo
         dp = n * (z * p0 - p1) / (z * z - 1._R8P)
         z1 = z ; z = z1 - p0 / dp
         if (abs(z - z1) <= 1.e-15_R8P) exit
      enddo
      s(n + 1 - i) = 0.5_R8P * (1._R8P + z) ! z descends with i: s ascends with n + 1 - i
      w(n + 1 - i) = 1._R8P / ((1._R8P - z * z) * dp * dp) ! half the [-1, 1] weight 2 / ((1 - z^2) P_n'^2)
   enddo
   endsubroutine gauss_legendre_01
endmodule fossil_swept_cells
