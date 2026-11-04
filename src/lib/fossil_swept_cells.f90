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
!< - otherwise (2-D rotation, 3-D): Gauss-Legendre at nq nodes, the apertures from a copy moved to
!<   each node pose (2-D: exact line crossings of every face; 3-D: `cut_cells_3d` face fractions);
!<   with the screw motion, `error_estimate` is the largest face difference between nq and 2 nq nodes.
!<
!< Check (`residual`): the volume balance of the rigid velocity field of the screw motion, a field
!< that moves with the body (Falcovitz's design requirement: uniform flow moving with the body is
!< exact). Per cell, residual = swept + int_0^1 sum_f (flux of the field out of the fluid part of
!< face f) ds; the field is linear in x, so a face flux is its fluid length (area) times the field at
!< the fluid centroid. It is round-off for the exact 2-D translation and the quadrature error
!< otherwise. Not computed (0) when the caller gives the node poses.
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_cut_cells_3d,       only : cut_cells_3d_t, cut_cells_3d
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t, pose_interpolate, pose_angle_2d, pose_log_3d
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P

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
   real(R8P)                 :: error_estimate = 0._R8P !< Max face |A(nq) - A(2 nq)| (screw motion).
endtype swept_cells_3d_t

type :: slist_t
   !< A list of breakpoints.
   real(R8P), allocatable :: s(:) !< Breakpoints.
endtype slist_t

contains
   subroutine swept_cells_2d(curve, pose0, pose1, x, y, sw, nq, poses, exact, status)
   !< Swept geometry of `curve` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1) (issue #41 step 5; see the module notes). `nq` (default 4) Gauss-Legendre nodes for the
   !< rotating case; `poses(nq)` the poses at the nodes (`gauss_legendre_01`) instead of the screw
   !< motion; `exact` (default .true.) uses the exact integration for a pure translation.
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   type(curve_object),     intent(in)            :: curve    !< Curve (body frame).
   type(rigid_pose_2d_t),  intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_2d_t),  intent(in)            :: pose1    !< Pose at the end.
   real(R8P),              intent(in)            :: x(:), y(:) !< Cell edges (strictly increasing).
   type(swept_cells_2d_t), intent(out)           :: sw       !< Swept geometry.
   integer(I4P),           intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_2d_t),  intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                intent(in),  optional :: exact    !< Exact translation path (default .true.).
   integer(I4P),           intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   type(curve_object)       :: c0, c1                !< Copies at the two poses.
   type(cut_cells_2d_t)     :: k0, k1                !< Their cut cells.
   real(R8P), allocatable   :: fu2(:,:), fv2(:,:), fl2(:,:), res(:,:) !< 2 nq face lengths, work, residual.
   integer(I4P)             :: nx, ny, n, i, j       !< Work.
   real(R8P)                :: d(2), dth             !< Displacement, rotation angle.
   logical                  :: translation, screw    !< Paths.

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
   if (present(exact)) translation = translation .and. exact
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

   subroutine swept_cells_3d(surface, pose0, pose1, x, y, z, sw, nq, poses, status)
   !< Swept geometry of `surface` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1), z(nz+1) (issue #41 step 5): volumes from `cut_cells_3d` at the two poses, face apertures
   !< by `nq`-point Gauss-Legendre (default 4) on copies moved to the screw poses or to `poses(nq)`.
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   type(surface_stl_object), intent(in)            :: surface  !< Surface (body frame).
   type(rigid_pose_t),       intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_t),       intent(in)            :: pose1    !< Pose at the end.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(swept_cells_3d_t),   intent(out)           :: sw       !< Swept geometry.
   integer(I4P),             intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_t),       intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   integer(I4P),             intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   type(surface_stl_object) :: s0                 !< Moved copy.
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
   s0 = surface ; call s0%move(pose0) ; call cut_cells_3d(s0, x, y, z, k0, st)
   s0 = surface ; call s0%move(pose1) ; call cut_cells_3d(s0, x, y, z, k1, st)
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
   call faces_gauss_3d(surface, pose0, pose1, x, y, z, n, sw%face_x, sw%face_y, sw%face_z, flux, poses)
   if (.not. present(poses)) then
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
