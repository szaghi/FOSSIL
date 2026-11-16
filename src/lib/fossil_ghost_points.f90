!< FOSSIL, ghost points of a finite-difference wall closure on a rectilinear block (#109 WP1 step 03).
module fossil_ghost_points
!< FOSSIL, ghost points of a finite-difference wall closure on a rectilinear block (#109 WP1 step 03).
!<
!< The solver passes its own point lines x(:), y(:)[, z(:)] (cell centres or nodes: any strictly increasing
!< coordinates, at least two per axis) and the stencil half-width `reach(a) >= 1` per axis. A point is solid when it is
!< inside the body by the exact perturbed node states of the cut cells (#69, #70): the point (x_i, y_j, z_k) is taken at
!< (x_i + eps1, y_j + eps2, z_k + eps3), so a point exactly on the wall belongs to the side the perturbation puts it
!< on, the low-side rule of the cut cells (a cell is (x_i, x_i+1]: the wall on the plane x = x_i leaves the point x_i on
!< the side of the cell above it). 3-D: the parities of the z lines (`cut_cells_3d_nodes`, `cut_cells_3d_nodes_device`);
!< 2-D: the parities of the rows (`node_flags`, `node_rows` + `node_parity`). On lines that are also the cut-cell edges
!< the flags are the cut's `node_in` bit for bit.
!<
!< A ghost point is a solid point with a fluid point within `reach(a)` steps along a grid line of axis a, in either
!< direction, inside the block. Per ghost: its indexes (i, j[, k]), coordinates, a direction mask (bit 2(a-1) set when
!< the -a side sees fluid within reach, bit 2a-1 for +a: the GHOST_* constants) and, per direction, the steps to the
!< first fluid point (1..reach(a), 0 when none within reach). The ghosts are listed in lexicographic order, i fastest.
!< Solid points that are not ghosts are not listed.
!<
!< Two paths, the same output bit for bit (integers and copied coordinates): the host loops, and `device=.true.`, the
!< node states on the device (3-D: copied back by `cut_cells_3d_nodes_device`, then mapped again; 2-D: they stay
!< there), then one thread per point counts (`ghost_mask_oac`), an exclusive prefix sum gives each ghost its slot
!< (`links_scan_dev`) and one thread per point fills it (`ghost_fill_oac`). Built without OpenACC the device loops run
!< on the host in order.
use fossil_curve_object,       only : curve_object
use fossil_cut_cells_2d,       only : node_flags, node_rows, node_parity
use fossil_cut_cells_3d,       only : cut_cells_3d_nodes, cut_cells_3d_nodes_device
use fossil_cut_graph,          only : links_scan_dev
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P

implicit none
private
public :: ghost_points_t, ghost_points_3d, ghost_points_2d, ghost_points_reset
public :: GHOST_XM, GHOST_XP, GHOST_YM, GHOST_YP, GHOST_ZM, GHOST_ZP
! the per-point kernels, shared by the host and the device paths (for the closure table of step 04)
public :: ghost_side_oac, ghost_mask_oac, ghost_fill_oac
! the 2-D ghosts of flat vertex arrays and the input check, for the closure table (#109 WP1 step 04; not re-exported)
public :: ghost_points_2d_flat, ghost_lines_check
! the 3-D ghosts of given node states, for the posed closure table (#109 WP1 step 06; not re-exported)
public :: ghost_points_nodes_3d

integer(I4P), parameter :: GHOST_XM = 1_I4P  !< Mask bit 0: fluid within reach on the -x side.
integer(I4P), parameter :: GHOST_XP = 2_I4P  !< Mask bit 1: on the +x side.
integer(I4P), parameter :: GHOST_YM = 4_I4P  !< Mask bit 2: on the -y side.
integer(I4P), parameter :: GHOST_YP = 8_I4P  !< Mask bit 3: on the +y side.
integer(I4P), parameter :: GHOST_ZM = 16_I4P !< Mask bit 4: on the -z side.
integer(I4P), parameter :: GHOST_ZP = 32_I4P !< Mask bit 5: on the +z side.

type :: ghost_points_t
   !< The ghost points of a block (2-D: the (u, v) plane of the curve, dims = 2, nk = 1).
   integer(I4P)              :: n = 0_I4P          !< Ghost points.
   integer(I4P)              :: dims = 0_I4P       !< Dimensions (2 or 3).
   integer(I4P)              :: ni = 0_I4P         !< Points along x (u).
   integer(I4P)              :: nj = 0_I4P         !< Points along y (v).
   integer(I4P)              :: nk = 0_I4P         !< Points along z (1 in 2-D).
   integer(I4P)              :: reach(3) = 0_I4P   !< Stencil half-width per axis (0 for z in 2-D).
   integer(I4P), allocatable :: ijk(:,:)           !< Point indexes [dims, n], lexicographic order (i fastest).
   real(R8P),    allocatable :: points(:,:)        !< Point coordinates [dims, n] (x(i), y(j)[, z(k)]).
   integer(I4P), allocatable :: mask(:)            !< Direction mask [n] (GHOST_* bits).
   integer(I4P), allocatable :: steps(:,:)         !< Steps to the first fluid point [2 dims, n]: row 2a-1 the -a side,
                                                   !< row 2a the +a side, 0 when none within reach.
endtype ghost_points_t

contains
   subroutine ghost_points_3d(surface, x, y, z, reach, gp, status, device)
   !< The ghost points of the surface on the point lines x, y, z with the stencil half-widths reach(3) (see the module
   !< header). `device=.true.` (default .false.): the device path (the surface needs its SAH BVH: `cut_cells_3d_nodes_device`;
   !< mapped by `surface%enter_device` or copied in), the same output bit for bit. A surface without facets has no ghost.
   !< `status = STATUS_INVALID_INPUT` on fewer than two points on an axis, lines not strictly increasing, size(reach) /= 3
   !< or a reach below 1, or (device) a surface without a SAH BVH; `gp` is then reset (n = 0).
   type(surface_stl_object), intent(in), target    :: surface          !< Surface.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Point lines (strictly increasing).
   integer(I4P),             intent(in)            :: reach(:)         !< Stencil half-width per axis [3], >= 1.
   type(ghost_points_t),     intent(inout)         :: gp               !< Ghost points (reset first).
   integer(I4P),             intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                  intent(in),  optional :: device           !< Device path (default .false.).
   logical,      allocatable :: nd(:,:,:)                              !< Point inside flags.
   real(R8P),    allocatable :: ex(:,:,:), ey(:,:,:), ez(:,:,:)        !< Edge fractions (not asked).
   integer(I4P)              :: st                                     !< Status.
   logical                   :: dev                                    !< Device path.

   call ghost_points_reset(gp)
   if (present(status)) status = STATUS_OK
   dev = .false. ; if (present(device)) dev = device
   st = lines_check(x, y, z, reach, 3_I4P)
   if (st == STATUS_OK) then
      if (surface%get_facets_number() < 1_I4P) then
         allocate(nd(size(x), size(y), size(z))) ; nd = .false.
      elseif (dev) then
         call cut_cells_3d_nodes_device(surface, x, y, z, nd, ex, ey, ez, status=st)
      else
         call cut_cells_3d_nodes(surface, x, y, z, nd, status=st)
      endif
   endif
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   gp%dims = 3_I4P ; gp%ni = size(x, kind=I4P) ; gp%nj = size(y, kind=I4P) ; gp%nk = size(z, kind=I4P)
   gp%reach = reach(1:3)
   if (dev) then
      !$acc data copyin(nd, x, y, z)
      call ghost_core_device(3_I4P, gp%ni, gp%nj, gp%nk, nd, x, y, z, gp%reach, gp)
      !$acc end data
   else
      call ghost_core_host(3_I4P, gp%ni, gp%nj, gp%nk, nd, x, y, z, gp%reach, gp)
   endif
   endsubroutine ghost_points_3d

   subroutine ghost_points_2d(curve, x, y, reach, gp, status, device)
   !< The ghost points of the curve on the point lines x (u), y (v) with the stencil half-widths reach(2) (see the module
   !< header; coordinates in the curve's (u, v) plane). `device=.true.` (default .false.): the node flags (`node_parity`
   !< per point over the rows of `node_rows`, built on the host), the classification and the fill on the device, the same
   !< output bit for bit. A curve without vertices has no ghost. `status = STATUS_INVALID_INPUT` on fewer than two points
   !< on an axis, lines not strictly increasing, size(reach) /= 2 or a reach below 1; `gp` is then reset (n = 0).
   type(curve_object),   intent(in)            :: curve      !< Curve.
   real(R8P),            intent(in)            :: x(:), y(:) !< Point lines (u, v; strictly increasing).
   integer(I4P),         intent(in)            :: reach(:)   !< Stencil half-width per axis [2], >= 1.
   type(ghost_points_t), intent(inout)         :: gp         !< Ghost points (reset first).
   integer(I4P),         intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,              intent(in),  optional :: device     !< Device path (default .false.).
   real(R8P),    allocatable :: u(:), v(:)                   !< Curve vertices.
   integer(I4P), allocatable :: nxt(:)                       !< Next vertex (0: none).
   real(R8P)                 :: uv(2)                        !< A vertex.
   integer(I4P)              :: nv, i                        !< Work.

   nv = curve%get_vertices_number()
   allocate(u(nv), v(nv), nxt(nv))
   do i = 1, nv
      uv = curve%get_vertex(i) ; u(i) = uv(1) ; v(i) = uv(2) ; nxt(i) = curve%get_next(i)
   enddo
   call ghost_points_2d_flat(u, v, nxt, x, y, reach, gp, status, device)
   endsubroutine ghost_points_2d

   subroutine ghost_points_2d_flat(u, v, nxt, x, y, reach, gp, status, device)
   !< `ghost_points_2d` on the vertices (u, v) and next-vertex links nxt of a curve (0: end of an open loop) instead of
   !< the curve object (#109 WP1 step 04: the closure table passes the vertices moved by a pose, as `curve%move` places
   !< them). Same output and statuses.
   real(R8P),            intent(in)            :: u(:), v(:) !< Vertices.
   integer(I4P),         intent(in)            :: nxt(:)     !< Next vertex (0: none).
   real(R8P),            intent(in)            :: x(:), y(:) !< Point lines (u, v; strictly increasing).
   integer(I4P),         intent(in)            :: reach(:)   !< Stencil half-width per axis [2], >= 1.
   type(ghost_points_t), intent(inout)         :: gp         !< Ghost points (reset first).
   integer(I4P),         intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,              intent(in),  optional :: device     !< Device path (default .false.).
   logical,      allocatable :: nd(:,:)                      !< Point inside flags.
   integer(I4P), allocatable :: rstart(:), rseg(:)           !< Segments crossing each row (device path).
   real(R8P)                 :: z(1)                         !< No z line.
   integer(I4P)              :: nv, i, j, nu, nr, ni, nj     !< Work.
   logical                   :: dev                          !< Device path.

   call ghost_points_reset(gp)
   if (present(status)) status = STATUS_OK
   dev = .false. ; if (present(device)) dev = device
   z = 0._R8P
   if (lines_check(x, y, z, reach, 2_I4P) /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   nv = size(u, kind=I4P)
   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P)
   gp%dims = 2_I4P ; gp%ni = ni ; gp%nj = nj ; gp%nk = 1_I4P
   gp%reach(1:2) = reach(1:2) ; gp%reach(3) = 0_I4P
   if (dev) then
      call node_rows(u, v, nxt, y, rstart, rseg)
      nu = nv ; nr = size(rseg, kind=I4P)
      allocate(nd(ni, nj))
      !$acc data copyin(u, v, nxt, rstart, rseg, x, y, z) create(nd)
      !$acc parallel loop gang vector collapse(2) default(present) present(u, v, nxt, rstart, rseg, x, y, nd)
      do j = 1_I4P, nj
         do i = 1_I4P, ni
            nd(i, j) = node_parity(x(i), y(j), rstart(j), rstart(j + 1) - 1_I4P, nr, rseg, nu, u, v, nxt)
         enddo
      enddo
      call ghost_core_device(2_I4P, ni, nj, 1_I4P, nd, x, y, z, gp%reach, gp)
      !$acc end data
   else
      call node_flags(u, v, nxt, x, y, nd)
      call ghost_core_host(2_I4P, ni, nj, 1_I4P, nd, x, y, z, gp%reach, gp)
   endif
   endsubroutine ghost_points_2d_flat

   subroutine ghost_points_nodes_3d(node_in, x, y, z, reach, gp)
   !< The 3-D ghost points of given point states node_in [size(x), size(y), size(z)] on the host (#109 WP1 step 06: the
   !< posed closure table computes the states on the refitted tree): what `ghost_points_3d` does after its node states,
   !< the same output. The lines and reach are not checked (`ghost_lines_check` first).
   logical,              intent(in)    :: node_in(:,:,:)   !< Point inside flags.
   real(R8P),            intent(in)    :: x(:), y(:), z(:) !< Point lines.
   integer(I4P),         intent(in)    :: reach(:)         !< Stencil half-width per axis [3].
   type(ghost_points_t), intent(inout) :: gp               !< Ghost points (reset first).

   call ghost_points_reset(gp)
   gp%dims = 3_I4P ; gp%ni = size(x, kind=I4P) ; gp%nj = size(y, kind=I4P) ; gp%nk = size(z, kind=I4P)
   gp%reach = reach(1:3)
   call ghost_core_host(3_I4P, gp%ni, gp%nj, gp%nk, node_in, x, y, z, gp%reach, gp)
   endsubroutine ghost_points_nodes_3d

   integer(I4P) function ghost_lines_check(x, y, z, reach, nd) result(st)
   !< The input check of the ghost points (#109 WP1 step 04): STATUS_OK if every line has at least two strictly
   !< increasing points and reach has nd entries >= 1 (z not checked in 2-D), else STATUS_INVALID_INPUT.
   real(R8P),    intent(in) :: x(:), y(:), z(:) !< Point lines.
   integer(I4P), intent(in) :: reach(:)         !< Stencil half-widths.
   integer(I4P), intent(in) :: nd               !< Dimensions.

   st = lines_check(x, y, z, reach, nd)
   endfunction ghost_lines_check

   subroutine ghost_points_reset(gp)
   !< Reset the ghost points (no ghost, arrays released).
   type(ghost_points_t), intent(inout) :: gp    !< Ghost points.
   type(ghost_points_t)                :: fresh !< Default state.

   gp = fresh
   endsubroutine ghost_points_reset

   ! private

   integer(I4P) function lines_check(x, y, z, reach, nd) result(st)
   !< STATUS_OK if every line has at least two strictly increasing points and reach has nd entries >= 1 (z not checked
   !< in 2-D), else STATUS_INVALID_INPUT.
   real(R8P),    intent(in) :: x(:), y(:), z(:) !< Point lines.
   integer(I4P), intent(in) :: reach(:)         !< Stencil half-widths.
   integer(I4P), intent(in) :: nd               !< Dimensions.

   st = STATUS_INVALID_INPUT
   if (size(reach) /= nd) return
   if (any(reach < 1_I4P)) return
   if (size(x) < 2 .or. size(y) < 2) return
   if (any(x(2:) <= x(:size(x) - 1)) .or. any(y(2:) <= y(:size(y) - 1))) return
   if (nd == 3_I4P) then
      if (size(z) < 2) return
      if (any(z(2:) <= z(:size(z) - 1))) return
   endif
   st = STATUS_OK
   endfunction lines_check

   subroutine ghost_core_host(nd, ni, nj, nk, node_in, x, y, z, r, gp)
   !< The ghost points from the point flags, on the host: a count pass, then a fill pass in lexicographic order, through
   !< the kernels of the device path.
   integer(I4P),         intent(in)    :: nd                  !< Dimensions.
   integer(I4P),         intent(in)    :: ni, nj, nk          !< Points per axis.
   logical,              intent(in)    :: node_in(ni, nj, nk) !< Point inside flags.
   real(R8P),            intent(in)    :: x(:), y(:), z(:)    !< Point lines.
   integer(I4P),         intent(in)    :: r(3)                !< Reach per axis.
   type(ghost_points_t), intent(inout) :: gp                  !< Ghost points (n and arrays written).
   integer(I4P) :: i, j, k, n, g                              !< Work.

   n = 0_I4P
   do k = 1_I4P, nk
      do j = 1_I4P, nj
         do i = 1_I4P, ni
            if (ghost_mask_oac(ni, nj, nk, node_in, i, j, k, nd, r(1), r(2), r(3)) /= 0_I4P) n = n + 1_I4P
         enddo
      enddo
   enddo
   allocate(gp%ijk(nd, n), gp%points(nd, n), gp%mask(n), gp%steps(2 * nd, n))
   gp%n = n
   g = 0_I4P
   do k = 1_I4P, nk
      do j = 1_I4P, nj
         do i = 1_I4P, ni
            if (ghost_mask_oac(ni, nj, nk, node_in, i, j, k, nd, r(1), r(2), r(3)) == 0_I4P) cycle
            g = g + 1_I4P
            call ghost_fill_oac(ni, nj, nk, node_in, i, j, k, nd, r(1), r(2), r(3), gp%mask(g), gp%steps(1, g))
            gp%ijk(1, g) = i ; gp%ijk(2, g) = j ; gp%points(1, g) = x(i) ; gp%points(2, g) = y(j)
            if (nd == 3_I4P) then
               gp%ijk(3, g) = k ; gp%points(3, g) = z(k)
            endif
         enddo
      enddo
   enddo
   endsubroutine ghost_core_host

   subroutine ghost_core_device(nd, ni, nj, nk, node_in, x, y, z, r, gp)
   !< The ghost points from the point flags, on the device: one thread per point counts its ghost (0 or 1), the exclusive
   !< prefix sum of the counts (`links_scan_dev`) gives each ghost its slot in lexicographic order, then one thread per
   !< point fills its slot; the outputs are copied back into gp. node_in, x, y, z must be present.
   integer(I4P),         intent(in)    :: nd                  !< Dimensions.
   integer(I4P),         intent(in)    :: ni, nj, nk          !< Points per axis.
   logical,              intent(in)    :: node_in(ni, nj, nk) !< Point inside flags (present).
   real(R8P),            intent(in)    :: x(:), y(:), z(:)    !< Point lines (present).
   integer(I4P),         intent(in)    :: r(3)                !< Reach per axis.
   type(ghost_points_t), intent(inout) :: gp                  !< Ghost points (n and arrays written).
   integer(I4P), allocatable :: cnt(:), off(:)                !< Ghost count and slot offset per point.
   integer(I4P), allocatable :: ijk(:,:), mask(:), steps(:,:) !< Outputs (built here).
   real(R8P),    allocatable :: pts(:,:)                      !< Outputs (built here).
   integer(I4P) :: i, j, k, p, g, n, np, r1, r2, r3           !< Work.

   np = ni * nj * nk ; r1 = r(1) ; r2 = r(2) ; r3 = r(3)
   allocate(cnt(np), off(np))
   !$acc data create(cnt, off)
   !$acc parallel loop gang vector collapse(3) default(present) present(node_in, cnt) private(p)
   do k = 1_I4P, nk
      do j = 1_I4P, nj
         do i = 1_I4P, ni
            p = i + ni * ((j - 1_I4P) + nj * (k - 1_I4P))
            cnt(p) = 0_I4P
            if (ghost_mask_oac(ni, nj, nk, node_in, i, j, k, nd, r1, r2, r3) /= 0_I4P) cnt(p) = 1_I4P
         enddo
      enddo
   enddo
   call links_scan_dev(np, cnt, off, n)
   allocate(ijk(nd, n), pts(nd, n), mask(n), steps(2 * nd, n))
   if (n > 0_I4P) then
      !$acc data copyout(ijk, pts, mask, steps)
      !$acc parallel loop gang vector collapse(3) default(present) &
      !$acc& present(node_in, x, y, z, cnt, off, ijk, pts, mask, steps) private(p, g)
      do k = 1_I4P, nk
         do j = 1_I4P, nj
            do i = 1_I4P, ni
               p = i + ni * ((j - 1_I4P) + nj * (k - 1_I4P))
               if (cnt(p) == 1_I4P) then
                  g = off(p) + 1_I4P
                  call ghost_fill_oac(ni, nj, nk, node_in, i, j, k, nd, r1, r2, r3, mask(g), steps(1, g))
                  ijk(1, g) = i ; ijk(2, g) = j ; pts(1, g) = x(i) ; pts(2, g) = y(j)
                  if (nd == 3_I4P) then
                     ijk(3, g) = k ; pts(3, g) = z(k)
                  endif
               endif
            enddo
         enddo
      enddo
      !$acc end data
   endif
   !$acc end data
   gp%n = n
   call move_alloc(ijk, gp%ijk) ; call move_alloc(pts, gp%points) ; call move_alloc(mask, gp%mask)
   call move_alloc(steps, gp%steps)
   endsubroutine ghost_core_device

   pure integer(I4P) function ghost_side_oac(ni, nj, nk, node_in, i, j, k, a, sgn, r) result(s)
   !< Kernel: steps from the point (i, j, k) to the first fluid point along axis a in the direction sgn (-1 or +1),
   !< 1..r, or 0 if the r points that way are solid or the block ends before a fluid one.
   !$acc routine seq
   integer(I4P), intent(in) :: ni, nj, nk          !< Points per axis.
   logical,      intent(in) :: node_in(ni, nj, nk) !< Point inside flags.
   integer(I4P), intent(in) :: i, j, k             !< Point.
   integer(I4P), intent(in) :: a, sgn, r           !< Axis, direction, reach.
   integer(I4P)             :: t, ii, jj, kk       !< Work.

   s = 0_I4P
   do t = 1_I4P, r
      ii = i ; jj = j ; kk = k
      if (a == 1_I4P) then
         ii = i + sgn * t
      elseif (a == 2_I4P) then
         jj = j + sgn * t
      else
         kk = k + sgn * t
      endif
      if (ii < 1_I4P .or. ii > ni .or. jj < 1_I4P .or. jj > nj .or. kk < 1_I4P .or. kk > nk) return
      if (.not. node_in(ii, jj, kk)) then
         s = t ; return
      endif
   enddo
   endfunction ghost_side_oac

   pure integer(I4P) function ghost_mask_oac(ni, nj, nk, node_in, i, j, k, nd, r1, r2, r3) result(m)
   !< Kernel: the direction mask of the point (i, j, k), 0 if it is fluid or no fluid point lies within reach (not a
   !< ghost).
   !$acc routine seq
   integer(I4P), intent(in) :: ni, nj, nk          !< Points per axis.
   logical,      intent(in) :: node_in(ni, nj, nk) !< Point inside flags.
   integer(I4P), intent(in) :: i, j, k             !< Point.
   integer(I4P), intent(in) :: nd                  !< Dimensions.
   integer(I4P), intent(in) :: r1, r2, r3          !< Reach per axis.
   integer(I4P)             :: a, r                !< Work.

   m = 0_I4P
   if (.not. node_in(i, j, k)) return
   do a = 1_I4P, nd
      r = r1 ; if (a == 2_I4P) r = r2 ; if (a == 3_I4P) r = r3
      if (ghost_side_oac(ni, nj, nk, node_in, i, j, k, a, -1_I4P, r) > 0_I4P) m = ibset(m, 2_I4P * (a - 1_I4P))
      if (ghost_side_oac(ni, nj, nk, node_in, i, j, k, a, 1_I4P, r) > 0_I4P) m = ibset(m, 2_I4P * a - 1_I4P)
   enddo
   endfunction ghost_mask_oac

   pure subroutine ghost_fill_oac(ni, nj, nk, node_in, i, j, k, nd, r1, r2, r3, mask, steps)
   !< Kernel: the direction mask and the steps to the first fluid point per direction of the point (i, j, k).
   !$acc routine seq
   integer(I4P), intent(in)  :: ni, nj, nk          !< Points per axis.
   logical,      intent(in)  :: node_in(ni, nj, nk) !< Point inside flags.
   integer(I4P), intent(in)  :: i, j, k             !< Point.
   integer(I4P), intent(in)  :: nd                  !< Dimensions.
   integer(I4P), intent(in)  :: r1, r2, r3          !< Reach per axis.
   integer(I4P), intent(out) :: mask                !< Direction mask.
   integer(I4P), intent(out) :: steps(2 * nd)       !< Steps per direction (-x, +x, -y, +y[, -z, +z]).
   integer(I4P)              :: a, r, sm, sp        !< Work.

   mask = 0_I4P
   do a = 1_I4P, nd
      r = r1 ; if (a == 2_I4P) r = r2 ; if (a == 3_I4P) r = r3
      sm = 0_I4P ; sp = 0_I4P
      if (node_in(i, j, k)) then
         sm = ghost_side_oac(ni, nj, nk, node_in, i, j, k, a, -1_I4P, r)
         sp = ghost_side_oac(ni, nj, nk, node_in, i, j, k, a, 1_I4P, r)
      endif
      steps(2_I4P * a - 1_I4P) = sm ; steps(2_I4P * a) = sp
      if (sm > 0_I4P) mask = ibset(mask, 2_I4P * (a - 1_I4P))
      if (sp > 0_I4P) mask = ibset(mask, 2_I4P * a - 1_I4P)
   enddo
   endsubroutine ghost_fill_oac
endmodule fossil_ghost_points
