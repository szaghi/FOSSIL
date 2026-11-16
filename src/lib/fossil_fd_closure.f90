!< FOSSIL, the closure table of a finite-difference wall closure on a rectilinear block (#109 WP1 step 04).
module fossil_fd_closure
!< FOSSIL, the closure table of a finite-difference wall closure on a rectilinear block (#109 WP1 step 04).
!<
!< One call gives, per ghost point of the solver's own point lines (`fossil_ghost_points`), every geometric ingredient
!< of a sharp-interface wall closure; the solver picks the closure (image line, two samples per ray, inverse
!< Lax-Wendroff, curvature condition, eikonal extension) and interpolates its own flow:
!<
!< - the ghost: indexes, coordinates, direction mask, steps to the first fluid point per direction (as `ghost_points_t`);
!< - the wall: signed distance (negative inside), foot point (closest point), pseudo-normal of the closest feature
!<   (outward, into the fluid), element (3-D facet, 2-D segment = its first vertex), feature (3-D FEATURE_FACE / _EDGE /
!<   _VERTEX; 2-D the closest vertex, 0 in a segment interior, as `query_many`'s `vertices`), body, patch;
!< - m probes along the normal: probe(:, k, g) = foot + lengths(k) normal, the same signed lengths for every ghost, in
!<   absolute units (+ outward; an image line of spacing l is lengths = [l, 2l, ...]); visible(k, g) when the open segment
!<   foot-probe crosses no element (the kernel of `probe_many`); probe_cell(a, k, g) the enclosing cell on the solver's
!<   lines, low-side rule: x(i) < p <= x(i+1) gives i, outside the block (or on its first line) 0;
!< - the wall curvature at the foot when `build_curvature` ran: 3-D kappa = (k1, k2), k1 >= k2, direction = (d1, d2)
!<   (`query_many(principal_curvatures=, principal_directions=)`); 2-D kappa(1) the curvature, direction(:, 1) the
!<   tangent (-n_v, n_u);
!< - the wall velocity at the foot of a rigid body when `linear_velocity` or `angular_velocity` is given:
!<   v = V + omega x (foot - c), c = `centre` (default: the pose's translation, the body origin in the world; the world
!<   origin without a pose). The rigid field is exact for any rigid motion at the instant of the pose. 2-D: omega is a
!<   scalar, omega x r = omega (-r_v, r_u);
!< - with `ilw=.true.`, the inverse Lax-Wendroff data per ghost and direction slot (slot 2a-1 the -a side, 2a the +a
!<   side, as `steps`): when the slot's steps s > 0, the crossing of that grid line with the wall between the last solid
!<   point (s - 1 steps away) and the first fluid point (s steps): ilw_t its coordinate along the axis, ilw_element the
!<   element, ilw_normal the pseudo-normal of the crossing point (a query of that point), ilw_kappa / ilw_direction the
!<   curvature there if built. Of the exact crossings between the two points (an odd number on a watertight wall, the
!<   node states being the parities of the same crossings), the one nearest the fluid point (ties by element index, the
!<   order of `line_crossings`). The crossings are the cut cells' (exact, `line_crossing_before` against the perturbed
!<   points: the low-side rule), so on the cut-cell lines the crossing lies on the cell edge the cut cells cut. A slot
!<   with no crossing (steps 0, or an open wall) has element 0 and zeros.
!<
!< Poses (`pose=`): the body is at the pose, everything is in the world frame. The ghosts and the ILW crossings are those
!< of the surface moved by the pose, the geometry the cut cells of the moved body see, without moving it (3-D, #109 WP1
!< step 06: the facet vertices placed by `to_world` as `move` places them, and the SAH BVH of the given surface refitted
!< to them as `move` refits it, both built for the call, on the device on the device path; 2-D: the vertices placed by
!< `to_world`, as `curve%move` places them); every query (foot, normal, probes, visibility, curvature, the ILW normals)
!< runs in the body frame on the given geometry and is mapped back, as `query_many(pose=)`. An octree surface (host path
!< only) takes a moved copy, `move` (the octree is rebuilt there).
!<
!< Two paths, the same table bit for bit. Host: `ghost_points_3d/_2d`, `query_many`, `probe_many`, `line_crossings` and
!< the ILW selection kernel. Device (`device=.true.`; the surface or curve mapped by `enter_device`, and its curvature if
!< built): every stage on the device in one pipeline (node states, ghost count, scan and fill, body-frame points,
!< queries, wall data, probes and visibility, ILW crossings, their queries), the same kernels; no round trip between the
!< stages. `resident=.true.` (implies the device path) leaves the table mapped: the resident pattern of the cut cells
!< (`on_device`, `device_key`, `enter_device` / `update_host` / `exit_device`; `ib_closure_reset`; every writer takes the
!< table intent(inout) and releases it first). Built without OpenACC the device loops run on the host in order.
use fossil_aabb_node_object,   only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, bvh_bmax_x_oac, &
                                      bvh_bmax_y_oac, bvh_bmax_z_oac, bvh_left_child_oac, bvh_right_child_oac, &
                                      bvh_payload_first_oac, bvh_payload_count_oac, bvh_node_refit_oac
use fossil_aabb_tree_object,   only : facet_distance_payload, AABB_TREE_SAH_BVH
use fossil_curve_object,       only : curve_object
use fossil_cut_cells_2d,       only : node_rows, node_parity
use fossil_cut_cells_3d,       only : cut_cells_3d_nodes_tree, cut_cells_3d_nodes_tree_present
use fossil_cut_graph,          only : links_scan_dev
use fossil_ghost_points,       only : ghost_points_t, ghost_points_3d, ghost_points_2d_flat, ghost_lines_check, &
                                      ghost_mask_oac, ghost_fill_oac, ghost_points_nodes_3d
use fossil_line_query,         only : bvh_line_count_oac, bvh_line_crossings_oac, line_crossing_before, &
                                      line_crossing_before_2d, seg_line_hit_oac
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t
use fossil_segment_visibility, only : bvh_segment_blocked_oac
use fossil_shape_operator,     only : foot_principal_oac
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use iso_c_binding,             only : c_intptr_t, c_loc
use iso_fortran_env,           only : stderr => error_unit, int64
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
private
public :: ib_closure_t, fd_closure_3d, fd_closure_2d, ib_closure_reset, ib_closure_mismatches
! the per-item kernels, shared by the two paths (public for the tests)
public :: fd_cell_index_oac, ilw_select_oac, ilw_item_2d_oac

type :: ib_closure_t
   !< The closure table of a block (see the module header). dims = 3 or 2 (2-D: the curve's (u, v) plane, nk = 1). Every
   !< array is allocated after a call (zero-size when not asked: kappa and direction without curvature, velocity without
   !< velocities, the ilw_* without `ilw`, the probes with m = 0).
   integer(I4P)              :: n = 0_I4P              !< Ghost points.
   integer(I4P)              :: dims = 0_I4P           !< Dimensions (2 or 3).
   integer(I4P)              :: m = 0_I4P              !< Probes per ghost.
   integer(I4P)              :: ni = 0_I4P             !< Points along x (u).
   integer(I4P)              :: nj = 0_I4P             !< Points along y (v).
   integer(I4P)              :: nk = 0_I4P             !< Points along z (1 in 2-D).
   integer(I4P)              :: reach(3) = 0_I4P       !< Stencil half-width per axis (0 for z in 2-D).
   logical                   :: with_curvature = .false. !< kappa and direction filled (build_curvature ran).
   logical                   :: with_velocity = .false.  !< velocity filled.
   logical                   :: with_ilw = .false.       !< ilw_* filled.
   integer(I4P), allocatable :: ijk(:,:)               !< Ghost indexes [dims, n], lexicographic order (i fastest).
   real(R8P),    allocatable :: points(:,:)            !< Ghost coordinates [dims, n].
   integer(I4P), allocatable :: mask(:)                !< Direction mask [n] (GHOST_* bits).
   integer(I4P), allocatable :: steps(:,:)             !< Steps to the first fluid point [2 dims, n] (0: none in reach).
   real(R8P),    allocatable :: distance(:)            !< Signed distance [n] (negative inside).
   real(R8P),    allocatable :: foot(:,:)              !< Foot point [dims, n].
   real(R8P),    allocatable :: normal(:,:)            !< Pseudo-normal at the foot [dims, n] (outward).
   integer(I4P), allocatable :: element(:)             !< Closest facet / segment [n].
   integer(I4P), allocatable :: feature(:)             !< 3-D FEATURE_*; 2-D closest vertex (0: segment interior) [n].
   integer(I4P), allocatable :: body(:)                !< Body [n].
   integer(I4P), allocatable :: patch(:)               !< Patch [n].
   real(R8P),    allocatable :: lengths(:)             !< Probe lengths [m] (+ outward, absolute).
   real(R8P),    allocatable :: probe(:,:,:)           !< Probe points [dims, m, n].
   logical,      allocatable :: visible(:,:)           !< Segment foot-probe crosses no element [m, n].
   integer(I4P), allocatable :: probe_cell(:,:,:)      !< Enclosing cell per axis [dims, m, n] (0: outside).
   real(R8P),    allocatable :: kappa(:,:)             !< Principal curvatures [dims-1, n] (2-D: the curvature).
   real(R8P),    allocatable :: direction(:,:,:)       !< Principal directions [dims, dims-1, n] (2-D: the tangent).
   real(R8P),    allocatable :: velocity(:,:)          !< Wall velocity at the foot [dims, n].
   real(R8P),    allocatable :: ilw_t(:,:)             !< ILW crossing coordinate along the slot's axis [2 dims, n].
   integer(I4P), allocatable :: ilw_element(:,:)       !< ILW crossing element [2 dims, n] (0: none).
   real(R8P),    allocatable :: ilw_normal(:,:,:)      !< Pseudo-normal at the crossing [dims, 2 dims, n].
   real(R8P),    allocatable :: ilw_kappa(:,:,:)       !< Curvature at the crossing [dims-1, 2 dims, n].
   real(R8P),    allocatable :: ilw_direction(:,:,:,:) !< Directions at the crossing [dims, dims-1, 2 dims, n].
   logical                   :: on_device = .false.    !< The arrays are mapped on the device (resident).
   integer(c_intptr_t)       :: device_key = 0         !< Host address of points when mapped (a copy has another one).
   contains
      procedure, pass(self) :: enter_device !< Map the table to the device (copyin).
      procedure, pass(self) :: update_host  !< Copy the mapped table back to the host.
      procedure, pass(self) :: exit_device  !< Release the device copies.
endtype ib_closure_t

integer(I4P), parameter :: ILW_BUF = 16_I4P !< Private crossing buffer of an ILW slot (more: the exact-size path).

interface mis_i
   module procedure mis_i1, mis_i2, mis_i3
endinterface mis_i
interface mis_r
   module procedure mis_r1, mis_r2, mis_r3, mis_r4
endinterface mis_r

contains
   subroutine fd_closure_3d(surface, x, y, z, reach, closure, status, lengths, pose, linear_velocity, angular_velocity, &
                            centre, ilw, device, resident)
   !< The closure table of the surface on the point lines x, y, z with the stencil half-widths reach(3) (see the module
   !< header). `lengths(m)`: probe lengths (absolute, + outward; default none). `pose`: the body at the pose.
   !< `linear_velocity(3)`, `angular_velocity(3)`, `centre(3)`: the rigid wall velocity. `ilw` (default .false.): the
   !< inverse Lax-Wendroff crossings. `device` (default .false.): the device path (`enter_device` first, and after
   !< `build_curvature` if built); `resident` (default .false., implies device): leave the table mapped.
   !< `status = STATUS_INVALID_INPUT` on bad lines or reach (`ghost_points_3d`), a bad size of the velocities or the
   !< centre, or (device) a surface with facets not mapped, or its curvature built and not mapped (a stderr message);
   !< the table is then reset (n = 0). A surface without facets gives no ghost.
   type(surface_stl_object), intent(in), target    :: surface             !< Surface.
   real(R8P),                intent(in)            :: x(:), y(:), z(:)    !< Point lines (strictly increasing).
   integer(I4P),             intent(in)            :: reach(:)            !< Stencil half-width per axis [3], >= 1.
   type(ib_closure_t),       intent(inout)         :: closure             !< Table (released and reset first).
   integer(I4P),             intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                intent(in),  optional :: lengths(:)          !< Probe lengths [m] (+ outward, absolute).
   type(rigid_pose_t),       intent(in),  optional :: pose                !< Pose of the body.
   real(R8P),                intent(in),  optional :: linear_velocity(:)  !< Velocity of the centre [3].
   real(R8P),                intent(in),  optional :: angular_velocity(:) !< Angular velocity [3].
   real(R8P),                intent(in),  optional :: centre(:)           !< Centre of the rotation [3].
   logical,                  intent(in),  optional :: ilw                 !< ILW crossings (default .false.).
   logical,                  intent(in),  optional :: device              !< Device path (default .false.).
   logical,                  intent(in),  optional :: resident            !< Leave the table mapped (default .false.).
   type(surface_stl_object), target  :: sm                    !< The surface moved by the pose (octree only).
   type(surface_stl_object), pointer :: geo                   !< Geometry of the ghosts and crossings.
   real(R8P)                         :: lv(3), av(3), c(3)    !< Rigid velocity.
   real(R8P)                         :: pr(3, 3), pt(3)       !< Pose.
   real(R8P), allocatable            :: ll(:)                 !< Probe lengths.
   integer(I4P)                      :: st, m                 !< Status, probes.
   logical                           :: dev, res, wv, wi, wc, posed, tree !< Options; posed SAH BVH refitted.

   call ib_closure_reset(closure)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   wi = .false. ; if (present(ilw)) wi = ilw
   posed = present(pose)
   st = ghost_lines_check(x, y, z, reach, 3_I4P)
   call velocity_inputs(3_I4P, st, linear_velocity, angular_velocity, centre, lv, av, c, wv)
   if (present(lengths)) then
      allocate(ll(size(lengths))) ; ll = lengths
   else
      allocate(ll(0))
   endif
   m = size(ll, kind=I4P)
   pr = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ; pt = 0._R8P
   if (posed) then
      pr = pose%r ; pt = pose%t
      if (.not. present(centre)) c = pose%t
   endif
   wc = surface%is_curvature_built()
   if (st == STATUS_OK .and. dev .and. surface%get_facets_number() > 0_I4P) then
      if (.not. surface%is_on_device) then
         write(stderr, '(A)') 'fossil fd_closure_3d: device path refused -- call surface%enter_device first'
         st = STATUS_INVALID_INPUT
      elseif (wc .and. .not. surface%is_curvature_on_device()) then
         write(stderr, '(A)') 'fossil fd_closure_3d: device path refused -- build_curvature, then enter_device'
         st = STATUS_INVALID_INPUT
      endif
   endif
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   geo => surface ; tree = .false.
   if (posed .and. surface%get_facets_number() > 0_I4P) then
      tree = surface%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(surface%aabb%node) .and. &
             allocated(surface%aabb%payload) .and. allocated(surface%tri_vertices)
      if (.not. tree) then ! an octree: a moved copy (move rebuilds it)
         sm = surface ; sm%is_on_device = .false. ; call sm%move(pose) ; geo => sm
      endif
   endif
   if (dev .and. surface%get_facets_number() > 0_I4P) then
      call closure_3d_device(surface, geo, tree, x, y, z, reach, ll, posed, pr, pt, wc, wv, lv, av, c, wi, res, closure)
   else
      call closure_3d_host(surface, geo, tree, x, y, z, reach, ll, pose, pr, pt, wc, wv, lv, av, c, wi, closure)
   endif
   endsubroutine fd_closure_3d

   subroutine fd_closure_2d(curve, x, y, reach, closure, status, lengths, pose, linear_velocity, angular_velocity, centre, &
                            ilw, device, resident)
   !< The closure table of the curve on the point lines x (u), y (v) with the stencil half-widths reach(2), in the curve's
   !< (u, v) plane (see the module header). `angular_velocity` is a scalar (rad per unit time, counter-clockwise in
   !< (u, v)), `linear_velocity(2)`, `centre(2)`. The other arguments and the statuses as `fd_closure_3d` (device: the
   !< curve mapped by `enter_device`, after `build_curvature` if built).
   type(curve_object),    intent(in), target    :: curve              !< Curve.
   real(R8P),             intent(in)            :: x(:), y(:)         !< Point lines (u, v; strictly increasing).
   integer(I4P),          intent(in)            :: reach(:)           !< Stencil half-width per axis [2], >= 1.
   type(ib_closure_t),    intent(inout)         :: closure            !< Table (released and reset first).
   integer(I4P),          intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),             intent(in),  optional :: lengths(:)         !< Probe lengths [m] (+ outward, absolute).
   type(rigid_pose_2d_t), intent(in),  optional :: pose               !< Pose of the curve.
   real(R8P),             intent(in),  optional :: linear_velocity(:) !< Velocity of the centre [2].
   real(R8P),             intent(in),  optional :: angular_velocity   !< Angular velocity (scalar).
   real(R8P),             intent(in),  optional :: centre(:)          !< Centre of the rotation [2].
   logical,               intent(in),  optional :: ilw                !< ILW crossings (default .false.).
   logical,               intent(in),  optional :: device             !< Device path (default .false.).
   logical,               intent(in),  optional :: resident           !< Leave the table mapped (default .false.).
   real(R8P),    allocatable :: u(:), v(:), ll(:)                     !< Vertices (placed), probe lengths.
   integer(I4P), allocatable :: nxt(:)                                !< Next vertex.
   real(R8P)                 :: lv(3), av(3), c(3), uv(2), z(1)       !< Rigid velocity, a vertex, no z line.
   real(R8P)                 :: pr(2, 2), pt(2), w1, w2               !< Pose, placed vertex.
   integer(I4P)              :: st, m, nv, i                          !< Work.
   logical                   :: dev, res, wv, wi, wc, posed           !< Options.

   call ib_closure_reset(closure)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   wi = .false. ; if (present(ilw)) wi = ilw
   posed = present(pose)
   z = 0._R8P
   st = ghost_lines_check(x, y, z, reach, 2_I4P)
   if (present(angular_velocity)) then
      call velocity_inputs(2_I4P, st, linear_velocity, [0._R8P, 0._R8P, angular_velocity], centre, lv, av, c, wv)
   else
      call velocity_inputs(2_I4P, st, linear_velocity, centre=centre, lv=lv, av=av, c=c, wv=wv)
   endif
   if (present(lengths)) then
      allocate(ll(size(lengths))) ; ll = lengths
   else
      allocate(ll(0))
   endif
   m = size(ll, kind=I4P)
   pr = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) ; pt = 0._R8P
   if (posed) then
      pr = pose%r ; pt = pose%t
      if (.not. present(centre)) c(1:2) = pose%t
   endif
   nv = curve%get_vertices_number()
   wc = curve%is_curvature_built()
   if (st == STATUS_OK .and. dev .and. nv > 0_I4P) then
      if (.not. curve%get_is_on_device()) then
         write(stderr, '(A)') 'fossil fd_closure_2d: device path refused -- call curve%enter_device first'
         st = STATUS_INVALID_INPUT
      elseif (wc .and. .not. curve%is_curvature_on_device()) then
         write(stderr, '(A)') 'fossil fd_closure_2d: device path refused -- build_curvature, then enter_device'
         st = STATUS_INVALID_INPUT
      endif
   endif
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(u(nv), v(nv), nxt(nv))
   do i = 1_I4P, nv ! the vertices, placed by the pose as `curve%move` places them
      uv = curve%get_vertex(i) ; nxt(i) = curve%get_next(i)
      if (posed) then
         call to_world_2d_oac(pr, pt, uv(1), uv(2), w1, w2) ; uv(1) = w1 ; uv(2) = w2
      endif
      u(i) = uv(1) ; v(i) = uv(2)
   enddo
   if (dev .and. nv > 0_I4P) then
      call closure_2d_device(curve, u, v, nxt, x, y, reach, ll, posed, pr, pt, wc, wv, lv, av(3), c, wi, res, closure)
   else
      call closure_2d_host(curve, u, v, nxt, x, y, reach, ll, pose, wc, wv, lv, av(3), c, wi, closure)
   endif
   endsubroutine fd_closure_2d

   subroutine ib_closure_reset(closure)
   !< Release a resident table from the device (`exit_device`), then reset it to the default (no ghost, arrays released):
   !< what `intent(out)` did, without freeing arrays that are still mapped.
   type(ib_closure_t), intent(inout) :: closure !< Table.
   type(ib_closure_t)                :: fresh   !< Default.

   if (closure%on_device) call closure%exit_device
   closure = fresh
   endsubroutine ib_closure_reset

   function ib_closure_mismatches(a, b) result(nmis)
   !< Number of differing entries of two tables, every field, reals and logicals by bits (the twin checks of the tests and
   !< benches). Arrays of different shapes count their larger size. Read on the host: `update_host` a resident table first.
   type(ib_closure_t), intent(in) :: a, b !< Tables.
   integer(I4P)                   :: nmis !< Mismatches.

   nmis = 0_I4P
   if (a%n /= b%n) nmis = nmis + 1_I4P
   if (a%dims /= b%dims) nmis = nmis + 1_I4P
   if (a%m /= b%m) nmis = nmis + 1_I4P
   if (a%ni /= b%ni .or. a%nj /= b%nj .or. a%nk /= b%nk) nmis = nmis + 1_I4P
   if (any(a%reach /= b%reach)) nmis = nmis + 1_I4P
   if ((a%with_curvature .neqv. b%with_curvature) .or. (a%with_velocity .neqv. b%with_velocity) .or. &
       (a%with_ilw .neqv. b%with_ilw)) nmis = nmis + 1_I4P
   nmis = nmis + mis_i(a%ijk, b%ijk) + mis_r(a%points, b%points) + mis_i(a%mask, b%mask) + mis_i(a%steps, b%steps)
   nmis = nmis + mis_r(a%distance, b%distance) + mis_r(a%foot, b%foot) + mis_r(a%normal, b%normal)
   nmis = nmis + mis_i(a%element, b%element) + mis_i(a%feature, b%feature) + mis_i(a%body, b%body)
   nmis = nmis + mis_i(a%patch, b%patch) + mis_r(a%lengths, b%lengths) + mis_r(a%probe, b%probe)
   nmis = nmis + mis_l(a%visible, b%visible) + mis_i(a%probe_cell, b%probe_cell) + mis_r(a%kappa, b%kappa)
   nmis = nmis + mis_r(a%direction, b%direction) + mis_r(a%velocity, b%velocity) + mis_r(a%ilw_t, b%ilw_t)
   nmis = nmis + mis_i(a%ilw_element, b%ilw_element) + mis_r(a%ilw_normal, b%ilw_normal)
   nmis = nmis + mis_r(a%ilw_kappa, b%ilw_kappa) + mis_r(a%ilw_direction, b%ilw_direction)
   endfunction ib_closure_mismatches

   ! the resident pattern

   subroutine enter_device(self)
   !< Map a host table to the device (copyin). Member arrays only, never `self` (#20). No-op if already mapped or if there
   !< is no ghost; a copy of a resident table is mapped as a new table.
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.

   if (is_mapped(self)) return
   self%on_device = .false. ; self%device_key = 0
   if (self%n < 1_I4P .or. .not. allocated(self%points)) return
   call groups(self, ga, gb, gc, gv, gi, gk)
   !$acc enter data copyin(self%ijk, self%points, self%mask, self%steps, self%distance, self%foot, self%normal, &
   !$acc&                  self%element, self%feature, self%body, self%patch) if(ga)
   !$acc enter data copyin(self%lengths, self%probe, self%visible, self%probe_cell) if(gb)
   !$acc enter data copyin(self%kappa, self%direction) if(gc)
   !$acc enter data copyin(self%velocity) if(gv)
   !$acc enter data copyin(self%ilw_t, self%ilw_element, self%ilw_normal) if(gi)
   !$acc enter data copyin(self%ilw_kappa, self%ilw_direction) if(gk)
   call mark_resident(self)
   endsubroutine enter_device

   subroutine update_host(self)
   !< Copy the mapped table back to the host; it stays mapped. No-op if not mapped; on a copy of a resident table a warning
   !< on stderr (its arrays are not on the device: update the original before copying it).
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.

   if (.not. self%on_device) return
   if (.not. is_mapped(self)) then
      write(stderr, '(A)') 'fossil ib_closure_t%update_host: a copy of a resident table is not on the device; '// &
                           'call update_host on the original before copying it'
      self%on_device = .false. ; self%device_key = 0
      return
   endif
   call groups(self, ga, gb, gc, gv, gi, gk)
   !$acc update host(self%ijk, self%points, self%mask, self%steps, self%distance, self%foot, self%normal, &
   !$acc&            self%element, self%feature, self%body, self%patch) if(ga)
   !$acc update host(self%lengths, self%probe, self%visible, self%probe_cell) if(gb)
   !$acc update host(self%kappa, self%direction) if(gc)
   !$acc update host(self%velocity) if(gv)
   !$acc update host(self%ilw_t, self%ilw_element, self%ilw_normal) if(gi)
   !$acc update host(self%ilw_kappa, self%ilw_direction) if(gk)
   endsubroutine update_host

   subroutine exit_device(self)
   !< Release the device copies; the host arrays are kept as they are (call `update_host` first to keep the device
   !< values). No-op if not mapped; a copy of a resident table only drops the flag.
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.

   if (.not. self%on_device) return
   if (is_mapped(self)) then
      call groups(self, ga, gb, gc, gv, gi, gk)
      !$acc exit data delete(self%ijk, self%points, self%mask, self%steps, self%distance, self%foot, self%normal, &
      !$acc&                 self%element, self%feature, self%body, self%patch) if(ga)
      !$acc exit data delete(self%lengths, self%probe, self%visible, self%probe_cell) if(gb)
      !$acc exit data delete(self%kappa, self%direction) if(gc)
      !$acc exit data delete(self%velocity) if(gv)
      !$acc exit data delete(self%ilw_t, self%ilw_element, self%ilw_normal) if(gi)
      !$acc exit data delete(self%ilw_kappa, self%ilw_direction) if(gk)
   endif
   self%on_device = .false. ; self%device_key = 0
   endsubroutine exit_device

   ! private: the two paths

   subroutine closure_3d_host(surface, geo, tree, x, y, z, reach, ll, pose, pr, pt, wc, wv, lv, av, c, wi, cl)
   !< The host path of `fd_closure_3d`: `ghost_points_3d` on geo (posed SAH BVH: the node states on the refitted tree,
   !< `posed_tree_3d`), `query_many` and `probe_many` on the surface (pose=), the ILW crossings by `line_crossings` on geo
   !< (posed: the BVH kernel on the refitted tree) and `ilw_select_oac`, their queries by `query_many`.
   type(surface_stl_object), intent(in), target   :: surface          !< Surface (queries).
   type(surface_stl_object), intent(in)           :: geo              !< Geometry of the ghosts and crossings.
   logical,                  intent(in)           :: tree             !< Posed SAH BVH: ghosts and crossings on the refit.
   real(R8P),                intent(in)           :: x(:), y(:), z(:) !< Point lines.
   integer(I4P),             intent(in)           :: reach(:)         !< Reach.
   real(R8P),                intent(in)           :: ll(:)            !< Probe lengths.
   type(rigid_pose_t),       intent(in), optional :: pose             !< Pose.
   real(R8P),                intent(in)           :: pr(3, 3), pt(3)  !< Pose (matrix, translation).
   logical,                  intent(in)           :: wc, wv, wi       !< Curvature, velocity, ILW.
   real(R8P),                intent(in)           :: lv(3), av(3), c(3) !< Rigid velocity.
   type(ib_closure_t),       intent(inout)        :: cl               !< Table.
   type(ghost_points_t)          :: gp                                !< Ghosts.
   type(aabb_node_object), allocatable :: nodem(:)                    !< Refitted tree (posed).
   real(R8P),        allocatable :: tvm(:,:)                          !< Moved vertices (posed).
   logical,          allocatable :: nd(:,:,:)                         !< Point inside flags (posed).
   type(vector_R8P), allocatable :: pts(:)                            !< Query points.
   real(R8P),        allocatable :: lens(:,:), cp(:,:), qd(:), qn(:,:), qk(:,:), qdir(:,:,:) !< Work.
   real(R8P),        allocatable :: ts(:)                             !< Crossings.
   integer(I4P),     allocatable :: ids(:), slot(:,:)                 !< Crossing facets, slots of the ILW queries.
   logical,          allocatable :: es(:)                             !< Entering flags.
   real(R8P)                     :: apos, bpos, elo, ehi, t0, t1      !< Slot window.
   integer(I4P)                  :: n, m, g, k, a, d, nc, best, sgn, axis, nq, q !< Work.

   if (tree) then
      allocate(nodem(0:ubound(surface%aabb%node, 1)), tvm(9, surface%get_facets_number()))
      call posed_tree_3d(surface, pr, pt, .false., nodem, tvm)
      allocate(nd(size(x), size(y), size(z)))
      call cut_cells_3d_nodes_tree(nodem, surface%aabb%payload, tvm, x, y, z, nd)
      call ghost_points_nodes_3d(nd, x, y, z, reach, gp)
   else
      call ghost_points_3d(geo, x, y, z, reach, gp)
   endif
   n = gp%n ; m = size(ll, kind=I4P)
   call closure_arrays(cl, 3_I4P, n, m, wc, wv, wi)
   cl%ni = size(x, kind=I4P) ; cl%nj = size(y, kind=I4P) ; cl%nk = size(z, kind=I4P) ; cl%reach = reach(1:3)
   cl%lengths = ll
   if (n < 1_I4P) return
   cl%ijk = gp%ijk ; cl%points = gp%points ; cl%mask = gp%mask ; cl%steps = gp%steps
   allocate(pts(n))
   do g = 1_I4P, n
      pts(g)%x = cl%points(1, g) ; pts(g)%y = cl%points(2, g) ; pts(g)%z = cl%points(3, g)
   enddo
   if (wc) then
      call surface%query_many(points=pts, distances=cl%distance, closest_points=cl%foot, facets=cl%element, &
                              features=cl%feature, normals=cl%normal, bodies=cl%body, patches=cl%patch, &
                              principal_curvatures=cl%kappa, principal_directions=cl%direction, pose=pose)
   else
      call surface%query_many(points=pts, distances=cl%distance, closest_points=cl%foot, facets=cl%element, &
                              features=cl%feature, normals=cl%normal, bodies=cl%body, patches=cl%patch, pose=pose)
   endif
   if (m > 0_I4P) then
      allocate(lens(m, n))
      do g = 1_I4P, n
         lens(:, g) = ll
      enddo
      call surface%probe_many(points=pts, lengths=lens, probe_points=cl%probe, visible=cl%visible, pose=pose)
      do g = 1_I4P, n
         do k = 1_I4P, m
            cl%probe_cell(1, k, g) = fd_cell_index_oac(cl%ni, x, cl%probe(1, k, g))
            cl%probe_cell(2, k, g) = fd_cell_index_oac(cl%nj, y, cl%probe(2, k, g))
            cl%probe_cell(3, k, g) = fd_cell_index_oac(cl%nk, z, cl%probe(3, k, g))
         enddo
      enddo
   endif
   if (wv) then
      do g = 1_I4P, n
         call rigid_velocity_3d_oac(lv, av, c, cl%foot(1, g), cl%foot(2, g), cl%foot(3, g), cl%velocity(1, g), &
                                    cl%velocity(2, g), cl%velocity(3, g))
      enddo
   endif
   if (.not. wi) return
   cl%ilw_t = 0._R8P ; cl%ilw_element = 0_I4P ; cl%ilw_normal = 0._R8P
   if (wc) then
      cl%ilw_kappa = 0._R8P ; cl%ilw_direction = 0._R8P
   endif
   allocate(cp(3, 6 * n), slot(2, 6 * n))
   nq = 0_I4P
   do g = 1_I4P, n
      do d = 1_I4P, 6_I4P
         if (cl%steps(d, g) < 1_I4P) cycle
         call ilw_slot_3d_oac(cl%ni, cl%nj, cl%nk, x, y, z, cl%ijk(1, g), cl%ijk(2, g), cl%ijk(3, g), d, cl%steps(d, g), &
                              axis, apos, bpos, elo, ehi, sgn, t0, t1)
         if (tree) then
            call tree_crossings(nodem, surface%aabb%payload, tvm, axis, apos, bpos, t0, t1, ts, ids, es, nc)
            call ilw_select_oac(nc, ts, ids, es, tvm, axis, apos, bpos, elo, ehi, sgn, best)
         else
            call geo%line_crossings(axis, apos, bpos, t0, t1, ts, ids, es, nc)
            call ilw_select_oac(nc, ts, ids, es, geo%tri_vertices, axis, apos, bpos, elo, ehi, sgn, best)
         endif
         if (best < 1_I4P) cycle
         cl%ilw_t(d, g) = ts(best) ; cl%ilw_element(d, g) = ids(best)
         nq = nq + 1_I4P ; slot(1, nq) = d ; slot(2, nq) = g
         call crossing_point_3d(axis, apos, bpos, ts(best), cp(1, nq), cp(2, nq), cp(3, nq))
      enddo
   enddo
   if (nq < 1_I4P) return
   deallocate(pts) ; allocate(pts(nq), qd(nq), qn(3, nq))
   do q = 1_I4P, nq
      pts(q)%x = cp(1, q) ; pts(q)%y = cp(2, q) ; pts(q)%z = cp(3, q)
   enddo
   if (wc) then
      allocate(qk(2, nq), qdir(3, 2, nq))
      call surface%query_many(points=pts, distances=qd, normals=qn, principal_curvatures=qk, principal_directions=qdir, &
                              pose=pose)
   else
      call surface%query_many(points=pts, distances=qd, normals=qn, pose=pose)
   endif
   do q = 1_I4P, nq
      d = slot(1, q) ; g = slot(2, q)
      cl%ilw_normal(:, d, g) = qn(:, q)
      if (wc) then
         cl%ilw_kappa(:, d, g) = qk(:, q)
         do a = 1, 2
            cl%ilw_direction(:, a, d, g) = qdir(:, a, q)
         enddo
      endif
   enddo
   endsubroutine closure_3d_host

   subroutine closure_3d_device(surface, geo, tree, x, y, z, reach, ll, posed, pr, pt, wc, wv, lv, av, c, wi, res, cl)
   !< The device path of `fd_closure_3d`: one pipeline, every stage on the device (see the module header). Posed (tree):
   !< the moved vertices and the refitted tree are built on the device (`posed_tree_3d`), nothing is copied in.
   type(surface_stl_object), intent(in), target :: surface          !< Surface (queries; mapped).
   type(surface_stl_object), intent(in), target :: geo              !< Geometry of the ghosts and crossings.
   logical,                  intent(in)         :: tree             !< Posed SAH BVH: ghosts and crossings on the refit.
   real(R8P),                intent(in)         :: x(:), y(:), z(:) !< Point lines.
   integer(I4P),             intent(in)         :: reach(:)         !< Reach.
   real(R8P),                intent(in)         :: ll(:)            !< Probe lengths.
   logical,                  intent(in)         :: posed            !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   logical,                  intent(in)         :: wc, wv, wi, res  !< Curvature, velocity, ILW, resident.
   real(R8P),                intent(in)         :: lv(3), av(3), c(3) !< Rigid velocity.
   type(ib_closure_t),       intent(inout)      :: cl               !< Table.
   type(ib_closure_t)                    :: w        !< The table built on the device.
   type(aabb_node_object),       pointer :: gnode(:) !< Alias of geo%aabb%node.
   type(facet_distance_payload), pointer :: gpay(:)  !< Alias of geo%aabb%payload.
   real(R8P),                    pointer :: gtv(:,:) !< Alias of geo%tri_vertices.
   type(aabb_node_object), allocatable, target :: nodem(:) !< Refitted tree (posed).
   real(R8P),              allocatable, target :: tvm(:,:) !< Moved vertices (posed).
   logical,      allocatable :: nd(:,:,:)            !< Point inside flags.
   integer(I4P), allocatable :: cnt(:), off(:), fpatch(:), fidx(:) !< Ghost counts, offsets, facet patches, indexes.
   real(R8P),    allocatable :: qp(:,:), cpb(:,:), nmb(:,:) !< Body-frame points, feet, normals.
   real(R8P)                 :: tol, dbx, dby, dbz   !< Visibility tolerance.
   type(vector_R8P)          :: bmin, bmax           !< Surface box.
   integer(I4P)              :: ni, nj, nk, np, n, m, f, r(3), nf !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nk = size(z, kind=I4P) ; np = ni * nj * nk
   m = size(ll, kind=I4P) ; r = reach(1:3) ; nf = surface%get_facets_number()
   bmin = surface%get_bmin() ; bmax = surface%get_bmax()
   dbx = bmax%x - bmin%x ; dby = bmax%y - bmin%y ; dbz = bmax%z - bmin%z
   tol = 1.e-9_R8P * sqrt(dbx * dbx + dby * dby + dbz * dbz) ! the tolerance of probe_many
   allocate(fpatch(nf))
   do f = 1_I4P, nf
      fpatch(f) = surface%facet_patch(f)
   enddo
   gnode => geo%aabb%node ; gpay => geo%aabb%payload ; gtv => geo%tri_vertices
   if (tree) then ! the tree and vertices of the moved surface, built on the device
      allocate(nodem(0:ubound(surface%aabb%node, 1)), tvm(9, nf))
      !$acc enter data create(nodem, tvm)
      call posed_tree_3d(surface, pr, pt, .true., nodem, tvm)
      gnode => nodem ; gtv => tvm
   endif
   allocate(nd(ni, nj, nk), cnt(np), off(np))
   !$acc data copyin(x, y, z, gnode, gpay, gtv, fpatch) create(nd, cnt, off)
   call cut_cells_3d_nodes_tree_present(gnode, gpay, gtv, x, y, z, nd)
   call ghost_count_dev(3_I4P, ni, nj, nk, nd, r, cnt)
   call links_scan_dev(np, cnt, off, n)
   call closure_arrays(w, 3_I4P, n, m, wc, wv, wi)
   w%ni = ni ; w%nj = nj ; w%nk = nk ; w%reach = r ; w%lengths = ll
   call map_create(w)
   if (n > 0_I4P) then
      call ghost_fill_dev(3_I4P, ni, nj, nk, nd, x, y, z, r, cnt, off, n, w%ijk, w%points, w%mask, w%steps)
      allocate(qp(3, n), cpb(3, n), nmb(3, n), fidx(n))
      !$acc data create(qp, cpb, nmb, fidx)
      call body_points_3d_dev(n, w%points, posed, pr, pt, qp)
      call surface%query_points_present(n, qp, .true., .false., w%distance, cpb, nmb, w%element, w%feature, fidx, w%body)
      call wall_3d_dev(n, nf, posed, pr, pt, cpb, nmb, w%element, fpatch, w%foot, w%normal, w%patch)
      if (wv) call velocity_3d_dev(n, lv, av, c, w%foot, w%velocity)
      if (wc) call curvature_3d_dev(surface, n, posed, pr, cpb, nmb, w%element, w%kappa, w%direction)
      if (m > 0_I4P) call probes_3d_dev(surface, n, m, posed, pr, pt, cpb, nmb, w%lengths, tol, ni, nj, nk, x, y, z, &
                                        w%probe, w%visible, w%probe_cell)
      if (wi) call ilw_3d_dev(surface, n, ni, nj, nk, x, y, z, w%ijk, w%steps, gnode, gpay, gtv, posed, pr, pt, wc, &
                              w%ilw_t, w%ilw_element, w%ilw_normal, size(w%ilw_kappa, 3, kind=I4P), w%ilw_kappa, &
                              w%ilw_direction)
      !$acc end data
   endif
   !$acc end data
   if (tree) then
      !$acc exit data delete(nodem, tvm)
   endif
   call finish_device(w, res, cl)
   endsubroutine closure_3d_device

   subroutine closure_2d_host(curve, u, v, nxt, x, y, reach, ll, pose, wc, wv, lv, aw, c, wi, cl)
   !< The host path of `fd_closure_2d`: `ghost_points_2d_flat` on the placed vertices, `query_many` and `probe_many` on
   !< the curve (pose=), the ILW crossings by `ilw_item_2d_oac` on the placed vertices, their queries by `query_many`.
   type(curve_object),    intent(in)           :: curve       !< Curve (queries).
   real(R8P),             intent(in)           :: u(:), v(:)  !< Vertices (placed).
   integer(I4P),          intent(in)           :: nxt(:)      !< Next vertex.
   real(R8P),             intent(in)           :: x(:), y(:)  !< Point lines.
   integer(I4P),          intent(in)           :: reach(:)    !< Reach.
   real(R8P),             intent(in)           :: ll(:)       !< Probe lengths.
   type(rigid_pose_2d_t), intent(in), optional :: pose        !< Pose.
   logical,               intent(in)           :: wc, wv, wi  !< Curvature, velocity, ILW.
   real(R8P),             intent(in)           :: lv(3), aw, c(3) !< Rigid velocity (angular: scalar).
   type(ib_closure_t),    intent(inout)        :: cl          !< Table.
   type(ghost_points_t)      :: gp                            !< Ghosts.
   real(R8P),    allocatable :: lens(:,:), cp(:,:), qd(:), qn(:,:), qk(:,:), qdir(:,:,:) !< Work.
   integer(I4P), allocatable :: rstart(:), rseg(:), cstart(:), cseg(:), slot(:,:) !< Row and column buckets, ILW slots.
   real(R8P)                 :: apos, elo, ehi, bt    !< Slot window, crossing.
   integer(I4P)              :: n, m, g, k, d, sgn, axis, nq, q, bs, nr, nc !< Work.

   call ghost_points_2d_flat(u, v, nxt, x, y, reach, gp)
   n = gp%n ; m = size(ll, kind=I4P)
   call closure_arrays(cl, 2_I4P, n, m, wc, wv, wi)
   cl%ni = size(x, kind=I4P) ; cl%nj = size(y, kind=I4P) ; cl%nk = 1_I4P ; cl%reach(1:2) = reach(1:2)
   cl%lengths = ll
   if (n < 1_I4P) return
   cl%ijk = gp%ijk ; cl%points = gp%points ; cl%mask = gp%mask ; cl%steps = gp%steps
   if (wc) then
      call curve%query_many(points=cl%points, distances=cl%distance, closest_points=cl%foot, segments=cl%element, &
                            vertices=cl%feature, normals=cl%normal, bodies=cl%body, patches=cl%patch, &
                            principal_curvatures=cl%kappa, principal_directions=cl%direction, pose=pose)
   else
      call curve%query_many(points=cl%points, distances=cl%distance, closest_points=cl%foot, segments=cl%element, &
                            vertices=cl%feature, normals=cl%normal, bodies=cl%body, patches=cl%patch, pose=pose)
   endif
   if (m > 0_I4P) then
      allocate(lens(m, n))
      do g = 1_I4P, n
         lens(:, g) = ll
      enddo
      call curve%probe_many(points=cl%points, lengths=lens, probe_points=cl%probe, visible=cl%visible, pose=pose)
      do g = 1_I4P, n
         do k = 1_I4P, m
            cl%probe_cell(1, k, g) = fd_cell_index_oac(cl%ni, x, cl%probe(1, k, g))
            cl%probe_cell(2, k, g) = fd_cell_index_oac(cl%nj, y, cl%probe(2, k, g))
         enddo
      enddo
   endif
   if (wv) then
      do g = 1_I4P, n
         call rigid_velocity_2d_oac(lv, aw, c, cl%foot(1, g), cl%foot(2, g), cl%velocity(1, g), cl%velocity(2, g))
      enddo
   endif
   if (.not. wi) return
   cl%ilw_t = 0._R8P ; cl%ilw_element = 0_I4P ; cl%ilw_normal = 0._R8P
   if (wc) then
      cl%ilw_kappa = 0._R8P ; cl%ilw_direction = 0._R8P
   endif
   call node_rows(u, v, nxt, y, rstart, rseg)  ! segments crossing each row (lines along u)
   call node_rows(v, u, nxt, x, cstart, cseg)  ! segments crossing each column (lines along v)
   nr = size(rseg, kind=I4P) ; nc = size(cseg, kind=I4P)
   allocate(cp(2, 4 * n), slot(2, 4 * n))
   nq = 0_I4P
   do g = 1_I4P, n
      do d = 1_I4P, 4_I4P
         if (cl%steps(d, g) < 1_I4P) cycle
         call ilw_slot_2d_oac(cl%ni, cl%nj, x, y, cl%ijk(1, g), cl%ijk(2, g), d, cl%steps(d, g), axis, apos, elo, ehi, sgn)
         if (axis == 1_I4P) then
            call ilw_item_2d_oac(size(u, kind=I4P), u, v, nxt, rstart(cl%ijk(2, g)), rstart(cl%ijk(2, g) + 1) - 1_I4P, nr, &
                                 rseg, axis, apos, elo, ehi, sgn, bt, bs)
         else
            call ilw_item_2d_oac(size(u, kind=I4P), u, v, nxt, cstart(cl%ijk(1, g)), cstart(cl%ijk(1, g) + 1) - 1_I4P, nc, &
                                 cseg, axis, apos, elo, ehi, sgn, bt, bs)
         endif
         if (bs < 1_I4P) cycle
         cl%ilw_t(d, g) = bt ; cl%ilw_element(d, g) = bs
         nq = nq + 1_I4P ; slot(1, nq) = d ; slot(2, nq) = g
         call crossing_point_2d(axis, apos, bt, cp(1, nq), cp(2, nq))
      enddo
   enddo
   if (nq < 1_I4P) return
   allocate(qd(nq), qn(2, nq))
   if (wc) then
      allocate(qk(1, nq), qdir(2, 1, nq))
      call curve%query_many(points=cp(:, 1:nq), distances=qd, normals=qn, principal_curvatures=qk, &
                            principal_directions=qdir, pose=pose)
   else
      call curve%query_many(points=cp(:, 1:nq), distances=qd, normals=qn, pose=pose)
   endif
   do q = 1_I4P, nq
      d = slot(1, q) ; g = slot(2, q)
      cl%ilw_normal(:, d, g) = qn(:, q)
      if (wc) then
         cl%ilw_kappa(1, d, g) = qk(1, q) ; cl%ilw_direction(:, 1, d, g) = qdir(:, 1, q)
      endif
   enddo
   endsubroutine closure_2d_host

   subroutine closure_2d_device(curve, u, v, nxt, x, y, reach, ll, posed, pr, pt, wc, wv, lv, aw, c, wi, res, cl)
   !< The device path of `fd_closure_2d`: one pipeline, every stage on the device (see the module header).
   type(curve_object),  intent(in), target :: curve          !< Curve (queries; mapped).
   real(R8P),           intent(in)         :: u(:), v(:)     !< Vertices (placed).
   integer(I4P),        intent(in)         :: nxt(:)         !< Next vertex.
   real(R8P),           intent(in)         :: x(:), y(:)     !< Point lines.
   integer(I4P),        intent(in)         :: reach(:)       !< Reach.
   real(R8P),           intent(in)         :: ll(:)          !< Probe lengths.
   logical,             intent(in)         :: posed          !< A pose.
   real(R8P),           intent(in)         :: pr(2, 2), pt(2) !< Pose.
   logical,             intent(in)         :: wc, wv, wi, res !< Curvature, velocity, ILW, resident.
   real(R8P),           intent(in)         :: lv(3), aw, c(3) !< Rigid velocity (angular: scalar).
   type(ib_closure_t),  intent(inout)      :: cl             !< Table.
   type(ib_closure_t)        :: w                            !< The table built on the device.
   logical,      allocatable :: nd(:,:,:)                    !< Point inside flags.
   integer(I4P), allocatable :: cnt(:), off(:), rstart(:), rseg(:), spatch(:), sgl(:), vxl(:) !< Work.
   real(R8P),    allocatable :: qp(:,:), cpb(:,:), nmb(:,:), pk(:) !< Body-frame points, feet, normals, curvatures.
   real(R8P)                 :: z(1)                         !< No z line.
   integer(I4P)              :: ni, nj, np, n, m, nr, nu, i, j, r(3) !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; np = ni * nj ; z = 0._R8P
   m = size(ll, kind=I4P) ; r(1:2) = reach(1:2) ; r(3) = 0_I4P ; nu = size(u, kind=I4P)
   allocate(spatch(nu))
   do i = 1_I4P, nu
      spatch(i) = curve%segment_patch(i)
   enddo
   call node_rows(u, v, nxt, y, rstart, rseg)
   nr = size(rseg, kind=I4P)
   allocate(nd(ni, nj, 1), cnt(np), off(np))
   !$acc data copyin(u, v, nxt, x, y, z, rstart, rseg, spatch) create(nd, cnt, off)
   !$acc parallel loop gang vector collapse(2) default(present) present(u, v, nxt, rstart, rseg, x, y, nd)
   do j = 1_I4P, nj
      do i = 1_I4P, ni
         nd(i, j, 1) = node_parity(x(i), y(j), rstart(j), rstart(j + 1) - 1_I4P, nr, rseg, nu, u, v, nxt)
      enddo
   enddo
   call ghost_count_dev(2_I4P, ni, nj, 1_I4P, nd, r, cnt)
   call links_scan_dev(np, cnt, off, n)
   call closure_arrays(w, 2_I4P, n, m, wc, wv, wi)
   w%ni = ni ; w%nj = nj ; w%nk = 1_I4P ; w%reach = r ; w%lengths = ll
   call map_create(w)
   if (n > 0_I4P) then
      call ghost_fill_dev(2_I4P, ni, nj, 1_I4P, nd, x, y, z, r, cnt, off, n, w%ijk, w%points, w%mask, w%steps)
      allocate(qp(2, n), cpb(2, n), nmb(2, n), sgl(n), vxl(n), pk(n))
      !$acc data create(qp, cpb, nmb, sgl, vxl, pk)
      call body_points_2d_dev(n, w%points, posed, pr, pt, qp)
      call curve%query_present(n, qp, .true., .false., w%distance, cpb, sgl, vxl, nmb, w%body)
      call wall_2d_dev(n, nu, posed, pr, pt, cpb, nmb, sgl, vxl, spatch, w%foot, w%normal, w%element, w%feature, w%patch)
      if (wv) call velocity_2d_dev(n, lv, aw, c, w%foot, w%velocity)
      if (wc) then
         call curve%curvature_present(n, sgl, cpb, pk)
         call curvature_2d_dev(n, posed, pr, nmb, pk, w%kappa, w%direction)
      endif
      if (m > 0_I4P) call probes_2d_dev(curve, n, m, posed, pr, pt, cpb, nmb, sgl, vxl, w%lengths, ni, nj, x, y, w%probe, &
                                        w%visible, w%probe_cell)
      if (wi) call ilw_2d_dev(curve, n, ni, nj, x, y, w%ijk, w%steps, u, v, nxt, posed, pr, pt, wc, w%ilw_t, &
                              w%ilw_element, w%ilw_normal, size(w%ilw_kappa, 3, kind=I4P), w%ilw_kappa, w%ilw_direction)
      !$acc end data
   endif
   !$acc end data
   call finish_device(w, res, cl)
   endsubroutine closure_2d_device

   subroutine posed_tree_3d(surface, pr, pt, dev, nodem, tvm)
   !< The vertices and the SAH BVH of the surface moved by the pose, without moving the surface (#109 WP1 step 06): tvm the
   !< facet vertices placed by `to_world` (the operations of `move`), nodem the tree of the surface with every box
   !< refitted to tvm as `refit` does it (a leaf the min/max of its facets' vertices, an internal node the union of its
   !< children's boxes, level by level from the deepest; children and payload slices kept). min and max are exact, so
   !< these are the vertices and boxes of the moved copy `move` makes, bit for bit, and every traversal on them is the
   !< traversal on that copy. dev: nodem and tvm present, written on the device (the surface mapped by `enter_device`);
   !< else on the host. The node depths are computed on the host from the given tree (integers, O(nodes)).
   type(surface_stl_object), intent(in), target :: surface         !< Surface (body frame).
   real(R8P),                intent(in)         :: pr(3, 3), pt(3) !< Pose.
   logical,                  intent(in)         :: dev             !< On the device.
   type(aabb_node_object),   intent(inout)      :: nodem(0:)       !< Refitted tree (same size as the surface's).
   real(R8P),                intent(inout)      :: tvm(:,:)        !< Moved vertices [9, facets].
   type(aabb_node_object),       pointer :: node_p(:)              !< Alias of surface%aabb%node.
   type(facet_distance_payload), pointer :: pay_p(:)               !< Alias of surface%aabb%payload.
   real(R8P),                    pointer :: tv_p(:,:)              !< Alias of surface%tri_vertices.
   integer(I4P), allocatable :: depth(:), lstart(:), lnode(:)      !< Node depths, internal nodes by depth (CSR).
   real(R8P)    :: r(3, 3), t(3), lox, loy, loz, hix, hiy, hiz      !< Pose copies, box.
   integer(I4P) :: nn, nf, n, lc, rc, f, k, maxd, lev, i, i0, i1    !< Work.

   node_p => surface%aabb%node ; pay_p => surface%aabb%payload ; tv_p => surface%tri_vertices
   nn = size(nodem, kind=I4P) ; nf = size(tvm, 2, kind=I4P) ; r = pr ; t = pt
   allocate(depth(0:nn - 1_I4P)) ; depth = 0_I4P
   do n = 0_I4P, nn - 1_I4P ! children have larger indices than their parent
      lc = node_p(n)%get_left_child() ; rc = node_p(n)%get_right_child()
      if (lc > 0_I4P) depth(lc) = depth(n) + 1_I4P
      if (rc > 0_I4P) depth(rc) = depth(n) + 1_I4P
   enddo
   maxd = maxval(depth)
   allocate(lstart(0:maxd + 1_I4P)) ; lstart = 0_I4P
   do n = 0_I4P, nn - 1_I4P
      if (node_p(n)%get_left_child() > 0_I4P .or. node_p(n)%get_right_child() > 0_I4P) &
         lstart(depth(n) + 1_I4P) = lstart(depth(n) + 1_I4P) + 1_I4P
   enddo
   do lev = 1_I4P, maxd + 1_I4P
      lstart(lev) = lstart(lev) + lstart(lev - 1_I4P)
   enddo
   allocate(lnode(max(1_I4P, lstart(maxd + 1_I4P))))
   block
      integer(I4P), allocatable :: cur(:) !< Cursor per depth.
      allocate(cur(0:maxd)) ; cur = lstart(0:maxd)
      do n = 0_I4P, nn - 1_I4P
         if (node_p(n)%get_left_child() > 0_I4P .or. node_p(n)%get_right_child() > 0_I4P) then
            cur(depth(n)) = cur(depth(n)) + 1_I4P ; lnode(cur(depth(n))) = n
         endif
      enddo
   endblock
   !$acc data copyin(r, t, lnode) if(dev)
   !$acc parallel loop gang vector default(present) present(tv_p, tvm) if(dev)
   do f = 1_I4P, nf
      call to_world_3d_oac(r, t, tv_p(1, f), tv_p(2, f), tv_p(3, f), tvm(1, f), tvm(2, f), tvm(3, f))
      call to_world_3d_oac(r, t, tv_p(4, f), tv_p(5, f), tv_p(6, f), tvm(4, f), tvm(5, f), tvm(6, f))
      call to_world_3d_oac(r, t, tv_p(7, f), tv_p(8, f), tv_p(9, f), tvm(7, f), tvm(8, f), tvm(9, f))
   enddo
   ! every node: children and payload slice; a leaf its box (an internal node is boxed below)
   !$acc parallel loop gang vector default(present) present(node_p, pay_p, tvm, nodem) &
   !$acc&  private(lox, loy, loz, hix, hiy, hiz, k, f) if(dev)
   do n = 0_I4P, nn - 1_I4P
      lox = huge(1._R8P) ; loy = huge(1._R8P) ; loz = huge(1._R8P) ; hix = -lox ; hiy = -loy ; hiz = -loz
      if (bvh_left_child_oac(node_p(n)) == 0_I4P .and. bvh_right_child_oac(node_p(n)) == 0_I4P) then
         do k = 0_I4P, bvh_payload_count_oac(node_p(n)) - 1_I4P
            f = pay_p(bvh_payload_first_oac(node_p(n)) + k)%facet_id
            lox = min(lox, min(tvm(1, f), tvm(4, f), tvm(7, f)))
            loy = min(loy, min(tvm(2, f), tvm(5, f), tvm(8, f)))
            loz = min(loz, min(tvm(3, f), tvm(6, f), tvm(9, f)))
            hix = max(hix, max(tvm(1, f), tvm(4, f), tvm(7, f)))
            hiy = max(hiy, max(tvm(2, f), tvm(5, f), tvm(8, f)))
            hiz = max(hiz, max(tvm(3, f), tvm(6, f), tvm(9, f)))
         enddo
      endif
      call bvh_node_refit_oac(node_p(n), lox, loy, loz, hix, hiy, hiz, nodem(n))
   enddo
   do lev = maxd - 1_I4P, 0_I4P, -1_I4P ! the internal nodes of one depth: their children are boxed
      i0 = lstart(lev) + 1_I4P ; i1 = lstart(lev + 1_I4P)
      !$acc parallel loop gang vector default(present) present(node_p, nodem, lnode) &
      !$acc&  private(n, lc, rc, lox, loy, loz, hix, hiy, hiz) if(dev)
      do i = i0, i1
         n = lnode(i) ; lc = bvh_left_child_oac(node_p(n)) ; rc = bvh_right_child_oac(node_p(n))
         lox = min(bvh_bmin_x_oac(nodem(lc)), bvh_bmin_x_oac(nodem(rc)))
         loy = min(bvh_bmin_y_oac(nodem(lc)), bvh_bmin_y_oac(nodem(rc)))
         loz = min(bvh_bmin_z_oac(nodem(lc)), bvh_bmin_z_oac(nodem(rc)))
         hix = max(bvh_bmax_x_oac(nodem(lc)), bvh_bmax_x_oac(nodem(rc)))
         hiy = max(bvh_bmax_y_oac(nodem(lc)), bvh_bmax_y_oac(nodem(rc)))
         hiz = max(bvh_bmax_z_oac(nodem(lc)), bvh_bmax_z_oac(nodem(rc)))
         call bvh_node_refit_oac(node_p(n), lox, loy, loz, hix, hiy, hiz, nodem(n))
      enddo
   enddo
   !$acc end data
   endsubroutine posed_tree_3d

   subroutine tree_crossings(node, payload, tv, axis, a, b, t0, t1, ts, ids, es, n)
   !< The sorted crossings of an axis-aligned line with a SAH BVH given as arrays, `line_crossings` (its BVH branch: the
   !< same kernel, the same growing buffer) on the refitted tree of the posed host path.
   type(aabb_node_object),       intent(in)  :: node(0:)   !< Tree nodes.
   type(facet_distance_payload), intent(in)  :: payload(:) !< Payload.
   real(R8P),                    intent(in)  :: tv(:,:)    !< Exact vertices [9, facets].
   integer(I4P),                 intent(in)  :: axis       !< Axis.
   real(R8P),                    intent(in)  :: a, b       !< Line position.
   real(R8P),                    intent(in)  :: t0, t1     !< Parameter range.
   real(R8P),    allocatable,    intent(out) :: ts(:)      !< Parameters [n].
   integer(I4P), allocatable,    intent(out) :: ids(:)     !< Facets [n].
   logical,      allocatable,    intent(out) :: es(:)      !< Entering flags [n].
   integer(I4P),                 intent(out) :: n          !< Crossings.
   integer(I4P) :: nmax !< Buffer size.
   logical      :: ovf  !< Overflow.

   nmax = 64_I4P
   do
      if (allocated(ts)) deallocate(ts, ids, es)
      allocate(ts(nmax), ids(nmax), es(nmax))
      call bvh_line_crossings_oac(node, payload, tv, axis, a, b, t0, t1, ts, ids, es, n, ovf)
      if (.not. ovf) exit
      nmax = 4_I4P * nmax
   enddo
   ts = ts(1:n) ; ids = ids(1:n) ; es = es(1:n)
   endsubroutine tree_crossings

   ! private: device stages (every array present; built without OpenACC the loops run on the host)

   subroutine ghost_count_dev(nd, ni, nj, nk, node_in, r, cnt)
   !< One thread per point: 1 if it is a ghost (`ghost_mask_oac`), else 0.
   integer(I4P), intent(in)    :: nd                  !< Dimensions.
   integer(I4P), intent(in)    :: ni, nj, nk          !< Points per axis.
   logical,      intent(in)    :: node_in(ni, nj, nk) !< Point inside flags (present).
   integer(I4P), intent(in)    :: r(3)                !< Reach.
   integer(I4P), intent(inout) :: cnt(ni * nj * nk)   !< Counts (present).
   integer(I4P) :: i, j, k, p, r1, r2, r3             !< Work.

   r1 = r(1) ; r2 = r(2) ; r3 = r(3)
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
   endsubroutine ghost_count_dev

   subroutine ghost_fill_dev(nd, ni, nj, nk, node_in, x, y, z, r, cnt, off, n, ijk, pts, mask, steps)
   !< One thread per point fills its ghost slot (`ghost_fill_oac`), the order of `ghost_points_3d/_2d`.
   integer(I4P), intent(in)    :: nd                  !< Dimensions.
   integer(I4P), intent(in)    :: ni, nj, nk          !< Points per axis.
   logical,      intent(in)    :: node_in(ni, nj, nk) !< Point inside flags (present).
   real(R8P),    intent(in)    :: x(:), y(:), z(:)    !< Point lines (present).
   integer(I4P), intent(in)    :: r(3)                !< Reach.
   integer(I4P), intent(in)    :: cnt(ni * nj * nk)   !< Counts (present).
   integer(I4P), intent(in)    :: off(ni * nj * nk)   !< Offsets (present).
   integer(I4P), intent(in)    :: n                   !< Ghosts.
   integer(I4P), intent(inout) :: ijk(nd, n)          !< Indexes (present).
   real(R8P),    intent(inout) :: pts(nd, n)          !< Coordinates (present).
   integer(I4P), intent(inout) :: mask(n)             !< Masks (present).
   integer(I4P), intent(inout) :: steps(2 * nd, n)    !< Steps (present).
   integer(I4P) :: i, j, k, p, g, r1, r2, r3          !< Work.

   r1 = r(1) ; r2 = r(2) ; r3 = r(3)
   !$acc parallel loop gang vector collapse(3) default(present) present(node_in, x, y, z, cnt, off, ijk, pts, mask, steps) &
   !$acc&  private(p, g)
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
   endsubroutine ghost_fill_dev

   subroutine body_points_3d_dev(n, pts, posed, pr, pt, qp)
   !< The query points in the body frame (`to_body`'s arithmetic), or copied without a pose.
   integer(I4P), intent(in)    :: n               !< Points.
   real(R8P),    intent(in)    :: pts(3, n)       !< World points (present).
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(3, 3), pt(3) !< Pose.
   real(R8P),    intent(inout) :: qp(3, n)        !< Body points (present).
   real(R8P)    :: r(3, 3), t(3)                  !< Pose copies.
   integer(I4P) :: g                              !< Counter.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(pts, qp) copyin(r, t)
   do g = 1_I4P, n
      if (posed) then
         call to_body_3d_oac(r, t, pts(1, g), pts(2, g), pts(3, g), qp(1, g), qp(2, g), qp(3, g))
      else
         qp(1, g) = pts(1, g) ; qp(2, g) = pts(2, g) ; qp(3, g) = pts(3, g)
      endif
   enddo
   endsubroutine body_points_3d_dev

   subroutine body_points_2d_dev(n, pts, posed, pr, pt, qp)
   !< 2-D twin of `body_points_3d_dev`.
   integer(I4P), intent(in)    :: n               !< Points.
   real(R8P),    intent(in)    :: pts(2, n)       !< World points (present).
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),    intent(inout) :: qp(2, n)        !< Body points (present).
   real(R8P)    :: r(2, 2), t(2)                  !< Pose copies.
   integer(I4P) :: g                              !< Counter.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(pts, qp) copyin(r, t)
   do g = 1_I4P, n
      if (posed) then
         call to_body_2d_oac(r, t, pts(1, g), pts(2, g), qp(1, g), qp(2, g))
      else
         qp(1, g) = pts(1, g) ; qp(2, g) = pts(2, g)
      endif
   enddo
   endsubroutine body_points_2d_dev

   subroutine wall_3d_dev(n, nf, posed, pr, pt, cpb, nmb, elem, fpatch, foot, normal, patch)
   !< Per ghost: the foot and the normal in the world (`to_world`, `rotate`), the patch of the facet.
   integer(I4P), intent(in)    :: n, nf           !< Ghosts, facets.
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(3, 3), pt(3) !< Pose.
   real(R8P),    intent(in)    :: cpb(3, n)       !< Body-frame feet (present).
   real(R8P),    intent(in)    :: nmb(3, n)       !< Body-frame normals (present).
   integer(I4P), intent(in)    :: elem(n)         !< Facets (present).
   integer(I4P), intent(in)    :: fpatch(nf)      !< Facet patches (present).
   real(R8P),    intent(inout) :: foot(3, n)      !< Feet (present).
   real(R8P),    intent(inout) :: normal(3, n)    !< Normals (present).
   integer(I4P), intent(inout) :: patch(n)        !< Patches (present).
   real(R8P)    :: r(3, 3), t(3)                  !< Copies.
   integer(I4P) :: g, f                           !< Work.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(cpb, nmb, elem, fpatch, foot, normal, patch) copyin(r, t) &
   !$acc&  private(f)
   do g = 1_I4P, n
      if (posed) then
         call to_world_3d_oac(r, t, cpb(1, g), cpb(2, g), cpb(3, g), foot(1, g), foot(2, g), foot(3, g))
         call rotate_3d_oac(r, nmb(1, g), nmb(2, g), nmb(3, g), normal(1, g), normal(2, g), normal(3, g))
      else
         foot(1, g) = cpb(1, g) ; foot(2, g) = cpb(2, g) ; foot(3, g) = cpb(3, g)
         normal(1, g) = nmb(1, g) ; normal(2, g) = nmb(2, g) ; normal(3, g) = nmb(3, g)
      endif
      f = elem(g) ; patch(g) = 0_I4P
      if (f >= 1_I4P .and. f <= nf) patch(g) = fpatch(f)
   enddo
   endsubroutine wall_3d_dev

   subroutine velocity_3d_dev(n, lv, av, c, foot, vel)
   !< Per ghost: the rigid velocity at the foot (`rigid_velocity_3d_oac`).
   integer(I4P), intent(in)    :: n                  !< Ghosts.
   real(R8P),    intent(in)    :: lv(3), av(3), c(3) !< Rigid velocity.
   real(R8P),    intent(in)    :: foot(3, n)         !< Feet (present).
   real(R8P),    intent(inout) :: vel(3, n)          !< Velocities (present).
   real(R8P)    :: vl(3), va(3), vc(3)               !< Copies.
   integer(I4P) :: g                                 !< Counter.

   vl = lv ; va = av ; vc = c
   !$acc parallel loop gang vector default(present) present(foot, vel) copyin(vl, va, vc)
   do g = 1_I4P, n
      call rigid_velocity_3d_oac(vl, va, vc, foot(1, g), foot(2, g), foot(3, g), vel(1, g), vel(2, g), vel(3, g))
   enddo
   endsubroutine velocity_3d_dev

   subroutine wall_2d_dev(n, nu, posed, pr, pt, cpb, nmb, sgl, vxl, spatch, foot, normal, elem, feat, patch)
   !< 2-D twin of `wall_3d_dev` (element = segment, feature = closest vertex).
   integer(I4P), intent(in)    :: n, nu           !< Ghosts, segments.
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),    intent(in)    :: cpb(2, n)       !< Body-frame feet (present).
   real(R8P),    intent(in)    :: nmb(2, n)       !< Body-frame normals (present).
   integer(I4P), intent(in)    :: sgl(n), vxl(n)  !< Segments, vertices (present).
   integer(I4P), intent(in)    :: spatch(nu)      !< Segment patches (present).
   real(R8P),    intent(inout) :: foot(2, n)      !< Feet (present).
   real(R8P),    intent(inout) :: normal(2, n)    !< Normals (present).
   integer(I4P), intent(inout) :: elem(n), feat(n), patch(n) !< Segments, vertices, patches (present).
   real(R8P)    :: r(2, 2), t(2)                  !< Copies.
   integer(I4P) :: g, s                           !< Work.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(cpb, nmb, sgl, vxl, spatch, foot, normal, elem, feat, patch) &
   !$acc&  copyin(r, t) private(s)
   do g = 1_I4P, n
      if (posed) then
         call to_world_2d_oac(r, t, cpb(1, g), cpb(2, g), foot(1, g), foot(2, g))
         call rotate_2d_oac(r, nmb(1, g), nmb(2, g), normal(1, g), normal(2, g))
      else
         foot(1, g) = cpb(1, g) ; foot(2, g) = cpb(2, g)
         normal(1, g) = nmb(1, g) ; normal(2, g) = nmb(2, g)
      endif
      s = sgl(g) ; elem(g) = s ; feat(g) = vxl(g) ; patch(g) = 0_I4P
      if (s >= 1_I4P .and. s <= nu) patch(g) = spatch(s)
   enddo
   endsubroutine wall_2d_dev

   subroutine velocity_2d_dev(n, lv, aw, c, foot, vel)
   !< Per ghost: the in-plane rigid velocity at the foot (`rigid_velocity_2d_oac`).
   integer(I4P), intent(in)    :: n               !< Ghosts.
   real(R8P),    intent(in)    :: lv(3), aw, c(3) !< Rigid velocity.
   real(R8P),    intent(in)    :: foot(2, n)      !< Feet (present).
   real(R8P),    intent(inout) :: vel(2, n)       !< Velocities (present).
   real(R8P)    :: vl(3), va, vc(3)               !< Copies.
   integer(I4P) :: g                              !< Counter.

   vl = lv ; va = aw ; vc = c
   !$acc parallel loop gang vector default(present) present(foot, vel) copyin(vl, vc)
   do g = 1_I4P, n
      call rigid_velocity_2d_oac(vl, va, vc, foot(1, g), foot(2, g), vel(1, g), vel(2, g))
   enddo
   endsubroutine velocity_2d_dev

   subroutine curvature_3d_dev(surface, n, posed, pr, cpb, nmb, elem, kappa, dirn)
   !< Per point: the shape operator at the body-frame foot (`foot_principal_oac`, the routine of `query_many`), the
   !< directions rotated to the world with a pose.
   type(surface_stl_object), intent(in), target :: surface    !< Surface (curvature mapped).
   integer(I4P),             intent(in)         :: n          !< Points.
   logical,                  intent(in)         :: posed      !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3)   !< Pose rotation.
   real(R8P),                intent(in)         :: cpb(3, n)  !< Body-frame feet (present).
   real(R8P),                intent(in)         :: nmb(3, n)  !< Body-frame normals (present).
   integer(I4P),             intent(in)         :: elem(n)    !< Facets (present).
   real(R8P),                intent(inout)      :: kappa(2, n) !< Curvatures (present).
   real(R8P),                intent(inout)      :: dirn(3, 2, n) !< Directions (present).
   real(R8P),    pointer :: ck_p(:,:), cd_p(:,:,:), tv_p(:,:) !< Aliases.
   integer(I4P), pointer :: cv_p(:,:)                         !< Alias.
   real(R8P)             :: r(3, 3), k1, k2, d1x, d1y, d1z, d2x, d2y, d2z !< Work.
   integer(I4P)          :: g                                 !< Counter.

   r = pr
   ck_p => surface%curv_k ; cd_p => surface%curv_d ; cv_p => surface%curv_vid ; tv_p => surface%tri_vertices
   !$acc parallel loop gang vector default(present) present(ck_p, cd_p, cv_p, tv_p, cpb, nmb, elem, kappa, dirn) &
   !$acc&  copyin(r) private(k1, k2, d1x, d1y, d1z, d2x, d2y, d2z)
   do g = 1_I4P, n
      call foot_principal_oac(tv_p, cv_p, ck_p, cd_p, elem(g), cpb(1, g), cpb(2, g), cpb(3, g), nmb(1, g), nmb(2, g), &
                              nmb(3, g), k1, k2, d1x, d1y, d1z, d2x, d2y, d2z)
      kappa(1, g) = k1 ; kappa(2, g) = k2
      if (posed) then
         call rotate_3d_oac(r, d1x, d1y, d1z, dirn(1, 1, g), dirn(2, 1, g), dirn(3, 1, g))
         call rotate_3d_oac(r, d2x, d2y, d2z, dirn(1, 2, g), dirn(2, 2, g), dirn(3, 2, g))
      else
         dirn(1, 1, g) = d1x ; dirn(2, 1, g) = d1y ; dirn(3, 1, g) = d1z
         dirn(1, 2, g) = d2x ; dirn(2, 2, g) = d2y ; dirn(3, 2, g) = d2z
      endif
   enddo
   endsubroutine curvature_3d_dev

   subroutine curvature_2d_dev(n, posed, pr, nmb, pk, kappa, dirn)
   !< Per point: the curvature and the tangent (-n_v, n_u) of the body-frame normal, rotated to the world with a pose.
   integer(I4P), intent(in)    :: n             !< Points.
   logical,      intent(in)    :: posed         !< A pose.
   real(R8P),    intent(in)    :: pr(2, 2)      !< Pose rotation.
   real(R8P),    intent(in)    :: nmb(2, n)     !< Body-frame normals (present).
   real(R8P),    intent(in)    :: pk(n)         !< Curvatures (present).
   real(R8P),    intent(inout) :: kappa(1, n)   !< Curvatures (present).
   real(R8P),    intent(inout) :: dirn(2, 1, n) !< Tangents (present).
   real(R8P)    :: r(2, 2), tu, tv              !< Work.
   integer(I4P) :: g                            !< Counter.

   r = pr
   !$acc parallel loop gang vector default(present) present(nmb, pk, kappa, dirn) copyin(r) private(tu, tv)
   do g = 1_I4P, n
      kappa(1, g) = pk(g)
      tu = -nmb(2, g) ; tv = nmb(1, g)
      if (posed) then
         call rotate_2d_oac(r, tu, tv, dirn(1, 1, g), dirn(2, 1, g))
      else
         dirn(1, 1, g) = tu ; dirn(2, 1, g) = tv
      endif
   enddo
   endsubroutine curvature_2d_dev

   subroutine probes_3d_dev(surface, n, m, posed, pr, pt, cpb, nmb, ll, tol, ni, nj, nk, x, y, z, probe, vis, pcell)
   !< Per (ghost, probe): the probe in the body frame (foot + l n), its visibility (`bvh_segment_blocked_oac` from the
   !< foot, the kernel of `probe_many`), the probe in the world, its enclosing cell on the lines.
   type(surface_stl_object), intent(in), target :: surface         !< Surface (mapped).
   integer(I4P),             intent(in)         :: n, m            !< Ghosts, probes per ghost.
   logical,                  intent(in)         :: posed           !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3) !< Pose.
   real(R8P),                intent(in)         :: cpb(3, n)       !< Body-frame feet (present).
   real(R8P),                intent(in)         :: nmb(3, n)       !< Body-frame normals (present).
   real(R8P),                intent(in)         :: ll(m)           !< Lengths (present).
   real(R8P),                intent(in)         :: tol             !< Visibility tolerance.
   integer(I4P),             intent(in)         :: ni, nj, nk      !< Points per axis.
   real(R8P),                intent(in)         :: x(ni), y(nj), z(nk) !< Point lines (present).
   real(R8P),                intent(inout)      :: probe(3, m, n)  !< Probes (present).
   logical,                  intent(inout)      :: vis(m, n)       !< Visibility (present).
   integer(I4P),             intent(inout)      :: pcell(3, m, n)  !< Cells (present).
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias of surface%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias of surface%aabb%payload.
   real(R8P),                    pointer :: tv_p(:,:)    !< Alias of surface%tri_vertices.
   real(R8P)    :: r(3, 3), t(3), cx, cy, cz, nx, ny, nz, px, py, pz, lk, wx, wy, wz !< Work.
   integer(I4P) :: g, k                                  !< Counters.

   r = pr ; t = pt
   node_p => surface%aabb%node ; payload_p => surface%aabb%payload ; tv_p => surface%tri_vertices
   !$acc parallel loop gang vector collapse(2) default(present) present(node_p, payload_p, tv_p, cpb, nmb, ll, x, y, z, &
   !$acc&  probe, vis, pcell) copyin(r, t) private(cx, cy, cz, nx, ny, nz, px, py, pz, lk, wx, wy, wz)
   do g = 1_I4P, n
      do k = 1_I4P, m
         cx = cpb(1, g) ; cy = cpb(2, g) ; cz = cpb(3, g)
         nx = nmb(1, g) ; ny = nmb(2, g) ; nz = nmb(3, g)
         lk = ll(k)
         px = cx + lk * nx ; py = cy + lk * ny ; pz = cz + lk * nz
         vis(k, g) = .not. bvh_segment_blocked_oac(node_p, payload_p, tv_p, cx, cy, cz, px, py, pz, tol)
         if (posed) then
            call to_world_3d_oac(r, t, px, py, pz, wx, wy, wz)
         else
            wx = px ; wy = py ; wz = pz
         endif
         probe(1, k, g) = wx ; probe(2, k, g) = wy ; probe(3, k, g) = wz
         pcell(1, k, g) = fd_cell_index_oac(ni, x, wx)
         pcell(2, k, g) = fd_cell_index_oac(nj, y, wy)
         pcell(3, k, g) = fd_cell_index_oac(nk, z, wz)
      enddo
   enddo
   endsubroutine probes_3d_dev

   subroutine probes_2d_dev(curve, n, m, posed, pr, pt, cpb, nmb, sgl, vxl, ll, ni, nj, x, y, probe, vis, pcell)
   !< 2-D twin of `probes_3d_dev`: the body-frame probes and visibility by `curve%probes_present`, then the world probes
   !< and the cells.
   type(curve_object), intent(in)    :: curve           !< Curve (mapped).
   integer(I4P),       intent(in)    :: n, m            !< Ghosts, probes per ghost.
   logical,            intent(in)    :: posed           !< A pose.
   real(R8P),          intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),          intent(in)    :: cpb(2, n)       !< Body-frame feet (present).
   real(R8P),          intent(in)    :: nmb(2, n)       !< Body-frame normals (present).
   integer(I4P),       intent(in)    :: sgl(n), vxl(n)  !< Foot segments, vertices (present).
   real(R8P),          intent(in)    :: ll(m)           !< Lengths (present).
   integer(I4P),       intent(in)    :: ni, nj          !< Points per axis.
   real(R8P),          intent(in)    :: x(ni), y(nj)    !< Point lines (present).
   real(R8P),          intent(inout) :: probe(2, m, n)  !< Probes (present).
   logical,            intent(inout) :: vis(m, n)       !< Visibility (present).
   integer(I4P),       intent(inout) :: pcell(2, m, n)  !< Cells (present).
   real(R8P), allocatable :: ppb(:,:,:)                 !< Body-frame probes.
   real(R8P)    :: r(2, 2), t(2), wu, wv                !< Work.
   integer(I4P) :: g, k                                 !< Counters.

   r = pr ; t = pt
   allocate(ppb(2, m, n))
   !$acc data create(ppb)
   call curve%probes_present(n, m, cpb, nmb, sgl, vxl, ll, ppb, vis)
   !$acc parallel loop gang vector collapse(2) default(present) present(ppb, x, y, probe, pcell) copyin(r, t) &
   !$acc&  private(wu, wv)
   do g = 1_I4P, n
      do k = 1_I4P, m
         if (posed) then
            call to_world_2d_oac(r, t, ppb(1, k, g), ppb(2, k, g), wu, wv)
         else
            wu = ppb(1, k, g) ; wv = ppb(2, k, g)
         endif
         probe(1, k, g) = wu ; probe(2, k, g) = wv
         pcell(1, k, g) = fd_cell_index_oac(ni, x, wu)
         pcell(2, k, g) = fd_cell_index_oac(nj, y, wv)
      enddo
   enddo
   !$acc end data
   endsubroutine probes_2d_dev

   subroutine ilw_3d_overflow_dev(n, ns, ni, nj, nk, x, y, z, ijk, steps, gnode, gpay, gtv, cnt, ilw_t, ilw_el)
   !< The ILW slots whose crossings overflowed the private buffer of `ilw_3d_dev` (cnt = -1) on exact-size buffers: their
   !< crossings counted (`bvh_line_count_oac`), a scan gives the buffers, `ilw_item_3d_oac` sorts and selects; cnt becomes
   !< 1 where a crossing was found, else 0 (the other slots keep theirs). Every array present.
   integer(I4P),                 intent(in)    :: n, ns, ni, nj, nk   !< Ghosts, slots, points per axis.
   real(R8P),                    intent(in)    :: x(ni), y(nj), z(nk) !< Point lines (present).
   integer(I4P),                 intent(in)    :: ijk(3, n)           !< Ghost indexes (present).
   integer(I4P),                 intent(in)    :: steps(6, n)         !< Steps (present).
   type(aabb_node_object),       intent(in)    :: gnode(0:)           !< Geometry tree (present).
   type(facet_distance_payload), intent(in)    :: gpay(:)             !< Geometry payload (present).
   real(R8P),                    intent(in)    :: gtv(:,:)            !< Geometry vertices (present).
   integer(I4P),                 intent(inout) :: cnt(ns)             !< Slot flags (present).
   real(R8P),                    intent(inout) :: ilw_t(6, n)         !< Crossing coordinates (present).
   integer(I4P),                 intent(inout) :: ilw_el(6, n)        !< Crossing facets (present).
   integer(I4P), allocatable :: oc(:), off(:), ids(:) !< Counts of the overflowed slots, offsets, buffer: facets.
   real(R8P),    allocatable :: ts(:)                 !< Buffer: parameters.
   logical,      allocatable :: es(:)                 !< Buffer: entering flags.
   real(R8P)    :: apos, bpos, elo, ehi, t0, t1, bt   !< Work.
   integer(I4P) :: l, g, d, s, axis, sgn, tot, bf     !< Work.

   allocate(oc(ns), off(ns))
   !$acc data create(oc, off)
   !$acc parallel loop gang vector default(present) present(x, y, z, ijk, steps, gnode, gpay, gtv, cnt, oc) &
   !$acc&  private(g, d, s, axis, apos, bpos, elo, ehi, sgn, t0, t1)
   do l = 1_I4P, ns
      oc(l) = 0_I4P
      if (cnt(l) < 0_I4P) then
         g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; s = steps(d, g)
         call ilw_slot_3d_oac(ni, nj, nk, x, y, z, ijk(1, g), ijk(2, g), ijk(3, g), d, s, axis, apos, bpos, elo, ehi, sgn, &
                              t0, t1)
         oc(l) = bvh_line_count_oac(gnode, gpay, gtv, axis, apos, bpos, t0, t1)
      endif
   enddo
   call links_scan_dev(ns, oc, off, tot)
   allocate(ts(tot + 1), ids(tot + 1), es(tot + 1))
   !$acc data create(ts, ids, es)
   !$acc parallel loop gang vector default(present) present(x, y, z, ijk, steps, gnode, gpay, gtv, cnt, oc, off, ts, ids, &
   !$acc&  es, ilw_t, ilw_el) private(g, d, s, axis, apos, bpos, elo, ehi, sgn, t0, t1, bt, bf)
   do l = 1_I4P, ns
      if (cnt(l) < 0_I4P) then
         g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; s = steps(d, g)
         call ilw_slot_3d_oac(ni, nj, nk, x, y, z, ijk(1, g), ijk(2, g), ijk(3, g), d, s, axis, apos, bpos, elo, ehi, sgn, &
                              t0, t1)
         call ilw_item_3d_oac(gnode, gpay, gtv, axis, apos, bpos, t0, t1, elo, ehi, sgn, oc(l), ts(off(l) + 1), &
                              ids(off(l) + 1), es(off(l) + 1), bt, bf)
         ilw_t(d, g) = bt ; ilw_el(d, g) = bf
         cnt(l) = 0_I4P ; if (bf > 0_I4P) cnt(l) = 1_I4P
      endif
   enddo
   !$acc end data
   !$acc end data
   endsubroutine ilw_3d_overflow_dev

   subroutine ilw_3d_dev(surface, n, ni, nj, nk, x, y, z, ijk, steps, gnode, gpay, gtv, posed, pr, pt, wc, ilw_t, ilw_el, &
                         ilw_n, nk_, ilw_k, ilw_d)
   !< The ILW crossings on the device: one thread per slot gathers the sorted crossings of its window into a private buffer
   !< of ILW_BUF entries and selects (`bvh_line_crossings_oac`, `ilw_select_oac`: one traversal, #109 WP1 step 06); the
   !< slots whose crossings overflow it are redone by the exact-size path (counted by `bvh_line_count_oac`, a scan gives
   !< the buffers, `ilw_item_3d_oac`): either way the same crossings in the same order, so the same selection. Then the
   !< found crossings are compacted (a scan), queried in the body frame (`query_points_present`) and scattered back
   !< (normal, curvature).
   type(surface_stl_object),     intent(in), target :: surface         !< Surface (queries; mapped).
   integer(I4P),                 intent(in)         :: n, ni, nj, nk   !< Ghosts, points per axis.
   real(R8P),                    intent(in)         :: x(ni), y(nj), z(nk) !< Point lines (present).
   integer(I4P),                 intent(in)         :: ijk(3, n)       !< Ghost indexes (present).
   integer(I4P),                 intent(in)         :: steps(6, n)     !< Steps (present).
   type(aabb_node_object),       intent(in)         :: gnode(0:)       !< Geometry tree (present).
   type(facet_distance_payload), intent(in)         :: gpay(:)         !< Geometry payload (present).
   real(R8P),                    intent(in)         :: gtv(:,:)        !< Geometry vertices (present).
   logical,                      intent(in)         :: posed           !< A pose.
   real(R8P),                    intent(in)         :: pr(3, 3), pt(3) !< Pose.
   logical,                      intent(in)         :: wc              !< Curvature.
   integer(I4P),                 intent(in)         :: nk_             !< Last extent of ilw_k, ilw_d (n with wc, else 0).
   real(R8P),                    intent(inout)      :: ilw_t(6, n)     !< Crossing coordinates (present).
   integer(I4P),                 intent(inout)      :: ilw_el(6, n)    !< Crossing facets (present).
   real(R8P),                    intent(inout)      :: ilw_n(3, 6, n)  !< Normals (present).
   real(R8P),                    intent(inout)      :: ilw_k(2, 6, nk_) !< Curvatures (present if wc).
   real(R8P),                    intent(inout)      :: ilw_d(3, 2, 6, nk_) !< Directions (present if wc).
   integer(I4P), allocatable :: cnt(:), off(:), qf(:), qft(:), qfi(:), qb(:) !< Work.
   real(R8P),    allocatable :: cq(:,:), qd(:), qc(:,:), qn(:,:), qk(:,:), qdir(:,:,:) !< Work.
   real(R8P)    :: r(3, 3), t(3), apos, bpos, elo, ehi, t0, t1, bt, cx, cy, cz !< Work.
   real(R8P)    :: tsl(ILW_BUF)                        !< Private buffer: parameters.
   integer(I4P) :: idsl(ILW_BUF)                       !< Private buffer: facets.
   logical      :: esl(ILW_BUF), ovf                   !< Private buffer: entering flags; overflow.
   integer(I4P) :: ns, l, g, d, s, axis, sgn, nq, q, bf, nf, best, nov !< Work.

   r = pr ; t = pt ; ns = 6_I4P * n
   allocate(cnt(ns), off(ns))
   !$acc data create(cnt, off) copyin(r, t)
   ! one traversal per slot into the private buffer; cnt = 1 if a crossing was found, -1 if the buffer overflowed
   !$acc parallel loop gang vector default(present) present(x, y, z, ijk, steps, gnode, gpay, gtv, cnt, ilw_t, ilw_el, &
   !$acc&  ilw_n) private(g, d, s, axis, apos, bpos, elo, ehi, sgn, t0, t1, bt, bf, tsl, idsl, esl, nf, ovf, best)
   do l = 1_I4P, ns
      g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; s = steps(d, g)
      bt = 0._R8P ; bf = 0_I4P ; cnt(l) = 0_I4P
      if (s > 0_I4P) then
         call ilw_slot_3d_oac(ni, nj, nk, x, y, z, ijk(1, g), ijk(2, g), ijk(3, g), d, s, axis, apos, bpos, elo, ehi, sgn, &
                              t0, t1)
         call bvh_line_crossings_oac(gnode, gpay, gtv, axis, apos, bpos, t0, t1, tsl, idsl, esl, nf, ovf)
         if (ovf) then
            cnt(l) = -1_I4P
         else
            call ilw_select_oac(nf, tsl, idsl, esl, gtv, axis, apos, bpos, elo, ehi, sgn, best)
            if (best > 0_I4P) then
               bt = tsl(best) ; bf = idsl(best)
            endif
         endif
      endif
      ilw_t(d, g) = bt ; ilw_el(d, g) = bf
      ilw_n(1, d, g) = 0._R8P ; ilw_n(2, d, g) = 0._R8P ; ilw_n(3, d, g) = 0._R8P
      if (bf > 0_I4P) cnt(l) = 1_I4P
   enddo
   nov = 0_I4P
   !$acc parallel loop gang vector default(present) present(cnt) reduction(+:nov)
   do l = 1_I4P, ns
      if (cnt(l) < 0_I4P) nov = nov + 1_I4P
   enddo
   if (nov > 0_I4P) call ilw_3d_overflow_dev(n, ns, ni, nj, nk, x, y, z, ijk, steps, gnode, gpay, gtv, cnt, ilw_t, ilw_el)
   if (wc) then
      !$acc parallel loop gang vector collapse(2) default(present) present(ilw_k, ilw_d)
      do g = 1_I4P, n
         do d = 1_I4P, 6_I4P
            ilw_k(1, d, g) = 0._R8P ; ilw_k(2, d, g) = 0._R8P
            ilw_d(1, 1, d, g) = 0._R8P ; ilw_d(2, 1, d, g) = 0._R8P ; ilw_d(3, 1, d, g) = 0._R8P
            ilw_d(1, 2, d, g) = 0._R8P ; ilw_d(2, 2, d, g) = 0._R8P ; ilw_d(3, 2, d, g) = 0._R8P
         enddo
      enddo
   endif
   call links_scan_dev(ns, cnt, off, nq)
   if (nq > 0_I4P) then
      allocate(cq(3, nq), qd(nq), qc(3, nq), qn(3, nq), qf(nq), qft(nq), qfi(nq), qb(nq))
      if (wc) then
         allocate(qk(2, nq), qdir(3, 2, nq))
      else
         allocate(qk(2, 1), qdir(3, 2, 1))
      endif
      !$acc data create(cq, qd, qc, qn, qf, qft, qfi, qb, qk, qdir)
      !$acc parallel loop gang vector default(present) present(x, y, z, ijk, steps, cnt, off, ilw_t, cq) &
      !$acc&  private(g, d, s, axis, apos, bpos, elo, ehi, sgn, t0, t1, q, cx, cy, cz)
      do l = 1_I4P, ns
         if (cnt(l) /= 1_I4P) cycle
         g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; s = steps(d, g) ; q = off(l) + 1_I4P
         call ilw_slot_3d_oac(ni, nj, nk, x, y, z, ijk(1, g), ijk(2, g), ijk(3, g), d, s, axis, apos, bpos, elo, ehi, sgn, &
                              t0, t1)
         call crossing_point_3d(axis, apos, bpos, ilw_t(d, g), cx, cy, cz)
         if (posed) then
            call to_body_3d_oac(r, t, cx, cy, cz, cq(1, q), cq(2, q), cq(3, q))
         else
            cq(1, q) = cx ; cq(2, q) = cy ; cq(3, q) = cz
         endif
      enddo
      call surface%query_points_present(nq, cq, .true., .false., qd, qc, qn, qf, qft, qfi, qb)
      if (wc) call curvature_3d_dev(surface, nq, posed, pr, qc, qn, qf, qk, qdir)
      !$acc parallel loop gang vector default(present) present(cnt, off, qn, ilw_n) private(g, d, q)
      do l = 1_I4P, ns
         if (cnt(l) /= 1_I4P) cycle
         g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; q = off(l) + 1_I4P
         if (posed) then
            call rotate_3d_oac(r, qn(1, q), qn(2, q), qn(3, q), ilw_n(1, d, g), ilw_n(2, d, g), ilw_n(3, d, g))
         else
            ilw_n(1, d, g) = qn(1, q) ; ilw_n(2, d, g) = qn(2, q) ; ilw_n(3, d, g) = qn(3, q)
         endif
      enddo
      if (wc) then
         !$acc parallel loop gang vector default(present) present(cnt, off, qk, qdir, ilw_k, ilw_d) private(g, d, q)
         do l = 1_I4P, ns
            if (cnt(l) /= 1_I4P) cycle
            g = (l - 1_I4P) / 6_I4P + 1_I4P ; d = mod(l - 1_I4P, 6_I4P) + 1_I4P ; q = off(l) + 1_I4P
            ilw_k(1, d, g) = qk(1, q) ; ilw_k(2, d, g) = qk(2, q)
            ilw_d(1, 1, d, g) = qdir(1, 1, q) ; ilw_d(2, 1, d, g) = qdir(2, 1, q) ; ilw_d(3, 1, d, g) = qdir(3, 1, q)
            ilw_d(1, 2, d, g) = qdir(1, 2, q) ; ilw_d(2, 2, d, g) = qdir(2, 2, q) ; ilw_d(3, 2, d, g) = qdir(3, 2, q)
         enddo
      endif
      !$acc end data
   endif
   !$acc end data
   endsubroutine ilw_3d_dev

   subroutine ilw_2d_dev(curve, n, ni, nj, x, y, ijk, steps, u, v, nxt, posed, pr, pt, wc, ilw_t, ilw_el, ilw_n, nk_, &
                         ilw_k, ilw_d)
   !< 2-D twin of `ilw_3d_dev`: one thread per slot walks the row or column bucket of its line (`ilw_item_2d_oac`, no
   !< buffer), then the found crossings are compacted, queried in the body frame (`curve%query_present`) and scattered.
   type(curve_object), intent(in)    :: curve           !< Curve (queries; mapped).
   integer(I4P),       intent(in)    :: n, ni, nj       !< Ghosts, points per axis.
   real(R8P),          intent(in)    :: x(ni), y(nj)    !< Point lines (present).
   integer(I4P),       intent(in)    :: ijk(2, n)       !< Ghost indexes (present).
   integer(I4P),       intent(in)    :: steps(4, n)     !< Steps (present).
   real(R8P),          intent(in)    :: u(:), v(:)      !< Vertices, placed (present).
   integer(I4P),       intent(in)    :: nxt(:)          !< Next vertex (present).
   logical,            intent(in)    :: posed           !< A pose.
   real(R8P),          intent(in)    :: pr(2, 2), pt(2) !< Pose.
   logical,            intent(in)    :: wc              !< Curvature.
   integer(I4P),       intent(in)    :: nk_             !< Last extent of ilw_k, ilw_d (n with wc, else 0).
   real(R8P),          intent(inout) :: ilw_t(4, n)     !< Crossing coordinates (present).
   integer(I4P),       intent(inout) :: ilw_el(4, n)    !< Crossing segments (present).
   real(R8P),          intent(inout) :: ilw_n(2, 4, n)  !< Normals (present).
   real(R8P),          intent(inout) :: ilw_k(1, 4, nk_) !< Curvatures (present if wc).
   real(R8P),          intent(inout) :: ilw_d(2, 1, 4, nk_) !< Tangents (present if wc).
   integer(I4P), allocatable :: rstart(:), rseg(:), cstart(:), cseg(:), cnt(:), off(:), qs(:), qv(:), qb(:) !< Work.
   real(R8P),    allocatable :: cq(:,:), qd(:), qc(:,:), qn(:,:), pk(:) !< Work.
   real(R8P)    :: r(2, 2), t(2), apos, elo, ehi, bt, cu, cv, tu, tv !< Work.
   integer(I4P) :: ns, l, g, d, s, axis, sgn, nq, q, bs, nu, nr, nc, k0, k1 !< Work.

   r = pr ; t = pt ; ns = 4_I4P * n ; nu = size(u, kind=I4P)
   call node_rows(u, v, nxt, y, rstart, rseg)
   call node_rows(v, u, nxt, x, cstart, cseg)
   nr = size(rseg, kind=I4P) ; nc = size(cseg, kind=I4P)
   allocate(cnt(ns), off(ns))
   !$acc data copyin(rstart, rseg, cstart, cseg, r, t) create(cnt, off)
   !$acc parallel loop gang vector default(present) present(x, y, ijk, steps, u, v, nxt, rstart, rseg, cstart, cseg, cnt, &
   !$acc&  ilw_t, ilw_el, ilw_n) private(g, d, s, axis, apos, elo, ehi, sgn, bt, bs, k0, k1)
   do l = 1_I4P, ns
      g = (l - 1_I4P) / 4_I4P + 1_I4P ; d = mod(l - 1_I4P, 4_I4P) + 1_I4P ; s = steps(d, g)
      bt = 0._R8P ; bs = 0_I4P
      if (s > 0_I4P) then
         call ilw_slot_2d_oac(ni, nj, x, y, ijk(1, g), ijk(2, g), d, s, axis, apos, elo, ehi, sgn)
         if (axis == 1_I4P) then
            k0 = rstart(ijk(2, g)) ; k1 = rstart(ijk(2, g) + 1) - 1_I4P
            call ilw_item_2d_oac(nu, u, v, nxt, k0, k1, nr, rseg, axis, apos, elo, ehi, sgn, bt, bs)
         else
            k0 = cstart(ijk(1, g)) ; k1 = cstart(ijk(1, g) + 1) - 1_I4P
            call ilw_item_2d_oac(nu, u, v, nxt, k0, k1, nc, cseg, axis, apos, elo, ehi, sgn, bt, bs)
         endif
      endif
      ilw_t(d, g) = bt ; ilw_el(d, g) = bs
      ilw_n(1, d, g) = 0._R8P ; ilw_n(2, d, g) = 0._R8P
      cnt(l) = 0_I4P ; if (bs > 0_I4P) cnt(l) = 1_I4P
   enddo
   if (wc) then
      !$acc parallel loop gang vector collapse(2) default(present) present(ilw_k, ilw_d)
      do g = 1_I4P, n
         do d = 1_I4P, 4_I4P
            ilw_k(1, d, g) = 0._R8P ; ilw_d(1, 1, d, g) = 0._R8P ; ilw_d(2, 1, d, g) = 0._R8P
         enddo
      enddo
   endif
   call links_scan_dev(ns, cnt, off, nq)
   if (nq > 0_I4P) then
      allocate(cq(2, nq), qd(nq), qc(2, nq), qn(2, nq), qs(nq), qv(nq), qb(nq), pk(nq))
      !$acc data create(cq, qd, qc, qn, qs, qv, qb, pk)
      !$acc parallel loop gang vector default(present) present(x, y, ijk, steps, cnt, off, ilw_t, cq) &
      !$acc&  private(g, d, s, axis, apos, elo, ehi, sgn, q, cu, cv)
      do l = 1_I4P, ns
         if (cnt(l) /= 1_I4P) cycle
         g = (l - 1_I4P) / 4_I4P + 1_I4P ; d = mod(l - 1_I4P, 4_I4P) + 1_I4P ; s = steps(d, g) ; q = off(l) + 1_I4P
         call ilw_slot_2d_oac(ni, nj, x, y, ijk(1, g), ijk(2, g), d, s, axis, apos, elo, ehi, sgn)
         call crossing_point_2d(axis, apos, ilw_t(d, g), cu, cv)
         if (posed) then
            call to_body_2d_oac(r, t, cu, cv, cq(1, q), cq(2, q))
         else
            cq(1, q) = cu ; cq(2, q) = cv
         endif
      enddo
      call curve%query_present(nq, cq, .true., .false., qd, qc, qs, qv, qn, qb)
      if (wc) call curve%curvature_present(nq, qs, qc, pk)
      !$acc parallel loop gang vector default(present) present(cnt, off, qn, ilw_n) private(g, d, q)
      do l = 1_I4P, ns
         if (cnt(l) /= 1_I4P) cycle
         g = (l - 1_I4P) / 4_I4P + 1_I4P ; d = mod(l - 1_I4P, 4_I4P) + 1_I4P ; q = off(l) + 1_I4P
         if (posed) then
            call rotate_2d_oac(r, qn(1, q), qn(2, q), ilw_n(1, d, g), ilw_n(2, d, g))
         else
            ilw_n(1, d, g) = qn(1, q) ; ilw_n(2, d, g) = qn(2, q)
         endif
      enddo
      if (wc) then
         !$acc parallel loop gang vector default(present) present(cnt, off, qn, pk, ilw_k, ilw_d) private(g, d, q, tu, tv)
         do l = 1_I4P, ns
            if (cnt(l) /= 1_I4P) cycle
            g = (l - 1_I4P) / 4_I4P + 1_I4P ; d = mod(l - 1_I4P, 4_I4P) + 1_I4P ; q = off(l) + 1_I4P
            ilw_k(1, d, g) = pk(q)
            tu = -qn(2, q) ; tv = qn(1, q)
            if (posed) then
               call rotate_2d_oac(r, tu, tv, ilw_d(1, 1, d, g), ilw_d(2, 1, d, g))
            else
               ilw_d(1, 1, d, g) = tu ; ilw_d(2, 1, d, g) = tv
            endif
         enddo
      endif
      !$acc end data
   endif
   !$acc end data
   endsubroutine ilw_2d_dev

   ! private: the table's arrays

   subroutine closure_arrays(cl, dims, n, m, wc, wv, wi)
   !< Allocate every array of the table (zero-size the ones not asked) and set its sizes and flags.
   type(ib_closure_t), intent(inout) :: cl         !< Table.
   integer(I4P),       intent(in)    :: dims, n, m !< Dimensions, ghosts, probes.
   logical,            intent(in)    :: wc, wv, wi !< Curvature, velocity, ILW.
   integer(I4P) :: nc, nvl, nil                    !< Sizes.

   cl%n = n ; cl%dims = dims ; cl%m = m
   cl%with_curvature = wc ; cl%with_velocity = wv ; cl%with_ilw = wi
   nc = 0_I4P ; if (wc) nc = n
   nvl = 0_I4P ; if (wv) nvl = n
   nil = 0_I4P ; if (wi) nil = n
   allocate(cl%ijk(dims, n), cl%points(dims, n), cl%mask(n), cl%steps(2 * dims, n), cl%distance(n), cl%foot(dims, n), &
            cl%normal(dims, n), cl%element(n), cl%feature(n), cl%body(n), cl%patch(n), cl%lengths(m), cl%probe(dims, m, n), &
            cl%visible(m, n), cl%probe_cell(dims, m, n), cl%kappa(dims - 1, nc), cl%direction(dims, dims - 1, nc), &
            cl%velocity(dims, nvl), cl%ilw_t(2 * dims, nil), cl%ilw_element(2 * dims, nil), &
            cl%ilw_normal(dims, 2 * dims, nil))
   if (wc) then
      allocate(cl%ilw_kappa(dims - 1, 2 * dims, nil), cl%ilw_direction(dims, dims - 1, 2 * dims, nil))
   else
      allocate(cl%ilw_kappa(dims - 1, 2 * dims, 0), cl%ilw_direction(dims, dims - 1, 2 * dims, 0))
   endif
   endsubroutine closure_arrays

   subroutine map_create(w)
   !< Create the device copies of the arrays of a table built on the device (unstructured data on a local table: its
   !< arrays are later moved into the caller's, the present table keyed on the unchanged host addresses).
   type(ib_closure_t), intent(inout) :: w !< Table.

   logical :: ga, gb, gc, gv, gi, gk !< Non-empty groups (zero-size arrays are never mapped).

   call groups(w, ga, gb, gc, gv, gi, gk)
   !$acc enter data create(w%ijk, w%points, w%mask, w%steps, w%distance, w%foot, w%normal, w%element, w%feature, w%body, &
   !$acc&                  w%patch) if(ga)
   !$acc enter data create(w%lengths, w%probe, w%visible, w%probe_cell) if(gb)
   !$acc update device(w%lengths) if(gb)
   !$acc enter data create(w%kappa, w%direction) if(gc)
   !$acc enter data create(w%velocity) if(gv)
   !$acc enter data create(w%ilw_t, w%ilw_element, w%ilw_normal) if(gi)
   !$acc enter data create(w%ilw_kappa, w%ilw_direction) if(gk)
   endsubroutine map_create

   subroutine finish_device(w, res, cl)
   !< End of a device path: resident, the arrays move into the caller's table, flagged; else they are copied back and
   !< released first.
   type(ib_closure_t), intent(inout) :: w   !< Table built on the device.
   logical,            intent(in)    :: res !< Resident.
   type(ib_closure_t), intent(inout) :: cl  !< Caller's table.

   logical :: ga, gb, gc, gv, gi, gk !< Non-empty groups.

   if (.not. res .or. w%n < 1_I4P) then
      call groups(w, ga, gb, gc, gv, gi, gk)
      !$acc exit data copyout(w%ijk, w%points, w%mask, w%steps, w%distance, w%foot, w%normal, w%element, w%feature, &
      !$acc&                  w%body, w%patch) if(ga)
      !$acc exit data copyout(w%lengths, w%probe, w%visible, w%probe_cell) if(gb)
      !$acc exit data copyout(w%kappa, w%direction) if(gc)
      !$acc exit data copyout(w%velocity) if(gv)
      !$acc exit data copyout(w%ilw_t, w%ilw_element, w%ilw_normal) if(gi)
      !$acc exit data copyout(w%ilw_kappa, w%ilw_direction) if(gk)
   endif
   cl%n = w%n ; cl%dims = w%dims ; cl%m = w%m ; cl%ni = w%ni ; cl%nj = w%nj ; cl%nk = w%nk ; cl%reach = w%reach
   cl%with_curvature = w%with_curvature ; cl%with_velocity = w%with_velocity ; cl%with_ilw = w%with_ilw
   call move_alloc(w%ijk, cl%ijk) ; call move_alloc(w%points, cl%points) ; call move_alloc(w%mask, cl%mask)
   call move_alloc(w%steps, cl%steps) ; call move_alloc(w%distance, cl%distance) ; call move_alloc(w%foot, cl%foot)
   call move_alloc(w%normal, cl%normal) ; call move_alloc(w%element, cl%element) ; call move_alloc(w%feature, cl%feature)
   call move_alloc(w%body, cl%body) ; call move_alloc(w%patch, cl%patch) ; call move_alloc(w%lengths, cl%lengths)
   call move_alloc(w%probe, cl%probe) ; call move_alloc(w%visible, cl%visible)
   call move_alloc(w%probe_cell, cl%probe_cell) ; call move_alloc(w%kappa, cl%kappa)
   call move_alloc(w%direction, cl%direction) ; call move_alloc(w%velocity, cl%velocity)
   call move_alloc(w%ilw_t, cl%ilw_t) ; call move_alloc(w%ilw_element, cl%ilw_element)
   call move_alloc(w%ilw_normal, cl%ilw_normal) ; call move_alloc(w%ilw_kappa, cl%ilw_kappa)
   call move_alloc(w%ilw_direction, cl%ilw_direction)
   if (res .and. cl%n > 0_I4P) call mark_resident(cl)
   endsubroutine finish_device

   subroutine groups(c, ga, gb, gc, gv, gi, gk)
   !< The non-empty array groups of a table (the device copies exist for these only): ghosts and wall, probes, curvature,
   !< velocity, ILW, ILW curvature.
   type(ib_closure_t), intent(in)  :: c                      !< Table.
   logical,            intent(out) :: ga, gb, gc, gv, gi, gk !< Groups.

   ga = c%n > 0_I4P
   gb = ga .and. c%m > 0_I4P
   gc = ga .and. c%with_curvature
   gv = ga .and. c%with_velocity
   gi = ga .and. c%with_ilw
   gk = gi .and. c%with_curvature
   endsubroutine groups

   function host_key(c) result(key)
   !< Host address of points (the key of a mapped table; 0 if points is not allocated or empty).
   type(ib_closure_t), intent(in), target :: c   !< Table.
   integer(c_intptr_t)                    :: key !< Address.

   key = 0
   if (allocated(c%points)) then
      if (size(c%points) > 0) key = transfer(c_loc(c%points), key)
   endif
   endfunction host_key

   subroutine mark_resident(c)
   !< Flag the arrays of `c` as mapped, keyed by the host address of points.
   type(ib_closure_t), intent(inout) :: c !< Table.

   c%on_device = .true. ; c%device_key = host_key(c)
   endsubroutine mark_resident

   logical function is_mapped(c)
   !< The flag is set and the arrays are the ones mapped (not a copy of a resident table).
   type(ib_closure_t), intent(in) :: c !< Table.

   is_mapped = .false.
   if (c%on_device) is_mapped = c%device_key == host_key(c)
   endfunction is_mapped

   subroutine velocity_inputs(dims, st, linear_velocity, angular_velocity, centre, lv, av, c, wv)
   !< Resolve the rigid velocity inputs (defaults 0) and check their sizes (dims, 3, dims); st = STATUS_INVALID_INPUT if
   !< one is wrong.
   integer(I4P), intent(in)           :: dims                !< Dimensions.
   integer(I4P), intent(inout)        :: st                  !< Status.
   real(R8P),    intent(in), optional :: linear_velocity(:)  !< Linear velocity.
   real(R8P),    intent(in), optional :: angular_velocity(:) !< Angular velocity [3] (2-D: [0, 0, omega]).
   real(R8P),    intent(in), optional :: centre(:)           !< Centre.
   real(R8P),    intent(out)          :: lv(3), av(3), c(3)  !< Resolved.
   logical,      intent(out)          :: wv                  !< Velocities asked.

   lv = 0._R8P ; av = 0._R8P ; c = 0._R8P
   wv = present(linear_velocity) .or. present(angular_velocity)
   if (present(linear_velocity)) then
      if (size(linear_velocity) /= dims) then
         st = STATUS_INVALID_INPUT
      else
         lv(1:dims) = linear_velocity
      endif
   endif
   if (present(angular_velocity)) then
      if (size(angular_velocity) /= 3) then
         st = STATUS_INVALID_INPUT
      else
         av = angular_velocity
      endif
   endif
   if (present(centre)) then
      if (size(centre) /= dims) then
         st = STATUS_INVALID_INPUT
      else
         c(1:dims) = centre
      endif
   endif
   endsubroutine velocity_inputs

   integer(I4P) function mis_i1(a, b) result(k)
   !< Differing entries of two integer rank-1 arrays (shapes differ: the larger size).
   integer(I4P), intent(in) :: a(:), b(:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a /= b)
   endif
   endfunction mis_i1

   integer(I4P) function mis_i2(a, b) result(k)
   !< Differing entries of two integer rank-2 arrays (shapes differ: the larger size).
   integer(I4P), intent(in) :: a(:,:), b(:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a /= b)
   endif
   endfunction mis_i2

   integer(I4P) function mis_i3(a, b) result(k)
   !< Differing entries of two integer rank-3 arrays (shapes differ: the larger size).
   integer(I4P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a /= b)
   endif
   endfunction mis_i3

   integer(I4P) function mis_r1(a, b) result(k)
   !< Differing entries of two real rank-1 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:), b(:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(a) /= bits(b))
   endif
   endfunction mis_r1

   integer(I4P) function mis_r2(a, b) result(k)
   !< Differing entries of two real rank-2 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:,:), b(:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(reshape(a, [size(a)])) /= bits(reshape(b, [size(b)])))
   endif
   endfunction mis_r2

   integer(I4P) function mis_r3(a, b) result(k)
   !< Differing entries of two real rank-3 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(reshape(a, [size(a)])) /= bits(reshape(b, [size(b)])))
   endif
   endfunction mis_r3

   integer(I4P) function mis_r4(a, b) result(k)
   !< Differing entries of two real rank-4 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:,:,:,:), b(:,:,:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(reshape(a, [size(a)])) /= bits(reshape(b, [size(b)])))
   endif
   endfunction mis_r4

   integer(I4P) function mis_l(a, b) result(k)
   !< Differing entries of two logical rank-2 arrays (sizes differ: the larger size).
   logical, intent(in) :: a(:,:), b(:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a .neqv. b)
   endif
   endfunction mis_l

   pure function bits(a) result(b)
   !< The bit patterns of reals.
   real(R8P), intent(in) :: a(:)        !< Reals.
   integer(int64)        :: b(size(a))  !< Bits.

   b = transfer(a, 0_int64, size(a))
   endfunction bits

   ! private kernels (`!$acc routine seq`, shared by the two paths)

   pure integer(I4P) function fd_cell_index_oac(n, e, c) result(idx)
   !< Kernel: the cell of the coordinate c on the lines e(1:n) by the low-side rule, e(i) < c <= e(i+1) gives i; 0 if c <= e(1),
   !< c > e(n) or c is NaN (binary search).
   !$acc routine seq
   integer(I4P), intent(in) :: n    !< Points.
   real(R8P),    intent(in) :: e(n) !< Lines (strictly increasing).
   real(R8P),    intent(in) :: c    !< Coordinate.
   integer(I4P)             :: lo, hi, mid !< Work.

   idx = 0_I4P
   if (n < 2_I4P) return
   if (.not. (c > e(1) .and. c <= e(n))) return
   lo = 1_I4P ; hi = n
   do while (hi - lo > 1_I4P)
      mid = (lo + hi) / 2_I4P
      if (c > e(mid)) then
         lo = mid
      else
         hi = mid
      endif
   enddo
   idx = lo
   endfunction fd_cell_index_oac

   pure subroutine ilw_slot_3d_oac(ni, nj, nk, x, y, z, i, j, k, d, s, axis, apos, bpos, elo, ehi, sgn, t0, t1)
   !< Kernel: the line of ILW slot d (2a-1: -a, 2a: +a) of the ghost (i, j, k) with s steps to the fluid: its axis, its
   !< position (a, b) in cyclic order, the coordinates elo < ehi of the last solid and the first fluid points, the walk
   !< direction sgn (the fluid at ehi if +1, at elo if -1) and the crossing window [elo - h, ehi + h], h = ehi - elo.
   !$acc routine seq
   integer(I4P), intent(in)  :: ni, nj, nk          !< Points per axis.
   real(R8P),    intent(in)  :: x(ni), y(nj), z(nk) !< Point lines.
   integer(I4P), intent(in)  :: i, j, k, d, s       !< Ghost, slot, steps.
   integer(I4P), intent(out) :: axis, sgn           !< Axis, direction.
   real(R8P),    intent(out) :: apos, bpos, elo, ehi, t0, t1 !< Line position, ends, window.
   integer(I4P)              :: p0, p1, lo, hi      !< Point indexes.
   real(R8P)                 :: h                   !< Spacing.

   axis = (d + 1_I4P) / 2_I4P ; sgn = 1_I4P ; if (mod(d, 2_I4P) == 1_I4P) sgn = -1_I4P
   if (axis == 1_I4P) then
      p0 = i + sgn * (s - 1_I4P) ; p1 = i + sgn * s
   elseif (axis == 2_I4P) then
      p0 = j + sgn * (s - 1_I4P) ; p1 = j + sgn * s
   else
      p0 = k + sgn * (s - 1_I4P) ; p1 = k + sgn * s
   endif
   lo = min(p0, p1) ; hi = max(p0, p1)
   if (axis == 1_I4P) then
      elo = x(lo) ; ehi = x(hi) ; apos = y(j) ; bpos = z(k)
   elseif (axis == 2_I4P) then
      elo = y(lo) ; ehi = y(hi) ; apos = z(k) ; bpos = x(i)
   else
      elo = z(lo) ; ehi = z(hi) ; apos = x(i) ; bpos = y(j)
   endif
   h = ehi - elo ; t0 = elo - h ; t1 = ehi + h
   endsubroutine ilw_slot_3d_oac

   pure subroutine ilw_slot_2d_oac(ni, nj, x, y, i, j, d, s, axis, apos, elo, ehi, sgn)
   !< Kernel: 2-D twin of `ilw_slot_3d_oac` (axis 1: the line along u at v = y(j); axis 2: along v at u = x(i)).
   !$acc routine seq
   integer(I4P), intent(in)  :: ni, nj        !< Points per axis.
   real(R8P),    intent(in)  :: x(ni), y(nj)  !< Point lines.
   integer(I4P), intent(in)  :: i, j, d, s    !< Ghost, slot, steps.
   integer(I4P), intent(out) :: axis, sgn     !< Axis, direction.
   real(R8P),    intent(out) :: apos, elo, ehi !< Line position, ends.
   integer(I4P)              :: p0, p1        !< Point indexes.

   axis = (d + 1_I4P) / 2_I4P ; sgn = 1_I4P ; if (mod(d, 2_I4P) == 1_I4P) sgn = -1_I4P
   if (axis == 1_I4P) then
      p0 = i + sgn * (s - 1_I4P) ; p1 = i + sgn * s
      elo = x(min(p0, p1)) ; ehi = x(max(p0, p1)) ; apos = y(j)
   else
      p0 = j + sgn * (s - 1_I4P) ; p1 = j + sgn * s
      elo = y(min(p0, p1)) ; ehi = y(max(p0, p1)) ; apos = x(i)
   endif
   endsubroutine ilw_slot_2d_oac

   pure subroutine ilw_select_oac(nc, ts, ids, es, tv, axis, apos, bpos, elo, ehi, sgn, best)
   !< Kernel: of the sorted crossings ts(1:nc) (facets ids, entering flags es) of the line along `axis` at (apos, bpos),
   !< the ones exactly between the perturbed points elo and ehi (`line_crossing_before`: before ehi, not before elo), the
   !< one nearest the fluid point: the last if sgn = +1 (fluid at ehi), the first if -1. best = its index, 0 if none.
   !$acc routine seq
   integer(I4P), intent(in)  :: nc          !< Crossings.
   real(R8P),    intent(in)  :: ts(nc)      !< Parameters (sorted: the order the selection relies on).
   integer(I4P), intent(in)  :: ids(nc)     !< Facets.
   logical,      intent(in)  :: es(nc)      !< Entering flags.
   real(R8P),    intent(in)  :: tv(:,:)     !< Exact facet vertices [9, facets].
   integer(I4P), intent(in)  :: axis, sgn   !< Axis, direction.
   real(R8P),    intent(in)  :: apos, bpos  !< Line position.
   real(R8P),    intent(in)  :: elo, ehi    !< Ends.
   integer(I4P), intent(out) :: best        !< Selected crossing (0: none).
   real(R8P)                 :: tri(3, 3)   !< Triangle.
   integer(I4P)              :: c, vv, ww   !< Counters.

   best = 0_I4P
   do c = 1_I4P, nc
      do vv = 1, 3
         do ww = 1, 3
            tri(ww, vv) = tv(3 * (vv - 1) + ww, ids(c))
         enddo
      enddo
      if (.not. line_crossing_before(tri, axis, apos, bpos, ehi, es(c))) cycle
      if (line_crossing_before(tri, axis, apos, bpos, elo, es(c))) cycle
      if (sgn > 0_I4P) then
         best = c
      elseif (best == 0_I4P) then
         best = c
      endif
   enddo
   endsubroutine ilw_select_oac

   pure subroutine ilw_item_3d_oac(node, payload, tv, axis, apos, bpos, t0, t1, elo, ehi, sgn, nc, ts, ids, es, bt, bf)
   !< Kernel: the sorted crossings of the slot's window into its buffers (nc = their count, `bvh_line_count_oac`), then
   !< `ilw_select_oac`: bt, bf = the crossing's coordinate and facet (0, 0 if none).
   !$acc routine seq
   type(aabb_node_object),       intent(in)    :: node(0:)   !< Tree nodes.
   type(facet_distance_payload), intent(in)    :: payload(:) !< Payload.
   real(R8P),                    intent(in)    :: tv(:,:)    !< Exact vertices [9, facets].
   integer(I4P),                 intent(in)    :: axis, sgn  !< Axis, direction.
   real(R8P),                    intent(in)    :: apos, bpos !< Line position.
   real(R8P),                    intent(in)    :: t0, t1     !< Window.
   real(R8P),                    intent(in)    :: elo, ehi   !< Ends.
   integer(I4P),                 intent(in)    :: nc         !< Buffer size.
   real(R8P),                    intent(inout) :: ts(nc)     !< Buffer: parameters.
   integer(I4P),                 intent(inout) :: ids(nc)    !< Buffer: facets.
   logical,                      intent(inout) :: es(nc)     !< Buffer: entering flags.
   real(R8P),                    intent(out)   :: bt         !< Crossing coordinate.
   integer(I4P),                 intent(out)   :: bf         !< Crossing facet.
   integer(I4P) :: nfound, best !< Work.
   logical      :: ovf          !< Overflow (never: nc is the count).

   bt = 0._R8P ; bf = 0_I4P
   if (nc < 1_I4P) return
   call bvh_line_crossings_oac(node, payload, tv, axis, apos, bpos, t0, t1, ts, ids, es, nfound, ovf)
   call ilw_select_oac(nfound, ts, ids, es, tv, axis, apos, bpos, elo, ehi, sgn, best)
   if (best > 0_I4P) then
      bt = ts(best) ; bf = ids(best)
   endif
   endsubroutine ilw_item_3d_oac

   pure subroutine ilw_item_2d_oac(nu, u, v, nxt, k0, k1, nb, bseg, axis, apos, elo, ehi, sgn, bt, bs)
   !< Kernel: the crossing of the 2-D ILW slot: of the segments bseg(k0:k1) crossing the line (a row or column bucket of
   !< `node_rows`), the ones exactly between the perturbed points elo and ehi (`line_crossing_before_2d`), the one nearest
   !< the fluid point (largest coordinate if sgn = +1, ties by the larger segment; smallest if -1, ties by the smaller):
   !< the order of `ilw_select_oac`. bt, bs = its coordinate and segment (0, 0 if none).
   !$acc routine seq
   integer(I4P), intent(in)  :: nu               !< Vertices.
   real(R8P),    intent(in)  :: u(nu), v(nu)     !< Vertices.
   integer(I4P), intent(in)  :: nxt(nu)          !< Next vertex.
   integer(I4P), intent(in)  :: k0, k1, nb       !< Bucket range, bucket array size.
   integer(I4P), intent(in)  :: bseg(nb)         !< Bucket segments (increasing).
   integer(I4P), intent(in)  :: axis, sgn        !< Axis (1: along u, 2: along v), direction.
   real(R8P),    intent(in)  :: apos, elo, ehi   !< Line position, ends.
   real(R8P),    intent(out) :: bt               !< Crossing coordinate.
   integer(I4P), intent(out) :: bs               !< Crossing segment.
   real(R8P)                 :: t                !< Parameter.
   integer(I4P)              :: kk, s, w         !< Work.
   logical                   :: hit, enter, take !< Work.

   bt = 0._R8P ; bs = 0_I4P
   do kk = k0, k1
      s = bseg(kk) ; w = nxt(s)
      if (w < 1_I4P) cycle
      call seg_line_hit_oac(u(s), v(s), u(w), v(w), axis, apos, hit, t, enter)
      if (.not. hit) cycle
      if (.not. line_crossing_before_2d(u(s), v(s), u(w), v(w), axis, apos, ehi, enter)) cycle
      if (line_crossing_before_2d(u(s), v(s), u(w), v(w), axis, apos, elo, enter)) cycle
      if (bs == 0_I4P) then
         take = .true.
      elseif (sgn > 0_I4P) then
         take = t > bt .or. (t == bt .and. s > bs)
      else
         take = t < bt .or. (t == bt .and. s < bs)
      endif
      if (take) then
         bt = t ; bs = s
      endif
   enddo
   endsubroutine ilw_item_2d_oac

   pure subroutine crossing_point_3d(axis, apos, bpos, t, cx, cy, cz)
   !< Kernel: the point of the line along `axis` at (apos, bpos) (cyclic order) with coordinate t.
   !$acc routine seq
   integer(I4P), intent(in)  :: axis          !< Axis.
   real(R8P),    intent(in)  :: apos, bpos, t !< Line position, coordinate.
   real(R8P),    intent(out) :: cx, cy, cz    !< Point.

   if (axis == 1_I4P) then
      cx = t ; cy = apos ; cz = bpos
   elseif (axis == 2_I4P) then
      cx = bpos ; cy = t ; cz = apos
   else
      cx = apos ; cy = bpos ; cz = t
   endif
   endsubroutine crossing_point_3d

   pure subroutine crossing_point_2d(axis, apos, t, cu, cv)
   !< Kernel: the point of the line along `axis` (1: u at v = apos; 2: v at u = apos) with coordinate t.
   !$acc routine seq
   integer(I4P), intent(in)  :: axis    !< Axis.
   real(R8P),    intent(in)  :: apos, t !< Line position, coordinate.
   real(R8P),    intent(out) :: cu, cv  !< Point.

   if (axis == 1_I4P) then
      cu = t ; cv = apos
   else
      cu = apos ; cv = t
   endif
   endsubroutine crossing_point_2d

   pure subroutine rigid_velocity_3d_oac(lv, av, c, px, py, pz, vx, vy, vz)
   !< Kernel: the rigid velocity V + omega x (p - c).
   !$acc routine seq
   real(R8P), intent(in)  :: lv(3), av(3), c(3) !< V, omega, centre.
   real(R8P), intent(in)  :: px, py, pz         !< Point.
   real(R8P), intent(out) :: vx, vy, vz         !< Velocity.
   real(R8P)              :: rx, ry, rz         !< p - c.

   rx = px - c(1) ; ry = py - c(2) ; rz = pz - c(3)
   vx = lv(1) + (av(2) * rz - av(3) * ry)
   vy = lv(2) + (av(3) * rx - av(1) * rz)
   vz = lv(3) + (av(1) * ry - av(2) * rx)
   endsubroutine rigid_velocity_3d_oac

   pure subroutine rigid_velocity_2d_oac(lv, aw, c, pu, pv, vu, vv)
   !< Kernel: the in-plane rigid velocity V + omega (-(p_v - c_v), p_u - c_u).
   !$acc routine seq
   real(R8P), intent(in)  :: lv(3), aw, c(3) !< V (first two), omega, centre (first two).
   real(R8P), intent(in)  :: pu, pv          !< Point.
   real(R8P), intent(out) :: vu, vv          !< Velocity.

   vu = lv(1) - aw * (pv - c(2))
   vv = lv(2) + aw * (pu - c(1))
   endsubroutine rigid_velocity_2d_oac

   pure subroutine to_body_3d_oac(r, t, x, y, z, bx, by, bz)
   !< Kernel: `to_body` (R^T (x - t)) on scalars, the same operations in the same order.
   !$acc routine seq
   real(R8P), intent(in)  :: r(3, 3), t(3) !< Pose.
   real(R8P), intent(in)  :: x, y, z       !< World point.
   real(R8P), intent(out) :: bx, by, bz    !< Body point.
   real(R8P)              :: dx, dy, dz    !< x - t.

   dx = x - t(1) ; dy = y - t(2) ; dz = z - t(3)
   bx = r(1, 1) * dx + r(2, 1) * dy + r(3, 1) * dz
   by = r(1, 2) * dx + r(2, 2) * dy + r(3, 2) * dz
   bz = r(1, 3) * dx + r(2, 3) * dy + r(3, 3) * dz
   endsubroutine to_body_3d_oac

   pure subroutine rotate_3d_oac(r, x, y, z, wx, wy, wz)
   !< Kernel: `rotate` (R v) on scalars.
   !$acc routine seq
   real(R8P), intent(in)  :: r(3, 3)    !< Rotation.
   real(R8P), intent(in)  :: x, y, z    !< Vector.
   real(R8P), intent(out) :: wx, wy, wz !< Rotated vector.

   wx = r(1, 1) * x + r(1, 2) * y + r(1, 3) * z
   wy = r(2, 1) * x + r(2, 2) * y + r(2, 3) * z
   wz = r(3, 1) * x + r(3, 2) * y + r(3, 3) * z
   endsubroutine rotate_3d_oac

   pure subroutine to_world_3d_oac(r, t, x, y, z, wx, wy, wz)
   !< Kernel: `to_world` (R x + t) on scalars.
   !$acc routine seq
   real(R8P), intent(in)  :: r(3, 3), t(3) !< Pose.
   real(R8P), intent(in)  :: x, y, z       !< Body point.
   real(R8P), intent(out) :: wx, wy, wz    !< World point.

   call rotate_3d_oac(r, x, y, z, wx, wy, wz)
   wx = wx + t(1) ; wy = wy + t(2) ; wz = wz + t(3)
   endsubroutine to_world_3d_oac

   pure subroutine to_body_2d_oac(r, t, x, y, bx, by)
   !< Kernel: 2-D `to_body` on scalars.
   !$acc routine seq
   real(R8P), intent(in)  :: r(2, 2), t(2) !< Pose.
   real(R8P), intent(in)  :: x, y          !< World point.
   real(R8P), intent(out) :: bx, by        !< Body point.
   real(R8P)              :: dx, dy        !< x - t.

   dx = x - t(1) ; dy = y - t(2)
   bx = r(1, 1) * dx + r(2, 1) * dy
   by = r(1, 2) * dx + r(2, 2) * dy
   endsubroutine to_body_2d_oac

   pure subroutine rotate_2d_oac(r, x, y, wx, wy)
   !< Kernel: 2-D `rotate` on scalars.
   !$acc routine seq
   real(R8P), intent(in)  :: r(2, 2) !< Rotation.
   real(R8P), intent(in)  :: x, y    !< Vector.
   real(R8P), intent(out) :: wx, wy  !< Rotated vector.

   wx = r(1, 1) * x + r(1, 2) * y
   wy = r(2, 1) * x + r(2, 2) * y
   endsubroutine rotate_2d_oac

   pure subroutine to_world_2d_oac(r, t, x, y, wx, wy)
   !< Kernel: 2-D `to_world` on scalars.
   !$acc routine seq
   real(R8P), intent(in)  :: r(2, 2), t(2) !< Pose.
   real(R8P), intent(in)  :: x, y          !< Body point.
   real(R8P), intent(out) :: wx, wy        !< World point.

   call rotate_2d_oac(r, x, y, wx, wy)
   wx = wx + t(1) ; wy = wy + t(2)
   endsubroutine to_world_2d_oac
endmodule fossil_fd_closure
