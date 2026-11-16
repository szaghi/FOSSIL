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
use fossil_cut_cells_2d,       only : node_rows, node_parity, node_flags
use fossil_cut_cells_3d,       only : cut_cells_3d_nodes, cut_cells_3d_nodes_tree, cut_cells_3d_nodes_tree_present
use fossil_cut_graph,          only : links_scan_dev
use fossil_ghost_points,       only : ghost_points_t, ghost_points_3d, ghost_points_2d_flat, ghost_lines_check, &
                                      ghost_mask_oac, ghost_fill_oac, ghost_points_nodes_3d, ghost_points_nodes_2d
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
public :: fd_closure_3d_bodies, fd_closure_2d_bodies
! the per-item kernels, shared by the two paths (public for the tests and for `fossil_surface_loads`)
public :: fd_cell_index_oac, ilw_select_oac, ilw_item_2d_oac
public :: to_world_3d_oac, rotate_3d_oac, to_world_2d_oac, rotate_2d_oac
public :: to_body_3d_oac, to_body_2d_oac ! for `fossil_field_transfer`
! the several-bodies stages, for the load samples of several bodies (`fossil_surface_loads`, #109 WP3 step 02)
public :: blocked_3d, blocked_2d, merge_clear, body_points_3d_dev, body_points_2d_dev, world_points_3d_dev, &
          world_points_2d_dev

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
   logical                   :: with_bodies = .false.    !< The several-bodies fields filled (`fd_closure_*_bodies`).
   integer(I4P)              :: bodies = 0_I4P           !< Surfaces (curves) of a bodies table (0: single body).
   integer(I4P)              :: overlap_points = 0_I4P   !< Points inside two or more surfaces (bodies table).
   real(R8P)                 :: shared_distance = 0._R8P !< Clearance below which a ghost is shared (bodies table).
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
   ! several bodies (`fd_closure_*_bodies`; zero-size in a single-body table): `body` is then the surface (curve) index
   integer(I4P), allocatable :: sub_body(:)            !< The owner's facet_body (shell) [n].
   real(R8P),    allocatable :: clearance(:)           !< Unsigned distance foot - nearest other surface [n] (huge: none).
   integer(I4P), allocatable :: second_body(:)         !< That nearest other surface [n] (0: none).
   real(R8P),    allocatable :: second_foot(:,:)       !< Its closest point to the foot [dims, n].
   logical,      allocatable :: shared(:)              !< clearance < shared_distance [n].
   logical,      allocatable :: overlap(:)             !< The ghost point lies inside two or more surfaces [n].
   integer(I4P), allocatable :: ilw_body(:,:)          !< Surface of the ILW crossing [2 dims, n] (0: none).
   logical                   :: on_device = .false.    !< The arrays are mapped on the device (resident).
   integer(c_intptr_t)       :: device_key = 0         !< Host address of points when mapped (a copy has another one).
   contains
      procedure, pass(self) :: enter_device !< Map the table to the device (copyin).
      procedure, pass(self) :: update_host  !< Copy the mapped table back to the host.
      procedure, pass(self) :: exit_device  !< Release the device copies.
endtype ib_closure_t

integer(I4P), parameter :: ILW_BUF = 16_I4P !< Private crossing buffer of an ILW slot (more: the exact-size path).

type :: posed_tree_t
   !< One body's refitted tree and moved vertices (the several-bodies host path).
   type(aabb_node_object), allocatable :: node(:) !< Refitted tree [0:].
   real(R8P),              allocatable :: tv(:,:) !< Moved vertices [9, facets].
endtype posed_tree_t

type :: placed_curve_t
   !< One curve's vertices placed by its pose, links and segment patches (the several-bodies 2-D paths).
   real(R8P),    allocatable :: u(:), v(:)  !< Vertices (placed).
   integer(I4P), allocatable :: nxt(:)      !< Next vertex.
   integer(I4P), allocatable :: spatch(:)   !< Segment patches.
endtype placed_curve_t

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

   subroutine fd_closure_3d_bodies(surfaces, x, y, z, reach, closure, status, poses, lengths, linear_velocities, &
                                   angular_velocities, centres, ilw, shared_distance, device, resident)
   !< The closure table of several bodies (#109 WP3 step 01; see the module header, "Several bodies"): surfaces(b) at
   !< poses(b) (absent: no body posed), body b's rigid velocity linear_velocities(:, b), angular_velocities(:, b) about
   !< centres(:, b) ([3, nb] each; the centre defaults as in `fd_closure_3d`, per body). `shared_distance`: the clearance
   !< below which a ghost is shared (default: the largest reach(a) times the largest spacing along axis a). The other
   !< arguments as `fd_closure_3d`. Every non-empty surface needs its SAH BVH (the default tree); the device path every
   !< non-empty surface mapped by `enter_device` (after `build_curvature` if built). Curvature is filled when every
   !< non-empty surface has it built. `status = STATUS_INVALID_INPUT` on bad lines or reach, no surface, size(poses) /= nb,
   !< a bad shape of a velocity array, a surface with an octree, or (device) an unmapped surface or curvature (a stderr
   !< message); the table is then reset (n = 0).
   type(surface_stl_object), intent(in), target    :: surfaces(:)             !< Surfaces (the bodies).
   real(R8P),                intent(in)            :: x(:), y(:), z(:)        !< Point lines (strictly increasing).
   integer(I4P),             intent(in)            :: reach(:)                !< Stencil half-width per axis [3], >= 1.
   type(ib_closure_t),       intent(inout)         :: closure                 !< Table (released and reset first).
   integer(I4P),             intent(out), optional :: status                  !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_t),       intent(in),  optional :: poses(:)                !< Poses of the bodies [nb].
   real(R8P),                intent(in),  optional :: lengths(:)              !< Probe lengths [m] (+ outward, absolute).
   real(R8P),                intent(in),  optional :: linear_velocities(:,:)  !< Velocities of the centres [3, nb].
   real(R8P),                intent(in),  optional :: angular_velocities(:,:) !< Angular velocities [3, nb].
   real(R8P),                intent(in),  optional :: centres(:,:)            !< Centres of the rotations [3, nb].
   logical,                  intent(in),  optional :: ilw                     !< ILW crossings (default .false.).
   real(R8P),                intent(in),  optional :: shared_distance         !< Shared clearance (default: see above).
   logical,                  intent(in),  optional :: device                  !< Device path (default .false.).
   logical,                  intent(in),  optional :: resident                !< Leave the table mapped (default .false.).
   real(R8P),    allocatable :: ll(:), lvs(:,:), avs(:,:), cs(:,:), prs(:,:,:), pts(:,:), tols(:) !< Per-body inputs.
   integer(I4P), allocatable :: nfs(:)                                       !< Facets per body.
   real(R8P)                 :: sd, dbx, dby, dbz                            !< Shared distance, box.
   type(vector_R8P)          :: bmin, bmax                                   !< Box.
   integer(I4P)              :: st, nb, b                                    !< Work.
   logical                   :: dev, res, wv, wi, wc, posed, some            !< Options.

   call ib_closure_reset(closure)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   wi = .false. ; if (present(ilw)) wi = ilw
   posed = present(poses)
   nb = size(surfaces, kind=I4P)
   st = ghost_lines_check(x, y, z, reach, 3_I4P)
   if (nb < 1_I4P) st = STATUS_INVALID_INPUT
   if (posed) then
      if (size(poses) /= nb) st = STATUS_INVALID_INPUT
   endif
   allocate(lvs(3, max(nb, 1_I4P)), avs(3, max(nb, 1_I4P)), cs(3, max(nb, 1_I4P)))
   call bodies_velocities(3_I4P, nb, st, linear_velocities, angular_velocities, centres, lvs, avs, cs, wv)
   if (present(lengths)) then
      allocate(ll(size(lengths))) ; ll = lengths
   else
      allocate(ll(0))
   endif
   allocate(prs(3, 3, max(nb, 1_I4P)), pts(3, max(nb, 1_I4P)), tols(max(nb, 1_I4P)), nfs(max(nb, 1_I4P)))
   wc = .true. ; some = .false.
   do b = 1_I4P, nb
      prs(:, :, b) = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3])
      pts(:, b) = 0._R8P
      if (posed .and. st == STATUS_OK) then
         prs(:, :, b) = poses(b)%r ; pts(:, b) = poses(b)%t
         if (.not. present(centres)) cs(:, b) = poses(b)%t
      endif
      nfs(b) = surfaces(b)%get_facets_number()
      bmin = surfaces(b)%get_bmin() ; bmax = surfaces(b)%get_bmax()
      dbx = bmax%x - bmin%x ; dby = bmax%y - bmin%y ; dbz = bmax%z - bmin%z
      tols(b) = 1.e-9_R8P * sqrt(dbx * dbx + dby * dby + dbz * dbz) ! the tolerance of probe_many
      if (nfs(b) < 1_I4P) cycle
      some = .true. ; wc = wc .and. surfaces(b)%is_curvature_built()
   enddo
   wc = wc .and. some
   do b = 1_I4P, nb
      if (nfs(b) < 1_I4P .or. st /= STATUS_OK) cycle
      if (.not. (surfaces(b)%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(surfaces(b)%aabb%node) .and. &
                 allocated(surfaces(b)%aabb%payload) .and. allocated(surfaces(b)%tri_vertices))) then
         write(stderr, '(A,I0,A)') 'fossil fd_closure_3d_bodies: refused -- surface ', b, ' has no SAH BVH'
         st = STATUS_INVALID_INPUT
      elseif (dev .and. .not. surfaces(b)%is_on_device) then
         write(stderr, '(A,I0,A)') 'fossil fd_closure_3d_bodies: device path refused -- call surface%enter_device '// &
                                   'first (surface ', b, ')'
         st = STATUS_INVALID_INPUT
      elseif (dev .and. wc .and. .not. surfaces(b)%is_curvature_on_device()) then
         write(stderr, '(A,I0,A)') 'fossil fd_closure_3d_bodies: device path refused -- build_curvature, then '// &
                                   'enter_device (surface ', b, ')'
         st = STATUS_INVALID_INPUT
      endif
   enddo
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   sd = default_shared(x, y, z, reach, 3_I4P) ; if (present(shared_distance)) sd = shared_distance
   if (dev .and. some) then
      call bodies_3d_device(surfaces, nfs, x, y, z, reach, ll, posed, prs, pts, tols, wc, wv, lvs, avs, cs, wi, sd, res, &
                            closure)
   elseif (posed) then
      call bodies_3d_host(surfaces, nfs, x, y, z, reach, ll, posed, prs, pts, tols, wc, wv, lvs, avs, cs, wi, sd, closure, &
                          poses)
   else
      call bodies_3d_host(surfaces, nfs, x, y, z, reach, ll, posed, prs, pts, tols, wc, wv, lvs, avs, cs, wi, sd, closure)
   endif
   endsubroutine fd_closure_3d_bodies

   subroutine fd_closure_2d_bodies(curves, x, y, reach, closure, status, poses, lengths, linear_velocities, &
                                   angular_velocities, centres, ilw, shared_distance, device, resident)
   !< 2-D twin of `fd_closure_3d_bodies` (#109 WP3 step 01): curves(b) at poses(b), in the curves' common (u, v) plane;
   !< `linear_velocities` [2, nb], `angular_velocities` [nb] (scalars, counter-clockwise), `centres` [2, nb]. The device
   !< path needs every non-empty curve mapped by `enter_device` (after `build_curvature` if built). The other arguments
   !< and the statuses as `fd_closure_3d_bodies` (no tree requirement: every curve has its segment BVH).
   type(curve_object),    intent(in), target    :: curves(:)              !< Curves (the bodies).
   real(R8P),             intent(in)            :: x(:), y(:)             !< Point lines (u, v; strictly increasing).
   integer(I4P),          intent(in)            :: reach(:)               !< Stencil half-width per axis [2], >= 1.
   type(ib_closure_t),    intent(inout)         :: closure                !< Table (released and reset first).
   integer(I4P),          intent(out), optional :: status                 !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_2d_t), intent(in),  optional :: poses(:)               !< Poses of the curves [nb].
   real(R8P),             intent(in),  optional :: lengths(:)             !< Probe lengths [m] (+ outward, absolute).
   real(R8P),             intent(in),  optional :: linear_velocities(:,:) !< Velocities of the centres [2, nb].
   real(R8P),             intent(in),  optional :: angular_velocities(:)  !< Angular velocities [nb] (scalars).
   real(R8P),             intent(in),  optional :: centres(:,:)           !< Centres of the rotations [2, nb].
   logical,               intent(in),  optional :: ilw                    !< ILW crossings (default .false.).
   real(R8P),             intent(in),  optional :: shared_distance        !< Shared clearance (default: as 3-D).
   logical,               intent(in),  optional :: device                 !< Device path (default .false.).
   logical,               intent(in),  optional :: resident               !< Leave the table mapped (default .false.).
   type(placed_curve_t), allocatable :: pc(:)                           !< Placed vertices per curve.
   real(R8P),    allocatable :: ll(:), lvs(:,:), avs(:,:), cs(:,:), prs(:,:,:), pts(:,:), a3(:,:) !< Per-body inputs.
   integer(I4P), allocatable :: nvs(:)                                  !< Vertices per curve.
   real(R8P)                 :: sd, uv(2), w1, w2, z(1)                 !< Shared distance, a vertex, no z line.
   integer(I4P)              :: st, nb, b, i                            !< Work.
   logical                   :: dev, res, wv, wi, wc, posed, some       !< Options.

   call ib_closure_reset(closure)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   wi = .false. ; if (present(ilw)) wi = ilw
   posed = present(poses)
   nb = size(curves, kind=I4P)
   z = 0._R8P
   st = ghost_lines_check(x, y, z, reach, 2_I4P)
   if (nb < 1_I4P) st = STATUS_INVALID_INPUT
   if (posed) then
      if (size(poses) /= nb) st = STATUS_INVALID_INPUT
   endif
   allocate(lvs(3, max(nb, 1_I4P)), avs(3, max(nb, 1_I4P)), cs(3, max(nb, 1_I4P)))
   if (present(angular_velocities)) then
      if (size(angular_velocities) /= nb) st = STATUS_INVALID_INPUT
      allocate(a3(3, max(nb, 1_I4P))) ; a3 = 0._R8P
      if (st == STATUS_OK) a3(3, 1:nb) = angular_velocities
      call bodies_velocities(2_I4P, nb, st, linear_velocities, a3, centres, lvs, avs, cs, wv)
   else
      call bodies_velocities(2_I4P, nb, st, linear_velocities, centres=centres, lvs=lvs, avs=avs, cs=cs, wv=wv)
   endif
   if (present(lengths)) then
      allocate(ll(size(lengths))) ; ll = lengths
   else
      allocate(ll(0))
   endif
   allocate(prs(2, 2, max(nb, 1_I4P)), pts(2, max(nb, 1_I4P)), nvs(max(nb, 1_I4P)), pc(max(nb, 1_I4P)))
   wc = .true. ; some = .false.
   do b = 1_I4P, nb
      prs(:, :, b) = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) ; pts(:, b) = 0._R8P
      if (posed .and. st == STATUS_OK) then
         prs(:, :, b) = poses(b)%r ; pts(:, b) = poses(b)%t
         if (.not. present(centres)) cs(1:2, b) = poses(b)%t
      endif
      nvs(b) = curves(b)%get_vertices_number()
      if (nvs(b) < 1_I4P) cycle
      some = .true. ; wc = wc .and. curves(b)%is_curvature_built()
      if (dev .and. st == STATUS_OK .and. .not. curves(b)%get_is_on_device()) then
         write(stderr, '(A,I0,A)') 'fossil fd_closure_2d_bodies: device path refused -- call curve%enter_device first '// &
                                   '(curve ', b, ')'
         st = STATUS_INVALID_INPUT
      endif
   enddo
   wc = wc .and. some
   do b = 1_I4P, nb
      if (nvs(b) < 1_I4P .or. .not. (dev .and. wc) .or. st /= STATUS_OK) cycle
      if (.not. curves(b)%is_curvature_on_device()) then
         write(stderr, '(A,I0,A)') 'fossil fd_closure_2d_bodies: device path refused -- build_curvature, then '// &
                                   'enter_device (curve ', b, ')'
         st = STATUS_INVALID_INPUT
      endif
   enddo
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   do b = 1_I4P, nb ! the vertices, placed by the pose as `curve%move` places them
      allocate(pc(b)%u(nvs(b)), pc(b)%v(nvs(b)), pc(b)%nxt(nvs(b)), pc(b)%spatch(nvs(b)))
      do i = 1_I4P, nvs(b)
         uv = curves(b)%get_vertex(i) ; pc(b)%nxt(i) = curves(b)%get_next(i) ; pc(b)%spatch(i) = curves(b)%segment_patch(i)
         if (posed) then
            call to_world_2d_oac(prs(:, :, b), pts(:, b), uv(1), uv(2), w1, w2) ; uv(1) = w1 ; uv(2) = w2
         endif
         pc(b)%u(i) = uv(1) ; pc(b)%v(i) = uv(2)
      enddo
   enddo
   sd = default_shared(x, y, z, reach, 2_I4P) ; if (present(shared_distance)) sd = shared_distance
   if (dev .and. some) then
      call bodies_2d_device(curves, nvs, pc, x, y, reach, ll, posed, prs, pts, wc, wv, lvs, avs, cs, wi, sd, res, closure)
   elseif (posed) then
      call bodies_2d_host(curves, nvs, pc, x, y, reach, ll, posed, prs, pts, wc, wv, lvs, avs, cs, wi, sd, closure, poses)
   else
      call bodies_2d_host(curves, nvs, pc, x, y, reach, ll, posed, prs, pts, wc, wv, lvs, avs, cs, wi, sd, closure)
   endif
   endsubroutine fd_closure_2d_bodies

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
   if ((a%with_bodies .neqv. b%with_bodies) .or. a%bodies /= b%bodies .or. a%overlap_points /= b%overlap_points .or. &
       transfer(a%shared_distance, 0_int64) /= transfer(b%shared_distance, 0_int64)) nmis = nmis + 1_I4P
   nmis = nmis + mis_i(a%sub_body, b%sub_body) + mis_r(a%clearance, b%clearance) + mis_i(a%second_body, b%second_body)
   nmis = nmis + mis_r(a%second_foot, b%second_foot) + mis_l1(a%shared, b%shared) + mis_l1(a%overlap, b%overlap)
   nmis = nmis + mis_i(a%ilw_body, b%ilw_body)
   endfunction ib_closure_mismatches

   ! the resident pattern

   subroutine enter_device(self)
   !< Map a host table to the device (copyin). Member arrays only, never `self` (#20). No-op if already mapped or if there
   !< is no ghost; a copy of a resident table is mapped as a new table.
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.
   logical :: gs, gsi                         !< Non-empty bodies groups.

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
   call groups_bodies(self, gs, gsi)
   !$acc enter data copyin(self%sub_body, self%clearance, self%second_body, self%second_foot, self%shared, &
   !$acc&                  self%overlap) if(gs)
   !$acc enter data copyin(self%ilw_body) if(gsi)
   call mark_resident(self)
   endsubroutine enter_device

   subroutine update_host(self)
   !< Copy the mapped table back to the host; it stays mapped. No-op if not mapped; on a copy of a resident table a warning
   !< on stderr (its arrays are not on the device: update the original before copying it).
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.
   logical :: gs, gsi                         !< Non-empty bodies groups.

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
   call groups_bodies(self, gs, gsi)
   !$acc update host(self%sub_body, self%clearance, self%second_body, self%second_foot, self%shared, self%overlap) if(gs)
   !$acc update host(self%ilw_body) if(gsi)
   endsubroutine update_host

   subroutine exit_device(self)
   !< Release the device copies; the host arrays are kept as they are (call `update_host` first to keep the device
   !< values). No-op if not mapped; a copy of a resident table only drops the flag.
   class(ib_closure_t), intent(inout) :: self !< Table.
   logical :: ga, gb, gc, gv, gi, gk          !< Non-empty groups.
   logical :: gs, gsi                         !< Non-empty bodies groups.

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
      call groups_bodies(self, gs, gsi)
      !$acc exit data delete(self%sub_body, self%clearance, self%second_body, self%second_foot, self%shared, &
      !$acc&                 self%overlap) if(gs)
      !$acc exit data delete(self%ilw_body) if(gsi)
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
   type(ghost_points_t)                :: gp                          !< Ghosts.
   type(aabb_node_object), allocatable :: nodem(:)                    !< Refitted tree (posed).
   real(R8P),              allocatable :: tvm(:,:)                    !< Moved vertices (posed).
   logical,                allocatable :: nd(:,:,:)                   !< Point inside flags (posed).
   integer(I4P)                        :: n, m                        !< Ghosts, probes.

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
   call fill_3d_host(surface, geo, tree, nodem, tvm, x, y, z, ll, pose, wc, wv, lv, av, c, wi, cl)
   endsubroutine closure_3d_host

   subroutine fill_3d_host(surface, geo, tree, nodem, tvm, x, y, z, ll, pose, wc, wv, lv, av, c, wi, cl)
   !< The host path of the table at given ghosts (cl%n > 0, its ijk, points, mask and steps set; every array allocated by
   !< `closure_arrays`): `query_many` and `probe_many` on the surface (pose=), the ILW crossings by `line_crossings` on
   !< geo (posed: the BVH kernel on the refitted tree nodem, tvm) and `ilw_select_oac`, their queries by `query_many`.
   type(surface_stl_object),            intent(in), target   :: surface          !< Surface (queries).
   type(surface_stl_object),            intent(in)           :: geo              !< Geometry of the crossings.
   logical,                             intent(in)           :: tree             !< Posed SAH BVH: crossings on nodem, tvm.
   type(aabb_node_object), allocatable, intent(in)           :: nodem(:)         !< Refitted tree (tree).
   real(R8P),              allocatable, intent(in)           :: tvm(:,:)         !< Moved vertices (tree).
   real(R8P),                           intent(in)           :: x(:), y(:), z(:) !< Point lines.
   real(R8P),                           intent(in)           :: ll(:)            !< Probe lengths.
   type(rigid_pose_t),                  intent(in), optional :: pose             !< Pose.
   logical,                             intent(in)           :: wc, wv, wi       !< Curvature, velocity, ILW.
   real(R8P),                           intent(in)           :: lv(3), av(3), c(3) !< Rigid velocity.
   type(ib_closure_t),                  intent(inout)        :: cl               !< Table.
   type(vector_R8P), allocatable :: pts(:)                            !< Query points.
   real(R8P),        allocatable :: lens(:,:), cp(:,:), qd(:), qn(:,:), qk(:,:), qdir(:,:,:) !< Work.
   real(R8P),        allocatable :: ts(:)                             !< Crossings.
   integer(I4P),     allocatable :: ids(:), slot(:,:)                 !< Crossing facets, slots of the ILW queries.
   logical,          allocatable :: es(:)                             !< Entering flags.
   real(R8P)                     :: apos, bpos, elo, ehi, t0, t1      !< Slot window.
   integer(I4P)                  :: n, m, g, k, a, d, nc, best, sgn, axis, nq, q !< Work.

   n = cl%n ; m = size(ll, kind=I4P)
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
   endsubroutine fill_3d_host

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
   integer(I4P), allocatable :: cnt(:), off(:), fpatch(:) !< Ghost counts, offsets, facet patches.
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
      call fill_3d_dev(surface, gnode, gpay, gtv, fpatch, tol, posed, pr, pt, wc, wv, lv, av, c, wi, x, y, z, w)
   endif
   !$acc end data
   if (tree) then
      !$acc exit data delete(nodem, tvm)
   endif
   call finish_device(w, res, cl)
   endsubroutine closure_3d_device

   subroutine fill_3d_dev(surface, gnode, gpay, gtv, fpatch, tol, posed, pr, pt, wc, wv, lv, av, c, wi, x, y, z, w)
   !< The device path of the table at given ghosts: w%n > 0, its arrays mapped (`map_create`), ijk, points, mask and steps
   !< filled on the device; x, y, z, gnode, gpay, gtv and fpatch present. The body-frame points, the queries, the wall
   !< data, the velocity, the curvature, the probes and the ILW crossings, every stage on the device.
   type(surface_stl_object),     intent(in), target :: surface          !< Surface (queries; mapped).
   type(aabb_node_object),       intent(in)         :: gnode(0:)        !< Geometry tree (present).
   type(facet_distance_payload), intent(in)         :: gpay(:)          !< Geometry payload (present).
   real(R8P),                    intent(in)         :: gtv(:,:)         !< Geometry vertices (present).
   integer(I4P),                 intent(in)         :: fpatch(:)        !< Facet patches (present).
   real(R8P),                    intent(in)         :: tol              !< Visibility tolerance.
   logical,                      intent(in)         :: posed            !< A pose.
   real(R8P),                    intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   logical,                      intent(in)         :: wc, wv, wi       !< Curvature, velocity, ILW.
   real(R8P),                    intent(in)         :: lv(3), av(3), c(3) !< Rigid velocity.
   real(R8P),                    intent(in)         :: x(:), y(:), z(:) !< Point lines (present).
   type(ib_closure_t),           intent(inout)      :: w                !< Table (mapped).
   integer(I4P), allocatable :: fidx(:)              !< Local feature indexes.
   real(R8P),    allocatable :: qp(:,:), cpb(:,:), nmb(:,:) !< Body-frame points, feet, normals.
   integer(I4P)              :: ni, nj, nk, n, m, nf !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nk = size(z, kind=I4P)
   n = w%n ; m = w%m ; nf = size(fpatch, kind=I4P)
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
   endsubroutine fill_3d_dev

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
   integer(I4P)              :: n, m                          !< Ghosts, probes.

   call ghost_points_2d_flat(u, v, nxt, x, y, reach, gp)
   n = gp%n ; m = size(ll, kind=I4P)
   call closure_arrays(cl, 2_I4P, n, m, wc, wv, wi)
   cl%ni = size(x, kind=I4P) ; cl%nj = size(y, kind=I4P) ; cl%nk = 1_I4P ; cl%reach(1:2) = reach(1:2)
   cl%lengths = ll
   if (n < 1_I4P) return
   cl%ijk = gp%ijk ; cl%points = gp%points ; cl%mask = gp%mask ; cl%steps = gp%steps
   call fill_2d_host(curve, u, v, nxt, x, y, ll, pose, wc, wv, lv, aw, c, wi, cl)
   endsubroutine closure_2d_host

   subroutine fill_2d_host(curve, u, v, nxt, x, y, ll, pose, wc, wv, lv, aw, c, wi, cl)
   !< 2-D twin of `fill_3d_host`: `query_many` and `probe_many` on the curve (pose=), the ILW crossings by
   !< `ilw_item_2d_oac` on the placed vertices, their queries by `query_many`.
   type(curve_object),    intent(in)           :: curve       !< Curve (queries).
   real(R8P),             intent(in)           :: u(:), v(:)  !< Vertices (placed).
   integer(I4P),          intent(in)           :: nxt(:)      !< Next vertex.
   real(R8P),             intent(in)           :: x(:), y(:)  !< Point lines.
   real(R8P),             intent(in)           :: ll(:)       !< Probe lengths.
   type(rigid_pose_2d_t), intent(in), optional :: pose        !< Pose.
   logical,               intent(in)           :: wc, wv, wi  !< Curvature, velocity, ILW.
   real(R8P),             intent(in)           :: lv(3), aw, c(3) !< Rigid velocity (angular: scalar).
   type(ib_closure_t),    intent(inout)        :: cl          !< Table.
   real(R8P),    allocatable :: lens(:,:), cp(:,:), qd(:), qn(:,:), qk(:,:), qdir(:,:,:) !< Work.
   integer(I4P), allocatable :: rstart(:), rseg(:), cstart(:), cseg(:), slot(:,:) !< Row and column buckets, ILW slots.
   real(R8P)                 :: apos, elo, ehi, bt    !< Slot window, crossing.
   integer(I4P)              :: n, m, g, k, d, sgn, axis, nq, q, bs, nr, nc !< Work.

   n = cl%n ; m = size(ll, kind=I4P)
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
   endsubroutine fill_2d_host

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
   integer(I4P), allocatable :: cnt(:), off(:), rstart(:), rseg(:), spatch(:) !< Work.
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
      call fill_2d_dev(curve, u, v, nxt, spatch, posed, pr, pt, wc, wv, lv, aw, c, wi, x, y, w)
   endif
   !$acc end data
   call finish_device(w, res, cl)
   endsubroutine closure_2d_device

   subroutine fill_2d_dev(curve, u, v, nxt, spatch, posed, pr, pt, wc, wv, lv, aw, c, wi, x, y, w)
   !< 2-D twin of `fill_3d_dev`: w%n > 0, mapped, its ghosts filled on the device; u, v, nxt, spatch, x, y present.
   type(curve_object), intent(in), target :: curve          !< Curve (queries; mapped).
   real(R8P),          intent(in)         :: u(:), v(:)     !< Vertices, placed (present).
   integer(I4P),       intent(in)         :: nxt(:)         !< Next vertex (present).
   integer(I4P),       intent(in)         :: spatch(:)      !< Segment patches (present).
   logical,            intent(in)         :: posed          !< A pose.
   real(R8P),          intent(in)         :: pr(2, 2), pt(2) !< Pose.
   logical,            intent(in)         :: wc, wv, wi     !< Curvature, velocity, ILW.
   real(R8P),          intent(in)         :: lv(3), aw, c(3) !< Rigid velocity (angular: scalar).
   real(R8P),          intent(in)         :: x(:), y(:)     !< Point lines (present).
   type(ib_closure_t), intent(inout)      :: w              !< Table (mapped).
   integer(I4P), allocatable :: sgl(:), vxl(:)              !< Foot segments, vertices.
   real(R8P),    allocatable :: qp(:,:), cpb(:,:), nmb(:,:), pk(:) !< Body-frame points, feet, normals, curvatures.
   integer(I4P)              :: ni, nj, n, m, nu            !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; n = w%n ; m = w%m ; nu = size(u, kind=I4P)
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
   endsubroutine fill_2d_dev

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

   ! private: several bodies (#109 WP3 step 01)

   subroutine bodies_3d_host(surfaces, nfs, x, y, z, reach, ll, posed, prs, pts, tols, wc, wv, lvs, avs, cs, wi, sd, cl, &
                             poses)
   !< The host path of `fd_closure_3d_bodies`: per body its node states (posed: on its refitted tree), their union and the
   !< ghosts; per body the single-body table at those ghosts (`fill_3d_host`) merged in body order (`merge_ghosts`,
   !< `merge_ilw`); the clearances (`query_many` of the owners' feet per other body, `merge_clear`), the probes' visibility
   !< against the other bodies (`blocked_3d`), the flags. The merge and visibility kernels are the device path's.
   type(surface_stl_object), intent(in), target   :: surfaces(:)        !< Surfaces.
   integer(I4P),             intent(in)           :: nfs(:)             !< Facets per body.
   real(R8P),                intent(in)           :: x(:), y(:), z(:)   !< Point lines.
   integer(I4P),             intent(in)           :: reach(:)           !< Reach.
   real(R8P),                intent(in)           :: ll(:)              !< Probe lengths.
   logical,                  intent(in)           :: posed              !< Poses given.
   real(R8P),                intent(in)           :: prs(:,:,:), pts(:,:) !< Poses [3, 3, nb], [3, nb].
   real(R8P),                intent(in)           :: tols(:)            !< Visibility tolerances [nb].
   logical,                  intent(in)           :: wc, wv, wi         !< Curvature, velocity, ILW.
   real(R8P),                intent(in)           :: lvs(:,:), avs(:,:), cs(:,:) !< Rigid velocities [3, nb].
   real(R8P),                intent(in)           :: sd                 !< Shared distance.
   type(ib_closure_t),       intent(inout)        :: cl                 !< Table.
   type(rigid_pose_t),       intent(in), optional :: poses(:)           !< Poses.
   type(posed_tree_t), allocatable :: pg(:)                             !< Refitted trees (posed).
   type(ghost_points_t)            :: gp                                !< Ghosts.
   type(ib_closure_t)              :: tb                                !< One body's table at the ghosts.
   logical,          allocatable   :: nd(:,:,:)                         !< Point inside flags.
   integer(I4P),     allocatable   :: ncnt(:,:,:)                       !< Bodies containing each point.
   type(vector_R8P), allocatable   :: fp(:)                             !< Owners' feet.
   real(R8P),        allocatable   :: qd(:), qc(:,:)                    !< Clearance queries.
   integer(I4P)                    :: nb, b, n, m, ni, nj, nk, g        !< Work.

   nb = size(nfs, kind=I4P) ; ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nk = size(z, kind=I4P)
   m = size(ll, kind=I4P)
   allocate(pg(nb), ncnt(ni, nj, nk)) ; ncnt = 0_I4P
   do b = 1_I4P, nb
      if (nfs(b) < 1_I4P) cycle
      if (posed) then
         allocate(pg(b)%node(0:ubound(surfaces(b)%aabb%node, 1)), pg(b)%tv(9, nfs(b)), nd(ni, nj, nk))
         call posed_tree_3d(surfaces(b), prs(:, :, b), pts(:, b), .false., pg(b)%node, pg(b)%tv)
         call cut_cells_3d_nodes_tree(pg(b)%node, surfaces(b)%aabb%payload, pg(b)%tv, x, y, z, nd)
      else
         call cut_cells_3d_nodes(surfaces(b), x, y, z, nd)
      endif
      where (nd) ncnt = ncnt + 1_I4P
      deallocate(nd)
   enddo
   allocate(nd(ni, nj, nk)) ; nd = ncnt > 0_I4P
   call ghost_points_nodes_3d(nd, x, y, z, reach, gp)
   n = gp%n
   call closure_arrays(cl, 3_I4P, n, m, wc, wv, wi) ; call closure_arrays_bodies(cl, nb, sd)
   cl%ni = ni ; cl%nj = nj ; cl%nk = nk ; cl%reach = reach(1:3) ; cl%lengths = ll
   cl%overlap_points = count(ncnt >= 2_I4P, kind=I4P)
   if (n < 1_I4P) return
   cl%ijk = gp%ijk ; cl%points = gp%points ; cl%mask = gp%mask ; cl%steps = gp%steps
   call bodies_init(.false., 3_I4P, n, size(cl%ilw_body, 2, kind=I4P), size(cl%ilw_kappa, 3, kind=I4P), cl%body, &
                    cl%distance, cl%clearance, cl%second_body, cl%second_foot, cl%ilw_t, cl%ilw_element, cl%ilw_normal, &
                    cl%ilw_body, cl%ilw_kappa, cl%ilw_direction)
   do b = 1_I4P, nb
      if (nfs(b) < 1_I4P) cycle
      call closure_arrays(tb, 3_I4P, n, m, wc, wv, wi)
      tb%ni = ni ; tb%nj = nj ; tb%nk = nk ; tb%reach = reach(1:3) ; tb%lengths = ll
      tb%ijk = cl%ijk ; tb%points = cl%points ; tb%mask = cl%mask ; tb%steps = cl%steps
      if (posed) then
         call fill_3d_host(surfaces(b), surfaces(b), .true., pg(b)%node, pg(b)%tv, x, y, z, ll, poses(b), wc, wv, &
                           lvs(:, b), avs(:, b), cs(:, b), wi, tb)
      else
         call fill_3d_host(surfaces(b), surfaces(b), .false., pg(b)%node, pg(b)%tv, x, y, z, ll, wc=wc, wv=wv, &
                           lv=lvs(:, b), av=avs(:, b), c=cs(:, b), wi=wi, cl=tb)
      endif
      call merge_table(.false., b, tb, cl)
      call ib_closure_reset(tb)
   enddo
   if (nb > 1_I4P) then
      allocate(fp(n), qd(n), qc(3, n))
      do g = 1_I4P, n
         fp(g)%x = cl%foot(1, g) ; fp(g)%y = cl%foot(2, g) ; fp(g)%z = cl%foot(3, g)
      enddo
      do b = 1_I4P, nb
         if (nfs(b) < 1_I4P) cycle
         if (posed) then
            call surfaces(b)%query_many(points=fp, distances=qd, closest_points=qc, pose=poses(b))
         else
            call surfaces(b)%query_many(points=fp, distances=qd, closest_points=qc)
         endif
         call merge_clear(.false., 3_I4P, n, b, qd, qc, cl%body, cl%clearance, cl%second_body, cl%second_foot)
         if (m > 0_I4P) call blocked_3d(.false., surfaces(b), b, n, m, posed, prs(:, :, b), pts(:, b), nb, tols, cl%body, &
                                        cl%foot, cl%probe, cl%visible)
      enddo
   endif
   call bodies_flags(.false., 3_I4P, n, ni, nj, nk, cl%ijk, ncnt, cl%clearance, sd, cl%shared, cl%overlap)
   endsubroutine bodies_3d_host

   subroutine bodies_3d_device(surfaces, nfs, x, y, z, reach, ll, posed, prs, pts, tols, wc, wv, lvs, avs, cs, wi, sd, &
                               res, cl)
   !< The device path of `fd_closure_3d_bodies`: the stages of the host path with the device kernels of `fd_closure_3d`,
   !< one launch set per body (`nodes_body_3d_dev`, `fill_body_3d_dev`, `clear_body_3d_dev`), merged in body order; no
   !< round trip between the stages.
   type(surface_stl_object), intent(in), target :: surfaces(:)          !< Surfaces (mapped).
   integer(I4P),             intent(in)         :: nfs(:)               !< Facets per body.
   real(R8P),                intent(in)         :: x(:), y(:), z(:)     !< Point lines.
   integer(I4P),             intent(in)         :: reach(:)             !< Reach.
   real(R8P),                intent(in)         :: ll(:)                !< Probe lengths.
   logical,                  intent(in)         :: posed                !< Poses given.
   real(R8P),                intent(in)         :: prs(:,:,:), pts(:,:) !< Poses.
   real(R8P),                intent(in)         :: tols(:)              !< Visibility tolerances [nb].
   logical,                  intent(in)         :: wc, wv, wi, res      !< Curvature, velocity, ILW, resident.
   real(R8P),                intent(in)         :: lvs(:,:), avs(:,:), cs(:,:) !< Rigid velocities.
   real(R8P),                intent(in)         :: sd                   !< Shared distance.
   type(ib_closure_t),       intent(inout)      :: cl                   !< Table.
   type(ib_closure_t)        :: w                                       !< The table built on the device.
   logical,      allocatable :: nd(:,:,:)                               !< Point inside flags.
   integer(I4P), allocatable :: ncnt(:,:,:), cnt(:), off(:)             !< Bodies per point, ghost counts, offsets.
   real(R8P),    allocatable :: tl(:)                                   !< Tolerances.
   integer(I4P)              :: nb, b, ni, nj, nk, np, n, m, r(3), op   !< Work.

   nb = size(nfs, kind=I4P) ; ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nk = size(z, kind=I4P)
   np = ni * nj * nk ; m = size(ll, kind=I4P) ; r = reach(1:3)
   allocate(nd(ni, nj, nk), ncnt(ni, nj, nk), cnt(np), off(np), tl(nb)) ; tl = tols(1:nb)
   !$acc data copyin(x, y, z, tl) create(nd, ncnt, cnt, off)
   call zero_counts_dev(np, ncnt)
   do b = 1_I4P, nb
      if (nfs(b) < 1_I4P) cycle
      call nodes_body_3d_dev(surfaces(b), posed, prs(:, :, b), pts(:, b), x, y, z, nd, ncnt)
   enddo
   call union_dev(np, ncnt, nd, op)
   call ghost_count_dev(3_I4P, ni, nj, nk, nd, r, cnt)
   call links_scan_dev(np, cnt, off, n)
   call closure_arrays(w, 3_I4P, n, m, wc, wv, wi) ; call closure_arrays_bodies(w, nb, sd)
   w%ni = ni ; w%nj = nj ; w%nk = nk ; w%reach = r ; w%lengths = ll ; w%overlap_points = op
   call map_create(w)
   if (n > 0_I4P) then
      call ghost_fill_dev(3_I4P, ni, nj, nk, nd, x, y, z, r, cnt, off, n, w%ijk, w%points, w%mask, w%steps)
      call bodies_init(.true., 3_I4P, n, size(w%ilw_body, 2, kind=I4P), size(w%ilw_kappa, 3, kind=I4P), w%body, &
                       w%distance, w%clearance, w%second_body, w%second_foot, w%ilw_t, w%ilw_element, w%ilw_normal, &
                       w%ilw_body, w%ilw_kappa, w%ilw_direction)
      do b = 1_I4P, nb
         if (nfs(b) < 1_I4P) cycle
         call fill_body_3d_dev(surfaces(b), b, posed, prs(:, :, b), pts(:, b), tl(b), wc, wv, lvs(:, b), avs(:, b), &
                               cs(:, b), wi, x, y, z, w)
      enddo
      if (nb > 1_I4P) then
         do b = 1_I4P, nb
            if (nfs(b) < 1_I4P) cycle
            call clear_body_3d_dev(surfaces(b), b, posed, prs(:, :, b), pts(:, b), w)
            if (m > 0_I4P) call blocked_3d(.true., surfaces(b), b, n, m, posed, prs(:, :, b), pts(:, b), nb, tl, w%body, &
                                           w%foot, w%probe, w%visible)
         enddo
      endif
      call bodies_flags(.true., 3_I4P, n, ni, nj, nk, w%ijk, ncnt, w%clearance, sd, w%shared, w%overlap)
   endif
   !$acc end data
   call finish_device(w, res, cl)
   endsubroutine bodies_3d_device

   subroutine nodes_body_3d_dev(surface, posed, pr, pt, x, y, z, nd, ncnt)
   !< One body's node states on the device (on its refitted tree if posed, built on the device), added to the counts.
   type(surface_stl_object), intent(in), target :: surface          !< Surface (mapped).
   logical,                  intent(in)         :: posed            !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   real(R8P),                intent(in)         :: x(:), y(:), z(:) !< Point lines (present).
   logical,                  intent(inout)      :: nd(size(x), size(y), size(z))   !< Point inside flags (present).
   integer(I4P),             intent(inout)      :: ncnt(size(x), size(y), size(z)) !< Bodies per point (present).
   type(aabb_node_object),       pointer :: gnode(:)               !< Geometry tree.
   type(facet_distance_payload), pointer :: gpay(:)                !< Geometry payload.
   real(R8P),                    pointer :: gtv(:,:)               !< Geometry vertices.
   type(aabb_node_object), allocatable, target :: nodem(:)         !< Refitted tree (posed).
   real(R8P),              allocatable, target :: tvm(:,:)         !< Moved vertices (posed).
   integer(I4P) :: i, j, k, ni, nj, nk                             !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nk = size(z, kind=I4P)
   gnode => surface%aabb%node ; gpay => surface%aabb%payload ; gtv => surface%tri_vertices
   if (posed) then
      allocate(nodem(0:ubound(surface%aabb%node, 1)), tvm(9, surface%get_facets_number()))
      !$acc enter data create(nodem, tvm)
      call posed_tree_3d(surface, pr, pt, .true., nodem, tvm)
      gnode => nodem ; gtv => tvm
   endif
   !$acc data copyin(gnode, gpay, gtv)
   call cut_cells_3d_nodes_tree_present(gnode, gpay, gtv, x, y, z, nd)
   !$acc parallel loop gang vector collapse(3) default(present) present(nd, ncnt)
   do k = 1_I4P, nk
      do j = 1_I4P, nj
         do i = 1_I4P, ni
            if (nd(i, j, k)) ncnt(i, j, k) = ncnt(i, j, k) + 1_I4P
         enddo
      enddo
   enddo
   !$acc end data
   if (posed) then
      !$acc exit data delete(nodem, tvm)
   endif
   endsubroutine nodes_body_3d_dev

   subroutine fill_body_3d_dev(surface, b, posed, pr, pt, tol, wc, wv, lv, av, c, wi, x, y, z, w)
   !< One body's table at the ghosts of w on the device (`fill_3d_dev` on a work table whose ghosts are copied on the
   !< device; the refitted tree built on the device if posed), merged into w (`merge_table`), the work table released.
   type(surface_stl_object), intent(in), target :: surface          !< Surface (mapped).
   integer(I4P),             intent(in)         :: b                !< Body index.
   logical,                  intent(in)         :: posed            !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   real(R8P),                intent(in)         :: tol              !< Visibility tolerance.
   logical,                  intent(in)         :: wc, wv, wi       !< Curvature, velocity, ILW.
   real(R8P),                intent(in)         :: lv(3), av(3), c(3) !< Rigid velocity.
   real(R8P),                intent(in)         :: x(:), y(:), z(:) !< Point lines (present).
   type(ib_closure_t),       intent(inout)      :: w                !< Table (mapped, ghosts filled).
   type(ib_closure_t)                    :: tb                     !< The body's work table.
   type(aabb_node_object),       pointer :: gnode(:)               !< Geometry tree.
   type(facet_distance_payload), pointer :: gpay(:)                !< Geometry payload.
   real(R8P),                    pointer :: gtv(:,:)               !< Geometry vertices.
   type(aabb_node_object), allocatable, target :: nodem(:)         !< Refitted tree (posed).
   real(R8P),              allocatable, target :: tvm(:,:)         !< Moved vertices (posed).
   integer(I4P),           allocatable         :: fpatch(:)        !< Facet patches.
   integer(I4P) :: nf, f                                           !< Work.

   nf = surface%get_facets_number()
   allocate(fpatch(nf))
   do f = 1_I4P, nf
      fpatch(f) = surface%facet_patch(f)
   enddo
   gnode => surface%aabb%node ; gpay => surface%aabb%payload ; gtv => surface%tri_vertices
   if (posed) then
      allocate(nodem(0:ubound(surface%aabb%node, 1)), tvm(9, nf))
      !$acc enter data create(nodem, tvm)
      call posed_tree_3d(surface, pr, pt, .true., nodem, tvm)
      gnode => nodem ; gtv => tvm
   endif
   !$acc data copyin(gnode, gpay, gtv, fpatch)
   call closure_arrays(tb, 3_I4P, w%n, w%m, wc, wv, wi)
   tb%ni = w%ni ; tb%nj = w%nj ; tb%nk = w%nk ; tb%reach = w%reach ; tb%lengths = w%lengths
   call map_create(tb)
   call copy_ghosts(.true., 3_I4P, w%n, w%ijk, w%points, w%mask, w%steps, tb%ijk, tb%points, tb%mask, tb%steps)
   call fill_3d_dev(surface, gnode, gpay, gtv, fpatch, tol, posed, pr, pt, wc, wv, lv, av, c, wi, x, y, z, tb)
   call merge_table(.true., b, tb, w)
   call map_delete(tb)
   !$acc end data
   if (posed) then
      !$acc exit data delete(nodem, tvm)
   endif
   endsubroutine fill_body_3d_dev

   subroutine clear_body_3d_dev(surface, b, posed, pr, pt, w)
   !< The clearance of every ghost's foot to body b on the device: the feet in b's frame, `query_points_present`, the
   !< closest points back in the world, `merge_clear` (the ghosts body b owns are skipped).
   type(surface_stl_object), intent(in), target :: surface         !< Surface (mapped).
   integer(I4P),             intent(in)         :: b               !< Body index.
   logical,                  intent(in)         :: posed           !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3) !< Pose.
   type(ib_closure_t),       intent(inout)      :: w               !< Table (mapped).
   real(R8P),    allocatable :: qp(:,:), qd(:), qc(:,:), qn(:,:), qw(:,:) !< Work.
   integer(I4P), allocatable :: qf(:), qft(:), qfi(:), qb(:)        !< Work.
   integer(I4P) :: n                                                !< Ghosts.

   n = w%n
   allocate(qp(3, n), qd(n), qc(3, n), qn(3, n), qw(3, n), qf(n), qft(n), qfi(n), qb(n))
   !$acc data create(qp, qd, qc, qn, qw, qf, qft, qfi, qb)
   call body_points_3d_dev(n, w%foot, posed, pr, pt, qp)
   call surface%query_points_present(n, qp, .true., .false., qd, qc, qn, qf, qft, qfi, qb)
   call world_points_3d_dev(n, posed, pr, pt, qc, qw)
   call merge_clear(.true., 3_I4P, n, b, qd, qw, w%body, w%clearance, w%second_body, w%second_foot)
   !$acc end data
   endsubroutine clear_body_3d_dev

   subroutine bodies_2d_host(curves, nvs, pc, x, y, reach, ll, posed, prs, pts, wc, wv, lvs, avs, cs, wi, sd, cl, poses)
   !< The host path of `fd_closure_2d_bodies` (see `bodies_3d_host`): node flags per curve on its placed vertices, the
   !< union, the ghosts, per curve `fill_2d_host` merged in curve order, the clearances, the visibility, the flags.
   type(curve_object),    intent(in), target   :: curves(:)            !< Curves.
   integer(I4P),          intent(in)           :: nvs(:)               !< Vertices per curve.
   type(placed_curve_t),  intent(in)           :: pc(:)                !< Placed vertices.
   real(R8P),             intent(in)           :: x(:), y(:)           !< Point lines.
   integer(I4P),          intent(in)           :: reach(:)             !< Reach.
   real(R8P),             intent(in)           :: ll(:)                !< Probe lengths.
   logical,               intent(in)           :: posed                !< Poses given.
   real(R8P),             intent(in)           :: prs(:,:,:), pts(:,:) !< Poses [2, 2, nb], [2, nb].
   logical,               intent(in)           :: wc, wv, wi           !< Curvature, velocity, ILW.
   real(R8P),             intent(in)           :: lvs(:,:), avs(:,:), cs(:,:) !< Rigid velocities [3, nb].
   real(R8P),             intent(in)           :: sd                   !< Shared distance.
   type(ib_closure_t),    intent(inout)        :: cl                   !< Table.
   type(rigid_pose_2d_t), intent(in), optional :: poses(:)             !< Poses.
   type(ghost_points_t)      :: gp                                     !< Ghosts.
   type(ib_closure_t)        :: tb                                     !< One curve's table at the ghosts.
   logical,      allocatable :: nd(:,:)                                !< Point inside flags.
   integer(I4P), allocatable :: ncnt(:,:,:)                            !< Curves containing each point.
   real(R8P),    allocatable :: qd(:), qc(:,:)                         !< Clearance queries.
   integer(I4P)              :: nb, b, n, m, ni, nj                    !< Work.

   nb = size(nvs, kind=I4P) ; ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; m = size(ll, kind=I4P)
   allocate(ncnt(ni, nj, 1)) ; ncnt = 0_I4P
   do b = 1_I4P, nb
      if (nvs(b) < 1_I4P) cycle
      call node_flags(pc(b)%u, pc(b)%v, pc(b)%nxt, x, y, nd)
      where (nd) ncnt(:, :, 1) = ncnt(:, :, 1) + 1_I4P
      deallocate(nd)
   enddo
   allocate(nd(ni, nj)) ; nd = ncnt(:, :, 1) > 0_I4P
   call ghost_points_nodes_2d(nd, x, y, reach, gp)
   n = gp%n
   call closure_arrays(cl, 2_I4P, n, m, wc, wv, wi) ; call closure_arrays_bodies(cl, nb, sd)
   cl%ni = ni ; cl%nj = nj ; cl%nk = 1_I4P ; cl%reach(1:2) = reach(1:2) ; cl%lengths = ll
   cl%overlap_points = count(ncnt >= 2_I4P, kind=I4P)
   if (n < 1_I4P) return
   cl%ijk = gp%ijk ; cl%points = gp%points ; cl%mask = gp%mask ; cl%steps = gp%steps
   call bodies_init(.false., 2_I4P, n, size(cl%ilw_body, 2, kind=I4P), size(cl%ilw_kappa, 3, kind=I4P), cl%body, &
                    cl%distance, cl%clearance, cl%second_body, cl%second_foot, cl%ilw_t, cl%ilw_element, cl%ilw_normal, &
                    cl%ilw_body, cl%ilw_kappa, cl%ilw_direction)
   do b = 1_I4P, nb
      if (nvs(b) < 1_I4P) cycle
      call closure_arrays(tb, 2_I4P, n, m, wc, wv, wi)
      tb%ni = ni ; tb%nj = nj ; tb%nk = 1_I4P ; tb%reach(1:2) = reach(1:2) ; tb%lengths = ll
      tb%ijk = cl%ijk ; tb%points = cl%points ; tb%mask = cl%mask ; tb%steps = cl%steps
      if (posed) then
         call fill_2d_host(curves(b), pc(b)%u, pc(b)%v, pc(b)%nxt, x, y, ll, poses(b), wc, wv, lvs(:, b), avs(3, b), &
                           cs(:, b), wi, tb)
      else
         call fill_2d_host(curves(b), pc(b)%u, pc(b)%v, pc(b)%nxt, x, y, ll, wc=wc, wv=wv, lv=lvs(:, b), aw=avs(3, b), &
                           c=cs(:, b), wi=wi, cl=tb)
      endif
      call merge_table(.false., b, tb, cl)
      call ib_closure_reset(tb)
   enddo
   if (nb > 1_I4P) then
      allocate(qd(n), qc(2, n))
      do b = 1_I4P, nb
         if (nvs(b) < 1_I4P) cycle
         if (posed) then
            call curves(b)%query_many(points=cl%foot, distances=qd, closest_points=qc, pose=poses(b))
         else
            call curves(b)%query_many(points=cl%foot, distances=qd, closest_points=qc)
         endif
         call merge_clear(.false., 2_I4P, n, b, qd, qc, cl%body, cl%clearance, cl%second_body, cl%second_foot)
         if (m > 0_I4P) call blocked_2d(.false., curves(b), b, n, m, posed, prs(:, :, b), pts(:, b), cl%body, cl%foot, &
                                        cl%probe, cl%visible)
      enddo
   endif
   call bodies_flags(.false., 2_I4P, n, ni, nj, 1_I4P, cl%ijk, ncnt, cl%clearance, sd, cl%shared, cl%overlap)
   endsubroutine bodies_2d_host

   subroutine bodies_2d_device(curves, nvs, pc, x, y, reach, ll, posed, prs, pts, wc, wv, lvs, avs, cs, wi, sd, res, cl)
   !< The device path of `fd_closure_2d_bodies` (see `bodies_3d_device`).
   type(curve_object),   intent(in), target :: curves(:)            !< Curves (mapped).
   integer(I4P),         intent(in)         :: nvs(:)               !< Vertices per curve.
   type(placed_curve_t), intent(in)         :: pc(:)                !< Placed vertices.
   real(R8P),            intent(in)         :: x(:), y(:)           !< Point lines.
   integer(I4P),         intent(in)         :: reach(:)             !< Reach.
   real(R8P),            intent(in)         :: ll(:)                !< Probe lengths.
   logical,              intent(in)         :: posed                !< Poses given.
   real(R8P),            intent(in)         :: prs(:,:,:), pts(:,:) !< Poses.
   logical,              intent(in)         :: wc, wv, wi, res      !< Curvature, velocity, ILW, resident.
   real(R8P),            intent(in)         :: lvs(:,:), avs(:,:), cs(:,:) !< Rigid velocities.
   real(R8P),            intent(in)         :: sd                   !< Shared distance.
   type(ib_closure_t),   intent(inout)      :: cl                   !< Table.
   type(ib_closure_t)        :: w                                   !< The table built on the device.
   logical,      allocatable :: nd(:,:,:)                           !< Point inside flags.
   integer(I4P), allocatable :: ncnt(:,:,:), cnt(:), off(:)         !< Curves per point, ghost counts, offsets.
   real(R8P)                 :: z(1)                                !< No z line.
   integer(I4P)              :: nb, b, ni, nj, np, n, m, r(3), op   !< Work.

   nb = size(nvs, kind=I4P) ; ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; np = ni * nj ; z = 0._R8P
   m = size(ll, kind=I4P) ; r(1:2) = reach(1:2) ; r(3) = 0_I4P
   allocate(nd(ni, nj, 1), ncnt(ni, nj, 1), cnt(np), off(np))
   !$acc data copyin(x, y, z) create(nd, ncnt, cnt, off)
   call zero_counts_dev(np, ncnt)
   do b = 1_I4P, nb
      if (nvs(b) < 1_I4P) cycle
      call nodes_body_2d_dev(pc(b)%u, pc(b)%v, pc(b)%nxt, x, y, ncnt)
   enddo
   call union_dev(np, ncnt, nd, op)
   call ghost_count_dev(2_I4P, ni, nj, 1_I4P, nd, r, cnt)
   call links_scan_dev(np, cnt, off, n)
   call closure_arrays(w, 2_I4P, n, m, wc, wv, wi) ; call closure_arrays_bodies(w, nb, sd)
   w%ni = ni ; w%nj = nj ; w%nk = 1_I4P ; w%reach = r ; w%lengths = ll ; w%overlap_points = op
   call map_create(w)
   if (n > 0_I4P) then
      call ghost_fill_dev(2_I4P, ni, nj, 1_I4P, nd, x, y, z, r, cnt, off, n, w%ijk, w%points, w%mask, w%steps)
      call bodies_init(.true., 2_I4P, n, size(w%ilw_body, 2, kind=I4P), size(w%ilw_kappa, 3, kind=I4P), w%body, &
                       w%distance, w%clearance, w%second_body, w%second_foot, w%ilw_t, w%ilw_element, w%ilw_normal, &
                       w%ilw_body, w%ilw_kappa, w%ilw_direction)
      do b = 1_I4P, nb
         if (nvs(b) < 1_I4P) cycle
         call fill_body_2d_dev(curves(b), b, pc(b), posed, prs(:, :, b), pts(:, b), wc, wv, lvs(:, b), avs(3, b), &
                               cs(:, b), wi, x, y, w)
      enddo
      if (nb > 1_I4P) then
         do b = 1_I4P, nb
            if (nvs(b) < 1_I4P) cycle
            call clear_body_2d_dev(curves(b), b, posed, prs(:, :, b), pts(:, b), w)
            if (m > 0_I4P) call blocked_2d(.true., curves(b), b, n, m, posed, prs(:, :, b), pts(:, b), w%body, w%foot, &
                                           w%probe, w%visible)
         enddo
      endif
      call bodies_flags(.true., 2_I4P, n, ni, nj, 1_I4P, w%ijk, ncnt, w%clearance, sd, w%shared, w%overlap)
   endif
   !$acc end data
   call finish_device(w, res, cl)
   endsubroutine bodies_2d_device

   subroutine nodes_body_2d_dev(u, v, nxt, x, y, ncnt)
   !< One curve's node flags on the device (`node_parity` per point over the rows of `node_rows`), added to the counts.
   real(R8P),    intent(in)    :: u(:), v(:)  !< Vertices (placed).
   integer(I4P), intent(in)    :: nxt(:)      !< Next vertex.
   real(R8P),    intent(in)    :: x(:), y(:)  !< Point lines (present).
   integer(I4P), intent(inout) :: ncnt(size(x), size(y), 1) !< Curves per point (present).
   integer(I4P), allocatable :: rstart(:), rseg(:) !< Rows.
   integer(I4P) :: i, j, ni, nj, nu, nr      !< Work.

   ni = size(x, kind=I4P) ; nj = size(y, kind=I4P) ; nu = size(u, kind=I4P)
   call node_rows(u, v, nxt, y, rstart, rseg)
   nr = size(rseg, kind=I4P)
   !$acc data copyin(u, v, nxt, rstart, rseg)
   !$acc parallel loop gang vector collapse(2) default(present) present(u, v, nxt, rstart, rseg, x, y, ncnt)
   do j = 1_I4P, nj
      do i = 1_I4P, ni
         if (node_parity(x(i), y(j), rstart(j), rstart(j + 1) - 1_I4P, nr, rseg, nu, u, v, nxt)) &
            ncnt(i, j, 1) = ncnt(i, j, 1) + 1_I4P
      enddo
   enddo
   !$acc end data
   endsubroutine nodes_body_2d_dev

   subroutine fill_body_2d_dev(curve, b, pc, posed, pr, pt, wc, wv, lv, aw, c, wi, x, y, w)
   !< 2-D twin of `fill_body_3d_dev`.
   type(curve_object),   intent(in), target :: curve           !< Curve (mapped).
   integer(I4P),         intent(in)         :: b               !< Curve index.
   type(placed_curve_t), intent(in)         :: pc              !< Placed vertices.
   logical,              intent(in)         :: posed           !< A pose.
   real(R8P),            intent(in)         :: pr(2, 2), pt(2) !< Pose.
   logical,              intent(in)         :: wc, wv, wi      !< Curvature, velocity, ILW.
   real(R8P),            intent(in)         :: lv(3), aw, c(3) !< Rigid velocity (angular: scalar).
   real(R8P),            intent(in)         :: x(:), y(:)      !< Point lines (present).
   type(ib_closure_t),   intent(inout)      :: w               !< Table (mapped, ghosts filled).
   type(ib_closure_t) :: tb                                    !< The curve's work table.

   !$acc data copyin(pc%u, pc%v, pc%nxt, pc%spatch)
   call closure_arrays(tb, 2_I4P, w%n, w%m, wc, wv, wi)
   tb%ni = w%ni ; tb%nj = w%nj ; tb%nk = w%nk ; tb%reach = w%reach ; tb%lengths = w%lengths
   call map_create(tb)
   call copy_ghosts(.true., 2_I4P, w%n, w%ijk, w%points, w%mask, w%steps, tb%ijk, tb%points, tb%mask, tb%steps)
   call fill_2d_dev(curve, pc%u, pc%v, pc%nxt, pc%spatch, posed, pr, pt, wc, wv, lv, aw, c, wi, x, y, tb)
   call merge_table(.true., b, tb, w)
   call map_delete(tb)
   !$acc end data
   endsubroutine fill_body_2d_dev

   subroutine clear_body_2d_dev(curve, b, posed, pr, pt, w)
   !< 2-D twin of `clear_body_3d_dev`.
   type(curve_object), intent(in), target :: curve           !< Curve (mapped).
   integer(I4P),       intent(in)         :: b               !< Curve index.
   logical,            intent(in)         :: posed           !< A pose.
   real(R8P),          intent(in)         :: pr(2, 2), pt(2) !< Pose.
   type(ib_closure_t), intent(inout)      :: w               !< Table (mapped).
   real(R8P),    allocatable :: qp(:,:), qd(:), qc(:,:), qn(:,:), qw(:,:) !< Work.
   integer(I4P), allocatable :: qs(:), qv(:), qb(:)                 !< Work.
   integer(I4P) :: n                                                !< Ghosts.

   n = w%n
   allocate(qp(2, n), qd(n), qc(2, n), qn(2, n), qw(2, n), qs(n), qv(n), qb(n))
   !$acc data create(qp, qd, qc, qn, qw, qs, qv, qb)
   call body_points_2d_dev(n, w%foot, posed, pr, pt, qp)
   call curve%query_present(n, qp, .true., .false., qd, qc, qs, qv, qn, qb)
   call world_points_2d_dev(n, posed, pr, pt, qc, qw)
   call merge_clear(.true., 2_I4P, n, b, qd, qw, w%body, w%clearance, w%second_body, w%second_foot)
   !$acc end data
   endsubroutine clear_body_2d_dev

   subroutine merge_table(dev, b, tb, w)
   !< Merge body b's table tb at the ghosts of w into w: `merge_ghosts`, and `merge_ilw` with ILW (arrays present if dev).
   logical,            intent(in)    :: dev !< On the device.
   integer(I4P),       intent(in)    :: b   !< Body index.
   type(ib_closure_t), intent(in)    :: tb  !< Body b's table.
   type(ib_closure_t), intent(inout) :: w   !< Table.

   call merge_ghosts(dev, w%dims, w%n, w%m, size(w%kappa, 2, kind=I4P), size(w%velocity, 2, kind=I4P), b, tb%distance, &
                     tb%foot, tb%normal, tb%element, tb%feature, tb%body, tb%patch, tb%probe, tb%visible, tb%probe_cell, &
                     tb%kappa, tb%direction, tb%velocity, w%distance, w%foot, w%normal, w%element, w%feature, w%body, &
                     w%sub_body, w%patch, w%probe, w%visible, w%probe_cell, w%kappa, w%direction, w%velocity)
   if (w%with_ilw) call merge_ilw(dev, w%dims, w%n, size(w%ilw_kappa, 3, kind=I4P), b, tb%ilw_t, tb%ilw_element, &
                                  tb%ilw_normal, tb%ilw_kappa, tb%ilw_direction, w%ilw_t, w%ilw_element, w%ilw_normal, &
                                  w%ilw_kappa, w%ilw_direction, w%ilw_body)
   endsubroutine merge_table

   subroutine merge_ghosts(dev, dims, n, m, nc, nvl, b, td, tf, tn, te, tft, tbd, tp, tpr, tvi, tpc, tk, tdir, tvel, &
                           d, f, nrm, e, ft, own, sub, p, pr, vi, pc, k, dir, vel)
   !< Per ghost: body b becomes the owner when the ghost has none yet or |its distance| is below the owner's (strict: ties
   !< keep the lower body index); then every per-ghost field of the owner is copied from b's table (body = b, sub_body =
   !< b's facet_body). Copies and comparisons only, so the host and device merges are equal bit for bit. One kernel per
   !< non-empty group (a zero-size array is never mapped), the owners and distances last: every kernel decides on the
   !< same values.
   logical,      intent(in)    :: dev                    !< On the device (arrays present).
   integer(I4P), intent(in)    :: dims, n, m, nc, nvl, b !< Dimensions, ghosts, probes, curvature and velocity sizes, body.
   real(R8P),    intent(in)    :: td(n), tf(dims, n), tn(dims, n) !< b: distances, feet, normals.
   integer(I4P), intent(in)    :: te(n), tft(n), tbd(n), tp(n)    !< b: elements, features, facet_bodies, patches.
   real(R8P),    intent(in)    :: tpr(dims, m, n)        !< b: probes.
   logical,      intent(in)    :: tvi(m, n)              !< b: visibility.
   integer(I4P), intent(in)    :: tpc(dims, m, n)        !< b: probe cells.
   real(R8P),    intent(in)    :: tk(dims - 1, nc), tdir(dims, dims - 1, nc), tvel(dims, nvl) !< b: curvature, velocity.
   real(R8P),    intent(inout) :: d(n), f(dims, n), nrm(dims, n) !< Distances, feet, normals.
   integer(I4P), intent(inout) :: e(n), ft(n), own(n), sub(n), p(n) !< Elements, features, owners, sub-bodies, patches.
   real(R8P),    intent(inout) :: pr(dims, m, n)         !< Probes.
   logical,      intent(inout) :: vi(m, n)               !< Visibility.
   integer(I4P), intent(inout) :: pc(dims, m, n)         !< Probe cells.
   real(R8P),    intent(inout) :: k(dims - 1, nc), dir(dims, dims - 1, nc), vel(dims, nvl) !< Curvature, velocity.
   integer(I4P) :: g, a, kk, q                           !< Counters.

   if (m > 0_I4P) then
      !$acc parallel loop gang vector default(present) present(td, own, d, tpr, tvi, tpc, pr, vi, pc) private(a, kk) if(dev)
      do g = 1_I4P, n
         if (own(g) == 0_I4P .or. abs(td(g)) < abs(d(g))) then
            do kk = 1_I4P, m
               vi(kk, g) = tvi(kk, g)
               do a = 1_I4P, dims
                  pr(a, kk, g) = tpr(a, kk, g) ; pc(a, kk, g) = tpc(a, kk, g)
               enddo
            enddo
         endif
      enddo
   endif
   if (nc > 0_I4P) then
      !$acc parallel loop gang vector default(present) present(td, own, d, tk, tdir, k, dir) private(a, q) if(dev)
      do g = 1_I4P, n
         if (own(g) == 0_I4P .or. abs(td(g)) < abs(d(g))) then
            do q = 1_I4P, dims - 1_I4P
               k(q, g) = tk(q, g)
               do a = 1_I4P, dims
                  dir(a, q, g) = tdir(a, q, g)
               enddo
            enddo
         endif
      enddo
   endif
   if (nvl > 0_I4P) then
      !$acc parallel loop gang vector default(present) present(td, own, d, tvel, vel) private(a) if(dev)
      do g = 1_I4P, n
         if (own(g) == 0_I4P .or. abs(td(g)) < abs(d(g))) then
            do a = 1_I4P, dims
               vel(a, g) = tvel(a, g)
            enddo
         endif
      enddo
   endif
   !$acc parallel loop gang vector default(present) present(td, tf, tn, te, tft, tbd, tp, d, f, nrm, e, ft, own, sub, p) &
   !$acc&  private(a) if(dev)
   do g = 1_I4P, n
      if (own(g) == 0_I4P .or. abs(td(g)) < abs(d(g))) then
         own(g) = b ; d(g) = td(g) ; e(g) = te(g) ; ft(g) = tft(g) ; sub(g) = tbd(g) ; p(g) = tp(g)
         do a = 1_I4P, dims
            f(a, g) = tf(a, g) ; nrm(a, g) = tn(a, g)
         enddo
      endif
   enddo
   endsubroutine merge_ghosts

   subroutine merge_ilw(dev, dims, n, nc, b, tt, tel, tnrm, tk, tdir, it, iel, inrm, ik, idir, ib)
   !< Per ILW slot: body b's crossing replaces the slot's when the slot has none yet or b's lies strictly nearer the fluid
   !< point (a larger coordinate on a + slot, a smaller one on a - slot; ties keep the lower body index); its element,
   !< normal and curvature are copied, ilw_body = b. The curvature kernel first (it decides on the unchanged slots).
   logical,      intent(in)    :: dev                     !< On the device (arrays present).
   integer(I4P), intent(in)    :: dims, n, nc, b          !< Dimensions, ghosts, curvature size, body.
   real(R8P),    intent(in)    :: tt(2 * dims, n)         !< b: crossing coordinates.
   integer(I4P), intent(in)    :: tel(2 * dims, n)        !< b: crossing elements.
   real(R8P),    intent(in)    :: tnrm(dims, 2 * dims, n) !< b: normals.
   real(R8P),    intent(in)    :: tk(dims - 1, 2 * dims, nc), tdir(dims, dims - 1, 2 * dims, nc) !< b: curvature.
   real(R8P),    intent(inout) :: it(2 * dims, n)         !< Crossing coordinates.
   integer(I4P), intent(inout) :: iel(2 * dims, n)        !< Crossing elements.
   real(R8P),    intent(inout) :: inrm(dims, 2 * dims, n) !< Normals.
   real(R8P),    intent(inout) :: ik(dims - 1, 2 * dims, nc), idir(dims, dims - 1, 2 * dims, nc) !< Curvature.
   integer(I4P), intent(inout) :: ib(2 * dims, n)         !< Crossing bodies.
   integer(I4P) :: g, l, a, q                             !< Counters.

   if (nc > 0_I4P) then
      !$acc parallel loop gang vector collapse(2) default(present) present(tt, tel, it, ib, tk, tdir, ik, idir) &
      !$acc&  private(a, q) if(dev)
      do g = 1_I4P, n
         do l = 1_I4P, 2_I4P * dims
            if (ilw_take_oac(l, tel(l, g), ib(l, g), tt(l, g), it(l, g))) then
               do q = 1_I4P, dims - 1_I4P
                  ik(q, l, g) = tk(q, l, g)
                  do a = 1_I4P, dims
                     idir(a, q, l, g) = tdir(a, q, l, g)
                  enddo
               enddo
            endif
         enddo
      enddo
   endif
   !$acc parallel loop gang vector collapse(2) default(present) present(tt, tel, tnrm, it, iel, inrm, ib) private(a) if(dev)
   do g = 1_I4P, n
      do l = 1_I4P, 2_I4P * dims
         if (ilw_take_oac(l, tel(l, g), ib(l, g), tt(l, g), it(l, g))) then
            ib(l, g) = b ; it(l, g) = tt(l, g) ; iel(l, g) = tel(l, g)
            do a = 1_I4P, dims
               inrm(a, l, g) = tnrm(a, l, g)
            enddo
         endif
      enddo
   enddo
   endsubroutine merge_ilw

   pure logical function ilw_take_oac(l, tel, ib, tt, it) result(take)
   !< Kernel: a body's crossing (element tel, coordinate tt) replaces the crossing of ILW slot l (body ib, coordinate it).
   !$acc routine seq
   integer(I4P), intent(in) :: l, tel, ib !< Slot, candidate element, current body.
   real(R8P),    intent(in) :: tt, it     !< Candidate and current coordinates.

   take = .false.
   if (tel < 1_I4P) return
   if (ib == 0_I4P) then
      take = .true.
   elseif (mod(l, 2_I4P) == 0_I4P) then
      take = tt > it
   else
      take = tt < it
   endif
   endfunction ilw_take_oac

   subroutine merge_clear(dev, dims, n, b, qd, qc, own, clr, sb, sf)
   !< Per ghost not owned by body b: |qd| (the unsigned distance of the owner's foot to b) replaces the clearance when
   !< the ghost has no second body yet or it is strictly smaller (ties keep the lower body index); second_foot = qc.
   logical,      intent(in)    :: dev          !< On the device (arrays present).
   integer(I4P), intent(in)    :: dims, n, b   !< Dimensions, ghosts, body.
   real(R8P),    intent(in)    :: qd(n)        !< Distances of the feet to b.
   real(R8P),    intent(in)    :: qc(dims, n)  !< Closest points of b (world).
   integer(I4P), intent(in)    :: own(n)       !< Owners.
   real(R8P),    intent(inout) :: clr(n)       !< Clearances.
   integer(I4P), intent(inout) :: sb(n)        !< Second bodies.
   real(R8P),    intent(inout) :: sf(dims, n)  !< Second feet.
   real(R8P)    :: dd                          !< Distance.
   integer(I4P) :: g, a                        !< Counters.

   !$acc parallel loop gang vector default(present) present(qd, qc, own, clr, sb, sf) private(dd, a) if(dev)
   do g = 1_I4P, n
      if (own(g) /= b) then
         dd = abs(qd(g))
         if (sb(g) == 0_I4P .or. dd < clr(g)) then
            clr(g) = dd ; sb(g) = b
            do a = 1_I4P, dims
               sf(a, g) = qc(a, g)
            enddo
         endif
      endif
   enddo
   endsubroutine merge_clear

   subroutine blocked_3d(dev, surface, b, n, m, posed, pr, pt, nb, tols, own, foot, probe, vis)
   !< The probes' visibility against body b (B4): for every probe of a ghost body b does not own, still visible, the open
   !< segment foot - probe (world) moved to b's frame (`to_body`) is tested against b's SAH BVH (`bvh_segment_blocked_oac`,
   !< the kernel of `probe_many`, with the owner's tolerance); blocked: visible = .false.
   logical,                  intent(in)         :: dev              !< On the device (arrays present).
   type(surface_stl_object), intent(in), target :: surface          !< Body b (mapped if dev).
   integer(I4P),             intent(in)         :: b, n, m, nb      !< Body, ghosts, probes, bodies.
   logical,                  intent(in)         :: posed            !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< b's pose.
   real(R8P),                intent(in)         :: tols(nb)         !< Visibility tolerances (present if dev).
   integer(I4P),             intent(in)         :: own(n)           !< Owners (present if dev).
   real(R8P),                intent(in)         :: foot(3, n)       !< Feet (present if dev).
   real(R8P),                intent(in)         :: probe(3, m, n)   !< Probes (present if dev).
   logical,                  intent(inout)      :: vis(m, n)        !< Visibility (present if dev).
   type(aabb_node_object),       pointer :: node_p(:)              !< Alias of surface%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:)           !< Alias of surface%aabb%payload.
   real(R8P),                    pointer :: tv_p(:,:)              !< Alias of surface%tri_vertices.
   real(R8P)    :: r(3, 3), t(3), ax, ay, az, bx, by, bz           !< Work.
   integer(I4P) :: g, k, o                                         !< Work.

   r = pr ; t = pt
   node_p => surface%aabb%node ; payload_p => surface%aabb%payload ; tv_p => surface%tri_vertices
   !$acc parallel loop gang vector collapse(2) default(present) present(node_p, payload_p, tv_p, tols, own, foot, probe, &
   !$acc&  vis) copyin(r, t) private(ax, ay, az, bx, by, bz, o) if(dev)
   do g = 1_I4P, n
      do k = 1_I4P, m
         o = own(g)
         if (o /= b .and. vis(k, g)) then
            if (posed) then
               call to_body_3d_oac(r, t, foot(1, g), foot(2, g), foot(3, g), ax, ay, az)
               call to_body_3d_oac(r, t, probe(1, k, g), probe(2, k, g), probe(3, k, g), bx, by, bz)
            else
               ax = foot(1, g) ; ay = foot(2, g) ; az = foot(3, g)
               bx = probe(1, k, g) ; by = probe(2, k, g) ; bz = probe(3, k, g)
            endif
            if (bvh_segment_blocked_oac(node_p, payload_p, tv_p, ax, ay, az, bx, by, bz, tols(o))) vis(k, g) = .false.
         endif
      enddo
   enddo
   endsubroutine blocked_3d

   subroutine blocked_2d(dev, curve, b, n, m, posed, pr, pt, own, foot, probe, vis)
   !< 2-D twin of `blocked_3d`: the closed segment foot - probe in curve b's frame against all its segments
   !< (`curve%segments_blocked_present`).
   logical,            intent(in)         :: dev             !< On the device (arrays present).
   type(curve_object), intent(in), target :: curve           !< Curve b (mapped if dev).
   integer(I4P),       intent(in)         :: b, n, m         !< Curve, ghosts, probes.
   logical,            intent(in)         :: posed           !< A pose.
   real(R8P),          intent(in)         :: pr(2, 2), pt(2) !< b's pose.
   integer(I4P),       intent(in)         :: own(n)          !< Owners (present if dev).
   real(R8P),          intent(in)         :: foot(2, n)      !< Feet (present if dev).
   real(R8P),          intent(in)         :: probe(2, m, n)  !< Probes (present if dev).
   logical,            intent(inout)      :: vis(m, n)       !< Visibility (present if dev).
   real(R8P), allocatable :: au(:), av(:), bu(:), bv(:)      !< Segments in b's frame.
   logical,   allocatable :: act(:), blk(:)                  !< Tested, blocked.
   real(R8P)    :: r(2, 2), t(2)                             !< Pose copies.
   integer(I4P) :: g, k, l, nm                               !< Work.

   r = pr ; t = pt ; nm = n * m
   allocate(au(nm), av(nm), bu(nm), bv(nm), act(nm), blk(nm))
   !$acc data create(au, av, bu, bv, act, blk) copyin(r, t) if(dev)
   !$acc parallel loop gang vector collapse(2) default(present) present(own, foot, probe, vis, au, av, bu, bv, act, r, t) &
   !$acc&  private(l) if(dev)
   do g = 1_I4P, n
      do k = 1_I4P, m
         l = k + m * (g - 1_I4P)
         act(l) = own(g) /= b .and. vis(k, g)
         if (posed) then
            call to_body_2d_oac(r, t, foot(1, g), foot(2, g), au(l), av(l))
            call to_body_2d_oac(r, t, probe(1, k, g), probe(2, k, g), bu(l), bv(l))
         else
            au(l) = foot(1, g) ; av(l) = foot(2, g) ; bu(l) = probe(1, k, g) ; bv(l) = probe(2, k, g)
         endif
      enddo
   enddo
   call curve%segments_blocked_present(nm, au, av, bu, bv, act, blk, dev)
   !$acc parallel loop gang vector collapse(2) default(present) present(vis, blk) private(l) if(dev)
   do g = 1_I4P, n
      do k = 1_I4P, m
         l = k + m * (g - 1_I4P)
         if (blk(l)) vis(k, g) = .false.
      enddo
   enddo
   !$acc end data
   endsubroutine blocked_2d

   subroutine bodies_init(dev, dims, n, nil, ncl, own, d, clr, sb, sf, it, iel, inrm, ib, ik, idir)
   !< Before the merges: no owner, distance 0, clearance huge, no second body; every ILW slot empty (zeros).
   logical,      intent(in)    :: dev                     !< On the device (arrays present).
   integer(I4P), intent(in)    :: dims, n, nil, ncl       !< Dimensions, ghosts, ILW and ILW curvature sizes.
   integer(I4P), intent(inout) :: own(n), sb(n)           !< Owners, second bodies.
   real(R8P),    intent(inout) :: d(n), clr(n), sf(dims, n) !< Distances, clearances, second feet.
   real(R8P),    intent(inout) :: it(2 * dims, nil)       !< Crossing coordinates.
   integer(I4P), intent(inout) :: iel(2 * dims, nil)      !< Crossing elements.
   real(R8P),    intent(inout) :: inrm(dims, 2 * dims, nil) !< Normals.
   integer(I4P), intent(inout) :: ib(2 * dims, nil)       !< Crossing bodies.
   real(R8P),    intent(inout) :: ik(dims - 1, 2 * dims, ncl), idir(dims, dims - 1, 2 * dims, ncl) !< Curvature.
   integer(I4P) :: g, l, a, q                             !< Counters.

   !$acc parallel loop gang vector default(present) present(own, sb, d, clr, sf) private(a) if(dev)
   do g = 1_I4P, n
      own(g) = 0_I4P ; sb(g) = 0_I4P ; d(g) = 0._R8P ; clr(g) = huge(1._R8P)
      do a = 1_I4P, dims
         sf(a, g) = 0._R8P
      enddo
   enddo
   if (nil > 0_I4P) then
      !$acc parallel loop gang vector collapse(2) default(present) present(it, iel, inrm, ib) private(a) if(dev)
      do g = 1_I4P, nil
         do l = 1_I4P, 2_I4P * dims
            it(l, g) = 0._R8P ; iel(l, g) = 0_I4P ; ib(l, g) = 0_I4P
            do a = 1_I4P, dims
               inrm(a, l, g) = 0._R8P
            enddo
         enddo
      enddo
   endif
   if (ncl > 0_I4P) then
      !$acc parallel loop gang vector collapse(2) default(present) present(ik, idir) private(a, q) if(dev)
      do g = 1_I4P, ncl
         do l = 1_I4P, 2_I4P * dims
            do q = 1_I4P, dims - 1_I4P
               ik(q, l, g) = 0._R8P
               do a = 1_I4P, dims
                  idir(a, q, l, g) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine bodies_init

   subroutine bodies_flags(dev, dims, n, ni, nj, nk, ijk, ncnt, clr, sd, shared, overlap)
   !< Per ghost: shared = clearance < sd, overlap = its point lies inside two or more bodies.
   logical,      intent(in)    :: dev                !< On the device (arrays present).
   integer(I4P), intent(in)    :: dims, n, ni, nj, nk !< Dimensions, ghosts, points per axis.
   integer(I4P), intent(in)    :: ijk(dims, n)       !< Ghost indexes.
   integer(I4P), intent(in)    :: ncnt(ni, nj, nk)   !< Bodies per point.
   real(R8P),    intent(in)    :: clr(n)             !< Clearances.
   real(R8P),    intent(in)    :: sd                 !< Shared distance.
   logical,      intent(inout) :: shared(n), overlap(n) !< Flags.
   integer(I4P) :: g, k                              !< Work.

   !$acc parallel loop gang vector default(present) present(ijk, ncnt, clr, shared, overlap) private(k) if(dev)
   do g = 1_I4P, n
      k = 1_I4P ; if (dims == 3_I4P) k = ijk(3, g)
      overlap(g) = ncnt(ijk(1, g), ijk(2, g), k) >= 2_I4P
      shared(g) = clr(g) < sd
   enddo
   endsubroutine bodies_flags

   subroutine copy_ghosts(dev, dims, n, ijk, pts, mask, steps, ijk_d, pts_d, mask_d, steps_d)
   !< Copy the ghosts of a table into a work table (arrays present if dev).
   logical,      intent(in)    :: dev                                   !< On the device.
   integer(I4P), intent(in)    :: dims, n                               !< Dimensions, ghosts.
   integer(I4P), intent(in)    :: ijk(dims, n), mask(n), steps(2 * dims, n) !< Source indexes, masks, steps.
   real(R8P),    intent(in)    :: pts(dims, n)                          !< Source points.
   integer(I4P), intent(inout) :: ijk_d(dims, n), mask_d(n), steps_d(2 * dims, n) !< Copies.
   real(R8P),    intent(inout) :: pts_d(dims, n)                        !< Copy.
   integer(I4P) :: g, a                                                 !< Counters.

   !$acc parallel loop gang vector default(present) present(ijk, pts, mask, steps, ijk_d, pts_d, mask_d, steps_d) &
   !$acc&  private(a) if(dev)
   do g = 1_I4P, n
      mask_d(g) = mask(g)
      do a = 1_I4P, dims
         ijk_d(a, g) = ijk(a, g) ; pts_d(a, g) = pts(a, g)
      enddo
      do a = 1_I4P, 2_I4P * dims
         steps_d(a, g) = steps(a, g)
      enddo
   enddo
   endsubroutine copy_ghosts

   subroutine zero_counts_dev(np, ncnt)
   !< Zero the bodies-per-point counts (present).
   integer(I4P), intent(in)    :: np       !< Points.
   integer(I4P), intent(inout) :: ncnt(np) !< Counts (present).
   integer(I4P) :: p                       !< Counter.

   !$acc parallel loop gang vector default(present) present(ncnt)
   do p = 1_I4P, np
      ncnt(p) = 0_I4P
   enddo
   endsubroutine zero_counts_dev

   subroutine union_dev(np, ncnt, nd, op)
   !< The union of the bodies' states (inside if inside one body) and the number of points inside two or more (present).
   integer(I4P), intent(in)    :: np       !< Points.
   integer(I4P), intent(in)    :: ncnt(np) !< Bodies per point (present).
   logical,      intent(inout) :: nd(np)   !< Point inside flags (present).
   integer(I4P), intent(out)   :: op       !< Points inside two or more bodies.
   integer(I4P) :: p, c                    !< Work.

   c = 0_I4P
   !$acc parallel loop gang vector default(present) present(ncnt, nd) reduction(+:c)
   do p = 1_I4P, np
      nd(p) = ncnt(p) > 0_I4P
      if (ncnt(p) >= 2_I4P) c = c + 1_I4P
   enddo
   op = c
   endsubroutine union_dev

   subroutine world_points_3d_dev(n, posed, pr, pt, qc, qw)
   !< Body-frame points to the world (`to_world`), or copied without a pose (present).
   integer(I4P), intent(in)    :: n               !< Points.
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(3, 3), pt(3) !< Pose.
   real(R8P),    intent(in)    :: qc(3, n)        !< Body points (present).
   real(R8P),    intent(inout) :: qw(3, n)        !< World points (present).
   real(R8P)    :: r(3, 3), t(3)                  !< Pose copies.
   integer(I4P) :: g                              !< Counter.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(qc, qw) copyin(r, t)
   do g = 1_I4P, n
      if (posed) then
         call to_world_3d_oac(r, t, qc(1, g), qc(2, g), qc(3, g), qw(1, g), qw(2, g), qw(3, g))
      else
         qw(1, g) = qc(1, g) ; qw(2, g) = qc(2, g) ; qw(3, g) = qc(3, g)
      endif
   enddo
   endsubroutine world_points_3d_dev

   subroutine world_points_2d_dev(n, posed, pr, pt, qc, qw)
   !< 2-D twin of `world_points_3d_dev`.
   integer(I4P), intent(in)    :: n               !< Points.
   logical,      intent(in)    :: posed           !< A pose.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),    intent(in)    :: qc(2, n)        !< Body points (present).
   real(R8P),    intent(inout) :: qw(2, n)        !< World points (present).
   real(R8P)    :: r(2, 2), t(2)                  !< Pose copies.
   integer(I4P) :: g                              !< Counter.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(qc, qw) copyin(r, t)
   do g = 1_I4P, n
      if (posed) then
         call to_world_2d_oac(r, t, qc(1, g), qc(2, g), qw(1, g), qw(2, g))
      else
         qw(1, g) = qc(1, g) ; qw(2, g) = qc(2, g)
      endif
   enddo
   endsubroutine world_points_2d_dev

   subroutine bodies_velocities(dims, nb, st, linear_velocities, angular_velocities, centres, lvs, avs, cs, wv)
   !< Resolve the per-body rigid velocity inputs (defaults 0) and check their shapes ([dims, nb], [3, nb], [dims, nb]);
   !< st = STATUS_INVALID_INPUT if one is wrong.
   integer(I4P), intent(in)           :: dims, nb                   !< Dimensions, bodies.
   integer(I4P), intent(inout)        :: st                         !< Status.
   real(R8P),    intent(in), optional :: linear_velocities(:,:)     !< Linear velocities.
   real(R8P),    intent(in), optional :: angular_velocities(:,:)    !< Angular velocities [3, nb] (2-D: [0, 0, omega]).
   real(R8P),    intent(in), optional :: centres(:,:)               !< Centres.
   real(R8P),    intent(inout)        :: lvs(:,:), avs(:,:), cs(:,:) !< Resolved [3, >= nb].
   logical,      intent(out)          :: wv                         !< Velocities asked.

   lvs = 0._R8P ; avs = 0._R8P ; cs = 0._R8P
   wv = present(linear_velocities) .or. present(angular_velocities)
   if (present(linear_velocities)) then
      if (size(linear_velocities, 1) /= dims .or. size(linear_velocities, 2) /= nb) then
         st = STATUS_INVALID_INPUT
      else
         lvs(1:dims, 1:nb) = linear_velocities
      endif
   endif
   if (present(angular_velocities)) then
      if (size(angular_velocities, 1) /= 3 .or. size(angular_velocities, 2) /= nb) then
         st = STATUS_INVALID_INPUT
      else
         avs(:, 1:nb) = angular_velocities
      endif
   endif
   if (present(centres)) then
      if (size(centres, 1) /= dims .or. size(centres, 2) /= nb) then
         st = STATUS_INVALID_INPUT
      else
         cs(1:dims, 1:nb) = centres
      endif
   endif
   endsubroutine bodies_velocities

   pure function default_shared(x, y, z, reach, dims) result(sd)
   !< The default shared distance: the largest reach(a) times the largest spacing along axis a.
   real(R8P),    intent(in) :: x(:), y(:), z(:) !< Point lines.
   integer(I4P), intent(in) :: reach(:)         !< Reach.
   integer(I4P), intent(in) :: dims             !< Dimensions.
   real(R8P)                :: sd               !< Shared distance.

   sd = max(reach(1) * spacing_max(x), reach(2) * spacing_max(y))
   if (dims == 3_I4P) sd = max(sd, reach(3) * spacing_max(z))
   endfunction default_shared

   pure function spacing_max(e) result(h)
   !< The largest spacing of a line (0 with fewer than two points).
   real(R8P), intent(in) :: e(:) !< Line.
   real(R8P)                :: h !< Spacing.
   integer(I4P)             :: i !< Counter.

   h = 0._R8P
   do i = 2_I4P, size(e, kind=I4P)
      h = max(h, e(i) - e(i - 1_I4P))
   enddo
   endfunction spacing_max

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
   allocate(cl%sub_body(0), cl%clearance(0), cl%second_body(0), cl%second_foot(dims, 0), cl%shared(0), cl%overlap(0), &
            cl%ilw_body(2 * dims, 0))
   endsubroutine closure_arrays

   subroutine closure_arrays_bodies(cl, nb, sd)
   !< After `closure_arrays`: the several-bodies fields sized for the n ghosts (`ilw_body` with ILW), the flag and scalars.
   type(ib_closure_t), intent(inout) :: cl !< Table.
   integer(I4P),       intent(in)    :: nb !< Surfaces (curves).
   real(R8P),          intent(in)    :: sd !< Shared distance.
   integer(I4P) :: n, nil                  !< Sizes.

   n = cl%n ; nil = 0_I4P ; if (cl%with_ilw) nil = n
   cl%with_bodies = .true. ; cl%bodies = nb ; cl%shared_distance = sd ; cl%overlap_points = 0_I4P
   deallocate(cl%sub_body, cl%clearance, cl%second_body, cl%second_foot, cl%shared, cl%overlap, cl%ilw_body)
   allocate(cl%sub_body(n), cl%clearance(n), cl%second_body(n), cl%second_foot(cl%dims, n), cl%shared(n), cl%overlap(n), &
            cl%ilw_body(2 * cl%dims, nil))
   endsubroutine closure_arrays_bodies

   subroutine map_create(w)
   !< Create the device copies of the arrays of a table built on the device (unstructured data on a local table: its
   !< arrays are later moved into the caller's, the present table keyed on the unchanged host addresses).
   type(ib_closure_t), intent(inout) :: w !< Table.

   logical :: ga, gb, gc, gv, gi, gk !< Non-empty groups (zero-size arrays are never mapped).
   logical :: gs, gsi                !< Non-empty bodies groups.

   call groups(w, ga, gb, gc, gv, gi, gk)
   !$acc enter data create(w%ijk, w%points, w%mask, w%steps, w%distance, w%foot, w%normal, w%element, w%feature, w%body, &
   !$acc&                  w%patch) if(ga)
   !$acc enter data create(w%lengths, w%probe, w%visible, w%probe_cell) if(gb)
   !$acc update device(w%lengths) if(gb)
   !$acc enter data create(w%kappa, w%direction) if(gc)
   !$acc enter data create(w%velocity) if(gv)
   !$acc enter data create(w%ilw_t, w%ilw_element, w%ilw_normal) if(gi)
   !$acc enter data create(w%ilw_kappa, w%ilw_direction) if(gk)
   call groups_bodies(w, gs, gsi)
   !$acc enter data create(w%sub_body, w%clearance, w%second_body, w%second_foot, w%shared, w%overlap) if(gs)
   !$acc enter data create(w%ilw_body) if(gsi)
   endsubroutine map_create

   subroutine map_delete(w)
   !< Release the device copies of a table built on the device without copying them back (the per-body work tables of
   !< the several-bodies device path).
   type(ib_closure_t), intent(inout) :: w !< Table.
   logical :: ga, gb, gc, gv, gi, gk     !< Non-empty groups.
   logical :: gs, gsi                    !< Non-empty bodies groups.

   call groups(w, ga, gb, gc, gv, gi, gk)
   !$acc exit data delete(w%ijk, w%points, w%mask, w%steps, w%distance, w%foot, w%normal, w%element, w%feature, w%body, &
   !$acc&                 w%patch) if(ga)
   !$acc exit data delete(w%lengths, w%probe, w%visible, w%probe_cell) if(gb)
   !$acc exit data delete(w%kappa, w%direction) if(gc)
   !$acc exit data delete(w%velocity) if(gv)
   !$acc exit data delete(w%ilw_t, w%ilw_element, w%ilw_normal) if(gi)
   !$acc exit data delete(w%ilw_kappa, w%ilw_direction) if(gk)
   call groups_bodies(w, gs, gsi)
   !$acc exit data delete(w%sub_body, w%clearance, w%second_body, w%second_foot, w%shared, w%overlap) if(gs)
   !$acc exit data delete(w%ilw_body) if(gsi)
   endsubroutine map_delete

   subroutine finish_device(w, res, cl)
   !< End of a device path: resident, the arrays move into the caller's table, flagged; else they are copied back and
   !< released first.
   type(ib_closure_t), intent(inout) :: w   !< Table built on the device.
   logical,            intent(in)    :: res !< Resident.
   type(ib_closure_t), intent(inout) :: cl  !< Caller's table.

   logical :: ga, gb, gc, gv, gi, gk !< Non-empty groups.
   logical :: gs, gsi                !< Non-empty bodies groups.

   if (.not. res .or. w%n < 1_I4P) then
      call groups_bodies(w, gs, gsi)
      !$acc exit data copyout(w%sub_body, w%clearance, w%second_body, w%second_foot, w%shared, w%overlap) if(gs)
      !$acc exit data copyout(w%ilw_body) if(gsi)
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
   cl%with_bodies = w%with_bodies ; cl%bodies = w%bodies ; cl%overlap_points = w%overlap_points
   cl%shared_distance = w%shared_distance
   call move_alloc(w%sub_body, cl%sub_body) ; call move_alloc(w%clearance, cl%clearance)
   call move_alloc(w%second_body, cl%second_body) ; call move_alloc(w%second_foot, cl%second_foot)
   call move_alloc(w%shared, cl%shared) ; call move_alloc(w%overlap, cl%overlap) ; call move_alloc(w%ilw_body, cl%ilw_body)
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

   subroutine groups_bodies(c, gs, gsi)
   !< The non-empty several-bodies groups of a table: the per-ghost fields, the ILW crossing bodies.
   type(ib_closure_t), intent(in)  :: c       !< Table.
   logical,            intent(out) :: gs, gsi !< Groups.

   gs = c%n > 0_I4P .and. c%with_bodies
   gsi = gs .and. c%with_ilw
   endsubroutine groups_bodies

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

   integer(I4P) function mis_l1(a, b) result(k)
   !< Differing entries of two logical rank-1 arrays (sizes differ: the larger size).
   logical, intent(in) :: a(:), b(:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a .neqv. b)
   endif
   endfunction mis_l1

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
