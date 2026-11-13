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
use fossil_curve_object,       only : curve_object, cut_cells_2d_t, cut_cells_2d_blocks_t, cut_cells_device_many, &
                                     place_for_crossings
use fossil_curve_slice,        only : slice_surface, slice_planes_t, slice_planes_prepare, slice_surface_plane, slice_triangles
use fossil_cut_cells_2d,       only : cut_cell_2d_pieces, cut_cell_2d_measures
use fossil_cut_cells_3d,       only : cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, cut_cells_3d_sections, &
                                     cut_cell_3d_pieces, cut_cell_3d_totals
use fossil_facet_object,       only : facet_object
use fossil_piece_map,          only : piece_map_2d_t, piece_map_3d_t, TRACK_UNTRUSTED, track_components, track_status_of
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t, pose_interpolate, pose_angle_2d, pose_log_3d, to_world, &
                                     to_body
use fossil_line_query,         only : line_crossing_before
use fossil_self_intersection,  only : intersection_pair_t
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, I8P, R8P
use, intrinsic :: iso_fortran_env, only : stderr => error_unit

implicit none
private
public :: swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, swept_cells_3d, gauss_legendre_01
public :: swept_cells_2d_device
public :: swept_cells_3d_device

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

integer(I4P), parameter :: SWEPT_CHUNK_CELLS = 2_I4P**22 !< Section cells per device launch set, and face sums per
                                                      !< plane group held on the host (`swept_cells_3d_device`).
real(R8P), parameter :: EVENT_GAP = 1e-10_R8P !< Piece map (issue #102): an interval between two events narrower than
                                             !< this is not evaluated (its events act as one; `map_cell`).
real(R8P), parameter :: EVENT_RES = 64._R8P * epsilon(1._R8P) !< Piece map (issue #104): the round-off of an
                                                             !< explicit event instant, in units of its formula's magnitudes
                                                             !< over its denominator (see `node_breakpoints`).
real(R8P), parameter :: FLOOR_H = 1e-6_R8P !< Pieces shorter than this converge to 64 eps of the measure times FLOOR_H
                                         !< (issue #69: slices within 1e-10 of a face carry ~1e-7 relative noise).
! A piece whose Gauss values miss that target at 64 nodes is still converged if the last change is within 64 eps of
! the round-off of ONE evaluation (issues #84, #86): 2-D, the coordinate magnitude xm = max |grid| + |t0| + |dt| (a
! crossing's round-off); 3-D, xm times the face side. Near an instant when a wall side lies along the face's line
! (2-D) or a body face on the grid plane (3-D) (a side or face starting on it, a corner leaving it tangentially),
! the crossings of a moved copy carry round-off amplified by 1 / angle, and the values stop agreeing at a floor
! that does not shrink with the piece. Anything a piece hides below that bound (noise, or an event within the
! round-off of a zero zone) moves the face's average by at most 64 eps xm (times the side in 3-D).

type :: rotation_planes_t
   !< The planes of one axis that `faces_rotation_3d` registers in device mode (issue #100 step 4): position, events
   !< and the convergence tolerances `one_plane` computed, for `rotation_axis_device`.
   integer(I4P)              :: n = 0_I4P   !< Planes.
   integer(I4P)              :: nbp = 0_I4P !< Breakpoints stored.
   integer(I4P), allocatable :: id(:)       !< Plane index along the axis [n].
   real(R8P),    allocatable :: lev(:)      !< Plane positions [n].
   integer(I4P), allocatable :: bfirst(:)   !< Events of plane p: bp(bfirst(p) : bfirst(p + 1) - 1) [n + 1].
   real(R8P),    allocatable :: bp(:)       !< Events (sorted per plane; capacity >= nbp).
   real(R8P),    allocatable :: tol(:,:)    !< sca, sco, sna, sno of `one_plane` [4, n].
endtype rotation_planes_t

type :: slist_t
   !< A list of breakpoints.
   real(R8P), allocatable :: s(:) !< Breakpoints.
endtype slist_t

type :: events_2d_t
   !< Events of one grid line (wall vertices crossing it) or of one grid node (it crossing wall segments) in the exact
   !< 2-D rotation (issue #94): they do not depend on the face, so every face on the line, or ending at the node, reuses
   !< them (the same operations, done once).
   logical                :: done = .false. !< Computed.
   logical                :: cert = .true.  !< Every `isolate_zeros` certified.
   integer(I4P)           :: n = 0_I4P      !< Events.
   real(R8P), allocatable :: z(:)           !< Instants [n] (a node: those inside a wall segment).
   real(R8P), allocatable :: w(:)           !< Line events: the coordinate along the line at the instant [n].
   logical,   allocatable :: u(:)           !< Not certified by a sign change (`isolate_zeros`, issue #102) [n].
   real(R8P), allocatable :: r(:)           !< Localisation half-width (`isolate_zeros`, issue #104) [n].
endtype events_2d_t

type :: rotation_2d_t
   !< State of `faces_rotation_2d` shared by its module-level helpers (issue #94).
   integer(I4P)                   :: nv = 0_I4P    !< Vertices.
   real(R8P),         allocatable :: vu(:), vv(:)  !< Vertices (body frame).
   integer(I4P),      allocatable :: nxt(:)        !< Next vertex (0: none).
   real(R8P)                      :: th0 = 0._R8P  !< Angle at pose0.
   real(R8P)                      :: dth = 0._R8P  !< Rotation over the step.
   real(R8P)                      :: t0(2) = 0._R8P, dt(2) = 0._R8P !< Translation at pose0, its change.
   real(R8P)                      :: xm = 0._R8P   !< Coordinate magnitude (round-off scale of a crossing).
   type(rigid_pose_2d_t)          :: pose0, pose1  !< End poses.
   type(curve_object)             :: cq            !< Work copy, placed at each Gauss node (`place_for_crossings`).
   logical                        :: exact = .true. !< Every breakpoint certified.
endtype rotation_2d_t

integer(I4P), parameter :: MAP_TRANSLATION = 1_I4P !< Piece map on the exact translation path (issue #102).
integer(I4P), parameter :: MAP_ROTATION    = 2_I4P !< Piece map on the exact rotation path.
integer(I4P), parameter :: MAP_GAUSS       = 3_I4P !< Piece map on the Gauss path (`exact=.false.`).

type :: map_ctx_t
   !< State of the piece map construction shared by its module-level helpers (issue #102).
   integer(I4P)              :: path = MAP_TRANSLATION !< MAP_TRANSLATION, MAP_ROTATION or MAP_GAUSS.
   integer(I4P)              :: nq = 4_I4P   !< Gauss nodes (MAP_GAUSS).
   integer(I4P)              :: nv = 0_I4P   !< Vertices.
   real(R8P),    allocatable :: bu(:), bv(:) !< Body vertices.
   real(R8P),    allocatable :: u0(:), v0(:) !< Vertices of the pose0 copy (`move`, as the end-pose cut).
   real(R8P),    allocatable :: u1(:), v1(:) !< Vertices of the pose1 copy.
   integer(I4P), allocatable :: nxt(:)       !< Next vertex (0: none).
   integer(I4P), allocatable :: loc(:)       !< Work: local index of a vertex in the current cell (where stamp is set).
   integer(I4P), allocatable :: stamp(:)     !< Work: the last cell that listed a vertex (its tick).
   integer(I4P)              :: tick = 0_I4P !< Work: a stamp per local curve.
   type(rigid_pose_2d_t)     :: pose0, pose1 !< End poses.
   real(R8P)                 :: d(2) = 0._R8P   !< Translation over the step (pose1%t - pose0%t).
   real(R8P)                 :: th0 = 0._R8P    !< Angle at pose0.
   real(R8P)                 :: dth = 0._R8P    !< Rotation over the step.
   real(R8P)                 :: xm = 0._R8P     !< Coordinate magnitude (the round-off scale of a crossing).
   real(R8P)                 :: bmin0(2) = 0._R8P !< Box minimum of the pose0 copy.
   real(R8P)                 :: bmin(2) = 0._R8P  !< Box minimum of the body.
   logical                   :: exact = .true.  !< Every multi-track integration converged.
   integer(I4P)              :: skipped = 0_I4P !< Quadrature nodes or pieces skipped in narrow intervals (`map_node`).
   integer(I4P)              :: contacts = 0_I4P !< Cells untrusted for an uncertified contact (`map_cell`).
   integer(I4P)              :: resolved = 0_I4P !< Cells with a state narrower than EVENT_GAP evaluated (issue #104).
endtype map_ctx_t

type :: cell_eval_t
   !< The pieces of one cell at one instant with their boundary features (issue #102, `map_eval`).
   logical                   :: closed = .false. !< The cell boundary closed.
   integer(I4P)              :: np = 0_I4P       !< Pieces.
   integer(I4P)              :: nk = 0_I4P       !< Features.
   real(R8P),    allocatable :: parea(:)         !< Piece areas [np].
   integer(I4P), allocatable :: key(:)           !< Features, increasing: fluid corners -4..-1, wall segments [nk].
   integer(I4P), allocatable :: kpiece(:)        !< Piece next to each feature [nk].
   real(R8P),    allocatable :: plen(:,:)        !< Fluid length of the faces (-u, +u, -v, +v) per piece [4, np].
   real(R8P),    allocatable :: pmom(:,:)        !< First moment of the faces' fluid parts per piece [4, np].
   logical                   :: fl(4) = .false.  !< The corner flags it was evaluated with.
endtype cell_eval_t

type :: map3_ctx_t
   !< State of the 3-D piece map construction shared by its module-level helpers (issue #102 step 2).
   integer(I4P)              :: path = MAP_TRANSLATION !< MAP_TRANSLATION, MAP_ROTATION or MAP_GAUSS.
   integer(I4P)              :: nq = 4_I4P      !< Gauss nodes (MAP_GAUSS).
   integer(I4P)              :: nf = 0_I4P      !< Facets.
   integer(I4P)              :: nvx = 0_I4P     !< Mesh vertices (distinct body-frame points).
   integer(I4P)              :: ned = 0_I4P     !< Mesh edges.
   logical                   :: rot = .false.   !< The step rotates: screw events, evaluations on posed facets.
   real(R8P),    allocatable :: tb(:,:,:)       !< Body-frame facets [3, 3, nf].
   real(R8P),    allocatable :: t0(:,:,:)       !< Facets of the pose0 copy (`move`, as the end-pose cut) [3, 3, nf].
   real(R8P),    allocatable :: t1(:,:,:)       !< Facets of the pose1 copy [3, 3, nf].
   real(R8P),    allocatable :: w(:,:,:)        !< R0 q of every facet vertex (screw events) [3, 3, nf].
   real(R8P),    allocatable :: vlo(:,:,:)      !< Box of every facet vertex's path over the step, low [3, 3, nf].
   real(R8P),    allocatable :: vhi(:,:,:)      !< Box of every facet vertex's path over the step, high [3, 3, nf].
   integer(I4P), allocatable :: fv(:,:)         !< Mesh vertex of every facet vertex [3, nf].
   integer(I4P), allocatable :: fe(:,:)         !< Mesh edge of every facet edge (m -> m + 1) [3, nf].
   integer(I4P), allocatable :: vstamp(:)       !< Work: the last cell that visited a mesh vertex (its tick).
   integer(I4P), allocatable :: estamp(:)       !< Work: the last cell that visited a mesh edge (its tick).
   integer(I4P)              :: tick = 0_I4P    !< Work: a stamp per cell.
   type(rigid_pose_t)        :: pose0, pose1    !< End poses.
   real(R8P)                 :: d(3) = 0._R8P   !< Translation over the step (pose1%t - pose0%t).
   real(R8P)                 :: kax(3) = 0._R8P !< Rotation axis of the step (unit).
   real(R8P)                 :: ang = 0._R8P    !< Rotation angle of the step.
   real(R8P)                 :: om(3) = 0._R8P  !< Rotation vector ang kax.
   real(R8P)                 :: xm = 0._R8P     !< Coordinate magnitude (the round-off scale of a crossing).
   real(R8P)                 :: tol = 0._R8P    !< Position tolerance of an event check, 1e-12 (max |grid| + 1).
   logical                   :: exact = .true.  !< Every multi-track integration converged.
   integer(I4P)              :: skipped = 0_I4P !< Quadrature nodes or pieces skipped in narrow intervals (`map3_node`).
   integer(I4P)              :: contacts = 0_I4P !< Cells untrusted for an uncertified contact (`map3_cell`).
   integer(I4P)              :: resolved = 0_I4P !< Cells with a state narrower than EVENT_GAP evaluated (issue #104).
   integer(I4P)              :: cs0(8) = 1_I4P  !< Current cell: corner states at pose0 (k0's node flags).
   logical                   :: known(8) = .false. !< Current cell: the corner's states follow from its events.
   integer(I4P)              :: ncn = 0_I4P     !< Current cell: corner events.
   real(R8P),    allocatable :: ctime(:)        !< Current cell: instants of the corner events [ncn].
   real(R8P),    allocatable :: ccor(:)         !< Current cell: their corners (1..8, as reals) [ncn].
endtype map3_ctx_t

type :: cell_eval3_t
   !< The pieces of one 3-D cell at one instant with their boundary features (issue #102 step 2, `map3_eval`).
   logical                   :: closed = .false. !< Trusted (every face section closed, the pieces counted).
   integer(I4P)              :: np = 0_I4P       !< Pieces.
   integer(I4P)              :: nk = 0_I4P       !< Features.
   integer(I4P), allocatable :: key(:)           !< Features, increasing: fluid corners -8..-1, facets (global ids) [nk].
   integer(I4P), allocatable :: kpiece(:)        !< Piece next to each feature [nk].
   real(R8P),    allocatable :: pvol(:)          !< Piece volumes (measures, several pieces) [np].
   real(R8P),    allocatable :: pa(:,:)          !< Fluid area of the faces (-x, +x, -y, +y, -z, +z) per piece [6, np].
   real(R8P),    allocatable :: pm(:,:,:)        !< int of the two in-face coordinates (cyclic) over them [2, 6, np].
endtype cell_eval3_t

contains
   subroutine swept_cells_2d(curve, pose0, pose1, x, y, sw, nq, poses, exact, status, check_input, map, refine)
   !< Swept geometry of `curve` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1) (issue #41 step 5; see the module notes). `nq` (default 4) Gauss-Legendre nodes for the
   !< rotating case; `poses(nq)` the poses at the nodes (`gauss_legendre_01`) instead of the screw
   !< motion; `exact` (default .true.) uses the exact integration for a pure translation and, since issue #69,
   !< for the screw motion with a rotation (certified breakpoints, Gauss doubled to round-off on every piece;
   !< `sw%exact` false if a piece did not converge at 64 nodes, to the round-off of one evaluation: see FLOOR_H).
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   !< `check_input` (default .false., issue #44): look for self-intersections of the curve first
   !< (`self_intersections`, one BVH pass); if any, `sw%self_intersecting` and a warning on stderr: the
   !< residual of the cells around a knot is not round-off (see the module notes). The result is
   !< computed anyway.
   !< `map` (issue #102): the piece map of the step, the tracks of every cell (`piece_map_2d_t`): the space-time
   !< components of the cell's fluid over the step, linking the pieces of the cut cells at pose0 to those at pose1, with
   !< per-track volumes, face apertures, wall and residual (see `piece_map_2d_build`). It needs the screw motion: with
   !< `poses`, `status = STATUS_INVALID_INPUT` and nothing is computed. Without `map` every other output is unchanged.
   !< On the Gauss path of the screw motion `sw%error_estimate` is the larger of the nq-vs-2nq face change and the
   !< largest |residual| over the largest rigid-field speed (issue #104): a face's Gauss error enters its cells' balance
   !< times the speed, and the face change alone understated it a hundredfold at 2 cells per step.
   !< `refine` (default .false., issue #104; screw motion only: with `poses`, `status = STATUS_INVALID_INPUT`): on the
   !< Gauss path every cell with an event inside the step integrates its faces per event interval on its own wall
   !< segments (the piece map's machinery, `piece_map_2d_build(all_meas=)`), so the residual drops to round-off; the
   !< other faces keep their Gauss values (no kink). Measured (256^2, 0.5-2 cells per step): 1.1-2.2 times the Gauss
   !< call, against 15-45 for the exact rotation path. Without `refine` every output but `error_estimate` is unchanged.
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
   type(piece_map_2d_t),   intent(out), optional :: map      !< Piece map of the step (issue #102).
   logical,                intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).

   call swept_cells_2d_core(curve, pose0, pose1, x, y, sw, .false., nq, poses, exact, status, check_input, map, refine)
   endsubroutine swept_cells_2d

   subroutine swept_cells_2d_device(curve, pose0, pose1, x, y, sw, nq, poses, exact, status, check_input, refine)
   !< `swept_cells_2d` with the cut cells of the two end poses on the device (issue #100 step 4): the same `sw`, bit for
   !< bit. Both moved copies go as the two blocks of one `cut_cells_device_many` call (one grid, cell arrays only, no
   !< pieces or intervals), left resident; only `alpha` comes back (`update host`), then the device copies are
   !< released. The face apertures, the wall terms and the residual stay on the host (`swept_cells_2d`'s routines):
   !< they are 94-99.9 % of a call at 256^2 and 1024^2 (step 4.0: `collect` and the node breakpoints of the
   !< translation, the line crossings of the Gauss path, `place_for_crossings` and the events of the rotation), so the
   !< end poses are the whole device part. Arguments, defaults and refusals as `swept_cells_2d`. Built without OpenACC
   !< it runs on the host. With `refine` the end-pose cuts are made on the host too (the refinement needs their pieces).
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
   logical,                intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).

   call swept_cells_2d_core(curve, pose0, pose1, x, y, sw, .true., nq, poses, exact, status, check_input, refine=refine)
   endsubroutine swept_cells_2d_device

   subroutine swept_cells_2d_core(curve, pose0, pose1, x, y, sw, device, nq, poses, exact, status, check_input, map, &
                                  refine)
   !< `swept_cells_2d` (device = .false.) and `swept_cells_2d_device` (device = .true.: the end-pose cut cells on the
   !< device); everything else is the same code. The piece map (`map`, host path only) is built last, from the end-pose
   !< cuts.
   type(curve_object),     intent(in)            :: curve    !< Curve (body frame).
   type(rigid_pose_2d_t),  intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_2d_t),  intent(in)            :: pose1    !< Pose at the end.
   real(R8P),              intent(in)            :: x(:), y(:) !< Cell edges (strictly increasing).
   type(swept_cells_2d_t), intent(out)           :: sw       !< Swept geometry.
   logical,                intent(in)            :: device   !< End-pose cut cells on the device.
   integer(I4P),           intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_2d_t),  intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                intent(in),  optional :: exact    !< Exact translation path (default .true.).
   integer(I4P),           intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                intent(in),  optional :: check_input !< Check the curve for self-intersections.
   type(piece_map_2d_t),   intent(out), optional :: map      !< Piece map of the step (issue #102).
   logical,                intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).
   integer(I4P), allocatable :: pairs(:,:)             !< Crossing pairs (check_input).
   type(piece_map_2d_t)     :: mp                    !< Piece map of the refinement (issue #104).
   logical,      allocatable :: meas_c(:), ru(:,:), rv(:,:) !< Cells integrated, faces replaced (refine).
   real(R8P)                :: vmax                  !< Largest rigid-field speed (error_estimate).
   logical                  :: rfn, flows            !< Refine; face outflows needed.
   type(curve_object)       :: ce(2)                 !< Copies at the two poses.
   type(cut_cells_2d_t)     :: k0, k1                !< Their cut cells (host).
   type(cut_cells_2d_blocks_t) :: kb                 !< Their cut cells (device, one call).
   integer(I4P)             :: st                    !< Status.
   real(R8P), allocatable   :: fu2(:,:), fv2(:,:), fl2(:,:), res(:,:) !< 2 nq face lengths, work, residual.
   real(R8P), allocatable   :: fou(:,:), fov(:,:)    !< Face outflows in +u / +v, time-integrated (the map only).
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
   rfn = .false. ; if (present(refine)) rfn = refine
   if (present(poses)) then
      if (size(poses) /= n .or. present(map) .or. rfn) then
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
   ce(1) = curve ; call ce(1)%move(pose0)
   ce(2) = curve ; call ce(2)%move(pose1)
   if (device) then ! both poses as the blocks of one call, cell arrays only; only alpha comes back
      call cut_cells_device_many(ce, x, y, kb, st, intervals=.false., pieces=.false., resident=.true.)
      !$acc update host(kb%alpha)
      sw%alpha0 = kb%alpha(:, :, 1) ; sw%alpha1 = kb%alpha(:, :, 2)
      call kb%exit_device
      if (rfn) then ! the refinement needs the end cuts' pieces: on the host, equal to the device's bit for bit
         call ce(1)%cut_cells(x, y, k0) ; call ce(2)%cut_cells(x, y, k1)
      endif
   else
      call ce(1)%cut_cells(x, y, k0) ; call ce(2)%cut_cells(x, y, k1)
      sw%alpha0 = k0%alpha ; sw%alpha1 = k1%alpha
   endif
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
   flows = present(map) .or. rfn
   if (flows) allocate(fou(nx + 1, ny), fov(nx, ny + 1))
   if (translation) then
      sw%exact = .true.
      call faces_translation_2d(ce(1), d, x, y, sw%face_u, sw%face_v, sw%evaluations)
      if (present(map)) then ! the field is uniform: the outflow is d times the fluid length
         fou = d(1) * sw%face_u ; fov = d(2) * sw%face_v
      endif
      do j = 1_I4P, ny ! the field is uniform (d): only the fluid lengths enter
         do i = 1_I4P, nx
            res(i, j) = d(1) * (sw%face_u(i + 1, j) - sw%face_u(i, j)) + d(2) * (sw%face_v(i, j + 1) - sw%face_v(i, j))
         enddo
      enddo
   elseif (screw .and. rotation_exact .and. dth /= 0._R8P) then
      ! exact rotation (issue #69): certified breakpoints per face, Gauss doubled to round-off on every piece
      if (present(map)) then
         call faces_rotation_2d(curve, pose0, pose1, x, y, sw%face_u, sw%face_v, res, sw%evaluations, sw%exact, fou, fov)
      else
         call faces_rotation_2d(curve, pose0, pose1, x, y, sw%face_u, sw%face_v, res, sw%evaluations, sw%exact)
      endif
   else
      if (flows) then ! no node poses with a map or refine
         call faces_gauss_2d(curve, pose0, pose1, x, y, n, sw%face_u, sw%face_v, res, fou=fou, fov=fov)
      else
         call faces_gauss_2d(curve, pose0, pose1, x, y, n, sw%face_u, sw%face_v, res, poses)
      endif
      sw%evaluations = n
      if (screw) then
         allocate(fu2(nx + 1, ny), fv2(nx, ny + 1), fl2(nx, ny))
         call faces_gauss_2d(curve, pose0, pose1, x, y, 2_I4P * n, fu2, fv2, fl2)
         if (rfn) then ! event cells per event interval (issue #104)
            call piece_map_2d_build(curve, ce(1), ce(2), pose0, pose1, x, y, k0, k1, sw, fou, fov, MAP_GAUSS, n, mp, &
                                    all_meas=.true., measured=meas_c)
            allocate(ru(nx + 1, ny), rv(nx, ny + 1))
            call refine_faces_2d(mp, meas_c, nx, ny, sw%face_u, sw%face_v, fou, fov, ru, rv)
            do j = 1_I4P, ny
               do i = 1_I4P, nx
                  res(i, j) = fou(i + 1, j) - fou(i, j) + fov(i, j + 1) - fov(i, j)
               enddo
            enddo
            sw%error_estimate = max(maxval(abs(fu2 - sw%face_u), mask=.not. ru), &
                                    maxval(abs(fv2 - sw%face_v), mask=.not. rv), 0._R8P)
         else
            sw%error_estimate = max(maxval(abs(fu2 - sw%face_u)), maxval(abs(fv2 - sw%face_v)))
         endif
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
   if (screw .and. .not. translation .and. .not. (rotation_exact .and. dth /= 0._R8P)) then ! the Gauss path (#104)
      vmax = abs(dth) * (maxval(abs([x, y])) + norm2(pose0%t) + norm2(d)) + norm2(d)
      if (vmax > 0._R8P) sw%error_estimate = max(sw%error_estimate, maxval(abs(sw%residual)) / vmax)
   endif
   if (rfn .and. present(map)) then
      call refine_map_2d(mp, meas_c, sw, fou, fov)
      map = mp
   elseif (present(map)) then
      if (translation) then
         call piece_map_2d_build(curve, ce(1), ce(2), pose0, pose1, x, y, k0, k1, sw, fou, fov, MAP_TRANSLATION, n, map)
      elseif (screw .and. rotation_exact .and. dth /= 0._R8P) then
         call piece_map_2d_build(curve, ce(1), ce(2), pose0, pose1, x, y, k0, k1, sw, fou, fov, MAP_ROTATION, n, map)
      else
         call piece_map_2d_build(curve, ce(1), ce(2), pose0, pose1, x, y, k0, k1, sw, fou, fov, MAP_GAUSS, n, map)
      endif
   endif
   endsubroutine swept_cells_2d_core

   subroutine faces_rotation_2d(curve, pose0, pose1, x, y, fu, fv, flux, nev, exact, fou, fov)
   !< Exact time-averaged fluid lengths and rigid-field outflows of every face under the screw motion of
   !< `pose_interpolate` with a rotation (issue #69; see the moving-bodies reference page, "Exact rotation"). A face
   !< outside the box swept by the body's bounding circle is fluid throughout (closed form). Every other face gets its
   !< breakpoints, the zeros of g(s) = (p0 + p1 s) cos phi + (q0 + q1 s) sin phi + k0 + k1 s isolated with the bound
   !< on |g''| (`isolate_zeros`): a wall vertex crossing the face's line inside its span, an end node crossing a wall
   !< segment inside it. On every piece between breakpoints the fluid length and the outflow are analytic: Gauss-
   !< Legendre doubled from 2 nodes until two values agree to round-off (at 64, to that of one evaluation: FLOOR_H;
   !< else exact = .false.). Each Gauss node places one work copy (`place_for_crossings`: the vertices and boxes of
   !< `move`, bit for bit) and the events of a grid line or node are computed once for all its faces (issue #94).
   type(curve_object),    intent(in)  :: curve        !< Curve (body frame).
   type(rigid_pose_2d_t), intent(in)  :: pose0, pose1 !< End poses.
   real(R8P),             intent(in)  :: x(:), y(:)   !< Cell edges.
   real(R8P),             intent(out) :: fu(:,:), fv(:,:) !< Time-averaged fluid lengths.
   real(R8P),             intent(out) :: flux(:,:)    !< Time integral of the rigid-field outflow per cell.
   integer(I4P),          intent(out) :: nev          !< Largest number of evaluations per face.
   logical,               intent(out) :: exact        !< Every piece converged.
   real(R8P), optional,   intent(out) :: fou(:,:), fov(:,:) !< Time integrals of the outflows in +u / +v (issue #102).
   type(rotation_2d_t)            :: ctx
   type(events_2d_t), allocatable :: line_u(:), line_v(:), node(:,:) !< Events of the lines u = x_i, v = y_j, the nodes.
   real(R8P),    allocatable :: ou(:,:), ov(:,:), bp(:)
   real(R8P)                 :: th0, dth, t0(2), dt(2), cen(2), rho, cbox(2, 2), uv(2), ss, cc(2), cmax
   integer(I4P)              :: nx, ny, nv, i, j, k, nb, ne
   logical                   :: ok

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   nv = curve%get_vertices_number()
   allocate(ctx%vu(nv), ctx%vv(nv), ctx%nxt(nv), ou(nx + 1, ny), ov(nx, ny + 1))
   do k = 1_I4P, nv
      uv = curve%get_vertex(k) ; ctx%vu(k) = uv(1) ; ctx%vv(k) = uv(2) ; ctx%nxt(k) = curve%get_next(k)
   enddo
   th0 = pose_angle_2d(pose0) ; dth = pose_angle_2d(pose1) - th0 ; dth = atan2(sin(dth), cos(dth))
   t0 = pose0%t ; dt = pose1%t - pose0%t
   ctx%nv = nv ; ctx%th0 = th0 ; ctx%dth = dth ; ctx%t0 = t0 ; ctx%dt = dt ; ctx%pose0 = pose0 ; ctx%pose1 = pose1
   ctx%xm = maxval(abs([x, y])) + norm2(t0) + norm2(dt) ! coordinate magnitude: the round-off scale of a crossing
   ctx%cq = curve ; ctx%exact = .true.
   allocate(line_u(nx + 1), line_v(ny + 1), node(nx + 1, ny + 1))
   ! the box swept by the body's bounding circle: centre path C(s) = R(th(s)) cen + t(s), |C''| <= dth^2 |cen|
   cen = 0.5_R8P * ([minval(ctx%vu), minval(ctx%vv)] + [maxval(ctx%vu), maxval(ctx%vv)])
   rho = sqrt(maxval((ctx%vu - cen(1))**2 + (ctx%vv - cen(2))**2))
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
         call face_events_2d(ctx, 1_I4P, x(i), y(j), y(j + 1), line_u(i), node(i, j), node(i, j + 1), bp, nb)
         call integrate_face_2d(ctx, curve, 1_I4P, x(i), y(j), y(j + 1), bp, nb, fu(i, j), ou(i, j), ne, ok)
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
         call face_events_2d(ctx, 2_I4P, y(j), x(i), x(i + 1), line_v(j), node(i, j), node(i + 1, j), bp, nb)
         call integrate_face_2d(ctx, curve, 2_I4P, y(j), x(i), x(i + 1), bp, nb, fv(i, j), ov(i, j), ne, ok)
         nev = max(nev, ne) ; exact = exact .and. ok
      enddo
   enddo
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         flux(i, j) = ou(i + 1, j) - ou(i, j) + ov(i, j + 1) - ov(i, j)
      enddo
   enddo
   exact = exact .and. ctx%exact
   if (present(fou)) fou = ou
   if (present(fov)) fov = ov
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
   endsubroutine faces_rotation_2d

   subroutine face_events_2d(ctx, nax, a, lo, hi, line, n1, n2, b, m)
   !< Breakpoints in (0, 1) of the face on the line normal to `nax` (1: u = a, 2: v = a), span [lo, hi], from the events
   !< of its line (filtered to the span) and of its two end nodes (lo, then hi), each computed at its first use.
   type(rotation_2d_t),    intent(inout) :: ctx      !< State.
   integer(I4P),           intent(in)    :: nax      !< Normal axis of the face.
   real(R8P),              intent(in)    :: a, lo, hi !< Line, span.
   type(events_2d_t),      intent(inout) :: line     !< Events of the line.
   type(events_2d_t),      intent(inout) :: n1, n2   !< Events of the end nodes at lo and at hi.
   real(R8P), allocatable, intent(out)   :: b(:)     !< Breakpoints.
   integer(I4P),           intent(out)   :: m        !< Their number.
   real(R8P), allocatable                :: tmp(:)   !< Work.
   integer(I4P)                          :: k        !< Counter.

   if (.not. line%done) call line_events_2d(ctx, nax, a, line)
   if (nax == 1_I4P) then
      if (.not. n1%done) call node_events_2d(ctx, [a, lo], n1)
      if (.not. n2%done) call node_events_2d(ctx, [a, hi], n2)
   else
      if (.not. n1%done) call node_events_2d(ctx, [lo, a], n1)
      if (.not. n2%done) call node_events_2d(ctx, [hi, a], n2)
   endif
   if (.not. (line%cert .and. n1%cert .and. n2%cert)) ctx%exact = .false.
   allocate(tmp(line%n + n1%n + n2%n))
   m = 0_I4P
   do k = 1_I4P, line%n
      if (line%w(k) >= lo .and. line%w(k) <= hi) then
         m = m + 1_I4P ; tmp(m) = line%z(k)
      endif
   enddo
   do k = 1_I4P, n1%n
      m = m + 1_I4P ; tmp(m) = n1%z(k)
   enddo
   do k = 1_I4P, n2%n
      m = m + 1_I4P ; tmp(m) = n2%z(k)
   enddo
   call sort_unique(tmp, m)
   call merge_close(tmp, m)
   b = tmp(1:m)
   endsubroutine face_events_2d

   subroutine line_events_2d(ctx, nax, a, ev)
   !< Wall vertices crossing the line normal to `nax` at a: X_nax(s) = a, with the other coordinate at the instant.
   type(rotation_2d_t), intent(in)    :: ctx  !< State.
   integer(I4P),        intent(in)    :: nax  !< Normal axis of the line.
   real(R8P),           intent(in)    :: a    !< Line position.
   type(events_2d_t),   intent(inout) :: ev   !< Events.
   real(R8P), allocatable             :: z(:) !< Zeros.
   logical,   allocatable             :: un(:) !< Zeros not certified by a sign change.
   real(R8P), allocatable             :: wd(:) !< Their widths.
   real(R8P)                          :: co(6), xs(2) !< Coefficients, position.
   integer(I4P)                       :: q, nz, kk !< Counters.
   logical                            :: cert !< Certified.

   ev%done = .true. ; ev%cert = .true. ; ev%n = 0_I4P
   if (allocated(ev%u)) deallocate(ev%u)
   if (allocated(ev%r)) deallocate(ev%r)
   allocate(ev%z(4), ev%w(4), ev%u(4), ev%r(4))
   do q = 1_I4P, ctx%nv
      if (nax == 1_I4P) then
         co = [ctx%vu(q), 0._R8P, -ctx%vv(q), 0._R8P, ctx%t0(1) - a, ctx%dt(1)]
      else
         co = [ctx%vv(q), 0._R8P, ctx%vu(q), 0._R8P, ctx%t0(2) - a, ctx%dt(2)]
      endif
      call isolate_zeros(co, ctx%th0, ctx%dth, abs(ctx%vu(q)) + abs(ctx%vv(q)) + norm2(ctx%t0) + abs(a) + norm2(ctx%dt), &
                         z, nz, cert, un, wd)
      if (.not. cert) ev%cert = .false.
      do kk = 1_I4P, nz
         xs = rot2(ctx%th0 + z(kk) * ctx%dth, [ctx%vu(q), ctx%vv(q)]) + ctx%t0 + z(kk) * ctx%dt
         call push_event(ev, z(kk), xs(3 - nax), un(kk), wd(kk))
      enddo
   enddo
   endsubroutine line_events_2d

   subroutine node_events_2d(ctx, node, ev)
   !< The grid node crossing wall segments: orient2d(X_q, X_q', N) = (R d) x (N - t(s)) - d x q, the instants at which
   !< the node lies inside the segment.
   type(rotation_2d_t), intent(in)    :: ctx     !< State.
   real(R8P),           intent(in)    :: node(2) !< Node.
   type(events_2d_t),   intent(inout) :: ev      !< Events.
   real(R8P), allocatable             :: z(:)    !< Zeros.
   logical,   allocatable             :: un(:)   !< Zeros not certified by a sign change.
   real(R8P), allocatable             :: wd(:)   !< Their widths.
   real(R8P)                          :: co(6), w0(2), w1(2), d(2), xs(2), r(2), lam !< Work.
   integer(I4P)                       :: q, nz, kk !< Counters.
   logical                            :: cert    !< Certified.

   ev%done = .true. ; ev%cert = .true. ; ev%n = 0_I4P
   if (allocated(ev%u)) deallocate(ev%u)
   if (allocated(ev%r)) deallocate(ev%r)
   allocate(ev%z(4), ev%w(4), ev%u(4), ev%r(4))
   w0 = node - ctx%t0 ; w1 = -ctx%dt
   do q = 1_I4P, ctx%nv
      if (ctx%nxt(q) == 0_I4P) cycle
      d = [ctx%vu(ctx%nxt(q)) - ctx%vu(q), ctx%vv(ctx%nxt(q)) - ctx%vv(q)]
      co = [d(1) * w0(2) - d(2) * w0(1), d(1) * w1(2) - d(2) * w1(1), &
            -(d(1) * w0(1) + d(2) * w0(2)), -(d(1) * w1(1) + d(2) * w1(2)), &
            -(d(1) * ctx%vv(q) - d(2) * ctx%vu(q)), 0._R8P]
      call isolate_zeros(co, ctx%th0, ctx%dth, norm2(d) * (norm2(w0) + norm2(w1) + abs(ctx%vu(q)) + abs(ctx%vv(q))), z, &
                         nz, cert, un, wd)
      if (.not. cert) ev%cert = .false.
      do kk = 1_I4P, nz
         xs = rot2(ctx%th0 + z(kk) * ctx%dth, [ctx%vu(q), ctx%vv(q)]) + ctx%t0 + z(kk) * ctx%dt
         r = rot2(ctx%th0 + z(kk) * ctx%dth, d)
         lam = dot_product(node - xs, r) / dot_product(r, r)
         if (lam >= 0._R8P .and. lam <= 1._R8P) call push_event(ev, z(kk), 0._R8P, un(kk), wd(kk))
      enddo
   enddo
   endsubroutine node_events_2d

   pure subroutine push_event(ev, z, w, u, r)
   !< Append one event (geometric growth, no self-referencing constructor).
   type(events_2d_t), intent(inout) :: ev        !< Events.
   real(R8P),         intent(in)    :: z, w      !< Instant, coordinate along the line.
   logical,           intent(in)    :: u         !< Not certified by a sign change.
   real(R8P),         intent(in)    :: r         !< Localisation half-width (issue #104).
   real(R8P), allocatable           :: tz(:), tw(:), tr(:) !< Grown buffers.
   logical,   allocatable           :: tu(:)     !< Grown buffer.

   if (ev%n == size(ev%z, kind=I4P)) then
      allocate(tz(2 * ev%n), tw(2 * ev%n), tu(2 * ev%n), tr(2 * ev%n))
      tz(1:ev%n) = ev%z(1:ev%n) ; tw(1:ev%n) = ev%w(1:ev%n) ; tu(1:ev%n) = ev%u(1:ev%n) ; tr(1:ev%n) = ev%r(1:ev%n)
      call move_alloc(tz, ev%z) ; call move_alloc(tw, ev%w) ; call move_alloc(tu, ev%u) ; call move_alloc(tr, ev%r)
   endif
   ev%n = ev%n + 1_I4P ; ev%z(ev%n) = z ; ev%w(ev%n) = w ; ev%u(ev%n) = u ; ev%r(ev%n) = r
   endsubroutine push_event

   subroutine integrate_face_2d(ctx, curve, nax, a, lo, hi, b, m, avg, out, ne, conv)
   !< Time averages over [0, 1] of the fluid length and of the rigid-field outflow of one face, piece by piece.
   type(rotation_2d_t), intent(inout) :: ctx      !< State.
   type(curve_object),  intent(in)    :: curve    !< Curve (body frame).
   integer(I4P),        intent(in)    :: nax      !< Normal axis of the face.
   real(R8P),           intent(in)    :: a, lo, hi !< Line, span.
   real(R8P),           intent(in)    :: b(:)     !< Breakpoints.
   integer(I4P),        intent(in)    :: m        !< Their number.
   real(R8P),           intent(out)   :: avg, out !< Time averages.
   integer(I4P),        intent(out)   :: ne       !< Evaluations.
   logical,             intent(out)   :: conv     !< Every piece converged.
   real(R8P)                          :: sa, sb, l1, o1, l2, o2, scl, sco, xm, sno !< Work.
   integer(I4P)                       :: q, ng    !< Work.

   avg = 0._R8P ; out = 0._R8P ; ne = 0_I4P ; conv = .true.
   xm = ctx%xm
   scl = hi - lo ; sco = scl * (abs(ctx%dth) * xm + norm2(ctx%dt))
   sno = xm * (abs(ctx%dth) * xm + norm2(ctx%dt))
   do q = 0_I4P, m
      sa = 0._R8P ; if (q > 0_I4P) sa = b(q)
      sb = 1._R8P ; if (q < m) sb = b(q + 1)
      if (sb <= sa) cycle
      ng = 2_I4P
      call piece_2d(ctx, curve, nax, a, lo, hi, sa, sb, ng, l1, o1) ; ne = ne + ng
      do
         ng = 2_I4P * ng
         call piece_2d(ctx, curve, nax, a, lo, hi, sa, sb, ng, l2, o2) ; ne = ne + ng
         if (abs(l2 - l1) <= 64._R8P * epsilon(1._R8P) * scl * max(sb - sa, FLOOR_H) .and. &
             abs(o2 - o1) <= 64._R8P * epsilon(1._R8P) * sco * max(sb - sa, FLOOR_H)) exit
         ! at 64 nodes: the change against the round-off of one evaluation (length: xm; issue #86, see FLOOR_H)
         if (ng >= 64_I4P) then
            conv = max(abs(l2 - l1) / xm, abs(o2 - o1) / sno) <= 64._R8P * epsilon(1._R8P) ; exit
         endif
         l1 = l2 ; o1 = o2
      enddo
      avg = avg + l2 ; out = out + o2
   enddo
   endsubroutine integrate_face_2d

   subroutine piece_2d(ctx, curve, nax, a, lo, hi, sa, sb, ng, l, o)
   !< ng-point Gauss-Legendre over [sa, sb] of the fluid length and the outflow of one face: at each node the work copy
   !< takes the vertices and tree boxes of `curve` moved to the node pose (`place_for_crossings`, as `move`).
   type(rotation_2d_t), intent(inout) :: ctx      !< State.
   type(curve_object),  intent(in)    :: curve    !< Curve (body frame).
   integer(I4P),        intent(in)    :: nax, ng  !< Normal axis of the face, nodes.
   real(R8P),           intent(in)    :: a, lo, hi, sa, sb !< Line, span, piece.
   real(R8P),           intent(out)   :: l, o     !< Fluid length and outflow integrated over the piece.
   real(R8P), allocatable             :: sq(:), wq(:) !< Gauss nodes, weights.
   real(R8P)                          :: ll, mm, sn, tq(2), dth, dt(2) !< Work.
   integer(I4P)                       :: q        !< Node.

   dth = ctx%dth ; dt = ctx%dt
   call gauss_legendre_01(ng, sq, wq)
   l = 0._R8P ; o = 0._R8P
   do q = 1_I4P, ng
      sn = sa + (sb - sa) * sq(q)
      call place_for_crossings(ctx%cq, curve, pose_interpolate(ctx%pose0, ctx%pose1, sn))
      tq = ctx%t0 + sn * dt
      if (nax == 1_I4P) then ! face u = a along v
         call face_measure(ctx%cq, 2_I4P, a, lo, hi, ll, mm)
         o = o + (sb - sa) * wq(q) * (-dth * (mm - tq(2) * ll) + dt(1) * ll)
      else                   ! face v = a along u
         call face_measure(ctx%cq, 1_I4P, a, lo, hi, ll, mm)
         o = o + (sb - sa) * wq(q) * (dth * (mm - tq(1) * ll) + dt(2) * ll)
      endif
      l = l + (sb - sa) * wq(q) * ll
   enddo
   endsubroutine piece_2d

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

   subroutine isolate_zeros(co, phi0, dphi, scale, z, nz, certified, uncertain, width)
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
   !< `uncertain` (issue #102, the piece map): per zero, true if no sign change certifies it (a zero zone's end, the
   !< middle of an interval narrower than GAP, or of one the stack cannot split) and |g| exceeds 4 err on both sides of
   !< it: a contact within the round-off of g inside the step (a tangency, or two crossings too close to separate),
   !< whose element may or may not change side. A zero with g at round-off all the way to an end of the step is a
   !< contact at the end pose, which the end cut decides.
   !< `width` (issue #104, the piece map): per zero, the half-width of the instants it may stand for in double precision:
   !< a sign change refined by bisection, the final bracket's half-width plus err over the slope of g there; any other
   !< zero, half the interval it stands for (a zero zone's ends: the whole zone). Two zeros farther apart than the sum
   !< of their widths are distinct crossings, whatever their gap.
   real(R8P),              intent(in)  :: co(6), phi0, dphi
   real(R8P),              intent(in)  :: scale !< Magnitude of the terms that make up the coefficients.
   real(R8P), allocatable, intent(out) :: z(:)
   integer(I4P),           intent(out) :: nz
   logical,                intent(out) :: certified
   logical, allocatable, intent(out), optional :: uncertain(:) !< Per zero: not certified by a sign change [nz].
   real(R8P), allocatable, intent(out), optional :: width(:) !< Per zero: localisation half-width (issue #104) [nz].
   integer(I4P), parameter             :: NSTACK = 4096, BUDGET = 4096
   real(R8P),    parameter             :: GAP = 1e-10_R8P !< Uncertified intervals stop here (see above).
   real(R8P)                           :: big_m, err, st(2, NSTACK), sa, sb, ga, gb, sm, gm, lo, hi, glo
   integer(I4P)                        :: top, used
   logical,      allocatable           :: uf(:) !< Flags of the zeros (with `uncertain`; at most 2 per interval examined).
   real(R8P),    allocatable           :: wf(:) !< Widths of the zeros (with `width`).
   real(R8P)                           :: slp   !< Slope of g over a bracket.
   integer(I4P)                        :: kz, j !< Work (`uncertain`).
   logical                             :: lft, rgt !< |g| above err before and after an uncertified zero.

   allocate(z(0)) ; nz = 0_I4P ; certified = .true. ; used = 0_I4P
   if (present(uncertain)) then
      allocate(uf(2_I4P * BUDGET + 2_I4P), uncertain(0))
   endif
   if (present(width)) then
      allocate(wf(2_I4P * BUDGET + 2_I4P), width(0))
   endif
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
         lo = sa ; hi = sb ; glo = ga ; slp = abs(gb - ga) / (sb - sa)
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
         call add(0.5_R8P * (lo + hi), .false., 0.5_R8P * (hi - lo) + err / max(slp, tiny(1._R8P)))
         cycle
      endif
      if (ga * gb > 0._R8P .and. min(abs(ga), abs(gb)) > big_m * (sb - sa)**2 / 8._R8P + err) cycle ! no zero
      if (abs(ga) <= err .and. abs(gb) <= err .and. big_m * (sb - sa)**2 / 8._R8P <= err) then
         ! a zero zone: |g| at round-off over the whole interval, the event undetermined in double precision
         call add(sa, .true., sb - sa) ; call add(sb, .true., sb - sa) ; cycle
      endif
      if (sb - sa < GAP) then
         call add(0.5_R8P * (sa + sb), .true., 0.5_R8P * (sb - sa)) ; cycle
      endif
      if (top + 2_I4P > NSTACK) then ! cannot happen for |g''| bounded; keep a breakpoint rather than drop one
         call add(0.5_R8P * (sa + sb), .true., 0.5_R8P * (sb - sa)) ; cycle
      endif
      sm = 0.5_R8P * (sa + sb)
      top = top + 1_I4P ; st(:, top) = [sm, sb]
      top = top + 1_I4P ; st(:, top) = [sa, sm]
   enddo
   if (present(uncertain)) then
      ! a contact inside the step has |g| clearly above err (4 err) on both sides of it (samples at z 2^-j towards 0, and
      ! towards 1); otherwise g stays near round-off up to an end: the contact is at the end pose, decided by the end cut
      do kz = 1_I4P, nz
         if (.not. uf(kz)) cycle
         lft = .false. ; rgt = .false.
         do j = 1_I4P, 30_I4P
            if (abs(g(z(kz) * (1._R8P - 0.5_R8P**j))) > 4._R8P * err) lft = .true.
            if (abs(g(z(kz) + (1._R8P - z(kz)) * (1._R8P - 0.5_R8P**j))) > 4._R8P * err) rgt = .true.
         enddo
         uf(kz) = lft .and. rgt
      enddo
      uncertain = uf(1:nz)
   endif
   if (present(width)) width = wf(1:nz)
   contains
      function g(s) result(v)
      real(R8P), intent(in) :: s
      real(R8P)             :: v, ph

      ph = phi0 + s * dphi
      v = (co(1) + co(2) * s) * cos(ph) + (co(3) + co(4) * s) * sin(ph) + co(5) + co(6) * s
      endfunction g

      subroutine add(s, u, w)
      real(R8P), intent(in) :: s
      logical,   intent(in) :: u !< Not certified by a sign change.
      real(R8P), intent(in) :: w !< Localisation half-width.

      if (s > GAP .and. s < 1._R8P - GAP) then
         z = [z, s] ; nz = nz + 1_I4P
         if (allocated(uf)) uf(nz) = u
         if (allocated(wf)) wf(nz) = w
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

   pure subroutine node_breakpoints(p, d, vu, vv, nxt, b, bw)
   !< Instants s in (0, 1) when the point p - s d lies on a segment of the curve; with `bw` (issue #104) the half-width
   !< of each instant's round-off (`EVENT_RES` times the magnitudes of Cramer's numerator and determinant over |det|).
   real(R8P),              intent(in)  :: p(2), d(2)  !< Node, displacement.
   real(R8P),              intent(in)  :: vu(:), vv(:) !< Vertices.
   integer(I4P),           intent(in)  :: nxt(:)      !< Segment ends.
   real(R8P), allocatable, intent(out) :: b(:)        !< Instants.
   real(R8P), allocatable, intent(out), optional :: bw(:) !< Their round-off half-widths.
   real(R8P), allocatable              :: tmp(:), tmw(:) !< Work.
   real(R8P)                           :: ex, ey, ax, ay, det, s, lam, plo(2), phi(2), xx, dd !< Work.
   integer(I4P)                        :: k, m        !< Work.

   allocate(tmp(size(vu)), tmw(size(vu)))
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
         xx = max(abs(vu(k)), abs(vu(nxt(k))), abs(vv(k)), abs(vv(nxt(k))), abs(p(1)), abs(p(2)))
         dd = max(abs(d(1)), abs(d(2)))
         tmw(m) = EVENT_RES * (4._R8P * xx * xx + 2._R8P * xx * dd * s) / abs(det)
      endif
   enddo
   b = tmp(1:m)
   if (present(bw)) bw = tmw(1:m)
   endsubroutine node_breakpoints

   pure subroutine sort_unique_w(a, w, m)
   !< Sort a(1:m) ascending (insertion), carrying w, and drop exact duplicates of a, keeping the larger w (issue #104).
   real(R8P),    intent(inout) :: a(:), w(:) !< Values, their half-widths.
   integer(I4P), intent(inout) :: m          !< Count.
   real(R8P)                   :: v, vw      !< Work.
   integer(I4P)                :: i, j, k    !< Work.

   do i = 2_I4P, m
      v = a(i) ; vw = w(i) ; j = i - 1_I4P
      do while (j >= 1_I4P)
         if (a(j) <= v) exit
         a(j + 1) = a(j) ; w(j + 1) = w(j) ; j = j - 1_I4P
      enddo
      a(j + 1) = v ; w(j + 1) = vw
   enddo
   k = min(m, 1_I4P)
   do i = 2_I4P, m
      if (a(i) /= a(k)) then
         k = k + 1_I4P ; a(k) = a(i) ; w(k) = w(i)
      else
         w(k) = max(w(k), w(i))
      endif
   enddo
   m = k
   endsubroutine sort_unique_w

   pure subroutine push_pair(a, w, n, z, r)
   !< Append the instant z with half-width r to a(1:n), w(1:n) (geometric growth; issue #104).
   real(R8P), allocatable, intent(inout) :: a(:), w(:) !< Instants, half-widths.
   integer(I4P),           intent(inout) :: n          !< Count.
   real(R8P),              intent(in)    :: z, r       !< Instant, half-width.
   real(R8P), allocatable                :: t(:)       !< Grown buffer.

   if (n == size(a, kind=I4P)) then
      allocate(t(2_I4P * n + 1_I4P)) ; t(1:n) = a(1:n) ; call move_alloc(t, a)
   endif
   if (n >= size(w, kind=I4P)) then
      allocate(t(size(a))) ; t(1:n) = w(1:n) ; call move_alloc(t, w)
   endif
   n = n + 1_I4P ; a(n) = z ; w(n) = r
   endsubroutine push_pair

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

   subroutine faces_gauss_2d(curve, pose0, pose1, x, y, n, fu, fv, flux, poses, fou, fov)
   !< Time-averaged fluid lengths by n-point Gauss-Legendre, and the time integral of the outflow of the
   !< rigid field of the screw motion per cell (`flux`, zero if `poses` given).
   type(curve_object),    intent(in)           :: curve      !< Curve (body frame).
   type(rigid_pose_2d_t), intent(in)           :: pose0, pose1 !< End poses.
   real(R8P),             intent(in)           :: x(:), y(:) !< Cell edges.
   integer(I4P),          intent(in)           :: n          !< Nodes.
   real(R8P),             intent(out)          :: fu(:,:), fv(:,:) !< Fluid lengths.
   real(R8P),             intent(out)          :: flux(:,:)  !< Outflow of the rigid field [nx, ny].
   type(rigid_pose_2d_t), intent(in), optional :: poses(:)   !< Node poses.
   real(R8P), intent(out), optional :: fou(:,:), fov(:,:)    !< Time integrals of the outflows in +u / +v (issue #102).
   type(curve_object)       :: cq                  !< Copy at a node.
   type(rigid_pose_2d_t)    :: pq                  !< Node pose.
   real(R8P), allocatable   :: sq(:), wq(:), ou(:,:), ov(:,:) !< Nodes, weights, face outflows.
   real(R8P)                :: l, m, d(2), dth, tq(2) !< Work.
   integer(I4P)             :: nx, ny, q, i, j     !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   call gauss_legendre_01(n, sq, wq)
   allocate(ou(nx + 1, ny), ov(nx, ny + 1))
   fu = 0._R8P ; fv = 0._R8P ; flux = 0._R8P
   if (present(fou)) fou = 0._R8P
   if (present(fov)) fov = 0._R8P
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
      if (present(fou)) fou = fou + wq(q) * ou
      if (present(fov)) fov = fov + wq(q) * ov
   enddo
   endsubroutine faces_gauss_2d

   subroutine piece_map_2d_build(curve, c0, c1, pose0, pose1, x, y, k0, k1, sw, fou, fov, path, nq, map, all_meas, measured)
   !< The piece map of the step (issue #102; `piece_map_2d_t`), from the end-pose cuts k0, k1 of the copies c0, c1.
   !<
   !< 1. Candidates: the cells that a segment can meet during the step (its endpoints' paths, boxed per sub-step with
   !<    the curvature bound of the screw motion), each with its local segments (CSR, increasing order).
   !< 2. Topology per cell (`map_cell`): the events of its four faces from its local segments (a wall vertex crossing a
   !<    face line inside the face, a corner crossing a wall segment: the translation formulas of `faces_translation_2d`,
   !<    or `line_events_2d` / `node_events_2d` certified by `isolate_zeros`); a cell without events whose two end cuts
   !<    have the same single piece (or none) takes the trivial track. Otherwise its pieces are evaluated at pose0 (the
   !<    end-pose copy, k0's corner flags: they must equal k0's pieces bit for bit, which fixes the numbering), at the
   !<    middle of every interval between events, and at pose1, and pieces of consecutive evaluations sharing a boundary
   !<    feature (a wall segment or a fluid corner) are one track (`track_components`).
   !< 3. Measures: one track takes the cell's `sw` values and face outflows (fou, fov); several tracks integrate their
   !<    face parts on the cell's event intervals (`map_measures`).
   !< 4. CSR fill in cell order.
   !< `all_meas` (issue #104, `refine`): every evaluated cell integrates its tracks per event interval, one track included
   !< (a single-track cell whose integration fails keeps the `sw` values); `measured` flags those cells (cell order).
   type(curve_object),     intent(in)  :: curve        !< Curve (body frame).
   type(curve_object),     intent(in)  :: c0, c1       !< Copies at pose0 and pose1 (`move`).
   type(rigid_pose_2d_t),  intent(in)  :: pose0, pose1 !< End poses.
   real(R8P),              intent(in)  :: x(:), y(:)   !< Cell edges.
   type(cut_cells_2d_t),   intent(in)  :: k0, k1       !< Cut cells of c0 and c1.
   type(swept_cells_2d_t), intent(in)  :: sw           !< Swept geometry of the step.
   real(R8P),              intent(in)  :: fou(:,:)     !< Time-integrated rigid-field outflow of the faces u = x_i (+u).
   real(R8P),              intent(in)  :: fov(:,:)     !< Time-integrated rigid-field outflow of the faces v = y_j (+v).
   integer(I4P),           intent(in)  :: path         !< MAP_TRANSLATION, MAP_ROTATION or MAP_GAUSS.
   integer(I4P),           intent(in)  :: nq           !< Gauss nodes (MAP_GAUSS).
   type(piece_map_2d_t),   intent(out) :: map          !< Piece map.
   logical,      optional,   intent(in)  :: all_meas   !< Integrate every evaluated cell (issue #104).
   logical, allocatable, optional, intent(out) :: measured(:) !< Cells integrated per event interval [nx ny].
   type(map_ctx_t)           :: mc                     !< State.
   integer(I4P), allocatable :: cstart(:), cseg(:)     !< Local segments of every cell.
   integer(I4P), allocatable :: ntr(:), ckind(:)       !< Tracks and kind (0 trivial, 1 evaluated, -1 untrusted) per cell.
   integer(I4P), allocatable :: otr(:), ntn(:)         !< Track (in its cell) of every old and new piece.
   integer(I4P), allocatable :: moff(:), oc(:), nc(:), fo(:), fn(:) !< Work.
   real(R8P),    allocatable :: meas(:,:), mloc(:,:)   !< Per-track measures of the multi-track cells.
   integer(I4P), allocatable :: so(:), sn(:)           !< Saved tracks of a cell's pieces (all_meas).
   integer(I4P)              :: nx, ny, ncell, c, i, j, nt, kind, nev, tot, t, p, a0, a1, b0, b1 !< Work.
   logical                   :: allm                   !< all_meas.

   allm = .false. ; if (present(all_meas)) allm = all_meas
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; ncell = nx * ny
   map%nx = nx ; map%ny = ny ; map%np0 = k0%np ; map%np1 = k1%np
   call map_context(curve, c0, c1, pose0, pose1, x, y, path, nq, mc)
   call map_candidates(mc, x, y, cstart, cseg)
   allocate(ntr(ncell), ckind(ncell), otr(k0%np), ntn(k1%np), moff(ncell))
   ! topology of every cell
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         a0 = k0%first_piece(i, j) ; a1 = a0 + k0%pieces(i, j) - 1_I4P
         b0 = k1%first_piece(i, j) ; b1 = b0 + k1%pieces(i, j) - 1_I4P
         call map_cell(mc, curve, c0, i, j, x, y, k0, k1, cseg(cstart(c):cstart(c + 1) - 1_I4P), .false., kind, nt, &
                       otr(a0:a1), ntn(b0:b1), mloc, nev)
         ckind(c) = kind ; ntr(c) = nt
         map%evaluations = map%evaluations + nev
         if (nev > 0_I4P) map%evaluated_cells = map%evaluated_cells + 1_I4P
      enddo
   enddo
   ! the multi-track cells again, with their per-track measures
   tot = 0_I4P ; moff = 0_I4P
   do c = 1_I4P, ncell
      if (ckind(c) == 1_I4P .and. (ntr(c) >= 2_I4P .or. allm)) then
         moff(c) = tot + 1_I4P ; tot = tot + ntr(c)
      endif
   enddo
   allocate(meas(8, tot))
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         if (moff(c) == 0_I4P) cycle
         a0 = k0%first_piece(i, j) ; a1 = a0 + k0%pieces(i, j) - 1_I4P
         b0 = k1%first_piece(i, j) ; b1 = b0 + k1%pieces(i, j) - 1_I4P
         so = otr(a0:a1) ; sn = ntn(b0:b1)
         call map_cell(mc, curve, c0, i, j, x, y, k0, k1, cseg(cstart(c):cstart(c + 1) - 1_I4P), .true., kind, nt, &
                       otr(a0:a1), ntn(b0:b1), mloc, nev)
         if (kind == 1_I4P .and. nt == ntr(c)) then
            meas(:, moff(c):moff(c) + nt - 1_I4P) = mloc(:, 1:nt)
         elseif (ntr(c) == 1_I4P .and. allm) then ! one track (all_meas): keep it with the sw values
            otr(a0:a1) = so ; ntn(b0:b1) = sn ; moff(c) = 0_I4P
         else
            ckind(c) = -1_I4P ; moff(c) = 0_I4P
         endif
      enddo
   enddo
   if (present(measured)) measured = moff > 0_I4P
   ! untrusted cells: one track with every piece
   allocate(map%tracks(nx, ny), map%first_track(nx, ny))
   nt = 0_I4P
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         if (ckind(c) == -1_I4P) then
            ntr(c) = 1_I4P
            a0 = k0%first_piece(i, j) ; a1 = a0 + k0%pieces(i, j) - 1_I4P
            b0 = k1%first_piece(i, j) ; b1 = b0 + k1%pieces(i, j) - 1_I4P
            otr(a0:a1) = 1_I4P ; ntn(b0:b1) = 1_I4P
            map%untrusted_cells = map%untrusted_cells + 1_I4P
         endif
         if (ntr(c) >= 2_I4P) map%multi_track_cells = map%multi_track_cells + 1_I4P
         map%tracks(i, j) = ntr(c) ; map%first_track(i, j) = nt + 1_I4P ; nt = nt + ntr(c)
      enddo
   enddo
   map%nt = nt
   allocate(map%track_cell(2, nt), map%status(nt), map%old_start(nt + 1), map%new_start(nt + 1), map%old_piece(k0%np), &
            map%new_piece(k1%np), map%old_track(k0%np), map%new_track(k1%np), map%volume0(nt), map%volume1(nt), &
            map%faces(4, nt), map%wall(2, nt), map%residual(nt), map%outflow(4, nt), oc(nt), nc(nt), fo(nt), fn(nt))
   map%volume0 = 0._R8P ; map%volume1 = 0._R8P ; oc = 0_I4P ; nc = 0_I4P
   do p = 1_I4P, k0%np ! pieces in order: volumes summed in piece order
      t = map%first_track(k0%piece_cell(1, p), k0%piece_cell(2, p)) + otr(p) - 1_I4P
      map%old_track(p) = t ; oc(t) = oc(t) + 1_I4P ; map%volume0(t) = map%volume0(t) + k0%piece_area(p)
   enddo
   do p = 1_I4P, k1%np
      t = map%first_track(k1%piece_cell(1, p), k1%piece_cell(2, p)) + ntn(p) - 1_I4P
      map%new_track(p) = t ; nc(t) = nc(t) + 1_I4P ; map%volume1(t) = map%volume1(t) + k1%piece_area(p)
   enddo
   map%old_start(1) = 1_I4P ; map%new_start(1) = 1_I4P
   do t = 1_I4P, nt
      map%old_start(t + 1) = map%old_start(t) + oc(t) ; map%new_start(t + 1) = map%new_start(t) + nc(t)
      fo(t) = map%old_start(t) ; fn(t) = map%new_start(t)
   enddo
   do p = 1_I4P, k0%np
      t = map%old_track(p) ; map%old_piece(fo(t)) = p ; fo(t) = fo(t) + 1_I4P
   enddo
   do p = 1_I4P, k1%np
      t = map%new_track(p) ; map%new_piece(fn(t)) = p ; fn(t) = fn(t) + 1_I4P
   enddo
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         do t = map%first_track(i, j), map%first_track(i, j) + ntr(c) - 1_I4P
            map%track_cell(:, t) = [i, j]
            if (ckind(c) == -1_I4P) then
               map%status(t) = TRACK_UNTRUSTED
            else
               map%status(t) = track_status_of(oc(t), nc(t))
            endif
            map%n_status(map%status(t)) = map%n_status(map%status(t)) + 1_I4P
            if (moff(c) > 0_I4P) then ! several tracks: integrated per track
               p = moff(c) + t - map%first_track(i, j)
               map%faces(:, t) = meas(1:4, p)
               map%outflow(:, t) = [-meas(5, p), meas(6, p), -meas(7, p), meas(8, p)]
               map%wall(:, t) = -[meas(2, p) - meas(1, p), meas(4, p) - meas(3, p)]
               map%residual(t) = map%volume1(t) - map%volume0(t) + (meas(6, p) - meas(5, p) + meas(8, p) - meas(7, p))
            else ! one track: the cell's values
               map%faces(:, t) = [sw%face_u(i, j), sw%face_u(i + 1, j), sw%face_v(i, j), sw%face_v(i, j + 1)]
               map%outflow(:, t) = [-fou(i, j), fou(i + 1, j), -fov(i, j), fov(i, j + 1)]
               map%wall(:, t) = sw%wall(:, i, j)
               map%residual(t) = sw%residual(i, j)
            endif
         enddo
      enddo
   enddo
   map%exact = mc%exact ; map%skipped_nodes = mc%skipped
   map%contact_cells = mc%contacts ; map%narrow_cells = mc%resolved
   endsubroutine piece_map_2d_build

   subroutine refine_faces_2d(map, measured, nx, ny, fu, fv, fou, fov, ru, rv)
   !< The Gauss path's faces next to the wall integrated per event interval (issue #104, `refine`): every cell the map
   !< integrated (`measured`, `piece_map_2d_build(all_meas=)`) gives its four faces and their outflows as the sums
   !< over its tracks; a face between two such cells takes the value of the low one. Every face with an event inside the
   !< step belongs to such a cell (an event on a face is an event of both its cells); the others keep their Gauss
   !< values, accurate for an integrand without kinks. ru, rv flag the faces replaced.
   type(piece_map_2d_t), intent(in)    :: map              !< Piece map (all_meas).
   logical,              intent(in)    :: measured(:)      !< Cells integrated per event interval [nx ny].
   integer(I4P),         intent(in)    :: nx, ny           !< Cells.
   real(R8P),            intent(inout) :: fu(:,:), fv(:,:) !< Time-averaged fluid lengths.
   real(R8P),            intent(inout) :: fou(:,:), fov(:,:) !< Time-integrated outflows in +u / +v.
   logical,              intent(out)   :: ru(:,:), rv(:,:) !< Faces replaced.
   real(R8P)                           :: fl(4), ol(4)     !< Sums over a cell's tracks.
   integer(I4P)                        :: i, j, c, t       !< Work.

   ru = .false. ; rv = .false.
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         if (.not. measured(c)) cycle
         fl = 0._R8P ; ol = 0._R8P
         do t = map%first_track(i, j), map%first_track(i, j) + map%tracks(i, j) - 1_I4P
            fl = fl + map%faces(:, t) ; ol = ol + map%outflow(:, t)
         enddo
         fu(i + 1, j) = fl(2) ; fou(i + 1, j) = ol(2) ; ru(i + 1, j) = .true.
         fv(i, j + 1) = fl(4) ; fov(i, j + 1) = ol(4) ; rv(i, j + 1) = .true.
         if (i == 1_I4P) then
            fu(i, j) = fl(1) ; fou(i, j) = -ol(1) ; ru(i, j) = .true.
         elseif (.not. measured(c - 1_I4P)) then
            fu(i, j) = fl(1) ; fou(i, j) = -ol(1) ; ru(i, j) = .true.
         endif
         if (j == 1_I4P) then
            fv(i, j) = fl(3) ; fov(i, j) = -ol(3) ; rv(i, j) = .true.
         elseif (.not. measured(c - nx)) then
            fv(i, j) = fl(3) ; fov(i, j) = -ol(3) ; rv(i, j) = .true.
         endif
      enddo
   enddo
   endsubroutine refine_faces_2d

   subroutine refine_map_2d(map, measured, sw, fou, fov)
   !< The tracks of the cells `refine_faces_2d` did not integrate take the final `sw` values (they were built from the
   !< Gauss faces): faces, outflows, wall, residual (issue #104).
   type(piece_map_2d_t),   intent(inout) :: map          !< Piece map.
   logical,                intent(in)    :: measured(:)  !< Cells integrated per event interval.
   type(swept_cells_2d_t), intent(in)    :: sw           !< Swept geometry (refined).
   real(R8P),              intent(in)    :: fou(:,:), fov(:,:) !< Time-integrated outflows in +u / +v (refined).
   integer(I4P)                          :: i, j, t      !< Work.

   do j = 1_I4P, map%ny
      do i = 1_I4P, map%nx
         if (measured(i + (j - 1_I4P) * map%nx)) cycle
         do t = map%first_track(i, j), map%first_track(i, j) + map%tracks(i, j) - 1_I4P
            map%faces(:, t) = [sw%face_u(i, j), sw%face_u(i + 1, j), sw%face_v(i, j), sw%face_v(i, j + 1)]
            map%outflow(:, t) = [-fou(i, j), fou(i + 1, j), -fov(i, j), fov(i, j + 1)]
            map%wall(:, t) = sw%wall(:, i, j) ; map%residual(t) = sw%residual(i, j)
         enddo
      enddo
   enddo
   endsubroutine refine_map_2d

   subroutine refine_faces_3d(map, measured, nx, ny, nz, fx, fy, fz, fox, foy, foz, rx, ry, rz)
   !< The 3-D twin of `refine_faces_2d`: six faces per integrated cell, the low cell's value on a shared face.
   type(piece_map_3d_t), intent(in)    :: map                       !< Piece map (all_meas).
   logical,              intent(in)    :: measured(:)               !< Cells integrated per event interval [nx ny nz].
   integer(I4P),         intent(in)    :: nx, ny, nz                !< Cells.
   real(R8P),            intent(inout) :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Time-averaged fluid areas.
   real(R8P),            intent(inout) :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Time-integrated outflows in +axis.
   logical,              intent(out)   :: rx(:,:,:), ry(:,:,:), rz(:,:,:) !< Faces replaced.
   real(R8P)                           :: fl(6), ol(6)              !< Sums over a cell's tracks.
   integer(I4P)                        :: i, j, k, c, t             !< Work.

   rx = .false. ; ry = .false. ; rz = .false.
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            if (.not. measured(c)) cycle
            fl = 0._R8P ; ol = 0._R8P
            do t = map%first_track(i, j, k), map%first_track(i, j, k) + map%tracks(i, j, k) - 1_I4P
               fl = fl + map%faces(:, t) ; ol = ol + map%outflow(:, t)
            enddo
            fx(i + 1, j, k) = fl(2) ; fox(i + 1, j, k) = ol(2) ; rx(i + 1, j, k) = .true.
            fy(i, j + 1, k) = fl(4) ; foy(i, j + 1, k) = ol(4) ; ry(i, j + 1, k) = .true.
            fz(i, j, k + 1) = fl(6) ; foz(i, j, k + 1) = ol(6) ; rz(i, j, k + 1) = .true.
            if (i == 1_I4P) then
               fx(i, j, k) = fl(1) ; fox(i, j, k) = -ol(1) ; rx(i, j, k) = .true.
            elseif (.not. measured(c - 1_I4P)) then
               fx(i, j, k) = fl(1) ; fox(i, j, k) = -ol(1) ; rx(i, j, k) = .true.
            endif
            if (j == 1_I4P) then
               fy(i, j, k) = fl(3) ; foy(i, j, k) = -ol(3) ; ry(i, j, k) = .true.
            elseif (.not. measured(c - nx)) then
               fy(i, j, k) = fl(3) ; foy(i, j, k) = -ol(3) ; ry(i, j, k) = .true.
            endif
            if (k == 1_I4P) then
               fz(i, j, k) = fl(5) ; foz(i, j, k) = -ol(5) ; rz(i, j, k) = .true.
            elseif (.not. measured(c - nx * ny)) then
               fz(i, j, k) = fl(5) ; foz(i, j, k) = -ol(5) ; rz(i, j, k) = .true.
            endif
         enddo
      enddo
   enddo
   endsubroutine refine_faces_3d

   subroutine refine_map_3d(map, measured, sw, fox, foy, foz)
   !< The 3-D twin of `refine_map_2d`.
   type(piece_map_3d_t),   intent(inout) :: map          !< Piece map.
   logical,                intent(in)    :: measured(:)  !< Cells integrated per event interval.
   type(swept_cells_3d_t), intent(in)    :: sw           !< Swept geometry (refined).
   real(R8P),              intent(in)    :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Outflows in +axis (refined).
   integer(I4P)                          :: i, j, k, t   !< Work.

   do k = 1_I4P, map%nz
      do j = 1_I4P, map%ny
         do i = 1_I4P, map%nx
            if (measured(i + map%nx * ((j - 1_I4P) + map%ny * (k - 1_I4P)))) cycle
            do t = map%first_track(i, j, k), map%first_track(i, j, k) + map%tracks(i, j, k) - 1_I4P
               map%faces(:, t) = [sw%face_x(i, j, k), sw%face_x(i + 1, j, k), sw%face_y(i, j, k), sw%face_y(i, j + 1, k), &
                                  sw%face_z(i, j, k), sw%face_z(i, j, k + 1)]
               map%outflow(:, t) = [-fox(i, j, k), fox(i + 1, j, k), -foy(i, j, k), foy(i, j + 1, k), -foz(i, j, k), &
                                    foz(i, j, k + 1)]
               map%wall(:, t) = sw%wall(:, i, j, k) ; map%residual(t) = sw%residual(i, j, k)
            enddo
         enddo
      enddo
   enddo
   endsubroutine refine_map_3d

   subroutine map_context(curve, c0, c1, pose0, pose1, x, y, path, nq, mc)
   !< The state of the piece map construction: vertices (body, pose0 and pose1 copies), motion, work arrays.
   type(curve_object),    intent(in)  :: curve        !< Curve (body frame).
   type(curve_object),    intent(in)  :: c0, c1       !< Copies at pose0 and pose1.
   type(rigid_pose_2d_t), intent(in)  :: pose0, pose1 !< End poses.
   real(R8P),             intent(in)  :: x(:), y(:)   !< Cell edges.
   integer(I4P),          intent(in)  :: path, nq     !< Path, Gauss nodes.
   type(map_ctx_t),       intent(out) :: mc           !< State.
   real(R8P)                          :: uv(2)        !< Vertex.
   integer(I4P)                       :: g            !< Vertex.

   mc%path = path ; mc%nq = nq ; mc%nv = curve%get_vertices_number()
   allocate(mc%bu(mc%nv), mc%bv(mc%nv), mc%u0(mc%nv), mc%v0(mc%nv), mc%u1(mc%nv), mc%v1(mc%nv), mc%nxt(mc%nv), &
            mc%loc(mc%nv), mc%stamp(mc%nv))
   do g = 1_I4P, mc%nv
      uv = curve%get_vertex(g) ; mc%bu(g) = uv(1) ; mc%bv(g) = uv(2)
      uv = c0%get_vertex(g) ; mc%u0(g) = uv(1) ; mc%v0(g) = uv(2)
      uv = c1%get_vertex(g) ; mc%u1(g) = uv(1) ; mc%v1(g) = uv(2)
      mc%nxt(g) = curve%get_next(g)
   enddo
   mc%loc = 0_I4P ; mc%stamp = 0_I4P ; mc%tick = 0_I4P
   mc%pose0 = pose0 ; mc%pose1 = pose1 ; mc%d = pose1%t - pose0%t
   mc%th0 = pose_angle_2d(pose0) ; mc%dth = pose_angle_2d(pose1) - mc%th0 ; mc%dth = atan2(sin(mc%dth), cos(mc%dth))
   mc%xm = max(maxval(abs(x)), maxval(abs(y))) + sqrt(pose0%t(1)**2 + pose0%t(2)**2) + sqrt(mc%d(1)**2 + mc%d(2)**2)
   mc%bmin0 = c0%get_bmin() ; mc%bmin = curve%get_bmin() ; mc%exact = .true.
   endsubroutine map_context

   subroutine map_candidates(mc, x, y, cstart, cseg)
   !< The cells a segment can meet during the step, with their segments (CSR over the cells c = i + (j - 1) nx, the
   !< segments of a cell in increasing order): the segment boxed per sub-step (`map_segment_box`), the cells meeting a
   !< box (closed) listed once per segment. Count, then fill.
   type(map_ctx_t),           intent(in)  :: mc         !< State.
   real(R8P),                 intent(in)  :: x(:), y(:) !< Cell edges.
   integer(I4P), allocatable, intent(out) :: cstart(:)  !< [ncell + 1].
   integer(I4P), allocatable, intent(out) :: cseg(:)    !< Segments.
   integer(I4P), allocatable :: cstamp(:), fill(:)      !< Last segment per cell, fill pointers.
   real(R8P)                 :: lo(2), hi(2), rho, h, l !< Box, work.
   integer(I4P)              :: nx, ny, nsub, pass, q, k, i1, i2, j1, j2, ii, jj, c !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   nsub = 1_I4P
   if (mc%path /= MAP_TRANSLATION .and. mc%nv > 0_I4P) then ! sub-steps no longer than the smallest cell side
      rho = sqrt(maxval(mc%bu**2 + mc%bv**2))
      l = abs(mc%dth) * rho + sqrt(mc%d(1)**2 + mc%d(2)**2)
      h = min(minval(x(2:) - x(:nx)), minval(y(2:) - y(:ny)))
      nsub = int(min(64._R8P, max(1._R8P, real(ceiling(l / h), R8P))), I4P)
   endif
   allocate(cstart(nx * ny + 1), cstamp(nx * ny), fill(nx * ny))
   cstart = 0_I4P
   do pass = 1_I4P, 2_I4P
      cstamp = 0_I4P
      if (pass == 2_I4P) then
         k = 1_I4P
         do c = 1_I4P, nx * ny
            fill(c) = k ; k = k + cstart(c) ; cstart(c) = fill(c)
         enddo
         cstart(nx * ny + 1) = k
         allocate(cseg(k - 1_I4P))
      endif
      do q = 1_I4P, mc%nv
         if (mc%nxt(q) == 0_I4P) cycle
         do k = 1_I4P, nsub
            call map_segment_box(mc, q, k, nsub, lo, hi)
            call index_range(x, lo(1), hi(1), i1, i2) ; call index_range(y, lo(2), hi(2), j1, j2)
            i1 = max(1_I4P, i1 - 1_I4P) ; i2 = min(nx, i2) ; j1 = max(1_I4P, j1 - 1_I4P) ; j2 = min(ny, j2)
            do jj = j1, j2
               do ii = i1, i2
                  c = ii + (jj - 1_I4P) * nx
                  if (cstamp(c) == q) cycle
                  cstamp(c) = q
                  if (pass == 1_I4P) then
                     cstart(c) = cstart(c) + 1_I4P
                  else
                     cseg(fill(c)) = q ; fill(c) = fill(c) + 1_I4P
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine map_candidates

   subroutine map_segment_box(mc, q, k, nsub, lo, hi)
   !< A box holding segment q over sub-step k of nsub: translation, the box of its ends at pose0 and moved by the
   !< displacement (the segment sweeps their convex hull); screw motion, the box of its ends at the two sub-step
   !< instants grown by the bound on the arcs, dth^2 |p| h^2 / 8 (|X''| = dth^2 |p| for X = R(th) p + t). Grown by 64
   !< eps of the coordinate magnitude in every case.
   type(map_ctx_t), intent(in)  :: mc           !< State.
   integer(I4P),    intent(in)  :: q, k, nsub   !< Segment, sub-step, sub-steps.
   real(R8P),       intent(out) :: lo(2), hi(2) !< Box.
   real(R8P)                    :: sa, sb, pad, pa(2), pb(2) !< Work.
   integer(I4P)                 :: w, e, g      !< Work.

   w = mc%nxt(q)
   if (mc%path == MAP_TRANSLATION) then
      lo(1) = min(mc%u0(q), mc%u0(w), mc%u0(q) + mc%d(1), mc%u0(w) + mc%d(1))
      hi(1) = max(mc%u0(q), mc%u0(w), mc%u0(q) + mc%d(1), mc%u0(w) + mc%d(1))
      lo(2) = min(mc%v0(q), mc%v0(w), mc%v0(q) + mc%d(2), mc%v0(w) + mc%d(2))
      hi(2) = max(mc%v0(q), mc%v0(w), mc%v0(q) + mc%d(2), mc%v0(w) + mc%d(2))
      pad = 0._R8P
   else
      sa = real(k - 1_I4P, R8P) / real(nsub, R8P) ; sb = real(k, R8P) / real(nsub, R8P)
      lo = huge(1._R8P) ; hi = -huge(1._R8P) ; pad = 0._R8P
      do e = 1_I4P, 2_I4P
         g = q ; if (e == 2_I4P) g = w
         pa = rot2(mc%th0 + sa * mc%dth, [mc%bu(g), mc%bv(g)]) + mc%pose0%t + sa * mc%d
         pb = rot2(mc%th0 + sb * mc%dth, [mc%bu(g), mc%bv(g)]) + mc%pose0%t + sb * mc%d
         lo = min(lo, pa, pb) ; hi = max(hi, pa, pb)
         pad = max(pad, mc%dth**2 * sqrt(mc%bu(g)**2 + mc%bv(g)**2) * (sb - sa)**2 / 8._R8P)
      enddo
   endif
   pad = pad + 64._R8P * epsilon(1._R8P) * mc%xm
   lo = lo - pad ; hi = hi + pad
   endsubroutine map_segment_box

   subroutine map_cell(mc, curve, c0, i, j, x, y, k0, k1, lseg, with_meas, kind, nt, otr, ntn, meas, nev)
   !< The tracks of cell (i, j) (see `piece_map_2d_build`, step 2): kind 0 (trivial: no evaluation), 1 (evaluated) or
   !< -1 (untrusted); nt tracks; the track of each of its old pieces (otr, in k0's order) and new pieces (ntn); with
   !< `with_meas` and nt >= 1, the per-track time integrals meas(1:4) of the faces' fluid lengths and meas(5:8) of
   !< their rigid-field outflows (`map_measures`; kind -1 if a face part could not be assigned).
   !<
   !< The events split [0, 1] into intervals; one narrower than EVENT_GAP (two events computed apart by round-off, or
   !< closer than the precision of an evaluation) is not evaluated: its two events act as one. The others are
   !< evaluated at their middle, at least EVENT_GAP / 2 from any event. A cell with an event that no sign change
   !< certifies (issue #102: a contact within the round-off of its function under a rotation, a corner grazing a
   !< segment or a vertex grazing a face line at the turning point of its arc) is untrusted: the piece next to it may
   !< shrink to a point and come back as another track, or not, undecidably in double precision.
   type(map_ctx_t),           intent(inout) :: mc         !< State.
   type(curve_object),        intent(in)    :: curve, c0  !< Curve (body frame), copy at pose0.
   integer(I4P),              intent(in)    :: i, j       !< Cell.
   real(R8P),                 intent(in)    :: x(:), y(:) !< Cell edges.
   type(cut_cells_2d_t),      intent(in)    :: k0, k1     !< End-pose cuts.
   integer(I4P),              intent(in)    :: lseg(:)    !< Local segments (increasing).
   logical,                   intent(in)    :: with_meas  !< Integrate the per-track measures.
   integer(I4P),              intent(out)   :: kind       !< 0 trivial, 1 evaluated, -1 untrusted.
   integer(I4P),              intent(out)   :: nt         !< Tracks.
   integer(I4P),              intent(out)   :: otr(:)     !< Track of each old piece of the cell.
   integer(I4P),              intent(out)   :: ntn(:)     !< Track of each new piece of the cell.
   real(R8P),    allocatable, intent(out)   :: meas(:,:)  !< Per-track measures [8, nt] (with_meas, nt >= 1).
   integer(I4P),              intent(out)   :: nev        !< Evaluations.
   type(cell_eval_t), allocatable :: ev(:)                 !< Evaluations: pose0, the wide intervals' middles, pose1.
   integer(I4P),      allocatable :: lg(:), lnx(:), nstart(:), lstart(:), lkey(:), lpiece(:), par(:), comp(:) !< Work.
   integer(I4P),      allocatable :: eref(:), ealt(:)      !< Per interval: its evaluation, the one after it (narrow).
   real(R8P),         allocatable :: b(:), bw(:)           !< Events, their half-widths.
   logical,           allocatable :: wide(:), small(:)     !< Per interval: evaluated; narrower than EVENT_GAP.
   real(R8P)                      :: x0, x1, y0, y1, sa, sb, cn(2, 4) !< Cell, interval, corners.
   logical                        :: fl(4)                 !< Corner flags.
   integer(I4P)                   :: np0, np1, nl, nb, ne, e, m, k, nk !< Work.
   logical                        :: cert, ok, unc         !< Certified events, measures assigned, uncertified contact.

   kind = 0_I4P ; nt = 0_I4P ; nev = 0_I4P
   x0 = x(i) ; x1 = x(i + 1) ; y0 = y(j) ; y1 = y(j + 1)
   np0 = k0%pieces(i, j) ; np1 = k1%pieces(i, j)
   if (.not. (k0%closed(i, j) .and. k1%closed(i, j))) then
      kind = -1_I4P ; return
   endif
   if (size(lseg) == 0) then ! no wall near the cell during the step: the same state throughout
      if (np0 == np1 .and. np0 <= 1_I4P) then
         nt = np0 ; otr = 1_I4P ; ntn = 1_I4P
      else
         kind = -1_I4P
      endif
      return
   endif
   call map_local(mc, lseg, nl, lg, lnx)
   call map_events(mc, x0, x1, y0, y1, nl, lg, lnx, b, bw, nb, cert, unc)
   if (.not. cert) then
      kind = -1_I4P ; return
   endif
   if (unc) then ! an uncertified contact: a state the evaluations cannot decide (issue #102)
      if (.not. with_meas) mc%contacts = mc%contacts + 1_I4P
      kind = -1_I4P ; return
   endif
   if (nb == 0_I4P .and. np0 == np1 .and. np0 <= 1_I4P) then ! one state throughout, at most one piece
      nt = np0 ; otr = 1_I4P ; ntn = 1_I4P ; return
   endif
   kind = 1_I4P
   cn = reshape([x0, y0, x1, y0, x1, y1, x0, y1], [2, 4])
   ! the intervals: a wide one has its evaluation (2 .. ne - 1); a narrow one refers to those before and after it
   allocate(eref(nb + 1), ealt(nb + 1), wide(nb + 1), small(nb + 1))
   ne = 1_I4P
   do m = 1_I4P, nb + 1_I4P
      call map_interval(b, nb, m, sa, sb)
      wide(m) = interval_wide(b, bw, nb, m) ; small(m) = sb - sa < EVENT_GAP
      if (wide(m)) then
         ne = ne + 1_I4P ; eref(m) = ne ; ealt(m) = ne
      else
         eref(m) = ne ; ealt(m) = 0_I4P
      endif
   enddo
   if (.not. with_meas .and. any(wide .and. small)) mc%resolved = mc%resolved + 1_I4P
   ne = ne + 1_I4P ; nev = ne
   e = ne
   do m = nb + 1_I4P, 1_I4P, -1_I4P
      if (ealt(m) == 0_I4P) then
         ealt(m) = e
      else
         e = eref(m)
      endif
   enddo
   allocate(ev(ne))
   fl = [k0%node_in(i, j), k0%node_in(i + 1, j), k0%node_in(i + 1, j + 1), k0%node_in(i, j + 1)]
   call map_eval(mc, 0_I4P, 0._R8P, x0, x1, y0, y1, lg, lnx, fl, ev(1))
   e = 1_I4P
   do m = 1_I4P, nb + 1_I4P
      call map_interval(b, nb, m, sa, sb)
      if (.not. wide(m)) cycle
      e = e + 1_I4P
      do k = 1_I4P, 4_I4P
         fl(k) = map_corner_in(mc, curve, c0, 0.5_R8P * (sa + sb), cn(:, k))
      enddo
      call map_eval(mc, 2_I4P, 0.5_R8P * (sa + sb), x0, x1, y0, y1, lg, lnx, fl, ev(e))
   enddo
   fl = [k1%node_in(i, j), k1%node_in(i + 1, j), k1%node_in(i + 1, j + 1), k1%node_in(i, j + 1)]
   call map_eval(mc, 1_I4P, 1._R8P, x0, x1, y0, y1, lg, lnx, fl, ev(ne))
   ! every evaluation closed; the end evaluations are the end-pose cuts' pieces, bit for bit
   do e = 1_I4P, ne
      if (.not. ev(e)%closed) kind = -1_I4P
   enddo
   if (ev(1)%np /= np0 .or. ev(ne)%np /= np1) kind = -1_I4P
   if (kind == -1_I4P) return
   do k = 1_I4P, np0
      if (ev(1)%parea(k) /= k0%piece_area(k0%first_piece(i, j) + k - 1_I4P)) kind = -1_I4P
   enddo
   do k = 1_I4P, np1
      if (ev(ne)%parea(k) /= k1%piece_area(k1%first_piece(i, j) + k - 1_I4P)) kind = -1_I4P
   enddo
   if (kind == -1_I4P) return
   ! tracks: pieces of consecutive evaluations sharing a feature
   allocate(nstart(ne + 1), lstart(ne + 1))
   nstart(1) = 1_I4P ; lstart(1) = 1_I4P
   do e = 1_I4P, ne
      nstart(e + 1) = nstart(e) + ev(e)%np ; lstart(e + 1) = lstart(e) + ev(e)%nk
   enddo
   nk = lstart(ne + 1) - 1_I4P
   allocate(lkey(nk), lpiece(nk), par(nstart(ne + 1) - 1_I4P), comp(nstart(ne + 1) - 1_I4P))
   do e = 1_I4P, ne
      lkey(lstart(e):lstart(e + 1) - 1_I4P) = ev(e)%key ; lpiece(lstart(e):lstart(e + 1) - 1_I4P) = ev(e)%kpiece
   enddo
   call track_components(ne, nstart, lstart, nk, lkey, lpiece, par, comp, nt)
   otr = comp(nstart(1):nstart(2) - 1_I4P) ; ntn = comp(nstart(ne):nstart(ne + 1) - 1_I4P)
   if (.not. (with_meas .and. nt >= 1_I4P)) return
   call map_measures(mc, curve, c0, x0, x1, y0, y1, lg, lnx, b, nb, eref, ealt, small, ev, nstart, comp, nt, meas, ok)
   if (.not. ok) kind = -1_I4P
   endsubroutine map_cell

   subroutine map_local(mc, lseg, nl, lg, lnx)
   !< The local curve of a cell: the end vertices of its segments in increasing order (lg, as the whole curve's indexes)
   !< and the next vertex of each in that numbering (0 if its segment is not local), so that `walk_all` walks the local
   !< segments in the whole curve's order.
   type(map_ctx_t),           intent(inout) :: mc      !< State.
   integer(I4P),              intent(in)    :: lseg(:) !< Local segments.
   integer(I4P),              intent(out)   :: nl      !< Local vertices.
   integer(I4P), allocatable, intent(out)   :: lg(:)   !< Their indexes in the curve.
   integer(I4P), allocatable, intent(out)   :: lnx(:)  !< Next local vertex.
   integer(I4P), allocatable                :: tmp(:)  !< Work.
   integer(I4P)                             :: k, e, g, a, h !< Work.

   mc%tick = mc%tick + 1_I4P
   allocate(tmp(2 * size(lseg))) ; nl = 0_I4P
   do k = 1_I4P, size(lseg, kind=I4P)
      do e = 1_I4P, 2_I4P
         g = lseg(k) ; if (e == 2_I4P) g = mc%nxt(lseg(k))
         if (mc%stamp(g) == mc%tick) cycle
         mc%stamp(g) = mc%tick ; nl = nl + 1_I4P ; tmp(nl) = g
      enddo
   enddo
   do k = 2_I4P, nl ! insertion sort (a few vertices)
      a = tmp(k) ; h = k - 1_I4P
      do while (h >= 1_I4P)
         if (tmp(h) <= a) exit
         tmp(h + 1_I4P) = tmp(h) ; h = h - 1_I4P
      enddo
      tmp(h + 1_I4P) = a
   enddo
   lg = tmp(1:nl)
   allocate(lnx(nl)) ; lnx = 0_I4P
   do k = 1_I4P, nl
      mc%loc(lg(k)) = k
   enddo
   do k = 1_I4P, size(lseg, kind=I4P)
      lnx(mc%loc(lseg(k))) = mc%loc(mc%nxt(lseg(k)))
   enddo
   endsubroutine map_local

   subroutine map_events(mc, x0, x1, y0, y1, nl, lg, lnx, b, bw, nb, cert, unc)
   !< The events of the cell [x0, x1] x [y0, y1] in (0, 1), sorted, each once: a local vertex crossing one of its four
   !< face lines inside the face, a corner crossing a local segment. Without rotation, the formulas of
   !< `faces_translation_2d` on the pose0 copy (`collect`, `node_breakpoints`); with a rotation, `line_events_2d` and
   !< `node_events_2d` on the local curve (`isolate_zeros`; cert false if one is not certified; unc true if an event of
   !< the cell is a contact within the round-off of its function, not certified by a sign change). bw (issue #104):
   !< each event's half-width, the instants it may stand for in double precision (`EVENT_RES` on the explicit formulas,
   !< `isolate_zeros`' width on the screw ones; an exact duplicate keeps the larger).
   type(map_ctx_t),        intent(in)  :: mc             !< State.
   real(R8P),              intent(in)  :: x0, x1, y0, y1 !< Cell.
   integer(I4P),           intent(in)  :: nl             !< Local vertices.
   integer(I4P),           intent(in)  :: lg(:), lnx(:)  !< Local vertices, next.
   real(R8P), allocatable, intent(out) :: b(:)           !< Events.
   real(R8P), allocatable, intent(out) :: bw(:)          !< Their half-widths (issue #104).
   integer(I4P),           intent(out) :: nb             !< Their number.
   logical,                intent(out) :: cert           !< Every event certified.
   logical,                intent(out) :: unc            !< An event is an uncertified contact (tangency).
   type(rotation_2d_t)                 :: lc             !< Local rotation state.
   type(events_2d_t)                   :: ev             !< Events of a line or a node.
   real(R8P), allocatable              :: tmp(:), tw(:), bn(:), bnw(:) !< Work.
   real(R8P)                           :: cn(2, 4), lin(4), lo(4), hi(4), sk, w !< Corners, face lines, spans.
   integer(I4P)                        :: k, f, g, q, nax !< Work.

   cert = .true. ; unc = .false. ; nb = 0_I4P
   cn = reshape([x0, y0, x1, y0, x1, y1, x0, y1], [2, 4])
   lin = [x0, x1, y0, y1] ; lo = [y0, y0, x0, x0] ; hi = [y1, y1, x1, x1]
   if (mc%dth == 0._R8P) then
      allocate(tmp(8_I4P * nl + 4_I4P), tw(8_I4P * nl + 4_I4P)) ! <= 4 line events per vertex, <= 1 node event per
                                                                ! segment and corner
      do k = 1_I4P, nl
         g = lg(k)
         do f = 1_I4P, 4_I4P
            nax = 1_I4P ; if (f > 2_I4P) nax = 2_I4P
            if (mc%d(nax) == 0._R8P) cycle
            if (nax == 1_I4P) then
               sk = (lin(f) - mc%u0(g)) / mc%d(1) ; w = mc%v0(g) + sk * mc%d(2)
            else
               sk = (lin(f) - mc%v0(g)) / mc%d(2) ; w = mc%u0(g) + sk * mc%d(1)
            endif
            if (sk > 0._R8P .and. sk < 1._R8P .and. w >= lo(f) .and. w <= hi(f)) then
               nb = nb + 1_I4P ; tmp(nb) = sk
               if (nax == 1_I4P) then
                  tw(nb) = EVENT_RES * (sk + (abs(lin(f)) + abs(mc%u0(g))) / abs(mc%d(1)))
               else
                  tw(nb) = EVENT_RES * (sk + (abs(lin(f)) + abs(mc%v0(g))) / abs(mc%d(2)))
               endif
            endif
         enddo
      enddo
      do f = 1_I4P, 4_I4P
         call node_breakpoints(cn(:, f), mc%d, mc%u0(lg), mc%v0(lg), lnx, bn, bnw)
         do q = 1_I4P, size(bn, kind=I4P)
            nb = nb + 1_I4P ; tmp(nb) = bn(q) ; tw(nb) = bnw(q)
         enddo
      enddo
   else
      lc%nv = nl ; lc%vu = mc%bu(lg) ; lc%vv = mc%bv(lg) ; lc%nxt = lnx
      lc%th0 = mc%th0 ; lc%dth = mc%dth ; lc%t0 = mc%pose0%t ; lc%dt = mc%d
      allocate(tmp(16), tw(16))
      do f = 1_I4P, 4_I4P
         nax = 1_I4P ; if (f > 2_I4P) nax = 2_I4P
         if (allocated(ev%z)) deallocate(ev%z, ev%w)
         call line_events_2d(lc, nax, lin(f), ev)
         if (.not. ev%cert) cert = .false.
         do q = 1_I4P, ev%n
            if (ev%w(q) >= lo(f) .and. ev%w(q) <= hi(f)) then
               call push_pair(tmp, tw, nb, ev%z(q), ev%r(q))
               if (ev%u(q)) unc = .true.
            endif
         enddo
      enddo
      do f = 1_I4P, 4_I4P
         if (allocated(ev%z)) deallocate(ev%z, ev%w)
         call node_events_2d(lc, cn(:, f), ev)
         if (.not. ev%cert) cert = .false.
         do q = 1_I4P, ev%n
            call push_pair(tmp, tw, nb, ev%z(q), ev%r(q))
            if (ev%u(q)) unc = .true.
         enddo
      enddo
   endif
   call sort_unique_w(tmp, tw, nb)
   b = tmp(1:nb) ; bw = tw(1:nb)
   endsubroutine map_events

   pure subroutine push_value(a, n, z)
   !< Append z to a(1:n) (geometric growth, no self-referencing constructor).
   real(R8P), allocatable, intent(inout) :: a(:) !< Values.
   integer(I4P),           intent(inout) :: n    !< Count.
   real(R8P),              intent(in)    :: z    !< Value.
   real(R8P), allocatable                :: t(:) !< Grown buffer.

   if (n == size(a, kind=I4P)) then
      allocate(t(2_I4P * n + 1_I4P)) ; t(1:n) = a(1:n) ; call move_alloc(t, a)
   endif
   n = n + 1_I4P ; a(n) = z
   endsubroutine push_value

   function map_corner_in(mc, curve, c0, s, p) result(inside)
   !< The corner p inside the solid at instant s: the parity of the curve's crossings left of p on its horizontal line
   !< (`line_crossings`, exact with the symbolic tie-break), p taken to the frame of the curve that the path moves:
   !< the pose0 copy shifted back by s d on the translation path, the body frame (`to_body` of the screw pose) else.
   type(map_ctx_t),    intent(in) :: mc        !< State.
   type(curve_object), intent(in) :: curve, c0 !< Curve (body frame), copy at pose0.
   real(R8P),          intent(in) :: s, p(2)   !< Instant, corner.
   logical                        :: inside    !< Inside the solid.
   real(R8P),    allocatable      :: t(:)      !< Crossings.
   integer(I4P), allocatable      :: g(:)      !< Segments.
   logical,      allocatable      :: e(:)      !< Entering flags.
   real(R8P)                      :: q(2)      !< Point in the curve's frame.
   integer(I4P)                   :: n         !< Crossings.

   if (mc%path == MAP_TRANSLATION) then
      q = p - s * mc%d
      call c0%line_crossings(1_I4P, q(2), min(q(1), mc%bmin0(1)) - 1._R8P, q(1), t, g, e, n)
   else
      q = to_body(pose_interpolate(mc%pose0, mc%pose1, s), p)
      call curve%line_crossings(1_I4P, q(2), min(q(1), mc%bmin(1)) - 1._R8P, q(1), t, g, e, n)
   endif
   inside = mod(n, 2_I4P) == 1_I4P
   endfunction map_corner_in

   subroutine map_eval(mc, kind, s, x0, x1, y0, y1, lg, lnx, flags, e)
   !< The pieces of the cell at one instant with their features (`cut_cell_2d_pieces` on the local curve): kind 0 the
   !< pose0 copy, 1 the pose1 copy, 2 the instant s (the pose0 copy moved by s d on the translation path, `to_world` of
   !< the screw pose else, as `place_for_crossings`). Features: the fluid corners (keys -4..-1) then the wall segments
   !< (their indexes in the curve), increasing.
   type(map_ctx_t),   intent(in)  :: mc             !< State.
   integer(I4P),      intent(in)  :: kind           !< 0 pose0, 1 pose1, 2 instant s.
   real(R8P),         intent(in)  :: s              !< Instant (kind 2).
   real(R8P),         intent(in)  :: x0, x1, y0, y1 !< Cell.
   integer(I4P),      intent(in)  :: lg(:), lnx(:)  !< Local vertices, next.
   logical,           intent(in)  :: flags(4)       !< Corner flags.
   type(cell_eval_t), intent(out) :: e              !< Pieces and features.
   real(R8P),    allocatable :: lu(:), lv(:), parea(:) !< Placed vertices, areas.
   integer(I4P), allocatable :: wseg(:), wpiece(:)  !< Wall pieces.
   type(rigid_pose_2d_t)     :: ps                  !< Pose at s.
   real(R8P)                 :: w(2)                !< Vertex.
   integer(I4P)              :: nl, k, nw, cp(4)    !< Work.

   nl = size(lg, kind=I4P)
   allocate(lu(nl), lv(nl))
   select case (kind)
   case (0_I4P)
      lu = mc%u0(lg) ; lv = mc%v0(lg)
   case (1_I4P)
      lu = mc%u1(lg) ; lv = mc%v1(lg)
   case default
      if (mc%path == MAP_TRANSLATION) then
         do k = 1_I4P, nl
            lu(k) = mc%u0(lg(k)) + s * mc%d(1) ; lv(k) = mc%v0(lg(k)) + s * mc%d(2)
         enddo
      else
         ps = pose_interpolate(mc%pose0, mc%pose1, s)
         do k = 1_I4P, nl
            w = to_world(ps, [mc%bu(lg(k)), mc%bv(lg(k))]) ; lu(k) = w(1) ; lv(k) = w(2)
         enddo
      endif
   endselect
   call cut_cell_2d_pieces(lu, lv, lnx, x0, x1, y0, y1, flags, e%closed, e%np, parea, nw, wseg, wpiece, cp, e%plen, &
                           e%pmom)
   call move_alloc(parea, e%parea) ; e%fl = flags
   e%nk = count(cp > 0_I4P, kind=I4P) + count(wpiece > 0_I4P, kind=I4P)
   allocate(e%key(e%nk), e%kpiece(e%nk))
   e%nk = 0_I4P
   do k = 4_I4P, 1_I4P, -1_I4P
      if (cp(k) == 0_I4P) cycle
      e%nk = e%nk + 1_I4P ; e%key(e%nk) = -k ; e%kpiece(e%nk) = cp(k)
   enddo
   do k = 1_I4P, nw ! walk order: increasing segments, one wall piece each
      if (wpiece(k) == 0_I4P) cycle
      e%nk = e%nk + 1_I4P ; e%key(e%nk) = lg(wseg(k)) ; e%kpiece(e%nk) = wpiece(k)
   enddo
   endsubroutine map_eval

   subroutine map_measures(mc, curve, c0, x0, x1, y0, y1, lg, lnx, b, nb, eref, ealt, small, ev, nstart, comp, nt, meas, &
                           ok)
   !< Per-track time integrals over the step of the fluid lengths of the cell's four faces (meas(1:4), -u, +u, -v, +v)
   !< and of the rigid-field outflows through them (meas(5:8)), interval by interval between the cell's events: inside
   !< an interval the pieces keep their features, so each piece of an evaluation at a quadrature node belongs to the
   !< track of the interval's evaluation that holds its features (`map_node`). Translation: two-point Gauss per
   !< interval (exact: the lengths are linear in s). Rotation: Gauss doubled from 2 nodes per interval until every
   !< track's values agree to round-off (the criteria of `integrate_face_2d`; mc%exact false if not at 64). The Gauss
   !< path integrates the same way, by its motion (two-point per interval without rotation, doubled with one): its nq
   !< nodes over the step would miss a track living between two events, or weight it heavily (issue #102 step 3:
   !< per-track residuals up to 0.39 of a cell). ok false if a piece at a node cannot be assigned.
   type(map_ctx_t),        intent(inout) :: mc             !< State.
   type(curve_object),     intent(in)    :: curve, c0      !< Curve (body frame), copy at pose0.
   real(R8P),              intent(in)    :: x0, x1, y0, y1 !< Cell.
   integer(I4P),           intent(in)    :: lg(:), lnx(:)  !< Local vertices, next.
   real(R8P),              intent(in)    :: b(:)           !< Events.
   integer(I4P),           intent(in)    :: nb             !< Their number.
   integer(I4P),           intent(in)    :: eref(:), ealt(:) !< Evaluations of each interval (see `map_cell`).
   logical,                intent(in)    :: small(:)       !< Intervals narrower than EVENT_GAP (nodes may be skipped).
   type(cell_eval_t),      intent(in)    :: ev(:)          !< Evaluations.
   integer(I4P),           intent(in)    :: nstart(:)      !< First node of each evaluation.
   integer(I4P),           intent(in)    :: comp(:)        !< Track of each node.
   integer(I4P),           intent(in)    :: nt             !< Tracks.
   real(R8P), allocatable, intent(out)   :: meas(:,:)      !< Measures [8, nt].
   logical,                intent(out)   :: ok             !< Every piece assigned.
   real(R8P), parameter   :: G = 0.5_R8P / sqrt(3._R8P)    !< Two-point Gauss offset on a unit interval.
   real(R8P), allocatable :: a1(:,:), a2(:,:), sq(:), wq(:) !< Work.
   real(R8P)              :: sa, sb, h, scl(4), sco(4), sno, vel, tol !< Work.
   integer(I4P)           :: m, ng, q, t, f, r, a      !< Work.
   logical                :: conv, rot                     !< Converged, a rotation.

   allocate(meas(8, nt)) ; meas = 0._R8P ; ok = .true.
   rot = mc%dth /= 0._R8P
   if (.not. rot) then
      do m = 1_I4P, nb + 1_I4P
         call map_interval(b, nb, m, sa, sb) ; h = sb - sa
         if (h <= 0._R8P) cycle
         r = eref(m) ; a = ealt(m)
         call map_node(mc, curve, c0, 0.5_R8P * (sa + sb) - G * h, 0.5_R8P * h, x0, x1, y0, y1, lg, lnx, ev(r), nstart(r), &
                       ev(a), nstart(a), r /= a .or. small(m), comp, meas, ok)
         call map_node(mc, curve, c0, 0.5_R8P * (sa + sb) + G * h, 0.5_R8P * h, x0, x1, y0, y1, lg, lnx, ev(r), nstart(r), &
                       ev(a), nstart(a), r /= a .or. small(m), comp, meas, ok)
      enddo
   else
      vel = abs(mc%dth) * mc%xm + sqrt(mc%d(1)**2 + mc%d(2)**2)
      scl = [y1 - y0, y1 - y0, x1 - x0, x1 - x0] ; sco = scl * vel ; sno = mc%xm * vel
      allocate(a1(8, nt), a2(8, nt))
      do m = 1_I4P, nb + 1_I4P
         call map_interval(b, nb, m, sa, sb) ; h = sb - sa
         if (h <= 0._R8P) cycle
         r = eref(m) ; a = ealt(m)
         ng = 2_I4P ; a1 = 0._R8P
         call gauss_legendre_01(ng, sq, wq)
         do q = 1_I4P, ng
            call map_node(mc, curve, c0, sa + h * sq(q), h * wq(q), x0, x1, y0, y1, lg, lnx, ev(r), nstart(r), ev(a), &
                          nstart(a), r /= a .or. small(m), comp, a1, ok)
         enddo
         do
            ng = 2_I4P * ng ; a2 = 0._R8P
            call gauss_legendre_01(ng, sq, wq)
            do q = 1_I4P, ng
               call map_node(mc, curve, c0, sa + h * sq(q), h * wq(q), x0, x1, y0, y1, lg, lnx, ev(r), nstart(r), ev(a), &
                             nstart(a), r /= a .or. small(m), comp, a2, ok)
            enddo
            if (.not. ok) exit
            tol = 64._R8P * epsilon(1._R8P) * max(h, FLOOR_H)
            conv = .true.
            do t = 1_I4P, nt
               do f = 1_I4P, 4_I4P
                  if (abs(a2(f, t) - a1(f, t)) > tol * scl(f) .or. abs(a2(4 + f, t) - a1(4 + f, t)) > tol * sco(f)) &
                     conv = .false.
               enddo
            enddo
            if (conv) exit
            if (ng >= 64_I4P) then ! the round-off of one evaluation (see FLOOR_H)
               do t = 1_I4P, nt
                  do f = 1_I4P, 4_I4P
                     if (max(abs(a2(f, t) - a1(f, t)) / mc%xm, abs(a2(4 + f, t) - a1(4 + f, t)) / sno) > &
                         64._R8P * epsilon(1._R8P)) mc%exact = .false.
                  enddo
               enddo
               exit
            endif
            a1 = a2
         enddo
         meas = meas + a2
      enddo
   endif
   endsubroutine map_measures

   pure logical function interval_wide(b, bw, nb, m) result(wide)
   !< Interval m (see `map_interval`) holds a state of its own (issue #104): at least EVENT_GAP wide, or wider than the
   !< sum of its two events' half-widths (`bw`; 0 at the step's ends), so the two are distinct crossings and the state
   !< between them is resolved at the interval's middle.
   real(R8P),    intent(in) :: b(:), bw(:) !< Events, their half-widths.
   integer(I4P), intent(in) :: nb, m       !< Events, interval.
   real(R8P)                :: sa, sb, wa, wb !< Ends, their half-widths.

   call map_interval(b, nb, m, sa, sb)
   wa = 0._R8P ; if (m > 1_I4P) wa = bw(m - 1_I4P)
   wb = 0._R8P ; if (m <= nb) wb = bw(m)
   wide = sb - sa >= EVENT_GAP .or. sb - sa > wa + wb
   endfunction interval_wide

   pure subroutine map_interval(b, nb, m, sa, sb)
   !< Interval m (1 .. nb + 1) between the events b(1:nb) of a cell, in [0, 1].
   real(R8P),    intent(in)  :: b(:)   !< Events.
   integer(I4P), intent(in)  :: nb, m  !< Events, interval.
   real(R8P),    intent(out) :: sa, sb !< Ends.

   sa = 0._R8P ; if (m > 1_I4P) sa = b(m - 1_I4P)
   sb = 1._R8P ; if (m <= nb) sb = b(m)
   endsubroutine map_interval

   subroutine map_node(mc, curve, c0, s, w, x0, x1, y0, y1, lg, lnx, evr, nr, eva, na, narrow, comp, acc, ok)
   !< One quadrature node (instant s, weight w): the pieces at s (corner flags of the interval's evaluation, computed
   !< at s in a narrow interval), each added to a track
   !< with w times its face lengths and rigid-field outflows. A piece takes the track of the pieces of the interval's
   !< evaluation evr (first node nr) that hold its features (wall segments, fluid corners); a narrow interval also tries
   !< the evaluation after it, eva (na). The features must point to one track; a piece none of whose features is known
   !< (born within round-off of a nearby event) must be negligible (face lengths below 1e-9 of the cell side), else ok
   !< is false. In a narrow interval (narrow: width below EVENT_GAP, the node within that of an event) a node whose
   !< cell does not close, or a piece that cannot be assigned, is skipped (mc%skipped): its weight is below EVENT_GAP.
   !< Outflow through a face u = a: -dth (m - t_v(s) l) + d_u l; through a face v = a: dth (m - t_u(s) l) + d_v l (l the
   !< fluid length, m its first moment along the face, t(s) = t0 + s d: as `faces_rotation_2d`).
   type(map_ctx_t),    intent(inout) :: mc             !< State.
   type(curve_object), intent(in)    :: curve, c0      !< Curve (body frame), copy at pose0.
   real(R8P),          intent(in)    :: s, w           !< Instant, weight.
   real(R8P),          intent(in)    :: x0, x1, y0, y1 !< Cell.
   integer(I4P),       intent(in)    :: lg(:), lnx(:)  !< Local vertices, next.
   type(cell_eval_t),  intent(in)    :: evr, eva       !< The interval's evaluation, the one after it.
   integer(I4P),       intent(in)    :: nr, na         !< Their first nodes.
   logical,            intent(in)    :: narrow         !< A narrow interval (evr and eva differ; issue #102).
   integer(I4P),       intent(in)    :: comp(:)        !< Track of each node.
   real(R8P),          intent(inout) :: acc(:,:)       !< Sums [8, nt].
   logical,            intent(inout) :: ok             !< Every piece assigned.
   type(cell_eval_t)                 :: e              !< Pieces at s.
   real(R8P)                         :: tq(2), l, mm, o, cn(2, 4) !< Work.
   logical                           :: fl(4)          !< Corner flags at s.
   integer(I4P)                      :: p, t, f, k     !< Work.

   if (.not. ok) return
   if (na == nr) then ! a wide interval: no corner crosses a segment inside it, its evaluation's flags hold
      fl = evr%fl
   else
      cn = reshape([x0, y0, x1, y0, x1, y1, x0, y1], [2, 4])
      do k = 1_I4P, 4_I4P
         fl(k) = map_corner_in(mc, curve, c0, s, cn(:, k))
      enddo
   endif
   call map_eval(mc, 2_I4P, s, x0, x1, y0, y1, lg, lnx, fl, e)
   if (.not. e%closed) then
      if (narrow) then
         mc%skipped = mc%skipped + 1_I4P
      else
         ok = .false.
      endif
      return
   endif
   tq = mc%pose0%t + s * mc%d
   do p = 1_I4P, e%np
      t = piece_track(e, p, evr, nr, comp)
      if (t == 0_I4P .and. narrow) t = piece_track(e, p, eva, na, comp)
      if (t <= 0_I4P) then ! no known feature (negligible, or a missed event) or two tracks
         if (t == 0_I4P .and. maxval(e%plen(:, p)) <= 1e-9_R8P * max(x1 - x0, y1 - y0)) cycle
         if (narrow) then
            mc%skipped = mc%skipped + 1_I4P ; cycle
         endif
         ok = .false. ; return
      endif
      do f = 1_I4P, 4_I4P
         l = e%plen(f, p) ; mm = e%pmom(f, p)
         if (f <= 2_I4P) then
            o = -mc%dth * (mm - tq(2) * l) + mc%d(1) * l
         else
            o = mc%dth * (mm - tq(1) * l) + mc%d(2) * l
         endif
         acc(f, t) = acc(f, t) + w * l ; acc(4 + f, t) = acc(4 + f, t) + w * o
      enddo
   enddo
   endsubroutine map_node

   pure integer(I4P) function piece_track(e, p, evr, nr, comp) result(t)
   !< The track of piece p of e through its features found in evr (first node nr): 0 if none is found, -1 if they point
   !< to two tracks.
   type(cell_eval_t), intent(in) :: e, evr  !< Evaluation, reference evaluation.
   integer(I4P),      intent(in) :: p, nr   !< Piece, first node of evr.
   integer(I4P),      intent(in) :: comp(:) !< Track of each node.
   integer(I4P)                  :: k, q, lo, hi, mid, tk !< Work.

   t = 0_I4P
   do k = 1_I4P, e%nk
      if (e%kpiece(k) /= p) cycle
      q = 0_I4P ; lo = 1_I4P ; hi = evr%nk ! the feature in evr (keys increasing)
      do while (lo <= hi)
         mid = (lo + hi) / 2_I4P
         if (evr%key(mid) == e%key(k)) then
            q = evr%kpiece(mid) ; exit
         elseif (evr%key(mid) < e%key(k)) then
            lo = mid + 1_I4P
         else
            hi = mid - 1_I4P
         endif
      enddo
      if (q == 0_I4P) cycle
      tk = comp(nr + q - 1_I4P)
      if (t == 0_I4P) then
         t = tk
      elseif (t /= tk) then
         t = -1_I4P ; return
      endif
   enddo
   endfunction piece_track

   subroutine piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, path, nq, map, all_meas, measured)
   !< The 3-D piece map of the step (issue #102 step 2; `piece_map_3d_t`), from the end-pose cuts k0, k1 of the copies
   !< of `surface` moved to pose0 (s0) and pose1; the 3-D twin of `piece_map_2d_build`.
   !<
   !< 1. Candidates: the cells a facet can meet during the step (the box of its vertices' paths: the segment of the
   !<    translation, or 17 samples of the screw path grown by its curvature bound), each with its facets (CSR,
   !<    increasing facet order).
   !< 2. Topology per cell (`map3_cell`): the events of the cell from its facets (`map3_events`: a mesh vertex crossing a
   !<    face plane inside the face, a mesh edge crossing a cell edge, a corner crossing a facet; explicit on a
   !<    translation, `isolate_zeros` under a rotation); a cell without events whose two end cuts have the same single
   !<    piece (or none) takes the trivial track. Otherwise its pieces are evaluated at pose0 (the pose0 copy, k0's
   !<    corner states and fractions: they must equal k0's pieces bit for bit, which fixes the numbering), at the middle
   !<    of every interval between events, and at pose1, and pieces of consecutive evaluations sharing a boundary
   !<    feature (a facet or a fluid corner) are one track (`track_components`).
   !< 3. Measures: one track takes the cell's `sw` values and the face outflows of the call (fox, foy, foz: the time
   !<    integral of the rigid field's flux in +axis through every grid face); several tracks integrate their face
   !<    parts on the cell's event intervals with the exact paths' rules, on the Gauss path too (`map3_measures`).
   !< 4. CSR fill in cell order.
   type(surface_stl_object), intent(in)  :: surface        !< Surface (body frame).
   type(surface_stl_object), intent(in)  :: s0             !< Copy at pose0 (`move`, the end-pose cut's).
   type(rigid_pose_t),       intent(in)  :: pose0, pose1   !< End poses.
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(in)  :: k0, k1         !< Cut cells at pose0 and pose1.
   type(swept_cells_3d_t),   intent(in)  :: sw             !< Swept geometry of the step.
   real(R8P),                intent(in)  :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Face outflows in +axis.
   integer(I4P),             intent(in)  :: path           !< MAP_TRANSLATION, MAP_ROTATION or MAP_GAUSS.
   integer(I4P),             intent(in)  :: nq             !< Gauss nodes (MAP_GAUSS).
   type(piece_map_3d_t),     intent(out) :: map            !< Piece map.
   logical,      optional,   intent(in)  :: all_meas     !< Integrate every evaluated cell (issue #104).
   logical, allocatable, optional, intent(out) :: measured(:) !< Cells integrated per event interval [nx ny nz].
   type(map3_ctx_t)          :: mc                         !< State.
   integer(I4P), allocatable :: cstart(:), cfac(:)         !< Facets of every cell.
   integer(I4P), allocatable :: ntr(:), ckind(:)           !< Tracks and kind (0 trivial, 1 evaluated, -1 untrusted).
   integer(I4P), allocatable :: otr(:), ntn(:)             !< Track (in its cell) of every old and new piece.
   integer(I4P), allocatable :: moff(:), oc(:), nc(:), fo(:), fn(:) !< Work.
   real(R8P),    allocatable :: meas(:,:), mloc(:,:)       !< Per-track measures of the multi-track cells.
   integer(I4P), allocatable :: so(:), sn(:)               !< Saved tracks of a cell's pieces (all_meas).
   integer(I4P)              :: nx, ny, nz, ncell, c, i, j, k, nt, kind, nev, tot, t, p, a0, a1, b0, b1 !< Work.
   logical                   :: allm                       !< all_meas.

   allm = .false. ; if (present(all_meas)) allm = all_meas
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   ncell = nx * ny * nz
   map%nx = nx ; map%ny = ny ; map%nz = nz ; map%np0 = k0%np ; map%np1 = k1%np
   call map3_context(surface, s0, pose0, pose1, x, y, z, path, nq, mc)
   call map3_candidates(mc, x, y, z, cstart, cfac)
   allocate(ntr(ncell), ckind(ncell), otr(k0%np), ntn(k1%np), moff(ncell))
   ! topology of every cell
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            a0 = k0%first_piece(i, j, k) ; a1 = a0 + max(0_I4P, k0%pieces(i, j, k)) - 1_I4P
            b0 = k1%first_piece(i, j, k) ; b1 = b0 + max(0_I4P, k1%pieces(i, j, k)) - 1_I4P
            call map3_cell(mc, surface, s0, i, j, k, x, y, z, k0, k1, cfac(cstart(c):cstart(c + 1) - 1_I4P), .false., &
                           kind, nt, otr(a0:a1), ntn(b0:b1), mloc, nev)
            ckind(c) = kind ; ntr(c) = nt
            map%evaluations = map%evaluations + nev
            if (nev > 0_I4P) map%evaluated_cells = map%evaluated_cells + 1_I4P
         enddo
      enddo
   enddo
   ! the multi-track cells again, with their per-track measures
   tot = 0_I4P ; moff = 0_I4P
   do c = 1_I4P, ncell
      if (ckind(c) == 1_I4P .and. (ntr(c) >= 2_I4P .or. allm)) then
         moff(c) = tot + 1_I4P ; tot = tot + ntr(c)
      endif
   enddo
   allocate(meas(12, tot))
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            if (moff(c) == 0_I4P) cycle
            a0 = k0%first_piece(i, j, k) ; a1 = a0 + max(0_I4P, k0%pieces(i, j, k)) - 1_I4P
            b0 = k1%first_piece(i, j, k) ; b1 = b0 + max(0_I4P, k1%pieces(i, j, k)) - 1_I4P
            so = otr(a0:a1) ; sn = ntn(b0:b1)
            call map3_cell(mc, surface, s0, i, j, k, x, y, z, k0, k1, cfac(cstart(c):cstart(c + 1) - 1_I4P), .true., &
                           kind, nt, otr(a0:a1), ntn(b0:b1), mloc, nev)
            if (kind == 1_I4P .and. nt == ntr(c)) then
               meas(:, moff(c):moff(c) + nt - 1_I4P) = mloc(:, 1:nt)
            elseif (ntr(c) == 1_I4P .and. allm) then ! one track (all_meas): keep it with the sw values
               otr(a0:a1) = so ; ntn(b0:b1) = sn ; moff(c) = 0_I4P
            else
               ckind(c) = -1_I4P ; moff(c) = 0_I4P
            endif
         enddo
      enddo
   enddo
   if (present(measured)) measured = moff > 0_I4P
   ! untrusted cells: one track with every piece
   allocate(map%tracks(nx, ny, nz), map%first_track(nx, ny, nz))
   nt = 0_I4P
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            if (ckind(c) == -1_I4P) then
               ntr(c) = 1_I4P
               a0 = k0%first_piece(i, j, k) ; a1 = a0 + max(0_I4P, k0%pieces(i, j, k)) - 1_I4P
               b0 = k1%first_piece(i, j, k) ; b1 = b0 + max(0_I4P, k1%pieces(i, j, k)) - 1_I4P
               otr(a0:a1) = 1_I4P ; ntn(b0:b1) = 1_I4P
               map%untrusted_cells = map%untrusted_cells + 1_I4P
            endif
            if (ntr(c) >= 2_I4P) map%multi_track_cells = map%multi_track_cells + 1_I4P
            map%tracks(i, j, k) = ntr(c) ; map%first_track(i, j, k) = nt + 1_I4P ; nt = nt + ntr(c)
         enddo
      enddo
   enddo
   map%nt = nt
   allocate(map%track_cell(3, nt), map%status(nt), map%old_start(nt + 1), map%new_start(nt + 1), map%old_piece(k0%np), &
            map%new_piece(k1%np), map%old_track(k0%np), map%new_track(k1%np), map%volume0(nt), map%volume1(nt), &
            map%faces(6, nt), map%outflow(6, nt), map%wall(3, nt), map%residual(nt), oc(nt), nc(nt), fo(nt), fn(nt))
   map%volume0 = 0._R8P ; map%volume1 = 0._R8P ; oc = 0_I4P ; nc = 0_I4P
   do p = 1_I4P, k0%np ! pieces in order: volumes summed in piece order
      t = map%first_track(k0%piece_cell(1, p), k0%piece_cell(2, p), k0%piece_cell(3, p)) + otr(p) - 1_I4P
      map%old_track(p) = t ; oc(t) = oc(t) + 1_I4P ; map%volume0(t) = map%volume0(t) + k0%piece_volume(p)
   enddo
   do p = 1_I4P, k1%np
      t = map%first_track(k1%piece_cell(1, p), k1%piece_cell(2, p), k1%piece_cell(3, p)) + ntn(p) - 1_I4P
      map%new_track(p) = t ; nc(t) = nc(t) + 1_I4P ; map%volume1(t) = map%volume1(t) + k1%piece_volume(p)
   enddo
   map%old_start(1) = 1_I4P ; map%new_start(1) = 1_I4P
   do t = 1_I4P, nt
      map%old_start(t + 1) = map%old_start(t) + oc(t) ; map%new_start(t + 1) = map%new_start(t) + nc(t)
      fo(t) = map%old_start(t) ; fn(t) = map%new_start(t)
   enddo
   do p = 1_I4P, k0%np
      t = map%old_track(p) ; map%old_piece(fo(t)) = p ; fo(t) = fo(t) + 1_I4P
   enddo
   do p = 1_I4P, k1%np
      t = map%new_track(p) ; map%new_piece(fn(t)) = p ; fn(t) = fn(t) + 1_I4P
   enddo
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            do t = map%first_track(i, j, k), map%first_track(i, j, k) + ntr(c) - 1_I4P
               map%track_cell(:, t) = [i, j, k]
               if (ckind(c) == -1_I4P) then
                  map%status(t) = TRACK_UNTRUSTED
               else
                  map%status(t) = track_status_of(oc(t), nc(t))
               endif
               map%n_status(map%status(t)) = map%n_status(map%status(t)) + 1_I4P
               if (moff(c) > 0_I4P) then ! several tracks: integrated per track
                  p = moff(c) + t - map%first_track(i, j, k)
                  map%faces(:, t) = meas(1:6, p) ; map%outflow(:, t) = meas(7:12, p)
                  map%wall(:, t) = -[meas(2, p) - meas(1, p), meas(4, p) - meas(3, p), meas(6, p) - meas(5, p)]
                  map%residual(t) = map%volume1(t) - map%volume0(t) + (meas(7, p) + meas(8, p) + meas(9, p) + &
                                    meas(10, p) + meas(11, p) + meas(12, p))
               else ! one track: the cell's values
                  map%faces(:, t) = [sw%face_x(i, j, k), sw%face_x(i + 1, j, k), sw%face_y(i, j, k), sw%face_y(i, j + 1, k), &
                                     sw%face_z(i, j, k), sw%face_z(i, j, k + 1)]
                  map%outflow(:, t) = [-fox(i, j, k), fox(i + 1, j, k), -foy(i, j, k), foy(i, j + 1, k), -foz(i, j, k), &
                                       foz(i, j, k + 1)]
                  map%wall(:, t) = sw%wall(:, i, j, k)
                  map%residual(t) = sw%residual(i, j, k)
               endif
            enddo
         enddo
      enddo
   enddo
   map%exact = mc%exact ; map%skipped_nodes = mc%skipped
   map%contact_cells = mc%contacts ; map%narrow_cells = mc%resolved
   endsubroutine piece_map_3d_build

   subroutine map3_context(surface, s0, pose0, pose1, x, y, z, path, nq, mc)
   !< The state of the 3-D piece map construction: facets (body, pose0 and pose1 copies), motion, vertex path boxes, mesh
   !< vertex and edge ids (exact equality of the body-frame points), work stamps.
   type(surface_stl_object), intent(in)  :: surface, s0       !< Surface (body frame), copy at pose0.
   type(rigid_pose_t),       intent(in)  :: pose0, pose1      !< End poses.
   real(R8P),                intent(in)  :: x(:), y(:), z(:)  !< Cell edges.
   integer(I4P),             intent(in)  :: path, nq          !< Path, Gauss nodes.
   type(map3_ctx_t),         intent(out) :: mc                !< State.
   type(surface_stl_object)    :: s1                          !< Copy at pose1.
   type(facet_object), pointer :: fp                          !< Facet.
   real(R8P),    allocatable   :: key(:,:)                    !< Sort keys.
   integer(I4P), allocatable   :: ord(:)                      !< Sorted order.
   real(R8P)                   :: xs(3), mrg, pad             !< Work.
   integer(I4P)                :: f, m, q, n, id              !< Work.

   mc%path = path ; mc%nq = nq ; mc%nf = surface%get_facets_number()
   allocate(mc%tb(3, 3, mc%nf), mc%t0(3, 3, mc%nf), mc%t1(3, 3, mc%nf), mc%w(3, 3, mc%nf), mc%vlo(3, 3, mc%nf), &
            mc%vhi(3, 3, mc%nf), mc%fv(3, mc%nf), mc%fe(3, mc%nf))
   s1 = surface ; call s1%move(pose1)
   do f = 1_I4P, mc%nf
      fp => surface%facet_at(f)
      do m = 1_I4P, 3_I4P
         mc%tb(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
      fp => s0%facet_at(f)
      do m = 1_I4P, 3_I4P
         mc%t0(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
      fp => s1%facet_at(f)
      do m = 1_I4P, 3_I4P
         mc%t1(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
   enddo
   mc%pose0 = pose0 ; mc%pose1 = pose1 ; mc%d = pose1%t - pose0%t
   call pose_log_3d(pose0, pose1, mc%kax, mc%ang) ; mc%om = mc%ang * mc%kax
   mc%rot = path /= MAP_TRANSLATION .and. .not. all(pose1%r == pose0%r)
   mc%xm = max(maxval(abs(x)), maxval(abs(y)), maxval(abs(z))) + norm2(pose0%t) + norm2(mc%d)
   mc%tol = 1e-12_R8P * (max(maxval(abs(x)), maxval(abs(y)), maxval(abs(z))) + 1._R8P)
   pad = 64._R8P * epsilon(1._R8P) * mc%xm + mc%tol
   do f = 1_I4P, mc%nf ! the box of every vertex path
      do m = 1_I4P, 3_I4P
         mc%w(:, m, f) = matmul(pose0%r, mc%tb(:, m, f))
         if (mc%rot) then ! 17 samples, |X''| <= ang^2 |w| (as `faces_rotation_3d`)
            mc%vlo(:, m, f) = huge(1._R8P) ; mc%vhi(:, m, f) = -huge(1._R8P)
            do q = 0_I4P, 16_I4P
               xs = screw_at(mc, mc%w(:, m, f), real(q, R8P) / 16._R8P)
               mc%vlo(:, m, f) = min(mc%vlo(:, m, f), xs) ; mc%vhi(:, m, f) = max(mc%vhi(:, m, f), xs)
            enddo
            mrg = mc%ang**2 * norm2(mc%w(:, m, f)) / (8._R8P * 16._R8P**2) + 1e-12_R8P * (norm2(mc%w(:, m, f)) + &
                  norm2(pose0%t) + norm2(mc%d) + 1._R8P)
         else
            mc%vlo(:, m, f) = min(mc%t0(:, m, f), mc%t0(:, m, f) + mc%d)
            mc%vhi(:, m, f) = max(mc%t0(:, m, f), mc%t0(:, m, f) + mc%d)
            mrg = 0._R8P
         endif
         mc%vlo(:, m, f) = mc%vlo(:, m, f) - mrg - pad ; mc%vhi(:, m, f) = mc%vhi(:, m, f) + mrg + pad
      enddo
   enddo
   ! mesh vertices: equal body-frame points (sorted lexicographically); mesh edges: equal vertex pairs
   n = 3_I4P * mc%nf
   allocate(key(3, n))
   do f = 1_I4P, mc%nf
      do m = 1_I4P, 3_I4P
         key(:, 3_I4P * (f - 1_I4P) + m) = mc%tb(:, m, f)
      enddo
   enddo
   call sort_lex(key, ord)
   id = 0_I4P
   do q = 1_I4P, n
      if (q == 1_I4P) then
         id = 1_I4P
      elseif (any(key(:, ord(q)) /= key(:, ord(q - 1_I4P)))) then
         id = id + 1_I4P
      endif
      mc%fv(mod(ord(q) - 1_I4P, 3_I4P) + 1_I4P, (ord(q) - 1_I4P) / 3_I4P + 1_I4P) = id
   enddo
   mc%nvx = id
   key = 0._R8P
   do f = 1_I4P, mc%nf
      do m = 1_I4P, 3_I4P
         q = 3_I4P * (f - 1_I4P) + m
         key(1, q) = real(min(mc%fv(m, f), mc%fv(mod(m, 3_I4P) + 1_I4P, f)), R8P)
         key(2, q) = real(max(mc%fv(m, f), mc%fv(mod(m, 3_I4P) + 1_I4P, f)), R8P)
      enddo
   enddo
   call sort_lex(key, ord)
   id = 0_I4P
   do q = 1_I4P, n
      if (q == 1_I4P) then
         id = 1_I4P
      elseif (any(key(:, ord(q)) /= key(:, ord(q - 1_I4P)))) then
         id = id + 1_I4P
      endif
      mc%fe(mod(ord(q) - 1_I4P, 3_I4P) + 1_I4P, (ord(q) - 1_I4P) / 3_I4P + 1_I4P) = id
   enddo
   mc%ned = id
   allocate(mc%vstamp(max(1_I4P, mc%nvx)), mc%estamp(max(1_I4P, mc%ned)))
   mc%vstamp = 0_I4P ; mc%estamp = 0_I4P ; mc%tick = 0_I4P ; mc%exact = .true. ; mc%skipped = 0_I4P
   endsubroutine map3_context

   subroutine sort_lex(key, ord)
   !< Stable bottom-up merge sort of the columns of key [3, n] in lexicographic order: key(:, ord(1)) <= ...
   real(R8P),                 intent(in)  :: key(:,:) !< Keys.
   integer(I4P), allocatable, intent(out) :: ord(:)   !< Order.
   integer(I4P), allocatable :: tmp(:)                !< Work.
   integer(I4P)              :: n, wdt, lo, mid, hi, a, b, o, r !< Work.
   logical                   :: le                    !< key(:, a) <= key(:, b).

   n = size(key, 2, kind=I4P)
   allocate(ord(n), tmp(n))
   do a = 1_I4P, n
      ord(a) = a
   enddo
   wdt = 1_I4P
   do while (wdt < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + wdt, n + 1_I4P) ; hi = min(lo + 2_I4P * wdt, n + 1_I4P)
         a = lo ; b = mid ; o = lo
         do while (a < mid .and. b < hi)
            le = .true.
            do r = 1_I4P, size(key, 1, kind=I4P)
               if (key(r, ord(a)) < key(r, ord(b))) exit
               if (key(r, ord(a)) > key(r, ord(b))) then
                  le = .false. ; exit
               endif
            enddo
            if (le) then
               tmp(o) = ord(a) ; a = a + 1_I4P
            else
               tmp(o) = ord(b) ; b = b + 1_I4P
            endif
            o = o + 1_I4P
         enddo
         do while (a < mid)
            tmp(o) = ord(a) ; a = a + 1_I4P ; o = o + 1_I4P
         enddo
         do while (b < hi)
            tmp(o) = ord(b) ; b = b + 1_I4P ; o = o + 1_I4P
         enddo
         lo = hi
      enddo
      ord = tmp ; wdt = 2_I4P * wdt
   enddo
   endsubroutine sort_lex

   pure function screw_at(mc, wv, s) result(xv)
   !< Position at instant s of the body point with R0 q = wv under the screw motion of the step (as `faces_rotation_3d`).
   type(map3_ctx_t), intent(in) :: mc    !< State.
   real(R8P),        intent(in) :: wv(3) !< R0 q.
   real(R8P),        intent(in) :: s     !< Instant.
   real(R8P)                    :: xv(3), cv(3) !< Position, work.

   cv = mc%kax * dot_product(mc%kax, wv)
   xv = (wv - cv) * cos(s * mc%ang) + cross(mc%kax, wv) * sin(s * mc%ang) + cv + mc%pose0%t + s * mc%d
   endfunction screw_at

   pure subroutine rodrigues_k(kax, wv, a3, b3, c3)
   !< Rot(phi) wv = a3 cos phi + b3 sin phi + c3 about the unit axis kax.
   real(R8P), intent(in)  :: kax(3), wv(3)  !< Axis, vector.
   real(R8P), intent(out) :: a3(3), b3(3), c3(3) !< Coefficients.

   c3 = kax * dot_product(kax, wv) ; a3 = wv - c3 ; b3 = cross(kax, wv)
   endsubroutine rodrigues_k

   subroutine map3_candidates(mc, x, y, z, cstart, cfac)
   !< The cells a facet can meet during the step, with their facets (CSR over the cells c = i + nx ((j - 1) + ny (k - 1)),
   !< the facets of a cell in increasing order): the box of the facet's vertex paths, the cells meeting it (closed).
   !< Count, then fill.
   type(map3_ctx_t),          intent(in)  :: mc               !< State.
   real(R8P),                 intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P), allocatable, intent(out) :: cstart(:)        !< [ncell + 1].
   integer(I4P), allocatable, intent(out) :: cfac(:)          !< Facets.
   integer(I4P), allocatable :: fill(:)                       !< Fill pointers.
   real(R8P)                 :: lo(3), hi(3)                  !< Box.
   integer(I4P)              :: nx, ny, nz, ncell, pass, f, q, i1, i2, j1, j2, l1, l2, ii, jj, kk, c !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   ncell = nx * ny * nz
   allocate(cstart(ncell + 1), fill(ncell))
   cstart = 0_I4P
   do pass = 1_I4P, 2_I4P
      if (pass == 2_I4P) then
         q = 1_I4P
         do c = 1_I4P, ncell
            fill(c) = q ; q = q + cstart(c) ; cstart(c) = fill(c)
         enddo
         cstart(ncell + 1) = q
         allocate(cfac(q - 1_I4P))
      endif
      do f = 1_I4P, mc%nf
         lo = min(mc%vlo(:, 1, f), mc%vlo(:, 2, f), mc%vlo(:, 3, f))
         hi = max(mc%vhi(:, 1, f), mc%vhi(:, 2, f), mc%vhi(:, 3, f))
         if (.not. (all(lo <= hi))) then ! a NaN coordinate: every cell
            lo = -huge(1._R8P) ; hi = huge(1._R8P)
         endif
         call index_range(x, lo(1), hi(1), i1, i2) ; call index_range(y, lo(2), hi(2), j1, j2)
         call index_range(z, lo(3), hi(3), l1, l2)
         i1 = max(1_I4P, i1 - 1_I4P) ; i2 = min(nx, i2) ; j1 = max(1_I4P, j1 - 1_I4P) ; j2 = min(ny, j2)
         l1 = max(1_I4P, l1 - 1_I4P) ; l2 = min(nz, l2)
         do kk = l1, l2
            do jj = j1, j2
               do ii = i1, i2
                  c = ii + nx * ((jj - 1_I4P) + ny * (kk - 1_I4P))
                  if (pass == 1_I4P) then
                     cstart(c) = cstart(c) + 1_I4P
                  else
                     cfac(fill(c)) = f ; fill(c) = fill(c) + 1_I4P
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine map3_candidates

   subroutine map3_cell(mc, surface, s0, i, j, k, x, y, z, k0, k1, lf, with_meas, kind, nt, otr, ntn, meas, nev)
   !< The tracks of cell (i, j, k) (see `piece_map_3d_build`, step 2): kind 0 (trivial: no evaluation), 1 (evaluated) or
   !< -1 (untrusted); nt tracks; the track of each of its old pieces (otr, in k0's order) and new pieces (ntn); with
   !< `with_meas` and nt >= 1, the per-track time integrals meas(1:6) of the faces' fluid areas and meas(7:12) of their
   !< outward rigid-field outflows (`map3_measures`; kind -1 if a face part could not be assigned). Intervals between
   !< events narrower than EVENT_GAP are not evaluated (their events act as one), as `map_cell`; a cell with an
   !< uncertified contact (a corner grazing a facet, a vertex grazing a face plane or an edge a cell edge, under a
   !< rotation) is untrusted.
   type(map3_ctx_t),          intent(inout) :: mc               !< State.
   type(surface_stl_object),  intent(in)    :: surface, s0      !< Surface (body frame), copy at pose0.
   integer(I4P),              intent(in)    :: i, j, k          !< Cell.
   real(R8P),                 intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),      intent(in)    :: k0, k1           !< End-pose cuts.
   integer(I4P),              intent(in)    :: lf(:)            !< Facets of the cell (increasing).
   logical,                   intent(in)    :: with_meas        !< Integrate the per-track measures.
   integer(I4P),              intent(out)   :: kind             !< 0 trivial, 1 evaluated, -1 untrusted.
   integer(I4P),              intent(out)   :: nt               !< Tracks.
   integer(I4P),              intent(out)   :: otr(:)           !< Track of each old piece of the cell.
   integer(I4P),              intent(out)   :: ntn(:)           !< Track of each new piece of the cell.
   real(R8P),    allocatable, intent(out)   :: meas(:,:)        !< Per-track measures [12, nt] (with_meas, nt >= 1).
   integer(I4P),              intent(out)   :: nev              !< Evaluations.
   type(cell_eval3_t), allocatable :: ev(:)                     !< Evaluations: pose0, the wide intervals' middles, pose1.
   integer(I4P),       allocatable :: nstart(:), lstart(:), lkey(:), lpiece(:), par(:), comp(:) !< Work.
   integer(I4P),       allocatable :: eref(:), ealt(:)          !< Per interval: its evaluation, the one after it (narrow).
   real(R8P),          allocatable :: b(:), bw(:)               !< Events, their half-widths.
   logical,            allocatable :: wide(:), small(:)         !< Per interval: evaluated; narrower than EVENT_GAP.
   real(R8P)                       :: c0(3), c1(3), sa, sb, fr(6) !< Cell, interval, face fractions.
   integer(I4P)                    :: np0, np1, nb, ne, e, m, q, nk, cs(8) !< Work.
   logical                         :: cert, unc                 !< Certified events, uncertified contact.

   kind = 0_I4P ; nt = 0_I4P ; nev = 0_I4P
   c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
   np0 = k0%pieces(i, j, k) ; np1 = k1%pieces(i, j, k)
   if (np0 < 0_I4P .or. np1 < 0_I4P) then
      kind = -1_I4P ; return
   endif
   if (size(lf) == 0) then ! no wall near the cell during the step: the same state throughout
      if (np0 == np1 .and. np0 <= 1_I4P) then
         nt = np0 ; otr = 1_I4P ; ntn = 1_I4P
      else
         kind = -1_I4P
      endif
      return
   endif
   call map3_events(mc, c0, c1, lf, b, bw, nb, cert, unc)
   if (.not. cert) then
      kind = -1_I4P ; return
   endif
   if (unc) then ! an uncertified contact: a state the evaluations cannot decide (issue #102)
      if (.not. with_meas) mc%contacts = mc%contacts + 1_I4P
      kind = -1_I4P ; return
   endif
   if (nb == 0_I4P .and. np0 == np1 .and. np0 <= 1_I4P) then ! one state throughout, at most one piece
      nt = np0 ; otr = 1_I4P ; ntn = 1_I4P ; return
   endif
   kind = 1_I4P
   ! the intervals: a wide one has its evaluation (2 .. ne - 1); a narrow one refers to those before and after it
   allocate(eref(nb + 1), ealt(nb + 1), wide(nb + 1), small(nb + 1))
   ne = 1_I4P
   do m = 1_I4P, nb + 1_I4P
      call map_interval(b, nb, m, sa, sb)
      wide(m) = interval_wide(b, bw, nb, m) ; small(m) = sb - sa < EVENT_GAP
      if (wide(m)) then
         ne = ne + 1_I4P ; eref(m) = ne ; ealt(m) = ne
      else
         eref(m) = ne ; ealt(m) = 0_I4P
      endif
   enddo
   if (.not. with_meas .and. any(wide .and. small)) mc%resolved = mc%resolved + 1_I4P
   ne = ne + 1_I4P ; nev = ne
   e = ne
   do m = nb + 1_I4P, 1_I4P, -1_I4P
      if (ealt(m) == 0_I4P) then
         ealt(m) = e
      else
         e = eref(m)
      endif
   enddo
   allocate(ev(ne))
   call node_states(k0%node_in, i, j, k, mc%cs0) ; call node_states(k1%node_in, i, j, k, cs)
   call map3_corner_known(mc, cs)
   call node_states(k0%node_in, i, j, k, cs)
   fr = [k0%face_x(i, j, k), k0%face_x(i + 1, j, k), k0%face_y(i, j, k), k0%face_y(i, j + 1, k), k0%face_z(i, j, k), &
         k0%face_z(i, j, k + 1)]
   call map3_eval(mc, surface, s0, 0_I4P, 0._R8P, c0, c1, lf, cs, k0%alpha(i, j, k), fr, np0 > 1_I4P, ev(1))
   e = 1_I4P
   do m = 1_I4P, nb + 1_I4P
      call map_interval(b, nb, m, sa, sb)
      if (.not. wide(m)) cycle
      e = e + 1_I4P
      call map3_eval(mc, surface, s0, 2_I4P, 0.5_R8P * (sa + sb), c0, c1, lf, cs, 0._R8P, fr, .false., ev(e))
   enddo
   call node_states(k1%node_in, i, j, k, cs)
   fr = [k1%face_x(i, j, k), k1%face_x(i + 1, j, k), k1%face_y(i, j, k), k1%face_y(i, j + 1, k), k1%face_z(i, j, k), &
         k1%face_z(i, j, k + 1)]
   call map3_eval(mc, surface, s0, 1_I4P, 1._R8P, c0, c1, lf, cs, k1%alpha(i, j, k), fr, np1 > 1_I4P, ev(ne))
   ! every evaluation trusted; the end evaluations are the end-pose cuts' pieces, bit for bit
   do e = 1_I4P, ne
      if (.not. ev(e)%closed) kind = -1_I4P
   enddo
   if (kind == -1_I4P) return
   if (ev(1)%np /= np0 .or. ev(ne)%np /= np1) kind = -1_I4P
   if (kind == -1_I4P) return
   if (np0 > 1_I4P) then
      do q = 1_I4P, np0
         if (ev(1)%pvol(q) /= k0%piece_volume(k0%first_piece(i, j, k) + q - 1_I4P)) kind = -1_I4P
      enddo
   endif
   if (np1 > 1_I4P) then
      do q = 1_I4P, np1
         if (ev(ne)%pvol(q) /= k1%piece_volume(k1%first_piece(i, j, k) + q - 1_I4P)) kind = -1_I4P
      enddo
   endif
   if (kind == -1_I4P) return
   ! tracks: pieces of consecutive evaluations sharing a feature
   allocate(nstart(ne + 1), lstart(ne + 1))
   nstart(1) = 1_I4P ; lstart(1) = 1_I4P
   do e = 1_I4P, ne
      nstart(e + 1) = nstart(e) + ev(e)%np ; lstart(e + 1) = lstart(e) + ev(e)%nk
   enddo
   nk = lstart(ne + 1) - 1_I4P
   allocate(lkey(nk), lpiece(nk), par(nstart(ne + 1) - 1_I4P), comp(nstart(ne + 1) - 1_I4P))
   do e = 1_I4P, ne
      lkey(lstart(e):lstart(e + 1) - 1_I4P) = ev(e)%key ; lpiece(lstart(e):lstart(e + 1) - 1_I4P) = ev(e)%kpiece
   enddo
   call track_components(ne, nstart, lstart, nk, lkey, lpiece, par, comp, nt)
   otr = comp(nstart(1):nstart(2) - 1_I4P) ; ntn = comp(nstart(ne):nstart(ne + 1) - 1_I4P)
   if (.not. (with_meas .and. nt >= 1_I4P)) return
   call map3_measures(mc, surface, s0, c0, c1, lf, b, nb, eref, ealt, small, ev, nstart, comp, nt, meas, cert)
   if (.not. cert) kind = -1_I4P
   endsubroutine map3_cell

   pure subroutine node_states(node_in, i, j, k, cs)
   !< The corner states of cell (i, j, k) from the exact node flags of a cut (0 solid, 1 fluid; corner 1 + bx + 2 by + 4 bz).
   logical,      intent(in)  :: node_in(:,:,:) !< Grid node inside the solid.
   integer(I4P), intent(in)  :: i, j, k        !< Cell.
   integer(I4P), intent(out) :: cs(8)          !< Corner states.
   integer(I4P)              :: q              !< Corner.

   do q = 1_I4P, 8_I4P
      cs(q) = merge(0_I4P, 1_I4P, node_in(i + mod(q - 1_I4P, 2_I4P), j + mod((q - 1_I4P) / 2_I4P, 2_I4P), k + (q - 1_I4P) / 4_I4P))
   enddo
   endsubroutine node_states

   subroutine map3_events(mc, c0, c1, lf, b, bw, nb, cert, unc)
   !< The events of the cell [c0, c1] in (0, 1), sorted, each once, from its facets lf: a mesh vertex crossing one of its
   !< six face planes inside the face, a mesh edge crossing one of its 12 edges, one of its 8 corners crossing a facet;
   !< each mesh vertex and edge once per cell. Without rotation, explicit (the formulas of `plane_breakpoints` on the
   !< pose0 copy); with a rotation, `isolate_zeros` on the screw forms of `plane_events` (cert false if one is not
   !< certified; unc true if a kept event is not certified by a sign change: an uncertified contact). The position
   !< checks are inclusive by mc%tol: an extra event only adds an evaluation. bw (issue #104): each event's half-width
   !< (`EVENT_RES` on the explicit formulas, `isolate_zeros`' width on the screw ones; a duplicate keeps the larger).
   type(map3_ctx_t),       intent(inout) :: mc           !< State.
   real(R8P),              intent(in)    :: c0(3), c1(3) !< Cell.
   integer(I4P),           intent(in)    :: lf(:)        !< Facets.
   real(R8P), allocatable, intent(out)   :: b(:)         !< Events.
   real(R8P), allocatable, intent(out)   :: bw(:)        !< Their half-widths.
   integer(I4P),           intent(out)   :: nb           !< Their number.
   logical,                intent(out)   :: cert         !< Every event certified.
   logical,                intent(out)   :: unc          !< An event is an uncertified contact.
   real(R8P), allocatable :: tmp(:), tw(:)               !< Work.
   real(R8P)              :: p(3), lo(3), hi(3)          !< Corner, facet box.
   integer(I4P)           :: q, f, m, cn, nb0, k, nc     !< Work.

   cert = .true. ; unc = .false. ; nb = 0_I4P ; mc%tick = mc%tick + 1_I4P
   allocate(tmp(16), tw(16))
   if (allocated(mc%ctime)) deallocate(mc%ctime, mc%ccor)
   allocate(mc%ctime(8), mc%ccor(8)) ; mc%ncn = 0_I4P
   do q = 1_I4P, size(lf, kind=I4P)
      f = lf(q)
      lo = min(mc%vlo(:, 1, f), mc%vlo(:, 2, f), mc%vlo(:, 3, f))
      hi = max(mc%vhi(:, 1, f), mc%vhi(:, 2, f), mc%vhi(:, 3, f))
      do cn = 1_I4P, 8_I4P
         p = corner_point(c0, c1, cn)
         if (any(p < lo) .or. any(p > hi)) cycle
         nb0 = nb
         call event_corner(mc, f, p, tmp, tw, nb, cert, unc)
         do k = nb0 + 1_I4P, nb ! the corner's own list, for its states within the step (`map3_corner_known`)
            nc = mc%ncn
            call push_value(mc%ctime, nc, tmp(k))
            call push_value(mc%ccor, mc%ncn, real(cn, R8P))
         enddo
      enddo
      do m = 1_I4P, 3_I4P
         if (mc%vstamp(mc%fv(m, f)) /= mc%tick) then
            mc%vstamp(mc%fv(m, f)) = mc%tick
            call event_vertex(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
         endif
         if (mc%estamp(mc%fe(m, f)) /= mc%tick) then
            mc%estamp(mc%fe(m, f)) = mc%tick
            call event_edge(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
         endif
      enddo
   enddo
   call sort_unique_w(tmp, tw, nb)
   b = tmp(1:nb) ; bw = tw(1:nb)
   endsubroutine map3_events

   pure function corner_point(c0, c1, cn) result(p)
   !< Corner cn (1 + bx + 2 by + 4 bz) of the box [c0, c1].
   real(R8P),    intent(in) :: c0(3), c1(3) !< Box.
   integer(I4P), intent(in) :: cn           !< Corner.
   real(R8P)                :: p(3)         !< Point.

   p(1) = merge(c1(1), c0(1), mod(cn - 1_I4P, 2_I4P) == 1_I4P)
   p(2) = merge(c1(2), c0(2), mod((cn - 1_I4P) / 2_I4P, 2_I4P) == 1_I4P)
   p(3) = merge(c1(3), c0(3), (cn - 1_I4P) / 4_I4P == 1_I4P)
   endfunction corner_point

   subroutine event_vertex(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
   !< Events of vertex m of facet f crossing a face plane of the cell [c0, c1] inside the face.
   type(map3_ctx_t),       intent(in)    :: mc           !< State.
   integer(I4P),           intent(in)    :: f, m         !< Facet, vertex.
   real(R8P),              intent(in)    :: c0(3), c1(3) !< Cell.
   real(R8P), allocatable, intent(inout) :: tmp(:)       !< Events.
   real(R8P), allocatable, intent(inout) :: tw(:)        !< Their half-widths (issue #104).
   integer(I4P),           intent(inout) :: nb           !< Their number.
   logical,                intent(inout) :: cert         !< Certified.
   logical,                intent(inout) :: unc          !< An uncertified contact kept.
   real(R8P), allocatable :: zz(:), wd(:)                !< Zeros, their widths.
   logical,   allocatable :: un(:)                       !< Zeros not certified by a sign change.
   real(R8P)              :: lev, s, pp(3), a3(3), b3(3), c3(3), co(6) !< Work.
   integer(I4P)           :: a, bb, cc, side, nz, q      !< Work.
   logical                :: ce                          !< Certified.

   do a = 1_I4P, 3_I4P
      bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
      if (mc%vhi(bb, m, f) < c0(bb) .or. mc%vlo(bb, m, f) > c1(bb) .or. mc%vhi(cc, m, f) < c0(cc) .or. &
          mc%vlo(cc, m, f) > c1(cc)) cycle
      do side = 0_I4P, 1_I4P
         lev = merge(c1(a), c0(a), side == 1_I4P)
         if (lev < mc%vlo(a, m, f) .or. lev > mc%vhi(a, m, f)) cycle
         if (.not. mc%rot) then
            if (mc%d(a) == 0._R8P) cycle
            s = (lev - mc%t0(a, m, f)) / mc%d(a)
            if (s <= 0._R8P .or. s >= 1._R8P) cycle
            pp = mc%t0(:, m, f) + s * mc%d
            if (in_rect(pp, c0, c1, bb, cc, mc%tol)) &
               call push_pair(tmp, tw, nb, s, EVENT_RES * (s + (abs(lev) + abs(mc%t0(a, m, f))) / abs(mc%d(a))))
         else
            call rodrigues_k(mc%kax, mc%w(:, m, f), a3, b3, c3)
            co = [a3(a), 0._R8P, b3(a), 0._R8P, c3(a) + mc%pose0%t(a) - lev, mc%d(a)]
            call isolate_zeros(co, 0._R8P, mc%ang, norm2(mc%w(:, m, f)) + norm2(mc%pose0%t) + abs(lev) + norm2(mc%d), zz, &
                               nz, ce, un, wd)
            if (.not. ce) cert = .false.
            do q = 1_I4P, nz
               pp = screw_at(mc, mc%w(:, m, f), zz(q))
               if (in_rect(pp, c0, c1, bb, cc, mc%tol)) then
                  call push_pair(tmp, tw, nb, zz(q), wd(q))
                  if (un(q)) unc = .true.
               endif
            enddo
         endif
      enddo
   enddo
   endsubroutine event_vertex

   pure logical function in_rect(pt, c0, c1, bb, cc, tol)
   !< The point within the box [c0, c1] along the axes bb and cc, inclusive by tol.
   real(R8P),    intent(in) :: pt(3), c0(3), c1(3) !< Point, box.
   integer(I4P), intent(in) :: bb, cc              !< Axes.
   real(R8P),    intent(in) :: tol                 !< Tolerance.

   in_rect = pt(bb) >= c0(bb) - tol .and. pt(bb) <= c1(bb) + tol .and. pt(cc) >= c0(cc) - tol .and. pt(cc) <= c1(cc) + tol
   endfunction in_rect

   subroutine event_edge(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
   !< Events of edge m (vertex m to m + 1) of facet f crossing one of the 12 edges of the cell [c0, c1].
   type(map3_ctx_t),       intent(in)    :: mc           !< State.
   integer(I4P),           intent(in)    :: f, m         !< Facet, edge.
   real(R8P),              intent(in)    :: c0(3), c1(3) !< Cell.
   real(R8P), allocatable, intent(inout) :: tmp(:)       !< Events.
   real(R8P), allocatable, intent(inout) :: tw(:)        !< Their half-widths (issue #104).
   integer(I4P),           intent(inout) :: nb           !< Their number.
   logical,                intent(inout) :: cert         !< Certified.
   logical,                intent(inout) :: unc          !< An uncertified contact kept.
   real(R8P), allocatable :: zz(:), wd(:)                !< Zeros, their widths.
   logical,   allocatable :: un(:)                       !< Zeros not certified by a sign change.
   real(R8P)              :: lo(3), hi(3), pb, pc, p(3), e(3), det, rb, rc, lam, s, pa, pp(3), el(3), d3(3), wxd(3) !< Work.
   real(R8P)              :: xx                          !< Coordinate magnitude (event width).
   real(R8P)              :: a3(3), b3(3), c3(3), ae(3), be(3), ce3(3), f0(3), f1(3), co(6), xq(3), dd(3), fx, ox !< Work.
   integer(I4P)           :: m2, g, a, bb, cc, kk, nz, q, ax, oa !< Work.
   logical                :: ce                          !< Certified.

   m2 = mod(m, 3_I4P) + 1_I4P
   lo = min(mc%vlo(:, m, f), mc%vlo(:, m2, f)) ; hi = max(mc%vhi(:, m, f), mc%vhi(:, m2, f))
   do g = 1_I4P, 12_I4P
      a = (g - 1_I4P) / 4_I4P + 1_I4P ; kk = mod(g - 1_I4P, 4_I4P)
      bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
      pb = merge(c1(bb), c0(bb), mod(kk, 2_I4P) == 1_I4P) ; pc = merge(c1(cc), c0(cc), kk / 2_I4P == 1_I4P)
      if (pb < lo(bb) .or. pb > hi(bb) .or. pc < lo(cc) .or. pc > hi(cc) .or. c1(a) < lo(a) .or. c0(a) > hi(a)) cycle
      if (.not. mc%rot) then ! p + lam e + s d on the line (b, c) = (pb, pc): Cramer
         p = mc%t0(:, m, f) ; e = mc%t0(:, m2, f) - p
         det = e(bb) * mc%d(cc) - mc%d(bb) * e(cc)
         if (det == 0._R8P) cycle
         rb = pb - p(bb) ; rc = pc - p(cc)
         lam = (rb * mc%d(cc) - mc%d(bb) * rc) / det ; s = (e(bb) * rc - rb * e(cc)) / det
         if (s <= 0._R8P .or. s >= 1._R8P .or. lam < -1e-12_R8P .or. lam > 1._R8P + 1e-12_R8P) cycle
         pa = p(a) + lam * e(a) + s * mc%d(a)
         xx = max(maxval(abs(p)), maxval(abs(mc%t0(:, m2, f))), abs(pb), abs(pc))
         if (pa >= c0(a) - mc%tol .and. pa <= c1(a) + mc%tol) &
            call push_pair(tmp, tw, nb, s, EVENT_RES * (4._R8P * xx * xx + 2._R8P * xx * maxval(abs(mc%d)) * s) / abs(det))
      else ! the forms of `plane_events` (b)
         pp(a) = c0(a) ; pp(bb) = pb ; pp(cc) = pc
         el = 0._R8P ; el(a) = 1._R8P
         d3 = mc%w(:, m2, f) - mc%w(:, m, f) ; wxd = cross(mc%w(:, m, f), d3)
         call rodrigues_k(mc%kax, d3, a3, b3, c3) ; call rodrigues_k(mc%kax, el, ae, be, ce3)
         f0 = cross(el, pp - mc%pose0%t) ; f1 = cross(el, -mc%d)
         co = [dot_product(a3, f0) - dot_product(ae, wxd), dot_product(a3, f1), dot_product(b3, f0) + dot_product(be, wxd), &
               dot_product(b3, f1), dot_product(c3, f0) - dot_product(ce3, wxd), dot_product(c3, f1)]
         call isolate_zeros(co, 0._R8P, mc%ang, norm2(d3) * (norm2(pp) + norm2(mc%pose0%t) + norm2(mc%d) + &
                            norm2(mc%w(:, m, f))), zz, nz, ce, un, wd)
         if (.not. ce) cert = .false.
         do q = 1_I4P, nz
            xq = screw_at(mc, mc%w(:, m, f), zz(q)) ; dd = screw_at(mc, mc%w(:, m2, f), zz(q)) - xq
            if (abs(dd(bb)) >= abs(dd(cc))) then
               ax = bb ; fx = pb ; oa = cc ; ox = pc
            else
               ax = cc ; fx = pc ; oa = bb ; ox = pb
            endif
            if (dd(ax) == 0._R8P) cycle
            lam = (fx - xq(ax)) / dd(ax)
            if (lam < -1e-12_R8P .or. lam > 1._R8P + 1e-12_R8P) cycle
            if (abs(xq(oa) + lam * dd(oa) - ox) > 1e3_R8P * mc%tol) cycle
            pa = xq(a) + lam * dd(a)
            if (pa >= c0(a) - mc%tol .and. pa <= c1(a) + mc%tol) then
               call push_pair(tmp, tw, nb, zz(q), wd(q))
               if (un(q)) unc = .true.
            endif
         enddo
      endif
   enddo
   endsubroutine event_edge

   subroutine event_corner(mc, f, p, tmp, tw, nb, cert, unc)
   !< Events of the corner p crossing facet f.
   type(map3_ctx_t),       intent(in)    :: mc     !< State.
   integer(I4P),           intent(in)    :: f      !< Facet.
   real(R8P),              intent(in)    :: p(3)   !< Corner.
   real(R8P), allocatable, intent(inout) :: tmp(:) !< Events.
   real(R8P), allocatable, intent(inout) :: tw(:)  !< Their half-widths (issue #104).
   integer(I4P),           intent(inout) :: nb     !< Their number.
   logical,                intent(inout) :: cert   !< Certified.
   logical,                intent(inout) :: unc    !< An uncertified contact kept.
   real(R8P), allocatable :: zz(:), wd(:)          !< Zeros, their widths.
   logical,   allocatable :: un(:)                 !< Zeros not certified by a sign change.
   real(R8P)              :: e1(3), e2(3), r(3), det, s, be, ga, nn(3), a3(3), b3(3), c3(3), f0(3), f1(3), co(6) !< Work.
   real(R8P)              :: xx                    !< Coordinate magnitude (event width).
   real(R8P)              :: xa(3), xb(3), xc(3), el(3), bary(3) !< Work.
   integer(I4P)           :: nz, q                 !< Work.
   logical                :: ce                    !< Certified.

   if (.not. mc%rot) then ! w1 + be e1 + ga e2 + s d = p (as `plane_breakpoints`)
      e1 = mc%t0(:, 2, f) - mc%t0(:, 1, f) ; e2 = mc%t0(:, 3, f) - mc%t0(:, 1, f)
      det = dot_product(e1, cross(e2, mc%d))
      if (det == 0._R8P) return
      r = p - mc%t0(:, 1, f)
      s = dot_product(e1, cross(e2, r)) / det
      if (s <= 0._R8P .or. s >= 1._R8P) return
      be = dot_product(r, cross(e2, mc%d)) / det ; ga = dot_product(e1, cross(r, mc%d)) / det
      xx = max(maxval(abs(mc%t0(:, :, f))), maxval(abs(p)))
      if (be >= -1e-12_R8P .and. ga >= -1e-12_R8P .and. be + ga <= 1._R8P + 1e-12_R8P) &
         call push_pair(tmp, tw, nb, s, EVENT_RES * (8._R8P * xx**3 + 4._R8P * xx**2 * maxval(abs(mc%d)) * s) / abs(det))
   else ! the forms of `plane_events` (c)
      nn = cross(mc%w(:, 2, f) - mc%w(:, 1, f), mc%w(:, 3, f) - mc%w(:, 1, f))
      call rodrigues_k(mc%kax, nn, a3, b3, c3)
      f0 = p - mc%pose0%t ; f1 = -mc%d
      co = [dot_product(a3, f0), dot_product(a3, f1), dot_product(b3, f0), dot_product(b3, f1), &
            dot_product(c3, f0) - dot_product(nn, mc%w(:, 1, f)), dot_product(c3, f1)]
      call isolate_zeros(co, 0._R8P, mc%ang, norm2(nn) * (norm2(p) + norm2(mc%pose0%t) + norm2(mc%d) + &
                         norm2(mc%w(:, 1, f))), zz, nz, ce, un, wd)
      if (.not. ce) cert = .false.
      do q = 1_I4P, nz
         xa = screw_at(mc, mc%w(:, 1, f), zz(q)) ; xb = screw_at(mc, mc%w(:, 2, f), zz(q))
         xc = screw_at(mc, mc%w(:, 3, f), zz(q))
         el = cross(xb - xa, xc - xa)
         if (dot_product(el, el) <= 0._R8P) cycle
         bary = [dot_product(cross(xc - xb, p - xb), el), dot_product(cross(xa - xc, p - xc), el), &
                 dot_product(cross(xb - xa, p - xa), el)] / dot_product(el, el)
         if (all(bary >= -1e-12_R8P)) then
            call push_pair(tmp, tw, nb, zz(q), wd(q))
            if (un(q)) unc = .true.
         endif
      enddo
   endif
   endsubroutine event_corner

   function map3_corner_in(mc, surface, s0, q) result(inside)
   !< The point q inside the solid, by the exact parity of the crossings below it on its z line (`line_crossings` and
   !< `line_crossing_before`, the perturbation of the cut cells' node flags), q given in the frame of the surface the
   !< path moves: the pose0 copy (q shifted back by s d) without rotation, the body frame (`to_body`) with one.
   type(map3_ctx_t),         intent(in) :: mc          !< State.
   type(surface_stl_object), intent(in) :: surface, s0 !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in) :: q(3)        !< Point.
   logical                              :: inside      !< Inside the solid.
   real(R8P),    allocatable :: t(:)                   !< Crossings.
   integer(I4P), allocatable :: fc(:)                  !< Facets.
   logical,      allocatable :: en(:)                  !< Entering flags.
   integer(I4P)              :: n, c, cnt              !< Work.

   if (mc%rot) then
      call surface%line_crossings(3_I4P, q(1), q(2), -huge(1._R8P), huge(1._R8P), t, fc, en, n)
   else
      call s0%line_crossings(3_I4P, q(1), q(2), -huge(1._R8P), huge(1._R8P), t, fc, en, n)
   endif
   cnt = 0_I4P
   do c = 1_I4P, n
      if (mc%rot) then
         if (line_crossing_before(reshape(surface%tri_vertices(:, fc(c)), [3, 3]), 3_I4P, q(1), q(2), q(3), en(c))) &
            cnt = cnt + 1_I4P
      else
         if (line_crossing_before(reshape(s0%tri_vertices(:, fc(c)), [3, 3]), 3_I4P, q(1), q(2), q(3), en(c))) &
            cnt = cnt + 1_I4P
      endif
   enddo
   inside = mod(cnt, 2_I4P) == 1_I4P
   endfunction map3_corner_in

   subroutine map3_corner_known(mc, cs1)
   !< Which corners of the current cell have states that follow from their events (the flip rule): the corner's events
   !< (it crossing a facet, `map3_events`) at least 1e-9 apart, so that each is one transversal crossing, and their
   !< parity taking the pose0 state (mc%cs0) to the pose1 state cs1. A corner passing through a mesh edge or vertex (two
   !< facets at once), touching a facet, or whose count disagrees takes the exact query (`map3_corner_in`) instead.
   type(map3_ctx_t), intent(inout) :: mc      !< State.
   integer(I4P),     intent(in)    :: cs1(8)  !< Corner states at pose1.
   real(R8P)                       :: t(64)   !< Instants of one corner.
   integer(I4P)                    :: q, k, n !< Work.

   do q = 1_I4P, 8_I4P
      n = 0_I4P ; mc%known(q) = .true.
      do k = 1_I4P, mc%ncn
         if (nint(mc%ccor(k), I4P) /= q) cycle
         if (n == size(t, kind=I4P)) then
            mc%known(q) = .false. ; exit
         endif
         n = n + 1_I4P ; t(n) = mc%ctime(k)
      enddo
      if (.not. mc%known(q)) cycle
      call sort_unique(t, n)
      do k = 2_I4P, n
         if (t(k) - t(k - 1_I4P) <= 1e-9_R8P) mc%known(q) = .false.
      enddo
      if (count(nint(mc%ccor(1:mc%ncn), I4P) == q) /= n) mc%known(q) = .false. ! exact repeats merged by sort_unique
      if (mod(mc%cs0(q) + n, 2_I4P) /= cs1(q)) mc%known(q) = .false.
   enddo
   endsubroutine map3_corner_known

   function map3_corner_state(mc, surface, s0, s, c0, c1, b0, b1, q, known) result(st)
   !< The state (0 solid, 1 fluid) of corner q of the cell [c0, c1] at instant s: by the flip rule from its events when
   !< known (`map3_corner_known`), else by the exact query `map3_corner_in` (the cell [b0, b1] shifted back on the
   !< pose0 copy without rotation, the corner taken to the body frame with one).
   type(map3_ctx_t),         intent(in) :: mc           !< State.
   type(surface_stl_object), intent(in) :: surface, s0  !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in) :: s            !< Instant.
   real(R8P),                intent(in) :: c0(3), c1(3) !< Cell.
   real(R8P),                intent(in) :: b0(3), b1(3) !< Cell shifted back (no rotation).
   integer(I4P),             intent(in) :: q            !< Corner.
   logical,                  intent(in) :: known        !< Its states follow from its events.
   integer(I4P)                         :: st           !< State.
   integer(I4P)                         :: k, n         !< Work.

   if (known) then
      n = 0_I4P
      do k = 1_I4P, mc%ncn
         if (nint(mc%ccor(k), I4P) == q .and. mc%ctime(k) < s) n = n + 1_I4P
      enddo
      st = mod(mc%cs0(q) + n, 2_I4P)
   elseif (mc%rot) then
      st = merge(0_I4P, 1_I4P, map3_corner_in(mc, surface, s0, to_body(pose_interpolate(mc%pose0, mc%pose1, s), &
                                                                         corner_point(c0, c1, q))))
   else
      st = merge(0_I4P, 1_I4P, map3_corner_in(mc, surface, s0, corner_point(b0, b1, q)))
   endif
   endfunction map3_corner_state

   subroutine map3_eval(mc, surface, s0, kind, s, c0, c1, lf, cs_in, alpha_in, faces_in, measures, e)
   !< The pieces of the cell [c0, c1] at one instant with their features (`cut_cell_3d_pieces` on its facets lf): kind 0
   !< the pose0 copy, 1 the pose1 copy (corner states cs_in, fractions alpha_in, faces_in of the end cut), 2 the instant
   !< s: without rotation the pose0 copy and the cell shifted back by s d, else the facets placed by `to_world` of the
   !< screw pose (as `slice_surface(pose=)`); corner states by `map3_corner_in`, the face fractions (needed by a
   !< multi-piece split, and with `measures`) from the six face sections (`map3_sections`), alpha from them
   !< (`cut_cell_3d_totals`). Without `measures` the kernel runs first with fractions of -1, which fail only a
   !< multi-piece split (`split_closes`): the sections are computed only then. Features: the fluid corners (keys
   !< -8..-1) then the facets (global ids), increasing. With `measures`: the face areas and first moments of every piece
   !< (one piece: the sections'), and the volumes of several pieces.
   type(map3_ctx_t),         intent(in)  :: mc          !< State.
   type(surface_stl_object), intent(in)  :: surface, s0 !< Surface (body frame), copy at pose0.
   integer(I4P),             intent(in)  :: kind        !< 0 pose0, 1 pose1, 2 instant s.
   real(R8P),                intent(in)  :: s           !< Instant (kind 2).
   real(R8P),                intent(in)  :: c0(3), c1(3) !< Cell.
   integer(I4P),             intent(in)  :: lf(:)       !< Facets (increasing).
   integer(I4P),             intent(in)  :: cs_in(8)    !< Corner states (kinds 0, 1).
   real(R8P),                intent(in)  :: alpha_in, faces_in(6) !< Cell and face fractions (kinds 0, 1).
   logical,                  intent(in)  :: measures    !< Piece measures.
   type(cell_eval3_t),       intent(out) :: e           !< Pieces and features.
   real(R8P),    allocatable :: tl(:,:,:), pvol(:), pfc(:,:), pfm(:,:,:) !< Facets, kernel measures.
   integer(I4P), allocatable :: apiece(:)               !< Piece of each facet.
   type(rigid_pose_t)        :: ps                      !< Pose at s.
   real(R8P)                 :: b0(3), b1(3), fa(6), fm(2, 6), fr(6), al, h(3), ar(6), xc(3) !< Work.
   integer(I4P)              :: n, q, m, np, cs(8), cp(8), f, a, bb, cc, p !< Work.
   logical                   :: secs, ok, cl            !< Sections computed, sections closed, faces closed.

   n = size(lf, kind=I4P)
   allocate(tl(3, 3, n))
   b0 = c0 ; b1 = c1 ; cs = cs_in
   select case (kind)
   case (0_I4P)
      tl = mc%t0(:, :, lf)
   case (1_I4P)
      tl = mc%t1(:, :, lf)
   case default
      if (.not. mc%rot) then
         tl = mc%t0(:, :, lf) ; b0 = c0 - s * mc%d ; b1 = c1 - s * mc%d
      else
         ps = pose_interpolate(mc%pose0, mc%pose1, s)
         do q = 1_I4P, n
            do m = 1_I4P, 3_I4P
               tl(:, m, q) = to_world(ps, mc%tb(:, m, lf(q)))
            enddo
         enddo
      endif
      do q = 1_I4P, 8_I4P
         cs(q) = map3_corner_state(mc, surface, s0, s, c0, c1, b0, b1, q, mc%known(q))
      enddo
   endselect
   h = b1 - b0 ; ar = [h(2) * h(3), h(2) * h(3), h(3) * h(1), h(3) * h(1), h(1) * h(2), h(1) * h(2)]
   secs = .false. ; fa = 0._R8P ; fm = 0._R8P
   if (kind /= 2_I4P) then
      fr = faces_in ; al = alpha_in
   elseif (measures) then
      call map3_sections(tl, b0, b1, cs, ok, fa, fm)
      if (.not. ok) return
      secs = .true. ; fr = fa / ar
      call cut_cell_3d_totals(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), fr(2), fr(4), fr(6), al, cl)
   else ! fractions of -1 make a multi-piece split untrusted (`split_closes`): the sections only if one comes
      fr = -1._R8P ; al = -1._R8P
   endif
   call cut_cell_3d_pieces(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, measures, np, apiece, cp, pvol, pfc, pfm)
   if (np < 0_I4P .and. kind == 2_I4P .and. .not. secs) then
      call map3_sections(tl, b0, b1, cs, ok, fa, fm)
      if (ok) then
         secs = .true. ; fr = fa / ar
         call cut_cell_3d_totals(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), fr(2), fr(4), fr(6), al, cl)
         call cut_cell_3d_pieces(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, measures, np, apiece, cp, pvol, &
                                 pfc, pfm)
      endif
   endif
   if (np < 0_I4P .and. kind == 2_I4P .and. any(mc%known)) then ! the corner states by the exact queries, and again
      do q = 1_I4P, 8_I4P
         cs(q) = map3_corner_state(mc, surface, s0, s, c0, c1, b0, b1, q, .false.)
      enddo
      call map3_sections(tl, b0, b1, cs, ok, fa, fm)
      if (.not. ok) return
      secs = .true. ; fr = fa / ar
      call cut_cell_3d_totals(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), fr(2), fr(4), fr(6), al, cl)
      call cut_cell_3d_pieces(tl, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, measures, np, apiece, cp, pvol, pfc, &
                              pfm)
   endif
   if (np < 0_I4P) return
   e%closed = .true. ; e%np = np
   e%nk = count(cp > 0_I4P, kind=I4P) + count(apiece > 0_I4P, kind=I4P)
   allocate(e%key(e%nk), e%kpiece(e%nk))
   e%nk = 0_I4P
   do q = 8_I4P, 1_I4P, -1_I4P
      if (cp(q) == 0_I4P) cycle
      e%nk = e%nk + 1_I4P ; e%key(e%nk) = -q ; e%kpiece(e%nk) = cp(q)
   enddo
   do q = 1_I4P, n
      if (apiece(q) == 0_I4P) cycle
      e%nk = e%nk + 1_I4P ; e%key(e%nk) = lf(q) ; e%kpiece(e%nk) = apiece(q)
   enddo
   if (.not. measures) return
   allocate(e%pa(6, np), e%pm(2, 6, np), e%pvol(np))
   e%pvol = 0._R8P
   if (np == 1_I4P) then
      e%pa(:, 1) = fa ; e%pm(:, :, 1) = fm
   elseif (np > 1_I4P) then
      xc = 0.5_R8P * (b0 + b1)
      do p = 1_I4P, np
         e%pvol(p) = pvol(p)
         do f = 1_I4P, 6_I4P
            a = (f + 1_I4P) / 2_I4P ; bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
            e%pa(f, p) = pfc(f, p) * ar(f)
            e%pm(1, f, p) = xc(bb) * e%pa(f, p) + h(bb) * pfm(2, f, p)
            e%pm(2, f, p) = xc(cc) * e%pa(f, p) + h(cc) * pfm(3, f, p)
         enddo
      enddo
   endif
   endsubroutine map3_eval

   subroutine map3_sections(tl, b0, b1, cs, ok, fa, fm)
   !< The fluid area and the first moments of the two in-face coordinates (cyclic) of the six faces of the cell [b0, b1]
   !< cut by the facets tl: the section of each face plane (`slice_triangles`, the cut cells' on-plane rule) measured on
   !< the face with the cell's corner states (`cut_cell_2d_measures`). ok false if a face does not close.
   real(R8P),    intent(in)  :: tl(:,:,:)    !< Facets.
   real(R8P),    intent(in)  :: b0(3), b1(3) !< Cell.
   integer(I4P), intent(in)  :: cs(8)        !< Corner states (0 solid, 1 fluid).
   logical,      intent(out) :: ok           !< Every face closed.
   real(R8P),    intent(out) :: fa(6)        !< Fluid areas (-x, +x, -y, +y, -z, +z).
   real(R8P),    intent(out) :: fm(2, 6)     !< int u, int v over them.
   real(R8P),    allocatable :: u(:), v(:)   !< Section.
   integer(I4P), allocatable :: nxt(:)       !< Next vertex.
   integer(I4P), parameter   :: BU(4) = [0_I4P, 1_I4P, 1_I4P, 0_I4P], BV(4) = [0_I4P, 0_I4P, 1_I4P, 1_I4P] !< Face corners.
   real(R8P)                 :: lev, al, cen(2), ar !< Work.
   integer(I4P)              :: f, a, bb, cc, q, bit(3) !< Work.
   logical                   :: fl(4), cl   !< Corner flags, closed.

   ok = .true. ; fa = 0._R8P ; fm = 0._R8P
   do f = 1_I4P, 6_I4P
      a = (f + 1_I4P) / 2_I4P ; bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
      lev = merge(b1(a), b0(a), mod(f, 2_I4P) == 0_I4P)
      do q = 1_I4P, 4_I4P
         bit(a) = merge(1_I4P, 0_I4P, mod(f, 2_I4P) == 0_I4P) ; bit(bb) = BU(q) ; bit(cc) = BV(q)
         fl(q) = cs(1_I4P + bit(1) + 2_I4P * bit(2) + 4_I4P * bit(3)) == 0_I4P
      enddo
      call slice_triangles(tl, a, lev, u, v, nxt)
      call cut_cell_2d_measures(u, v, nxt, b0(bb), b1(bb), b0(cc), b1(cc), fl, cl, al, cen)
      if (.not. cl) then
         ok = .false. ; return
      endif
      ar = al * ((b1(bb) - b0(bb)) * (b1(cc) - b0(cc)))
      fa(f) = ar ; fm(1, f) = ar * cen(1) ; fm(2, f) = ar * cen(2)
   enddo
   endsubroutine map3_sections

   subroutine map3_measures(mc, surface, s0, c0, c1, lf, b, nb, eref, ealt, small, ev, nstart, comp, nt, meas, ok)
   !< Per-track time integrals over the step of the fluid areas of the cell's six faces (meas(1:6)) and of the outward
   !< rigid-field outflows through them (meas(7:12)), interval by interval between the cell's events, on every path:
   !< without rotation two-point Gauss per interval (exact: the areas are quadratic in s), with one Gauss doubled from 2
   !< nodes per interval until every track's values agree to round-off (the criteria of `one_plane`; mc%exact false if
   !< not at 64). A track exists only between events, so its integrand jumps there: the step's nq Gauss nodes of the
   !< Gauss path (`exact=.false.`) would give a short-lived track no node, or one of a large weight (2-D, NACA 0012,
   !< nq = 4: a per-track residual of 0.39 of a cell where the cell's is 1e-6). On the Gauss path the sums over a cell's
   !< tracks therefore differ from the cell's Gauss values in `sw` by the quadrature error of these (the per-track
   !< values are the better ones). ok false if a piece at a node cannot be assigned.
   type(map3_ctx_t),         intent(inout) :: mc           !< State.
   type(surface_stl_object), intent(in)    :: surface, s0  !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: c0(3), c1(3) !< Cell.
   integer(I4P),             intent(in)    :: lf(:)        !< Facets.
   real(R8P),                intent(in)    :: b(:)         !< Events.
   integer(I4P),             intent(in)    :: nb           !< Their number.
   integer(I4P),             intent(in)    :: eref(:), ealt(:) !< Evaluations of each interval (see `map3_cell`).
   logical,                  intent(in)    :: small(:)     !< Intervals narrower than EVENT_GAP (nodes may be skipped).
   type(cell_eval3_t),       intent(in)    :: ev(:)        !< Evaluations.
   integer(I4P),             intent(in)    :: nstart(:)    !< First node of each evaluation.
   integer(I4P),             intent(in)    :: comp(:)      !< Track of each node.
   integer(I4P),             intent(in)    :: nt           !< Tracks.
   real(R8P), allocatable,   intent(out)   :: meas(:,:)    !< Measures [12, nt].
   logical,                  intent(out)   :: ok           !< Every piece assigned.
   real(R8P), parameter   :: G = 0.5_R8P / sqrt(3._R8P)    !< Two-point Gauss offset on a unit interval.
   real(R8P), allocatable :: a1(:,:), a2(:,:), sq(:), wq(:) !< Work.
   real(R8P)              :: sa, sb, h, hc(3), scl(6), sco(6), sna, sno, vel, tol !< Work.
   integer(I4P)           :: m, ng, q, t, f, r, a      !< Work.
   logical                :: conv                      !< Converged.

   allocate(meas(12, nt)) ; meas = 0._R8P ; ok = .true.
   if (.not. mc%rot) then
      do m = 1_I4P, nb + 1_I4P
         call map_interval(b, nb, m, sa, sb) ; h = sb - sa
         if (h <= 0._R8P) cycle
         r = eref(m) ; a = ealt(m)
         call map3_node(mc, surface, s0, 0.5_R8P * (sa + sb) - G * h, 0.5_R8P * h, c0, c1, lf, ev(r), nstart(r), ev(a), &
                        nstart(a), r /= a .or. small(m), comp, meas, ok)
         call map3_node(mc, surface, s0, 0.5_R8P * (sa + sb) + G * h, 0.5_R8P * h, c0, c1, lf, ev(r), nstart(r), ev(a), &
                        nstart(a), r /= a .or. small(m), comp, meas, ok)
      enddo
   else
      vel = norm2(mc%om) * mc%xm + norm2(mc%d)
      hc = c1 - c0
      scl = [hc(2) * hc(3), hc(2) * hc(3), hc(3) * hc(1), hc(3) * hc(1), hc(1) * hc(2), hc(1) * hc(2)]
      sco = scl * vel ; sna = mc%xm * maxval(hc) ; sno = sna * vel
      allocate(a1(12, nt), a2(12, nt))
      do m = 1_I4P, nb + 1_I4P
         call map_interval(b, nb, m, sa, sb) ; h = sb - sa
         if (h <= 0._R8P) cycle
         r = eref(m) ; a = ealt(m)
         ng = 2_I4P ; a1 = 0._R8P
         call gauss_legendre_01(ng, sq, wq)
         do q = 1_I4P, ng
            call map3_node(mc, surface, s0, sa + h * sq(q), h * wq(q), c0, c1, lf, ev(r), nstart(r), ev(a), nstart(a), &
                           r /= a .or. small(m), comp, a1, ok)
         enddo
         do
            ng = 2_I4P * ng ; a2 = 0._R8P
            call gauss_legendre_01(ng, sq, wq)
            do q = 1_I4P, ng
               call map3_node(mc, surface, s0, sa + h * sq(q), h * wq(q), c0, c1, lf, ev(r), nstart(r), ev(a), nstart(a), &
                              r /= a .or. small(m), comp, a2, ok)
            enddo
            if (.not. ok) exit
            tol = 64._R8P * epsilon(1._R8P) * max(h, FLOOR_H)
            conv = .true.
            do t = 1_I4P, nt
               do f = 1_I4P, 6_I4P
                  if (abs(a2(f, t) - a1(f, t)) > tol * scl(f) .or. abs(a2(6 + f, t) - a1(6 + f, t)) > tol * sco(f)) &
                     conv = .false.
               enddo
            enddo
            if (conv) exit
            if (ng >= 64_I4P) then ! the round-off of one evaluation (see FLOOR_H)
               do t = 1_I4P, nt
                  do f = 1_I4P, 6_I4P
                     if (max(abs(a2(f, t) - a1(f, t)) / sna, abs(a2(6 + f, t) - a1(6 + f, t)) / sno) > &
                         64._R8P * epsilon(1._R8P)) mc%exact = .false.
                  enddo
               enddo
               exit
            endif
            a1 = a2
         enddo
         meas = meas + a2
      enddo
   endif
   endsubroutine map3_measures

   subroutine map3_node(mc, surface, s0, s, w, c0, c1, lf, evr, nr, eva, na, narrow, comp, acc, ok)
   !< One quadrature node (instant s, weight w): the pieces at s with their measures, each added to a track with w times
   !< its face areas and outward rigid-field outflows; a piece takes the track of the pieces of the interval's evaluation
   !< evr (first node nr) that hold its features, a narrow interval also tries the one after it, eva (na); as `map_node`
   !< (a negligible piece without a known feature is dropped; in a narrow interval a node or piece that cannot be
   !< assigned is skipped, mc%skipped). Outflow of a face normal to a, (a, b, c) cyclic: the field
   !< om x (x - t(s)) + d, int (om_b (x_c - t_c) - om_c (x_b - t_b) + d_a) dA from the area and the first moments (as
   !< `faces_gauss_3d` / `plane_piece`); without rotation d_a A.
   type(map3_ctx_t),         intent(inout) :: mc           !< State.
   type(surface_stl_object), intent(in)    :: surface, s0  !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: s, w         !< Instant, weight.
   real(R8P),                intent(in)    :: c0(3), c1(3) !< Cell.
   integer(I4P),             intent(in)    :: lf(:)        !< Facets.
   type(cell_eval3_t),       intent(in)    :: evr, eva     !< The interval's evaluation, the one after it.
   integer(I4P),             intent(in)    :: nr, na       !< Their first nodes.
   logical,                  intent(in)    :: narrow       !< A narrow interval (evr and eva differ; issue #102).
   integer(I4P),             intent(in)    :: comp(:)      !< Track of each node.
   real(R8P),                intent(inout) :: acc(:,:)     !< Sums [12, nt].
   logical,                  intent(inout) :: ok           !< Every piece assigned.
   type(cell_eval3_t)                      :: e            !< Pieces at s.
   real(R8P)                               :: tq(3), ar, o, fz(6) !< Work.
   integer(I4P)                            :: p, t, f, a, bb, cc, cs(8) !< Work.

   if (.not. ok) return
   cs = 0_I4P ; fz = 0._R8P
   call map3_eval(mc, surface, s0, 2_I4P, s, c0, c1, lf, cs, 0._R8P, fz, .true., e)
   if (.not. e%closed) then
      if (narrow) then
         mc%skipped = mc%skipped + 1_I4P
      else
         ok = .false.
      endif
      return
   endif
   tq = mc%pose0%t + s * mc%d
   do p = 1_I4P, e%np
      t = piece_track3(e, p, evr, nr, comp)
      if (t == 0_I4P .and. narrow) t = piece_track3(e, p, eva, na, comp)
      if (t <= 0_I4P) then ! no known feature (negligible, or a missed event) or two tracks
         if (t == 0_I4P .and. maxval(e%pa(:, p)) <= 1e-9_R8P * maxval(c1 - c0)**2) cycle
         if (narrow) then
            mc%skipped = mc%skipped + 1_I4P ; cycle
         endif
         ok = .false. ; return
      endif
      do f = 1_I4P, 6_I4P
         a = (f + 1_I4P) / 2_I4P ; bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
         ar = e%pa(f, p)
         if (mc%rot) then
            o = mc%om(bb) * (e%pm(2, f, p) - tq(cc) * ar) - mc%om(cc) * (e%pm(1, f, p) - tq(bb) * ar) + mc%d(a) * ar
         else
            o = mc%d(a) * ar
         endif
         if (mod(f, 2_I4P) == 1_I4P) o = -o
         acc(f, t) = acc(f, t) + w * ar ; acc(6 + f, t) = acc(6 + f, t) + w * o
      enddo
   enddo
   endsubroutine map3_node

   pure integer(I4P) function piece_track3(e, p, evr, nr, comp) result(t)
   !< The track of piece p of e through its features found in evr (first node nr): 0 if none is found, -1 if they point
   !< to two tracks (`piece_track` in 3-D).
   type(cell_eval3_t), intent(in) :: e, evr  !< Evaluation, reference evaluation.
   integer(I4P),       intent(in) :: p, nr   !< Piece, first node of evr.
   integer(I4P),       intent(in) :: comp(:) !< Track of each node.
   integer(I4P)                   :: k, q, lo, hi, mid, tk !< Work.

   t = 0_I4P
   do k = 1_I4P, e%nk
      if (e%kpiece(k) /= p) cycle
      q = 0_I4P ; lo = 1_I4P ; hi = evr%nk ! the feature in evr (keys increasing)
      do while (lo <= hi)
         mid = (lo + hi) / 2_I4P
         if (evr%key(mid) == e%key(k)) then
            q = evr%kpiece(mid) ; exit
         elseif (evr%key(mid) < e%key(k)) then
            lo = mid + 1_I4P
         else
            hi = mid - 1_I4P
         endif
      enddo
      if (q == 0_I4P) cycle
      tk = comp(nr + q - 1_I4P)
      if (t == 0_I4P) then
         t = tk
      elseif (t /= tk) then
         t = -1_I4P ; return
      endif
   enddo
   endfunction piece_track3

   subroutine swept_cells_3d(surface, pose0, pose1, x, y, z, sw, nq, poses, exact, status, check_input, map, refine)
   !< Swept geometry of `surface` between `pose0` and `pose1` on the block with cell edges x(nx+1),
   !< y(ny+1), z(nz+1) (issue #41 step 5): volumes from `cut_cells_3d` at the two poses; face apertures
   !< exact for a pure translation and, since issue #69, for the screw motion with a rotation (`exact`, default
   !< .true.; see the module notes and the moving-bodies page), else by `nq`-point Gauss-Legendre (default 4) on
   !< copies moved to the screw poses or to `poses(nq)`.
   !< `status = STATUS_INVALID_INPUT` on bad edges, nq < 1 or a poses size other than nq.
   !< `map` (issue #102 step 2): the piece map of the step, the tracks of every cell (`piece_map_3d_t`): the space-time
   !< components of the cell's fluid over the step, linking the pieces of the cut cells at pose0 to those at pose1, with
   !< per-track volumes, face apertures and outflows, wall and residual (see `piece_map_3d_build`). It needs the screw
   !< motion: with `poses`, `status = STATUS_INVALID_INPUT` and nothing is computed. Without `map` every other output is
   !< unchanged.
   !< `error_estimate` and `refine` (issue #104): as `swept_cells_2d` (the residual is a volume, the face change an area).
   !< Measured with `refine` (32^3, 0.5-2 cells per step): 2.7-22 times the Gauss call, about the exact path's cost.
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
   type(piece_map_3d_t),     intent(out), optional :: map      !< Piece map of the step (issue #102 step 2).
   logical,                  intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).

   call swept_3d_run(surface, pose0, pose1, x, y, z, sw, .false., SWEPT_CHUNK_CELLS, nq, poses, exact, status, check_input, &
                     map, refine)
   endsubroutine swept_cells_3d

   subroutine swept_cells_3d_device(surface, pose0, pose1, x, y, z, sw, nq, poses, exact, status, chunk_cells, refine)
   !< `swept_cells_3d` on the device (OpenACC, issue #100 step 4): the same `sw`, bit for bit, on every path.
   !< - End poses: `cut_cells_3d_device(..., alpha_only=.true.)` on the moved copies `swept_cells_3d` cuts.
   !< - Exact translation: per axis, the host computes every plane's breakpoints (`plane_breakpoints`, unchanged) and
   !<   the Gauss nodes of `plane_average`; the node sections come from one plane sweep (`slice_planes_prepare` on the
   !<   node levels, `slice_surface_plane`), their cut cells on the shifted grids from `cut_cells_device_many` (a grid per
   !<   block), and a device kernel adds 0.5 h alpha per face over the plane's nodes in node order
   !<   (`translation_axis_device`).
   !< - Exact rotation: the host computes the events of every plane (`plane_events`, unchanged); the Gauss doubling of
   !<   `one_plane` runs in rounds over all the pieces of a group of planes, each round one batch of posed sections
   !<   (`slice_surface_plane` with `pose`, the index built on the swept spans of the facets) cut on the fixed grid; a
   !<   device kernel sums each (piece, node count) over its nodes in node order, the host makes the convergence
   !<   decisions and the piece sums in piece order (`rotation_axis_device`).
   !< - Gauss poses: each node's faces by `cut_cells_3d_sections(..., device=.true., pose=)` (the faces of `cut_cells_3d`
   !<   on the moved copy, sliced at the pose without copying), the node sums on the host as in `faces_gauss_3d`.
   !< - `flux`, the wall terms and the residual: the host loops of `swept_cells_3d`.
   !< The surface needs a SAH BVH (as `cut_cells_3d_device`); `status = STATUS_INVALID_INPUT` also when the device cut
   !< cells refuse it. No `check_input` (call `swept_cells_3d` or `find_self_intersections` for that).
   !< `chunk_cells` (default 2^22) bounds the section cells of one device launch set and the face sums of one plane
   !< group of the rotation held on the host: smaller uses less memory (and more launches), with the same bits.
   !< `refine` (issue #104): the end poses' full cuts (`cut_cells_3d_device` without `alpha_only`), then the host
   !< refinement of `swept_cells_3d`.
   type(surface_stl_object), intent(in)            :: surface  !< Surface (body frame).
   type(rigid_pose_t),       intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_t),       intent(in)            :: pose1    !< Pose at the end.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(swept_cells_3d_t),   intent(out)           :: sw       !< Swept geometry.
   integer(I4P),             intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_t),       intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                  intent(in),  optional :: exact    !< Exact paths (default .true.).
   integer(I4P),             intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),             intent(in),  optional :: chunk_cells !< Section cells per launch set (default 2^22).
   logical,                  intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).
   integer(I4P)                                    :: ck       !< Chunk.

   ck = SWEPT_CHUNK_CELLS ; if (present(chunk_cells)) ck = max(1_I4P, chunk_cells)
   call swept_3d_run(surface, pose0, pose1, x, y, z, sw, .true., ck, nq, poses, exact, status, refine=refine)
   endsubroutine swept_cells_3d_device

   subroutine swept_3d_run(surface, pose0, pose1, x, y, z, sw, dev, chunk, nq, poses, exact, status, check_input, map, &
                           refine)
   !< The work of `swept_cells_3d` (dev = .false.) and `swept_cells_3d_device` (dev = .true.): the same operations in
   !< the same order, the cut cells and face sums on the device with dev. The piece map (`map`, host path only) is built
   !< last, from the end-pose cuts and the face outflows (an extra output of the face routines, their arithmetic as is).
   type(surface_stl_object), intent(in)            :: surface  !< Surface (body frame).
   type(rigid_pose_t),       intent(in)            :: pose0    !< Pose at the start.
   type(rigid_pose_t),       intent(in)            :: pose1    !< Pose at the end.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(swept_cells_3d_t),   intent(out)           :: sw       !< Swept geometry.
   logical,                  intent(in)            :: dev      !< On the device.
   integer(I4P),             intent(in)            :: chunk    !< Section cells per launch set (dev).
   integer(I4P),             intent(in),  optional :: nq       !< Gauss nodes (default 4).
   type(rigid_pose_t),       intent(in),  optional :: poses(:) !< Poses at the Gauss nodes.
   logical,                  intent(in),  optional :: exact    !< Exact translation path (default .true.).
   logical,                  intent(in),  optional :: check_input !< Check the surface for self-intersections.
   integer(I4P),             intent(out), optional :: status   !< STATUS_OK or STATUS_INVALID_INPUT.
   type(piece_map_3d_t),     intent(out), optional :: map      !< Piece map of the step (issue #102 step 2).
   logical,                  intent(in),  optional :: refine   !< Gauss path: event cells per event interval (#104).
   type(surface_stl_object) :: s0                 !< Moved copy.
   type(piece_map_3d_t)     :: mp                 !< Piece map of the refinement (issue #104).
   logical,     allocatable :: meas_c(:), rx(:,:,:), ry(:,:,:), rz(:,:,:) !< Cells integrated, faces replaced.
   real(R8P)                :: vmax, ax(3), ang   !< Largest rigid-field speed; rotation of the step.
   logical                  :: rfn, flows         !< Refine; face outflows needed.
   type(intersection_pair_t), allocatable :: pairs(:) !< Crossing facet pairs (check_input).
   real(R8P)                :: d(3)               !< Displacement.
   logical                  :: translation        !< Exact path.
   logical                  :: rotation_exact     !< Exact rotation path (issue #69).
   type(cut_cells_3d_t)     :: k0, k1             !< Cut cells at the end poses.
   real(R8P), allocatable   :: fx2(:,:,:), fy2(:,:,:), fz2(:,:,:), flux(:,:,:), fl2(:,:,:) !< Work.
   real(R8P), allocatable   :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Face outflows in +axis (map).
   integer(I4P)             :: nx, ny, nz, n, i, j, k, st !< Work.
   logical                  :: want_map           !< The piece map asked (host path).

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
   rfn = .false. ; if (present(refine)) rfn = refine
   if (present(poses)) then
      if (size(poses) /= n .or. present(map) .or. rfn) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   want_map = present(map) .and. .not. dev
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
   if (dev) then ! alpha alone (issue #100 step 4), the same copies
      s0 = surface ; call s0%move(pose1) ; call cut_cells_3d_device(s0, x, y, z, k1, st, alpha_only=.not. rfn)
      if (st == STATUS_OK) then
         s0 = surface ; call s0%move(pose0) ; call cut_cells_3d_device(s0, x, y, z, k0, st, alpha_only=.not. rfn)
      endif
      if (st /= STATUS_OK) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   else
      s0 = surface ; call s0%move(pose1) ; call cut_cells_3d(s0, x, y, z, k1, st)
      s0 = surface ; call s0%move(pose0) ; call cut_cells_3d(s0, x, y, z, k0, st)
   endif
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
   flows = want_map .or. rfn
   if (flows) allocate(fox(nx + 1, ny, nz), foy(nx, ny + 1, nz), foz(nx, ny, nz + 1))
   if (rotation_exact) then
      ! exact rotation (issue #69): certified breakpoints per grid plane, Gauss doubled to round-off on every piece
      if (want_map) then
         call faces_rotation_3d(surface, pose0, pose1, x, y, z, sw%face_x, sw%face_y, sw%face_z, flux, sw%evaluations, &
                                sw%exact, dev, chunk, fox, foy, foz)
      else
         call faces_rotation_3d(surface, pose0, pose1, x, y, z, sw%face_x, sw%face_y, sw%face_z, flux, sw%evaluations, &
                                sw%exact, dev, chunk)
      endif
   elseif (translation) then
      sw%exact = .true.
      d = pose1%t - pose0%t
      if (dev) then
         call faces_translation_3d_device(s0, d, x, y, z, sw%face_x, sw%face_y, sw%face_z, sw%evaluations, chunk)
      else
         call faces_translation_3d(s0, d, x, y, z, sw%face_x, sw%face_y, sw%face_z, sw%evaluations)
      endif
      do k = 1_I4P, nz ! the field is uniform (d): only the fluid areas enter
         do j = 1_I4P, ny
            do i = 1_I4P, nx
               flux(i, j, k) = d(1) * (sw%face_x(i + 1, j, k) - sw%face_x(i, j, k)) + &
                               d(2) * (sw%face_y(i, j + 1, k) - sw%face_y(i, j, k)) + &
                               d(3) * (sw%face_z(i, j, k + 1) - sw%face_z(i, j, k))
            enddo
         enddo
      enddo
      if (want_map) then ! the uniform field's outflows in +axis
         fox = d(1) * sw%face_x ; foy = d(2) * sw%face_y ; foz = d(3) * sw%face_z
      endif
   else
      if (flows) then
         call faces_gauss_3d(surface, pose0, pose1, x, y, z, n, sw%face_x, sw%face_y, sw%face_z, flux, dev, poses, fox, foy, &
                             foz)
      else
         call faces_gauss_3d(surface, pose0, pose1, x, y, z, n, sw%face_x, sw%face_y, sw%face_z, flux, dev, poses)
      endif
      sw%evaluations = n
   endif
   if (.not. present(poses) .and. .not. translation .and. .not. rotation_exact) then
      allocate(fx2(nx + 1, ny, nz), fy2(nx, ny + 1, nz), fz2(nx, ny, nz + 1), fl2(nx, ny, nz))
      call faces_gauss_3d(surface, pose0, pose1, x, y, z, 2_I4P * n, fx2, fy2, fz2, fl2, dev)
      if (rfn) then ! event cells per event interval (issue #104), on the host for both paths
         call piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, MAP_GAUSS, n, mp, &
                                 all_meas=.true., measured=meas_c)
         allocate(rx(nx + 1, ny, nz), ry(nx, ny + 1, nz), rz(nx, ny, nz + 1))
         call refine_faces_3d(mp, meas_c, nx, ny, nz, sw%face_x, sw%face_y, sw%face_z, fox, foy, foz, rx, ry, rz)
         do k = 1_I4P, nz
            do j = 1_I4P, ny
               do i = 1_I4P, nx
                  flux(i, j, k) = fox(i + 1, j, k) - fox(i, j, k) + foy(i, j + 1, k) - foy(i, j, k) + &
                                  foz(i, j, k + 1) - foz(i, j, k)
               enddo
            enddo
         enddo
         sw%error_estimate = max(maxval(abs(fx2 - sw%face_x), mask=.not. rx), maxval(abs(fy2 - sw%face_y), mask=.not. ry), &
                                 maxval(abs(fz2 - sw%face_z), mask=.not. rz), 0._R8P)
      else
         sw%error_estimate = max(maxval(abs(fx2 - sw%face_x)), maxval(abs(fy2 - sw%face_y)), maxval(abs(fz2 - sw%face_z)))
      endif
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
   if (.not. present(poses) .and. .not. translation .and. .not. rotation_exact) then ! the Gauss path (issue #104)
      call pose_log_3d(pose0, pose1, ax, ang)
      d = pose1%t - pose0%t
      vmax = abs(ang) * (maxval(abs([x, y, z])) + norm2(pose0%t) + norm2(d)) + norm2(d)
      if (vmax > 0._R8P) sw%error_estimate = max(sw%error_estimate, maxval(abs(sw%residual)) / vmax)
   endif
   if (rfn .and. want_map) then
      call refine_map_3d(mp, meas_c, sw, fox, foy, foz)
      map = mp
   elseif (want_map) then
      if (translation) then
         call piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, MAP_TRANSLATION, n, map)
      elseif (rotation_exact) then
         call piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, MAP_ROTATION, n, map)
      else
         call piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, MAP_GAUSS, n, map)
      endif
   endif
   endsubroutine swept_3d_run

   subroutine faces_rotation_3d(surface, pose0, pose1, x, y, z, fx, fy, fz, flux, nev, exact, device, chunk, fox, foy, foz)
   !< Exact time-averaged fluid face areas and rigid-field outflows under the screw motion with a rotation (issue #69;
   !< moving-bodies page, "Exact rotation"), one grid plane at a time. The body point w = R0 q moves to
   !< X(s) = Rot(s angle) w + t0 + s dt about the fixed unit axis k (Rodrigues: Rot(phi) w = A cos phi + B sin phi + C,
   !< C = k (k . w), A = w - C, B = k x w). A plane's breakpoints are the zeros (`isolate_zeros`) of: a mesh vertex
   !< on the plane, a mesh edge meeting one of its grid lines (det[D, e_l, N - t(s)] - det[d, Rot^T e_l, w]), one of
   !< its grid nodes on a triangle ((Rot n0) . (N - t(s)) - n0 . w_a), each checked to lie on the plane's grid and
   !< prefiltered by the box the element sweeps. A plane outside the box the body sweeps is fluid throughout (closed
   !< form). On every piece, Gauss-Legendre doubled from 2 nodes until the plane's face areas and outflows change by
   !< round-off (at 64 nodes, to that of one evaluation: FLOOR_H; else exact = .false.); one evaluation slices a
   !< moved copy (`slice_surface` with `on_plane_low` and the 2-D cut cells, as `cut_cells_3d`).
   !< With `device` (issue #100 step 4) `one_plane` only registers each plane's events and tolerances (after the same
   !< host work), and `rotation_axis_device` evaluates the registered planes of an axis in batches (see there).
   type(surface_stl_object), intent(in)  :: surface        !< Surface (body frame).
   type(rigid_pose_t),       intent(in)  :: pose0, pose1   !< End poses.
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   real(R8P),                intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Time-averaged fluid areas.
   real(R8P),                intent(out) :: flux(:,:,:)    !< Time integral of the rigid-field outflow per cell.
   integer(I4P),             intent(out) :: nev            !< Largest number of evaluations per plane.
   logical,                  intent(out) :: exact          !< Every piece converged.
   logical,      optional,   intent(in)  :: device         !< The section cut cells and face sums on the device.
   integer(I4P), optional,   intent(in)  :: chunk          !< Section cells per launch set (device).
   real(R8P),    optional,   intent(out) :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Time-integrated outflows in +axis per
                                                           !< face (issue #102: the piece map's single-track cells).
   type(facet_object), pointer :: f
   real(R8P), allocatable      :: w(:,:,:), vlo(:,:,:), vhi(:,:,:), ox(:,:,:), oy(:,:,:), oz(:,:,:), bp(:)
   real(R8P), allocatable      :: area(:,:), outf(:,:)
   real(R8P)                   :: kax(3), ang, t0(3), dt(3), om(3), blo(3), bhi(3), xs(3), mrg, g(3, 3)
   integer(I4P)                :: nf, nx, ny, nz, i, j, k, q, c, ne, nb
   logical                     :: ok
   logical                     :: dev                  !< Device mode.
   type(rotation_planes_t)     :: rp                   !< Registered planes of the current axis (device mode).
   real(R8P), allocatable      :: ra(:,:,:), ro(:,:,:) !< Their face areas and outflows (device mode).
   integer(I4P), allocatable   :: rne(:)               !< Their evaluations (device mode).
   integer(I4P)                :: n0                   !< Registered planes before a plane (device mode).
   integer(I4P)                :: ck                   !< Section cells per launch set (device mode).

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
   dev = .false. ; if (present(device)) dev = device
   ck = SWEPT_CHUNK_CELLS ; if (present(chunk)) ck = chunk
   g(:, 1) = [minval(x), minval(y), minval(z)] ; g(:, 2) = [maxval(x), maxval(y), maxval(z)]
   if (dev) call rotation_planes_reset(rp, nx + 1_I4P)
   do i = 1_I4P, nx + 1_I4P ! planes x = x_i: (u, v) = (y, z)
      n0 = rp%n
      call one_plane(1_I4P, x(i), y, z, area, outf)
      if (rp%n > n0) rp%id(rp%n) = i
      fx(i, :, :) = area ; ox(i, :, :) = outf
   enddo
   if (dev) then
      call rotation_axis_device(surface, pose0, pose1, t0, dt, om, 1_I4P, y, z, nf, w, vlo, vhi, rp, ck, ra, ro, rne, &
                                exact)
      do q = 1_I4P, rp%n
         fx(rp%id(q), :, :) = ra(:, :, q) ; ox(rp%id(q), :, :) = ro(:, :, q) ; nev = max(nev, rne(q))
      enddo
      call rotation_planes_reset(rp, ny + 1_I4P)
   endif
   do j = 1_I4P, ny + 1_I4P ! planes y = y_j: (u, v) = (z, x)
      n0 = rp%n
      call one_plane(2_I4P, y(j), z, x, area, outf)
      if (rp%n > n0) rp%id(rp%n) = j
      do k = 1_I4P, nz
         do i = 1_I4P, nx
            fy(i, j, k) = area(k, i) ; oy(i, j, k) = outf(k, i)
         enddo
      enddo
   enddo
   if (dev) then
      call rotation_axis_device(surface, pose0, pose1, t0, dt, om, 2_I4P, z, x, nf, w, vlo, vhi, rp, ck, ra, ro, rne, &
                                exact)
      do q = 1_I4P, rp%n
         do k = 1_I4P, nz
            do i = 1_I4P, nx
               fy(i, rp%id(q), k) = ra(k, i, q) ; oy(i, rp%id(q), k) = ro(k, i, q)
            enddo
         enddo
         nev = max(nev, rne(q))
      enddo
      call rotation_planes_reset(rp, nz + 1_I4P)
   endif
   do k = 1_I4P, nz + 1_I4P ! planes z = z_k: (u, v) = (x, y)
      n0 = rp%n
      call one_plane(3_I4P, z(k), x, y, area, outf)
      if (rp%n > n0) rp%id(rp%n) = k
      fz(:, :, k) = area ; oz(:, :, k) = outf
   enddo
   if (dev) then
      call rotation_axis_device(surface, pose0, pose1, t0, dt, om, 3_I4P, x, y, nf, w, vlo, vhi, rp, ck, ra, ro, rne, &
                                exact)
      do q = 1_I4P, rp%n
         fz(:, :, rp%id(q)) = ra(:, :, q) ; oz(:, :, rp%id(q)) = ro(:, :, q) ; nev = max(nev, rne(q))
      enddo
   endif
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            flux(i, j, k) = ox(i + 1, j, k) - ox(i, j, k) + oy(i, j + 1, k) - oy(i, j, k) + oz(i, j, k + 1) - oz(i, j, k)
         enddo
      enddo
   enddo
   if (present(fox)) fox = ox
   if (present(foy)) foy = oy
   if (present(foz)) foz = oz
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
      real(R8P)                           :: sa, sb, sca, sco, xm, sna, sno
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
      xm = maxval(abs([x, y, z])) + norm2(t0) + norm2(dt) ! coordinate magnitude: the round-off scale of a crossing
      sco = sca * (norm2(om) * xm + norm2(dt))
      sna = xm * max(maxval(gu(2:) - gu(:nu)), maxval(gv(2:) - gv(:nv))) ! round-off of one face area
      sno = sna * (norm2(om) * xm + norm2(dt))
      if (dev) then ! issue #100 step 4: evaluated with the other planes of the axis (`rotation_axis_device`)
         call rotation_planes_add(rp, lev, bp, nb, [sca, sco, sna, sno])
         return
      endif
      do m = 0_I4P, nb
         sa = 0._R8P ; if (m > 0_I4P) sa = bp(m)
         sb = 1._R8P ; if (m < nb) sb = bp(m + 1)
         if (sb <= sa) cycle
         ng = 2_I4P
         call plane_piece(surface, pose0, pose1, t0, dt, om, na, lev, gu, gv, iu, iv, sa, sb, ng, a1, o1) ; ne = ne + ng
         do
            ng = 2_I4P * ng
            call plane_piece(surface, pose0, pose1, t0, dt, om, na, lev, gu, gv, iu, iv, sa, sb, ng, a2, o2) ; ne = ne + ng
            if (maxval(abs(a2 - a1)) <= 64._R8P * epsilon(1._R8P) * sca * max(sb - sa, FLOOR_H) .and. &
                maxval(abs(o2 - o1)) <= 64._R8P * epsilon(1._R8P) * sco * max(sb - sa, FLOOR_H)) exit
            ! at 64 nodes: the change against the round-off of one evaluation (area: sna; issue #84, see FLOOR_H)
            if (ng >= 64_I4P) then
               if (max(maxval(abs(a2 - a1)) / sna, maxval(abs(o2 - o1)) / sno) > 64._R8P * epsilon(1._R8P)) &
                  exact = .false.
               exit
            endif
            a1 = a2 ; o1 = o2
         enddo
         ar = ar + a2 ; ou = ou + o2
      enddo
      nev = max(nev, ne)
      endsubroutine one_plane

      subroutine plane_events(na, lev, gu, gv, iu, iv, b, m)
      !< Breakpoints in (0, 1) of the plane normal to na at lev with grid gu x gv.
      integer(I4P),           intent(in)  :: na, iu, iv
      real(R8P),              intent(in)  :: lev, gu(:), gv(:)
      real(R8P), allocatable, intent(out) :: b(:)
      integer(I4P),           intent(out) :: m
      real(R8P), allocatable              :: tmp(:), zz(:)
      real(R8P)                           :: co(6), a3(3), b3(3), c3(3), ae(3), be(3), ce(3), d3(3), f0(3), f1(3)
      real(R8P)                           :: nn(3), xq(3), dd(3), lam, pp(3), tol, el(3), ebox(3, 2), tb(3, 2), bary(3)
      real(R8P)                           :: xa(3), xb(3), xc(3), wxd(3), lval
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
                  ! not merge(gu(li), gv(li), ...): merge evaluates both, out of bounds when the sizes differ
                  if (ldir == 1) then
                     lval = gu(li)
                  else
                     lval = gv(li)
                  endif
                  pp = 0._R8P ; pp(na) = lev ; pp(lfix) = lval
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
                     if (abs(pp(lfix) - lval) > tol * 1e3_R8P) cycle
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

   subroutine plane_piece(surface, pose0, pose1, t0, dt, om, na, lev, gu, gv, iu, iv, sa, sb, ng, ar, ou)
   !< ng-point Gauss-Legendre over [sa, sb] of the fluid areas and outflows of one plane (`faces_rotation_3d`). Each node
   !< slices the body-frame surface at the node pose (`slice_surface(pose=)`, issue #94): the section of the copy moved
   !< there, bit for bit, without copying and moving the surface (vertex pool, tree, pseudo-normals) per node.
   type(surface_stl_object), intent(in)  :: surface      !< Surface (body frame).
   type(rigid_pose_t),       intent(in)  :: pose0, pose1 !< End poses.
   real(R8P),                intent(in)  :: t0(3), dt(3) !< Translation at s = 0, its change over the step.
   real(R8P),                intent(in)  :: om(3)        !< Rotation vector of the step (angle times unit axis).
   integer(I4P),             intent(in)  :: na, iu, iv   !< Plane normal axis, in-plane axes (cyclic).
   integer(I4P),             intent(in)  :: ng           !< Gauss nodes.
   real(R8P),                intent(in)  :: lev          !< Plane position.
   real(R8P),                intent(in)  :: gu(:), gv(:) !< In-plane cell edges.
   real(R8P),                intent(in)  :: sa, sb       !< Piece.
   real(R8P), allocatable,   intent(out) :: ar(:,:)      !< Fluid areas integrated over the piece.
   real(R8P), allocatable,   intent(out) :: ou(:,:)      !< Outflows integrated over the piece.
   type(curve_object)                    :: sec          !< Section.
   type(cut_cells_2d_t)                  :: c2           !< Its cut cells.
   real(R8P), allocatable                :: sq(:), wq(:) !< Gauss nodes, weights.
   real(R8P)                             :: sn, tq(3), a, cw(3) !< Work.
   integer(I4P)                          :: q, ju, jv, st !< Work.

   allocate(ar(size(gu) - 1, size(gv) - 1), ou(size(gu) - 1, size(gv) - 1)) ; ar = 0._R8P ; ou = 0._R8P
   call gauss_legendre_01(ng, sq, wq)
   do q = 1_I4P, ng
      sn = sa + (sb - sa) * sq(q)
      tq = t0 + sn * dt
      call slice_surface(surface=surface, normal_axis=na, level=lev, curve=sec, on_plane_low=.true., &
                         pose=pose_interpolate(pose0, pose1, sn))
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

   subroutine faces_gauss_3d(surface, pose0, pose1, x, y, z, n, fx, fy, fz, flux, device, poses, fox, foy, foz)
   !< Time-averaged fluid face areas by n-point Gauss-Legendre (`cut_cells_3d` at each node pose), and
   !< the time integral of the outflow of the rigid field of the screw motion per cell. With `device` (issue #100 step
   !< 4) each node's faces come from `cut_cells_3d_sections(surface, ..., device=.true., pose=pq)`: the faces of
   !< `cut_cells_3d` on the moved copy, the same bits, without the copy; the node sums are the same host loops.
   type(surface_stl_object), intent(in)           :: surface      !< Surface (body frame).
   type(rigid_pose_t),       intent(in)           :: pose0, pose1 !< End poses.
   real(R8P),                intent(in)           :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),             intent(in)           :: n            !< Nodes.
   real(R8P),                intent(out)          :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Fluid areas.
   real(R8P),                intent(out)          :: flux(:,:,:)  !< Outflow of the rigid field.
   logical,                  intent(in)           :: device       !< The faces by the device sections.
   type(rigid_pose_t),       intent(in), optional :: poses(:)     !< Node poses.
   real(R8P),                intent(out), optional :: fox(:,:,:), foy(:,:,:), foz(:,:,:) !< Time-integrated outflows in
                                                                  !< +axis per face (issue #102, the piece map).
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
   if (present(fox)) fox = 0._R8P
   if (present(foy)) foy = 0._R8P
   if (present(foz)) foz = 0._R8P
   d = pose1%t - pose0%t
   call pose_log_3d(pose0, pose1, axis, angle) ; om = angle * axis
   do q = 1_I4P, n
      if (present(poses)) then
         pq = poses(q)
      else
         pq = pose_interpolate(pose0, pose1, sq(q))
      endif
      tq = pose0%t + sq(q) * d
      if (device) then ! the faces of the moved copy, sliced at the pose (no copy, no move)
         call cut_cells_3d_sections(surface, x, y, z, kc, device=.true., status=st, pose=pq)
      else
         sq_ = surface ; call sq_%move(pq) ; call cut_cells_3d(sq_, x, y, z, kc, st)
      endif
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
      if (present(fox)) fox = fox + wq(q) * ox
      if (present(foy)) foy = foy + wq(q) * oy
      if (present(foz)) foz = foz + wq(q) * oz
   enddo
   endsubroutine faces_gauss_3d

   subroutine rotation_planes_reset(rp, nmax)
   !< Empty the registered planes of `faces_rotation_3d` (device mode) for an axis of at most nmax planes.
   type(rotation_planes_t), intent(inout) :: rp   !< Registered planes.
   integer(I4P),            intent(in)    :: nmax !< Planes of the axis.

   if (allocated(rp%id)) deallocate(rp%id, rp%lev, rp%bfirst, rp%bp, rp%tol)
   allocate(rp%id(nmax), rp%lev(nmax), rp%bfirst(nmax + 1_I4P), rp%bp(64), rp%tol(4, nmax))
   rp%n = 0_I4P ; rp%nbp = 0_I4P ; rp%bfirst(1) = 1_I4P
   endsubroutine rotation_planes_reset

   subroutine rotation_planes_add(rp, lev, bp, nb, tol)
   !< Register a plane of `faces_rotation_3d` (device mode): its position, events bp(1:nb) and tolerances.
   type(rotation_planes_t), intent(inout) :: rp     !< Registered planes.
   real(R8P),               intent(in)    :: lev    !< Plane position.
   real(R8P),               intent(in)    :: bp(:)  !< Events.
   integer(I4P),            intent(in)    :: nb     !< Their number.
   real(R8P),               intent(in)    :: tol(4) !< sca, sco, sna, sno.
   real(R8P), allocatable                 :: tmp(:) !< Grown storage.

   if (rp%nbp + nb > size(rp%bp, kind=I4P)) then ! grown through a new array (no self-referencing growth)
      allocate(tmp(max(2_I4P * size(rp%bp, kind=I4P), rp%nbp + nb)))
      tmp(1:rp%nbp) = rp%bp(1:rp%nbp)
      call move_alloc(tmp, rp%bp)
   endif
   rp%n = rp%n + 1_I4P
   rp%lev(rp%n) = lev ; rp%tol(:, rp%n) = tol
   rp%bp(rp%nbp + 1_I4P:rp%nbp + nb) = bp(1:nb)
   rp%nbp = rp%nbp + nb
   rp%bfirst(rp%n + 1_I4P) = rp%nbp + 1_I4P
   endsubroutine rotation_planes_add

   subroutine rotation_axis_device(surface, pose0, pose1, t0, dt, om, na, gu, gv, nf, w, vlo, vhi, rp, ck, ra, ro, rne, &
                                   exact)
   !< The planes of one axis registered by `faces_rotation_3d` in device mode (issue #100 step 4): the face areas
   !< ra(nu, nv, p) and outflows ro(nu, nv, p) of plane p, and its evaluations rne(p), as `one_plane` computes them.
   !<
   !< Posed plane sweep: every node of plane p slices the body-frame surface at the plane through one index of the axis,
   !< built on the swept span of each facet along na: [min vlo, max vhi] over its vertices, grown by 1e-9 of its
   !< magnitude (vlo/vhi bound the path of each vertex over the whole step, `faces_rotation_3d`, with a margin of 1e-12 of
   !< the coordinate magnitude, far above the round-off by which `to_world(pose_interpolate(...))` differs from that
   !< path). A facet `slice_surface(pose=)` keeps at a node has the level in its posed span, hence in that span: the
   !< index lists a superset of those facets in increasing order, and `slice_surface_plane(pose=)` builds the same
   !< curve bit for bit. A facet with a NaN coordinate or bound is listed on every plane.
   !<
   !< The pieces of a group of planes (in plane order, piece order: those of `one_plane`) are then evaluated by
   !< `rotation_pieces_device`, and each plane's areas and outflows are summed over its pieces in piece order on the
   !< host (`ar = ar + a2` of `one_plane`). A group holds at most ck face values per piece array (at least one plane).
   !< ck also bounds the section cells per launch set.
   type(surface_stl_object), intent(in)    :: surface         !< Surface (body frame).
   type(rigid_pose_t),       intent(in)    :: pose0, pose1    !< End poses.
   real(R8P),                intent(in)    :: t0(3), dt(3), om(3) !< Translation at s = 0, its change, rotation vector.
   integer(I4P),             intent(in)    :: na              !< Plane normal axis.
   real(R8P),                intent(in)    :: gu(:), gv(:)    !< In-plane cell edges (cyclic order).
   integer(I4P),             intent(in)    :: nf              !< Facets.
   real(R8P),                intent(in)    :: w(:,:,:)        !< R0 times the body vertices [3, 3, nf].
   real(R8P),                intent(in)    :: vlo(:,:,:), vhi(:,:,:) !< Boxes of the vertex paths [3, 3, nf].
   type(rotation_planes_t),  intent(in)    :: rp              !< Registered planes.
   integer(I4P),             intent(in)    :: ck              !< Section cells per launch set, face values per group.
   real(R8P), allocatable,   intent(out)   :: ra(:,:,:)       !< Face areas [nu, nv, rp%n].
   real(R8P), allocatable,   intent(out)   :: ro(:,:,:)       !< Face outflows [nu, nv, rp%n].
   integer(I4P), allocatable, intent(out)  :: rne(:)          !< Evaluations per plane [rp%n].
   logical,                  intent(inout) :: exact           !< Cleared if a piece misses the round-off at 64 nodes.
   type(slice_planes_t)      :: idx                           !< Posed plane sweep index.
   real(R8P),    allocatable :: sp(:,:)                       !< Swept spans of the facets [2, nf].
   integer(I4P), allocatable :: kp(:), pk(:)                  !< Plane of each piece, pieces of each plane [npl + 1].
   real(R8P),    allocatable :: ksa(:), ksb(:)                !< Piece ends.
   real(R8P),    allocatable :: fa(:,:,:), fo(:,:,:)          !< Piece sums of a group.
   real(R8P)                 :: lo, hi, mg, sa, sb            !< Work.
   integer(I4P)              :: nu, nv, npl, npc, p, m, nb, k, c, i, p0, p1, k0, k1 !< Work.
   logical                   :: bad                           !< A NaN in a facet.

   nu = size(gu, kind=I4P) - 1_I4P ; nv = size(gv, kind=I4P) - 1_I4P ; npl = rp%n
   allocate(ra(nu, nv, npl), ro(nu, nv, npl), rne(npl))
   ra = 0._R8P ; ro = 0._R8P ; rne = 0_I4P
   if (npl == 0_I4P) return
   allocate(sp(2, max(1_I4P, nf)))
   do i = 1_I4P, nf
      lo = vlo(na, 1, i) ; hi = vhi(na, 1, i)
      do c = 2, 3
         if (vlo(na, c, i) < lo) lo = vlo(na, c, i)
         if (vhi(na, c, i) > hi) hi = vhi(na, c, i)
      enddo
      bad = any(w(:, :, i) /= w(:, :, i)) .or. any(vlo(na, :, i) /= vlo(na, :, i)) .or. any(vhi(na, :, i) /= vhi(na, :, i))
      if (bad .or. .not. (abs(lo) <= huge(1._R8P) .and. abs(hi) <= huge(1._R8P))) then
         sp(1, i) = -huge(1._R8P) ; sp(2, i) = huge(1._R8P)
      else
         mg = 1e-9_R8P * (abs(lo) + abs(hi))
         sp(1, i) = lo - mg ; sp(2, i) = hi + mg
      endif
   enddo
   call slice_planes_prepare(surface, na, rp%lev(1:npl), idx, on_plane_low=.true., spans=sp(:, 1:nf))
   ! the pieces of `one_plane`, in plane order and piece order
   allocate(pk(npl + 1_I4P))
   npc = 0_I4P
   do p = 1_I4P, npl
      pk(p) = npc + 1_I4P
      nb = rp%bfirst(p + 1_I4P) - rp%bfirst(p)
      do m = 0_I4P, nb
         sa = 0._R8P ; if (m > 0_I4P) sa = rp%bp(rp%bfirst(p) + m - 1_I4P)
         sb = 1._R8P ; if (m < nb) sb = rp%bp(rp%bfirst(p) + m)
         if (sb <= sa) cycle
         npc = npc + 1_I4P
      enddo
   enddo
   pk(npl + 1_I4P) = npc + 1_I4P
   allocate(kp(npc), ksa(npc), ksb(npc))
   k = 0_I4P
   do p = 1_I4P, npl
      nb = rp%bfirst(p + 1_I4P) - rp%bfirst(p)
      do m = 0_I4P, nb
         sa = 0._R8P ; if (m > 0_I4P) sa = rp%bp(rp%bfirst(p) + m - 1_I4P)
         sb = 1._R8P ; if (m < nb) sb = rp%bp(rp%bfirst(p) + m)
         if (sb <= sa) cycle
         k = k + 1_I4P ; kp(k) = p ; ksa(k) = sa ; ksb(k) = sb
      enddo
   enddo
   ! groups of whole planes, then the plane sums in piece order
   p0 = 1_I4P
   do while (p0 <= npl)
      p1 = p0
      do while (p1 < npl)
         if (int(pk(p1 + 2_I4P) - pk(p0), I8P) * nu * nv > ck) exit
         p1 = p1 + 1_I4P
      enddo
      k0 = pk(p0) ; k1 = pk(p1 + 1_I4P) - 1_I4P
      if (k1 >= k0) then
         call rotation_pieces_device(surface, idx, pose0, pose1, t0, dt, om, na, gu, gv, rp, ck, k1 - k0 + 1_I4P, &
                                     kp(k0:k1), ksa(k0:k1), ksb(k0:k1), fa, fo, rne, exact)
         do k = k0, k1
            p = kp(k)
            ra(:, :, p) = ra(:, :, p) + fa(:, :, k - k0 + 1_I4P) ; ro(:, :, p) = ro(:, :, p) + fo(:, :, k - k0 + 1_I4P)
         enddo
      endif
      p0 = p1 + 1_I4P
   enddo
   endsubroutine rotation_axis_device

   subroutine rotation_pieces_device(surface, idx, pose0, pose1, t0, dt, om, na, gu, gv, rp, ck, npc, kp, ksa, ksb, fa, &
                                     fo, rne, exact)
   !< The Gauss doubling of `one_plane` on npc pieces at once (issue #100 step 4): round 1 evaluates every piece at 2
   !< and 4 nodes, each later round the live pieces at twice their last count, through one batch per round
   !< (`rotation_groups_device`). The decisions per piece are those of `one_plane` on the same values (the same
   !< `maxval` tests on whole arrays, the same tolerances), so every piece ends at the node count, with the sums, the
   !< evaluations and the exactness flag of the one-piece loop. fa, fo: the final sums of each piece (a2, o2).
   type(surface_stl_object), intent(in)    :: surface         !< Surface (body frame).
   type(slice_planes_t),     intent(in)    :: idx             !< Posed plane sweep index of the axis.
   type(rigid_pose_t),       intent(in)    :: pose0, pose1    !< End poses.
   real(R8P),                intent(in)    :: t0(3), dt(3), om(3) !< Translation at s = 0, its change, rotation vector.
   integer(I4P),             intent(in)    :: na              !< Plane normal axis.
   real(R8P),                intent(in)    :: gu(:), gv(:)    !< In-plane cell edges.
   type(rotation_planes_t),  intent(in)    :: rp              !< Registered planes.
   integer(I4P),             intent(in)    :: ck              !< Section cells per launch set.
   integer(I4P),             intent(in)    :: npc             !< Pieces.
   integer(I4P),             intent(in)    :: kp(:)           !< Plane of each piece [npc].
   real(R8P),                intent(in)    :: ksa(:), ksb(:)  !< Piece ends [npc].
   real(R8P), allocatable,   intent(out)   :: fa(:,:,:)       !< Areas of each piece [nu, nv, npc].
   real(R8P), allocatable,   intent(out)   :: fo(:,:,:)       !< Outflows of each piece [nu, nv, npc].
   integer(I4P),             intent(inout) :: rne(:)          !< Evaluations per plane.
   logical,                  intent(inout) :: exact           !< Cleared if a piece misses the round-off at 64 nodes.
   real(R8P),    allocatable :: pa(:,:,:), po(:,:,:)          !< Last level of each piece (then its final sums).
   real(R8P),    allocatable :: sa(:,:,:), so(:,:,:)          !< Sums of the round's groups.
   real(R8P),    allocatable :: a1(:,:), o1(:,:), a2(:,:), o2(:,:) !< One piece's levels (as in `one_plane`).
   integer(I4P), allocatable :: ng(:), gk(:), gng(:)          !< Node count per piece; group piece, node count.
   logical,      allocatable :: live(:)                       !< Piece still doubling.
   real(R8P)                 :: sca, sco, sna, sno            !< Tolerances of the piece's plane.
   integer(I4P)              :: nu, nv, ngr, g, k, p          !< Work.
   logical                   :: first                         !< Round 1.

   nu = size(gu, kind=I4P) - 1_I4P ; nv = size(gv, kind=I4P) - 1_I4P
   allocate(pa(nu, nv, npc), po(nu, nv, npc), ng(npc), live(npc), a1(nu, nv), o1(nu, nv), a2(nu, nv), o2(nu, nv))
   ng = 4_I4P ; live = .true. ; first = .true.
   do
      ngr = count(live, kind=I4P) ; if (first) ngr = 2_I4P * ngr
      if (ngr == 0_I4P) exit
      allocate(gk(ngr), gng(ngr))
      g = 0_I4P
      do k = 1_I4P, npc
         if (.not. live(k)) cycle
         if (first) then
            g = g + 1_I4P ; gk(g) = k ; gng(g) = 2_I4P
         endif
         g = g + 1_I4P ; gk(g) = k ; gng(g) = ng(k)
      enddo
      call rotation_groups_device(surface, idx, pose0, pose1, t0, dt, om, na, gu, gv, rp%lev, ck, ngr, kp(gk), ksa(gk), &
                                  ksb(gk), gng, sa, so)
      g = 0_I4P
      do k = 1_I4P, npc
         if (.not. live(k)) cycle
         p = kp(k)
         sca = rp%tol(1, p) ; sco = rp%tol(2, p) ; sna = rp%tol(3, p) ; sno = rp%tol(4, p)
         if (first) then ! a1: the 2-node level
            g = g + 1_I4P
            pa(:, :, k) = sa(:, :, g) ; po(:, :, k) = so(:, :, g) ; rne(p) = rne(p) + 2_I4P
         endif
         g = g + 1_I4P
         rne(p) = rne(p) + ng(k)
         a1 = pa(:, :, k) ; o1 = po(:, :, k) ; a2 = sa(:, :, g) ; o2 = so(:, :, g)
         pa(:, :, k) = a2 ; po(:, :, k) = o2
         if (maxval(abs(a2 - a1)) <= 64._R8P * epsilon(1._R8P) * sca * max(ksb(k) - ksa(k), FLOOR_H) .and. &
             maxval(abs(o2 - o1)) <= 64._R8P * epsilon(1._R8P) * sco * max(ksb(k) - ksa(k), FLOOR_H)) then
            live(k) = .false.
         elseif (ng(k) >= 64_I4P) then
            if (max(maxval(abs(a2 - a1)) / sna, maxval(abs(o2 - o1)) / sno) > 64._R8P * epsilon(1._R8P)) exact = .false.
            live(k) = .false.
         else
            ng(k) = 2_I4P * ng(k)
         endif
      enddo
      deallocate(gk, gng)
      first = .false.
   enddo
   call move_alloc(pa, fa) ; call move_alloc(po, fo)
   endsubroutine rotation_pieces_device

   subroutine rotation_groups_device(surface, idx, pose0, pose1, t0, dt, om, na, gu, gv, lev, ck, ngr, gpl, gsa, gsb, &
                                     gng, sa, so)
   !< One round of `rotation_pieces_device`: group g is the gng(g)-node Gauss rule of `plane_piece` on the piece
   !< [gsa(g), gsb(g)] of plane gpl(g); sa(:, :, g), so(:, :, g) are its sums ar, ou of `plane_piece`, bit for bit.
   !<
   !< The host lists the nodes in group order, each with what `plane_piece` computes for it (sn = sa + (sb - sa) sq,
   !< the weight (sb - sa) wq, tq = t0 + sn dt, the pose `pose_interpolate(pose0, pose1, sn)`), and slices them through
   !< the posed index (`slice_surface_plane(pose=)`, equal to `slice_surface(pose=)`) in chunks of at most ck section
   !< cells (at least one node); each chunk is cut in one `cut_cells_device_many` (the fixed in-plane grid,
   !< resident: alpha and the centroids are those of `curve%cut_cells`) and `rotation_reduce_dev` adds the chunk's nodes
   !< to their groups. The sums live on the device across the chunks, each updated by one thread per face, its nodes in
   !< node order: the additions of `plane_piece` in its order.
   type(surface_stl_object), intent(in)  :: surface        !< Surface (body frame).
   type(slice_planes_t),     intent(in)  :: idx            !< Posed plane sweep index of the axis.
   type(rigid_pose_t),       intent(in)  :: pose0, pose1   !< End poses.
   real(R8P),                intent(in)  :: t0(3), dt(3), om(3) !< Translation at s = 0, its change, rotation vector.
   integer(I4P),             intent(in)  :: na             !< Plane normal axis.
   real(R8P),                intent(in)  :: gu(:), gv(:)   !< In-plane cell edges.
   real(R8P),                intent(in)  :: lev(:)         !< Plane positions (index planes).
   integer(I4P),             intent(in)  :: ck             !< Section cells per launch set.
   integer(I4P),             intent(in)  :: ngr            !< Groups.
   integer(I4P),             intent(in)  :: gpl(:)         !< Plane of each group [ngr].
   real(R8P),                intent(in)  :: gsa(:), gsb(:) !< Piece of each group [ngr].
   integer(I4P),             intent(in)  :: gng(:)         !< Gauss nodes of each group [ngr].
   real(R8P), allocatable,   intent(out) :: sa(:,:,:)      !< Area sums [nu, nv, ngr].
   real(R8P), allocatable,   intent(out) :: so(:,:,:)      !< Outflow sums [nu, nv, ngr].
   type(curve_object),   allocatable :: secs(:)            !< Sections of a chunk.
   type(cut_cells_2d_blocks_t)       :: cb                 !< Their cut cells.
   real(R8P),    allocatable :: sq(:), wq(:)               !< Gauss rule.
   real(R8P),    allocatable :: nsn(:), nwf(:), ntq(:,:), nlv(:) !< Node instants, weights, translations, levels.
   integer(I4P), allocatable :: npl(:), ngp(:), gfirst(:)  !< Node plane, node group; nodes of each group [ngr + 1].
   real(R8P),    allocatable :: sal(:,:,:), sol(:,:,:)     !< Sums (built here).
   real(R8P),    allocatable :: gul(:), gvl(:)             !< Contiguous copies of the cell edges (mapped).
   real(R8P)                 :: om1, om2, om3, dtn         !< Rotation vector, dt(na).
   integer(I4P)              :: nu, nv, nn, g, q, b, last, c0, c1, kc, st !< Work.

   nu = size(gu, kind=I4P) - 1_I4P ; nv = size(gv, kind=I4P) - 1_I4P
   nn = sum(gng(1:ngr))
   allocate(nsn(nn), nwf(nn), ntq(3, nn), nlv(nn), npl(nn), ngp(nn), gfirst(ngr + 1_I4P))
   b = 0_I4P ; last = 0_I4P
   do g = 1_I4P, ngr
      if (gng(g) /= last) then
         call gauss_legendre_01(gng(g), sq, wq) ; last = gng(g)
      endif
      gfirst(g) = b + 1_I4P
      do q = 1_I4P, gng(g)
         b = b + 1_I4P
         nsn(b) = gsa(g) + (gsb(g) - gsa(g)) * sq(q)
         nwf(b) = (gsb(g) - gsa(g)) * wq(q)
         ntq(:, b) = t0 + nsn(b) * dt
         nlv(b) = lev(gpl(g)) ; npl(b) = gpl(g) ; ngp(b) = g
      enddo
   enddo
   gfirst(ngr + 1_I4P) = nn + 1_I4P
   allocate(sal(nu, nv, ngr), sol(nu, nv, ngr))
   sal = 0._R8P ; sol = 0._R8P
   om1 = om(1) ; om2 = om(2) ; om3 = om(3) ; dtn = dt(na)
   kc = max(1_I4P, ck / (nu * nv))
   gul = gu ; gvl = gv
   !$acc data copyin(gul, gvl, gfirst, nlv, nwf, ntq) copy(sal, sol)
   c0 = 1_I4P
   do while (c0 <= nn)
      c1 = min(nn, c0 + kc - 1_I4P)
      allocate(secs(c1 - c0 + 1_I4P))
      do b = c0, c1
         call slice_surface_plane(surface, idx, npl(b), secs(b - c0 + 1_I4P), pose=pose_interpolate(pose0, pose1, nsn(b)))
      enddo
      call cut_cells_device_many(secs, gul, gvl, cb, st, intervals=.false., pieces=.false., resident=.true.)
      deallocate(secs)
      call rotation_reduce_dev(nu, nv, c1 - c0 + 1_I4P, cb%alpha, cb%centroid, gul, gvl, na, om1, om2, om3, dtn, c0, &
                               ngp(c0), ngp(c1), ngr, gfirst, nn, nlv, nwf, ntq, sal, sol)
      call cb%exit_device
      c0 = c1 + 1_I4P
   enddo
   !$acc end data
   call move_alloc(sal, sa) ; call move_alloc(sol, so)
   endsubroutine rotation_groups_device

   subroutine rotation_reduce_dev(nu, nv, nbk, alpha, cen, gu, gv, na, om1, om2, om3, dtn, c0, g0, g1, ngr, gfirst, nn, &
                                  nlv, nwf, ntq, sa, so)
   !< Device kernel (issue #100 step 4): the chunk of nodes c0 .. c0 + nbk - 1 (blocks 1 .. nbk of the section cut
   !< cells) added to the sums of groups g0 .. g1, one thread per (group, face), the group's nodes of the chunk in node
   !< order: per node `plane_piece`'s a = alpha du dv, ar = ar + w a, ou = ou + w a (cross_k(om, cw - tq, na) + dt(na))
   !< with w = (sb - sa) wq and cw the fluid centroid on the plane (the component na of the cross product written out
   !< as `cross` computes it). Every array must be present on the device.
   integer(I4P), intent(in)    :: nu, nv, nbk, na      !< Face cells, chunk blocks, plane normal axis.
   real(R8P),    intent(in)    :: alpha(nu, nv, nbk)   !< Section fluid fractions.
   real(R8P),    intent(in)    :: cen(2, nu, nv, nbk)  !< Section fluid centroids.
   real(R8P),    intent(in)    :: gu(nu + 1), gv(nv + 1) !< In-plane cell edges.
   real(R8P),    intent(in)    :: om1, om2, om3, dtn   !< Rotation vector, dt(na).
   integer(I4P), intent(in)    :: c0, g0, g1, ngr, nn  !< First node of the chunk, its groups, groups, nodes.
   integer(I4P), intent(in)    :: gfirst(ngr + 1)      !< Nodes of each group.
   real(R8P),    intent(in)    :: nlv(nn), nwf(nn), ntq(3, nn) !< Node levels, weights, translations.
   real(R8P),    intent(inout) :: sa(nu, nv, ngr), so(nu, nv, ngr) !< Sums.
   real(R8P)                   :: ra, ro, a, wa, v1, v2, v3, cr !< Work.
   integer(I4P)                :: g, ju, jv, b, b0, b1, l !< Work.

   !$acc parallel loop gang vector collapse(3) default(present) private(ra, ro, a, wa, v1, v2, v3, cr, b0, b1, l) &
   !$acc& present(alpha, cen, gu, gv, gfirst, nlv, nwf, ntq, sa, so)
   do g = g0, g1
      do jv = 1_I4P, nv
         do ju = 1_I4P, nu
            b0 = max(gfirst(g), c0) ; b1 = min(gfirst(g + 1) - 1_I4P, c0 + nbk - 1_I4P)
            ra = sa(ju, jv, g) ; ro = so(ju, jv, g)
            !$acc loop seq
            do b = b0, b1
               l = b - c0 + 1_I4P
               a = alpha(ju, jv, l) * (gu(ju + 1) - gu(ju)) * (gv(jv + 1) - gv(jv))
               select case (na) ! v = cw - tq with cw(na) = level, cw(iu) = centroid 1, cw(iv) = centroid 2
               case (1_I4P)
                  v1 = nlv(b) - ntq(1, b) ; v2 = cen(1, ju, jv, l) - ntq(2, b) ; v3 = cen(2, ju, jv, l) - ntq(3, b)
                  cr = om2 * v3 - om3 * v2
               case (2_I4P)
                  v2 = nlv(b) - ntq(2, b) ; v3 = cen(1, ju, jv, l) - ntq(3, b) ; v1 = cen(2, ju, jv, l) - ntq(1, b)
                  cr = om3 * v1 - om1 * v3
               case default
                  v3 = nlv(b) - ntq(3, b) ; v1 = cen(1, ju, jv, l) - ntq(1, b) ; v2 = cen(2, ju, jv, l) - ntq(2, b)
                  cr = om1 * v2 - om2 * v1
               endselect
               wa = nwf(b) * a
               ra = ra + wa
               ro = ro + wa * (cr + dtn)
            enddo
            sa(ju, jv, g) = ra ; so(ju, jv, g) = ro
         enddo
      enddo
   enddo
   endsubroutine rotation_reduce_dev

   subroutine faces_translation_3d_device(s0, d, x, y, z, fx, fy, fz, nev, ck)
   !< `faces_translation_3d` with the face sums on the device (issue #100 step 4): per axis the breakpoints of every
   !< plane (`plane_breakpoints`, unchanged), then `translation_axis_device` for all the planes of the axis.
   type(surface_stl_object), intent(in)  :: s0        !< Surface at pose0.
   real(R8P),                intent(in)  :: d(3)      !< Displacement over the step.
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   real(R8P),                intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:) !< Fluid areas.
   integer(I4P),             intent(out) :: nev       !< Largest number of evaluations per plane.
   integer(I4P),             intent(in)  :: ck        !< Section cells per launch set.
   type(facet_object), pointer :: fp                  !< Facet.
   real(R8P), allocatable      :: tri(:,:,:)          !< Triangles [3, 3, facets].
   real(R8P), allocatable      :: ue(:), ve(:), bp(:), lv(:), bps(:), tmp(:), area(:,:,:) !< Work.
   integer(I4P), allocatable   :: bf(:)               !< Breakpoints of each plane [np + 1].
   real(R8P)                   :: dl(3)               !< Displacement in the plane frame.
   integer(I4P)                :: nf, f, m, ax, a1, a2, ip, np, nb, nbp !< Work.

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
      case (1_I4P) ; np = size(x, kind=I4P) ; ue = y ; ve = z ; lv = x
      case (2_I4P) ; np = size(y, kind=I4P) ; ue = z ; ve = x ; lv = y
      case default ; np = size(z, kind=I4P) ; ue = x ; ve = y ; lv = z
      endselect
      dl = [d(ax), d(a1), d(a2)]
      if (allocated(bf)) deallocate(bf, bps)
      allocate(bf(np + 1_I4P), bps(64))
      bf(1) = 1_I4P ; nbp = 0_I4P
      do ip = 1_I4P, np
         call plane_breakpoints(tri, [ax, a1, a2], lv(ip), dl, ue, ve, bp, nb)
         if (nbp + nb > size(bps, kind=I4P)) then ! grown through a new array (no self-referencing growth)
            allocate(tmp(max(2_I4P * size(bps, kind=I4P), nbp + nb)))
            tmp(1:nbp) = bps(1:nbp)
            call move_alloc(tmp, bps)
         endif
         bps(nbp + 1_I4P:nbp + nb) = bp(1:nb) ; nbp = nbp + nb ; bf(ip + 1_I4P) = nbp + 1_I4P
         nev = max(nev, 2_I4P * (nb + 1_I4P))
      enddo
      call translation_axis_device(s0, ax, lv, dl, ue, ve, bf, bps, ck, area)
      do ip = 1_I4P, np
         select case (ax)
         case (1_I4P) ; fx(ip, :, :) = area(:, :, ip)
         case (2_I4P) ; fy(:, ip, :) = transpose(area(:, :, ip))
         case default ; fz(:, :, ip) = area(:, :, ip)
         endselect
      enddo
   enddo
   endsubroutine faces_translation_3d_device

   subroutine translation_axis_device(s0, ax, lv, dl, ue, ve, bf, bp, ck, area)
   !< The faces of the planes x_ax = lv(p) of `faces_translation_3d` (issue #100 step 4): area(:, :, p) equals
   !< `plane_average` of plane p bit for bit.
   !<
   !< The host lists the nodes of `plane_average` in its order (plane by plane, piece by piece, the two Gauss points
   !< sg = (sa + sb) / 2 -+ G h of every piece with h > 0): level lv(p) - sg dl(1), weight 0.5 h, grid ue - sg dl(2),
   !< ve - sg dl(3). The node sections come from a plane sweep per chunk on its node levels (`slice_planes_prepare`,
   !< `slice_surface_plane`: the curves of `slice_surface`), cut in chunks of at most ck cells (one node at least) by
   !< `cut_cells_device_many` with a grid per block (the cut cells of `curve%cut_cells` on the shifted grid); per chunk
   !< `translation_reduce_dev` adds 0.5 h alpha of its nodes to their planes' faces, one thread per face, in node order
   !< (`area = area + 0.5 h alpha` of `plane_average`). The sums stay on the device across the chunks; the host then
   !< scales them by the face sizes as `plane_average` does.
   type(surface_stl_object), intent(in)  :: s0           !< Surface at pose0.
   integer(I4P),             intent(in)  :: ax           !< Plane normal axis.
   real(R8P),                intent(in)  :: lv(:)        !< Plane positions [np].
   real(R8P),                intent(in)  :: dl(3)        !< Displacement in the plane frame.
   real(R8P),                intent(in)  :: ue(:), ve(:) !< In-plane grid edges.
   integer(I4P),             intent(in)  :: bf(:)        !< Breakpoints of plane p: bp(bf(p) : bf(p + 1) - 1) [np + 1].
   real(R8P),                intent(in)  :: bp(:)        !< Breakpoints.
   integer(I4P),             intent(in)  :: ck           !< Section cells per launch set.
   real(R8P), allocatable,   intent(out) :: area(:,:,:)  !< Time-averaged fluid areas [nu, nv, np].
   real(R8P), parameter      :: G = 0.5_R8P / sqrt(3._R8P) !< Two-point Gauss offset on a unit piece (`plane_average`).
   type(slice_planes_t)      :: idx                      !< Plane sweep index of the node levels.
   type(curve_object),   allocatable :: secs(:)          !< Sections of a chunk.
   type(cut_cells_2d_blocks_t)       :: cb               !< Their cut cells.
   real(R8P),    allocatable :: nlev(:), nsg(:), nhw(:)  !< Node levels, instants, weights.
   integer(I4P), allocatable :: npl(:), pfirst(:)        !< Node plane; nodes of each plane [np + 1].
   real(R8P),    allocatable :: acc(:,:,:), xg(:,:), yg(:,:) !< Sums; grids of a chunk.
   real(R8P)                 :: sa, sb, h, sg            !< Work.
   integer(I4P)              :: np, nu, nv, nn, p, m, q, g2, b, c0, c1, kc, st, j, k !< Work.

   np = size(lv, kind=I4P) ; nu = size(ue, kind=I4P) - 1_I4P ; nv = size(ve, kind=I4P) - 1_I4P
   nn = 0_I4P
   do p = 1_I4P, np
      m = bf(p + 1_I4P) - bf(p)
      do q = 0_I4P, m
         sa = 0._R8P ; if (q > 0_I4P) sa = bp(bf(p) + q - 1_I4P)
         sb = 1._R8P ; if (q < m) sb = bp(bf(p) + q)
         if (sb - sa <= 0._R8P) cycle
         nn = nn + 2_I4P
      enddo
   enddo
   allocate(nlev(nn), nsg(nn), nhw(nn), npl(nn), pfirst(np + 1_I4P))
   b = 0_I4P
   do p = 1_I4P, np
      pfirst(p) = b + 1_I4P
      m = bf(p + 1_I4P) - bf(p)
      do q = 0_I4P, m
         sa = 0._R8P ; if (q > 0_I4P) sa = bp(bf(p) + q - 1_I4P)
         sb = 1._R8P ; if (q < m) sb = bp(bf(p) + q)
         h = sb - sa
         if (h <= 0._R8P) cycle
         do g2 = -1_I4P, 1_I4P, 2_I4P
            sg = 0.5_R8P * (sa + sb) + g2 * G * h
            b = b + 1_I4P
            nsg(b) = sg ; nlev(b) = lv(p) - sg * dl(1) ; nhw(b) = 0.5_R8P * h ; npl(b) = p
         enddo
      enddo
   enddo
   pfirst(np + 1_I4P) = nn + 1_I4P
   allocate(acc(nu, nv, np)) ; acc = 0._R8P
   kc = max(1_I4P, ck / (nu * nv))
   !$acc data copyin(nhw, pfirst) copy(acc)
   c0 = 1_I4P
   do while (c0 <= nn)
      c1 = min(nn, c0 + kc - 1_I4P)
      allocate(secs(c1 - c0 + 1_I4P), xg(nu + 1_I4P, c1 - c0 + 1_I4P), yg(nv + 1_I4P, c1 - c0 + 1_I4P))
      call slice_planes_prepare(s0, ax, nlev(c0:c1), idx, on_plane_low=.true.) ! the chunk's levels (index memory bound)
      do b = c0, c1
         call slice_surface_plane(s0, idx, b - c0 + 1_I4P, secs(b - c0 + 1_I4P))
         xg(:, b - c0 + 1_I4P) = ue - nsg(b) * dl(2) ; yg(:, b - c0 + 1_I4P) = ve - nsg(b) * dl(3)
      enddo
      call cut_cells_device_many(secs, xg, yg, cb, st, intervals=.false., pieces=.false., resident=.true.)
      deallocate(secs, xg, yg)
      call translation_reduce_dev(nu, nv, c1 - c0 + 1_I4P, cb%alpha, c0, npl(c0), npl(c1), np, pfirst, nn, nhw, acc)
      call cb%exit_device
      c0 = c1 + 1_I4P
   enddo
   !$acc end data
   do p = 1_I4P, np ! the face sizes, as `plane_average`
      do k = 1_I4P, nv
         do j = 1_I4P, nu
            acc(j, k, p) = acc(j, k, p) * (ue(j + 1) - ue(j)) * (ve(k + 1) - ve(k))
         enddo
      enddo
   enddo
   call move_alloc(acc, area)
   endsubroutine translation_axis_device

   subroutine translation_reduce_dev(nu, nv, nbk, alpha, c0, p0, p1, np, pfirst, nn, hw, acc)
   !< Device kernel (issue #100 step 4): the chunk of nodes c0 .. c0 + nbk - 1 (blocks 1 .. nbk of the section cut
   !< cells) added to the face sums of planes p0 .. p1, one thread per (plane, face), the plane's nodes of the chunk in
   !< node order: acc = acc + hw alpha, hw = 0.5 h (`plane_average`). Every array must be present on the device.
   integer(I4P), intent(in)    :: nu, nv, nbk          !< Face cells, chunk blocks.
   real(R8P),    intent(in)    :: alpha(nu, nv, nbk)   !< Section fluid fractions.
   integer(I4P), intent(in)    :: c0, p0, p1, np, nn   !< First node of the chunk, its planes, planes, nodes.
   integer(I4P), intent(in)    :: pfirst(np + 1)       !< Nodes of each plane.
   real(R8P),    intent(in)    :: hw(nn)               !< Node weights 0.5 h.
   real(R8P),    intent(inout) :: acc(nu, nv, np)      !< Sums.
   real(R8P)                   :: sm                   !< Work.
   integer(I4P)                :: p, j, k, b, b0, b1   !< Work.

   !$acc parallel loop gang vector collapse(3) default(present) private(sm, b0, b1) present(alpha, pfirst, hw, acc)
   do p = p0, p1
      do k = 1_I4P, nv
         do j = 1_I4P, nu
            b0 = max(pfirst(p), c0) ; b1 = min(pfirst(p + 1) - 1_I4P, c0 + nbk - 1_I4P)
            sm = acc(j, k, p)
            !$acc loop seq
            do b = b0, b1
               sm = sm + hw(b) * alpha(j, k, b - c0 + 1_I4P)
            enddo
            acc(j, k, p) = sm
         enddo
      enddo
   enddo
   endsubroutine translation_reduce_dev

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
