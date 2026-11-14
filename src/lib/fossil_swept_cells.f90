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
use fossil_curve_slice,        only : slice_surface, slice_planes_t, slice_planes_prepare, slice_surface_plane
use fossil_cut_cells_2d,       only : cut_cell_2d_walk, cut_cell_2d_one_sizes, cut_cell_2d_one_pieces
use fossil_cut_cells_3d,       only : cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, cut_cells_3d_sections, &
                                     cut_cell_3d_one_sizes, cut_cell_3d_one_count, cut_cell_3d_one_pieces, cut_cell_3d_one_totals, &
                                     cut_cell_3d_face_sizes, cut_cell_3d_face_count, cut_cell_3d_face_measures
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
                                             !< this is evaluated only if wider than its events' round-off widths
                                             !< (issue #104, `interval_wide`), else its events act as one; its
                                             !< quadrature nodes may be skipped (`map_piece_assign`).
real(R8P), parameter :: EVENT_RES = 64._R8P * epsilon(1._R8P) !< Piece map (issue #104): the round-off of an
                                                             !< explicit event instant, in units of its formula's magnitudes
                                                             !< over its denominator (see `node_breakpoints`).
real(R8P), parameter :: IZ_GAP = 1e-10_R8P !< `isolate_zeros`: uncertified intervals stop at this width (see there).
integer(I4P), parameter :: IZ_STACK = 40_I4P !< `isolate_zeros`: stack of intervals. An interval is halved only if it is
                                             !< at least IZ_GAP wide; halves of [0, 1] are exact dyadic intervals of
                                             !< width 2^-d, so a halved one has d <= 33 (2^-33 >= IZ_GAP > 2^-34). The
                                             !< left half is examined first, so the stack holds at most one pending right
                                             !< half per level 1..d plus the two new halves: d + 2 <= 35 entries.
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
   type(rigid_pose_2d_t)     :: pose0, pose1 !< End poses.
   real(R8P)                 :: d(2) = 0._R8P   !< Translation over the step (pose1%t - pose0%t).
   real(R8P)                 :: th0 = 0._R8P    !< Angle at pose0.
   real(R8P)                 :: dth = 0._R8P    !< Rotation over the step.
   real(R8P)                 :: xm = 0._R8P     !< Coordinate magnitude (the round-off scale of a crossing).
   real(R8P)                 :: bmin0(2) = 0._R8P !< Box minimum of the pose0 copy.
   real(R8P)                 :: bmin(2) = 0._R8P  !< Box minimum of the body.
   logical                   :: exact = .true.  !< Every multi-track integration converged.
   integer(I4P)              :: skipped = 0_I4P !< Quadrature nodes or pieces skipped in narrow intervals (`map_measures_2d`).
endtype map_ctx_t

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
   type(rigid_pose_t)        :: pose0, pose1    !< End poses.
   real(R8P)                 :: d(3) = 0._R8P   !< Translation over the step (pose1%t - pose0%t).
   real(R8P)                 :: kax(3) = 0._R8P !< Rotation axis of the step (unit).
   real(R8P)                 :: ang = 0._R8P    !< Rotation angle of the step.
   real(R8P)                 :: om(3) = 0._R8P  !< Rotation vector ang kax.
   real(R8P)                 :: xm = 0._R8P     !< Coordinate magnitude (the round-off scale of a crossing).
   real(R8P)                 :: tol = 0._R8P    !< Position tolerance of an event check, 1e-12 (max |grid| + 1).
   logical                   :: exact = .true.  !< Every multi-track integration converged.
   integer(I4P)              :: skipped = 0_I4P !< Quadrature nodes or pieces skipped in narrow intervals (`map_measures_3d`).
endtype map3_ctx_t

integer(I4P), parameter :: PLAN_UNTRUSTED = -1_I4P !< Map plan (issue #105 step 1c): the cell is untrusted.
integer(I4P), parameter :: PLAN_TRIVIAL   = 0_I4P  !< Map plan: the cell keeps one state, its trivial track (no evaluation).
integer(I4P), parameter :: PLAN_EVAL      = 1_I4P  !< Map plan: the cell is evaluated.

type :: map_plan_t
   !< The host phase of the piece map (issue #105 step 1c), as flat CSR arrays over the cells c (2-D c = i + (j - 1) nx,
   !< 3-D c = i + nx ((j - 1) + ny (k - 1))): every cell's status and, for an evaluated cell (PLAN_EVAL), its local
   !< elements, events, interval layout, evaluation list and (3-D) corner events; a pose table, one placed pose per
   !< distinct instant of a middle evaluation (on the paths that place by a pose: 2-D off the translation path, 3-D
   !< with a rotation), computed by `pose_interpolate` as before; the Gauss node sets of the measures. Built once per
   !< step by `map_plan_2d` / `map_plan_3d` (count per cell, scan, fill; a cell's slices depend on that cell only), read
   !< by the evaluations (`map_evals_2d` / `map_evals_3d`, issue #105 step 1d) and the measures. A slice of cell c is
   !< x(xstart(c) : xstart(c + 1) - 1), empty unless PLAN_EVAL.
   integer(I4P)              :: ncell = 0_I4P    !< Cells.
   integer(I4P)              :: contacts = 0_I4P !< Cells untrusted for an uncertified contact (issue #102).
   integer(I4P)              :: resolved = 0_I4P !< Evaluated cells with a state narrower than EVENT_GAP (issue #104).
   integer(I4P), allocatable :: status(:)        !< PLAN_UNTRUSTED, PLAN_TRIVIAL or PLAN_EVAL [ncell].
   integer(I4P), allocatable :: lstart(:)        !< Local elements [ncell + 1].
   integer(I4P), allocatable :: lel(:)           !< 2-D: local vertices (curve indexes, increasing); 3-D: facets (increasing).
   integer(I4P), allocatable :: lnext(:)         !< 2-D: the next local vertex in the cell's numbering (0: none).
   integer(I4P), allocatable :: bstart(:)        !< Events [ncell + 1].
   real(R8P),    allocatable :: bs(:)            !< Event instants in (0, 1), sorted, distinct.
   real(R8P),    allocatable :: bw(:)            !< Their round-off half-widths (issue #104).
   integer(I4P), allocatable :: istart(:)        !< Intervals between the events (events + 1 per evaluated cell) [ncell + 1].
   logical,      allocatable :: wide(:)          !< The interval is evaluated (`interval_wide`).
   logical,      allocatable :: small(:)         !< The interval is narrower than EVENT_GAP.
   integer(I4P), allocatable :: eref(:)          !< Its evaluation (index in the cell's list).
   integer(I4P), allocatable :: ealt(:)          !< The evaluation after it (= eref for a wide interval).
   integer(I4P), allocatable :: vstart(:)        !< Evaluations: pose0, the wide intervals' middles, pose1 [ncell + 1].
   real(R8P),    allocatable :: vs(:)            !< Their instants (0, 0.5 (sa + sb), 1).
   integer(I4P), allocatable :: vkind(:)         !< 0 the pose0 copy, 1 the pose1 copy, 2 the instant vs.
   integer(I4P), allocatable :: vpose(:)         !< Pose of a kind-2 evaluation in the table (0: placed by s d, or kind 0/1).
   integer(I4P)              :: npose = 0_I4P    !< Poses.
   real(R8P),    allocatable :: pose_s(:)        !< Their instants, increasing, distinct [npose].
   real(R8P),    allocatable :: pose_r(:,:,:)    !< Rotations [dim, dim, npose].
   real(R8P),    allocatable :: pose_t(:,:)      !< Translations [dim, npose].
   integer(I4P), allocatable :: cnstart(:)       !< Corner events (`map_corner_known`) [ncell + 1].
   real(R8P),    allocatable :: cntime(:)        !< Their instants (in the cell's event order).
   integer(I4P), allocatable :: cncor(:)         !< Their corners (1..8).
   real(R8P)                 :: gq_s(126) = 0._R8P !< Gauss-Legendre nodes on [0, 1], n = 2, 4, ..., 64: set n at n - 1 .. 2n - 2.
   real(R8P)                 :: gq_w(126) = 0._R8P !< Their weights.
endtype map_plan_t

type :: plan_cell_t
   !< One candidate cell of the map plan before the scan (`map_plan_2d`, `map_plan_3d`).
   integer(I4P)              :: status = PLAN_UNTRUSTED !< PLAN_UNTRUSTED, PLAN_TRIVIAL or PLAN_EVAL.
   logical                   :: contact = .false.  !< Untrusted for an uncertified contact.
   logical                   :: resolved = .false. !< A wide interval narrower than EVENT_GAP.
   integer(I4P)              :: ne = 0_I4P         !< Evaluations.
   integer(I4P), allocatable :: lg(:), lnx(:)      !< 2-D local vertices and next.
   real(R8P),    allocatable :: b(:), bw(:)        !< Events, half-widths.
   real(R8P),    allocatable :: ct(:), cc(:)       !< Corner events: instants, corners (as reals).
endtype plan_cell_t

integer(I8P), parameter :: MAP_CHUNK_BYTES = 2_I8P**19 !< Map evaluations (issue #105 step 1d): workspace bytes per chunk
                                                       !< of items (one item at least). On the host a chunk that stays in
                                                       !< cache (64 MB chunks made the 3-D map 15-34 % slower); a device
                                                       !< launch would take far larger ones. No output depends on it.

type :: map_evals_t
   !< The evaluations of the piece map (issue #105 step 1d) as flat arrays over the items v, the plan's evaluation slots
   !< (`map_plan_t`: v in vstart(c) .. vstart(c + 1) - 1 for cell c, in its evaluation order). Built by count, scan and
   !< fill: the sizes pass (`map_item_2d_count`, `map_item_3d_count`: the placed geometry and the count kernels give
   !< each item its workspace and the capacities of its outputs), the evaluation pass in chunks of items (one workspace
   !< per item; `map_item_2d_eval`, `map_item_3d_eval`: closed, np, nk and the outputs in capacity slots), then the
   !< scans of np and nk and the fill of the keys (CSR, `kstart`). Read by the tracks (`map_tracks`) and the measures.
   integer(I4P)              :: n = 0_I4P        !< Items.
   logical,      allocatable :: closed(:)        !< The cell boundary closed (3-D: the pieces trusted) [n].
   integer(I4P), allocatable :: np(:), nk(:)     !< Pieces, features (0 if not closed) [n].
   integer(I4P), allocatable :: pstart(:)        !< First piece (node of the tracks) of each item [n + 1].
   integer(I4P), allocatable :: kstart(:)        !< First feature of each item [n + 1].
   integer(I4P), allocatable :: key(:), kpiece(:) !< Features (increasing per item; corners first) and their pieces.
   integer(I4P), allocatable :: mstart(:)        !< Piece measures of each item (capacity slots) [n + 1].
   real(R8P),    allocatable :: pm(:)            !< 2-D piece areas (end items), 3-D piece volumes (end items, several).
   logical,      allocatable :: fl(:,:)          !< 2-D corner flags (inside the solid) [4, n].
   integer(I4P), allocatable :: cs(:,:)          !< 3-D corner states (0 solid, 1 fluid) [8, n].
endtype map_evals_t

integer(I4P), parameter :: PAIR_ACTIVE  = 0_I4P !< Map measures (issue #105 step 1e): the pair takes the next round.
integer(I4P), parameter :: PAIR_DONE    = 1_I4P !< Converged (or at 64 nodes, or without rotation after its round).
integer(I4P), parameter :: PAIR_FAILED  = 2_I4P !< A node could not be assigned: its cell is not measured.
integer(I4P), parameter :: PAIR_DROPPED = 3_I4P !< After a failed pair of its cell (not integrated).

type :: map_pairs_t
   !< The (cell, interval) pairs of the measures of the piece map (issue #105 step 1e): every interval of positive width
   !< of every measured cell, in cell order then interval order (measured cell k = cell(k) has the pairs kstart(k) ..
   !< kstart(k + 1) - 1), with its state over the rounds of `map_measures_2d` / `map_measures_3d`.
   integer(I4P)              :: n = 0_I4P        !< Pairs.
   integer(I4P)              :: nk = 0_I4P       !< Measured cells.
   integer(I4P), allocatable :: cell(:)          !< Cell of each measured cell [nk].
   integer(I4P), allocatable :: kstart(:)        !< Its pairs [nk + 1].
   integer(I4P), allocatable :: k(:)             !< Measured cell of each pair [n].
   real(R8P),    allocatable :: sa(:), sb(:)     !< Ends of its interval [n].
   integer(I4P), allocatable :: vr(:), va(:)     !< Its evaluation and the one after it (items of `map_evals_t`) [n].
   logical,      allocatable :: narrow(:)        !< vr and va differ, or the interval is narrower than EVENT_GAP [n].
   integer(I4P), allocatable :: col(:)           !< First column of its sums (ntr of its cell) [n].
   integer(I4P), allocatable :: state(:)         !< PAIR_ACTIVE, PAIR_DONE, PAIR_FAILED or PAIR_DROPPED [n].
   integer(I4P), allocatable :: skip(:)          !< Nodes or pieces skipped (narrow) [n].
   logical,      allocatable :: inexact(:)       !< Not converged to round-off at 64 nodes [n].
endtype map_pairs_t

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
   !< refined by bisection to one ulp. An interval narrower than IZ_GAP = 1e-10 without a certificate (a zero at an end of
   !< the step, a tangency, g within round-off of zero) gives one breakpoint at its middle: a kink within IZ_GAP of a
   !< breakpoint changes the Gauss value of its piece by ~slope IZ_GAP^2, far below round-off. Zeros within IZ_GAP of 0 or
   !< 1, or of each other, are merged by the callers (`merge_close`). An interval on which |g| is at round-off
   !< throughout (end values and the curvature bound below `err`: a zero zone, e.g. a mesh face starting on a grid
   !< plane) gets breakpoints at its two ends only. A g with every coefficient at
   !< round-off of `scale` (identically zero: the element never leaves the face's line or plane) has no event. At most 4096
   !< intervals are examined: beyond, `certified` is false and the caller must not claim exactness. Its helpers are the
   !< module procedures `iz_g` (g) and `iz_add` (zeros grown by doubling, trimmed at the end; issue #105 step 1a); the
   !< interval stack has IZ_STACK entries, its depth bound; were it full, the interval would give an uncertified middle.
   !< `uncertain` (issue #102, the piece map): per zero, true if no sign change certifies it (a zero zone's end, the
   !< middle of an interval narrower than IZ_GAP, or of one the stack cannot split) and |g| exceeds 4 err on both sides of
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
   integer(I4P), parameter             :: BUDGET = 4096
   real(R8P)                           :: big_m, err, st(2, IZ_STACK), sa, sb, ga, gb, sm, gm, lo, hi, glo
   integer(I4P)                        :: top, used
   logical,      allocatable           :: uf(:) !< Flags of the zeros (capacity of z; returned with `uncertain`).
   real(R8P),    allocatable           :: wf(:) !< Widths of the zeros (capacity of z; returned with `width`).
   real(R8P),    allocatable           :: tz(:) !< Trim of z.
   real(R8P)                           :: slp   !< Slope of g over a bracket.
   integer(I4P)                        :: kz, j !< Work (`uncertain`).
   logical                             :: lft, rgt !< |g| above err before and after an uncertified zero.

   allocate(z(0), uf(0), wf(0)) ; nz = 0_I4P ; certified = .true. ; used = 0_I4P
   if (present(uncertain)) allocate(uncertain(0))
   if (present(width)) allocate(width(0))
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
      ga = iz_g(co, phi0, dphi, sa) ; gb = iz_g(co, phi0, dphi, sb)
      if (ga * gb < 0._R8P .and. abs(ga) > err .and. abs(gb) > err .and. big_m * (sb - sa)**2 / 8._R8P < err) then
         ! one sign change, g nearly linear here: bisection to one ulp
         lo = sa ; hi = sb ; glo = ga ; slp = abs(gb - ga) / (sb - sa)
         do while (nearest(lo, 1._R8P) < hi)
            sm = 0.5_R8P * (lo + hi)
            if (sm <= lo .or. sm >= hi) exit
            gm = iz_g(co, phi0, dphi, sm)
            if (gm == 0._R8P) then
               lo = sm ; hi = sm ; exit
            endif
            if ((gm > 0._R8P) .eqv. (glo > 0._R8P)) then
               lo = sm ; glo = gm
            else
               hi = sm
            endif
         enddo
         call iz_add(0.5_R8P * (lo + hi), .false., 0.5_R8P * (hi - lo) + err / max(slp, tiny(1._R8P)), z, uf, wf, &
                     nz)
         cycle
      endif
      if (ga * gb > 0._R8P .and. min(abs(ga), abs(gb)) > big_m * (sb - sa)**2 / 8._R8P + err) cycle ! no zero
      if (abs(ga) <= err .and. abs(gb) <= err .and. big_m * (sb - sa)**2 / 8._R8P <= err) then
         ! a zero zone: |g| at round-off over the whole interval, the event undetermined in double precision
         call iz_add(sa, .true., sb - sa, z, uf, wf, nz) ; call iz_add(sb, .true., sb - sa, z, uf, wf, nz) ; cycle
      endif
      if (sb - sa < IZ_GAP) then
         call iz_add(0.5_R8P * (sa + sb), .true., 0.5_R8P * (sb - sa), z, uf, wf, nz) ; cycle
      endif
      if (top + 2_I4P > IZ_STACK) then ! cannot happen (see IZ_STACK); keep a breakpoint rather than drop one
         call iz_add(0.5_R8P * (sa + sb), .true., 0.5_R8P * (sb - sa), z, uf, wf, nz) ; cycle
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
            if (abs(iz_g(co, phi0, dphi, z(kz) * (1._R8P - 0.5_R8P**j))) > 4._R8P * err) lft = .true.
            if (abs(iz_g(co, phi0, dphi, z(kz) + (1._R8P - z(kz)) * (1._R8P - 0.5_R8P**j))) > 4._R8P * err) &
               rgt = .true.
         enddo
         uf(kz) = lft .and. rgt
      enddo
      uncertain = uf(1:nz)
   endif
   if (present(width)) width = wf(1:nz)
   if (size(z, kind=I4P) /= nz) then
      allocate(tz(nz)) ; tz = z(1:nz) ; call move_alloc(tz, z)
   endif
   endsubroutine isolate_zeros

   pure function iz_g(co, phi0, dphi, s) result(v)
   !< The function of `isolate_zeros`: g(s) = (co(1) + co(2) s) cos phi + (co(3) + co(4) s) sin phi + co(5) + co(6) s,
   !< phi = phi0 + s dphi.
   real(R8P), intent(in) :: co(6)      !< Coefficients [p0, p1, q0, q1, k0, k1].
   real(R8P), intent(in) :: phi0, dphi !< Angle at s = 0, its change over the step.
   real(R8P), intent(in) :: s          !< Instant.
   real(R8P)             :: v          !< g(s).
   real(R8P)             :: ph         !< Angle.

   ph = phi0 + s * dphi
   v = (co(1) + co(2) * s) * cos(ph) + (co(3) + co(4) * s) * sin(ph) + co(5) + co(6) * s
   endfunction iz_g

   pure subroutine iz_add(s, u, w, z, uf, wf, nz)
   !< Append a zero of `isolate_zeros` (kept only inside (IZ_GAP, 1 - IZ_GAP)) with its flag and width; z, uf and wf
   !< share one capacity, doubled when full (at least 8; the caller trims z to nz).
   real(R8P),              intent(in)    :: s     !< Zero.
   logical,                intent(in)    :: u     !< Not certified by a sign change.
   real(R8P),              intent(in)    :: w     !< Localisation half-width.
   real(R8P), allocatable, intent(inout) :: z(:)  !< Zeros [capacity].
   logical,   allocatable, intent(inout) :: uf(:) !< Flags [capacity].
   real(R8P), allocatable, intent(inout) :: wf(:) !< Widths [capacity].
   integer(I4P),           intent(inout) :: nz    !< Zeros so far.
   real(R8P), allocatable                :: tz(:), tw(:) !< Grown z, wf.
   logical,   allocatable                :: tu(:) !< Grown uf.
   integer(I4P)                          :: cap   !< New capacity.

   if (.not. (s > IZ_GAP .and. s < 1._R8P - IZ_GAP)) return
   if (nz >= size(z, kind=I4P)) then
      cap = max(8_I4P, 2_I4P * size(z, kind=I4P))
      allocate(tz(cap), tu(cap), tw(cap))
      tz(1:nz) = z(1:nz) ; tu(1:nz) = uf(1:nz) ; tw(1:nz) = wf(1:nz)
      call move_alloc(tz, z) ; call move_alloc(tu, uf) ; call move_alloc(tw, wf)
   endif
   nz = nz + 1_I4P ; z(nz) = s ; uf(nz) = u ; wf(nz) = w
   endsubroutine iz_add

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
   !<    The host phase (`map_plan_2d`, issue #105 step 1c) then gives every cell its status and, if it is evaluated,
   !<    its local curve, events, interval layout and evaluation instants, as flat CSR arrays, with a pose table.
   !< 2. Topology (from the plan): the events of a cell's four faces from its local segments (a wall vertex crossing a
   !<    face line inside the face, a corner crossing a wall segment: the translation formulas of
   !<    `faces_translation_2d`, or `line_events_2d` / `node_events_2d` certified by `isolate_zeros`); a cell without
   !<    events whose two end cuts have the same single piece (or none) takes the trivial track. Otherwise its pieces
   !<    are evaluated at pose0 (the end-pose copy, k0's corner flags: they must equal k0's pieces bit for bit, which
   !<    fixes the numbering), at the middle of every interval between events (corner flags by the flip rule where it
   !<    holds, `map_flags_2d`), and at pose1: every (cell, evaluation) one item of the flat evaluations
   !<    (`map_evals_2d`, issue #105 step 1d). The end evaluations are checked (`map_end_check`), and pieces of
   !<    consecutive evaluations sharing a boundary feature (a wall segment or a fluid corner) are one track
   !<    (`map_tracks`, `track_components`).
   !< 3. Measures: one track takes the cell's `sw` values and face outflows (fou, fov); several tracks integrate their
   !<    face parts on the cell's event intervals (`map_measures_2d`: by rounds over flat (cell, interval, node) items,
   !<    issue #105 step 1e).
   !< 4. Assembly (issue #105 step 1f): tracks per cell and their scan (`map_asm_tracks`), the old and new pieces per
   !<    track by count / scan / fill per cell (`map_asm_pieces`), the statuses and measures per cell (`map_asm_status`,
   !<    `map_asm_measures_*`); the same bits as the former loops in cell and piece order.
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
   type(map_plan_t)          :: plan                   !< Host phase (issue #105 step 1c).
   type(map_evals_t)         :: ep                     !< Evaluations (issue #105 step 1d).
   integer(I4P), allocatable :: cstart(:), cseg(:)     !< Local segments of every cell.
   integer(I4P), allocatable :: ntr(:), ckind(:)       !< Tracks and kind (0 trivial, 1 evaluated, -1 untrusted) per cell.
   integer(I4P), allocatable :: otr(:), ntn(:)         !< Track (in its cell) of every old and new piece.
   integer(I4P), allocatable :: comp(:)                !< Track (in its cell) of every piece of every evaluation.
   integer(I4P), allocatable :: tw(:)                  !< Work of `map_tracks`.
   integer(I4P), allocatable :: moff(:)                !< First measure of every measured cell's tracks.
   real(R8P),    allocatable :: meas(:,:)              !< Per-track measures of the measured cells.
   logical,      allocatable :: cok(:)                 !< Measured cells whose measures are assigned.
   integer(I4P)              :: nx, ny, ncell, c, i, j, nt, tot, t, a0, a1, b0, b1, v0, ne, p0, nn, q0, nk !< Work.
   logical                   :: allm                   !< all_meas.

   allm = .false. ; if (present(all_meas)) allm = all_meas
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; ncell = nx * ny
   map%nx = nx ; map%ny = ny ; map%np0 = k0%np ; map%np1 = k1%np
   call map_context(curve, c0, c1, pose0, pose1, x, y, path, nq, mc)
   call map_candidates(mc, x, y, cstart, cseg)
   call map_plan_2d(mc, x, y, k0, k1, cstart, cseg, plan)
   deallocate(cstart, cseg)
   allocate(ntr(ncell), ckind(ncell), otr(k0%np), ntn(k1%np), moff(ncell))
   ! the evaluations of every evaluated cell, flat (issue #105 step 1d), then the tracks of every cell
   call map_flags_2d(mc, plan, curve, c0, x, y, k0, k1, ep)
   call map_evals_2d(mc, plan, x, y, ep)
   allocate(comp(ep%pstart(ep%n + 1) - 1_I4P), tw(map_tracks_work(plan, ep)))
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         c = i + (j - 1_I4P) * nx
         a0 = k0%first_piece(i, j) ; a1 = a0 + k0%pieces(i, j) - 1_I4P
         b0 = k1%first_piece(i, j) ; b1 = b0 + k1%pieces(i, j) - 1_I4P
         ntr(c) = 0_I4P
         select case (plan%status(c))
         case (PLAN_UNTRUSTED)
            ckind(c) = -1_I4P
         case (PLAN_TRIVIAL)
            ckind(c) = 0_I4P ; ntr(c) = k0%pieces(i, j) ; otr(a0:a1) = 1_I4P ; ntn(b0:b1) = 1_I4P
         case default
            v0 = plan%vstart(c) ; ne = plan%vstart(c + 1) - v0
            map%evaluations = map%evaluations + ne ; map%evaluated_cells = map%evaluated_cells + 1_I4P
            ckind(c) = map_end_check(ep, v0, ne, k0%pieces(i, j), k1%pieces(i, j), k0%piece_area(a0:a1), &
                                     k1%piece_area(b0:b1))
            if (ckind(c) == 1_I4P) then ! tracks: pieces of consecutive evaluations sharing a feature
               p0 = ep%pstart(v0) ; nn = ep%pstart(v0 + ne) - p0 ; q0 = ep%kstart(v0) ; nk = ep%kstart(v0 + ne) - q0
               call map_tracks(ne, ep%pstart(v0:v0 + ne), ep%kstart(v0:v0 + ne), nk, ep%key(q0:q0 + nk - 1_I4P), &
                               ep%kpiece(q0:q0 + nk - 1_I4P), nn, tw, comp(p0:p0 + nn - 1_I4P), ntr(c))
               otr(a0:a1) = comp(ep%pstart(v0):ep%pstart(v0 + 1_I4P) - 1_I4P)
               ntn(b0:b1) = comp(ep%pstart(v0 + ne - 1_I4P):ep%pstart(v0 + ne) - 1_I4P)
            endif
         endselect
      enddo
   enddo
   ! the per-track measures of the multi-track cells (all evaluated cells with all_meas), from the flat evaluations
   tot = 0_I4P ; moff = 0_I4P
   do c = 1_I4P, ncell
      if (ckind(c) == 1_I4P .and. (ntr(c) >= 2_I4P .or. allm)) then
         moff(c) = tot + 1_I4P ; tot = tot + ntr(c)
      endif
   enddo
   allocate(meas(8, tot), cok(ncell)) ; meas = 0._R8P
   call map_measures_2d(mc, plan, ep, comp, curve, c0, x, y, ntr, moff, meas, cok)
   do c = 1_I4P, ncell
      if (moff(c) == 0_I4P .or. ntr(c) == 0_I4P) cycle ! not measured, or no track (solid throughout): measured
      if (cok(c)) cycle
      if (ntr(c) == 1_I4P .and. allm) then ! one track (all_meas): keep it with the sw values
         moff(c) = 0_I4P
      else
         ckind(c) = -1_I4P ; moff(c) = 0_I4P
      endif
   enddo
   if (present(measured)) measured = moff > 0_I4P
   ! assembly (issue #105 step 1f): tracks per cell, scan; pieces per track, scan, fill; status; measures
   allocate(map%tracks(nx, ny), map%first_track(nx, ny))
   call map_asm_tracks(ncell, ckind, ntr, map%first_track, nt, map%untrusted_cells, map%multi_track_cells)
   map%tracks = reshape(ntr, [nx, ny]) ; map%nt = nt
   allocate(map%track_cell(2, nt), map%status(nt), map%old_start(nt + 1), map%new_start(nt + 1), map%old_piece(k0%np), &
            map%new_piece(k1%np), map%old_track(k0%np), map%new_track(k1%np), map%volume0(nt), map%volume1(nt), &
            map%faces(4, nt), map%wall(2, nt), map%residual(nt), map%outflow(4, nt))
   call map_asm_pieces(ncell, ckind, ntr, map%first_track, k0%first_piece, k0%pieces, k0%np, otr, k0%piece_area, nt, &
                       map%old_start, map%old_piece, map%old_track, map%volume0)
   call map_asm_pieces(ncell, ckind, ntr, map%first_track, k1%first_piece, k1%pieces, k1%np, ntn, k1%piece_area, nt, &
                       map%new_start, map%new_piece, map%new_track, map%volume1)
   !$omp parallel do default(shared) private(c, i, j) schedule(static)
   do c = 1_I4P, ncell
      i = 1_I4P + mod(c - 1_I4P, nx) ; j = 1_I4P + (c - 1_I4P) / nx
      call map_asm_status(2_I4P, [i, j], ckind(c) == -1_I4P, map%first_track(i, j), ntr(c), nt, map%old_start, &
                          map%new_start, map%status, map%track_cell)
      call map_asm_measures_2d(map%first_track(i, j), ntr(c), moff(c), size(meas, 2, kind=I4P), meas, &
                               [sw%face_u(i, j), sw%face_u(i + 1, j), sw%face_v(i, j), sw%face_v(i, j + 1)], &
                               [-fou(i, j), fou(i + 1, j), -fov(i, j), fov(i, j + 1)], sw%wall(:, i, j), sw%residual(i, j), &
                               nt, map%volume0, map%volume1, map%faces, map%outflow, map%wall, map%residual)
   enddo
   !$omp end parallel do
   do t = 1_I4P, nt
      map%n_status(map%status(t)) = map%n_status(map%status(t)) + 1_I4P
   enddo
   map%exact = mc%exact ; map%skipped_nodes = mc%skipped
   map%contact_cells = plan%contacts ; map%narrow_cells = plan%resolved
   endsubroutine piece_map_2d_build

   subroutine map_asm_tracks(ncell, ckind, ntr, ft, nt, nun, nmul)
   !< Assembly of the map (issue #105 step 1f), tracks per cell: an untrusted cell (ckind -1) gets one track with every
   !< piece, then the exclusive scan of the counts gives every cell its first track (cell order); counts of the untrusted
   !< and multi-track cells.
   integer(I4P), intent(in)    :: ncell       !< Cells.
   integer(I4P), intent(in)    :: ckind(ncell) !< Kind per cell (0 trivial, 1 evaluated, -1 untrusted).
   integer(I4P), intent(inout) :: ntr(ncell)  !< Tracks per cell (untrusted: set to 1).
   integer(I4P), intent(out)   :: ft(ncell)   !< First track of every cell.
   integer(I4P), intent(out)   :: nt          !< Tracks.
   integer(I4P), intent(inout) :: nun, nmul   !< Untrusted cells, multi-track cells.
   integer(I4P), allocatable   :: w(:)        !< Counts -> starts [ncell + 1].
   integer(I4P)                :: c, nu, nm   !< Work.

   allocate(w(ncell + 1)) ; nu = 0_I4P ; nm = 0_I4P
   !$omp parallel do default(shared) private(c) reduction(+:nu, nm) schedule(static)
   do c = 1_I4P, ncell
      if (ckind(c) == -1_I4P) then
         ntr(c) = 1_I4P ; nu = nu + 1_I4P
      endif
      if (ntr(c) >= 2_I4P) nm = nm + 1_I4P
      w(c) = ntr(c)
   enddo
   !$omp end parallel do
   call csr_scan(w)
   ft = w(1:ncell) ; nt = w(ncell + 1) - 1_I4P ; nun = nun + nu ; nmul = nmul + nm
   endsubroutine map_asm_tracks

   subroutine map_asm_pieces(ncell, ckind, ntr, ft, fp, pc, np, ltr, pv, nt, start, piece, trk, vol)
   !< Assembly of the map (issue #105 step 1f), the pieces of one end (old or new) per track: the count pass per cell
   !< (`map_asm_count`: every piece's track, every track's count), the scan of the counts into the CSR starts, the fill
   !< pass per cell (`map_asm_fill`: the pieces of every track in increasing order, their volumes summed in that order).
   !< A cell's pieces (fp(c) .. fp(c) + pc(c) - 1) and tracks are contiguous, so every pass writes only its cell's
   !< entries, and the order of every sum is that of the former global loop over the pieces: the same bits.
   integer(I4P), intent(in)    :: ncell          !< Cells.
   integer(I4P), intent(in)    :: ckind(ncell)   !< Kind per cell (-1 untrusted: one track with every piece).
   integer(I4P), intent(in)    :: ntr(ncell)     !< Tracks per cell.
   integer(I4P), intent(in)    :: ft(ncell)      !< First track per cell.
   integer(I4P), intent(in)    :: fp(ncell)      !< First piece per cell.
   integer(I4P), intent(in)    :: pc(ncell)      !< Pieces per cell (-1: untrusted cut cell, no piece).
   integer(I4P), intent(in)    :: np             !< Pieces.
   integer(I4P), intent(in)    :: ltr(np)        !< Track in its cell of every piece (not read in untrusted cells).
   real(R8P),    intent(in)    :: pv(np)         !< Volume (area) of every piece.
   integer(I4P), intent(in)    :: nt             !< Tracks.
   integer(I4P), intent(out)   :: start(nt + 1)  !< CSR starts of the pieces of every track.
   integer(I4P), intent(out)   :: piece(np)      !< Pieces of every track.
   integer(I4P), intent(out)   :: trk(np)        !< Track of every piece.
   real(R8P),    intent(out)   :: vol(nt)        !< Volume of every track.
   integer(I4P)                :: c              !< Work.

   !$omp parallel do default(shared) private(c) schedule(static)
   do c = 1_I4P, ncell
      call map_asm_count(ckind(c) == -1_I4P, ft(c), ntr(c), fp(c), max(0_I4P, pc(c)), np, ltr, nt, trk, start)
   enddo
   !$omp end parallel do
   call csr_scan(start)
   !$omp parallel do default(shared) private(c) schedule(static)
   do c = 1_I4P, ncell
      call map_asm_fill(ft(c), ntr(c), fp(c), max(0_I4P, pc(c)), np, trk, pv, nt, start, piece, vol)
   enddo
   !$omp end parallel do
   endsubroutine map_asm_pieces

   pure subroutine map_asm_count(untrusted, ft, ntc, a0, na, np, ltr, nt, trk, cnt)
   !< Kernel (issue #105 step 1f), count pass of one cell: the track of each of its pieces (its first track + its track in
   !< the cell - 1; the first track for an untrusted cell) and the pieces of each of its tracks (cnt(t), before the scan).
   !$acc routine seq
   logical,      intent(in)    :: untrusted   !< Untrusted cell: one track with every piece.
   integer(I4P), intent(in)    :: ft, ntc     !< First track, tracks of the cell.
   integer(I4P), intent(in)    :: a0, na      !< First piece, pieces of the cell.
   integer(I4P), intent(in)    :: np          !< Pieces.
   integer(I4P), intent(in)    :: ltr(np)     !< Track in its cell of every piece.
   integer(I4P), intent(in)    :: nt          !< Tracks.
   integer(I4P), intent(inout) :: trk(np)     !< Track of every piece.
   integer(I4P), intent(inout) :: cnt(nt + 1) !< Pieces of every track.
   integer(I4P)                :: p, t        !< Work.

   do t = ft, ft + ntc - 1_I4P
      cnt(t) = 0_I4P
   enddo
   do p = a0, a0 + na - 1_I4P
      if (untrusted) then
         t = ft
      else
         t = ft + ltr(p) - 1_I4P
      endif
      trk(p) = t ; cnt(t) = cnt(t) + 1_I4P
   enddo
   endsubroutine map_asm_count

   pure subroutine map_asm_fill(ft, ntc, a0, na, np, trk, pv, nt, start, piece, vol)
   !< Kernel (issue #105 step 1f), fill pass of one cell: for each of its tracks, its pieces in increasing order into the
   !< track's CSR slots and their volumes summed in that order (no cursor array: a cell has few tracks and pieces).
   !$acc routine seq
   integer(I4P), intent(in)    :: ft, ntc       !< First track, tracks of the cell.
   integer(I4P), intent(in)    :: a0, na        !< First piece, pieces of the cell.
   integer(I4P), intent(in)    :: np            !< Pieces.
   integer(I4P), intent(in)    :: trk(np)       !< Track of every piece.
   real(R8P),    intent(in)    :: pv(np)        !< Volume of every piece.
   integer(I4P), intent(in)    :: nt            !< Tracks.
   integer(I4P), intent(in)    :: start(nt + 1) !< CSR starts.
   integer(I4P), intent(inout) :: piece(np)     !< Pieces of every track.
   real(R8P),    intent(inout) :: vol(nt)       !< Volume of every track.
   integer(I4P)                :: p, t, q       !< Work.
   real(R8P)                   :: v             !< Sum.

   do t = ft, ft + ntc - 1_I4P
      q = start(t) ; v = 0._R8P
      do p = a0, a0 + na - 1_I4P
         if (trk(p) /= t) cycle
         piece(q) = p ; q = q + 1_I4P ; v = v + pv(p)
      enddo
      vol(t) = v
   enddo
   endsubroutine map_asm_fill

   pure subroutine map_asm_status(dim, ijk, untrusted, ft, ntc, nt, ostart, nstart, status, track_cell)
   !< Kernel (issue #105 step 1f): the status and cell of every track of one cell (TRACK_UNTRUSTED in an untrusted cell,
   !< else `track_status_of` its old and new piece counts).
   !$acc routine seq
   integer(I4P), intent(in)    :: dim                 !< 2 or 3.
   integer(I4P), intent(in)    :: ijk(dim)            !< Cell.
   logical,      intent(in)    :: untrusted           !< Untrusted cell.
   integer(I4P), intent(in)    :: ft, ntc             !< First track, tracks of the cell.
   integer(I4P), intent(in)    :: nt                  !< Tracks.
   integer(I4P), intent(in)    :: ostart(nt + 1)      !< CSR starts of the old pieces.
   integer(I4P), intent(in)    :: nstart(nt + 1)      !< CSR starts of the new pieces.
   integer(I4P), intent(inout) :: status(nt)          !< Status of every track.
   integer(I4P), intent(inout) :: track_cell(dim, nt) !< Cell of every track.
   integer(I4P)                :: t, a                !< Work.

   do t = ft, ft + ntc - 1_I4P
      do a = 1_I4P, dim
         track_cell(a, t) = ijk(a)
      enddo
      if (untrusted) then
         status(t) = TRACK_UNTRUSTED
      else
         status(t) = track_status_of(ostart(t + 1) - ostart(t), nstart(t + 1) - nstart(t))
      endif
   enddo
   endsubroutine map_asm_status

   pure subroutine map_asm_measures_2d(ft, ntc, mo, nm, meas, cf, co, cw, cr, nt, v0, v1, faces, outflow, wall, residual)
   !< Kernel (issue #105 step 1f): faces, outflows, wall and residual of every track of one 2-D cell, from the per-track
   !< measures (mo > 0: meas(1:4) face lengths, meas(5:8) outflows in +axis of the faces -u, +u, -v, +v) or, one track,
   !< from the cell's values (cf faces, co outward outflows, cw wall, cr residual of `sw`).
   !$acc routine seq
   integer(I4P), intent(in)    :: ft, ntc        !< First track, tracks of the cell.
   integer(I4P), intent(in)    :: mo             !< First measure of the cell's tracks (0: the cell's values).
   integer(I4P), intent(in)    :: nm             !< Measured tracks.
   real(R8P),    intent(in)    :: meas(8, nm)    !< Per-track measures.
   real(R8P),    intent(in)    :: cf(4), co(4)   !< The cell's faces and outward outflows.
   real(R8P),    intent(in)    :: cw(2), cr      !< The cell's wall and residual.
   integer(I4P), intent(in)    :: nt             !< Tracks.
   real(R8P),    intent(in)    :: v0(nt), v1(nt) !< Track volumes at pose0 and pose1.
   real(R8P),    intent(inout) :: faces(4, nt), outflow(4, nt), wall(2, nt), residual(nt) !< Track values.
   integer(I4P)                :: t, p, a        !< Work.

   do t = ft, ft + ntc - 1_I4P
      if (mo > 0_I4P) then ! several tracks: integrated per track
         p = mo + t - ft
         do a = 1_I4P, 4_I4P
            faces(a, t) = meas(a, p)
         enddo
         outflow(1, t) = -meas(5, p) ; outflow(2, t) = meas(6, p) ; outflow(3, t) = -meas(7, p)
         outflow(4, t) = meas(8, p)
         wall(1, t) = -(meas(2, p) - meas(1, p)) ; wall(2, t) = -(meas(4, p) - meas(3, p))
         residual(t) = v1(t) - v0(t) + (meas(6, p) - meas(5, p) + meas(8, p) - meas(7, p))
      else ! one track: the cell's values
         do a = 1_I4P, 4_I4P
            faces(a, t) = cf(a) ; outflow(a, t) = co(a)
         enddo
         wall(1, t) = cw(1) ; wall(2, t) = cw(2) ; residual(t) = cr
      endif
   enddo
   endsubroutine map_asm_measures_2d

   pure subroutine map_asm_measures_3d(ft, ntc, mo, nm, meas, cf, co, cw, cr, nt, v0, v1, faces, outflow, wall, residual)
   !< Kernel (issue #105 step 1f), the 3-D twin of `map_asm_measures_2d`: meas(1:6) face areas, meas(7:12) outward
   !< outflows of the faces -x, +x, -y, +y, -z, +z.
   !$acc routine seq
   integer(I4P), intent(in)    :: ft, ntc        !< First track, tracks of the cell.
   integer(I4P), intent(in)    :: mo             !< First measure of the cell's tracks (0: the cell's values).
   integer(I4P), intent(in)    :: nm             !< Measured tracks.
   real(R8P),    intent(in)    :: meas(12, nm)   !< Per-track measures.
   real(R8P),    intent(in)    :: cf(6), co(6)   !< The cell's faces and outward outflows.
   real(R8P),    intent(in)    :: cw(3), cr      !< The cell's wall and residual.
   integer(I4P), intent(in)    :: nt             !< Tracks.
   real(R8P),    intent(in)    :: v0(nt), v1(nt) !< Track volumes at pose0 and pose1.
   real(R8P),    intent(inout) :: faces(6, nt), outflow(6, nt), wall(3, nt), residual(nt) !< Track values.
   integer(I4P)                :: t, p, a        !< Work.

   do t = ft, ft + ntc - 1_I4P
      if (mo > 0_I4P) then ! several tracks: integrated per track
         p = mo + t - ft
         do a = 1_I4P, 6_I4P
            faces(a, t) = meas(a, p) ; outflow(a, t) = meas(6 + a, p)
         enddo
         wall(1, t) = -(meas(2, p) - meas(1, p)) ; wall(2, t) = -(meas(4, p) - meas(3, p))
         wall(3, t) = -(meas(6, p) - meas(5, p))
         residual(t) = v1(t) - v0(t) + (meas(7, p) + meas(8, p) + meas(9, p) + meas(10, p) + meas(11, p) + meas(12, p))
      else ! one track: the cell's values
         do a = 1_I4P, 6_I4P
            faces(a, t) = cf(a) ; outflow(a, t) = co(a)
         enddo
         do a = 1_I4P, 3_I4P
            wall(a, t) = cw(a)
         enddo
         residual(t) = cr
      endif
   enddo
   endsubroutine map_asm_measures_3d

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
   allocate(mc%bu(mc%nv), mc%bv(mc%nv), mc%u0(mc%nv), mc%v0(mc%nv), mc%u1(mc%nv), mc%v1(mc%nv), mc%nxt(mc%nv))
   do g = 1_I4P, mc%nv
      uv = curve%get_vertex(g) ; mc%bu(g) = uv(1) ; mc%bv(g) = uv(2)
      uv = c0%get_vertex(g) ; mc%u0(g) = uv(1) ; mc%v0(g) = uv(2)
      uv = c1%get_vertex(g) ; mc%u1(g) = uv(1) ; mc%v1(g) = uv(2)
      mc%nxt(g) = curve%get_next(g)
   enddo
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

   subroutine map_plan_2d(mc, x, y, k0, k1, cstart, cseg, plan)
   !< The host phase of the 2-D piece map (issue #105 step 1c; `map_plan_t`). Per cell, the checks of the map that
   !< come before its evaluations: both end cuts closed, a wall near the cell during the step, its events certified and
   !< free of uncertified contacts, more than one state. An evaluated cell also gets its local curve (`map_local`, a
   !< per-cell deduplication), its events (`map_events`), its interval layout and its evaluation list (`plan_layout`).
   !< Then the pose table (`plan_instants`, off the translation path). The candidate cells are independent: each fills
   !< its own `plan_cell_t` (an OpenMP loop when built with OpenMP), then `plan_fill` scans and copies them in cell order,
   !< so the plan does not depend on the threads.
   type(map_ctx_t),      intent(in)  :: mc                 !< State.
   real(R8P),            intent(in)  :: x(:), y(:)         !< Cell edges.
   type(cut_cells_2d_t), intent(in)  :: k0, k1             !< End-pose cuts.
   integer(I4P),         intent(in)  :: cstart(:), cseg(:) !< Local segments of every cell (`map_candidates`).
   type(map_plan_t),     intent(out) :: plan               !< Plan.
   type(plan_cell_t), allocatable :: pc(:)                 !< Candidate cells.
   integer(I4P),      allocatable :: cand(:)               !< Their cells.
   type(rigid_pose_2d_t)          :: ps                    !< Pose.
   integer(I4P)                   :: nx, ny, ncell, nc, ic, c, i, j, np0, np1, nb, k !< Work.
   logical                        :: cert, unc             !< Events certified, an uncertified contact.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; ncell = nx * ny
   plan%ncell = ncell ; allocate(plan%status(ncell))
   nc = 0_I4P
   do c = 1_I4P, ncell
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      np0 = k0%pieces(i, j) ; np1 = k1%pieces(i, j)
      if (.not. (k0%closed(i, j) .and. k1%closed(i, j))) then
         plan%status(c) = PLAN_UNTRUSTED
      elseif (cstart(c + 1) == cstart(c)) then ! no wall near the cell during the step: the same state throughout
         plan%status(c) = merge(PLAN_TRIVIAL, PLAN_UNTRUSTED, np0 == np1 .and. np0 <= 1_I4P)
      else
         plan%status(c) = PLAN_EVAL ; nc = nc + 1_I4P
      endif
   enddo
   allocate(cand(nc), pc(nc))
   nc = 0_I4P
   do c = 1_I4P, ncell
      if (plan%status(c) /= PLAN_EVAL) cycle
      nc = nc + 1_I4P ; cand(nc) = c
   enddo
   !$omp parallel do default(shared) private(ic, c, i, j, nb, cert, unc) schedule(dynamic, 4)
   do ic = 1_I4P, nc ! independent cells, each into its own pc(ic)
      c = cand(ic) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      call map_local(mc%nxt, cseg(cstart(c):cstart(c + 1) - 1_I4P), pc(ic)%lg, pc(ic)%lnx)
      call map_events(mc, x(i), x(i + 1), y(j), y(j + 1), size(pc(ic)%lg, kind=I4P), pc(ic)%lg, pc(ic)%lnx, pc(ic)%b, &
                      pc(ic)%bw, nb, cert, unc, pc(ic)%ct, pc(ic)%cc)
      call plan_classify(pc(ic), cert, unc, k0%pieces(i, j), k1%pieces(i, j))
   enddo
   !$omp end parallel do
   call plan_fill(plan, cand, pc, 2_I4P)
   if (mc%path /= MAP_TRANSLATION) then ! the middle evaluations placed by the screw pose (`map_place_2d`)
      call plan_instants(plan, 2_I4P)
      do k = 1_I4P, plan%npose
         ps = pose_interpolate(mc%pose0, mc%pose1, plan%pose_s(k))
         plan%pose_r(:, :, k) = ps%r ; plan%pose_t(:, k) = ps%t
      enddo
   endif
   endsubroutine map_plan_2d

   subroutine plan_classify(pc, cert, unc, np0, np1)
   !< The status of a candidate cell from its events (the map's checks, in their order): untrusted if an event is not
   !< certified or is an uncertified contact; trivial without events and with the same single piece (or none) at both
   !< ends; evaluated otherwise, with its evaluations counted and its narrow states flagged (`plan_layout`). Only an
   !< evaluated cell keeps its local elements and events.
   type(plan_cell_t), intent(inout) :: pc        !< Candidate cell (events in pc%b).
   logical,           intent(in)    :: cert, unc !< Events certified, an uncertified contact.
   integer(I4P),      intent(in)    :: np0, np1  !< Pieces at the two ends.
   logical,      allocatable        :: wide(:), small(:) !< Interval layout.
   integer(I4P), allocatable        :: eref(:), ealt(:)  !< Interval layout.
   integer(I4P)                     :: nb        !< Events.

   nb = size(pc%b, kind=I4P)
   if (.not. cert) then
      pc%status = PLAN_UNTRUSTED
   elseif (unc) then ! an uncertified contact: a state the evaluations cannot decide (issue #102)
      pc%status = PLAN_UNTRUSTED ; pc%contact = .true.
   elseif (nb == 0_I4P .and. np0 == np1 .and. np0 <= 1_I4P) then ! one state throughout, at most one piece
      pc%status = PLAN_TRIVIAL
   else
      pc%status = PLAN_EVAL
      allocate(wide(nb + 1), small(nb + 1), eref(nb + 1), ealt(nb + 1))
      call plan_layout(pc%b, pc%bw, nb, wide, small, eref, ealt, pc%ne)
      pc%resolved = any(wide .and. small)
   endif
   if (pc%status /= PLAN_EVAL) then
      if (allocated(pc%lg)) deallocate(pc%lg)
      if (allocated(pc%lnx)) deallocate(pc%lnx)
      if (allocated(pc%ct)) deallocate(pc%ct, pc%cc)
      deallocate(pc%b, pc%bw)
   endif
   endsubroutine plan_classify

   pure subroutine plan_layout(b, bw, nb, wide, small, eref, ealt, ne)
   !< The interval layout of an evaluated cell (issue #102): the events b(1:nb) split [0, 1] into nb + 1 intervals; a
   !< wide one (`interval_wide`) has an evaluation of its own at its middle (2 .. ne - 1, ne the evaluations with pose0
   !< and pose1); a narrow one refers to the evaluations before (eref) and after (ealt) it.
   real(R8P),    intent(in)  :: b(:), bw(:)     !< Events, their half-widths.
   integer(I4P), intent(in)  :: nb              !< Events.
   logical,      intent(out) :: wide(:), small(:) !< Per interval: evaluated; narrower than EVENT_GAP [nb + 1].
   integer(I4P), intent(out) :: eref(:), ealt(:)  !< Per interval: its evaluation, the one after it [nb + 1].
   integer(I4P), intent(out) :: ne              !< Evaluations.
   real(R8P)                 :: sa, sb          !< Interval.
   integer(I4P)              :: m, e            !< Work.

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
   ne = ne + 1_I4P
   e = ne
   do m = nb + 1_I4P, 1_I4P, -1_I4P
      if (ealt(m) == 0_I4P) then
         ealt(m) = e
      else
         e = eref(m)
      endif
   enddo
   endsubroutine plan_layout

   subroutine plan_fill(plan, cand, pc, dim)
   !< The flat arrays of the plan from its candidate cells (`map_plan_2d`, `map_plan_3d`): the counts per cell, an
   !< exclusive scan into the CSR starts, then each evaluated cell's slices (local elements, events, interval layout,
   !< evaluation list: pose0, the wide intervals' middles 0.5 (sa + sb), pose1; 3-D corner events); the counters; the
   !< Gauss node sets n = 2, 4, ..., 64 of the measures.
   type(map_plan_t),  intent(inout) :: plan    !< Plan (status set).
   integer(I4P),      intent(in)    :: cand(:) !< Candidate cells.
   type(plan_cell_t), intent(in)    :: pc(:)   !< Their host phase.
   integer(I4P),      intent(in)    :: dim     !< 2 or 3.
   real(R8P), allocatable           :: sq(:), wq(:) !< Gauss nodes, weights.
   real(R8P)                        :: sa, sb  !< Interval.
   integer(I4P)                     :: ncell, ic, c, l0, b0, i0, v0, k0, nl, nb, ncn, m, e, n, ne !< Work.

   ncell = plan%ncell
   allocate(plan%lstart(ncell + 1), plan%bstart(ncell + 1), plan%istart(ncell + 1), plan%vstart(ncell + 1), &
            plan%cnstart(ncell + 1))
   plan%lstart = 0_I4P ; plan%bstart = 0_I4P ; plan%istart = 0_I4P ; plan%vstart = 0_I4P ; plan%cnstart = 0_I4P
   plan%contacts = 0_I4P ; plan%resolved = 0_I4P
   do ic = 1_I4P, size(cand, kind=I4P)
      c = cand(ic) ; plan%status(c) = pc(ic)%status
      if (pc(ic)%contact) plan%contacts = plan%contacts + 1_I4P
      if (pc(ic)%status /= PLAN_EVAL) cycle
      if (pc(ic)%resolved) plan%resolved = plan%resolved + 1_I4P
      plan%lstart(c) = size(pc(ic)%lg, kind=I4P) ; plan%bstart(c) = size(pc(ic)%b, kind=I4P)
      plan%istart(c) = plan%bstart(c) + 1_I4P ; plan%vstart(c) = pc(ic)%ne
      if (allocated(pc(ic)%ct)) plan%cnstart(c) = size(pc(ic)%ct, kind=I4P)
   enddo
   call csr_scan(plan%lstart) ; call csr_scan(plan%bstart) ; call csr_scan(plan%istart) ; call csr_scan(plan%vstart)
   call csr_scan(plan%cnstart)
   allocate(plan%lel(plan%lstart(ncell + 1) - 1_I4P), plan%bs(plan%bstart(ncell + 1) - 1_I4P), &
            plan%bw(plan%bstart(ncell + 1) - 1_I4P), plan%wide(plan%istart(ncell + 1) - 1_I4P), &
            plan%small(plan%istart(ncell + 1) - 1_I4P), plan%eref(plan%istart(ncell + 1) - 1_I4P), &
            plan%ealt(plan%istart(ncell + 1) - 1_I4P), plan%vs(plan%vstart(ncell + 1) - 1_I4P), &
            plan%vkind(plan%vstart(ncell + 1) - 1_I4P), plan%vpose(plan%vstart(ncell + 1) - 1_I4P), &
            plan%cntime(plan%cnstart(ncell + 1) - 1_I4P), plan%cncor(plan%cnstart(ncell + 1) - 1_I4P))
   if (dim == 2_I4P) allocate(plan%lnext(plan%lstart(ncell + 1) - 1_I4P))
   plan%vpose = 0_I4P
   do ic = 1_I4P, size(cand, kind=I4P)
      if (pc(ic)%status /= PLAN_EVAL) cycle
      c = cand(ic)
      l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
      b0 = plan%bstart(c) ; nb = plan%bstart(c + 1) - b0
      i0 = plan%istart(c) ; v0 = plan%vstart(c) ; k0 = plan%cnstart(c) ; ncn = plan%cnstart(c + 1) - k0
      plan%lel(l0:l0 + nl - 1_I4P) = pc(ic)%lg
      if (dim == 2_I4P) plan%lnext(l0:l0 + nl - 1_I4P) = pc(ic)%lnx
      plan%bs(b0:b0 + nb - 1_I4P) = pc(ic)%b ; plan%bw(b0:b0 + nb - 1_I4P) = pc(ic)%bw
      call plan_layout(plan%bs(b0:b0 + nb - 1_I4P), plan%bw(b0:b0 + nb - 1_I4P), nb, plan%wide(i0:i0 + nb), &
                       plan%small(i0:i0 + nb), plan%eref(i0:i0 + nb), plan%ealt(i0:i0 + nb), ne)
      plan%vkind(v0) = 0_I4P ; plan%vs(v0) = 0._R8P
      e = v0
      do m = 1_I4P, nb + 1_I4P
         if (.not. plan%wide(i0 + m - 1_I4P)) cycle
         call map_interval(plan%bs(b0:b0 + nb - 1_I4P), nb, m, sa, sb)
         e = e + 1_I4P ; plan%vkind(e) = 2_I4P ; plan%vs(e) = 0.5_R8P * (sa + sb)
      enddo
      plan%vkind(v0 + ne - 1_I4P) = 1_I4P ; plan%vs(v0 + ne - 1_I4P) = 1._R8P
      if (ncn > 0_I4P) then
         plan%cntime(k0:k0 + ncn - 1_I4P) = pc(ic)%ct ; plan%cncor(k0:k0 + ncn - 1_I4P) = nint(pc(ic)%cc, I4P)
      endif
   enddo
   do m = 1_I4P, 6_I4P
      n = 2_I4P**m
      call gauss_legendre_01(n, sq, wq)
      plan%gq_s(n - 1_I4P:2_I4P * n - 2_I4P) = sq ; plan%gq_w(n - 1_I4P:2_I4P * n - 2_I4P) = wq
   enddo
   plan%npose = 0_I4P
   allocate(plan%pose_s(0), plan%pose_r(dim, dim, 0), plan%pose_t(dim, 0))
   endsubroutine plan_fill

   pure subroutine csr_scan(a)
   !< Counts a(1:n) to CSR starts in place: a(1) = 1, a(k + 1) = a(k) + count k (a(n + 1) on input is ignored).
   integer(I4P), intent(inout) :: a(:) !< Counts [n + 1] -> starts.
   integer(I4P)                :: k, s, t !< Work.

   s = 1_I4P
   do k = 1_I4P, size(a, kind=I4P) - 1_I4P
      t = a(k) ; a(k) = s ; s = s + t
   enddo
   a(size(a)) = s
   endsubroutine csr_scan

   subroutine plan_instants(plan, dim)
   !< The pose table's instants: the distinct instants of the kind-2 evaluations, increasing (exact equality), each
   !< evaluation pointing to its own (vpose). The caller places the poses.
   type(map_plan_t), intent(inout) :: plan  !< Plan (evaluations filled).
   integer(I4P),     intent(in)    :: dim   !< 2 or 3.
   real(R8P),    allocatable       :: key(:,:) !< Instants.
   integer(I4P), allocatable       :: idx(:), ord(:) !< Their evaluations, sorted order.
   integer(I4P)                    :: n, v, r !< Work.

   n = count(plan%vkind == 2_I4P, kind=I4P)
   allocate(key(1, n), idx(n))
   n = 0_I4P
   do v = 1_I4P, size(plan%vkind, kind=I4P)
      if (plan%vkind(v) /= 2_I4P) cycle
      n = n + 1_I4P ; key(1, n) = plan%vs(v) ; idx(n) = v
   enddo
   call sort_lex(key, ord)
   plan%npose = 0_I4P
   do r = 1_I4P, n
      if (r == 1_I4P) then
         plan%npose = 1_I4P
      elseif (key(1, ord(r)) /= key(1, ord(r - 1_I4P))) then
         plan%npose = plan%npose + 1_I4P
      endif
      plan%vpose(idx(ord(r))) = plan%npose
   enddo
   deallocate(plan%pose_s, plan%pose_r, plan%pose_t)
   allocate(plan%pose_s(plan%npose), plan%pose_r(dim, dim, plan%npose), plan%pose_t(dim, plan%npose))
   do r = 1_I4P, n
      plan%pose_s(plan%vpose(idx(ord(r)))) = key(1, ord(r))
   enddo
   endsubroutine plan_instants

   pure function plan_pose_2d(plan, k) result(ps)
   !< Pose k of the plan's table (the identity for k = 0).
   type(map_plan_t), intent(in) :: plan !< Plan.
   integer(I4P),     intent(in) :: k    !< Pose (0: none).
   type(rigid_pose_2d_t)        :: ps   !< Pose.

   if (k > 0_I4P) then
      ps%r = plan%pose_r(:, :, k) ; ps%t = plan%pose_t(:, k)
   endif
   endfunction plan_pose_2d

   subroutine map_flags_2d(mc, plan, curve, c0, x, y, k0, k1, ep)
   !< The corner flags of every evaluation (issue #105 step 1d; `map_evals_t`): pose0 and pose1 take the end cuts'
   !< node flags (k0, k1); a middle one, at s, the flip rule (`map_corner_known`, as the 3-D map): a corner whose events
   !< are each one transversal crossing of a wall segment (at least 1e-9 apart) and whose parity takes k0's flag to
   !< k1's has the flag of pose0 flipped once per event before s; any other corner takes the exact query
   !< `map_corner_in` (the screw pose of the plan's table off the translation path). Host.
   type(map_ctx_t),      intent(in)    :: mc         !< State.
   type(map_plan_t),     intent(in)    :: plan       !< Host phase.
   type(curve_object),   intent(in)    :: curve, c0  !< Curve (body frame), copy at pose0.
   real(R8P),            intent(in)    :: x(:), y(:) !< Cell edges.
   type(cut_cells_2d_t), intent(in)    :: k0, k1     !< End-pose cuts.
   type(map_evals_t),    intent(inout) :: ep         !< Evaluations (fl set).
   real(R8P)                           :: cn(2, 4), s !< Corners, instant.
   logical                             :: f0(4), f1(4), known(4) !< End flags, flip rule.
   integer(I4P)                        :: nx, c, i, j, v, q, k, n, n0, n1, cs0(4), cs1(4) !< Work.

   nx = size(x, kind=I4P) - 1_I4P
   ep%n = plan%vstart(plan%ncell + 1) - 1_I4P
   allocate(ep%fl(4, ep%n))
   do c = 1_I4P, plan%ncell
      if (plan%status(c) /= PLAN_EVAL) cycle
      i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
      cn = reshape([x(i), y(j), x(i + 1), y(j), x(i + 1), y(j + 1), x(i), y(j + 1)], [2, 4])
      f0 = [k0%node_in(i, j), k0%node_in(i + 1, j), k0%node_in(i + 1, j + 1), k0%node_in(i, j + 1)]
      f1 = [k1%node_in(i, j), k1%node_in(i + 1, j), k1%node_in(i + 1, j + 1), k1%node_in(i, j + 1)]
      cs0 = merge(0_I4P, 1_I4P, f0) ; cs1 = merge(0_I4P, 1_I4P, f1)
      n0 = plan%cnstart(c) ; n1 = plan%cnstart(c + 1) - 1_I4P
      call map_corner_known(4_I4P, plan%cntime(n0:n1), plan%cncor(n0:n1), cs0, cs1, known)
      do v = plan%vstart(c), plan%vstart(c + 1) - 1_I4P
         select case (plan%vkind(v))
         case (0_I4P)
            ep%fl(:, v) = f0
         case (1_I4P)
            ep%fl(:, v) = f1
         case default
            s = plan%vs(v)
            do q = 1_I4P, 4_I4P
               if (known(q)) then
                  n = 0_I4P
                  do k = n0, n1
                     if (plan%cncor(k) == q .and. plan%cntime(k) < s) n = n + 1_I4P
                  enddo
                  ep%fl(q, v) = mod(cs0(q) + n, 2_I4P) == 0_I4P
               else
                  ep%fl(q, v) = map_corner_in(mc, curve, c0, s, plan_pose_2d(plan, plan%vpose(v)), cn(:, q))
               endif
            enddo
         endselect
      enddo
   enddo
   endsubroutine map_flags_2d

   subroutine map_evals_2d(mc, plan, x, y, ep)
   !< The evaluations of the 2-D map (issue #105 step 1d; `map_evals_t`, corner flags set by `map_flags_2d`): every
   !< item's pieces and features on its cell's local curve, placed at its instant (pose0, pose1 or s: `map_place_2d`).
   !< Sizes pass (`map_item_2d_count`: the walk's counts), scans of the workspaces and of the output capacities, the
   !< evaluation pass (`map_item_2d_eval`) in chunks of at most MAP_CHUNK_BYTES of workspace, then the scans of the
   !< pieces and features and the fill of the features. Bit for bit the per-cell `map_eval` of issue #102 (the kernels of
   !< `cut_cell_2d_pieces`).
   type(map_ctx_t),   intent(in)    :: mc         !< State.
   type(map_plan_t),  intent(in)    :: plan       !< Host phase.
   real(R8P),         intent(in)    :: x(:), y(:) !< Cell edges.
   type(map_evals_t), intent(inout) :: ep         !< Evaluations (fl set; the rest built).
   integer(I4P), allocatable        :: vc(:), nxs(:), nps(:), wr(:), wi(:), kcap(:), keyc(:), kpc(:) !< Work.
   real(R8P),    allocatable        :: rw(:)      !< Workspace of a chunk.
   integer(I4P), allocatable        :: iw(:)      !< Workspace of a chunk.
   integer(I4P), allocatable        :: ro(:), io(:) !< Workspace offsets in a chunk.
   real(R8P)                        :: pr(2, 2), pt(2) !< Pose of an item.
   integer(I4P)                     :: nx, n, c, i, j, v, a, b, l0, nl, mp, mk !< Work.
   logical                          :: tr         !< Translation path.

   nx = size(x, kind=I4P) - 1_I4P ; n = ep%n ; tr = mc%path == MAP_TRANSLATION
   allocate(vc(n), nxs(n), nps(n), wr(n), wi(n), kcap(n + 1), ro(n + 1), io(n + 1))
   do c = 1_I4P, plan%ncell
      do v = plan%vstart(c), plan%vstart(c + 1) - 1_I4P
         vc(v) = c ; wr(v) = 2_I4P * (plan%lstart(c + 1) - plan%lstart(c)) ; wi(v) = 0_I4P
      enddo
   enddo
   ! sizes pass: the placed local curve, the walk's counts
   a = 1_I4P
   do while (a <= n)
      call map_chunk(a, n, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, c, i, j, l0, nl, pr, pt) schedule(dynamic, 16)
      do v = a, b
         c = vc(v) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call plan_pose_rt(plan, plan%vpose(v), 2_I4P, pr, pt)
         call map_item_2d_count(nl, plan%lel(l0), plan%lnext(l0), mc%nv, mc%u0, mc%v0, mc%u1, mc%v1, mc%bu, mc%bv, &
                                plan%vkind(v), tr, plan%vs(v), mc%d, pr, pt, x(i), x(i + 1), y(j), y(j + 1), wr(v), &
                                rw(ro(v)), nxs(v), nps(v))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   ! workspaces and capacities: features, piece areas (the end evaluations' only)
   allocate(ep%mstart(n + 1))
   do v = 1_I4P, n
      c = vc(v) ; nl = plan%lstart(c + 1) - plan%lstart(c)
      call map_item_2d_sizes(nl, nxs(v), nps(v), wr(v), wi(v), mp, mk)
      kcap(v) = mk ; ep%mstart(v) = merge(0_I4P, mp, plan%vkind(v) == 2_I4P)
   enddo
   call csr_scan(kcap) ; call csr_scan(ep%mstart)
   allocate(ep%closed(n), ep%np(n), ep%nk(n), keyc(kcap(n + 1) - 1_I4P), kpc(kcap(n + 1) - 1_I4P), &
            ep%pm(ep%mstart(n + 1))) ! + 1: every slot start inside
   ! evaluation pass: one workspace per item
   a = 1_I4P
   do while (a <= n)
      call map_chunk(a, n, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, c, i, j, l0, nl, pr, pt) schedule(dynamic, 16)
      do v = a, b
         c = vc(v) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call plan_pose_rt(plan, plan%vpose(v), 2_I4P, pr, pt)
         call map_item_2d_eval(nl, plan%lel(l0), plan%lnext(l0), mc%nv, mc%u0, mc%v0, mc%u1, mc%v1, mc%bu, mc%bv, &
                               plan%vkind(v), tr, plan%vs(v), mc%d, pr, pt, x(i), x(i + 1), y(j), y(j + 1), ep%fl(:, v), &
                               nxs(v), nps(v), wr(v), rw(ro(v)), wi(v), iw(io(v)), ep%closed(v), ep%np(v), ep%nk(v), &
                               kcap(v + 1) - kcap(v), keyc(kcap(v)), kpc(kcap(v)), ep%mstart(v + 1) - ep%mstart(v), &
                               ep%pm(ep%mstart(v)))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   call map_evals_fill(ep, kcap, keyc, kpc)
   endsubroutine map_evals_2d

   subroutine map_evals_fill(ep, kcap, keyc, kpc)
   !< The CSR of the evaluations' pieces (nodes of the tracks) and features from their counts (scans), then the features
   !< copied from their capacity slots (issue #105 step 1d).
   type(map_evals_t), intent(inout) :: ep              !< Evaluations (closed, np, nk set).
   integer(I4P),      intent(in)    :: kcap(:)         !< Capacity slots of the features [n + 1].
   integer(I4P),      intent(in)    :: keyc(:), kpc(:) !< Features and pieces in their slots.
   integer(I4P)                     :: v, k            !< Work.

   allocate(ep%pstart(ep%n + 1), ep%kstart(ep%n + 1))
   ep%pstart(1:ep%n) = ep%np ; ep%kstart(1:ep%n) = ep%nk
   call csr_scan(ep%pstart) ; call csr_scan(ep%kstart)
   allocate(ep%key(ep%kstart(ep%n + 1) - 1_I4P), ep%kpiece(ep%kstart(ep%n + 1) - 1_I4P))
   do v = 1_I4P, ep%n
      do k = 0_I4P, ep%nk(v) - 1_I4P
         ep%key(ep%kstart(v) + k) = keyc(kcap(v) + k) ; ep%kpiece(ep%kstart(v) + k) = kpc(kcap(v) + k)
      enddo
   enddo
   endsubroutine map_evals_fill

   subroutine map_chunk(a, n, wr, wi, b, ro, io, rw, iw)
   !< The chunk of items a .. b whose workspaces (wr reals, wi integers each) add up to at most MAP_CHUNK_BYTES (one
   !< item at least), their offsets in the chunk's workspace (ro, io) and the workspace (grown if needed, issue #105
   !< step 1d).
   integer(I4P),              intent(in)    :: a, n           !< First item, items.
   integer(I4P),              intent(in)    :: wr(:), wi(:)   !< Workspace of every item [n].
   integer(I4P),              intent(out)   :: b              !< Last item of the chunk.
   integer(I4P),              intent(inout) :: ro(:), io(:)   !< Offsets of the chunk's items [n + 1].
   real(R8P),    allocatable, intent(inout) :: rw(:)          !< Real workspace.
   integer(I4P), allocatable, intent(inout) :: iw(:)          !< Integer workspace.
   integer(I8P)                             :: tot            !< Bytes.
   integer(I4P)                             :: v              !< Work.

   b = a ; tot = 8_I8P * wr(a) + 4_I8P * wi(a)
   do while (b < n)
      if (tot + 8_I8P * wr(b + 1) + 4_I8P * wi(b + 1) > MAP_CHUNK_BYTES) exit
      b = b + 1_I4P ; tot = tot + 8_I8P * wr(b) + 4_I8P * wi(b)
   enddo
   ro(a) = 1_I4P ; io(a) = 1_I4P
   do v = a, b
      ro(v + 1) = ro(v) + wr(v) ; io(v + 1) = io(v) + wi(v)
   enddo
   if (allocated(rw)) then
      if (size(rw, kind=I4P) < ro(b + 1)) deallocate(rw)
   endif
   if (.not. allocated(rw)) allocate(rw(max(1_I4P, ro(b + 1))))
   if (allocated(iw)) then
      if (size(iw, kind=I4P) < io(b + 1)) deallocate(iw)
   endif
   if (.not. allocated(iw)) allocate(iw(max(1_I4P, io(b + 1))))
   endsubroutine map_chunk

   pure subroutine plan_pose_rt(plan, k, dim, r, t)
   !< Rotation and translation of pose k of the plan's table (the identity for k = 0), as plain arrays.
   type(map_plan_t), intent(in)  :: plan          !< Plan.
   integer(I4P),     intent(in)  :: k, dim        !< Pose (0: none), 2 or 3.
   real(R8P),        intent(out) :: r(dim, dim), t(dim) !< Rotation, translation.
   integer(I4P)                  :: a             !< Work.

   if (k > 0_I4P) then
      r = plan%pose_r(:, :, k) ; t = plan%pose_t(:, k)
   else
      r = 0._R8P ; t = 0._R8P
      do a = 1_I4P, dim
         r(a, a) = 1._R8P
      enddo
   endif
   endsubroutine plan_pose_rt

   pure subroutine map_place_2d(nl, lg, nv, u0, v0, u1, v1, bu, bv, kind, tr, s, d, pr, pt, lu, lv)
   !< Kernel (issue #105 step 1d): the local vertices lg of an evaluation, placed: kind 0 the pose0 copy, 1 the pose1
   !< copy, 2 the instant s (the pose0 copy moved by s d on the translation path, else the body vertices by the pose
   !< (pr, pt) of the plan's table, R q + t as `to_world`). As the per-cell `map_eval` of issue #102, bit for bit.
   !$acc routine seq
   integer(I4P), intent(in)  :: nl                 !< Local vertices.
   integer(I4P), intent(in)  :: lg(nl)             !< Their curve indexes.
   integer(I4P), intent(in)  :: nv                 !< Curve vertices.
   real(R8P),    intent(in)  :: u0(nv), v0(nv)     !< Pose0 copy.
   real(R8P),    intent(in)  :: u1(nv), v1(nv)     !< Pose1 copy.
   real(R8P),    intent(in)  :: bu(nv), bv(nv)     !< Body frame.
   integer(I4P), intent(in)  :: kind               !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)  :: tr                 !< Translation path.
   real(R8P),    intent(in)  :: s, d(2)            !< Instant, translation over the step.
   real(R8P),    intent(in)  :: pr(2, 2), pt(2)    !< Pose at s (kind 2 off the translation path).
   real(R8P),    intent(out) :: lu(nl), lv(nl)     !< Placed vertices.
   real(R8P)                 :: w1, w2             !< Rotated vertex.
   integer(I4P)              :: k, g               !< Work.

   do k = 1_I4P, nl
      g = lg(k)
      select case (kind)
      case (0_I4P)
         lu(k) = u0(g) ; lv(k) = v0(g)
      case (1_I4P)
         lu(k) = u1(g) ; lv(k) = v1(g)
      case default
         if (tr) then
            lu(k) = u0(g) + s * d(1) ; lv(k) = v0(g) + s * d(2)
         else
            w1 = pr(1, 1) * bu(g) + pr(1, 2) * bv(g) ; w2 = pr(2, 1) * bu(g) + pr(2, 2) * bv(g)
            lu(k) = w1 + pt(1) ; lv(k) = w2 + pt(2)
         endif
      endselect
   enddo
   endsubroutine map_place_2d

   pure subroutine map_item_2d_sizes(nl, nxs, nps, nrw, niw, mp, mk)
   !< Workspace of a 2-D evaluation (`map_item_2d_eval`) with nl local vertices and the walk's counts nxs, nps (issue
   !< #105 step 1d): reals, the placed vertices (2 nl), the one-cell workspace (`cut_cell_2d_one_sizes`), the piece
   !< outputs (areas, face lengths, first moments: 9 mp); integers, the one-cell workspace and the wall pieces' segments
   !< and pieces (2 nps, + 1 so that their starts stay inside). mp = 4 + nxs + nps pieces at most, mk = 4 + nps
   !< features at most.
   !$acc routine seq
   integer(I4P), intent(in)  :: nl, nxs, nps     !< Local vertices, crossings, wall pieces.
   integer(I4P), intent(out) :: nrw, niw, mp, mk !< Reals, integers, piece and feature capacities.
   integer(I4P)              :: nr1, ni1         !< One-cell workspace.

   call cut_cell_2d_one_sizes(nxs, nps, nr1, ni1)
   mp = 4_I4P + nxs + nps ; mk = 4_I4P + nps
   nrw = 2_I4P * nl + nr1 + 9_I4P * mp ; niw = ni1 + 2_I4P * nps + 1_I4P
   endsubroutine map_item_2d_sizes

   pure subroutine map_item_2d_layout(nl, nxs, nps, nr1, ni1, mp, pa, pl, pq)
   !< Layout of a 2-D evaluation's real workspace (`map_item_2d_sizes`; issue #105 step 1e): the one-cell workspace
   !< (nr1 reals, ni1 integers) after the placed vertices (2 nl), then the piece outputs of capacity mp: areas from pa,
   !< face lengths [4, mp] from pl, first moments [4, mp] from pq.
   !$acc routine seq
   integer(I4P), intent(in)  :: nl, nxs, nps   !< Local vertices, crossings, wall pieces.
   integer(I4P), intent(out) :: nr1, ni1, mp   !< One-cell workspace, piece capacity.
   integer(I4P), intent(out) :: pa, pl, pq     !< Starts of the piece outputs.

   call cut_cell_2d_one_sizes(nxs, nps, nr1, ni1)
   mp = 4_I4P + nxs + nps
   pa = 2_I4P * nl + nr1 + 1_I4P ; pl = pa + mp ; pq = pl + 4_I4P * mp
   endsubroutine map_item_2d_layout

   pure subroutine map_item_2d_count(nl, lg, lnx, nv, u0, v0, u1, v1, bu, bv, kind, tr, s, d, pr, pt, x0, x1, y0, y1, nrw, &
                                     rw, nxs, nps)
   !< Kernel (issue #105 step 1d), sizes pass of a 2-D evaluation: the placed local curve (`map_place_2d`, in rw) and the
   !< count pass of its walk over the cell (`cut_cell_2d_walk`): nxs crossings, nps wall pieces.
   !$acc routine seq
   integer(I4P), intent(in)    :: nl                 !< Local vertices.
   integer(I4P), intent(in)    :: lg(nl), lnx(nl)    !< Their curve indexes, next local vertex.
   integer(I4P), intent(in)    :: nv                 !< Curve vertices.
   real(R8P),    intent(in)    :: u0(nv), v0(nv), u1(nv), v1(nv), bu(nv), bv(nv) !< Pose0 copy, pose1 copy, body.
   integer(I4P), intent(in)    :: kind               !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)    :: tr                 !< Translation path.
   real(R8P),    intent(in)    :: s, d(2), pr(2, 2), pt(2) !< Instant, translation, pose at s.
   real(R8P),    intent(in)    :: x0, x1, y0, y1     !< Cell.
   integer(I4P), intent(in)    :: nrw                !< Workspace (>= 2 nl).
   real(R8P),    intent(inout) :: rw(nrw)            !< Workspace.
   integer(I4P), intent(out)   :: nxs, nps           !< Crossings, wall pieces.
   integer(I4P)                :: idum(2, 1)         !< Work.
   real(R8P)                   :: rdum(4, 1)         !< Work.

   call map_place_2d(nl, lg, nv, u0, v0, u1, v1, bu, bv, kind, tr, s, d, pr, pt, rw(1), rw(nl + 1_I4P))
   call cut_cell_2d_walk(nl, rw(1), rw(nl + 1_I4P), lnx, x0, x1, y0, y1, .false., 0_I4P, 0_I4P, nxs, nps, idum, rdum, &
                         idum, rdum, idum, rdum, idum)
   endsubroutine map_item_2d_count

   subroutine map_item_2d_eval(nl, lg, lnx, nv, u0, v0, u1, v1, bu, bv, kind, tr, s, d, pr, pt, x0, x1, y0, y1, fl, nxs, &
                               nps, nrw, rw, niw, iw, closed, np, nk, mk, key, kpiece, mm, parea)
   !< Kernel (issue #105 step 1d): one 2-D evaluation on its workspace (`map_item_2d_sizes`): the placed local curve,
   !< the cell's pieces with the corner flags fl (`cut_cell_2d_one_pieces`), then its features, increasing: the fluid
   !< corners (keys -4..-1) then the wall segments (their curve indexes), each with its piece; the piece areas into
   !< parea (capacity mm, the end evaluations'). The per-cell `map_eval` of issue #102 bit for bit. Not closed: np = nk = 0.
   !$acc routine seq
   integer(I4P), intent(in)    :: nl                 !< Local vertices.
   integer(I4P), intent(in)    :: lg(nl), lnx(nl)    !< Their curve indexes, next local vertex.
   integer(I4P), intent(in)    :: nv                 !< Curve vertices.
   real(R8P),    intent(in)    :: u0(nv), v0(nv), u1(nv), v1(nv), bu(nv), bv(nv) !< Pose0 copy, pose1 copy, body.
   integer(I4P), intent(in)    :: kind               !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)    :: tr                 !< Translation path.
   real(R8P),    intent(in)    :: s, d(2), pr(2, 2), pt(2) !< Instant, translation, pose at s.
   real(R8P),    intent(in)    :: x0, x1, y0, y1     !< Cell.
   logical,      intent(in)    :: fl(4)              !< Corner flags.
   integer(I4P), intent(in)    :: nxs, nps           !< The walk's counts (`map_item_2d_count`).
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes (`map_item_2d_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   logical,      intent(out)   :: closed             !< The cell boundary closed.
   integer(I4P), intent(out)   :: np, nk             !< Pieces, features.
   integer(I4P), intent(in)    :: mk                 !< Feature capacity.
   integer(I4P), intent(inout) :: key(mk), kpiece(mk) !< Features, their pieces.
   integer(I4P), intent(in)    :: mm                 !< Piece area capacity (0: none kept).
   real(R8P),    intent(inout) :: parea(mm)          !< Piece areas.
   integer(I4P)                :: nr1, ni1, mp, pa, pl, pq, k, cp(4) !< Work.

   call map_item_2d_layout(nl, nxs, nps, nr1, ni1, mp, pa, pl, pq)
   call map_place_2d(nl, lg, nv, u0, v0, u1, v1, bu, bv, kind, tr, s, d, pr, pt, rw(1), rw(nl + 1_I4P))
   call cut_cell_2d_one_pieces(nl, rw(1), rw(nl + 1_I4P), lnx, x0, x1, y0, y1, fl(1), fl(2), fl(3), fl(4), nxs, nps, nr1, &
                               rw(2_I4P * nl + 1_I4P), ni1, iw(1), closed, np, mp, rw(pa), iw(ni1 + 1_I4P), &
                               iw(ni1 + nps + 1_I4P), cp, rw(pl), rw(pq))
   nk = 0_I4P
   if (.not. closed) then
      np = 0_I4P ; return
   endif
   do k = 4_I4P, 1_I4P, -1_I4P
      if (cp(k) == 0_I4P) cycle
      nk = nk + 1_I4P ; key(nk) = -k ; kpiece(nk) = cp(k)
   enddo
   do k = 1_I4P, nps ! walk order: increasing segments, one wall piece each
      if (iw(ni1 + nps + k) == 0_I4P) cycle
      nk = nk + 1_I4P ; key(nk) = lg(iw(ni1 + k)) ; kpiece(nk) = iw(ni1 + nps + k)
   enddo
   do k = 1_I4P, min(np, mm)
      parea(k) = rw(pa + k - 1_I4P)
   enddo
   endsubroutine map_item_2d_eval

   integer(I4P) function map_end_check(ep, v0, ne, np0, np1, ref0, ref1) result(kind)
   !< A cell's evaluations v0 .. v0 + ne - 1 checked (issue #105 step 1d): every one closed, the end evaluations the
   !< end-pose cuts' pieces, bit for bit (counts, then the measures they keep: 2-D the areas,
   !< 3-D the volumes of several pieces). 1 passed, -1 untrusted.
   type(map_evals_t), intent(in) :: ep         !< Evaluations.
   integer(I4P),      intent(in) :: v0, ne     !< First evaluation, evaluations.
   integer(I4P),      intent(in) :: np0, np1   !< Pieces of the end cuts.
   real(R8P),         intent(in) :: ref0(:), ref1(:) !< Their measures.
   integer(I4P)                  :: e, k, v1   !< Work.

   kind = 1_I4P ; v1 = v0 + ne - 1_I4P
   do e = v0, v1
      if (.not. ep%closed(e)) kind = -1_I4P
   enddo
   if (ep%np(v0) /= np0 .or. ep%np(v1) /= np1) kind = -1_I4P
   if (kind == -1_I4P) return
   if (ep%mstart(v0 + 1) > ep%mstart(v0)) then
      do k = 1_I4P, np0
         if (ep%pm(ep%mstart(v0) + k - 1_I4P) /= ref0(k)) kind = -1_I4P
      enddo
   endif
   if (ep%mstart(v1 + 1) > ep%mstart(v1)) then
      do k = 1_I4P, np1
         if (ep%pm(ep%mstart(v1) + k - 1_I4P) /= ref1(k)) kind = -1_I4P
      enddo
   endif
   endfunction map_end_check

   integer(I4P) function map_tracks_work(plan, ep) result(nw)
   !< The workspace of `map_tracks` for the largest evaluated cell.
   type(map_plan_t),  intent(in) :: plan !< Plan.
   type(map_evals_t), intent(in) :: ep   !< Evaluations.
   integer(I4P)                  :: c, ne !< Work.

   nw = 1_I4P
   do c = 1_I4P, plan%ncell
      ne = plan%vstart(c + 1) - plan%vstart(c)
      if (ne > 0_I4P) nw = max(nw, 2_I4P * ne + 3_I4P + ep%pstart(plan%vstart(c + 1)) - ep%pstart(plan%vstart(c)))
   enddo
   endfunction map_tracks_work

   subroutine map_tracks(ne, pst, kst, nk, key, kpiece, nn, w, comp, nt)
   !< Kernel (issue #105 step 1d): the tracks of one cell (`track_components`) from its slices of the flat evaluations:
   !< pst, kst the CSR starts of its ne evaluations' pieces and features (global), key, kpiece its features, comp(nn)
   !< the track of each of its pieces (nodes).
   !$acc routine seq
   integer(I4P), intent(in)    :: ne                 !< Evaluations.
   integer(I4P), intent(in)    :: pst(ne + 1)        !< First piece of each evaluation.
   integer(I4P), intent(in)    :: kst(ne + 1)        !< First feature of each evaluation.
   integer(I4P), intent(in)    :: nk                 !< Features.
   integer(I4P), intent(in)    :: key(nk), kpiece(nk) !< Features, their pieces.
   integer(I4P), intent(in)    :: nn                 !< Pieces (nodes).
   integer(I4P), intent(inout) :: w(2 * ne + 3 + nn) !< Work: local starts, union-find.
   integer(I4P), intent(inout) :: comp(nn)           !< Track of each node.
   integer(I4P), intent(out)   :: nt                 !< Tracks.
   integer(I4P)                :: e                  !< Work.

   do e = 1_I4P, ne + 1_I4P
      w(e) = pst(e) - pst(1) + 1_I4P ; w(ne + 1_I4P + e) = kst(e) - kst(1) + 1_I4P
   enddo
   call track_components(ne, w(1), w(ne + 2_I4P), nk, key, kpiece, w(2_I4P * ne + 3_I4P), comp, nt)
   endsubroutine map_tracks

   pure subroutine map_local(nxt, lseg, lg, lnx)
   !< The local curve of a cell: the end vertices of its segments in increasing order (lg, as the whole curve's indexes)
   !< and the next vertex of each in that numbering (0 if its segment is not local), so that `walk_all` walks the local
   !< segments in the whole curve's order. A per-cell deduplication (sort, then unique; issue #105 step 1c).
   integer(I4P),              intent(in)  :: nxt(:)  !< Next vertex of every curve vertex.
   integer(I4P),              intent(in)  :: lseg(:) !< Local segments.
   integer(I4P), allocatable, intent(out) :: lg(:)   !< Their end vertices, increasing, distinct.
   integer(I4P), allocatable, intent(out) :: lnx(:)  !< Next local vertex.
   integer(I4P), allocatable              :: tmp(:)  !< Work.
   integer(I4P)                           :: k, n, a, h, nl !< Work.

   n = 2_I4P * size(lseg, kind=I4P)
   allocate(tmp(n))
   do k = 1_I4P, size(lseg, kind=I4P)
      tmp(2_I4P * k - 1_I4P) = lseg(k) ; tmp(2_I4P * k) = nxt(lseg(k))
   enddo
   do k = 2_I4P, n ! insertion sort (a few vertices)
      a = tmp(k) ; h = k - 1_I4P
      do while (h >= 1_I4P)
         if (tmp(h) <= a) exit
         tmp(h + 1_I4P) = tmp(h) ; h = h - 1_I4P
      enddo
      tmp(h + 1_I4P) = a
   enddo
   nl = min(n, 1_I4P)
   do k = 2_I4P, n
      if (tmp(k) /= tmp(nl)) then
         nl = nl + 1_I4P ; tmp(nl) = tmp(k)
      endif
   enddo
   lg = tmp(1:nl)
   allocate(lnx(nl)) ; lnx = 0_I4P
   do k = 1_I4P, size(lseg, kind=I4P)
      lnx(local_index(lg, lseg(k))) = local_index(lg, nxt(lseg(k)))
   enddo
   endsubroutine map_local

   pure integer(I4P) function local_index(lg, g) result(k)
   !< The position of g in the increasing list lg (bisection; g is in lg).
   integer(I4P), intent(in) :: lg(:) !< Increasing, distinct.
   integer(I4P), intent(in) :: g     !< Value.
   integer(I4P)             :: lo, hi, mid !< Work.

   lo = 1_I4P ; hi = size(lg, kind=I4P) ; k = 0_I4P
   do while (lo <= hi)
      mid = (lo + hi) / 2_I4P
      if (lg(mid) == g) then
         k = mid ; return
      elseif (lg(mid) < g) then
         lo = mid + 1_I4P
      else
         hi = mid - 1_I4P
      endif
   enddo
   endfunction local_index

   subroutine map_events(mc, x0, x1, y0, y1, nl, lg, lnx, b, bw, nb, cert, unc, ct, cc)
   !< The events of the cell [x0, x1] x [y0, y1] in (0, 1), sorted, each once: a local vertex crossing one of its four
   !< face lines inside the face, a corner crossing a local segment. Without rotation, the formulas of
   !< `faces_translation_2d` on the pose0 copy (`collect`, `node_breakpoints`); with a rotation, `line_events_2d` and
   !< `node_events_2d` on the local curve (`isolate_zeros`; cert false if one is not certified; unc true if an event of
   !< the cell is a contact within the round-off of its function, not certified by a sign change). bw (issue #104):
   !< each event's half-width, the instants it may stand for in double precision (`EVENT_RES` on the explicit formulas,
   !< `isolate_zeros`' width on the screw ones; an exact duplicate keeps the larger). ct, cc (issue #105 step 1d): the
   !< corner events, every one (no deduplication), instants and corners (1..4: (x0,y0), (x1,y0), (x1,y1), (x0,y1), as reals), for
   !< `map_corner_known`.
   type(map_ctx_t),        intent(in)  :: mc             !< State.
   real(R8P),              intent(in)  :: x0, x1, y0, y1 !< Cell.
   integer(I4P),           intent(in)  :: nl             !< Local vertices.
   integer(I4P),           intent(in)  :: lg(:), lnx(:)  !< Local vertices, next.
   real(R8P), allocatable, intent(out) :: b(:)           !< Events.
   real(R8P), allocatable, intent(out) :: bw(:)          !< Their half-widths (issue #104).
   integer(I4P),           intent(out) :: nb             !< Their number.
   logical,                intent(out) :: cert           !< Every event certified.
   logical,                intent(out) :: unc            !< An event is an uncertified contact (tangency).
   real(R8P), allocatable, intent(out) :: ct(:)          !< Corner events: instants.
   real(R8P), allocatable, intent(out) :: cc(:)          !< Their corners (1..4, as reals).
   type(rotation_2d_t)                 :: lc             !< Local rotation state.
   type(events_2d_t)                   :: ev             !< Events of a line or a node.
   real(R8P), allocatable              :: tmp(:), tw(:), bn(:), bnw(:) !< Work.
   real(R8P)                           :: cn(2, 4), lin(4), lo(4), hi(4), sk, w !< Corners, face lines, spans.
   integer(I4P)                        :: k, f, g, q, nax, ncn, nc !< Work.

   cert = .true. ; unc = .false. ; nb = 0_I4P
   allocate(ct(4), cc(4)) ; ncn = 0_I4P
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
            nc = ncn ; call push_value(ct, nc, bn(q)) ; call push_value(cc, ncn, real(f, R8P))
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
            nc = ncn ; call push_value(ct, nc, ev%z(q)) ; call push_value(cc, ncn, real(f, R8P))
         enddo
      enddo
   endif
   call sort_unique_w(tmp, tw, nb)
   b = tmp(1:nb) ; bw = tw(1:nb)
   tmp = ct(1:ncn) ; call move_alloc(tmp, ct)
   tw = cc(1:ncn) ; call move_alloc(tw, cc)
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

   function map_corner_in(mc, curve, c0, s, ps, p) result(inside)
   !< The corner p inside the solid at instant s: the parity of the curve's crossings left of p on its horizontal line
   !< (`line_crossings`, exact with the symbolic tie-break), p taken to the frame of the curve that the path moves:
   !< the pose0 copy shifted back by s d on the translation path, the body frame (`to_body` of the screw pose ps at s,
   !< `pose_interpolate`) else.
   type(map_ctx_t),       intent(in) :: mc        !< State.
   type(curve_object),    intent(in) :: curve, c0 !< Curve (body frame), copy at pose0.
   real(R8P),             intent(in) :: s, p(2)   !< Instant, corner.
   type(rigid_pose_2d_t), intent(in) :: ps        !< Screw pose at s (off the translation path).
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
      q = to_body(ps, p)
      call curve%line_crossings(1_I4P, q(2), min(q(1), mc%bmin(1)) - 1._R8P, q(1), t, g, e, n)
   endif
   inside = mod(n, 2_I4P) == 1_I4P
   endfunction map_corner_in

   subroutine map_pairs_build(plan, ntr, moff, rot, pr, ncol)
   !< The (cell, interval) pairs of the measures (issue #105 step 1e; `map_pairs_t`): every interval of positive width of
   !< every measured cell (moff > 0, at least one track), in cell order then interval order, all active. Columns of the
   !< sums: with a rotation each pair has its own ntr columns (ncol in all), else the pairs of a cell add into the cell's
   !< columns of the measures (moff, ncol = 0).
   type(map_plan_t),  intent(in)  :: plan        !< Host phase.
   integer(I4P),      intent(in)  :: ntr(:)      !< Tracks of every cell.
   integer(I4P),      intent(in)  :: moff(:)     !< First column of every cell's measures (0: not measured).
   logical,           intent(in)  :: rot         !< The step rotates (doubling rounds).
   type(map_pairs_t), intent(out) :: pr          !< Pairs.
   integer(I4P),      intent(out) :: ncol        !< Columns of the pairs' own sums.
   real(R8P)                      :: sa, sb      !< Interval.
   integer(I4P)                   :: c, k, m, nb, b0, i0, n, r, a !< Work.

   pr%nk = 0_I4P ; n = 0_I4P
   do c = 1_I4P, plan%ncell
      if (moff(c) == 0_I4P .or. ntr(c) == 0_I4P) cycle
      pr%nk = pr%nk + 1_I4P ; n = n + plan%istart(c + 1) - plan%istart(c)
   enddo
   allocate(pr%cell(pr%nk), pr%kstart(pr%nk + 1), pr%k(n), pr%sa(n), pr%sb(n), pr%vr(n), pr%va(n), pr%narrow(n), &
            pr%col(n), pr%state(n), pr%skip(n), pr%inexact(n))
   k = 0_I4P ; n = 0_I4P ; ncol = 0_I4P
   do c = 1_I4P, plan%ncell
      if (moff(c) == 0_I4P .or. ntr(c) == 0_I4P) cycle
      k = k + 1_I4P ; pr%cell(k) = c ; pr%kstart(k) = n + 1_I4P
      b0 = plan%bstart(c) ; nb = plan%bstart(c + 1) - b0 ; i0 = plan%istart(c)
      do m = 1_I4P, nb + 1_I4P
         call map_interval(plan%bs(b0:b0 + nb - 1_I4P), nb, m, sa, sb)
         if (sb - sa <= 0._R8P) cycle
         n = n + 1_I4P ; pr%k(n) = k ; pr%sa(n) = sa ; pr%sb(n) = sb
         r = plan%eref(i0 + m - 1_I4P) ; a = plan%ealt(i0 + m - 1_I4P)
         pr%vr(n) = plan%vstart(c) + r - 1_I4P ; pr%va(n) = plan%vstart(c) + a - 1_I4P
         pr%narrow(n) = r /= a .or. plan%small(i0 + m - 1_I4P)
         if (rot) then
            pr%col(n) = ncol + 1_I4P ; ncol = ncol + ntr(c)
         else
            pr%col(n) = moff(c)
         endif
      enddo
   enddo
   pr%kstart(pr%nk + 1) = n + 1_I4P ; pr%n = n
   pr%state = PAIR_ACTIVE ; pr%skip = 0_I4P ; pr%inexact = .false.
   endsubroutine map_pairs_build

   subroutine map_round_nodes(pr, rot, ng, gq_s, gq_w, nn, nj, ns, nw)
   !< The nodes of one round of the measures (issue #105 step 1e), active pair after active pair, node after node:
   !< without rotation the two-point Gauss nodes of the pair's interval (the only round), with one the ng Gauss nodes of
   !< the round (`map_plan_t` sets). As the per-cell `map_measures` of issue #102 placed them, bit for bit.
   type(map_pairs_t),         intent(in)  :: pr          !< Pairs.
   logical,                   intent(in)  :: rot         !< The step rotates.
   integer(I4P),              intent(in)  :: ng          !< Nodes per pair of the round (rotation).
   real(R8P),                 intent(in)  :: gq_s(:), gq_w(:) !< Gauss node sets.
   integer(I4P),              intent(out) :: nn          !< Nodes.
   integer(I4P), allocatable, intent(out) :: nj(:)       !< Pair of each node.
   real(R8P),    allocatable, intent(out) :: ns(:), nw(:) !< Instants, weights.
   real(R8P), parameter :: G = 0.5_R8P / sqrt(3._R8P)  !< Two-point Gauss offset on a unit interval.
   real(R8P)            :: sa, sb, h                     !< Interval.
   integer(I4P)         :: jp, q, npp                    !< Work.

   npp = merge(ng, 2_I4P, rot)
   nn = npp * count(pr%state == PAIR_ACTIVE, kind=I4P)
   allocate(nj(nn), ns(nn), nw(nn))
   nn = 0_I4P
   do jp = 1_I4P, pr%n
      if (pr%state(jp) /= PAIR_ACTIVE) cycle
      sa = pr%sa(jp) ; sb = pr%sb(jp) ; h = sb - sa
      if (.not. rot) then
         nn = nn + 1_I4P ; nj(nn) = jp ; ns(nn) = 0.5_R8P * (sa + sb) - G * h ; nw(nn) = 0.5_R8P * h
         nn = nn + 1_I4P ; nj(nn) = jp ; ns(nn) = 0.5_R8P * (sa + sb) + G * h ; nw(nn) = 0.5_R8P * h
      else
         do q = ng - 1_I4P, 2_I4P * ng - 2_I4P
            nn = nn + 1_I4P ; nj(nn) = jp ; ns(nn) = sa + h * gq_s(q) ; nw(nn) = h * gq_w(q)
         enddo
      endif
   enddo
   endsubroutine map_round_nodes

   subroutine map_round_zero(pr, ntr, cur)
   !< The sums of the active pairs zeroed before a doubling round (issue #105 step 1e).
   type(map_pairs_t), intent(in)    :: pr       !< Pairs.
   integer(I4P),      intent(in)    :: ntr(:)   !< Tracks of every cell.
   real(R8P),         intent(inout) :: cur(:,:) !< Sums of the pairs.
   integer(I4P)                     :: jp, nt   !< Work.

   do jp = 1_I4P, pr%n
      if (pr%state(jp) /= PAIR_ACTIVE) cycle
      nt = ntr(pr%cell(pr%k(jp)))
      cur(:, pr%col(jp):pr%col(jp) + nt - 1_I4P) = 0._R8P
   enddo
   endsubroutine map_round_zero

   subroutine map_node_outcome(pr, jp, fail, nskip, ok)
   !< The outcome of one node for its pair (issue #105 step 1e), nodes taken in order: a node of a pair that is no
   !< longer active is ignored; a failed node fails its pair (the per-cell `map_node` of issue #102 stopped the cell
   !< there); else its skips are counted and ok says its pieces are to be added.
   type(map_pairs_t), intent(inout) :: pr     !< Pairs.
   integer(I4P),      intent(in)    :: jp     !< Pair.
   logical,           intent(in)    :: fail   !< The node failed.
   integer(I4P),      intent(in)    :: nskip  !< Its skips.
   logical,           intent(out)   :: ok     !< Add its pieces.

   ok = .false.
   if (pr%state(jp) /= PAIR_ACTIVE) return
   if (fail) then
      pr%state(jp) = PAIR_FAILED ; return
   endif
   pr%skip(jp) = pr%skip(jp) + nskip ; ok = .true.
   endsubroutine map_node_outcome

   pure subroutine map_node_reduce(nv, np, ptk, pv, nc, acc)
   !< Kernel (issue #105 step 1e): a node's values added to the sums of their tracks, in piece order (the summation order
   !< of the per-cell `map_node` of issue #102): acc(:, ptk(p)) += pv(:, p) for every piece with a track.
   !$acc routine seq
   integer(I4P), intent(in)    :: nv          !< Values per track.
   integer(I4P), intent(in)    :: np          !< Pieces.
   integer(I4P), intent(in)    :: ptk(np)     !< Track of each piece (0: none).
   real(R8P),    intent(in)    :: pv(nv, np)  !< Values of each piece.
   integer(I4P), intent(in)    :: nc          !< Tracks.
   real(R8P),    intent(inout) :: acc(nv, nc) !< Sums.
   integer(I4P)                :: p, f, t     !< Work.

   do p = 1_I4P, np
      t = ptk(p)
      if (t <= 0_I4P) cycle
      do f = 1_I4P, nv
         acc(f, t) = acc(f, t) + pv(f, p)
      enddo
   enddo
   endsubroutine map_node_reduce

   pure subroutine map_round_test(nf, nt, a1, a2, h, scl, sco, sna, sno, conv, inexact)
   !< Kernel (issue #105 step 1e): the convergence test of a doubling round (the criteria of `integrate_face_2d` /
   !< `one_plane`): every track's nf face measures (a1, a2: the rounds n / 2 and n) agree within 64 eps max(h, FLOOR_H)
   !< times the face scale scl, the nf outflows within that times sco; inexact (meaningful at 64 nodes): a change above
   !< 64 eps of the round-off scale sna (measures), sno (outflows).
   !$acc routine seq
   integer(I4P), intent(in)  :: nf, nt                 !< Faces, tracks.
   real(R8P),    intent(in)  :: a1(2 * nf, nt), a2(2 * nf, nt) !< Sums of the two rounds.
   real(R8P),    intent(in)  :: h                      !< Interval width.
   real(R8P),    intent(in)  :: scl(nf), sco(nf)       !< Scales of the measures and outflows.
   real(R8P),    intent(in)  :: sna, sno               !< Round-off scales.
   logical,      intent(out) :: conv, inexact          !< Converged; not converged to round-off.
   real(R8P)                 :: tol                    !< Tolerance.
   integer(I4P)              :: t, f                   !< Work.

   tol = 64._R8P * epsilon(1._R8P) * max(h, FLOOR_H)
   conv = .true. ; inexact = .false.
   do t = 1_I4P, nt
      do f = 1_I4P, nf
         if (abs(a2(f, t) - a1(f, t)) > tol * scl(f) .or. abs(a2(nf + f, t) - a1(nf + f, t)) > tol * sco(f)) &
            conv = .false.
      enddo
   enddo
   if (conv) return
   do t = 1_I4P, nt
      do f = 1_I4P, nf
         if (max(abs(a2(f, t) - a1(f, t)) / sna, abs(a2(nf + f, t) - a1(nf + f, t)) / sno) > 64._R8P * epsilon(1._R8P)) &
            inexact = .true.
      enddo
   enddo
   endsubroutine map_round_test

   subroutine map_round_close(pr, ntr, round, ng, nf, cur, prv, ksc, kso, ksa, ksn)
   !< The end of a doubling round (issue #105 step 1e), per active pair: the first round's sums kept, a later one's tested
   !< against them (`map_round_test`): converged, or at 64 nodes (inexact recorded), the pair is done, else its sums are
   !< kept for the next round. Then the active pairs of a cell after its first failed pair are dropped (the per-cell
   !< `map_node` of issue #102 stopped the cell there).
   type(map_pairs_t), intent(inout) :: pr                 !< Pairs.
   integer(I4P),      intent(in)    :: ntr(:)             !< Tracks of every cell.
   integer(I4P),      intent(in)    :: round, ng, nf      !< Round, its nodes per pair, faces.
   real(R8P),         intent(in)    :: cur(:,:)           !< Sums of the round.
   real(R8P),         intent(inout) :: prv(:,:)           !< Sums of the round before.
   real(R8P),         intent(in)    :: ksc(:,:), kso(:,:) !< Face scales per measured cell [nf, nk].
   real(R8P),         intent(in)    :: ksa(:), ksn(:)     !< Round-off scales per measured cell [nk].
   integer(I4P)                     :: jp, k, nt, c0      !< Work.
   logical                          :: conv, inexact, failed !< Work.

   do jp = 1_I4P, pr%n
      if (pr%state(jp) /= PAIR_ACTIVE) cycle
      k = pr%k(jp) ; nt = ntr(pr%cell(k)) ; c0 = pr%col(jp)
      if (round > 1_I4P) then
         call map_round_test(nf, nt, prv(:, c0:c0 + nt - 1_I4P), cur(:, c0:c0 + nt - 1_I4P), pr%sb(jp) - pr%sa(jp), &
                             ksc(:, k), kso(:, k), ksa(k), ksn(k), conv, inexact)
         if (conv) then
            pr%state(jp) = PAIR_DONE ; cycle
         elseif (ng >= 64_I4P) then ! the round-off of one evaluation (see FLOOR_H)
            pr%state(jp) = PAIR_DONE ; pr%inexact(jp) = inexact ; cycle
         endif
      endif
      prv(:, c0:c0 + nt - 1_I4P) = cur(:, c0:c0 + nt - 1_I4P)
   enddo
   do k = 1_I4P, pr%nk
      failed = .false.
      do jp = pr%kstart(k), pr%kstart(k + 1) - 1_I4P
         if (failed .and. pr%state(jp) == PAIR_ACTIVE) pr%state(jp) = PAIR_DROPPED
         if (pr%state(jp) == PAIR_FAILED) failed = .true.
      enddo
   enddo
   endsubroutine map_round_close

   subroutine map_measures_finish(pr, rot, ntr, moff, cur, meas, cok, skipped, exact)
   !< The measures of every measured cell from its pairs (issue #105 step 1e), as the per-cell `map_measures` of issue
   !< #102 ended: a cell whose pair failed is not ok; the skips and the inexact flags count up to its first failed pair
   !< (the later pairs did not run); with a rotation the measures are the sums of its pairs, added in interval order to
   !< zero (without one the pairs added into them directly).
   type(map_pairs_t), intent(in)    :: pr        !< Pairs.
   logical,           intent(in)    :: rot       !< The step rotates.
   integer(I4P),      intent(in)    :: ntr(:), moff(:) !< Tracks, first measure column of every cell.
   real(R8P),         intent(in)    :: cur(:,:)  !< Sums of the pairs (rotation).
   real(R8P),         intent(inout) :: meas(:,:) !< Measures (zero on input with a rotation).
   logical,           intent(inout) :: cok(:)    !< Measured cells ok.
   integer(I4P),      intent(inout) :: skipped   !< Skipped nodes or pieces.
   logical,           intent(inout) :: exact     !< Every integration converged.
   integer(I4P)                     :: k, c, nt, jp, m0, c0 !< Work.
   logical                          :: ok        !< Work.

   do k = 1_I4P, pr%nk
      c = pr%cell(k) ; nt = ntr(c) ; m0 = moff(c) ; ok = .true.
      do jp = pr%kstart(k), pr%kstart(k + 1) - 1_I4P
         if (pr%state(jp) == PAIR_FAILED) then
            ok = .false. ; exit
         endif
         skipped = skipped + pr%skip(jp)
         if (pr%inexact(jp)) exact = .false.
      enddo
      cok(c) = ok
      if (.not. (ok .and. rot)) cycle
      do jp = pr%kstart(k), pr%kstart(k + 1) - 1_I4P
         c0 = pr%col(jp)
         meas(:, m0:m0 + nt - 1_I4P) = meas(:, m0:m0 + nt - 1_I4P) + cur(:, c0:c0 + nt - 1_I4P)
      enddo
   enddo
   endsubroutine map_measures_finish

   subroutine node_instants(n, s, pid, us, nu)
   !< The pose table of a round of nodes (issue #105 step 1e, as `plan_instants`): the distinct instants s, increasing
   !< (exact equality), and the index of each node's.
   integer(I4P),           intent(in)  :: n       !< Nodes.
   real(R8P),              intent(in)  :: s(:)    !< Their instants.
   integer(I4P),           intent(out) :: pid(:)  !< Pose of each node.
   real(R8P), allocatable, intent(out) :: us(:)   !< Distinct instants.
   integer(I4P),           intent(out) :: nu      !< Their number.
   real(R8P),    allocatable :: key(:,:)          !< Instants.
   integer(I4P), allocatable :: ord(:)            !< Sorted order.
   integer(I4P)              :: r                 !< Work.

   allocate(key(1, n)) ; key(1, :) = s(1:n)
   call sort_lex(key, ord)
   nu = 0_I4P
   do r = 1_I4P, n
      if (r == 1_I4P) then
         nu = 1_I4P
      elseif (key(1, ord(r)) /= key(1, ord(r - 1_I4P))) then
         nu = nu + 1_I4P
      endif
      pid(ord(r)) = nu
   enddo
   allocate(us(nu))
   do r = 1_I4P, n
      us(pid(ord(r))) = key(1, ord(r))
   enddo
   endsubroutine node_instants

   pure subroutine node_pose_rt(qr, qt, k, dim, r, t)
   !< Rotation and translation of pose k of a round's table (the identity for k = 0), as plain arrays.
   real(R8P),    intent(in)  :: qr(:,:,:), qt(:,:) !< Table.
   integer(I4P), intent(in)  :: k, dim             !< Pose (0: none), 2 or 3.
   real(R8P),    intent(out) :: r(dim, dim), t(dim) !< Rotation, translation.
   integer(I4P)              :: a                  !< Work.

   if (k > 0_I4P) then
      r = qr(:, :, k) ; t = qt(:, k)
   else
      r = 0._R8P ; t = 0._R8P
      do a = 1_I4P, dim
         r(a, a) = 1._R8P
      enddo
   endif
   endsubroutine node_pose_rt

   pure integer(I4P) function piece_track(nk, ekey, ekpiece, p, nkr, rkey, rkpiece, nrp, rcomp) result(t)
   !< The track of piece p of an evaluation (features ekey, ekpiece) through its features found in a reference
   !< evaluation (rkey, rkpiece, increasing keys; rcomp the tracks of its pieces): 0 if none is found, -1 if they point to
   !< two tracks. 2-D and 3-D (issue #105 step 1d: the reference from the flat evaluations; step 1e: explicit shapes).
   !$acc routine seq
   integer(I4P), intent(in) :: nk                    !< Features of the evaluation.
   integer(I4P), intent(in) :: ekey(nk), ekpiece(nk) !< Its features, their pieces.
   integer(I4P), intent(in) :: p                     !< Piece.
   integer(I4P), intent(in) :: nkr                   !< Features of the reference.
   integer(I4P), intent(in) :: rkey(nkr), rkpiece(nkr) !< Its features, their pieces.
   integer(I4P), intent(in) :: nrp                   !< Pieces of the reference.
   integer(I4P), intent(in) :: rcomp(nrp)            !< Their tracks.
   integer(I4P)             :: k, q, lo, hi, mid, tk !< Work.

   t = 0_I4P
   do k = 1_I4P, nk
      if (ekpiece(k) /= p) cycle
      q = 0_I4P ; lo = 1_I4P ; hi = nkr ! the feature in the reference (keys increasing)
      do while (lo <= hi)
         mid = (lo + hi) / 2_I4P
         if (rkey(mid) == ekey(k)) then
            q = rkpiece(mid) ; exit
         elseif (rkey(mid) < ekey(k)) then
            lo = mid + 1_I4P
         else
            hi = mid - 1_I4P
         endif
      enddo
      if (q == 0_I4P) cycle
      tk = rcomp(q)
      if (t == 0_I4P) then
         t = tk
      elseif (t /= tk) then
         t = -1_I4P ; return
      endif
   enddo
   endfunction piece_track

   pure subroutine map_piece_assign(nk, ekey, ekpiece, p, narrow, negligible, nkr, rkey, rkpiece, nrp, rcomp, nka, akey, &
                                    akpiece, npa, acomp, t, act)
   !< Kernel (issue #105 step 1e): piece p of a node (features ekey, ekpiece) assigned as the per-cell `map_node` /
   !< `map3_node` of issue #102 did: the track of its features in the interval's evaluation r (`piece_track`), in a
   !< narrow interval else in the one after it, a. act: 1 added to track t; 0 dropped (no known feature, negligible:
   !< born within round-off of an event); 2 skipped (no track, or two, in a narrow interval); 3 the node fails.
   !$acc routine seq
   integer(I4P), intent(in)  :: nk                    !< Features of the node.
   integer(I4P), intent(in)  :: ekey(nk), ekpiece(nk) !< Its features, their pieces.
   integer(I4P), intent(in)  :: p                     !< Piece.
   logical,      intent(in)  :: narrow                !< A narrow interval.
   logical,      intent(in)  :: negligible            !< The piece's face measures below 1e-9 of the cell's.
   integer(I4P), intent(in)  :: nkr                   !< Features of r.
   integer(I4P), intent(in)  :: rkey(nkr), rkpiece(nkr) !< Its features, their pieces.
   integer(I4P), intent(in)  :: nrp                   !< Pieces of r.
   integer(I4P), intent(in)  :: rcomp(nrp)            !< Their tracks.
   integer(I4P), intent(in)  :: nka                   !< Features of a.
   integer(I4P), intent(in)  :: akey(nka), akpiece(nka) !< Its features, their pieces.
   integer(I4P), intent(in)  :: npa                   !< Pieces of a.
   integer(I4P), intent(in)  :: acomp(npa)            !< Their tracks.
   integer(I4P), intent(out) :: t, act                !< Track, action.

   t = piece_track(nk, ekey, ekpiece, p, nkr, rkey, rkpiece, nrp, rcomp)
   if (t == 0_I4P .and. narrow) t = piece_track(nk, ekey, ekpiece, p, nka, akey, akpiece, npa, acomp)
   act = 1_I4P
   if (t <= 0_I4P) then ! no known feature (negligible, or a missed event) or two tracks
      if (t == 0_I4P .and. negligible) then
         act = 0_I4P
      elseif (narrow) then
         act = 2_I4P
      else
         act = 3_I4P
      endif
   endif
   endsubroutine map_piece_assign

   subroutine map_measures_2d(mc, plan, ep, comp, curve, c0, x, y, ntr, moff, meas, cok)
   !< Per-track time integrals over the step of the fluid lengths of the cell's four faces (meas(1:4), -u, +u, -v, +v)
   !< and of the rigid-field outflows through them (meas(5:8)), for every measured cell (moff), interval by interval
   !< between the cell's events: inside an interval the pieces keep their features, so each piece of an evaluation at a
   !< quadrature node belongs to the track of the interval's evaluation that holds its features (`map_piece_assign`).
   !< Translation: two-point Gauss per interval (exact: the lengths are linear in s). Rotation: Gauss doubled from 2
   !< nodes per interval until every track's values agree to round-off (the criteria of `integrate_face_2d`; mc%exact
   !< false if not at 64). The Gauss path integrates the same way, by its motion (two-point per interval without
   !< rotation, doubled with one): its nq nodes over the step would miss a track living between two events, or weight
   !< it heavily (issue #102 step 3: per-track residuals up to 0.39 of a cell). cok false if a piece at a node cannot be
   !< assigned.
   !< By rounds over flat items (issue #105 step 1e): the (cell, interval) pairs (`map_pairs_t`), per round the nodes
   !< of the active pairs (`map_round_nodes`: one round of two nodes per pair without rotation, rounds of 2, 4, ..., 64
   !< with one), each node one item (`map_nodes_2d`), the round's test per pair (`map_round_close`), then the measures
   !< (`map_measures_finish`). Summation order of the per-cell form: within a pair node after node, piece after piece,
   !< then pair after pair (without rotation one running sum per cell): bit for bit.
   type(map_ctx_t),    intent(inout) :: mc         !< State (skipped, exact).
   type(map_plan_t),   intent(in)    :: plan       !< Host phase.
   type(map_evals_t),  intent(in)    :: ep         !< Evaluations.
   integer(I4P),       intent(in)    :: comp(:)    !< Track of every piece of every evaluation.
   type(curve_object), intent(in)    :: curve, c0  !< Curve (body frame), copy at pose0.
   real(R8P),          intent(in)    :: x(:), y(:) !< Cell edges.
   integer(I4P),       intent(in)    :: ntr(:)     !< Tracks of every cell.
   integer(I4P),       intent(in)    :: moff(:)    !< First column of every measured cell's measures (0: not measured).
   real(R8P),          intent(inout) :: meas(:,:)  !< Measures [8, *] (zero on input).
   logical,            intent(out)   :: cok(:)     !< Measured cells ok [nx ny].
   type(map_pairs_t)         :: pr                 !< Pairs.
   real(R8P),    allocatable :: cur(:,:), prv(:,:) !< Sums of the pairs (rotation).
   real(R8P),    allocatable :: ksc(:,:), kso(:,:), ksa(:), ksn(:) !< Scales of the test per measured cell.
   real(R8P),    allocatable :: ns(:), nw(:)       !< Nodes of a round: instants, weights.
   integer(I4P), allocatable :: nj(:)              !< Nodes of a round: pairs.
   real(R8P)                 :: vel                !< Speed bound.
   integer(I4P)              :: nx, k, c, i, j, round, ng, nn, ncol !< Work.
   logical                   :: rot                !< A rotation.

   nx = size(x, kind=I4P) - 1_I4P ; rot = mc%dth /= 0._R8P ; cok = .true.
   call map_pairs_build(plan, ntr, moff, rot, pr, ncol)
   allocate(cur(8, ncol), prv(8, ncol), ksc(4, pr%nk), kso(4, pr%nk), ksa(pr%nk), ksn(pr%nk))
   if (rot) then
      vel = abs(mc%dth) * mc%xm + sqrt(mc%d(1)**2 + mc%d(2)**2)
      do k = 1_I4P, pr%nk
         c = pr%cell(k) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         ksc(:, k) = [y(j + 1) - y(j), y(j + 1) - y(j), x(i + 1) - x(i), x(i + 1) - x(i)] ; kso(:, k) = ksc(:, k) * vel
         ksa(k) = mc%xm ; ksn(k) = mc%xm * vel
      enddo
   endif
   do round = 1_I4P, merge(6_I4P, 1_I4P, rot)
      ng = 2_I4P**round
      call map_round_nodes(pr, rot, ng, plan%gq_s, plan%gq_w, nn, nj, ns, nw)
      if (nn == 0_I4P) exit
      if (rot) then
         call map_round_zero(pr, ntr, cur)
         call map_nodes_2d(mc, plan, ep, comp, curve, c0, x, y, ntr, pr, nn, nj, ns, nw, cur)
         call map_round_close(pr, ntr, round, ng, 4_I4P, cur, prv, ksc, kso, ksa, ksn)
      else
         call map_nodes_2d(mc, plan, ep, comp, curve, c0, x, y, ntr, pr, nn, nj, ns, nw, meas)
      endif
   enddo
   call map_measures_finish(pr, rot, ntr, moff, cur, meas, cok, mc%skipped, mc%exact)
   endsubroutine map_measures_2d

   subroutine map_nodes_2d(mc, plan, ep, comp, curve, c0, x, y, ntr, pr, nn, nj, ns, nw, acc)
   !< One round of 2-D nodes (issue #105 step 1e). Host: the poses at the round's distinct instants off the translation
   !< path (`pose_interpolate`, `node_instants`), the corner flags of every node (its interval's evaluation's, unless
   !< that evaluation and the one after it differ in their first piece: then the exact queries `map_corner_in` at s, as
   !< the per-cell `map_node` of issue #102). Then the sizes pass (`map_item_2d_count`) and, in chunks of at most
   !< MAP_CHUNK_BYTES of workspace, the nodes (`map_node_2d_eval`) and their sums into their pairs' columns of acc, node
   !< after node (`map_node_outcome`, `map_node_reduce`).
   type(map_ctx_t),    intent(in)    :: mc         !< State.
   type(map_plan_t),   intent(in)    :: plan       !< Host phase.
   type(map_evals_t),  intent(in)    :: ep         !< Evaluations.
   integer(I4P),       intent(in)    :: comp(:)    !< Track of every piece of every evaluation.
   type(curve_object), intent(in)    :: curve, c0  !< Curve (body frame), copy at pose0.
   real(R8P),          intent(in)    :: x(:), y(:) !< Cell edges.
   integer(I4P),       intent(in)    :: ntr(:)     !< Tracks of every cell.
   type(map_pairs_t),  intent(inout) :: pr         !< Pairs.
   integer(I4P),       intent(in)    :: nn         !< Nodes.
   integer(I4P),       intent(in)    :: nj(:)      !< Pair of each node.
   real(R8P),          intent(in)    :: ns(:), nw(:) !< Instants, weights.
   real(R8P),          intent(inout) :: acc(:,:)   !< Sums [8, *].
   logical,      allocatable :: fl(:,:), nf(:)     !< Corner flags, failed.
   integer(I4P), allocatable :: pid(:), nxs(:), nps(:), wr(:), wi(:), ro(:), io(:), np(:), nsk(:), iw(:) !< Work.
   real(R8P),    allocatable :: us(:), qr(:,:,:), qt(:,:), rw(:) !< Pose table, workspace.
   type(rigid_pose_2d_t)     :: ps                 !< Pose.
   real(R8P)                 :: cn(2, 4), prr(2, 2), ptt(2), tq(2) !< Work.
   integer(I4P)              :: nx, v, jp, c, i, j, l0, nl, vr, va, a, b, k, q, npose, nr, ni, mp, mk !< Work.
   logical                   :: tr, ok             !< Translation path, node added.

   nx = size(x, kind=I4P) - 1_I4P ; tr = mc%path == MAP_TRANSLATION
   allocate(fl(4, nn), nf(nn), pid(nn), nxs(nn), nps(nn), wr(nn), wi(nn), ro(nn + 1), io(nn + 1), np(nn), nsk(nn))
   pid = 0_I4P ; npose = 0_I4P
   if (.not. tr) then ! the screw poses at the round's instants
      call node_instants(nn, ns, pid, us, npose)
      allocate(qr(2, 2, npose), qt(2, npose))
      do k = 1_I4P, npose
         ps = pose_interpolate(mc%pose0, mc%pose1, us(k)) ; qr(:, :, k) = ps%r ; qt(:, k) = ps%t
      enddo
   else
      allocate(qr(2, 2, 0), qt(2, 0))
   endif
   do v = 1_I4P, nn ! corner flags (host: the exact queries)
      jp = nj(v) ; c = pr%cell(pr%k(jp)) ; vr = pr%vr(jp) ; va = pr%va(jp)
      if (ep%pstart(va) == ep%pstart(vr)) then ! a wide interval: no corner crosses a segment inside it
         fl(:, v) = ep%fl(:, vr)
      else
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         cn = reshape([x(i), y(j), x(i + 1), y(j), x(i + 1), y(j + 1), x(i), y(j + 1)], [2, 4])
         if (pid(v) > 0_I4P) then
            ps%r = qr(:, :, pid(v)) ; ps%t = qt(:, pid(v))
         endif
         do q = 1_I4P, 4_I4P
            fl(q, v) = map_corner_in(mc, curve, c0, ns(v), ps, cn(:, q))
         enddo
      endif
      wr(v) = 2_I4P * (plan%lstart(c + 1) - plan%lstart(c)) ; wi(v) = 0_I4P
   enddo
   ! sizes pass: the placed local curve, the walk's counts
   a = 1_I4P
   do while (a <= nn)
      call map_chunk(a, nn, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, c, i, j, l0, nl, prr, ptt) schedule(dynamic, 16)
      do v = a, b
         c = pr%cell(pr%k(nj(v))) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call node_pose_rt(qr, qt, pid(v), 2_I4P, prr, ptt)
         call map_item_2d_count(nl, plan%lel(l0), plan%lnext(l0), mc%nv, mc%u0, mc%v0, mc%u1, mc%v1, mc%bu, mc%bv, 2_I4P, &
                                tr, ns(v), mc%d, prr, ptt, x(i), x(i + 1), y(j), y(j + 1), wr(v), rw(ro(v)), nxs(v), nps(v))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   do v = 1_I4P, nn
      c = pr%cell(pr%k(nj(v)))
      call map_node_2d_sizes(plan%lstart(c + 1) - plan%lstart(c), nxs(v), nps(v), wr(v), wi(v), mp)
   enddo
   ! the nodes, chunk by chunk, each chunk then added to its pairs' sums in node order
   a = 1_I4P
   do while (a <= nn)
      call map_chunk(a, nn, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, jp, c, i, j, l0, nl, vr, va, prr, ptt, tq) schedule(dynamic, 16)
      do v = a, b
         jp = nj(v) ; c = pr%cell(pr%k(jp)) ; i = mod(c - 1_I4P, nx) + 1_I4P ; j = (c - 1_I4P) / nx + 1_I4P
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0 ; vr = pr%vr(jp) ; va = pr%va(jp)
         call node_pose_rt(qr, qt, pid(v), 2_I4P, prr, ptt)
         tq = mc%pose0%t + ns(v) * mc%d
         call map_node_2d_eval(nl, plan%lel(l0), plan%lnext(l0), mc%nv, mc%u0, mc%v0, mc%u1, mc%v1, mc%bu, mc%bv, tr, ns(v), &
                               mc%d, prr, ptt, x(i), x(i + 1), y(j), y(j + 1), fl(:, v), nxs(v), nps(v), nw(v), tq, mc%dth, &
                               pr%narrow(jp), ep%nk(vr), ep%key(ep%kstart(vr):ep%kstart(vr + 1) - 1_I4P), &
                               ep%kpiece(ep%kstart(vr):ep%kstart(vr + 1) - 1_I4P), ep%np(vr), &
                               comp(ep%pstart(vr):ep%pstart(vr + 1) - 1_I4P), ep%nk(va), &
                               ep%key(ep%kstart(va):ep%kstart(va + 1) - 1_I4P), &
                               ep%kpiece(ep%kstart(va):ep%kstart(va + 1) - 1_I4P), ep%np(va), &
                               comp(ep%pstart(va):ep%pstart(va + 1) - 1_I4P), wr(v), rw(ro(v)), wi(v), iw(io(v)), np(v), &
                               nf(v), nsk(v))
      enddo
      !$omp end parallel do
      do v = a, b
         jp = nj(v)
         call map_node_outcome(pr, jp, nf(v), nsk(v), ok)
         if (.not. ok) cycle
         c = pr%cell(pr%k(jp))
         call map_item_2d_sizes(plan%lstart(c + 1) - plan%lstart(c), nxs(v), nps(v), nr, ni, mp, mk)
         call map_node_reduce(8_I4P, np(v), iw(io(v) + ni + 2_I4P * mk), rw(ro(v) + nr), ntr(c), &
                              acc(:, pr%col(jp):pr%col(jp) + ntr(c) - 1_I4P))
      enddo
      a = b + 1_I4P
   enddo
   endsubroutine map_nodes_2d

   pure subroutine map_node_2d_sizes(nl, nxs, nps, nrw, niw, mp)
   !< Workspace of a 2-D node (`map_node_2d_eval`, issue #105 step 1e): the evaluation's (`map_item_2d_sizes`: nr reals,
   !< ni integers), then its outputs: reals, the 8 values of every piece (8 mp); integers, the features and their pieces
   !< (2 mk), then the track of every piece (mp).
   !$acc routine seq
   integer(I4P), intent(in)  :: nl, nxs, nps !< Local vertices, crossings, wall pieces.
   integer(I4P), intent(out) :: nrw, niw, mp !< Reals, integers, piece capacity.
   integer(I4P)              :: nr, ni, mk   !< Evaluation workspace, feature capacity.

   call map_item_2d_sizes(nl, nxs, nps, nr, ni, mp, mk)
   nrw = nr + 8_I4P * mp ; niw = ni + 2_I4P * mk + mp
   endsubroutine map_node_2d_sizes

   subroutine map_node_2d_eval(nl, lg, lnx, nv, u0, v0, u1, v1, bu, bv, tr, s, d, pr, pt, x0, x1, y0, y1, fl, nxs, nps, w, &
                               tq, dth, narrow, nkr, rkey, rkpiece, nrp, rcomp, nka, akey, akpiece, npa, acomp, nrw, rw, niw, &
                               iw, np, fail, nskip)
   !< Kernel (issue #105 step 1e): one quadrature node of the 2-D measures (instant s, weight w), the per-cell
   !< `map_node` of issue #102 without the exact corner queries (the corner flags fl come with the node): the pieces at
   !< s (`map_item_2d_eval`, kind 2), each assigned to a track (`map_piece_assign`: the interval's evaluation r, and in
   !< a narrow interval the one after it, a), with its 8 values: w times the fluid lengths of the faces (-u, +u, -v, +v)
   !< and w times the rigid-field outflows through them, -dth (m - tq_v l) + d_u l through a face u = a, dth (m - tq_u
   !< l) + d_v l through a face v = a (l the fluid length, m its first moment along the face, tq = t0 + s d the
   !< translation at s). Outputs after the evaluation's workspace (`map_node_2d_sizes`): the values pv(8, mp) and the
   !< track of every piece (0: none). fail: outside a narrow interval the cell does not close or a piece has no track;
   !< nskip: in a narrow interval the node (1) or the pieces skipped.
   !$acc routine seq
   integer(I4P), intent(in)    :: nl                 !< Local vertices.
   integer(I4P), intent(in)    :: lg(nl), lnx(nl)    !< Their curve indexes, next local vertex.
   integer(I4P), intent(in)    :: nv                 !< Curve vertices.
   real(R8P),    intent(in)    :: u0(nv), v0(nv), u1(nv), v1(nv), bu(nv), bv(nv) !< Pose0 copy, pose1 copy, body.
   logical,      intent(in)    :: tr                 !< Translation path.
   real(R8P),    intent(in)    :: s, d(2), pr(2, 2), pt(2) !< Instant, translation over the step, pose at s.
   real(R8P),    intent(in)    :: x0, x1, y0, y1     !< Cell.
   logical,      intent(in)    :: fl(4)              !< Corner flags at s.
   integer(I4P), intent(in)    :: nxs, nps           !< The walk's counts (`map_item_2d_count`).
   real(R8P),    intent(in)    :: w, tq(2), dth      !< Weight, translation at s, rotation over the step.
   logical,      intent(in)    :: narrow             !< A narrow interval.
   integer(I4P), intent(in)    :: nkr                !< Features of r.
   integer(I4P), intent(in)    :: rkey(nkr), rkpiece(nkr) !< Its features, their pieces.
   integer(I4P), intent(in)    :: nrp                !< Pieces of r.
   integer(I4P), intent(in)    :: rcomp(nrp)         !< Their tracks.
   integer(I4P), intent(in)    :: nka                !< Features of a.
   integer(I4P), intent(in)    :: akey(nka), akpiece(nka) !< Its features, their pieces.
   integer(I4P), intent(in)    :: npa                !< Pieces of a.
   integer(I4P), intent(in)    :: acomp(npa)         !< Their tracks.
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes (`map_node_2d_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   integer(I4P), intent(out)   :: np                 !< Pieces.
   logical,      intent(out)   :: fail               !< The node cannot be assigned.
   integer(I4P), intent(out)   :: nskip              !< Skipped (narrow).
   real(R8P)                   :: rdum(1), l, mm, o, big !< Work.
   integer(I4P)                :: nr, ni, mp, mk, nr1, ni1, pa, pl, pq, nk, p, f, t, act, b !< Work.
   logical                     :: closed             !< The cell closed.

   call map_item_2d_sizes(nl, nxs, nps, nr, ni, mp, mk)
   call map_item_2d_layout(nl, nxs, nps, nr1, ni1, mp, pa, pl, pq)
   fail = .false. ; nskip = 0_I4P
   call map_item_2d_eval(nl, lg, lnx, nv, u0, v0, u1, v1, bu, bv, 2_I4P, tr, s, d, pr, pt, x0, x1, y0, y1, fl, nxs, nps, nr, &
                         rw, ni, iw, closed, np, nk, mk, iw(ni + 1_I4P), iw(ni + mk + 1_I4P), 0_I4P, rdum)
   if (.not. closed) then
      np = 0_I4P
      if (narrow) then
         nskip = 1_I4P
      else
         fail = .true.
      endif
      return
   endif
   do p = 1_I4P, np
      b = 4_I4P * (p - 1_I4P) ; iw(ni + 2_I4P * mk + p) = 0_I4P
      big = rw(pl + b)
      do f = 2_I4P, 4_I4P
         if (rw(pl + b + f - 1_I4P) > big) big = rw(pl + b + f - 1_I4P)
      enddo
      call map_piece_assign(nk, iw(ni + 1_I4P), iw(ni + mk + 1_I4P), p, narrow, big <= 1e-9_R8P * max(x1 - x0, y1 - y0), &
                            nkr, rkey, rkpiece, nrp, rcomp, nka, akey, akpiece, npa, acomp, t, act)
      if (act == 3_I4P) then
         fail = .true. ; return
      elseif (act == 2_I4P) then
         nskip = nskip + 1_I4P ; cycle
      elseif (act == 0_I4P) then
         cycle
      endif
      iw(ni + 2_I4P * mk + p) = t
      do f = 1_I4P, 4_I4P
         l = rw(pl + b + f - 1_I4P) ; mm = rw(pq + b + f - 1_I4P)
         if (f <= 2_I4P) then
            o = -dth * (mm - tq(2) * l) + d(1) * l
         else
            o = dth * (mm - tq(1) * l) + d(2) * l
         endif
         rw(nr + 2_I4P * b + f) = w * l ; rw(nr + 2_I4P * b + 4_I4P + f) = w * o
      enddo
   enddo
   endsubroutine map_node_2d_eval

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

   subroutine piece_map_3d_build(surface, s0, pose0, pose1, x, y, z, k0, k1, sw, fox, foy, foz, path, nq, map, all_meas, measured)
   !< The 3-D piece map of the step (issue #102 step 2; `piece_map_3d_t`), from the end-pose cuts k0, k1 of the copies
   !< of `surface` moved to pose0 (s0) and pose1; the 3-D twin of `piece_map_2d_build`.
   !<
   !< 1. Candidates: the cells a facet can meet during the step (the box of its vertices' paths: the segment of the
   !<    translation, or 17 samples of the screw path grown by its curvature bound), each with its facets (CSR,
   !<    increasing facet order). The host phase (`map_plan_3d`, issue #105 step 1c) then gives every cell its status
   !<    and, if it is evaluated, its facets, events, corner events, interval layout and evaluation instants, as flat CSR
   !<    arrays, with a pose table.
   !< 2. Topology (from the plan): the events of a cell from its facets (`map3_events`: a mesh vertex crossing a face
   !<    plane inside the face, a mesh edge crossing a cell edge, a corner crossing a facet; explicit on a translation,
   !<    `isolate_zeros` under a rotation); a cell without events whose two end cuts have the same single piece (or
   !<    none) takes the trivial track. Otherwise its pieces are evaluated at pose0 (the pose0 copy, k0's corner states
   !<    and fractions: they must equal k0's pieces bit for bit, which fixes the numbering), at the middle of every
   !<    interval between events, and at pose1: every (cell, evaluation) one item of the flat evaluations
   !<    (`map_evals_3d`, issue #105 step 1d). The end evaluations are checked (`map_end_check`), and pieces of
   !<    consecutive evaluations sharing a boundary feature (a facet or a fluid corner) are one track (`map_tracks`).
   !< 3. Measures: one track takes the cell's `sw` values and the face outflows of the call (fox, foy, foz: the time
   !<    integral of the rigid field's flux in +axis through every grid face); several tracks integrate their face
   !<    parts on the cell's event intervals with the exact paths' rules, on the Gauss path too (`map_measures_3d`: by
   !<    rounds over flat (cell, interval, node) items, issue #105 step 1e).
   !< 4. Assembly (issue #105 step 1f): tracks per cell and their scan (`map_asm_tracks`), the old and new pieces per
   !<    track by count / scan / fill per cell (`map_asm_pieces`), the statuses and measures per cell (`map_asm_status`,
   !<    `map_asm_measures_*`); the same bits as the former loops in cell and piece order.
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
   type(map_plan_t)          :: plan                       !< Host phase (issue #105 step 1c).
   type(map_evals_t)         :: ep                         !< Evaluations (issue #105 step 1d).
   integer(I4P), allocatable :: cstart(:), cfac(:)         !< Facets of every cell.
   integer(I4P), allocatable :: ntr(:), ckind(:)           !< Tracks and kind (0 trivial, 1 evaluated, -1 untrusted).
   integer(I4P), allocatable :: otr(:), ntn(:)             !< Track (in its cell) of every old and new piece.
   integer(I4P), allocatable :: comp(:)                    !< Track (in its cell) of every piece of every evaluation.
   integer(I4P), allocatable :: tw(:)                      !< Work of `map_tracks`.
   integer(I4P), allocatable :: moff(:)                    !< First measure of every measured cell's tracks.
   real(R8P),    allocatable :: meas(:,:)                  !< Per-track measures of the measured cells.
   logical,      allocatable :: cok(:)                     !< Measured cells whose measures are assigned.
   integer(I4P)              :: nx, ny, nz, ncell, c, i, j, k, nt, tot, t, a0, a1, b0, b1, v0, ne, p0, nn, q0, nk !< Work.
   logical                   :: allm                       !< all_meas.

   allm = .false. ; if (present(all_meas)) allm = all_meas
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   ncell = nx * ny * nz
   map%nx = nx ; map%ny = ny ; map%nz = nz ; map%np0 = k0%np ; map%np1 = k1%np
   call map3_context(surface, s0, pose0, pose1, x, y, z, path, nq, mc)
   call map3_candidates(mc, x, y, z, cstart, cfac)
   call map_plan_3d(mc, x, y, z, k0, k1, cstart, cfac, plan)
   deallocate(cstart, cfac)
   allocate(ntr(ncell), ckind(ncell), otr(k0%np), ntn(k1%np), moff(ncell))
   ! the evaluations of every evaluated cell, flat (issue #105 step 1d), then the tracks of every cell
   call map_evals_3d(mc, plan, surface, s0, x, y, z, k0, k1, ep)
   allocate(comp(ep%pstart(ep%n + 1) - 1_I4P), tw(map_tracks_work(plan, ep)))
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            a0 = k0%first_piece(i, j, k) ; a1 = a0 + max(0_I4P, k0%pieces(i, j, k)) - 1_I4P
            b0 = k1%first_piece(i, j, k) ; b1 = b0 + max(0_I4P, k1%pieces(i, j, k)) - 1_I4P
            ntr(c) = 0_I4P
            select case (plan%status(c))
            case (PLAN_UNTRUSTED)
               ckind(c) = -1_I4P
            case (PLAN_TRIVIAL)
               ckind(c) = 0_I4P ; ntr(c) = k0%pieces(i, j, k) ; otr(a0:a1) = 1_I4P ; ntn(b0:b1) = 1_I4P
            case default
               v0 = plan%vstart(c) ; ne = plan%vstart(c + 1) - v0
               map%evaluations = map%evaluations + ne ; map%evaluated_cells = map%evaluated_cells + 1_I4P
               ckind(c) = map_end_check(ep, v0, ne, k0%pieces(i, j, k), k1%pieces(i, j, k), k0%piece_volume(a0:a1), &
                                        k1%piece_volume(b0:b1))
               if (ckind(c) == 1_I4P) then ! tracks: pieces of consecutive evaluations sharing a feature
                  p0 = ep%pstart(v0) ; nn = ep%pstart(v0 + ne) - p0 ; q0 = ep%kstart(v0) ; nk = ep%kstart(v0 + ne) - q0
                  call map_tracks(ne, ep%pstart(v0:v0 + ne), ep%kstart(v0:v0 + ne), nk, ep%key(q0:q0 + nk - 1_I4P), &
                                  ep%kpiece(q0:q0 + nk - 1_I4P), nn, tw, comp(p0:p0 + nn - 1_I4P), ntr(c))
                  otr(a0:a1) = comp(ep%pstart(v0):ep%pstart(v0 + 1_I4P) - 1_I4P)
                  ntn(b0:b1) = comp(ep%pstart(v0 + ne - 1_I4P):ep%pstart(v0 + ne) - 1_I4P)
               endif
            endselect
         enddo
      enddo
   enddo
   ! the per-track measures of the multi-track cells (all evaluated cells with all_meas), from the flat evaluations
   tot = 0_I4P ; moff = 0_I4P
   do c = 1_I4P, ncell
      if (ckind(c) == 1_I4P .and. (ntr(c) >= 2_I4P .or. allm)) then
         moff(c) = tot + 1_I4P ; tot = tot + ntr(c)
      endif
   enddo
   allocate(meas(12, tot), cok(ncell)) ; meas = 0._R8P
   call map_measures_3d(mc, plan, ep, comp, surface, s0, x, y, z, k0, k1, ntr, moff, meas, cok)
   do c = 1_I4P, ncell
      if (moff(c) == 0_I4P .or. ntr(c) == 0_I4P) cycle ! not measured, or no track (solid throughout): measured
      if (cok(c)) cycle
      if (ntr(c) == 1_I4P .and. allm) then ! one track (all_meas): keep it with the sw values
         moff(c) = 0_I4P
      else
         ckind(c) = -1_I4P ; moff(c) = 0_I4P
      endif
   enddo
   if (present(measured)) measured = moff > 0_I4P
   ! assembly (issue #105 step 1f): tracks per cell, scan; pieces per track, scan, fill; status; measures
   allocate(map%tracks(nx, ny, nz), map%first_track(nx, ny, nz))
   call map_asm_tracks(ncell, ckind, ntr, map%first_track, nt, map%untrusted_cells, map%multi_track_cells)
   map%tracks = reshape(ntr, [nx, ny, nz]) ; map%nt = nt
   allocate(map%track_cell(3, nt), map%status(nt), map%old_start(nt + 1), map%new_start(nt + 1), map%old_piece(k0%np), &
            map%new_piece(k1%np), map%old_track(k0%np), map%new_track(k1%np), map%volume0(nt), map%volume1(nt), &
            map%faces(6, nt), map%outflow(6, nt), map%wall(3, nt), map%residual(nt))
   call map_asm_pieces(ncell, ckind, ntr, map%first_track, k0%first_piece, k0%pieces, k0%np, otr, k0%piece_volume, nt, &
                       map%old_start, map%old_piece, map%old_track, map%volume0)
   call map_asm_pieces(ncell, ckind, ntr, map%first_track, k1%first_piece, k1%pieces, k1%np, ntn, k1%piece_volume, nt, &
                       map%new_start, map%new_piece, map%new_track, map%volume1)
   !$omp parallel do default(shared) private(c, i, j, k) schedule(static)
   do c = 1_I4P, ncell
      i = 1_I4P + mod(c - 1_I4P, nx) ; j = 1_I4P + mod((c - 1_I4P) / nx, ny) ; k = 1_I4P + (c - 1_I4P) / (nx * ny)
      call map_asm_status(3_I4P, [i, j, k], ckind(c) == -1_I4P, map%first_track(i, j, k), ntr(c), nt, map%old_start, &
                          map%new_start, map%status, map%track_cell)
      call map_asm_measures_3d(map%first_track(i, j, k), ntr(c), moff(c), size(meas, 2, kind=I4P), meas, &
                               [sw%face_x(i, j, k), sw%face_x(i + 1, j, k), sw%face_y(i, j, k), sw%face_y(i, j + 1, k), &
                                sw%face_z(i, j, k), sw%face_z(i, j, k + 1)], &
                               [-fox(i, j, k), fox(i + 1, j, k), -foy(i, j, k), foy(i, j + 1, k), -foz(i, j, k), &
                                foz(i, j, k + 1)], sw%wall(:, i, j, k), sw%residual(i, j, k), &
                               nt, map%volume0, map%volume1, map%faces, map%outflow, map%wall, map%residual)
   enddo
   !$omp end parallel do
   do t = 1_I4P, nt
      map%n_status(map%status(t)) = map%n_status(map%status(t)) + 1_I4P
   enddo
   map%exact = mc%exact ; map%skipped_nodes = mc%skipped
   map%contact_cells = plan%contacts ; map%narrow_cells = plan%resolved
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
   mc%exact = .true. ; mc%skipped = 0_I4P
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

   subroutine map_plan_3d(mc, x, y, z, k0, k1, cstart, cfac, plan)
   !< The host phase of the 3-D piece map (issue #105 step 1c; `map_plan_t`), the twin of `map_plan_2d`: per cell the
   !< checks of the map that come before its evaluations (both end cuts trusted, a wall near the cell during the
   !< step, its events certified and free of uncertified contacts, more than one state); an evaluated cell also gets its
   !< facets, its events and corner events (`map3_events`, a per-cell deduplication of mesh vertices and edges), its
   !< interval layout and its evaluation list; then the pose table (`plan_instants`, with a rotation).
   type(map3_ctx_t),     intent(in)  :: mc                 !< State.
   real(R8P),            intent(in)  :: x(:), y(:), z(:)   !< Cell edges.
   type(cut_cells_3d_t), intent(in)  :: k0, k1             !< End-pose cuts.
   integer(I4P),         intent(in)  :: cstart(:), cfac(:) !< Facets of every cell (`map3_candidates`).
   type(map_plan_t),     intent(out) :: plan               !< Plan.
   type(plan_cell_t), allocatable :: pc(:)                 !< Candidate cells.
   integer(I4P),      allocatable :: cand(:)               !< Their cells.
   type(rigid_pose_t)             :: ps                    !< Pose.
   real(R8P)                      :: b0(3), b1(3)          !< Cell.
   integer(I4P)                   :: nx, ny, nz, ncell, nc, ic, c, i, j, k, np0, np1, nb !< Work.
   logical                        :: cert, unc             !< Events certified, an uncertified contact.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   ncell = nx * ny * nz
   plan%ncell = ncell ; allocate(plan%status(ncell))
   nc = 0_I4P
   do c = 1_I4P, ncell
      call cell_ijk(c, nx, ny, i, j, k)
      np0 = k0%pieces(i, j, k) ; np1 = k1%pieces(i, j, k)
      if (np0 < 0_I4P .or. np1 < 0_I4P) then
         plan%status(c) = PLAN_UNTRUSTED
      elseif (cstart(c + 1) == cstart(c)) then ! no wall near the cell during the step: the same state throughout
         plan%status(c) = merge(PLAN_TRIVIAL, PLAN_UNTRUSTED, np0 == np1 .and. np0 <= 1_I4P)
      else
         plan%status(c) = PLAN_EVAL ; nc = nc + 1_I4P
      endif
   enddo
   allocate(cand(nc), pc(nc))
   nc = 0_I4P
   do c = 1_I4P, ncell
      if (plan%status(c) /= PLAN_EVAL) cycle
      nc = nc + 1_I4P ; cand(nc) = c
   enddo
   !$omp parallel do default(shared) private(ic, c, i, j, k, b0, b1, nb, cert, unc) schedule(dynamic, 4)
   do ic = 1_I4P, nc ! independent cells, each into its own pc(ic)
      c = cand(ic) ; call cell_ijk(c, nx, ny, i, j, k)
      b0 = [x(i), y(j), z(k)] ; b1 = [x(i + 1), y(j + 1), z(k + 1)]
      pc(ic)%lg = cfac(cstart(c):cstart(c + 1) - 1_I4P)
      call map3_events(mc, b0, b1, pc(ic)%lg, pc(ic)%b, pc(ic)%bw, nb, cert, unc, pc(ic)%ct, pc(ic)%cc)
      call plan_classify(pc(ic), cert, unc, k0%pieces(i, j, k), k1%pieces(i, j, k))
   enddo
   !$omp end parallel do
   call plan_fill(plan, cand, pc, 3_I4P)
   if (mc%rot) then ! the middle evaluations placed by the screw pose (`map_place_3d`)
      call plan_instants(plan, 3_I4P)
      do k = 1_I4P, plan%npose
         ps = pose_interpolate(mc%pose0, mc%pose1, plan%pose_s(k))
         plan%pose_r(:, :, k) = ps%r ; plan%pose_t(:, k) = ps%t
      enddo
   endif
   endsubroutine map_plan_3d

   pure subroutine cell_ijk(c, nx, ny, i, j, k)
   !< The indexes of cell c = i + nx ((j - 1) + ny (k - 1)).
   integer(I4P), intent(in)  :: c, nx, ny !< Cell, cells along x and y.
   integer(I4P), intent(out) :: i, j, k   !< Indexes.

   i = mod(c - 1_I4P, nx) + 1_I4P ; j = mod((c - 1_I4P) / nx, ny) + 1_I4P ; k = (c - 1_I4P) / (nx * ny) + 1_I4P
   endsubroutine cell_ijk

   pure function plan_pose_3d(plan, k) result(ps)
   !< Pose k of the plan's table (the identity for k = 0).
   type(map_plan_t), intent(in) :: plan !< Plan.
   integer(I4P),     intent(in) :: k    !< Pose (0: none).
   type(rigid_pose_t)           :: ps   !< Pose.

   if (k > 0_I4P) then
      ps%r = plan%pose_r(:, :, k) ; ps%t = plan%pose_t(:, k)
   endif
   endfunction plan_pose_3d

   subroutine map_states_3d(mc, plan, surface, s0, x, y, z, k0, k1, ep, anyk)
   !< The corner states of every 3-D evaluation (issue #105 step 1d; `map_evals_t`): pose0 and pose1 take the end cuts'
   !< node states, a middle one the flip rule where it holds and the exact query `map3_corner_in` elsewhere
   !< (`map3_corner_state`, as the per-cell `map3_eval` of issue #102); anyk: a corner of the item's cell is known by
   !< the flip rule (a failed middle evaluation is then retried with the exact states). Host.
   type(map3_ctx_t),         intent(in)    :: mc               !< State.
   type(map_plan_t),         intent(in)    :: plan             !< Host phase.
   type(surface_stl_object), intent(in)    :: surface, s0      !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(in)    :: k0, k1           !< End-pose cuts.
   type(map_evals_t),        intent(inout) :: ep               !< Evaluations (cs set).
   logical, allocatable,     intent(out)   :: anyk(:)          !< A corner of the item's cell known [n].
   real(R8P)                               :: c0(3), c1(3), b0(3), b1(3), s !< Cell, shifted cell, instant.
   integer(I4P)                            :: nx, ny, c, i, j, k, v, q, n0, n1, cs0(8), cs1(8) !< Work.
   logical                                 :: known(8)         !< Flip rule.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   ep%n = plan%vstart(plan%ncell + 1) - 1_I4P
   allocate(ep%cs(8, ep%n), anyk(ep%n))
   do c = 1_I4P, plan%ncell
      if (plan%status(c) /= PLAN_EVAL) cycle
      call cell_ijk(c, nx, ny, i, j, k)
      c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
      call node_states(k0%node_in, i, j, k, cs0) ; call node_states(k1%node_in, i, j, k, cs1)
      n0 = plan%cnstart(c) ; n1 = plan%cnstart(c + 1) - 1_I4P
      call map_corner_known(8_I4P, plan%cntime(n0:n1), plan%cncor(n0:n1), cs0, cs1, known)
      do v = plan%vstart(c), plan%vstart(c + 1) - 1_I4P
         anyk(v) = any(known)
         select case (plan%vkind(v))
         case (0_I4P)
            ep%cs(:, v) = cs0
         case (1_I4P)
            ep%cs(:, v) = cs1
         case default
            s = plan%vs(v) ; b0 = c0 ; b1 = c1
            if (.not. mc%rot) then
               b0 = c0 - s * mc%d ; b1 = c1 - s * mc%d
            endif
            do q = 1_I4P, 8_I4P
               ep%cs(q, v) = map3_corner_state(mc, surface, s0, s, plan_pose_3d(plan, plan%vpose(v)), c0, c1, b0, b1, q, &
                                               known(q), cs0, plan%cntime(n0:n1), plan%cncor(n0:n1))
            enddo
         endselect
      enddo
   enddo
   endsubroutine map_states_3d

   subroutine map_evals_3d(mc, plan, surface, s0, x, y, z, k0, k1, ep)
   !< The evaluations of the 3-D map (issue #105 step 1d; `map_evals_t`, the twin of `map_evals_2d`): every item's
   !< pieces and features from its cell's facets placed at its instant (`map_place_3d`). Corner states by
   !< `map_states_3d`; sizes pass (`map_item_3d_count`: the piercings); the evaluation pass (`map_item_3d_eval`) in
   !< chunks; a middle evaluation that fails with a corner known by the flip rule is evaluated again with the exact
   !< corner states (the last retry of the per-cell `map3_eval` of issue #102); then the CSR. Bit for bit
   !< that `map3_eval`.
   type(map3_ctx_t),         intent(in)    :: mc               !< State.
   type(map_plan_t),         intent(in)    :: plan             !< Host phase.
   type(surface_stl_object), intent(in)    :: surface, s0      !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(in)    :: k0, k1           !< End-pose cuts.
   type(map_evals_t),        intent(inout) :: ep               !< Evaluations.
   integer(I4P), allocatable :: vc(:), npc(:), wr(:), wi(:), kcap(:), keyc(:), kpc(:), mode(:), ro(:), io(:) !< Work.
   integer(I4P), allocatable :: rl(:)                          !< Items to evaluate again.
   real(R8P),    allocatable :: rw(:), al(:), fr(:,:)          !< Workspace, end fractions.
   integer(I4P), allocatable :: iw(:)                          !< Workspace.
   logical,      allocatable :: anyk(:), vm(:)                 !< A corner known, piece measures.
   real(R8P)                 :: pr(3, 3), pt(3), c0(3), c1(3), b0(3), b1(3), s !< Work.
   integer(I4P)              :: nx, ny, n, c, i, j, k, v, a, b, l0, nl, mk, nr, q, r, np0, np1 !< Work.
   integer(I4P)              :: csz(8)                         !< Unused pose0 states (exact queries).

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   call map_states_3d(mc, plan, surface, s0, x, y, z, k0, k1, ep, anyk)
   n = ep%n ; csz = 0_I4P
   allocate(vc(n), npc(n), wr(n), wi(n), kcap(n + 1), mode(n), ro(n + 1), io(n + 1), al(n), fr(6, n), vm(n))
   do c = 1_I4P, plan%ncell
      if (plan%vstart(c + 1) == plan%vstart(c)) cycle
      call cell_ijk(c, nx, ny, i, j, k)
      np0 = k0%pieces(i, j, k) ; np1 = k1%pieces(i, j, k)
      do v = plan%vstart(c), plan%vstart(c + 1) - 1_I4P
         vc(v) = c ; nl = plan%lstart(c + 1) - plan%lstart(c) ; wr(v) = 9_I4P * nl ; wi(v) = nl
         select case (plan%vkind(v))
         case (1_I4P)
            mode(v) = 0_I4P ; al(v) = k1%alpha(i, j, k) ; vm(v) = np1 > 1_I4P
            fr(:, v) = [k1%face_x(i, j, k), k1%face_x(i + 1, j, k), k1%face_y(i, j, k), k1%face_y(i, j + 1, k), &
                        k1%face_z(i, j, k), k1%face_z(i, j, k + 1)]
         case (0_I4P)
            mode(v) = 0_I4P ; al(v) = k0%alpha(i, j, k) ; vm(v) = np0 > 1_I4P
            fr(:, v) = [k0%face_x(i, j, k), k0%face_x(i + 1, j, k), k0%face_y(i, j, k), k0%face_y(i, j + 1, k), &
                        k0%face_z(i, j, k), k0%face_z(i, j, k + 1)]
         case default
            mode(v) = 1_I4P ; al(v) = 0._R8P ; vm(v) = .false. ; fr(:, v) = 0._R8P
         endselect
      enddo
   enddo
   ! sizes pass: the placed facets, the piercings
   a = 1_I4P
   do while (a <= n)
      call map_chunk(a, n, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, c, i, j, k, l0, nl, pr, pt, c0, c1) schedule(dynamic, 16)
      do v = a, b
         c = vc(v) ; call cell_ijk(c, nx, ny, i, j, k)
         c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call plan_pose_rt(plan, plan%vpose(v), 3_I4P, pr, pt)
         call map_item_3d_count(nl, plan%lel(l0), mc%nf, mc%t0, mc%t1, mc%tb, plan%vkind(v), mc%rot, plan%vs(v), mc%d, pr, &
                                pt, c0, c1, wr(v), rw(ro(v)), wi(v), iw(io(v)), npc(v))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   ! workspaces and capacities: features, piece volumes (the end evaluations with several pieces)
   allocate(ep%mstart(n + 1))
   do v = 1_I4P, n
      nl = plan%lstart(vc(v) + 1) - plan%lstart(vc(v))
      call map_item_3d_sizes(nl, npc(v), plan%vkind(v), vm(v), wr(v), wi(v), mk, ep%mstart(v))
      kcap(v) = mk
   enddo
   call csr_scan(kcap) ; call csr_scan(ep%mstart)
   allocate(ep%closed(n), ep%np(n), ep%nk(n), keyc(kcap(n + 1) - 1_I4P), kpc(kcap(n + 1) - 1_I4P), &
            ep%pm(ep%mstart(n + 1))) ! + 1: every slot start inside
   ! evaluation pass, then the failed middle evaluations with a known corner again, with the exact corner states
   do r = 1_I4P, 2_I4P
      if (r == 1_I4P) then
         rl = [(v, v = 1_I4P, n)]
      else
         nr = 0_I4P
         do v = 1_I4P, n
            if (mode(v) == 1_I4P .and. .not. ep%closed(v) .and. anyk(v)) nr = nr + 1_I4P
         enddo
         deallocate(rl) ; allocate(rl(nr)) ; nr = 0_I4P
         do v = 1_I4P, n
            if (mode(v) == 1_I4P .and. .not. ep%closed(v) .and. anyk(v)) then
               nr = nr + 1_I4P ; rl(nr) = v ; mode(v) = 2_I4P
               c = vc(v) ; call cell_ijk(c, nx, ny, i, j, k)
               c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
               s = plan%vs(v) ; b0 = c0 ; b1 = c1
               if (.not. mc%rot) then
                  b0 = c0 - s * mc%d ; b1 = c1 - s * mc%d
               endif
               do q = 1_I4P, 8_I4P
                  ep%cs(q, v) = map3_corner_state(mc, surface, s0, s, plan_pose_3d(plan, plan%vpose(v)), c0, c1, b0, b1, &
                                                  q, .false., csz, plan%cntime(1:0), plan%cncor(1:0))
               enddo
            endif
         enddo
      endif
      call map_evals_3d_pass(mc, plan, x, y, z, rl, vc, mode, al, fr, vm, npc, wr, wi, kcap, keyc, kpc, ep)
   enddo
   call map_evals_fill(ep, kcap, keyc, kpc)
   endsubroutine map_evals_3d

   subroutine map_evals_3d_pass(mc, plan, x, y, z, rl, vc, mode, al, fr, vm, npc, wr, wi, kcap, keyc, kpc, ep)
   !< The evaluation pass of `map_evals_3d` over the items rl, in chunks of at most MAP_CHUNK_BYTES of workspace.
   type(map3_ctx_t),  intent(in)    :: mc               !< State.
   type(map_plan_t),  intent(in)    :: plan             !< Host phase.
   real(R8P),         intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),      intent(in)    :: rl(:)            !< Items (increasing).
   integer(I4P),      intent(in)    :: vc(:), mode(:)   !< Cell and mode of every item.
   real(R8P),         intent(in)    :: al(:), fr(:,:)   !< End fractions of every item.
   logical,           intent(in)    :: vm(:)            !< Piece measures of every item.
   integer(I4P),      intent(in)    :: npc(:), wr(:), wi(:), kcap(:) !< Piercings, workspace, feature slots.
   integer(I4P),      intent(inout) :: keyc(*), kpc(*)  !< Features in their slots.
   type(map_evals_t), intent(inout) :: ep               !< Evaluations.
   integer(I4P), allocatable        :: rr(:), ri(:), ro(:), io(:), iw(:) !< Work.
   real(R8P),    allocatable        :: rw(:)            !< Workspace.
   real(R8P)                        :: pr(3, 3), pt(3), c0(3), c1(3) !< Work.
   integer(I4P)                     :: nx, ny, m, a, b, h, v, c, i, j, k, l0, nl !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; m = size(rl, kind=I4P)
   allocate(rr(m), ri(m), ro(m + 1), io(m + 1))
   do h = 1_I4P, m
      rr(h) = wr(rl(h)) ; ri(h) = wi(rl(h))
   enddo
   a = 1_I4P
   do while (a <= m)
      call map_chunk(a, m, rr, ri, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(h, v, c, i, j, k, l0, nl, pr, pt, c0, c1) schedule(dynamic, 16)
      do h = a, b
         v = rl(h) ; c = vc(v) ; call cell_ijk(c, nx, ny, i, j, k)
         c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call plan_pose_rt(plan, plan%vpose(v), 3_I4P, pr, pt)
         call map_item_3d_eval(nl, plan%lel(l0), mc%nf, mc%t0, mc%t1, mc%tb, plan%vkind(v), mc%rot, plan%vs(v), mc%d, pr, &
                               pt, c0, c1, ep%cs(:, v), al(v), fr(:, v), vm(v), mode(v), npc(v), wr(v), rw(ro(h)), wi(v), &
                               iw(io(h)), ep%closed(v), ep%np(v), ep%nk(v), kcap(v + 1) - kcap(v), keyc(kcap(v)), &
                               kpc(kcap(v)), ep%mstart(v + 1) - ep%mstart(v), ep%pm(ep%mstart(v)))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   endsubroutine map_evals_3d_pass

   pure subroutine map_place_3d(n, lf, nf, t0, t1, tb, kind, rot, s, d, pr, pt, c0, c1, tl, b0, b1)
   !< Kernel (issue #105 step 1d): the facets lf of a 3-D evaluation, placed, and its cell: kind 0 the pose0 copy, 1 the
   !< pose1 copy, 2 the instant s (without rotation the pose0 copy and the cell shifted back by s d, else the body
   !< facets by the pose (pr, pt) of the plan's table, R q + t as `to_world`). As the per-cell `map3_eval` of issue
   !< #102, bit for bit.
   !$acc routine seq
   integer(I4P), intent(in)  :: n                    !< Facets.
   integer(I4P), intent(in)  :: lf(n)                !< Their ids.
   integer(I4P), intent(in)  :: nf                   !< Surface facets.
   real(R8P),    intent(in)  :: t0(3, 3, nf), t1(3, 3, nf), tb(3, 3, nf) !< Pose0 copy, pose1 copy, body frame.
   integer(I4P), intent(in)  :: kind                 !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)  :: rot                  !< The step rotates.
   real(R8P),    intent(in)  :: s, d(3)              !< Instant, translation over the step.
   real(R8P),    intent(in)  :: pr(3, 3), pt(3)      !< Pose at s (kind 2 with a rotation).
   real(R8P),    intent(in)  :: c0(3), c1(3)         !< Cell.
   real(R8P),    intent(out) :: tl(3, 3, n)          !< Placed facets.
   real(R8P),    intent(out) :: b0(3), b1(3)         !< Cell (shifted back on the translation path).
   real(R8P)                 :: w(3)                 !< Rotated vertex.
   integer(I4P)              :: q, m, a              !< Work.

   do a = 1_I4P, 3_I4P
      b0(a) = c0(a) ; b1(a) = c1(a)
   enddo
   do q = 1_I4P, n
      do m = 1_I4P, 3_I4P
         select case (kind)
         case (0_I4P)
            do a = 1_I4P, 3_I4P
               tl(a, m, q) = t0(a, m, lf(q))
            enddo
         case (1_I4P)
            do a = 1_I4P, 3_I4P
               tl(a, m, q) = t1(a, m, lf(q))
            enddo
         case default
            if (.not. rot) then
               do a = 1_I4P, 3_I4P
                  tl(a, m, q) = t0(a, m, lf(q))
               enddo
            else
               do a = 1_I4P, 3_I4P
                  w(a) = pr(a, 1) * tb(1, m, lf(q)) + pr(a, 2) * tb(2, m, lf(q)) + pr(a, 3) * tb(3, m, lf(q))
               enddo
               do a = 1_I4P, 3_I4P
                  tl(a, m, q) = w(a) + pt(a)
               enddo
            endif
         endselect
      enddo
   enddo
   if (kind == 2_I4P .and. .not. rot) then
      do a = 1_I4P, 3_I4P
         b0(a) = c0(a) - s * d(a) ; b1(a) = c1(a) - s * d(a)
      enddo
   endif
   endsubroutine map_place_3d

   pure subroutine map_item_3d_sizes(n, npc, kind, measures, nrw, niw, mk, mm)
   !< Workspace of a 3-D evaluation (`map_item_3d_eval`) with n facets and npc piercings (issue #105 step 1d): reals, the
   !< placed facets (9 n) and a region that holds either the one-cell workspace (`cut_cell_3d_one_sizes`) with its piece
   !< outputs (volumes, face fractions, face moments: 25 npo) or, for a middle evaluation (kind 2), the face sections
   !< (`cut_cell_3d_face_sizes` at the bound of n segments: at most 4 n crossings and n wall pieces), + 1; integers, the
   !< facets' pieces (n) and a region for the one-cell or face workspace or the candidate list (n). mk = 8 + n features
   !< at most; mm = npo piece volumes kept (with measures, the end evaluations with several pieces).
   !$acc routine seq
   integer(I4P), intent(in)  :: n, npc            !< Facets, piercings.
   integer(I4P), intent(in)  :: kind              !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)  :: measures          !< Piece measures.
   integer(I4P), intent(out) :: nrw, niw, mk, mm  !< Reals, integers, feature and volume capacities.
   integer(I4P)              :: nrp, nip, npo, nrf, nif, rr, ii !< Work.

   call cut_cell_3d_one_sizes(n, npc, measures, nrp, nip, npo)
   rr = nrp + 25_I4P * npo ; ii = max(nip, n)
   if (kind == 2_I4P) then
      call cut_cell_3d_face_sizes(n, 4_I4P * n, n, nrf, nif)
      rr = max(rr, nrf) ; ii = max(ii, nif)
   endif
   nrw = 9_I4P * n + rr + 1_I4P ; niw = n + ii ; mk = 8_I4P + n ; mm = npo
   endsubroutine map_item_3d_sizes

   subroutine map_item_3d_count(n, lf, nf, t0, t1, tb, kind, rot, s, d, pr, pt, c0, c1, nrw, rw, niw, iw, npc)
   !< Kernel (issue #105 step 1d), sizes pass of a 3-D evaluation: the placed facets (`map_place_3d`, in rw) and the count
   !< pass of the cell's piercings (`cut_cell_3d_one_count`).
   !$acc routine seq
   integer(I4P), intent(in)    :: n                  !< Facets.
   integer(I4P), intent(in)    :: lf(n)              !< Their ids.
   integer(I4P), intent(in)    :: nf                 !< Surface facets.
   real(R8P),    intent(in)    :: t0(3, 3, nf), t1(3, 3, nf), tb(3, 3, nf) !< Pose0 copy, pose1 copy, body frame.
   integer(I4P), intent(in)    :: kind               !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)    :: rot                !< The step rotates.
   real(R8P),    intent(in)    :: s, d(3), pr(3, 3), pt(3) !< Instant, translation, pose at s.
   real(R8P),    intent(in)    :: c0(3), c1(3)       !< Cell.
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes (>= 9 n, >= n).
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   integer(I4P), intent(out)   :: npc                !< Piercings.
   real(R8P)                   :: b0(3), b1(3)       !< Cell (shifted).

   call map_place_3d(n, lf, nf, t0, t1, tb, kind, rot, s, d, pr, pt, c0, c1, rw, b0, b1)
   call cut_cell_3d_one_count(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), niw, iw, npc)
   endsubroutine map_item_3d_count

   subroutine map_item_3d_eval(n, lf, nf, t0, t1, tb, kind, rot, s, d, pr, pt, c0, c1, cs, al_in, fr_in, measures, mode, &
                               npc, nrw, rw, niw, iw, closed, np, nk, mk, key, kpiece, mm, pvol)
   !< Kernel (issue #105 step 1d): one 3-D evaluation on its workspace (`map_item_3d_sizes`), the per-cell `map3_eval`
   !< of issue #102 without the exact corner queries (the corner states cs come with the item). Mode 0 (an end
   !< evaluation): the end cut's fractions al_in, fr_in. Mode 1 (a middle one): fractions of -1, which fail only a
   !< multi-piece split; on a failure the six face sections (`map3_sections_ws`), alpha from them
   !< (`cut_cell_3d_one_totals`) and the pieces again. Mode 2 (a middle one again, with the exact corner states): the
   !< sections first. Features, increasing: the fluid corners (keys -8..-1) then the facets (ids), each with its piece;
   !< with measures and several pieces the piece volumes into pvol (capacity mm). Not trusted: closed false, np = nk =
   !< 0.
   !$acc routine seq
   integer(I4P), intent(in)    :: n                  !< Facets.
   integer(I4P), intent(in)    :: lf(n)              !< Their ids.
   integer(I4P), intent(in)    :: nf                 !< Surface facets.
   real(R8P),    intent(in)    :: t0(3, 3, nf), t1(3, 3, nf), tb(3, 3, nf) !< Pose0 copy, pose1 copy, body frame.
   integer(I4P), intent(in)    :: kind               !< 0 pose0, 1 pose1, 2 instant s.
   logical,      intent(in)    :: rot                !< The step rotates.
   real(R8P),    intent(in)    :: s, d(3), pr(3, 3), pt(3) !< Instant, translation, pose at s.
   real(R8P),    intent(in)    :: c0(3), c1(3)       !< Cell.
   integer(I4P), intent(in)    :: cs(8)              !< Corner states (0 solid, 1 fluid).
   real(R8P),    intent(in)    :: al_in, fr_in(6)    !< Cell and face fractions (mode 0).
   logical,      intent(in)    :: measures           !< Piece measures (mode 0).
   integer(I4P), intent(in)    :: mode               !< 0 end, 1 middle, 2 middle with the exact corner states.
   integer(I4P), intent(in)    :: npc                !< Piercings (`map_item_3d_count`).
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes (`map_item_3d_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   logical,      intent(out)   :: closed             !< Trusted.
   integer(I4P), intent(out)   :: np, nk             !< Pieces, features.
   integer(I4P), intent(in)    :: mk                 !< Feature capacity.
   integer(I4P), intent(inout) :: key(mk), kpiece(mk) !< Features, their pieces.
   integer(I4P), intent(in)    :: mm                 !< Volume capacity (0: none kept).
   real(R8P),    intent(inout) :: pvol(mm)           !< Piece volumes.
   real(R8P)                   :: b0(3), b1(3), h(3), ar(6), fa(6), fm(2, 6), fr(6), al !< Work.
   integer(I4P)                :: nrp, nip, npo, r0, i0, q, cp(8) !< Work.
   logical                     :: ok, cl, mp         !< Sections closed, faces closed, piece measures.

   call cut_cell_3d_one_sizes(n, npc, measures .and. mode == 0_I4P, nrp, nip, npo)
   r0 = 9_I4P * n + 1_I4P ; i0 = n + 1_I4P ; mp = measures .and. mode == 0_I4P
   call map_place_3d(n, lf, nf, t0, t1, tb, kind, rot, s, d, pr, pt, c0, c1, rw, b0, b1)
   do q = 1_I4P, 3_I4P
      h(q) = b1(q) - b0(q)
   enddo
   ar(1) = h(2) * h(3) ; ar(2) = h(2) * h(3) ; ar(3) = h(3) * h(1) ; ar(4) = h(3) * h(1) ; ar(5) = h(1) * h(2)
   ar(6) = h(1) * h(2)
   np = -1_I4P ; ok = .true.
   if (mode == 0_I4P) then
      do q = 1_I4P, 6_I4P
         fr(q) = fr_in(q)
      enddo
      al = al_in
   elseif (mode == 1_I4P) then
      do q = 1_I4P, 6_I4P
         fr(q) = -1._R8P
      enddo
      al = -1._R8P
   endif
   if (mode /= 2_I4P) call cut_cell_3d_one_pieces(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, mp, npc, nrp, &
                                                  rw(r0), nip, iw(i0), np, iw, cp, npo, rw(r0 + nrp), rw(r0 + nrp + npo), &
                                                  rw(r0 + nrp + 7_I4P * npo))
   if (np < 0_I4P .and. mode /= 0_I4P) then
      call map3_sections_ws(n, rw, b0, b1, cs, nrw - r0 + 1_I4P, rw(r0), niw - n, iw(i0), ok, fa, fm)
      if (ok) then
         do q = 1_I4P, 6_I4P
            fr(q) = fa(q) / ar(q)
         enddo
         call cut_cell_3d_one_totals(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), fr(2), fr(4), fr(6), niw - n, &
                                     iw(i0), al, cl)
         call cut_cell_3d_one_pieces(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, .false., npc, nrp, &
                                     rw(r0), nip, iw(i0), np, iw, cp, npo, rw(r0 + nrp), rw(r0 + nrp + npo), &
                                     rw(r0 + nrp + 7_I4P * npo))
      endif
   endif
   closed = np >= 0_I4P ; nk = 0_I4P
   if (.not. closed) then
      np = 0_I4P ; return
   endif
   do q = 8_I4P, 1_I4P, -1_I4P
      if (cp(q) == 0_I4P) cycle
      nk = nk + 1_I4P ; key(nk) = -q ; kpiece(nk) = cp(q)
   enddo
   do q = 1_I4P, n
      if (iw(q) == 0_I4P) cycle
      nk = nk + 1_I4P ; key(nk) = lf(q) ; kpiece(nk) = iw(q)
   enddo
   if (mp .and. np > 1_I4P) then
      do q = 1_I4P, min(np, mm)
         pvol(q) = rw(r0 + nrp + q - 1_I4P)
      enddo
   endif
   endsubroutine map_item_3d_eval

   subroutine map3_sections_ws(n, tl, b0, b1, cs, nrw, rw, niw, iw, ok, fa, fm)
   !< Kernel (issue #105 step 1d): the fluid areas of the six faces of the cell [b0, b1] cut by the facets tl, on a
   !< workspace (at least `cut_cell_3d_face_sizes(n, 4 n, n)`): the per-cell `map3_sections` of issue #102 bit for bit
   !< (each face's count pass, then its measures with the counts; a face section has at most n segments, so at most 4 n
   !< crossings and n wall pieces).
   !$acc routine seq
   integer(I4P), intent(in)    :: n                  !< Facets.
   real(R8P),    intent(in)    :: tl(3, 3, n)        !< Facets.
   real(R8P),    intent(in)    :: b0(3), b1(3)       !< Cell.
   integer(I4P), intent(in)    :: cs(8)              !< Corner states (0 solid, 1 fluid).
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes.
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   logical,      intent(out)   :: ok                 !< Every face closed.
   real(R8P),    intent(out)   :: fa(6)              !< Fluid areas (-x, +x, -y, +y, -z, +z).
   real(R8P),    intent(out)   :: fm(2, 6)           !< int u, int v over them (in-face coordinates in cyclic order).
   real(R8P)                   :: lev, al, cen(2), ar !< Work.
   integer(I4P)                :: f, a, bb, cc, q, bit(3), nv, nxs, nps, nr0, ni0, nr1, ni1 !< Work.
   logical                     :: fl(4), cl          !< Corner flags, closed.

   ok = .true.
   do f = 1_I4P, 6_I4P
      fa(f) = 0._R8P ; fm(1, f) = 0._R8P ; fm(2, f) = 0._R8P
   enddo
   call cut_cell_3d_face_sizes(n, 0_I4P, 0_I4P, nr0, ni0)
   do f = 1_I4P, 6_I4P
      a = (f + 1_I4P) / 2_I4P ; bb = mod(a, 3_I4P) + 1_I4P ; cc = mod(a + 1_I4P, 3_I4P) + 1_I4P
      if (mod(f, 2_I4P) == 0_I4P) then
         lev = b1(a) ; bit(a) = 1_I4P
      else
         lev = b0(a) ; bit(a) = 0_I4P
      endif
      do q = 1_I4P, 4_I4P
         bit(bb) = 0_I4P ; if (q == 2_I4P .or. q == 3_I4P) bit(bb) = 1_I4P
         bit(cc) = 0_I4P ; if (q >= 3_I4P) bit(cc) = 1_I4P
         fl(q) = cs(1_I4P + bit(1) + 2_I4P * bit(2) + 4_I4P * bit(3)) == 0_I4P
      enddo
      call cut_cell_3d_face_count(n, tl, a, lev, b0(bb), b1(bb), b0(cc), b1(cc), nr0, rw, ni0, iw, nv, nxs, nps)
      if (nxs > 4_I4P * n .or. nps > n) then ! cannot happen (see above)
         ok = .false. ; return
      endif
      call cut_cell_3d_face_sizes(n, nxs, nps, nr1, ni1)
      call cut_cell_3d_face_measures(n, tl, a, lev, b0(bb), b1(bb), b0(cc), b1(cc), fl(1), fl(2), fl(3), fl(4), nxs, nps, &
                                     nr1, rw, ni1, iw, cl, al, cen)
      if (.not. cl) then
         ok = .false. ; return
      endif
      ar = al * ((b1(bb) - b0(bb)) * (b1(cc) - b0(cc)))
      fa(f) = ar ; fm(1, f) = ar * cen(1) ; fm(2, f) = ar * cen(2)
   enddo
   endsubroutine map3_sections_ws

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

   subroutine map3_events(mc, c0, c1, lf, b, bw, nb, cert, unc, ct, cc)
   !< The events of the cell [c0, c1] in (0, 1), sorted, each once, from its facets lf: a mesh vertex crossing one of its
   !< six face planes inside the face, a mesh edge crossing one of its 12 edges, one of its 8 corners crossing a facet;
   !< each mesh vertex and edge once per cell (at its first facet vertex, edge, in the order of lf: a per-cell
   !< deduplication, `first_in_order`, issue #105 step 1c). Without rotation, explicit (the formulas of
   !< `plane_breakpoints` on the pose0 copy); with a rotation, `isolate_zeros` on the screw forms of `plane_events` (cert
   !< false if one is not certified; unc true if a kept event is not certified by a sign change: an uncertified contact).
   !< The position checks are inclusive by mc%tol: an extra event only adds an evaluation. bw (issue #104): each event's
   !< half-width (`EVENT_RES` on the explicit formulas, `isolate_zeros`' width on the screw ones; a duplicate keeps the
   !< larger). ct, cc: the corner events, instants and corners (as reals), for `map_corner_known`.
   type(map3_ctx_t),       intent(in)  :: mc           !< State.
   real(R8P),              intent(in)  :: c0(3), c1(3) !< Cell.
   integer(I4P),           intent(in)  :: lf(:)        !< Facets.
   real(R8P), allocatable, intent(out) :: b(:)         !< Events.
   real(R8P), allocatable, intent(out) :: bw(:)        !< Their half-widths.
   integer(I4P),           intent(out) :: nb           !< Their number.
   logical,                intent(out) :: cert         !< Every event certified.
   logical,                intent(out) :: unc          !< An event is an uncertified contact.
   real(R8P), allocatable, intent(out) :: ct(:)        !< Corner events: instants (in event order).
   real(R8P), allocatable, intent(out) :: cc(:)        !< Their corners (1..8, as reals).
   real(R8P),    allocatable :: tmp(:), tw(:)           !< Work.
   integer(I4P)              :: ids(3 * size(lf))       !< Mesh vertex or edge of every facet vertex (q, m), visit order.
   logical                   :: fv1(3 * size(lf))       !< First visit of the mesh vertex of facet vertex (q, m).
   logical                   :: fe1(3 * size(lf))       !< First visit of the mesh edge of facet edge (q, m).
   real(R8P)                 :: p(3), lo(3), hi(3)      !< Corner, facet box.
   integer(I4P)              :: q, f, m, cn, nb0, k, nc, ncn, n !< Work.

   cert = .true. ; unc = .false. ; nb = 0_I4P
   n = size(lf, kind=I4P)
   allocate(tmp(16), tw(16), ct(8), cc(8)) ; ncn = 0_I4P
   do q = 1_I4P, n
      ids(3_I4P * (q - 1_I4P) + 1_I4P:3_I4P * q) = mc%fv(:, lf(q))
   enddo
   call first_in_order(ids, fv1)
   do q = 1_I4P, n
      ids(3_I4P * (q - 1_I4P) + 1_I4P:3_I4P * q) = mc%fe(:, lf(q))
   enddo
   call first_in_order(ids, fe1)
   do q = 1_I4P, n
      f = lf(q)
      lo = min(mc%vlo(:, 1, f), mc%vlo(:, 2, f), mc%vlo(:, 3, f))
      hi = max(mc%vhi(:, 1, f), mc%vhi(:, 2, f), mc%vhi(:, 3, f))
      do cn = 1_I4P, 8_I4P
         p = corner_point(c0, c1, cn)
         if (any(p < lo) .or. any(p > hi)) cycle
         nb0 = nb
         call event_corner(mc, f, p, tmp, tw, nb, cert, unc)
         do k = nb0 + 1_I4P, nb ! the corner's own list, for its states within the step (`map_corner_known`)
            nc = ncn
            call push_value(ct, nc, tmp(k))
            call push_value(cc, ncn, real(cn, R8P))
         enddo
      enddo
      do m = 1_I4P, 3_I4P
         if (fv1(3_I4P * (q - 1_I4P) + m)) call event_vertex(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
         if (fe1(3_I4P * (q - 1_I4P) + m)) call event_edge(mc, f, m, c0, c1, tmp, tw, nb, cert, unc)
      enddo
   enddo
   call sort_unique_w(tmp, tw, nb)
   b = tmp(1:nb) ; bw = tw(1:nb)
   tmp = ct(1:ncn) ; call move_alloc(tmp, ct)
   tw = cc(1:ncn) ; call move_alloc(tw, cc)
   endsubroutine map3_events

   pure subroutine first_in_order(ids, first)
   !< first(r): ids(r) differs from every ids(r') with r' < r (the positions sorted by id with a stable insertion sort,
   !< then the head of each run of equal ids: its smallest position). A few dozen ids per cell.
   integer(I4P), intent(in)  :: ids(:)   !< Ids in visit order [n].
   logical,      intent(out) :: first(:) !< First occurrences [n].
   integer(I4P)              :: ord(size(ids)) !< Positions sorted by id.
   integer(I4P)              :: r, h, a  !< Work.

   do r = 1_I4P, size(ids, kind=I4P)
      a = r ; h = r - 1_I4P
      do while (h >= 1_I4P)
         if (ids(ord(h)) <= ids(a)) exit
         ord(h + 1_I4P) = ord(h) ; h = h - 1_I4P
      enddo
      ord(h + 1_I4P) = a
   enddo
   do r = 1_I4P, size(ids, kind=I4P)
      if (r == 1_I4P) then
         first(ord(r)) = .true.
      else
         first(ord(r)) = ids(ord(r)) /= ids(ord(r - 1_I4P))
      endif
   enddo
   endsubroutine first_in_order

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

   pure subroutine map_corner_known(nq, ct, cc, cs0, cs1, known)
   !< Which corners of a cell have states that follow from their events (the flip rule): the corner's events (it
   !< crossing a facet, `map3_events`, or a wall segment, `map_events`: instants ct, corners cc) at least 1e-9 apart, so
   !< that each is one transversal crossing, and their parity taking the pose0 state cs0 to the pose1 state cs1. A corner
   !< passing through a mesh edge or vertex (two facets or segments at once), touching one, or whose count disagrees
   !< takes the exact query (`map3_corner_in`, `map_corner_in`) instead. 3-D: 8 corners; 2-D (issue #105 step 1d): 4.
   integer(I4P), intent(in)  :: nq        !< Corners.
   real(R8P),    intent(in)  :: ct(:)     !< The cell's corner events: instants.
   integer(I4P), intent(in)  :: cc(:)     !< Their corners.
   integer(I4P), intent(in)  :: cs0(nq)   !< Corner states at pose0 (0 solid, 1 fluid).
   integer(I4P), intent(in)  :: cs1(nq)   !< Corner states at pose1.
   logical,      intent(out) :: known(nq) !< The corner's states follow from its events.
   real(R8P)                 :: t(64)     !< Instants of one corner.
   integer(I4P)              :: q, k, n   !< Work.

   do q = 1_I4P, nq
      n = 0_I4P ; known(q) = .true.
      do k = 1_I4P, size(ct, kind=I4P)
         if (cc(k) /= q) cycle
         if (n == size(t, kind=I4P)) then
            known(q) = .false. ; exit
         endif
         n = n + 1_I4P ; t(n) = ct(k)
      enddo
      if (.not. known(q)) cycle
      call sort_unique(t, n)
      do k = 2_I4P, n
         if (t(k) - t(k - 1_I4P) <= 1e-9_R8P) known(q) = .false.
      enddo
      if (count(cc == q) /= n) known(q) = .false. ! exact repeats merged by sort_unique
      if (mod(cs0(q) + n, 2_I4P) /= cs1(q)) known(q) = .false.
   enddo
   endsubroutine map_corner_known

   function map3_corner_state(mc, surface, s0, s, ps, c0, c1, b0, b1, q, known, cs0, ct, cc) result(st)
   !< The state (0 solid, 1 fluid) of corner q of the cell [c0, c1] at instant s: by the flip rule from its events (ct,
   !< cc) and its pose0 state cs0 when known (`map_corner_known`), else by the exact query `map3_corner_in` (the cell
   !< [b0, b1] shifted back on the pose0 copy without rotation, the corner taken to the body frame by the screw pose ps
   !< at s with one).
   type(map3_ctx_t),         intent(in) :: mc           !< State.
   type(surface_stl_object), intent(in) :: surface, s0  !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in) :: s            !< Instant.
   type(rigid_pose_t),       intent(in) :: ps           !< Screw pose at s (with a rotation).
   real(R8P),                intent(in) :: c0(3), c1(3) !< Cell.
   real(R8P),                intent(in) :: b0(3), b1(3) !< Cell shifted back (no rotation).
   integer(I4P),             intent(in) :: q            !< Corner.
   logical,                  intent(in) :: known        !< Its states follow from its events.
   integer(I4P),             intent(in) :: cs0(8)       !< Corner states at pose0.
   real(R8P),                intent(in) :: ct(:)        !< The cell's corner events: instants.
   integer(I4P),             intent(in) :: cc(:)        !< Their corners.
   integer(I4P)                         :: st           !< State.
   integer(I4P)                         :: k, n         !< Work.

   if (known) then
      n = 0_I4P
      do k = 1_I4P, size(ct, kind=I4P)
         if (cc(k) == q .and. ct(k) < s) n = n + 1_I4P
      enddo
      st = mod(cs0(q) + n, 2_I4P)
   elseif (mc%rot) then
      st = merge(0_I4P, 1_I4P, map3_corner_in(mc, surface, s0, to_body(ps, corner_point(c0, c1, q))))
   else
      st = merge(0_I4P, 1_I4P, map3_corner_in(mc, surface, s0, corner_point(b0, b1, q)))
   endif
   endfunction map3_corner_state

   subroutine map_measures_3d(mc, plan, ep, comp, surface, s0, x, y, z, k0, k1, ntr, moff, meas, cok)
   !< Per-track time integrals over the step of the fluid areas of the cell's six faces (meas(1:6)) and of the outward
   !< rigid-field outflows through them (meas(7:12)), for every measured cell (moff), interval by interval between the
   !< cell's events, on every path: without rotation two-point Gauss per interval (exact: the areas are quadratic in s),
   !< with one Gauss doubled from 2 nodes per interval until every track's values agree to round-off (the criteria of
   !< `one_plane`; mc%exact false if not at 64). A track exists only between events, so its integrand jumps there: the
   !< step's nq Gauss nodes of the Gauss path (`exact=.false.`) would give a short-lived track no node, or one of a large
   !< weight (2-D, NACA 0012, nq = 4: a per-track residual of 0.39 of a cell where the cell's is 1e-6). On the Gauss path
   !< the sums over a cell's tracks therefore differ from the cell's Gauss values in `sw` by the quadrature error of
   !< these (the per-track values are the better ones). cok false if a piece at a node cannot be assigned.
   !< By rounds over flat items (issue #105 step 1e), as `map_measures_2d`; the nodes by `map_nodes_3d`.
   type(map3_ctx_t),         intent(inout) :: mc               !< State (skipped, exact).
   type(map_plan_t),         intent(in)    :: plan             !< Host phase.
   type(map_evals_t),        intent(in)    :: ep               !< Evaluations.
   integer(I4P),             intent(in)    :: comp(:)          !< Track of every piece of every evaluation.
   type(surface_stl_object), intent(in)    :: surface, s0      !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(in)    :: k0, k1           !< End-pose cuts (corner states).
   integer(I4P),             intent(in)    :: ntr(:)           !< Tracks of every cell.
   integer(I4P),             intent(in)    :: moff(:)          !< First column of every measured cell's measures (0: not).
   real(R8P),                intent(inout) :: meas(:,:)        !< Measures [12, *] (zero on input).
   logical,                  intent(out)   :: cok(:)           !< Measured cells ok [nx ny nz].
   type(map_pairs_t)         :: pr                             !< Pairs.
   real(R8P),    allocatable :: cur(:,:), prv(:,:)             !< Sums of the pairs (rotation).
   real(R8P),    allocatable :: ksc(:,:), kso(:,:), ksa(:), ksn(:) !< Scales of the test per measured cell.
   integer(I4P), allocatable :: kcs0(:,:)                      !< Corner states at pose0 per measured cell.
   logical,      allocatable :: kknown(:,:)                    !< Corners known by the flip rule per measured cell.
   real(R8P),    allocatable :: ns(:), nw(:)                   !< Nodes of a round: instants, weights.
   integer(I4P), allocatable :: nj(:)                          !< Nodes of a round: pairs.
   real(R8P)                 :: vel, hc(3)                     !< Speed bound, cell.
   integer(I4P)              :: nx, ny, k, c, i, j, kk, n0, n1, round, ng, nn, ncol, cs1(8) !< Work.
   logical                   :: rot                            !< A rotation.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; rot = mc%rot ; cok = .true.
   call map_pairs_build(plan, ntr, moff, rot, pr, ncol)
   allocate(cur(12, ncol), prv(12, ncol), ksc(6, pr%nk), kso(6, pr%nk), ksa(pr%nk), ksn(pr%nk), kcs0(8, pr%nk), &
            kknown(8, pr%nk))
   if (rot) vel = norm2(mc%om) * mc%xm + norm2(mc%d)
   do k = 1_I4P, pr%nk
      c = pr%cell(k) ; call cell_ijk(c, nx, ny, i, j, kk)
      n0 = plan%cnstart(c) ; n1 = plan%cnstart(c + 1) - 1_I4P
      call node_states(k0%node_in, i, j, kk, kcs0(:, k)) ; call node_states(k1%node_in, i, j, kk, cs1)
      call map_corner_known(8_I4P, plan%cntime(n0:n1), plan%cncor(n0:n1), kcs0(:, k), cs1, kknown(:, k))
      if (rot) then
         hc = [x(i + 1), y(j + 1), z(kk + 1)] - [x(i), y(j), z(kk)]
         ksc(:, k) = [hc(2) * hc(3), hc(2) * hc(3), hc(3) * hc(1), hc(3) * hc(1), hc(1) * hc(2), hc(1) * hc(2)]
         kso(:, k) = ksc(:, k) * vel ; ksa(k) = mc%xm * maxval(hc) ; ksn(k) = ksa(k) * vel
      endif
   enddo
   do round = 1_I4P, merge(6_I4P, 1_I4P, rot)
      ng = 2_I4P**round
      call map_round_nodes(pr, rot, ng, plan%gq_s, plan%gq_w, nn, nj, ns, nw)
      if (nn == 0_I4P) exit
      if (rot) then
         call map_round_zero(pr, ntr, cur)
         call map_nodes_3d(mc, plan, ep, comp, surface, s0, x, y, z, ntr, kcs0, kknown, pr, nn, nj, ns, nw, cur)
         call map_round_close(pr, ntr, round, ng, 6_I4P, cur, prv, ksc, kso, ksa, ksn)
      else
         call map_nodes_3d(mc, plan, ep, comp, surface, s0, x, y, z, ntr, kcs0, kknown, pr, nn, nj, ns, nw, meas)
      endif
   enddo
   call map_measures_finish(pr, rot, ntr, moff, cur, meas, cok, mc%skipped, mc%exact)
   endsubroutine map_measures_3d

   subroutine map_nodes_3d(mc, plan, ep, comp, surface, s0, x, y, z, ntr, kcs0, kknown, pr, nn, nj, ns, nw, acc)
   !< One round of 3-D nodes (issue #105 step 1e). Host: the poses at the round's distinct instants with a rotation
   !< (`pose_interpolate`, `node_instants`), the corner states of every node (`map3_corner_state`: the flip rule where
   !< it holds, else the exact query). Then the sizes pass (`map_item_3d_count`) and, in chunks of at most
   !< MAP_CHUNK_BYTES of workspace, the nodes (`map_node_3d_eval`); a node whose pieces failed with a corner of its cell
   !< known by the flip rule is evaluated again with the exact corner states (the per-cell `map3_eval`'s retry); then
   !< the chunk's sums into their pairs' columns of acc, node after node (`map_node_outcome`, `map_node_reduce`).
   type(map3_ctx_t),         intent(in)    :: mc               !< State.
   type(map_plan_t),         intent(in)    :: plan             !< Host phase.
   type(map_evals_t),        intent(in)    :: ep               !< Evaluations.
   integer(I4P),             intent(in)    :: comp(:)          !< Track of every piece of every evaluation.
   type(surface_stl_object), intent(in)    :: surface, s0      !< Surface (body frame), copy at pose0.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),             intent(in)    :: ntr(:)           !< Tracks of every cell.
   integer(I4P),             intent(in)    :: kcs0(:,:)        !< Corner states at pose0 per measured cell.
   logical,                  intent(in)    :: kknown(:,:)      !< Corners known by the flip rule per measured cell.
   type(map_pairs_t),        intent(inout) :: pr               !< Pairs.
   integer(I4P),             intent(in)    :: nn               !< Nodes.
   integer(I4P),             intent(in)    :: nj(:)            !< Pair of each node.
   real(R8P),                intent(in)    :: ns(:), nw(:)     !< Instants, weights.
   real(R8P),                intent(inout) :: acc(:,:)         !< Sums [12, *].
   integer(I4P), allocatable :: cs(:,:), pid(:), npc(:), wr(:), wi(:), ro(:), io(:), np(:), nsk(:), iw(:), rl(:) !< Work.
   logical,      allocatable :: nf(:), rt(:)                   !< Failed, to retry.
   real(R8P),    allocatable :: us(:), qr(:,:,:), qt(:,:), rw(:) !< Pose table, workspace.
   type(rigid_pose_t)        :: ps                             !< Pose.
   real(R8P)                 :: prr(3, 3), ptt(3), c0(3), c1(3), b0(3), b1(3), tq(3) !< Work.
   integer(I4P)              :: nx, ny, v, jp, k, c, i, j, kk, l0, nl, n0, n1, q, a, b, h, m, pass, npose, nr, ni, mk, npo, &
                                nrt !< Work.
   integer(I4P)              :: csz(8)                         !< Unused pose0 states (exact queries).
   logical                   :: ok                             !< Node added.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; csz = 0_I4P
   allocate(cs(8, nn), pid(nn), npc(nn), wr(nn), wi(nn), ro(nn + 1), io(nn + 1), np(nn), nsk(nn), nf(nn), rt(nn), rl(nn))
   pid = 0_I4P ; npose = 0_I4P
   if (mc%rot) then ! the screw poses at the round's instants
      call node_instants(nn, ns, pid, us, npose)
      allocate(qr(3, 3, npose), qt(3, npose))
      do k = 1_I4P, npose
         ps = pose_interpolate(mc%pose0, mc%pose1, us(k)) ; qr(:, :, k) = ps%r ; qt(:, k) = ps%t
      enddo
   else
      allocate(qr(3, 3, 0), qt(3, 0))
   endif
   do v = 1_I4P, nn ! corner states (host: the flip rule or the exact queries)
      jp = nj(v) ; k = pr%k(jp) ; c = pr%cell(k) ; call cell_ijk(c, nx, ny, i, j, kk)
      c0 = [x(i), y(j), z(kk)] ; c1 = [x(i + 1), y(j + 1), z(kk + 1)] ; b0 = c0 ; b1 = c1
      if (.not. mc%rot) then
         b0 = c0 - ns(v) * mc%d ; b1 = c1 - ns(v) * mc%d
      endif
      if (pid(v) > 0_I4P) then
         ps%r = qr(:, :, pid(v)) ; ps%t = qt(:, pid(v))
      endif
      n0 = plan%cnstart(c) ; n1 = plan%cnstart(c + 1) - 1_I4P
      do q = 1_I4P, 8_I4P
         cs(q, v) = map3_corner_state(mc, surface, s0, ns(v), ps, c0, c1, b0, b1, q, kknown(q, k), kcs0(:, k), &
                                      plan%cntime(n0:n1), plan%cncor(n0:n1))
      enddo
      nl = plan%lstart(c + 1) - plan%lstart(c) ; wr(v) = 9_I4P * nl ; wi(v) = nl
   enddo
   ! sizes pass: the placed facets, the piercings
   a = 1_I4P
   do while (a <= nn)
      call map_chunk(a, nn, wr, wi, b, ro, io, rw, iw)
      !$omp parallel do default(shared) private(v, c, i, j, kk, l0, nl, prr, ptt, c0, c1) schedule(dynamic, 16)
      do v = a, b
         c = pr%cell(pr%k(nj(v))) ; call cell_ijk(c, nx, ny, i, j, kk)
         c0 = [x(i), y(j), z(kk)] ; c1 = [x(i + 1), y(j + 1), z(kk + 1)]
         l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
         call node_pose_rt(qr, qt, pid(v), 3_I4P, prr, ptt)
         call map_item_3d_count(nl, plan%lel(l0), mc%nf, mc%t0, mc%t1, mc%tb, 2_I4P, mc%rot, ns(v), mc%d, prr, ptt, c0, c1, &
                                wr(v), rw(ro(v)), wi(v), iw(io(v)), npc(v))
      enddo
      !$omp end parallel do
      a = b + 1_I4P
   enddo
   do v = 1_I4P, nn
      c = pr%cell(pr%k(nj(v)))
      call map_node_3d_sizes(plan%lstart(c + 1) - plan%lstart(c), npc(v), wr(v), wi(v), npo)
   enddo
   ! the nodes, chunk by chunk: evaluated, the failures with a known corner again, then the chunk's sums in node order
   a = 1_I4P
   do while (a <= nn)
      call map_chunk(a, nn, wr, wi, b, ro, io, rw, iw)
      do v = a, b
         rl(v - a + 1_I4P) = v
      enddo
      nrt = b - a + 1_I4P
      do pass = 1_I4P, 2_I4P
         !$omp parallel do default(shared) private(h, v, jp, c, i, j, kk, l0, nl, prr, ptt, c0, c1, tq) schedule(dynamic, 16)
         do h = 1_I4P, nrt
            v = rl(h) ; jp = nj(v) ; c = pr%cell(pr%k(jp)) ; call cell_ijk(c, nx, ny, i, j, kk)
            c0 = [x(i), y(j), z(kk)] ; c1 = [x(i + 1), y(j + 1), z(kk + 1)]
            l0 = plan%lstart(c) ; nl = plan%lstart(c + 1) - l0
            call node_pose_rt(qr, qt, pid(v), 3_I4P, prr, ptt)
            tq = mc%pose0%t + ns(v) * mc%d
            call map_node_3d_eval(nl, plan%lel(l0), mc%nf, mc%t0, mc%t1, mc%tb, mc%rot, ns(v), mc%d, prr, ptt, c0, c1, &
                                  cs(:, v), npc(v), nw(v), tq, mc%om, pr%narrow(jp), ep%nk(pr%vr(jp)), &
                                  ep%key(ep%kstart(pr%vr(jp)):ep%kstart(pr%vr(jp) + 1) - 1_I4P), &
                                  ep%kpiece(ep%kstart(pr%vr(jp)):ep%kstart(pr%vr(jp) + 1) - 1_I4P), ep%np(pr%vr(jp)), &
                                  comp(ep%pstart(pr%vr(jp)):ep%pstart(pr%vr(jp) + 1) - 1_I4P), ep%nk(pr%va(jp)), &
                                  ep%key(ep%kstart(pr%va(jp)):ep%kstart(pr%va(jp) + 1) - 1_I4P), &
                                  ep%kpiece(ep%kstart(pr%va(jp)):ep%kstart(pr%va(jp) + 1) - 1_I4P), ep%np(pr%va(jp)), &
                                  comp(ep%pstart(pr%va(jp)):ep%pstart(pr%va(jp) + 1) - 1_I4P), wr(v), rw(ro(v)), wi(v), &
                                  iw(io(v)), np(v), nf(v), nsk(v), rt(v))
         enddo
         !$omp end parallel do
         if (pass == 2_I4P) exit
         m = nrt ; nrt = 0_I4P
         do h = 1_I4P, m ! the pieces failed with a known corner: once more, with the exact corner states (host queries)
            v = rl(h) ; k = pr%k(nj(v))
            if (.not. (rt(v) .and. any(kknown(:, k)))) cycle
            c = pr%cell(k) ; call cell_ijk(c, nx, ny, i, j, kk)
            c0 = [x(i), y(j), z(kk)] ; c1 = [x(i + 1), y(j + 1), z(kk + 1)] ; b0 = c0 ; b1 = c1
            if (.not. mc%rot) then
               b0 = c0 - ns(v) * mc%d ; b1 = c1 - ns(v) * mc%d
            endif
            if (pid(v) > 0_I4P) then
               ps%r = qr(:, :, pid(v)) ; ps%t = qt(:, pid(v))
            endif
            do q = 1_I4P, 8_I4P
               cs(q, v) = map3_corner_state(mc, surface, s0, ns(v), ps, c0, c1, b0, b1, q, .false., csz, plan%cntime(1:0), &
                                            plan%cncor(1:0))
            enddo
            nrt = nrt + 1_I4P ; rl(nrt) = v
         enddo
         if (nrt == 0_I4P) exit
      enddo
      do v = a, b
         jp = nj(v)
         call map_node_outcome(pr, jp, nf(v), nsk(v), ok)
         if (.not. ok) cycle
         c = pr%cell(pr%k(jp))
         call map_item_3d_sizes(plan%lstart(c + 1) - plan%lstart(c), npc(v), 2_I4P, .true., nr, ni, mk, npo)
         call map_node_reduce(12_I4P, np(v), iw(io(v) + ni + 2_I4P * mk), rw(ro(v) + nr), ntr(c), &
                              acc(:, pr%col(jp):pr%col(jp) + ntr(c) - 1_I4P))
      enddo
      a = b + 1_I4P
   enddo
   endsubroutine map_nodes_3d

   pure subroutine map_node_3d_sizes(n, npc, nrw, niw, npo)
   !< Workspace of a 3-D node (`map_node_3d_eval`, issue #105 step 1e): a middle evaluation's with measures
   !< (`map_item_3d_sizes`: nr reals, ni integers, mk features, npo pieces), then its outputs: reals, the 12 values of
   !< every piece (12 npo); integers, the features and their pieces (2 mk), then the track of every piece (npo).
   !$acc routine seq
   integer(I4P), intent(in)  :: n, npc        !< Facets, piercings.
   integer(I4P), intent(out) :: nrw, niw, npo !< Reals, integers, piece capacity.
   integer(I4P)              :: nr, ni, mk    !< Evaluation workspace, feature capacity.

   call map_item_3d_sizes(n, npc, 2_I4P, .true., nr, ni, mk, npo)
   nrw = nr + 12_I4P * npo ; niw = ni + 2_I4P * mk + npo
   endsubroutine map_node_3d_sizes

   subroutine map_node_3d_eval(n, lf, nf, t0, t1, tb, rot, s, d, pr, pt, c0, c1, cs, npc, w, tq, om, narrow, nkr, rkey, &
                               rkpiece, nrp, rcomp, nka, akey, akpiece, npa, acomp, nrw, rw, niw, iw, np, fail, nskip, retry)
   !< Kernel (issue #105 step 1e): one quadrature node of the 3-D measures (instant s, weight w), the per-cell
   !< `map3_node` of issue #102 without the exact corner queries (the corner states cs come with the node): the facets
   !< placed at s (`map_place_3d`), the six face sections with their first moments (`map3_sections_ws`), alpha from them
   !< (`cut_cell_3d_one_totals`), the pieces with their measures (`cut_cell_3d_one_pieces`; one piece takes the
   !< sections'), then per piece its track (`map_piece_assign`) and its 12 values: w times the fluid areas of the faces
   !< (-x, +x, -y, +y, -z, +z) and w times the outward rigid-field outflows through them, of a face normal to a ((a, b,
   !< c) cyclic) the field om x (x - tq) + d: int (om_b (x_c - tq_c) - om_c (x_b - tq_b) + d_a) dA from the area and the
   !< first moments, d_a A without rotation. Outputs after the evaluation's workspace (`map_node_3d_sizes`): the values
   !< pv(12, npo) and the track of every piece (0: none). fail, nskip as `map_node_2d_eval`; retry: the sections closed
   !< but the pieces failed (the per-cell `map3_eval` retried those with the exact corner states).
   !$acc routine seq
   integer(I4P), intent(in)    :: n                  !< Facets.
   integer(I4P), intent(in)    :: lf(n)              !< Their ids.
   integer(I4P), intent(in)    :: nf                 !< Surface facets.
   real(R8P),    intent(in)    :: t0(3, 3, nf), t1(3, 3, nf), tb(3, 3, nf) !< Pose0 copy, pose1 copy, body frame.
   logical,      intent(in)    :: rot                !< The step rotates.
   real(R8P),    intent(in)    :: s, d(3), pr(3, 3), pt(3) !< Instant, translation over the step, pose at s.
   real(R8P),    intent(in)    :: c0(3), c1(3)       !< Cell.
   integer(I4P), intent(in)    :: cs(8)              !< Corner states (0 solid, 1 fluid).
   integer(I4P), intent(in)    :: npc                !< Piercings (`map_item_3d_count`).
   real(R8P),    intent(in)    :: w, tq(3), om(3)    !< Weight, translation at s, rotation vector of the step.
   logical,      intent(in)    :: narrow             !< A narrow interval.
   integer(I4P), intent(in)    :: nkr                !< Features of r.
   integer(I4P), intent(in)    :: rkey(nkr), rkpiece(nkr) !< Its features, their pieces.
   integer(I4P), intent(in)    :: nrp                !< Pieces of r.
   integer(I4P), intent(in)    :: rcomp(nrp)         !< Their tracks.
   integer(I4P), intent(in)    :: nka                !< Features of a.
   integer(I4P), intent(in)    :: akey(nka), akpiece(nka) !< Its features, their pieces.
   integer(I4P), intent(in)    :: npa                !< Pieces of a.
   integer(I4P), intent(in)    :: acomp(npa)         !< Their tracks.
   integer(I4P), intent(in)    :: nrw, niw           !< Workspace sizes (`map_node_3d_sizes`).
   real(R8P),    intent(inout) :: rw(nrw)            !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)            !< Integer workspace.
   integer(I4P), intent(out)   :: np                 !< Pieces.
   logical,      intent(out)   :: fail               !< The node cannot be assigned.
   integer(I4P), intent(out)   :: nskip              !< Skipped (narrow).
   logical,      intent(out)   :: retry              !< The pieces failed after the sections closed.
   real(R8P)    :: b0(3), b1(3), h(3), ar(6), fa(6), fm(2, 6), fr(6), al, xc(3), pa(6), pm(2, 6), big, hm, o !< Work.
   integer(I4P) :: nr, ni, mk, npo, nrp1, nip, npo1, r0, i0, q, cp(8), nk, p, f, ax, bb, ec, t, act, npx, rv, rf, rm !< Work.
   logical      :: ok, cl                            !< Sections closed, faces closed.

   call map_item_3d_sizes(n, npc, 2_I4P, .true., nr, ni, mk, npo)
   call cut_cell_3d_one_sizes(n, npc, .true., nrp1, nip, npo1)
   r0 = 9_I4P * n + 1_I4P ; i0 = n + 1_I4P
   fail = .false. ; nskip = 0_I4P ; retry = .false. ; np = 0_I4P ; npx = -1_I4P
   call map_place_3d(n, lf, nf, t0, t1, tb, 2_I4P, rot, s, d, pr, pt, c0, c1, rw, b0, b1)
   do q = 1_I4P, 3_I4P
      h(q) = b1(q) - b0(q)
   enddo
   ar(1) = h(2) * h(3) ; ar(2) = h(2) * h(3) ; ar(3) = h(3) * h(1) ; ar(4) = h(3) * h(1) ; ar(5) = h(1) * h(2)
   ar(6) = h(1) * h(2)
   call map3_sections_ws(n, rw, b0, b1, cs, nr - r0 + 1_I4P, rw(r0), ni - n, iw(i0), ok, fa, fm)
   if (ok) then
      do q = 1_I4P, 6_I4P
         fr(q) = fa(q) / ar(q)
      enddo
      call cut_cell_3d_one_totals(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), fr(2), fr(4), fr(6), ni - n, iw(i0), al, cl)
      call cut_cell_3d_one_pieces(n, rw, b0(1), b1(1), b0(2), b1(2), b0(3), b1(3), cs, al, fr, .true., npc, nrp1, rw(r0), &
                                  nip, iw(i0), npx, iw, cp, npo1, rw(r0 + nrp1), rw(r0 + nrp1 + npo1), &
                                  rw(r0 + nrp1 + 7_I4P * npo1))
      retry = npx < 0_I4P
   endif
   if (npx < 0_I4P) then
      if (narrow) then
         nskip = 1_I4P
      else
         fail = .true.
      endif
      return
   endif
   np = npx ; nk = 0_I4P
   do q = 8_I4P, 1_I4P, -1_I4P
      if (cp(q) == 0_I4P) cycle
      nk = nk + 1_I4P ; iw(ni + nk) = -q ; iw(ni + mk + nk) = cp(q)
   enddo
   do q = 1_I4P, n
      if (iw(q) == 0_I4P) cycle
      nk = nk + 1_I4P ; iw(ni + nk) = lf(q) ; iw(ni + mk + nk) = iw(q)
   enddo
   hm = c1(1) - c0(1)
   do q = 2_I4P, 3_I4P
      if (c1(q) - c0(q) > hm) hm = c1(q) - c0(q)
   enddo
   do q = 1_I4P, 3_I4P
      xc(q) = 0.5_R8P * (b0(q) + b1(q))
   enddo
   rv = r0 + nrp1 ; rf = rv + npo1 ; rm = rf + 6_I4P * npo1 ! piece volumes, face fractions, face moments
   do p = 1_I4P, np
      iw(ni + 2_I4P * mk + p) = 0_I4P
      do f = 1_I4P, 6_I4P
         if (np == 1_I4P) then
            pa(f) = fa(f) ; pm(1, f) = fm(1, f) ; pm(2, f) = fm(2, f)
         else
            ax = (f + 1_I4P) / 2_I4P ; bb = mod(ax, 3_I4P) + 1_I4P ; ec = mod(ax + 1_I4P, 3_I4P) + 1_I4P
            pa(f) = rw(rf + 6_I4P * (p - 1_I4P) + f - 1_I4P) * ar(f)
            pm(1, f) = xc(bb) * pa(f) + h(bb) * rw(rm + 18_I4P * (p - 1_I4P) + 3_I4P * (f - 1_I4P) + 1_I4P)
            pm(2, f) = xc(ec) * pa(f) + h(ec) * rw(rm + 18_I4P * (p - 1_I4P) + 3_I4P * (f - 1_I4P) + 2_I4P)
         endif
      enddo
      big = pa(1)
      do f = 2_I4P, 6_I4P
         if (pa(f) > big) big = pa(f)
      enddo
      call map_piece_assign(nk, iw(ni + 1_I4P), iw(ni + mk + 1_I4P), p, narrow, big <= 1e-9_R8P * hm**2, nkr, rkey, rkpiece, &
                            nrp, rcomp, nka, akey, akpiece, npa, acomp, t, act)
      if (act == 3_I4P) then
         fail = .true. ; return
      elseif (act == 2_I4P) then
         nskip = nskip + 1_I4P ; cycle
      elseif (act == 0_I4P) then
         cycle
      endif
      iw(ni + 2_I4P * mk + p) = t
      do f = 1_I4P, 6_I4P
         ax = (f + 1_I4P) / 2_I4P ; bb = mod(ax, 3_I4P) + 1_I4P ; ec = mod(ax + 1_I4P, 3_I4P) + 1_I4P
         if (rot) then
            o = om(bb) * (pm(2, f) - tq(ec) * pa(f)) - om(ec) * (pm(1, f) - tq(bb) * pa(f)) + d(ax) * pa(f)
         else
            o = d(ax) * pa(f)
         endif
         if (mod(f, 2_I4P) == 1_I4P) o = -o
         rw(nr + 12_I4P * (p - 1_I4P) + f) = w * pa(f) ; rw(nr + 12_I4P * (p - 1_I4P) + 6_I4P + f) = w * o
      enddo
   enddo
   endsubroutine map_node_3d_eval

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
