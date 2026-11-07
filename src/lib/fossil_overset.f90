!< FOSSIL, overset assembly kernels on curvilinear blocks (issue #72).
module fossil_overset
!< FOSSIL, overset assembly kernels on curvilinear blocks (issue #72).
!<
!< Design: docs/guide/advanced/overset.md. FOSSIL gives kernels and measures; the assembly policy stays with the
!< caller. Points are the vertices of a set of blocks (cell-centred data: the blocks of the cell centres).
!<
!< Cutters (`wall_cutter`, `wall_cutter_2d`): the wall faces of the blocks (faces covered by TOPO_BOUNDARY records with
!< the given name) as a surface (curve) oriented with the block interior, the fluid, on the outer side. A bilinear
!< face is split along its shorter diagonal; `gap` reports the largest distance between a face's bilinear centre and
!< its two triangles. Overlapping wall patches of two blocks are trimmed: the elements of the coarser patch (larger
!< mean element size over the overlap) whose centroid lies within the tolerance of the finer patch are dropped, so the
!< elements straddling the edge of the finer patch split between kept and dropped: a residual strip of over- and
!< under-coverage about one coarse element wide, whose effects on the winding number partly cancel.
!< Two blocks whose patches overlap share a wall (`shared`).
!<
!< Hole flags (`hole_flags`, `hole_flags_2d`): a vertex is a hole when the signed distance to the cutter, signed by the
!< fast winding number (negative inside), is below `offset` (default 0: inside; a vertex on the wall is a field
!< point). `cuts(b)` selects the blocks the cutter cuts (default all). device=.true. runs the same kernels on the GPU
!< (the cutter is entered and left around the call if it is not on the device already), bit-identical.
!<
!< Hole map (`hole_map`, Chiu & Meakin / PEGASUS 5): a Cartesian map of n elements per direction over the cutter's
!< box enlarged by one element; the elements the wall meets are cut (exact BOX_CUT), the others are painted outside
!< from the corner elements through face neighbours, the rest inside. A vertex in an uncut element takes its status;
!< in a cut element it looks along a grid axis (+x, -x, +y, ... in turn) to the first uncut element of the column
!< and takes that element's status, flipped by the parity of the exact wall crossings in between (`line_crossings`,
!< #40): exact. pegasus=.true. instead takes the status of a face-neighbour element only if no crossing lies in
!< between (the PEGASUS 5 line of sight), and leaves the vertex a field point otherwise (`unresolved` counts them):
!< inside a body thinner than the map no such neighbour exists and the hole is not cut.
!< X-rays (`xray_flags`, Meakin): the exact crossings of the wall along the lines parallel to z through a raster of
!< n(1) x n(2) points over the cutter's box; a vertex takes the side of its nearest ray at its own z (inside if the
!< last crossing below it enters the solid): exact on the rays, wrong only near the wall between rays.
!<
!< Donor candidates (`donor_candidates`): every vertex of every block located (WP2's robust path; build_search
!< first) in every other block the pair matrix allows; every candidate other than LOC_OUTSIDE is kept, with its
!< quality (donor volume, receptor volume = mean volume of the cells around the vertex, their ratio, the certified min
!< det J of the donor over its volume). device=.true. runs the location on the GPU (blocks entered once,
!< enter_device), bit-identical. The pair matrix (`overlap_matrix`): PAIR_OVERLAP if a boundary vertex of one block
!< lies in a cell of the other, or the middle vertex of one does (embedded); PAIR_SHARED where the cutter found a
!< shared wall; PAIR_NONE otherwise.
!<
!< Classification (`classify`, `classify_2d`; Chesshire & Henshaw, Lee & Baeder). Holes (from the cutters) are 0.
!< A vertex wants a donor if it lies on an overlap boundary (faces covered by TOPO_BOUNDARY records named `overlap`;
!< mandatory), or by the rule: SELECT_EXPLICIT, within `layers` index rings of a hole (the fringe); SELECT_IMPLICIT,
!< when a candidate's donor cell is smaller than its own cells (ratio < 1). Vertices within `immune` index layers of a
!< wall face never receive. A candidate is valid when every node of its width-q donor stencil (WP2 stencils, nonzero
!< weights) is no hole, and, with explicit interpolation (default), no receptor either; candidates are tried by
!< increasing donor volume. Validity depends on the neighbours' classes, so the classification is iterated
!< (Jacobi: every block from the previous state, hence independent of the block order) to a fixed point. A mandatory
!< vertex is a receptor point whether or not it finds a donor (CMPGRD): it never serves as an explicit stencil node,
!< which also keeps the iteration from oscillating; without a valid donor it is an orphan. A desired one without a
!< valid donor stays a field point. With implicit interpolation,
!< receptors in each other's stencils (circular interpolation, Sherer & Scott) are counted.
!<
!< Wall projection (`wall_projection`, PEGASUS 5 / PROGRD): for a receptor block a and a donor block b sharing a curved
!< wall discretized twice, every wall vertex of a moves onto the closest point of b's wall (its cutter) when that is
!< within `max_distance`; the shift is carried up the wall-normal grid line, scaled by 1 - s / max_distance (s: the
!< distance from the wall vertex along the line) down to zero. The projected point is moved 1e-9 max_distance along
!< the wall normal into the fluid, so the rounding of the foot point never leaves it outside the donor block's exact
!< containment test. The shifted points serve the donor search only
!< (`donor_candidates(..., shifts=)`); the grid never moves.

use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P
use fossil_facet_object,       only : facet_object
use fossil_surface_stl_object, only : surface_stl_object, SIGN_FWN => SIGN_FAST_WINDING_NUMBER
use fossil_curve_object,       only : curve_object
use fossil_curvilinear_block,  only : curvilinear_block_t, curvilinear_block_2d_t, TOPO_BOUNDARY
use fossil_predicates,         only : orient2d_xy
use fossil_box_query,          only : BOX_CUT
use fossil_block_locate,       only : LOC_INSIDE, LOC_OUTSIDE
use fossil_block_geometry,     only : hex_volume, quad_area, hex_validity, quad_validity
use fossil_block_stencil,      only : STENCIL_FAILED

implicit none
private
public :: wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, block_real_t
public :: hole_map, hole_map_2d, xray_flags, xray_flags_2d
public :: donor_set_t, overlap_matrix, overlap_matrix_2d, donor_candidates, donor_candidates_2d
public :: PAIR_NONE, PAIR_OVERLAP, PAIR_SHARED
public :: overset_mask_t, classify, classify_2d, SELECT_EXPLICIT, SELECT_IMPLICIT
public :: fringe_levels, fringe_levels_2d, overlap_width, overlap_width_2d
public :: ORPHAN_NONE, ORPHAN_NO_CANDIDATE, ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR, ORPHAN_BOUNDARY, ORPHAN_WALL
public :: FRINGE_INNER_1, FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
public :: wall_projection, wall_projection_2d

integer(I4P), parameter :: ORPHAN_NONE           = 0_I4P !< Not an orphan.
integer(I4P), parameter :: ORPHAN_NO_CANDIDATE   = 1_I4P !< No donor cell anywhere: insufficient overlap.
integer(I4P), parameter :: ORPHAN_STENCIL_HOLE   = 2_I4P !< Every candidate's stencil hits a hole.
integer(I4P), parameter :: ORPHAN_RECEPTOR_DONOR = 3_I4P !< The stencils hold receptors (explicit interpolation).
integer(I4P), parameter :: ORPHAN_BOUNDARY       = 4_I4P !< The stencils fail near a block boundary.
integer(I4P), parameter :: ORPHAN_WALL           = 5_I4P !< A wall vertex outside every donor (curved-wall mismatch).
integer(I4P), parameter :: FRINGE_INNER_1        = 1_I4P !< Shares a cell with a hole.
integer(I4P), parameter :: FRINGE_INNER_2        = 2_I4P !< Shares a cell with an INNER 1 vertex.
integer(I4P), parameter :: FRINGE_OUTER_1        = 11_I4P !< On an overlap boundary.
integer(I4P), parameter :: FRINGE_OUTER_2        = 12_I4P !< Shares a cell with an OUTER 1 vertex.

integer(I4P), parameter :: SELECT_EXPLICIT = 1_I4P !< Holes from cutters; the hole fringe and overlap boundaries receive.
integer(I4P), parameter :: SELECT_IMPLICIT = 2_I4P !< Lee & Baeder: a point receives when a smaller donor cell exists.

integer(I4P), parameter :: PAIR_NONE    = 0_I4P !< The blocks do not overlap.
integer(I4P), parameter :: PAIR_OVERLAP = 1_I4P !< They overlap (and may cut each other).
integer(I4P), parameter :: PAIR_SHARED  = 2_I4P !< They share a wall: overlapping without mutual cutting.

type :: block_logical_t
   !< Logical values per vertex of a block (linear vertex index, i fastest).
   logical, allocatable :: v(:) !< Values.
endtype block_logical_t

type :: block_real_t
   !< Real values per vertex of a block (linear vertex index, i fastest).
   real(R8P), allocatable :: v(:) !< Values.
endtype block_real_t

type :: donor_set_t
   !< Donor candidates of the vertices of one receptor block, CSR by vertex (linear index, i fastest).
   integer(I4P), allocatable :: start(:)           !< First candidate of each vertex [nv + 1].
   integer(I4P), allocatable :: block(:)           !< Donor block [nc].
   integer(I4P), allocatable :: cell(:,:)          !< Donor cell [d, nc].
   real(R8P),    allocatable :: xi(:,:)            !< Local coordinates [d, nc].
   integer(I4P), allocatable :: status(:)          !< LOC_* (never LOC_OUTSIDE) [nc].
   real(R8P),    allocatable :: donor_volume(:)    !< |volume| of the donor cell [nc].
   real(R8P),    allocatable :: ratio(:)           !< donor_volume / receptor_volume [nc].
   real(R8P),    allocatable :: min_jacobian(:)    !< Certified min det J of the donor cell over its volume [nc].
   real(R8P),    allocatable :: receptor_volume(:) !< Mean |volume| of the cells around the vertex [nv].
   real(R8P),    allocatable :: wall_distance(:)   !< Distance to the wall [nv] (if a wall is given).
endtype donor_set_t

type :: overset_mask_t
   !< Classification of the vertices of one block (linear vertex index, i fastest).
   integer(I4P), allocatable :: iota(:)   !< 1 field, 0 hole, -k receptor from block k.
   integer(I4P), allocatable :: donor(:)  !< Chosen candidate (index into the block's donor_set_t), 0 if none.
   integer(I4P), allocatable :: orphan(:) !< ORPHAN_* (mandatory vertices without a valid donor).
   integer(I4P), allocatable :: fringe(:) !< FRINGE_* level (0: none), from fringe_levels.
endtype overset_mask_t

type :: cls_t
   !< Work descriptor of one block for the classification.
   integer(I4P)              :: n(3) = 1_I4P     !< Vertices per direction.
   logical,      allocatable :: hole(:)          !< Holes [nv].
   logical,      allocatable :: mandatory(:)     !< Overlap-boundary vertices [nv].
   logical,      allocatable :: immune(:)        !< Never receptors [nv].
   logical,      allocatable :: wall(:)          !< On a wall face [nv].
   logical,      allocatable :: desired(:)       !< Want a donor by the selection rule [nv].
   integer(I4P), allocatable :: order(:)         !< Candidates of each vertex in preference order [nc].
   integer(I4P), allocatable :: nodes(:,:)       !< Donor stencil nodes (linear, 0 unused) [q^d, nc].
endtype cls_t

type :: ray_t
   !< Exact crossings of one X-ray.
   integer(I4P)              :: n = 0_I4P !< Crossings.
   real(R8P),    allocatable :: t(:)      !< Parameters, sorted.
   logical,      allocatable :: ent(:)    !< Entering the solid.
endtype ray_t

type :: tri_list_t
   !< Wall triangles of one block.
   integer(I4P)           :: n = 0_I4P      !< Triangles.
   real(R8P), allocatable :: v(:,:,:)       !< Vertices [3, 3, n].
   logical,   allocatable :: keep(:)        !< Kept after trimming [n].
endtype tri_list_t

type :: seg_list_t
   !< Wall polylines of one block (2-D): segments with their polyline.
   integer(I4P)              :: n = 0_I4P   !< Segments.
   real(R8P),    allocatable :: a(:,:), b(:,:) !< End points [2, n] (solid on the left of a -> b).
   integer(I4P), allocatable :: line(:)     !< Polyline of each segment [n].
   logical,      allocatable :: keep(:)     !< Kept after trimming [n].
endtype seg_list_t

contains
   ! ---------------------------------------------------------------------------------------------------------- cutters
   subroutine wall_cutter(blocks, cutter, wall_name, tolerance, gap, shared, owner, status)
   !< The wall cutter of a set of 3-D blocks (module notes). owner(f): the block of facet f; shared(a, b): blocks a
   !< and b share a wall. tolerance: trimming distance (default 0.1 of the square root of the element area).
   type(curvilinear_block_t), intent(in)            :: blocks(:)   !< Blocks.
   type(surface_stl_object),  intent(inout)         :: cutter      !< Cutter.
   character(*),              intent(in),  optional :: wall_name   !< Record name of the walls (default 'wall').
   real(R8P),                 intent(in),  optional :: tolerance   !< Trimming distance.
   real(R8P),                 intent(out), optional :: gap         !< Largest bilinear-to-triangle distance.
   logical,      allocatable, intent(out), optional :: shared(:,:) !< Shared walls [nb, nb].
   integer(I4P), allocatable, intent(out), optional :: owner(:)    !< Block of each facet.
   integer(I4P),              intent(out), optional :: status      !< 0, or 1 if no wall face was found.
   type(tri_list_t), allocatable :: tl(:)          !< Triangles per block.
   type(facet_object), allocatable :: f(:)         !< Facets.
   logical,  allocatable         :: sh(:,:)        !< Shared walls.
   integer(I4P), allocatable     :: own(:)         !< Owners.
   character(len=:), allocatable :: wname          !< Wall name.
   real(R8P)                     :: gmax           !< Gap.
   integer(I4P)                  :: b, t, n, k     !< Counters.

   wname = 'wall' ; if (present(wall_name)) wname = wall_name
   allocate(tl(size(blocks))) ; gmax = 0._R8P
   do b = 1, size(blocks, kind=I4P)
      call block_wall_triangles(blocks(b), wname, tl(b), gmax)
   enddo
   call trim_triangles(tl, tolerance, sh)
   n = 0_I4P
   do b = 1, size(tl, kind=I4P)
      n = n + count(tl(b)%keep)
   enddo
   if (present(status)) status = merge(0_I4P, 1_I4P, n > 0)
   allocate(f(n), own(n)) ; n = 0_I4P
   do b = 1, size(tl, kind=I4P)
      do t = 1, tl(b)%n
         if (.not. tl(b)%keep(t)) cycle
         n = n + 1_I4P ; own(n) = b
         do k = 1, 3
            f(n)%vertex(k) = vector_R8P(tl(b)%v(1, k, t), tl(b)%v(2, k, t), tl(b)%v(3, k, t))
         enddo
      enddo
   enddo
   if (n > 0) then
      call cutter%adopt_facets(facets=f)
      do b = 1, size(blocks, kind=I4P)
         if (any(own == b)) call cutter%set_patch(facets=pack([(k, k = 1, n)], own == b), name=blocks(b)%name)
      enddo
   endif
   if (present(gap)) gap = gmax
   if (present(shared)) call move_alloc(sh, shared)
   if (present(owner)) call move_alloc(own, owner)
   endsubroutine wall_cutter

   subroutine block_wall_triangles(blk, wname, tl, gmax)
   !< Wall triangles of one block, oriented with the block interior on the outer side, shorter diagonal.
   type(curvilinear_block_t), intent(in)    :: blk   !< Block.
   character(*),              intent(in)    :: wname !< Wall name.
   type(tri_list_t),          intent(inout) :: tl    !< Triangles.
   real(R8P),                 intent(inout) :: gmax  !< Gap.
   real(R8P), allocatable                   :: v(:,:,:), w(:,:,:) !< Work.
   integer(I4P)                             :: r, nd, d1, d2, m1, m2, p(3), q(3), nv(3), s, cap !< Counters.
   real(R8P)                                :: c(3, 4), xin(3), tri(3, 3, 2), nrm(3), cen(3) !< Work.
   integer(I4P)                             :: k                 !< Triangle.

   nv = [blk%ni, blk%nj, blk%nk] ; cap = 0_I4P ; allocate(v(3, 3, 64))
   do r = 1, blk%nrec
      if (blk%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blk%rec(r)%name) /= wname) cycle
      nd = (blk%rec(r)%face + 1_I4P) / 2_I4P ; d1 = mod(nd, 3_I4P) + 1_I4P ; d2 = mod(nd + 1_I4P, 3_I4P) + 1_I4P
      s = merge(1_I4P, -1_I4P, mod(blk%rec(r)%face, 2_I4P) == 1_I4P)
      do m2 = blk%rec(r)%lo(d2), blk%rec(r)%hi(d2) - 1_I4P
         do m1 = blk%rec(r)%lo(d1), blk%rec(r)%hi(d1) - 1_I4P
            p = blk%rec(r)%lo ; p(d1) = m1 ; p(d2) = m2
            q = p ; c(:, 1) = pt3(blk, q)
            q = p ; q(d1) = q(d1) + 1_I4P ; c(:, 2) = pt3(blk, q)
            q = p ; q(d1) = q(d1) + 1_I4P ; q(d2) = q(d2) + 1_I4P ; c(:, 3) = pt3(blk, q)
            q = p ; q(d2) = q(d2) + 1_I4P ; c(:, 4) = pt3(blk, q)
            q = p ; q(nd) = q(nd) + s ; q(d1) = min(q(d1) + 1_I4P, nv(d1))
            if (q(nd) < 1 .or. q(nd) > nv(nd)) cycle
            xin = pt3(blk, q)
            gmax = max(gmax, 0.25_R8P * norm2(c(:, 1) + c(:, 3) - c(:, 2) - c(:, 4)))
            if (norm2(c(:, 3) - c(:, 1)) <= norm2(c(:, 4) - c(:, 2))) then
               tri(:, :, 1) = c(:, [1, 2, 3]) ; tri(:, :, 2) = c(:, [1, 3, 4])
            else
               tri(:, :, 1) = c(:, [1, 2, 4]) ; tri(:, :, 2) = c(:, [2, 3, 4])
            endif
            do k = 1, 2
               nrm = cross(tri(:, 2, k) - tri(:, 1, k), tri(:, 3, k) - tri(:, 1, k))
               if (norm2(nrm) == 0._R8P) cycle ! collapsed (a pole)
               cen = sum(tri(:, :, k), dim=2) / 3._R8P
               if (dot_product(nrm, xin - cen) < 0._R8P) tri(:, :, k) = tri(:, [1, 3, 2], k)
               cap = cap + 1_I4P
               if (cap > size(v, 3)) then ! grow by copy (no self-referencing reshape)
                  allocate(w(3, 3, 2 * size(v, 3))) ; w(:, :, 1:cap - 1) = v(:, :, 1:cap - 1) ; call move_alloc(w, v)
               endif
               v(:, :, cap) = tri(:, :, k)
            enddo
         enddo
      enddo
   enddo
   tl%n = cap ; tl%v = v(:, :, 1:cap) ; allocate(tl%keep(cap)) ; tl%keep = .true.
   endsubroutine block_wall_triangles

   subroutine trim_triangles(tl, tolerance, sh)
   !< Trim overlapping patches of block pairs (module notes); sh(a, b): the patches of a and b overlap.
   type(tri_list_t),     intent(inout)        :: tl(:)     !< Triangles per block.
   real(R8P),            intent(in), optional :: tolerance !< Trimming distance.
   logical, allocatable, intent(out)          :: sh(:,:)   !< Shared walls.
   type(surface_stl_object), allocatable      :: sb(:)     !< Surface per block.
   type(facet_object),       allocatable      :: f(:)      !< Facets.
   logical,                  allocatable      :: near(:)   !< Triangles of a near b.
   integer(I4P)                               :: a, b, t, k, nb !< Counters.
   real(R8P)                                  :: mean_a, mean_b !< Mean sizes over the overlap.
   logical,                  allocatable      :: nearb(:)  !< Triangles of b near a.

   nb = size(tl, kind=I4P)
   allocate(sh(nb, nb), sb(nb)) ; sh = .false.
   do b = 1, nb
      if (tl(b)%n == 0) cycle
      allocate(f(tl(b)%n))
      do t = 1, tl(b)%n
         do k = 1, 3
            f(t)%vertex(k) = vector_R8P(tl(b)%v(1, k, t), tl(b)%v(2, k, t), tl(b)%v(3, k, t))
         enddo
      enddo
      call sb(b)%adopt_facets(facets=f)
   enddo
   do a = 1, nb
      do b = a + 1, nb
         if (tl(a)%n == 0 .or. tl(b)%n == 0) cycle
         call near_set(tl(a), sb(b), tolerance, near) ; call near_set(tl(b), sb(a), tolerance, nearb)
         if (.not. (any(near) .or. any(nearb))) cycle
         sh(a, b) = .true. ; sh(b, a) = .true.
         mean_a = mean_size(tl(a), near) ; mean_b = mean_size(tl(b), nearb)
         if (mean_a >= mean_b) then ! a coarser (ties: the higher index wins)
            tl(a)%keep = tl(a)%keep .and. .not. near
         else
            tl(b)%keep = tl(b)%keep .and. .not. nearb
         endif
      enddo
   enddo
   endsubroutine trim_triangles

   subroutine near_set(ta, sbs, tolerance, near)
   !< near(t): the centroid of triangle t lies within the tolerance of the surface sbs.
   type(tri_list_t),         intent(in)           :: ta        !< Triangles.
   type(surface_stl_object), intent(in)           :: sbs       !< Other surface.
   real(R8P),                intent(in), optional :: tolerance !< Distance.
   logical, allocatable,     intent(out)          :: near(:)   !< Flags [n].
   type(vector_R8P), allocatable                  :: pts(:)    !< Vertices.
   real(R8P),        allocatable                  :: d(:)      !< Distances.
   integer(I4P)                                   :: t         !< Counter.
   real(R8P)                                      :: tol       !< Tolerance.

   allocate(pts(ta%n), d(ta%n), near(ta%n))
   do t = 1, ta%n
      pts(t) = vector_R8P(sum(ta%v(1, :, t)) / 3._R8P, sum(ta%v(2, :, t)) / 3._R8P, sum(ta%v(3, :, t)) / 3._R8P)
   enddo
   call sbs%distance_many(points=pts, distances=d, is_signed=.false., is_square_root=.true.)
   do t = 1, ta%n
      tol = 0.1_R8P * sqrt(0.5_R8P * norm2(cross(ta%v(:, 2, t) - ta%v(:, 1, t), ta%v(:, 3, t) - ta%v(:, 1, t))))
      if (present(tolerance)) tol = tolerance
      near(t) = d(t) <= tol
   enddo
   endsubroutine near_set

   pure real(R8P) function mean_size(ta, mask) result(m)
   !< Mean square root of the triangle area over the masked triangles (huge if none).
   type(tri_list_t), intent(in) :: ta      !< Triangles.
   logical,          intent(in) :: mask(:) !< Mask.
   integer(I4P)                 :: t       !< Counter.

   m = 0._R8P
   do t = 1, ta%n
      if (mask(t)) m = m + sqrt(0.5_R8P * norm2(cross(ta%v(:, 2, t) - ta%v(:, 1, t), ta%v(:, 3, t) - ta%v(:, 1, t))))
   enddo
   if (count(mask) > 0) then
      m = m / count(mask)
   else
      m = huge(1._R8P)
   endif
   endfunction mean_size

   subroutine wall_cutter_2d(blocks, cutter, wall_name, tolerance, shared, owner, status)
   !< The wall cutter of a set of 2-D blocks: polylines along the wall edges, solid on the left (module notes).
   !< A polyline whose ends coincide is a closed loop; trimming splits polylines into open pieces.
   type(curvilinear_block_2d_t), intent(in)            :: blocks(:)   !< Blocks.
   type(curve_object),           intent(inout)         :: cutter      !< Cutter.
   character(*),                 intent(in),  optional :: wall_name   !< Record name of the walls (default 'wall').
   real(R8P),                    intent(in),  optional :: tolerance   !< Trimming distance (default 0.1 segment).
   logical,      allocatable,    intent(out), optional :: shared(:,:) !< Shared walls [nb, nb].
   integer(I4P), allocatable,    intent(out), optional :: owner(:)    !< Block of each loop.
   integer(I4P),                 intent(out), optional :: status      !< 0, or 1 if no wall edge was found.
   type(seg_list_t), allocatable :: sl(:)          !< Segments per block.
   logical,          allocatable :: sh(:,:)        !< Shared walls.
   real(R8P),        allocatable :: u(:), v(:)     !< Loop vertices.
   integer(I4P),     allocatable :: sizes(:), own(:) !< Loop sizes, owners.
   logical,          allocatable :: closed(:)      !< Loop closed.
   character(len=:), allocatable :: wname          !< Wall name.
   integer(I4P)                  :: b, s, s0, nl   !< Counters.

   wname = 'wall' ; if (present(wall_name)) wname = wall_name
   allocate(sl(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      call block_wall_segments(blocks(b), wname, sl(b))
   enddo
   call trim_segments(sl, tolerance, sh)
   allocate(u(0), v(0), sizes(0), own(0), closed(0))
   do b = 1, size(sl, kind=I4P)
      s = 1_I4P
      do while (s <= sl(b)%n)
         if (.not. sl(b)%keep(s)) then
            s = s + 1_I4P ; cycle
         endif
         s0 = s ! a run of kept segments of the same polyline
         do while (s < sl(b)%n)
            if (.not. (sl(b)%keep(s + 1) .and. sl(b)%line(s + 1) == sl(b)%line(s0))) exit
            s = s + 1_I4P
         enddo
         nl = s - s0 + 1_I4P
         if (all(sl(b)%a(:, s0) == sl(b)%b(:, s))) then ! closed loop
            u = [u, sl(b)%a(1, s0:s)] ; v = [v, sl(b)%a(2, s0:s)] ; sizes = [sizes, nl] ; closed = [closed, .true.]
         else
            u = [u, sl(b)%a(1, s0:s), sl(b)%b(1, s)] ; v = [v, sl(b)%a(2, s0:s), sl(b)%b(2, s)]
            sizes = [sizes, nl + 1_I4P] ; closed = [closed, .false.]
         endif
         own = [own, b] ; s = s + 1_I4P
      enddo
   enddo
   if (present(status)) status = merge(0_I4P, 1_I4P, size(sizes) > 0)
   if (size(sizes) > 0) call cutter%adopt_loops(u=u, v=v, loop_sizes=sizes, closed=closed)
   if (present(shared)) call move_alloc(sh, shared)
   if (present(owner)) call move_alloc(own, owner)
   endsubroutine wall_cutter_2d

   subroutine block_wall_segments(blk, wname, sl)
   !< Wall segments of one 2-D block, one polyline per record, solid on the left (the block interior on the right).
   type(curvilinear_block_2d_t), intent(in)    :: blk   !< Block.
   character(*),                 intent(in)    :: wname !< Wall name.
   type(seg_list_t),             intent(inout) :: sl    !< Segments.
   real(R8P),    allocatable                   :: a(:,:), b(:,:) !< Work.
   integer(I4P), allocatable                   :: line(:) !< Work.
   integer(I4P)                                :: r, nd, d1, s, m, n0, n, p(2), q(2), nv(2), nr !< Counters.
   real(R8P)                                   :: o   !< Orientation of the first segment.

   nv = [blk%ni, blk%nj] ; nr = 0_I4P ; n = 0_I4P
   do r = 1, blk%nrec ! capacity: the wall edges
      if (blk%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blk%rec(r)%name) /= wname) cycle
      nd = (blk%rec(r)%face + 1_I4P) / 2_I4P ; d1 = 3_I4P - nd
      n = n + max(blk%rec(r)%hi(d1) - blk%rec(r)%lo(d1), 0_I4P)
   enddo
   allocate(a(2, n), b(2, n), line(n)) ; n = 0_I4P
   do r = 1, blk%nrec
      if (blk%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blk%rec(r)%name) /= wname) cycle
      nr = nr + 1_I4P
      nd = (blk%rec(r)%face + 1_I4P) / 2_I4P ; d1 = 3_I4P - nd
      s = merge(1_I4P, -1_I4P, mod(blk%rec(r)%face, 2_I4P) == 1_I4P)
      n0 = n
      do m = blk%rec(r)%lo(d1), blk%rec(r)%hi(d1) - 1_I4P
         p = blk%rec(r)%lo(1:2) ; p(d1) = m ; q = p ; q(d1) = m + 1_I4P
         if (blk%u(p(1), p(2)) == blk%u(q(1), q(2)) .and. blk%v(p(1), p(2)) == blk%v(q(1), q(2))) cycle
         n = n + 1_I4P
         a(:, n) = [blk%u(p(1), p(2)), blk%v(p(1), p(2))] ; b(:, n) = [blk%u(q(1), q(2)), blk%v(q(1), q(2))]
         line(n) = nr
      enddo
      if (n == n0) cycle
      ! orientation from the first segment and the interior vertex next to its start
      p = blk%rec(r)%lo(1:2) ; p(d1) = blk%rec(r)%lo(d1) ; q = p ; q(nd) = q(nd) + s
      o = orient2d_xy(a(1, n0 + 1), a(2, n0 + 1), b(1, n0 + 1), b(2, n0 + 1), blk%u(q(1), q(2)), blk%v(q(1), q(2)))
      if (o > 0._R8P) then ! interior on the left: reverse the polyline
         a(:, n0 + 1:n) = a(:, n:n0 + 1:-1) ; b(:, n0 + 1:n) = b(:, n:n0 + 1:-1)
         call swap_ends(a(:, n0 + 1:n), b(:, n0 + 1:n))
      endif
   enddo
   sl%n = n ; sl%a = a(:, 1:n) ; sl%b = b(:, 1:n) ; sl%line = line(1:n) ; allocate(sl%keep(n)) ; sl%keep = .true.
   endsubroutine block_wall_segments

   pure subroutine swap_ends(a, b)
   !< Swap the end points of every segment.
   real(R8P), intent(inout) :: a(:,:), b(:,:) !< Ends.
   real(R8P)                :: t(size(a, 1), size(a, 2)) !< Work.

   t = a ; a = b ; b = t
   endsubroutine swap_ends

   subroutine trim_segments(sl, tolerance, sh)
   !< Trim overlapping polylines of block pairs (as trim_triangles).
   type(seg_list_t),     intent(inout)        :: sl(:)     !< Segments per block.
   real(R8P),            intent(in), optional :: tolerance !< Trimming distance.
   logical, allocatable, intent(out)          :: sh(:,:)   !< Shared walls.
   type(curve_object), allocatable            :: cb(:)     !< Curve per block.
   logical,            allocatable            :: na(:), nbb(:) !< Near flags.
   integer(I4P)                               :: a, b, nb  !< Counters.
   real(R8P)                                  :: ma, mb    !< Mean lengths.

   nb = size(sl, kind=I4P)
   allocate(sh(nb, nb), cb(nb)) ; sh = .false.
   do b = 1, nb
      if (sl(b)%n == 0) cycle
      call segment_curve(sl(b), cb(b))
   enddo
   do a = 1, nb
      do b = a + 1, nb
         if (sl(a)%n == 0 .or. sl(b)%n == 0) cycle
         call near_segs(sl(a), cb(b), tolerance, na) ; call near_segs(sl(b), cb(a), tolerance, nbb)
         if (.not. (any(na) .or. any(nbb))) cycle
         sh(a, b) = .true. ; sh(b, a) = .true.
         ma = mean_len(sl(a), na) ; mb = mean_len(sl(b), nbb)
         if (ma >= mb) then
            sl(a)%keep = sl(a)%keep .and. .not. na
         else
            sl(b)%keep = sl(b)%keep .and. .not. nbb
         endif
      enddo
   enddo
   endsubroutine trim_segments

   subroutine segment_curve(sl, c)
   !< The curve of a block's segments, each an open 2-vertex loop.
   type(seg_list_t),   intent(in)    :: sl !< Segments.
   type(curve_object), intent(inout) :: c  !< Curve.
   real(R8P)                         :: u(2 * sl%n), v(2 * sl%n) !< Interleaved ends.
   integer(I4P)                      :: s  !< Counter.

   do s = 1, sl%n
      u(2 * s - 1) = sl%a(1, s) ; u(2 * s) = sl%b(1, s) ; v(2 * s - 1) = sl%a(2, s) ; v(2 * s) = sl%b(2, s)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[(2_I4P, s = 1, sl%n)], closed=[(.false., s = 1, sl%n)])
   endsubroutine segment_curve

   subroutine near_segs(sa, cb, tolerance, near)
   !< near(s): the midpoint of segment s lies within the tolerance of the curve cb.
   type(seg_list_t),     intent(in)           :: sa        !< Segments.
   type(curve_object),   intent(in)           :: cb        !< Other curve.
   real(R8P),            intent(in), optional :: tolerance !< Distance.
   logical, allocatable, intent(out)          :: near(:)   !< Flags.
   real(R8P), allocatable                     :: p(:,:), d(:) !< Ends, distances.
   integer(I4P)                               :: s         !< Counter.
   real(R8P)                                  :: tol       !< Tolerance.

   allocate(p(2, sa%n), d(sa%n), near(sa%n))
   p = 0.5_R8P * (sa%a + sa%b)
   call cb%distance_many(points=p, distances=d, is_signed=.false.)
   do s = 1, sa%n
      tol = 0.1_R8P * norm2(sa%b(:, s) - sa%a(:, s)) ; if (present(tolerance)) tol = tolerance
      near(s) = abs(d(s)) <= tol
   enddo
   endsubroutine near_segs

   pure real(R8P) function mean_len(sa, mask) result(m)
   !< Mean length of the masked segments (huge if none).
   type(seg_list_t), intent(in) :: sa      !< Segments.
   logical,          intent(in) :: mask(:) !< Mask.
   integer(I4P)                 :: s       !< Counter.

   m = 0._R8P
   do s = 1, sa%n
      if (mask(s)) m = m + norm2(sa%b(:, s) - sa%a(:, s))
   enddo
   if (count(mask) > 0) then
      m = m / count(mask)
   else
      m = huge(1._R8P)
   endif
   endfunction mean_len

   ! ------------------------------------------------------------------------------------------------------ hole flags
   subroutine hole_flags(blocks, cutter, holes, offset, cuts, device, distances)
   !< Hole flags of every vertex of the 3-D blocks against the cutter (module notes).
   type(curvilinear_block_t),          intent(in)            :: blocks(:)    !< Blocks.
   type(surface_stl_object),           intent(inout)         :: cutter       !< Cutter.
   type(block_logical_t), allocatable, intent(out)           :: holes(:)     !< Holes per block.
   real(R8P),                          intent(in),  optional :: offset       !< Blank below this signed distance.
   logical,                            intent(in),  optional :: cuts(:)      !< Blocks cut (default all).
   logical,                            intent(in),  optional :: device       !< On the GPU.
   type(block_real_t),    allocatable, intent(out), optional :: distances(:) !< Signed distances per block.
   type(vector_R8P), allocatable :: p(:)  !< Points.
   real(R8P),        allocatable :: d(:)  !< Distances.
   real(R8P)                     :: off   !< Offset.
   logical                       :: dev, entered !< Device flags.
   integer(I4P)                  :: b, n, i, j, k !< Counters.

   off = 0._R8P ; if (present(offset)) off = offset
   dev = .false. ; if (present(device)) dev = device
   entered = .false.
   if (dev .and. .not. cutter%is_on_device) then
      call cutter%enter_device ; entered = .true.
   endif
   allocate(holes(size(blocks)))
   if (present(distances)) allocate(distances(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      n = blocks(b)%ni * blocks(b)%nj * blocks(b)%nk
      allocate(holes(b)%v(n)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      allocate(p(n), d(n)) ; n = 0_I4P
      do k = 1, blocks(b)%nk
         do j = 1, blocks(b)%nj
            do i = 1, blocks(b)%ni
               n = n + 1_I4P ; p(n) = vector_R8P(blocks(b)%x(i, j, k), blocks(b)%y(i, j, k), blocks(b)%z(i, j, k))
            enddo
         enddo
      enddo
      if (dev) then
         call cutter%distance_many_device(points=p, distances=d, is_signed=.true., sign_algorithm=SIGN_FWN, &
                                          is_square_root=.true.)
      else
         call cutter%distance_many(points=p, distances=d, is_signed=.true., sign_algorithm=SIGN_FWN, &
                                   is_square_root=.true.)
      endif
      holes(b)%v = d < off
      if (present(distances)) call move_alloc(d, distances(b)%v)
      if (allocated(d)) deallocate(d)
      deallocate(p)
   enddo
   if (entered) call cutter%exit_device
   endsubroutine hole_flags

   subroutine hole_flags_2d(blocks, cutter, holes, offset, cuts, device, distances)
   !< Hole flags of every vertex of the 2-D blocks against the cutter (module notes).
   type(curvilinear_block_2d_t),       intent(in)            :: blocks(:)    !< Blocks.
   type(curve_object),                 intent(inout)         :: cutter       !< Cutter.
   type(block_logical_t), allocatable, intent(out)           :: holes(:)     !< Holes per block.
   real(R8P),                          intent(in),  optional :: offset       !< Blank below this signed distance.
   logical,                            intent(in),  optional :: cuts(:)      !< Blocks cut (default all).
   logical,                            intent(in),  optional :: device       !< On the GPU.
   type(block_real_t),    allocatable, intent(out), optional :: distances(:) !< Signed distances per block.
   real(R8P),        allocatable :: p(:,:), d(:) !< Points, distances.
   real(R8P)                     :: off   !< Offset.
   logical                       :: dev, entered !< Device flags.
   integer(I4P)                  :: b, n, i, j !< Counters.

   off = 0._R8P ; if (present(offset)) off = offset
   dev = .false. ; if (present(device)) dev = device
   entered = .false.
   if (dev .and. .not. cutter%get_is_on_device()) then
      call cutter%enter_device ; entered = .true.
   endif
   allocate(holes(size(blocks)))
   if (present(distances)) allocate(distances(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      n = blocks(b)%ni * blocks(b)%nj
      allocate(holes(b)%v(n)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      allocate(p(2, n), d(n)) ; n = 0_I4P
      do j = 1, blocks(b)%nj
         do i = 1, blocks(b)%ni
            n = n + 1_I4P ; p(:, n) = [blocks(b)%u(i, j), blocks(b)%v(i, j)]
         enddo
      enddo
      if (dev) then
         call cutter%distance_many_device(points=p, distances=d, is_signed=.true., sign_algorithm=SIGN_FWN)
      else
         call cutter%distance_many(points=p, distances=d, is_signed=.true., sign_algorithm=SIGN_FWN)
      endif
      holes(b)%v = d < off
      if (present(distances)) call move_alloc(d, distances(b)%v)
      if (allocated(d)) deallocate(d)
      deallocate(p)
   enddo
   if (entered) call cutter%exit_device
   endsubroutine hole_flags_2d

   ! -------------------------------------------------------------------------------------- hole map and object X-rays
   subroutine hole_map(blocks, cutter, holes, n, cuts, pegasus, unresolved)
   !< Hole flags from a Cartesian hole map (module notes).
   type(curvilinear_block_t),          intent(in)            :: blocks(:)  !< Blocks.
   type(surface_stl_object),           intent(in)            :: cutter     !< Cutter.
   type(block_logical_t), allocatable, intent(out)           :: holes(:)   !< Holes per block.
   integer(I4P),                       intent(in)            :: n(3)       !< Map elements per direction.
   logical,                            intent(in),  optional :: cuts(:)    !< Blocks cut (default all).
   logical,                            intent(in),  optional :: pegasus    !< PEGASUS 5 line of sight.
   integer(I4P),                       intent(out), optional :: unresolved !< Vertices left unresolved (pegasus).
   real(R8P),    allocatable :: lo(:,:), hi(:,:)  !< Element boxes.
   integer(I4P), allocatable :: codes(:)          !< BOX_* codes.
   integer(I4P), allocatable :: state(:,:,:)      !< 0 cut, 1 outside, 2 inside.
   real(R8P)                 :: b0(3), b1(3), h(3), p(3) !< Map box, spacing, point.
   integer(I4P)              :: i, j, k, e, b, nv, m, ijk(3), nu !< Counters.
   logical                   :: peg, ins, ok   !< Flags.
   type(vector_R8P)          :: bmin, bmax     !< Cutter box.

   peg = .false. ; if (present(pegasus)) peg = pegasus
   bmin = cutter%get_bmin() ; bmax = cutter%get_bmax()
   h = ([bmax%x, bmax%y, bmax%z] - [bmin%x, bmin%y, bmin%z]) / (n - 2)
   b0 = [bmin%x, bmin%y, bmin%z] - h ; b1 = b0 + h * n
   allocate(lo(3, product(n)), hi(3, product(n)), codes(product(n)), state(n(1), n(2), n(3))) ; e = 0_I4P
   do k = 1, n(3)
      do j = 1, n(2)
         do i = 1, n(1)
            e = e + 1_I4P ; lo(:, e) = b0 + h * [i - 1, j - 1, k - 1] ; hi(:, e) = b0 + h * [i, j, k]
         enddo
      enddo
   enddo
   call cutter%box_classify_many(lo, hi, codes)
   state = reshape(merge(0_I4P, 2_I4P, codes == BOX_CUT), n)
   call paint_outside(state)
   allocate(holes(size(blocks))) ; nu = 0_I4P
   do b = 1, size(blocks, kind=I4P)
      nv = blocks(b)%ni * blocks(b)%nj * blocks(b)%nk
      allocate(holes(b)%v(nv)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      m = 0_I4P
      do k = 1, blocks(b)%nk
         do j = 1, blocks(b)%nj
            do i = 1, blocks(b)%ni
               m = m + 1_I4P ; p = [blocks(b)%x(i, j, k), blocks(b)%y(i, j, k), blocks(b)%z(i, j, k)]
               if (any(p < b0) .or. any(p >= b1)) cycle
               ijk = min(int((p - b0) / h) + 1_I4P, n)
               if (state(ijk(1), ijk(2), ijk(3)) /= 0) then
                  holes(b)%v(m) = state(ijk(1), ijk(2), ijk(3)) == 2
               elseif (cutter%distance(point=vector_R8P(p(1), p(2), p(3)), is_signed=.false., is_square_root=.true.) &
                       == 0._R8P) then
                  holes(b)%v(m) = .false. ! on the wall: a field point
               else
                  call line_of_sight_3d(cutter, state, b0, h, p, ijk, peg, ins, ok)
                  holes(b)%v(m) = ins .and. ok
                  if (.not. ok) nu = nu + 1_I4P
               endif
            enddo
         enddo
      enddo
   enddo
   if (present(unresolved)) unresolved = nu
   endsubroutine hole_map

   subroutine paint_outside(state)
   !< Paint state 2 (uncut) elements face-connected to the corner elements as 1 (outside), breadth first.
   integer(I4P), intent(inout) :: state(:,:,:) !< 0 cut, 1 outside, 2 inside (on entry: uncut).
   integer(I4P), allocatable   :: queue(:,:)   !< Queue.
   integer(I4P)                :: n(3), head, tail, c(3), d, s, q(3), a, e !< Work.

   n = shape(state) ; allocate(queue(3, product(n))) ; head = 1_I4P ; tail = 0_I4P
   do e = 0, 7
      c = merge(n, [1, 1, 1], [btest(e, 0), btest(e, 1), btest(e, 2)])
      if (state(c(1), c(2), c(3)) == 2) then
         state(c(1), c(2), c(3)) = 1 ; tail = tail + 1_I4P ; queue(:, tail) = c
      endif
   enddo
   do while (head <= tail)
      c = queue(:, head) ; head = head + 1_I4P
      do d = 1, 3
         do s = -1, 1, 2
            q = c ; q(d) = q(d) + s
            if (q(d) < 1 .or. q(d) > n(d)) cycle
            a = state(q(1), q(2), q(3))
            if (a /= 2) cycle
            state(q(1), q(2), q(3)) = 1 ; tail = tail + 1_I4P ; queue(:, tail) = q
         enddo
      enddo
   enddo
   endsubroutine paint_outside

   subroutine line_of_sight_3d(cutter, state, b0, h, p, ijk, peg, ins, ok)
   !< Status of the vertex p in the cut element ijk by a line of sight along a grid axis (module notes).
   type(surface_stl_object), intent(in)  :: cutter       !< Cutter.
   integer(I4P),             intent(in)  :: state(:,:,:) !< Element states.
   real(R8P),                intent(in)  :: b0(3), h(3), p(3) !< Map origin, spacing, point.
   integer(I4P),             intent(in)  :: ijk(3)       !< Element of p.
   logical,                  intent(in)  :: peg          !< PEGASUS 5 rule.
   logical,                  intent(out) :: ins, ok      !< Inside, resolved.
   real(R8P),    allocatable             :: t(:)         !< Crossings.
   integer(I4P), allocatable             :: fac(:)       !< Facets.
   logical,      allocatable             :: ent(:)       !< Entering flags.
   integer(I4P)                          :: d, s, q(3), nc, na, nb, steps !< Work.
   real(R8P)                             :: tend         !< End of the segment along the axis.

   ins = .false. ; ok = .false.
   do d = 1, 3
      na = mod(d, 3_I4P) + 1_I4P ; nb = mod(d + 1_I4P, 3_I4P) + 1_I4P
      do s = 1, -1, -2
         q = ijk ; steps = 0_I4P
         do
            q(d) = q(d) + s ; steps = steps + 1_I4P
            if (q(d) < 1 .or. q(d) > size(state, d)) exit
            if (state(q(1), q(2), q(3)) /= 0) exit
            if (peg) exit
         enddo
         if (q(d) < 1 .or. q(d) > size(state, d)) cycle
         if (state(q(1), q(2), q(3)) == 0) cycle
         tend = b0(d) + h(d) * (q(d) - 0.5_R8P)
         call cutter%line_crossings(d, p(na), p(nb), min(p(d), tend), max(p(d), tend), t, fac, ent, nc)
         if (any(t == p(d))) then ! on the wall: a field point
            ok = .true. ; ins = .false. ; return
         endif
         if (peg .and. nc > 0) cycle
         ins = (state(q(1), q(2), q(3)) == 2) .neqv. (mod(nc, 2_I4P) == 1)
         ok = .true. ; return
      enddo
   enddo
   endsubroutine line_of_sight_3d

   subroutine xray_flags(blocks, cutter, holes, n, cuts)
   !< Hole flags from object X-rays along z (module notes).
   type(curvilinear_block_t),          intent(in)           :: blocks(:) !< Blocks.
   type(surface_stl_object),           intent(in)           :: cutter    !< Cutter.
   type(block_logical_t), allocatable, intent(out)          :: holes(:)  !< Holes per block.
   integer(I4P),                       intent(in)           :: n(2)      !< Rays in x and y.
   logical,                            intent(in), optional :: cuts(:)   !< Blocks cut (default all).
   type(ray_t), allocatable :: rays(:,:)   !< Crossings per ray.
   integer(I4P), allocatable :: fac(:)     !< Facets.
   real(R8P)                :: b0(3), b1(3), h(2), p(3) !< Box, spacing, point.
   integer(I4P)             :: i, j, k, b, m, r(2), c, nv !< Counters.
   type(vector_R8P)         :: bmin, bmax  !< Cutter box.

   bmin = cutter%get_bmin() ; bmax = cutter%get_bmax()
   b0 = [bmin%x, bmin%y, bmin%z] ; b1 = [bmax%x, bmax%y, bmax%z] ; h = (b1(1:2) - b0(1:2)) / (n - 1)
   allocate(rays(n(1), n(2)))
   do j = 1, n(2)
      do i = 1, n(1)
         call cutter%line_crossings(3_I4P, b0(1) + h(1) * (i - 1), b0(2) + h(2) * (j - 1), b0(3) - 1._R8P, &
                                    b1(3) + 1._R8P, rays(i, j)%t, fac, rays(i, j)%ent, rays(i, j)%n)
      enddo
   enddo
   allocate(holes(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      nv = blocks(b)%ni * blocks(b)%nj * blocks(b)%nk
      allocate(holes(b)%v(nv)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      m = 0_I4P
      do k = 1, blocks(b)%nk
         do j = 1, blocks(b)%nj
            do i = 1, blocks(b)%ni
               m = m + 1_I4P ; p = [blocks(b)%x(i, j, k), blocks(b)%y(i, j, k), blocks(b)%z(i, j, k)]
               if (any(p < b0) .or. any(p > b1)) cycle
               r = nint((p(1:2) - b0(1:2)) / h) + 1_I4P
               c = count(rays(r(1), r(2))%t(1:rays(r(1), r(2))%n) < p(3))
               if (c > 0) holes(b)%v(m) = rays(r(1), r(2))%ent(c)
            enddo
         enddo
      enddo
   enddo
   endsubroutine xray_flags

   subroutine hole_map_2d(blocks, cutter, holes, n, cuts, pegasus, unresolved)
   !< Hole flags from a 2-D Cartesian hole map (as hole_map).
   type(curvilinear_block_2d_t),       intent(in)            :: blocks(:)  !< Blocks.
   type(curve_object),                 intent(in)            :: cutter     !< Cutter.
   type(block_logical_t), allocatable, intent(out)           :: holes(:)   !< Holes per block.
   integer(I4P),                       intent(in)            :: n(2)       !< Map elements per direction.
   logical,                            intent(in),  optional :: cuts(:)    !< Blocks cut (default all).
   logical,                            intent(in),  optional :: pegasus    !< PEGASUS 5 line of sight.
   integer(I4P),                       intent(out), optional :: unresolved !< Vertices left unresolved.
   real(R8P),    allocatable :: lo(:,:), hi(:,:)  !< Element boxes.
   integer(I4P), allocatable :: codes(:)          !< BOX_* codes.
   integer(I4P), allocatable :: state(:,:,:)      !< 0 cut, 1 outside, 2 inside.
   real(R8P)                 :: b0(2), b1(2), h(2), p(2), bb(2) !< Map box, spacing, point.
   integer(I4P)              :: i, j, e, b, nv, m, ij(2), nu !< Counters.
   logical                   :: peg, ins, ok   !< Flags.

   peg = .false. ; if (present(pegasus)) peg = pegasus
   bb = cutter%get_bmin() ; b0 = bb ; bb = cutter%get_bmax()
   h = (bb - b0) / (n - 2) ; b0 = b0 - h ; b1 = b0 + h * n
   allocate(lo(2, product(n)), hi(2, product(n)), codes(product(n)), state(n(1), n(2), 1)) ; e = 0_I4P
   do j = 1, n(2)
      do i = 1, n(1)
         e = e + 1_I4P ; lo(:, e) = b0 + h * [i - 1, j - 1] ; hi(:, e) = b0 + h * [i, j]
      enddo
   enddo
   call cutter%box_classify_many(lo, hi, codes)
   state = reshape(merge(0_I4P, 2_I4P, codes == BOX_CUT), [n(1), n(2), 1])
   call paint_outside(state)
   allocate(holes(size(blocks))) ; nu = 0_I4P
   do b = 1, size(blocks, kind=I4P)
      nv = blocks(b)%ni * blocks(b)%nj
      allocate(holes(b)%v(nv)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      m = 0_I4P
      do j = 1, blocks(b)%nj
         do i = 1, blocks(b)%ni
            m = m + 1_I4P ; p = [blocks(b)%u(i, j), blocks(b)%v(i, j)]
            if (any(p < b0) .or. any(p >= b1)) cycle
            ij = min(int((p - b0) / h) + 1_I4P, n)
            if (state(ij(1), ij(2), 1) /= 0) then
               holes(b)%v(m) = state(ij(1), ij(2), 1) == 2
            elseif (on_curve(cutter, p)) then
               holes(b)%v(m) = .false. ! on the wall: a field point
            else
               call line_of_sight_2d(cutter, state(:,:,1), b0, h, p, ij, peg, ins, ok)
               holes(b)%v(m) = ins .and. ok
               if (.not. ok) nu = nu + 1_I4P
            endif
         enddo
      enddo
   enddo
   if (present(unresolved)) unresolved = nu
   endsubroutine hole_map_2d

   logical function on_curve(cutter, p) result(yes)
   !< The point is exactly on the curve (distance 0).
   type(curve_object), intent(in) :: cutter !< Curve.
   real(R8P),          intent(in) :: p(2)   !< Point.
   real(R8P)                      :: d(1)   !< Distance.

   call cutter%distance_many(points=reshape(p, [2, 1]), distances=d, is_signed=.false.)
   yes = d(1) == 0._R8P
   endfunction on_curve

   subroutine line_of_sight_2d(cutter, state, b0, h, p, ij, peg, ins, ok)
   !< Status of the vertex p in the cut element ij by a line of sight along a grid axis (as line_of_sight_3d).
   type(curve_object), intent(in)  :: cutter     !< Cutter.
   integer(I4P),       intent(in)  :: state(:,:) !< Element states.
   real(R8P),          intent(in)  :: b0(2), h(2), p(2) !< Map origin, spacing, point.
   integer(I4P),       intent(in)  :: ij(2)      !< Element of p.
   logical,            intent(in)  :: peg        !< PEGASUS 5 rule.
   logical,            intent(out) :: ins, ok    !< Inside, resolved.
   real(R8P),    allocatable       :: t(:)       !< Crossings.
   integer(I4P), allocatable       :: seg(:)     !< Segments.
   logical,      allocatable       :: ent(:)     !< Entering flags.
   integer(I4P)                    :: d, s, q(2), nc !< Work.
   real(R8P)                       :: tend       !< End of the segment along the axis.

   ins = .false. ; ok = .false.
   do d = 1, 2
      do s = 1, -1, -2
         q = ij
         do
            q(d) = q(d) + s
            if (q(d) < 1 .or. q(d) > size(state, d)) exit
            if (state(q(1), q(2)) /= 0) exit
            if (peg) exit
         enddo
         if (q(d) < 1 .or. q(d) > size(state, d)) cycle
         if (state(q(1), q(2)) == 0) cycle
         tend = b0(d) + h(d) * (q(d) - 0.5_R8P)
         call cutter%line_crossings(d, p(3 - d), min(p(d), tend), max(p(d), tend), t, seg, ent, nc)
         if (any(t == p(d))) then
            ok = .true. ; ins = .false. ; return
         endif
         if (peg .and. nc > 0) cycle
         ins = (state(q(1), q(2)) == 2) .neqv. (mod(nc, 2_I4P) == 1)
         ok = .true. ; return
      enddo
   enddo
   endsubroutine line_of_sight_2d

   subroutine xray_flags_2d(blocks, cutter, holes, n, cuts)
   !< Hole flags from object X-rays along v through n points in u (as xray_flags).
   type(curvilinear_block_2d_t),       intent(in)           :: blocks(:) !< Blocks.
   type(curve_object),                 intent(in)           :: cutter    !< Cutter.
   type(block_logical_t), allocatable, intent(out)          :: holes(:)  !< Holes per block.
   integer(I4P),                       intent(in)           :: n         !< Rays.
   logical,                            intent(in), optional :: cuts(:)   !< Blocks cut (default all).
   type(ray_t), allocatable :: rays(:)     !< Crossings per ray.
   integer(I4P), allocatable :: seg(:)     !< Segments.
   real(R8P)                :: b0(2), b1(2), h, p(2) !< Box, spacing, point.
   integer(I4P)             :: i, j, b, m, r, c, nv !< Counters.

   b0 = cutter%get_bmin() ; b1 = cutter%get_bmax() ; h = (b1(1) - b0(1)) / (n - 1)
   allocate(rays(n))
   do i = 1, n
      call cutter%line_crossings(2_I4P, b0(1) + h * (i - 1), b0(2) - 1._R8P, b1(2) + 1._R8P, rays(i)%t, seg, &
                                 rays(i)%ent, rays(i)%n)
   enddo
   allocate(holes(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      nv = blocks(b)%ni * blocks(b)%nj
      allocate(holes(b)%v(nv)) ; holes(b)%v = .false.
      if (present(cuts)) then
         if (.not. cuts(b)) cycle
      endif
      m = 0_I4P
      do j = 1, blocks(b)%nj
         do i = 1, blocks(b)%ni
            m = m + 1_I4P ; p = [blocks(b)%u(i, j), blocks(b)%v(i, j)]
            if (any(p < b0) .or. any(p > b1)) cycle
            r = nint((p(1) - b0(1)) / h) + 1_I4P
            c = count(rays(r)%t(1:rays(r)%n) < p(2))
            if (c > 0) holes(b)%v(m) = rays(r)%ent(c)
         enddo
      enddo
   enddo
   endsubroutine xray_flags_2d

   ! ------------------------------------------------------------------------------------------------ donor candidates
   subroutine overlap_matrix(blocks, pair, shared)
   !< The three-valued pair matrix of a set of 3-D blocks (build_search done) (module notes).
   type(curvilinear_block_t), intent(in)           :: blocks(:)   !< Blocks.
   integer(I4P), allocatable, intent(out)          :: pair(:,:)   !< PAIR_* [nb, nb].
   logical,                   intent(in), optional :: shared(:,:) !< Shared walls (wall_cutter).
   real(R8P),    allocatable :: px(:), py(:), pz(:), xi(:,:) !< Points, local coordinates.
   integer(I4P), allocatable :: cell(:,:), st(:) !< Results.
   integer(I4P)              :: a, b, nb         !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(pair(nb, nb)) ; pair = PAIR_NONE
   do a = 1, nb
      call boundary_points_3d(blocks(a), px, py, pz)
      do b = 1, nb
         if (a == b) cycle
         call blocks(b)%locate_many(px, py, pz, cell, xi, st)
         if (any(st == LOC_INSIDE)) then
            pair(a, b) = PAIR_OVERLAP ; pair(b, a) = PAIR_OVERLAP
         endif
      enddo
   enddo
   if (present(shared)) then
      where (shared) pair = PAIR_SHARED
   endif
   do a = 1, nb
      pair(a, a) = PAIR_NONE
   enddo
   endsubroutine overlap_matrix

   subroutine boundary_points_3d(blk, px, py, pz)
   !< The boundary vertices of a block and its middle vertex.
   type(curvilinear_block_t), intent(in)  :: blk !< Block.
   real(R8P), allocatable,    intent(out) :: px(:), py(:), pz(:) !< Points.
   integer(I4P)                           :: i, j, k, n !< Counters.

   allocate(px(blk%ni * blk%nj * blk%nk), py(blk%ni * blk%nj * blk%nk), pz(blk%ni * blk%nj * blk%nk)) ; n = 0_I4P
   do k = 1, blk%nk
      do j = 1, blk%nj
         do i = 1, blk%ni
            if (.not. (i == 1 .or. j == 1 .or. k == 1 .or. i == blk%ni .or. j == blk%nj .or. k == blk%nk .or. &
                (i == (blk%ni + 1) / 2 .and. j == (blk%nj + 1) / 2 .and. k == (blk%nk + 1) / 2))) cycle
            n = n + 1_I4P ; px(n) = blk%x(i, j, k) ; py(n) = blk%y(i, j, k) ; pz(n) = blk%z(i, j, k)
         enddo
      enddo
   enddo
   px = px(1:n) ; py = py(1:n) ; pz = pz(1:n)
   endsubroutine boundary_points_3d

   subroutine donor_candidates(blocks, donors, pair, wall, device, shifts)
   !< Donor candidates of every vertex of every 3-D block (build_search done) (module notes).
   type(curvilinear_block_t), target,  intent(inout)        :: blocks(:) !< Blocks.
   type(donor_set_t),    allocatable,  intent(out)          :: donors(:) !< Candidates per receptor block.
   integer(I4P),                       intent(in), optional :: pair(:,:) !< PAIR_* (default: every pair).
   type(surface_stl_object),           intent(in), optional :: wall      !< Wall (for the wall distances).
   logical,                            intent(in), optional :: device    !< Location on the GPU.
   type(block_real_t),                 intent(in), optional :: shifts(:,:) !< Search shifts [nb, nb] (wall_projection):
                                                                       !< shifts(a, b)%v(3 nv) for receptors of a in b.
   real(R8P),    allocatable :: px(:), py(:), pz(:), xi(:,:), vol(:,:,:), cand_xi(:,:,:), qx(:), qy(:), qz(:) !< Work.
   integer(I4P), allocatable :: cell(:,:), st(:), cand_cell(:,:,:), cand_st(:,:), cnt(:) !< Work.
   real(R8P),    allocatable :: dw(:) !< Wall distances.
   integer(I4P), allocatable :: bw(:) !< Wall bodies.
   type(vector_R8P), allocatable :: pv(:) !< Points.
   logical                   :: dev  !< Device.
   integer(I4P)              :: a, b, nb, nv, v, c, k, i, j, kk, nc, slot !< Counters.
   real(R8P)                 :: jm, o !< Validity bound, orientation.
   integer(I4P)              :: cs    !< Validity status.

   dev = .false. ; if (present(device)) dev = device
   nb = size(blocks, kind=I4P) ; allocate(donors(nb))
   do a = 1, nb
      nv = blocks(a)%ni * blocks(a)%nj * blocks(a)%nk
      call vertex_points_3d(blocks(a), px, py, pz)
      allocate(cand_cell(3, nb, nv), cand_xi(3, nb, nv), cand_st(nb, nv), cnt(nv)) ; cand_st = LOC_OUTSIDE
      do b = 1, nb
         if (b == a) cycle
         if (present(pair)) then
            if (pair(a, b) == PAIR_NONE) cycle
         endif
         qx = px ; qy = py ; qz = pz
         if (present(shifts)) then
            if (allocated(shifts(a, b)%v)) then
               qx = px + shifts(a, b)%v(1:3 * nv:3) ; qy = py + shifts(a, b)%v(2:3 * nv:3)
               qz = pz + shifts(a, b)%v(3:3 * nv:3)
            endif
         endif
         if (dev) then
            call blocks(b)%locate_many_device(qx, qy, qz, cell, xi, st)
         else
            call blocks(b)%locate_many(qx, qy, qz, cell, xi, st)
         endif
         cand_cell(:, b, :) = cell ; cand_xi(:, b, :) = xi ; cand_st(b, :) = st
      enddo
      cnt = count(cand_st /= LOC_OUTSIDE, dim=1) ; nc = sum(cnt)
      allocate(donors(a)%start(nv + 1), donors(a)%block(nc), donors(a)%cell(3, nc), donors(a)%xi(3, nc), &
               donors(a)%status(nc), donors(a)%donor_volume(nc), donors(a)%ratio(nc), donors(a)%min_jacobian(nc))
      call blocks(a)%geometry(volume=vol)
      donors(a)%receptor_volume = vertex_volumes_3d(abs(vol))
      donors(a)%start(1) = 1_I4P ; slot = 0_I4P
      do v = 1, nv
         do b = 1, nb
            if (cand_st(b, v) == LOC_OUTSIDE) cycle
            slot = slot + 1_I4P ; c = slot
            donors(a)%block(c) = b ; donors(a)%cell(:, c) = cand_cell(:, b, v) ; donors(a)%xi(:, c) = cand_xi(:, b, v)
            donors(a)%status(c) = cand_st(b, v)
            i = cand_cell(1, b, v) ; j = cand_cell(2, b, v) ; kk = cand_cell(3, b, v)
            donors(a)%donor_volume(c) = abs(hex_volume(blocks(b)%corners(i, j, kk)))
            donors(a)%ratio(c) = donors(a)%donor_volume(c) / max(donors(a)%receptor_volume(v), tiny(1._R8P))
            o = blocks(b)%orient
            call hex_validity(blocks(b)%corners(i, j, kk), o, cs, jm)
            donors(a)%min_jacobian(c) = jm / max(donors(a)%donor_volume(c), tiny(1._R8P))
         enddo
         donors(a)%start(v + 1) = slot + 1_I4P
      enddo
      if (present(wall)) then
         allocate(pv(nv), dw(nv), bw(nv))
         do k = 1, nv
            pv(k) = vector_R8P(px(k), py(k), pz(k))
         enddo
         call wall%wall_distance_many(points=pv, distances=dw, bodies=bw)
         call move_alloc(dw, donors(a)%wall_distance) ; deallocate(pv, bw)
      endif
      deallocate(cand_cell, cand_xi, cand_st, cnt)
   enddo
   endsubroutine donor_candidates

   subroutine vertex_points_3d(blk, px, py, pz)
   !< Every vertex of a block, linear index (i fastest).
   type(curvilinear_block_t), intent(in)  :: blk !< Block.
   real(R8P), allocatable,    intent(out) :: px(:), py(:), pz(:) !< Points.

   px = reshape(blk%x, [size(blk%x)]) ; py = reshape(blk%y, [size(blk%y)]) ; pz = reshape(blk%z, [size(blk%z)])
   endsubroutine vertex_points_3d

   pure function vertex_volumes_3d(vol) result(vv)
   !< Mean of the cell volumes around each vertex.
   real(R8P), intent(in)  :: vol(:,:,:) !< Cell volumes [ni-1, nj-1, nk-1].
   real(R8P), allocatable :: vv(:)      !< Per vertex [ni nj nk].
   real(R8P), allocatable :: acc(:,:,:), cnt(:,:,:) !< Work.
   integer(I4P)           :: n(3), a, b, c !< Counters.

   n = shape(vol) + 1
   allocate(acc(n(1), n(2), n(3)), cnt(n(1), n(2), n(3))) ; acc = 0._R8P ; cnt = 0._R8P
   do c = 0, 1
      do b = 0, 1
         do a = 0, 1
            acc(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b, 1 + c:n(3) - 1 + c) = &
               acc(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b, 1 + c:n(3) - 1 + c) + vol
            cnt(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b, 1 + c:n(3) - 1 + c) = &
               cnt(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b, 1 + c:n(3) - 1 + c) + 1._R8P
         enddo
      enddo
   enddo
   vv = reshape(acc / cnt, [product(n)])
   endfunction vertex_volumes_3d

   subroutine overlap_matrix_2d(blocks, pair, shared)
   !< The pair matrix of a set of 2-D blocks (as overlap_matrix).
   type(curvilinear_block_2d_t), intent(in)           :: blocks(:)   !< Blocks.
   integer(I4P), allocatable,    intent(out)          :: pair(:,:)   !< PAIR_* [nb, nb].
   logical,                      intent(in), optional :: shared(:,:) !< Shared walls.
   real(R8P),    allocatable :: pu(:), pv(:), xi(:,:) !< Points, local coordinates.
   integer(I4P), allocatable :: cell(:,:), st(:) !< Results.
   integer(I4P)              :: a, b, nb, i, j, n !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(pair(nb, nb)) ; pair = PAIR_NONE
   do a = 1, nb
      allocate(pu(blocks(a)%ni * blocks(a)%nj), pv(blocks(a)%ni * blocks(a)%nj)) ; n = 0_I4P
      do j = 1, blocks(a)%nj
         do i = 1, blocks(a)%ni
            if (.not. (i == 1 .or. j == 1 .or. i == blocks(a)%ni .or. j == blocks(a)%nj .or. &
                (i == (blocks(a)%ni + 1) / 2 .and. j == (blocks(a)%nj + 1) / 2))) cycle
            n = n + 1_I4P ; pu(n) = blocks(a)%u(i, j) ; pv(n) = blocks(a)%v(i, j)
         enddo
      enddo
      do b = 1, nb
         if (a == b) cycle
         call blocks(b)%locate_many(pu(1:n), pv(1:n), cell, xi, st)
         if (any(st == LOC_INSIDE)) then
            pair(a, b) = PAIR_OVERLAP ; pair(b, a) = PAIR_OVERLAP
         endif
      enddo
      deallocate(pu, pv)
   enddo
   if (present(shared)) then
      where (shared) pair = PAIR_SHARED
   endif
   do a = 1, nb
      pair(a, a) = PAIR_NONE
   enddo
   endsubroutine overlap_matrix_2d

   subroutine donor_candidates_2d(blocks, donors, pair, wall, device, shifts)
   !< Donor candidates of every vertex of every 2-D block (as donor_candidates).
   type(curvilinear_block_2d_t), target, intent(inout)        :: blocks(:) !< Blocks.
   type(donor_set_t),    allocatable,    intent(out)          :: donors(:) !< Candidates per receptor block.
   integer(I4P),                         intent(in), optional :: pair(:,:) !< PAIR_* (default: every pair).
   type(curve_object),                   intent(in), optional :: wall      !< Wall (for the wall distances).
   logical,                              intent(in), optional :: device    !< Location on the GPU.
   type(block_real_t),                   intent(in), optional :: shifts(:,:) !< Search shifts [nb, nb]: (2 nv).
   real(R8P),    allocatable :: pu(:), pv(:), xi(:,:), ar(:,:), cand_xi(:,:,:), pts(:,:), dw(:), qu(:), qv(:) !< Work.
   integer(I4P), allocatable :: cell(:,:), st(:), cand_cell(:,:,:), cand_st(:,:), cnt(:), bw(:) !< Work.
   logical                   :: dev  !< Device.
   integer(I4P)              :: a, b, nb, nv, v, c, i, j, nc, slot, cs !< Counters.
   real(R8P)                 :: jm   !< Validity bound.

   dev = .false. ; if (present(device)) dev = device
   nb = size(blocks, kind=I4P) ; allocate(donors(nb))
   do a = 1, nb
      nv = blocks(a)%ni * blocks(a)%nj
      pu = reshape(blocks(a)%u, [nv]) ; pv = reshape(blocks(a)%v, [nv])
      allocate(cand_cell(2, nb, nv), cand_xi(2, nb, nv), cand_st(nb, nv), cnt(nv)) ; cand_st = LOC_OUTSIDE
      do b = 1, nb
         if (b == a) cycle
         if (present(pair)) then
            if (pair(a, b) == PAIR_NONE) cycle
         endif
         qu = pu ; qv = pv
         if (present(shifts)) then
            if (allocated(shifts(a, b)%v)) then
               qu = pu + shifts(a, b)%v(1:2 * nv:2) ; qv = pv + shifts(a, b)%v(2:2 * nv:2)
            endif
         endif
         if (dev) then
            call blocks(b)%locate_many_device(qu, qv, cell, xi, st)
         else
            call blocks(b)%locate_many(qu, qv, cell, xi, st)
         endif
         cand_cell(:, b, :) = cell ; cand_xi(:, b, :) = xi ; cand_st(b, :) = st
      enddo
      cnt = count(cand_st /= LOC_OUTSIDE, dim=1) ; nc = sum(cnt)
      allocate(donors(a)%start(nv + 1), donors(a)%block(nc), donors(a)%cell(2, nc), donors(a)%xi(2, nc), &
               donors(a)%status(nc), donors(a)%donor_volume(nc), donors(a)%ratio(nc), donors(a)%min_jacobian(nc))
      call blocks(a)%geometry(area=ar)
      donors(a)%receptor_volume = vertex_volumes_2d(abs(ar))
      donors(a)%start(1) = 1_I4P ; slot = 0_I4P
      do v = 1, nv
         do b = 1, nb
            if (cand_st(b, v) == LOC_OUTSIDE) cycle
            slot = slot + 1_I4P ; c = slot
            donors(a)%block(c) = b ; donors(a)%cell(:, c) = cand_cell(:, b, v) ; donors(a)%xi(:, c) = cand_xi(:, b, v)
            donors(a)%status(c) = cand_st(b, v)
            i = cand_cell(1, b, v) ; j = cand_cell(2, b, v)
            donors(a)%donor_volume(c) = abs(quad_area(blocks(b)%corners(i, j)))
            donors(a)%ratio(c) = donors(a)%donor_volume(c) / max(donors(a)%receptor_volume(v), tiny(1._R8P))
            call quad_validity(blocks(b)%corners(i, j), blocks(b)%orient, cs, jm)
            donors(a)%min_jacobian(c) = jm / max(donors(a)%donor_volume(c), tiny(1._R8P))
         enddo
         donors(a)%start(v + 1) = slot + 1_I4P
      enddo
      if (present(wall)) then
         allocate(pts(2, nv), dw(nv), bw(nv)) ; pts(1, :) = pu ; pts(2, :) = pv
         call wall%wall_distance_many(points=pts, distances=dw, bodies=bw)
         call move_alloc(dw, donors(a)%wall_distance) ; deallocate(pts, bw)
      endif
      deallocate(cand_cell, cand_xi, cand_st, cnt)
   enddo
   endsubroutine donor_candidates_2d

   pure function vertex_volumes_2d(ar) result(vv)
   !< Mean of the cell areas around each vertex.
   real(R8P), intent(in)  :: ar(:,:) !< Cell areas [ni-1, nj-1].
   real(R8P), allocatable :: vv(:)   !< Per vertex [ni nj].
   real(R8P), allocatable :: acc(:,:), cnt(:,:) !< Work.
   integer(I4P)           :: n(2), a, b !< Counters.

   n = shape(ar) + 1
   allocate(acc(n(1), n(2)), cnt(n(1), n(2))) ; acc = 0._R8P ; cnt = 0._R8P
   do b = 0, 1
      do a = 0, 1
         acc(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b) = acc(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b) + ar
         cnt(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b) = cnt(1 + a:n(1) - 1 + a, 1 + b:n(2) - 1 + b) + 1._R8P
      enddo
   enddo
   vv = reshape(acc / cnt, [product(n)])
   endfunction vertex_volumes_2d

   ! -------------------------------------------------------------------------------------------------- classification
   subroutine classify(blocks, donors, masks, rule, q, holes, explicit_interpolation, layers, immune, overlap_name, &
                       wall_name, iterations, circular)
   !< Classify the vertices of a set of 3-D blocks (module notes).
   type(curvilinear_block_t),          intent(in)            :: blocks(:)   !< Blocks.
   type(donor_set_t),                  intent(in)            :: donors(:)   !< Candidates (donor_candidates).
   type(overset_mask_t), allocatable,  intent(out)           :: masks(:)    !< Classification per block.
   integer(I4P),                       intent(in)            :: rule        !< SELECT_*.
   integer(I4P),                       intent(in)            :: q           !< Donor stencil width.
   type(block_logical_t),              intent(in),  optional :: holes(:)    !< Holes per block (default none).
   logical,                            intent(in),  optional :: explicit_interpolation !< Default .true.
   integer(I4P),                       intent(in),  optional :: layers      !< Fringe rings (explicit rule, default 1).
   integer(I4P),                       intent(in),  optional :: immune      !< Wall-normal immune layers (default 0).
   character(*),                       intent(in),  optional :: overlap_name !< Overlap boundary name ('overlap').
   character(*),                       intent(in),  optional :: wall_name   !< Wall name ('wall').
   integer(I4P),                       intent(out), optional :: iterations  !< Fixed-point iterations.
   integer(I4P),                       intent(out), optional :: circular    !< Circular receptor pairs.
   type(cls_t), allocatable :: cl(:) !< Descriptors.
   integer(I4P), allocatable :: nodes(:,:,:), sst(:), cell(:,:) !< Stencils.
   real(R8P),    allocatable :: w(:,:), xi(:,:) !< Weights.
   integer(I4P)              :: a, b, c, nb, nc, k, m !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(cl(nb))
   do a = 1, nb ! every block's sizes first: the stencil nodes of block b are linearised with cl(b)%n
      cl(a)%n = [blocks(a)%ni, blocks(a)%nj, blocks(a)%nk]
   enddo
   do a = 1, nb
      call face_flags_3d(blocks(a), cl(a), overlap_name, wall_name, immune)
      if (present(holes)) then
         cl(a)%hole = holes(a)%v
      else
         allocate(cl(a)%hole(product(cl(a)%n))) ; cl(a)%hole = .false.
      endif
      nc = size(donors(a)%block, kind=I4P)
      allocate(cl(a)%nodes(q ** 3, nc)) ; cl(a)%nodes = 0_I4P
      do b = 1, nb ! stencils of the candidates in block b
         m = count(donors(a)%block == b)
         if (m == 0) cycle
         cell = donors(a)%cell(:, pack([(c, c = 1, nc)], donors(a)%block == b))
         xi = donors(a)%xi(:, pack([(c, c = 1, nc)], donors(a)%block == b))
         call blocks(b)%stencils(cell, xi, q, nodes, w, sst)
         k = 0_I4P
         do c = 1, nc
            if (donors(a)%block(c) /= b) cycle
            k = k + 1_I4P
            if (sst(k) == STENCIL_FAILED) then
               cl(a)%nodes(:, c) = -1_I4P
            else
               where (w(:, k) /= 0._R8P) cl(a)%nodes(:, c) = nodes(1, :, k) + cl(b)%n(1) * &
                  (nodes(2, :, k) - 1_I4P + cl(b)%n(2) * (nodes(3, :, k) - 1_I4P))
            endif
         enddo
      enddo
   enddo
   call classify_core(cl, donors, masks, rule, explicit_interpolation, layers, iterations, circular)
   endsubroutine classify

   subroutine classify_2d(blocks, donors, masks, rule, q, holes, explicit_interpolation, layers, immune, overlap_name, &
                          wall_name, iterations, circular)
   !< Classify the vertices of a set of 2-D blocks (as classify).
   type(curvilinear_block_2d_t),       intent(in)            :: blocks(:)   !< Blocks.
   type(donor_set_t),                  intent(in)            :: donors(:)   !< Candidates (donor_candidates_2d).
   type(overset_mask_t), allocatable,  intent(out)           :: masks(:)    !< Classification per block.
   integer(I4P),                       intent(in)            :: rule        !< SELECT_*.
   integer(I4P),                       intent(in)            :: q           !< Donor stencil width.
   type(block_logical_t),              intent(in),  optional :: holes(:)    !< Holes per block (default none).
   logical,                            intent(in),  optional :: explicit_interpolation !< Default .true.
   integer(I4P),                       intent(in),  optional :: layers      !< Fringe rings (explicit rule, default 1).
   integer(I4P),                       intent(in),  optional :: immune      !< Wall-normal immune layers (default 0).
   character(*),                       intent(in),  optional :: overlap_name !< Overlap boundary name ('overlap').
   character(*),                       intent(in),  optional :: wall_name   !< Wall name ('wall').
   integer(I4P),                       intent(out), optional :: iterations  !< Fixed-point iterations.
   integer(I4P),                       intent(out), optional :: circular    !< Circular receptor pairs.
   type(cls_t), allocatable :: cl(:) !< Descriptors.
   integer(I4P), allocatable :: nodes(:,:,:), sst(:), cell(:,:) !< Stencils.
   real(R8P),    allocatable :: w(:,:), xi(:,:) !< Weights.
   integer(I4P)              :: a, b, c, nb, nc, k, m !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(cl(nb))
   do a = 1, nb
      cl(a)%n = [blocks(a)%ni, blocks(a)%nj, 1_I4P]
   enddo
   do a = 1, nb
      call face_flags_2d(blocks(a), cl(a), overlap_name, wall_name, immune)
      if (present(holes)) then
         cl(a)%hole = holes(a)%v
      else
         allocate(cl(a)%hole(product(cl(a)%n))) ; cl(a)%hole = .false.
      endif
      nc = size(donors(a)%block, kind=I4P)
      allocate(cl(a)%nodes(q ** 2, nc)) ; cl(a)%nodes = 0_I4P
      do b = 1, nb
         m = count(donors(a)%block == b)
         if (m == 0) cycle
         cell = donors(a)%cell(:, pack([(c, c = 1, nc)], donors(a)%block == b))
         xi = donors(a)%xi(:, pack([(c, c = 1, nc)], donors(a)%block == b))
         call blocks(b)%stencils(cell, xi, q, nodes, w, sst)
         k = 0_I4P
         do c = 1, nc
            if (donors(a)%block(c) /= b) cycle
            k = k + 1_I4P
            if (sst(k) == STENCIL_FAILED) then
               cl(a)%nodes(:, c) = -1_I4P
            else
               where (w(:, k) /= 0._R8P) cl(a)%nodes(:, c) = nodes(1, :, k) + cl(b)%n(1) * (nodes(2, :, k) - 1_I4P)
            endif
         enddo
      enddo
   enddo
   call classify_core(cl, donors, masks, rule, explicit_interpolation, layers, iterations, circular)
   endsubroutine classify_2d

   subroutine face_flags_3d(blk, cl, overlap_name, wall_name, immune)
   !< Mandatory (overlap boundary) and immune (near-wall) vertices of a 3-D block.
   type(curvilinear_block_t), intent(in)           :: blk          !< Block.
   type(cls_t),               intent(inout)        :: cl           !< Descriptor.
   character(*),              intent(in), optional :: overlap_name !< Overlap name.
   character(*),              intent(in), optional :: wall_name    !< Wall name.
   integer(I4P),              intent(in), optional :: immune       !< Immune layers.
   integer(I4P)                                    :: r, nimm      !< Counters.

   nimm = 0_I4P ; if (present(immune)) nimm = immune
   allocate(cl%mandatory(product(cl%n)), cl%immune(product(cl%n)), cl%wall(product(cl%n)))
   cl%mandatory = .false. ; cl%immune = .false. ; cl%wall = .false.
   do r = 1, blk%nrec
      if (blk%rec(r)%kind /= TOPO_BOUNDARY) cycle
      if (trim(blk%rec(r)%name) == name_or(wall_name, 'wall')) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, cl%wall)
      if (trim(blk%rec(r)%name) == name_or(overlap_name, 'overlap')) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, cl%mandatory)
      if (trim(blk%rec(r)%name) == name_or(wall_name, 'wall') .and. nimm > 0) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, nimm, cl%immune)
   enddo
   endsubroutine face_flags_3d

   subroutine face_flags_2d(blk, cl, overlap_name, wall_name, immune)
   !< Mandatory (overlap boundary) and immune (near-wall) vertices of a 2-D block.
   type(curvilinear_block_2d_t), intent(in)           :: blk          !< Block.
   type(cls_t),                  intent(inout)        :: cl           !< Descriptor.
   character(*),                 intent(in), optional :: overlap_name !< Overlap name.
   character(*),                 intent(in), optional :: wall_name    !< Wall name.
   integer(I4P),                 intent(in), optional :: immune       !< Immune layers.
   integer(I4P)                                       :: r, nimm      !< Counters.

   nimm = 0_I4P ; if (present(immune)) nimm = immune
   allocate(cl%mandatory(product(cl%n)), cl%immune(product(cl%n)), cl%wall(product(cl%n)))
   cl%mandatory = .false. ; cl%immune = .false. ; cl%wall = .false.
   do r = 1, blk%nrec
      if (blk%rec(r)%kind /= TOPO_BOUNDARY) cycle
      if (trim(blk%rec(r)%name) == name_or(wall_name, 'wall')) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, cl%wall)
      if (trim(blk%rec(r)%name) == name_or(overlap_name, 'overlap')) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, cl%mandatory)
      if (trim(blk%rec(r)%name) == name_or(wall_name, 'wall') .and. nimm > 0) &
         call mark_range(cl%n, blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, nimm, cl%immune)
   enddo
   endsubroutine face_flags_2d

   pure function name_or(name, default) result(n)
   !< The name if present, else the default.
   character(*), intent(in), optional :: name    !< Name.
   character(*), intent(in)           :: default !< Default.
   character(len=:), allocatable      :: n       !< Result.

   n = default ; if (present(name)) n = trim(name)
   endfunction name_or

   pure subroutine mark_range(n, face, lo, hi, depth, flag)
   !< Mark the vertices of a face range and of the next depth - 1 layers into the block.
   integer(I4P), intent(in)    :: n(3), face, lo(3), hi(3), depth !< Sizes, face, range, layers.
   logical,      intent(inout) :: flag(:)                         !< Flags [nv].
   integer(I4P)                :: l(3), h(3), nd, i, j, k         !< Work.

   nd = (face + 1_I4P) / 2_I4P ; l = lo ; h = hi
   if (mod(face, 2_I4P) == 1_I4P) then
      h(nd) = min(lo(nd) + depth - 1_I4P, n(nd))
   else
      l(nd) = max(hi(nd) - depth + 1_I4P, 1_I4P)
   endif
   do k = l(3), h(3)
      do j = l(2), h(2)
         do i = l(1), h(1)
            flag(i + n(1) * (j - 1_I4P + n(2) * (k - 1_I4P))) = .true.
         enddo
      enddo
   enddo
   endsubroutine mark_range

   subroutine classify_core(cl, donors, masks, rule, explicit_interpolation, layers, iterations, circular)
   !< The fixed point of the classification on the descriptors (module notes).
   type(cls_t),                       intent(inout)         :: cl(:)      !< Descriptors.
   type(donor_set_t),                 intent(in)            :: donors(:)  !< Candidates.
   type(overset_mask_t), allocatable, intent(out)           :: masks(:)   !< Classification.
   integer(I4P),                      intent(in)            :: rule       !< SELECT_*.
   logical,                           intent(in),  optional :: explicit_interpolation !< Default .true.
   integer(I4P),                      intent(in),  optional :: layers     !< Fringe rings.
   integer(I4P),                      intent(out), optional :: iterations !< Iterations.
   integer(I4P),                      intent(out), optional :: circular   !< Circular pairs.
   type(overset_mask_t), allocatable :: old(:) !< Previous state.
   integer(I4P), allocatable         :: key(:) !< Sort keys.
   logical                           :: expl, changed, valid !< Flags.
   integer(I4P)                      :: nb, a, v, c, j, it, nl, b, u, nv, ncirc, cu !< Counters.
   integer(I4P), parameter           :: MAX_IT = 50_I4P !< Iteration cap.

   expl = .true. ; if (present(explicit_interpolation)) expl = explicit_interpolation
   nl = 1_I4P ; if (present(layers)) nl = layers
   nb = size(cl, kind=I4P) ; allocate(masks(nb))
   do a = 1, nb
      nv = product(cl(a)%n)
      allocate(masks(a)%iota(nv), masks(a)%donor(nv), masks(a)%orphan(nv), masks(a)%fringe(nv))
      masks(a)%iota = merge(0_I4P, 1_I4P, cl(a)%hole) ; masks(a)%donor = 0_I4P ; masks(a)%orphan = 0_I4P
      masks(a)%fringe = 0_I4P
      allocate(cl(a)%desired(nv)) ; cl(a)%desired = .false.
      if (rule == SELECT_EXPLICIT) then
         call hole_rings(cl(a)%n, cl(a)%hole, nl, cl(a)%desired)
      else
         do v = 1, nv
            cl(a)%desired(v) = any(donors(a)%ratio(donors(a)%start(v):donors(a)%start(v + 1) - 1) < 1._R8P)
         enddo
      endif
      ! candidates in preference order: increasing donor volume (stable)
      allocate(cl(a)%order(size(donors(a)%block)))
      do v = 1, nv
         key = [(c, c = donors(a)%start(v), donors(a)%start(v + 1) - 1)]
         call sort_by_volume(key, donors(a)%donor_volume)
         cl(a)%order(donors(a)%start(v):donors(a)%start(v + 1) - 1) = key
      enddo
   enddo
   do it = 1, MAX_IT
      old = masks ; changed = .false.
      do a = 1, nb
         do v = 1, product(cl(a)%n)
            if (cl(a)%hole(v)) cycle
            masks(a)%iota(v) = 1_I4P ; masks(a)%donor(v) = 0_I4P ; masks(a)%orphan(v) = 0_I4P
            if (cl(a)%immune(v) .or. .not. (cl(a)%mandatory(v) .or. cl(a)%desired(v))) cycle
            do j = donors(a)%start(v), donors(a)%start(v + 1) - 1
               c = cl(a)%order(j) ; b = donors(a)%block(c)
               if (donors(a)%status(c) /= LOC_INSIDE .or. cl(a)%nodes(1, c) < 0) cycle
               valid = .true.
               do u = 1, size(cl(a)%nodes, 1, I4P)
                  if (cl(a)%nodes(u, c) == 0) cycle
                  if (old(b)%iota(cl(a)%nodes(u, c)) == 0 .or. (expl .and. (old(b)%iota(cl(a)%nodes(u, c)) < 0 .or. &
                      cl(b)%mandatory(cl(a)%nodes(u, c))))) then
                     valid = .false. ; exit
                  endif
               enddo
               if (valid) then
                  masks(a)%iota(v) = -b ; masks(a)%donor(v) = c ; exit
               endif
            enddo
            if (masks(a)%iota(v) == 1 .and. cl(a)%mandatory(v)) masks(a)%orphan(v) = orphan_cause(cl, old, donors, a, v, expl)
         enddo
         changed = changed .or. any(masks(a)%iota /= old(a)%iota)
      enddo
      if (.not. changed) exit
   enddo
   if (present(iterations)) iterations = it
   if (present(circular)) then ! receptors whose stencils hold each other
      ncirc = 0_I4P
      do a = 1, nb
         do v = 1, product(cl(a)%n)
            if (masks(a)%iota(v) >= 0) cycle
            c = masks(a)%donor(v) ; b = -masks(a)%iota(v)
            do u = 1, size(cl(a)%nodes, 1, I4P)
               if (cl(a)%nodes(u, c) <= 0) cycle
               if (masks(b)%iota(cl(a)%nodes(u, c)) /= -a) cycle
               cu = masks(b)%donor(cl(a)%nodes(u, c))
               if (any(cl(b)%nodes(:, cu) == v)) ncirc = ncirc + 1_I4P
            enddo
         enddo
      enddo
      circular = ncirc / 2_I4P
   endif
   endsubroutine classify_core

   integer(I4P) function orphan_cause(cl, old, donors, a, v, expl) result(cause)
   !< Why the mandatory vertex v of block a has no valid donor (the first cause found, in the order of the codes).
   type(cls_t),          intent(in) :: cl(:)     !< Descriptors.
   type(overset_mask_t), intent(in) :: old(:)    !< Classes the validity was tested against.
   type(donor_set_t),    intent(in) :: donors(:) !< Candidates.
   integer(I4P),         intent(in) :: a, v      !< Block, vertex.
   logical,              intent(in) :: expl      !< Explicit interpolation.
   integer(I4P)                     :: c, b, u, nd !< Counters.
   logical                          :: hole, rec, fail, any_inside !< Causes seen.

   hole = .false. ; rec = .false. ; fail = .false. ; any_inside = .false.
   do c = donors(a)%start(v), donors(a)%start(v + 1) - 1
      if (donors(a)%status(c) /= LOC_INSIDE) cycle
      any_inside = .true. ; b = donors(a)%block(c)
      if (cl(a)%nodes(1, c) < 0) then
         fail = .true. ; cycle
      endif
      do u = 1, size(cl(a)%nodes, 1, I4P)
         nd = cl(a)%nodes(u, c)
         if (nd == 0) cycle
         if (old(b)%iota(nd) == 0) hole = .true.
         if (expl .and. (old(b)%iota(nd) < 0 .or. cl(b)%mandatory(nd))) rec = .true.
      enddo
   enddo
   if (.not. any_inside) then
      cause = merge(ORPHAN_WALL, ORPHAN_NO_CANDIDATE, cl(a)%wall(v))
   elseif (hole) then
      cause = ORPHAN_STENCIL_HOLE
   elseif (rec) then
      cause = ORPHAN_RECEPTOR_DONOR
   else
      cause = ORPHAN_BOUNDARY
   endif
   endfunction orphan_cause

   subroutine fringe_levels(blocks, masks, levels)
   !< Fringe levels of the vertices of 3-D blocks (SUGGAR's node rules): INNER 1 shares a cell with a hole, INNER 2
   !< with an INNER 1 vertex; OUTER 1 is a mandatory receptor (overlap boundary: iota < 0 or orphan on it), OUTER 2
   !< shares a cell with an OUTER 1 vertex. levels (default 2) caps the depth.
   type(curvilinear_block_t), intent(in)           :: blocks(:) !< Blocks.
   type(overset_mask_t),      intent(inout)        :: masks(:)  !< Classification (fringe set here).
   integer(I4P),              intent(in), optional :: levels    !< Depth (1 or 2).
   integer(I4P)                                    :: a, nl     !< Counters.

   nl = 2_I4P ; if (present(levels)) nl = levels
   do a = 1, size(blocks, kind=I4P)
      call fringe_core([blocks(a)%ni, blocks(a)%nj, blocks(a)%nk], outer_flags_3d(blocks(a)), masks(a), nl)
   enddo
   endsubroutine fringe_levels

   subroutine fringe_levels_2d(blocks, masks, levels)
   !< Fringe levels of the vertices of 2-D blocks (as fringe_levels).
   type(curvilinear_block_2d_t), intent(in)           :: blocks(:) !< Blocks.
   type(overset_mask_t),         intent(inout)        :: masks(:)  !< Classification.
   integer(I4P),                 intent(in), optional :: levels    !< Depth (1 or 2).
   integer(I4P)                                       :: a, nl     !< Counters.

   nl = 2_I4P ; if (present(levels)) nl = levels
   do a = 1, size(blocks, kind=I4P)
      call fringe_core([blocks(a)%ni, blocks(a)%nj, 1_I4P], outer_flags_2d(blocks(a)), masks(a), nl)
   enddo
   endsubroutine fringe_levels_2d

   function outer_flags_3d(blk) result(f)
   !< Overlap-boundary vertices of a 3-D block (records named 'overlap').
   type(curvilinear_block_t), intent(in) :: blk  !< Block.
   logical, allocatable                  :: f(:) !< Flags.
   integer(I4P)                          :: r    !< Counter.

   allocate(f(blk%ni * blk%nj * blk%nk)) ; f = .false.
   do r = 1, blk%nrec
      if (blk%rec(r)%kind == TOPO_BOUNDARY .and. trim(blk%rec(r)%name) == 'overlap') &
         call mark_range([blk%ni, blk%nj, blk%nk], blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, f)
   enddo
   endfunction outer_flags_3d

   function outer_flags_2d(blk) result(f)
   !< Overlap-boundary vertices of a 2-D block (records named 'overlap').
   type(curvilinear_block_2d_t), intent(in) :: blk  !< Block.
   logical, allocatable                     :: f(:) !< Flags.
   integer(I4P)                             :: r    !< Counter.

   allocate(f(blk%ni * blk%nj)) ; f = .false.
   do r = 1, blk%nrec
      if (blk%rec(r)%kind == TOPO_BOUNDARY .and. trim(blk%rec(r)%name) == 'overlap') &
         call mark_range([blk%ni, blk%nj, 1_I4P], blk%rec(r)%face, blk%rec(r)%lo, blk%rec(r)%hi, 1_I4P, f)
   enddo
   endfunction outer_flags_2d

   subroutine fringe_core(n, outer, mask, nl)
   !< Node fringe rules on one block.
   integer(I4P),         intent(in)    :: n(3)     !< Vertices per direction.
   logical,              intent(in)    :: outer(:) !< Overlap-boundary vertices.
   type(overset_mask_t), intent(inout) :: mask     !< Classification.
   integer(I4P),         intent(in)    :: nl       !< Depth.
   logical, allocatable                :: s1(:), s2(:) !< Seeds, rings.

   mask%fringe = 0_I4P
   s1 = mask%iota == 0 ; call cell_ring(n, s1, s2)
   where (s2 .and. mask%iota /= 0) mask%fringe = FRINGE_INNER_1
   if (nl > 1) then
      s1 = mask%fringe == FRINGE_INNER_1 ; call cell_ring(n, s1, s2)
      where (s2 .and. mask%iota /= 0 .and. mask%fringe == 0) mask%fringe = FRINGE_INNER_2
   endif
   where (outer .and. mask%iota /= 0 .and. mask%fringe == 0) mask%fringe = FRINGE_OUTER_1
   if (nl > 1) then
      s1 = outer ; call cell_ring(n, s1, s2)
      where (s2 .and. mask%iota /= 0 .and. mask%fringe == 0) mask%fringe = FRINGE_OUTER_2
   endif
   endsubroutine fringe_core

   pure subroutine cell_ring(n, seed, ring)
   !< The vertices sharing a cell with a seed vertex (seeds included).
   integer(I4P),         intent(in)  :: n(3)    !< Vertices per direction.
   logical,              intent(in)  :: seed(:) !< Seeds.
   logical, allocatable, intent(out) :: ring(:) !< Ring.

   allocate(ring(size(seed)))
   call hole_rings(n, seed, 1_I4P, ring)
   ring = ring .or. seed
   endsubroutine cell_ring

   subroutine overlap_width(blocks, masks, d)
   !< Chesshire's overlap width per block pair: the smallest distance between a receptor of a from b and a receptor
   !< of b from a (huge if either set is empty).
   type(curvilinear_block_t), intent(in)  :: blocks(:) !< Blocks.
   type(overset_mask_t),      intent(in)  :: masks(:)  !< Classification.
   real(R8P), allocatable,    intent(out) :: d(:,:)    !< Widths [nb, nb].
   type(block_real_t), allocatable        :: x(:), y(:), z(:) !< Flat coordinates.
   integer(I4P)                           :: a, nb     !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(x(nb), y(nb), z(nb))
   do a = 1, nb
      x(a)%v = reshape(blocks(a)%x, [size(blocks(a)%x)]) ; y(a)%v = reshape(blocks(a)%y, [size(blocks(a)%y)])
      z(a)%v = reshape(blocks(a)%z, [size(blocks(a)%z)])
   enddo
   call width_core(x, y, z, masks, d)
   endsubroutine overlap_width

   subroutine overlap_width_2d(blocks, masks, d)
   !< Chesshire's overlap width per pair of 2-D blocks (as overlap_width).
   type(curvilinear_block_2d_t), intent(in)  :: blocks(:) !< Blocks.
   type(overset_mask_t),         intent(in)  :: masks(:)  !< Classification.
   real(R8P), allocatable,       intent(out) :: d(:,:)    !< Widths [nb, nb].
   type(block_real_t), allocatable           :: x(:), y(:), z(:) !< Flat coordinates.
   integer(I4P)                              :: a, nb     !< Counters.

   nb = size(blocks, kind=I4P) ; allocate(x(nb), y(nb), z(nb))
   do a = 1, nb
      x(a)%v = reshape(blocks(a)%u, [size(blocks(a)%u)]) ; y(a)%v = reshape(blocks(a)%v, [size(blocks(a)%v)])
      allocate(z(a)%v(size(x(a)%v))) ; z(a)%v = 0._R8P
   enddo
   call width_core(x, y, z, masks, d)
   endsubroutine overlap_width_2d

   subroutine width_core(x, y, z, masks, d)
   !< The overlap widths from flat coordinates.
   type(block_real_t),     intent(in)  :: x(:), y(:), z(:) !< Coordinates per block.
   type(overset_mask_t),   intent(in)  :: masks(:)  !< Classification.
   real(R8P), allocatable, intent(out) :: d(:,:)    !< Widths.
   integer(I4P)                        :: a, b, v, u, nb !< Counters.

   nb = size(masks, kind=I4P) ; allocate(d(nb, nb)) ; d = huge(1._R8P)
   do a = 1, nb
      do b = a + 1, nb
         do v = 1, size(masks(a)%iota, kind=I4P)
            if (masks(a)%iota(v) /= -b) cycle
            do u = 1, size(masks(b)%iota, kind=I4P)
               if (masks(b)%iota(u) /= -a) cycle
               d(a, b) = min(d(a, b), norm2([x(a)%v(v) - x(b)%v(u), y(a)%v(v) - y(b)%v(u), z(a)%v(v) - z(b)%v(u)]))
            enddo
         enddo
         d(b, a) = d(a, b)
      enddo
   enddo
   endsubroutine width_core

   pure subroutine hole_rings(n, hole, nl, ring)
   !< Vertices within nl index rings (Chebyshev) of a hole, holes excluded.
   integer(I4P), intent(in)  :: n(3), nl  !< Sizes, rings.
   logical,      intent(in)  :: hole(:)   !< Holes.
   logical,      intent(out) :: ring(:)   !< Ring vertices.
   integer(I4P)              :: i, j, k, di, dj, dk, ii, jj, kk, dz !< Counters.

   ring = .false. ; dz = merge(0_I4P, nl, n(3) == 1)
   do k = 1, n(3)
      do j = 1, n(2)
         do i = 1, n(1)
            if (.not. hole(i + n(1) * (j - 1_I4P + n(2) * (k - 1_I4P)))) cycle
            do dk = -dz, dz
               do dj = -nl, nl
                  do di = -nl, nl
                     ii = i + di ; jj = j + dj ; kk = k + dk
                     if (ii < 1 .or. jj < 1 .or. kk < 1 .or. ii > n(1) .or. jj > n(2) .or. kk > n(3)) cycle
                     ring(ii + n(1) * (jj - 1_I4P + n(2) * (kk - 1_I4P))) = .true.
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   ring = ring .and. .not. hole
   endsubroutine hole_rings

   pure subroutine sort_by_volume(key, vol)
   !< Stable insertion sort of candidate indices by increasing donor volume.
   integer(I4P), intent(inout) :: key(:) !< Candidates.
   real(R8P),    intent(in)    :: vol(:) !< Donor volumes.
   integer(I4P)                :: i, j, t !< Work.

   do i = 2, size(key)
      t = key(i) ; j = i - 1
      do while (j >= 1)
         if (vol(key(j)) <= vol(t)) exit
         key(j + 1) = key(j) ; j = j - 1
      enddo
      key(j + 1) = t
   enddo
   endsubroutine sort_by_volume

   ! ------------------------------------------------------------------------------------------------- wall projection
   subroutine wall_projection(blocks, a, b, shift, max_distance, wall_name)
   !< Search shifts of the vertices of block a for donors in block b (module notes); shift(3 nv) (zero off the walls).
   type(curvilinear_block_t), intent(in)           :: blocks(:)    !< Blocks.
   integer(I4P),              intent(in)           :: a, b         !< Receptor and donor blocks.
   real(R8P), allocatable,    intent(out)          :: shift(:)     !< Shifts [3 nv].
   real(R8P),                 intent(in)           :: max_distance !< Projection and decay distance.
   character(*),              intent(in), optional :: wall_name    !< Wall name ('wall').
   type(surface_stl_object)      :: wb      !< Wall of b.
   type(vector_R8P), allocatable :: p(:)    !< Wall vertices of a.
   real(R8P),        allocatable :: d(:), cp(:,:), nr(:,:) !< Distances, closest points, normals.
   integer(I4P),     allocatable :: idx(:,:) !< Wall vertex index and face.
   integer(I4P)                  :: r, nd, s, i, j, k, m, n, st, w(3), q(3), nv(3), lin !< Counters.
   real(R8P)                     :: sh(3), x0(3), x(3), dist, f !< Work.

   nv = [blocks(a)%ni, blocks(a)%nj, blocks(a)%nk]
   allocate(shift(3 * product(nv))) ; shift = 0._R8P
   call wall_cutter(blocks(b:b), wb, wall_name=wall_name, status=st)
   if (st /= 0) return
   n = 0_I4P ; allocate(p(product(nv)), idx(4, product(nv)))
   do r = 1, blocks(a)%nrec
      if (blocks(a)%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blocks(a)%rec(r)%name) /= name_or(wall_name, 'wall')) cycle
      do k = blocks(a)%rec(r)%lo(3), blocks(a)%rec(r)%hi(3)
         do j = blocks(a)%rec(r)%lo(2), blocks(a)%rec(r)%hi(2)
            do i = blocks(a)%rec(r)%lo(1), blocks(a)%rec(r)%hi(1)
               n = n + 1_I4P ; p(n) = vector_R8P(blocks(a)%x(i, j, k), blocks(a)%y(i, j, k), blocks(a)%z(i, j, k))
               idx(:, n) = [i, j, k, blocks(a)%rec(r)%face]
            enddo
         enddo
      enddo
   enddo
   if (n == 0) return
   allocate(d(n), cp(3, n), nr(3, n))
   call wb%query_many(points=p(1:n), distances=d, closest_points=cp, normals=nr, is_signed=.false.)
   do m = 1, n
      if (abs(d(m)) > max_distance) cycle
      w = idx(1:3, m) ; nd = (idx(4, m) + 1_I4P) / 2_I4P ; s = merge(1_I4P, -1_I4P, mod(idx(4, m), 2_I4P) == 1)
      x0 = [p(m)%x, p(m)%y, p(m)%z]
      sh = cp(:, m) + 1e-9_R8P * max_distance * nr(:, m) / max(norm2(nr(:, m)), tiny(1._R8P)) - x0
      q = w
      do
         x = [blocks(a)%x(q(1), q(2), q(3)), blocks(a)%y(q(1), q(2), q(3)), blocks(a)%z(q(1), q(2), q(3))]
         dist = norm2(x - x0) ; f = 1._R8P - dist / max_distance
         if (f <= 0._R8P) exit
         lin = q(1) + nv(1) * (q(2) - 1_I4P + nv(2) * (q(3) - 1_I4P))
         shift(3 * lin - 2:3 * lin) = f * sh
         q(nd) = q(nd) + s
         if (q(nd) < 1 .or. q(nd) > nv(nd)) exit
      enddo
   enddo
   endsubroutine wall_projection

   subroutine wall_projection_2d(blocks, a, b, shift, max_distance, wall_name)
   !< Search shifts of the vertices of the 2-D block a for donors in block b (as wall_projection); shift(2 nv).
   type(curvilinear_block_2d_t), intent(in)           :: blocks(:)    !< Blocks.
   integer(I4P),                 intent(in)           :: a, b         !< Receptor and donor blocks.
   real(R8P), allocatable,       intent(out)          :: shift(:)     !< Shifts [2 nv].
   real(R8P),                    intent(in)           :: max_distance !< Projection and decay distance.
   character(*),                 intent(in), optional :: wall_name    !< Wall name ('wall').
   type(curve_object)        :: wb      !< Wall of b.
   real(R8P),    allocatable :: p(:,:), d(:), cp(:,:), nr(:,:) !< Wall vertices, distances, closest points, normals.
   integer(I4P), allocatable :: idx(:,:) !< Wall vertex index and face.
   integer(I4P)              :: r, nd, s, i, j, m, n, st, q(2), nv(2), lin !< Counters.
   real(R8P)                 :: sh(2), x0(2), x(2), dist, f !< Work.

   nv = [blocks(a)%ni, blocks(a)%nj]
   allocate(shift(2 * product(nv))) ; shift = 0._R8P
   call wall_cutter_2d(blocks(b:b), wb, wall_name=wall_name, status=st)
   if (st /= 0) return
   n = 0_I4P ; allocate(p(2, product(nv)), idx(3, product(nv)))
   do r = 1, blocks(a)%nrec
      if (blocks(a)%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blocks(a)%rec(r)%name) /= name_or(wall_name, 'wall')) cycle
      do j = blocks(a)%rec(r)%lo(2), blocks(a)%rec(r)%hi(2)
         do i = blocks(a)%rec(r)%lo(1), blocks(a)%rec(r)%hi(1)
            n = n + 1_I4P ; p(:, n) = [blocks(a)%u(i, j), blocks(a)%v(i, j)] ; idx(:, n) = [i, j, blocks(a)%rec(r)%face]
         enddo
      enddo
   enddo
   if (n == 0) return
   allocate(d(n), cp(2, n), nr(2, n))
   call wb%query_many(points=p(:, 1:n), distances=d, closest_points=cp, normals=nr, is_signed=.false.)
   do m = 1, n
      if (abs(d(m)) > max_distance) cycle
      nd = (idx(3, m) + 1_I4P) / 2_I4P ; s = merge(1_I4P, -1_I4P, mod(idx(3, m), 2_I4P) == 1)
      x0 = p(:, m) ; q = idx(1:2, m)
      sh = cp(:, m) + 1e-9_R8P * max_distance * nr(:, m) / max(norm2(nr(:, m)), tiny(1._R8P)) - x0
      do
         x = [blocks(a)%u(q(1), q(2)), blocks(a)%v(q(1), q(2))]
         dist = norm2(x - x0) ; f = 1._R8P - dist / max_distance
         if (f <= 0._R8P) exit
         lin = q(1) + nv(1) * (q(2) - 1_I4P)
         shift(2 * lin - 1:2 * lin) = f * sh
         q(nd) = q(nd) + s
         if (q(nd) < 1 .or. q(nd) > nv(nd)) exit
      enddo
   enddo
   endsubroutine wall_projection_2d

   ! ------------------------------------------------------------------------------------------------------- utilities
   pure function pt3(blk, p) result(x)
   !< Vertex p of a 3-D block.
   type(curvilinear_block_t), intent(in) :: blk  !< Block.
   integer(I4P),              intent(in) :: p(3) !< Index.
   real(R8P)                             :: x(3) !< Coordinates.

   x = [blk%x(p(1), p(2), p(3)), blk%y(p(1), p(2), p(3)), blk%z(p(1), p(2), p(3))]
   endfunction pt3

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross
endmodule fossil_overset
