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

use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P
use fossil_facet_object,       only : facet_object
use fossil_surface_stl_object, only : surface_stl_object, SIGN_FWN => SIGN_FAST_WINDING_NUMBER
use fossil_curve_object,       only : curve_object
use fossil_curvilinear_block,  only : curvilinear_block_t, curvilinear_block_2d_t, TOPO_BOUNDARY
use fossil_predicates,         only : orient2d_xy
use fossil_box_query,          only : BOX_CUT

implicit none
private
public :: wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, block_real_t
public :: hole_map, hole_map_2d, xray_flags, xray_flags_2d

type :: block_logical_t
   !< Logical values per vertex of a block (linear vertex index, i fastest).
   logical, allocatable :: v(:) !< Values.
endtype block_logical_t

type :: block_real_t
   !< Real values per vertex of a block (linear vertex index, i fastest).
   real(R8P), allocatable :: v(:) !< Values.
endtype block_real_t

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
