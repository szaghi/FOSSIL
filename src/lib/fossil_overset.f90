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
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P
use fossil_facet_object,       only : facet_object
use fossil_surface_stl_object, only : surface_stl_object
use fossil_curve_object,       only : curve_object
use fossil_curvilinear_block,  only : curvilinear_block_t, curvilinear_block_2d_t, TOPO_BOUNDARY
use fossil_predicates,         only : orient2d_xy

implicit none
private
public :: wall_cutter, wall_cutter_2d

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
   real(R8P), allocatable                   :: v(:,:,:) !< Work.
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
               if (cap > size(v, 3)) v = reshape(v, [3, 3, 2 * size(v, 3)], pad=v)
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

   nv = [blk%ni, blk%nj] ; allocate(a(2, 0), b(2, 0), line(0)) ; nr = 0_I4P
   do r = 1, blk%nrec
      if (blk%rec(r)%kind /= TOPO_BOUNDARY .or. trim(blk%rec(r)%name) /= wname) cycle
      nr = nr + 1_I4P
      nd = (blk%rec(r)%face + 1_I4P) / 2_I4P ; d1 = 3_I4P - nd
      s = merge(1_I4P, -1_I4P, mod(blk%rec(r)%face, 2_I4P) == 1_I4P)
      n0 = size(line, kind=I4P)
      do m = blk%rec(r)%lo(d1), blk%rec(r)%hi(d1) - 1_I4P
         p = blk%rec(r)%lo(1:2) ; p(d1) = m ; q = p ; q(d1) = m + 1_I4P
         if (blk%u(p(1), p(2)) == blk%u(q(1), q(2)) .and. blk%v(p(1), p(2)) == blk%v(q(1), q(2))) cycle
         a = reshape([a, [blk%u(p(1), p(2)), blk%v(p(1), p(2))]], [2, size(line) + 1])
         b = reshape([b, [blk%u(q(1), q(2)), blk%v(q(1), q(2))]], [2, size(line) + 1])
         line = [line, nr]
      enddo
      n = size(line, kind=I4P)
      if (n == n0) cycle
      ! orientation from the first segment and the interior vertex next to its start
      p = blk%rec(r)%lo(1:2) ; p(d1) = blk%rec(r)%lo(d1) ; q = p ; q(nd) = q(nd) + s
      o = orient2d_xy(a(1, n0 + 1), a(2, n0 + 1), b(1, n0 + 1), b(2, n0 + 1), blk%u(q(1), q(2)), blk%v(q(1), q(2)))
      if (o > 0._R8P) then ! interior on the left: reverse the polyline
         a(:, n0 + 1:n) = a(:, n:n0 + 1:-1) ; b(:, n0 + 1:n) = b(:, n:n0 + 1:-1)
         call swap_ends(a(:, n0 + 1:n), b(:, n0 + 1:n))
      endif
   enddo
   sl%n = size(line, kind=I4P) ; sl%a = a ; sl%b = b ; sl%line = line ; allocate(sl%keep(sl%n)) ; sl%keep = .true.
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
