!< Alpha wrapping: a watertight, 2-manifold surrogate of any triangle soup (issue #18 §1.6, issue #50).
!<
!< Reference: Portaneri, Hemmer, Birdal, Mandad & Alliez, *Alpha Wrapping with an Offset* (SIGGRAPH 2022). CGAL
!< package: Alpha_wrap_3.
!<
!< `awrap_run` (behind `surface%alpha_wrap`) computes the wrap on a cell grid (issue #50): the boundary of the space a
!< ball of diameter α carves from outside, laid on the ε-offset surface of the input, extracted by Marching Cubes on a
!< well-composed binary lattice (closed, 2-manifold), with an exact inclusion check and a refinement loop. See its
!< header. (The first, octree pipeline of issue #18, an α-voxel shell, was removed by issue #54.)
module fossil_alpha_wrap
!< Alpha wrapping — watertight surrogate from broken triangle-soup input.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_marching_cubes, only : EDGE_VERTEX, VOX, MC_EDGE_TABLE, MC_TRI_TABLE
use fossil_tritri_exact, only : pair_classify, PAIR_SEPARATED, PAIR_CROSSING
use penf, only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none
private

public :: awrap_run
public :: AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE
public :: AWRAP_STATUS_NOT_CONVERGED
public :: AWRAP_DEFAULT_MAX_OUTER_ITERATIONS

integer(I4P), parameter :: AWRAP_STATUS_OK         = 0_I4P
integer(I4P), parameter :: AWRAP_STATUS_BAD_INPUT  = 1_I4P  !< Empty facet list, or alpha <= 0.
integer(I4P), parameter :: AWRAP_STATUS_DEGENERATE = 2_I4P  !< Alpha >= bbox diagonal: the wrap is empty.
integer(I4P), parameter :: AWRAP_STATUS_NOT_CONVERGED = 3_I4P  !< The finest grid level allowed still has a defect.

integer(I4P), parameter :: AWRAP_DEFAULT_MAX_OUTER_ITERATIONS = 5_I4P    !< Grid levels of awrap_run (cell size halved each).
! the wrap of awrap_run (issue #50)
real(R8P),    parameter :: AWRAP_CELLS_PER_ALPHA = 4._R8P    !< First level: cell size alpha/4.
real(R8P),    parameter :: AWRAP_MAX_CELLS = 8388608._R8P    !< Refinement stops before a grid of more cells (2^23).
integer(I4P), parameter :: AWRAP_CULL_BLOCK = 8_I4P          !< Distance culling block, cells per side.
real(R8P),    parameter :: AWRAP_EDT_FAR = 1.e20_R8P         !< No reached cell (squared index distance).
integer(I4P), parameter :: AWRAP_TRACE_STEPS = 64_I4P        !< Sphere-tracing steps to the offset surface.
real(R8P),    parameter :: AWRAP_TRACE_TOL = 1.e-9_R8P       !< Sphere tracing stops within this fraction of offset.
integer(I4P), parameter :: AWRAP_BISECTION_STEPS = 50_I4P    !< Bisection steps to the end of the carving balls.
real(R8P),    parameter :: AWRAP_FLOAT_ON = 1.e-13_R8P     !< float_pair: a vertex this close (relative) is on the plane.
real(R8P),    parameter :: AWRAP_FLOAT_BAND = 1.e-8_R8P    !< float_pair: below this (relative), the exact test decides.
integer(I4P), parameter :: PLANE_EXACT = 0_I4P, PLANE_APART = 1_I4P, PLANE_COPLANAR = 2_I4P !< plane_relation codes.
integer(I4P), parameter :: AWRAP_PUSH_ROUNDS = 6_I4P         !< Rounds of the repair pass (each halves the move).

type :: awrap_grid_t
   !< Cell grid of one level of awrap_run (issue #50). Cell (i, j, k) has centre o + (i - 1/2, j - 1/2, k - 1/2) h.
   type(vector_R8P)          :: o                  !< Lower corner of the grid.
   real(R8P)                 :: h    = 0._R8P      !< Cell size.
   real(R8P)                 :: rho  = 0._R8P      !< Radius of the carving ball (alpha/2).
   real(R8P)                 :: eps  = 0._R8P      !< Offset.
   integer(I4P)              :: n(3) = 0_I4P       !< Cells per axis.
   integer(I4P)              :: nvertex = 0_I4P    !< Wrap vertices.
   real(R8P),    allocatable :: u(:,:,:)           !< Unsigned distance at the centres (a lower bound where culled).
   logical,      allocatable :: reached(:,:,:)     !< Free and reached from outside (R).
   real(R8P),    allocatable :: dsq(:,:,:)         !< Power distance to the reached centres (index units), <= 0: covered.
   integer(I4P), allocatable :: site(:,:,:)        !< Linear index of the minimising reached cell.
   logical,      allocatable :: carved(:,:,:)      !< Carved (outside the wrap).
   integer(I4P), allocatable :: bin_first(:)       !< CSR offsets of the facets binned per cell.
   integer(I4P), allocatable :: bin_facet(:)       !< Binned facet ids.
endtype awrap_grid_t

contains

   pure subroutine compute_facet_bbox(facet, bmin, bmax)
   !< Tight axis-aligned bbox of all vertices in the facet array.
   type(facet_object), intent(in)  :: facet(:)
   type(vector_R8P),   intent(out) :: bmin, bmax
   integer(I4P)                    :: f, v

   bmin = facet(1)%vertex(1)
   bmax = facet(1)%vertex(1)
   do f = 1_I4P, size(facet, kind=I4P)
      do v = 1_I4P, 3_I4P
         bmin%x = min(bmin%x, facet(f)%vertex(v)%x)
         bmin%y = min(bmin%y, facet(f)%vertex(v)%y)
         bmin%z = min(bmin%z, facet(f)%vertex(v)%z)
         bmax%x = max(bmax%x, facet(f)%vertex(v)%x)
         bmax%y = max(bmax%y, facet(f)%vertex(v)%y)
         bmax%z = max(bmax%z, facet(f)%vertex(v)%z)
      enddo
   enddo
   endsubroutine compute_facet_bbox

   subroutine awrap_run(input_facet, input_tree, alpha, offset, max_iterations, wrapped_facets, status)
   !< Alpha wrap with an offset (issue #50): the boundary of the space an α-ball carves from outside, laid on the
   !< ε-offset surface of the input (Portaneri et al. 2022, on a cell grid instead of a Delaunay triangulation).
   !<
   !< One level, cell size h (α/4 first):
   !<   1. Unsigned distance u to the input at the cell centres (blocks far from the input are culled).
   !<   2. Free cells: u > α/2 + ε, where the carving ball (diameter α) fits outside the offset surface. A flood fill
   !<      from the grid corner over free cells gives the reached set R: the space the ball reaches from outside.
   !<   3. Carved cells: covered by the carving balls (within u(r) - ε of a reached centre r: an exact power-distance
   !<      transform, see grid_carve), and not touching the input (u > h√3/2). Every other cell is solid, so the input
   !<      always lies in solid cells, even a sheet carved on both sides: this is the offset-surface barrier, the two
   !<      sides of a wall are separate surfaces that never collapse onto each other.
   !<   4. The carved set is made well-composed (no 2-D checkerboard, no 3-D diagonal pair) by turning carved cells
   !<      solid; Marching Cubes on a well-composed binary lattice is a closed 2-manifold.
   !<   5. Marching Cubes on the cell centres. The vertex of a carved-solid edge is computed once per edge: the first
   !<      point from the carved centre where u = ε (sphere tracing: the offset surface). An edge that never comes
   !<      within ε of the input: to a solid cell touching the input, its midpoint projected onto the offset surface;
   !<      to one that does not (a cap over a gap narrower than α), the point where the carving balls end.
   !<   6. Repair: a wrap triangle that meets an input facet (exact `pair_classify` against the facets binned around
   !<      it), or crosses a wrap triangle near it (exact, with the adjacency rule of find_self_intersections), has
   !<      its vertices moved halfway back to their binary (edge-midpoint) positions, where the Marching Cubes surface
   !<      is embedded, for up to AWRAP_PUSH_ROUNDS rounds; a vertex that ends closer than ε to the input is projected
   !<      back onto the offset surface (issue #53). A chord across a convex edge or corner of the input, its
   !<      vertices on the offset surface, cuts the input when the cells are larger than ε; vertices converging on a
   !<      concave crease of the offset surface fold the triangles around them.
   !<
   !< Refinement loop: a level that leaves a triangle meeting the input, or crossing the wrap, is redone with
   !< h/2, up to `max_iterations` levels (default AWRAP_DEFAULT_MAX_OUTER_ITERATIONS) and AWRAP_MAX_CELLS cells.
   !< Status OK: the wrap encloses the input without touching it and does not cross itself; NOT_CONVERGED: the finest
   !< level still has a defect (its wrap is returned). The wrap is closed and 2-manifold in both cases. No vertex is
   !< closer than ε to the input: they lie on the ε-offset surface, except where the repair pass moved them outwards
   !< (convex edges and corners, concave creases) and on the caps (up to about α/2 + ε from the input).
   type(facet_object),              intent(in)            :: input_facet(:) !< Input facets (any topology).
   type(aabb_tree_object),          intent(in)            :: input_tree     !< AABB tree over input_facet.
   real(R8P),                       intent(in)            :: alpha          !< Diameter of the carving ball.
   real(R8P),                       intent(in)            :: offset         !< Offset ε of the wrap from the input.
   integer(I4P),                    intent(in),  optional :: max_iterations !< Refinement levels.
   type(facet_object), allocatable, intent(out)           :: wrapped_facets(:) !< Wrap.
   integer(I4P),                    intent(out), optional :: status         !< Status code.
   type(facet_object), allocatable                        :: trial(:)
   type(vector_R8P)                                       :: bmin, bmax
   real(R8P)                                              :: bdiag, h
   integer(I4P)                                           :: levels, level
   logical                                                :: contact

   if (present(status)) status = AWRAP_STATUS_OK
   if (size(input_facet, kind=I4P) == 0_I4P .or. alpha <= 0._R8P .or. offset <= 0._R8P) then
      allocate(wrapped_facets(0))
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif
   call compute_facet_bbox(facet=input_facet, bmin=bmin, bmax=bmax)
   bdiag = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2)
   if (alpha >= bdiag) then
      allocate(wrapped_facets(0))
      if (present(status)) status = AWRAP_STATUS_DEGENERATE
      return
   endif
   levels = AWRAP_DEFAULT_MAX_OUTER_ITERATIONS
   if (present(max_iterations)) levels = max(1_I4P, max_iterations)
   h = alpha / AWRAP_CELLS_PER_ALPHA
   do level = 1_I4P, levels
      if (level > 1_I4P .and. grid_cells(bmin=bmin, bmax=bmax, pad=grid_pad(alpha, offset, h), h=h) > AWRAP_MAX_CELLS) exit
      call wrap_level(input_facet=input_facet, input_tree=input_tree, bmin=bmin, bmax=bmax, alpha=alpha, &
                      offset=offset, h=h, wrap=trial, contact=contact)
      call move_alloc(from=trial, to=wrapped_facets)
      if (.not. contact) return
      h = 0.5_R8P * h
   enddo
   if (present(status)) status = AWRAP_STATUS_NOT_CONVERGED
   endsubroutine awrap_run

   pure function grid_pad(alpha, offset, h) result(pad)
   !< Padding of the grid around the input box: the two outermost cell layers are free and reached.
   real(R8P), intent(in) :: alpha, offset, h
   real(R8P)             :: pad

   pad = 0.5_R8P * alpha + offset + 2.5_R8P * h
   endfunction grid_pad

   pure function grid_cells(bmin, bmax, pad, h) result(n)
   !< Number of cells of the grid of size h over the padded box (real, to avoid overflow).
   type(vector_R8P), intent(in) :: bmin, bmax
   real(R8P),        intent(in) :: pad, h
   real(R8P)                    :: n

   n = real(ceiling((bmax%x - bmin%x + 2._R8P * pad) / h), R8P) * real(ceiling((bmax%y - bmin%y + 2._R8P * pad) / h), R8P) &
       * real(ceiling((bmax%z - bmin%z + 2._R8P * pad) / h), R8P)
   endfunction grid_cells

   subroutine wrap_level(input_facet, input_tree, bmin, bmax, alpha, offset, h, wrap, contact)
   !< One level of the wrap on cells of size h (see awrap_run).
   type(facet_object),              intent(in)  :: input_facet(:) !< Input facets.
   type(aabb_tree_object),          intent(in)  :: input_tree     !< AABB tree over input_facet.
   type(vector_R8P),                intent(in)  :: bmin, bmax     !< Input box.
   real(R8P),                       intent(in)  :: alpha, offset  !< Ball diameter, offset.
   real(R8P),                       intent(in)  :: h              !< Cell size.
   type(facet_object), allocatable, intent(out) :: wrap(:)        !< Wrap triangles.
   logical,                         intent(out) :: contact        !< A wrap triangle touches the input.
   type(awrap_grid_t)                           :: g
   real(R8P)                                    :: pad

   pad = grid_pad(alpha, offset, h)
   g%h = h
   g%rho = 0.5_R8P * alpha
   g%eps = offset
   g%o = vector_R8P(bmin%x - pad, bmin%y - pad, bmin%z - pad)
   g%n(1) = ceiling((bmax%x - bmin%x + 2._R8P * pad) / h)
   g%n(2) = ceiling((bmax%y - bmin%y + 2._R8P * pad) / h)
   g%n(3) = ceiling((bmax%z - bmin%z + 2._R8P * pad) / h)
   call grid_distances(g=g, input_facet=input_facet, input_tree=input_tree)
   call grid_reach(g=g)
   call grid_carve(g=g)
   call grid_well_composed(g=g)
   call grid_bin_facets(g=g, input_facet=input_facet)
   call grid_extract(g=g, input_facet=input_facet, input_tree=input_tree, wrap=wrap, contact=contact)
   endsubroutine wrap_level

   function udist(input_facet, input_tree, p) result(u)
   !< Unsigned distance from p to the input.
   type(facet_object),     intent(in) :: input_facet(:) !< Input facets.
   type(aabb_tree_object), intent(in) :: input_tree     !< AABB tree over input_facet.
   type(vector_R8P),       intent(in) :: p              !< Point.
   real(R8P)                          :: u
   real(R8P)                          :: d2
   integer(I4P)                       :: cf, cr

   call input_tree%distance_tree_with_region(facet=input_facet, point=p, distance=d2, closest_facet=cf, closest_region=cr)
   u = sqrt(max(d2, 0._R8P))
   endfunction udist

   pure function cell_centre(o, h, i, j, k) result(p)
   !< Centre of cell (i, j, k) of the grid with lower corner o and cell size h. The helpers take the grid's scalars,
   !< not the grid: gfortran finalises a derived-type dummy with allocatable components on every call.
   type(vector_R8P), intent(in) :: o
   real(R8P),        intent(in) :: h
   integer(I4P),     intent(in) :: i, j, k
   type(vector_R8P)             :: p

   p%x = o%x + (real(i, R8P) - 0.5_R8P) * h
   p%y = o%y + (real(j, R8P) - 0.5_R8P) * h
   p%z = o%z + (real(k, R8P) - 0.5_R8P) * h
   endfunction cell_centre

   pure function lin(n, i, j, k) result(l)
   !< Linear index of cell (i, j, k) of a grid of n cells per axis.
   integer(I4P), intent(in) :: n(3)
   integer(I4P), intent(in) :: i, j, k
   integer(I4P)             :: l

   l = i + n(1) * ((j - 1_I4P) + n(2) * (k - 1_I4P))
   endfunction lin

   pure subroutine unlin(n, l, i, j, k)
   !< Cell (i, j, k) of linear index l of a grid of n cells per axis.
   integer(I4P), intent(in)  :: n(3)
   integer(I4P), intent(in)  :: l
   integer(I4P), intent(out) :: i, j, k
   integer(I4P)              :: r

   r = l - 1_I4P
   i = mod(r, n(1)) + 1_I4P ; r = r / n(1)
   j = mod(r, n(2)) + 1_I4P
   k = r / n(2) + 1_I4P
   endsubroutine unlin

   subroutine grid_distances(g, input_facet, input_tree)
   !< Unsigned distance at the cell centres. A block of AWRAP_CULL_BLOCK^3 cells whose centres are all farther than
   !< the thresholds used (free: rho + eps, touching: h√3/2) gets the lower bound of its centre's distance.
   type(awrap_grid_t),     intent(inout) :: g              !< Grid.
   type(facet_object),     intent(in)    :: input_facet(:) !< Input facets.
   type(aabb_tree_object), intent(in)    :: input_tree     !< AABB tree over input_facet.
   integer(I4P)                          :: i0, j0, k0, i1, j1, k1, i, j, k
   real(R8P)                             :: ub, rb, thr
   type(vector_R8P)                      :: c0, c1, cb

   allocate(g%u(g%n(1), g%n(2), g%n(3)))
   thr = max(g%rho + g%eps, 0.5_R8P * sqrt(3._R8P) * g%h)
   do k0 = 1_I4P, g%n(3), AWRAP_CULL_BLOCK
      k1 = min(k0 + AWRAP_CULL_BLOCK - 1_I4P, g%n(3))
      do j0 = 1_I4P, g%n(2), AWRAP_CULL_BLOCK
         j1 = min(j0 + AWRAP_CULL_BLOCK - 1_I4P, g%n(2))
         do i0 = 1_I4P, g%n(1), AWRAP_CULL_BLOCK
            i1 = min(i0 + AWRAP_CULL_BLOCK - 1_I4P, g%n(1))
            c0 = cell_centre(g%o, g%h, i0, j0, k0) ; c1 = cell_centre(g%o, g%h, i1, j1, k1)
            cb = vector_R8P(0.5_R8P * (c0%x + c1%x), 0.5_R8P * (c0%y + c1%y), 0.5_R8P * (c0%z + c1%z))
            rb = 0.5_R8P * sqrt((c1%x - c0%x)**2 + (c1%y - c0%y)**2 + (c1%z - c0%z)**2)
            ub = udist(input_facet, input_tree, cb)
            if (ub - rb > thr) then
               g%u(i0:i1, j0:j1, k0:k1) = ub - rb
               cycle
            endif
            do k = k0, k1
               do j = j0, j1
                  do i = i0, i1
                     g%u(i, j, k) = udist(input_facet, input_tree, cell_centre(g%o, g%h, i, j, k))
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine grid_distances

   subroutine grid_reach(g)
   !< Flood fill from cell (1, 1, 1) over the free cells (u > rho + eps), 6-connected: the reached set R.
   type(awrap_grid_t), intent(inout) :: g !< Grid.
   integer(I4P), allocatable         :: queue(:)
   integer(I4P)                      :: head, tail, i, j, k, a, ni, nj, nk
   integer(I4P), parameter           :: STEP(3, 6) = reshape([-1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1], &
                                                             [3, 6])
   real(R8P)                         :: rfree

   rfree = g%rho + g%eps
   allocate(g%reached(g%n(1), g%n(2), g%n(3)))
   g%reached = .false.
   allocate(queue(g%n(1) * g%n(2) * g%n(3)))
   head = 1_I4P ; tail = 1_I4P
   queue(1) = lin(g%n, 1_I4P, 1_I4P, 1_I4P)
   g%reached(1, 1, 1) = .true.
   do while (head <= tail)
      call unlin(g%n, queue(head), i, j, k)
      head = head + 1_I4P
      do a = 1_I4P, 6_I4P
         ni = i + STEP(1, a) ; nj = j + STEP(2, a) ; nk = k + STEP(3, a)
         if (ni < 1_I4P .or. nj < 1_I4P .or. nk < 1_I4P) cycle
         if (ni > g%n(1) .or. nj > g%n(2) .or. nk > g%n(3)) cycle
         if (g%reached(ni, nj, nk)) cycle
         if (g%u(ni, nj, nk) <= rfree) cycle
         g%reached(ni, nj, nk) = .true.
         tail = tail + 1_I4P
         queue(tail) = lin(g%n, ni, nj, nk)
      enddo
   enddo
   endsubroutine grid_reach

   subroutine grid_carve(g)
   !< Carved cells: covered by the carving balls and not touching the input (u > h√3/2). Every point within
   !< u(r) - eps of a reached centre r is covered by a ball of diameter alpha whose centre is reachable from r inside
   !< the free space (the centre on the segment from the point to r, alpha/2 from the point, has u > alpha/2 + eps),
   !< so a cell c is carved if min over r of |c - r|^2 - (u(r) - eps)^2 <= 0: an exact power-distance transform,
   !< separable like the squared EDT (Felzenszwalb & Huttenlocher 2012) with f(r) = -(u(r) - eps)^2, which also keeps
   !< the minimising centre (site).
   type(awrap_grid_t), intent(inout) :: g !< Grid.
   real(R8P), allocatable            :: f(:), d(:)
   integer(I4P), allocatable         :: arg(:), col(:)
   integer(I4P)                      :: i, j, k, m
   real(R8P)                         :: rtouch

   allocate(g%dsq(g%n(1), g%n(2), g%n(3)), g%site(g%n(1), g%n(2), g%n(3)))
   m = maxval(g%n)
   allocate(f(m), d(m), arg(m), col(m))
   ! x pass
   do k = 1_I4P, g%n(3)
      do j = 1_I4P, g%n(2)
         do i = 1_I4P, g%n(1)
            f(i) = AWRAP_EDT_FAR
            if (g%reached(i, j, k)) f(i) = -((g%u(i, j, k) - g%eps) / g%h)**2
         enddo
         call edt_1d(f=f(1:g%n(1)), d=d(1:g%n(1)), arg=arg(1:g%n(1)))
         do i = 1_I4P, g%n(1)
            g%dsq(i, j, k) = d(i) ; g%site(i, j, k) = lin(g%n, arg(i), j, k)
         enddo
      enddo
   enddo
   ! y pass
   do k = 1_I4P, g%n(3)
      do i = 1_I4P, g%n(1)
         f(1:g%n(2)) = g%dsq(i, :, k) ; col(1:g%n(2)) = g%site(i, :, k)
         call edt_1d(f=f(1:g%n(2)), d=d(1:g%n(2)), arg=arg(1:g%n(2)))
         do j = 1_I4P, g%n(2)
            g%dsq(i, j, k) = d(j) ; g%site(i, j, k) = col(arg(j))
         enddo
      enddo
   enddo
   ! z pass
   do j = 1_I4P, g%n(2)
      do i = 1_I4P, g%n(1)
         f(1:g%n(3)) = g%dsq(i, j, :) ; col(1:g%n(3)) = g%site(i, j, :)
         call edt_1d(f=f(1:g%n(3)), d=d(1:g%n(3)), arg=arg(1:g%n(3)))
         do k = 1_I4P, g%n(3)
            g%dsq(i, j, k) = d(k) ; g%site(i, j, k) = col(arg(k))
         enddo
      enddo
   enddo
   rtouch = 0.5_R8P * sqrt(3._R8P) * g%h
   allocate(g%carved(g%n(1), g%n(2), g%n(3)))
   g%carved = (g%dsq <= 0._R8P) .and. (g%u > rtouch)
   endsubroutine grid_carve

   pure subroutine edt_1d(f, d, arg)
   !< 1-D squared distance transform d(q) = min_p (q - p)^2 + f(p), with the minimiser arg(q) (lower envelope of
   !< parabolas, Felzenszwalb & Huttenlocher). f >= AWRAP_EDT_FAR marks no site.
   real(R8P),    intent(in)  :: f(:)   !< Sampled function.
   real(R8P),    intent(out) :: d(:)   !< Transform.
   integer(I4P), intent(out) :: arg(:) !< Minimiser.
   integer(I4P)              :: v(size(f)), n, k, q
   real(R8P)                 :: z(size(f) + 1), s

   n = size(f, kind=I4P)
   k = 1_I4P ; v(1) = 1_I4P ; z(1) = -huge(1._R8P) ; z(2) = huge(1._R8P)
   do q = 2_I4P, n
      do
         s = ((f(q) + real(q, R8P)**2) - (f(v(k)) + real(v(k), R8P)**2)) / (2._R8P * real(q - v(k), R8P))
         if (s > z(k)) exit
         k = k - 1_I4P
      enddo
      k = k + 1_I4P
      v(k) = q ; z(k) = s ; z(k + 1) = huge(1._R8P)
   enddo
   k = 1_I4P
   do q = 1_I4P, n
      do while (z(k + 1) < real(q, R8P))
         k = k + 1_I4P
      enddo
      d(q) = real(q - v(k), R8P)**2 + f(v(k)) ; arg(q) = v(k)
   enddo
   endsubroutine edt_1d

   subroutine grid_well_composed(g)
   !< Turn carved cells solid until the carved set is well-composed: no 2x2 square of cells with a carved diagonal
   !< and a solid one, no 2x2x2 block whose only two carved (or only two solid) cells are opposite corners. Monotone:
   !< carved cells only decrease, so it ends; the outer layers (free, far from the input) never change.
   type(awrap_grid_t), intent(inout) :: g !< Grid.
   integer(I4P)                      :: i, j, k, c, nc, a, b
   logical                           :: changed, q(0:7)
   integer(I4P), parameter           :: OPP(0:7) = [6, 7, 4, 5, 2, 3, 0, 1] ! opposite corner, VOX numbering

   do
      changed = .false.
      do k = 1_I4P, g%n(3) - 1_I4P
         do j = 1_I4P, g%n(2) - 1_I4P
            do i = 1_I4P, g%n(1) - 1_I4P
               ! the 2x2 squares of the block's three lower faces (the upper ones are the next blocks' lower faces,
               ! except on the last layer: those squares lie in the outer, all-carved layers)
               call square(i, j, k, i + 1_I4P, j, k, i, j + 1_I4P, k, i + 1_I4P, j + 1_I4P, k)
               call square(i, j, k, i + 1_I4P, j, k, i, j, k + 1_I4P, i + 1_I4P, j, k + 1_I4P)
               call square(i, j, k, i, j + 1_I4P, k, i, j, k + 1_I4P, i, j + 1_I4P, k + 1_I4P)
               do c = 0_I4P, 7_I4P
                  q(c) = g%carved(i + VOX(0, c), j + VOX(1, c), k + VOX(2, c))
               enddo
               nc = count(q)
               if (nc /= 2_I4P .and. nc /= 6_I4P) cycle
               do c = 0_I4P, 7_I4P
                  if (q(c) .eqv. (nc == 6_I4P)) cycle ! the minority corners
                  if (q(OPP(c)) .eqv. q(c)) then
                     if (nc == 2_I4P) then
                        g%carved(i + VOX(0, c), j + VOX(1, c), k + VOX(2, c)) = .false.
                        a = OPP(c)
                        g%carved(i + VOX(0, a), j + VOX(1, a), k + VOX(2, a)) = .false.
                     else
                        b = ieor(c, 1_I4P) ! the x-neighbour of the solid corner: carved, turned solid
                        g%carved(i + VOX(0, b), j + VOX(1, b), k + VOX(2, b)) = .false.
                     endif
                     changed = .true.
                  endif
                  exit
               enddo
            enddo
         enddo
      enddo
      if (.not. changed) exit
   enddo
   contains
      subroutine square(i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4)
      !< Cells 1-4 of a square: 1 and 4 are one diagonal, 2 and 3 the other.
      integer(I4P), intent(in) :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
      logical                  :: c1, c2, c3, c4

      c1 = g%carved(i1, j1, k1) ; c2 = g%carved(i2, j2, k2) ; c3 = g%carved(i3, j3, k3) ; c4 = g%carved(i4, j4, k4)
      if (c1 .and. c4 .and. .not. c2 .and. .not. c3) then
         g%carved(i1, j1, k1) = .false. ; g%carved(i4, j4, k4) = .false. ; changed = .true.
      elseif (c2 .and. c3 .and. .not. c1 .and. .not. c4) then
         g%carved(i2, j2, k2) = .false. ; g%carved(i3, j3, k3) = .false. ; changed = .true.
      endif
      endsubroutine square
   endsubroutine grid_well_composed

   subroutine grid_bin_facets(g, input_facet)
   !< Bin the input facets into the cells their bounding box overlaps (CSR lists: g%bin_first, g%bin_facet).
   type(awrap_grid_t), intent(inout) :: g              !< Grid.
   type(facet_object), intent(in)    :: input_facet(:) !< Input facets.
   integer(I4P)                      :: f, i, j, k, lo(3), hi(3), ncell, pass, l
   integer(I4P), allocatable         :: fill(:)

   ncell = g%n(1) * g%n(2) * g%n(3)
   allocate(g%bin_first(ncell + 1_I4P))
   g%bin_first = 0_I4P
   do pass = 1_I4P, 2_I4P
      do f = 1_I4P, size(input_facet, kind=I4P)
         call facet_cells(input_facet(f), lo, hi)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  l = lin(g%n, i, j, k)
                  if (pass == 1_I4P) then
                     g%bin_first(l + 1_I4P) = g%bin_first(l + 1_I4P) + 1_I4P
                  else
                     fill(l) = fill(l) + 1_I4P
                     g%bin_facet(fill(l)) = f
                  endif
               enddo
            enddo
         enddo
      enddo
      if (pass == 1_I4P) then
         g%bin_first(1) = 1_I4P
         do l = 1_I4P, ncell
            g%bin_first(l + 1_I4P) = g%bin_first(l + 1_I4P) + g%bin_first(l)
         enddo
         allocate(g%bin_facet(g%bin_first(ncell + 1_I4P) - 1_I4P))
         allocate(fill(ncell))
         fill = g%bin_first(1:ncell) - 1_I4P
      endif
   enddo
   contains
      subroutine facet_cells(fc, lo, hi)
      !< Cell index range of a facet's bounding box.
      type(facet_object), intent(in)  :: fc
      integer(I4P),       intent(out) :: lo(3), hi(3)
      real(R8P)                       :: vmin(3), vmax(3), o(3)
      integer(I4P)                    :: v, a

      vmin = huge(1._R8P) ; vmax = -huge(1._R8P)
      do v = 1_I4P, 3_I4P
         vmin = min(vmin, [fc%vertex(v)%x, fc%vertex(v)%y, fc%vertex(v)%z])
         vmax = max(vmax, [fc%vertex(v)%x, fc%vertex(v)%y, fc%vertex(v)%z])
      enddo
      o = [g%o%x, g%o%y, g%o%z]
      do a = 1_I4P, 3_I4P
         lo(a) = max(1_I4P, min(g%n(a), floor((vmin(a) - o(a)) / g%h, I4P) + 1_I4P))
         hi(a) = max(1_I4P, min(g%n(a), floor((vmax(a) - o(a)) / g%h, I4P) + 1_I4P))
      enddo
      endsubroutine facet_cells
   endsubroutine grid_bin_facets

   subroutine grid_extract(g, input_facet, input_tree, wrap, contact)
   !< Marching Cubes on the cell centres, solid = inside (MC bit set); one vertex per carved-solid edge, computed once
   !< (open-addressing hash on the edge key). Then the repair pass of awrap_run (step 6). `contact` reports a triangle
   !< still meeting the input or crossing the wrap after AWRAP_PUSH_ROUNDS rounds.
   type(awrap_grid_t),              intent(inout) :: g              !< Grid.
   type(facet_object),              intent(in)    :: input_facet(:) !< Input facets.
   type(aabb_tree_object),          intent(in)    :: input_tree     !< AABB tree over input_facet.
   type(facet_object), allocatable, intent(out)   :: wrap(:)        !< Wrap triangles.
   logical,                         intent(out)   :: contact        !< A defect is left.
   integer(I8P), allocatable                      :: key(:)
   type(vector_R8P), allocatable                  :: point(:), mid(:)
   integer(I4P), allocatable                      :: slot_vertex(:), tv(:,:), tcube(:), itmp(:,:), ctmp(:)
   integer(I4P), allocatable                      :: cube_first(:), cube_tri(:), stamp(:)
   real(R8P), allocatable                         :: tbox(:,:)      !< Per triangle: min x, y, z, max x, y, z.
   logical, allocatable                           :: tdeg(:)        !< Per triangle: two coincident vertices.
   logical, allocatable                           :: defect(:), push(:), moved(:)
   integer(I4P)                                   :: ev(0:11), nedge, nslot, i, j, k, a, c, e, t, cube_case, mask
   integer(I4P)                                   :: c0, c1, ci(3), lo(3), idx, ntri, cap, round, v, t1
   integer(I8P)                                   :: kk
   type(vector_R8P)                               :: tri(3), cp
   real(R8P)                                      :: u

   contact = .false.
   ! count the carved-solid edges, size the hash
   nedge = 0_I4P
   do k = 1_I4P, g%n(3)
      do j = 1_I4P, g%n(2)
         do i = 1_I4P, g%n(1)
            if (i < g%n(1)) then ; if (g%carved(i, j, k) .neqv. g%carved(i + 1, j, k)) nedge = nedge + 1_I4P ; endif
            if (j < g%n(2)) then ; if (g%carved(i, j, k) .neqv. g%carved(i, j + 1, k)) nedge = nedge + 1_I4P ; endif
            if (k < g%n(3)) then ; if (g%carved(i, j, k) .neqv. g%carved(i, j, k + 1)) nedge = nedge + 1_I4P ; endif
         enddo
      enddo
   enddo
   nslot = 64_I4P
   do while (nslot < 2_I4P * nedge)
      nslot = 2_I4P * nslot
   enddo
   allocate(key(0:nslot - 1), slot_vertex(0:nslot - 1), point(max(1_I4P, nedge)), mid(max(1_I4P, nedge)))
   key = -1_I8P
   g%nvertex = 0_I4P
   cap = max(1024_I4P, 2_I4P * nedge) ; ntri = 0_I4P
   allocate(tv(3, cap), tcube(cap))
   ! marching cubes: vertices and triangles (vertex ids)
   do k = 1_I4P, g%n(3) - 1_I4P
      do j = 1_I4P, g%n(2) - 1_I4P
         do i = 1_I4P, g%n(1) - 1_I4P
            cube_case = 0_I4P
            do c = 0_I4P, 7_I4P
               if (.not. g%carved(i + VOX(0, c), j + VOX(1, c), k + VOX(2, c))) cube_case = ibset(cube_case, c)
            enddo
            mask = MC_EDGE_TABLE(cube_case)
            if (mask == 0_I4P) cycle
            do e = 0_I4P, 11_I4P
               if (.not. btest(mask, e)) cycle
               c0 = EDGE_VERTEX(1, e) ; c1 = EDGE_VERTEX(2, e)
               if (sum(VOX(:, c0)) > sum(VOX(:, c1))) then ; c0 = EDGE_VERTEX(2, e) ; c1 = EDGE_VERTEX(1, e) ; endif
               lo = [i + VOX(0, c0), j + VOX(1, c0), k + VOX(2, c0)]
               a = maxloc(abs(VOX(:, c1) - VOX(:, c0)), dim=1)
               kk = 3_I8P * int(lin(g%n, lo(1), lo(2), lo(3)) - 1_I4P, I8P) + int(a - 1_I4P, I8P)
               idx = lookup(kk)
               if (idx == 0_I4P) then
                  ci = lo ; ci(a) = ci(a) + 1_I4P
                  g%nvertex = g%nvertex + 1_I4P
                  if (g%carved(lo(1), lo(2), lo(3))) then
                     point(g%nvertex) = edge_point(g, input_facet, input_tree, lo, ci)
                  else
                     point(g%nvertex) = edge_point(g, input_facet, input_tree, ci, lo)
                  endif
                  mid(g%nvertex) = 0.5_R8P * (cell_centre(g%o, g%h, lo(1), lo(2), lo(3)) + &
                                              cell_centre(g%o, g%h, ci(1), ci(2), ci(3)))
                  call insert(kk, g%nvertex)
                  idx = g%nvertex
               endif
               ev(e) = idx
            enddo
            t = 0_I4P
            do
               if (t > 13_I4P) exit
               if (MC_TRI_TABLE(t, cube_case) < 0_I4P) exit
               if (ntri == cap) then
                  cap = 2_I4P * cap
                  allocate(itmp(3, cap), ctmp(cap))
                  itmp(:, 1:ntri) = tv(:, 1:ntri) ; ctmp(1:ntri) = tcube(1:ntri)
                  call move_alloc(from=itmp, to=tv) ; call move_alloc(from=ctmp, to=tcube)
               endif
               ntri = ntri + 1_I4P
               ! Bourke's listing reversed (as in extract_isosurface): outward normals for FOSSIL's convention
               tv(1, ntri) = ev(MC_TRI_TABLE(t, cube_case))
               tv(3, ntri) = ev(MC_TRI_TABLE(t + 1_I4P, cube_case))
               tv(2, ntri) = ev(MC_TRI_TABLE(t + 2_I4P, cube_case))
               tcube(ntri) = lin(g%n, i, j, k)
               t = t + 3_I4P
            enddo
         enddo
      enddo
   enddo
   ! repair pass: a triangle that meets the input, or crosses a wrap triangle near it, has its vertices (and the other
   ! triangle's) moved halfway back to their binary positions, where the Marching Cubes surface is embedded
   allocate(defect(ntri), push(g%nvertex), moved(g%nvertex), stamp(ntri), tbox(6, ntri), tdeg(ntri))
   defect = .false.
   moved = .true.
   do round = 0_I4P, AWRAP_PUSH_ROUNDS
      call cube_lists
      stamp = 0_I4P
      push = .false.
      do t = 1_I4P, ntri
         if (.not. any(moved(tv(:, t)))) cycle
         defect(t) = .false.
         if (degenerate(point(tv(:, t)))) cycle
         if (touches_input(tcube(t), point(tv(:, t)))) then
            defect(t) = .true.
         else
            t1 = crossed_triangle(t)
            if (t1 > 0_I4P) then
               defect(t) = .true.
               push(tv(:, t1)) = .true.
            endif
         endif
         if (defect(t)) push(tv(:, t)) = .true.
      enddo
      if (.not. any(push) .or. round == AWRAP_PUSH_ROUNDS) exit
      do v = 1_I4P, g%nvertex
         if (.not. push(v)) cycle
         point(v) = 0.5_R8P * (point(v) + mid(v))
         ! back onto the offset surface if the move went closer to the input than eps (issue #53)
         call closest_point(input_facet, input_tree, point(v), cp, u)
         if (u < g%eps .and. u > 0._R8P) point(v) = cp + (g%eps / u) * (point(v) - cp)
      enddo
      moved = push
   enddo
   contact = any(defect)
   ! the facets (a triangle with coincident vertices has no area: dropped, both sides of its only edge go with it)
   allocate(wrap(ntri))
   c = 0_I4P
   do t = 1_I4P, ntri
      if (tv(1, t) == tv(2, t) .or. tv(2, t) == tv(3, t) .or. tv(3, t) == tv(1, t)) cycle
      tri = point(tv(:, t))
      if (degenerate(tri)) cycle
      c = c + 1_I4P
      wrap(c)%vertex = tri
      call wrap(c)%compute_metrix
   enddo
   if (c < ntri) wrap = wrap(1:c)
   contains
      function slot_of(kk) result(s)
      !< First probe slot of a key (multiplicative hash).
      integer(I8P), intent(in) :: kk
      integer(I4P)             :: s

      s = int(iand(ishft(kk * 2654435761_I8P, -7), int(nslot - 1_I4P, I8P)), I4P)
      endfunction slot_of

      function lookup(kk) result(v)
      !< Vertex of a key, 0 if absent.
      integer(I8P), intent(in) :: kk
      integer(I4P)             :: v, s

      v = 0_I4P
      s = slot_of(kk)
      do while (key(s) /= -1_I8P)
         if (key(s) == kk) then ; v = slot_vertex(s) ; return ; endif
         s = iand(s + 1_I4P, nslot - 1_I4P)
      enddo
      endfunction lookup

      subroutine insert(kk, v)
      !< Insert a key (absent).
      integer(I8P), intent(in) :: kk
      integer(I4P), intent(in) :: v
      integer(I4P)             :: s

      s = slot_of(kk)
      do while (key(s) /= -1_I8P)
         s = iand(s + 1_I4P, nslot - 1_I4P)
      enddo
      key(s) = kk ; slot_vertex(s) = v
      endsubroutine insert

      subroutine cube_lists
      !< Triangles per dual cube, by bounding box (CSR over the linear index of the cube's lower cell): a triangle is
      !< listed in every dual cube its box overlaps, so the cubes overlapping a triangle's box list every triangle
      !< that can meet it. Rebuilt every repair round (the vertices move).
      integer(I4P)              :: ncell, l, tt, pass, lo(3), hi(3), ii, jj, kk3
      integer(I4P), allocatable :: fill(:)

      ncell = g%n(1) * g%n(2) * g%n(3)
      if (allocated(cube_first)) deallocate(cube_first)
      if (allocated(cube_tri)) deallocate(cube_tri)
      allocate(cube_first(ncell + 1_I4P), fill(ncell))
      cube_first = 0_I4P
      do tt = 1_I4P, ntri
         tri = point(tv(:, tt))
         tdeg(tt) = degenerate(tri)
         tbox(:, tt) = [min(tri(1)%x, tri(2)%x, tri(3)%x), min(tri(1)%y, tri(2)%y, tri(3)%y), &
                        min(tri(1)%z, tri(2)%z, tri(3)%z), max(tri(1)%x, tri(2)%x, tri(3)%x), &
                        max(tri(1)%y, tri(2)%y, tri(3)%y), max(tri(1)%z, tri(2)%z, tri(3)%z)]
      enddo
      do pass = 1_I4P, 2_I4P
         do tt = 1_I4P, ntri
            call dual_range(tbox(:, tt), lo, hi)
            do kk3 = lo(3), hi(3)
               do jj = lo(2), hi(2)
                  do ii = lo(1), hi(1)
                     l = lin(g%n, ii, jj, kk3)
                     if (pass == 1_I4P) then
                        cube_first(l + 1_I4P) = cube_first(l + 1_I4P) + 1_I4P
                     else
                        fill(l) = fill(l) + 1_I4P
                        cube_tri(fill(l)) = tt
                     endif
                  enddo
               enddo
            enddo
         enddo
         if (pass == 1_I4P) then
            cube_first(1) = 1_I4P
            do l = 1_I4P, ncell
               cube_first(l + 1_I4P) = cube_first(l + 1_I4P) + cube_first(l)
            enddo
            allocate(cube_tri(max(1_I4P, cube_first(ncell + 1_I4P) - 1_I4P)))
            fill = cube_first(1:ncell) - 1_I4P
         endif
      enddo
      endsubroutine cube_lists

      subroutine dual_range(box, lo, hi)
      !< Range of the dual cubes (by lower cell) a bounding box overlaps: cube i spans the centres i, i+1.
      real(R8P),    intent(in)  :: box(6)
      integer(I4P), intent(out) :: lo(3), hi(3)
      real(R8P)                 :: o(3)
      integer(I4P)              :: a

      o = [g%o%x, g%o%y, g%o%z]
      do a = 1_I4P, 3_I4P
         lo(a) = max(1_I4P, min(g%n(a) - 1_I4P, floor((box(a) - o(a)) / g%h - 0.5_R8P, I4P) + 1_I4P))
         hi(a) = max(1_I4P, min(g%n(a) - 1_I4P, floor((box(a + 3) - o(a)) / g%h - 0.5_R8P, I4P) + 1_I4P))
      enddo
      endsubroutine dual_range

      function crossed_triangle(t0) result(t1)
      !< A wrap triangle listed in the dual cubes t0's box overlaps that t0 crosses (0 if none): the two share a segment of
      !< positive length (exact `pair_classify`). Two or three vertices in common (bit for bit) is adjacency, as in
      !< find_self_intersections (issue #46); a (near-)coplanar pair crosses if it overlaps with positive area.
      integer(I4P), intent(in) :: t0
      integer(I4P)             :: t1
      integer(I4P)             :: ii, jj, kk3, l, b, kind, ns, x, y, tc, rel, lo(3), hi(3)
      type(vector_R8P)         :: ta(3), tb(3), p, q

      t1 = 0_I4P
      ta = point(tv(:, t0))
      call dual_range(tbox(:, t0), lo, hi)
      do kk3 = lo(3), hi(3)
         do jj = lo(2), hi(2)
            do ii = lo(1), hi(1)
               l = lin(g%n, ii, jj, kk3)
               do b = cube_first(l), cube_first(l + 1_I4P) - 1_I4P
                  tc = cube_tri(b)
                  if (tc == t0 .or. stamp(tc) == t0) cycle
                  stamp(tc) = t0
                  if (tdeg(tc)) cycle
                  if (any(tbox(1:3, tc) > tbox(4:6, t0)) .or. any(tbox(1:3, t0) > tbox(4:6, tc))) cycle
                  tb = point(tv(:, tc))
                  ns = 0_I4P
                  do x = 1_I4P, 3_I4P
                     do y = 1_I4P, 3_I4P
                        if (same_xyz(ta(x), tb(y))) ns = ns + 1_I4P
                     enddo
                  enddo
                  if (ns >= 2_I4P) cycle
                  rel = plane_relation(ta, tb)
                  if (rel == PLANE_APART) cycle
                  if (rel == PLANE_COPLANAR) then
                     if (coplanar_overlap(ta, tb)) then ; t1 = tc ; return ; endif
                     cycle
                  endif
                  if (plane_relation(tb, ta) == PLANE_APART) cycle
                  kind = float_pair(ta, tb)
                  if (kind == 1_I4P) then ; t1 = tc ; return ; endif
                  if (kind == 0_I4P .or. kind == 2_I4P) cycle ! a point contact is not a crossing
                  call pair_classify(ta=ta, tb=tb, kind=kind, p=p, q=q)
                  if (kind == PAIR_CROSSING) then
                     t1 = tc
                     return
                  endif
               enddo
            enddo
         enddo
      enddo
      endfunction crossed_triangle

      function touches_input(cube, tri) result(yes)
      !< Does the triangle of dual cube `cube` (linear index of its lower cell) meet an input facet binned in its 8
      !< cells (exact orient3d)? The pushed vertices stay within the cells' reach of their facets: the bins of the
      !< 8 cells hold every facet the dual cube's triangles can meet, a pushed one included only while it stays in it,
      !< so the bins of the 27 cells around the lower one are searched.
      integer(I4P),     intent(in) :: cube
      type(vector_R8P), intent(in) :: tri(3)
      logical                      :: yes
      integer(I4P)                 :: ii, jj, kk3, di, dj, dk, l, b, kind
      type(vector_R8P)             :: p, q

      yes = .false.
      call unlin(g%n, cube, ii, jj, kk3)
      do dk = -1_I4P, 2_I4P
         do dj = -1_I4P, 2_I4P
            do di = -1_I4P, 2_I4P
               if (ii + di < 1_I4P .or. jj + dj < 1_I4P .or. kk3 + dk < 1_I4P) cycle
               if (ii + di > g%n(1) .or. jj + dj > g%n(2) .or. kk3 + dk > g%n(3)) cycle
               l = lin(g%n, ii + di, jj + dj, kk3 + dk)
               do b = g%bin_first(l), g%bin_first(l + 1_I4P) - 1_I4P
                  if (.not. boxes_overlap(tri, input_facet(g%bin_facet(b))%vertex)) cycle
                  if (plane_relation(tri, input_facet(g%bin_facet(b))%vertex) == PLANE_APART) cycle
                  if (plane_relation(input_facet(g%bin_facet(b))%vertex, tri) == PLANE_APART) cycle
                  kind = float_pair(tri, input_facet(g%bin_facet(b))%vertex)
                  if (kind == 1_I4P) then ; yes = .true. ; return ; endif
                  if (kind == 0_I4P) cycle
                  call pair_classify(ta=tri, tb=input_facet(g%bin_facet(b))%vertex, kind=kind, p=p, q=q)
                  if (kind /= PAIR_SEPARATED) then ; yes = .true. ; return ; endif
               enddo
            enddo
         enddo
      enddo
      endfunction touches_input
   endsubroutine grid_extract

   function edge_point(g, input_facet, input_tree, cc, cs) result(p)
   !< Wrap vertex of the edge from the carved cell cc to the solid cell cs: the first point from cc where u = eps
   !< (sphere tracing, which never passes the first root: on cc's side of any wall; near-tangent edges, where it
   !< converges slowly, end with a projection of the last point onto the offset surface). If the edge never comes
   !< within eps: cs touches the input, the edge midpoint projected onto the offset surface along its closest-point
   !< direction (the midpoint is on cc's side: the edge did not cross the input); otherwise (a cap over a gap narrower
   !< than alpha) where the union of the carving balls ends (bisection on the distance to the nearest reached
   !< centres), else the midpoint.
   type(awrap_grid_t),     intent(in) :: g              !< Grid.
   type(facet_object),     intent(in) :: input_facet(:) !< Input facets.
   type(aabb_tree_object), intent(in) :: input_tree     !< AABB tree over input_facet.
   integer(I4P),           intent(in) :: cc(3), cs(3)   !< Carved and solid cells.
   type(vector_R8P)                   :: p              !< Vertex.
   type(vector_R8P)                   :: xc, xs, xm, cp
   real(R8P)                          :: t, gap, ta, tb, tm, us
   integer(I4P)                       :: it

   xc = cell_centre(g%o, g%h, cc(1), cc(2), cc(3)) ; xs = cell_centre(g%o, g%h, cs(1), cs(2), cs(3))
   t = 0._R8P
   do it = 1_I4P, AWRAP_TRACE_STEPS
      p = along(t)
      gap = udist(input_facet, input_tree, p) - g%eps
      if (gap <= AWRAP_TRACE_TOL * g%eps) return
      t = t + gap / g%h
      if (t >= 1._R8P) exit
   enddo
   if (t < 1._R8P) then ! grazing: stopped short of the offset surface by gap, on the carved side: projected
      call closest_point(input_facet, input_tree, p, cp, us)
      p%x = cp%x + (g%eps / us) * (p%x - cp%x)
      p%y = cp%y + (g%eps / us) * (p%y - cp%y)
      p%z = cp%z + (g%eps / us) * (p%z - cp%z)
      return
   endif
   if (g%u(cs(1), cs(2), cs(3)) <= 0.5_R8P * sqrt(3._R8P) * g%h) then
      ! the solid cell touches the input but the edge never comes within eps: its midpoint, on cc's side (the edge
      ! did not cross the input), projected onto the offset surface
      xm = along(0.5_R8P)
      call closest_point(input_facet, input_tree, xm, cp, us)
      p%x = cp%x + (g%eps / us) * (xm%x - cp%x)
      p%y = cp%y + (g%eps / us) * (xm%y - cp%y)
      p%z = cp%z + (g%eps / us) * (xm%z - cp%z)
      return
   endif
   if (ball_gap(xc) <= 0._R8P .and. ball_gap(xs) > 0._R8P) then
      ta = 0._R8P ; tb = 1._R8P
      do it = 1_I4P, AWRAP_BISECTION_STEPS
         tm = 0.5_R8P * (ta + tb)
         if (ball_gap(along(tm)) <= 0._R8P) then ; ta = tm ; else ; tb = tm ; endif
      enddo
      p = along(0.5_R8P * (ta + tb))
   else
      p = along(0.5_R8P)
   endif
   contains
      function along(s) result(x)
      !< Point at parameter s from xc to xs.
      real(R8P), intent(in) :: s
      type(vector_R8P)      :: x

      x%x = xc%x + s * (xs%x - xc%x) ; x%y = xc%y + s * (xs%y - xc%y) ; x%z = xc%z + s * (xs%z - xc%z)
      endfunction along

      function ball_gap(x) result(b)
      !< Signed gap from x to the carving balls (distance to a reached centre r minus u(r) - eps), over the minimising
      !< sites of both cells and their face neighbours.
      type(vector_R8P), intent(in) :: x
      real(R8P)                    :: b
      integer(I4P)                 :: w, a, sgn, nb(3), si, sj, sk
      type(vector_R8P)             :: r

      b = huge(1._R8P)
      do w = 1_I4P, 2_I4P
         do a = 0_I4P, 3_I4P
            do sgn = -1_I4P, 1_I4P, 2_I4P
               if (a == 0_I4P .and. sgn == 1_I4P) cycle
               nb = merge(cc, cs, w == 1_I4P)
               if (a > 0_I4P) nb(a) = nb(a) + sgn
               if (any(nb < 1_I4P) .or. any(nb > g%n)) cycle
               if (g%dsq(nb(1), nb(2), nb(3)) >= AWRAP_EDT_FAR) cycle
               call unlin(g%n, g%site(nb(1), nb(2), nb(3)), si, sj, sk)
               r = cell_centre(g%o, g%h, si, sj, sk)
               b = min(b, sqrt((x%x - r%x)**2 + (x%y - r%y)**2 + (x%z - r%z)**2) - (g%u(si, sj, sk) - g%eps))
            enddo
         enddo
      enddo
      endfunction ball_gap
   endfunction edge_point

   subroutine closest_point(input_facet, input_tree, x, cp, u)
   !< Closest input point cp to x, and the distance u.
   type(facet_object),     intent(in)  :: input_facet(:) !< Input facets.
   type(aabb_tree_object), intent(in)  :: input_tree     !< AABB tree over input_facet.
   type(vector_R8P),       intent(in)  :: x              !< Point.
   type(vector_R8P),       intent(out) :: cp             !< Closest point.
   real(R8P),              intent(out) :: u              !< Distance.
   real(R8P)                           :: d2
   integer(I4P)                        :: cf, cr

   call input_tree%distance_tree_with_region(facet=input_facet, point=x, distance=d2, closest_facet=cf, closest_region=cr)
   call input_facet(cf)%compute_distance_with_region(point=x, distance=d2, closest=cp, region=cr)
   u = sqrt(max(d2, 0._R8P))
   endsubroutine closest_point

   pure function boxes_overlap(ta, tb) result(yes)
   !< Do the bounding boxes of two triangles overlap (closed)?
   type(vector_R8P), intent(in) :: ta(3), tb(3)
   logical                      :: yes

   yes = max(ta(1)%x, ta(2)%x, ta(3)%x) >= min(tb(1)%x, tb(2)%x, tb(3)%x) .and. &
         max(tb(1)%x, tb(2)%x, tb(3)%x) >= min(ta(1)%x, ta(2)%x, ta(3)%x) .and. &
         max(ta(1)%y, ta(2)%y, ta(3)%y) >= min(tb(1)%y, tb(2)%y, tb(3)%y) .and. &
         max(tb(1)%y, tb(2)%y, tb(3)%y) >= min(ta(1)%y, ta(2)%y, ta(3)%y) .and. &
         max(ta(1)%z, ta(2)%z, ta(3)%z) >= min(tb(1)%z, tb(2)%z, tb(3)%z) .and. &
         max(tb(1)%z, tb(2)%z, tb(3)%z) >= min(ta(1)%z, ta(2)%z, ta(3)%z)
   endfunction boxes_overlap

   pure function plane_relation(ta, tb) result(rel)
   !< Floating-point filter before the exact test of two wrap triangles, on the vertices of tb not shared with ta (bit
   !< for bit) against the plane of ta, with a relative tolerance of 1e-9: PLANE_APART, all on one side (the two meet
   !< at most at the shared vertex); PLANE_COPLANAR, all within the tolerance; PLANE_EXACT otherwise. Saves the exact
   !< orient3d expansions on the fan around a vertex and on the coplanar neighbours of a flat face.
   type(vector_R8P), intent(in) :: ta(3), tb(3)
   integer(I4P)                 :: rel
   real(R8P)                    :: n(3), d(3), tol
   logical                      :: own(3)
   integer(I4P)                 :: i

   n = normal3(ta)
   tol = 0._R8P
   do i = 1_I4P, 3_I4P
      own(i) = .not. (same_xyz(tb(i), ta(1)) .or. same_xyz(tb(i), ta(2)) .or. same_xyz(tb(i), ta(3)))
      d(i) = n(1) * (tb(i)%x - ta(1)%x) + n(2) * (tb(i)%y - ta(1)%y) + n(3) * (tb(i)%z - ta(1)%z)
      tol = max(tol, abs(tb(i)%x - ta(1)%x) + abs(tb(i)%y - ta(1)%y) + abs(tb(i)%z - ta(1)%z))
   enddo
   tol = 1.e-9_R8P * norm2(n) * tol
   rel = PLANE_EXACT
   if (all(abs(d) <= tol .or. .not. own)) then
      rel = PLANE_COPLANAR
   elseif (all(d > tol .or. .not. own) .or. all(d < -tol .or. .not. own)) then
      rel = PLANE_APART
   endif
   endfunction plane_relation

   pure function float_pair(ta, tb) result(res)
   !< Floating-point crossing test of two triangles that meet each other's planes (a common vertex, bit for bit, or a
   !< vertex within a relative AWRAP_FLOAT_ON of a plane lies on it): each meets the other's plane in a segment on the
   !< planes' common line, and they cross iff the two segments overlap with positive length. 1: crossing, 0: not,
   !< 2: they meet at a single point at most (a section is one point: a vertex on the other plane, the rest on one
   !< side), -1: not decided here, left to the exact pair_classify (a vertex within AWRAP_FLOAT_BAND of a plane but not
   !< on it, a section of three points, or, without a common vertex, an overlap within AWRAP_FLOAT_BAND of zero).
   type(vector_R8P), intent(in) :: ta(3), tb(3)
   integer(I4P)                 :: res
   real(R8P)                    :: xa(3, 2), xb(3, 2), l(3), sa(2), sb(2), ov, scale
   logical                      :: ok, shared
   integer(I4P)                 :: m1

   res = -1_I4P
   call section(ta, tb, xa, ok, shared, m1) ; if (.not. ok) then ; if (m1 == 1_I4P) res = 2_I4P ; return ; endif
   call section(tb, ta, xb, ok, shared, m1) ; if (.not. ok) then ; if (m1 == 1_I4P) res = 2_I4P ; return ; endif
   l = cross_arr(normal3(ta), normal3(tb))
   if (norm2(l) == 0._R8P) return
   sa = [dot_product(xa(:, 1), l), dot_product(xa(:, 2), l)]
   sb = [dot_product(xb(:, 1), l), dot_product(xb(:, 2), l)]
   ov = min(maxval(sa), maxval(sb)) - max(minval(sa), minval(sb))
   scale = AWRAP_FLOAT_BAND * max(maxval(sa) - minval(sa), maxval(sb) - minval(sb))
   if (ov > scale) then
      res = 1_I4P
   elseif (shared .or. ov < -scale) then
      res = 0_I4P ! both sections start at the common vertex (the same bits): no overlap beyond it is a touch there
   endif
   contains
      pure subroutine section(t, o, x, ok, shared, m)
      !< The two points where t meets o's plane: its vertices on the plane (a common vertex, bit for bit, or within the
      !< tolerance), and the crossings of its edges whose ends are strictly on either side.
      type(vector_R8P), intent(in)  :: t(3), o(3)
      real(R8P),        intent(out) :: x(3, 2)
      logical,          intent(out) :: ok, shared
      integer(I4P),     intent(out) :: m
      real(R8P)                     :: n(3), d(3), tol, v(3, 3)
      logical                       :: on(3)
      integer(I4P)                  :: i, j

      n = normal3(o)
      tol = 0._R8P
      do i = 1_I4P, 3_I4P
         v(:, i) = [t(i)%x, t(i)%y, t(i)%z]
         on(i) = same_xyz(t(i), o(1)) .or. same_xyz(t(i), o(2)) .or. same_xyz(t(i), o(3))
         d(i) = n(1) * (t(i)%x - o(1)%x) + n(2) * (t(i)%y - o(1)%y) + n(3) * (t(i)%z - o(1)%z)
         tol = max(tol, abs(t(i)%x - o(1)%x) + abs(t(i)%y - o(1)%y) + abs(t(i)%z - o(1)%z))
      enddo
      shared = any(on)
      ok = .false.
      m = -1_I4P
      if (any(.not. on .and. abs(d) > AWRAP_FLOAT_ON * norm2(n) * tol .and. abs(d) <= AWRAP_FLOAT_BAND * norm2(n) * tol)) return
      on = on .or. abs(d) <= AWRAP_FLOAT_ON * norm2(n) * tol
      m = 0_I4P
      do i = 1_I4P, 3_I4P
         if (on(i)) then
            m = m + 1_I4P ; if (m > 2_I4P) return
            x(:, m) = v(:, i)
         endif
      enddo
      do i = 1_I4P, 3_I4P
         j = mod(i, 3_I4P) + 1_I4P
         if (on(i) .or. on(j)) cycle
         if (d(i) * d(j) >= 0._R8P) cycle
         m = m + 1_I4P ; if (m > 2_I4P) return
         x(:, m) = v(:, i) + (d(i) / (d(i) - d(j))) * (v(:, j) - v(:, i))
      enddo
      ok = (m == 2_I4P)
      endsubroutine section
   endfunction float_pair

   pure function cross_arr(a, b) result(c)
   !< Cross product of two arrays.
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: c(3)

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross_arr

   pure function coplanar_overlap(ta, tb) result(yes)
   !< Do two (near-)coplanar triangles overlap with positive area (a fold of the wrap within a flat region)? 2-D
   !< separating axes in the plane of ta's dominant normal component: separated if, for an edge of either triangle,
   !< the other's vertices are all on its outer side or on it (relative tolerance 1e-9).
   type(vector_R8P), intent(in) :: ta(3), tb(3)
   logical                      :: yes
   real(R8P)                    :: n(3), a(2, 3), b(2, 3)
   integer(I4P)                 :: drop

   n = normal3(ta)
   drop = maxloc(abs(n), dim=1)
   a = project2(ta, drop) ; b = project2(tb, drop)
   yes = .not. (separates(a, b) .or. separates(b, a))
   contains
      pure function project2(t, drop) result(q)
      type(vector_R8P), intent(in) :: t(3)
      integer(I4P),     intent(in) :: drop
      real(R8P)                    :: q(2, 3)
      integer(I4P)                 :: i

      do i = 1_I4P, 3_I4P
         select case (drop)
         case (1_I4P) ; q(:, i) = [t(i)%y, t(i)%z]
         case (2_I4P) ; q(:, i) = [t(i)%z, t(i)%x]
         case default ; q(:, i) = [t(i)%x, t(i)%y]
         endselect
      enddo
      endfunction project2

      pure function separates(p, q) result(sep)
      !< Is there an edge of p with all of q on its outer side or on it?
      real(R8P), intent(in) :: p(2, 3), q(2, 3)
      logical               :: sep
      real(R8P)             :: e(2), s(3), inner, tol
      integer(I4P)          :: i, j, k

      sep = .false.
      do i = 1_I4P, 3_I4P
         j = mod(i, 3_I4P) + 1_I4P ; k = mod(j, 3_I4P) + 1_I4P
         e = p(:, j) - p(:, i)
         inner = e(1) * (p(2, k) - p(2, i)) - e(2) * (p(1, k) - p(1, i))
         if (inner == 0._R8P) cycle
         tol = 1.e-9_R8P * abs(inner)
         do k = 1_I4P, 3_I4P
            s(k) = sign(1._R8P, inner) * (e(1) * (q(2, k) - p(2, i)) - e(2) * (q(1, k) - p(1, i)))
         enddo
         if (all(s <= tol)) then
            sep = .true.
            return
         endif
      enddo
      endfunction separates
   endfunction coplanar_overlap

   pure function normal3(t) result(n)
   !< Unnormalised normal of a triangle.
   type(vector_R8P), intent(in) :: t(3)
   real(R8P)                    :: n(3)

   n = [(t(2)%y - t(1)%y) * (t(3)%z - t(1)%z) - (t(2)%z - t(1)%z) * (t(3)%y - t(1)%y), &
        (t(2)%z - t(1)%z) * (t(3)%x - t(1)%x) - (t(2)%x - t(1)%x) * (t(3)%z - t(1)%z), &
        (t(2)%x - t(1)%x) * (t(3)%y - t(1)%y) - (t(2)%y - t(1)%y) * (t(3)%x - t(1)%x)]
   endfunction normal3

   pure function degenerate(tri) result(yes)
   !< Has the triangle two coincident vertices (bit for bit)?
   type(vector_R8P), intent(in) :: tri(3)
   logical                      :: yes

   yes = same_xyz(tri(1), tri(2)) .or. same_xyz(tri(2), tri(3)) .or. same_xyz(tri(3), tri(1))
   endfunction degenerate

   pure function same_xyz(a, b) result(same)
   !< True if a and b are the same point, bit for bit.
   type(vector_R8P), intent(in) :: a, b
   logical                      :: same

   same = a%x == b%x .and. a%y == b%y .and. a%z == b%z
   endfunction same_xyz
endmodule fossil_alpha_wrap
