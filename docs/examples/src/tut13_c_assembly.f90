!run tut13c tut13_c_assembly
!fig tut13-holes
!fig tut13-donors
!fig tut13-fringes
!< Tutorial 13, step 3: the O-grid of step 1 assembled with a Cartesian background: the wall cutter, the background
!< holes (hole_flags_2d, checked against hole_map_2d), the donor candidates, the explicit classification, the fringe
!< levels, the overlap width and the connectivity file. Figures: the cutter and the holes it blanks; the fringe levels
!< of both grids.
program tut13_c_assembly
!region background
use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY
use fossil,     only : curve_object, wall_cutter_2d, hole_flags_2d, hole_map_2d, block_logical_t, donor_set_t, &
                       donor_candidates_2d, overset_mask_t, classify_2d, fringe_levels_2d, overlap_width_2d, &
                       connectivity_t, build_connectivity, save_connectivity, SELECT_EXPLICIT, FRINGE_INNER_1, &
                       FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
use penf,       only : I4P, R8P
!endregion background
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_label, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
!region background

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33
type(curvilinear_block_2d_t)       :: blocks(2)              ! 1: the O-grid, 2: the background
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(block_logical_t), allocatable :: holes(:), holes_map(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)               :: conn
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ)
real(R8P),             allocatable :: d(:,:)
integer(I4P)                       :: i, j, b, iterations, circular, un
character(200)                     :: line
!endregion background
real(R8P),    parameter            :: LIFT = 1.e-3_R8P
type(vector_R8P), allocatable      :: pa(:), pb(:), pts(:)
character(100)                     :: caption
integer(I4P)                       :: f, m
!region background

call make_ogrid(u, v)                                        ! the O-grid of step 1, with its records
call blocks(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(blocks(1:1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, MJ                                                 ! a Cartesian background, spacing 0.15
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
do b = 1, 2
   call blocks(b)%build_search
enddo
!endregion background
!region holes
call wall_cutter_2d(blocks(1:1), cutter)                     ! the cutter: the edges of the faces named 'wall'
print '(A,I0,A)', 'cutter: ', cutter%get_vertices_number(), ' vertices, one closed loop'
call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.]) ! the background only: the O-grid's wall is the cutter
print '(A,I0,A,I0)', 'background holes ', count(holes(2)%v), ' of ', size(holes(2)%v)
call hole_map_2d(blocks, cutter, holes_map, [64, 64], cuts=[.false., .true.])  ! a Cartesian hole map instead
print '(A,L1)', 'hole map gives the same holes: ', all(holes_map(2)%v .eqv. holes(2)%v)
!endregion holes
!region classify
call donor_candidates_2d(blocks, donors)                     ! every vertex located in the other block
print '(A,I0,A,I0)', 'donor candidates: O-grid ', size(donors(1)%block), ', background ', size(donors(2)%block)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, iterations=iterations, &
                 circular=circular)                          ! 3 x 3 donor stencils, explicit interpolation
do b = 1, 2
   print '(2X,A,A,I0,A,I0,A,I0,A,I0)', blocks(b)%name(1:10), ': field ', count(masks(b)%iota == 1), ', holes ', &
         count(masks(b)%iota == 0), ', receptors ', count(masks(b)%iota < 0), ', orphans ', count(masks(b)%orphan > 0)
enddo
print '(A,I0,A,I0)', 'fixed point in ', iterations, ' iterations, circular pairs ', circular
!endregion classify
!region fringes
call fringe_levels_2d(blocks, masks)                         ! INNER 1, 2 around the holes; OUTER 1, 2 at the overlap
print '(A,4(1X,I0))', 'fringe INNER 1, INNER 2, OUTER 1, OUTER 2:', &
      sum([(count(masks(b)%fringe == FRINGE_INNER_1), b = 1, 2)]), sum([(count(masks(b)%fringe == FRINGE_INNER_2), b = 1, 2)]), &
      sum([(count(masks(b)%fringe == FRINGE_OUTER_1), b = 1, 2)]), sum([(count(masks(b)%fringe == FRINGE_OUTER_2), b = 1, 2)])
call overlap_width_2d(blocks, masks, d)                      ! Chesshire's width between the two fringes
print '(A,F7.4)', 'overlap width ', d(1, 2)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)   ! the receptor list a solver reads
call save_connectivity(conn, 'tut13-connectivity.txt')
open(newunit=un, file='tut13-connectivity.txt', action='read')
do i = 1, 5                                                  ! the header and the first receptor
   read(un, '(A)') line
   print '(A)', trim(line)
enddo
close(un)
!endregion fringes

! the cutter and the holes it blanks
write(caption, '(A,I0,A)') 'the cutter (yellow) blanks ', count(holes(2)%v), ' background vertices (ember)'
call fig_open(f, 'tut13-holes', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(0._R8P, 0._R8P), 1.25_R8P)
call background_lines('box', 0._R8P)
call cutter_segments
call fig_segments(f, 'accent', pa, pb)
pts = pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)], holes(2)%v)
call fig_points(f, 'outside', pts)
call fig_close(f)

! two receptors and their donor stencils
call fig_open(f, 'tut13-donors', 'receptors: background teal, O-grid yellow; two 3 x 3 donor stencils', &
              -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(0._R8P, 0.98_R8P), 0.62_R8P)
call background_lines('box', 0._R8P)
call ogrid_lines('ghost', LIFT, 1_I4P)
call cutter_segments
call fig_segments(f, 'accent', pa, pb)
pts = pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)], masks(2)%iota < 0)
call fig_points(f, 'probe', pts)
pts = pack([((pt(u(i, j), v(i, j), 2 * LIFT), i = 1, NI), j = 1, NJ)], masks(1)%iota < 0)
call fig_points(f, 'foot', pts)
m = findloc(masks(2)%iota(21:MI * MJ:MI) < 0, .true., dim=1, back=.true.)
call donor_stencil(2_I4P, 21_I4P + (m - 1_I4P) * MI)         ! the top background receptor of column 21 (x = 0)
call donor_stencil(1_I4P, 13_I4P + (NJ - 1_I4P) * NI)         ! O-grid vertex (13, 17), on its overlap boundary
call fig_close(f)

! the fringe levels of both grids
call fig_open(f, 'tut13-fringes', 'fringe_levels_2d: INNER 1 ember, INNER 2 yellow; OUTER 1 teal, OUTER 2 sepia', &
              -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(0.45_R8P, 0.45_R8P), 1.3_R8P)
call background_lines('box', 0._R8P)
call ogrid_lines('ghost', LIFT, NJ - 2_I4P)
call cutter_segments
call fig_segments(f, 'accent', pa, pb)
call fig_points(f, 'ghost', mpts(2, 0))                     ! the holes
call fig_points(f, 'patch', mpts(2, FRINGE_INNER_1))
call fig_points(f, 'foot', mpts(2, FRINGE_INNER_2))
call fig_points(f, 'probe', mpts(1, FRINGE_OUTER_1))
call fig_points(f, 'body', mpts(1, FRINGE_OUTER_2))
call fig_close(f)
contains
   subroutine make_ogrid(u, v)
   !< The O-grid of step 1: elliptic coordinates about the ellipse of semi-axes 1, 0.4, layers clustered at the wall.
   real(R8P), intent(out) :: u(NI, NJ), v(NI, NJ) !< Vertices.
   real(R8P), parameter   :: PI = acos(-1._R8P)   !< Pi.
   real(R8P)              :: c, mu0, mu, s        !< Coordinates.
   integer(I4P)           :: i, j                 !< Counters.

   c = sqrt(1._R8P - 0.4_R8P**2) ; mu0 = atanh(0.4_R8P)
   do j = 1, NJ
      s = real(j - 1, R8P) / (NJ - 1)
      mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)
      do i = 1, NI
         u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
         v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
      enddo
   enddo
   u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)
   endsubroutine make_ogrid

   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt

   function mpts(blk, level) result(p)
   !< The vertices of block blk with fringe level `level` (0: the holes), lifted over the grids.
   integer(I4P), intent(in)      :: blk, level !< Block, level.
   type(vector_R8P), allocatable :: p(:)       !< Points.
   type(vector_R8P), allocatable :: all(:)     !< Every vertex.
   integer(I4P)                  :: i, j       !< Counters.

   if (blk == 1) then
      all = [((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, NJ)]
   else
      all = [((pt(x(i, j), y(i, j), 3 * LIFT), i = 1, MI), j = 1, MJ)]
   endif
   if (level == 0) then
      p = pack(all, masks(blk)%iota == 0)
   else
      p = pack(all, masks(blk)%fringe == level)
   endif
   endfunction mpts

   subroutine donor_stencil(rb, rv)
   !< The receptor rv of block rb (linear vertex index) and the 3 x 3 stencil of its chosen donor.
   integer(I4P), intent(in)      :: rb, rv            !< Receptor block, vertex.
   integer(I4P), allocatable     :: nodes(:,:,:), sst(:) !< Stencil.
   real(R8P),    allocatable     :: wt(:,:)           !< Weights.
   type(vector_R8P), allocatable :: sa(:), sb(:), sp(:) !< Stencil lattice and nodes.
   integer(I4P)                  :: c, db, k, l, n    !< Candidate, donor block, counters.

   if (masks(rb)%iota(rv) >= 0) error stop 'not a receptor'
   c = masks(rb)%donor(rv) ; db = donors(rb)%block(c)
   call blocks(db)%stencils(donors(rb)%cell(:, c:c), donors(rb)%xi(:, c:c), 3_I4P, nodes, wt, sst)
   allocate(sa(12), sb(12), sp(9)) ; n = 0
   do k = 1, 9
      sp(k) = vtx(db, nodes(1, k, 1), nodes(2, k, 1))
   enddo
   do l = 0, 2
      do k = 1, 2
         n = n + 1 ; sa(n) = sp(3 * l + k) ; sb(n) = sp(3 * l + k + 1)
         n = n + 1 ; sa(n) = sp(3 * (k - 1) + l + 1) ; sb(n) = sp(3 * k + l + 1)
      enddo
   enddo
   call fig_segments(f, 'patch', sa, sb)
   call fig_points(f, 'patch', sp)
   if (rb == 1) then
      sp(1) = vtx(1_I4P, mod(rv - 1, NI) + 1, (rv - 1) / NI + 1)
      call fig_label(f, 'foot', sp(1) + 0.08_R8P * ex_R8P + 0.03_R8P * ey_R8P, 'O-grid receptor')
   else
      sp(1) = vtx(2_I4P, mod(rv - 1, MI) + 1, (rv - 1) / MI + 1)
      call fig_label(f, 'probe', sp(1) + 0.08_R8P * ex_R8P + 0.05_R8P * ey_R8P, 'background receptor')
   endif
   call fig_points(f, 'inside', [sp(1) + LIFT * ez_R8P])
   endsubroutine donor_stencil

   function vtx(blk, i, j) result(q)
   !< Vertex (i, j) of block blk, lifted over the grids.
   integer(I4P), intent(in) :: blk, i, j !< Block, vertex.
   type(vector_R8P)         :: q         !< Point.

   if (blk == 1) then
      q = pt(u(i, j), v(i, j), 4 * LIFT)
   else
      q = pt(x(i, j), y(i, j), 4 * LIFT)
   endif
   endfunction vtx

   subroutine cutter_segments
   !< The segments of the cutter into pa, pb.
   real(R8P)    :: uv(2) !< Vertex.
   integer(I4P) :: s     !< Counter.

   if (allocated(pa)) deallocate(pa, pb)
   allocate(pa(cutter%get_vertices_number()), pb(cutter%get_vertices_number())) ; m = 0
   do s = 1, cutter%get_vertices_number()
      if (cutter%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = cutter%get_vertex(s) ; pa(m) = pt(uv(1), uv(2), 2 * LIFT)
      uv = cutter%get_vertex(cutter%get_next(s)) ; pb(m) = pt(uv(1), uv(2), 2 * LIFT)
   enddo
   pa = pa(1:m) ; pb = pb(1:m)
   endsubroutine cutter_segments

   subroutine background_lines(style, dz)
   !< The grid lines of the background.
   character(*), intent(in)      :: style      !< Style.
   real(R8P),    intent(in)      :: dz         !< Lift.
   type(vector_R8P), allocatable :: a(:), b(:) !< Segments.
   integer(I4P)                  :: i, j, k    !< Counters.

   allocate(a(2 * MI * MJ), b(2 * MI * MJ)) ; k = 0
   do j = 1, MJ
      do i = 1, MI
         if (i < MI) then
            k = k + 1 ; a(k) = pt(x(i, j), y(i, j), dz) ; b(k) = pt(x(i + 1, j), y(i + 1, j), dz)
         endif
         if (j < MJ) then
            k = k + 1 ; a(k) = pt(x(i, j), y(i, j), dz) ; b(k) = pt(x(i, j + 1), y(i, j + 1), dz)
         endif
      enddo
   enddo
   call fig_segments(f, style, a(1:k), b(1:k))
   endsubroutine background_lines

   subroutine ogrid_lines(style, dz, jmin)
   !< The grid lines of the O-grid from the ring jmin out.
   character(*), intent(in)      :: style      !< Style.
   real(R8P),    intent(in)      :: dz         !< Lift.
   integer(I4P), intent(in)      :: jmin       !< First ring.
   type(vector_R8P), allocatable :: a(:), b(:) !< Segments.
   integer(I4P)                  :: i, j, k    !< Counters.

   allocate(a(2 * NI * NJ), b(2 * NI * NJ)) ; k = 0
   do j = jmin, NJ
      do i = 1, NI
         if (i < NI) then
            k = k + 1 ; a(k) = pt(u(i, j), v(i, j), dz) ; b(k) = pt(u(i + 1, j), v(i + 1, j), dz)
         endif
         if (j < NJ) then
            k = k + 1 ; a(k) = pt(u(i, j), v(i, j), dz) ; b(k) = pt(u(i, j + 1), v(i, j + 1), dz)
         endif
      enddo
   enddo
   call fig_segments(f, style, a(1:k), b(1:k))
   endsubroutine ogrid_lines
endprogram tut13_c_assembly
