!run fig13 fig_tut13
!fig tut13-block
!fig tut13-overset
!< Figures of tutorial 13. First the O-grid around the elliptic cylinder: the spiral of points the chapter locates
!< (inside the block, outside it), and the 4 x 4 stencil of a point next to the seam, which takes vertices from both
!< sides of it. Then the assembly of the O-grid cut to 11 cell layers with the Cartesian background: the background
!< holes, the background receptors around them, the O-grid's outer receptors and its orphans.
program fig_tut13
use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY, LOC_INSIDE, &
                       curve_object, wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, &
                       overset_mask_t, classify_2d, SELECT_EXPLICIT
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17, NT = 12  ! O-grid vertices; layers of the thin O-grid + 1
real(R8P),    parameter      :: PI = acos(-1._R8P), A = 1._R8P, B = 0.4_R8P, LIFT = 1.e-3_R8P
type(curvilinear_block_2d_t) :: blocks(2)
type(topo_record_t)          :: rec
type(curve_object)           :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
real(R8P)                    :: u(NI, NJ), v(NI, NJ), x(41, 33), y(41, 33), c, mu0, mu, s
real(R8P),    allocatable    :: pu(:), pv(:), xi(:,:), w(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:), nodes(:,:,:), sst(:)
type(vector_R8P), allocatable :: pa(:), pb(:)
character(100)               :: caption
integer(I4P)                 :: i, j, n, f, m

c = sqrt(A**2 - B**2) ; mu0 = atanh(B / A)
do j = 1, NJ
   s = real(j - 1, R8P) / (NJ - 1)
   mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)
   do i = 1, NI
      u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
      v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
   enddo
enddo
u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)

! the block, the located spiral, one stencil across the seam
call blocks(1)%initialize(u, v)
call detect_topology_2d(blocks(1:1))
call blocks(1)%build_search
n = 200
allocate(pu(n + 1), pv(n + 1))
do i = 1, n
   s = real(i - 1, R8P) / (n - 1)
   pu(i) = (0.9_R8P + 0.9_R8P * s) * cos(14 * PI * s) ; pv(i) = (0.3_R8P + 1.1_R8P * s) * sin(14 * PI * s)
enddo
pu(n + 1) = 1.3_R8P ; pv(n + 1) = -0.03_R8P                   ! next to the seam, below it
call blocks(1)%locate_many(pu, pv, cell, xi, st)
call blocks(1)%stencils(cell(:, n + 1:n + 1), xi(:, n + 1:n + 1), 4_I4P, nodes, w, sst)
write(caption, '(A,I0,A,I0,A)') 'O-grid: ', count(st(1:n) == LOC_INSIDE), ' of ', n, &
                                 ' points located; a 4 x 4 stencil across the seam'
call fig_open(f, 'tut13-block', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
call grid_lines(blocks(1), 'body', 0._R8P)
call fig_points(f, 'inside', pack([(pt(pu(i), pv(i), 2 * LIFT), i = 1, n)], st(1:n) == LOC_INSIDE))
call fig_points(f, 'outside', pack([(pt(pu(i), pv(i), 2 * LIFT), i = 1, n)], st(1:n) /= LOC_INSIDE))
call fig_segments(f, 'patch', [pt(u(1, 1), v(1, 1), 3 * LIFT)], [pt(u(1, NJ), v(1, NJ), 3 * LIFT)])
allocate(pa(24), pb(24)) ; m = 0                              ! the stencil's lattice: nodes k = i + 4 (j - 1)
do j = 0, 3
   do i = 1, 3
      m = m + 1 ; pa(m) = node(4 * j + i) ; pb(m) = node(4 * j + i + 1)
      m = m + 1 ; pa(m) = node(4 * (i - 1) + j + 1) ; pb(m) = node(4 * i + j + 1)
   enddo
enddo
call fig_segments(f, 'accent', pa, pb)
call fig_points(f, 'foot', [(node(i), i = 1, 16)])
call fig_points(f, 'hole', [pt(pu(n + 1), pv(n + 1), 5 * LIFT)])
call fig_label(f, 'patch', pt(1.62_R8P, 0.12_R8P, 5 * LIFT), 'seam')
deallocate(pa, pb)
call fig_close(f)

! the assembly of the thin O-grid with the background
call blocks(1)%destroy
call blocks(1)%initialize(u(:, 1:NT), v(:, 1:NT))
call detect_topology_2d(blocks(1:1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, NT, 1] ; rec%hi = [NI, NT, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, 33
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
call donor_candidates_2d(blocks, donors)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'O-grid of ', NT - 1, ' layers: ', count(masks(1)%iota < 0), ' receptors, ', &
                                     count(masks(1)%orphan > 0), ' orphans, ', count(masks(2)%iota == 0), ' holes'
call fig_open(f, 'tut13-overset', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(0._R8P, 0._R8P), 1.3_R8P)
allocate(pa(2 * 40 * 33), pb(2 * 40 * 33)) ; m = 0
do j = 1, 33
   do i = 1, 40
      m = m + 1 ; pa(m) = pt(x(i, j), y(i, j)) ; pb(m) = pt(x(i + 1, j), y(i + 1, j))
   enddo
enddo
do j = 1, 32
   do i = 1, 41
      m = m + 1 ; pa(m) = pt(x(i, j), y(i, j)) ; pb(m) = pt(x(i, j + 1), y(i, j + 1))
   enddo
enddo
call fig_segments(f, 'box', pa(1:m), pb(1:m))
call grid_lines(blocks(1), 'body', LIFT)
call fig_points(f, 'ghost', pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, 41), j = 1, 33)], masks(2)%iota == 0))
call fig_points(f, 'probe', pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, 41), j = 1, 33)], masks(2)%iota < 0))
call fig_points(f, 'foot', pack([((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, NT)], masks(1)%iota < 0))
call fig_points(f, 'outside', pack([((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, NT)], masks(1)%orphan > 0))
call fig_close(f)
contains
   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt

   function node(k) result(q)
   !< Node k of the stencil of the last point, lifted over the grid.
   integer(I4P), intent(in) :: k !< Node.
   type(vector_R8P)         :: q !< Point.

   q = pt(u(nodes(1, k, 1), nodes(2, k, 1)), v(nodes(1, k, 1), nodes(2, k, 1)), 3 * LIFT)
   endfunction node

   subroutine grid_lines(blk, style, dz)
   !< The grid lines of a 2-D block.
   type(curvilinear_block_2d_t), intent(in) :: blk   !< Block.
   character(*),                 intent(in) :: style !< Style.
   real(R8P),                    intent(in) :: dz    !< Lift.
   type(vector_R8P), allocatable            :: a(:), b(:) !< Segments.
   integer(I4P)                             :: i, j, k    !< Counters.

   allocate(a(2 * blk%ni * blk%nj), b(2 * blk%ni * blk%nj)) ; k = 0
   do j = 1, blk%nj
      do i = 1, blk%ni
         if (i < blk%ni) then
            k = k + 1 ; a(k) = pt(blk%u(i, j), blk%v(i, j), dz) ; b(k) = pt(blk%u(i + 1, j), blk%v(i + 1, j), dz)
         endif
         if (j < blk%nj) then
            k = k + 1 ; a(k) = pt(blk%u(i, j), blk%v(i, j), dz) ; b(k) = pt(blk%u(i, j + 1), blk%v(i, j + 1), dz)
         endif
      enddo
   enddo
   call fig_segments(f, style, a(1:k), b(1:k))
   endsubroutine grid_lines
endprogram fig_tut13
