!< Tutorial 13, step 2: points located in the O-grid of step 1 (build_search, locate_many, the walk) and interpolated
!< from it by Lagrange stencils of width 2, 3, 4. Figures: the spiral of points, inside the block or not, with the cells
!< that hold them; a 4 x 4 stencil next to the seam, which takes vertices from both sides of it.
program tut13_b_locate
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, LOC_INSIDE, LOC_OUTSIDE, STENCIL_CENTRED
use penf,       only : I4P, R8P
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_quad, fig_label, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17, N = 200
real(R8P),    parameter      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: ogrid(1)
real(R8P)                    :: u(NI, NJ), v(NI, NJ), pu(N), pv(N), s
real(R8P),    allocatable    :: xi(:,:), xi_walk(:,:)
integer(I4P), allocatable    :: cell(:,:), cell_walk(:,:), st(:), st_walk(:)
integer(I4P)                 :: i, fallbacks, jumps
real(R8P),    allocatable    :: w(:,:)
integer(I4P), allocatable    :: nodes(:,:,:), sst(:)
real(R8P)                    :: err
integer(I4P)                 :: q, k
real(R8P),    parameter      :: LIFT = 1.e-3_R8P
type(vector_R8P), allocatable :: pa(:), pb(:)
real(R8P)                    :: su(1), sv(1)
character(100)               :: caption
integer(I4P)                 :: f, j, m

call make_ogrid(u, v)                                        ! the O-grid of step 1
call ogrid(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(ogrid)                               ! the seam records
call ogrid(1)%build_search                                   ! validity, handedness, the cell-box tree
do i = 1, N                                                  ! points on a spiral around the body
   s = real(i - 1, R8P) / (N - 1)
   pu(i) = (0.9_R8P + 0.9_R8P * s) * cos(14 * PI * s) ; pv(i) = (0.3_R8P + 1.1_R8P * s) * sin(14 * PI * s)
enddo
call ogrid(1)%locate_many(pu, pv, cell, xi, st)              ! robust: the tree, exact containment
print '(A,I0,A,I0,A,I0)', 'located ', count(st == LOC_INSIDE), ', outside ', count(st == LOC_OUTSIDE), ' of ', N
print '(A,I0,A,2I3,A,2F7.4)', 'point ', 100, ': cell', cell(:, 100), ', xi', xi(:, 100)
call ogrid(1)%locate_many(pu, pv, cell_walk, xi_walk, st_walk, walk=.true., jumps=jumps, fallbacks=fallbacks)
print '(A,L1,A,I0,A,I0)', 'walk equal to the robust path: ', all(cell_walk == cell) .and. all(xi_walk == xi) .and. &
      all(st_walk == st), ', cells crossed ', jumps, ', fallbacks ', fallbacks
do q = 2, 4
   call ogrid(1)%stencils(cell, xi, q, nodes, w, sst)        ! q x q Lagrange stencils on the vertices
   err = 0._R8P
   do i = 1, N
      if (st(i) /= LOC_INSIDE) cycle
      err = max(err, abs(sum([(w(k, i) * fun(u(nodes(1, k, i), nodes(2, k, i)), v(nodes(1, k, i), nodes(2, k, i))), &
                                k = 1, q**2)]) - fun(pu(i), pv(i))))
   enddo
   print '(A,I0,A,I0,A,ES9.2)', 'q = ', q, ': centred stencils ', count(sst == STENCIL_CENTRED .and. st == LOC_INSIDE), &
         ', max error of sin(2x) cos(y) ', err
enddo

! the located spiral and the cells holding its points
write(caption, '(A,I0,A,I0,A)') 'locate_many: ', count(st == LOC_INSIDE), ' of ', N, &
                                ' points inside (teal), their cells shaded; outside ember'
call fig_open(f, 'tut13-locate', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
do i = 1, N
   if (st(i) /= LOC_INSIDE) cycle
   call fig_quad(f, 'frac', pt(u(cell(1, i), cell(2, i)), v(cell(1, i), cell(2, i))), &
                 pt(u(cell(1, i) + 1, cell(2, i)), v(cell(1, i) + 1, cell(2, i))), &
                 pt(u(cell(1, i) + 1, cell(2, i) + 1), v(cell(1, i) + 1, cell(2, i) + 1)), &
                 pt(u(cell(1, i), cell(2, i) + 1), v(cell(1, i), cell(2, i) + 1)), value=0.45_R8P)
enddo
call grid_lines('body', LIFT)
call fig_segments(f, 'ghost', [(pt(pu(i), pv(i), 2 * LIFT), i = 1, N - 1)], [(pt(pu(i), pv(i), 2 * LIFT), i = 2, N)])
call fig_points(f, 'inside', pack([(pt(pu(i), pv(i), 3 * LIFT), i = 1, N)], st == LOC_INSIDE))
call fig_points(f, 'outside', pack([(pt(pu(i), pv(i), 3 * LIFT), i = 1, N)], st /= LOC_INSIDE))
call fig_close(f)

! a 4 x 4 stencil next to the seam
su = 1.3_R8P ; sv = -0.03_R8P                                ! just below the seam
call ogrid(1)%locate_many(su, sv, cell, xi, st)
call ogrid(1)%stencils(cell, xi, 4_I4P, nodes, w, sst)
write(caption, '(A,2I3,A)') '4 x 4 stencil of a point in cell', cell(:, 1), ': columns 47, 48, seam, 2'
call fig_open(f, 'tut13-stencil', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(1.35_R8P, 0._R8P), 0.5_R8P)
call grid_lines('body', 0._R8P)
call fig_segments(f, 'patch', [pt(u(1, 1), v(1, 1), LIFT)], [pt(u(1, NJ), v(1, NJ), LIFT)])
allocate(pa(24), pb(24)) ; m = 0                              ! the stencil's lattice: nodes k = i + 4 (j - 1)
do j = 0, 3
   do i = 1, 3
      m = m + 1 ; pa(m) = node(4 * j + i) ; pb(m) = node(4 * j + i + 1)
      m = m + 1 ; pa(m) = node(4 * (i - 1) + j + 1) ; pb(m) = node(4 * i + j + 1)
   enddo
enddo
call fig_segments(f, 'accent', pa, pb)
call fig_points(f, 'foot', [(node(i), i = 1, 16)])
call fig_points(f, 'hole', [pt(su(1), sv(1), 5 * LIFT)])
call fig_label(f, 'patch', pt(1.6_R8P, 0.04_R8P, 5 * LIFT), 'seam')
call fig_close(f)
contains
   pure function fun(px, py) result(val)
   !< The interpolated function.
   real(R8P), intent(in) :: px, py !< Point.
   real(R8P)             :: val    !< Value.

   val = sin(2 * px) * cos(py)
   endfunction fun

   subroutine make_ogrid(u, v)
   !< The O-grid of step 1: elliptic coordinates about the ellipse of semi-axes 1, 0.4, layers clustered at the wall.
   real(R8P), intent(out) :: u(NI, NJ), v(NI, NJ) !< Vertices.
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

   function node(k) result(q)
   !< Node k of the stencil, lifted over the grid.
   integer(I4P), intent(in) :: k !< Node.
   type(vector_R8P)         :: q !< Point.

   q = pt(u(nodes(1, k, 1), nodes(2, k, 1)), v(nodes(1, k, 1), nodes(2, k, 1)), 3 * LIFT)
   endfunction node

   subroutine grid_lines(style, dz)
   !< The grid lines of the O-grid.
   character(*), intent(in)      :: style      !< Style.
   real(R8P),    intent(in)      :: dz         !< Lift.
   type(vector_R8P), allocatable :: a(:), b(:) !< Segments.
   integer(I4P)                  :: i, j, k    !< Counters.

   allocate(a(2 * NI * NJ), b(2 * NI * NJ)) ; k = 0
   do j = 1, NJ
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
   endsubroutine grid_lines
endprogram tut13_b_locate
