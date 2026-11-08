!run tut13 tut13_blocks_overset
!< Tutorial 13: an O-grid around an elliptic cylinder (topology, validity, quality, location, stencils), then its
!< overset assembly with a Cartesian background (wall cutter, holes, donors, classification, fringes, connectivity).
program tut13_blocks_overset
!region block
use fossil, only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_ABUT, TOPO_BOUNDARY, &
                   CELL_VALID, CELL_INVALID, LOC_INSIDE, LOC_OUTSIDE, STENCIL_CENTRED
use fossil, only : curve_object, wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, &
                   overset_mask_t, classify_2d, fringe_levels_2d, overlap_width_2d, connectivity_t, &
                   build_connectivity, save_connectivity, SELECT_EXPLICIT, FRINGE_INNER_1, FRINGE_INNER_2, &
                   FRINGE_OUTER_1, FRINGE_OUTER_2, ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17            ! O-grid vertices: around, out from the wall
real(R8P),    parameter      :: PI = acos(-1._R8P)
real(R8P),    parameter      :: A = 1._R8P, B = 0.4_R8P      ! semi-axes of the elliptic cylinder
type(curvilinear_block_2d_t) :: blocks(2), folded, thin(2)
type(topo_record_t)          :: rec
real(R8P)                    :: u(NI, NJ), v(NI, NJ), x(41, 33), y(41, 33), c, mu0, mu, s
real(R8P),    allocatable    :: area(:,:), skew(:,:), ortho(:,:), aspect(:,:), xi(:,:), xi_walk(:,:), w(:,:)
real(R8P),    allocatable    :: pu(:), pv(:), d(:,:)
integer(I4P), allocatable    :: valid(:,:), cell(:,:), cell_walk(:,:), st(:), st_walk(:), nodes(:,:,:), sst(:)
integer(I4P)                 :: i, j, k, r, q, n, fallbacks, un

c = sqrt(A**2 - B**2) ; mu0 = atanh(B / A)                   ! elliptic coordinates: the wall is mu = mu0
do j = 1, NJ
   s = real(j - 1, R8P) / (NJ - 1)
   mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)   ! clustered at the wall
   do i = 1, NI
      u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
      v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
   enddo
enddo
u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)                          ! the seam i = 1 = NI, written exactly
call blocks(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(blocks(1:1))                         ! the seam: one TOPO_ABUT record per side
do r = 1, blocks(1)%nrec
   print '(A,I0,A,L1,A,I0,A,2I3)', 'record ', r, ': abut ', blocks(1)%rec(r)%kind == TOPO_ABUT, ', face ', &
         blocks(1)%rec(r)%face, ', transform', blocks(1)%rec(r)%transform(1:2)
enddo
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)                               ! face 3 (-j): the wall
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)                               ! face 4 (+j): the boundary that receives
!endregion block
!region validity
call blocks(1)%geometry(area=area)
print '(A,F10.6,A,F10.6)', 'O-grid area ', sum(area), ', annulus between the ellipses ', &
      PI * c**2 * (cosh(mu) * sinh(mu) - cosh(mu0) * sinh(mu0))
call blocks(1)%validity(valid)
print '(A,I0,A,I0)', 'valid cells ', count(valid == CELL_VALID), ' of ', size(valid)
call blocks(1)%quality(skewness=skew, orthogonality=ortho, aspect=aspect)
print '(A,F6.3,A,F6.3,A,F6.2)', 'max skewness ', maxval(skew), ', min orthogonality ', minval(ortho), &
      ', max aspect ', maxval(aspect)
folded = blocks(1)
folded%u(13, 5) = folded%u(13, 7) ; folded%v(13, 5) = folded%v(13, 7)  ! a vertex dragged two layers out
call folded%validity(valid)
print '(A,I0,A)', 'folded copy: ', count(valid == CELL_INVALID), ' invalid cells'
do j = 1, NJ - 1
   do i = 1, NI - 1
      if (valid(i, j) == CELL_INVALID) print '(A,2I3)', '  invalid cell', i, j
   enddo
enddo
!endregion validity
!region locate
call blocks(1)%build_search                                  ! validity, orientation, the cell-box tree
n = 200
allocate(pu(n), pv(n))
do i = 1, n                                                  ! points on a spiral around the body
   s = real(i - 1, R8P) / (n - 1)
   pu(i) = (0.9_R8P + 0.9_R8P * s) * cos(14 * PI * s) ; pv(i) = (0.3_R8P + 1.1_R8P * s) * sin(14 * PI * s)
enddo
call blocks(1)%locate_many(pu, pv, cell, xi, st)             ! robust: the cell-box tree, exact containment
print '(A,I0,A,I0,A,I0)', 'located ', count(st == LOC_INSIDE), ', outside ', count(st == LOC_OUTSIDE), ' of ', n
call blocks(1)%locate_many(pu, pv, cell_walk, xi_walk, st_walk, walk=.true., fallbacks=fallbacks)
print '(A,L1,A,I0)', 'walk equal to the robust path: ', all(cell_walk == cell) .and. all(xi_walk == xi) .and. &
      all(st_walk == st), ', fallbacks ', fallbacks
!endregion locate
!region stencils
do q = 2, 4
   call blocks(1)%stencils(cell, xi, q, nodes, w, sst)       ! q x q Lagrange stencils, order-matched xi
   s = 0._R8P
   do i = 1, n
      if (st(i) /= LOC_INSIDE) cycle
      s = max(s, abs(sum([(w(r, i) * f(blocks(1)%u(nodes(1, r, i), nodes(2, r, i)), &
                                       blocks(1)%v(nodes(1, r, i), nodes(2, r, i))), r = 1, q**2)]) - f(pu(i), pv(i))))
   enddo
   print '(A,I0,A,I0,A,ES9.2)', 'q = ', q, ': centred stencils ', count(sst == STENCIL_CENTRED .and. st == LOC_INSIDE), &
         ', max error of sin(2x) cos(y) ', s
enddo
!endregion stencils
!region assembly
do j = 1, 33                                                 ! a Cartesian background, spacing 0.15
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(2)%build_search
block
type(curve_object)                 :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)               :: conn
integer(I4P)                       :: iterations, circular
character(200)                     :: line

call wall_cutter_2d(blocks(1:1), wall)                       ! the cutter: the O-grid's wall edges
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.]) ! cut the background only
print '(A,I0,A,I0)', 'background holes ', count(holes(2)%v), ' of ', size(holes(2)%v)
call donor_candidates_2d(blocks, donors)                     ! every vertex located in the other block
print '(A,I0,A,I0)', 'donor candidates: O-grid ', size(donors(1)%block), ', background ', size(donors(2)%block)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, iterations=iterations, &
                 circular=circular)
call fringe_levels_2d(blocks, masks)
do r = 1, 2
   print '(A,A,A,I0,A,I0,A,I0,A,I0)', '  ', blocks(r)%name, ': field ', count(masks(r)%iota == 1), ', holes ', &
         count(masks(r)%iota == 0), ', receptors ', count(masks(r)%iota < 0), ', orphans ', count(masks(r)%orphan > 0)
enddo
print '(A,I0,A,I0)', 'fixed point in ', iterations, ' iterations, circular pairs ', circular
print '(A,4(1X,I0))', 'fringe INNER 1, INNER 2, OUTER 1, OUTER 2:', &
      sum([(count(masks(r)%fringe == FRINGE_INNER_1), r = 1, 2)]), sum([(count(masks(r)%fringe == FRINGE_INNER_2), r = 1, 2)]), &
      sum([(count(masks(r)%fringe == FRINGE_OUTER_1), r = 1, 2)]), sum([(count(masks(r)%fringe == FRINGE_OUTER_2), r = 1, 2)])
call overlap_width_2d(blocks, masks, d)
print '(A,F7.4)', 'overlap width ', d(1, 2)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
call save_connectivity(conn, 'tut13-connectivity.txt')
open(newunit=un, file='tut13-connectivity.txt', action='read')
do r = 1, 5                                                  ! the header and the first receptor
   read(un, '(A)') line
   print '(A)', trim(line)
enddo
close(un)
!endregion assembly
!region orphans
do k = 10, 14, 2                                             ! the first k - 1 cell layers of the O-grid only
   call thin(1)%destroy
   call thin(1)%initialize(u(:, 1:k), v(:, 1:k), name='thin')
   call detect_topology_2d(thin(1:1))
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
   call thin(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, k, 1] ; rec%hi = [NI, k, 1] ; rec%name = 'overlap'
   call thin(1)%add_record(rec)
   thin(2) = blocks(2)
   call thin(1)%build_search
   call wall_cutter_2d(thin(1:1), wall)
   call hole_flags_2d(thin, wall, holes, cuts=[.false., .true.])
   call donor_candidates_2d(thin, donors)
   call classify_2d(thin, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', 'O-grid of ', k - 1, ' layers: receptors ', count(masks(1)%iota < 0), &
         ', orphans ', count(masks(1)%orphan > 0), ' (stencil on a hole ', &
         count(masks(1)%orphan == ORPHAN_STENCIL_HOLE), ', receptor in the stencil ', &
         count(masks(1)%orphan == ORPHAN_RECEPTOR_DONOR), ')'
enddo
endblock
!endregion orphans
contains
   pure function f(px, py) result(val)
   !< The interpolated function.
   real(R8P), intent(in) :: px, py !< Point.
   real(R8P)             :: val    !< Value.

   val = sin(2 * px) * cos(py)
   endfunction f
endprogram tut13_blocks_overset
