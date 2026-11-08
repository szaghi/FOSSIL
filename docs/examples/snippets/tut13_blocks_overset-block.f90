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
