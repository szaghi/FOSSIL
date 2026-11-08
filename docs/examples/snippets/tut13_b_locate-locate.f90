use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, LOC_INSIDE, LOC_OUTSIDE, STENCIL_CENTRED
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17, N = 200
real(R8P),    parameter      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: ogrid(1)
real(R8P)                    :: u(NI, NJ), v(NI, NJ), pu(N), pv(N), s
real(R8P),    allocatable    :: xi(:,:), xi_walk(:,:)
integer(I4P), allocatable    :: cell(:,:), cell_walk(:,:), st(:), st_walk(:)
integer(I4P)                 :: i, fallbacks, jumps

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
