use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_ABUT, TOPO_BOUNDARY, &
                       CELL_VALID, CELL_INVALID
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17            ! O-grid vertices: around the body, out from the wall
real(R8P),    parameter      :: PI = acos(-1._R8P)
real(R8P),    parameter      :: A = 1._R8P, B = 0.4_R8P      ! semi-axes of the elliptic cylinder
type(curvilinear_block_2d_t) :: ogrid(1), folded       ! a set of one block
type(topo_record_t)          :: rec
real(R8P)                    :: u(NI, NJ), v(NI, NJ), c, mu0, mu, s
integer(I4P)                 :: i, j, r

c = sqrt(A**2 - B**2) ; mu0 = atanh(B / A)                   ! elliptic coordinates: the wall is mu = mu0
do j = 1, NJ
   s = real(j - 1, R8P) / (NJ - 1)
   mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)   ! layers clustered at the wall
   do i = 1, NI
      u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
      v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
   enddo
enddo
u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)                          ! the seam i = 1 = NI, written exactly
call ogrid(1)%initialize(u, v, name='ellipse')
