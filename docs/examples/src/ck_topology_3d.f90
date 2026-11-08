!run ck_topology_3d ck_topology_3d
program ck_topology_3d
use fossil, only : curvilinear_block_t, detect_topology, TOPO_SINGULAR
use penf,   only : I4P, R8P
implicit none
character(9), parameter   :: KIND(4) = ['boundary ', 'abut     ', 'periodic ', 'singular ']
real(R8P),    parameter   :: PI = acos(-1._R8P)
type(curvilinear_block_t) :: shell(1)
real(R8P)                 :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph
integer(I4P)              :: i, j, k

do k = 1, 3 ; do j = 1, 13 ; do i = 1, 25      ! a spherical shell 1 <= r <= 1.5: i azimuth, j polar angle, k radius
   r = 1._R8P + 0.25_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:)                ! the azimuth seam, written exactly
x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P             ! the poles, collapsed exactly
call shell(1)%initialize(x, y, z)
call detect_topology(shell)                                 ! exact coincidence
do i = 1, shell(1)%nrec
   associate(rec => shell(1)%rec(i))
      if (rec%kind == TOPO_SINGULAR) then
         print '(A,I0,3A,I0)', 'face ', rec%face, ': ', KIND(rec%kind), 'collapsed along direction ', rec%collapsed
      else
         print '(A,I0,3A,I0,A,3I3)', 'face ', rec%face, ': ', KIND(rec%kind), 'face ', rec%donor_face, &
               ', transform', rec%transform
      endif
   endassociate
enddo
endprogram ck_topology_3d
