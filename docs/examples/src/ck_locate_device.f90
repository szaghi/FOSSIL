!run ck_locate_device ck_locate_device
program ck_locate_device
use fossil, only : curvilinear_block_t, LOC_INSIDE
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
integer(I4P), parameter   :: N = 4000
type(curvilinear_block_t) :: shell
integer(I4P), allocatable :: cell_dev(:,:), st_dev(:), cell_host(:,:), st_host(:)
real(R8P),    allocatable :: xi_dev(:,:), xi_host(:,:)
real(R8P)                 :: x(49, 25, 9), y(49, 25, 9), z(49, 25, 9), r, th, ph, px(N), py(N), pz(N)
integer(I4P)              :: i, j, k

do k = 1, 9 ; do j = 1, 25 ; do i = 1, 49                   ! a spherical shell 1 <= r <= 2, 48 x 24 x 8 cells
   r = 1._R8P + 0.125_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 48 ; ph = PI * (j - 1) / 24
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(49,:,:) = x(1,:,:) ; y(49,:,:) = y(1,:,:) ; x(:,[1,25],:) = 0._R8P ; y(:,[1,25],:) = 0._R8P
call shell%initialize(x, y, z)
call shell%build_search
do i = 1, N                                                 ! points on a spiral through the shell and beyond
   r = 0.8_R8P + 1.4_R8P * (i - 1) / (N - 1) ; th = 0.37_R8P * i ; ph = PI * (i - 0.5_R8P) / N
   px(i) = r * sin(ph) * cos(th) ; py(i) = r * sin(ph) * sin(th) ; pz(i) = r * cos(ph)
enddo
call shell%enter_device                                     ! copy the search arrays to the device, once
call shell%locate_many_device(px, py, pz, cell_dev, xi_dev, st_dev)
call shell%exit_device
call shell%locate_many(px, py, pz, cell_host, xi_host, st_host)
print '(A,I0,A,I0)', 'points ', N, ', inside a cell ', count(st_dev == LOC_INSIDE)
print '(A,I0,A,I0)', 'device vs host mismatches: cell or status ', &
                     count(any(cell_dev /= cell_host, 1) .or. st_dev /= st_host), ', xi ', count(xi_dev /= xi_host)
endprogram ck_locate_device
