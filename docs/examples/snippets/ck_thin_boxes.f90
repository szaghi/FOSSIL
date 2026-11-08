program ck_thin_boxes
use fossil, only : surface_stl_object, THICKNESS_MEDIAL, STATUS_OK
use penf,   only : I4P, R8P
implicit none
integer(I4P), parameter  :: NX = 40, NY = 8           ! one layer of cells of side h over the wing section
real(R8P),    parameter  :: H = 1._R8P / 32._R8P      ! chord / 32
type(surface_stl_object) :: wing
real(R8P), allocatable   :: thickness(:), gap(:)
real(R8P)                :: lo(3, NX * NY), hi(3, NX * NY)
logical                  :: touched(NX * NY)
integer(I4P)             :: i, j, c, status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL, status=status)
call wing%build_thin(thickness=thickness, gap=gap, status=status)   ! local size per facet, in the tree
if (status /= STATUS_OK) error stop 'build_thin failed'
c = 0
do j = 1, NY ; do i = 1, NX
   c = c + 1 ; lo(:, c) = [-0.3_R8P + H * (i - 1), -0.125_R8P + H * (j - 1), 0.5_R8P] ; hi(:, c) = lo(:, c) + H
enddo ; enddo
call wing%box_thin_many(lo=lo, hi=hi, factor=1._R8P, touched=touched, status=status)  ! thinner than 1 h
print '(A,I0,A,I0)', 'cells ', NX * NY, ', meeting a wall thinner than h ', count(touched)
print '(A,F7.4,A,F7.4)', 'their x range ', minval(lo(1, :), touched), ' to', maxval(hi(1, :), touched)
endprogram ck_thin_boxes
