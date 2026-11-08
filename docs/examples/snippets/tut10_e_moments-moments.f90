use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_coarsen_2d, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)        :: section
type(cut_cells_2d_t)      :: cut, coarse, direct
type(cut_links_t)         :: links, coarse_links
real(R8P)                 :: x(49), h, xc(2), qx
integer(I4P), allocatable :: fine_to_coarse(:)
integer(I4P)              :: i, k, p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]
call section%cut_cells(x=x, y=x, cut=cut, moments_order=2_I4P, boundaries=.true., quadrature=.true.)
p = cut%first_piece(23, 41)                          ! the first piece of the tail cell (23, 41), step 2
xc = 0.5_R8P * [x(23) + x(24), x(41) + x(42)]        ! the cell centre: moments in X = (x - xc) / h, Y = (y - yc) / h
print '(A,6F9.5)', 'moments / h^2 (1, X, Y, X^2, XY, Y^2):', cut%piece_moments(:, p) / h**2
print '(A,2F9.5,A,2F9.5)', '  centroid from the moments ', xc + h * cut%piece_moments(2:3, p) / cut%piece_moments(1, p), &
                           ', piece_centroid ', cut%piece_centroid(:, p)
