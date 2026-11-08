!run figck-thin-boxes fig_ck_thin_boxes
!fig ck-thin-boxes
!< Figure of the cookbook recipe "Refine the cells across a thin wall": the layer of ck_thin_boxes (40 x 8 cells of
!< side chord / 32) over the NACA 0012 wing, seen along the span, with the wing's section in ember; the cells box_thin_many
!< flags dark, the others light.
program fig_ck_thin_boxes
use fossil,     only : surface_stl_object, curve_object, slice_surface, THICKNESS_MEDIAL
use fossil_fig, only : fig_open, fig_segments, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter  :: NX = 40, NY = 8
real(R8P),    parameter  :: H = 1._R8P / 32._R8P
type(surface_stl_object) :: wing
type(curve_object)       :: section
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P), allocatable   :: thickness(:), gap(:)
real(R8P)                :: lo(3, NX * NY), hi(3, NX * NY), z, uv(2)
logical                  :: touched(NX * NY)
integer(I4P)             :: i, j, c, u, n
character(90)            :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL)
call wing%build_thin(thickness=thickness, gap=gap)
c = 0
do j = 1, NY ; do i = 1, NX
   c = c + 1 ; lo(:, c) = [-0.3_R8P + H * (i - 1), -0.125_R8P + H * (j - 1), 0.5_R8P] ; hi(:, c) = lo(:, c) + H
enddo ; enddo
call wing%box_thin_many(lo=lo, hi=hi, factor=1._R8P, touched=touched)
write(caption, '(A,I0,A,I0,A)') 'box_thin_many(factor=1): ', count(touched), ' of ', NX * NY, &
                                ' cells meet a wall thinner than h = chord/32'
call fig_open(u, 'ck-thin-boxes', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
z = 1.01_R8P
do c = 1, NX * NY
   call fig_quad(u, 'frac', vector_R8P(lo(1, c), lo(2, c), z), vector_R8P(hi(1, c), lo(2, c), z), &
                 vector_R8P(hi(1, c), hi(2, c), z), vector_R8P(lo(1, c), hi(2, c), z), &
                 value=merge(0.9_R8P, 0.05_R8P, touched(c)))
enddo
call slice_surface(wing, 3_I4P, 0.5_R8P, section)
n = section%get_vertices_number() ; allocate(a(n), b(n))
do i = 1, n
   uv = section%get_vertex(i) ; a(i) = vector_R8P(uv(1), uv(2), 1.02_R8P)
   uv = section%get_vertex(section%get_next(i)) ; b(i) = vector_R8P(uv(1), uv(2), 1.02_R8P)
enddo
call fig_segments(u, 'patch', a, b)
call fig_close(u)
endprogram fig_ck_thin_boxes
