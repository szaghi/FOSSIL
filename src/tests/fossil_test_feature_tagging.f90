!< FOSSIL test: feature tagging for AMR (issue #61, roadmap #58 WP3, step 2).
!<
!< 1. Before build_features: box_feature_many and tag_block(feature_tags=) give STATUS_INVALID_INPUT.
!< 2. cube-01 and a block of 14^3 cells over [-0.2, 1.2]^3: the cells box_feature_many marks touched are exactly those of
!<    a brute force over every feature segment with the same predicate (segment_box_overlap, box_tolerance).
!< 3. tag_block(feature_band=2, feature_tags=): equal to cell_bands of the brute-force touched set; the surface tags are
!<    those of tag_block without the feature arguments.
!< 4. box_feature_many_device after enter_device equals the CPU (host build: the device code path, directives inert).
!< 5. After analyze (a geometry change) the tree is dropped: STATUS_INVALID_INPUT again.
!< 6. Sanity: the cell at the centre of a face is not touched, a cell on a cube edge is.
program fossil_test_feature_tagging

use fossil,             only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use fossil_box_query,   only : box_tolerance, cell_bands
use fossil_feature_bvh, only : segment_box_overlap
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,               only : I4P, R8P
use vecfor,             only : vector_R8P

implicit none

integer(I4P), parameter           :: N = 14
type(surface_stl_object)          :: s
type(vertex_pool_object), pointer :: pool
type(vector_R8P)                  :: pa, pb, bmin, bmax
real(R8P)                         :: x(N + 1), lo(3, N**3), hi(3, N**3), tol, seg(6)
logical                           :: touched(N**3), brute(N**3), dev(N**3), ok, all_ok
integer(I4P)                      :: tags(N, N, N), tags0(N, N, N), ftags(N, N, N), ref(N, N, N)
integer(I4P), allocatable         :: edges(:,:)
integer(I4P)                      :: status, i, j, k, c, e

all_ok = .true.
x = [(-0.2_R8P + 1.4_R8P * (i - 1) / N, i = 1, N + 1)]
c = 0
do k = 1, N ; do j = 1, N ; do i = 1, N
   c = c + 1 ; lo(:, c) = [x(i), x(j), x(k)] ; hi(:, c) = [x(i + 1), x(j + 1), x(k + 1)]
enddo ; enddo ; enddo
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)

! ---- 1. not built
call s%box_feature_many(lo, hi, touched, status=status)
ok = status == STATUS_INVALID_INPUT
call s%tag_block(x, x, x, 2_I4P, tags, status=status, feature_band=2_I4P, feature_tags=ftags)
ok = ok .and. status == STATUS_INVALID_INPUT
print '(A,L1)', '1. refused before build_features: ', ok
all_ok = all_ok .and. ok

! ---- 2. against brute force
call s%build_features(status=status)
call s%box_feature_many(lo, hi, touched, status=k)
call s%feature_edges(edges=edges)
pool => s%get_vertex_pool()
bmin = s%get_bmin() ; bmax = s%get_bmax()
tol = box_tolerance([bmin%x, bmin%y, bmin%z], [bmax%x, bmax%y, bmax%z])
brute = .false.
do c = 1, N**3
   do e = 1, size(edges, 2)
      pa = pool%coord(edges(1, e)) ; pb = pool%coord(edges(2, e))
      seg = [pa%x, pa%y, pa%z, pb%x, pb%y, pb%z]
      if (segment_box_overlap(lo(:, c), hi(:, c), seg(1:3), seg(4:6), tol)) then
         brute(c) = .true. ; exit
      endif
   enddo
enddo
ok = status == STATUS_OK .and. k == STATUS_OK .and. all(touched .eqv. brute) .and. count(touched) > 0
print '(A,I0,A,I0,A,L1)', '2. touched cells ', count(touched), ' of ', N**3, ', equal to brute force: ', ok
all_ok = all_ok .and. ok

! ---- 3. feature tags
call s%tag_block(x, x, x, 2_I4P, tags0, status=status)
call s%tag_block(x, x, x, 2_I4P, tags, status=k, feature_band=2_I4P, feature_tags=ftags)
call cell_bands(reshape(brute, [N, N, N]), 2_I4P, ref)
ok = status == STATUS_OK .and. k == STATUS_OK .and. all(ftags == ref) .and. all(tags == tags0)
print '(A,I0,A,I0,A,L1)', '3. feature tags: band cells ', count(ftags <= 2), ', zero ', count(ftags == 0), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. device path
call s%enter_device
call s%box_feature_many_device(lo, hi, dev, status=status)
call s%exit_device
ok = status == STATUS_OK .and. all(dev .eqv. touched)
print '(A,L1)', '4. device path equals the CPU: ', ok
all_ok = all_ok .and. ok

! ---- 6. sanity
c = cell_of([0.5_R8P, 0.5_R8P, 1.0_R8P + 0.01_R8P])            ! a face centre
ok = .not. touched(c)
c = cell_of([1.0_R8P + 0.01_R8P, 0.5_R8P, 1.0_R8P + 0.01_R8P])  ! on the edge x = 1, z = 1
ok = ok .and. touched(c)
print '(A,L1)', '6. face centre not touched, cube edge touched: ', ok
all_ok = all_ok .and. ok

! ---- 5. dropped by analyze
call s%analyze
call s%box_feature_many(lo, hi, touched, status=status)
ok = status == STATUS_INVALID_INPUT
print '(A,L1)', '5. dropped by analyze: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   pure function cell_of(p) result(cc)
   !< Linear index of the cell containing p.
   real(R8P), intent(in) :: p(3)
   integer(I4P)          :: cc, ii(3)

   ii = int((p + 0.2_R8P) / (1.4_R8P / N)) + 1
   cc = ii(1) + N * (ii(2) - 1) + N * N * (ii(3) - 1)
   endfunction cell_of
endprogram fossil_test_feature_tagging
