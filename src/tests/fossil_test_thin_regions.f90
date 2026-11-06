!< FOSSIL test: under-resolution report and AMR criterion (issue #62, roadmap #58 WP4, step 3).
!<
!< A plate [0, 2] x [0, 1] x [0, 0.5] with a fin [0.9, 0.95] x [0, 1] x [0.3, 1] (boolean union), medial thickness.
!< 1. thin_regions at h = 0.1, factor 1: every fin side facet (|n_x| > 0.5, z > 0.55) lies in region 1, no facet below
!<    the plate top (centroid z < 0.45) is thin, and region 1 has at least the area of the two fin sides (2 x 0.5 x 1);
!<    the regions come sorted by decreasing area. As h decreases the thin area never grows; below min(local size) / factor
!<    nothing is reported.
!< 2. Before build_thin: box_thin_many and tag_block(thin_tags=) give STATUS_INVALID_INPUT.
!< 3. box_thin_many on a 16^3 block over the bounding box (factor 2, so the 0.07 cells see the fin) equals a brute force
!<    over every facet with the same predicate; tag_block(thin_band=1, thin_tags=) equals cell_bands of it.
!< 4. box_thin_many_device equals the CPU (host build: the device code path).
!< 5. After move (a rotation and a translation) the thin data is kept and still equals the brute force on the moved
!<    facets; after analyze it is dropped (STATUS_INVALID_INPUT).
program fossil_test_thin_regions

use fossil,                only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle, STATUS_OK, &
                                  STATUS_INVALID_INPUT, THICKNESS_MEDIAL
use fossil_boolean,        only : BOOL_UNION, BOOL_STATUS_OK
use fossil_box_predicates, only : box_triangle_overlap
use fossil_box_query,      only : box_tolerance, cell_bands
use penf,                  only : I4P, R8P
use vecfor,                only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

integer(I4P), parameter     :: N = 16
type(surface_stl_object)    :: s, fin
type(facet_object), pointer :: fp
type(rigid_pose_t)          :: pose
type(vector_R8P)            :: bmin, bmax
real(R8P), allocatable      :: th(:), area(:), lo(:,:), hi(:,:)
integer(I4P), allocatable   :: region(:)
logical, allocatable        :: touched(:), brute(:), dev(:)
integer(I4P)                :: tags(N, N, N), ttags(N, N, N), ref(N, N, N)
real(R8P)                   :: h, prev, x(N + 1), y(N + 1), z(N + 1)
integer(I4P)                :: status, i, k
logical                     :: ok, all_ok

all_ok = .true.
call box(s, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P])
call box(fin, [0.9_R8P, 0._R8P, 0.3_R8P], [0.95_R8P, 1._R8P, 1._R8P])
call s%boolean(fin, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK .or. .not. s%is_watertight()) error stop 'plate and fin: union failed'
call s%thickness(thickness=th, method=THICKNESS_MEDIAL, status=status)

! ---- 1. thin regions
call s%thin_regions(thickness=th, h=0.1_R8P, factor=1._R8P, region=region, region_area=area, status=status)
ok = status == STATUS_OK .and. size(area) >= 1
if (ok) ok = area(1) >= 1._R8P - 1e-12_R8P .and. all(area(:size(area) - 1) >= area(2:))
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   if (abs(fp%normal%x) > 0.5_R8P .and. fp%centroid%z > 0.55_R8P) ok = ok .and. region(i) == 1
   if (fp%centroid%z < 0.45_R8P) ok = ok .and. region(i) == 0
enddo
print '(A,I0,A,F7.4,A,L1)', '1. h = 0.1: regions ', size(area), ', region 1 area ', area(1), '  ok ', ok
prev = huge(1._R8P)
h = 0.2_R8P
do while (h > 0.5_R8P * minval(th))
   call s%thin_regions(thickness=th, h=h, factor=1._R8P, region=region, region_area=area)
   ok = ok .and. sum(area) <= prev
   prev = sum(area) ; h = 0.8_R8P * h
enddo
ok = ok .and. size(area) == 0
print '(A,ES10.3,A,L1)', '   h sweep: thin area never grows, none below min(local size) = ', minval(th), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. not built
bmin = s%get_bmin() ; bmax = s%get_bmax()
x = [(bmin%x - 0.05_R8P + (bmax%x - bmin%x + 0.1_R8P) * (i - 1) / N, i = 1, N + 1)]
y = [(bmin%y - 0.05_R8P + (bmax%y - bmin%y + 0.1_R8P) * (i - 1) / N, i = 1, N + 1)]
z = [(bmin%z - 0.05_R8P + (bmax%z - bmin%z + 0.1_R8P) * (i - 1) / N, i = 1, N + 1)]
call cells(lo, hi)
allocate(touched(N**3), brute(N**3), dev(N**3))
call s%box_thin_many(lo, hi, 2._R8P, touched, status=status)
ok = status == STATUS_INVALID_INPUT
call s%tag_block(x, y, z, 1_I4P, tags, status=status, thin_factor=2._R8P, thin_band=1_I4P, thin_tags=ttags)
ok = ok .and. status == STATUS_INVALID_INPUT
print '(A,L1)', '2. refused before build_thin: ', ok
all_ok = all_ok .and. ok

! ---- 3. against brute force
call s%build_thin(thickness=th, status=status)
call s%box_thin_many(lo, hi, 2._R8P, touched, status=k)
call brute_force(brute)
call s%tag_block(x, y, z, 1_I4P, tags, status=i, thin_factor=2._R8P, thin_band=1_I4P, thin_tags=ttags)
call cell_bands(reshape(brute, [N, N, N]), 1_I4P, ref)
ok = status == STATUS_OK .and. k == STATUS_OK .and. i == STATUS_OK .and. all(touched .eqv. brute) .and. &
     count(touched) > 0 .and. all(ttags == ref)
print '(A,I0,A,I0,A,L1)', '3. thin cells ', count(touched), ' of ', N**3, ', equal to brute force, tags too: ', ok
all_ok = all_ok .and. ok

! ---- 4. device path
call s%enter_device
call s%box_thin_many_device(lo, hi, 2._R8P, dev, status=status)
call s%exit_device
ok = status == STATUS_OK .and. all(dev .eqv. touched)
print '(A,L1)', '4. device path equals the CPU: ', ok
all_ok = all_ok .and. ok

! ---- 5. move keeps, analyze drops
pose = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.3_R8P, [0.1_R8P, -0.2_R8P, 0.05_R8P])
call s%move(pose)
call s%box_thin_many(lo, hi, 2._R8P, touched, status=status)
call brute_force(brute)
ok = status == STATUS_OK .and. all(touched .eqv. brute)
call s%analyze
call s%box_thin_many(lo, hi, 2._R8P, touched, status=status)
ok = ok .and. status == STATUS_INVALID_INPUT
print '(A,L1)', '5. kept by move (equal to brute force), dropped by analyze: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine box(b, blo, bhi)
   !< The box [blo, bhi] from cube-01 ([0, 1]^3), resized and moved.
   type(surface_stl_object), intent(inout) :: b
   real(R8P),                intent(in)    :: blo(3), bhi(3)

   call b%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
   call b%resize(factor=(bhi(1) - blo(1)) * ex_R8P + (bhi(2) - blo(2)) * ey_R8P + (bhi(3) - blo(3)) * ez_R8P)
   call b%translate(delta=blo(1) * ex_R8P + blo(2) * ey_R8P + blo(3) * ez_R8P)
   call b%analyze
   endsubroutine box

   subroutine cells(clo, chi)
   !< The N^3 cell boxes of the block.
   real(R8P), allocatable, intent(out) :: clo(:,:), chi(:,:)
   integer(I4P)                        :: ii, jj, kk, c

   allocate(clo(3, N**3), chi(3, N**3))
   c = 0
   do kk = 1, N ; do jj = 1, N ; do ii = 1, N
      c = c + 1 ; clo(:, c) = [x(ii), y(jj), z(kk)] ; chi(:, c) = [x(ii + 1), y(jj + 1), z(kk + 1)]
   enddo ; enddo ; enddo
   endsubroutine cells

   subroutine brute_force(bt)
   !< Every cell against every facet thinner than 2 h_cell, the predicate of box_thin_many.
   logical, intent(out) :: bt(:)
   type(vector_R8P)     :: b0, b1
   real(R8P)            :: tol
   integer(I4P)         :: c, f

   b0 = s%get_bmin() ; b1 = s%get_bmax()
   tol = box_tolerance([b0%x, b0%y, b0%z], [b1%x, b1%y, b1%z])
   bt = .false.
   do c = 1, N**3
      do f = 1, s%get_facets_number()
         if (th(f) >= 2._R8P * maxval(hi(:, c) - lo(:, c))) cycle
         fp => s%facet_at(f)
         if (box_triangle_overlap(lo(:, c), hi(:, c), [fp%vertex(1)%x, fp%vertex(1)%y, fp%vertex(1)%z], &
                                  [fp%vertex(2)%x, fp%vertex(2)%y, fp%vertex(2)%z], &
                                  [fp%vertex(3)%x, fp%vertex(3)%y, fp%vertex(3)%z], tol)) then
            bt(c) = .true. ; exit
         endif
      enddo
   enddo
   endsubroutine brute_force
endprogram fossil_test_thin_regions
