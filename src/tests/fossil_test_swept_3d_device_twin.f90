!< FOSSIL test: the device path of the 3-D swept geometry equals the host path bit for bit (issue #100 step 4).
program fossil_test_swept_3d_device_twin
!< FOSSIL test: the device path of the 3-D swept geometry equals the host path bit for bit (issue #100 step 4).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, batching,
!< chunking and ordering bugs of `swept_cells_3d_device` (node lists, the posed plane sweep, the Gauss doubling in
!< rounds, the face sums across chunks), not device codegen (the distance bench does that on the GPU, `swept cells
!< 3-D device vs serial mismatches`). Every output of `swept_cells_3d_t` is compared bitwise
!< (`swept_cells_3d_mismatches`), with the default chunk and with chunk_cells = 37 (a few section blocks per launch
!< set, groups and planes split across chunks), on 4^3 blocks:
!< - exact translation: the unit cube (faces through no grid plane), and moving parallel to the x planes (d_x = 0:
!<   every node of an x plane at the same level), the cone tilted;
!< - exact rotation: the unit cube turning 0.3 rad about (1, 1, 1) through its centre while translating, the L-prism
!<   0.5 rad about a tilted off-centre axis;
!< - Gauss poses: the cube rotation with exact = .false. (nq = 2, the 2 nq error estimate), explicit poses (nq = 3);
!< - refusals: decreasing edges (both paths), an octree surface (the device path only).
use fossil,                   only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, &
                                     swept_cells_3d_t, swept_cells_3d, swept_cells_3d_device, gauss_legendre_01, &
                                     STATUS_OK, STATUS_INVALID_INPUT, AABB_TREE_OCTREE
use fossil_swept_3d_compare,  only : swept_cells_3d_mismatches
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable       :: b(:)
type(surface_stl_object)        :: s
type(rigid_pose_t)              :: p0, p1
type(rigid_pose_t), allocatable :: pq(:)
type(swept_cells_3d_t)          :: sh, sd
real(R8P),          allocatable :: sq(:), wq(:)
real(R8P)                       :: e(5), c(3)
integer(I4P)                    :: i, st, st2, total
logical                         :: ok

ok = .true. ; total = 0_I4P
call bench_bodies(b)
p0 = pose_identity()
call to_surface(b(1)%soup, s) ! unit cube [0, 1]^3
e = [(-0.23_R8P + 0.45_R8P * i, i = 0, 4)]
p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0.31_R8P, 0.17_R8P, 0.23_R8P])
call compare('cube translation', .true.)
p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0._R8P, 0.17_R8P, -0.23_R8P])
call compare('cube translation parallel to the x planes', .true.)
call to_surface(b(7)%soup, s) ! cone, at a rotated pose at both ends
p0 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.4_R8P, [0.6_R8P, 0.7_R8P, 0.5_R8P])
p1 = p0 ; p1%t = p1%t + [0.21_R8P, -0.13_R8P, 0.17_R8P]
e = [(-0.6_R8P + 0.675_R8P * i, i = 0, 4)]
call compare('tilted cone translation', .true.)
p0 = pose_identity()
call to_surface(b(1)%soup, s)
e = [(-0.25_R8P + 0.4_R8P * i, i = 0, 4)]
c = 0.5_R8P
p1 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P)
p1%t = c - matmul(p1%r, c) + [0.13_R8P, -0.07_R8P, 0.05_R8P]
call compare('cube rotation (exact)', .true.)
call compare_gauss('cube rotation (Gauss nq = 2)')
call gauss_legendre_01(3_I4P, sq, wq)
allocate(pq(3))
do i = 1, 3
   pq(i) = pose_interpolate(p0, p1, sq(i))
enddo
call swept_cells_3d(s, p0, p1, e, e, e, sh, nq=3_I4P, poses=pq, status=st)
call swept_cells_3d_device(s, p0, p1, e, e, e, sd, nq=3_I4P, poses=pq, status=st2)
call report('cube, explicit poses (nq = 3)', swept_cells_3d_mismatches(sh, sd) + abs(st - st2) + abs(st))
call to_surface(b(8)%soup, s) ! L-prism [0, 2]^2 x [0, 1]
e = [(-1.1_R8P + 0.9_R8P * i, i = 0, 4)]
c = [0.7_R8P, 0.9_R8P, 0.4_R8P]
p1 = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0.5_R8P)
p1%t = c - matmul(p1%r, c) + [0.05_R8P, 0.1_R8P, -0.08_R8P]
call compare('L-prism rotation (exact)', .false.)
! refusals
call swept_cells_3d(s, p0, p1, e(5:1:-1), e, e, sh, status=st)
call swept_cells_3d_device(s, p0, p1, e(5:1:-1), e, e, sd, status=st2)
call report('decreasing edges refused', merge(0_I4P, 1_I4P, st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT))
call s%analyze(aabb_tree_kind=AABB_TREE_OCTREE)
call swept_cells_3d_device(s, p0, p1, e, e, e, sd, status=st2)
call report('octree surface refused by the device path', merge(0_I4P, 1_I4P, st2 == STATUS_INVALID_INPUT))
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine compare(label, chunked)
   !< swept_cells_3d against swept_cells_3d_device (default chunk; with `chunked` also chunk_cells = 37) on s, e^3.
   character(*), intent(in) :: label   !< Case.
   logical,      intent(in) :: chunked !< Also the small chunks.
   integer(I4P)             :: n       !< Mismatches.

   call swept_cells_3d(s, p0, p1, e, e, e, sh, status=st)
   call swept_cells_3d_device(s, p0, p1, e, e, e, sd, status=st2)
   n = swept_cells_3d_mismatches(sh, sd) + abs(st - st2) + abs(st)
   if (chunked) then
      call swept_cells_3d_device(s, p0, p1, e, e, e, sd, status=st2, chunk_cells=37_I4P)
      n = n + swept_cells_3d_mismatches(sh, sd) + abs(st - st2)
   endif
   call report(label, n)
   endsubroutine compare

   subroutine compare_gauss(label)
   !< The Gauss path (exact = .false., nq = 2, with the 4-node error estimate), both chunk sizes.
   character(*), intent(in) :: label !< Case.
   integer(I4P)             :: n     !< Mismatches.

   call swept_cells_3d(s, p0, p1, e, e, e, sh, nq=2_I4P, exact=.false., status=st)
   call swept_cells_3d_device(s, p0, p1, e, e, e, sd, nq=2_I4P, exact=.false., status=st2)
   n = swept_cells_3d_mismatches(sh, sd) + abs(st - st2) + abs(st)
   call report(label, n)
   endsubroutine compare_gauss

   subroutine report(label, n)
   !< One line; a mismatch fails the test.
   character(*), intent(in) :: label !< Case.
   integer(I4P), intent(in) :: n     !< Mismatches.

   print '(A,A,I0)', label, ': mismatches ', n
   total = total + n
   if (n /= 0_I4P) ok = .false.
   endsubroutine report
endprogram fossil_test_swept_3d_device_twin
