!< FOSSIL test: the device probes and images of a curve equal the host ones bit for bit, and the shared visibility
!< kernel equals the exact brute-force definition (WP1 step 01).
program fossil_test_probe_device_twin_2d
!< FOSSIL test: the device probes and images of a curve equal the host ones bit for bit, and the shared visibility
!< kernel equals the exact brute-force definition (WP1 step 01).
!<
!< Built without OpenACC (this CI build) the device loops run on the host (the curve bench runs them on the GPU, `probe
!< device vs serial mismatches`).
!< 1. An ellipse (160 segments) with a square hole and a square island beside it, 400 points in the box, 4 lengths per
!<    point (sqrt(2) dx, -dx, 0.3, -0.5): the visibility of `probe_many` (the kernel `curve_segment_blocked_oac`) equals
!<    the issue #37 definition evaluated here by brute force (the closed segment foot-probe meets no segment but those
!<    at the foot, exact `segments_intersect_2d`) for every probe and image, with visible and blocked probes; the
!<    one-length form equals the m-length form column by column; `probe_many_device` / `image_points_device` equal
!<    the host bitwise (probe points, foot points, normals, patches, visibility, images), posed too.
!< 2. Two unit squares 0.2 apart: the probes of 0.1, -0.3 and 0.5 from (0.9, 0.4) are visible, visible and blocked on
!<    both paths; the image of (0.7, 0.4) is blocked.
!< 3. Refusals: a device call on a curve not on the device, or with mismatched sizes, is STATUS_INVALID_INPUT.
use fossil,                only : curve_object, rigid_pose_2d_t, pose_2d, STATUS_OK, STATUS_INVALID_INPUT
use fossil_box_predicates, only : segments_intersect_2d
use penf,                  only : I4P, I8P, R8P

implicit none

integer(I4P), parameter :: N = 400_I4P, M = 4_I4P, NE = 160_I4P
type(curve_object)      :: c, two
type(rigid_pose_2d_t)   :: pose
real(R8P), allocatable  :: u(:), v(:)
integer(I4P), allocatable :: nxt(:)
real(R8P)               :: p(2, N), pw(2, N), lengths(M, N), dx, d(N)
real(R8P)               :: pp(2, M, N), pd(2, M, N), ft(2, N), fd(2, N), nm(2, N), nd(2, N), p1(2, N)
real(R8P)               :: img(2, N), imd(2, N), fi(2, N), fid(2, N)
logical                 :: vis(M, N), vid(M, N), v1(N), ivis(N), ivid(N)
integer(I4P)            :: pat(N), pad(N), ipat(N), ipad(N), sg(N), vx(N), i, k, st, nblocked, nref
integer(I8P)            :: seed
logical                 :: ok(6)

ok = .true.
allocate(u(NE + 8), v(NE + 8))
do i = 1, NE ! ellipse, counter-clockwise
   u(i) = 0.6_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / NE) ; v(i) = 0.35_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / NE)
enddo
u(NE + 1:NE + 4) = [-0.2_R8P, -0.2_R8P, 0.2_R8P, 0.2_R8P] ; v(NE + 1:NE + 4) = [-0.1_R8P, 0.1_R8P, 0.1_R8P, -0.1_R8P] ! hole, CW
u(NE + 5:NE + 8) = [0.8_R8P, 1.0_R8P, 1.0_R8P, 0.8_R8P] ; v(NE + 5:NE + 8) = [-0.1_R8P, -0.1_R8P, 0.1_R8P, 0.1_R8P] ! island
call c%adopt_loops(u=u, v=v, loop_sizes=[NE, 4_I4P, 4_I4P])
allocate(nxt(NE + 8))
nxt = [(i + 1, i = 1, NE + 8)] ; nxt(NE) = 1 ; nxt(NE + 4) = NE + 1 ; nxt(NE + 8) = NE + 5
dx = 0.02_R8P
seed = 109_I8P
do i = 1, N
   p(:, i) = [-0.8_R8P + 1.9_R8P * rnd(), -0.5_R8P + rnd()]
   lengths(:, i) = [sqrt(2._R8P) * dx, -dx, 0.3_R8P, -0.5_R8P]
enddo

! 1. host kernel = the brute-force definition; one length = m lengths; device = host
call c%probe_many(points=p, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft, normals=nm, patches=pat, &
                  status=st)
call c%query_many(points=p, distances=d, segments=sg, vertices=vx)
ok(1) = st == STATUS_OK
nblocked = count(.not. vis) ; nref = 0
do i = 1, N
   do k = 1, M
      if (vis(k, i) .neqv. brute_visible(ft(:, i), pp(:, k, i), sg(i), vx(i))) nref = nref + 1
   enddo
enddo
call c%image_points(points=p, images=img, visible=ivis, foot_points=fi, patches=ipat)
do i = 1, N
   if (ivis(i) .neqv. brute_visible(fi(:, i), img(:, i), sg(i), vx(i))) nref = nref + 1
enddo
ok(1) = ok(1) .and. nref == 0 .and. nblocked > 0 .and. nblocked < M * N
print '(A,I0,A,I0,A,I0,A)', 'ellipse + hole + island: ', nref, ' visibility mismatches vs brute force (', nblocked, &
      ' of ', M * N, ' probes blocked)'
do k = 1, M
   call c%probe_many(points=p, lengths=lengths(k, :), probe_points=p1, visible=v1)
   ok(2) = ok(2) .and. same(p1, pp(:, k, :)) .and. all(v1 .eqv. vis(k, :))
enddo
print '(A,L1)', 'one length per point = m lengths, column by column: ', ok(2)
call c%enter_device
call c%probe_many_device(points=p, lengths=lengths, probe_points=pd, visible=vid, foot_points=fd, normals=nd, patches=pad, &
                         status=st)
call c%image_points_device(points=p, images=imd, visible=ivid, foot_points=fid, patches=ipad)
ok(3) = st == STATUS_OK .and. same(reshape(pp, [2, M * N]), reshape(pd, [2, M * N])) .and. all(vis .eqv. vid) .and. &
        same(ft, fd) .and. same(nm, nd) .and. all(pat == pad) .and. same(img, imd) .and. all(ivis .eqv. ivid) .and. &
        same(fi, fid) .and. all(ipat == ipad)
call c%probe_many_device(points=p, lengths=lengths(1, :), probe_points=p1, visible=v1)
ok(3) = ok(3) .and. same(p1, pp(:, 1, :)) .and. all(v1 .eqv. vis(1, :))
print '(A,L1)', 'device probes and images = host, bitwise: ', ok(3)
pose = pose_2d(0.47_R8P, translation=[0.13_R8P, -0.09_R8P])
do i = 1, N
   pw(:, i) = matmul(pose%r, p(:, i)) + pose%t
enddo
call c%probe_many(points=pw, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft, normals=nm, pose=pose)
call c%probe_many_device(points=pw, lengths=lengths, probe_points=pd, visible=vid, foot_points=fd, normals=nd, pose=pose)
call c%image_points(points=pw, images=img, visible=ivis, foot_points=fi, pose=pose)
call c%image_points_device(points=pw, images=imd, visible=ivid, foot_points=fid, pose=pose)
ok(4) = same(reshape(pp, [2, M * N]), reshape(pd, [2, M * N])) .and. all(vis .eqv. vid) .and. same(ft, fd) .and. &
        same(nm, nd) .and. same(img, imd) .and. all(ivis .eqv. ivid) .and. same(fi, fid) .and. count(.not. vis) > 0
print '(A,L1)', 'posed: device = host, bitwise: ', ok(4)
call c%exit_device

! 2. two squares 0.2 apart: a blocked probe and a blocked image on both paths
call two%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P,  1.2_R8P, 2.2_R8P, 2.2_R8P, 1.2_R8P], &
                     v=[0._R8P, 0._R8P, 1._R8P, 1._R8P,  0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
call two%enter_device
call two%probe_many_device(points=reshape([0.9_R8P, 0.4_R8P], [2, 1]), lengths=reshape([0.1_R8P, -0.3_R8P, 0.5_R8P], [3, 1]), &
                           probe_points=pd(:, 1:3, 1:1), visible=vid(1:3, 1:1))
call two%probe_many(points=reshape([0.9_R8P, 0.4_R8P], [2, 1]), lengths=reshape([0.1_R8P, -0.3_R8P, 0.5_R8P], [3, 1]), &
                    probe_points=pp(:, 1:3, 1:1), visible=vis(1:3, 1:1))
call two%image_points_device(points=reshape([0.9_R8P, 0.4_R8P, 0.7_R8P, 0.4_R8P], [2, 2]), images=imd(:, 1:2), &
                             visible=ivid(1:2))
ok(5) = all(vid(1:3, 1) .eqv. [.true., .true., .false.]) .and. all(vis(1:3, 1) .eqv. vid(1:3, 1)) .and. &
        same(pp(:, 1:3, 1), pd(:, 1:3, 1)) .and. ivid(1) .and. .not. ivid(2) .and. &
        all(abs(imd(:, 2) - [1.3_R8P, 0.4_R8P]) < 1.e-15_R8P)
print '(A,L1)', 'two squares: blocked probe and image, device = host: ', ok(5)
call two%exit_device

! 3. refusals
call two%probe_many_device(points=reshape([0.9_R8P, 0.4_R8P], [2, 1]), lengths=[0.1_R8P], probe_points=p1(:, 1:1), &
                           visible=v1(1:1), status=st)
ok(6) = st == STATUS_INVALID_INPUT
call two%enter_device
call two%probe_many_device(points=reshape([0.9_R8P, 0.4_R8P], [2, 1]), lengths=reshape([0.1_R8P, 0.2_R8P], [2, 1]), &
                           probe_points=pd(:, 1:1, 1:1), visible=vid(1:2, 1:1), status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call two%image_points_device(points=reshape([0.9_R8P, 0.4_R8P], [2, 1]), images=imd(:, 1:2), status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call two%exit_device
print '(A,L1)', 'refusals (not on the device, size mismatches): ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function brute_visible(a, b, s0, v0) result(clear)
   !< The issue #37 definition by brute force: the closed segment a-b meets no curve segment but those at the foot (the
   !< closest segment s0 when v0 = 0, else the two at the closest vertex v0).
   real(R8P),    intent(in) :: a(2), b(2) !< Foot point, probe.
   integer(I4P), intent(in) :: s0, v0     !< Closest segment and vertex.
   logical                  :: clear      !< Result.
   integer(I4P)             :: j          !< Counter.

   clear = .true.
   do j = 1, size(u, kind=I4P)
      if (v0 == 0_I4P .and. j == s0) cycle
      if (v0 > 0_I4P .and. (j == v0 .or. nxt(j) == v0)) cycle
      if (segments_intersect_2d(a(1), a(2), b(1), b(2), u(j), v(j), u(nxt(j)), v(nxt(j)))) clear = .false.
   enddo
   endfunction brute_visible

   pure function same(a, b) result(yes)
   !< Bitwise equality of two rank-2 arrays.
   real(R8P), intent(in) :: a(:,:), b(:,:) !< Arrays.
   logical               :: yes            !< Result.

   yes = all(transfer(a, [0_I8P]) == transfer(b, [0_I8P]))
   endfunction same

   function rnd() result(r)
   !< Deterministic uniform in (0, 1) (Park-Miller, the generator of `fossil_test_cut_links_3d`).
   real(R8P) :: r !< Value.

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_probe_device_twin_2d
