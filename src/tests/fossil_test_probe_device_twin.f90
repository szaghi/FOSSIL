!< FOSSIL test: the device probes and images equal the host ones bit for bit, and the shared visibility kernel equals the
!< former ray-query visibility (WP1 step 01).
program fossil_test_probe_device_twin
!< FOSSIL test: the device probes and images equal the host ones bit for bit, and the shared visibility kernel equals the
!< former ray-query visibility (WP1 step 01).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver bugs of
!< `probe_many_device` / `image_points_device` (sizes, the foot points of `query_many_device`, the per-probe loop, the
!< pose mapping), not device codegen (the distance bench does that on the GPU, `probe device vs serial mismatches`).
!< 1. The dragon (`dragon.stl`, open), 300 points in its box, 4 lengths per point (sqrt(2) dx, sqrt(3) dx, -dx and a
!<    long one that often crosses the body): the visibility of `probe_many` (the kernel `bvh_segment_blocked_oac`)
!<    equals the issue #37 definition evaluated here by `intersect_ray_all` (a crossing at tol < t <= |length|) for
!<    every probe and image, with both visible and blocked probes present; the one-length form equals the m-length
!<    form column by column; `probe_many_device` / `image_points_device` equal the host bitwise (probe points, foot
!<    points, normals, patches, visibility, images), posed too (a rotation and a translation).
!< 2. Two unit cubes 0.2 apart: the probes of 0.1, -0.3 and 0.5 from (0.9, 0.55, 0.4) are visible, visible and blocked
!<    (the 0.5 probe crosses the second cube) on both paths; the image of (0.7, 0.55, 0.4) is blocked.
!< 3. Refusals: a device call on a surface not on the device, or with mismatched sizes, is STATUS_INVALID_INPUT.
use fossil, only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle, STATUS_OK, &
                   STATUS_INVALID_INPUT
use fossil_ray_query, only : ray_hit_t
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none

integer(I4P), parameter         :: N = 300_I4P, M = 4_I4P
type(surface_stl_object)        :: s, two, cube
type(facet_object), allocatable :: f(:), g(:)
type(facet_object), pointer     :: fp
type(vector_R8P)                :: p(N), pb(N), bmn, bmx
type(rigid_pose_t)              :: pose
real(R8P)                       :: lengths(M, N), lo(3), hi(3), diag, dx, tol
real(R8P)                       :: pp(3, M, N), pd(3, M, N), ft(3, N), fd(3, N), nm(3, N), nd(3, N), p1(3, N)
real(R8P)                       :: img(3, N), imd(3, N), fi(3, N), fid(3, N)
logical                         :: vis(M, N), vid(M, N), v1(N), ivis(N), ivid(N)
integer(I4P)                    :: pat(N), pad(N), ipat(N), ipad(N), i, k, st, nblocked, nref, nf
integer(I8P)                    :: seed
logical                         :: ok(6)

ok = .true.
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
bmn = s%get_bmin() ; bmx = s%get_bmax()
lo = [bmn%x, bmn%y, bmn%z] ; hi = [bmx%x, bmx%y, bmx%z]
diag = norm2(hi - lo) ; dx = 0.02_R8P * diag ; tol = 1.e-9_R8P * diag
seed = 109_I8P
do i = 1, N
   p(i) = vector_R8P(lo(1) + rnd() * (hi(1) - lo(1)), lo(2) + rnd() * (hi(2) - lo(2)), lo(3) + rnd() * (hi(3) - lo(3)))
   lengths(:, i) = [sqrt(2._R8P) * dx, sqrt(3._R8P) * dx, -dx, (0.1_R8P + 0.3_R8P * rnd()) * diag]
enddo

! 1. host kernel = the issue #37 ray-query definition; one length = m lengths; device = host
call s%probe_many(points=p, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft, normals=nm, patches=pat, &
                  status=st)
ok(1) = st == STATUS_OK
nblocked = count(.not. vis) ; nref = 0
do i = 1, N
   do k = 1, M
      if (vis(k, i) .neqv. ray_visible(ft(:, i), pp(:, k, i))) nref = nref + 1
   enddo
enddo
call s%image_points(points=p, images=img, visible=ivis, foot_points=fi, patches=ipat)
do i = 1, N
   if (ivis(i) .neqv. ray_visible(fi(:, i), img(:, i))) nref = nref + 1
enddo
ok(1) = ok(1) .and. nref == 0 .and. nblocked > 0 .and. nblocked < M * N
print '(A,I0,A,I0,A,I0,A)', 'dragon: ', nref, ' visibility mismatches vs the ray query (', nblocked, ' of ', M * N, &
      ' probes blocked)'
do k = 1, M
   call s%probe_many(points=p, lengths=lengths(k, :), probe_points=p1, visible=v1)
   ok(2) = ok(2) .and. same2(p1, pp(:, k, :)) .and. all(v1 .eqv. vis(k, :))
enddo
print '(A,L1)', 'dragon: one length per point = m lengths, column by column: ', ok(2)
call s%enter_device
call s%probe_many_device(points=p, lengths=lengths, probe_points=pd, visible=vid, foot_points=fd, normals=nd, patches=pad, &
                         status=st)
call s%image_points_device(points=p, images=imd, visible=ivid, foot_points=fid, patches=ipad)
ok(3) = st == STATUS_OK .and. same3(pp, pd) .and. all(vis .eqv. vid) .and. same2(ft, fd) .and. same2(nm, nd) .and. &
        all(pat == pad) .and. same2(img, imd) .and. all(ivis .eqv. ivid) .and. same2(fi, fid) .and. all(ipat == ipad)
call s%probe_many_device(points=p, lengths=lengths(1, :), probe_points=p1, visible=v1)
ok(3) = ok(3) .and. same2(p1, pp(:, 1, :)) .and. all(v1 .eqv. vis(1, :))
print '(A,L1)', 'dragon: device probes and images = host, bitwise: ', ok(3)
pose = pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P, translation=[0.11_R8P, -0.07_R8P, 0.05_R8P])
do i = 1, N ! world points whose body-frame images are p, perturbed
   pb(i) = world(pose, p(i))
enddo
call s%probe_many(points=pb, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft, normals=nm, pose=pose)
call s%probe_many_device(points=pb, lengths=lengths, probe_points=pd, visible=vid, foot_points=fd, normals=nd, pose=pose)
call s%image_points(points=pb, images=img, visible=ivis, foot_points=fi, pose=pose)
call s%image_points_device(points=pb, images=imd, visible=ivid, foot_points=fid, pose=pose)
ok(4) = same3(pp, pd) .and. all(vis .eqv. vid) .and. same2(ft, fd) .and. same2(nm, nd) .and. same2(img, imd) .and. &
        all(ivis .eqv. ivid) .and. same2(fi, fid) .and. count(.not. vis) > 0
print '(A,L1)', 'dragon posed: device = host, bitwise: ', ok(4)
call s%exit_device

! 2. two cubes 0.2 apart: a blocked probe and a blocked image on both paths
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
nf = cube%get_facets_number()
allocate(f(nf), g(nf))
do i = 1, nf
   fp => cube%facet_at(i)
   f(i) = fp ; f(i)%vertex_id = 0_I4P
   g(i) = fp ; g(i)%vertex_id = 0_I4P
   do k = 1, 3
      g(i)%vertex(k)%x = g(i)%vertex(k)%x + 1.2_R8P
   enddo
   call g(i)%compute_metrix
enddo
f = [f, g]
call two%adopt_facets(facets=f)
call two%enter_device
call two%probe_many_device(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P)], &
                           lengths=reshape([0.1_R8P, -0.3_R8P, 0.5_R8P], [3, 1]), probe_points=pd(:, 1:3, 1:1), &
                           visible=vid(1:3, 1:1))
call two%probe_many(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P)], lengths=reshape([0.1_R8P, -0.3_R8P, 0.5_R8P], [3, 1]), &
                    probe_points=pp(:, 1:3, 1:1), visible=vis(1:3, 1:1))
call two%image_points_device(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P), vector_R8P(0.7_R8P, 0.55_R8P, 0.4_R8P)], &
                             images=imd(:, 1:2), visible=ivid(1:2))
ok(5) = all(vid(1:3, 1) .eqv. [.true., .true., .false.]) .and. all(vis(1:3, 1) .eqv. vid(1:3, 1)) .and. &
        same3(pp(:, 1:3, 1:1), pd(:, 1:3, 1:1)) .and. ivid(1) .and. .not. ivid(2) .and. &
        all(abs(imd(:, 2) - [1.3_R8P, 0.55_R8P, 0.4_R8P]) < 1.e-15_R8P)
print '(A,L1)', 'two cubes: blocked probe and image, device = host: ', ok(5)
call two%exit_device

! 3. refusals
call two%probe_many_device(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P)], lengths=[0.1_R8P], probe_points=p1(:, 1:1), &
                           visible=v1(1:1), status=st)
ok(6) = st == STATUS_INVALID_INPUT
call two%enter_device
call two%probe_many_device(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P)], lengths=reshape([0.1_R8P, 0.2_R8P], [2, 1]), &
                           probe_points=pd(:, 1:1, 1:1), visible=vid(1:2, 1:1), status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call two%image_points_device(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P)], images=imd(:, 1:2), status=st)
ok(6) = ok(6) .and. st == STATUS_INVALID_INPUT
call two%exit_device
print '(A,L1)', 'refusals (not on the device, size mismatches): ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function ray_visible(a, b) result(clear)
   !< The issue #37 definition: no crossing of the ray from a along (b - a)/|b - a| at tol < t <= |b - a|.
   real(R8P), intent(in)        :: a(3), b(3) !< Foot point, probe.
   logical                      :: clear      !< Result.
   type(ray_hit_t), allocatable :: hits(:)    !< Ray hits.
   type(vector_R8P)             :: va, vb, dir !< Work.
   real(R8P)                    :: len        !< Segment length.
   integer(I4P)                 :: h          !< Counter.

   clear = .true.
   va = vector_R8P(a(1), a(2), a(3)) ; vb = vector_R8P(b(1), b(2), b(3))
   len = sqrt((vb%x - va%x) * (vb%x - va%x) + (vb%y - va%y) * (vb%y - va%y) + (vb%z - va%z) * (vb%z - va%z))
   if (len <= tol) return
   dir = (vb - va) / len
   call s%intersect_ray_all(ray_origin=va, ray_direction=dir, hits=hits)
   do h = 1, size(hits)
      if (hits(h)%t > tol .and. hits(h)%t <= len) clear = .false.
   enddo
   endfunction ray_visible

   function world(pz, q) result(w)
   !< A world point near the body point q: q placed by the pose (its inverse of to_body), rounded.
   type(rigid_pose_t), intent(in) :: pz !< Pose.
   type(vector_R8P),   intent(in) :: q  !< Body point.
   type(vector_R8P)               :: w  !< World point.
   real(R8P)                      :: x(3) !< Work.

   x = matmul(pz%r, [q%x, q%y, q%z]) + pz%t
   w = vector_R8P(x(1), x(2), x(3))
   endfunction world

   pure function same2(a, b) result(yes)
   !< Bitwise equality of two rank-2 arrays.
   real(R8P), intent(in) :: a(:,:), b(:,:) !< Arrays.
   logical               :: yes            !< Result.

   yes = all(transfer(a, [0_I8P]) == transfer(b, [0_I8P]))
   endfunction same2

   pure function same3(a, b) result(yes)
   !< Bitwise equality of two rank-3 arrays.
   real(R8P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.
   logical               :: yes                !< Result.

   yes = all(transfer(a, [0_I8P]) == transfer(b, [0_I8P]))
   endfunction same3

   function rnd() result(r)
   !< Deterministic uniform in (0, 1) (Park-Miller, the generator of `fossil_test_cut_links_3d`).
   real(R8P) :: r !< Value.

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_probe_device_twin
