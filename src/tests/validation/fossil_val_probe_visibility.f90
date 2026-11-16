!< VALIDATION (WP1 step 01 follow-up): the visibility of `probe_many` / `image_points` (the shared any-hit kernels of
!< `fossil_segment_visibility`) equals the issue #37 definition it replaced, on many meshes, run by
!< scripts/run_validation.sh at -O2.
program fossil_val_probe_visibility
!< VALIDATION (WP1 step 01 follow-up): the visibility of `probe_many` / `image_points` (the shared any-hit kernels of
!< `fossil_segment_visibility`) equals the issue #37 definition it replaced, on many meshes.
!<
!< 3-D (kernel `bvh_segment_blocked_oac`): the reference is the former host path `segment_clear`, rewritten here: a
!< ray from the foot point a along (b - a) / |b - a| (`intersect_ray_all`), blocked by a hit at tol < t <= |b - a|,
!< tol = 1e-9 of the box diagonal; clear when |b - a| <= tol. Meshes: bunny, dragon-fine, dragon (open), cube-01,
!< cube-05, naca0012-binary, cube_disconnect, a thin plate (cube-01 flattened to 1e-3 of its height) and a narrow
!< gap (two cube-01 copies 1e-3 of their width apart).
!< 2-D (kernel `curve_segment_blocked_oac`): the reference is brute force over every segment, the closed segment a-b
!< meets no segment but those at the foot (the closest segment, or the two at the closest vertex), exact
!< `segments_intersect_2d` (the former box-filtered `path_clear` tested a subset holding every meeting segment).
!< Curves: an ellipse with a square hole and a square island, a NACA 0012 (closed trailing edge, 240 vertices), two
!< thin plates 1e-3 apart.
!< Per mesh NP points (fixed Park-Miller seed) in the box enlarged by 10 % of its diagonal, every fifth with one
!< coordinate (every tenth with two) rounded to a coarse grid (points on facet planes and lines of the axis-aligned
!< bodies); per point 4 probes (sqrt(2) dx, -dx with dx = 0.02 diagonal, uniform lengths in +-0.4 and +-1.2 diagonal:
!< the long inward ones cross the far side of a convex body) and the image point. Passes iff every mismatch count is 0
!< and each mesh has blocked and visible probes.
use fossil,                only : surface_stl_object, facet_object, curve_object
use fossil_box_predicates, only : segments_intersect_2d
use fossil_ray_query,      only : ray_hit_t
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

integer(I4P), parameter         :: NP = 20000_I4P !< Points per mesh.
integer(I4P), parameter         :: M = 4_I4P      !< Probes per point.
character(len=*), parameter     :: FILES(7) = [character(len=32) :: 'src/tests/bunny.stl', 'src/tests/dragon-fine.stl', &
                                               'src/tests/dragon.stl', 'src/tests/cube-01.stl', 'src/tests/cube-05.stl', &
                                               'src/tests/naca0012-binary.stl', 'src/tests/cube_disconnect.stl']
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(vector_R8P)                :: bmn, bmx
real(R8P), allocatable          :: u(:), v(:)
integer(I4P), allocatable       :: nxt(:)
real(R8P)                       :: zmin, xmax, w, x, yt
integer(I4P)                    :: i, k, nf, total, mism
integer(I8P)                    :: seed
logical                         :: ok

ok = .true. ; total = 0_I4P ; mism = 0_I4P
seed = 20261115_I8P

do i = 1, size(FILES)
   call s%destroy
   call s%load_from_file(file_name=trim(FILES(i)), guess_format=.true.)
   call check_surface(trim(FILES(i)(11:)))
enddo

! thin plate: cube-01 flattened to 1e-3 of its height
call s%destroy
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
bmn = s%get_bmin() ; zmin = bmn%z ; bmx = s%get_bmax()
nf = s%get_facets_number()
allocate(f(nf))
do i = 1, nf
   fp => s%facet_at(i)
   f(i) = fp ; f(i)%vertex_id = 0_I4P
   do k = 1, 3
      f(i)%vertex(k)%z = zmin + 1.e-3_R8P * (f(i)%vertex(k)%z - zmin)
   enddo
   call f(i)%compute_metrix
enddo
call s%destroy
call s%adopt_facets(facets=f)
call check_surface('thin plate (cube-01, 1e-3 thick)')

! narrow gap: two cube-01 copies 1e-3 of their width apart along x
call s%destroy
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
bmn = s%get_bmin() ; bmx = s%get_bmax() ; xmax = bmx%x - bmn%x
allocate(f(2 * nf))
do i = 1, nf
   fp => s%facet_at(i)
   f(i) = fp ; f(i)%vertex_id = 0_I4P
   f(nf + i) = fp ; f(nf + i)%vertex_id = 0_I4P
   do k = 1, 3
      f(nf + i)%vertex(k)%x = f(nf + i)%vertex(k)%x + 1.001_R8P * xmax
   enddo
   call f(nf + i)%compute_metrix
enddo
call s%destroy
call s%adopt_facets(facets=f)
call check_surface('narrow gap (two cube-01, 1e-3 apart)')

! 2-D: ellipse (160 segments, CCW) + square hole (CW) + square island (CCW)
allocate(u(168), v(168), nxt(168))
do i = 1, 160
   u(i) = 0.6_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / 160) ; v(i) = 0.35_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / 160)
enddo
u(161:164) = [-0.2_R8P, -0.2_R8P, 0.2_R8P, 0.2_R8P] ; v(161:164) = [-0.1_R8P, 0.1_R8P, 0.1_R8P, -0.1_R8P]
u(165:168) = [0.8_R8P, 1.0_R8P, 1.0_R8P, 0.8_R8P] ; v(165:168) = [-0.1_R8P, -0.1_R8P, 0.1_R8P, 0.1_R8P]
call set_next([160_I4P, 4_I4P, 4_I4P])
call check_curve('ellipse + hole + island', [160_I4P, 4_I4P, 4_I4P])

! 2-D: NACA 0012, closed trailing edge, 240 vertices CCW (upper surface first, from the trailing edge)
deallocate(u, v, nxt) ; allocate(u(240), v(240), nxt(240))
do i = 1, 240
   w = 2._R8P * acos(-1._R8P) * (i - 1) / 240
   x = 0.5_R8P * (1._R8P + cos(w))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   u(i) = x ; v(i) = sign(1._R8P, sin(w)) * yt
   if (i == 1 .or. i == 121) v(i) = 0._R8P
enddo
call set_next([240_I4P])
call check_curve('NACA 0012', [240_I4P])

! 2-D: two thin plates [0,1] x [0,0.002] and [0,1] x [0.003,0.005] (CCW)
deallocate(u, v, nxt) ; allocate(u(8), v(8), nxt(8))
u = [0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P]
v = [0._R8P, 0._R8P, 0.002_R8P, 0.002_R8P, 0.003_R8P, 0.003_R8P, 0.005_R8P, 0.005_R8P]
call set_next([4_I4P, 4_I4P])
call check_curve('two thin plates 1e-3 apart', [4_I4P, 4_I4P])

print '(A,I0,A,I0)', 'all meshes: probes and images ', total, ', visibility mismatches ', mism
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine check_surface(name)
   !< NP points around s: probe/image visibility of the library against the ray-query reference.
   character(len=*), intent(in) :: name !< Mesh label.
   type(vector_R8P), allocatable :: p(:) !< Points.
   real(R8P), allocatable :: lengths(:,:), pp(:,:,:), ft(:,:), img(:,:), fi(:,:) !< Work.
   logical, allocatable   :: vis(:,:), ivis(:) !< Visibility.
   real(R8P)              :: lo(3), hi(3), diag, dx, tol !< Box, scales.
   integer(I4P)           :: i, k, nref, nb, nprobe !< Counters.

   allocate(p(NP), lengths(M, NP), pp(3, M, NP), ft(3, NP), img(3, NP), fi(3, NP), vis(M, NP), ivis(NP))
   bmn = s%get_bmin() ; bmx = s%get_bmax()
   lo = [bmn%x, bmn%y, bmn%z] ; hi = [bmx%x, bmx%y, bmx%z]
   diag = norm2(hi - lo) ; dx = 0.02_R8P * diag ; tol = 1.e-9_R8P * diag
   lo = lo - 0.1_R8P * diag ; hi = hi + 0.1_R8P * diag
   do i = 1, NP
      p(i) = vector_R8P(lo(1) + rnd() * (hi(1) - lo(1)), lo(2) + rnd() * (hi(2) - lo(2)), lo(3) + rnd() * (hi(3) - lo(3)))
      if (mod(i, 5_I4P) == 0_I4P) p(i)%y = grid(p(i)%y, lo(2), hi(2))
      if (mod(i, 10_I4P) == 0_I4P) p(i)%z = grid(p(i)%z, lo(3), hi(3))
      lengths(1:2, i) = [sqrt(2._R8P) * dx, -dx]
      lengths(3, i) = (rnd() - 0.5_R8P) * 0.8_R8P * diag
      lengths(4, i) = (rnd() - 0.5_R8P) * 2.4_R8P * diag
   enddo
   call s%probe_many(points=p, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft)
   call s%image_points(points=p, images=img, visible=ivis, foot_points=fi)
   nref = 0_I4P ; nb = count(.not. vis) + count(.not. ivis) ; nprobe = (M + 1_I4P) * NP
   do i = 1, NP
      do k = 1, M
         if (vis(k, i) .neqv. ray_visible(ft(:, i), pp(:, k, i), tol)) nref = nref + 1_I4P
      enddo
      if (ivis(i) .neqv. ray_visible(fi(:, i), img(:, i), tol)) nref = nref + 1_I4P
   enddo
   call report(name, s%get_facets_number(), nprobe, nb, nref)
   endsubroutine check_surface

   function ray_visible(a, b, tol) result(clear)
   !< The former `segment_clear`: no hit of the ray from a along (b - a) / |b - a| at tol < t <= |b - a|.
   real(R8P), intent(in)        :: a(3), b(3) !< Foot point, probe.
   real(R8P), intent(in)        :: tol        !< Hits at t <= tol (the facets at a) ignored.
   logical                      :: clear      !< Result.
   type(ray_hit_t), allocatable :: hits(:)    !< Ray hits.
   type(vector_R8P)             :: va, vb, d, dir !< Work.
   real(R8P)                    :: len        !< Segment length.
   integer(I4P)                 :: h          !< Counter.

   clear = .true.
   va = vector_R8P(a(1), a(2), a(3)) ; vb = vector_R8P(b(1), b(2), b(3))
   d = vb - va
   len = sqrt(d%x * d%x + d%y * d%y + d%z * d%z)
   if (len <= tol) return
   dir = d / len
   call s%intersect_ray_all(ray_origin=va, ray_direction=dir, hits=hits)
   do h = 1, size(hits)
      if (hits(h)%t > tol .and. hits(h)%t <= len) then
         clear = .false.
         return
      endif
   enddo
   endfunction ray_visible

   subroutine set_next(sizes)
   !< nxt(j): the vertex after j in its closed loop.
   integer(I4P), intent(in) :: sizes(:) !< Loop sizes.
   integer(I4P)             :: l, j, j0 !< Counters.

   j0 = 0_I4P
   do l = 1, size(sizes)
      do j = j0 + 1, j0 + sizes(l)
         nxt(j) = j + 1_I4P
      enddo
      nxt(j0 + sizes(l)) = j0 + 1_I4P
      j0 = j0 + sizes(l)
   enddo
   endsubroutine set_next

   subroutine check_curve(name, sizes)
   !< NP points around the curve (u, v): probe/image visibility of the library against brute force.
   character(len=*), intent(in) :: name     !< Curve label.
   integer(I4P),     intent(in) :: sizes(:) !< Loop sizes.
   type(curve_object)     :: c              !< Curve.
   real(R8P), allocatable :: p(:,:), lengths(:,:), pp(:,:,:), ft(:,:), img(:,:), fi(:,:), d(:) !< Work.
   integer(I4P), allocatable :: sg(:), vx(:) !< Closest segments and vertices.
   logical, allocatable   :: vis(:,:), ivis(:) !< Visibility.
   real(R8P)              :: lo(2), hi(2), diag, dx !< Box, scales.
   integer(I4P)           :: i, k, nref, nb, nprobe !< Counters.

   call c%adopt_loops(u=u, v=v, loop_sizes=sizes)
   allocate(p(2, NP), lengths(M, NP), pp(2, M, NP), ft(2, NP), img(2, NP), fi(2, NP), d(NP), sg(NP), vx(NP), &
            vis(M, NP), ivis(NP))
   lo = [minval(u), minval(v)] ; hi = [maxval(u), maxval(v)]
   diag = norm2(hi - lo) ; dx = 0.02_R8P * diag
   lo = lo - 0.1_R8P * diag ; hi = hi + 0.1_R8P * diag
   do i = 1, NP
      p(:, i) = [lo(1) + rnd() * (hi(1) - lo(1)), lo(2) + rnd() * (hi(2) - lo(2))]
      if (mod(i, 5_I4P) == 0_I4P) p(2, i) = grid(p(2, i), lo(2), hi(2))
      if (mod(i, 10_I4P) == 0_I4P) p(1, i) = grid(p(1, i), lo(1), hi(1))
      lengths(1:2, i) = [sqrt(2._R8P) * dx, -dx]
      lengths(3, i) = (rnd() - 0.5_R8P) * 0.8_R8P * diag
      lengths(4, i) = (rnd() - 0.5_R8P) * 2.4_R8P * diag
   enddo
   call c%probe_many(points=p, lengths=lengths, probe_points=pp, visible=vis, foot_points=ft)
   call c%image_points(points=p, images=img, visible=ivis, foot_points=fi)
   call c%query_many(points=p, distances=d, segments=sg, vertices=vx)
   nref = 0_I4P ; nb = count(.not. vis) + count(.not. ivis) ; nprobe = (M + 1_I4P) * NP
   do i = 1, NP
      do k = 1, M
         if (vis(k, i) .neqv. brute_visible(ft(:, i), pp(:, k, i), sg(i), vx(i))) nref = nref + 1_I4P
      enddo
      if (ivis(i) .neqv. brute_visible(fi(:, i), img(:, i), sg(i), vx(i))) nref = nref + 1_I4P
   enddo
   call report(name, size(u, kind=I4P), nprobe, nb, nref)
   endsubroutine check_curve

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
      if (segments_intersect_2d(a(1), a(2), b(1), b(2), u(j), v(j), u(nxt(j)), v(nxt(j)))) then
         clear = .false.
         return
      endif
   enddo
   endfunction brute_visible

   subroutine report(name, nel, nprobe, nb, nref)
   !< One line per mesh; updates the totals and the pass flag.
   character(len=*), intent(in) :: name   !< Mesh label.
   integer(I4P),     intent(in) :: nel    !< Facets or segments.
   integer(I4P),     intent(in) :: nprobe !< Probes + images.
   integer(I4P),     intent(in) :: nb     !< Blocked.
   integer(I4P),     intent(in) :: nref   !< Mismatches.

   print '(A,A,I0,A,I0,A,I0,A,I0)', name, ': elements ', nel, ', probes and images ', nprobe, ', blocked ', nb, &
         ', visibility mismatches ', nref
   total = total + nprobe ; mism = mism + nref
   ok = ok .and. nref == 0_I4P .and. nb > 0_I4P .and. nb < nprobe
   endsubroutine report

   pure function grid(y, lo, hi) result(g)
   !< y rounded to a multiple of the power of 2 just below (hi - lo) / 16: on the facet planes and lines of bodies with
   !< dyadic coordinates (the cubes, the thin plates).
   real(R8P), intent(in) :: y, lo, hi !< Value, range.
   real(R8P)             :: g         !< Rounded value.
   real(R8P)             :: h         !< Grid step (a power of 2).

   h = 2._R8P**floor(log((hi - lo) / 16._R8P) / log(2._R8P))
   g = real(nint(y / h, I8P), R8P) * h
   endfunction grid

   function rnd() result(r)
   !< Deterministic uniform in (0, 1) (Park-Miller, the generator of `fossil_test_cut_links_3d`).
   real(R8P) :: r !< Value.

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_val_probe_visibility
