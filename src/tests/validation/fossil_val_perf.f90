!< FOSSIL validation: the performance bench (issue #94, roadmap #92 WP1).
!<
!< One line per kernel: group, kernel, problem size, items, t_min and t_median [s] of N timed repetitions after one
!< warm-up run, t_min per item [s], and a 64-bit hash (16 hex digits) of the bit patterns of every output the kernel
!< produces. The hash makes "bit-identical" a check, not an eye test: `scripts/run_perf.sh --compare A B` fails when a
!< hash differs. Wall clock by `system_clock` (int64 counts); single thread (validation-gnu has no -fopenmp, and
!< run_perf.sh sets OMP_NUM_THREADS=1).
!<
!< Arguments (any order): `reps=N` (default 5, or the environment variable FOSSIL_PERF_REPS; run_validation.sh sets 1),
!< `commit=SHA` and `host=NAME` for the header (run_perf.sh passes both), `only=g1,g2` to run some groups (profiling;
!< a baseline runs them all).
!<
!< The check: every kernel ran, its hash is nonzero, and every repetition (warm-up included) gave the same hash.
!< Baselines: src/tests/validation/perf/NAME.txt, written by scripts/run_perf.sh NAME (machine-specific).
module fossil_val_perf_kit
!< Hash, timer and report of the performance bench, and its kernels.
use fossil,                   only : surface_stl_object, curve_object, facet_object, rigid_pose_t, rigid_pose_2d_t, &
                                     pose_identity, pose_from_axis_angle, pose_2d, cut_cells_2d_t, cut_cells_3d_t, &
                                     cut_cells_3d, cut_links_t, swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, &
                                     swept_cells_3d, curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, &
                                     TOPO_BOUNDARY, wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, &
                                     donor_set_t, donor_candidates, donor_candidates_2d, overset_mask_t, classify, &
                                     classify_2d, SELECT_EXPLICIT
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none
private
public :: hash_t, hash_add, tic, toc, report, header_line, perf_ok
public :: k_distance_many, k_query_many, k_box_classify, k_curve_distance, k_body_cubes, k_body_large
public :: k_cut_2d, k_cut_3d, k_swept_2d, k_swept_3d, k_toi_3d, k_toi_cubes, k_toi_2d, k_remesh, k_decimate
public :: k_overset_2d, k_overset_3d, ellipse, scale_to_unit, centre, want, perf_only

integer(I8P), parameter :: P1 = 2147483647_I8P, P2 = 2147483629_I8P !< Primes < 2^31.
integer(I8P), parameter :: B1 = 1000003_I8P, B2 = 65599_I8P          !< Multipliers (h * B < 2^51: no overflow).
real(R8P),    parameter :: PI = acos(-1._R8P)                         !< Pi.
logical                 :: perf_ok = .true.                           !< All checks so far.
character(256)          :: perf_only = ''                             !< Groups to run, comma separated (all if empty).

type :: hash_t
   !< Two polynomial hashes mod primes below 2^31 over the 16-bit chunks of 64-bit words (no integer overflow).
   integer(I8P) :: a = 1_I8P !< Hash mod P1.
   integer(I8P) :: b = 7_I8P !< Hash mod P2.
endtype hash_t

interface hash_add
   module procedure hash_r0, hash_r1, hash_r2, hash_r3, hash_r4, hash_i0, hash_i1, hash_i2, hash_i3, hash_l0, hash_l1, &
                    hash_l2, hash_l3
endinterface hash_add

contains
   ! ------------------------------------------------------------------------------------------------- hash, timer, report
   pure subroutine hash_word(h, w)
   !< Add a 64-bit word: its four 16-bit chunks, low first.
   type(hash_t), intent(inout) :: h !< Hash.
   integer(I8P), intent(in)    :: w !< Word.
   integer(I8P)                :: c !< Chunk.
   integer(I4P)                :: k !< Counter.

   do k = 0, 3
      c = ibits(w, 16 * k, 16)
      h%a = mod(h%a * B1 + c + 1_I8P, P1) ; h%b = mod(h%b * B2 + c + 1_I8P, P2)
   enddo
   endsubroutine hash_word

   pure subroutine hash_r0(h, x)
   !< Add a real.
   type(hash_t), intent(inout) :: h !< Hash.
   real(R8P),    intent(in)    :: x !< Value.

   call hash_word(h, transfer(x, 0_I8P))
   endsubroutine hash_r0

   pure subroutine hash_r1(h, x)
   !< Add a real array, its size first.
   type(hash_t), intent(inout) :: h    !< Hash.
   real(R8P),    intent(in)    :: x(:) !< Values.
   integer(I4P)                :: i    !< Counter.

   call hash_word(h, int(size(x), I8P))
   do i = 1, size(x)
      call hash_word(h, transfer(x(i), 0_I8P))
   enddo
   endsubroutine hash_r1

   pure subroutine hash_r2(h, x)
   !< Add a rank-2 real array.
   type(hash_t), intent(inout) :: h      !< Hash.
   real(R8P),    intent(in)    :: x(:,:) !< Values.

   call hash_r1(h, reshape(x, [size(x)]))
   endsubroutine hash_r2

   pure subroutine hash_r3(h, x)
   !< Add a rank-3 real array.
   type(hash_t), intent(inout) :: h        !< Hash.
   real(R8P),    intent(in)    :: x(:,:,:) !< Values.

   call hash_r1(h, reshape(x, [size(x)]))
   endsubroutine hash_r3

   pure subroutine hash_r4(h, x)
   !< Add a rank-4 real array.
   type(hash_t), intent(inout) :: h          !< Hash.
   real(R8P),    intent(in)    :: x(:,:,:,:) !< Values.

   call hash_r1(h, reshape(x, [size(x)]))
   endsubroutine hash_r4

   pure subroutine hash_i0(h, x)
   !< Add an integer.
   type(hash_t), intent(inout) :: h !< Hash.
   integer(I4P), intent(in)    :: x !< Value.

   call hash_word(h, int(x, I8P))
   endsubroutine hash_i0

   pure subroutine hash_i1(h, x)
   !< Add an integer array, its size first.
   type(hash_t), intent(inout) :: h    !< Hash.
   integer(I4P), intent(in)    :: x(:) !< Values.
   integer(I4P)                :: i    !< Counter.

   call hash_word(h, int(size(x), I8P))
   do i = 1, size(x)
      call hash_word(h, int(x(i), I8P))
   enddo
   endsubroutine hash_i1

   pure subroutine hash_i2(h, x)
   !< Add a rank-2 integer array.
   type(hash_t), intent(inout) :: h      !< Hash.
   integer(I4P), intent(in)    :: x(:,:) !< Values.

   call hash_i1(h, reshape(x, [size(x)]))
   endsubroutine hash_i2

   pure subroutine hash_i3(h, x)
   !< Add a rank-3 integer array.
   type(hash_t), intent(inout) :: h        !< Hash.
   integer(I4P), intent(in)    :: x(:,:,:) !< Values.

   call hash_i1(h, reshape(x, [size(x)]))
   endsubroutine hash_i3

   pure subroutine hash_l0(h, x)
   !< Add a logical.
   type(hash_t), intent(inout) :: h !< Hash.
   logical,      intent(in)    :: x !< Value.

   call hash_word(h, merge(1_I8P, 0_I8P, x))
   endsubroutine hash_l0

   pure subroutine hash_l1(h, x)
   !< Add a logical array, its size first.
   type(hash_t), intent(inout) :: h    !< Hash.
   logical,      intent(in)    :: x(:) !< Values.

   call hash_i1(h, merge(1_I4P, 0_I4P, x))
   endsubroutine hash_l1

   pure subroutine hash_l2(h, x)
   !< Add a rank-2 logical array.
   type(hash_t), intent(inout) :: h      !< Hash.
   logical,      intent(in)    :: x(:,:) !< Values.

   call hash_l1(h, reshape(x, [size(x)]))
   endsubroutine hash_l2

   pure subroutine hash_l3(h, x)
   !< Add a rank-3 logical array.
   type(hash_t), intent(inout) :: h        !< Hash.
   logical,      intent(in)    :: x(:,:,:) !< Values.

   call hash_l1(h, reshape(x, [size(x)]))
   endsubroutine hash_l3

   function hash_hex(h) result(s)
   !< 16 hex digits.
   type(hash_t), intent(in) :: h !< Hash.
   character(16)            :: s !< Digits.

   write(s, '(Z8.8,Z8.8)') h%a, h%b
   endfunction hash_hex

   function tic() result(c)
   !< Clock count now.
   integer(I8P) :: c !< Count.

   call system_clock(c)
   endfunction tic

   function toc(c0) result(t)
   !< Seconds since the count c0.
   integer(I8P), intent(in) :: c0 !< Start count.
   real(R8P)                :: t  !< Seconds.
   integer(I8P)             :: c, rate !< Clock.

   call system_clock(c, rate)
   t = real(c - c0, R8P) / real(rate, R8P)
   endfunction toc

   subroutine header_line
   !< The column titles, in the columns of `report`.
   character(12) :: cg = '# group' !< Group column.
   character(28) :: ck = 'kernel'  !< Kernel column.
   character(52) :: cp = 'problem size' !< Problem column.

   print '(A,1X,A,1X,A,1X,A10,3(1X,A11),1X,A)', cg, ck, cp, 'items', 't_min[s]', 't_median[s]', 't_min/item', 'hash'
   endsubroutine header_line

   subroutine report(group, kernel, problem, items, times, hashes)
   !< Print one kernel line; check the hashes (nonzero, equal over the repetitions, warm-up = hashes(0) included).
   character(*), intent(in) :: group, kernel, problem !< Labels.
   integer(I8P), intent(in) :: items                  !< Items (calls, cells, points) per repetition.
   real(R8P),    intent(in) :: times(:)               !< Timed repetitions.
   type(hash_t), intent(in) :: hashes(0:)             !< Hashes, warm-up first.
   real(R8P), allocatable   :: s(:)                   !< Sorted times.
   real(R8P)                :: w, tmed                !< Work, median.
   integer(I4P)             :: i, j, n                !< Counters.
   logical                  :: same                   !< Hashes agree.
   character(12)            :: cg                     !< Group column.
   character(28)            :: ck                     !< Kernel column.
   character(52)            :: cp                     !< Problem column.

   s = times ; n = size(s)
   do i = 2, n
      w = s(i) ; j = i - 1
      do while (j >= 1)
         if (s(j) <= w) exit
         s(j + 1) = s(j) ; j = j - 1
      enddo
      s(j + 1) = w
   enddo
   if (mod(n, 2) == 1) then
      tmed = s((n + 1) / 2)
   else
      tmed = 0.5_R8P * (s(n / 2) + s(n / 2 + 1))
   endif
   same = .true.
   do i = 1, ubound(hashes, 1)
      same = same .and. hashes(i)%a == hashes(0)%a .and. hashes(i)%b == hashes(0)%b
   enddo
   cg = group ; ck = kernel ; cp = problem ! left-justified fixed columns
   print '(A,1X,A,1X,A,1X,I10,3(1X,ES11.4),1X,A16)', cg, ck, cp, items, s(1), tmed, &
         s(1) / real(max(items, 1_I8P), R8P), hash_hex(hashes(0))
   if (.not. same) print '(A)', '#  FAIL: the hash differs between repetitions of ' // trim(kernel)
   if (hashes(0)%a == 0_I8P .and. hashes(0)%b == 0_I8P) print '(A)', '#  FAIL: zero hash of ' // trim(kernel)
   perf_ok = perf_ok .and. same .and. .not. (hashes(0)%a == 0_I8P .and. hashes(0)%b == 0_I8P)
   endsubroutine report

   ! ---------------------------------------------------------------------------------------------------------- helpers
   subroutine grid_points(s, m, pts, lo, hi)
   !< The m^3 cell centres of the bounding box of s grown by 10 %, and the cells as boxes.
   type(surface_stl_object),      intent(in)  :: s        !< Surface.
   integer(I4P),                  intent(in)  :: m        !< Cells per direction.
   type(vector_R8P), allocatable, intent(out) :: pts(:)   !< Cell centres [m^3].
   real(R8P),        allocatable, intent(out) :: lo(:,:), hi(:,:) !< Cell boxes [3, m^3].
   type(vector_R8P)                           :: bmin, bmax !< Box.
   real(R8P)                                  :: o(3), h(3) !< Origin, spacing.
   integer(I4P)                               :: i, j, k, l !< Counters.

   bmin = s%get_bmin() ; bmax = s%get_bmax()
   h = 1.2_R8P * [bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z] / m
   o = [bmin%x, bmin%y, bmin%z] - 0.1_R8P * [bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z]
   allocate(pts(m ** 3), lo(3, m ** 3), hi(3, m ** 3))
   l = 0
   do k = 1, m ; do j = 1, m ; do i = 1, m
      l = l + 1
      lo(:, l) = o + h * [i - 1, j - 1, k - 1] ; hi(:, l) = lo(:, l) + h
      pts(l) = vector_R8P(lo(1, l) + 0.5_R8P * h(1), lo(2, l) + 0.5_R8P * h(2), lo(3, l) + 0.5_R8P * h(3))
   enddo ; enddo ; enddo
   endsubroutine grid_points

   subroutine ellipse(c, n, ra, rb, cx, cy)
   !< Closed ellipse of n segments, counter-clockwise (solid inside).
   type(curve_object), intent(inout) :: c              !< Curve.
   integer(I4P),       intent(in)    :: n              !< Segments.
   real(R8P),          intent(in)    :: ra, rb, cx, cy !< Semi-axes, centre.
   real(R8P)                         :: u(n), v(n)     !< Vertices.
   integer(I4P)                      :: i              !< Counter.

   do i = 1, n
      u(i) = cx + ra * cos(2 * PI * (i - 1) / n) ; v(i) = cy + rb * sin(2 * PI * (i - 1) / n)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[n])
   endsubroutine ellipse

   function want(group) result(yes)
   !< The group is selected (`only=` argument, comma-separated groups; all by default).
   character(*), intent(in) :: group !< Group.
   logical                  :: yes   !< Selected.

   yes = len_trim(perf_only) == 0 .or. index(',' // trim(perf_only) // ',', ',' // group // ',') > 0
   endfunction want

   subroutine centre(sf)
   !< Bounding-box centre to the origin.
   type(surface_stl_object), intent(inout) :: sf     !< Body.
   type(vector_R8P)                        :: lo, hi !< Box.

   lo = sf%get_bmin() ; hi = sf%get_bmax()
   call sf%translate(delta=-0.5_R8P * (lo + hi)) ; call sf%analyze
   endsubroutine centre

   subroutine scale_to_unit(sf)
   !< Centre at the origin, unit bounding-box diagonal.
   type(surface_stl_object), intent(inout) :: sf     !< Body.
   type(vector_R8P)                        :: lo, hi !< Box.
   real(R8P)                               :: f      !< Factor.

   lo = sf%get_bmin() ; hi = sf%get_bmax()
   call sf%translate(delta=-0.5_R8P * (lo + hi))
   f = 1._R8P / norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call sf%resize(factor=vector_R8P(f, f, f)) ; call sf%analyze
   endsubroutine scale_to_unit

   subroutine hash_surface(h, s)
   !< Add the facets' vertices of a surface.
   type(hash_t),             intent(inout) :: h !< Hash.
   type(surface_stl_object), intent(in)    :: s !< Surface.
   type(facet_object), pointer             :: f !< Facet.
   integer(I4P)                            :: i, k !< Counters.

   call hash_add(h, s%get_facets_number())
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      do k = 1, 3
         call hash_add(h, [f%vertex(k)%x, f%vertex(k)%y, f%vertex(k)%z])
      enddo
   enddo
   call hash_add(h, s%get_volume())
   endsubroutine hash_surface

   pure subroutine prng(seed, r)
   !< Park-Miller draws in (0, 1), fixed seed (a subroutine: one draw per statement would be a side effect).
   integer(I8P), intent(inout) :: seed !< State.
   real(R8P),    intent(out)   :: r(:) !< Draws.
   integer(I4P)                :: i    !< Counter.

   do i = 1, size(r)
      seed = mod(seed * 48271_I8P, 2147483647_I8P)
      r(i) = real(seed, R8P) / 2147483647._R8P
   enddo
   endsubroutine prng

   ! ---------------------------------------------------------------------------------------------------------- queries
   subroutine k_distance_many(s, label, m, nrep)
   !< Signed distance batch (pseudo-normal, square root) at the m^3 cell centres.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Surface label.
   integer(I4P),             intent(in) :: m, nrep !< Grid, repetitions.
   type(vector_R8P), allocatable        :: pts(:)  !< Points.
   real(R8P),        allocatable        :: lo(:,:), hi(:,:), d(:) !< Boxes, distances.
   real(R8P)                            :: t(nrep) !< Times.
   type(hash_t)                         :: h(0:nrep) !< Hashes.
   integer(I8P)                         :: c0      !< Clock.
   integer(I4P)                         :: r       !< Counter.
   character(64)                        :: sz      !< Size.

   call grid_points(s, m, pts, lo, hi) ; allocate(d(size(pts)))
   do r = 0, nrep
      c0 = tic()
      call s%distance_many(points=pts, distances=d, is_signed=.true., is_square_root=.true.)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), d)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, s%get_facets_number(), ' facets, ', m, '^3 pts'
   call report('queries', 'distance_many', sz, int(size(pts), I8P), t, h)
   endsubroutine k_distance_many

   subroutine k_query_many(s, label, m, nrep)
   !< Rich query batch (closest point, facet, feature, normal, body) at the m^3 cell centres.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Surface label.
   integer(I4P),             intent(in) :: m, nrep !< Grid, repetitions.
   type(vector_R8P), allocatable        :: pts(:)  !< Points.
   real(R8P),        allocatable        :: lo(:,:), hi(:,:), d(:), cp(:,:), nr(:,:) !< Outputs.
   integer(I4P),     allocatable        :: fc(:), ft(:), fi(:), bd(:) !< Outputs.
   real(R8P)                            :: t(nrep) !< Times.
   type(hash_t)                         :: h(0:nrep) !< Hashes.
   integer(I8P)                         :: c0      !< Clock.
   integer(I4P)                         :: r, n    !< Counters.
   character(64)                        :: sz      !< Size.

   call grid_points(s, m, pts, lo, hi) ; n = size(pts)
   allocate(d(n), cp(3, n), nr(3, n), fc(n), ft(n), fi(n), bd(n))
   do r = 0, nrep
      c0 = tic()
      call s%query_many(points=pts, distances=d, closest_points=cp, facets=fc, features=ft, feature_indexes=fi, &
                        normals=nr, bodies=bd)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), d) ; call hash_add(h(r), cp) ; call hash_add(h(r), nr) ; call hash_add(h(r), fc)
      call hash_add(h(r), ft) ; call hash_add(h(r), fi) ; call hash_add(h(r), bd)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, s%get_facets_number(), ' facets, ', m, '^3 pts'
   call report('queries', 'query_many', sz, int(n, I8P), t, h)
   endsubroutine k_query_many

   subroutine k_box_classify(s, label, m, nrep)
   !< Box classification (default FWN centre sign, near radius = one cell) of the m^3 cells.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Surface label.
   integer(I4P),             intent(in) :: m, nrep !< Grid, repetitions.
   type(vector_R8P), allocatable        :: pts(:)  !< Points.
   real(R8P),        allocatable        :: lo(:,:), hi(:,:) !< Boxes.
   integer(I4P),     allocatable        :: codes(:) !< Codes.
   logical,          allocatable        :: inside(:) !< Centre signs.
   real(R8P)                            :: t(nrep) !< Times.
   type(hash_t)                         :: h(0:nrep) !< Hashes.
   integer(I8P)                         :: c0      !< Clock.
   integer(I4P)                         :: r, n    !< Counters.
   character(64)                        :: sz      !< Size.

   call grid_points(s, m, pts, lo, hi) ; n = size(pts) ; allocate(codes(n), inside(n))
   do r = 0, nrep
      c0 = tic()
      call s%box_classify_many(lo=lo, hi=hi, codes=codes, radius=hi(1, 1) - lo(1, 1), inside=inside)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), codes) ; call hash_add(h(r), inside)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, s%get_facets_number(), ' facets, ', m, '^3 boxes'
   call report('queries', 'box_classify_many', sz, int(n, I8P), t, h)
   endsubroutine k_box_classify

   subroutine k_curve_distance(nseg, m, nrep)
   !< 2-D signed distance batch: an ellipse of nseg segments, m^2 points.
   integer(I4P), intent(in) :: nseg, m, nrep  !< Segments, grid, repetitions.
   type(curve_object)       :: c              !< Curve.
   real(R8P), allocatable   :: p(:,:), d(:)   !< Points, distances.
   real(R8P)                :: t(nrep)        !< Times.
   type(hash_t)             :: h(0:nrep)      !< Hashes.
   integer(I8P)             :: c0             !< Clock.
   integer(I4P)             :: r, i, j        !< Counters.
   character(64)            :: sz             !< Size.

   call ellipse(c, nseg, 0.4_R8P, 0.25_R8P, 0.5_R8P, 0.5_R8P)
   allocate(p(2, m * m), d(m * m))
   do j = 1, m ; do i = 1, m
      p(:, i + (j - 1) * m) = [(i - 0.5_R8P) / m, (j - 0.5_R8P) / m]
   enddo ; enddo
   do r = 0, nrep
      c0 = tic()
      call c%distance_many(points=p, distances=d)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), d)
   enddo
   write(sz, '(A,I0,A,I0,A)') 'ellipse ', nseg, ' segments, ', m, '^2 pts'
   call report('queries', 'curve.distance_many', sz, int(m * m, I8P), t, h)
   endsubroutine k_curve_distance

   ! ---------------------------------------------------------------------------------------------------- body distance
   subroutine k_body_cubes(ncalls, nrep)
   !< distance_to of two 12-facet cubes, ncalls poses of the first.
   integer(I4P), intent(in)     :: ncalls, nrep !< Calls, repetitions.
   type(surface_stl_object)     :: a, b         !< Cubes.
   type(rigid_pose_t), allocatable :: ps(:)     !< Poses.
   real(R8P),    allocatable    :: d(:), pa(:,:), pb(:,:) !< Outputs.
   integer(I4P), allocatable    :: fa(:), fb(:) !< Outputs.
   real(R8P)                    :: t(nrep)      !< Times.
   type(hash_t)                 :: h(0:nrep)    !< Hashes.
   integer(I8P)                 :: c0, seed     !< Clock, seed.
   real(R8P)                    :: w(6)         !< Draws.
   integer(I4P)                 :: r, k         !< Counters.
   character(64)                :: sz           !< Size.

   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%translate(x=1.8_R8P, y=0.3_R8P, z=0.1_R8P) ; call b%analyze
   allocate(ps(ncalls), d(ncalls), pa(3, ncalls), pb(3, ncalls), fa(ncalls), fb(ncalls))
   seed = 20261109_I8P
   do k = 1, ncalls
      call prng(seed, w)
      ps(k) = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [0.5_R8P * w(5), 0.2_R8P * w(6), 0._R8P])
   enddo
   do r = 0, nrep
      c0 = tic()
      do k = 1, ncalls
         call a%distance_to(b, d(k), point_self=pa(:, k), point_other=pb(:, k), facet_self=fa(k), facet_other=fb(k), &
                            pose_self=ps(k))
      enddo
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), d) ; call hash_add(h(r), pa) ; call hash_add(h(r), pb) ; call hash_add(h(r), fa)
      call hash_add(h(r), fb)
   enddo
   write(sz, '(A,I0,A)') 'cube vs cube, 12 facets each, ', ncalls, ' poses'
   call report('body', 'distance_to.cubes', sz, int(ncalls, I8P), t, h)
   endsubroutine k_body_cubes

   subroutine k_body_large(s, label, ncalls, nrep)
   !< distance_to of a surface and a translated copy, ncalls poses of the copy.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Label.
   integer(I4P),             intent(in) :: ncalls, nrep !< Calls, repetitions.
   type(surface_stl_object)     :: b            !< Copy.
   type(rigid_pose_t), allocatable :: ps(:)     !< Poses.
   type(vector_R8P)             :: lo, hi       !< Box.
   real(R8P),    allocatable    :: d(:), pa(:,:), pb(:,:) !< Outputs.
   integer(I4P), allocatable    :: fa(:), fb(:) !< Outputs.
   real(R8P)                    :: t(nrep), ext(3) !< Times, extent.
   type(hash_t)                 :: h(0:nrep)    !< Hashes.
   integer(I8P)                 :: c0, seed     !< Clock, seed.
   real(R8P)                    :: w(6)         !< Draws.
   integer(I4P)                 :: r, k         !< Counters.
   character(64)                :: sz           !< Size.

   b = s
   lo = s%get_bmin() ; hi = s%get_bmax() ; ext = [hi%x - lo%x, hi%y - lo%y, hi%z - lo%z]
   allocate(ps(ncalls), d(ncalls), pa(3, ncalls), pb(3, ncalls), fa(ncalls), fb(ncalls))
   seed = 20261110_I8P
   do k = 1, ncalls
      call prng(seed, w(1:5))
      ps(k) = pose_from_axis_angle(w(1:3) - 0.5_R8P, 0.5_R8P * w(4), [0.8_R8P * ext(1), 0.2_R8P * ext(2) * w(5), 0._R8P])
   enddo
   do r = 0, nrep
      c0 = tic()
      do k = 1, ncalls
         call s%distance_to(b, d(k), point_self=pa(:, k), point_other=pb(:, k), facet_self=fa(k), facet_other=fb(k), &
                            pose_other=ps(k))
      enddo
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), d) ; call hash_add(h(r), pa) ; call hash_add(h(r), pb) ; call hash_add(h(r), fa)
      call hash_add(h(r), fb)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, s%get_facets_number(), ' facets vs copy, ', ncalls, ' poses'
   call report('body', 'distance_to.large', sz, int(ncalls, I8P), t, h)
   endsubroutine k_body_large

   ! ------------------------------------------------------------------------------------------------------- cut cells
   subroutine k_cut_2d(nseg, m, nrep)
   !< 2-D cut cells of an ellipse of nseg segments on m^2 cells (per cut cell too).
   integer(I4P), intent(in) :: nseg, m, nrep !< Segments, cells, repetitions.
   type(curve_object)       :: c             !< Curve.
   type(cut_cells_2d_t)     :: cut           !< Cut cells.
   real(R8P), allocatable   :: x(:)          !< Cell edges.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0, ncut      !< Clock, cut cells.
   integer(I4P)             :: r, i          !< Counters.
   character(64)            :: sz            !< Size.

   call ellipse(c, nseg, 0.4_R8P, 0.25_R8P, 0.5_R8P, 0.5_R8P)
   x = [(0.0031_R8P + (1._R8P - 0.0062_R8P) * i / m, i = 0, m)]
   do r = 0, nrep
      c0 = tic()
      call c%cut_cells(x=x, y=x, cut=cut)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), cut%alpha) ; call hash_add(h(r), cut%aperture) ; call hash_add(h(r), cut%centroid)
      call hash_add(h(r), cut%wet_length) ; call hash_add(h(r), cut%wet_normal) ; call hash_add(h(r), cut%pieces)
      call hash_add(h(r), cut%piece_area) ; call hash_add(h(r), cut%piece_centroid)
      call hash_add(h(r), cut%piece_iv_lo) ; call hash_add(h(r), cut%piece_iv_hi)
   enddo
   ncut = count(cut%wet_length > 0._R8P)
   write(sz, '(A,I0,A,I0,A,I0,A)') 'ellipse ', nseg, ' segments, ', m, '^2 cells (', ncut, ' cut)'
   call report('cut', 'cut_cells_2d', sz, ncut, t, h)
   endsubroutine k_cut_2d

   subroutine k_cut_3d(s, label, m, with_links, nrep)
   !< 3-D cut cells (pieces) of a closed body on m^3 cells over its box grown by 10 %; optionally the links.
   type(surface_stl_object), intent(in) :: s     !< Body.
   character(*),             intent(in) :: label !< Label.
   integer(I4P),             intent(in) :: m, nrep !< Cells, repetitions.
   logical,                  intent(in) :: with_links !< Links too.
   type(cut_cells_3d_t)     :: cut           !< Cut cells.
   type(cut_links_t)        :: lk            !< Links.
   type(vector_R8P)         :: lo, hi        !< Box.
   real(R8P), allocatable   :: x(:), y(:), z(:) !< Cell edges.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0, ncut      !< Clock, cut cells.
   integer(I4P)             :: r, i          !< Counters.
   character(64)            :: sz            !< Size.

   lo = s%get_bmin() ; hi = s%get_bmax()
   x = [(lo%x - 0.1_R8P * (hi%x - lo%x) + 1.2_R8P * (hi%x - lo%x) * i / m, i = 0, m)]
   y = [(lo%y - 0.1_R8P * (hi%y - lo%y) + 1.2_R8P * (hi%y - lo%y) * i / m, i = 0, m)]
   z = [(lo%z - 0.1_R8P * (hi%z - lo%z) + 1.2_R8P * (hi%z - lo%z) * i / m, i = 0, m)]
   do r = 0, nrep
      c0 = tic()
      if (with_links) then
         call cut_cells_3d(s, x, y, z, cut, links=lk)
      else
         call cut_cells_3d(s, x, y, z, cut)
      endif
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), cut%alpha) ; call hash_add(h(r), cut%face_x) ; call hash_add(h(r), cut%face_y)
      call hash_add(h(r), cut%face_z) ; call hash_add(h(r), cut%centroid) ; call hash_add(h(r), cut%wet_area)
      call hash_add(h(r), cut%wet_normal) ; call hash_add(h(r), cut%pieces) ; call hash_add(h(r), cut%piece_volume)
      call hash_add(h(r), cut%piece_centroid) ; call hash_add(h(r), cut%piece_face) ; call hash_add(h(r), cut%node_in)
      if (with_links) then
         call hash_add(h(r), lk%start) ; call hash_add(h(r), lk%piece) ; call hash_add(h(r), lk%face)
         call hash_add(h(r), lk%area) ; call hash_add(h(r), lk%unmatched)
      endif
   enddo
   ncut = count(cut%wet_area > 0._R8P)
   write(sz, '(A,1X,I0,A,I0,A,I0,A)') label, s%get_facets_number(), ' facets, ', m, '^3 cells (', ncut, ' cut)'
   if (with_links) then
      call report('cut', 'cut_cells_3d.links', sz, ncut, t, h)
   else
      call report('cut', 'cut_cells_3d', sz, ncut, t, h)
   endif
   endsubroutine k_cut_3d

   ! ----------------------------------------------------------------------------------------------------- swept cells
   subroutine k_swept_2d(nseg, m, rotation, nrep)
   !< One step of swept_cells_2d of an ellipse of nseg segments on m^2 cells: a translation or a rotation.
   integer(I4P), intent(in) :: nseg, m, nrep !< Segments, cells, repetitions.
   logical,      intent(in) :: rotation      !< Rotation step, else translation.
   type(curve_object)       :: c             !< Curve.
   type(swept_cells_2d_t)   :: sw            !< Swept geometry.
   type(rigid_pose_2d_t)    :: p0, p1        !< Poses.
   real(R8P), allocatable   :: x(:)          !< Cell edges.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r, i          !< Counters.
   character(64)            :: sz            !< Size.

   call ellipse(c, nseg, 0.4_R8P, 0.25_R8P, 0._R8P, 0._R8P)
   x = [(-0.5_R8P + 0.0013_R8P + (1._R8P - 0.0026_R8P) * i / m, i = 0, m)]
   p0 = pose_2d(angle=0._R8P)
   if (rotation) then
      p1 = pose_2d(angle=PI / 18._R8P)
   else
      p1 = pose_2d(angle=0._R8P, translation=[1.7_R8P / m, 0.6_R8P / m])
   endif
   do r = 0, nrep
      c0 = tic()
      call swept_cells_2d(c, p0, p1, x, x, sw)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), sw%alpha0) ; call hash_add(h(r), sw%alpha1) ; call hash_add(h(r), sw%swept)
      call hash_add(h(r), sw%face_u) ; call hash_add(h(r), sw%face_v) ; call hash_add(h(r), sw%wall)
      call hash_add(h(r), sw%residual) ; call hash_add(h(r), sw%exact) ; call hash_add(h(r), sw%evaluations)
   enddo
   write(sz, '(A,I0,A,I0,A)') 'ellipse ', nseg, ' segments, ', m, '^2 cells, 1 step'
   if (rotation) then
      call report('swept', 'swept_cells_2d.rotation', sz, 1_I8P, t, h)
   else
      call report('swept', 'swept_cells_2d.translation', sz, 1_I8P, t, h)
   endif
   endsubroutine k_swept_2d

   subroutine k_swept_3d(s, label, m, rotation, nrep)
   !< One step of swept_cells_3d of a body (centred at the origin) on m^3 cells: a translation or a rotation.
   type(surface_stl_object), intent(in) :: s     !< Body.
   character(*),             intent(in) :: label !< Label.
   integer(I4P),             intent(in) :: m, nrep !< Cells, repetitions.
   logical,                  intent(in) :: rotation !< Rotation step, else translation.
   type(swept_cells_3d_t)   :: sw            !< Swept geometry.
   type(rigid_pose_t)       :: p0, p1        !< Poses.
   real(R8P), allocatable   :: x(:)          !< Cell edges.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r, i          !< Counters.
   character(64)            :: sz            !< Size.

   x = [(-0.8_R8P + 0.0017_R8P + (1.6_R8P - 0.0034_R8P) * i / m, i = 0, m)]
   p0 = pose_identity()
   if (rotation) then
      p1 = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0.3_R8P) ! about an axis through (0.2, 0.1, 0)
      p1%t = [0.2_R8P, 0.1_R8P, 0._R8P] - matmul(p1%r, [0.2_R8P, 0.1_R8P, 0._R8P])
   else
      p1 = pose_identity() ; p1%t = [1.7_R8P, 0.6_R8P, -0.4_R8P] * (1.6_R8P / m)
   endif
   do r = 0, nrep
      c0 = tic()
      call swept_cells_3d(s, p0, p1, x, x, x, sw)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), sw%alpha0) ; call hash_add(h(r), sw%alpha1) ; call hash_add(h(r), sw%swept)
      call hash_add(h(r), sw%face_x) ; call hash_add(h(r), sw%face_y) ; call hash_add(h(r), sw%face_z)
      call hash_add(h(r), sw%wall) ; call hash_add(h(r), sw%residual) ; call hash_add(h(r), sw%exact)
      call hash_add(h(r), sw%evaluations)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, s%get_facets_number(), ' facets, ', m, '^3 cells, 1 step'
   if (rotation) then
      call report('swept', 'swept_cells_3d.rotation', sz, 1_I8P, t, h)
   else
      call report('swept', 'swept_cells_3d.translation', sz, 1_I8P, t, h)
   endif
   endsubroutine k_swept_3d

   ! -------------------------------------------------------------------------------------------------------------- CCD
   subroutine k_toi_3d(a, b, nsteps, nrep)
   !< time_of_impact (default, adaptive) of body a moving across body b over nsteps fixed-seed steps (fossil_val_ccd 1).
   type(surface_stl_object), intent(in) :: a, b   !< Moving body, fixed body.
   integer(I4P),             intent(in) :: nsteps, nrep !< Steps, repetitions.
   type(rigid_pose_t), allocatable :: q0(:), q1(:) !< Poses of a.
   type(rigid_pose_t)       :: b0            !< Pose of b.
   real(R8P), allocatable   :: toi(:), d(:)  !< Outputs.
   integer(I4P), allocatable :: st(:), it(:) !< Outputs.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0, seed      !< Clock, seed.
   real(R8P)                :: w(6)          !< Draws.
   integer(I4P)             :: r, k          !< Counters.
   character(64)            :: sz            !< Size.

   allocate(q0(nsteps), q1(nsteps), toi(nsteps), d(nsteps), st(nsteps), it(nsteps))
   b0 = pose_identity() ; seed = 20261107_I8P
   do k = 1, nsteps
      call prng(seed, w)
      q0(k) = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [-1.2_R8P, 0.4_R8P * (w(5:6) - 0.5_R8P)])
      call prng(seed, w)
      q1(k) = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [1.2_R8P, 0.4_R8P * (w(5:6) - 0.5_R8P)])
   enddo
   do r = 0, nrep
      c0 = tic()
      do k = 1, nsteps
         call b%time_of_impact(a, q0(k), q1(k), b0, b0, toi(k), st(k), tolerance=1e-6_R8P, distance=d(k), &
                               iterations=it(k))
      enddo
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), toi) ; call hash_add(h(r), d) ; call hash_add(h(r), st) ; call hash_add(h(r), it)
   enddo
   write(sz, '(A,I0,A,I0,A,I0,A)') 'bunny ', a%get_facets_number(), ' vs dragon ', b%get_facets_number(), ', ', &
         sum(it), ' steps'
   call report('ccd', 'time_of_impact', sz, int(nsteps, I8P), t, h)
   endsubroutine k_toi_3d

   subroutine k_toi_cubes(nrep)
   !< time_of_impact, plain bound, of a cube sliding past another at a gap of 16 tol: the 1000-step cap.
   integer(I4P), intent(in) :: nrep          !< Repetitions.
   type(surface_stl_object) :: a, b          !< Cubes.
   type(rigid_pose_t)       :: a0, a1, b0    !< Poses.
   real(R8P)                :: toi, d        !< Outputs.
   integer(I4P)             :: st, it        !< Outputs.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r             !< Counter.
   character(64)            :: sz            !< Size.

   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%translate(x=2._R8P, y=1._R8P + 16e-6_R8P, z=0._R8P) ; call b%analyze
   a0 = pose_identity() ; a1 = a0 ; a1%t = [4._R8P, 0._R8P, 0._R8P] ; b0 = pose_identity()
   do r = 0, nrep
      c0 = tic()
      call a%time_of_impact(b, a0, a1, b0, b0, toi, st, tolerance=1e-6_R8P, distance=d, iterations=it, &
                            directional=.false.)
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), toi) ; call hash_add(h(r), d) ; call hash_add(h(r), st) ; call hash_add(h(r), it)
   enddo
   write(sz, '(A,I0,A)') 'sliding cubes, gap 16 tol, plain, ', it, ' steps'
   call report('ccd', 'time_of_impact.sliding', sz, int(it, I8P), t, h)
   endsubroutine k_toi_cubes

   subroutine k_toi_2d(nseg, nsteps, nrep)
   !< curve time_of_impact (default) of an ellipse turning and translating across another, nsteps fixed-seed steps.
   integer(I4P), intent(in) :: nseg, nsteps, nrep !< Segments, steps, repetitions.
   type(curve_object)       :: a, b          !< Curves.
   type(rigid_pose_2d_t), allocatable :: q0(:), q1(:) !< Poses of a.
   type(rigid_pose_2d_t)    :: b0            !< Pose of b.
   real(R8P), allocatable   :: toi(:), d(:)  !< Outputs.
   integer(I4P), allocatable :: st(:), it(:) !< Outputs.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0, seed      !< Clock, seed.
   real(R8P)                :: w(6)          !< Draws.
   integer(I4P)             :: r, k          !< Counters.
   character(64)            :: sz            !< Size.

   call ellipse(a, nseg, 0.5_R8P, 0.2_R8P, 0._R8P, 0._R8P) ; call ellipse(b, nseg, 0.4_R8P, 0.3_R8P, 0._R8P, 0._R8P)
   allocate(q0(nsteps), q1(nsteps), toi(nsteps), d(nsteps), st(nsteps), it(nsteps))
   b0 = pose_2d(0._R8P) ; seed = 20261108_I8P
   do k = 1, nsteps
      call prng(seed, w(1:4))
      q0(k) = pose_2d(6.28_R8P * w(1), [-1.5_R8P, 0.4_R8P * (w(2) - 0.5_R8P)])
      q1(k) = pose_2d(6.28_R8P * w(3), [1.5_R8P, 0.4_R8P * (w(4) - 0.5_R8P)])
   enddo
   do r = 0, nrep
      c0 = tic()
      do k = 1, nsteps
         call a%time_of_impact(b, q0(k), q1(k), b0, b0, toi(k), st(k), tolerance=1e-6_R8P, distance=d(k), &
                               iterations=it(k))
      enddo
      if (r > 0) t(r) = toc(c0)
      call hash_add(h(r), toi) ; call hash_add(h(r), d) ; call hash_add(h(r), st) ; call hash_add(h(r), it)
   enddo
   write(sz, '(A,I0,A,I0,A)') 'ellipses ', nseg, ' segments, ', sum(it), ' steps'
   call report('ccd', 'curve.time_of_impact', sz, int(nsteps, I8P), t, h)
   endsubroutine k_toi_2d

   ! ------------------------------------------------------------------------------------------------------- processing
   subroutine k_remesh(file_name, label, nrep)
   !< isotropic_remesh with the default parameters (median edge length, 5 iterations), a fresh load per repetition.
   character(*), intent(in) :: file_name, label !< Input, label.
   integer(I4P), intent(in) :: nrep          !< Repetitions.
   type(surface_stl_object) :: s             !< Surface.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r, st, n0     !< Counters, status.
   character(64)            :: sz            !< Size.

   do r = 0, nrep
      call s%destroy ; call s%load_from_file(file_name=file_name, guess_format=.true.) ; n0 = s%get_facets_number()
      c0 = tic()
      call s%isotropic_remesh(status=st)
      if (r > 0) t(r) = toc(c0)
      call hash_surface(h(r), s) ; call hash_add(h(r), st)
   enddo
   write(sz, '(A,1X,I0,A,I0,A)') label, n0, ' facets -> ', s%get_facets_number(), ', defaults'
   call report('processing', 'isotropic_remesh', sz, int(n0, I8P), t, h)
   endsubroutine k_remesh

   subroutine k_decimate(file_name, label, nrep)
   !< decimate to a quarter of the facets (the docs recipe), a fresh load per repetition.
   character(*), intent(in) :: file_name, label !< Input, label.
   integer(I4P), intent(in) :: nrep          !< Repetitions.
   type(surface_stl_object) :: s             !< Surface.
   real(R8P)                :: t(nrep)       !< Times.
   type(hash_t)             :: h(0:nrep)     !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r, st, n0     !< Counters, status.
   character(64)            :: sz            !< Size.

   do r = 0, nrep
      call s%destroy ; call s%load_from_file(file_name=file_name, guess_format=.true.) ; n0 = s%get_facets_number()
      c0 = tic()
      call s%decimate(target_facets=n0 / 4, status=st)
      if (r > 0) t(r) = toc(c0)
      call hash_surface(h(r), s) ; call hash_add(h(r), st)
   enddo
   write(sz, '(A,1X,I0,A,I0)') label, n0, ' facets -> ', s%get_facets_number()
   call report('processing', 'decimate', sz, int(n0, I8P), t, h)
   endsubroutine k_decimate

   ! ---------------------------------------------------------------------------------------------------------- overset
   subroutine cylinder_assembly(m, blocks)
   !< The cylinder O-grid (32m x 12m cells, wall inside, overlap outside) in a (40m)^2 Cartesian background.
   integer(I4P),                 intent(in)    :: m         !< Refinement.
   type(curvilinear_block_2d_t), intent(inout) :: blocks(2) !< Blocks.
   type(topo_record_t)    :: rec                            !< Record.
   real(R8P), allocatable :: u(:,:), v(:,:), x(:,:), y(:,:) !< Vertices.
   integer(I4P)           :: ni, nj, nb, i, j               !< Counters.
   real(R8P)              :: rr                             !< Radius.

   ni = 32 * m + 1 ; nj = 12 * m + 1 ; nb = 40 * m + 1
   allocate(u(ni, nj), v(ni, nj), x(nb, nb), y(nb, nb))
   do j = 1, nj
      rr = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = rr * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = rr * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   call blocks(1)%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   do j = 1, nb
      do i = 1, nb
         x(i, j) = -3._R8P + 6._R8P * (i - 1) / (nb - 1) ; y(i, j) = -3._R8P + 6._R8P * (j - 1) / (nb - 1)
      enddo
   enddo
   call blocks(2)%initialize(x, y)
   endsubroutine cylinder_assembly

   subroutine sphere_assembly(blocks)
   !< The spherical-shell O-grid (64 x 32 x 16 cells, poles collapsed) in a 49^3 background (fossil_val_overset 2).
   type(curvilinear_block_t), intent(inout) :: blocks(2) !< Blocks.
   type(topo_record_t)    :: rec                         !< Record.
   real(R8P), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)           :: i, j, k, na, np, nr, nb     !< Counters.
   real(R8P)              :: rr, th, ph                  !< Coordinates.

   na = 64 ; np = 32 ; nr = 16 ; nb = 49
   allocate(x(na + 1, np + 1, nr + 1), y(na + 1, np + 1, nr + 1), z(na + 1, np + 1, nr + 1))
   do k = 1, nr + 1
      rr = 1._R8P + (k - 1._R8P) / nr
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            x(i, j, k) = rr * sin(ph) * cos(th) ; y(i, j, k) = rr * sin(ph) * sin(th) ; z(i, j, k) = rr * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            x(:, j, k) = 0._R8P ; y(:, j, k) = 0._R8P ; z(:, j, k) = merge(rr, -rr, j == 1)
         endif
      enddo
   enddo
   x(na + 1,:,:) = x(1,:,:) ; y(na + 1,:,:) = y(1,:,:) ; z(na + 1,:,:) = z(1,:,:)
   call blocks(1)%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [na + 1, np + 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 6 ; rec%lo = [1, 1, nr + 1] ; rec%hi = [na + 1, np + 1, nr + 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   deallocate(x, y, z) ; allocate(x(nb, nb, nb), y(nb, nb, nb), z(nb, nb, nb))
   do k = 1, nb
      do j = 1, nb
         do i = 1, nb
            x(i, j, k) = -2.6_R8P + 5.2_R8P * (i - 1) / (nb - 1) ; y(i, j, k) = -2.6_R8P + 5.2_R8P * (j - 1) / (nb - 1)
            z(i, j, k) = -2.6_R8P + 5.2_R8P * (k - 1) / (nb - 1)
         enddo
      enddo
   enddo
   call blocks(2)%initialize(x, y, z)
   endsubroutine sphere_assembly

   subroutine hash_overset(h, holes, donors, masks)
   !< Add the holes, the donor candidates and the classification.
   type(hash_t),          intent(inout) :: h          !< Hash.
   type(block_logical_t), intent(in)    :: holes(:)   !< Holes.
   type(donor_set_t),     intent(in)    :: donors(:)  !< Candidates.
   type(overset_mask_t),  intent(in)    :: masks(:)   !< Classes.
   integer(I4P)                         :: b          !< Counter.

   do b = 1, size(holes)
      call hash_add(h, holes(b)%v)
   enddo
   do b = 1, size(donors)
      call hash_add(h, donors(b)%start) ; call hash_add(h, donors(b)%block) ; call hash_add(h, donors(b)%cell)
      call hash_add(h, donors(b)%xi) ; call hash_add(h, donors(b)%status) ; call hash_add(h, donors(b)%ratio)
      call hash_add(h, donors(b)%min_jacobian)
   enddo
   do b = 1, size(masks)
      call hash_add(h, masks(b)%iota) ; call hash_add(h, masks(b)%donor) ; call hash_add(h, masks(b)%orphan)
   enddo
   endsubroutine hash_overset

   subroutine k_overset_2d(m, nrep)
   !< 2-D explicit assembly of the cylinder case (q = 3): search + wall cutter + holes, donor candidates, classify.
   integer(I4P), intent(in) :: m, nrep        !< Refinement, repetitions.
   type(curvilinear_block_2d_t)       :: blocks(2) !< Blocks.
   type(curve_object)                 :: wall      !< Cutter.
   type(block_logical_t), allocatable :: holes(:)  !< Holes.
   type(donor_set_t),     allocatable :: donors(:) !< Candidates.
   type(overset_mask_t),  allocatable :: masks(:)  !< Classes.
   real(R8P)                :: t(nrep, 3)    !< Times.
   type(hash_t)             :: h(0:nrep, 3)  !< Hashes.
   integer(I8P)             :: c0            !< Clock.
   integer(I4P)             :: r, k          !< Counters.
   character(64)            :: sz            !< Size.

   do r = 0, nrep
      call cylinder_assembly(m, blocks) ; call wall%destroy
      c0 = tic()
      call blocks(1)%build_search ; call blocks(2)%build_search
      call wall_cutter_2d(blocks(1:1), wall)
      call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
      if (r > 0) t(r, 1) = toc(c0)
      c0 = tic()
      call donor_candidates_2d(blocks, donors)
      if (r > 0) t(r, 2) = toc(c0)
      c0 = tic()
      call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, layers=2_I4P)
      if (r > 0) t(r, 3) = toc(c0)
      do k = 1, size(holes)
         call hash_add(h(r, 1), holes(k)%v)
      enddo
      call hash_overset(h(r, 2), holes(1:0), donors, masks(1:0))
      call hash_overset(h(r, 3), holes(1:0), donors(1:0), masks)
   enddo
   write(sz, '(A,I0,A,I0,A,I0,A)') 'cylinder O ', 32 * m, 'x', 12 * m, ' + ', 40 * m, '^2 cells'
   call report('overset', 'overset_2d.holes', sz, int((32 * m + 1) * (12 * m + 1) + (40 * m + 1) ** 2, I8P), t(:, 1), h(:, 1))
   call report('overset', 'overset_2d.donors', sz, int((32 * m + 1) * (12 * m + 1) + (40 * m + 1) ** 2, I8P), t(:, 2), &
               h(:, 2))
   call report('overset', 'overset_2d.classify', sz, int((32 * m + 1) * (12 * m + 1) + (40 * m + 1) ** 2, I8P), t(:, 3), &
               h(:, 3))
   endsubroutine k_overset_2d

   subroutine k_overset_3d(nrep)
   !< 3-D explicit assembly of the sphere case (q = 3): search + wall cutter + holes, donor candidates, classify.
   integer(I4P), intent(in) :: nrep           !< Repetitions.
   type(curvilinear_block_t)          :: blocks(2) !< Blocks.
   type(surface_stl_object)           :: wall      !< Cutter.
   type(block_logical_t), allocatable :: holes(:)  !< Holes.
   type(donor_set_t),     allocatable :: donors(:) !< Candidates.
   type(overset_mask_t),  allocatable :: masks(:)  !< Classes.
   real(R8P)                :: t(nrep, 3)    !< Times.
   type(hash_t)             :: h(0:nrep, 3)  !< Hashes.
   integer(I8P)             :: c0, nv        !< Clock, vertices.
   integer(I4P)             :: r, k          !< Counters.
   character(64)            :: sz            !< Size.

   do r = 0, nrep
      call sphere_assembly(blocks) ; call wall%destroy
      c0 = tic()
      call blocks(1)%build_search ; call blocks(2)%build_search
      call wall_cutter(blocks(1:1), wall)
      call hole_flags(blocks, wall, holes, cuts=[.false., .true.])
      if (r > 0) t(r, 1) = toc(c0)
      c0 = tic()
      call donor_candidates(blocks, donors)
      if (r > 0) t(r, 2) = toc(c0)
      c0 = tic()
      call classify(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, layers=2_I4P)
      if (r > 0) t(r, 3) = toc(c0)
      do k = 1, size(holes)
         call hash_add(h(r, 1), holes(k)%v)
      enddo
      call hash_overset(h(r, 2), holes(1:0), donors, masks(1:0))
      call hash_overset(h(r, 3), holes(1:0), donors(1:0), masks)
   enddo
   nv = 65_I8P * 33 * 17 + 49_I8P ** 3
   sz = 'sphere shell 64x32x16 + 48^3 cells'
   call report('overset', 'overset_3d.holes', sz, nv, t(:, 1), h(:, 1))
   call report('overset', 'overset_3d.donors', sz, nv, t(:, 2), h(:, 2))
   call report('overset', 'overset_3d.classify', sz, nv, t(:, 3), h(:, 3))
   endsubroutine k_overset_3d
endmodule fossil_val_perf_kit

program fossil_val_perf
!< FOSSIL validation: the performance bench (issue #94).
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use fossil,                   only : surface_stl_object
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use fossil_val_perf_kit
use penf,                     only : I4P, I8P, R8P

implicit none

integer(I4P), parameter   :: SW2_SEG = 256, SW2_M = 128, SW2_MR = 64 !< Swept 2-D: segments, cells (translation, rotation).
integer(I4P), parameter   :: SW3_M = 16, SW3_MR = 8 !< Swept 3-D: cells (translation, rotation).
type(surface_stl_object)  :: dragon_fine, sphere, ccd_a, ccd_b, cone
type(body_t), allocatable :: bb(:)
character(256)            :: arg, commit, host
character(8)              :: date
character(10)             :: time
integer(I4P)              :: nrep, i, st, u
integer(I8P)              :: c0
real(R8P)                 :: ttotal

nrep = 5 ; commit = 'unknown' ; host = ''
call get_environment_variable('FOSSIL_PERF_REPS', arg, status=st)
if (st == 0 .and. len_trim(arg) > 0) read(arg, *) nrep
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (index(arg, 'reps=') == 1) read(arg(6:), *) nrep
   if (index(arg, 'commit=') == 1) commit = arg(8:)
   if (index(arg, 'host=') == 1) host = arg(6:)
   if (index(arg, 'only=') == 1) perf_only = arg(6:)
enddo
nrep = max(nrep, 1_I4P)
if (len_trim(host) == 0) call get_environment_variable('HOSTNAME', host)
if (len_trim(host) == 0) then
   open(newunit=u, file='/etc/hostname', action='read', status='old', iostat=st)
   if (st == 0) then
      read(u, '(A)', iostat=st) host ; close(u)
   endif
endif
call date_and_time(date, time)
print '(A)', '# FOSSIL performance bench (issue #94): fossil_val_perf'
print '(A)', '# commit:      ' // trim(commit)
print '(A)', '# host:        ' // trim(host)
print '(A)', '# compiler:    ' // compiler_version()
print '(A)', '# flags:       ' // compiler_options()
print '(A)', '# date:        ' // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // ' ' // time(1:2) // ':' // time(3:4)
print '(A,I0,A)', '# repetitions: ', nrep, ' timed after one warm-up; wall clock (system_clock int64), one thread'
call header_line
c0 = tic()

if (want('queries')) then
   call dragon_fine%load_from_file(file_name='src/tests/dragon-fine.stl', guess_format=.true.)
   call k_distance_many(dragon_fine, 'dragon-fine', 32_I4P, nrep)
   call k_query_many(dragon_fine, 'dragon-fine', 32_I4P, nrep)
   call k_box_classify(dragon_fine, 'dragon-fine', 32_I4P, nrep)
   call k_curve_distance(4096_I4P, 512_I4P, nrep)
endif
if (want('body')) then
   if (dragon_fine%get_facets_number() == 0) &
      call dragon_fine%load_from_file(file_name='src/tests/dragon-fine.stl', guess_format=.true.)
   call k_body_cubes(20000_I4P, nrep)
   call k_body_large(dragon_fine, 'dragon-fine', 20_I4P, nrep)
endif
if (want('cut')) then
   call mc_body(BODY_SPHERE, 24_I4P, sphere)
   call k_cut_2d(1024_I4P, 512_I4P, nrep)
   call k_cut_3d(sphere, 'mc sphere', 48_I4P, .false., nrep)
   call k_cut_3d(sphere, 'mc sphere', 48_I4P, .true., nrep)
endif
if (want('swept')) then
   call bench_bodies(bb) ; call to_surface(bb(7)%soup, cone) ; call centre(cone)
   call k_swept_2d(SW2_SEG, SW2_M, .false., nrep)
   call k_swept_2d(SW2_SEG, SW2_MR, .true., nrep)
   call k_swept_3d(cone, 'cone', SW3_M, .false., nrep)
   call k_swept_3d(cone, 'cone', SW3_MR, .true., nrep)
endif
if (want('ccd')) then
   call ccd_a%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
   call ccd_a%decimate(target_facets=5000_I4P)
   call ccd_b%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
   call ccd_b%decimate(target_facets=2000_I4P)
   call scale_to_unit(ccd_a) ; call scale_to_unit(ccd_b)
   call k_toi_3d(ccd_a, ccd_b, 8_I4P, nrep)
   call k_toi_cubes(nrep)
   call k_toi_2d(256_I4P, 16_I4P, nrep)
endif
if (want('processing')) then
   call k_remesh('src/tests/dragon.stl', 'dragon', nrep)
   call k_decimate('src/tests/dragon.stl', 'dragon', nrep)
endif
if (want('overset')) then
   call k_overset_2d(16_I4P, nrep)
   call k_overset_3d(nrep)
endif

ttotal = toc(c0)
print '(A,F8.1,A)', '# total: ', ttotal, ' s'
print '(A,L1)', 'Are all tests passed? ', perf_ok
if (.not. perf_ok) error stop 1
endprogram fossil_val_perf
