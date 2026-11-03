!< FOSSIL test: box classification on the trees (issue #38 step 2).
program fossil_test_box_classify
!< FOSSIL test: box classification on the trees (issue #38 steps 2 and 3).
!<
!< 3-D (`surface%box_classify_many`, `box_intersects_surface`):
!< - cube-02 and bunny, random boxes (sizes 1-30% and 1-5% of the diagonal) with a near radius: the BVH codes
!<   equal a brute-force reference over all facets (same predicate, cut/near/centre sign); all four
!<   codes occur; `inside` equals `is_point_inside` at the centre for every box not cut;
!< - cube-02: every box not cut has its 8 corners on the side of its centre (closed surface);
!< - cube-02 on the octree (brute-force fallback) gives the BVH codes;
!< - touching: a box sharing the face plane of the cube is cut; moved off by 1e-6 it is not, it is
!<   near with radius 1e-5 and outside with radius 0;
!< - refusals: negative radius, size mismatch.
!< 2-D (`curve%box_classify_many`, `box_intersects_curve`):
!< - NACA 0012 and an annulus with an island, random boxes: codes equal a brute-force reference
!<   over all segments (exact cut, distance for near), all codes occur, corner consistency;
!< - touching: a box sharing the annulus outer side is cut, one ulp off it is not.
!< Device (`box_classify_many_device`, step 3; on the host in tests-gnu): codes and `inside` equal the
!< CPU ones on cube-02 and bunny (pseudo-normal and fast-winding-number signs) and on the NACA 0012
!< and the annulus (both signs); use before `enter_device` is refused. With no `sign_algorithm` the codes
!< and `inside` equal the explicit fast-winding-number ones, CPU and device (the default, issue #39).

use fossil,                only : surface_stl_object, curve_object, facet_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, &
                                  BOX_OUTSIDE, AABB_TREE_OCTREE, STATUS_INVALID_INPUT, STATUS_OK, &
                                  SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use fossil_box_predicates, only : box_triangle_overlap, box_segment_overlap_2d, box_segment_distance2_2d
use fossil_box_query,      only : box_tolerance
use penf,                  only : I4P, R8P
use vecfor,                only : vector_R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube, bunny, cube_oct
type(curve_object)       :: c
real(R8P)                :: u(200), v(200), b0(3), b1(3)
integer(I4P)             :: k, st, st2, code(1)
logical                  :: ok(10)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call cube_oct%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)

ok(1) = check3(cube, 2000, 0.3_R8P, .true., 'cube-02')
ok(2) = check3(bunny, 400, 0.05_R8P, .false., 'bunny')
ok(3) = same_codes(cube, cube_oct, 1000)
print '(A,L1)', 'cube-02 octree (brute force) = BVH codes: ', ok(3)

! touching / separated from the face x = bmax of the cube
call bbox(cube, b0, b1)
ok(4) = cube%box_intersects_surface([b1(1), 0.5_R8P * (b0(2) + b1(2)), 0.5_R8P * (b0(3) + b1(3))], &
                                    [b1(1) + 0.1_R8P, 0.6_R8P * b1(2) + 0.4_R8P * b0(2), 0.6_R8P * b1(3) + 0.4_R8P * b0(3)])
ok(4) = ok(4) .and. .not. cube%box_intersects_surface([b1(1) + 1.e-6_R8P, b0(2), b0(3)], [b1(1) + 0.1_R8P, b1(2), b1(3)])
call cube%box_classify_many(lo=reshape([b1(1) + 1.e-6_R8P, b0(2), b0(3)], [3, 1]), &
                            hi=reshape([b1(1) + 0.1_R8P, b1(2), b1(3)], [3, 1]), codes=code, radius=1.e-5_R8P)
ok(4) = ok(4) .and. code(1) == BOX_NEAR
call cube%box_classify_many(lo=reshape([b1(1) + 1.e-6_R8P, b0(2), b0(3)], [3, 1]), &
                            hi=reshape([b1(1) + 0.1_R8P, b1(2), b1(3)], [3, 1]), codes=code)
ok(4) = ok(4) .and. code(1) == BOX_OUTSIDE
call cube%box_classify_many(lo=reshape(b0, [3, 1]), hi=reshape(b1, [3, 1]), codes=code, radius=-1._R8P, status=st)
call cube%box_classify_many(lo=reshape([b0, b0], [3, 2]), hi=reshape(b1, [3, 1]), codes=code, status=st2)
ok(4) = ok(4) .and. st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT
print '(A,L1)', '3-D touching, separated, refusals: ', ok(4)

! 2-D
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(5) = check2(c, [-0.3_R8P, -0.3_R8P], [1.3_R8P, 0.3_R8P], 'NACA 0012')
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P,  2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P,  2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
ok(6) = check2(c, [-1._R8P, -1._R8P], [7._R8P, 7._R8P], 'annulus with an island')
ok(7) = c%box_intersects_curve([6._R8P, 1._R8P], [7._R8P, 2._R8P]) .and. &
        .not. c%box_intersects_curve([6._R8P + spacing(6._R8P), 1._R8P], [7._R8P, 2._R8P]) .and. &
        c%box_intersects_curve([2.5_R8P, 2.5_R8P], [3.5_R8P, 4._R8P]) .and. &
        .not. c%box_intersects_curve([2.5_R8P, 2.5_R8P], [3.5_R8P, 4._R8P - spacing(4._R8P)])
print '(A,L1)', '2-D touching and one-ulp separated: ', ok(7)
call c%box_classify_many(lo=reshape([0._R8P, 0._R8P], [2, 1]), hi=reshape([1._R8P, 1._R8P], [2, 1]), codes=code, &
                         radius=-1._R8P, status=st)
ok(8) = st == STATUS_INVALID_INPUT
print '(A,L1)', '2-D refusal: ', ok(8)

! device paths (step 3)
call cube_oct%box_classify_many_device(lo=reshape(b0, [3, 1]), hi=reshape(b1, [3, 1]), codes=code, status=st)
call cube%enter_device
call bunny%enter_device
ok(9) = st == STATUS_INVALID_INPUT .and. dev3(cube, 0.3_R8P) .and. dev3(bunny, 0.05_R8P)
print '(A,L1)', '3-D device = CPU codes (pseudo-normal, FWN), refusal before enter_device: ', ok(9)
call c%box_classify_many_device(lo=reshape([0._R8P, 0._R8P], [2, 1]), hi=reshape([1._R8P, 1._R8P], [2, 1]), codes=code, &
                                status=st)
ok(10) = st == STATUS_INVALID_INPUT .and. dev2(c, [-1._R8P, -1._R8P], [7._R8P, 7._R8P])
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(10) = ok(10) .and. dev2(c, [-0.3_R8P, -0.3_R8P], [1.3_R8P, 0.3_R8P])
print '(A,L1)', '2-D device = CPU codes (pseudo-normal, FWN), refusal before enter_device: ', ok(10)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check3(s, n, smax, corners, label) result(good)
   !< BVH classification vs the brute-force reference on n random boxes.
   type(surface_stl_object), intent(in) :: s       !< Surface.
   integer,                  intent(in) :: n       !< Boxes.
   real(R8P),                intent(in) :: smax    !< Largest box size (fraction of the diagonal).
   logical,                  intent(in) :: corners !< Also check the corner consistency (closed surface).
   character(*),             intent(in) :: label   !< Label.
   logical                              :: good    !< Result.
   real(R8P)        :: lo(3, n), hi(3, n), r(6), bmin(3), bmax(3), diag, rad, tol !< Work.
   integer(I4P)     :: codes(n), ref, m, mism, cnt(4), cbad, ibad, j, nseed !< Work.
   logical          :: ins(n), cin                 !< Work.

   call bbox(s, bmin, bmax)
   diag = norm2(bmax - bmin) ; rad = 0.1_R8P * smax * diag ; tol = box_tolerance(bmin, bmax)
   call random_seed(size=nseed)
   call random_seed(put=[(31 + m, m = 1, nseed)])
   do m = 1, n
      call random_number(r)
      lo(:, m) = bmin - 0.2_R8P * (bmax - bmin) + 1.4_R8P * r(1:3) * (bmax - bmin)
      hi(:, m) = lo(:, m) + (0.01_R8P + (smax - 0.01_R8P) * r(4:6)) * diag
   enddo
   call s%box_classify_many(lo=lo, hi=hi, codes=codes, radius=rad, inside=ins)
   mism = 0 ; cnt = 0 ; cbad = 0 ; ibad = 0
   do m = 1, n
      ref = brute3(s, lo(:, m), hi(:, m), tol, rad + tol)
      cin = s%is_point_inside(vector_R8P(0.5_R8P * (lo(1, m) + hi(1, m)), 0.5_R8P * (lo(2, m) + hi(2, m)), &
                                         0.5_R8P * (lo(3, m) + hi(3, m))))
      if (ref == 0) then
         ref = BOX_OUTSIDE ; if (cin) ref = BOX_INSIDE
      endif
      if (ref /= codes(m)) mism = mism + 1
      if (codes(m) /= BOX_CUT .and. (ins(m) .neqv. cin)) ibad = ibad + 1
      if (s%box_intersects_surface(lo(:, m), hi(:, m)) .neqv. (codes(m) == BOX_CUT)) mism = mism + 1
      cnt(codes(m)) = cnt(codes(m)) + 1
      if (corners .and. codes(m) /= BOX_CUT) then
         do j = 0, 7
            if (s%is_point_inside(vector_R8P(merge(hi(1, m), lo(1, m), btest(j, 0)), merge(hi(2, m), lo(2, m), btest(j, 1)), &
                                             merge(hi(3, m), lo(3, m), btest(j, 2)))) .neqv. cin) cbad = cbad + 1
         enddo
      endif
   enddo
   good = mism == 0 .and. ibad == 0 .and. cbad == 0 .and. all(cnt > 0)
   print '(A,I0,A,4(I0,1X),A,I0,A,I0,A,I0)', label//': ', n, ' boxes, cut/near/in/out ', cnt, &
         ' mismatches vs brute ', mism, ' inside errors ', ibad, ' corner errors ', cbad
   endfunction check3

   function brute3(s, lo, hi, tol, reach) result(code)
   !< 1 cut, 2 near, 0 none, over all facets.
   type(surface_stl_object), intent(in) :: s            !< Surface.
   real(R8P),                intent(in) :: lo(3), hi(3) !< Box.
   real(R8P),                intent(in) :: tol, reach   !< Inflations.
   integer(I4P)                         :: code         !< Result.
   type(facet_object), pointer          :: f            !< Facet.
   real(R8P)                            :: v(3, 3)      !< Vertices.
   integer(I4P)                         :: i, j         !< Counters.

   code = 0
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      do j = 1, 3
         v(:, j) = [f%vertex(j)%x, f%vertex(j)%y, f%vertex(j)%z]
      enddo
      if (box_triangle_overlap(lo, hi, v(:, 1), v(:, 2), v(:, 3), tol)) then
         code = BOX_CUT
         return
      endif
      if (box_triangle_overlap(lo, hi, v(:, 1), v(:, 2), v(:, 3), reach)) code = BOX_NEAR
   enddo
   endfunction brute3

   function same_codes(s1, s2, n) result(same)
   !< Codes of two surfaces equal on n random boxes.
   type(surface_stl_object), intent(in) :: s1, s2 !< Surfaces.
   integer,                  intent(in) :: n      !< Boxes.
   logical                              :: same   !< Result.
   real(R8P)    :: lo(3, n), hi(3, n), r(6), bmin(3), bmax(3) !< Work.
   integer(I4P) :: c1(n), c2(n), m !< Work.

   call bbox(s1, bmin, bmax)
   do m = 1, n
      call random_number(r)
      lo(:, m) = bmin - 0.2_R8P * (bmax - bmin) + 1.4_R8P * r(1:3) * (bmax - bmin)
      hi(:, m) = lo(:, m) + (0.01_R8P + 0.29_R8P * r(4:6)) * norm2(bmax - bmin)
   enddo
   call s1%box_classify_many(lo=lo, hi=hi, codes=c1, radius=0.05_R8P)
   call s2%box_classify_many(lo=lo, hi=hi, codes=c2, radius=0.05_R8P)
   same = all(c1 == c2)
   endfunction same_codes

   function check2(cur, blo, bhi, label) result(good)
   !< Curve classification vs the brute-force reference on random boxes.
   type(curve_object), intent(in) :: cur        !< Curve.
   real(R8P),          intent(in) :: blo(2), bhi(2) !< Sampling region.
   character(*),       intent(in) :: label      !< Label.
   logical                        :: good       !< Result.
   integer, parameter :: N = 3000 !< Boxes.
   real(R8P)    :: lo(2, N), hi(2, N), r(4), rad, d2, cc(2), a(2), b(2) !< Work.
   integer(I4P) :: codes(N), ref, m, mism, cnt(4), cbad, ibad, j, l, i1, i2, i, nseed, st !< Work.
   logical      :: ins(N), cin !< Work.

   rad = 0.02_R8P * norm2(bhi - blo)
   call random_seed(size=nseed)
   call random_seed(put=[(47 + m, m = 1, nseed)])
   do m = 1, N
      call random_number(r)
      lo(:, m) = blo + r(1:2) * (bhi - blo)
      hi(:, m) = lo(:, m) + (0.005_R8P + 0.15_R8P * r(3:4)) * (bhi - blo)
   enddo
   call cur%box_classify_many(lo=lo, hi=hi, codes=codes, radius=rad, inside=ins, status=st)
   mism = 0 ; cnt = 0 ; cbad = 0 ; ibad = 0
   if (st /= STATUS_OK) mism = mism + 1
   do m = 1, N
      ref = 0 ; d2 = huge(1._R8P)
      loops: do l = 1, cur%get_loops_number()
         call cur%get_loop_range(l, i1, i2)
         do i = i1, i2
            if (i == i2 .and. .not. cur%is_loop_closed(l)) cycle
            a = cur%get_vertex(i) ; b = cur%get_vertex(merge(i1, i + 1, i == i2))
            if (box_segment_overlap_2d(lo(:, m), hi(:, m), a, b)) then
               ref = BOX_CUT
               exit loops
            endif
            d2 = min(d2, box_segment_distance2_2d(lo(:, m), hi(:, m), a, b))
         enddo
      enddo loops
      if (ref == 0 .and. d2 <= rad * rad) ref = BOX_NEAR
      cc = 0.5_R8P * (lo(:, m) + hi(:, m))
      cin = cur%is_point_inside(cc)
      if (ref == 0) then
         ref = BOX_OUTSIDE ; if (cin) ref = BOX_INSIDE
      endif
      if (ref /= codes(m)) mism = mism + 1
      if (codes(m) /= BOX_CUT .and. (ins(m) .neqv. cin)) ibad = ibad + 1
      if (cur%box_intersects_curve(lo(:, m), hi(:, m)) .neqv. (codes(m) == BOX_CUT)) mism = mism + 1
      cnt(codes(m)) = cnt(codes(m)) + 1
      if (codes(m) /= BOX_CUT) then
         do j = 0, 3
            if (cur%is_point_inside([merge(hi(1, m), lo(1, m), btest(j, 0)), merge(hi(2, m), lo(2, m), btest(j, 1))]) &
                .neqv. cin) cbad = cbad + 1
         enddo
      endif
   enddo
   good = mism == 0 .and. ibad == 0 .and. cbad == 0 .and. all(cnt > 0)
   print '(A,I0,A,4(I0,1X),A,I0,A,I0,A,I0)', label//': ', N, ' boxes, cut/near/in/out ', cnt, &
         ' mismatches vs brute ', mism, ' inside errors ', ibad, ' corner errors ', cbad
   endfunction check2

   function dev3(s, smax) result(same)
   !< Device codes and inside flags equal the CPU ones, both device signs.
   type(surface_stl_object), intent(inout) :: s    !< Surface (on the device).
   real(R8P),                intent(in)    :: smax !< Largest box size (fraction of the diagonal).
   logical                                 :: same !< Result.
   integer, parameter :: N = 1500 !< Boxes.
   integer(I4P), parameter :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   real(R8P)    :: lo(3, N), hi(3, N), r(6), bmin(3), bmax(3), diag !< Work.
   integer(I4P) :: c1(N), c2(N), m, a, st !< Work.
   logical      :: i1(N), i2(N) !< Work.

   call bbox(s, bmin, bmax) ; diag = norm2(bmax - bmin)
   do m = 1, N
      call random_number(r)
      lo(:, m) = bmin - 0.2_R8P * (bmax - bmin) + 1.4_R8P * r(1:3) * (bmax - bmin)
      hi(:, m) = lo(:, m) + (0.01_R8P + (smax - 0.01_R8P) * r(4:6)) * diag
   enddo
   same = .true.
   do a = 1, 2
      call s%box_classify_many(lo=lo, hi=hi, codes=c1, radius=0.1_R8P * smax * diag, sign_algorithm=ALGO(a), inside=i1)
      call s%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=0.1_R8P * smax * diag, sign_algorithm=ALGO(a), &
                                      inside=i2, status=st)
      same = same .and. st == STATUS_OK .and. all(c1 == c2) .and. all(i1 .eqv. i2)
      print '(A,I0,A,I0,A,I0)', '   3-D sign ', ALGO(a), ': box device vs serial mismatches ', count(c1 /= c2), &
            '  inside mismatches ', count(i1 .neqv. i2)
   enddo
   ! default sign (issue #39): the fast winding number, CPU and device (c1, i1: explicit FWN on the CPU)
   call s%box_classify_many(lo=lo, hi=hi, codes=c2, radius=0.1_R8P * smax * diag, inside=i2)
   same = same .and. all(c2 == c1) .and. all(i2 .eqv. i1)
   call s%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=0.1_R8P * smax * diag, inside=i2, status=st)
   same = same .and. st == STATUS_OK .and. all(c2 == c1) .and. all(i2 .eqv. i1)
   print '(A,L1)', '   3-D default sign = fast winding number (CPU, device): ', all(c2 == c1) .and. all(i2 .eqv. i1)
   endfunction dev3

   function dev2(cur, blo, bhi) result(same)
   !< Device codes and inside flags equal the CPU ones, both signs.
   type(curve_object), intent(inout) :: cur            !< Curve.
   real(R8P),          intent(in)    :: blo(2), bhi(2) !< Sampling region.
   logical                           :: same           !< Result.
   integer, parameter :: N = 2000 !< Boxes.
   integer(I4P), parameter :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   real(R8P)    :: lo(2, N), hi(2, N), r(4), rad !< Work.
   integer(I4P) :: c1(N), c2(N), m, a, st !< Work.
   logical      :: i1(N), i2(N) !< Work.

   rad = 0.02_R8P * norm2(bhi - blo)
   do m = 1, N
      call random_number(r)
      lo(:, m) = blo + r(1:2) * (bhi - blo)
      hi(:, m) = lo(:, m) + (0.005_R8P + 0.15_R8P * r(3:4)) * (bhi - blo)
   enddo
   call cur%enter_device
   same = .true.
   do a = 1, 2
      call cur%box_classify_many(lo=lo, hi=hi, codes=c1, radius=rad, sign_algorithm=ALGO(a), inside=i1)
      call cur%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=rad, sign_algorithm=ALGO(a), inside=i2, status=st)
      same = same .and. st == STATUS_OK .and. all(c1 == c2) .and. all(i1 .eqv. i2)
      print '(A,I0,A,I0,A,I0)', '   2-D sign ', ALGO(a), ': box device vs serial mismatches ', count(c1 /= c2), &
            '  inside mismatches ', count(i1 .neqv. i2)
   enddo
   ! default sign (issue #39): the fast winding number, CPU and device (c1, i1: explicit FWN on the CPU)
   call cur%box_classify_many(lo=lo, hi=hi, codes=c2, radius=rad, inside=i2)
   same = same .and. all(c2 == c1) .and. all(i2 .eqv. i1)
   call cur%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=rad, inside=i2, status=st)
   same = same .and. st == STATUS_OK .and. all(c2 == c1) .and. all(i2 .eqv. i1)
   print '(A,L1)', '   2-D default sign = fast winding number (CPU, device): ', all(c2 == c1) .and. all(i2 .eqv. i1)
   call cur%exit_device
   endfunction dev2

   subroutine bbox(s, bmin, bmax)
   !< Bounding box of a surface as arrays.
   type(surface_stl_object), intent(in)  :: s                !< Surface.
   real(R8P),                intent(out) :: bmin(3), bmax(3) !< Box.
   type(vector_R8P)                      :: p                !< Work.

   p = s%get_bmin() ; bmin = [p%x, p%y, p%z]
   p = s%get_bmax() ; bmax = [p%x, p%y, p%z]
   endsubroutine bbox

   pure function thickness(x) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: x !< Chord position.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.126_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_box_classify
