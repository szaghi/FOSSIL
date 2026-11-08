!< FOSSIL test: wall cutters of curvilinear blocks (issue #72, step 2).
program fossil_test_overset_cutters
!< FOSSIL test: wall cutters of curvilinear blocks (issue #72, step 2).
!<
!< 2-D winding numbers by the exact sum.
!< 1. a 2-D O-grid around the unit cylinder, wall on face 3: one closed loop, solid on the left (winding 1 at the
!<    centre, 0 outside; signed distance negative inside);
!< 2. two 2-D blocks whose wall arcs overlap by 0.4 pi (fine and coarse): shared(1, 2); trimmed, the winding number at
!<    the centre is 1 within 0.02 (a residual strip one coarse element wide); untrimmed (tolerance 0) it counts the
!<    overlap twice (1.4);
!< 3. a Yin-Yang pair of spherical shells (Yang coarser): shared; trimmed, winding 1 within 0.06 at the centre (1.040
!<    since issue #91: the straddling Yang triangles stay, double cover; 1.014 with the centroid rule) and 0 within 0.02
!<    outside; untrimmed above 1.05; the gap between bilinear faces and triangles is reported (> 0, small);
!< 4. a background block without walls shares nothing; a block set without walls gives status 1;
!< 5. no gap (issue #91), 2-D: the annuli of overset.md's cutter figure (fine 0-200 deg, 24 cells; coarse 153-373 deg,
!<    11 cells): every angle of the coarse patch's range lies in the range of a kept segment of one patch or the other
!<    (exact: the midpoints between the sorted segment ends are tested); the remaining double cover measured (20 deg:
!<    193-200 and 0-13); strict=.false. (the centroid rule) leaves the gap 353-360 (7 deg);
!< 6. no gap, 3-D: the Yin-Yang pair; the Yang wall (coarse) sampled at 10 x 10 points per face, the ray from the
!<    centre through each sample hits a kept triangle of Yin or of Yang; double cover measured (both hit: 0.53 of the
!<    unit sphere's 4 pi); strict=.false. leaves samples uncovered (0.027).

use fossil_overset,           only : wall_cutter, wall_cutter_2d
use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY
use fossil_surface_stl_object, only : surface_stl_object
use fossil_curve_object,      only : curve_object
use fossil_facet_object,      only : facet_object
use fossil,                   only : SIGN_FAST_WINDING_NUMBER
use vecfor,                   only : vector_R8P
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(surface_stl_object)                  :: s
logical,      allocatable                 :: sh(:,:)
integer(I4P), allocatable                 :: own(:)
real(R8P)                                 :: w, gap, d(1)
integer(I4P)                              :: st
logical                                   :: all_ok
real(R8P)                                 :: unc, dbl

all_ok = .true.

! 1. O-grid
allocate(q(1))
call arc_block(q(1), 0._R8P, 2 * PI, 65_I4P)
call wall_cutter_2d(q, c, owner=own)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call c%distance_many(reshape([0.5_R8P, 0._R8P], [2, 1]), d, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
call check(abs(w - 1) < 1e-12_R8P .and. abs(c%winding_number([2._R8P, 0._R8P], use_brute_force=.true.)) < 1e-12_R8P &
           .and. d(1) < 0, &
           'O-grid wall: closed loop, solid on the left', w)

! 2. overlapping arcs
deallocate(q) ; allocate(q(2))
call arc_block(q(1), -0.2_R8P * PI, 1.4_R8P * PI, 49_I4P)
call arc_block(q(2), 0.8_R8P * PI, 1.4_R8P * PI, 23_I4P)
call wall_cutter_2d(q, c, shared=sh)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call check(sh(1, 2) .and. abs(w - 1) < 0.02_R8P, 'overlapping arcs: shared wall, trimmed winding 1', w)
call wall_cutter_2d(q, c, tolerance=0._R8P)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call check(w > 1.15_R8P, 'overlapping arcs untrimmed: overlap counted twice', w)

! 3. Yin-Yang
allocate(b(3))
call yin_block(b(1), 25_I4P, 9_I4P, .false.)
call yin_block(b(2), 17_I4P, 7_I4P, .true.)
call box_block(b(3))
call wall_cutter(b, s, gap=gap, shared=sh, status=st)
w = s%winding_number(vector_R8P(0._R8P, 0._R8P, 0._R8P))
call check(st == 0 .and. sh(1, 2) .and. .not. sh(1, 3) .and. .not. sh(2, 3) .and. abs(w - 1) < 0.06_R8P .and. &
           abs(s%winding_number(vector_R8P(3._R8P, 0.1_R8P, 0.2_R8P))) < 0.02_R8P, &
           'Yin-Yang: shared wall, trimmed winding 1 inside, 0 outside; box shares nothing', w)
call check(gap > 0._R8P .and. gap < 0.05_R8P, 'Yin-Yang: bilinear-to-triangle gap reported', gap)
call wall_cutter(b, s, tolerance=0._R8P)
w = s%winding_number(vector_R8P(0._R8P, 0._R8P, 0._R8P))
call check(w > 1.05_R8P, 'Yin-Yang untrimmed: overlap counted twice', w)

! 4. no walls
call wall_cutter(b(3:3), s, status=st)
call check(st == 1, 'no wall faces: status 1', 0._R8P)

! 5. no gap, 2-D
deallocate(q) ; allocate(q(2))
call arc_block(q(1), 0._R8P, 200 * PI / 180, 25_I4P)
call arc_block(q(2), 153 * PI / 180, 220 * PI / 180, 12_I4P)
call wall_cutter_2d(q, c)
call angle_cover(c, 153._R8P, 373._R8P, unc, dbl)
call check(unc == 0._R8P, 'cutter figure annuli: no gap in 153-373 deg (uncovered deg)', unc)
call check(abs(dbl - 20._R8P) < 1e-9_R8P, 'cutter figure annuli: double cover 20 deg (193-200, 0-13)', dbl)
call wall_cutter_2d(q, c, strict=.false.)
call angle_cover(c, 153._R8P, 373._R8P, unc, dbl)
call check(abs(unc - 7._R8P) < 1e-9_R8P, 'cutter figure annuli, centroid rule: gap 353-360 deg', unc)

! 6. no gap, 3-D
call wall_cutter(b, s, owner=own)
call ray_cover(s, own, b(2), 10_I4P, unc, dbl)
call check(unc == 0._R8P, 'Yin-Yang: no gap over the Yang wall (uncovered area)', unc)
print '(A,ES11.3)', '  Yin-Yang double-covered area (unit sphere, 4 pi total): ', dbl
call wall_cutter(b, s, owner=own, strict=.false.)
call ray_cover(s, own, b(2), 10_I4P, unc, dbl)
call check(unc > 0._R8P, 'Yin-Yang, centroid rule: gaps over the Yang wall (uncovered area)', unc)

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES11.3)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   subroutine angle_cover(cut, t0, t1, uncovered, double)
   !< Over the angles [t0, t1] (degrees, t1 - t0 < 360): uncovered, the measure of the angles in the range of no
   !< segment of the cutter; double, the measure of the angles in the range of two or more. Exact up to the round-off
   !< of atan2: both are constant between consecutive segment ends.
   type(curve_object), intent(in)  :: cut               !< Cutter.
   real(R8P),          intent(in)  :: t0, t1            !< Range.
   real(R8P),          intent(out) :: uncovered, double !< Measures.
   real(R8P), allocatable          :: lo(:), hi(:), e(:) !< Segment ranges, sorted ends.
   real(R8P)                       :: p(2), r(2), x, y, tmp, m !< Work.
   integer(I4P)                    :: j, n, i, k, hits  !< Counters.

   n = 0_I4P ; allocate(lo(cut%get_vertices_number()), hi(cut%get_vertices_number()))
   do j = 1, cut%get_vertices_number()
      if (cut%get_next(j) == 0) cycle
      p = cut%get_vertex(j) ; r = cut%get_vertex(cut%get_next(j))
      x = shift(atan2(p(2), p(1)) * 180 / PI, t0) ; y = shift(atan2(r(2), r(1)) * 180 / PI, t0)
      if (x > y) then
         tmp = x ; x = y ; y = tmp
      endif
      if (y - x > 180._R8P) cycle ! across the cut angle t0 - 63: outside the range here
      n = n + 1_I4P ; lo(n) = x ; hi(n) = y
   enddo
   e = [t0, t1, pack(lo(1:n), lo(1:n) > t0 .and. lo(1:n) < t1), pack(hi(1:n), hi(1:n) > t0 .and. hi(1:n) < t1)]
   do i = 2, size(e) ! insertion sort
      tmp = e(i) ; k = i - 1
      do while (k >= 1)
         if (e(k) <= tmp) exit
         e(k + 1) = e(k) ; k = k - 1
      enddo
      e(k + 1) = tmp
   enddo
   uncovered = 0._R8P ; double = 0._R8P
   do i = 1, size(e) - 1
      if (e(i + 1) <= e(i)) cycle
      m = 0.5_R8P * (e(i) + e(i + 1))
      hits = count(lo(1:n) <= m .and. m <= hi(1:n))
      if (hits == 0) uncovered = uncovered + (e(i + 1) - e(i))
      if (hits >= 2) double = double + (e(i + 1) - e(i))
   enddo
   endsubroutine angle_cover

   pure real(R8P) function shift(t, t0) result(r)
   !< The angle t (degrees) moved by whole turns into [t0 - 63, t0 + 297).
   real(R8P), intent(in) :: t, t0 !< Angles.

   r = modulo(t - (t0 - 63._R8P), 360._R8P) + (t0 - 63._R8P)
   endfunction shift

   subroutine ray_cover(cut, own, blk, m, uncovered, double)
   !< The wall (face k = 1) of blk sampled at m x m bilinear points per face; the ray from the origin through each
   !< sample is tested against the kept triangles of each owner. uncovered: the area (on the bilinear faces) of the
   !< samples hit by no triangle; double: of those hit by triangles of two owners.
   type(surface_stl_object),  intent(in)  :: cut               !< Cutter.
   integer(I4P),              intent(in)  :: own(:)            !< Owner of each facet.
   type(curvilinear_block_t), intent(in)  :: blk               !< Block of the sampled wall.
   integer(I4P),              intent(in)  :: m                 !< Samples per face direction.
   real(R8P),                 intent(out) :: uncovered, double !< Areas.
   real(R8P)                              :: c(3, 4), p(3), w, sa, ta, t(3, 3) !< Work.
   integer(I4P)                           :: i, j, ii, jj, f, k, owners !< Counters.
   logical                                :: hit(maxval(own))  !< Hit per owner.
   type(facet_object), pointer            :: fp                !< Facet.

   uncovered = 0._R8P ; double = 0._R8P
   do j = 1, blk%nj - 1
      do i = 1, blk%ni - 1
         c(:, 1) = [blk%x(i, j, 1), blk%y(i, j, 1), blk%z(i, j, 1)]
         c(:, 2) = [blk%x(i + 1, j, 1), blk%y(i + 1, j, 1), blk%z(i + 1, j, 1)]
         c(:, 3) = [blk%x(i + 1, j + 1, 1), blk%y(i + 1, j + 1, 1), blk%z(i + 1, j + 1, 1)]
         c(:, 4) = [blk%x(i, j + 1, 1), blk%y(i, j + 1, 1), blk%z(i, j + 1, 1)]
         w = 0.5_R8P * (norm2(cross3(c(:, 2) - c(:, 1), c(:, 3) - c(:, 1))) + &
                        norm2(cross3(c(:, 3) - c(:, 1), c(:, 4) - c(:, 1)))) / m**2
         do jj = 1, m
            do ii = 1, m
               sa = (ii - 0.5_R8P) / m ; ta = (jj - 0.5_R8P) / m
               p = (1 - sa) * (1 - ta) * c(:, 1) + sa * (1 - ta) * c(:, 2) + sa * ta * c(:, 3) + (1 - sa) * ta * c(:, 4)
               hit = .false.
               do f = 1, cut%get_facets_number()
                  if (hit(own(f))) cycle
                  fp => cut%facet_at(f)
                  do k = 1, 3
                     t(:, k) = [fp%vertex(k)%x, fp%vertex(k)%y, fp%vertex(k)%z]
                  enddo
                  hit(own(f)) = ray_hits(p, t)
               enddo
               owners = count(hit)
               if (owners == 0) uncovered = uncovered + w
               if (owners >= 2) double = double + w
            enddo
         enddo
      enddo
   enddo
   endsubroutine ray_cover

   pure logical function ray_hits(d, t) result(h)
   !< The ray from the origin along d meets the triangle t (Moller-Trumbore, closed with a 1e-12 slack).
   real(R8P), intent(in) :: d(3), t(3, 3) !< Direction, vertices.
   real(R8P)             :: e1(3), e2(3), pv(3), qv(3), det, u, v, s !< Work.

   e1 = t(:, 2) - t(:, 1) ; e2 = t(:, 3) - t(:, 1) ; pv = cross3(d, e2) ; det = dot_product(e1, pv) ; h = .false.
   if (abs(det) < tiny(1._R8P)) return
   u = dot_product(-t(:, 1), pv) / det ; qv = cross3(-t(:, 1), e1) ; v = dot_product(d, qv) / det
   s = dot_product(e2, qv) / det
   h = u >= -1e-12_R8P .and. v >= -1e-12_R8P .and. u + v <= 1 + 1e-12_R8P .and. s > 0
   endfunction ray_hits

   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   subroutine arc_block(blk, t0, span, ni)
   !< 2-D block over the arc [t0, t0 + span] of the annulus 1 <= r <= 1.5, wall on face 3; a full circle copies the seam.
   type(curvilinear_block_2d_t), intent(inout) :: blk    !< Block.
   real(R8P),                    intent(in)    :: t0, span !< Arc.
   integer(I4P),                 intent(in)    :: ni     !< Vertices along the arc.
   real(R8P)                                   :: u(ni, 5), v(ni, 5), t, r !< Vertices.
   integer(I4P)                                :: i, j   !< Counters.
   type(topo_record_t)                         :: rec    !< Wall record.

   do j = 1, 5
      r = 1._R8P + 0.125_R8P * (j - 1)
      do i = 1, ni
         t = t0 + span * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   if (span >= 2 * PI) then
      u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   endif
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine arc_block

   subroutine yin_block(blk, ni, nj, yang)
   !< Spherical shell 1 <= r <= 2 over |lon| <= 3 pi / 4 + 0.15, |lat| <= pi / 4 + 0.15; Yang = (-x, z, y) of Yin.
   type(curvilinear_block_t), intent(inout) :: blk    !< Block.
   integer(I4P),              intent(in)    :: ni, nj !< Vertices in longitude, latitude.
   logical,                   intent(in)    :: yang   !< Yang copy.
   real(R8P)                                :: x(ni, nj, 3), y(ni, nj, 3), z(ni, nj, 3), lo, la, r, p(3) !< Data.
   integer(I4P)                             :: i, j, k !< Counters.
   type(topo_record_t)                      :: rec    !< Wall record.

   do k = 1, 3
      r = 1._R8P + 0.5_R8P * (k - 1)
      do j = 1, nj
         la = (PI / 4 + 0.15_R8P) * (2._R8P * (j - 1) / (nj - 1) - 1)
         do i = 1, ni
            lo = (3 * PI / 4 + 0.15_R8P) * (2._R8P * (i - 1) / (ni - 1) - 1)
            p = r * [cos(la) * cos(lo), cos(la) * sin(lo), sin(la)]
            if (yang) p = [-p(1), p(3), p(2)]
            x(i, j, k) = p(1) ; y(i, j, k) = p(2) ; z(i, j, k) = p(3)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine yin_block

   subroutine box_block(blk)
   !< A Cartesian background block around the sphere, no wall.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(5, 5, 5), y(5, 5, 5), z(5, 5, 5) !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 5
      do j = 1, 5
         do i = 1, 5
            x(i, j, k) = -3 + 1.5_R8P * (i - 1) ; y(i, j, k) = -3 + 1.5_R8P * (j - 1) ; z(i, j, k) = -3 + 1.5_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   endsubroutine box_block
endprogram fossil_test_overset_cutters
