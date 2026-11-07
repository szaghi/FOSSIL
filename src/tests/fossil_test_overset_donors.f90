!< FOSSIL test: donor candidates, quality and the pair matrix (issue #72, step 5).
program fossil_test_overset_donors
!< FOSSIL test: donor candidates, quality and the pair matrix (issue #72, step 5).
!<
!< 1. 2-D O-grid around the unit cylinder in a Cartesian background: the pair matrix says they overlap; every
!<    background vertex's candidate in the O-grid equals a brute force over all O-grid cells (the first cell that
!<    contains it, exact orient2d), and the O-grid's candidates in the background likewise; candidate xi maps back to
!<    the vertex; the receptor volume of a uniform background vertex is h^2, the ratio donor / receptor;
!< 2. the wall distance of each background vertex equals the distance to the cutter polygon;
!< 3. two arcs sharing the cylinder wall are PAIR_SHARED; a far block is PAIR_NONE;
!< 4. 3-D: a cylindrical O-grid in a Cartesian background, candidates equal a brute force over all cells (Newton
!<    containment), xi maps back to 1e-12;
!< 5. device=.true. gives the same candidates, bit for bit (the blocks entered once on the device).

use fossil,         only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, &
                           wall_cutter_2d, LOC_INSIDE
use fossil_overset, only : donor_set_t, overlap_matrix, overlap_matrix_2d, donor_candidates, donor_candidates_2d, &
                           PAIR_NONE, PAIR_OVERLAP, PAIR_SHARED
use fossil_block_geometry, only : bilinear_point, trilinear_point
use fossil_block_locate,   only : newton_3d_oac, in_box3, XI_TOL
use fossil_predicates,     only : orient2d_xy
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(donor_set_t),            allocatable :: dn(:), dd(:)
integer(I4P),                 allocatable :: pair(:,:)
logical,                      allocatable :: sh(:,:)
real(R8P)                                 :: err, d(1)
integer(I4P)                              :: nbad
logical                                   :: all_ok

all_ok = .true.

! 1-2. 2-D
allocate(q(2))
call o_grid(q(1), 0._R8P, 2 * PI, 49_I4P) ; call background_2d(q(2))
call q(1)%build_search ; call q(2)%build_search
call overlap_matrix_2d(q, pair)
call check(pair(1, 2) == PAIR_OVERLAP .and. pair(2, 1) == PAIR_OVERLAP, '2-D: O-grid and background overlap', 0._R8P)
call wall_cutter_2d(q(1:1), c)
call donor_candidates_2d(q, dn, pair=pair, wall=c)
call brute_2d(q, dn, nbad, err)
call check(nbad == 0 .and. err <= 1e-13_R8P, '2-D: candidates = brute force, xi maps back', err)
err = maxval(abs(dn(2)%receptor_volume - 0.1_R8P ** 2))
call check(err <= 1e-15_R8P .and. all(abs(dn(2)%ratio - dn(2)%donor_volume / 0.1_R8P ** 2) <= 1e-12_R8P), &
           '2-D: receptor volume h^2, ratio donor / receptor', err)
call c%distance_many(reshape([q(2)%u(3, 7), q(2)%v(3, 7)], [2, 1]), d, is_signed=.false.)
call check(abs(dn(2)%wall_distance(3 + 41 * 6) - abs(d(1))) <= 1e-14_R8P, '2-D: wall distance', abs(d(1)))
call donor_candidates_2d(q, dd, pair=pair, device=.true.)
call check(same(dn, dd), '2-D: device candidates = CPU', 0._R8P)

! 3. shared and none
deallocate(q) ; allocate(q(3))
call o_grid(q(1), -0.2_R8P * PI, 1.4_R8P * PI, 41_I4P) ; call o_grid(q(2), 0.8_R8P * PI, 1.4_R8P * PI, 21_I4P)
call far_2d(q(3))
call q(1)%build_search ; call q(2)%build_search ; call q(3)%build_search
call wall_cutter_2d(q, c, shared=sh)
call overlap_matrix_2d(q, pair, shared=sh)
call check(pair(1, 2) == PAIR_SHARED .and. pair(1, 3) == PAIR_NONE .and. pair(2, 3) == PAIR_NONE, &
           'shared wall: PAIR_SHARED; a far block: PAIR_NONE', 0._R8P)

! 4-5. 3-D
allocate(b(2))
call cylinder_3d(b(1)) ; call background_3d(b(2))
call b(1)%build_search ; call b(2)%build_search
call overlap_matrix(b, pair)
call donor_candidates(b, dn, pair=pair)
call brute_3d(b, dn, nbad, err)
call check(pair(1, 2) == PAIR_OVERLAP .and. nbad == 0 .and. err <= 1e-12_R8P, &
           '3-D: candidates = brute force, xi maps back', err)
call donor_candidates(b, dd, pair=pair, device=.true.)
call check(same(dn, dd), '3-D: device candidates = CPU', 0._R8P)

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

   logical function same(a, b) result(yes)
   !< Two candidate sets equal bit for bit.
   type(donor_set_t), intent(in) :: a(:), b(:) !< Sets.
   integer(I4P)                  :: k          !< Counter.

   yes = .true.
   do k = 1, size(a)
      yes = yes .and. all(a(k)%start == b(k)%start) .and. all(a(k)%block == b(k)%block) .and. &
            all(a(k)%cell == b(k)%cell) .and. all(a(k)%xi == b(k)%xi) .and. all(a(k)%status == b(k)%status)
   enddo
   endfunction same

   subroutine brute_2d(q, dn, nbad, err)
   !< Each vertex: the first cell of the other block containing it (exact), against its candidate.
   type(curvilinear_block_2d_t), intent(in)  :: q(:)   !< Blocks.
   type(donor_set_t),            intent(in)  :: dn(:)  !< Candidates.
   integer(I4P),                 intent(out) :: nbad   !< Mismatches.
   real(R8P),                    intent(out) :: err    !< Map-back error.
   integer(I4P)                              :: a, o, v, i, j, ci, cj, first(2), cand !< Counters.
   real(R8P)                                 :: p(2), cc(2, 4), s !< Work.

   nbad = 0 ; err = 0._R8P
   do a = 1, 2
      o = 3 - a ; v = 0
      do j = 1, q(a)%nj
         do i = 1, q(a)%ni
            v = v + 1 ; p = [q(a)%u(i, j), q(a)%v(i, j)] ; first = 0
            outer: do cj = 1, q(o)%nj - 1
               do ci = 1, q(o)%ni - 1
                  cc = q(o)%corners(ci, cj) ; s = q(o)%orient
                  if (s * orient2d_xy(cc(1, 1), cc(2, 1), cc(1, 2), cc(2, 2), p(1), p(2)) < 0) cycle
                  if (s * orient2d_xy(cc(1, 2), cc(2, 2), cc(1, 4), cc(2, 4), p(1), p(2)) < 0) cycle
                  if (s * orient2d_xy(cc(1, 4), cc(2, 4), cc(1, 3), cc(2, 3), p(1), p(2)) < 0) cycle
                  if (s * orient2d_xy(cc(1, 3), cc(2, 3), cc(1, 1), cc(2, 1), p(1), p(2)) < 0) cycle
                  first = [ci, cj] ; exit outer
               enddo
            enddo outer
            cand = 0
            if (dn(a)%start(v + 1) > dn(a)%start(v)) cand = dn(a)%start(v)
            if (first(1) == 0) then
               if (cand /= 0) nbad = nbad + 1
            else
               if (cand == 0) then
                  nbad = nbad + 1 ; cycle
               endif
               if (any(dn(a)%cell(:, cand) /= first) .or. dn(a)%status(cand) /= LOC_INSIDE) nbad = nbad + 1
               err = max(err, maxval(abs(bilinear_point(q(o)%corners(first(1), first(2)), dn(a)%xi(:, cand)) - p)))
            endif
         enddo
      enddo
   enddo
   endsubroutine brute_2d

   subroutine brute_3d(b, dn, nbad, err)
   !< Each background vertex: the first O-grid cell containing it (Newton), against its candidate.
   type(curvilinear_block_t), intent(in)  :: b(:)   !< Blocks.
   type(donor_set_t),         intent(in)  :: dn(:)  !< Candidates.
   integer(I4P),              intent(out) :: nbad   !< Mismatches.
   real(R8P),                 intent(out) :: err    !< Map-back error.
   integer(I4P)                           :: v, i, j, k, ci, cj, ck, first(3), cand !< Counters.
   real(R8P)                              :: p(3), xi(3) !< Work.
   logical                                :: conv   !< Newton.

   nbad = 0 ; err = 0._R8P ; v = 0
   do k = 1, b(2)%nk
      do j = 1, b(2)%nj
         do i = 1, b(2)%ni
            v = v + 1 ; p = [b(2)%x(i, j, k), b(2)%y(i, j, k), b(2)%z(i, j, k)] ; first = 0
            outer: do ck = 1, b(1)%nk - 1
               do cj = 1, b(1)%nj - 1
                  do ci = 1, b(1)%ni - 1
                     call newton_3d_oac(b(1)%corners(ci, cj, ck), p, xi, conv)
                     if (conv .and. in_box3(xi, XI_TOL)) then
                        first = [ci, cj, ck] ; exit outer
                     endif
                  enddo
               enddo
            enddo outer
            cand = 0
            if (dn(2)%start(v + 1) > dn(2)%start(v)) cand = dn(2)%start(v)
            if (first(1) == 0) then
               if (cand /= 0) nbad = nbad + 1
            else
               if (cand == 0) then
                  nbad = nbad + 1 ; cycle
               endif
               if (any(dn(2)%cell(:, cand) /= first)) nbad = nbad + 1
               err = max(err, maxval(abs(trilinear_point(b(1)%corners(first(1), first(2), first(3)), dn(2)%xi(:, cand)) &
                                         - p)))
            endif
         enddo
      enddo
   enddo
   endsubroutine brute_3d

   subroutine o_grid(blk, t0, span, ni)
   !< Annulus 1 <= r <= 1.5 over [t0, t0 + span], wall on face 3 (seam copied on a full circle).
   type(curvilinear_block_2d_t), intent(inout) :: blk      !< Block.
   real(R8P),                    intent(in)    :: t0, span !< Arc.
   integer(I4P),                 intent(in)    :: ni       !< Vertices along the arc.
   real(R8P)                                   :: u(ni, 7), v(ni, 7), r, t !< Vertices.
   integer(I4P)                                :: i, j     !< Counters.
   type(topo_record_t)                         :: rec      !< Wall.

   do j = 1, 7
      r = 1._R8P + 0.5_R8P * (j - 1) / 6
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
   endsubroutine o_grid

   subroutine background_2d(blk)
   !< Cartesian block [-2.03, -2.03] + 0.1 (i - 1, j - 1), 41 x 41.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: u(41, 41), v(41, 41) !< Vertices.
   integer(I4P)                                :: i, j !< Counters.

   do j = 1, 41
      do i = 1, 41
         u(i, j) = -2.03_R8P + 0.1_R8P * (i - 1) ; v(i, j) = -2.03_R8P + 0.1_R8P * (j - 1)
      enddo
   enddo
   call blk%initialize(u, v)
   endsubroutine background_2d

   subroutine far_2d(blk)
   !< A block far away.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: u(3, 3), v(3, 3) !< Vertices.
   integer(I4P)                                :: i, j !< Counters.

   do j = 1, 3
      do i = 1, 3
         u(i, j) = 10._R8P + i ; v(i, j) = 10._R8P + j
      enddo
   enddo
   call blk%initialize(u, v)
   endsubroutine far_2d

   subroutine cylinder_3d(blk)
   !< Cylindrical O-grid 1 <= r <= 1.5, 0 <= z <= 1 (seam copied).
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(25, 5, 4), y(25, 5, 4), z(25, 5, 4), r !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 4
      do j = 1, 5
         r = 1._R8P + 0.125_R8P * (j - 1)
         do i = 1, 25
            x(i, j, k) = r * cos(2 * PI * (i - 1) / 24) ; y(i, j, k) = r * sin(2 * PI * (i - 1) / 24)
            z(i, j, k) = (k - 1) / 3._R8P
         enddo
      enddo
   enddo
   x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:)
   call blk%initialize(x, y, z)
   endsubroutine cylinder_3d

   subroutine background_3d(blk)
   !< Cartesian block around the cylinder, 13 x 13 x 6.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(13, 13, 6), y(13, 13, 6), z(13, 13, 6) !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 6
      do j = 1, 13
         do i = 1, 13
            x(i, j, k) = -1.83_R8P + 0.31_R8P * (i - 1) ; y(i, j, k) = -1.83_R8P + 0.31_R8P * (j - 1)
            z(i, j, k) = -0.13_R8P + 0.25_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   endsubroutine background_3d
endprogram fossil_test_overset_donors
