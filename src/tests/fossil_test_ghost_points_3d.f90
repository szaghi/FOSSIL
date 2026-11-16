!< FOSSIL test: ghost points of a finite-difference wall closure on a 3-D block (#109 WP1 step 03).
program fossil_test_ghost_points_3d
!< FOSSIL test: ghost points of a finite-difference wall closure on a 3-D block (#109 WP1 step 03).
!<
!< 1. Hand case, the low-side rule: the unit cube `cube-01.stl` on lines of step 1/4 through its faces (x = -1/2 .. 3/2).
!<    The perturbed points put x = 0 inside and x = 1 outside, so the solid points are {0, 1/4, 1/2, 3/4}^3 (64): with
!<    reach 1 the 56 on the boundary of that block are ghosts, with reach 2 all 64; masks and steps of chosen points.
!< 2. The point flags: `cut_cells_3d_nodes` equals `cut_cells_3d`'s `node_in` bit for bit on the same lines.
!< 3. Brute force: a marching-cubes sphere on generic lines (uniform and stretched, mixed reach), the flags from
!<    `is_point_inside` (no point on the wall) and the ghosts from a direct loop over all points, against
!<    `ghost_points_3d`: count, indexes, coordinates, masks, steps.
!< 4. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path bit for bit on every
!<    case above; an octree surface is refused by the device path (STATUS_INVALID_INPUT), accepted by the host path; an
!<    empty surface has no ghost on either path.
!< 5. Bad input: one point on an axis, lines not increasing, reach 0, size(reach) /= 3: STATUS_INVALID_INPUT, n = 0.
use fossil,          only : surface_stl_object, extract_isosurface, cut_cells_3d_t, cut_cells_3d, ghost_points_t, &
                            ghost_points_3d, STATUS_OK, STATUS_INVALID_INPUT, AABB_TREE_OCTREE, GHOST_XM, GHOST_XP, &
                            GHOST_YM, GHOST_ZM
use fossil_cut_cells_3d, only : cut_cells_3d_nodes
use fossil_facet_object, only : facet_object
use penf,            only : I4P, I8P, R8P
use vecfor,          only : vector_R8P

implicit none

type(surface_stl_object)        :: cube, sphere, oct, empty
type(facet_object), allocatable :: fc(:)
type(ghost_points_t)            :: gh, gd
type(cut_cells_3d_t)            :: cut
real(R8P),          allocatable :: values(:,:,:), x(:), y(:), z(:)
logical,            allocatable :: nd(:,:,:)
real(R8P)                       :: dx, px, py, pz
integer(I4P)                    :: i, j, k, st, g
logical                         :: ok

ok = .true.
! 1. the unit cube on lines through its faces
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call ghost_points_3d(cube, x, x, x, [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('cube reach 1: status', st == STATUS_OK)
call check('cube reach 1: 56 ghosts', gh%n == 56_I4P)
call check('cube reach 1: every ghost in {0, 1/4, 1/2, 3/4}^3', &
           all(gh%points >= 0._R8P) .and. all(gh%points <= 0.75_R8P))
g = find(gh, [3_I4P, 3_I4P, 3_I4P])                    ! the point (0, 0, 0): on three faces, inside by the rule
call check('cube (0,0,0): ghost, mask -x -y -z', g > 0)
if (g > 0) call check('cube (0,0,0): mask and steps', gh%mask(g) == ior(ior(GHOST_XM, GHOST_YM), GHOST_ZM) .and. &
                       all(gh%steps(:, g) == [1_I4P, 0_I4P, 1_I4P, 0_I4P, 1_I4P, 0_I4P]))
g = find(gh, [6_I4P, 5_I4P, 5_I4P])                    ! (3/4, 1/2, 1/2): x = 1 is outside by the rule
call check('cube (3/4,1/2,1/2): ghost, mask +x', g > 0)
if (g > 0) call check('cube (3/4,1/2,1/2): mask and steps', gh%mask(g) == GHOST_XP .and. &
                       all(gh%steps(:, g) == [0_I4P, 1_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P]))
call check('cube (1,1/2,1/2) is not a ghost (fluid)', find(gh, [7_I4P, 5_I4P, 5_I4P]) == 0)
call twin('cube reach 1', cube, x, x, x, [1_I4P, 1_I4P, 1_I4P], gh)
call ghost_points_3d(cube, x, x, x, [2_I4P, 2_I4P, 2_I4P], gh, status=st)
call check('cube reach 2: 64 ghosts', st == STATUS_OK .and. gh%n == 64_I4P)
g = find(gh, [4_I4P, 4_I4P, 4_I4P])                    ! (1/4, 1/4, 1/4): fluid 2 steps back, 3 ahead
if (g > 0) call check('cube (1/4,1/4,1/4) reach 2: steps', all(gh%steps(:, g) == [2_I4P, 0_I4P, 2_I4P, 0_I4P, 2_I4P, 0_I4P]))
call check('cube (1/4,1/4,1/4) reach 2: ghost', g > 0)
call twin('cube reach 2', cube, x, x, x, [2_I4P, 2_I4P, 2_I4P], gh)
call ghost_points_3d(cube, x, x, x, [1_I4P, 2_I4P, 3_I4P], gh, status=st) ! mixed: the reach along z covers all
call check('cube reach (1, 2, 3): 64 ghosts', st == STATUS_OK .and. gh%n == 64_I4P)
call twin('cube reach (1, 2, 3)', cube, x, x, x, [1_I4P, 2_I4P, 3_I4P], gh)
! 2. the point flags are the cut cells' node states
call cut_cells_3d_nodes(cube, x, x, x, nd, status=st)
call cut_cells_3d(cube, x, x, x, cut)
call check('cube: cut_cells_3d_nodes = cut%node_in', st == STATUS_OK .and. all(nd .eqv. cut%node_in))
call check('cube: 64 solid points', count(nd) == 64)
y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)]
call cut_cells_3d_nodes(cube, x, y, x, nd, status=st)
call cut_cells_3d(cube, x, y, x, cut)
call check('cube, other y lines: cut_cells_3d_nodes = cut%node_in', st == STATUS_OK .and. all(nd .eqv. cut%node_in))
! 3. a marching-cubes sphere of radius 1 on generic lines, against a brute force
allocate(values(17, 17, 17))
dx = 4._R8P / 16._R8P
do k = 1, 17 ; do j = 1, 17 ; do i = 1, 17
   px = -2._R8P + (i - 1) * dx ; py = -2._R8P + (j - 1) * dx ; pz = -2._R8P + (k - 1) * dx
   values(i, j, k) = sqrt(px**2 + py**2 + pz**2) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=values, bmin=vector_R8P(-2._R8P, -2._R8P, -2._R8P), &
                        bmax=vector_R8P(2._R8P, 2._R8P, 2._R8P), iso=0._R8P, surface=fc)
call sphere%adopt_facets(facets=fc)
x = [(-1.31_R8P + 0.1373_R8P * i, i = 0, 19)] ; y = [(-1.27_R8P + 0.1411_R8P * i, i = 0, 18)]
z = [(-1.23_R8P + 0.1297_R8P * i, i = 0, 19)]
call brute('sphere reach 1', sphere, x, y, z, [1_I4P, 1_I4P, 1_I4P])
call brute('sphere reach (3, 1, 2)', sphere, x, y, z, [3_I4P, 1_I4P, 2_I4P])
x = [(-1.29_R8P + 2.61_R8P * (real(i, R8P) / 21._R8P)**1.3_R8P, i = 0, 21)]           ! stretched
call brute('sphere stretched, reach 2', sphere, x, y, z, [2_I4P, 2_I4P, 2_I4P])
call brute('sphere stretched, reach 4', sphere, x, y, z, [4_I4P, 4_I4P, 4_I4P])
! 4. an octree surface: the device path needs the SAH BVH
call oct%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call ghost_points_3d(oct, x, x, x, [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('octree, host: 56 ghosts', st == STATUS_OK .and. gh%n == 56_I4P)
call ghost_points_3d(oct, x, x, x, [1_I4P, 1_I4P, 1_I4P], gd, status=st, device=.true.)
call check('octree, device: refused', st == STATUS_INVALID_INPUT .and. gd%n == 0_I4P)
! an empty surface: no ghost on either path
call ghost_points_3d(empty, x, x, x, [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('empty surface: no ghost', st == STATUS_OK .and. gh%n == 0_I4P .and. allocated(gh%ijk))
call ghost_points_3d(empty, x, x, x, [1_I4P, 1_I4P, 1_I4P], gd, status=st, device=.true.)
call check('empty surface, device: no ghost', st == STATUS_OK .and. gd%n == 0_I4P)
! 5. bad input
call ghost_points_3d(cube, x, x, [0.5_R8P], [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('one z point: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P .and. .not. allocated(gh%ijk))
call ghost_points_3d(cube, x, x, x(9:1:-1), [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('decreasing z: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
call ghost_points_3d(cube, x, x, x, [1_I4P, 0_I4P, 1_I4P], gh, status=st)
call check('reach 0: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
call ghost_points_3d(cube, x, x, x, [1_I4P, 1_I4P], gh, status=st, device=.true.)
call check('two reaches: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label, pass)
   !< Record one check.
   character(*), intent(in) :: label !< Check.
   logical,      intent(in) :: pass  !< Result.

   if (.not. pass) then
      ok = .false. ; print '(A)', 'FAILED: '//label
   endif
   endsubroutine check

   integer(I4P) function find(gp, ijk) result(g)
   !< Slot of the point ijk among the ghosts, 0 if not a ghost.
   type(ghost_points_t), intent(in) :: gp     !< Ghost points.
   integer(I4P),         intent(in) :: ijk(3) !< Point.
   integer(I4P)                     :: q      !< Work.

   g = 0_I4P
   do q = 1, gp%n
      if (all(gp%ijk(:, q) == ijk)) then
         g = q ; return
      endif
   enddo
   endfunction find

   integer(I4P) function mismatches(a, b) result(m)
   !< Differences between two ghost sets (coordinates by their bits).
   type(ghost_points_t), intent(in) :: a, b !< Ghost points.

   m = 0_I4P
   if (a%n /= b%n .or. a%dims /= b%dims .or. a%ni /= b%ni .or. a%nj /= b%nj .or. a%nk /= b%nk .or. &
       any(a%reach /= b%reach)) then
      m = 1_I4P ; return
   endif
   m = count(a%ijk /= b%ijk) + count(a%mask /= b%mask) + count(a%steps /= b%steps) + &
       count(transfer(a%points, 1_I8P, size(a%points)) /= transfer(b%points, 1_I8P, size(b%points)))
   endfunction mismatches

   subroutine twin(label, s, x, y, z, reach, host)
   !< The device path against the host result.
   character(*),             intent(in) :: label            !< Case.
   type(surface_stl_object), intent(in) :: s                !< Surface.
   real(R8P),                intent(in) :: x(:), y(:), z(:) !< Lines.
   integer(I4P),             intent(in) :: reach(3)         !< Reach.
   type(ghost_points_t),     intent(in) :: host             !< Host result.
   type(ghost_points_t)                 :: dev              !< Device result.
   integer(I4P)                         :: st2, m           !< Work.

   call ghost_points_3d(s, x, y, z, reach, dev, status=st2, device=.true.)
   m = mismatches(host, dev)
   call check(label//': device twin', st2 == STATUS_OK .and. m == 0_I4P)
   if (m /= 0_I4P) print '(A,I0)', '   mismatches ', m
   endsubroutine twin

   subroutine brute(label, s, x, y, z, reach)
   !< `ghost_points_3d` against a direct loop over all points with flags from `is_point_inside`; then the device twin.
   character(*),             intent(in) :: label            !< Case.
   type(surface_stl_object), intent(in) :: s                !< Surface.
   real(R8P),                intent(in) :: x(:), y(:), z(:) !< Lines.
   integer(I4P),             intent(in) :: reach(3)         !< Reach.
   type(ghost_points_t)                 :: gp               !< Module result.
   logical,      allocatable            :: ins(:,:,:)       !< Brute-force flags.
   integer(I4P)                         :: off(3, 6), ijk(3), q(3), stp(6), msk, n, i, j, k, d, t, st2, m !< Work.

   allocate(ins(size(x), size(y), size(z)))
   do k = 1, size(z) ; do j = 1, size(y) ; do i = 1, size(x)
      ins(i, j, k) = s%is_point_inside(vector_R8P(x(i), y(j), z(k)))
   enddo ; enddo ; enddo
   off = reshape([-1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1], [3, 6])
   call ghost_points_3d(s, x, y, z, reach, gp, status=st2)
   call check(label//': status', st2 == STATUS_OK)
   call cut_cells_3d_nodes(s, x, y, z, nd)
   call check(label//': flags = is_point_inside', all(nd .eqv. ins))
   n = 0_I4P ; m = 0_I4P
   do k = 1, size(z) ; do j = 1, size(y) ; do i = 1, size(x)
      if (.not. ins(i, j, k)) cycle
      msk = 0_I4P ; stp = 0_I4P
      do d = 1, 6
         do t = 1, reach((d + 1) / 2)
            q = [i, j, k] + t * off(:, d)
            if (any(q < 1) .or. q(1) > size(x) .or. q(2) > size(y) .or. q(3) > size(z)) exit
            if (.not. ins(q(1), q(2), q(3))) then
               stp(d) = t ; msk = msk + 2_I4P**(d - 1) ; exit
            endif
         enddo
      enddo
      if (msk == 0_I4P) cycle
      n = n + 1_I4P
      if (n > gp%n) cycle
      ijk = [i, j, k]
      if (any(gp%ijk(:, n) /= ijk) .or. gp%mask(n) /= msk .or. any(gp%steps(:, n) /= stp) .or. &
          any(gp%points(:, n) /= [x(i), y(j), z(k)])) m = m + 1_I4P
   enddo ; enddo ; enddo
   call check(label//': count = brute force', n == gp%n .and. n > 0_I4P)
   call check(label//': ghosts = brute force', m == 0_I4P)
   if (n /= gp%n .or. m /= 0_I4P) print '(A,3(I0,1X))', '   brute, module, mismatches ', n, gp%n, m
   call twin(label, s, x, y, z, reach, gp)
   endsubroutine brute
endprogram fossil_test_ghost_points_3d
