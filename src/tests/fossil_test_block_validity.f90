!< FOSSIL test: certified validity of the cells of curvilinear blocks (issue #71, step 4).
program fossil_test_block_validity
!< FOSSIL test: certified validity of the cells of curvilinear blocks (issue #71, step 4).
!<
!< 1. a folded cell with det J > 0 at all eight corners but negative on its boundary: the corner test calls it valid,
!<    `hex_validity` invalid with a negative bound;
!< 2. affine cells, the cells of a wavy O-grid (left-handed: orientation -1) and 200 randomly perturbed unit cubes
!<    (corners moved by up to 0.08) are valid, with a positive lower bound below the corner minimum;
!< 3. a flat cell (top = bottom) is invalid; a flat cell on a tilted plane is invalid or undetermined;
!< 4. a polar-axis cell (an edge collapsed onto the axis) is singular;
!< 5. statuses do not depend on how the cell is described: the 48 symmetries of the cube relabel the corners (the 24
!<    reflections with the orientation reversed) and give the same status;
!< 6. 2-D: a convex quad valid, a dart invalid, a collapsed edge singular, a flat quad invalid, exactly; the mirrored
!<    block reverses the orientation and keeps the statuses.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t
use fossil_block_geometry,    only : hex_validity, quad_validity, trilinear_jacobian, det3, CELL_VALID, CELL_SINGULAR, &
                                     CELL_INVALID, CELL_UNDETERMINED
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter         :: PI = acos(-1._R8P)
real(R8P), parameter         :: FOLD(3, 8) = reshape([0.13_R8P, 0.25_R8P, 0.29_R8P, 0.13_R8P, -0.27_R8P, 0.23_R8P, &
                                                      0.58_R8P, 1.25_R8P, -0.61_R8P, 0.76_R8P, 0.86_R8P, -0.18_R8P, &
                                                      0.47_R8P, 0.54_R8P, 1.59_R8P, 0.51_R8P, -0.56_R8P, 0.33_R8P, &
                                                      -0.11_R8P, 1.05_R8P, 0.39_R8P, 1.12_R8P, 1.11_R8P, 0.91_R8P], [3, 8])
type(curvilinear_block_t)    :: b
type(curvilinear_block_2d_t) :: q
real(R8P), allocatable       :: x(:,:,:), y(:,:,:), z(:,:,:), jm(:,:,:), u(:,:), v(:,:), jm2(:,:)
integer(I4P), allocatable    :: st(:,:,:), st2(:,:)
real(R8P)                    :: c(3, 8), cube(3, 8), jmin, cmin, r(3), c2(2, 4), h
integer(I4P)                 :: s, i, j, k, n, ni, nj, nk, nbad
logical                      :: ok, all_ok

all_ok = .true.

! 1. folded cell
cmin = corner_min(FOLD)
call hex_validity(FOLD, 1._R8P, s, jmin)
call check(cmin > 0._R8P, 'folded cell: every corner det J > 0', cmin)
call check(s == CELL_INVALID .and. jmin < 0._R8P, 'folded cell: certified invalid, negative det J found', jmin)

! 2. valid cells
do k = 0, 1
   do j = 0, 1
      do i = 0, 1
         cube(:, 1 + i + 2 * j + 4 * k) = [i, j, k] * 1._R8P
      enddo
   enddo
enddo
c = matmul(reshape([1.2_R8P, 0.1_R8P, -0.2_R8P, 0.3_R8P, 0.9_R8P, 0.15_R8P, -0.1_R8P, 0.25_R8P, 1.1_R8P], [3, 3]), cube)
call hex_validity(c, 1._R8P, s, jmin)
call check(s == CELL_VALID .and. abs(jmin - corner_min(c)) <= 1e-12_R8P, 'affine cell valid, bound = det A', jmin)
call random_seed(put=[(17 + 3 * i, i = 1, 64)])
nbad = 0
do n = 1, 200
   do i = 1, 8
      call random_number(r) ; c(:, i) = cube(:, i) + 0.16_R8P * (r - 0.5_R8P)
   enddo
   call hex_validity(c, 1._R8P, s, jmin)
   if (s /= CELL_VALID .or. jmin <= 0._R8P .or. jmin > corner_min(c)) nbad = nbad + 1
enddo
call check(nbad == 0, '200 perturbed cubes valid, 0 < bound <= corner minimum', real(nbad, R8P))
ni = 25 ; nj = 7 ; nk = 5
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         h = (1._R8P + 0.15_R8P * (j - 1)) * (1._R8P + 0.05_R8P * sin(3 * 2 * PI * (i - 1) / (ni - 1) + 0.4_R8P * k))
         x(i, j, k) = h * cos(2 * PI * (i - 1) / (ni - 1)) ; y(i, j, k) = h * sin(2 * PI * (i - 1) / (ni - 1))
         z(i, j, k) = 0.2_R8P * (k - 1) + 0.03_R8P * cos(2 * PI * (i - 1) / (ni - 1)) * j
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%validity(st, jm)
call check(b%orientation() < 0._R8P .and. all(st == CELL_VALID) .and. all(jm > 0._R8P), &
           'wavy O-grid: left-handed, every cell valid', minval(jm))

! 3. flat cells
c = cube ; c(3, 5:8) = 0._R8P
call hex_validity(c, 1._R8P, s, jmin)
call check(s == CELL_INVALID, 'flat cell: invalid', jmin)
c = cube ; c(3, 5:8) = 0._R8P
do i = 1, 8
   c(3, i) = 0.3_R8P * c(1, i) - 0.7_R8P * c(2, i)
enddo
call hex_validity(c, 1._R8P, s, jmin)
call check(s == CELL_INVALID .or. s == CELL_UNDETERMINED, 'flat cell on a tilted plane: invalid or undetermined', jmin)

! 4. polar axis cell: vertices (r, theta, z), r = 0 collapsed
do k = 0, 1
   do j = 0, 1
      do i = 0, 1
         h = 0.3_R8P * j ; cmin = 0.4_R8P * i
         c(:, 1 + i + 2 * j + 4 * k) = [h * cos(cmin), h * sin(cmin), 0.5_R8P * k]
      enddo
   enddo
enddo
call hex_validity(c, sign(1._R8P, det3(trilinear_jacobian(c, [0.5_R8P, 0.5_R8P, 0.5_R8P]))), s, jmin)
call check(s == CELL_SINGULAR .and. jmin == 0._R8P, 'polar-axis cell: singular', jmin)

! 5. relabelling
ok = same_under_symmetries(FOLD, CELL_INVALID)
c = cube
do i = 1, 8
   c(:, i) = c(:, i) + 0.07_R8P * [sin(3._R8P * i), cos(5._R8P * i), sin(7._R8P * i + 1)]
enddo
ok = ok .and. same_under_symmetries(c, CELL_VALID)
call check(ok, 'statuses identical under the 48 relabellings of the cube', 0._R8P)

! 6. 2-D
c2 = reshape([0._R8P, 0._R8P, 1._R8P, 0.1_R8P, 0.2_R8P, 1._R8P, 1.1_R8P, 0.9_R8P], [2, 4])
call quad_validity(c2, 1._R8P, s, jmin) ; ok = s == CELL_VALID
c2(:, 4) = [0.3_R8P, 0.3_R8P] ! dart: corner 4 pulled inside
call quad_validity(c2, 1._R8P, s, jmin) ; ok = ok .and. s == CELL_INVALID
c2(:, 4) = c2(:, 3) ! collapsed edge: a triangle
call quad_validity(c2, 1._R8P, s, jmin) ; ok = ok .and. s == CELL_SINGULAR
c2 = reshape([0._R8P, 0._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P, 3._R8P, 3._R8P], [2, 4]) ! collinear
call quad_validity(c2, 1._R8P, s, jmin) ; ok = ok .and. s == CELL_INVALID
call check(ok, '2-D: valid, dart, triangle, flat', 0._R8P)
ni = 9 ; nj = 5
allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   do i = 1, ni
      u(i, j) = (1 + 0.2_R8P * j) * cos(PI * (i - 1) / (ni - 1)) ; v(i, j) = (1 + 0.2_R8P * j) * sin(PI * (i - 1) / (ni - 1))
   enddo
enddo
call q%initialize(u, v) ; call q%validity(st2, jm2)
ok = all(st2 == CELL_VALID) .and. q%orientation() < 0._R8P
u(5, 3) = 0._R8P ; v(5, 3) = 0.9_R8P ! fold the cells around vertex (5, 3)
call q%initialize(u, v) ; call q%validity(st2)
ok = ok .and. count(st2 == CELL_INVALID) > 0
call q%initialize(-u, v) ; call q%validity(st2, jm2)
ok = ok .and. q%orientation() > 0._R8P .and. count(st2 == CELL_INVALID) > 0
call check(ok, '2-D block: valid half annulus, folded vertex caught, mirror keeps statuses', 0._R8P)

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

   function corner_min(c) result(m)
   !< Minimum det J over the eight corners.
   real(R8P), intent(in) :: c(3, 8) !< Corners.
   real(R8P)             :: m       !< Minimum.
   integer(I4P)          :: a, b, e !< Bits.

   m = huge(1._R8P)
   do e = 0, 1
      do b = 0, 1
         do a = 0, 1
            m = min(m, det3(trilinear_jacobian(c, [a, b, e] * 1._R8P)))
         enddo
      enddo
   enddo
   endfunction corner_min

   function same_under_symmetries(c, expected) result(same)
   !< The status of the cell relabelled by each of the 48 symmetries of the cube equals `expected`.
   real(R8P),    intent(in) :: c(3, 8)  !< Corners.
   integer(I4P), intent(in) :: expected !< Status.
   logical                  :: same     !< Result.
   integer(I4P)             :: perm(3, 6), p, sg, a, b, e, w(3), vv(3), st !< Work.
   real(R8P)                :: cc(3, 8), jj, mat(3, 3) !< Relabelled cell.

   perm = reshape([1, 2, 3, 2, 3, 1, 3, 1, 2, 2, 1, 3, 1, 3, 2, 3, 2, 1], [3, 6])
   same = .true.
   do p = 1, 6
      do sg = 0, 7
         mat = 0._R8P
         do a = 1, 3
            mat(perm(a, p), a) = merge(-1._R8P, 1._R8P, btest(sg, a - 1))
         enddo
         do e = 0, 1
            do b = 0, 1
               do a = 0, 1
                  vv = [a, b, e]
                  ! new label w = mat (v - 1/2) + 1/2
                  w = nint(matmul(mat, vv - 0.5_R8P) + 0.5_R8P)
                  cc(:, 1 + w(1) + 2 * w(2) + 4 * w(3)) = c(:, 1 + a + 2 * b + 4 * e)
               enddo
            enddo
         enddo
         call hex_validity(cc, sign(1._R8P, det3(mat)), st, jj)
         same = same .and. st == expected
      enddo
   enddo
   endfunction same_under_symmetries
endprogram fossil_test_block_validity
