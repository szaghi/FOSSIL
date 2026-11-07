!< FOSSIL test: the links between the cut-cell pieces (issue #70, step 3).
program fossil_test_cut_links
!< FOSSIL test: the links between the cut-cell pieces (issue #70, step 3).
!<
!< 2-D, a circle (64 segments) and a thin plate (two pieces per cut cell) on generic and grid-aligned grids:
!< 1. per piece and face, the links sum to the piece's aperture times the edge length (relative 1e-14), and a face on
!<    the block boundary has no link;
!< 2. symmetry: p -> q over a face exists exactly when q -> p does, with the same length;
!< 3. the per-piece apertures add up to the cell apertures;
!< 4. the thin plate: in the 5 cells it crosses (two pieces each) a piece above the plate never links to a piece below
!<    it among those cells (the cells at its ends hold one piece wrapping around the end).

use fossil, only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d
use penf,   only : I4P, R8P

implicit none

type(curve_object)   :: c
type(cut_cells_2d_t) :: cut
type(cut_links_t)    :: lk
real(R8P), allocatable :: x(:), y(:), u(:), v(:)
real(R8P)            :: emax(3)
integer(I4P)         :: i, n, nbad, nlinks
logical              :: ok

emax = 0._R8P ; nbad = 0 ; nlinks = 0
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call check2(c, x, y)
x = [(0.1_R8P * (i - 1), i = 1, 11)]
call check2(c, x, x)
! thin plate [0.2, 0.8] x [0.4995, 0.5005], counter-clockwise
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P], v=[0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P], loop_sizes=[4_I4P])
x = [(0.05_R8P + 0.1_R8P * (i - 1), i = 1, 10)]
call check2(c, x, x)
call plate2(c, x)
print '(A,I0,A,I0)', '2-D links ', nlinks, '  bad ', nbad
print '(A,3ES10.2)', 'max link-sum / symmetry / aperture-sum errors: ', emax
ok = nbad == 0 .and. nlinks > 0 .and. all(emax <= 1.e-14_R8P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check2(cv, xx, yy)
   !< Checks 1-3 on the block xx x yy.
   type(curve_object), intent(in) :: cv
   real(R8P),          intent(in) :: xx(:), yy(:)
   integer(I4P)                   :: p, f, k, q, kk, i0, j0, nx, ny
   real(R8P)                      :: s, len(4), ap(4)
   logical                        :: border, found

   call cv%cut_cells(x=xx, y=yy, cut=cut)
   call cut_links_2d(cut, lk)
   nx = size(xx) - 1 ; ny = size(yy) - 1
   nlinks = nlinks + size(lk%piece)
   do p = 1, cut%np
      i0 = cut%piece_cell(1, p) ; j0 = cut%piece_cell(2, p)
      len = [yy(j0 + 1) - yy(j0), yy(j0 + 1) - yy(j0), xx(i0 + 1) - xx(i0), xx(i0 + 1) - xx(i0)]
      do f = 1, 4
         s = 0._R8P
         do k = lk%start(p), lk%start(p + 1) - 1
            if (lk%face(k) == f) s = s + lk%area(k)
         enddo
         border = (f == 1 .and. i0 == 1) .or. (f == 2 .and. i0 == nx) .or. (f == 3 .and. j0 == 1) .or. (f == 4 .and. j0 == ny)
         if (border) then
            if (s /= 0._R8P) nbad = nbad + 1
         else
            emax(1) = max(emax(1), abs(s - cut%piece_aperture(f, p) * len(f)) / len(f))
         endif
      enddo
      do k = lk%start(p), lk%start(p + 1) - 1
         q = lk%piece(k) ; found = .false.
         do kk = lk%start(q), lk%start(q + 1) - 1
            if (lk%piece(kk) == p) then
               found = .true. ; emax(2) = max(emax(2), abs(lk%area(kk) - lk%area(k)) / len(lk%face(k)))
            endif
         enddo
         if (.not. found) nbad = nbad + 1
      enddo
   enddo
   do j0 = 1, ny
      do i0 = 1, nx
         if (cut%pieces(i0, j0) < 1) cycle
         ap = 0._R8P
         do p = cut%first_piece(i0, j0), cut%first_piece(i0, j0) + cut%pieces(i0, j0) - 1
            ap = ap + cut%piece_aperture(:, p)
         enddo
         emax(3) = max(emax(3), maxval(abs(ap - cut%aperture(:, i0, j0))))
      enddo
   enddo
   endsubroutine check2

   subroutine plate2(cv, xx)
   !< Check 4: in the cut cells of the plate the pieces above the plate never link to a piece below it.
   type(curve_object), intent(in) :: cv
   real(R8P),          intent(in) :: xx(:)
   integer(I4P)                   :: p, k, q, ntwo

   call cv%cut_cells(x=xx, y=xx, cut=cut)
   call cut_links_2d(cut, lk)
   ntwo = 0
   do p = 1, cut%np
      if (cut%pieces(cut%piece_cell(1, p), cut%piece_cell(2, p)) == 2) ntwo = ntwo + 1
      do k = lk%start(p), lk%start(p + 1) - 1
         q = lk%piece(k)
         if (two(p) .and. two(q) .and. &
             (cut%piece_centroid(2, p) - 0.5_R8P) * (cut%piece_centroid(2, q) - 0.5_R8P) < 0._R8P) nbad = nbad + 1
      enddo
   enddo
   if (ntwo /= 10) nbad = nbad + 1 ! 5 cells across the plate, two pieces each
   endsubroutine plate2

   logical function two(pp)
   !< Piece pp lies in a cell with two pieces.
   integer(I4P), intent(in) :: pp

   two = cut%pieces(cut%piece_cell(1, pp), cut%piece_cell(2, pp)) == 2
   endfunction two
endprogram fossil_test_cut_links
