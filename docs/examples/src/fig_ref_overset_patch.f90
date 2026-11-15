!run figref-overset-patch fig_ref_overset_patch
!fig ref-overset-patch
!fig ref-overset-degenerate
!< Figures of the reference page "Overset assembly kernels", section "Conservative interface: the patched interface".
!< 1. ref-overset-patch: the outer boundary (r = 2) of the cylinder O-grid of the Example cuts the Cartesian background
!<    of spacing 0.15 (overset_patch_2d): the holes inside the minor region, the cut cells shaded by their fluid
!<    fraction, the smallest faces (the boundary split at the background's grid lines) alternately gold and ember by
!<    minor face, each with a tick along sf_area (out of the background piece, into the O-grid cell).
!< 2. ref-overset-degenerate: a square ring whose outer boundary lies on grid lines (a 12 x 12 background of unit
!<    cells): the low-side rule puts the boundary into the cells on its solid side, which hold degenerate pieces of zero
!<    area; each is joined to its merge_into piece across the coincident edge.
program fig_ref_overset_patch
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_quad, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter      :: PI = acos(-1._R8P), H = 0.15_R8P, LIFT = 1.e-3_R8P
type(curvilinear_block_2d_t) :: blocks(1)
type(topo_record_t)          :: rec
type(overset_patch_2d_t)     :: patch
type(vector_R8P), allocatable :: a(:), b(:), c(:)
real(R8P)                    :: u(33, 13), v(33, 13), x(41), r, q(2), m(2), y(13), sq(2, 4)
integer(I4P)                 :: u_, i, j, k, p, st, f, ndeg, nmerge, ncut
character(150)               :: caption

! 1. the patched interface of the cylinder O-grid
do j = 1, 13
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
x = [(-3._R8P + H * (i - 1), i = 1, 41)]
call overset_patch_2d(blocks, [1_I4P], x, x, patch, status=st)
ncut = count(patch%cut%pieces > 0 .and. patch%cut%wet_length > 0._R8P)
print '(A,I0,A,I0,A,I0,A,I0)', 'status ', st, ', holes ', count(patch%hole), ', cut cells ', ncut, &
                               ', smallest faces ', patch%ns
write(caption, '(A,I0,A,I0,A,I0,A)') 'holes (teal) ', count(patch%hole), ', cut cells (by solid fraction) ', ncut, &
                                     ', smallest faces ', patch%ns, ' (gold/ember by minor face)'
call fig_open(u_, 'ref-overset-patch', trim(caption), -90._R8P, 90._R8P, width=820_I4P, height=620_I4P)
call fig_zoom(u_, pt([1.25_R8P, 1.25_R8P], 0._R8P), 1.05_R8P)
do j = 1, 40
   do i = 1, 40
      if (patch%hole(i, j)) then
         call fig_quad(u_, 'field', pt([x(i), x(j)], 0._R8P), pt([x(i + 1), x(j)], 0._R8P), &
                       pt([x(i + 1), x(j + 1)], 0._R8P), pt([x(i), x(j + 1)], 0._R8P), -1._R8P)
      elseif (patch%cut%alpha(i, j) < 1._R8P) then
         call fig_quad(u_, 'frac', pt([x(i), x(j)], 0._R8P), pt([x(i + 1), x(j)], 0._R8P), &
                       pt([x(i + 1), x(j + 1)], 0._R8P), pt([x(i), x(j + 1)], 0._R8P), 1._R8P - patch%cut%alpha(i, j))
      endif
   enddo
enddo
allocate(a(0), b(0))
do i = 1, 41 ! the background's grid lines
   a = [a, pt([x(i), x(1)], LIFT), pt([x(1), x(i)], LIFT)] ; b = [b, pt([x(i), x(41)], LIFT), pt([x(41), x(i)], LIFT)]
enddo
call fig_segments(u_, 'box', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do j = 1, 13 ! the O-grid's lines
   do i = 1, 32
      a = [a, pt([u(i, j), v(i, j)], LIFT)] ; b = [b, pt([u(i + 1, j), v(i + 1, j)], LIFT)]
   enddo
enddo
do i = 1, 32
   a = [a, pt([u(i, 1), v(i, 1)], LIFT)] ; b = [b, pt([u(i, 13), v(i, 13)], LIFT)]
enddo
call fig_segments(u_, 'ghost', a, b)
do k = 1, 2 ! the smallest faces, alternately by minor face, and their ticks along sf_area
   deallocate(a, b) ; allocate(a(0), b(0))
   do f = k, patch%nf, 2
      do j = patch%face_sf_start(f), patch%face_sf_start(f + 1) - 1
         p = patch%face_sf(j) ; i = patch%sf_seg(p)
         a = [a, pt(patch%cut%piece_seg(1:2, i), 3 * LIFT)] ; b = [b, pt(patch%cut%piece_seg(3:4, i), 3 * LIFT)]
      enddo
   enddo
   call fig_segments(u_, trim(merge('accent', 'patch ', k == 1)), a, b)
enddo
deallocate(a, b) ; allocate(a(0), b(0), c(0))
do p = 1, patch%ns
   q = patch%sf_area(:, p) / norm2(patch%sf_area(:, p))
   a = [a, pt(patch%sf_centroid(:, p), 4 * LIFT)] ; b = [b, pt(patch%sf_centroid(:, p) + 0.05_R8P * q, 4 * LIFT)]
enddo
call fig_segments(u_, 'probe', a, b)
call fig_label(u_, 'inside', pt([0.65_R8P, 1.15_R8P], 5 * LIFT), 'holes')
call fig_label(u_, 'accent', pt([1.55_R8P, 1.5_R8P], 5 * LIFT), 'minor boundary r = 2')
call fig_label(u_, 'box', pt([1.9_R8P, 2.05_R8P], 5 * LIFT), 'major grid')
call fig_close(u_)

! 2. degenerate pieces: a square ring on grid lines, unit cells
y = [(real(i - 1, R8P), i = 1, 13)]
sq = reshape([3._R8P, 3._R8P, 9._R8P, 3._R8P, 9._R8P, 9._R8P, 3._R8P, 9._R8P], [2, 4])
call ring(sq, blocks(1))
call overset_patch_2d(blocks, [1_I4P], y, y, patch, status=st)
ndeg = count(patch%degenerate) ; nmerge = count(patch%degenerate .and. patch%merge_into > 0)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'square on grid lines: status ', st, ', holes ', count(patch%hole), ', pieces ', &
      patch%cut%np, ', degenerate ', ndeg, ', with merge_into ', nmerge
write(caption, '(A,I0,A,I0,A)') 'degenerate pieces (ember) ', ndeg, ', ', nmerge, &
                                ' joined to merge_into (gold); holes teal'
call fig_open(u_, 'ref-overset-degenerate', trim(caption), -90._R8P, 90._R8P, width=640_I4P, height=600_I4P)
do j = 1, 12
   do i = 1, 12
      if (patch%hole(i, j)) call fig_quad(u_, 'field', pt([y(i), y(j)], 0._R8P), pt([y(i + 1), y(j)], 0._R8P), &
                                          pt([y(i + 1), y(j + 1)], 0._R8P), pt([y(i), y(j + 1)], 0._R8P), -1._R8P)
   enddo
enddo
deallocate(a, b) ; allocate(a(0), b(0))
do i = 1, 13
   a = [a, pt([y(i), y(1)], LIFT), pt([y(1), y(i)], LIFT)] ; b = [b, pt([y(i), y(13)], LIFT), pt([y(13), y(i)], LIFT)]
enddo
call fig_segments(u_, 'box', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do p = 1, patch%ns
   i = patch%sf_seg(p)
   a = [a, pt(patch%cut%piece_seg(1:2, i), 2 * LIFT)] ; b = [b, pt(patch%cut%piece_seg(3:4, i), 2 * LIFT)]
enddo
call fig_segments(u_, 'accent', a, b)
deallocate(a, b) ; allocate(a(0), b(0)) ; deallocate(c) ; allocate(c(0))
do p = 1, patch%cut%np
   if (.not. patch%degenerate(p)) cycle
   m = 0.5_R8P * (y(patch%cut%piece_cell(:, p)) + y(patch%cut%piece_cell(:, p) + 1))
   c = [c, pt(m, 3 * LIFT)]
   if (patch%merge_into(p) == 0) cycle
   q = 0.5_R8P * (y(patch%cut%piece_cell(:, patch%merge_into(p))) + y(patch%cut%piece_cell(:, patch%merge_into(p)) + 1))
   a = [a, pt(m, 3 * LIFT)] ; b = [b, pt(m + 0.8_R8P * (q - m), 3 * LIFT)]
enddo
call fig_segments(u_, 'accent', a, b)
call fig_points(u_, 'patch', c)
call fig_label(u_, 'inside', pt([5.4_R8P, 6._R8P], 5 * LIFT), 'holes')
call fig_close(u_)
contains
   function pt(xy, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: xy(2) !< In-plane coordinates.
   real(R8P), intent(in) :: dz    !< Lift.
   type(vector_R8P)      :: r     !< Point.

   r = xy(1) * ex_R8P + xy(2) * ey_R8P + dz * ez_R8P
   endfunction pt

   subroutine ring(corners, blk)
   !< A block of 2 layers between the square of the given corners (outer, named 'overlap', 6 vertices per side) and its
   !< copy shrunk by half towards its centre.
   real(R8P),                    intent(in)    :: corners(2, 4) !< Square corners, counter-clockwise.
   type(curvilinear_block_2d_t), intent(inout) :: blk           !< Block.
   real(R8P)                                   :: bu(25, 3), bv(25, 3), cc(2), s, t
   type(topo_record_t)                         :: rr
   integer(I4P)                                :: ii, jj, side, l

   cc = sum(corners, dim=2) / 4
   do side = 1, 4
      do l = 1, 6
         ii = (side - 1) * 6 + l ; t = (l - 1) / 6._R8P
         bu(ii, 3) = corners(1, side) + t * (corners(1, mod(side, 4) + 1) - corners(1, side))
         bv(ii, 3) = corners(2, side) + t * (corners(2, mod(side, 4) + 1) - corners(2, side))
      enddo
   enddo
   bu(25, 3) = bu(1, 3) ; bv(25, 3) = bv(1, 3)
   do jj = 1, 2
      s = 0.5_R8P + 0.25_R8P * (jj - 1)
      bu(:, jj) = cc(1) + s * (bu(:, 3) - cc(1)) ; bv(:, jj) = cc(2) + s * (bv(:, 3) - cc(2))
   enddo
   call blk%destroy
   call blk%initialize(bu, bv)
   rr%kind = TOPO_BOUNDARY ; rr%face = 4 ; rr%lo = [1, 3, 1] ; rr%hi = [25, 3, 1] ; rr%name = 'overlap'
   call blk%add_record(rr)
   endsubroutine ring
endprogram fig_ref_overset_patch
