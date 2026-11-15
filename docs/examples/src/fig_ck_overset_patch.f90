!run figck-overset-patch fig_ck_overset_patch
!fig ck-overset-patch
!< Figure of the cookbook recipe "Couple an overset grid conservatively": the cylinder O-grid over the Cartesian
!< background of spacing 0.15, twice. Left, A (overset_patch_2d): the background holes in teal, the cut cells by
!< their solid fraction, the smallest faces in gold. Right, B (overlap_weights_2d, drawn 6.6 to the right): the overlap
!< kept, each background cell shaded by the fraction of it the O-grid's outer ring of receptors covers.
program fig_ck_overset_patch
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d, &
                       overlap_weights_t, overlap_weights_2d
use fossil_fig, only : fig_open, fig_segments, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), H = 0.15_R8P, LIFT = 1.e-3_R8P, SHIFT = 6.6_R8P
type(curvilinear_block_2d_t)  :: blocks(1)
type(topo_record_t)           :: rec
type(overset_patch_2d_t)      :: patch
type(overlap_weights_t)       :: w
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: u(33, 13), v(33, 13), x(41), r, cov(40, 40), dx
integer(I4P)                  :: u_, i, j, k, p, st, cells(2, 32), panel

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
cells(1, :) = [(i, i = 1, 32)] ; cells(2, :) = 12
call overlap_weights_2d(blocks(1), cells, x, x, w, status=st)
cov = 0._R8P
do k = 1, w%nw
   cov(w%donor(1, k), w%donor(2, k)) = cov(w%donor(1, k), w%donor(2, k)) + w%volume(k) / H ** 2
enddo
print '(A,I0,A,I0,A,I0,A,F6.4)', 'A: holes ', count(patch%hole), ', smallest faces ', patch%ns, &
                                 '; B: background cells under the ring ', count(cov > 0._R8P), ', most covered ', maxval(cov)
call fig_open(u_, 'ck-overset-patch', 'A: holes, cut cells, smallest faces (gold);  B: the overlap kept, '// &
              'cells by the share the receptor ring covers', -90._R8P, 90._R8P, width=900_I4P, height=470_I4P)
do j = 1, 40
   do i = 1, 40
      if (patch%hole(i, j)) then
         call fig_quad(u_, 'field', pt([x(i), x(j)], 0._R8P), pt([x(i + 1), x(j)], 0._R8P), &
                       pt([x(i + 1), x(j + 1)], 0._R8P), pt([x(i), x(j + 1)], 0._R8P), -1._R8P)
      elseif (patch%cut%alpha(i, j) < 1._R8P) then
         call fig_quad(u_, 'frac', pt([x(i), x(j)], 0._R8P), pt([x(i + 1), x(j)], 0._R8P), &
                       pt([x(i + 1), x(j + 1)], 0._R8P), pt([x(i), x(j + 1)], 0._R8P), 1._R8P - patch%cut%alpha(i, j))
      endif
      if (cov(i, j) > 0._R8P) &
         call fig_quad(u_, 'frac', pt([x(i) + SHIFT, x(j)], 0._R8P), pt([x(i + 1) + SHIFT, x(j)], 0._R8P), &
                       pt([x(i + 1) + SHIFT, x(j + 1)], 0._R8P), pt([x(i) + SHIFT, x(j + 1)], 0._R8P), cov(i, j))
   enddo
enddo
do panel = 0, 1
   dx = panel * SHIFT
   allocate(a(0), b(0))
   do i = 1, 41 ! the background's lines
      a = [a, pt([x(i) + dx, x(1)], LIFT), pt([x(1) + dx, x(i)], LIFT)]
      b = [b, pt([x(i) + dx, x(41)], LIFT), pt([x(41) + dx, x(i)], LIFT)]
   enddo
   call fig_segments(u_, 'box', a, b)
   deallocate(a, b) ; allocate(a(0), b(0))
   do j = 1, 13, 12 ! the O-grid's wall and outer boundary
      do i = 1, 32
         a = [a, pt([u(i, j) + dx, v(i, j)], 2 * LIFT)] ; b = [b, pt([u(i + 1, j) + dx, v(i + 1, j)], 2 * LIFT)]
      enddo
   enddo
   if (panel == 1) then
      do i = 1, 32 ! the ring of receptors
         a = [a, pt([u(i, 12) + dx, v(i, 12)], 2 * LIFT)] ; b = [b, pt([u(i + 1, 12) + dx, v(i + 1, 12)], 2 * LIFT)]
         a = [a, pt([u(i, 12) + dx, v(i, 12)], 2 * LIFT)] ; b = [b, pt([u(i, 13) + dx, v(i, 13)], 2 * LIFT)]
      enddo
   endif
   call fig_segments(u_, trim(merge('patch', 'ghost', panel == 1)), a, b)
   deallocate(a, b)
enddo
allocate(a(0), b(0))
do p = 1, patch%ns
   i = patch%sf_seg(p)
   a = [a, pt(patch%cut%piece_seg(1:2, i), 3 * LIFT)] ; b = [b, pt(patch%cut%piece_seg(3:4, i), 3 * LIFT)]
enddo
call fig_segments(u_, 'accent', a, b)
call fig_label(u_, 'accent', pt([-2.9_R8P, 2.75_R8P], 5 * LIFT), 'A: patched interface')
call fig_label(u_, 'accent', pt([-2.9_R8P + SHIFT, 2.75_R8P], 5 * LIFT), 'B: intersection weights')
call fig_close(u_)
contains
   function pt(xy, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: xy(2) !< In-plane coordinates.
   real(R8P), intent(in) :: dz    !< Lift.
   type(vector_R8P)      :: r     !< Point.

   r = xy(1) * ex_R8P + xy(2) * ey_R8P + dz * ez_R8P
   endfunction pt
endprogram fig_ck_overset_patch
