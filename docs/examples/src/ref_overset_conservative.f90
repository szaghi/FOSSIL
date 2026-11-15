!run ref_overset_conservative ref_overset_conservative
program ref_overset_conservative
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d, &
                   overlap_weights_t, overlap_weights_2d
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P), H = 0.15_R8P
type(curvilinear_block_2d_t) :: blocks(1)
type(topo_record_t)          :: rec
type(overset_patch_2d_t)     :: patch
type(overlap_weights_t)      :: w
real(R8P)                    :: u(33, 13), v(33, 13), x(41), r, sw(2), ap(4), closure, cover, area, exact, q, ql
real(R8P)                    :: xc(2), lin_err, sum_err
integer(I4P)                 :: cells(2, 32), i, j, p, k, st, ncut

do j = 1, 13 ! the O-grid of the Example: the unit cylinder out to r = 2
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec) ! the outer boundary: the minor grid's side of the interface
x = [(-3._R8P + H * (i - 1), i = 1, 41)] ! the major grid: [-3, 3]^2, 40 x 40 cells

call overset_patch_2d(blocks, [1_I4P], x, x, patch, status=st) ! A: the patched interface (Wang 1995)
ncut = count(patch%cut%pieces > 0 .and. patch%cut%wet_length > 0._R8P)
print '(A,I0,A,I0,A,I0,A,I0)', 'A: status ', st, ', minor faces ', patch%nf, ', holes ', count(patch%hole), &
                               ', cut cells ', ncut
print '(A,I0,A,I0,A,I0)', '   pieces ', patch%cut%np, ', smallest faces ', patch%ns, ', degenerate ', &
                          count(patch%degenerate)
closure = 0._R8P
do p = 1, patch%cut%np ! every piece closes: sum_f W_f S_f n_f + its smallest faces = 0
   sw = 0._R8P
   do k = patch%piece_sf_start(p), patch%piece_sf_start(p + 1) - 1
      sw = sw + patch%sf_area(:, k)
   enddo
   ap = patch%cut%piece_aperture(:, p)
   closure = max(closure, norm2(H * [ap(2) - ap(1), ap(4) - ap(3)] + sw) / H)
enddo
cover = 0._R8P
do k = 1, patch%nf ! every minor face is the sum of its smallest faces
   cover = max(cover, norm2(patch%face_covered(:, k) - patch%face_vector(:, k)) / norm2(patch%face_vector(:, k)))
enddo
area = sum(patch%cut%piece_area)
exact = 36._R8P - 0.5_R8P * 32 * 4 * sin(2 * PI / 32) ! the box less the 32-gon of radius 2
print '(A,ES8.1,A,ES8.1,A,F12.8,A,ES8.1)', '   closure / h ', closure, ', faces ', cover, ', fluid area ', area, &
                                            ', vs box - polygon ', abs(area - exact) / exact
k = patch%face_sf_start(1)
print '(A,I0,A)', '   minor face 1 (cell (1, 12), +j): ', patch%face_sf_start(2) - k, ' smallest faces'
do k = patch%face_sf_start(1), patch%face_sf_start(2) - 1
   p = patch%sf_piece(patch%face_sf(k))
   print '(A,2I4,A,F8.5,A,2F9.5)', '     major cell', patch%cut%piece_cell(:, p), ', length', &
         norm2(patch%sf_area(:, patch%face_sf(k))), ', sf_area', patch%sf_area(:, patch%face_sf(k))
enddo

cells(1, :) = [(i, i = 1, 32)] ; cells(2, :) = 12 ! B: the outer ring of cells as receptors (Picard et al. 2021)
call overlap_weights_2d(blocks(1), cells, x, x, w, status=st)
sum_err = 0._R8P ; lin_err = 0._R8P
do k = 1, w%nr
   sum_err = max(sum_err, abs(sum(w%weight(w%start(k):w%start(k + 1) - 1)) - 1._R8P))
   ql = 0._R8P ! a linear field q = 1 + 2 x - y, exact gradient: q_c = sum w (q_c' + g . (x_cc' - x_c'))
   do j = w%start(k), w%start(k + 1) - 1
      xc = [x(w%donor(1, j)), x(w%donor(2, j))] + 0.5_R8P * H ! donor cell centre
      q = 1._R8P + 2 * xc(1) - xc(2)
      ql = ql + w%weight(j) * (q + dot_product([2._R8P, -1._R8P], w%centroid(:, j) - xc))
   enddo
   xc = centroid(cells(1, k)) ! the receptor's mean of q is q at its centroid
   lin_err = max(lin_err, abs(ql - (1._R8P + 2 * xc(1) - xc(2))))
enddo
print '(A,I0,A,I0,A,I0,A,L1)', 'B: status ', st, ', receptors ', w%nr, ', intersections ', w%nw, &
                               ', all covered ', all(w%covered == 1._R8P)
print '(A,ES8.1,A,ES8.1)', '   max |sum w - 1| ', sum_err, ', linear field error ', lin_err
print '(A,I0,A,F9.6)', '   receptor 1, cell (1, 12): donors ', w%start(2) - w%start(1), ', area ', w%receptor_volume(1)
do j = w%start(1), w%start(2) - 1
   print '(A,2I4,A,F8.5,A,2F9.5)', '     donor cell', w%donor(:, j), ', weight', w%weight(j), ', centroid', &
         w%centroid(:, j)
enddo
contains
   function centroid(i) result(c)
   !< Centroid of the receptor cell (i, 12), a quadrilateral (shoelace).
   integer(I4P), intent(in) :: i
   real(R8P)                :: c(2), pu(4), pv(4), a, cr
   integer(I4P)             :: l, m

   pu = [u(i, 12), u(i + 1, 12), u(i + 1, 13), u(i, 13)] ; pv = [v(i, 12), v(i + 1, 12), v(i + 1, 13), v(i, 13)]
   a = 0._R8P ; c = 0._R8P
   do l = 1, 4
      m = mod(l, 4) + 1 ; cr = pu(l) * pv(m) - pu(m) * pv(l)
      a = a + cr ; c = c + cr * [pu(l) + pu(m), pv(l) + pv(m)]
   enddo
   c = c / (3 * a)
   endfunction centroid
endprogram ref_overset_conservative
