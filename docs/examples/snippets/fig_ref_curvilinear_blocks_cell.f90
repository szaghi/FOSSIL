!< Figure of the reference page "Curvilinear structured blocks", section Cell geometry: one trilinear cell, a unit cube
!< whose top face is twisted and lifted at a corner, so its faces are curved (bilinear). Drawn: the faces (translucent,
!< sampled 8 x 8 to show their curvature), the 2 x 2 x 2 Gauss points of the volume rule mapped into the cell, and the
!< six outward face area vectors 1/2 (c - a) x (d - b) from the face centres (scaled by 0.5).
program fig_ref_curvilinear_blocks_cell
use fossil,     only : curvilinear_block_t
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter   :: NS = 8
real(R8P),    parameter   :: G1 = 0.5_R8P - 0.5_R8P / sqrt(3._R8P), G2 = 0.5_R8P + 0.5_R8P / sqrt(3._R8P)
type(curvilinear_block_t) :: b
real(R8P)                 :: x(2, 2, 2), y(2, 2, 2), z(2, 2, 2), gs(2), fv(3, 6), s, t
real(R8P), allocatable    :: vol(:,:,:), cen(:,:,:,:), fi(:,:,:,:), fj(:,:,:,:), fk(:,:,:,:)
type(vector_R8P)          :: gp(8), fc(6), fe(6)
integer(I4P)              :: a, c, d, f, m, w, la, lb

x(1,:,:) = 0._R8P ; x(2,:,:) = 1._R8P ; y(:,1,:) = 0._R8P ; y(:,2,:) = 1._R8P ; z(:,:,1) = 0._R8P ; z(:,:,2) = 1._R8P
x(:,:,2) = x(:,:,2) + reshape([0.25_R8P, 0.25_R8P, -0.25_R8P, -0.25_R8P], [2, 2])  ! top face turned by ~30 degrees
y(:,:,2) = y(:,:,2) + reshape([-0.25_R8P, 0.25_R8P, -0.25_R8P, 0.25_R8P], [2, 2])
z(2, 2, 2) = 1.5_R8P                                                                  ! and one top corner lifted
call b%initialize(x, y, z)
call b%geometry(volume=vol, centroid=cen, face_i=fi, face_j=fj, face_k=fk)
fv(:, 1) = -fi(:, 1, 1, 1) ; fv(:, 2) = fi(:, 2, 1, 1) ; fv(:, 3) = -fj(:, 1, 1, 1) ; fv(:, 4) = fj(:, 1, 2, 1)
fv(:, 5) = -fk(:, 1, 1, 1) ; fv(:, 6) = fk(:, 1, 1, 2)
print '(A,F12.9,A,3F9.5)', 'volume ', vol(1, 1, 1), ', centroid ', cen(:, 1, 1, 1)
do f = 1, 6
   print '(A,I0,A,3F9.5)', 'face ', f, ' outward area vector ', fv(:, f) + 0._R8P
enddo
print '(A,ES9.2)', 'largest component of the sum of the six ', maxval(abs(sum(fv, 2)))
gs = [G1, G2] ; m = 0
do c = 1, 2
   do d = 1, 2
      do a = 1, 2
         m = m + 1 ; gp(m) = tri(gs(a), gs(d), gs(c))
      enddo
   enddo
enddo
do f = 1, 6   ! face centres: the bilinear map at (1/2, 1/2)
   select case (f)
   case (1, 2) ; fc(f) = tri(real(f - 1, R8P), 0.5_R8P, 0.5_R8P)
   case (3, 4) ; fc(f) = tri(0.5_R8P, real(f - 3, R8P), 0.5_R8P)
   case default ; fc(f) = tri(0.5_R8P, 0.5_R8P, real(f - 5, R8P))
   endselect
   fe(f) = fc(f) + vector_R8P(0.5_R8P * fv(1, f), 0.5_R8P * fv(2, f), 0.5_R8P * fv(3, f))
enddo

call fig_open(w, 'ref-curvilinear-blocks-cell', 'trilinear cell: 2 x 2 x 2 Gauss points (gold), outward face '// &
              'vectors (teal), sum 0', -58._R8P, 22._R8P, width=640_I4P, height=520_I4P)
do f = 1, 6
   do lb = 0, NS - 1
      do la = 0, NS - 1
         s = real(la, R8P) / NS ; t = real(lb, R8P) / NS
         call fig_quad(w, 'ghost', fp(f, s, t), fp(f, s + 1._R8P / NS, t), fp(f, s + 1._R8P / NS, t + 1._R8P / NS), &
                       fp(f, s, t + 1._R8P / NS))
      enddo
   enddo
enddo
call fig_segments(w, 'body', [tri(0._R8P, 0._R8P, 0._R8P), tri(0._R8P, 1._R8P, 0._R8P), tri(1._R8P, 1._R8P, 0._R8P), &
                              tri(1._R8P, 0._R8P, 0._R8P), tri(0._R8P, 0._R8P, 1._R8P), tri(0._R8P, 1._R8P, 1._R8P), &
                              tri(1._R8P, 1._R8P, 1._R8P), tri(1._R8P, 0._R8P, 1._R8P), tri(0._R8P, 0._R8P, 0._R8P), &
                              tri(1._R8P, 0._R8P, 0._R8P), tri(1._R8P, 1._R8P, 0._R8P), tri(0._R8P, 1._R8P, 0._R8P)], &
                             [tri(0._R8P, 1._R8P, 0._R8P), tri(1._R8P, 1._R8P, 0._R8P), tri(1._R8P, 0._R8P, 0._R8P), &
                              tri(0._R8P, 0._R8P, 0._R8P), tri(0._R8P, 1._R8P, 1._R8P), tri(1._R8P, 1._R8P, 1._R8P), &
                              tri(1._R8P, 0._R8P, 1._R8P), tri(0._R8P, 0._R8P, 1._R8P), tri(0._R8P, 0._R8P, 1._R8P), &
                              tri(1._R8P, 0._R8P, 1._R8P), tri(1._R8P, 1._R8P, 1._R8P), tri(0._R8P, 1._R8P, 1._R8P)])
call fig_points(w, 'accent', gp)
call fig_segments(w, 'probe', fc, fe)
call fig_points(w, 'probe', fe)
do f = 1, 6
   call fig_label(w, 'probe', fe(f) + vector_R8P(0.04_R8P, 0.04_R8P, 0.04_R8P), char(48 + f))
enddo
call fig_close(w)
contains
   function tri(xi, eta, zeta) result(q)
   !< The trilinear map of the cell.
   real(R8P), intent(in) :: xi, eta, zeta !< Local coordinates.
   type(vector_R8P)      :: q             !< Point.
   real(R8P)             :: wi(2), wj(2), wk(2), ww
   integer(I4P)          :: i, j, k

   wi = [1._R8P - xi, xi] ; wj = [1._R8P - eta, eta] ; wk = [1._R8P - zeta, zeta]
   q = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   do k = 1, 2
      do j = 1, 2
         do i = 1, 2
            ww = wi(i) * wj(j) * wk(k)
            q = q + vector_R8P(ww * x(i, j, k), ww * y(i, j, k), ww * z(i, j, k))
         enddo
      enddo
   enddo
   endfunction tri

   function fp(face, s1, s2) result(q)
   !< Point (s1, s2) of a face of the cell.
   integer(I4P), intent(in) :: face   !< Face 1-6.
   real(R8P),    intent(in) :: s1, s2 !< Face coordinates.
   type(vector_R8P)         :: q      !< Point.

   select case (face)
   case (1, 2) ; q = tri(real(face - 1, R8P), s1, s2)
   case (3, 4) ; q = tri(s1, real(face - 3, R8P), s2)
   case default ; q = tri(s1, s2, real(face - 5, R8P))
   endselect
   endfunction fp
endprogram fig_ref_curvilinear_blocks_cell
