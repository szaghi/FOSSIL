!< Figure of the reference page "Curvilinear structured blocks", section Quality metrics: a 2-D wall block of 24 x 8 cells
!< over a Gaussian bump, vertical grid lines, the wall-normal spacing growing by 1.5 per layer. The same block three times,
!< coloured by the equiangle skewness (top), 1 - orthogonality (middle) and log10 of the aspect ratio (bottom), each
!< from cream (0) to sepia (its maximum, in the label).
program fig_ref_curvilinear_blocks_quality
use fossil,     only : curvilinear_block_2d_t
use fossil_fig, only : fig_open, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter      :: NI = 25, NJ = 9
real(R8P),    parameter      :: GAP = 2.3_R8P
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: u(NI, NJ), v(NI, NJ), yw, s(NJ), vmax(3)
real(R8P), allocatable       :: sk(:,:), orth(:,:), asp(:,:), mj(:,:), m(:,:,:)
integer(I4P)                 :: i, j, k, w, il(2)
character(60)                :: text(3)

s(1) = 0._R8P
do j = 2, NJ
   s(j) = s(j - 1) + 1.5_R8P ** (j - 2)
enddo
s = s / s(NJ)
do i = 1, NI
   yw = 0.6_R8P * exp(-((0.25_R8P * (i - 1) - 3._R8P) / 0.8_R8P) ** 2)
   do j = 1, NJ
      u(i, j) = 0.25_R8P * (i - 1) ; v(i, j) = yw + (1.6_R8P - yw) * s(j)
   enddo
enddo
call b%initialize(u, v)
call b%quality(sk, orth, asp, mj)
allocate(m(NI - 1, NJ - 1, 3))
m(:,:,1) = sk ; m(:,:,2) = 1._R8P - orth ; m(:,:,3) = log10(asp)
vmax = [maxval(m(:,:,1)), maxval(m(:,:,2)), maxval(m(:,:,3))]
il = maxloc(sk)
write(text(1), '(A,F5.3,A,I0,A,I0,A)') 'skewness, max ', maxval(sk), ' at cell (', il(1), ',', il(2), ')'
il = minloc(orth)
write(text(2), '(A,F5.3,A,I0,A,I0,A)') '1 - orthogonality, max ', 1._R8P - minval(orth), ' at cell (', il(1), ',', il(2), ')'
il = maxloc(asp)
write(text(3), '(A,F5.2,A,I0,A,I0,A)') 'log10 aspect ratio; aspect max ', maxval(asp), ' at cell (', il(1), ',', il(2), ')'
do k = 1, 3
   print '(A)', trim(text(k))
enddo
print '(A,F6.3)', 'min J / mean J, lowest ', minval(mj)

call fig_open(w, 'ref-curvilinear-blocks-quality', 'one wall block over a bump coloured by three metrics, '// &
              'cream 0 to sepia max', -90._R8P, 90._R8P, width=760_I4P, height=600_I4P)
do k = 1, 3
   do j = 1, NJ - 1
      do i = 1, NI - 1
         call fig_quad(w, 'frac', p(i, j, k), p(i + 1, j, k), p(i + 1, j + 1, k), p(i, j + 1, k), m(i, j, k) / vmax(k))
      enddo
   enddo
   call fig_label(w, 'accent', vector_R8P(0._R8P, 1.75_R8P - GAP * (k - 1), 0.01_R8P), trim(text(k)))
enddo
call fig_close(w)
contains
   function p(i, j, k) result(q)
   !< Vertex (i, j) of panel k.
   integer(I4P), intent(in) :: i, j, k !< Indices, panel.
   type(vector_R8P)         :: q       !< Point.

   q = vector_R8P(u(i, j), v(i, j) - GAP * (k - 1), 0._R8P)
   endfunction p
endprogram fig_ref_curvilinear_blocks_quality
