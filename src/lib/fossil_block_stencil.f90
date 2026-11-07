!< FOSSIL, Lagrange interpolation stencils on curvilinear blocks: 1-D kernels (issue #71).
module fossil_block_stencil
!< FOSSIL, Lagrange interpolation stencils on curvilinear blocks: 1-D kernels (issue #71).
!<
!< A stencil of width q in one direction takes the vertices start..start+q-1 around the index coordinate s = i + xi
!< of a located point (vertex m sits at s = m): centred (q even: i - q/2 + 1 .. i + q/2; q odd: around the nearest
!< vertex), shifted one-sided at a physical boundary, or reduced in width. The weights are the 1-D Lagrange basis at s,
!< so they reproduce the polynomials of degree q - 1 in s exactly; the tensor product does the same in 2-D and 3-D.
use penf, only : I4P, R8P

implicit none
private
public :: lagrange_weights, centred_start
public :: STENCIL_CENTRED, STENCIL_SHIFTED, STENCIL_REDUCED, STENCIL_FAILED, STENCIL_MAX

integer(I4P), parameter :: STENCIL_CENTRED = 1_I4P !< Centred in every direction (crossing records if needed).
integer(I4P), parameter :: STENCIL_SHIFTED = 2_I4P !< Shifted one-sided at a physical boundary.
integer(I4P), parameter :: STENCIL_REDUCED = 3_I4P !< Width reduced at a physical boundary.
integer(I4P), parameter :: STENCIL_FAILED  = 4_I4P !< No stencil (point not located, or the block is too small).
integer(I4P), parameter :: STENCIL_MAX     = 8_I4P !< Largest width.

contains
   pure subroutine lagrange_weights(q, start, s, w, dw)
   !< Lagrange weights (and their derivatives) at s of the nodes start..start+q-1.
   integer(I4P), intent(in)            :: q, start !< Width, first node.
   real(R8P),    intent(in)            :: s        !< Index coordinate.
   real(R8P),    intent(out)           :: w(:)     !< Weights [q].
   real(R8P),    intent(out), optional :: dw(:)    !< Derivatives [q].
   real(R8P)                           :: t(q), num, den, acc, term !< Work.
   integer(I4P)                        :: m, l, r  !< Counters.

   do m = 1, q
      t(m) = s - real(start + m - 1_I4P, R8P)
   enddo
   do m = 1, q
      num = 1._R8P ; den = 1._R8P
      do l = 1, q
         if (l == m) cycle
         num = num * t(l) ; den = den * real(m - l, R8P)
      enddo
      w(m) = num / den
      if (present(dw)) then
         acc = 0._R8P
         do r = 1, q
            if (r == m) cycle
            term = 1._R8P
            do l = 1, q
               if (l == m .or. l == r) cycle
               term = term * t(l)
            enddo
            acc = acc + term
         enddo
         dw(m) = acc / den
      endif
   enddo
   endsubroutine lagrange_weights

   pure integer(I4P) function centred_start(i, xi, q) result(start)
   !< First node of the centred stencil of width q around s = i + xi.
   integer(I4P), intent(in) :: i  !< Cell.
   real(R8P),    intent(in) :: xi !< Local coordinate.
   integer(I4P), intent(in) :: q  !< Width.

   if (mod(q, 2_I4P) == 0_I4P) then
      start = i - q / 2_I4P + 1_I4P
   else
      start = i + merge(1_I4P, 0_I4P, xi >= 0.5_R8P) - (q - 1_I4P) / 2_I4P
   endif
   endfunction centred_start
endmodule fossil_block_stencil
