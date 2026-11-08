!run figref-isotropic-remesh-torus fig_ref_isotropic_remesh_torus
!fig ref-isotropic-remesh-torus
!< Figure of the reference page "Isotropic remeshing", section "What it does": a torus (radii 0.7 and 0.3) by marching
!< cubes on 24^3 samples of the box [-1.1, 1.1]^2 x [-0.4, 0.4], cells 2.75 times flatter than wide, so its triangles
!< are stretched; then isotropic_remesh(target_length = 0.08, iterations = 3). Every facet is coloured by its longest
!< edge over L (cream 0, sepia 2 and above); the captions give the edge-length range of each mesh.
program fig_ref_isotropic_remesh_torus
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, REM_STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 24
real(R8P),    parameter         :: L = 0.08_R8P
type(surface_stl_object)        :: torus, remeshed
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(NG, NG, NG), p(3)
type(vector_R8P)                :: shift
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   p = [-1.1_R8P + 2.2_R8P / (NG - 1) * (i - 1), -1.1_R8P + 2.2_R8P / (NG - 1) * (j - 1), &
        -0.4_R8P + 0.8_R8P / (NG - 1) * (k - 1)]
   phi(i, j, k) = hypot(hypot(p(1), p(2)) - 0.7_R8P, p(3)) - 0.3_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.1_R8P, -1.1_R8P, -0.4_R8P), &
                        bmax=vector_R8P(1.1_R8P, 1.1_R8P, 0.4_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call torus%adopt_facets(facets=facets)
call torus%sanitize
remeshed = torus
call remeshed%isotropic_remesh(target_length=L, iterations=3_I4P, status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
write(caption, '(A,I0,A,F5.3,A,F5.3,A,I0,A,F5.3,A,F5.3)') 'input ', torus%get_facets_number(), ' facets, edges ', &
      torus%smallest_edge_len(), '-', torus%largest_edge_len(), '; remeshed ', remeshed%get_facets_number(), ', ', &
      remeshed%smallest_edge_len(), '-', remeshed%largest_edge_len()
print '(A)', trim(caption)
shift = vector_R8P(2.4_R8P * 0.866_R8P, 2.4_R8P * 0.5_R8P, 0._R8P)
call fig_open(u, 'ref-isotropic-remesh-torus', trim(caption), -60._R8P, 35._R8P, 720_I4P, 320_I4P)
call draw(torus, vector_R8P(0._R8P, 0._R8P, 0._R8P))
call draw(remeshed, shift)
call fig_label(u, 'accent', vector_R8P(0._R8P, -1.1_R8P, -0.5_R8P), 'marching cubes')
call fig_label(u, 'accent', vector_R8P(0._R8P, -1.1_R8P, -0.5_R8P) + shift, 'isotropic_remesh(L = 0.08)')
call fig_close(u)
contains
   subroutine draw(s, delta)
   !< The facets of s moved by delta, coloured by their longest edge over L.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   type(vector_R8P),         intent(in) :: delta !< Shift.
   type(facet_object), pointer          :: f
   type(vector_R8P)                     :: d(3)
   integer(I4P)                         :: m

   do m = 1, s%get_facets_number()
      f => s%facet_at(m)
      d(1) = f%vertex(2) - f%vertex(1) ; d(2) = f%vertex(3) - f%vertex(2) ; d(3) = f%vertex(1) - f%vertex(3)
      call fig_triangle(u, 'frac', f%vertex(1) + delta, f%vertex(2) + delta, f%vertex(3) + delta, &
                        min(1._R8P, max(d(1)%normL2(), d(2)%normL2(), d(3)%normL2()) / (2._R8P * L)))
   enddo
   endsubroutine draw
endprogram fig_ref_isotropic_remesh_torus
