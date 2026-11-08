!run figck-vertex-curvature fig_ck_vertex_curvature
!fig ck-vertex-curvature
!< Figure of the cookbook recipe "Curvature at every vertex": the marching-cubes unit sphere of ck_vertex_curvature
!< twice, coloured by H - 1 averaged over each facet's vertices: left the cotangent-Laplacian H (mean_curvature),
!< right the cubic-jet H (jet_curvature); teal below 1, ember above, saturating at +-0.3.
program fig_ck_vertex_curvature
use fossil,     only : surface_stl_object, facet_object, extract_isosurface
use fossil_fig, only : fig_open, fig_triangle, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: SCALE = 0.3_R8P
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
real(R8P),          allocatable :: h_cotan(:), h_jet(:), k_jet(:)
type(vector_R8P)                :: shift
real(R8P)                       :: phi(31, 31, 31)
integer(I4P)                    :: i, j, k, status, u

do k = 1, 31 ; do j = 1, 31 ; do i = 1, 31
   phi(i, j, k) = norm2(-1.53_R8P + 0.1_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.53_R8P, -1.53_R8P, -1.53_R8P), &
                        bmax=vector_R8P(1.47_R8P, 1.47_R8P, 1.47_R8P), iso=0._R8P, surface=facets, status=status)
call sphere%adopt_facets(facets=facets)
call sphere%mean_curvature(H=h_cotan)
call sphere%jet_curvature(H=h_jet, K=k_jet)
call fig_open(u, 'ck-vertex-curvature', 'H - 1 per facet: cotangent (left) and cubic jet (right), colour to +-0.3', &
              -90._R8P, 20._R8P, height=340_I4P)
shift = vector_R8P(2.4_R8P, 0._R8P, 0._R8P)
do i = 1, sphere%get_facets_number()
   f => sphere%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1) - shift * 0.5_R8P, f%vertex(2) - shift * 0.5_R8P, &
                     f%vertex(3) - shift * 0.5_R8P, value=(sum(h_cotan(f%vertex_id)) / 3._R8P - 1._R8P) / SCALE)
   call fig_triangle(u, 'field', f%vertex(1) + shift * 0.5_R8P, f%vertex(2) + shift * 0.5_R8P, &
                     f%vertex(3) + shift * 0.5_R8P, value=(sum(h_jet(f%vertex_id)) / 3._R8P - 1._R8P) / SCALE)
enddo
call fig_label(u, 'accent', vector_R8P(-1.5_R8P, 0._R8P, 1.25_R8P), 'mean_curvature')
call fig_label(u, 'accent', vector_R8P(0.9_R8P, 0._R8P, 1.25_R8P), 'jet_curvature')
call fig_close(u)
endprogram fig_ck_vertex_curvature
