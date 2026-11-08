!< Figures of the reference page "Mesh decimation" on a rounded box (half-sides 0.8, 0.45, 0.3, corner radius 0.15)
!< extracted by marching cubes on 40^3 samples: the input next to decimate(target_facets = input/16), where the flat
!< faces are left with a few large triangles and the rounded edges and corners keep the small ones; then the quadric
!< error of the input, stages 1-2 of the pipeline: every facet coloured by the cheapest collapse of its three edges,
!< m^T (Q_u + Q_v) m at the edge midpoint m with Q_v the sum of the plane quadrics p p^T of the facets around v
!< (log scale, cream = 0 on the flat faces, sepia = the most expensive). The heap collapses the cream edges first.
program fig_ref_decimate_quadric
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_triangle, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P

implicit none
integer(I4P), parameter         :: NG = 40
type(surface_stl_object)        :: box, coarse
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
real(R8P)                       :: phi(NG, NG, NG), p(3), q(3), plane(4)
real(R8P), allocatable          :: quadric(:,:,:), cost(:), value(:)
type(vector_R8P)                :: shift
character(120)                  :: caption
real(R8P)                       :: lo, hi, e
integer(I4P)                    :: u, i, j, k, status, nv, a, b

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   p = -1.1_R8P + 2.2_R8P / (NG - 1) * ([i, j, k] - 1)
   q = abs(p) - [0.8_R8P, 0.45_R8P, 0.3_R8P]
   phi(i, j, k) = norm2(max(q, 0._R8P)) + min(maxval(q), 0._R8P) - 0.15_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.1_R8P, -1.1_R8P, -1.1_R8P), &
                        bmax=vector_R8P(1.1_R8P, 1.1_R8P, 1.1_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call box%adopt_facets(facets=facets)
call box%sanitize
coarse = box
call coarse%decimate(target_facets=box%get_facets_number() / 16, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
write(caption, '(A,I0,A,I0,A,F6.3,A,F6.3)') 'input ', box%get_facets_number(), ' facets, decimated ', &
      coarse%get_facets_number(), '; volume ', box%get_volume(), ' -> ', coarse%get_volume()
print '(A)', trim(caption)
shift = 2.1_R8P * ex_R8P
call coarse%translate(delta=shift)
call fig_open(u, 'ref-decimate-box', trim(caption), -90._R8P, 35._R8P, 720_I4P, 300_I4P)
call fig_surface(u, 'body', box)
call fig_surface(u, 'body', coarse)
call fig_label(u, 'accent', vector_R8P(-0.9_R8P, -0.6_R8P, -0.6_R8P), 'marching cubes')
call fig_label(u, 'accent', vector_R8P(-0.9_R8P, -0.6_R8P, -0.6_R8P) + shift, 'decimate(target_facets = input/16)')
call fig_close(u)

nv = 0                                                               ! stage 1: the vertex quadrics
do i = 1, box%get_facets_number()
   f => box%facet_at(i)
   nv = max(nv, maxval(f%vertex_id))
enddo
allocate(quadric(4, 4, nv), cost(box%get_facets_number()), value(box%get_facets_number()))
quadric = 0._R8P
do i = 1, box%get_facets_number()
   f => box%facet_at(i)
   plane = [f%normal%x, f%normal%y, f%normal%z, -f%d]
   do j = 1, 3
      do k = 1, 4
         quadric(:, k, f%vertex_id(j)) = quadric(:, k, f%vertex_id(j)) + plane * plane(k)
      enddo
   enddo
enddo
do i = 1, box%get_facets_number()                                    ! stage 2: the cheapest edge of each facet
   f => box%facet_at(i)
   cost(i) = huge(1._R8P)
   do j = 1, 3
      a = f%vertex_id(j) ; b = f%vertex_id(mod(j, 3) + 1)
      p = 0.5_R8P * ([f%vertex(j)%x, f%vertex(j)%y, f%vertex(j)%z] + &
                     [f%vertex(mod(j, 3) + 1)%x, f%vertex(mod(j, 3) + 1)%y, f%vertex(mod(j, 3) + 1)%z])
      plane = [p, 1._R8P]
      e = dot_product(plane, matmul(quadric(:, :, a) + quadric(:, :, b), plane))
      cost(i) = min(cost(i), max(e, 0._R8P))
   enddo
enddo
lo = -10._R8P ; hi = log10(maxval(cost))
value = (log10(max(cost, 1.e-10_R8P)) - lo) / (hi - lo)
write(caption, '(A,I0,A,I0,A,ES8.1)') 'cheapest collapse: ', count(cost < 1.e-10_R8P), ' of ', &
      box%get_facets_number(), ' facets below 1e-10 (cream), max ', maxval(cost)
print '(A)', trim(caption)
call fig_open(u, 'ref-decimate-quadric', trim(caption), -60._R8P, 30._R8P)
do i = 1, box%get_facets_number()
   f => box%facet_at(i)
   call fig_triangle(u, 'frac', f%vertex(1), f%vertex(2), f%vertex(3), value(i))
enddo
call fig_close(u)
endprogram fig_ref_decimate_quadric
