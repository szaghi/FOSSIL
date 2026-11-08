!< Figure of the reference page "Body-fitted grid support": the projection target of the NACA 0012 wing with its rear
!< 30 % of chord given a patch of its own (a flap): the constraint curves (sharp edges and the patch boundaries) and the
!< corners where they meet.
program fig_ref_body_fitted_target
use fossil,              only : surface_stl_object, projection_target_t
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_triangle, fig_segments, fig_points, fig_label, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P
implicit none
type(surface_stl_object)      :: wing
type(projection_target_t)     :: target
type(facet_object), pointer   :: f
type(vector_R8P), allocatable :: a(:), b(:), c(:)
integer(I4P), allocatable     :: flap(:)
integer(I4P)                  :: i, k, n, u

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
allocate(flap(0))
do i = 1, wing%get_facets_number()
   f => wing%facet_at(i)
   if (min(f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x) >= 0.45_R8P) flap = [flap, i]
enddo
call wing%set_patch(flap, 'flap')
call wing%projection_target(target)
n = size(target%seg, 2)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'facets ', wing%get_facets_number(), ', patches ', target%np, &
                                    ', constraint curves ', target%nc, ', segments ', n, ', corners ', target%nk
allocate(a(10 * n), b(10 * n), c(target%nk))
do k = 1, n ! each segment drawn in 10 pieces, so that the depth sort keeps them over the long spanwise facets
   do i = 1, 10
      a(10 * (k - 1) + i) = vector_R8P(target%seg(1, k), target%seg(2, k), target%seg(3, k)) * ((10 - i + 1) / 10._R8P) &
                          + vector_R8P(target%seg(4, k), target%seg(5, k), target%seg(6, k)) * ((i - 1) / 10._R8P)
      b(10 * (k - 1) + i) = vector_R8P(target%seg(1, k), target%seg(2, k), target%seg(3, k)) * ((10 - i) / 10._R8P) &
                          + vector_R8P(target%seg(4, k), target%seg(5, k), target%seg(6, k)) * (i / 10._R8P)
   enddo
enddo
do k = 1, target%nk
   c(k) = vector_R8P(target%corner(1, k), target%corner(2, k), target%corner(3, k))
enddo
call fig_open(u, 'ref-body-fitted-target', 'two patches: constraint curves (gold), corners (teal)', 70._R8P, 25._R8P)
do i = 1, wing%get_facets_number()
   f => wing%facet_at(i)
   call fig_triangle(u, merge('patch', 'body ', wing%facet_patch(i) == 2), f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_segments(u, 'accent', a, b)
call fig_points(u, 'probe', c)
call fig_label(u, 'accent', vector_R8P(0.75_R8P, 0._R8P, 1.1_R8P), 'trailing-edge lines')
call fig_label(u, 'accent', vector_R8P(0.565_R8P, 0.03_R8P, 0.6_R8P), 'patch boundary')
call fig_label(u, 'accent', vector_R8P(0.7_R8P, 0.02_R8P, 0.3_R8P), 'flap')
call fig_close(u)
endprogram fig_ref_body_fitted_target
