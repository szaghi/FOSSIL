!run anim-ref-alpha-wrap-shrink anim_ref_alpha_wrap_shrink
!fig ref-alpha-wrap-dragon
!fig anim-ref-alpha-wrap-shrink
!< Figures of the reference page "Alpha wrapping" on the dragon fixture as loaded (not closed: 61 non-manifold edges):
!< the input with its non-manifold edges (gold) next to its wrap at alpha = diagonal/14 (default offset alpha/30);
!< then an animation: the wrap, drawn translucent over the dragon, at alpha = diagonal/4, /6, /9, /14, /20: as the ball
!< shrinks the wrap enters the gaps between the claws, the tail and the body and settles onto the input.
program anim_ref_alpha_wrap_shrink
use fossil,                    only : surface_stl_object, AWRAP_STATUS_OK, FEATURE_EDGE_NONMANIFOLD
use fossil_fig,                only : fig_open, fig_surface, fig_segments, fig_label, fig_anim, fig_frame, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P

implicit none
real(R8P),    parameter           :: DIVISOR(5) = [4._R8P, 6._R8P, 9._R8P, 14._R8P, 20._R8P]
type(surface_stl_object)          :: dragon, wrapped, kept
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:)
type(vector_R8P)                  :: e, shift
integer(I4P), allocatable         :: edges(:,:), kinds(:)
real(R8P)                         :: diagonal
character(120)                    :: caption
integer(I4P)                      :: u, i, k, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
e = dragon%get_bmax() - dragon%get_bmin()
diagonal = e%normL2()
call dragon%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
pool => dragon%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
call fig_open(u, 'anim-ref-alpha-wrap-shrink', '', -60._R8P, 22._R8P)
call fig_anim(u, 1.6_R8P)
call fig_surface(u, 'body', dragon)
do k = 1, size(DIVISOR)
   call dragon%alpha_wrap(alpha=diagonal / DIVISOR(k), wrapped=wrapped, status=status)
   if (status /= AWRAP_STATUS_OK) error stop 'alpha_wrap failed'
   write(caption, '(A,I0,A,F5.2,A,I0,A,F7.1)') 'alpha = diagonal/', nint(DIVISOR(k)), ' = ', diagonal / DIVISOR(k), &
         ': ', wrapped%get_facets_number(), ' facets, volume ', wrapped%get_volume()
   print '(A)', trim(caption)
   if (k == size(DIVISOR)) then
      call fig_frame(u, caption=trim(caption), seconds=3._R8P)
   else
      call fig_frame(u, caption=trim(caption))
   endif
   call fig_surface(u, 'ghost', wrapped)
   if (k == 4) kept = wrapped                                        ! diagonal/14, for the side-by-side figure
enddo
call fig_close(u)

write(caption, '(A,I0,A,I0,A,L1,A,L1)') 'input: ', count(kinds == FEATURE_EDGE_NONMANIFOLD), &
      ' non-manifold edges; wrap: ', kept%get_facets_number(), ' facets, watertight ', kept%is_watertight(), &
      ', manifold ', kept%is_manifold()
print '(A)', trim(caption)
shift = (1.15_R8P * e%x) * ex_R8P
call kept%translate(delta=shift)
call fig_open(u, 'ref-alpha-wrap-dragon', trim(caption), -90._R8P, 12._R8P, 720_I4P, 330_I4P)
call fig_surface(u, 'body', dragon)
call fig_segments(u, 'accent', pack(a, kinds == FEATURE_EDGE_NONMANIFOLD), pack(b, kinds == FEATURE_EDGE_NONMANIFOLD))
call fig_surface(u, 'body', kept)
call fig_label(u, 'accent', dragon%get_bmin(), 'dragon.stl')
call fig_label(u, 'accent', dragon%get_bmin() + shift, 'alpha_wrap(alpha = diagonal/14)')
call fig_close(u)
endprogram anim_ref_alpha_wrap_shrink
