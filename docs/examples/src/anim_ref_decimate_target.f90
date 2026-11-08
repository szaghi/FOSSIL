!run anim-ref-decimate-target anim_ref_decimate_target
!fig anim-ref-decimate-target
!< Animation of the reference page "Mesh decimation": the rounded box of the page's figures (marching cubes on 40^3
!< samples) decimated from its 6460 facets to targets halved frame by frame, each frame a fresh decimate of the input,
!< seen from the same camera: the flat faces empty first, the cylindrical edges are thinned along their axis, the
!< corners keep their triangles longest; the caption gives the facets and the volume change.
program anim_ref_decimate_target
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 40
type(surface_stl_object)        :: box, coarse
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(NG, NG, NG), p(3), q(3)
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status, target

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
call fig_open(u, 'anim-ref-decimate-target', '', -60._R8P, 30._R8P)
call fig_anim(u, 1.2_R8P)
do k = 0, 6
   target = box%get_facets_number() / 2**k
   coarse = box
   if (k > 0) then
      call coarse%decimate(target_facets=target, status=status)
      if (status /= DEC_STATUS_OK) error stop 'decimate failed'
   endif
   write(caption, '(A,I0,A,I0,A,ES9.2)') 'target_facets ', target, ': ', coarse%get_facets_number(), &
         ' facets, volume change ', coarse%get_volume() / box%get_volume() - 1._R8P
   print '(A,A,L1,L2)', trim(caption), '; watertight, manifold ', coarse%is_watertight(), coarse%is_manifold()
   if (k == 6) then
      call fig_frame(u, caption=trim(caption), seconds=3._R8P)
   else
      call fig_frame(u, caption=trim(caption))
   endif
   call fig_surface(u, 'body', coarse)
enddo
call fig_close(u)
endprogram anim_ref_decimate_target
