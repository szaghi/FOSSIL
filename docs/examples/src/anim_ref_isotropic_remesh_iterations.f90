!run anim-ref-isotropic-remesh-iterations anim_ref_isotropic_remesh_iterations
!fig ref-isotropic-remesh-passes
!fig anim-ref-isotropic-remesh-iterations
!< Figures of the reference page "Isotropic remeshing" on the page's example: a unit sphere by marching cubes on 25^3
!< samples, target L = 0.1. The first figure is what the passes of the first iteration act on: the edges longer than
!< 4L/3 (ember, split), shorter than 4L/5 (teal, collapsed) and in the band (grey), and the vertices of valence other
!< than 6 (gold, the flips drive them towards 6). The animation calls isotropic_remesh(iterations=1) five times and
!< draws the same classes after each call: the splits dominate (the input is coarser than L in places, and the
!< remesh does not coarsen), the long edges disappear, the short ones fall to a fraction.
program anim_ref_isotropic_remesh_iterations
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, REM_STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 25
real(R8P),    parameter         :: L = 0.1_R8P
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(NG, NG, NG)
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   phi(i, j, k) = norm2(-1.5_R8P + 3._R8P / (NG - 1) * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)

call fig_open(u, 'ref-isotropic-remesh-passes', '', -60._R8P, 25._R8P)
call fig_zoom(u, vector_R8P(0.45_R8P, -0.75_R8P, 0.45_R8P), 0.7_R8P)
call draw('input', valence=.true.)
call fig_close(u)

call fig_open(u, 'anim-ref-isotropic-remesh-iterations', '', -60._R8P, 25._R8P)
call fig_anim(u, 1.5_R8P)
call fig_frame(u)
call draw('input')
do k = 1, 5
   call sphere%isotropic_remesh(target_length=L, iterations=1_I4P, status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   write(caption, '(A,I0)') 'iteration ', k
   if (k == 5) then
      call fig_frame(u, seconds=3._R8P)
   else
      call fig_frame(u)
   endif
   call draw(trim(caption))
enddo
call fig_close(u)
contains
   subroutine draw(label, valence)
   !< The sphere, its edges by class (each edge once), optionally the vertices of valence other than 6; the caption.
   character(*), intent(in)           :: label   !< State.
   logical,      intent(in), optional :: valence !< Draw the irregular vertices.
   type(vector_R8P), allocatable      :: a(:), b(:), p(:)
   type(vector_R8P)                   :: d
   integer(I4P),     allocatable      :: cls(:), facets_at(:)
   type(facet_object), pointer        :: f
   real(R8P)                          :: length
   integer(I4P)                       :: m, e, nv

   allocate(a(0), b(0), cls(0))
   nv = 0
   do m = 1, sphere%get_facets_number()
      f => sphere%facet_at(m)
      nv = max(nv, maxval(f%vertex_id))
      do e = 1, 3
         if (f%vertex_id(e) > f%vertex_id(mod(e, 3) + 1)) cycle         ! each edge once
         d = f%vertex(mod(e, 3) + 1) - f%vertex(e) ; length = d%normL2()
         a = [a, f%vertex(e)] ; b = [b, f%vertex(mod(e, 3) + 1)]
         if (length > 4._R8P * L / 3._R8P) then
            cls = [cls, 1]
         elseif (length < 0.8_R8P * L) then
            cls = [cls, 2]
         else
            cls = [cls, 3]
         endif
      enddo
   enddo
   write(caption, '(A,A,I0,A,I0,A,I0,A,I0)') label, ': ', sphere%get_facets_number(), ' facets; edges > 4L/3 ', &
         count(cls == 1), ', < 4L/5 ', count(cls == 2), ', in band ', count(cls == 3)
   print '(A)', trim(caption)
   if (.not. present(valence)) write(u, '(A)') 'caption '//trim(caption)
   call fig_surface(u, 'body', sphere)
   call fig_segments(u, 'box', pack(a, cls == 3), pack(b, cls == 3))
   call fig_segments(u, 'outside', pack(a, cls == 1), pack(b, cls == 1))
   call fig_segments(u, 'probe', pack(a, cls == 2), pack(b, cls == 2))
   if (present(valence)) then
      allocate(facets_at(nv), p(nv))
      facets_at = 0
      do m = 1, sphere%get_facets_number()
         f => sphere%facet_at(m)
         do e = 1, 3
            facets_at(f%vertex_id(e)) = facets_at(f%vertex_id(e)) + 1 ; p(f%vertex_id(e)) = f%vertex(e)
         enddo
      enddo
      print '(A,I0,A,I0)', 'vertices ', nv, ', valence other than 6: ', count(facets_at /= 6)
      write(caption, '(A,I0,A,I0,A,I0,A)') 'split ', count(cls == 1), ' (ember), collapse ', count(cls == 2), &
            ' (teal), ', count(facets_at /= 6), ' vertices of valence not 6 (gold)'
      write(u, '(A)') 'caption '//trim(caption)
      call fig_points(u, 'accent', pack(p, facets_at /= 6))
   endif
   endsubroutine draw
endprogram anim_ref_isotropic_remesh_iterations
