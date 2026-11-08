!< Reference page "SDF segmentation": a plate with a thickness step (0.1, then 0.5) split into two labels, then the
!< dragon into four. The plate is built here, watertight and outward-oriented; the figure colours it by label.
program ref_sdf_segmentation
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use fossil,              only : surface_stl_object, SDF_STATUS_OK, SDF_LABEL_UNASSIGNED
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: NX = 16, NY = 4             ! cells of the plate's top and bottom faces
type(surface_stl_object)        :: plate, dragon
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
integer(I4P),       allocatable :: labels(:)
real(R8P),          allocatable :: sdf(:), xc(:)
integer(I4P)                    :: k, m, status, u

call build_plate(f)                                    ! [0, 4] x [0, 1], thickness 0.1 for x < 2, 0.5 for x > 2.25
call plate%adopt_facets(facets=f)
call plate%analyze
call plate%segment_sdf(facet_labels=labels, sdf=sdf, num_clusters=2_I4P, status=status)
if (status /= SDF_STATUS_OK) error stop 'segment_sdf failed'
allocate(xc(plate%get_facets_number()))
do k = 1, plate%get_facets_number()                    ! facet centroids along the plate
   fp => plate%facet_at(k)
   xc(k) = (fp%vertex(1)%x + fp%vertex(2)%x + fp%vertex(3)%x) / 3._R8P
enddo
print '(A,I0,A,L1,A,I0)', 'plate: ', plate%get_facets_number(), ' facets, watertight ', plate%is_watertight(), &
                          ', unassigned ', count(labels == SDF_LABEL_UNASSIGNED)
do k = 1, 2
   print '(A,I0,A,I3,A,F6.3,A,F6.3,A,F6.3,A,F6.3)', '  label ', k, ': ', count(labels == k), ' facets, sdf ', &
         minval(sdf, labels == k), ' to ', maxval(sdf, labels == k), ', centroids x ', minval(xc, labels == k), &
         ' to ', maxval(xc, labels == k)
enddo
call fig_open(u, 'ref-sdf-segmentation-plate', 'plate with a thickness step, num_clusters=2: label 1 teal, label 2 ember', &
              -60._R8P, 30._R8P)
do k = 1, plate%get_facets_number()
   fp => plate%facet_at(k)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=merge(-0.8_R8P, 0.8_R8P, labels(k) == 1))
enddo
call fig_close(u)
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%segment_sdf(facet_labels=labels, sdf=sdf, status=status)    ! defaults: 4 clusters, 30 rays, 120 degrees
if (status /= SDF_STATUS_OK) error stop 'segment_sdf failed'
print '(A,I0,A,I0)', 'dragon: ', size(labels), ' facets, unassigned ', count(labels == SDF_LABEL_UNASSIGNED)
do k = 1, 4
   print '(A,I0,A,I5,A,F6.3,A,F6.3)', '  label ', k, ': ', count(labels == k), ' facets, sdf ', &
         minval(sdf, labels == k), ' to ', maxval(sdf, labels == k)
enddo
contains
   subroutine build_plate(ff)
   !< The plate [0, 4] x [0, 1], z in [-t/2, t/2]: t = 0.1 up to x = 2, 0.5 from x = 2.25, linear in between.
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: x(0:NX), y(0:NY), t(0:NX), zb(0:NX), zt(0:NX)
   integer(I4P)                                 :: i, j

   x = [(4._R8P * i / NX, i = 0, NX)] ; y = [(1._R8P * j / NY, j = 0, NY)]
   t = merge(0.1_R8P, 0.5_R8P, x < 2.1_R8P) ; zb = -0.5_R8P * t ; zt = 0.5_R8P * t
   allocate(ff(4 * NX * NY + 4 * NX + 4 * NY))
   m = 0
   do j = 0, NY - 1
      do i = 0, NX - 1
         call quad(ff, [0._R8P, 0._R8P, 1._R8P], [x(i), y(j), zt(i)], [x(i + 1), y(j), zt(i + 1)], &
                   [x(i + 1), y(j + 1), zt(i + 1)], [x(i), y(j + 1), zt(i)])
         call quad(ff, [0._R8P, 0._R8P, -1._R8P], [x(i), y(j), zb(i)], [x(i + 1), y(j), zb(i + 1)], &
                   [x(i + 1), y(j + 1), zb(i + 1)], [x(i), y(j + 1), zb(i)])
      enddo
   enddo
   do i = 0, NX - 1
      call quad(ff, [0._R8P, -1._R8P, 0._R8P], [x(i), y(0), zb(i)], [x(i + 1), y(0), zb(i + 1)], &
                [x(i + 1), y(0), zt(i + 1)], [x(i), y(0), zt(i)])
      call quad(ff, [0._R8P, 1._R8P, 0._R8P], [x(i), y(NY), zb(i)], [x(i + 1), y(NY), zb(i + 1)], &
                [x(i + 1), y(NY), zt(i + 1)], [x(i), y(NY), zt(i)])
   enddo
   do j = 0, NY - 1
      call quad(ff, [-1._R8P, 0._R8P, 0._R8P], [x(0), y(j), zb(0)], [x(0), y(j + 1), zb(0)], &
                [x(0), y(j + 1), zt(0)], [x(0), y(j), zt(0)])
      call quad(ff, [1._R8P, 0._R8P, 0._R8P], [x(NX), y(j), zb(NX)], [x(NX), y(j + 1), zb(NX)], &
                [x(NX), y(j + 1), zt(NX)], [x(NX), y(j), zt(NX)])
   enddo
   endsubroutine build_plate

   subroutine quad(ff, out, p1, p2, p3, p4)
   !< Two facets of the quad p1 p2 p3 p4 at ff(m + 1:m + 2), oriented along `out`.
   type(facet_object), intent(inout) :: ff(:)
   real(R8P),          intent(in)    :: out(3), p1(3), p2(3), p3(3), p4(3)

   m = m + 1 ; call tri(ff(m), out, p1, p2, p3)
   m = m + 1 ; call tri(ff(m), out, p1, p3, p4)
   endsubroutine quad

   subroutine tri(ft, out, p1, p2, p3)
   !< A facet p1 p2 p3, its vertices ordered so that its normal points along `out`.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: out(3), p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), n(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(n, out) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri
endprogram ref_sdf_segmentation
