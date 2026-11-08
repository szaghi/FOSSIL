!< Figures of the reference page "SDF segmentation", on the dragon with the defaults (4 clusters, 30 rays, 120 degrees):
!< - the smoothed shape diameter of each facet, teal (0) to ember (the largest value);
!< - the four labels, teal (1, thinnest) to ember (4, thickest); the unassigned facets (label 0) in gold.
program fig_ref_sdf_segmentation
use fossil,              only : surface_stl_object, SDF_STATUS_OK, SDF_LABEL_UNASSIGNED
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: dragon
type(facet_object), pointer :: fp
integer(I4P), allocatable   :: labels(:)
real(R8P),    allocatable   :: sdf(:)
real(R8P)                   :: smax
character(90)               :: caption
integer(I4P)                :: u, i, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%segment_sdf(facet_labels=labels, sdf=sdf, status=status)
if (status /= SDF_STATUS_OK) error stop 'segment_sdf failed'
smax = maxval(sdf)
write(caption, '(A,F5.2,A)') 'smoothed shape diameter, teal 0 to ember ', smax, ' (the largest)'
call fig_open(u, 'ref-sdf-segmentation-sdf', trim(caption), -60._R8P, 22._R8P)
do i = 1, dragon%get_facets_number()
   fp => dragon%facet_at(i)
   if (labels(i) == SDF_LABEL_UNASSIGNED) then
      call fig_triangle(u, 'foot', fp%vertex(1), fp%vertex(2), fp%vertex(3))
   else
      call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), &
                        value=2._R8P * sdf(i) / smax - 1._R8P)
   endif
enddo
call fig_close(u)
write(caption, '(A,I0,A)') 'labels 1 (teal, thinnest) to 4 (ember, thickest); ', &
                           count(labels == SDF_LABEL_UNASSIGNED), ' unassigned (gold)'
call fig_open(u, 'ref-sdf-segmentation-labels', trim(caption), -60._R8P, 22._R8P)
do i = 1, dragon%get_facets_number()
   fp => dragon%facet_at(i)
   if (labels(i) == SDF_LABEL_UNASSIGNED) then
      call fig_triangle(u, 'foot', fp%vertex(1), fp%vertex(2), fp%vertex(3))
   else
      call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=(2._R8P * labels(i) - 5._R8P) / 3._R8P)
   endif
enddo
call fig_close(u)
endprogram fig_ref_sdf_segmentation
