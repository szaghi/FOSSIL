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
