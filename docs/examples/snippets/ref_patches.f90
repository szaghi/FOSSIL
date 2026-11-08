!< Run example of the reference page "Named patches" and its figure: a unit box split at mid-height into patches low
!< and high, united with a copy moved by (0.5, 0.5, 0.5) with patches b_low and b_high, then isotropically remeshed;
!< every facet is drawn in the colour of its patch (low teal, high stone, b_low gold, b_high ember): the patch
!< boundaries (the mid-height lines and the intersection curve) stay straight through the boolean and the remesh.
program ref_patches
use fossil, only : surface_stl_object, facet_object, curve_object, slice_surface, BOOL_UNION, BOOL_STATUS_OK, &
                   REM_STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
character(6), parameter     :: STYLE(4) = ['probe ', 'body  ', 'accent', 'patch ']
type(surface_stl_object)    :: a, b
type(facet_object), pointer :: f
type(curve_object)          :: section
character(96)               :: caption
integer(I4P)                :: status, u, i
call split_box(a, 'low', 'high', vector_R8P(0._R8P, 0._R8P, 0._R8P))         ! z < 0.5 and z > 0.5 of a unit box
call split_box(b, 'b_low', 'b_high', vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))   ! the same, moved by (0.5, 0.5, 0.5)
call a%boolean(other=b, op=BOOL_UNION, status=status)                       ! each sub-triangle keeps its patch
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
call areas('union', a)
call a%isotropic_remesh(target_length=0.15_R8P, iterations=3_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
call areas('remeshed', a)                                                     ! patch boundaries stay where they are
call slice_surface(surface=a, normal_axis=3_I4P, level=0.75_R8P, curve=section)   ! the section z = 0.75
call lengths(section)                                                         ! each segment has its facet's patch
write(caption, '(A,I0,A)') 'union of two split boxes, remeshed: ', a%get_facets_number(), &
                           ' facets, coloured by patch'
call fig_open(u, 'ref-patches-remeshed', trim(caption), -60._R8P, 25._R8P)
do i = 1, a%get_facets_number()
   f => a%facet_at(i)
   call fig_triangle(u, trim(STYLE(a%facet_patch(i) - 1)), f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_close(u)
contains
   subroutine split_box(box, lower, upper, shift)
   !< The unit box moved by shift, its facets below its mid-height in patch `lower`, the others in `upper`.
   type(surface_stl_object), intent(inout) :: box
   character(*),             intent(in)    :: lower, upper
   type(vector_R8P),         intent(in)    :: shift
   type(facet_object), pointer             :: f
   integer(I4P), allocatable               :: low(:), high(:)
   integer(I4P)                            :: i
   call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call box%translate(delta=shift) ; call box%analyze
   allocate(low(0), high(0))
   do i = 1, box%get_facets_number()
      f => box%facet_at(i)
      if (f%centroid%z < shift%z + 0.5_R8P) then ; low = [low, i] ; else ; high = [high, i] ; endif
   enddo
   call box%set_patch(facets=low, name=lower)
   call box%set_patch(facets=high, name=upper)
   endsubroutine split_box

   subroutine areas(label, s)
   !< Facets and area of every non-empty patch of s.
   character(*),             intent(in) :: label
   type(surface_stl_object), intent(in) :: s
   type(facet_object), pointer          :: f
   real(R8P)                            :: area
   integer(I4P)                         :: p, i, n
   print '(A,A,I0,A)', label, ': ', s%get_facets_number(), ' facets'
   do p = 1, s%get_patches_number()
      area = 0._R8P ; n = 0
      do i = 1, s%get_facets_number()
         if (s%facet_patch(i) /= p) cycle
         f => s%facet_at(i)
         area = area + f%area() ; n = n + 1
      enddo
      if (n > 0) print '(A,A8,A,I4,A,F11.8)', '  patch ', s%patch_name(p), ': facets ', n, ', area ', area
   enddo
   endsubroutine areas

   subroutine lengths(c)
   !< Length of every non-empty patch of the curve c.
   type(curve_object), intent(in) :: c
   real(R8P)                      :: length
   integer(I4P)                   :: p, i
   print '(A,I0,A)', 'section z = 0.75: ', c%get_vertices_number(), ' segments'
   do p = 1, c%get_patches_number()
      length = 0._R8P
      do i = 1, c%get_vertices_number()
         if (c%segment_patch(i) == p) length = length + norm2(c%get_vertex(c%get_next(i)) - c%get_vertex(i))
      enddo
      if (length > 0._R8P) print '(A,A8,A,F11.8)', '  patch ', c%patch_name(p), ': length ', length
   enddo
   endsubroutine lengths
endprogram ref_patches
