program ref_alpha_wrap
use fossil, only : surface_stl_object, facet_object, AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)                :: source, holed, wrapped
type(facet_object), pointer             :: facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
integer(I4P)                            :: status, i
call source%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube, 8 facets per face
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1))                                  ! drop facet 1: a hole on the face x = 0
do i = 1, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1)
enddo
call holed%adopt_facets(facets=holed_facets)
print '(A,I0,A,L1,A,L1)', 'input:   ', holed%get_facets_number(), ' facets, watertight ', holed%is_watertight(), &
                          ', volume ', holed%is_volume()
call holed%alpha_wrap(alpha=0.6_R8P, offset=0.02_R8P, wrapped=wrapped, status=status)   ! ball wider than the hole
call report('alpha 0.60')
call holed%alpha_wrap(alpha=0.15_R8P, offset=0.033_R8P, wrapped=wrapped, status=status) ! ball through the hole
call report('alpha 0.15')
call holed%alpha_wrap(alpha=-1._R8P, wrapped=wrapped, status=status)
print '(A,L1)', 'alpha < 0:   AWRAP_STATUS_BAD_INPUT  ', status == AWRAP_STATUS_BAD_INPUT
call holed%alpha_wrap(alpha=2._R8P, wrapped=wrapped, status=status)
print '(A,L1)', 'alpha > box: AWRAP_STATUS_DEGENERATE ', status == AWRAP_STATUS_DEGENERATE
contains
   subroutine report(label)
   character(*), intent(in) :: label
   print '(2A,L1,A,I0,A,L1,A,L1,A,F7.4)', label, ': ok ', status == AWRAP_STATUS_OK, ', ', wrapped%get_facets_number(), &
         ' facets, watertight ', wrapped%is_watertight(), ', manifold ', wrapped%is_manifold(), ', volume ', &
         wrapped%get_volume()
   endsubroutine report
endprogram ref_alpha_wrap
