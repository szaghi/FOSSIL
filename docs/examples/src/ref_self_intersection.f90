!run ref_self_intersection ref_self_intersection
!fig ref-self-intersection-cubes
!< Run example of the reference page "Self-intersection" and its figure: two unit cubes, the second moved by
!< (0.3, 0.4, 0.45), merged into one surface (two shells passing through each other); its crossing pairs found exactly;
!< then removed by resolve_self_intersections (issue #78). The figure draws the merged surface translucent with the
!< segments the detection reported (ember): the closed curves along which the two shells cross, the seam of the
!< union.
program ref_self_intersection
use fossil,     only : surface_stl_object, intersection_pair_t, SELFINT_CROSSING, SELFINT_SHEETS, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)               :: body, other
type(intersection_pair_t), allocatable :: pairs(:), crossings(:)
character(96)                          :: caption
integer(I4P)                           :: status, u, i
!region merge
call body%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call other%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call other%translate(delta=vector_R8P(0.3_R8P, 0.4_R8P, 0.45_R8P))
call body%merge_solids(other=other)                         ! one surface, two shells passing through each other
call report('merged:   ', body)
!endregion merge
!region detect
call body%find_self_intersections(pairs=crossings, status=status)
print '(A,I0,A,I0,A,I0)', 'pairs:     ', size(crossings), ', crossing ', count(crossings%kind == SELFINT_CROSSING), &
      ', sheets ', count(crossings%kind == SELFINT_SHEETS)
!endregion detect
write(caption, '(A,I0,A)') 'two cubes merged into one surface (translucent): ', size(crossings), &
                           ' crossing segments (ember)'
call fig_open(u, 'ref-self-intersection-cubes', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'ghost', body)
do i = 1, size(crossings)
   call fig_segment(u, 'hole', crossings(i)%p, crossings(i)%q)
enddo
call fig_close(u)
!region resolve
call body%resolve_self_intersections(status=status)         ! the boundary of the region wound at least once
if (status /= BOOL_STATUS_OK) error stop 'resolve_self_intersections failed'
call report('resolved: ', body)
!endregion resolve
contains
   subroutine report(label, s)
   !< Facets, closure, shells, volume and remaining self-intersections of s.
   character(*),             intent(in) :: label
   type(surface_stl_object), intent(in) :: s
   call s%find_self_intersections(pairs=pairs)
   print '(A,I3,A,L1,A,I0,A,F6.4,A,I0)', label, s%get_facets_number(), ' facets, watertight ', s%is_watertight(), &
         ', shells ', s%get_shells_number(), ', volume ', s%get_volume(), ', self-intersecting pairs ', size(pairs)
   endsubroutine report
endprogram ref_self_intersection
