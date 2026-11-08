!< Tutorial 2, third run: how far the repair moved the dragon. Certified brackets of the one-sided Hausdorff distance in
!< both directions; the figure marks the facets of the repaired dragon whose centroid is farther than 1e-3 from
!< the original.
program tut02_c_deviation
use fossil_fig,          only : fig_open, fig_surface, fig_triangle, fig_points, fig_close
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: repaired, original
type(facet_object), pointer :: f
real(R8P), allocatable      :: d(:)
real(R8P)                   :: lower, upper
integer(I4P)                :: i, u
character(96)               :: caption

call repaired%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
call repaired%deviation(original, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'repaired -> original ', lower, ' <= h <= ', upper
call original%deviation(repaired, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'original -> repaired ', lower, ' <= h <= ', upper
allocate(d(repaired%get_facets_number()))                                 ! each facet centroid's distance
do i = 1, repaired%get_facets_number()
   f => repaired%facet_at(i)
   d(i) = original%distance(point=f%centroid, is_square_root=.true.)
enddo
print '(A,I0,A,ES10.3)', 'facet centroids farther than 1e-3 from the original: ', count(d > 1.e-3_R8P), &
                         ', the farthest ', maxval(d)
write(caption, '(I0,A,F5.3,A)') count(d > 1.e-3_R8P), ' facets farther than 1e-3 from the original (ember), up to ', &
                                maxval(d)
call fig_open(u, 'tut02-deviation', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', repaired)
do i = 1, repaired%get_facets_number()
   if (d(i) <= 1.e-3_R8P) cycle
   f => repaired%facet_at(i)
   call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
   call fig_points(u, 'outside', [f%centroid])
enddo
call fig_close(u)
endprogram tut02_c_deviation
