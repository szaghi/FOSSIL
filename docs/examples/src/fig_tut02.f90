!run fig02 fig_tut02
!fig tut02-hole
!fig tut02-patch
!< Figures of tutorial 2: a hole and the patch fill_holes closes it with. The dragon is repaired (sanitize, then
!< fill_holes closes the thin slits it leaves), a round hole is cut on its flank (the facets whose centroid is within
!< 5 % of the diagonal of a point of the surface), and fill_holes(refine, fair) closes it.
program fig_tut02
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_zoom, fig_surface, fig_boundary, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P), parameter            :: DEG = 57.29577951308232_R8P
type(surface_stl_object)        :: surface
type(facet_object), pointer     :: f
type(facet_object), allocatable :: kept(:)
type(vector_R8P)                :: bmin, bmax, p, centre, normal
real(R8P)                       :: diag, radius, d, azimuth, elevation
integer(I4P)                    :: u, i, n, n0, filled

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
call surface%sanitize
call surface%fill_holes(holes_filled=filled)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
diag = norm2([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z]) ; radius = 0.05_R8P * diag
p%x = 0.5_R8P * (bmin%x + bmax%x) ; p%y = bmin%y - 0.5_R8P * diag ; p%z = bmin%z + 0.7_R8P * (bmax%z - bmin%z)
call surface%compute_distance(point=p, distance=d, closest_point=centre, normal=normal)
azimuth = atan2(normal%y, normal%x) * DEG                                        ! look at the hole along its normal
elevation = asin(max(-1._R8P, min(1._R8P, normal%z / norm2([normal%x, normal%y, normal%z])))) * DEG
n = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   if (norm2([f%centroid%x - centre%x, f%centroid%y - centre%y, f%centroid%z - centre%z]) > radius) n = n + 1
enddo
allocate(kept(n)) ; n = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   if (norm2([f%centroid%x - centre%x, f%centroid%y - centre%y, f%centroid%z - centre%z]) <= radius) cycle
   n = n + 1 ; kept(n) = f
enddo
call surface%adopt_facets(facets=kept)
call fig_open(u, 'tut02-hole', 'a hole: its boundary loop, the edges without a neighbour', azimuth, elevation)
call fig_zoom(u, centre, 2.2_R8P * radius)
call fig_surface(u, 'body', surface)
call fig_boundary(u, 'hole', surface)
call fig_close(u)
n0 = surface%get_facets_number()
call surface%fill_holes(holes_filled=filled, refine=.true., fair=.true., max_hole_size=0.5_R8P)
call fig_open(u, 'tut02-patch', 'fill_holes(refine, fair): the patch, refined and faired', azimuth, elevation)
call fig_zoom(u, centre, 2.2_R8P * radius)
call fig_surface(u, 'body', surface, last=n0)
call fig_surface(u, 'patch', surface, first=n0 + 1_I4P)
call fig_close(u)
endprogram fig_tut02
