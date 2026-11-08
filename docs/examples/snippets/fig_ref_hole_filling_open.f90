!< Figure of the reference page "Hole filling", section "What it does": a unit sphere (marching cubes on 33^3 samples)
!< with three small holes cut in it (the facets within 0.05, 0.07 and 0.085 of three points of the sphere) and its cap
!< z > 0.55 cut off: an opening, not a defect. Left the input with its boundary loops (red); right after fill_holes at
!< its defaults: the three holes are closed by patches (ember), the opening, larger than max_hole_size (0.1 of the box
!< diagonal), stays open.
program fig_ref_hole_filling_open
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ey_R8P

implicit none
integer(I4P), parameter         :: NG = 33
type(surface_stl_object)        :: sphere, cut
type(facet_object), allocatable :: facets(:), kept(:)
type(vector_R8P)                :: centre(3), shift
real(R8P),        parameter     :: RADIUS(3) = [0.05_R8P, 0.07_R8P, 0.085_R8P]
real(R8P)                       :: phi(NG, NG, NG)
character(120)                  :: caption
logical,          allocatable   :: keep(:)
integer(I4P)                    :: u, i, j, k, n, status, filled, skipped, open_edges

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   phi(i, j, k) = norm2(-1.5_R8P + 3._R8P / (NG - 1) * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
call sphere%sanitize
centre(1) = vector_R8P(0.6_R8P, -0.8_R8P, 0._R8P)
centre(2) = vector_R8P(0.25_R8P, -0.75_R8P, -0.61_R8P)
centre(3) = vector_R8P(-0.5_R8P, -0.85_R8P, 0.15_R8P)
n = sphere%get_facets_number()
allocate(keep(n))
do i = 1, n
   associate(f => sphere%facet_at(i))
   keep(i) = f%centroid%z <= 0.55_R8P
   do k = 1, 3
      if (norm2([f%centroid%x - centre(k)%x, f%centroid%y - centre(k)%y, f%centroid%z - centre(k)%z]) < RADIUS(k)) &
         keep(i) = .false.
   enddo
   endassociate
enddo
allocate(kept(count(keep)))
j = 0
do i = 1, n
   if (.not. keep(i)) cycle
   j = j + 1 ; kept(j) = sphere%facet_at(i)
enddo
call cut%adopt_facets(facets=kept)
n = cut%get_facets_number()
sphere = cut
call sphere%fill_holes(holes_filled=filled, holes_skipped=skipped, open_edges=open_edges)
write(caption, '(A,I0,A,I0,A,I0,A)') 'fill_holes: ', filled, ' holes closed by ', sphere%get_facets_number() - n, &
      ' facets (ember), ', skipped, ' loop over max_hole_size left open'
print '(A)', trim(caption)
shift = 2.6_R8P * vector_R8P(0.866_R8P, 0.5_R8P, 0._R8P)
call sphere%translate(delta=shift)
call sphere%analyze
call fig_open(u, 'ref-hole-filling-open', trim(caption), -60._R8P, 25._R8P, 720_I4P, 340_I4P)
call fig_surface(u, 'body', cut)
call fig_boundary(u, 'hole', cut)
call fig_surface(u, 'body', sphere, last=n)
call fig_surface(u, 'patch', sphere, first=n + 1)
call fig_boundary(u, 'hole', sphere)
call fig_label(u, 'accent', vector_R8P(0._R8P, 0._R8P, -1.25_R8P), 'input: 4 boundary loops')
call fig_label(u, 'accent', vector_R8P(0._R8P, 0._R8P, -1.25_R8P) + shift, 'after fill_holes')
call fig_close(u)
endprogram fig_ref_hole_filling_open
