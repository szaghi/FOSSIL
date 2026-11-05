!run figref-feature-edges fig_ref_feature_edges
!fig ref-feature-edges
!< Figure of the reference page "Feature edges": a machined part built with boolean (a 2 x 1 x 0.5 block, united with a
!< raised boss, a 24-sided round hole drilled through it), its feature edges at the default 30 degrees (ember) and its
!< corners (gold). The hole's rims are two closed polylines without corners: the hole's facets meet at
!< 15 degrees, below the threshold.
program fig_ref_feature_edges
use fossil,                    only : surface_stl_object, facet_object, STATUS_OK
use fossil_boolean,            only : BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_fig,                only : fig_open, fig_surface, fig_segments, fig_points, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter              :: PI = acos(-1._R8P)
type(surface_stl_object)          :: part, boss, hole
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:), c(:)
integer(I4P), allocatable         :: edges(:,:), kinds(:), corners(:), ps(:), pv(:)
character(96)                     :: caption
integer(I4P)                      :: u, i, status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P])
call box(boss, [1.4_R8P, 0.2_R8P, 0._R8P], [1.9_R8P, 0.8_R8P, 0.9_R8P])
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call cylinder(hole, [0.6_R8P, 0.5_R8P], 0.25_R8P, -0.1_R8P, 0.6_R8P, 24_I4P)
call part%boolean(hole, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
call part%feature_edges(edges=edges, kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0,A,I0,A,L1)', 'part: ', part%get_facets_number(), ' facets, ', size(kinds), ' feature edges, ', &
                                    size(corners), ' corners, ', size(ps) - 1, ' polylines; watertight ', part%is_watertight()
pool => part%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)), c(size(corners)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
do i = 1, size(corners)
   c(i) = pool%coord(corners(i))
enddo
write(caption, '(I0,A,I0,A,I0,A)') size(kinds), ' feature edges at 30 degrees, ', size(corners), ' corners, ', &
                                   size(ps) - 1, ' polylines'
call fig_open(u, 'ref-feature-edges', trim(caption), -55._R8P, 30._R8P)
call fig_surface(u, 'body', part)
call fig_segments(u, 'hole', a, b)
call fig_points(u, 'accent', c)
call fig_close(u)
contains
   subroutine box(s, lo, hi)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)

   call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call s%resize(factor=(hi(1) - lo(1)) * ex_R8P + (hi(2) - lo(2)) * ey_R8P + (hi(3) - lo(3)) * ez_R8P)
   call s%translate(delta=lo(1) * ex_R8P + lo(2) * ey_R8P + lo(3) * ez_R8P)
   call s%analyze
   endsubroutine box

   subroutine cylinder(s, centre, r, z0, z1, m)
   !< A regular m-gon of radius r about centre, extruded over [z0, z1], caps fanned, outward.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: centre(2), r, z0, z1
   integer(I4P),             intent(in)    :: m
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: x(m), y(m)
   integer(I4P)                            :: k, j, n

   x = [(centre(1) + r * cos(2 * PI * (k - 1) / m), k = 1, m)] ; y = [(centre(2) + r * sin(2 * PI * (k - 1) / m), k = 1, m)]
   allocate(f(2 * (m - 2) + 2 * m))
   n = 0
   do k = 2, m - 1
      n = n + 1 ; call tri(f(n), [x(1), y(1), z1], [x(k), y(k), z1], [x(k+1), y(k+1), z1])
      n = n + 1 ; call tri(f(n), [x(1), y(1), z0], [x(k+1), y(k+1), z0], [x(k), y(k), z0])
   enddo
   do k = 1, m
      j = mod(k, m) + 1
      n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z0], [x(j), y(j), z1])
      n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z1], [x(k), y(k), z1])
   enddo
   call s%destroy
   call s%adopt_facets(facets=f)
   endsubroutine cylinder

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_metrix
   endsubroutine tri
endprogram fig_ref_feature_edges
