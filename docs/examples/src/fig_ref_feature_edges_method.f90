!run figref-feature-edges-method fig_ref_feature_edges_method
!fig ref-feature-edges-method
!< Figure of the reference page "Feature edges", Method: a D-shaped block (a rectangle closed by a half-disc of 8
!< chords, extruded), its feature edges from feature_edges in ember and its corners in gold. The walls of the rounded
!< end meet at 22.5 degrees and the flat sides meet the round at 11.25: below the 30-degree threshold, not sharp (their
!< edges drawn faint). The rims meet the walls at 90 degrees: sharp. The rim vertices along the round have two feature
!< edges turning by 22.5 degrees, below the corner angle: no corners there; the four back vertices have three.
program fig_ref_feature_edges_method
use fossil,                    only : surface_stl_object, facet_object, STATUS_OK
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_fig,                only : fig_open, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
real(R8P),    parameter           :: PI = acos(-1._R8P)
integer(I4P), parameter           :: NARC = 8, NP = NARC + 3      ! profile vertices
type(surface_stl_object)          :: part
type(facet_object), allocatable   :: f(:)
type(vertex_pool_object), pointer :: pool
type(vector_R8P)                  :: p(NP), cen, wa(NP), wb(NP)
type(vector_R8P), allocatable     :: a(:), b(:), c(:)
integer(I4P),     allocatable     :: edges(:,:), kinds(:), corners(:), ps(:), pv(:)
real(R8P)                         :: t
integer(I4P)                      :: i, j, m, u, status
character(90)                     :: caption

p(1) = vector_R8P(-1._R8P, -0.5_R8P, 0._R8P)                        ! the profile, counter-clockwise
do i = 0, NARC
   t = -0.5_R8P * PI + PI * i / NARC
   p(i + 2) = vector_R8P(0.5_R8P * cos(t), 0.5_R8P * sin(t), 0._R8P)
enddo
p(NP) = vector_R8P(-1._R8P, 0.5_R8P, 0._R8P)
cen = vector_R8P(-0.3_R8P, 0._R8P, 0._R8P)
allocate(f(4 * NP))
m = 0
do i = 1, NP
   j = mod(i, NP) + 1
   call tri(cen, p(j), p(i))                                         ! bottom cap, normal -z
   call tri(cen + up(), p(i) + up(), p(j) + up())                    ! top cap, normal +z
   call tri(p(i), p(j), p(j) + up())                                 ! wall, outward
   call tri(p(i), p(j) + up(), p(i) + up())
enddo
call part%adopt_facets(facets=f) ; call part%analyze
call part%feature_edges(edges=edges, kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0,A,L1)', 'D-block: ', size(kinds), ' feature edges, ', size(corners), ' corners, ', &
                               size(ps) - 1, ' polylines; watertight ', part%is_watertight()
print '(A,F6.2,A,F6.2,A)', 'wall dihedrals: round ', 180._R8P / NARC, ' deg, flat-to-round ', 90._R8P / NARC, &
                           ' deg, rim 90 deg; threshold 30'
pool => part%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)), c(size(corners)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
do i = 1, size(corners)
   c(i) = pool%coord(corners(i))
enddo
do i = 1, NP                                                         ! the vertical wall edges, for reference
   wa(i) = p(i) ; wb(i) = p(i) + up()
enddo
write(caption, '(I0,A,I0,A,I0,A)') size(kinds), ' feature edges, ', size(corners), ' corners, ', size(ps) - 1, &
                                   ' polylines; threshold 30 degrees'
call fig_open(u, 'ref-feature-edges-method', trim(caption), -50._R8P, 35._R8P)
call fig_surface(u, 'body', part)
call fig_segments(u, 'box', wa, wb)
call fig_segments(u, 'patch', a, b)
call fig_points(u, 'foot', c)
call fig_label(u, 'patch', 0.5_R8P * (p(1) + p(1) + up()), '90 deg: sharp')
call fig_label(u, 'box', 0.5_R8P * (p(4) + p(4) + up()), '22.5 deg: not sharp')
call fig_label(u, 'box', 0.5_R8P * (p(2) + p(2) + up()), '11.25 deg')
call fig_label(u, 'foot', p(1) + up(), 'corner: 3 feature edges')
call fig_label(u, 'patch', p(5) + up(), 'turns 22.5 deg: no corner')
call fig_close(u)
contains
   function up() result(v)
   !< The extrusion vector.
   type(vector_R8P) :: v

   v = vector_R8P(0._R8P, 0._R8P, 1._R8P)
   endfunction up

   subroutine tri(p1, p2, p3)
   !< The facet p1 p2 p3 (counter-clockwise seen from outside), appended to f.
   type(vector_R8P), intent(in) :: p1, p2, p3

   m = m + 1
   f(m)%vertex(1) = p1 ; f(m)%vertex(2) = p2 ; f(m)%vertex(3) = p3
   call f(m)%compute_metrix
   endsubroutine tri
endprogram fig_ref_feature_edges_method
