!run fig12 fig_tut12
!fig tut12-part
!fig tut12-thickness
!< Figures of tutorial 12: the bracket of the chapter (a base, a round boss united, a bore drilled 0.06 from the x = 0
!< end) after isotropic_remesh. First its patches (base sepia, boss teal, bore ember) with its feature edges (gold) and
!< corners. Then, seen from the x = 0 end, each facet coloured by its ray thickness from teal (thin) to ember (0.5 and
!< more), the facets of the thin regions under h = 0.1 outlined in gold: the two faces of the web left by the bore.
program fig_tut12
use fossil,                    only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                                      REM_STATUS_OK, THICKNESS_RAY, STATUS_OK
use fossil_fig,                only : fig_open, fig_triangle, fig_segments, fig_points, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter              :: PI = acos(-1._R8P)
type(surface_stl_object)          :: part, boss, bore, mesh
type(vertex_pool_object), pointer :: pool
type(facet_object),       pointer :: f
type(vector_R8P), allocatable     :: a(:), b(:), c(:)
integer(I4P),     allocatable     :: edges(:,:), kinds(:), corners(:), region(:)
real(R8P),        allocatable     :: thickness(:), area(:)
character(100)                    :: caption
character(8)                      :: style
integer(I4P)                      :: u, i, k, n, status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss')
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore')
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
mesh = part
call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=6_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'

! patches, feature edges and corners
call mesh%feature_edges(edges=edges, kinds=kinds, corners=corners, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
pool => mesh%get_vertex_pool()
allocate(a(size(kinds)), b(size(kinds)), c(size(corners)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
do i = 1, size(corners)
   c(i) = pool%coord(corners(i))
enddo
write(caption, '(A,I0,A,I0,A,I0,A)') 'patches base, boss (teal), bore (ember); ', size(kinds), ' feature edges, ', &
                                      size(corners), ' corners'
call fig_open(u, 'tut12-part', trim(caption), -60._R8P, 30._R8P)
do i = 1, mesh%get_facets_number()
   f => mesh%facet_at(i)
   select case (mesh%patch_name(f%patch))
   case ('boss') ; style = 'probe'
   case ('bore') ; style = 'patch'
   case default  ; style = 'body'
   endselect
   call fig_triangle(u, trim(style), f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_segments(u, 'foot', a, b)
call fig_points(u, 'accent', c)
call fig_close(u)

! ray thickness and the thin regions under h = 0.1
call mesh%thickness(thickness=thickness, method=THICKNESS_RAY, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call mesh%thin_regions(thickness=thickness, h=0.1_R8P, factor=1._R8P, region=region, region_area=area)
n = count(region > 0)
deallocate(a, b) ; allocate(a(3 * n), b(3 * n))
k = 0
do i = 1, mesh%get_facets_number()
   if (region(i) == 0) cycle
   f => mesh%facet_at(i)
   a(k + 1) = f%vertex(1) ; b(k + 1) = f%vertex(2)
   a(k + 2) = f%vertex(2) ; b(k + 2) = f%vertex(3)
   a(k + 3) = f%vertex(3) ; b(k + 3) = f%vertex(1)
   k = k + 3
enddo
write(caption, '(A,F6.4,A,I0,A,I0,A)') 'ray thickness, teal thin to ember 0.5: thinnest ', minval(thickness), ', ', &
                                       n, ' facets under h = 0.1 in ', size(area), ' regions'
call fig_open(u, 'tut12-thickness', trim(caption), -150._R8P, 30._R8P)
do i = 1, mesh%get_facets_number()
   f => mesh%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), &
                     value=2._R8P * min(1._R8P, thickness(i) / 0.5_R8P) - 1._R8P)
enddo
call fig_segments(u, 'foot', a, b)
call fig_close(u)
contains
   subroutine box(s, lo, hi, name)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)
   character(*),             intent(in)    :: name
   integer(I4P)                            :: k

   call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call s%resize(factor=(hi(1) - lo(1)) * ex_R8P + (hi(2) - lo(2)) * ey_R8P + (hi(3) - lo(3)) * ez_R8P)
   call s%translate(delta=lo(1) * ex_R8P + lo(2) * ey_R8P + lo(3) * ez_R8P)
   call s%analyze
   call s%set_patch(facets=[(k, k = 1, s%get_facets_number())], name=name)
   endsubroutine box

   subroutine cylinder(s, c, r, z0, z1, m, name)
   !< A regular m-gon of radius r about c, extruded over [z0, z1], caps fanned, outward, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: c(2), r, z0, z1
   integer(I4P),             intent(in)    :: m
   character(*),             intent(in)    :: name
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: x(m), y(m)
   integer(I4P)                            :: k, j, n

   x = [(c(1) + r * cos(2 * PI * (k - 1) / m), k = 1, m)] ; y = [(c(2) + r * sin(2 * PI * (k - 1) / m), k = 1, m)]
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
   call s%set_patch(facets=[(k, k = 1, s%get_facets_number())], name=name)
   endsubroutine cylinder

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_metrix
   endsubroutine tri
endprogram fig_tut12
