!run anim12 anim_tut12_bracket
!fig anim-tut12-bracket
!< Animation of tutorial 12: the bracket built step by step. The base, the boss, their union; the bore, the difference;
!< isotropic_remesh after 2 and 6 iterations; the feature edges and corners; the ray thickness with the facets under
!< h = 0.1; the principal axes of inertia. Every number in a caption is computed here.
program anim_tut12_bracket
use fossil,                    only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                                      REM_STATUS_OK, THICKNESS_RAY, STATUS_OK
use fossil_fig,                only : fig_open, fig_surface, fig_triangle, fig_segment, fig_segments, fig_points, &
                                      fig_anim, fig_frame, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter              :: PI = acos(-1._R8P)
type(surface_stl_object)          :: part, boss, bore, mesh
type(vertex_pool_object), pointer :: pool
type(facet_object),       pointer :: g
type(vector_R8P), allocatable     :: a(:), b(:), c(:)
type(vector_R8P)                  :: centre_v, axis
integer(I4P),     allocatable     :: edges(:,:), kinds(:), corners(:), region(:)
real(R8P),        allocatable     :: thickness(:), area(:)
real(R8P)                         :: mass, centre(3), moments(3), axes(3,3)
character(100)                    :: caption
integer(I4P)                      :: u, i, k, status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss')
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore')
call fig_open(u, 'anim-tut12-bracket', 'the bracket, step by step', -60._R8P, 30._R8P)
call fig_anim(u, 1.4_R8P)
write(caption, '(A,I0,A)') 'box: the base, 2 x 1 x 0.5, ', part%get_facets_number(), ' facets'
call fig_frame(u, trim(caption))
call by_patch(part)
write(caption, '(A,I0,A)') 'cylinder: a 32-gon prism, the boss, ', boss%get_facets_number(), ' facets'
call fig_frame(u, trim(caption))
call fig_surface(u, 'ghost', part)
call fig_surface(u, 'probe', boss)
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
write(caption, '(A,I0,A,F7.5)') 'boolean(BOOL_UNION): ', part%get_facets_number(), ' facets, volume ', &
                                part%get_volume()
call fig_frame(u, trim(caption))
call by_patch(part)
call fig_frame(u, 'cylinder: the bore, 0.06 from the x = 0 end')
call fig_surface(u, 'ghost', part)
call fig_surface(u, 'patch', bore)
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
write(caption, '(A,I0,A,F7.5)') 'boolean(BOOL_DIFFERENCE): ', part%get_facets_number(), ' facets, volume ', &
                                part%get_volume()
call fig_frame(u, trim(caption), 2._R8P)
call by_patch(part)
do k = 2, 6, 4
   mesh = part
   call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=k, preserve_features=.true., status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   write(caption, '(A,I0,A,I0,A)') 'isotropic_remesh, ', k, ' iterations: ', mesh%get_facets_number(), ' facets'
   call fig_frame(u, trim(caption))
   call by_patch(mesh)
enddo
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
write(caption, '(A,I0,A,I0,A)') 'feature_edges: ', size(kinds), ' sharp edges, ', size(corners), ' corners'
call fig_frame(u, trim(caption))
call fig_surface(u, 'body', mesh)
call fig_segments(u, 'foot', a, b)
call fig_points(u, 'accent', c)
call mesh%thickness(thickness=thickness, method=THICKNESS_RAY, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call mesh%thin_regions(thickness=thickness, h=0.1_R8P, factor=1._R8P, region=region, region_area=area)
write(caption, '(A,I0,A,F6.4)') 'thickness (teal thin): ', count(region > 0), ' facets under 0.1, min ', &
                                minval(thickness)
call fig_frame(u, trim(caption))
do i = 1, mesh%get_facets_number()
   g => mesh%facet_at(i)
   call fig_triangle(u, 'field', g%vertex(1), g%vertex(2), g%vertex(3), &
                     value=2._R8P * min(1._R8P, thickness(i) / 0.5_R8P) - 1._R8P)
enddo
call mesh%mass_properties(mass=mass, center=centre, principal_moments=moments, principal_axes=axes, &
                          density=2700._R8P, status=status)
if (status /= STATUS_OK) error stop 'not a closed body'
centre_v = vector_R8P(centre(1), centre(2), centre(3))
write(caption, '(A,F7.1,A)') 'mass_properties: ', mass, ' kg in aluminium, principal axes'
call fig_frame(u, trim(caption), 2.5_R8P)
call fig_surface(u, 'ghost', mesh)
call fig_points(u, 'accent', [centre_v])
do k = 1, 3
   axis = sqrt(moments(k) / mass) * vector_R8P(axes(1, k), axes(2, k), axes(3, k))
   call fig_segment(u, 'accent', centre_v - axis, centre_v + axis)
enddo
call fig_close(u)
contains
   subroutine by_patch(s)
   !< The facets of s coloured by patch: base sepia, boss teal, bore ember.
   type(surface_stl_object), intent(in) :: s
   integer(I4P)                         :: j

   do j = 1, s%get_facets_number()
      g => s%facet_at(j)
      select case (s%patch_name(g%patch))
      case ('boss') ; call fig_triangle(u, 'probe', g%vertex(1), g%vertex(2), g%vertex(3))
      case ('bore') ; call fig_triangle(u, 'patch', g%vertex(1), g%vertex(2), g%vertex(3))
      case default  ; call fig_triangle(u, 'body', g%vertex(1), g%vertex(2), g%vertex(3))
      endselect
   enddo
   endsubroutine by_patch

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
endprogram anim_tut12_bracket
