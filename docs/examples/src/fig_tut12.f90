!run fig12 fig_tut12
!fig tut12-build
!fig tut12-remesh
!fig tut12-deviation
!fig tut12-part
!fig tut12-thickness
!fig tut12-mass
!< Figures of tutorial 12: the bracket of the chapter (a base, a round boss united, a bore drilled 0.06 from the x = 0
!< end; tut12_a_build) and its remesh towards edges of 0.2 (tut12_b_remesh). The boolean result by patch with the cut
!< curves; the boolean result next to the remesh; the remesh coloured by its distance to the part; patches with the
!< feature edges and corners; the ray thickness with the thin regions under h = 0.1; the principal axes of inertia.
program fig_tut12
use fossil,                    only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                                      REM_STATUS_OK, THICKNESS_RAY, STATUS_OK
use fossil_fig,                only : fig_open, fig_zoom, fig_surface, fig_triangle, fig_segment, fig_segments, &
                                      fig_points, fig_label, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter              :: PI = acos(-1._R8P)
type(surface_stl_object)          :: part, mesh, shifted
type(vertex_pool_object), pointer :: pool
type(facet_object),       pointer :: f, g
type(vector_R8P), allocatable     :: a(:), b(:), c(:), samples(:)
type(vector_R8P)                  :: centre_v, worst
integer(I4P),     allocatable     :: edges(:,:), kinds(:), corners(:), region(:)
real(R8P),        allocatable     :: thickness(:), area(:), d(:), dfacet(:)
real(R8P)                         :: mass, centre(3), moments(3), axes(3,3), dmax, r
character(120)                    :: caption
character(4)                      :: num
integer(I4P)                      :: u, i, k, e, n, status

call bracket(part, mesh)

! the boolean result by patch, with the cut curves: the edges between facets of two patches
n = 0
do i = 1, part%get_facets_number()
   f => part%facet_at(i)
   do e = 1, 3
      if (f%fcon_edge(e) > i) then
         if (part%facet_patch(f%fcon_edge(e)) /= f%patch) n = n + 1
      endif
   enddo
enddo
allocate(a(n), b(n)) ; n = 0
do i = 1, part%get_facets_number()
   f => part%facet_at(i)
   do e = 1, 3
      if (f%fcon_edge(e) > i) then
         if (part%facet_patch(f%fcon_edge(e)) /= f%patch) then
            n = n + 1 ; a(n) = f%vertex(e) ; b(n) = f%vertex(mod(e, 3) + 1)
         endif
      endif
   enddo
enddo
write(caption, '(A,I0,A,I0,A)') 'union, difference: ', part%get_facets_number(), &
                               ' facets; boss teal, bore ember; ', n, ' cut edges gold'
call fig_open(u, 'tut12-build', trim(caption), -60._R8P, 30._R8P)
call by_patch(part)
call fig_segments(u, 'foot', a, b)
call fig_close(u)

! the boolean result (left) and its remesh (right)
shifted = mesh
call shifted%translate(x=2.4_R8P, y=0._R8P, z=0._R8P)
write(caption, '(A,I0,A,I0,A)') 'isotropic_remesh, target 0.2: ', part%get_facets_number(), ' facets (left) to ', &
                                mesh%get_facets_number(), ' (right)'
call fig_open(u, 'tut12-remesh', trim(caption), -75._R8P, 35._R8P, height=330_I4P)
call by_patch(part)
call by_patch(shifted)
call fig_close(u)

! the distance from the remesh to the part, per facet: the largest of 7 samples (centroid, edge midpoints, inner points)
n = mesh%get_facets_number()
allocate(samples(7 * n), d(7 * n), dfacet(n))
do i = 1, n
   f => mesh%facet_at(i)
   samples(7 * i - 6) = (f%vertex(1) + f%vertex(2) + f%vertex(3)) / 3._R8P
   do k = 1, 3
      samples(7 * i - 6 + k) = 0.5_R8P * (f%vertex(k) + f%vertex(mod(k, 3) + 1))
      samples(7 * i - 3 + k) = (4._R8P * f%vertex(k) + f%vertex(mod(k, 3) + 1) + f%vertex(mod(k + 1, 3) + 1)) / 6._R8P
   enddo
enddo
call part%distance_many(points=samples, distances=d, is_signed=.false., is_square_root=.true.)
dfacet = [(maxval(d(7 * i - 6:7 * i)), i = 1, n)]
dmax = maxval(dfacet) ; k = maxloc(dfacet, dim=1)
f => mesh%facet_at(k) ; worst = (f%vertex(1) + f%vertex(2) + f%vertex(3)) / 3._R8P
write(caption, '(A,ES8.2)') 'remesh to part, per facet: cream 0 to ember ', dmax
call fig_open(u, 'tut12-deviation', trim(caption), -60._R8P, 30._R8P)
call fig_zoom(u, vector_R8P(1.5_R8P, 0.5_R8P, 0.6_R8P), 0.45_R8P)
do i = 1, n
   f => mesh%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), value=dfacet(i) / dmax)
enddo
call fig_points(u, 'accent', [worst])
call fig_label(u, 'accent', worst, 'largest')
call fig_close(u)

! patches, feature edges and corners
call mesh%feature_edges(edges=edges, kinds=kinds, corners=corners, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
pool => mesh%get_vertex_pool()
deallocate(a, b) ; allocate(a(size(kinds)), b(size(kinds)), c(size(corners)))
do i = 1, size(kinds)
   a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
enddo
do i = 1, size(corners)
   c(i) = pool%coord(corners(i))
enddo
write(caption, '(A,I0,A,I0,A)') 'patches base, boss (teal), bore (ember); ', size(kinds), ' feature edges, ', &
                                 size(corners), ' corners'
call fig_open(u, 'tut12-part', trim(caption), -60._R8P, 30._R8P)
call by_patch(mesh)
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
write(caption, '(A,F6.4,A,I0,A)') 'ray thickness, teal thin to ember 0.5: min ', minval(thickness), ', ', &
                                  n, ' facets < 0.1'
call fig_open(u, 'tut12-thickness', trim(caption), -150._R8P, 30._R8P)
do i = 1, mesh%get_facets_number()
   f => mesh%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), &
                     value=2._R8P * min(1._R8P, thickness(i) / 0.5_R8P) - 1._R8P)
enddo
call fig_segments(u, 'foot', a, b)
call fig_close(u)

! centre of mass and principal axes, each drawn over its radius of gyration sqrt(I / m) on both sides
call mesh%mass_properties(mass=mass, center=centre, principal_moments=moments, principal_axes=axes, &
                          density=2700._R8P, status=status)
if (status /= STATUS_OK) error stop 'not a closed body'
centre_v = vector_R8P(centre(1), centre(2), centre(3))
write(caption, '(A,F7.1,A,3F6.3,A)') 'mass ', mass, ' kg, centre', centre, ', principal axes'
call fig_open(u, 'tut12-mass', trim(caption), -60._R8P, 30._R8P)
call fig_surface(u, 'ghost', mesh)
call fig_points(u, 'accent', [centre_v])
do k = 1, 3
   r = sqrt(moments(k) / mass)
   a = [centre_v + r * vector_R8P(axes(1, k), axes(2, k), axes(3, k))]
   b = [centre_v - r * vector_R8P(axes(1, k), axes(2, k), axes(3, k))]
   call fig_segment(u, 'accent', b(1), a(1))
   write(num, '(I4)') nint(moments(k))
   call fig_label(u, 'accent', a(1), 'I'//achar(48 + k)//' = '//trim(adjustl(num)))
enddo
call fig_points(u, 'box', [vector_R8P(1._R8P, 0.5_R8P, 0.25_R8P)])
call fig_label(u, 'box', vector_R8P(1._R8P, 0.5_R8P, 0.25_R8P), 'base centre')
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

   subroutine bracket(part, mesh)
   !< The bracket of tut12_a_build (base, boss united, bore drilled) and, if asked, its remesh of tut12_b_remesh.
   type(surface_stl_object), intent(inout)           :: part
   type(surface_stl_object), intent(inout), optional :: mesh
   type(surface_stl_object)                          :: boss, bore
   integer(I4P)                                      :: st

   call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
   call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss')
   call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore')
   call part%boolean(boss, op=BOOL_UNION, status=st)
   if (st /= BOOL_STATUS_OK) error stop 'union failed'
   call part%boolean(bore, op=BOOL_DIFFERENCE, status=st)
   if (st /= BOOL_STATUS_OK) error stop 'difference failed'
   if (.not. present(mesh)) return
   mesh = part
   call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=6_I4P, preserve_features=.true., status=st)
   if (st /= REM_STATUS_OK) error stop 'remesh failed'
   endsubroutine bracket

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
