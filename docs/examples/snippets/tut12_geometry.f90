!< Tutorial 12: a machined bracket built by booleans, remeshed for a solver, checked by the certified deviation, then
!< analysed: patches, feature edges, thin walls, mass properties.
program tut12_geometry
use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, REM_STATUS_OK, &
                   DEV_STATUS_OK, FEATURE_EDGE_SHARP, THICKNESS_RAY, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: part, boss, bore, mesh
integer(I4P), allocatable :: kinds(:), corners(:), ps(:), pv(:), region(:)
real(R8P),    allocatable :: thickness(:), area(:), d_vertex(:)
real(R8P)                 :: lower, upper, mass, centre(3), moments(3), v_part
integer(I4P)              :: i, n, status, invalid

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')      ! a 2 x 1 x 0.5 base
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss') ! a round boss on it
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore') ! a bore, 0.06 from x = 0
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
v_part = part%get_volume()
print '(A,I0,A,F8.5,A,L1)', 'part: ', part%get_facets_number(), ' facets, volume ', v_part, &
                            ', watertight ', part%is_watertight()
mesh = part
call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=6_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
print '(A,I0,A,F7.4,A,L1)', 'remeshed: ', mesh%get_facets_number(), ' facets, longest edge ', &
                            mesh%largest_edge_len(), ', watertight ', mesh%is_watertight()
call mesh%deviation(part, lower=lower, upper=upper, per_vertex=d_vertex, status=status)  ! how far the mesh strays
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2,A,ES9.2)', 'h(mesh, part): ', lower, ' <= h <= ', upper, ', at the vertices ', maxval(d_vertex)
call part%deviation(mesh, lower=lower, upper=upper, tolerance=5.e-3_R8P, status=status) ! what of the part it misses
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2)', 'h(part, mesh): ', lower, ' <= h <= ', upper
print '(A,ES9.2)', 'volume change, relative: ', (mesh%get_volume() - v_part) / v_part
do i = 1, mesh%get_patches_number()
   n = count(mesh%facet_patch(facets()) == i)
   if (n > 0) print '(3A,I0,A)', 'patch ', mesh%patch_name(i), ': ', n, ' facets'   ! empty ones skipped
enddo
call part%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv)
print '(A,I0,A,I0,A,I0)', 'part:     feature edges ', size(kinds), ', corners ', size(corners), &
                          ', polylines ', size(ps) - 1
call mesh%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0,A,I0)', 'remeshed: feature edges ', size(kinds), ', corners ', size(corners), &
                               ', polylines ', size(ps) - 1, ', all sharp ', count(kinds == FEATURE_EDGE_SHARP)
call mesh%thickness(thickness=thickness, method=THICKNESS_RAY, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call mesh%thin_regions(thickness=thickness, h=0.1_R8P, factor=1._R8P, region=region, region_area=area)
print '(A,I0,A,I0,A,I0)', 'thinner than h = 0.1: ', count(region > 0), ' facets in ', size(area), &
                          ' regions; invalid ', invalid
do i = 1, size(area)
   print '(A,I0,A,I0,A,F7.4,A,F7.4)', '  region ', i, ': ', count(region == i), ' facets, area ', area(i), &
                                      ', thinnest ', minval(thickness, mask=region == i)
enddo
call mesh%mass_properties(mass=mass, center=centre, principal_moments=moments, density=2700._R8P, status=status)
if (status /= STATUS_OK) error stop 'not a closed body'
print '(A,F9.2,A,3F8.4)', 'mass ', mass, ', centre of mass', centre
print '(A,3F10.3)', 'principal moments', moments
contains
   function facets() result(ids)
   !< All facet ids of the mesh.
   integer(I4P), allocatable :: ids(:)

   ids = [(i, i = 1, mesh%get_facets_number())]
   endfunction facets

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
endprogram tut12_geometry
