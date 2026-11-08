!run tut12-features tut12_d_features
!< Tutorial 12, step 4: the patches the bracket of steps 1-2 carries, and its feature edges before and after the
!< remesh.
program tut12_d_features
use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, REM_STATUS_OK, &
                   FEATURE_EDGE_SHARP, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: part, mesh
integer(I4P), allocatable :: ids(:), kinds(:), corners(:), ps(:), pv(:)
integer(I4P)              :: i, n, status

call bracket(part, mesh)                                                    ! steps 1 and 2, in a helper below
!region patches
ids = [(i, i = 1, mesh%get_facets_number())]
do i = 1, mesh%get_patches_number()
   n = count(mesh%facet_patch(ids) == i)
   if (n > 0) print '(3A,I0,A)', 'patch ', mesh%patch_name(i), ': ', n, ' facets'   ! empty ones skipped
enddo
!endregion patches
!region features
call part%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv)
print '(A,I0,A,I0,A,I0)', 'part:     feature edges ', size(kinds), ', corners ', size(corners), &
                          ', polylines ', size(ps) - 1
call mesh%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0,A,I0)', 'remeshed: feature edges ', size(kinds), ', corners ', size(corners), &
                               ', polylines ', size(ps) - 1, ', all sharp ', count(kinds == FEATURE_EDGE_SHARP)
n = count([(pv(ps(i)) == pv(ps(i + 1) - 1), i = 1, size(ps) - 1)])         ! first vertex = last: a loop
print '(A,I0,A,I0)', '  closed loops without a corner ', n, ', chains between corners ', size(ps) - 1 - n
!endregion features
contains
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
endprogram tut12_d_features
