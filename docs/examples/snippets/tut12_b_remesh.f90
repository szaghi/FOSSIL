!< Tutorial 12, step 2: the bracket of step 1 remeshed towards edges of 0.2, iteration after iteration; then the same
!< remesh towards 0.4, which coarsens it (before issue #79 the collapse pass rejected almost every edge).
program tut12_b_remesh
use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, REM_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: part, mesh
integer(I4P)             :: it, n0, status

call bracket(part)                                                          ! step 1, in a helper below
print '(A,I5,A,F7.4,A,F7.4)', 'boolean result: ', part%get_facets_number(), ' facets, edges from ', &
                               part%smallest_edge_len(), ' to', part%largest_edge_len()
do it = 1, 6
   if (it == 3 .or. it == 5) cycle
   mesh = part                                                              ! the part stays the reference
   call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=it, preserve_features=.true., status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   print '(A,I0,A,I5,A,F7.4,A,F7.4,A,L1)', 'iterations ', it, ': ', mesh%get_facets_number(), ' facets, edges from ', &
         mesh%smallest_edge_len(), ' to', mesh%largest_edge_len(), ', watertight ', mesh%is_watertight()
enddo
n0 = mesh%get_facets_number()
call mesh%isotropic_remesh(target_length=0.4_R8P, iterations=6_I4P, preserve_features=.true., status=status)
print '(A,I0,A,I0,A)', 'target 0.4 on the remeshed part: ', n0, ' -> ', mesh%get_facets_number(), &
                       ' facets'
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
endprogram tut12_b_remesh
