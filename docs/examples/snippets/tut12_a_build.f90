!< Tutorial 12, step 1: a machined bracket built from primitives by exact booleans, each step checked against the
!< volume of the polygonal solids it should give.
program tut12_a_build
use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: part, boss, bore
real(R8P)                :: a_boss, a_bore, v_exact
integer(I4P)             :: status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')      ! a 2 x 1 x 0.5 base
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss') ! a round boss on it
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore') ! a bore, 0.06 from x = 0
call report('base', part)
call report('boss', boss)
call report('bore', bore)
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
a_boss = 16 * 0.35_R8P**2 * sin(2 * PI / 32)                       ! area of the 32-gon of radius 0.35
v_exact = 1._R8P + a_boss * 0.4_R8P                                 ! the boss stands 0.4 above the base
call report('base + boss', part)
print '(A,ES9.2)', '  volume - exact polygonal volume: ', part%get_volume() - v_exact
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
a_bore = 16 * 0.31_R8P**2 * sin(2 * PI / 32)
v_exact = v_exact - a_bore * 0.5_R8P                                ! the bore goes through the 0.5 of the base
call report('base + boss - bore', part)
print '(A,ES9.2)', '  volume - exact polygonal volume: ', part%get_volume() - v_exact
contains
   subroutine report(what, s)
   !< Facets, volume, watertightness and non-empty patches of a surface.
   character(*),             intent(in) :: what
   type(surface_stl_object), intent(in) :: s
   character(60)                        :: names
   character(19)                        :: label
   integer(I4P)                         :: p, k

   label = what ; names = ''
   do p = 1, s%get_patches_number()
      if (any(s%facet_patch([(k, k = 1, s%get_facets_number())]) == p)) names = trim(names)//' '//s%patch_name(p)
   enddo
   print '(2A,I4,A,F8.5,A,L1,2A)', label, ': ', s%get_facets_number(), ' facets, volume ', s%get_volume(), &
                                       ', watertight ', s%is_watertight(), ', patches', trim(names)
   endsubroutine report

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
endprogram tut12_a_build
