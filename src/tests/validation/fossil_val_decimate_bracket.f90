!< FOSSIL validation: decimate on a body built by booleans (issue #81).
!<
!< The bracket of tutorial 12 with 128-sided prisms: cube-01 resized to a 2 x 1 x 0.5 base, a boss prism united, a
!< bore prism subtracted (exact booleans, 3 patches, 1868 facets). Decimated to 1000 and to 600 facets (each a fresh
!< call on the input). Asserted: a valid surface (watertight, manifold, one shell) with OK or NO_PROGRESS.
!< Not asserted, open (#81 follow-up): both calls stop at NO_PROGRESS at 1582 facets having lost 36.6 % of the volume
!< (59.5 % at 1580 facets before #81). Not a topology defect (the results are valid, before and after #81): the patch
!< boundaries (3 rims of 128 vertices) are locked (issue #66: two boundary vertices never merge) and the face-centre fans of
!< cube-01 cannot be emptied without folding, so the cheap collapses run out; QEM has no cost ceiling, so the loop
!< goes on with costly collapses (box corners sheared along the box edges) until the heap is empty.

program fossil_val_decimate_bracket

use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                   DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

real(R8P),    parameter :: PI       = acos(-1._R8P)
integer(I4P), parameter :: M        = 128_I4P                        !< Sides of the prisms.
integer(I4P), parameter :: TARGET(2) = [1000_I4P, 600_I4P]
type(surface_stl_object) :: part, boss, bore, dec
real(R8P)                :: v0, dv
integer(I4P)             :: status, t
logical                  :: are_tests_passed(2)

are_tests_passed = .false.
call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, M, 'boss')
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, M, 'bore')
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
v0 = part%get_volume()
print '(A,I0,A,F10.6,A,L1)', 'bracket: ', part%get_facets_number(), ' facets, volume ', v0, ', watertight ', &
      part%is_watertight()

do t = 1, size(TARGET)
   dec = part
   call dec%decimate(target_facets=TARGET(t), status=status)
   dv = (dec%get_volume() - v0) / v0
   print '(A,I0,A,I0,A,I0,A,ES10.2,A,L1,A,L1,A,I0)', 'decimate to ', TARGET(t), ': status ', status, ', facets ', &
         dec%get_facets_number(), ', relative volume change ', dv, ', watertight ', dec%is_watertight(), &
         ', manifold ', dec%is_manifold(), ', shells ', dec%get_shells_number()
   are_tests_passed(t) = (status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. &
                         dec%is_watertight() .and. dec%is_manifold() .and. dec%get_shells_number() == 1_I4P
   if (status == DEC_STATUS_OK) are_tests_passed(t) = are_tests_passed(t) .and. dec%get_facets_number() <= TARGET(t)
enddo

print '(A,2L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   subroutine box(s, lo, hi, name)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)
   character(*),             intent(in)    :: name
   integer(I4P)                            :: k

   call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
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
endprogram fossil_val_decimate_bracket
