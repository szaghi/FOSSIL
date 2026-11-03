!< FOSSIL test: shells and bodies, 3-D surfaces and 2-D curves (issue #37 step 1).
program fossil_test_bodies
!< FOSSIL test: shells and bodies, 3-D surfaces and 2-D curves (issue #37 step 1).
!<
!< A shell is a set of edge-connected facets; a body is an outer shell (even nesting depth) with
!< the cavities it contains; an island inside a cavity is a body of its own. Same rule for curve
!< loops (holes join their outer loop).
!< 3-D, built from cube-02 (unit cube): one cube (1 shell, 1 body); two disjoint cubes (2, 2);
!< a hollow cube, outer [-1,2]^3 with the reversed unit cube as cavity (2 shells, 1 body); the same
!< with an island [0.4,0.6]^3 inside the cavity (3 shells, 2 bodies).
!< 2-D: an annulus (1 body); an annulus with an island in its hole (2 bodies, the hole in body 1,
!< the island in body 2); two disjoint squares (2); a square and an open polyline (2).

use fossil, only : surface_stl_object, facet_object, curve_object
use penf,   only : I4P, R8P

implicit none

type(surface_stl_object)        :: cube, s
type(curve_object)              :: c
type(facet_object), allocatable :: unit(:), f(:)
type(facet_object), pointer     :: fp
integer(I4P)                    :: i, n
logical                         :: ok(6)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
n = cube%get_facets_number()
allocate(unit(n))
do i = 1, n
   fp => cube%facet_at(i)
   unit(i) = fp ; unit(i)%vertex_id = 0_I4P
enddo

ok(1) = cube%get_shells_number() == 1 .and. cube%get_bodies_number() == 1 .and. cube%facet_body(1_I4P) == 1
print '(A,L1)', 'one cube: ', ok(1)

f = [unit, moved(unit, 1._R8P, 3._R8P, .false.)]
call s%adopt_facets(facets=f)
ok(2) = s%get_shells_number() == 2 .and. s%get_bodies_number() == 2 .and. s%facet_body(1_I4P) == 1 .and. &
        s%facet_body(n + 1_I4P) == 2
print '(A,L1)', 'two disjoint cubes: ', ok(2)

f = [moved(unit, 3._R8P, -1._R8P, .false.), moved(unit, 1._R8P, 0._R8P, .true.)]
call s%adopt_facets(facets=f)
ok(3) = s%get_shells_number() == 2 .and. s%get_bodies_number() == 1 .and. s%facet_body(n + 1_I4P) == 1
print '(A,L1,2I3)', 'hollow cube (cavity joins its body): ', ok(3), s%get_shells_number(), s%get_bodies_number()

f = [moved(unit, 3._R8P, -1._R8P, .false.), moved(unit, 1._R8P, 0._R8P, .true.), moved(unit, 0.2_R8P, 0.4_R8P, .false.)]
call s%adopt_facets(facets=f)
ok(4) = s%get_shells_number() == 3 .and. s%get_bodies_number() == 2 .and. s%facet_body(n + 1_I4P) == 1 .and. &
        s%facet_body(2_I4P * n + 1_I4P) == 2
print '(A,L1,2I3)', 'island in a cavity (own body): ', ok(4), s%get_shells_number(), s%get_bodies_number()

! 2-D
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
ok(5) = c%get_bodies_number() == 1 .and. c%loop_body(2_I4P) == 1
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P,  2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P,  2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
ok(5) = ok(5) .and. c%get_bodies_number() == 2 .and. c%loop_body(1_I4P) == 1 .and. c%loop_body(2_I4P) == 1 .and. &
        c%loop_body(3_I4P) == 2
print '(A,L1)', '2-D annulus, island in the hole: ', ok(5)
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P,  3._R8P, 4._R8P, 4._R8P, 3._R8P], &
                   v=[0._R8P, 0._R8P, 1._R8P, 1._R8P,  0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
ok(6) = c%get_bodies_number() == 2 .and. c%loop_body(2_I4P) == 2
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P,  0.2_R8P, 0.8_R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P,  0.5_R8P, 0.5_R8P], &
                   loop_sizes=[4_I4P, 2_I4P], closed=[.true., .false.])
ok(6) = ok(6) .and. c%get_bodies_number() == 2 .and. c%loop_body(2_I4P) == 2
print '(A,L1)', '2-D disjoint squares, open polyline: ', ok(6)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function moved(src, scale, shift, reverse) result(dst)
   !< Facets scaled by `scale` and shifted by `shift` on every axis, winding reversed if asked.
   type(facet_object), intent(in) :: src(:)  !< Source facets.
   real(R8P),          intent(in) :: scale   !< Scale.
   real(R8P),          intent(in) :: shift   !< Shift.
   logical,            intent(in) :: reverse !< Reverse the winding.
   type(facet_object)             :: dst(size(src)) !< Result.
   integer(I4P)                   :: k, j    !< Counters.

   dst = src
   do k = 1, size(src)
      do j = 1, 3
         dst(k)%vertex(j)%x = scale * src(k)%vertex(j)%x + shift
         dst(k)%vertex(j)%y = scale * src(k)%vertex(j)%y + shift
         dst(k)%vertex(j)%z = scale * src(k)%vertex(j)%z + shift
      enddo
      if (reverse) then
         dst(k)%vertex(2) = src(k)%vertex(3) ; dst(k)%vertex(3) = src(k)%vertex(2)
         do j = 2, 3
            dst(k)%vertex(j)%x = scale * dst(k)%vertex(j)%x + shift
            dst(k)%vertex(j)%y = scale * dst(k)%vertex(j)%y + shift
            dst(k)%vertex(j)%z = scale * dst(k)%vertex(j)%z + shift
         enddo
      endif
      dst(k)%vertex_id = 0_I4P
      call dst(k)%compute_metrix
   enddo
   endfunction moved
endprogram fossil_test_bodies
