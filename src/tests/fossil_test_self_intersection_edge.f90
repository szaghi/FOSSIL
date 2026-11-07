!< FOSSIL test: facets sharing an edge, decided exactly (issue #69, step 3).
!<
!< 1. Coplanar fold: A = (0,0,0), (1,0,0), (0,1,0) and B on the same side of the shared edge (0,0,0)-(1,0,0): one
!<    SELFINT_EDGE_FOLD; B on the other side (a flat manifold pair): none.
!< 2. Four facets around the edge (0,0,0)-(0,0,1), free vertices at 0, 90, 180 and 270 degrees. Two flat sheets
!<    crossing along it (0/180 and 90/270, each sheet with opposite directions on its two sides): the order of
!<    directions is +, +, -, -, two SELFINT_SHEETS pairs. Two bent sheets touching along it (0/90 and 180/270): +, -,
!<    +, -, none.
!< 3. Two unit cubes sharing the vertical edge x = y = 1: the second turned by 225 degrees about it overlaps the
!<    first in a 45-degree wedge (sheets crossing: reported); turned by 45 degrees it only touches the first along the
!<    edge (none).
!< Negative control (base 04e5bbf): pairs sharing an edge were skipped, so 1 and 2 report none. In 3 the overlapping
!< cubes also cross away from the edge (16 crossing pairs on the base), so 3 checks the SELFINT_SHEETS kind.
program fossil_test_self_intersection_edge

use fossil,              only : surface_stl_object, intersection_pair_t, SELFINT_EDGE_FOLD, SELFINT_SHEETS
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P, ez_R8P

implicit none

real(R8P), parameter                   :: PI = acos(-1._R8P)
real(R8P), parameter                   :: A(9) = [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P]
type(intersection_pair_t), allocatable :: found(:)
type(facet_object),        allocatable :: f(:)
logical                                :: ok(3)
integer(I4P)                           :: n1, n2

! 1. fold
allocate(f(2))
call set(f(1), A) ; call set(f(2), [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0.5_R8P, 0.5_R8P, 0._R8P])
call run(f, found)
ok(1) = size(found) == 1
if (ok(1)) ok(1) = found(1)%kind == SELFINT_EDGE_FOLD
call set(f(2), [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0.5_R8P, -0.5_R8P, 0._R8P])
call run(f, found)
ok(1) = ok(1) .and. size(found) == 0
print '(A,L1)', '1. coplanar fold: ', ok(1)

! 2. fans of four around the z edge: +1 = (p0, p1, w), -1 = (p1, p0, w)
deallocate(f) ; allocate(f(4))
call fan(f, [0._R8P, 90._R8P, 180._R8P, 270._R8P], [1, 1, -1, -1])   ! flat sheets 0/180 and 90/270
call run(f, found)
n1 = count(found%kind == SELFINT_SHEETS)
call fan(f, [0._R8P, 90._R8P, 180._R8P, 270._R8P], [1, -1, 1, -1])   ! bent sheets 0/90 and 180/270
call run(f, found)
n2 = size(found)
ok(2) = n1 == 2 .and. n2 == 0
print '(A,L1,A,I0,A,I0)', '2. fans: ', ok(2), '  crossing sheets ', n1, '  touching sheets ', n2

! 3. two cubes sharing the edge x = y = 1
n1 = cubes(225._R8P) ; n2 = cubes(45._R8P)
ok(3) = n1 > 0 .and. n2 == 0
print '(A,L1,A,I0,A,I0)', '3. cubes on an edge: ', ok(3), '  overlapping ', n1, '  touching ', n2

print '(A,L1)', 'Are all tests passed? ', all(ok)
contains
   subroutine set(t, v)
   !< Facet from 9 coordinates.
   type(facet_object), intent(inout) :: t
   real(R8P),          intent(in)    :: v(9)
   integer(I4P)                      :: k

   do k = 1, 3
      t%vertex(k) = vector_R8P(v(3 * k - 2), v(3 * k - 1), v(3 * k))
   enddo
   endsubroutine set

   subroutine fan(t, deg, dirs)
   !< Facets around the edge (0,0,0)-(0,0,1), free vertex at the angle deg(k) in the plane z = 0.5, direction dirs(k).
   type(facet_object), intent(inout) :: t(:)
   real(R8P),          intent(in)    :: deg(:)
   integer,            intent(in)    :: dirs(:)
   real(R8P)                         :: w(3)
   integer(I4P)                      :: k

   do k = 1, size(t)
      w = [cos(deg(k) * PI / 180._R8P), sin(deg(k) * PI / 180._R8P), 0.5_R8P]
      if (dirs(k) > 0) then
         call set(t(k), [0._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, w])
      else
         call set(t(k), [0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, w])
      endif
   enddo
   endsubroutine fan

   subroutine run(t, found)
   !< Self-intersections of the surface of facets t.
   type(facet_object),                     intent(in)  :: t(:)
   type(intersection_pair_t), allocatable, intent(out) :: found(:)
   type(surface_stl_object)                            :: s
   type(facet_object), allocatable                     :: c(:)
   integer(I4P)                                        :: k

   c = t
   do k = 1, size(c)
      c(k)%id = k
   enddo
   call s%adopt_facets(facets=c)
   call s%find_self_intersections(pairs=found)
   endsubroutine run

   function cubes(angle_deg) result(nsheets)
   !< Sheet pairs of cube-01 plus a copy turned by angle_deg about the z axis through its corner (0, 0), moved to (1, 1).
   real(R8P), intent(in)           :: angle_deg
   integer(I4P)                    :: nsheets
   type(surface_stl_object)        :: c1, c2
   type(facet_object), allocatable :: all_f(:)
   type(facet_object), pointer     :: p
   type(intersection_pair_t), allocatable :: pr(:)
   integer(I4P)                    :: k, n

   call c1%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
   call c2%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
   call c2%rotate(axis=ez_R8P, angle=angle_deg * PI / 180._R8P)
   call c2%translate(delta=vector_R8P(1._R8P, 1._R8P, 0._R8P))
   n = c1%get_facets_number()
   allocate(all_f(2 * n))
   do k = 1, n
      p => c1%facet_at(k) ; all_f(k) = p
      p => c2%facet_at(k) ; all_f(n + k) = p
   enddo
   call run(all_f, pr)
   nsheets = count(pr%kind == SELFINT_SHEETS)
   endfunction cubes
endprogram fossil_test_self_intersection_edge
