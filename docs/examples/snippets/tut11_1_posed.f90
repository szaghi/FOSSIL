!< Tutorial 11, step 1: the dragon of chapters 2-3 at a pose. The same query asked in two ways: on the body as built with
!< `pose=` (world points mapped into the body frame) and on a copy moved there (`move`, tree refitted). The figure: the
!< dragon at rest (grey) and at the pose (solid), a circle of probe points with the segments to their closest points.
program tut11_1_posed
use fossil,     only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_surface, fig_points, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                         ! cell centres per direction
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object)      :: surface, moved
type(rigid_pose_t)            :: p0, p1
type(vector_R8P), allocatable :: centres(:)
real(R8P),        allocatable :: d_posed(:), d_moved(:), foot_posed(:,:), foot_moved(:,:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, c

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)   ! chapters 2-3
h = 1.2_R8P / n
allocate(centres(n**3), d_posed(n**3), d_moved(n**3), foot_posed(3, n**3), foot_moved(3, n**3))
c = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   c = c + 1
   centres(c) = (-0.6_R8P + (i - 0.5_R8P) * h) * ex_R8P + (-0.6_R8P + (j - 0.5_R8P) * h) * ey_R8P + &
                (-0.6_R8P + (k - 0.5_R8P) * h) * ez_R8P
enddo ; enddo ; enddo
p0 = pose_identity()                                               ! the rest pose
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, &  ! 10 degrees about y
                          translation=[0.02_R8P, 0._R8P, 0._R8P])                ! and 0.02 along x
print '(A,3F8.4)', 'pose1: R row 1', p1%r(1, :)
print '(A,3F8.4)', '       t      ', p1%t
call surface%query_many(points=centres, distances=d_posed, closest_points=foot_posed, &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER, pose=p1)        ! the tree stays as built
moved = surface
call moved%move(p1)                                                        ! a moved copy, tree refitted
call moved%query_many(points=centres, distances=d_moved, closest_points=foot_moved, &
                      sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,ES9.2)', 'distances:      posed - moved, max ', maxval(abs(d_posed - d_moved))
print '(A,ES9.2)', 'closest points: posed - moved, max ', maxval(abs(foot_posed - foot_moved))
print '(A,I0,A,I0)', 'centres inside at pose1: ', count(d_posed < 0._R8P), ' of ', size(centres)
call figure
contains
   subroutine figure
   !< The dragon at rest and at pose1, seen along y; a circle of probes in the plane y = 0 and their feet (posed queries).
   type(vector_R8P) :: probes(16), feet(16)
   real(R8P)        :: d(16), cp(3, 16)
   integer(I4P)     :: u, m

   do m = 1, 16                                                   ! a circle of probes in the plane y = 0
      probes(m) = 0.62_R8P * cos(2 * pi * (m - 0.5_R8P) / 16) * ex_R8P + 0.62_R8P * sin(2 * pi * (m - 0.5_R8P) / 16) * ez_R8P
   enddo
   call surface%query_many(points=probes, distances=d, closest_points=cp, pose=p1)
   do m = 1, 16
      feet(m)%x = cp(1, m) ; feet(m)%y = cp(2, m) ; feet(m)%z = cp(3, m)
   enddo
   call fig_open(u, 'tut11-posed', 'pose1 = 10 deg about y + 0.02 x: rest grey, posed solid, probes and their feet', &
                 -90._R8P, 0._R8P)
   call fig_surface(u, 'ghost', surface)
   call fig_surface(u, 'body', moved)
   call fig_segments(u, 'foot', probes, feet)
   call fig_points(u, 'probe', probes)
   call fig_points(u, 'foot', feet)
   call fig_close(u)
   endsubroutine figure
endprogram tut11_1_posed
