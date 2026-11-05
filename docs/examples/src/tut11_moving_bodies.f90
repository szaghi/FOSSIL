!run tut11 tut11_moving_bodies
!< Tutorial 11: a moving dragon: posed queries, a moved copy, fresh and dead cells, swept geometry of the section.
program tut11_moving_bodies
!region pose
use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, &
                   pose_2d, swept_cells_2d_t, swept_cells_2d, inertia_to_world, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                        ! cell centres per direction
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object)      :: surface, moved, wrap
type(curve_object)            :: section
type(rigid_pose_t)            :: p0, p1
type(rigid_pose_2d_t)         :: q0, q1
type(swept_cells_2d_t)        :: sw
type(vector_R8P), allocatable :: centres(:)
real(R8P),        allocatable :: d_posed(:), d_moved(:)
logical,          allocatable :: fresh(:), dead(:)
real(R8P)                     :: h, x(49), mass, centre(3), moments(3), inertia(3,3), inertia1(3,3)
integer(I4P)                  :: i, j, k, c, checked, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / n
allocate(centres(n**3), d_posed(n**3), d_moved(n**3), fresh(n**3), dead(n**3))
c = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   c = c + 1
   centres(c) = (-0.6_R8P + (i - 0.5_R8P) * h) * ex_R8P + (-0.6_R8P + (j - 0.5_R8P) * h) * ey_R8P + &
                (-0.6_R8P + (k - 0.5_R8P) * h) * ez_R8P
enddo ; enddo ; enddo
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
!endregion pose
!region queries
call surface%distance_many(points=centres, distances=d_posed, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, pose=p1)    ! the tree stays as built
moved = surface
call moved%move(p1)                                                       ! a moved copy, tree refitted
call moved%distance_many(points=centres, distances=d_moved, is_signed=.true., is_square_root=.true., &
                         sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,L1)', 'posed = moved to 1e-14: ', maxval(abs(d_posed - d_moved)) < 1.e-14_R8P
call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells: fresh ', count(fresh), ', dead ', count(dead), ', second query ', checked, &
                               ' of ', size(centres)
!endregion queries
!region swept
call section%load_from_file(file_name='dragon-section.dat')
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.02_R8P, 0._R8P])            ! a pure translation, 0.8 of a cell
call swept_cells_2d(section, q0, q1, x, x, sw)
print '(A,L1,A,I0,A,I0)', 'section swept: exact ', sw%exact, ', cells with swept area ', count(sw%swept /= 0._R8P), &
                         ', cleared ', count(sw%cleared)
print '(A,L1,A,I0)', '  sum of swept = 0 to 1e-14: ', abs(sum(sw%swept)) < 1.e-14_R8P, &
                    ', cells with residual > 1e-14: ', count(abs(sw%residual) > 1.e-14_R8P)
print '(A,ES9.2,A,2I3,A)', '  largest residual ', maxval(abs(sw%residual)), ' in cell', maxloc(abs(sw%residual)), &
                           ': the knot of the section (chapter 7)'
!endregion swept
!region inertia
call surface%mass_properties(mass=mass, status=status)
print '(A,L1)', 'the placed dragon is a closed body: ', status == STATUS_OK               ! chapter 2: it is not
call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap, status=status)     ! a closed surrogate
call wrap%mass_properties(mass=mass, center=centre, inertia=inertia, principal_moments=moments, &
                          density=1000._R8P, status=status)                      ! a body of density 1000
if (status /= STATUS_OK) error stop 'the wrap is not closed'
print '(A,ES11.4,A,3F8.4)', 'mass ', mass, ', centre of mass', centre
print '(A,3ES11.4)', 'principal moments', moments
call wrap%move(p1)
call wrap%mass_properties(inertia=inertia1, density=1000._R8P)
print '(A,ES9.2)', 'moved wrap against inertia_to_world: ', maxval(abs(inertia1 - inertia_to_world(inertia, p1))) / &
                                                           maxval(abs(inertia))
!endregion inertia
endprogram tut11_moving_bodies
