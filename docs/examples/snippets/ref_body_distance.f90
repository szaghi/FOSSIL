!< Run example of the reference page "Body-body distance" and its figure: three unit cubes, the second 0.25 beyond the
!< first along x, the third 0.6 beyond it along y; proximity_pairs over the bodies of one surface holding all three and
!< proximity_pairs_of over the three surfaces, cutoff 0.62; the closest pair of every body pair by distance_to; and the
!< 2-D twin, a 64-gon of radius 0.5 against a unit square turned 45 degrees by a pose. The figure draws the three cubes
!< translucent and the closest pair of each body pair under the cutoff (gold), with its distance.
program ref_body_distance
use fossil,     only : surface_stl_object, curve_object, proximity_pairs_of, pose_2d, STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: cube(3), all
type(curve_object)        :: circle, square
integer(I4P), allocatable :: pairs(:, :)
real(R8P),    allocatable :: distances(:)
real(R8P)                 :: d, pa(3), pb(3), qa(2), qb(2), cu(64), cv(64)
character(16)             :: text
integer(I4P)              :: fa, fb, i, k, status
do k = 1, 3
   call cube(k)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
enddo
call cube(2)%translate(delta=vector_R8P(1.25_R8P, 0._R8P, 0._R8P)) ; call cube(2)%analyze   ! 0.25 from cube 1
call cube(3)%translate(delta=vector_R8P(0._R8P, 1.6_R8P, 0._R8P))  ; call cube(3)%analyze   ! 0.6 from cube 1
all = cube(1)
call all%merge_solids(other=cube(2)) ; call all%merge_solids(other=cube(3))   ! one surface, three bodies
call all%proximity_pairs(cutoff=0.62_R8P, pairs=pairs, distances=distances, status=status)
if (status /= STATUS_OK) error stop 'proximity_pairs failed'
print '(A,I0,A)', 'one surface, ', all%get_bodies_number(), ' bodies, cutoff 0.62:'
do i = 1, size(distances)
   print '(A,I0,A,I0,A,F7.4)', '  bodies ', pairs(1, i), ' and ', pairs(2, i), ': ', distances(i)
enddo
call proximity_pairs_of(surfaces=cube, cutoff=0.62_R8P, pairs=pairs, distances=distances)
print '(A)', 'three surfaces, cutoff 0.62:'
do i = 1, size(distances)
   print '(A,I0,A,I0,A,F7.4)', '  surfaces ', pairs(1, i), ' and ', pairs(2, i), ': ', distances(i)
enddo
call cube(1)%distance_to(cube(3), distance=d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb)
print '(A,F7.4,A,I0,A,I0)', 'cubes 1 and 3: ', d, ', facets ', fa, ' and ', fb
print '(A,3F6.2,A,3F6.2)', '  from', pa, ' to', pb
do i = 1, 64
   cu(i) = 0.5_R8P * cos(2 * PI * (i - 1) / 64) ; cv(i) = 0.5_R8P * sin(2 * PI * (i - 1) / 64)
enddo
call circle%adopt_loops(u=cu, v=cv, loop_sizes=[64_I4P])                       ! a 64-gon of radius 0.5
call square%adopt_loops(u=[-0.5_R8P, 0.5_R8P, 0.5_R8P, -0.5_R8P], v=[-0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P], &
                        loop_sizes=[4_I4P])                                    ! a unit square about the origin
call circle%distance_to(square, distance=d, point_self=qa, point_other=qb, &
                        pose_other=pose_2d(angle=PI / 4, translation=[1.5_R8P, 0._R8P]))   ! a corner towards the circle
print '(A,F8.5,A,F8.5,A,2F8.5,A,2F8.5)', '2-D: ', d, ' (exact ', 1._R8P - sqrt(0.5_R8P), ') from', qa, ' to', qb
call fig_open(k, name='ref-body-distance-pairs', caption='proximity_pairs, cutoff 0.62: the closest pair of each '// &
              'body pair (gold)', azimuth=-60._R8P, elevation=55._R8P)
do i = 1, 3
   call fig_surface(k, 'ghost', cube(i))
enddo
call draw_pair(1, 2) ; call draw_pair(1, 3)
call fig_close(k)
contains
   subroutine draw_pair(a, b)
   !< The closest pair of cubes a and b as a segment, labelled with its length.
   integer(I4P), intent(in) :: a, b
   call cube(a)%distance_to(cube(b), distance=d, point_self=pa, point_other=pb)
   call fig_segment(k, 'accent', vector_R8P(pa(1), pa(2), pa(3)), vector_R8P(pb(1), pb(2), pb(3)))
   write(text, '(F4.2)') d
   call fig_label(k, 'accent', vector_R8P(0.5_R8P * (pa(1) + pb(1)), 0.5_R8P * (pa(2) + pb(2)), 0.5_R8P * (pa(3) + pb(3))), &
                  trim(text))
   endsubroutine draw_pair
endprogram ref_body_distance
