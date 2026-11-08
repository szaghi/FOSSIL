!< Example of the reference page "2-D geometry (curves)": a curve built from vertex arrays with two defects (a hole
!< wound the wrong way, a loop that misses closing by 1e-7), checked, repaired, then queried.
program ref_2d_geometry
use fossil, only : curve_object, CURVE_VALID, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: c
real(R8P)            :: u(16), v(16), d, foot(2), p(2, 3)
integer(I4P)         :: k, l, defect, seg, vtx, gaps, reversed, still_open

u(1:4) = [0._R8P, 2._R8P, 2._R8P, 0._R8P] ; v(1:4) = [0._R8P, 0._R8P, 1._R8P, 1._R8P]  ! outer plate, counter-clockwise
do k = 0, 7                                                                          ! hole, an octagon: counter-clockwise,
   u(5 + k) = 0.6_R8P + 0.25_R8P * cos(PI * k / 4) ; v(5 + k) = 0.5_R8P + 0.25_R8P * sin(PI * k / 4) ! the wrong way
enddo
u(13:16) = [3._R8P, 4._R8P, 3._R8P, 3._R8P] ; v(13:16) = [0._R8P, 0._R8P, 1._R8P, 1.e-7_R8P] ! a triangle that misses
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 8_I4P, 4_I4P], closed=[.true., .true., .false.]) ! closing: open chain
defect = c%check_validity(loop=l)
print '(A,I0,A,I0)', 'before repair: defect ', defect, ' at loop ', l
call c%repair(tolerance=1.e-6_R8P, gaps_closed=gaps, loops_reversed=reversed, loops_open=still_open)
print '(A,I0,A,I0,A,I0)', 'repair: gaps closed ', gaps, ', loops reversed ', reversed, ', still open ', still_open
defect = c%check_validity(loop=l, check_self_intersections=.true.)
print '(A,L1,A,I0,A,I0,A,F8.5)', 'after repair: valid ', defect == CURVE_VALID, ', loops ', c%get_loops_number(), &
                                 ', bodies ', c%get_bodies_number(), ', area ', c%get_area()
do l = 1, c%get_loops_number()
   print '(A,I0,A,I0,A,F9.5)', '  loop ', l, ': nesting depth ', c%nesting_depth(l), ', signed area ', c%get_loop_area(l)
enddo
p(:, 1) = [1.5_R8P, 0.5_R8P] ; p(:, 2) = [0.6_R8P, 0.5_R8P] ; p(:, 3) = [2.3_R8P, 0.2_R8P]   ! solid, hole, outside
do k = 1, 3
   call c%compute_distance(p(:, k), distance=d, closest_point=foot, closest_segment=seg, closest_vertex=vtx)
   print '(A,2F5.2,A,F8.5,A,2F8.5,A,I0,A,I0)', 'p = ', p(:, k), ': d ', d, ', foot ', foot, ', segment ', seg, &
                                              ', vertex ', vtx
   print '(A,F8.5,A,F8.5,A,I0,A,L1,A,L1)', '  w fast ', c%winding_number(p(:, k)), ', w sum ', &
         c%winding_number(p(:, k), use_brute_force=.true.), ', w exact ', c%winding_number_exact(p(:, k)), &
         ', inside ', c%is_point_inside(p(:, k)), ', inside (fwn) ', &
         c%is_point_inside(p(:, k), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo
endprogram ref_2d_geometry
