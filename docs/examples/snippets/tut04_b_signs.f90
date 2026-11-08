!< Tutorial 4, second run: the four sign algorithms at three points: the origin (inside), a point beyond the head
!< (outside), and an outside point whose closest feature is a vertex with a wrong pseudo-normal on this open dragon.
!< The figure shows that point, its closest point and the vertex pseudo-normal, which points away from it.
program tut04_b_signs
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_close
use fossil,     only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_RAY_INTERSECTIONS, &
                       SIGN_SOLID_ANGLE
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: ALG(4) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_RAY_INTERSECTIONS, &
                                      SIGN_SOLID_ANGLE]
type(surface_stl_object) :: surface
type(vector_R8P)         :: p(3), foot, normal
real(R8P)                :: d
integer(I4P)             :: k, a, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
p(1) = 0._R8P * ex_R8P
p(2) = 0.6_R8P * ey_R8P
p(3) = 0.36_R8P * ex_R8P - 0.2_R8P * ey_R8P + 0.2_R8P * ez_R8P
print '(A)', '        point          pseudo-n.  winding n.  rays  solid angle'
do k = 1, 3
   print '(3F7.3,4L10)', p(k)%x, p(k)%y, p(k)%z, &
                         (surface%is_point_inside(point=p(k), sign_algorithm=ALG(a)), a = 1, 4)
enddo
call surface%compute_distance(point=p(3), distance=d, is_signed=.true., is_square_root=.true., closest_point=foot, &
                              normal=normal)
print '(A,F8.4,A,3F8.4)', 'point 3: pseudo-normal d ', d, ', vertex pseudo-normal ', normal%x, normal%y, normal%z
call fig_open(u, 'tut04-signs', 'point 3 (ember), its closest vertex (gold), the vertex pseudo-normal (teal)', &
              -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', surface)
call fig_segment(u, 'accent', p(3), foot)
call fig_segment(u, 'inside', foot, foot + 0.15_R8P * normal)
call fig_points(u, 'outside', [p(3)])
call fig_points(u, 'foot', [foot])
call fig_close(u)
endprogram tut04_b_signs
