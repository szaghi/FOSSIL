!run figref-moving-bodies fig_ref_moving_bodies
!fig ref-moving-bodies-posed
!fig ref-moving-bodies-swept
!< Figures of the reference page "Moving bodies":
!< - body-frame queries: the unit cube at rest (edges) and at a pose (translucent) (turned 35 degrees about z and lifted), a ring of
!<   world points queried with `query_many(pose=)` on the tree built at rest; each point joined to its closest point,
!<   which comes back in the world frame, on the posed cube (teal: outside, ember: inside);
!< - swept geometry over one step of a screw motion (0.35 rad about a tilted axis, plus a shift): the cells of the
!<   middle layer of a 12 x 12 x 3 block seen from above, shaded by their swept volume over the cell volume (ember:
!<   fluid gained, teal: fluid lost), with the sections of the cube at the two poses at mid-layer (dark: start,
!<   yellow: end).
program fig_ref_moving_bodies
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, pose_identity, pose_from_axis_angle, &
                       swept_cells_3d_t, swept_cells_3d, slice_surface
use fossil_fig, only : fig_open, fig_box, fig_surface, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NP = 40
real(R8P),    parameter       :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(surface_stl_object)      :: cube, moved
type(curve_object)            :: section
type(rigid_pose_t)            :: p0, p1
type(swept_cells_3d_t)        :: sw
type(vector_R8P)              :: points(NP)
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: d(NP), foot(3, NP), x(13), z(4), phi
integer(I4P)                  :: u, i, j, n, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)                   ! centred at the origin
call cube%analyze

! ---- body-frame queries: the tree stays as built at rest, the points and the results are in the world frame
p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=35._R8P * PI / 180._R8P, translation=[0.2_R8P, 0._R8P, 0.25_R8P])
do i = 1, NP
   phi = 2._R8P * PI * (i - 1) / NP
   points(i) = vec(0.2_R8P + 0.6_R8P * cos(phi), 0.6_R8P * sin(phi), 0.25_R8P)
enddo
call cube%query_many(points=points, distances=d, closest_points=foot, pose=p1)
print '(A,I0,A,I0,A,F7.4,A,F7.4)', 'posed queries: outside ', count(d > 0._R8P), ', inside ', count(d < 0._R8P), &
                                   ', distances ', minval(d), ' to ', maxval(d)
moved = cube
call moved%move(p1)
call fig_open(u, 'ref-moving-bodies-posed', 'query_many(pose=): world points, feet on the posed cube (teal out, ember in)', &
              -60._R8P, 62._R8P)
call fig_box(u, 'box', vec(-0.5_R8P, -0.5_R8P, -0.5_R8P), vec(0.5_R8P, 0.5_R8P, 0.5_R8P))   ! the cube at rest
call fig_surface(u, 'ghost', moved)
allocate(a(NP), b(NP))
do i = 1, NP
   a(i) = points(i) ; b(i) = vec(foot(1, i), foot(2, i), foot(3, i))
enddo
call fig_segments(u, 'foot', a, b)
call fig_points(u, 'foot', b)
call fig_points(u, 'probe', pack(points, d > 0._R8P))
call fig_points(u, 'outside', pack(points, d < 0._R8P))
call fig_close(u)

! ---- swept geometry of one step of a screw motion
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0.3_R8P, 0.2_R8P, 1._R8P], angle=0.35_R8P, translation=[0.15_R8P, 0.05_R8P, 0._R8P])
x = [(-1._R8P + i / 6._R8P, i = 0, 12)]                ! 12 x 12 x 3 cells; the faces x, y = +-0.5 start on grid planes
z = [-0.75_R8P, -0.25_R8P, 0.25_R8P, 0.75_R8P]
call swept_cells_3d(cube, p0, p1, x, x, z, sw, status=status)
if (.not. sw%exact) error stop 'swept_cells_3d: off the exact path (issue #84)'
print '(A,L1,A,ES9.2,A,ES9.2)', 'swept: exact ', sw%exact, ', sum of swept ', sum(sw%swept), &
                                 ', largest residual ', maxval(abs(sw%residual))
call fig_open(u, 'ref-moving-bodies-swept', &
              'swept volume / cell volume (ember gained, teal lost); start dark, end yellow', &
              -90._R8P, 90._R8P)
do j = 1, 12 ; do i = 1, 12
   call fig_quad(u, 'field', vec(x(i), x(j), 0._R8P), vec(x(i + 1), x(j), 0._R8P), vec(x(i + 1), x(j + 1), 0._R8P), &
                 vec(x(i), x(j + 1), 0._R8P), value=sw%swept(i, j, 2) / ((x(i + 1) - x(i)) * (x(j + 1) - x(j)) * 0.5_R8P))
enddo ; enddo
moved = cube
call outline(moved, 'cell')
call moved%move(p1)
call outline(moved, 'accent')
call fig_close(u)
contains
   pure function vec(x, y, z) result(v)
   !< A vector.
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: v

   v%x = x ; v%y = y ; v%z = z
   endfunction vec

   subroutine outline(s, style)
   !< The section of s by the plane z = 0 (the middle of the layer), drawn just above the cells.
   type(surface_stl_object), intent(in) :: s
   character(*),             intent(in) :: style
   type(vector_R8P), allocatable        :: a(:), b(:)
   real(R8P)                            :: uv(2)
   integer(I4P)                         :: i, n, nv

   call slice_surface(s, 3_I4P, 0._R8P, section)
   nv = section%get_vertices_number()
   allocate(a(nv), b(nv)) ; n = 0
   do i = 1, nv
      if (section%get_next(i) == 0_I4P) cycle
      n = n + 1
      uv = section%get_vertex(i) ; a(n) = vec(uv(1), uv(2), LIFT)
      uv = section%get_vertex(section%get_next(i)) ; b(n) = vec(uv(1), uv(2), LIFT)
   enddo
   call fig_segments(u, style, a(1:n), b(1:n))
   endsubroutine outline
endprogram fig_ref_moving_bodies
