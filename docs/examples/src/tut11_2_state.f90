!run tut11-state tut11_2_state
!fig tut11-move
!< Tutorial 11, step 2: the cells whose centre changes state over the step of step 1 (state_changes). The figure: the
!< plane of cell centres y = const with the most state changes, the dragon cut there at pose0 (grey) and at pose1
!< (yellow), the cells flagged fresh (solid to fluid, teal) and dead (fluid to solid, ember).
program tut11_2_state
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, pose_identity, pose_from_axis_angle, &
                       slice_surface, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object)      :: surface
type(rigid_pose_t)            :: p0, p1
type(vector_R8P), allocatable :: centres(:)
logical,          allocatable :: fresh(:), dead(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, c, checked

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / n
allocate(centres(n**3), fresh(n**3), dead(n**3))
c = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   c = c + 1
   centres(c) = (-0.6_R8P + (i - 0.5_R8P) * h) * ex_R8P + (-0.6_R8P + (j - 0.5_R8P) * h) * ey_R8P + &
                (-0.6_R8P + (k - 0.5_R8P) * h) * ez_R8P
enddo ; enddo ; enddo
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
!region state
call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
print '(A,I0,A,I0)', 'fresh (solid to fluid) ', count(fresh), ', dead (fluid to solid) ', count(dead)
print '(A,I0,A,I0,A)', 'second query needed for ', checked, ' of ', size(centres), ' centres'
!endregion state
call figure
contains
   subroutine figure
   !< The plane y = const of cell centres with the most state changes, seen along y.
   real(R8P),    parameter       :: LIFT = 1.e-3_R8P
   type(surface_stl_object)      :: moved
   type(curve_object)            :: cut
   type(vector_R8P), allocatable :: a(:), b(:)
   character(80)                 :: caption
   real(R8P)                     :: yc
   logical                       :: in_plane(n**3)
   integer(I4P)                  :: u, jp, status

   jp = 1
   do j = 2, n
      if (count((fresh .or. dead) .and. plane(j)) > count((fresh .or. dead) .and. plane(jp))) jp = j
   enddo
   yc = -0.6_R8P + (jp - 0.5_R8P) * h
   print '(A,F7.4,A,I0,A,I0)', 'plane y = ', yc, ': fresh ', count(fresh .and. plane(jp)), ', dead ', &
                               count(dead .and. plane(jp))
   write(caption, '(A,F6.3,A)') 'state_changes, y = ', yc, ': pose0 grey, pose1 yellow, fresh teal, dead ember'
   call fig_open(u, 'tut11-move', trim(caption), -90._R8P, 0._R8P)
   in_plane = plane(jp)
   do c = 1, size(centres)
      if (.not. (fresh(c) .or. dead(c)) .or. .not. in_plane(c)) cycle
      call fig_quad(u, 'field', corner(c, -1, -1), corner(c, 1, -1), corner(c, 1, 1), corner(c, -1, 1), &
                    value=merge(-1._R8P, 1._R8P, fresh(c)))
   enddo
   call slice_surface(surface, 2_I4P, yc, cut, status=status)
   call curve_segments(cut, yc - LIFT, a, b)
   call fig_segments(u, 'box', a, b)
   moved = surface
   call moved%move(p1)
   call slice_surface(moved, 2_I4P, yc, cut, status=status)
   call curve_segments(cut, yc - LIFT, a, b)
   call fig_segments(u, 'accent', a, b)
   call fig_close(u)
   endsubroutine figure

   function plane(jj) result(mask)
   !< The cell centres of the plane j = jj.
   integer(I4P), intent(in) :: jj         !< Plane.
   logical                  :: mask(n**3) !< Mask.

   mask = abs(centres%y - (-0.6_R8P + (jj - 0.5_R8P) * h)) < 0.25_R8P * h
   endfunction plane

   function corner(cc, si, sk) result(v)
   !< Corner (si, sk = -1, 1 along x, z) of the cell of centre cc, in its plane y = const.
   integer(I4P), intent(in) :: cc, si, sk !< Cell, corner.
   type(vector_R8P)         :: v          !< Corner.

   v = centres(cc) ; v%x = v%x + 0.5_R8P * si * h ; v%z = v%z + 0.5_R8P * sk * h
   endfunction corner

   subroutine curve_segments(curve, level, a, b)
   !< The segments of a section by a plane y = const, drawn at y = level.
   type(curve_object),            intent(in)  :: curve !< Curve.
   real(R8P),                     intent(in)  :: level !< Plane drawn.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Ends.
   real(R8P)                                  :: uv(2), p(3)
   integer(I4P)                               :: s, m, ax

   ax = curve%get_normal_axis()
   allocate(a(curve%get_vertices_number()), b(curve%get_vertices_number())) ; m = 0
   do s = 1, curve%get_vertices_number()
      if (curve%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = curve%get_vertex(s)
      p(ax) = level ; p(mod(ax, 3) + 1) = uv(1) ; p(mod(ax + 1, 3) + 1) = uv(2) ; a(m) = p(1) * ex_R8P + p(2) * ey_R8P + &
                                                                                       p(3) * ez_R8P
      uv = curve%get_vertex(curve%get_next(s))
      p(ax) = level ; p(mod(ax, 3) + 1) = uv(1) ; p(mod(ax + 1, 3) + 1) = uv(2) ; b(m) = p(1) * ex_R8P + p(2) * ey_R8P + &
                                                                                       p(3) * ez_R8P
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine curve_segments
endprogram tut11_2_state
