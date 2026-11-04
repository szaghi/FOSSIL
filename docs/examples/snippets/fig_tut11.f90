!< Figures of tutorial 11: the step of the chapter (a pitch of 10 degrees about y and 0.02 along x). First the plane of
!< cell centres y = const with the most state changes: the dragon cut there at pose0 (grey) and moved to pose1 (move,
!< yellow), the cells state_changes flags fresh (solid to fluid, teal) and dead (fluid to solid, ember). Then the
!< section translated by 0.02: each cell shaded by its swept area (swept_cells_2d, sw%swept / h^2: ember gains fluid,
!< teal loses it), the section at pose0 (grey) and moved to pose1 (yellow).
!< The figure programs run before the tutorial ones, so the placed dragon (chapters 2-3) and its section (chapter 7)
!< are rebuilt here with the same calls.
program fig_tut11
use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                       pose_from_axis_angle, pose_2d, swept_cells_2d_t, swept_cells_2d, slice_surface, &
                       SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
real(R8P),    parameter       :: LIFT = 1.e-3_R8P            ! the sections drawn just above the cells (seen from +z)
type(surface_stl_object)      :: surface, moved
type(curve_object)            :: section, shifted, cut
type(rigid_pose_t)            :: p0, p1
type(rigid_pose_2d_t)         :: q0, q1
type(swept_cells_2d_t)        :: sw
type(vector_R8P)              :: bmin, bmax, centre
type(vector_R8P), allocatable :: centres(:), a(:), b(:)
logical,          allocatable :: fresh(:), dead(:)
character(80)                 :: caption
real(R8P)                     :: h, x(49), scale, r, yc
integer(I4P)                  :: u, i, j, k, c, jp, checked, open_loops, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! chapter 2: repair
call surface%sanitize(fill_holes=.true.)
call surface%save_into_file(file_name='fig11-dragon-clean.stl')
call surface%load_from_file(file_name='fig11-dragon-clean.stl', guess_format=.true.) ! chapter 3: place
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig11-dragon-placed.stl')
call surface%load_from_file(file_name='fig11-dragon-placed.stl', guess_format=.true.)

h = 1.2_R8P / n                                                                     ! chapter 11: the step
allocate(centres(n**3), fresh(n**3), dead(n**3))
c = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   c = c + 1
   centres(c) = (-0.6_R8P + (i - 0.5_R8P) * h) * ex_R8P + (-0.6_R8P + (j - 0.5_R8P) * h) * ey_R8P + &
                (-0.6_R8P + (k - 0.5_R8P) * h) * ez_R8P
enddo ; enddo ; enddo
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
moved = surface
call moved%move(p1)
call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
jp = 1                                                                              ! the plane y = const of
do j = 2, n                                                                         ! cell centres with the most
   if (count((fresh .or. dead) .and. plane(j)) > count((fresh .or. dead) .and. plane(jp))) jp = j   ! state changes
enddo
yc = -0.6_R8P + (jp - 0.5_R8P) * h
write(caption, '(A,F6.3,A)') 'state_changes, y = ', yc, ': pose0 grey, pose1 yellow, fresh teal, dead ember'
call fig_open(u, 'tut11-move', trim(caption), -90._R8P, 0._R8P)
do c = 1, size(centres)
   if (.not. (fresh(c) .or. dead(c)) .or. abs(centres(c)%y - yc) > 0.25_R8P * h) cycle
   call fig_quad(u, 'field', cell_corner(c, -1, -1), cell_corner(c, 1, -1), cell_corner(c, 1, 1), &
                 cell_corner(c, -1, 1), value=merge(-1._R8P, 1._R8P, fresh(c)))      ! fresh teal, dead ember
enddo
call slice_surface(surface, 2_I4P, yc, cut, status=status)                          ! the body at pose0 and pose1
call curve_segments(cut, yc - LIFT)
call fig_segments(u, 'box', a, b)
call slice_surface(moved, 2_I4P, yc, cut, status=status)
call curve_segments(cut, yc - LIFT)
call fig_segments(u, 'accent', a, b)
call fig_close(u)

call slice_surface(surface, 3_I4P, 0._R8P, section, merge_collinear=.true., open_loops=open_loops, status=status)
call section%save_into_file(file_name='fig11-dragon-section.dat')                   ! chapter 7: the section
call section%load_from_file(file_name='fig11-dragon-section.dat')
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.02_R8P, 0._R8P])
call swept_cells_2d(section, q0, q1, x, x, sw)
call fig_open(u, 'tut11-swept', 'swept_cells_2d: swept area (ember gained, teal lost), pose0 grey, pose1 yellow', &
              0._R8P, 90._R8P)
call curve_segments(section, LIFT)
bmin%x = minval(a%x) ; bmin%y = minval(a%y) ; bmax%x = maxval(a%x) + 0.02_R8P ; bmax%y = maxval(a%y)
centre = 0.5_R8P * (bmin + bmax) ; centre%z = 0._R8P
r = max(0.5_R8P * (bmax%x - bmin%x), 0.3_R8P * (bmax%y - bmin%y)) + 0.03_R8P      ! x across, y along the card
call fig_zoom(u, centre, r)
do j = 1, 48 ; do i = 1, 48
   if (sw%swept(i, j) == 0._R8P) cycle
   call fig_quad(u, 'field', vec(x(i), x(j)), vec(x(i + 1), x(j)), vec(x(i + 1), x(j + 1)), vec(x(i), x(j + 1)), &
                 value=sw%swept(i, j) / 0.025_R8P**2)
enddo ; enddo
call fig_segments(u, 'box', a, b)
shifted = section
call shifted%move(q1)                                                                ! the section at pose1
call curve_segments(shifted, LIFT)
call fig_segments(u, 'accent', a, b)
call fig_close(u)
contains
   subroutine curve_segments(curve, level)
   !< The segments of a curve in its plane, at `level` along its normal axis, into a, b.
   type(curve_object), intent(in) :: curve !< Curve.
   real(R8P),          intent(in) :: level !< Position of the plane drawn.
   real(R8P)                      :: uv(2)
   integer(I4P)                   :: s, m

   if (allocated(a)) deallocate(a, b)
   allocate(a(curve%get_vertices_number()), b(curve%get_vertices_number())) ; m = 0
   do s = 1, curve%get_vertices_number()
      if (curve%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = curve%get_vertex(s) ; a(m) = in_plane(curve%get_normal_axis(), uv, level)
      uv = curve%get_vertex(curve%get_next(s)) ; b(m) = in_plane(curve%get_normal_axis(), uv, level)
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine curve_segments

   function in_plane(axis, uv, level) result(v)
   !< The point of the plane x(axis) = level at the in-plane coordinates uv (u, v: the other axes in cyclic order).
   integer(I4P), intent(in) :: axis  !< Normal axis.
   real(R8P),    intent(in) :: uv(2) !< In-plane coordinates.
   real(R8P),    intent(in) :: level !< Plane position.
   type(vector_R8P)         :: v     !< Point.
   real(R8P)                :: p(3)

   p(axis) = level ; p(mod(axis, 3) + 1) = uv(1) ; p(mod(axis + 1, 3) + 1) = uv(2)
   v%x = p(1) ; v%y = p(2) ; v%z = p(3)
   endfunction in_plane

   pure function plane(jj) result(mask)
   !< The cell centres of the plane j = jj.
   integer(I4P), intent(in) :: jj               !< Plane.
   logical                  :: mask(n**3)       !< Mask.

   mask = abs(centres%y - (-0.6_R8P + (jj - 0.5_R8P) * h)) < 0.25_R8P * h
   endfunction plane

   function cell_corner(cc, si, sk) result(v)
   !< Corner (si, sk = -1, 1 along x, z) of the cell of centre cc, in its plane y = const.
   integer(I4P), intent(in) :: cc, si, sk !< Cell, corner.
   type(vector_R8P)         :: v          !< Corner.

   v = centres(cc) ; v%x = v%x + 0.5_R8P * si * h ; v%z = v%z + 0.5_R8P * sk * h
   endfunction cell_corner

   function vec(p, q, z) result(v)
   !< The point (p, q, z), z = 0 by default.
   real(R8P), intent(in)           :: p, q !< In-plane coordinates.
   real(R8P), intent(in), optional :: z    !< Height.
   type(vector_R8P)                :: v    !< Point.

   v%x = p ; v%y = q ; v%z = 0._R8P
   if (present(z)) v%z = z
   endfunction vec
endprogram fig_tut11
