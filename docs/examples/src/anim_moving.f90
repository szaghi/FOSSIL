!run anim-moving anim_moving
!fig anim-moving
!< Animation of tutorial 11: the dragon passing through a plane of cells. The placed dragon of chapters 2-3 (rebuilt
!< here with the same calls, decimated to 700 facets for a light drawing) moves through the plane x = 0.01 of cell centres,
!< travelling 0.3 along x while rolling 45 degrees about x, its pose interpolated over 8 steps (pose_interpolate).
!< Each frame draws the body at that pose (move; the part on the viewer's side of the plane translucent, so the cells
!< show), its section by the plane (slice_surface, yellow) and the cells whose centre state_changes(pose0=previous,
!< pose1=current) flags: fresh (solid to fluid, teal) and dead (fluid to solid, ember).
program anim_moving
use fossil,              only : surface_stl_object, curve_object, rigid_pose_t, pose_from_axis_angle, pose_interpolate, &
                                slice_surface, SIGN_FAST_WINDING_NUMBER, DEC_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_triangle, fig_quad, fig_box, fig_segments, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: NY = 36, NZ = 30, NS = 8              ! cells of the plane, steps
real(R8P),    parameter       :: PI = 4._R8P * atan(1._R8P), H = 0.04_R8P, Y0 = -0.72_R8P, Z0 = -0.6_R8P
real(R8P),    parameter       :: XC = 0.01_R8P                          ! the plane of cell centres
type(surface_stl_object)      :: surface, moved
type(curve_object)            :: section
type(facet_object), pointer   :: f
type(rigid_pose_t)            :: start, finish, p0, p1
type(vector_R8P)              :: bmin, bmax, centres(NY * NZ)
type(vector_R8P), allocatable :: a(:), b(:)
logical                       :: fresh(NY * NZ), dead(NY * NZ)
character(80)                 :: caption
real(R8P)                     :: scale
integer(I4P)                  :: u, i, k, c, s, checked, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)            ! chapter 2: repair
call surface%sanitize(fill_holes=.true.)
call surface%save_into_file(file_name='anim-moving-dragon-clean.stl')
call surface%load_from_file(file_name='anim-moving-dragon-clean.stl', guess_format=.true.) ! chapter 3: place
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%decimate(target_facets=700_I4P, status=status)                          ! a light drawing
if (status /= DEC_STATUS_OK) print '(A,I0)', 'decimate status ', status
print '(A,I0,A)', 'dragon: ', surface%get_facets_number(), ' facets'

c = 0                                                                                 ! the cell centres of the
do k = 1, NZ ; do i = 1, NY                                                           ! plane x = XC
   c = c + 1 ; centres(c) = vector_R8P(XC, Y0 + (i - 0.5_R8P) * H, Z0 + (k - 0.5_R8P) * H)
enddo ; enddo
start = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[-0.15_R8P, 0._R8P, 0._R8P])
finish = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI / 4._R8P, translation=[0.15_R8P, 0._R8P, 0._R8P])

call fig_open(u, 'anim-moving', '', -25._R8P, 15._R8P)
call fig_anim(u, 0.7_R8P)
call fig_box(u, 'box', vector_R8P(XC, Y0, Z0), vector_R8P(XC, Y0 + NY * H, Z0 + NZ * H))   ! the plane of cells
p1 = start
do s = 0, NS
   p0 = p1 ; p1 = pose_interpolate(start, finish, real(s, R8P) / NS)
   moved = surface
   call moved%move(p1)
   if (s == 0) then
      fresh = .false. ; dead = .false.
      call fig_frame(u, caption='pose0: the dragon and the plane of cell centres x = 0.01', seconds=1.4_R8P)
   else
      call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                                 sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
      write(caption, '(A,I2,A,I0,A,I0,A,I0,A)') 'step ', s, ': state_changes: ', count(fresh), ' fresh, ', &
                                               count(dead), ' dead (', checked, ' checked)'
      print '(A)', trim(caption)
      call fig_frame(u, caption=trim(caption), seconds=merge(2.0_R8P, 0.7_R8P, s == NS))
   endif
   do i = 1, moved%get_facets_number()                       ! the part on the viewer's side of the plane (x > XC)
      f => moved%facet_at(i)                                 ! translucent, the cells behind it visible
      call fig_triangle(u, trim(merge('ghost', 'body ', f%centroid%x > XC)), f%vertex(1), f%vertex(2), f%vertex(3))
   enddo
   do c = 1, size(centres)
      if (.not. (fresh(c) .or. dead(c))) cycle
      call fig_quad(u, 'field', corner(c, -1, -1), corner(c, 1, -1), corner(c, 1, 1), corner(c, -1, 1), &
                    value=merge(-1._R8P, 1._R8P, fresh(c)))  ! fresh teal, dead ember
   enddo
   call slice_surface(moved, 1_I4P, XC, section, status=status)
   call curve_segments(section)
   call fig_segments(u, 'accent', a, b)
enddo
call fig_close(u)
contains
   function corner(cc, si, sk) result(v)
   !< Corner (si, sk = -1, 1 along y, z) of the cell of centre cc.
   integer(I4P), intent(in) :: cc, si, sk !< Cell, corner.
   type(vector_R8P)         :: v          !< Corner.

   v = centres(cc) ; v%y = v%y + 0.5_R8P * si * H ; v%z = v%z + 0.5_R8P * sk * H
   endfunction corner

   subroutine curve_segments(curve)
   !< The segments of a section of normal x, in the plane x = XC, into a, b.
   type(curve_object), intent(in) :: curve !< Section.
   real(R8P)                      :: uv(2)
   integer(I4P)                   :: j, m

   if (allocated(a)) deallocate(a, b)
   allocate(a(curve%get_vertices_number()), b(curve%get_vertices_number())) ; m = 0
   do j = 1, curve%get_vertices_number()
      if (curve%get_next(j) == 0_I4P) cycle
      m = m + 1
      uv = curve%get_vertex(j) ; a(m) = vector_R8P(XC, uv(1), uv(2))           ! normal x: u = y, v = z
      uv = curve%get_vertex(curve%get_next(j)) ; b(m) = vector_R8P(XC, uv(1), uv(2))
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine curve_segments
endprogram anim_moving
