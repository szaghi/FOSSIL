!run tut08-ghosts tut08_a_ghosts
!fig tut08-ghosts
!< Tutorial 8, first run: the ghost points of a 41^3 grid over the domain: the nodes inside the dragon within one grid
!< spacing of the wall. The figure shows the grid plane y = 0, seen along y: fluid nodes (grey), ghost points (ember),
!< the solid nodes deeper inside (teal), and the section of the dragon by the plane.
program tut08_a_ghosts
use fossil_fig, only : fig_open, fig_segment, fig_points, fig_close
!region ghosts
use fossil,     only : surface_stl_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 41                         ! grid points per direction
type(surface_stl_object)      :: surface
type(curve_object)            :: section
type(vector_R8P), allocatable :: points(:), ghosts(:)
real(R8P),        allocatable :: d(:)
logical,          allocatable :: plane(:)
real(R8P)                     :: h, uv(2)
integer(I4P)                  :: i, j, k, p, m, u, s
character(96)                 :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / (n - 1)
allocate(points(n**3), d(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call surface%distance_many(points=points, distances=d, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
ghosts = pack(points, d < 0._R8P .and. d > -h)                  ! inside, within one spacing of the wall
m = size(ghosts)
print '(A,I0,A,I0,A,I0)', 'grid points ', n**3, ', inside ', count(d < 0._R8P), ', ghost points ', m
!endregion ghosts
plane = abs(points%y) < 0.5_R8P * h .and. abs(points%x) < 0.32_R8P .and. abs(points%z) < 0.42_R8P   ! y = 0, near the body
write(caption, '(A,I0,A,I0,A)') 'plane y = 0: ', count(plane .and. d < 0._R8P .and. d > -h), ' ghost points (ember), ', &
                                count(plane .and. d <= -h), ' deeper solid nodes (teal)'
print '(A)', trim(caption)
call fig_open(u, 'tut08-ghosts', trim(caption), -90._R8P, 0._R8P)
call slice_surface(surface, 2_I4P, 0._R8P, section)                ! (u, v) = (z, x)
do s = 1, section%get_vertices_number()
   if (section%get_next(s) == 0_I4P) cycle
   uv = section%get_vertex(s)
   call fig_segment(u, 'body', uv(2) * ex_R8P + uv(1) * ez_R8P, xz(section%get_vertex(section%get_next(s))))
enddo
call fig_points(u, 'box', pack(points, plane .and. d >= 0._R8P))
call fig_points(u, 'outside', pack(points, plane .and. d < 0._R8P .and. d > -h))
call fig_points(u, 'inside', pack(points, plane .and. d <= -h))
call fig_close(u)
contains
   function xz(w) result(q)
   !< The point of the plane y = 0 with section coordinates (u, v) = (z, x).
   real(R8P), intent(in) :: w(2) !< (z, x).
   type(vector_R8P)      :: q    !< Point.

   q = w(2) * ex_R8P + w(1) * ez_R8P
   endfunction xz
endprogram tut08_a_ghosts
