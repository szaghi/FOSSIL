!run tut15-face tut15_2_face
!fig tut15-face
!< Tutorial 15, step 2: a 41 x 21 block face above the upper side of the wing projected onto it by project_block_face,
!< twice: 0.005 off the wing with tolerance 0.01 (larger than the last chordwise steps: vertices collapse onto the
!< trailing-edge corners), then 0.001 off with tolerance 0.0012. Animation: the root end of the trailing edge seen
!< from behind, the face before and after each projection, its vertices coloured by the constraint taken.
program tut15_2_face
!region face
use fossil,     only : surface_stl_object, projection_target_t, curvilinear_block_t, project_block_face, &
                       PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_label, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P

implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
real(R8P), parameter      :: OFFSET(2) = [0.005_R8P, 0.0005_R8P], TOLERANCE(2) = [0.01_R8P, 0.0014_R8P]
type(surface_stl_object)  :: wing
type(projection_target_t) :: target
type(curvilinear_block_t) :: block(2), block0(2)
real(R8P),    allocatable :: moved(:,:)
integer(I4P), allocatable :: kinds(:,:,:), k1(:,:)
integer(I4P)              :: t, i, status, u
character(8)              :: num

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%projection_target(target)
allocate(kinds(41, 21, 2))
block(1) = face_block(0._R8P)
print '(A,2F8.5)', 'last two chordwise steps at the trailing edge: ', block(1)%x(40, 1, 1) - block(1)%x(39, 1, 1), &
                                                                   block(1)%x(41, 1, 1) - block(1)%x(40, 1, 1)
do t = 1, 2
   block(t) = face_block(OFFSET(t))                ! layers k = 1, 2: face 5 is k = 1, the one near the wing
   block0(t) = block(t)
   call project_block_face(block(t), face=5_I4P, target=target, tolerance=TOLERANCE(t), moved=moved, &
                           kinds=k1, status=status)
   kinds(:, :, t) = k1
   print '(A,F7.4,A,F7.4,A,I0)', 'face ', OFFSET(t), ' off the wing, tolerance ', TOLERANCE(t), ': status ', status
   print '(A,I0,A,I0,A,I0,A,F7.4)', '  on corners ', count(kinds(:, :, t) == PROJ_CORNER), ', on curves ', &
         count(kinds(:, :, t) == PROJ_CURVE), ', on the patch ', count(kinds(:, :, t) == PROJ_PATCH), &
         ', largest move ', maxval(moved)
   print '(A,L1,A,L1,A,I0)', '  trailing-edge row on x = 0.75: ', all(block(t)%x(41, :, 1) == 0.75_R8P), &
         ', root and tip rows on z = 0, 1: ', all(block(t)%z(:, 1, 1) == 0._R8P) .and. all(block(t)%z(:, 21, 1) == 1._R8P), &
         ', face edges of zero length: ', zero_edges(block(t))
enddo
!endregion face

call wing%rotate(axis=ex_R8P, angle=PI / 2._R8P)              ! drawn with the span down the page: (x, -z, y)
call fig_open(u, 'tut15-face', '', -90._R8P, 90._R8P)          ! the root end of the trailing edge, from above
call fig_zoom(u, vector_R8P(0.7465_R8P, -0.0065_R8P, 0._R8P), 0.0085_R8P)
call fig_surface(u, 'ghost', wing)
call fig_anim(u, 2._R8P)
do t = 1, 2
   call fig_frame(u, caption=caption(t, 'before'))
   call grid_lines(u, 'box', block0(t))
   do i = 38, 41
      write(num, '(A,I0)') 'i=', i
      call fig_label(u, 'body', vector_R8P(block0(t)%x(i, 1, 1) - 0.0008_R8P, -0.0035_R8P - 0.0012_R8P * (i - 38), &
                     0.01_R8P), trim(num))
   enddo
   call fig_frame(u, caption=caption(t, 'projected'))
   call grid_lines(u, 'box', block(t))
   call fig_points(u, 'body', pack(points(block(t)), kinds(:, :, t) == PROJ_PATCH))
   call fig_points(u, 'accent', pack(points(block(t)), kinds(:, :, t) == PROJ_CURVE))
   call fig_points(u, 'probe', pack(points(block(t)), kinds(:, :, t) == PROJ_CORNER))
enddo
call fig_close(u)
contains
   function face_block(off) result(bl)
   !< A 41 x 21 x 2 block over the upper side, its face k = 1 `off` above the wing, cosine-clustered chordwise.
   real(R8P), intent(in)     :: off !< Distance of the face k = 1 from the wing.
   type(curvilinear_block_t) :: bl  !< Block.
   real(R8P)                 :: x(41, 21, 2), y(41, 21, 2), z(41, 21, 2), xc, yt
   integer(I4P)              :: i, j

   do j = 1, 21
      do i = 1, 41
         xc = 0.5_R8P * (1._R8P - cos(PI * (i - 1) / 40._R8P))
         yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - &
                         0.1015_R8P * xc**4)
         x(i, j, :) = xc - 0.25_R8P ; y(i, j, :) = [yt + off, yt + 0.1_R8P] ; z(i, j, :) = (j - 1) / 20._R8P
      enddo
   enddo
   call bl%initialize(x, y, z)
   endfunction face_block

   function zero_edges(bl) result(n)
   !< Edges of the face k = 1 whose two ends coincide.
   type(curvilinear_block_t), intent(in) :: bl !< Block.
   integer(I4P)                          :: n  !< Edges.
   integer(I4P)                          :: i, j

   n = 0
   do j = 1, 21
      do i = 1, 41
         if (i < 41) n = n + merge(1, 0, same(bl, i, j, i + 1, j))
         if (j < 21) n = n + merge(1, 0, same(bl, i, j, i, j + 1))
      enddo
   enddo
   endfunction zero_edges

   pure logical function same(bl, i, j, k, l)
   !< Vertices (i, j, 1) and (k, l, 1) coincide.
   type(curvilinear_block_t), intent(in) :: bl         !< Block.
   integer(I4P),              intent(in) :: i, j, k, l !< Vertices.

   same = bl%x(i, j, 1) == bl%x(k, l, 1) .and. bl%y(i, j, 1) == bl%y(k, l, 1) .and. bl%z(i, j, 1) == bl%z(k, l, 1)
   endfunction same

   function caption(t, what) result(c)
   !< Caption of a frame.
   integer(I4P), intent(in) :: t    !< Case.
   character(*), intent(in) :: what !< Before or projected.
   character(:), allocatable :: c   !< Caption.
   character(120)           :: buf

   write(buf, '(A,F7.4,A,F7.4,A,A)') 'face', OFFSET(t), ' off, tolerance', TOLERANCE(t), ': ', what
   if (what == 'projected') buf = trim(buf)//merge(' (3 vertices on one corner)', '                           ', t == 1)
   c = trim(buf)
   endfunction caption

   subroutine grid_lines(un, style, bl)
   !< The lines of the face k = 1 of a block.
   integer(I4P),              intent(in) :: un    !< Unit of the scene.
   character(*),              intent(in) :: style !< Style.
   type(curvilinear_block_t), intent(in) :: bl    !< Block.
   type(vector_R8P)                      :: a(1660), b(1660), p(41, 21)
   integer(I4P)                          :: i, j, n

   p = points(bl) ; n = 0
   do j = 1, 21
      do i = 1, 40
         n = n + 1 ; a(n) = p(i, j) ; b(n) = p(i + 1, j)
      enddo
   enddo
   do i = 1, 41
      do j = 1, 20
         n = n + 1 ; a(n) = p(i, j) ; b(n) = p(i, j + 1)
      enddo
   enddo
   call fig_segments(un, style, a(1:n), b(1:n))
   endsubroutine grid_lines

   function points(bl) result(p)
   !< The vertices of the face k = 1 of a block, drawn as (x, -z, y).
   type(curvilinear_block_t), intent(in) :: bl       !< Block.
   type(vector_R8P)                      :: p(41, 21) !< Points.
   integer(I4P)                          :: i, j

   do j = 1, 21
      do i = 1, 41
         p(i, j) = vector_R8P(bl%x(i, j, 1), -bl%z(i, j, 1), bl%y(i, j, 1))
      enddo
   enddo
   endfunction points
endprogram tut15_2_face
