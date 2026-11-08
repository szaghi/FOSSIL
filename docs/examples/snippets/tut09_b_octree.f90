!< Tutorial 9, step 2: an octree refined level by level around the placed dragon. Each level is one box_classify_many
!< batch; a box is refined when it is cut, or near within the diagonal of its children (Matsushita & Aoki 2021). The
!< animation: the cut boxes of levels 3 to 5 as voxels, side 0.15 to 0.0375, then the dragon they wrap.
program tut09_b_octree
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_quad, fig_surface, fig_close
use vecfor,     only : vector_R8P
use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter   :: LEVELS = 5                     ! finest side 1.2 / 2^5 = 0.0375
type(surface_stl_object)  :: surface
real(R8P),    allocatable :: lo(:,:), hi(:,:), next_lo(:,:), next_hi(:,:)
integer(I4P), allocatable :: codes(:)
integer(I4P)              :: level, n, c, i, j, leaves, leaves_in, leaves_out
real(R8P)                 :: side
character(96)             :: caption
integer(I4P)              :: u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
lo = reshape([-0.6_R8P, -0.6_R8P, -0.6_R8P], [3, 1]) ; hi = -lo          ! the root box, the domain [-0.6, 0.6]^3
leaves = 0 ; leaves_in = 0 ; leaves_out = 0
call fig_open(u, 'anim-tut09-octree', '', -60._R8P, 22._R8P)
call fig_anim(u, 1.4_R8P)
call fig_zoom(u, vector_R8P(0._R8P, 0._R8P, 0._R8P), 0.55_R8P)
print '(A)', 'level   side     boxes    cut   near inside outside'
do level = 0, LEVELS
   n = size(lo, 2) ; side = hi(1, 1) - lo(1, 1)
   allocate(codes(n))
   call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.5_R8P * side * sqrt(3._R8P))  ! child diagonal
   print '(I5,F8.4,I10,4I7)', level, side, n, count(codes == BOX_CUT), count(codes == BOX_NEAR), &
                              count(codes == BOX_INSIDE), count(codes == BOX_OUTSIDE)
   leaves_in = leaves_in + count(codes == BOX_INSIDE) ; leaves_out = leaves_out + count(codes == BOX_OUTSIDE)
   if (level == LEVELS) exit
   allocate(next_lo(3, 8 * n), next_hi(3, 8 * n))
   c = 0
   do i = 1, n
      if (codes(i) /= BOX_CUT .and. codes(i) /= BOX_NEAR) cycle            ! inside or outside: a leaf
      do j = 0, 7                                                          ! the eight children
         c = c + 1
         next_lo(:, c) = lo(:, i) + 0.5_R8P * side * [mod(j, 2), mod(j / 2, 2), j / 4]
         next_hi(:, c) = next_lo(:, c) + 0.5_R8P * side
      enddo
   enddo
   call draw_level
   leaves = leaves + n - c / 8
   lo = next_lo(:, 1:c) ; hi = next_hi(:, 1:c)
   deallocate(codes, next_lo, next_hi)
enddo
leaves = leaves + n
print '(A,I0,A,I0,A,I0,A,I0)', 'leaves ', leaves, ' (inside ', leaves_in, ', outside ', leaves_out, &
                               '); uniform grid of the finest side: ', 8**LEVELS
call draw_level
call fig_frame(u, caption='the placed dragon the finest cut boxes wrap', seconds=2._R8P)
call fig_surface(u, 'body', surface)
call fig_close(u)
contains
   subroutine draw_level
   !< One frame: the cut boxes of this level as voxels (only the faces not shared with another cut box).
   logical, allocatable :: cut(:,:,:)
   integer(I4P)         :: m, b, ijk(3), d, s, nb(3), ax, ay, az
   real(R8P)            :: p(3), g

   if (level < 3) return                                            ! levels 0-2: boxes larger than the dragon's parts
   m = nint(1.2_R8P / side) ; g = 0.06_R8P * side
   allocate(cut(0:m + 1, 0:m + 1, 0:m + 1)) ; cut = .false.
   do b = 1, n
      ijk = nint((lo(:, b) + 0.6_R8P) / side) + 1
      if (codes(b) == BOX_CUT) cut(ijk(1), ijk(2), ijk(3)) = .true.
   enddo
   if (level < LEVELS) then
      write(caption, '(A,I0,A,F6.4,A,I0,A,I0,A)') 'level ', level, ', side ', side, ': ', count(codes == BOX_CUT), &
            ' cut + ', count(codes == BOX_NEAR), ' near boxes refined'
   else
      write(caption, '(A,I0,A,F6.4,A,I0,A)') 'level ', level, ', side ', side, ': ', count(codes == BOX_CUT), &
            ' cut boxes, the finest level'
   endif
   call fig_frame(u, caption=trim(caption), seconds=merge(2._R8P, 1.4_R8P, level == LEVELS))
   do b = 1, n
      if (codes(b) /= BOX_CUT) cycle
      ijk = nint((lo(:, b) + 0.6_R8P) / side) + 1
      do d = 1, 3
         do s = -1, 1, 2
            nb = ijk ; nb(d) = nb(d) + s
            if (cut(nb(1), nb(2), nb(3))) cycle
            ax = d ; ay = mod(d, 3) + 1 ; az = mod(d + 1, 3) + 1
            p = lo(:, b) ; if (s > 0) p(ax) = hi(ax, b)
            p = shift(shift(p, ay, g), az, g)                                ! the face inset by g: the cells show
            call fig_quad(u, 'field', v(p), v(shift(p, ay, side - 2 * g)), &
                          v(shift(shift(p, ay, side - 2 * g), az, side - 2 * g)), v(shift(p, az, side - 2 * g)), value=0.55_R8P)
         enddo
      enddo
   enddo
   endsubroutine draw_level

   pure function shift(p, ax, d) result(q)
   real(R8P),    intent(in) :: p(3), d
   integer(I4P), intent(in) :: ax
   real(R8P)                :: q(3)

   q = p ; q(ax) = q(ax) + d
   endfunction shift

   pure function v(p) result(w)
   real(R8P), intent(in) :: p(3)
   type(vector_R8P)      :: w

   w%x = p(1) ; w%y = p(2) ; w%z = p(3)
   endfunction v
endprogram tut09_b_octree
