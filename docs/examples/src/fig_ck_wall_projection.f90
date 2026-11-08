!run figck-wall-projection fig_ck_wall_projection
!fig ck-wall-projection
!< Figure of the cookbook recipe "Find donors next to a curved wall that two blocks discretize differently": the
!< plane z = 0.5 of ck_wall_projection, seen along z and zoomed on the overlap at the wall: the fine block's grid
!< lines (faint), the coarse block's (stone); the fine wall vertices without a donor before the projection (ember),
!< outside the coarse chords, each joined to its shifted search point on the coarse wall (teal).
program fig_ck_wall_projection
use fossil, only : curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, donor_set_t, donor_candidates, block_real_t, &
                   wall_projection
use fossil_fig, only : fig_open, fig_zoom, fig_points, fig_segments, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter           :: PI = acos(-1._R8P)
type(curvilinear_block_t)      :: blocks(2)
type(donor_set_t), allocatable :: donors(:)
type(block_real_t)             :: shifts(2, 2)
type(vector_R8P), allocatable  :: miss(:), moved(:)
integer(I4P)                   :: i, v, u

call arc(blocks(1), 0._R8P, 0.9_R8P * PI, 41_I4P)            ! fine: 41 vertices along the wall
call arc(blocks(2), 0.62_R8P * PI, 1.02_R8P * PI, 5_I4P)     ! coarse: 5, overlapping the fine one's end
call blocks(1)%build_search ; call blocks(2)%build_search
call donor_candidates(blocks, donors)
allocate(miss(0), moved(0))
do i = 1, 41
   v = i + 41 * 8 + 41 * 9
   if (0.9_R8P * PI * (i - 1) / 40 <= 0.62_R8P * PI + 1e-9_R8P .or. donors(1)%start(v + 1) > donors(1)%start(v)) cycle
   miss = [miss, vector_R8P(blocks(1)%x(i, 9, 2), blocks(1)%y(i, 9, 2), 0.5_R8P)]
enddo
call wall_projection(blocks, 1_I4P, 2_I4P, shifts(1, 2)%v, max_distance=0.3_R8P)
do i = 1, 41
   v = i + 41 * 8 + 41 * 9
   if (0.9_R8P * PI * (i - 1) / 40 <= 0.62_R8P * PI + 1e-9_R8P .or. donors(1)%start(v + 1) > donors(1)%start(v)) cycle
   moved = [moved, vector_R8P(blocks(1)%x(i, 9, 2) + shifts(1, 2)%v(3 * v - 2), &
                              blocks(1)%y(i, 9, 2) + shifts(1, 2)%v(3 * v - 1), 0.5_R8P)]
enddo
call fig_open(u, 'ck-wall-projection', 'fine wall vertices outside the coarse chords (ember), shifted onto them (teal)', &
              0._R8P, 90._R8P, height=420_I4P)
call fig_zoom(u, vector_R8P(1.96_R8P * cos(0.77_R8P * PI), 1.96_R8P * sin(0.77_R8P * PI), 0.5_R8P), 0.17_R8P)
call lines(blocks(1), 'ghost')
call lines(blocks(2), 'body')
call fig_segments(u, 'probe', miss, moved)
call fig_points(u, 'hole', miss)
call fig_points(u, 'inside', moved)
call fig_close(u)
contains
   subroutine lines(blk, style)
   !< The grid lines of the layer k = 2.
   type(curvilinear_block_t), intent(in) :: blk
   character(*),              intent(in) :: style
   type(vector_R8P), allocatable         :: a(:), b(:)
   integer(I4P)                          :: i, j

   allocate(a(0), b(0))
   do j = 1, 9 ; do i = 1, blk%ni
      if (i < blk%ni) then
         a = [a, vector_R8P(blk%x(i, j, 2), blk%y(i, j, 2), 0.5_R8P)]
         b = [b, vector_R8P(blk%x(i + 1, j, 2), blk%y(i + 1, j, 2), 0.5_R8P)]
      endif
      if (j < 9) then
         a = [a, vector_R8P(blk%x(i, j, 2), blk%y(i, j, 2), 0.5_R8P)]
         b = [b, vector_R8P(blk%x(i, j + 1, 2), blk%y(i, j + 1, 2), 0.5_R8P)]
      endif
   enddo ; enddo
   call fig_segments(u, style, a, b)
   endsubroutine lines

   subroutine arc(blk, t0, t1, ni)
   !< A sector of the annulus 1 <= r <= 2 over [t0, t1], 0 <= z <= 1; the wall is the concave r = 2 (j = 9, face 4),
   !< discretized by its own chords.
   type(curvilinear_block_t), intent(inout) :: blk
   real(R8P),                 intent(in)    :: t0, t1
   integer(I4P),              intent(in)    :: ni
   real(R8P)                                :: x(ni, 9, 3), y(ni, 9, 3), z(ni, 9, 3), r, t
   integer(I4P)                             :: i, j, k
   type(topo_record_t)                      :: rec

   do k = 1, 3 ; do j = 1, 9 ; do i = 1, ni
      r = 2._R8P - ((9 - j) / 8._R8P) ** 2 ; t = t0 + (t1 - t0) * (i - 1) / (ni - 1)   ! clustered at the wall
      x(i, j, k) = r * cos(t) ; y(i, j, k) = r * sin(t) ; z(i, j, k) = 0.5_R8P * (k - 1)
   enddo ; enddo ; enddo
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%name = 'wall' ; rec%lo = [1, 9, 1] ; rec%hi = [ni, 9, 3]
   call blk%add_record(rec)
   endsubroutine arc

endprogram fig_ck_wall_projection
