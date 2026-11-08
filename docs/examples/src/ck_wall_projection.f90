!run ck_wall_projection ck_wall_projection
program ck_wall_projection
use fossil, only : curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, donor_set_t, donor_candidates, block_real_t, &
                   wall_projection
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter           :: PI = acos(-1._R8P)
type(curvilinear_block_t)      :: blocks(2)
type(donor_set_t), allocatable :: donors(:)
type(block_real_t)             :: shifts(2, 2)

call arc(blocks(1), 0._R8P, 0.9_R8P * PI, 41_I4P)            ! fine: 41 vertices along the wall
call arc(blocks(2), 0.62_R8P * PI, 1.02_R8P * PI, 5_I4P)     ! coarse: 5, overlapping the fine one's end
call blocks(1)%build_search ; call blocks(2)%build_search
call donor_candidates(blocks, donors)
print '(A,I0,A,I0)', 'without projection: fine wall vertices in the overlap with a donor ', found(donors(1)), &
                     ' of ', in_overlap()
call wall_projection(blocks, 1_I4P, 2_I4P, shifts(1, 2)%v, max_distance=0.3_R8P)  ! fine wall onto the coarse one
call donor_candidates(blocks, donors, shifts=shifts)
print '(A,I0,A,I0)', 'with projection:    fine wall vertices in the overlap with a donor ', found(donors(1)), &
                     ' of ', in_overlap()
contains
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

   integer(I4P) function in_overlap() result(n)
   !< Fine wall vertices (layer k = 2) strictly inside the coarse block's angular range.
   integer(I4P) :: i

   n = count([(0.9_R8P * PI * (i - 1) / 40 > 0.62_R8P * PI + 1e-9_R8P, i = 1, 41)])
   endfunction in_overlap

   integer(I4P) function found(d) result(n)
   !< Of those, the ones with a donor candidate.
   type(donor_set_t), intent(in) :: d
   integer(I4P)                  :: i, v

   n = 0
   do i = 1, 41
      if (0.9_R8P * PI * (i - 1) / 40 <= 0.62_R8P * PI + 1e-9_R8P) cycle
      v = i + 41 * 8 + 41 * 9                                ! vertex (i, 9, 2)
      if (d%start(v + 1) > d%start(v)) n = n + 1
   enddo
   endfunction found
endprogram ck_wall_projection
