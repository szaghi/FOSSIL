!< Tutorial 9, step 1: classify the 8^3 blocks of the domain against the placed dragon (box_classify_many), and ask one
!< block alone (box_intersects_surface). The figure: the dragon with its cut blocks (ember) and near blocks (teal).
program tut09_a_blocks
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use vecfor,     only : vector_R8P
use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nb = 8                      ! blocks per direction
type(surface_stl_object) :: surface
real(R8P)                :: lo(3, nb**3), hi(3, nb**3), hb
integer(I4P)             :: codes(nb**3), i, j, k, b, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
hb = 1.2_R8P / nb                                       ! the domain [-0.6, 0.6]^3 cut into 8^3 blocks of side 0.15
b = 0
do k = 1, nb ; do j = 1, nb ; do i = 1, nb
   b = b + 1
   lo(:, b) = -0.6_R8P + [i - 1, j - 1, k - 1] * hb
   hi(:, b) = lo(:, b) + hb
enddo ; enddo ; enddo
call surface%box_classify_many(lo=lo, hi=hi, codes=codes)                      ! cut, or inside / outside
print '(A,I0,A,I0,A,I0)', 'radius 0:    cut ', count(codes == BOX_CUT), ', inside ', count(codes == BOX_INSIDE), &
                          ', outside ', count(codes == BOX_OUTSIDE)
call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.05_R8P)     ! and near within 0.05
print '(A,I0,A,I0,A,I0,A,I0)', 'radius 0.05: cut ', count(codes == BOX_CUT), ', near ', count(codes == BOX_NEAR), &
                               ', inside ', count(codes == BOX_INSIDE), ', outside ', count(codes == BOX_OUTSIDE)
b = findloc(codes, BOX_CUT, dim=1)                                             ! the first cut block
print '(A,I0,A,L1)', 'block ', b, ' meets the dragon: ', surface%box_intersects_surface(lo(:, b), hi(:, b))
print '(A,L1)', 'its lower half meets it:  ', &
      surface%box_intersects_surface(lo(:, b), [hi(1, b), hi(2, b), lo(3, b) + 0.5_R8P * hb])
call fig_open(u, 'tut09-blocks', 'box_classify_many, 8^3 blocks: cut (ember), near within 0.05 (teal)', -60._R8P, 22._R8P)
call fig_surface(u, 'body', surface)
do b = 1, nb**3
   if (codes(b) == BOX_CUT) call fig_box(u, 'outside', vector_R8P(lo(1, b), lo(2, b), lo(3, b)), &
                                         vector_R8P(hi(1, b), hi(2, b), hi(3, b)))
   if (codes(b) == BOX_NEAR) call fig_box(u, 'probe', vector_R8P(lo(1, b), lo(2, b), lo(3, b)), &
                                          vector_R8P(hi(1, b), hi(2, b), hi(3, b)))
enddo
call fig_close(u)
endprogram tut09_a_blocks
