!run figck-overset-3d fig_ck_overset_3d
!fig ck-overset-3d
!< Figure of the cookbook recipe "Assemble a 3-D overset grid": the plane y = 0 of ck_overset_3d's assembly, seen
!< along y: the shell's grid lines in that plane (stone), the background vertices of the plane, holes in teal, receptors
!< in ember, field vertices faint; the shell's receptors (its outer ring, OUTER 1) in ember too.
program fig_ck_overset_3d
use fossil, only : surface_stl_object, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, block_logical_t, &
                   wall_cutter, hole_flags, donor_set_t, donor_candidates, overset_mask_t, classify, fringe_levels, &
                   SELECT_EXPLICIT, FRINGE_INNER_1, FRINGE_OUTER_1
use fossil_fig, only : fig_open, fig_points, fig_segments, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_t)          :: blocks(2)
type(topo_record_t)                :: rec
type(surface_stl_object)           :: cutter
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
real(R8P)                          :: x(25, 13, 6), y(25, 13, 6), z(25, 13, 6), r, th, ph, g(17)
real(R8P)                          :: bx(17, 17, 17), by(17, 17, 17), bz(17, 17, 17)
integer(I4P)                       :: i, j, k, n, u
type(vector_R8P), allocatable      :: bg(:), bgh(:), bgr(:), sa(:), sb(:), sr(:)
character(100)                     :: caption

do k = 1, 6 ; do j = 1, 13 ; do i = 1, 25                   ! block 1: a shell 1 <= r <= 2 around a sphere
   r = 1._R8P + 0.2_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call blocks(1)%initialize(x, y, z)
rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [25, 13, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)                              ! k = 1: the wall
rec%face = 6 ; rec%lo = [1, 1, 6] ; rec%hi = [25, 13, 6] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)                              ! k = 6: receives from the background
g = [(-2.4_R8P + 0.3_R8P * (i - 1), i = 1, 17)]             ! block 2: a Cartesian background, h = 0.3
do k = 1, 17 ; do j = 1, 17
   bx(:, j, k) = g ; by(:, j, k) = g(j) ; bz(:, j, k) = g(k)
enddo ; enddo
call blocks(2)%initialize(bx, by, bz)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter(blocks(1:1), cutter)
call hole_flags(blocks, cutter, holes, cuts=[.false., .true.])
call donor_candidates(blocks, donors)                       ! every vertex located in the other block
call classify(blocks, donors, masks, SELECT_EXPLICIT, 2_I4P, holes=holes)   ! stencils of width 2
call fringe_levels(blocks, masks)
allocate(bg(0), bgh(0), bgr(0), sa(0), sb(0), sr(0))
do k = 1, 17 ; do i = 1, 17                                 ! the background plane j = 9 (y = 0)
   n = i + 17 * 8 + 17 * 17 * (k - 1)
   if (masks(2)%iota(n) == 0) then ; bgh = [bgh, vector_R8P(bx(i, 9, k), 0._R8P, bz(i, 9, k))]
   elseif (masks(2)%iota(n) < 0) then ; bgr = [bgr, vector_R8P(bx(i, 9, k), 0._R8P, bz(i, 9, k))]
   else ; bg = [bg, vector_R8P(bx(i, 9, k), 0._R8P, bz(i, 9, k))] ; endif
enddo ; enddo
do i = 1, 13, 12                                            ! the shell's meridians th = 0 and th = pi
   do k = 1, 6 ; do j = 1, 13
      if (j < 13) then
         sa = [sa, vector_R8P(x(i, j, k), 0._R8P, z(i, j, k))] ; sb = [sb, vector_R8P(x(i, j + 1, k), 0._R8P, z(i, j + 1, k))]
      endif
      if (k < 6) then
         sa = [sa, vector_R8P(x(i, j, k), 0._R8P, z(i, j, k))] ; sb = [sb, vector_R8P(x(i, j, k + 1), 0._R8P, z(i, j, k + 1))]
      endif
      n = i + 25 * (j - 1) + 25 * 13 * (k - 1)
      if (masks(1)%iota(n) < 0) sr = [sr, vector_R8P(x(i, j, k), 0._R8P, z(i, j, k))]
   enddo ; enddo
enddo
write(caption, '(A,I0,A,I0,A)') 'plane y = 0: holes (teal, ', count(masks(2)%iota == 0), ' in all) and receptors (ember, ', &
                                count(masks(2)%iota < 0) + count(masks(1)%iota < 0), ')'
call fig_open(u, 'ck-overset-3d', trim(caption), -90._R8P, 0._R8P, height=460_I4P)
call fig_segments(u, 'body', sa, sb)
call fig_points(u, 'ghost', bg)
call fig_points(u, 'inside', bgh)
call fig_points(u, 'outside', bgr)
call fig_points(u, 'outside', sr)
call fig_close(u)
endprogram fig_ck_overset_3d
