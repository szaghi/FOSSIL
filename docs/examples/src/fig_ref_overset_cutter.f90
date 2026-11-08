!run figref-overset-cutter fig_ref_overset_cutter
!fig ref-overset-cutter
!< Figure of the reference page "Overset assembly kernels", section "Cutters": two partial annuli around the unit
!< cylinder, a fine one (24 cells over 0 to 200 degrees) and a coarse one (11 cells over 153 to 373 degrees), whose
!< wall patches overlap at both ends. wall_cutter_2d keeps the fine patch, drops the coarse segments it covers (both
!< ends and the midpoint within the tolerance: ghost) and keeps the others: double-covered strips at 193 to 200 and
!< 0 to 13 degrees, no gap (issue #91; the centroid rule, strict=.false., dropped the segment from 353 to 373 and left
!< 353 to 360 uncovered); every kept segment carries a tick towards the fluid (solid on the left of a -> b).
program fig_ref_overset_cutter
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d
use fossil_fig, only : fig_open, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter        :: PI = acos(-1._R8P), D2R = PI / 180._R8P, LIFT = 1.e-3_R8P
type(curvilinear_block_2d_t)   :: blocks(2)
type(curve_object)             :: cutter
logical,          allocatable  :: shared(:,:), kept(:)
integer(I4P),     allocatable  :: owner(:)
type(vector_R8P), allocatable  :: a(:), b(:), ta(:), tb(:)
real(R8P)                      :: p(2), q(2), mid(2), nrm(2), dbl, gap, lo
integer(I4P)                   :: u, i, j, l, first, last, nk, ns
character(120)                 :: caption

call annulus(blocks(1), 0._R8P, 200._R8P, 24_I4P, 1.5_R8P)
call annulus(blocks(2), 153._R8P, 373._R8P, 11_I4P, 1.7_R8P)
call wall_cutter_2d(blocks, cutter, shared=shared, owner=owner)
print '(A,I0,A,L1)', 'cutter loops: ', cutter%get_loops_number(), ', shared wall: ', shared(1, 2)
do l = 1, cutter%get_loops_number()
   call cutter%get_loop_range(l, first, last)
   print '(A,I0,A,I0,A,I0)', 'loop ', l, ' from block ', owner(l), ', vertices ', last - first + 1
enddo
! which of the coarse wall segments (block 2, j = 1) survive: a cutter segment with the same end points
allocate(kept(11)) ; kept = .false.
do i = 1, 11
   p = [blocks(2)%u(i, 1), blocks(2)%v(i, 1)] ; q = [blocks(2)%u(i + 1, 1), blocks(2)%v(i + 1, 1)]
   do j = 1, cutter%get_vertices_number()
      if (cutter%get_next(j) == 0) cycle
      if (all(abs(cutter%get_vertex(j) - q) < 1.e-12_R8P) .and. &
          all(abs(cutter%get_vertex(cutter%get_next(j)) - p) < 1.e-12_R8P)) kept(i) = .true.
      if (all(abs(cutter%get_vertex(j) - p) < 1.e-12_R8P) .and. &
          all(abs(cutter%get_vertex(cutter%get_next(j)) - q) < 1.e-12_R8P)) kept(i) = .true.
   enddo
enddo
print '(A,11L2)', 'coarse wall segments kept: ', kept
nk = count(kept)
! the coarse range 153-373 deg against the fine 0-200 (360-560): covered twice where a kept coarse segment lies over
! the fine patch, not at all where a dropped one lies off it
dbl = 0._R8P ; gap = 0._R8P
do i = 1, 11
   lo = 153._R8P + 20 * (i - 1)
   if (kept(i)) dbl = dbl + max(0._R8P, min(lo + 20, 200._R8P) - max(lo, 0._R8P)) &
                          + max(0._R8P, min(lo + 20, 560._R8P) - max(lo, 360._R8P))
   if (.not. kept(i)) gap = gap + 20 - max(0._R8P, min(lo + 20, 200._R8P) - max(lo, 0._R8P)) &
                                     - max(0._R8P, min(lo + 20, 560._R8P) - max(lo, 360._R8P))
enddo
print '(A,F4.1,A,F3.1,A)', 'double-covered: ', dbl, ' deg, uncovered: ', gap, ' deg'

write(caption, '(A,I0,A,I0,A)') 'fine patch kept (gold), coarse patch (ember): ', 11 - nk, &
                                ' dropped (grey, drawn inside), ', nk, ' kept (drawn outside)'
call fig_open(u, 'ref-overset-cutter', trim(caption), -90._R8P, 90._R8P, width=820_I4P, height=600_I4P)
call grid(blocks(1), 'box', 0._R8P)
call grid(blocks(2), 'cell', 0._R8P)
ns = 0
allocate(a(0), b(0), ta(0), tb(0))
do j = 1, cutter%get_vertices_number()
   if (cutter%get_next(j) == 0) cycle
   p = cutter%get_vertex(j) ; q = cutter%get_vertex(cutter%get_next(j))
   mid = 0.5_R8P * (p + q) ; nrm = [q(2) - p(2), p(1) - q(1)] / norm2(q - p) ! right of a -> b: the fluid
   ta = [ta, pt(mid, 3 * LIFT)] ; tb = [tb, pt(mid + 0.12_R8P * nrm, 3 * LIFT)]
enddo
do l = 1, cutter%get_loops_number()
   call cutter%get_loop_range(l, first, last)
   if (allocated(a)) deallocate(a, b)
   allocate(a(0), b(0))
   do j = first, last
      if (cutter%get_next(j) == 0) cycle
      a = [a, pt(cutter%get_vertex(j) * merge(1._R8P, 1.03_R8P, owner(l) == 1), (2 + owner(l)) * LIFT)]
      b = [b, pt(cutter%get_vertex(cutter%get_next(j)) * merge(1._R8P, 1.03_R8P, owner(l) == 1), (2 + owner(l)) * LIFT)]
   enddo
   call fig_segments(u, trim(merge('accent', 'patch ', owner(l) == 1)), a, b)
enddo
call fig_segments(u, 'probe', ta, tb)
deallocate(a, b) ; allocate(a(0), b(0))
do i = 1, 11
   if (kept(i)) cycle
   a = [a, pt([blocks(2)%u(i, 1), blocks(2)%v(i, 1)] * 0.95_R8P, 2 * LIFT)]
   b = [b, pt([blocks(2)%u(i + 1, 1), blocks(2)%v(i + 1, 1)] * 0.95_R8P, 2 * LIFT)]
enddo
call fig_segments(u, 'ghost', a, b)
call fig_label(u, 'accent', pt([0.2_R8P, 1.62_R8P], 5 * LIFT), 'fine, 0-200 deg')
call fig_label(u, 'patch', pt([-0.3_R8P, -1.9_R8P], 5 * LIFT), 'coarse, 153-373 deg')
call fig_label(u, 'ghost', pt(0.72_R8P * [cos(172 * D2R), sin(172 * D2R)], 5 * LIFT), 'dropped')
call fig_label(u, 'hole', pt(0.8_R8P * [cos(197 * D2R), sin(197 * D2R)], 5 * LIFT), 'double, 193-200')
call fig_label(u, 'hole', pt(1.12_R8P * [cos(6 * D2R), sin(6 * D2R)], 5 * LIFT), 'double, 0-13')
call fig_close(u)
contains
   subroutine annulus(blk, t0, t1, nc, rmax)
   !< A partial annulus from angle t0 to t1 (degrees), nc cells around, 4 cells out to rmax, wall on face 3.
   type(curvilinear_block_2d_t), intent(inout) :: blk         !< Block.
   real(R8P),                    intent(in)    :: t0, t1, rmax !< Angles, outer radius.
   integer(I4P),                 intent(in)    :: nc          !< Cells around.
   real(R8P)                                   :: bu(nc + 1, 5), bv(nc + 1, 5), t, r
   type(topo_record_t)                         :: rec
   integer(I4P)                                :: ii, jj

   do jj = 1, 5
      r = 1._R8P + (rmax - 1._R8P) * (jj - 1) / 4._R8P
      do ii = 1, nc + 1
         t = (t0 + (t1 - t0) * (ii - 1) / nc) * D2R
         bu(ii, jj) = r * cos(t) ; bv(ii, jj) = r * sin(t)
      enddo
   enddo
   call blk%initialize(bu, bv)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [nc + 1, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine annulus

   subroutine grid(blk, style, dz)
   !< The lines of a block.
   type(curvilinear_block_2d_t), intent(in) :: blk   !< Block.
   character(*),                 intent(in) :: style !< Style.
   real(R8P),                    intent(in) :: dz    !< Lift.
   type(vector_R8P), allocatable :: ga(:), gb(:)
   integer(I4P)                  :: ii, jj

   allocate(ga(0), gb(0))
   do jj = 1, blk%nj
      do ii = 1, blk%ni - 1
         ga = [ga, pt([blk%u(ii, jj), blk%v(ii, jj)], dz)] ; gb = [gb, pt([blk%u(ii + 1, jj), blk%v(ii + 1, jj)], dz)]
      enddo
   enddo
   do jj = 1, blk%nj - 1
      do ii = 1, blk%ni
         ga = [ga, pt([blk%u(ii, jj), blk%v(ii, jj)], dz)] ; gb = [gb, pt([blk%u(ii, jj + 1), blk%v(ii, jj + 1)], dz)]
      enddo
   enddo
   call fig_segments(u, style, ga, gb)
   endsubroutine grid

   function pt(xy, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: xy(2) !< In-plane coordinates.
   real(R8P), intent(in) :: dz    !< Lift.
   type(vector_R8P)      :: r     !< Point.

   r = xy(1) * ex_R8P + xy(2) * ey_R8P + dz * ez_R8P
   endfunction pt
endprogram fig_ref_overset_cutter
