!< Tutorial 13, step 1: an O-grid around an elliptic cylinder, its topology records (the seam found by
!< detect_topology_2d, the wall and the overlap boundary added by hand), its exact cell geometry, validity and quality,
!< and a folded copy whose two inverted cells validity flags. Figures: the O-grid with its named faces; its cells
!< shaded by aspect ratio; the folded copy around the dragged vertex.
program tut13_a_block
use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_ABUT, TOPO_BOUNDARY, &
                       CELL_VALID, CELL_INVALID
use penf,       only : I4P, R8P
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_quad, fig_label, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter      :: NI = 49, NJ = 17            ! O-grid vertices: around the body, out from the wall
real(R8P),    parameter      :: PI = acos(-1._R8P)
real(R8P),    parameter      :: A = 1._R8P, B = 0.4_R8P      ! semi-axes of the elliptic cylinder
type(curvilinear_block_2d_t) :: ogrid(1), folded       ! a set of one block
type(topo_record_t)          :: rec
real(R8P)                    :: u(NI, NJ), v(NI, NJ), c, mu0, mu, s
integer(I4P)                 :: i, j, r
real(R8P),    allocatable    :: area(:,:), skew(:,:), ortho(:,:), aspect(:,:)
integer(I4P), allocatable    :: valid(:,:)
real(R8P),    parameter      :: LIFT = 1.e-3_R8P
character(100)               :: caption
integer(I4P)                 :: f

c = sqrt(A**2 - B**2) ; mu0 = atanh(B / A)                   ! elliptic coordinates: the wall is mu = mu0
do j = 1, NJ
   s = real(j - 1, R8P) / (NJ - 1)
   mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)   ! layers clustered at the wall
   do i = 1, NI
      u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
      v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
   enddo
enddo
u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)                          ! the seam i = 1 = NI, written exactly
call ogrid(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(ogrid)                            ! the seam: one TOPO_ABUT record per side
do r = 1, ogrid(1)%nrec
   print '(A,I0,A,L1,A,I0,A,2I3)', 'record ', r, ': abut ', ogrid(1)%rec(r)%kind == TOPO_ABUT, ', face ', &
         ogrid(1)%rec(r)%face, ', transform', ogrid(1)%rec(r)%transform(1:2)
enddo
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call ogrid(1)%add_record(rec)                                ! face 3 (-j): the wall
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call ogrid(1)%add_record(rec)                                ! face 4 (+j): the boundary that receives
print '(A,I0,A)', 'records after add_record: ', ogrid(1)%nrec, ' (2 abut, wall, overlap)'
call ogrid(1)%geometry(area=area)                            ! exact cell areas
print '(A,F10.6,A,F10.6)', 'O-grid area ', sum(area), ', annulus between the ellipses ', &
      PI * c**2 * (cosh(mu) * sinh(mu) - cosh(mu0) * sinh(mu0))
call ogrid(1)%validity(valid)                                ! exact: CELL_VALID, CELL_INVALID, CELL_SINGULAR
print '(A,I0,A,I0)', 'valid cells ', count(valid == CELL_VALID), ' of ', size(valid)
call ogrid(1)%quality(skewness=skew, orthogonality=ortho, aspect=aspect)
print '(A,F6.3,A,F6.3,A,F6.2,A,F5.2)', 'max skewness ', maxval(skew), ', min orthogonality ', minval(ortho), &
      ', aspect from ', maxval(aspect), ' at the wall to ', minval(aspect)
folded = ogrid(1)
folded%u(13, 5) = folded%u(13, 7) ; folded%v(13, 5) = folded%v(13, 7)  ! a vertex dragged two layers out
call folded%validity(valid)
print '(A,I0,A)', 'folded copy: ', count(valid == CELL_INVALID), ' invalid cells'
do j = 1, NJ - 1
   do i = 1, NI - 1
      if (valid(i, j) == CELL_INVALID) print '(A,2I3)', '  invalid cell', i, j
   enddo
enddo

! the O-grid and its named faces
call fig_open(f, 'tut13-ogrid', 'O-grid of 48 x 16 cells: wall, overlap boundary and seam records', &
              -90._R8P, 90._R8P, height=460_I4P)
call grid_lines(ogrid(1), 'body', 0._R8P)
call fig_segments(f, 'accent', [(pt(u(i, 1), v(i, 1), LIFT), i = 1, NI - 1)], [(pt(u(i, 1), v(i, 1), LIFT), i = 2, NI)])
call fig_segments(f, 'probe', [(pt(u(i, NJ), v(i, NJ), LIFT), i = 1, NI - 1)], &
                  [(pt(u(i, NJ), v(i, NJ), LIFT), i = 2, NI)])
call fig_segments(f, 'patch', [pt(u(1, 1), v(1, 1), 2 * LIFT)], [pt(u(1, NJ), v(1, NJ), 2 * LIFT)])
call fig_label(f, 'accent', pt(-0.45_R8P, 0._R8P, 3 * LIFT), 'wall, face 3 (j = 1)')
call fig_label(f, 'probe', pt(-1.25_R8P, 1.25_R8P, 3 * LIFT), 'overlap, face 4 (j = 17)')
call fig_label(f, 'patch', pt(1.05_R8P, 0.1_R8P, 3 * LIFT), 'seam: faces 1, 2')
call fig_close(f)

! the cells shaded by aspect ratio
write(caption, '(A,F5.2,A,F4.2)') 'quality: aspect ratio from ', maxval(aspect), ' at the wall (dark) to ', minval(aspect)
call fig_open(f, 'tut13-quality', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(0.9_R8P, 0._R8P), 0.75_R8P)
do j = 1, NJ - 1
   do i = 1, NI - 1
      call fig_quad(f, 'frac', pt(u(i, j), v(i, j)), pt(u(i + 1, j), v(i + 1, j)), pt(u(i + 1, j + 1), v(i + 1, j + 1)), &
                    pt(u(i, j + 1), v(i, j + 1)), value=log(aspect(i, j) / minval(aspect)) / log(maxval(aspect) / minval(aspect)))
   enddo
enddo
call grid_lines(ogrid(1), 'box', LIFT)
call fig_close(f)

! the folded copy around the dragged vertex
call fig_open(f, 'tut13-folded', 'folded copy: cells (12, 5) and (13, 5) inverted, CELL_INVALID', &
              -90._R8P, 90._R8P, height=460_I4P)
call fig_zoom(f, pt(folded%u(13, 5), folded%v(13, 5)), 0.17_R8P)
do j = 1, NJ - 1
   do i = 1, NI - 1
      if (valid(i, j) /= CELL_INVALID) cycle
      call fig_quad(f, 'patch', pt(folded%u(i, j), folded%v(i, j)), pt(folded%u(i + 1, j), folded%v(i + 1, j)), &
                    pt(folded%u(i + 1, j + 1), folded%v(i + 1, j + 1)), pt(folded%u(i, j + 1), folded%v(i, j + 1)))
   enddo
enddo
call grid_lines(folded, 'body', LIFT)
call fig_points(f, 'ghost', [pt(u(13, 5), v(13, 5), 2 * LIFT)])
call fig_points(f, 'accent', [pt(folded%u(13, 5), folded%v(13, 5), 2 * LIFT)])
call fig_label(f, 'ghost', pt(u(13, 5), v(13, 5), 3 * LIFT), '(13, 5) before')
call fig_label(f, 'accent', pt(folded%u(13, 5), folded%v(13, 5), 3 * LIFT), '(13, 5) dragged')
call fig_close(f)
contains
   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt

   subroutine grid_lines(blk, style, dz)
   !< The grid lines of a 2-D block.
   type(curvilinear_block_2d_t), intent(in) :: blk        !< Block.
   character(*),                 intent(in) :: style      !< Style.
   real(R8P),                    intent(in) :: dz         !< Lift.
   type(vector_R8P), allocatable            :: a(:), b(:) !< Segments.
   integer(I4P)                             :: i, j, k    !< Counters.

   allocate(a(2 * blk%ni * blk%nj), b(2 * blk%ni * blk%nj)) ; k = 0
   do j = 1, blk%nj
      do i = 1, blk%ni
         if (i < blk%ni) then
            k = k + 1 ; a(k) = pt(blk%u(i, j), blk%v(i, j), dz) ; b(k) = pt(blk%u(i + 1, j), blk%v(i + 1, j), dz)
         endif
         if (j < blk%nj) then
            k = k + 1 ; a(k) = pt(blk%u(i, j), blk%v(i, j), dz) ; b(k) = pt(blk%u(i, j + 1), blk%v(i, j + 1), dz)
         endif
      enddo
   enddo
   call fig_segments(f, style, a(1:k), b(1:k))
   endsubroutine grid_lines
endprogram tut13_a_block
