!< Figure of the reference page "AMR tagging", section "How it works": `tag_block` on a 24 x 18 block around a square
!< ring with an island in its hole. Every cell carries its tag (the Chebyshev cell distance to the nearest cut cell,
!< capped at band + 1 = 3, negative inside); the uncut cells fall into 8-connected components, each signed by one
!< point query at the centre of its first cell (column-major order), outlined in gold.
program fig_ref_amr_tagging_bands
use fossil,     only : curve_object, STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: NX = 24, NY = 18, BAND = 2
type(curve_object)        :: body
type(vector_R8P), allocatable :: seeds(:)
real(R8P)                 :: x(NX + 1), y(NY + 1)
integer(I4P)              :: tags(NX, NY), comp(NX, NY), i, j, ncomp, status, u
character(8)              :: txt
character(120)            :: caption

! outer square (counter-clockwise), hole (clockwise), island (counter-clockwise), off the grid lines
call body%adopt_loops(u=[-1.03_R8P, 1.03_R8P, 1.03_R8P, -1.03_R8P,  -0.61_R8P, -0.61_R8P, 0.61_R8P, 0.61_R8P, &
                         -0.22_R8P, 0.22_R8P, 0.22_R8P, -0.22_R8P], &
                      v=[-0.83_R8P, -0.83_R8P, 0.83_R8P, 0.83_R8P,  -0.47_R8P, 0.47_R8P, 0.47_R8P, -0.47_R8P, &
                         -0.16_R8P, -0.16_R8P, 0.16_R8P, 0.16_R8P], loop_sizes=[4_I4P, 4_I4P, 4_I4P], status=status)
if (status /= STATUS_OK) error stop 'adopt_loops failed'
x = [(-1.5_R8P + 0.125_R8P * i, i = 0, NX)] ; y = [(-1.125_R8P + 0.125_R8P * i, i = 0, NY)]
call body%tag_block(x=x, y=y, band=BAND, tags=tags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
do i = -(BAND + 1), BAND + 1
   print '(A,I3,A,I4)', 'tag', i, ': cells', count(tags == i)
enddo
call components
print '(A,I0,A)', 'uncut cells in ', ncomp, ' components, one point query each'

write(caption, '(A,I0,A,I0,A)') 'tag_block, band ', BAND, ' (cap 3): ', ncomp, &
                                ' components of uncut cells, one point query each, in the outlined cell'
call fig_open(u, 'ref-amr-tagging-bands', trim(caption), -90._R8P, 90._R8P, width=820_I4P, height=660_I4P)
do j = 1, NY
   do i = 1, NX
      call fig_quad(u, 'field', p(x(i), y(j)), p(x(i + 1), y(j)), p(x(i + 1), y(j + 1)), p(x(i), y(j + 1)), &
                    sign(0.3_R8P + 0.7_R8P * (abs(tags(i, j)) - 1) / BAND, real(tags(i, j), R8P)) &
                    * merge(0._R8P, 1._R8P, tags(i, j) == 0))
      if (tags(i, j) == 0) cycle
      write(txt, '(I0)') tags(i, j)
      call fig_label(u, 'ink', p(x(i) + 0.012_R8P, y(j) + 0.012_R8P, 3.e-3_R8P), trim(txt))
   enddo
enddo
do i = 0, NX
   call fig_segment(u, 'box', p(x(i + 1), y(1), 1.e-3_R8P), p(x(i + 1), y(NY + 1), 1.e-3_R8P))
enddo
do j = 0, NY
   call fig_segment(u, 'box', p(x(1), y(j + 1), 1.e-3_R8P), p(x(NX + 1), y(j + 1), 1.e-3_R8P))
enddo
call draw_curve
call fig_segments(u, 'foot', [(seeds(i), i = 1, size(seeds))], [(seeds(4 * ((i - 1) / 4) + mod(i, 4) + 1), &
                  i = 1, size(seeds))])
call fig_close(u)
contains
   subroutine components
   !< The 8-connected components of the uncut cells (flood fill), and the centre of the first cell of each.
   integer(I4P) :: stk(2, NX * NY), sp, a, b, di, dj, ii, jj

   comp = 0 ; ncomp = 0 ; allocate(seeds(0))
   do jj = 1, NY
      do ii = 1, NX
         if (tags(ii, jj) == 0 .or. comp(ii, jj) /= 0) cycle
         ncomp = ncomp + 1 ; comp(ii, jj) = ncomp
         seeds = [seeds, p(x(ii), y(jj), 4.e-3_R8P), p(x(ii + 1), y(jj), 4.e-3_R8P), p(x(ii + 1), y(jj + 1), 4.e-3_R8P), &
                  p(x(ii), y(jj + 1), 4.e-3_R8P)]
         sp = 1 ; stk(:, 1) = [ii, jj]
         do while (sp > 0)
            a = stk(1, sp) ; b = stk(2, sp) ; sp = sp - 1
            do dj = -1, 1
               do di = -1, 1
                  if (a + di < 1 .or. a + di > NX .or. b + dj < 1 .or. b + dj > NY) cycle
                  if (tags(a + di, b + dj) == 0 .or. comp(a + di, b + dj) /= 0) cycle
                  comp(a + di, b + dj) = ncomp ; sp = sp + 1 ; stk(:, sp) = [a + di, b + dj]
               enddo
            enddo
         enddo
         print '(A,I0,A,I4,A,I3,A,I3,A)', 'component ', ncomp, ': ', count(comp == ncomp), ' cells, first cell (', &
               ii, ',', jj, '), ' // trim(merge('inside ', 'outside', tags(ii, jj) < 0))
      enddo
   enddo
   endsubroutine components

   subroutine draw_curve
   !< The loops, just above the plane.
   real(R8P)    :: a(2), b(2)
   integer(I4P) :: s

   do s = 1, body%get_vertices_number()
      a = body%get_vertex(s) ; b = body%get_vertex(body%get_next(s))
      call fig_segment(u, 'accent', p(a(1), a(2), 2.e-3_R8P), p(b(1), b(2), 2.e-3_R8P))
   enddo
   endsubroutine draw_curve

   function p(px, py, dz) result(q)
   !< The point (px, py, dz) of the figure, dz = 0 by default.
   real(R8P), intent(in)           :: px, py
   real(R8P), intent(in), optional :: dz
   type(vector_R8P)                :: q

   q = vector_R8P(px, py, 0._R8P)
   if (present(dz)) q%z = dz
   endfunction p
endprogram fig_ref_amr_tagging_bands
