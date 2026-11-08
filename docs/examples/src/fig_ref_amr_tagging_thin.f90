!run figref-amr-tagging-thin fig_ref_amr_tagging_thin
!fig ref-amr-tagging-thin
!< Figure of the reference page "AMR tagging", section "Refining thin walls and gaps": a box 1 x 1 x 0.6 and a plate
!< 0.04 thick standing 0.06 off its +x face. `thickness` and `build_thin`, then `tag_block(thin_factor=1)` on a block
!< of cells of side 0.1: the cells meeting the plate or the face across the gap are tagged; on a block of side 0.03
!< none is, since every local size is above the cell. `box_thin_many` on the coarse cells agrees.
program fig_ref_amr_tagging_thin
use fossil,              only : surface_stl_object, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_zoom, fig_surface, fig_box, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: NF = 6
type(surface_stl_object)        :: body
type(facet_object), allocatable :: f(:)
real(R8P),    allocatable       :: th(:), gp(:), x(:), y(:), z(:), lo(:,:), hi(:,:)
integer(I4P), allocatable       :: tags(:,:,:), ttags(:,:,:)
logical,      allocatable       :: touched(:)
integer(I4P)                    :: i, j, k, m, status, u, nx, ny, nz
character(120)                  :: caption

allocate(f(24 * NF * NF))
call box(f(1:12 * NF * NF),  [0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0.6_R8P])
call box(f(12 * NF * NF + 1:), [1.06_R8P, 0._R8P, 0._R8P], [1.1_R8P, 1._R8P, 0.6_R8P])
call body%adopt_facets(facets=f) ; call body%analyze
call body%thickness(thickness=th, gap=gp, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call body%build_thin(thickness=th, gap=gp, status=status)
if (status /= STATUS_OK) error stop 'build_thin failed'
print '(A,F6.3,A,F6.3)', 'local size min(thickness, gap): plate ', minval(min(th(12 * NF * NF + 1:), gp(12 * NF * NF + 1:))), &
                         ', box face across the gap ', minval(min(th(1:12 * NF * NF), gp(1:12 * NF * NF)))

call coarse_block(0.1_R8P)
print '(A,I0,A,I0,A,I0)', 'cells of side 0.10: ', nx * ny * nz, ', cut ', count(tags == 0), ', thin ', count(ttags == 0)
allocate(lo(3, nx * ny * nz), hi(3, nx * ny * nz), touched(nx * ny * nz))
m = 0
do k = 1, nz ; do j = 1, ny ; do i = 1, nx
   m = m + 1 ; lo(:, m) = [x(i), y(j), z(k)] ; hi(:, m) = [x(i + 1), y(j + 1), z(k + 1)]
enddo ; enddo ; enddo
call body%box_thin_many(lo=lo, hi=hi, factor=1._R8P, touched=touched, status=status)
print '(A,I0,A,L1)', 'box_thin_many on the same cells: ', count(touched), ' touched, same cells as the tags ', &
      all(touched .eqv. reshape(ttags == 0, [nx * ny * nz]))
write(caption, '(A,I0,A)') 'plate 0.04 (ember), gap 0.06; of the cells of side 0.1, ', count(ttags == 0), &
      ' meet a facet thinner than the cell (gold)'
call fig_open(u, 'ref-amr-tagging-thin', trim(caption), -78._R8P, 18._R8P, width=820_I4P, height=520_I4P)
call fig_zoom(u, vector_R8P(1.0_R8P, 0.3_R8P, 0.3_R8P), 0.55_R8P)
call fig_surface(u, 'body', body, last=12 * NF * NF)
call fig_surface(u, 'patch', body, first=12 * NF * NF + 1)
do k = 1, nz ; do j = 1, ny ; do i = 1, nx
   if (ttags(i, j, k) == 0) call fig_box(u, 'accent', vector_R8P(x(i), y(j), z(k)), vector_R8P(x(i + 1), y(j + 1), z(k + 1)))
enddo ; enddo ; enddo
call fig_close(u)
call coarse_block(0.03_R8P)
print '(A,I0,A,I0,A,I0)', 'cells of side 0.03: ', nx * ny * nz, ', cut ', count(tags == 0), ', thin ', count(ttags == 0)
contains
   subroutine coarse_block(h)
   !< A block of cells of side h over [-0.15, 1.25] x [-0.15, 1.15] x [-0.15, 0.75], tagged with band 2 and thin band 1.
   real(R8P), intent(in) :: h
   integer(I4P)          :: ii

   nx = nint(1.4_R8P / h) ; ny = nint(1.3_R8P / h) ; nz = nint(0.9_R8P / h)
   x = [(-0.15_R8P + h * ii, ii = 0, nx)] ; y = [(-0.15_R8P + h * ii, ii = 0, ny)] ; z = [(-0.15_R8P + h * ii, ii = 0, nz)]
   if (allocated(tags)) deallocate(tags, ttags)
   allocate(tags(nx, ny, nz), ttags(nx, ny, nz))
   call body%tag_block(x=x, y=y, z=z, band=2, tags=tags, thin_factor=1._R8P, thin_band=1, thin_tags=ttags, status=status)
   if (status /= STATUS_OK) error stop 'tag_block failed'
   endsubroutine coarse_block

   subroutine box(ff, a, b)
   !< The facets of the box [a, b], outward: each side split into NF x NF quads, two facets each.
   type(facet_object), intent(inout) :: ff(12 * NF * NF)
   real(R8P),          intent(in)    :: a(3), b(3)
   real(R8P)                         :: c(3, 8), p00(3), p10(3), p01(3), e1(3), e2(3)
   integer(I4P), parameter           :: q(4, 6) = reshape([1, 3, 4, 2,  5, 6, 8, 7,  1, 2, 6, 5,  3, 7, 8, 4, &
                                                           1, 5, 7, 3,  2, 4, 8, 6], [4, 6])
   integer(I4P)                      :: ii, s, iu, iv, mm

   do ii = 1, 8
      c(:, ii) = merge(b, a, [btest(ii - 1, 0), btest(ii - 1, 1), btest(ii - 1, 2)])
   enddo
   mm = 0
   do s = 1, 6
      p00 = c(:, q(1, s)) ; p10 = c(:, q(2, s)) ; p01 = c(:, q(4, s))
      e1 = (p10 - p00) / NF ; e2 = (p01 - p00) / NF
      do iv = 0, NF - 1
         do iu = 0, NF - 1
            mm = mm + 1
            call tri(ff(mm), p00 + iu * e1 + iv * e2, p00 + (iu + 1) * e1 + iv * e2, &
                     p00 + (iu + 1) * e1 + (iv + 1) * e2, 0.5_R8P * (a + b))
            mm = mm + 1
            call tri(ff(mm), p00 + iu * e1 + iv * e2, p00 + (iu + 1) * e1 + (iv + 1) * e2, &
                     p00 + iu * e1 + (iv + 1) * e2, 0.5_R8P * (a + b))
         enddo
      enddo
   enddo
   endsubroutine box

   subroutine tri(ft, p1, p2, p3, centre)
   !< A facet oriented away from the centre of its box.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3), centre(3)
   real(R8P)                         :: e1(3), e2(3), nn(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - centre) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri

endprogram fig_ref_amr_tagging_thin
