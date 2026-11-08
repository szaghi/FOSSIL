!< Figure of the reference page "Overset assembly kernels", section "Donor candidates": the cylinder's O-grid over the
!< background of spacing 0.15 and a refined patch of spacing 0.075 over [1.65, 2.4] x [-0.6, 0.6]. donor_candidates_2d keeps every
!< candidate: three vertices of the O-grid's outer ring inside the patch each have two, a background cell (brown) and
!< a patch cell (gold, in the patch drawn in teal), labelled with their donor / receptor volume ratios.
program fig_ref_overset_donors
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, donor_set_t, donor_candidates_2d
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter        :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
integer(I4P), parameter        :: RV(3) = [32_I4P, 1_I4P, 2_I4P]   ! O-grid vertices (i, 13): -11.25, 0, 11.25 deg
type(curvilinear_block_2d_t)   :: blocks(3)
type(topo_record_t)            :: rec
type(donor_set_t), allocatable :: donors(:)
real(R8P)                      :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), px(11, 17), py(11, 17), r
integer(I4P)                   :: i, j, k, c, n, uu, ncand(3)
character(60)                  :: t
character(120)                 :: caption

do j = 1, 13
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
do j = 1, 41
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
do j = 1, 17
   do i = 1, 11
      px(i, j) = 1.65_R8P + 0.075_R8P * (i - 1) ; py(i, j) = -0.6_R8P + 0.075_R8P * (j - 1)
   enddo
enddo
call blocks(3)%initialize(px, py)
do k = 1, 3
   call blocks(k)%build_search
enddo
call donor_candidates_2d(blocks, donors)
ncand = 0
do n = 1, 33 * 13
   k = donors(1)%start(n + 1) - donors(1)%start(n)
   if (k >= 1 .and. k <= 3) ncand(k) = ncand(k) + 1
enddo
print '(A,3I5)', 'O-grid vertices with 1, 2, 3 candidates: ', ncand
write(caption, '(A,I0,A,I0,A)') 'all candidates kept: ', ncand(1), ' O-grid vertices have one, ', ncand(2), &
                                ' two (background and patch)'
call fig_open(uu, 'ref-overset-donors', trim(caption), -90._R8P, 90._R8P, width=640_I4P, height=560_I4P)
call fig_zoom(uu, pt(1.85_R8P, 0._R8P, 0._R8P), 0.7_R8P)
call grid(x, y, 'box', 0._R8P)
call grid(px, py, 'probe', 0.5_R8P * LIFT)
call grid(u, v, 'body', LIFT)
do k = 1, 3
   n = RV(k) + 33 * 12
   t = ''
   do c = donors(1)%start(n), donors(1)%start(n + 1) - 1
      call cell_quad(donors(1)%block(c), donors(1)%cell(:, c))
      write(t, '(A,A,A,F4.2)') trim(t), trim(merge(' background ', ', patch     ', donors(1)%block(c) == 2)), ' ', &
                                  donors(1)%ratio(c)
      print '(A,I0,A,I0,A,2I3,A,2F7.3,A,F7.4)', 'O-grid vertex ', n, ': block ', donors(1)%block(c), ' cell', &
            donors(1)%cell(:, c), ' xi', donors(1)%xi(:, c), ' ratio ', donors(1)%ratio(c)
   enddo
   call fig_points(uu, 'hole', [pt(u(RV(k), 13), v(RV(k), 13), 6 * LIFT)])
   call fig_label(uu, 'accent', pt(u(RV(k), 13), v(RV(k), 13), 7 * LIFT), trim(adjustl(t)))
enddo
call fig_label(uu, 'accent', pt(1.2_R8P, 0.62_R8P, 7 * LIFT), 'donor / receptor volume per candidate')
call fig_close(uu)
contains
   subroutine cell_quad(b, cell)
   !< Fill a donor cell: background brown, patch teal.
   integer(I4P), intent(in) :: b, cell(2) !< Block, cell.
   integer(I4P)             :: a1, a2

   a1 = cell(1) ; a2 = cell(2)
   if (b == 2) then
      call fig_quad(uu, 'cell', pt(x(a1, a2), y(a1, a2), 2 * LIFT), pt(x(a1 + 1, a2), y(a1 + 1, a2), 2 * LIFT), &
                    pt(x(a1 + 1, a2 + 1), y(a1 + 1, a2 + 1), 2 * LIFT), pt(x(a1, a2 + 1), y(a1, a2 + 1), 2 * LIFT))
   else
      call fig_quad(uu, 'foot', pt(px(a1, a2), py(a1, a2), 3 * LIFT), pt(px(a1 + 1, a2), py(a1 + 1, a2), 3 * LIFT), &
                    pt(px(a1 + 1, a2 + 1), py(a1 + 1, a2 + 1), 3 * LIFT), pt(px(a1, a2 + 1), py(a1, a2 + 1), 3 * LIFT))
   endif
   endsubroutine cell_quad

   subroutine grid(gx, gy, style, dz)
   !< The lines of a block.
   real(R8P),    intent(in) :: gx(:,:), gy(:,:) !< Vertices.
   character(*), intent(in) :: style            !< Style.
   real(R8P),    intent(in) :: dz               !< Lift.
   type(vector_R8P), allocatable :: a(:), b(:)
   integer(I4P)             :: n1, n2, m, p, q

   n1 = size(gx, 1) ; n2 = size(gx, 2)
   allocate(a((n1 - 1) * n2 + n1 * (n2 - 1)), b((n1 - 1) * n2 + n1 * (n2 - 1))) ; m = 0
   do q = 1, n2
      do p = 1, n1 - 1
         m = m + 1 ; a(m) = pt(gx(p, q), gy(p, q), dz) ; b(m) = pt(gx(p + 1, q), gy(p + 1, q), dz)
      enddo
   enddo
   do q = 1, n2 - 1
      do p = 1, n1
         m = m + 1 ; a(m) = pt(gx(p, q), gy(p, q), dz) ; b(m) = pt(gx(p, q + 1), gy(p, q + 1), dz)
      enddo
   enddo
   call fig_segments(uu, style, a, b)
   endsubroutine grid

   function pt(qx, qy, dz) result(q)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: qx, qy, dz !< In-plane coordinates, lift.
   type(vector_R8P)      :: q          !< Point.

   q = qx * ex_R8P + qy * ey_R8P + dz * ez_R8P
   endfunction pt
endprogram fig_ref_overset_donors
