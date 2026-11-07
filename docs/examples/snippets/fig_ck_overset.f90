!< Figure of the cookbook recipe "Assemble an overset grid": the O-grid around the cylinder and the background, the
!< background holes, the receptors of both blocks.
program fig_ck_overset
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                       hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, &
                       SELECT_EXPLICIT
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter         :: PI = acos(-1._R8P), lift = 1.e-3_R8P
type(curvilinear_block_2d_t) :: blocks(2)
type(topo_record_t)          :: rec
type(curve_object)           :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(vector_R8P),      allocatable :: a(:), b(:)
real(R8P)                    :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                 :: i, j, m, uu

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
rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, 41
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
call donor_candidates_2d(blocks, donors)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)

call fig_open(uu, 'ck-overset', 'O-grid over a background: holes, the receptors around them, the O-grid receptors', &
              -90._R8P, 90._R8P, height=520_I4P)
allocate(a(40 * 41 * 2), b(40 * 41 * 2)) ; m = 0
do j = 1, 41
   do i = 1, 40
      m = m + 1 ; a(m) = pt(x(i, j), y(i, j)) ; b(m) = pt(x(i + 1, j), y(i + 1, j))
      m = m + 1 ; a(m) = pt(x(j, i), y(j, i)) ; b(m) = pt(x(j, i + 1), y(j, i + 1))
   enddo
enddo
call fig_segments(uu, 'box', a(1:m), b(1:m))
deallocate(a, b) ; allocate(a(32 * 13 + 33 * 12), b(32 * 13 + 33 * 12)) ; m = 0
do j = 1, 13
   do i = 1, 32
      m = m + 1 ; a(m) = pt(u(i, j), v(i, j), lift) ; b(m) = pt(u(i + 1, j), v(i + 1, j), lift)
   enddo
enddo
do j = 1, 12
   do i = 1, 33
      m = m + 1 ; a(m) = pt(u(i, j), v(i, j), lift) ; b(m) = pt(u(i, j + 1), v(i, j + 1), lift)
   enddo
enddo
call fig_segments(uu, 'body', a(1:m), b(1:m))
call fig_points(uu, 'hole', pack([((pt(x(i, j), y(i, j), 2 * lift), i = 1, 41), j = 1, 41)], masks(2)%iota == 0))
call fig_points(uu, 'probe', pack([((pt(x(i, j), y(i, j), 3 * lift), i = 1, 41), j = 1, 41)], masks(2)%iota < 0))
call fig_points(uu, 'foot', pack([((pt(u(i, j), v(i, j), 3 * lift), i = 1, 33), j = 1, 13)], masks(1)%iota < 0))
call fig_label(uu, 'accent', pt(2.2_R8P, 2.2_R8P, 4 * lift), 'O-grid receptors')
call fig_label(uu, 'accent', pt(-0.6_R8P, 0._R8P, 4 * lift), 'holes')
call fig_close(uu)
contains
   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt
endprogram fig_ck_overset
