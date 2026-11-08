!run tut13d tut13_d_orphans
!fig tut13-overset
!< Tutorial 13, step 4: the assembly of step 3 with the O-grid cut to fewer cell layers. The overlap narrows until the
!< donor stencils of the O-grid's outer vertices reach the background's holes or receptors: orphans, by cause. Figure:
!< the O-grid of 11 layers over the background, with the holes, the receptors of both grids and the orphans.
program tut13_d_orphans
!region orphans
use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY, curve_object, &
                       wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, &
                       overset_mask_t, classify_2d, SELECT_EXPLICIT, ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR
use penf,       only : I4P, R8P
!endregion orphans
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
!region orphans

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ)
integer(I4P)                       :: i, j, k
!endregion orphans
real(R8P),    parameter            :: LIFT = 1.e-3_R8P
type(vector_R8P), allocatable      :: pa(:), pb(:)
character(100)                     :: caption
integer(I4P)                       :: f, m
!region orphans

call make_ogrid(u, v)                                        ! the O-grid of step 1
do j = 1, MJ                                                 ! the background of step 3
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(2)%build_search
do k = 10, 16, 2                                             ! the first k - 1 cell layers of the O-grid only
   call blocks(1)%destroy
   call blocks(1)%initialize(u(:, 1:k), v(:, 1:k), name='thin')
   call detect_topology_2d(blocks(1:1))
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, k, 1] ; rec%hi = [NI, k, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   call blocks(1)%build_search
   call wall_cutter_2d(blocks(1:1), cutter)
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.])
   call donor_candidates_2d(blocks, donors)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   print '(A,I2,A,I2,A,I2,A,I2,A,I2,A)', 'O-grid of ', k - 1, ' layers: receptors ', count(masks(1)%iota < 0), &
         ', orphans ', count(masks(1)%orphan > 0), ' (stencil on a hole ', &
         count(masks(1)%orphan == ORPHAN_STENCIL_HOLE), ', receptor in the stencil ', &
         count(masks(1)%orphan == ORPHAN_RECEPTOR_DONOR), ')'
!endregion orphans
   if (k == 12) call draw                                    ! the figure: 11 layers
!region orphans
enddo
!endregion orphans
contains
   subroutine draw
   !< The assembly of the O-grid of 11 layers.
   write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'O-grid of 11 layers: ', count(masks(1)%iota < 0), ' receptors, ', &
                                       count(masks(1)%orphan > 0), ' orphans, ', count(masks(2)%iota == 0), ' holes'
   call fig_open(f, 'tut13-overset', trim(caption), -90._R8P, 90._R8P, height=460_I4P)
   call fig_zoom(f, pt(0._R8P, 0._R8P), 1.3_R8P)
   allocate(pa(2 * MI * MJ), pb(2 * MI * MJ)) ; m = 0
   do j = 1, MJ
      do i = 1, MI - 1
         m = m + 1 ; pa(m) = pt(x(i, j), y(i, j)) ; pb(m) = pt(x(i + 1, j), y(i + 1, j))
      enddo
   enddo
   do j = 1, MJ - 1
      do i = 1, MI
         m = m + 1 ; pa(m) = pt(x(i, j), y(i, j)) ; pb(m) = pt(x(i, j + 1), y(i, j + 1))
      enddo
   enddo
   call fig_segments(f, 'box', pa(1:m), pb(1:m))
   m = 0
   do j = 1, k
      do i = 1, NI
         if (i < NI) then
            m = m + 1 ; pa(m) = pt(u(i, j), v(i, j), LIFT) ; pb(m) = pt(u(i + 1, j), v(i + 1, j), LIFT)
         endif
         if (j < k) then
            m = m + 1 ; pa(m) = pt(u(i, j), v(i, j), LIFT) ; pb(m) = pt(u(i, j + 1), v(i, j + 1), LIFT)
         endif
      enddo
   enddo
   call fig_segments(f, 'body', pa(1:m), pb(1:m))
   call fig_points(f, 'ghost', pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)], masks(2)%iota == 0))
   call fig_points(f, 'probe', pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)], masks(2)%iota < 0))
   call fig_points(f, 'foot', pack([((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, k)], masks(1)%iota < 0 .and. &
                                   masks(1)%orphan == 0))
   call fig_points(f, 'outside', pack([((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, k)], masks(1)%orphan > 0))
   call fig_close(f)
   deallocate(pa, pb)
   endsubroutine draw

   subroutine make_ogrid(u, v)
   !< The O-grid of step 1: elliptic coordinates about the ellipse of semi-axes 1, 0.4, layers clustered at the wall.
   real(R8P), intent(out) :: u(NI, NJ), v(NI, NJ) !< Vertices.
   real(R8P), parameter   :: PI = acos(-1._R8P)   !< Pi.
   real(R8P)              :: c, mu0, mu, s        !< Coordinates.
   integer(I4P)           :: i, j                 !< Counters.

   c = sqrt(1._R8P - 0.4_R8P**2) ; mu0 = atanh(0.4_R8P)
   do j = 1, NJ
      s = real(j - 1, R8P) / (NJ - 1)
      mu = mu0 + 0.78_R8P * (exp(2._R8P * s) - 1._R8P) / (exp(2._R8P) - 1._R8P)
      do i = 1, NI
         u(i, j) = c * cosh(mu) * cos(2 * PI * (i - 1) / (NI - 1))
         v(i, j) = c * sinh(mu) * sin(2 * PI * (i - 1) / (NI - 1))
      enddo
   enddo
   u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)
   endsubroutine make_ogrid

   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt
endprogram tut13_d_orphans
