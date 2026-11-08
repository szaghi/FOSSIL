!< Animation of the overset reference page and of the landing page: an O-grid around a unit cylinder crossing a
!< Cartesian background. Each frame moves the O-grid, then recomputes the whole assembly with the fossil_overset 2-D
!< kernels: the cutter (wall_cutter_2d), the background holes (hole_flags_2d), the donor candidates
!< (donor_candidates_2d) and the classification (classify_2d, explicit rule, 3-point stencils). Holes are red dots, the
!< background receptors teal, the O-grid receptors yellow, the orphans ember squares: they appear where the O-grid's
!< outer boundary leaves the background and its receptors find no donor.
program anim_overset
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                       hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, &
                       SELECT_EXPLICIT
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter      :: NF = 16, NI = 33, NJ = 13, MI = 61, MJ = 33
real(R8P),    parameter      :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, H = 0.15_R8P
type(curvilinear_block_2d_t) :: blocks(2)
type(topo_record_t)          :: rec
type(curve_object)           :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(vector_R8P),      allocatable :: a(:), b(:), p(:)
real(R8P)                    :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ), r, cx, cy
integer(I4P)                 :: i, j, m, f, uu, nh, nrb, nro, nor
logical                      :: orphan(NI * NJ + MI * MJ)
character(96)                :: caption

do j = 1, MJ                                                  ! the background, 9 x 4.8
   do i = 1, MI
      x(i, j) = -4.5_R8P + H * (i - 1) ; y(i, j) = -2.4_R8P + H * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(2)%build_search

call fig_open(uu, 'anim-overset', '', -90._R8P, 90._R8P, width=720_I4P, height=440_I4P)
call fig_anim(uu, 0.45_R8P)
allocate(a(MI + MJ), b(MI + MJ)) ; m = 0                       ! the background is drawn in every frame
do j = 1, MJ
   m = m + 1 ; a(m) = pt(x(1, j), y(1, j)) ; b(m) = pt(x(MI, j), y(MI, j))
enddo
do i = 1, MI
   m = m + 1 ; a(m) = pt(x(i, 1), y(i, 1)) ; b(m) = pt(x(i, MJ), y(i, MJ))
enddo
call fig_segments(uu, 'box', a(1:m), b(1:m))
deallocate(a, b)

do f = 0, NF
   cx = -3.2_R8P + 6.4_R8P * f / NF ; cy = 0.5_R8P * sin(PI * f / NF)   ! the cylinder's centre
   do j = 1, NJ                                                ! the O-grid around it, radii 1 to 2
      r = 1._R8P + (j - 1) / real(NJ - 1, R8P)
      do i = 1, NI
         u(i, j) = cx + r * cos(2 * PI * (i - 1) / (NI - 1)) ; v(i, j) = cy + r * sin(2 * PI * (i - 1) / (NI - 1))
      enddo
   enddo
   u(NI,:) = u(1,:) ; v(NI,:) = v(1,:)
   call blocks(1)%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   call blocks(1)%build_search
   call wall_cutter_2d(blocks(1:1), wall)
   call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
   call donor_candidates_2d(blocks, donors)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   nh = count(masks(2)%iota == 0) ; nrb = count(masks(2)%iota < 0) ; nro = count(masks(1)%iota < 0)
   orphan = [masks(1)%orphan > 0, masks(2)%orphan > 0] ; nor = count(orphan)
   write(caption, '(A,I3,A,I3,A,I3,A,I3,A,I3)') 'frame ', f, ': holes ', nh, ', receptors ', nrb, ' + ', nro, &
                                               ', orphans ', nor
   print '(A)', trim(caption)
   call fig_frame(uu, caption=trim(caption))
   allocate(a(2 * NI * NJ), b(2 * NI * NJ)) ; m = 0
   do j = 1, NJ
      do i = 1, NI - 1
         m = m + 1 ; a(m) = pt(u(i, j), v(i, j), LIFT) ; b(m) = pt(u(i + 1, j), v(i + 1, j), LIFT)
      enddo
   enddo
   do i = 1, NI - 1                                            ! radial lines: straight
      m = m + 1 ; a(m) = pt(u(i, 1), v(i, 1), LIFT) ; b(m) = pt(u(i, NJ), v(i, NJ), LIFT)
   enddo
   call fig_segments(uu, 'body', a(1:m), b(1:m))
   deallocate(a, b)
   call fig_points(uu, 'hole', pack([((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)], masks(2)%iota == 0))
   call fig_points(uu, 'probe', pack([((pt(x(i, j), y(i, j), 3 * LIFT), i = 1, MI), j = 1, MJ)], masks(2)%iota < 0))
   call fig_points(uu, 'foot', pack([((pt(u(i, j), v(i, j), 3 * LIFT), i = 1, NI), j = 1, NJ)], masks(1)%iota < 0))
   if (nor > 0) then                                           ! orphans: a square around each
      p = pack([((pt(u(i, j), v(i, j), 4 * LIFT), i = 1, NI), j = 1, NJ), &
                ((pt(x(i, j), y(i, j), 4 * LIFT), i = 1, MI), j = 1, MJ)], orphan)
      allocate(a(4 * nor), b(4 * nor))
      do i = 1, nor
         a(4*i-3) = p(i) + pt(-0.06_R8P, -0.06_R8P) ; b(4*i-3) = p(i) + pt( 0.06_R8P, -0.06_R8P)
         a(4*i-2) = p(i) + pt( 0.06_R8P, -0.06_R8P) ; b(4*i-2) = p(i) + pt( 0.06_R8P,  0.06_R8P)
         a(4*i-1) = p(i) + pt( 0.06_R8P,  0.06_R8P) ; b(4*i-1) = p(i) + pt(-0.06_R8P,  0.06_R8P)
         a(4*i  ) = p(i) + pt(-0.06_R8P,  0.06_R8P) ; b(4*i  ) = p(i) + pt(-0.06_R8P, -0.06_R8P)
      enddo
      call fig_segments(uu, 'patch', a, b)
      call fig_points(uu, 'patch', p)
      deallocate(a, b)
   endif
enddo
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
endprogram anim_overset
