!< Tutorial 13, step 5: the O-grid of step 1 pitching and heaving over the background, assembled at every step with
!< poses: the O-grid, its search tree and its cutter stay in the body frame, built once; hole_flags_2d and
!< donor_candidates_2d take the pose. exposed_points lists the background vertices uncovered by the motion (holes or
!< receptors that become field points: they need values) and counts the holes that became field points without
!< receiving first. The animation shows one period: the O-grid, the holes, the receptors of both grids and the vertices
!< exposed at each step.
program tut13_e_moving
use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY, rigid_pose_2d_t, &
                       pose_2d, to_world, curve_object, wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, &
                       donor_candidates_2d, overset_mask_t, classify_2d, exposed_points, SELECT_EXPLICIT
use penf,       only : I4P, R8P
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_segments, fig_points, fig_close
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33, NF = 12   ! NF steps per period
real(R8P),    parameter            :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(rigid_pose_2d_t)              :: poses(2)
type(block_logical_t), allocatable :: holes(:), exposed(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:), old(:)
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ), angle
integer(I4P)                       :: i, j, k, violations, n
real(R8P),    parameter            :: LIFT = 1.e-3_R8P
type(vector_R8P), allocatable      :: a(:), b(:), w(:)
character(100)                     :: caption
integer(I4P)                       :: f, m

call make_ogrid(u, v)                                        ! the O-grid of step 1, in its body frame
call blocks(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(blocks(1:1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, MJ
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(1)%build_search ; call blocks(2)%build_search    ! once: the O-grid moves by its pose only
call wall_cutter_2d(blocks(1:1), cutter)                     ! once, in the body frame

call fig_open(f, 'tut13-moving', '', -90._R8P, 90._R8P, width=720_I4P, height=460_I4P)
call fig_anim(f, 0.5_R8P)
call fig_zoom(f, pt(0._R8P, 0._R8P), 1.95_R8P)
allocate(a(MI + MJ), b(MI + MJ)) ; m = 0                     ! the background, drawn in every frame
do j = 1, MJ
   m = m + 1 ; a(m) = pt(x(1, j), y(1, j)) ; b(m) = pt(x(MI, j), y(MI, j))
enddo
do i = 1, MI
   m = m + 1 ; a(m) = pt(x(i, 1), y(i, 1)) ; b(m) = pt(x(i, MJ), y(i, MJ))
enddo
call fig_segments(f, 'box', a(1:m), b(1:m))
print '(A)', 'step  pitch  holes  receptors (bg + O-grid)  orphans  exposed  holes to field'
do k = 0, NF
   angle = PI / 6._R8P * sin(2 * PI * mod(k, NF) / NF)     ! pitch 30 degrees, heave 0.3
   poses(1) = pose_2d(angle=angle, translation=[0._R8P, 0.3_R8P * sin(2 * PI * k / NF + PI / 2._R8P)])
   poses(2) = pose_2d(angle=0._R8P)                          ! the background stays
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   if (k == 0) old = masks
   call exposed_points(old, masks, exposed, violations)      ! uncovered since the previous step
   print '(I4,F7.1,I7,I8,A,I3,I15,I9,I12)', k, angle * 180._R8P / PI, count(masks(2)%iota == 0), &
         count(masks(2)%iota < 0), ' +', count(masks(1)%iota < 0), count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0), &
         count(exposed(2)%v), violations
   if (k < NF) then
      write(caption, '(A,I0,A,F6.1,A,I0,A,I0)') 'step ', k, ', pitch ', angle * 180._R8P / PI, ' deg: holes ', &
           count(masks(2)%iota == 0), ', exposed ', count(exposed(2)%v)
      call fig_frame(f, caption=trim(caption))
      if (allocated(a)) deallocate(a, b)
      allocate(a(2 * NI * NJ), b(2 * NI * NJ)) ; m = 0       ! the O-grid at this pose
      do j = 1, NJ
         do i = 1, NI
            if (i < NI) then
               m = m + 1 ; a(m) = wpt(i, j, LIFT) ; b(m) = wpt(i + 1, j, LIFT)
            endif
            if (j < NJ .and. mod(i - 1, 2) == 0) then
               m = m + 1 ; a(m) = wpt(i, j, LIFT) ; b(m) = wpt(i, j + 1, LIFT)
            endif
         enddo
      enddo
      call fig_segments(f, 'ghost', a(1:m), b(1:m))
      w = [((pt(x(i, j), y(i, j), 2 * LIFT), i = 1, MI), j = 1, MJ)]
      call fig_points(f, 'body', pack(w, masks(2)%iota == 0))
      call fig_points(f, 'probe', pack(w, masks(2)%iota < 0))
      call fig_points(f, 'patch', pack(w, exposed(2)%v))
      w = [((wpt(i, j, 2 * LIFT), i = 1, NI), j = 1, NJ)]
      call fig_points(f, 'foot', pack(w, masks(1)%iota < 0))
      call fig_points(f, 'hole', pack(w, masks(1)%orphan > 0))
   endif
   old = masks
enddo
n = 0
do k = 0, 2 * NF                                             ! the same period in steps half as long
   poses(1) = pose_2d(angle=PI / 6._R8P * sin(PI * mod(k, 2 * NF) / NF), &
                      translation=[0._R8P, 0.3_R8P * sin(PI * k / NF + PI / 2._R8P)])
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   if (k > 0) call exposed_points(old, masks, exposed, violations)
   if (k > 0) n = n + violations
   old = masks
enddo
print '(A,I0,A,I0)', 'steps per period ', 2 * NF, ': holes to field ', n
call fig_close(f)
contains
   subroutine make_ogrid(u, v)
   !< The O-grid of step 1: elliptic coordinates about the ellipse of semi-axes 1, 0.4, layers clustered at the wall.
   real(R8P), intent(out) :: u(NI, NJ), v(NI, NJ) !< Vertices.
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

   function wpt(i, j, dz) result(q)
   !< O-grid vertex (i, j) at the current pose.
   integer(I4P), intent(in) :: i, j !< Vertex.
   real(R8P),    intent(in) :: dz   !< Lift.
   type(vector_R8P)         :: q    !< Point.
   real(R8P)                :: p(2) !< World coordinates.

   p = to_world(poses(1), [u(i, j), v(i, j)])
   q = pt(p(1), p(2), dz)
   endfunction wpt
endprogram tut13_e_moving
