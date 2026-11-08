!< Figure of the reference page "Overset assembly kernels", section "Moving grids": the cylinder's O-grid moves with a
!< pose (vertices and cutter in the body frame, search built once, `previous=` walks from the last donors) over the
!< background of spacing h = 0.15, by 0.6 h (left) and by 3.6 h (right) in one step. Background vertices: holes and
!< receptors after the step; exposed points (holes or receptors before, field points after) in gold; the holes that
!< became field points without being receptors first (Meakin & Suhs violations) in ember.
program fig_ref_overset_exposed
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                       hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, &
                       SELECT_EXPLICIT, exposed_points, rigid_pose_2d_t, pose_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter               :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, H = 0.15_R8P, SH = 4.6_R8P
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: wall
type(donor_set_t),     allocatable :: d0(:), d1(:)
type(overset_mask_t),  allocatable :: m0(:), m1(:)
type(block_logical_t), allocatable :: exposed(:)
type(vector_R8P)                   :: p(41 * 41)
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r, step(2)
integer(I4P)                       :: i, j, uu, run, viol
character(120)                     :: caption

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
      x(i, j) = -3._R8P + H * (i - 1) ; y(i, j) = -3._R8P + H * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search   ! once: the O-grid stays in its body frame
call wall_cutter_2d(blocks(1:1), wall)
call assemble(0._R8P, m0, d0)
step = [0.6_R8P * H, 3.6_R8P * H]

caption = ''
do run = 1, 2
   call assemble(step(run), m1, d1, d0)
   call exposed_points(m0, m1, exposed, viol)
   print '(A,F5.1,A,I0,A,I0,A,I0,A,I0)', 'step ', step(run) / H, ' h: holes ', count(m1(2)%iota == 0), &
         ', receptors ', count(m1(2)%iota < 0), ', exposed ', count(exposed(2)%v), ', violations ', viol
   if (run == 2) caption = trim(caption)//';'
   write(caption, '(A,1X,F3.1,A,I0,A,I0,A)') trim(caption), step(run) / H, ' h: ', count(exposed(2)%v), &
                                             ' exposed, ', viol, ' holes to field'
enddo
call fig_open(uu, 'ref-overset-exposed', 'steps of '//trim(adjustl(caption)), -90._R8P, 90._R8P, width=900_I4P, &
              height=470_I4P)
do run = 1, 2
   call assemble(step(run), m1, d1, d0)
   call exposed_points(m0, m1, exposed, viol)
   p = [((pt(x(i, j) + SH * (run - 1), y(i, j), 2 * LIFT), i = 1, 41), j = 1, 41)]
   call fig_segments(uu, 'box', pack(p, [((abs(x(i, j)) <= 2.2_R8P .and. abs(y(i, j)) <= 2.2_R8P .and. &
                     i < 41, i = 1, 41), j = 1, 41)]), pack([((pt(x(i, j) + H + SH * (run - 1), y(i, j)), i = 1, 41), &
                     j = 1, 41)], [((abs(x(i, j)) <= 2.2_R8P .and. abs(y(i, j)) <= 2.2_R8P .and. i < 41, i = 1, 41), &
                     j = 1, 41)]))
   call fig_segments(uu, 'box', pack(p, [((abs(x(i, j)) <= 2.2_R8P .and. abs(y(i, j)) <= 2.2_R8P .and. &
                     j < 41, i = 1, 41), j = 1, 41)]), pack([((pt(x(i, j) + SH * (run - 1), y(i, j) + H), i = 1, 41), &
                     j = 1, 41)], [((abs(x(i, j)) <= 2.2_R8P .and. abs(y(i, j)) <= 2.2_R8P .and. j < 41, i = 1, 41), &
                     j = 1, 41)]))
   call circle(0._R8P + SH * (run - 1), 'ghost')
   call circle(step(run) + SH * (run - 1), 'hole')
   call fig_label(uu, 'ghost', pt(-1.45_R8P + SH * (run - 1), -1.15_R8P, 5 * LIFT), 'n')
   call fig_label(uu, 'hole', pt(0.75_R8P + step(run) + SH * (run - 1), -1.15_R8P, 5 * LIFT), 'n + 1')
   call fig_points(uu, 'cell', pack(p, m1(2)%iota == 0))
   call fig_points(uu, 'probe', pack(p, m1(2)%iota < 0))
   call fig_points(uu, 'foot', pack(p, exposed(2)%v .and. m0(2)%iota < 0))
   call fig_points(uu, 'patch', pack(p, exposed(2)%v .and. m0(2)%iota == 0))
   call fig_label(uu, 'accent', pt(-2.2_R8P + SH * (run - 1), 2.45_R8P), trim(merge('step 0.6 h', 'step 3.6 h', &
                  run == 1)))
enddo
call fig_close(uu)
contains
   subroutine assemble(dx, masks, donors, previous)
   !< Holes, candidates and the explicit classification with the O-grid at pose (0, (dx, 0)).
   real(R8P),                         intent(in)           :: dx          !< Shift of the O-grid.
   type(overset_mask_t), allocatable, intent(out)          :: masks(:)    !< Classification.
   type(donor_set_t),    allocatable, intent(out)          :: donors(:)   !< Candidates.
   type(donor_set_t),                 intent(in), optional :: previous(:) !< Last candidates (walk starts).
   type(rigid_pose_2d_t)              :: poses(2)
   type(block_logical_t), allocatable :: holes(:)

   poses(1) = pose_2d(angle=0._R8P, translation=[dx, 0._R8P]) ; poses(2) = pose_2d(angle=0._R8P)
   call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses, previous=previous)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   endsubroutine assemble

   subroutine circle(cx, style)
   !< The cylinder's wall (the O-grid's j = 1 line) centred at (cx, 0).
   real(R8P),    intent(in) :: cx    !< Centre.
   character(*), intent(in) :: style !< Style.

   call fig_segments(uu, style, [(pt(u(i, 1) + cx, v(i, 1), 3 * LIFT), i = 1, 32)], &
                     [(pt(u(i + 1, 1) + cx, v(i + 1, 1), 3 * LIFT), i = 1, 32)])
   endsubroutine circle

   function pt(qx, qy, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: qx, qy !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = qx * ex_R8P + qy * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt
endprogram fig_ref_overset_exposed
