!< FOSSIL, wall-spacing and quality report of curvilinear blocks (issue #74).
module fossil_block_report
!< FOSSIL, wall-spacing and quality report of curvilinear blocks (issue #74).
!<
!< Design: docs/guide/advanced/body-fitted.md. For every wall point of a block (the faces covered by TOPO_BOUNDARY records
!< named `wall`, 2-D: the sides), along the wall-normal grid line (the index normal to the wall face, inward):
!< - y1 along the grid line: |x_2 - x_1|, the first off-wall vertex; with cell_centred = .true. the points are the wall
!<   cells instead: y1 = |c - w|, c the first cell's centre (mean of its 8 vertices, 2-D 4) and w its wall face's
!<   centre (mean of 4 vertices, 2-D 2);
!< - y1 to the wall: the distance of that off-wall point (vertex or cell centre) to the wall cutter (the caller's
!<   `wall_cutter` / `wall_cutter_2d` surface or curve);
!< - y1+ = y1_wall u_tau / nu: u_tau given (a scalar or one per wall point), or estimated by `utau_flat_plate` from
!<   u_inf and the distance x - x_le along +x (flagged: such estimates can be off by about 300 % at a leading edge);
!< - the growth ratio: the largest h_{k+1} / h_k along the grid line, h_k = |x_{k+1} - x_k|.
!< Totals: min, max, mean, 95th percentile of y1+ (lower median rank), the fraction of wall points above target_y1p, the
!< largest growth ratio, and the block's quality (WP2 `quality`): max and 95th percentile of skewness, min
!< orthogonality, max aspect, max warp (3-D), min of the scaled certified min J, the cells not certified valid (scaled
!< min J <= 0: invalid or undetermined).
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P
use fossil_surface_stl_object, only : surface_stl_object
use fossil_curve_object,       only : curve_object
use fossil_curvilinear_block,  only : curvilinear_block_t, curvilinear_block_2d_t, TOPO_BOUNDARY

implicit none
private
public :: block_report_t, block_report, block_report_2d, save_block_report, utau_flat_plate

type :: block_report_t
   !< Wall spacing and quality of a block.
   integer(I4P)              :: nw = 0_I4P            !< Wall points.
   integer(I4P), allocatable :: point(:,:)            !< Wall vertex (or the low vertex of the wall cell) [3, nw].
   integer(I4P), allocatable :: face(:)               !< Wall face of each point [nw].
   real(R8P),    allocatable :: y1_line(:)            !< y1 along the grid line [nw].
   real(R8P),    allocatable :: y1_wall(:)            !< y1 to the wall [nw].
   real(R8P),    allocatable :: y1p(:)                !< y1+ [nw] (0 without u_tau).
   real(R8P),    allocatable :: growth(:)             !< Largest wall-normal growth ratio [nw].
   logical                   :: cell_centred = .false. !< Points are wall cells.
   logical                   :: utau_given = .false.  !< u_tau given by the caller.
   logical                   :: utau_estimated = .false. !< u_tau from the flat-plate correlation.
   real(R8P)                 :: target_y1p = 1._R8P   !< Target y1+.
   real(R8P)                 :: y1p_min = 0._R8P, y1p_max = 0._R8P, y1p_mean = 0._R8P, y1p_p95 = 0._R8P !< y1+ stats.
   real(R8P)                 :: above = 0._R8P        !< Fraction of wall points with y1+ > target.
   real(R8P)                 :: growth_max = 0._R8P   !< Largest growth ratio.
   real(R8P)                 :: skew_max = 0._R8P, skew_p95 = 0._R8P, orth_min = 0._R8P, aspect_max = 0._R8P !< Quality.
   real(R8P)                 :: warp_max = 0._R8P, minj_min = 0._R8P !< Quality.
   integer(I4P)              :: not_valid = 0_I4P     !< Cells not certified valid.
   integer(I4P)              :: cells = 0_I4P         !< Cells.
endtype block_report_t

contains
   pure function utau_flat_plate(u_inf, nu, x) result(ut)
   !< Friction velocity of a turbulent flat plate at distance x from its leading edge (estimate, gk Schlichting):
   !< C_f = 0.0576 Re_x^(-1/5), u_tau = u_inf sqrt(C_f / 2); x <= 0 gives 0.
   real(R8P), intent(in) :: u_inf, nu, x !< Free stream, kinematic viscosity, distance.
   real(R8P)             :: ut           !< u_tau.

   ut = 0._R8P
   if (x > 0._R8P .and. nu > 0._R8P .and. u_inf > 0._R8P) ut = u_inf * sqrt(0.5_R8P * 0.0576_R8P * (u_inf * x / nu)**(-0.2_R8P))
   endfunction utau_flat_plate

   subroutine block_report(block, cutter, nu, report, u_tau, u_tau_wall, u_inf, x_le, target_y1p, cell_centred, &
                           wall_name, status)
   !< Wall spacing and quality report of a 3-D block (see the module notes). status: 0, 1 (no wall face, or bad input).
   type(curvilinear_block_t), intent(in)            :: block          !< Block.
   type(surface_stl_object),  intent(in)            :: cutter         !< Wall cutter.
   real(R8P),                 intent(in)            :: nu             !< Kinematic viscosity.
   type(block_report_t),      intent(out)           :: report         !< Report.
   real(R8P),                 intent(in),  optional :: u_tau          !< Friction velocity (all wall points).
   real(R8P),                 intent(in),  optional :: u_tau_wall(:)  !< Friction velocity per wall point.
   real(R8P),                 intent(in),  optional :: u_inf, x_le    !< Free stream and leading-edge x (estimate).
   real(R8P),                 intent(in),  optional :: target_y1p     !< Target y1+ (default 1).
   logical,                   intent(in),  optional :: cell_centred   !< Wall cells and first cell centres.
   character(*),              intent(in),  optional :: wall_name      !< Record name (default 'wall').
   integer(I4P),              intent(out), optional :: status         !< 0, 1.
   real(R8P),    allocatable                        :: pw(:,:), po(:,:), gr(:), ut(:), d(:), sk(:,:,:), orth(:,:,:)
   real(R8P),    allocatable                        :: asp(:,:,:), wp(:,:,:), mj(:,:,:)
   integer(I4P), allocatable                        :: pt(:,:), fc(:)
   type(vector_R8P),          allocatable           :: q(:)
   character(:),              allocatable           :: wn
   integer(I4P)                                     :: r, a, s, i, j, k, n, i0(3), i1(3), m, nl, l, nc(3)
   real(R8P)                                        :: w(3), c(3), h0, h1, g

   wn = 'wall' ; if (present(wall_name)) wn = wall_name
   report%cell_centred = .false. ; if (present(cell_centred)) report%cell_centred = cell_centred
   if (present(target_y1p)) report%target_y1p = target_y1p
   if (present(status)) status = 1
   if (.not. allocated(block%x) .or. cutter%get_facets_number() == 0) return
   n = 0
   allocate(pw(3, 0), po(3, 0), gr(0), pt(3, 0), fc(0))
   do r = 1, block%nrec
      if (block%rec(r)%kind /= TOPO_BOUNDARY .or. trim(block%rec(r)%name) /= wn) cycle
      a = (block%rec(r)%face + 1) / 2 ; s = merge(1_I4P, -1_I4P, mod(block%rec(r)%face, 2) == 1)
      nl = block%nk
      if (a == 1) nl = block%ni
      if (a == 2) nl = block%nj
      i0 = block%rec(r)%lo ; i1 = block%rec(r)%hi
      if (report%cell_centred) i1 = i1 - 1 ; i1(a) = i0(a)
      do k = i0(3), i1(3)
         do j = i0(2), i1(2)
            do i = i0(1), i1(1)
               call wall_point(block, [i, j, k], a, s, report%cell_centred, w, c)
               g = 0._R8P ; h0 = 0._R8P
               do l = 1, nl - 1 ! growth along the vertex line through [i, j, k]
                  h1 = line_step(block, [i, j, k], a, s, l)
                  if (l > 1 .and. h0 > 0._R8P) g = max(g, h1 / h0)
                  h0 = h1
               enddo
               call append_point(pw, po, gr, pt, fc, n, w, c, g, [i, j, k], block%rec(r)%face)
            enddo
         enddo
      enddo
   enddo
   if (n == 0) return
   allocate(q(n), d(n))
   do m = 1, n
      q(m) = vector_R8P(po(1, m), po(2, m), po(3, m))
   enddo
   call cutter%query_many(points=q, distances=d, is_signed=.false.)
   report%nw = n ; report%point = pt(:, 1:n) ; report%face = fc(1:n) ; report%growth = gr(1:n)
   report%y1_line = [(norm2(po(:, m) - pw(:, m)), m = 1, n)] ; report%y1_wall = abs(d)
   call friction(report, n, pw, ut, u_tau, u_tau_wall, u_inf, x_le, nu)
   report%y1p = report%y1_wall * ut / nu
   call y1p_stats(report)
   call block%quality(skewness=sk, orthogonality=orth, aspect=asp, warp=wp, min_jacobian=mj)
   nc = shape(sk)
   report%cells = product(nc)
   report%skew_max = maxval(sk) ; report%skew_p95 = percentile(pack(sk, .true.), 0.95_R8P)
   report%orth_min = minval(orth) ; report%aspect_max = maxval(asp) ; report%warp_max = maxval(wp)
   report%minj_min = minval(mj) ; report%not_valid = count(.not. (mj > 0._R8P))
   if (present(status)) status = 0
   endsubroutine block_report

   subroutine block_report_2d(block, cutter, nu, report, u_tau, u_tau_wall, u_inf, x_le, target_y1p, cell_centred, &
                              wall_name, status)
   !< 2-D twin of block_report (the curve cutter; x along u).
   type(curvilinear_block_2d_t), intent(in)            :: block         !< Block.
   type(curve_object),           intent(in)            :: cutter        !< Wall cutter.
   real(R8P),                    intent(in)            :: nu            !< Kinematic viscosity.
   type(block_report_t),         intent(out)           :: report        !< Report.
   real(R8P),                    intent(in),  optional :: u_tau         !< Friction velocity (all wall points).
   real(R8P),                    intent(in),  optional :: u_tau_wall(:) !< Friction velocity per wall point.
   real(R8P),                    intent(in),  optional :: u_inf, x_le   !< Free stream and leading-edge u (estimate).
   real(R8P),                    intent(in),  optional :: target_y1p    !< Target y1+ (default 1).
   logical,                      intent(in),  optional :: cell_centred  !< Wall cells and first cell centres.
   character(*),                 intent(in),  optional :: wall_name     !< Record name (default 'wall').
   integer(I4P),                 intent(out), optional :: status        !< 0, 1.
   real(R8P),    allocatable                           :: pw(:,:), po(:,:), gr(:), ut(:), d(:), sk(:,:), orth(:,:)
   real(R8P),    allocatable                           :: asp(:,:), mj(:,:), q(:,:)
   integer(I4P), allocatable                           :: pt(:,:), fc(:)
   character(:), allocatable                           :: wn
   integer(I4P)                                        :: r, a, s, i, j, n, i0(3), i1(3), nl, l, m
   real(R8P)                                           :: w(3), c(3), h0, h1, g

   wn = 'wall' ; if (present(wall_name)) wn = wall_name
   report%cell_centred = .false. ; if (present(cell_centred)) report%cell_centred = cell_centred
   if (present(target_y1p)) report%target_y1p = target_y1p
   if (present(status)) status = 1
   if (.not. allocated(block%u)) return
   n = 0
   allocate(pw(3, 0), po(3, 0), gr(0), pt(3, 0), fc(0))
   do r = 1, block%nrec
      if (block%rec(r)%kind /= TOPO_BOUNDARY .or. trim(block%rec(r)%name) /= wn) cycle
      a = (block%rec(r)%face + 1) / 2 ; s = merge(1_I4P, -1_I4P, mod(block%rec(r)%face, 2) == 1)
      nl = merge(block%ni, block%nj, a == 1)
      i0 = block%rec(r)%lo ; i1 = block%rec(r)%hi
      if (report%cell_centred) i1(1:2) = i1(1:2) - 1 ; i1(a) = i0(a)
      do j = i0(2), i1(2)
         do i = i0(1), i1(1)
            call wall_point_2d(block, [i, j], a, s, report%cell_centred, w, c)
            g = 0._R8P ; h0 = 0._R8P
            do l = 1, nl - 1
               h1 = line_step_2d(block, [i, j], a, s, l)
               if (l > 1 .and. h0 > 0._R8P) g = max(g, h1 / h0)
               h0 = h1
            enddo
            call append_point(pw, po, gr, pt, fc, n, w, c, g, [i, j, 1_I4P], block%rec(r)%face)
         enddo
      enddo
   enddo
   if (n == 0) return
   allocate(q(2, n), d(n)) ; q = po(1:2, 1:n)
   call cutter%query_many(points=q, distances=d, is_signed=.false.)
   report%nw = n ; report%point = pt(:, 1:n) ; report%face = fc(1:n) ; report%growth = gr(1:n)
   report%y1_line = [(norm2(po(:, m) - pw(:, m)), m = 1, n)] ; report%y1_wall = abs(d)
   call friction(report, n, pw, ut, u_tau, u_tau_wall, u_inf, x_le, nu)
   report%y1p = report%y1_wall * ut / nu
   call y1p_stats(report)
   call block%quality(skewness=sk, orthogonality=orth, aspect=asp, min_jacobian=mj)
   report%cells = size(sk, kind=I4P)
   report%skew_max = maxval(sk) ; report%skew_p95 = percentile(pack(sk, .true.), 0.95_R8P)
   report%orth_min = minval(orth) ; report%aspect_max = maxval(asp) ; report%warp_max = 0._R8P
   report%minj_min = minval(mj) ; report%not_valid = count(.not. (mj > 0._R8P))
   if (present(status)) status = 0
   endsubroutine block_report_2d

   subroutine save_block_report(report, file_name, status)
   !< Write the report as plain text: the totals, then one line per wall point (indices, face, y1 line, y1 wall, y1+,
   !< growth).
   type(block_report_t), intent(in)            :: report    !< Report.
   character(*),         intent(in)            :: file_name !< File.
   integer(I4P),         intent(out), optional :: status    !< 0, or the iostat of the open.
   integer                                     :: u, ios, m

   open(newunit=u, file=file_name, status='replace', action='write', iostat=ios)
   if (present(status)) status = int(ios, I4P)
   if (ios /= 0) return
   write(u, '(A)') '# FOSSIL block report (issue #74)'
   write(u, '(A,I0,A,L1,A,L1,A,L1)') '# wall points ', report%nw, ', cell centred ', report%cell_centred, &
                                     ', u_tau given ', report%utau_given, ', u_tau estimated (flat plate) ', &
                                     report%utau_estimated
   write(u, '(A,4ES14.6,A,ES12.4,A,F8.4)') '# y1+ min max mean p95 ', report%y1p_min, report%y1p_max, report%y1p_mean, &
                                           report%y1p_p95, '  target ', report%target_y1p, '  fraction above ', report%above
   write(u, '(A,ES14.6)') '# growth ratio max ', report%growth_max
   write(u, '(A,I0,A,2ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,I0)') '# cells ', report%cells, &
      '  skewness max p95 ', report%skew_max, report%skew_p95, '  orthogonality min ', report%orth_min, &
      '  aspect max ', report%aspect_max, '  warp max ', report%warp_max, '  scaled min J min ', report%minj_min, &
      '  not certified valid ', report%not_valid
   write(u, '(A)') '# i j k face y1_line y1_wall y1+ growth'
   do m = 1, report%nw
      write(u, '(3I7,I3,4ES16.8)') report%point(:, m), report%face(m), report%y1_line(m), report%y1_wall(m), &
                                   report%y1p(m), report%growth(m)
   enddo
   close(u)
   endsubroutine save_block_report

   subroutine wall_point(block, p, a, s, cells, w, c)
   !< The wall point w and the first off-wall point c of the vertex (or wall cell with low vertex) p, direction a, step s.
   type(curvilinear_block_t), intent(in)  :: block  !< Block.
   integer(I4P),              intent(in)  :: p(3)   !< Vertex / low vertex of the wall cell.
   integer(I4P),              intent(in)  :: a, s   !< Normal direction, inward step.
   logical,                   intent(in)  :: cells  !< Cell centred.
   real(R8P),                 intent(out) :: w(3), c(3) !< Points.
   integer(I4P)                           :: q(3), e1, e2, m1, m2, k

   if (.not. cells) then
      w = vx(block, p) ; q = p ; q(a) = q(a) + s ; c = vx(block, q)
      return
   endif
   e1 = merge(2_I4P, 1_I4P, a == 1) ; e2 = merge(2_I4P, 3_I4P, a == 3)
   w = 0._R8P ; c = 0._R8P
   do m2 = 0, 1
      do m1 = 0, 1
         q = p ; q(e1) = q(e1) + m1 ; q(e2) = q(e2) + m2
         w = w + vx(block, q)
         do k = 0, 1
            c = c + vx(block, [q(1) + merge(s * k, 0_I4P, a == 1), q(2) + merge(s * k, 0_I4P, a == 2), &
                               q(3) + merge(s * k, 0_I4P, a == 3)])
         enddo
      enddo
   enddo
   w = w / 4 ; c = c / 8
   endsubroutine wall_point

   subroutine wall_point_2d(block, p, a, s, cells, w, c)
   !< 2-D twin of wall_point.
   type(curvilinear_block_2d_t), intent(in)  :: block  !< Block.
   integer(I4P),                 intent(in)  :: p(2)   !< Vertex / low vertex of the wall cell.
   integer(I4P),                 intent(in)  :: a, s   !< Normal direction, inward step.
   logical,                      intent(in)  :: cells  !< Cell centred.
   real(R8P),                    intent(out) :: w(3), c(3) !< Points.
   integer(I4P)                              :: q(2), e, m, k

   if (.not. cells) then
      w = [block%u(p(1), p(2)), block%v(p(1), p(2)), 0._R8P] ; q = p ; q(a) = q(a) + s
      c = [block%u(q(1), q(2)), block%v(q(1), q(2)), 0._R8P]
      return
   endif
   e = 3 - a ; w = 0._R8P ; c = 0._R8P
   do m = 0, 1
      q = p ; q(e) = q(e) + m
      w = w + [block%u(q(1), q(2)), block%v(q(1), q(2)), 0._R8P]
      do k = 0, 1
         c = c + [block%u(q(1) + merge(s * k, 0_I4P, a == 1), q(2) + merge(s * k, 0_I4P, a == 2)), &
                  block%v(q(1) + merge(s * k, 0_I4P, a == 1), q(2) + merge(s * k, 0_I4P, a == 2)), 0._R8P]
      enddo
   enddo
   w = w / 2 ; c = c / 4
   endsubroutine wall_point_2d

   pure function vx(block, p) result(x)
   !< Vertex p of a 3-D block.
   type(curvilinear_block_t), intent(in) :: block !< Block.
   integer(I4P),              intent(in) :: p(3)  !< Index.
   real(R8P)                             :: x(3)  !< Coordinates.

   x = [block%x(p(1), p(2), p(3)), block%y(p(1), p(2), p(3)), block%z(p(1), p(2), p(3))]
   endfunction vx

   pure function line_step(block, p, a, s, l) result(h)
   !< |x_{l+1} - x_l| along the line through p in direction a (step s), l counted from the wall.
   type(curvilinear_block_t), intent(in) :: block  !< Block.
   integer(I4P),              intent(in) :: p(3)   !< Wall vertex.
   integer(I4P),              intent(in) :: a, s, l !< Direction, step, layer.
   real(R8P)                             :: h      !< Step length.
   integer(I4P)                          :: q0(3), q1(3)

   q0 = p ; q0(a) = p(a) + s * (l - 1) ; q1 = p ; q1(a) = p(a) + s * l
   h = norm2(vx(block, q1) - vx(block, q0))
   endfunction line_step

   pure function line_step_2d(block, p, a, s, l) result(h)
   !< 2-D twin of line_step.
   type(curvilinear_block_2d_t), intent(in) :: block  !< Block.
   integer(I4P),                 intent(in) :: p(2)   !< Wall vertex.
   integer(I4P),                 intent(in) :: a, s, l !< Direction, step, layer.
   real(R8P)                                :: h      !< Step length.
   integer(I4P)                             :: q0(2), q1(2)

   q0 = p ; q0(a) = p(a) + s * (l - 1) ; q1 = p ; q1(a) = p(a) + s * l
   h = hypot(block%u(q1(1), q1(2)) - block%u(q0(1), q0(2)), block%v(q1(1), q1(2)) - block%v(q0(1), q0(2)))
   endfunction line_step_2d

   subroutine append_point(pw, po, gr, pt, fc, n, w, c, g, p, f)
   !< Append a wall point (growing the arrays by doubling, no self-referencing growth).
   real(R8P),    allocatable, intent(inout) :: pw(:,:), po(:,:), gr(:) !< Points, growth.
   integer(I4P), allocatable, intent(inout) :: pt(:,:), fc(:)          !< Indices, faces.
   integer(I4P),              intent(inout) :: n                       !< Count.
   real(R8P),                 intent(in)    :: w(3), c(3), g           !< Point data.
   integer(I4P),              intent(in)    :: p(3), f                 !< Index, face.
   real(R8P),    allocatable                :: t2(:,:), t1(:)
   integer(I4P), allocatable                :: i2(:,:), i1(:)
   integer(I4P)                             :: cap

   if (n == size(gr)) then
      cap = max(64_I4P, 2 * n)
      allocate(t2(3, cap)) ; t2(:, 1:n) = pw(:, 1:n) ; call move_alloc(t2, pw)
      allocate(t2(3, cap)) ; t2(:, 1:n) = po(:, 1:n) ; call move_alloc(t2, po)
      allocate(t1(cap)) ; t1(1:n) = gr(1:n) ; call move_alloc(t1, gr)
      allocate(i2(3, cap)) ; i2(:, 1:n) = pt(:, 1:n) ; call move_alloc(i2, pt)
      allocate(i1(cap)) ; i1(1:n) = fc(1:n) ; call move_alloc(i1, fc)
   endif
   n = n + 1 ; pw(:, n) = w ; po(:, n) = c ; gr(n) = g ; pt(:, n) = p ; fc(n) = f
   endsubroutine append_point

   subroutine friction(report, n, pw, ut, u_tau, u_tau_wall, u_inf, x_le, nu)
   !< u_tau per wall point: given, per point, or the flat-plate estimate (0 without any).
   type(block_report_t),   intent(inout)        :: report        !< Report (flags).
   integer(I4P),           intent(in)           :: n             !< Wall points.
   real(R8P),              intent(in)           :: pw(:,:)       !< Wall points [3, n].
   real(R8P), allocatable, intent(out)          :: ut(:)         !< u_tau [n].
   real(R8P),              intent(in), optional :: u_tau, u_tau_wall(:), u_inf, x_le !< Inputs.
   real(R8P),              intent(in)           :: nu            !< Kinematic viscosity.
   integer(I4P)                                 :: m
   real(R8P)                                    :: x0            !< Leading edge.

   allocate(ut(n)) ; ut = 0._R8P
   if (present(u_tau_wall)) then
      if (size(u_tau_wall) == n) then
         ut = u_tau_wall ; report%utau_given = .true. ; return
      endif
   endif
   if (present(u_tau)) then
      ut = u_tau ; report%utau_given = .true. ; return
   endif
   if (present(u_inf)) then
      x0 = 0._R8P
      if (present(x_le)) x0 = x_le
      do m = 1, n
         ut(m) = utau_flat_plate(u_inf, nu, pw(1, m) - x0)
      enddo
      report%utau_estimated = .true.
   endif
   endsubroutine friction

   subroutine y1p_stats(report)
   !< Totals of y1+ and of the growth ratio.
   type(block_report_t), intent(inout) :: report !< Report.

   report%y1p_min = minval(report%y1p) ; report%y1p_max = maxval(report%y1p)
   report%y1p_mean = sum(report%y1p) / report%nw ; report%y1p_p95 = percentile(report%y1p, 0.95_R8P)
   report%above = real(count(report%y1p > report%target_y1p), R8P) / report%nw
   report%growth_max = maxval(report%growth)
   endsubroutine y1p_stats

   function percentile(a, f) result(p)
   !< The value of rank ceiling(f n) of a (insertion-free: a sorted copy by merge passes).
   real(R8P), intent(in) :: a(:) !< Values.
   real(R8P), intent(in) :: f    !< Fraction.
   real(R8P)             :: p    !< Percentile.
   real(R8P), allocatable :: b(:), t(:)
   integer(I4P)          :: n, w, lo, mid, hi, i, j, k

   n = size(a, kind=I4P) ; p = 0._R8P
   if (n == 0) return
   b = a ; allocate(t(n)) ; w = 1
   do while (w < n)
      lo = 1
      do while (lo <= n)
         mid = min(lo + w - 1, n) ; hi = min(lo + 2 * w - 1, n) ; i = lo ; j = mid + 1 ; k = lo
         do while (i <= mid .and. j <= hi)
            if (b(j) < b(i)) then
               t(k) = b(j) ; j = j + 1
            else
               t(k) = b(i) ; i = i + 1
            endif
            k = k + 1
         enddo
         t(k:k + mid - i) = b(i:mid) ; k = k + mid - i + 1 ; t(k:k + hi - j) = b(j:hi)
         b(lo:hi) = t(lo:hi) ; lo = lo + 2 * w
      enddo
      w = 2 * w
   enddo
   p = b(max(1_I4P, min(n, ceiling(f * n, I4P))))
   endfunction percentile
endmodule fossil_block_report
