!< FOSSIL validation: the 2-D exact rotation by local integration against the per-face reference (issue #104).
program fossil_val_rotation_2d_local
!< FOSSIL validation: the 2-D exact rotation by local integration against the per-face reference (issue #104).
!<
!< Reference: `swept_cells_2d_reference` (the exact rotation path of issue #69 on every face, `faces_rotation_2d`:
!< certified breakpoints per face, Gauss doubled to round-off per piece). Against it:
!< - G: the Gauss path, nq = 4 (`exact=.false.`), the cost unit;
!< - a: the Gauss path with `refine=.true.` (the map's event cells integrated per event interval, the other faces Gauss);
!< - b: `swept_cells_2d` on the exact rotation path since issue #104 decision 2 (every candidate cell of the map
!<   integrated per event interval to round-off, the faces of no candidate cell in closed form, the rest by
!<   `faces_rotation_2d`: `sw%fallback_faces`).
!< Bodies: the 1024-gon ellipse and the 1022-vertex NACA 0012 of the curve bench, on the curve bench grid (n^2 cells over
!< [0.0031, 0.9969]^2 grown by 2 %). Motions about the block centre (0.5, 0.5) from a 3-degree pose: 0.5, 1 and 2 cells
!< per step at the body's extremity (the vertex farthest from the centre), 20 degrees, and a screw (1 cell of rotation
!< plus a translation (0.7, 0.4) h).
!< Per case: max |face difference| / h and max |outflow difference| / h^2 against the reference, the largest per-cell
!< residual / h^2 (closed cells), `exact`, the best of nrep times (one thread), for b its fallback faces and
!< evaluations, and with `map=` the evaluated cells and the largest difference between a cell's track sums and its
!< swept faces / outflows (round-off: a face between two measured cells takes the low cell's value).
!< Checks (each failure counted): every call STATUS_OK; the reference and b exact; b without fallback faces; b's faces
!< within 16 eps xm of the reference's and its outflows within 16 eps xm vmax (xm = max |grid| + |t0| + |dt|, vmax =
!< |dth| xm + |dt|); b with `map=` equal to b without bit for bit (faces, outflows, residual).
!< Usage: fossil_val_rotation_2d_local [n [body [motion [nrep]]]] (0 = all; default n = 256, all bodies and motions,
!< nrep = 1; n = 1024 is the step-0 table's second size, the reference dominating its time).
use fossil, only : curve_object, pose_2d, rigid_pose_2d_t, swept_cells_2d_t, swept_cells_2d, swept_cells_2d_reference, &
                   piece_map_2d_t, STATUS_OK
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)            !< Pi.
real(R8P), parameter      :: CEN(2) = [0.5_R8P, 0.5_R8P]   !< Block centre.
character(10), parameter  :: MNAME(5) = [character(10) :: '0.5 cell', '1 cell', '2 cells', '20 deg', 'screw'] !< Motions.
type(curve_object)        :: cb(2)                         !< Bodies.
character(10)             :: bname(2)                      !< Body names.
real(R8P),    allocatable :: uu(:), vv(:)                  !< Work.
character(32)             :: arg                           !< Argument.
integer(I4P)              :: n, bsel, msel, nrep, b, m, k, bad !< Selection, work.

n = 256_I4P ; bsel = 0_I4P ; msel = 0_I4P ; nrep = 1_I4P
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg) ; read(arg, *) n
endif
if (command_argument_count() >= 2) then
   call get_command_argument(2, arg) ; read(arg, *) bsel
endif
if (command_argument_count() >= 3) then
   call get_command_argument(3, arg) ; read(arg, *) msel
endif
if (command_argument_count() >= 4) then
   call get_command_argument(4, arg) ; read(arg, *) nrep
endif
allocate(uu(1024), vv(1024))
do k = 1, 1024
   uu(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 1024) ; vv(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 1024)
enddo
call cb(1)%adopt_loops(u=uu, v=vv, loop_sizes=[1024_I4P])
deallocate(uu, vv) ; allocate(uu(1022), vv(1022))
call naca0012(511_I4P, uu, vv)
uu = 0.15_R8P + 0.7_R8P * uu ; vv = 0.5_R8P + 0.7_R8P * vv
call cb(2)%adopt_loops(u=uu, v=vv, loop_sizes=[1022_I4P])
bname = [character(10) :: 'ellipse', 'NACA 0012']
print '(A,I0,A,I0)', 'exact 2-D rotation, local vs per-face reference (issue #104), n = ', n, ', nrep = ', nrep
print '(A)', 'units: faces / h, outflows / h^2, residuals / h^2 (cell area); times in s, best of nrep, one thread'
bad = 0_I4P
do b = 1, 2
   if (bsel /= 0 .and. b /= bsel) cycle
   do m = 1, 5
      if (msel /= 0 .and. m /= msel) cycle
      call run_case(b, m, bad)
   enddo
enddo
print '(A,I0)', 'failed checks: ', bad
print '(A,L1)', 'Are all tests passed? ', bad == 0_I4P

contains
   subroutine run_case(b, m, bad)
   !< One body and motion: the reference, G, a, b, then b with `map=`.
   integer(I4P), intent(in)    :: b, m  !< Body, motion.
   integer(I4P), intent(inout) :: bad   !< Failed checks.
   type(swept_cells_2d_t)      :: sr, sg, sa, sb, sm !< Reference, Gauss, refine, local, local with the map.
   type(piece_map_2d_t)        :: mp    !< Map of b.
   type(rigid_pose_2d_t)       :: p0, p1 !< End poses.
   real(R8P),    allocatable   :: x(:), y(:), our(:,:), ovr(:,:), oua(:,:), ova(:,:), oub(:,:), ovb(:,:), oum(:,:), &
                                  ovm(:,:) !< Grid, outflows.
   logical,      allocatable   :: cl(:,:) !< Closed cells.
   real(R8P)                   :: h, rmax, ang, th0, tr(2), tm(4), t1, dfa(2), dfb(2), rr(4), dmap(2), xm, vmax !< Work.
   integer(I4P)                :: r, st(5), nb       !< Work.

   call grid_edges([0.0031_R8P, 0.0031_R8P], [0.9969_R8P, 0.9969_R8P], n, x, y)
   h = x(2) - x(1)
   rmax = extremity(cb(b))
   th0 = 3._R8P * PI / 180._R8P ; tr = 0._R8P
   select case (m)
   case (1) ; ang = 0.5_R8P * h / rmax
   case (2) ; ang = h / rmax
   case (3) ; ang = 2._R8P * h / rmax
   case (4) ; ang = 20._R8P * PI / 180._R8P
   case default ; ang = h / rmax ; tr = [0.7_R8P, 0.4_R8P] * h
   endselect
   p0 = pose_2d(th0) ; p0%t = CEN - matmul(p0%r, CEN)
   p1 = pose_2d(th0 + ang) ; p1%t = CEN - matmul(p1%r, CEN) + tr
   xm = max(maxval(abs(x)), maxval(abs(y))) + norm2(p0%t) + norm2(p1%t - p0%t)
   vmax = abs(ang) * xm + norm2(tr)
   allocate(our(n + 1, n), ovr(n, n + 1), oua(n + 1, n), ova(n, n + 1), oub(n + 1, n), ovb(n, n + 1), oum(n + 1, n), &
            ovm(n, n + 1))
   tm = huge(1._R8P) ; nb = 0_I4P
   do r = 1, nrep
      t1 = clock() ; call swept_cells_2d_reference(cb(b), p0, p1, x, y, sr, status=st(1), outflow_u=our, outflow_v=ovr)
      tm(1) = min(tm(1), clock() - t1)
      t1 = clock() ; call swept_cells_2d(cb(b), p0, p1, x, y, sg, nq=4_I4P, exact=.false., status=st(2))
      tm(2) = min(tm(2), clock() - t1)
      t1 = clock()
      call swept_cells_2d(cb(b), p0, p1, x, y, sa, nq=4_I4P, exact=.false., refine=.true., status=st(3), outflow_u=oua, &
                          outflow_v=ova)
      tm(3) = min(tm(3), clock() - t1)
      t1 = clock() ; call swept_cells_2d(cb(b), p0, p1, x, y, sb, status=st(4), outflow_u=oub, outflow_v=ovb)
      tm(4) = min(tm(4), clock() - t1)
   enddo
   dfa(1) = max(maxval(abs(sa%face_u - sr%face_u)), maxval(abs(sa%face_v - sr%face_v))) / h
   dfa(2) = max(maxval(abs(oua - our)), maxval(abs(ova - ovr))) / h**2
   dfb(1) = max(maxval(abs(sb%face_u - sr%face_u)), maxval(abs(sb%face_v - sr%face_v))) / h
   dfb(2) = max(maxval(abs(oub - our)), maxval(abs(ovb - ovr))) / h**2
   ! b with the map (untimed)
   call swept_cells_2d(cb(b), p0, p1, x, y, sm, status=st(5), map=mp, outflow_u=oum, outflow_v=ovm)
   if (any(sm%face_u /= sb%face_u) .or. any(sm%face_v /= sb%face_v) .or. any(sm%residual /= sb%residual) .or. &
       any(oum /= oub) .or. any(ovm /= ovb)) nb = nb + 1_I4P
   dmap = map_vs_faces(mp, sm, oum, ovm) / [h, h**2]
   cl = sr%alpha0 >= 0._R8P .and. sr%alpha1 >= 0._R8P
   rr(1) = maxval(abs(sr%residual), mask=cl) / h**2 ; rr(2) = maxval(abs(sg%residual), mask=cl) / h**2
   rr(3) = maxval(abs(sa%residual), mask=cl) / h**2 ; rr(4) = maxval(abs(sb%residual), mask=cl) / h**2
   nb = nb + count(st /= STATUS_OK, kind=I4P)
   if (.not. sr%exact .or. .not. sb%exact .or. sb%fallback_faces /= 0_I4P) nb = nb + 1_I4P
   if (dfb(1) * h > 16._R8P * epsilon(1._R8P) * xm .or. dfb(2) * h**2 > 16._R8P * epsilon(1._R8P) * xm * vmax) &
      nb = nb + 1_I4P
   bad = bad + nb
   print '(A,A,A,A,ES9.2)', trim(bname(b)), ' ', trim(MNAME(m)), ': angle [deg] ', ang * 180._R8P / PI
   print '(A,4ES9.2,A,L1,A,L1)', '  residual ref/G/a/b ', rr, ', exact ref ', sr%exact, ' b ', sb%exact
   print '(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)', '  |dface| a ', dfa(1), ' b ', dfb(1), ', |doutflow| a ', dfa(2), ' b ', &
         dfb(2)
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '  b: fallback faces ', sb%fallback_faces, ', evaluations ', sb%evaluations, &
         ' (reference, per face: ', sr%evaluations, '), map evaluated cells ', mp%evaluated_cells, ', tracks ', mp%nt
   print '(A,4F9.4,A,F7.2,A,F7.2,A,F7.2)', '  times ref/G/a/b ', tm, ', ref/b ', tm(1) / tm(4), ' a/G ', tm(3) / tm(2), &
         ' b/G ', tm(4) / tm(2)
   print '(A,ES9.2,A,ES9.2,A,I0)', '  b map: max |track sums - faces| ', dmap(1), ', outflows ', dmap(2), &
         ', failed checks ', nb
   endsubroutine run_case

   function clock() result(t)
   !< Wall clock, seconds.
   real(R8P)    :: t        !< Seconds.
   integer(I8P) :: c, rate  !< Clock.

   call system_clock(c, rate)
   t = real(c, R8P) / real(rate, R8P)
   endfunction clock

   function map_vs_faces(mp, sw, fou, fov) result(d)
   !< Largest |sum over a cell's tracks - the cell's swept face| (faces, outflows) over the cells with tracks.
   type(piece_map_2d_t),   intent(in) :: mp         !< Map.
   type(swept_cells_2d_t), intent(in) :: sw         !< Swept geometry.
   real(R8P),              intent(in) :: fou(:,:), fov(:,:) !< Outflows.
   real(R8P)                          :: d(2), fl(4), ol(4) !< Result, sums.
   integer(I4P)                       :: i, j, t    !< Work.

   d = 0._R8P
   do j = 1, mp%ny
      do i = 1, mp%nx
         if (mp%tracks(i, j) == 0) cycle
         fl = 0._R8P ; ol = 0._R8P
         do t = mp%first_track(i, j), mp%first_track(i, j) + mp%tracks(i, j) - 1
            if (mp%status(t) < 0) cycle
            fl = fl + mp%faces(:, t) ; ol = ol + mp%outflow(:, t)
         enddo
         if (mp%status(mp%first_track(i, j)) < 0) cycle
         d(1) = max(d(1), maxval(abs(fl - [sw%face_u(i, j), sw%face_u(i + 1, j), sw%face_v(i, j), sw%face_v(i, j + 1)])))
         d(2) = max(d(2), maxval(abs(ol - [-fou(i, j), fou(i + 1, j), -fov(i, j), fov(i, j + 1)])))
      enddo
   enddo
   endfunction map_vs_faces

   function extremity(c) result(r)
   !< Largest distance of a vertex from the block centre.
   type(curve_object), intent(in) :: c  !< Curve.
   real(R8P)                      :: r, uv(2) !< Result, vertex.
   integer(I4P)                   :: k  !< Counter.

   r = 0._R8P
   do k = 1, c%get_vertices_number()
      uv = c%get_vertex(k) ; r = max(r, norm2(uv - CEN))
   enddo
   endfunction extremity

   subroutine grid_edges(lo, hi, n, x, y)
   !< n x n cells over the box [lo, hi] grown by 2 % (the grid of the curve bench).
   real(R8P),              intent(in)  :: lo(2), hi(2) !< Box.
   integer(I4P),           intent(in)  :: n            !< Cells per axis.
   real(R8P), allocatable, intent(out) :: x(:), y(:)   !< Cell edges.
   real(R8P)                           :: a(2), b(2)   !< Grown box.
   integer(I4P)                        :: k            !< Counter.

   a = lo - 0.02_R8P * (hi - lo) ; b = hi + 0.02_R8P * (hi - lo)
   allocate(x(n + 1), y(n + 1))
   do k = 1, n + 1
      x(k) = a(1) + (b(1) - a(1)) * (k - 1) / real(n, R8P) ; y(k) = a(2) + (b(2) - a(2)) * (k - 1) / real(n, R8P)
   enddo
   endsubroutine grid_edges

   subroutine naca0012(m, x, y)
   !< NACA 0012 with a closed trailing edge, 2m vertices, cosine spacing, counter-clockwise.
   integer(I4P), intent(in)  :: m          !< Points per side.
   real(R8P),    intent(out) :: x(:), y(:) !< Vertices [2m].
   real(R8P)                 :: xi         !< Abscissa.
   integer(I4P)              :: kk, j      !< Counters.

   j = 0
   do kk = m, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = half_thickness(xi)
   enddo
   do kk = 1, m - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = -half_thickness(xi)
   enddo
   endsubroutine naca0012

   pure function half_thickness(x) result(y)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction half_thickness
endprogram fossil_val_rotation_2d_local
