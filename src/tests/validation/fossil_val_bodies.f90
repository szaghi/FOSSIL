!< FOSSIL validation: the closure table, the load samples and the cut-cell loads of several bodies in one block (#109 WP3
!< step 05). A geometry study, no flow solve.
program fossil_val_bodies
!< FOSSIL validation: the closure table, the load samples and the cut-cell loads of several bodies in one block (#109 WP3
!< step 05). A geometry study, no flow solve.
!<
!< 1. Locality. Two bodies farther apart than twice the stencil reach (gap > 2 reach h), posed, with their own rigid
!<    velocities, curvature, probes and ILW crossings: every ghost of `fd_closure_3d_bodies` / `_2d_bodies` equals the
!<    ghost at the same indexes of the single-body table of its body (`fd_closure_3d` / `_2d` at the body's pose), every
!<    field by bits, and the ghosts are the union of the two; every load sample of `load_samples_3d_bodies` / `_2d_bodies`
!<    equals the slot of the single-body set of its body, by bits, in body order. 3-D: an icosphere (radius 1, 1 280
!<    facets) and an ellipsoid (semi-axes 1, 0.7, 0.5, 5 120 facets) at three pose pairs, h = 0.1; 2-D: a NACA 0012
!<    (800 segments) and an ellipse 0.3 x 0.12 (512 segments) at three pose pairs, h = 0.01; reach 2 per axis.
!< 2. Narrow gap. Two spheres of radius R (icospheres of 5 120 facets, R = 1, h = 0.1) / circles (1 024 segments,
!<    R = 0.5, h = 0.02) at gaps delta = 0.25, 0.5, 1, 2, 4 h, posed (rotated about their centres), reach 2 (default
!<    shared distance 2 h), probes h, 2 h, 3 h. The analytic clearance of a ghost is the distance from its (library) foot
!<    f to the true other sphere, |f - c_B| - R. The other body is a polyhedron inscribed in its sphere and containing the
!<    ball of radius R - s_B (s_B = R - least facet-plane distance, the sag), so the library clearance lies in
!<    [analytic, analytic + s_B]; the shared ghosts therefore number between the analytic counts N(d < sd - s_B) and
!<    N(d < sd), and a probe segment (foot to probe) is blocked if it meets the ball R - s_B of the other body and not
!<    if it misses the ball R (the own body is convex). The least clearance is at least the gap (the feet lie in their
!<    own ball).
!< 3. Loads through the union. In p = p0 + g . x the force on each closed body is -g V_k and its moment
!<    V_k (x_c,k - c) x (-g) (divergence theorem), whatever the other body. FV path: `cut_cells_3d_bodies(walls=,
!<    piece_bodies=)` / `cut_cells_2d_bodies`, the wall terms of body k selected by `piece_bodies%wall_body` (the
!<    pressure and its gradient zero on the other walls), `cut_loads_3d` with the gyration of the merged surface /
!<    `cut_loads_2d` with the gradient; FD path: `load_samples_3d_bodies` / `_2d_bodies`, `integrate_loads` with
!<    `per_body_force` / `per_body_moment`, rule orders 1, 2, 3. Errors in units of eps times the sum of the absolute
!<    terms of the body (as `fossil_val_surface_loads`). Cases: two posed spheres a quarter cell apart, a sphere and an
!<    ellipsoid half a cell apart, two unit cubes in face contact (aligned and shifted) on grids whose planes miss
!<    (h = 0.13) and contain (h = 0.125) the contact; 2-D twins (circles, a NACA 0012 and an ellipse, unit squares in
!<    edge contact). In a contact both bodies keep their wall terms on the contact face (opposite vector areas): each
!<    body's sum is the integral over its own closed surface.
!<
!< Pass criteria (`Are all tests passed? T` iff all hold; exit status 1 otherwise):
!<   1. locality: 0 ghost and 0 sample mismatches at every pose pair, the ghost count the sum of the single-body ones, and
!<      the least clearance of the table above 2 reach h;
!<   2. narrow gap: N(d < sd - s_B) <= shared <= N(d < sd) with every ghost consistent (shared and analytic < sd; not
!<      shared and analytic >= sd - s_B), the clearance error in [-1e-13, s_B + 1e-13], the least clearance >= delta -
!<      1e-13, and no probe against its bracket (blocked yet missing the ball R; visible yet meeting the ball R - s_B);
!<   3. loads: per body the force within 1 (eps S_F) on both paths and every rule order, the moment within 1 (eps S_M)
!<      where it is exact (FD from order 2, FV with the gradient (and the gyration)); every cut cell trusted;
!<   4. every call STATUS_OK.
use fossil,                  only : surface_stl_object, curve_object, facet_object, ib_closure_t, fd_closure_3d, &
                                    fd_closure_2d, fd_closure_3d_bodies, fd_closure_2d_bodies, load_samples_t, &
                                    load_samples_3d, load_samples_2d, load_samples_3d_bodies, load_samples_2d_bodies, &
                                    integrate_loads, cut_cells_3d_t, cut_walls_3d_t, cut_cells_2d_t, cut_walls_2d_t, &
                                    piece_bodies_t, cut_cells_3d_bodies, cut_cells_2d_bodies, cut_walls_gyration, &
                                    cut_loads_3d, cut_loads_2d, rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, &
                                    pose_2d, to_world, STATUS_OK
use fossil_thickness_bodies, only : make_sphere
use iso_fortran_env,         only : int64
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P), EPS = epsilon(1._R8P), TOL = 1.e-13_R8P
real(R8P),    parameter :: G3(3) = [0.3_R8P, -0.7_R8P, 1.1_R8P], C3(3) = [0.1_R8P, -0.2_R8P, 0.3_R8P], P0 = 2._R8P
real(R8P),    parameter :: G2(2) = [0.4_R8P, -0.9_R8P], C2(2) = [0.15_R8P, -0.25_R8P]
real(R8P),    parameter :: GAPS(5) = [0.25_R8P, 0.5_R8P, 1._R8P, 2._R8P, 4._R8P]
integer(I4P), parameter :: REACH = 2_I4P
integer(I4P)            :: nstat = 0_I4P
logical                 :: ok = .true.

call locality_3d
call locality_2d
call gap_3d
call gap_2d
call loads_3d
call loads_2d
call check('every call STATUS_OK', nstat == 0_I4P)
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   ! 1. locality

   subroutine locality_3d()
   !< Part 1, 3-D: an icosphere and an ellipsoid far apart at three pose pairs.
   type(surface_stl_object)        :: sph, ell, two(2)
   type(facet_object), allocatable :: fc(:)
   type(ib_closure_t)              :: cb, cs(2)
   type(load_samples_t)            :: lb, ls(2)
   type(rigid_pose_t)              :: pp(2, 3)
   real(R8P), allocatable          :: x(:), y(:)
   real(R8P)                       :: vl(3, 2), wl(3, 2), cl(3, 2), h
   integer(I4P)                    :: c, b, g, st, nbad, sbad, off
   integer(I4P), allocatable       :: idx(:,:,:,:)

   print '(A)', '1. Locality: bodies farther apart than 2 reach h; bodies table vs the single-body tables, bit for bit'
   print '(A)', '   3-D icosphere (1 280 facets) + ellipsoid (5 120 facets), h = 0.1, reach 2, probes h, 2 h, curvature, ILW'
   print '(A)', '   pose   ghosts  (body 1, body 2)  ghost mismatches  least clearance / (reach h)   samples  sample mismatches'
   h = 0.1_R8P
   call make_sphere(3_I4P, fc) ; call sph%adopt_facets(facets=fc) ; call sph%build_curvature(status=st) ; call stat(st)
   call make_sphere(4_I4P, fc) ; call ell%adopt_facets(facets=fc)
   call ell%resize(factor=vector_R8P(1._R8P, 0.7_R8P, 0.5_R8P)) ; call ell%analyze ; call ell%build_curvature(status=st)
   call stat(st)
   two(1) = sph ; two(2) = ell
   x = [(-2.65_R8P + h * g, g = 0, 53)] ; y = [(-1.45_R8P + h * g, g = 0, 29)]
   vl = reshape([0.1_R8P, -0.2_R8P, 0.3_R8P, -0.3_R8P, 0.05_R8P, 0.2_R8P], [3, 2])
   wl = reshape([0.5_R8P, -1._R8P, 2._R8P, -0.7_R8P, 0.4_R8P, 0.1_R8P], [3, 2])
   cl = reshape([0.01_R8P, 0.02_R8P, -0.03_R8P, 0.2_R8P, -0.1_R8P, 0.05_R8P], [3, 2])
   pp(1, 1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.2_R8P, [-1.3_R8P, 0.03_R8P, 0.01_R8P])
   pp(2, 1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], -0.3_R8P, [1.3_R8P, -0.02_R8P, 0.04_R8P])
   pp(1, 2) = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.7_R8P, [-1.27_R8P, 0.11_R8P, -0.07_R8P])
   pp(2, 2) = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], 1.1_R8P, [1.31_R8P, -0.13_R8P, 0.09_R8P])
   pp(1, 3) = pose_from_axis_angle([1._R8P, -1._R8P, 0._R8P], -0.4_R8P, [-1.25_R8P, -0.2_R8P, 0.15_R8P])
   pp(2, 3) = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 1._R8P, [1.25_R8P, 0.2_R8P, -0.1_R8P])
   do c = 1, 3
      call fd_closure_3d_bodies(two, x, y, y, [REACH, REACH, REACH], cb, status=st, poses=pp(:, c), &
                                lengths=[h, 2._R8P * h], linear_velocities=vl, angular_velocities=wl, centres=cl, ilw=.true.)
      call stat(st)
      do b = 1, 2
         call fd_closure_3d(two(b), x, y, y, [REACH, REACH, REACH], cs(b), status=st, lengths=[h, 2._R8P * h], &
                            pose=pp(b, c), linear_velocity=vl(:, b), angular_velocity=wl(:, b), centre=cl(:, b), ilw=.true.)
         call stat(st)
      enddo
      call ghost_index(cs, size(x, kind=I4P), size(y, kind=I4P), size(y, kind=I4P), idx)
      nbad = 0_I4P
      do g = 1, cb%n
         b = cb%body(g)
         if (.not. same_ghost(cb, g, cs(b), idx(cb%ijk(1, g), cb%ijk(2, g), cb%ijk(3, g), b), b)) nbad = nbad + 1_I4P
      enddo
      if (cb%n /= cs(1)%n + cs(2)%n .or. cb%overlap_points /= 0_I4P) nbad = nbad + 1_I4P
      call load_samples_3d_bodies(two, x, y, y, lb, status=st, poses=pp(:, c), order=3_I4P, lengths=[h, 2._R8P * h])
      call stat(st)
      sbad = 0_I4P ; off = 0_I4P
      do b = 1, 2
         call load_samples_3d(two(b), x, y, y, ls(b), status=st, order=3_I4P, lengths=[h, 2._R8P * h], pose=pp(b, c))
         call stat(st)
         if (.not. same_slice(lb, off, ls(b), b)) sbad = sbad + 1_I4P
         off = off + ls(b)%n
      enddo
      if (lb%n /= off) sbad = sbad + 1_I4P
      print '(3X,I4,I9,2X,"(",I6,",",I6,")",I14,F25.2,I14,I12)', c, cb%n, cs(1)%n, cs(2)%n, nbad, &
            minval(cb%clearance) / (REACH * h), lb%n, sbad
      call check('locality 3-D, pose pair '//trim(itoa(c))//': every ghost and sample = its single-body run, bit for bit', &
                 nbad == 0_I4P .and. sbad == 0_I4P .and. cb%n > 0_I4P .and. minval(cb%clearance) > 2 * REACH * h)
   enddo
   endsubroutine locality_3d

   subroutine locality_2d()
   !< Part 1, 2-D: a NACA 0012 and an ellipse far apart at three pose pairs.
   type(curve_object)     :: naca, ell, two(2)
   type(ib_closure_t)     :: cb, cs(2)
   type(load_samples_t)   :: lb, ls(2)
   type(rigid_pose_2d_t)  :: pp(2, 3)
   real(R8P), allocatable :: x(:), y(:), u(:), v(:)
   real(R8P)              :: vl(2, 2), wl(2), cl(2, 2), h
   integer(I4P)           :: c, b, g, st, nbad, sbad, off
   integer(I4P), allocatable :: idx(:,:,:,:)

   print '(A)', '   2-D NACA 0012 (800 segments) + ellipse 0.3 x 0.12 (512 segments), h = 0.01, reach 2, probes h, 2 h'
   print '(A)', '   pose   ghosts  (body 1, body 2)  ghost mismatches  least clearance / (reach h)   samples  sample mismatches'
   h = 0.01_R8P
   call naca0012(400_I4P, u, v)
   call naca%adopt_loops(u=u, v=v, loop_sizes=[size(u, kind=I4P)]) ; call naca%build_curvature(status=st) ; call stat(st)
   u = [(0.3_R8P * cos(2 * PI * g / 512), g = 0, 511)] ; v = [(0.12_R8P * sin(2 * PI * g / 512), g = 0, 511)]
   call ell%adopt_loops(u=u, v=v, loop_sizes=[512_I4P]) ; call ell%build_curvature(status=st) ; call stat(st)
   two(1) = naca ; two(2) = ell
   x = [(-0.253_R8P + h * g, g = 0, 200)] ; y = [(-0.447_R8P + h * g, g = 0, 105)]
   vl = reshape([0.1_R8P, -0.2_R8P, -0.3_R8P, 0.05_R8P], [2, 2]) ; wl = [2._R8P, -1.5_R8P]
   cl = reshape([0.25_R8P, 0._R8P, 1.4_R8P, 0.1_R8P], [2, 2])
   pp(1, 1) = pose_2d(0._R8P, [0._R8P, 0._R8P])          ; pp(2, 1) = pose_2d(0.3_R8P, [1.45_R8P, 0.02_R8P])
   pp(1, 2) = pose_2d(0.17_R8P, [0.01_R8P, 0.03_R8P])    ; pp(2, 2) = pose_2d(-0.6_R8P, [1.4_R8P, -0.1_R8P])
   pp(1, 3) = pose_2d(-0.05_R8P, [0.02_R8P, -0.08_R8P])  ; pp(2, 3) = pose_2d(-0.5_R8P, [0.4_R8P, 0.33_R8P])
   do c = 1, 3
      call fd_closure_2d_bodies(two, x, y, [REACH, REACH], cb, status=st, poses=pp(:, c), lengths=[h, 2._R8P * h], &
                                linear_velocities=vl, angular_velocities=wl, centres=cl, ilw=.true.)
      call stat(st)
      do b = 1, 2
         call fd_closure_2d(two(b), x, y, [REACH, REACH], cs(b), status=st, lengths=[h, 2._R8P * h], pose=pp(b, c), &
                            linear_velocity=vl(:, b), angular_velocity=wl(b), centre=cl(:, b), ilw=.true.)
         call stat(st)
      enddo
      call ghost_index(cs, size(x, kind=I4P), size(y, kind=I4P), 1_I4P, idx)
      nbad = 0_I4P
      do g = 1, cb%n
         b = cb%body(g)
         if (.not. same_ghost(cb, g, cs(b), idx(cb%ijk(1, g), cb%ijk(2, g), 1, b), b)) nbad = nbad + 1_I4P
      enddo
      if (cb%n /= cs(1)%n + cs(2)%n .or. cb%overlap_points /= 0_I4P) nbad = nbad + 1_I4P
      call load_samples_2d_bodies(two, x, y, lb, status=st, poses=pp(:, c), order=3_I4P, lengths=[h, 2._R8P * h])
      call stat(st)
      sbad = 0_I4P ; off = 0_I4P
      do b = 1, 2
         call load_samples_2d(two(b), x, y, ls(b), status=st, order=3_I4P, lengths=[h, 2._R8P * h], pose=pp(b, c))
         call stat(st)
         if (.not. same_slice(lb, off, ls(b), b)) sbad = sbad + 1_I4P
         off = off + ls(b)%n
      enddo
      if (lb%n /= off) sbad = sbad + 1_I4P
      print '(3X,I4,I9,2X,"(",I6,",",I6,")",I14,F25.2,I14,I12)', c, cb%n, cs(1)%n, cs(2)%n, nbad, &
            minval(cb%clearance) / (REACH * h), lb%n, sbad
      call check('locality 2-D, pose pair '//trim(itoa(c))//': every ghost and sample = its single-body run, bit for bit', &
                 nbad == 0_I4P .and. sbad == 0_I4P .and. cb%n > 0_I4P .and. minval(cb%clearance) > 2 * REACH * h)
   enddo
   endsubroutine locality_2d

   subroutine ghost_index(cs, ni, nj, nk, idx)
   !< idx(i, j, k, b) = the ghost of the single-body table cs(b) at (i, j, k), 0 if none.
   type(ib_closure_t),        intent(in)  :: cs(:)          !< Single-body tables.
   integer(I4P),              intent(in)  :: ni, nj, nk     !< Points per axis.
   integer(I4P), allocatable, intent(out) :: idx(:,:,:,:)   !< Index.
   integer(I4P) :: b, g, k                                  !< Counters.

   allocate(idx(ni, nj, nk, size(cs))) ; idx = 0_I4P
   do b = 1, size(cs)
      do g = 1, cs(b)%n
         k = 1_I4P ; if (cs(b)%dims == 3_I4P) k = cs(b)%ijk(3, g)
         idx(cs(b)%ijk(1, g), cs(b)%ijk(2, g), k, b) = g
      enddo
   enddo
   endsubroutine ghost_index

   logical function same_ghost(cb, g, cs, s, b) result(same)
   !< Ghost g of the bodies table cb equals ghost s of body b's single-body table cs, every field by bits (body = b,
   !< sub_body = cs's body; ilw_body = b where a crossing). 2-D: the meaningful slices of kappa and the directions.
   type(ib_closure_t), intent(in) :: cb, cs   !< Bodies and single-body tables.
   integer(I4P),       intent(in) :: g, s, b  !< Ghost, single-body ghost, body.

   same = .false.
   if (s < 1_I4P) return
   same = all(cb%ijk(:, g) == cs%ijk(:, s)) .and. eqr(cb%points(:, g), cs%points(:, s)) .and. cb%mask(g) == cs%mask(s) .and. &
          all(cb%steps(:, g) == cs%steps(:, s))
   same = same .and. eqr([cb%distance(g)], [cs%distance(s)]) .and. eqr(cb%foot(:, g), cs%foot(:, s)) .and. &
          eqr(cb%normal(:, g), cs%normal(:, s)) .and. cb%element(g) == cs%element(s) .and. &
          cb%feature(g) == cs%feature(s) .and. cb%patch(g) == cs%patch(s) .and. cb%body(g) == b .and. &
          cb%sub_body(g) == cs%body(s)
   same = same .and. eqr(flat(cb%probe(:, :, g)), flat(cs%probe(:, :, s))) .and. &
          all(cb%visible(:, g) .eqv. cs%visible(:, s)) .and. all(cb%probe_cell(:, :, g) == cs%probe_cell(:, :, s))
   same = same .and. eqr(cb%velocity(:, g), cs%velocity(:, s)) .and. eqr(cb%ilw_t(:, g), cs%ilw_t(:, s)) .and. &
          all(cb%ilw_element(:, g) == cs%ilw_element(:, s)) .and. &
          eqr(flat(cb%ilw_normal(:, :, g)), flat(cs%ilw_normal(:, :, s))) .and. &
          all((cb%ilw_body(:, g) == b) .eqv. (cs%ilw_element(:, s) > 0_I4P))
   if (cb%dims == 3_I4P) then
      same = same .and. eqr(cb%kappa(:, g), cs%kappa(:, s)) .and. eqr(flat(cb%direction(:, :, g)), flat(cs%direction(:, :, s)))
      same = same .and. eqr(flat(cb%ilw_kappa(:, :, g)), flat(cs%ilw_kappa(:, :, s))) .and. &
             eqr(reshape(cb%ilw_direction(:, :, :, g), [size(cb%ilw_direction(:, :, :, g))]), &
                 reshape(cs%ilw_direction(:, :, :, s), [size(cs%ilw_direction(:, :, :, s))]))
   else
      same = same .and. eqr(cb%kappa(:, g), cs%kappa(:, s)) .and. eqr(cb%direction(:, 1, g), cs%direction(:, 1, s)) .and. &
             eqr(cb%ilw_kappa(1, :, g), cs%ilw_kappa(1, :, s)) .and. &
             eqr(flat(cb%ilw_direction(:, 1, :, g)), flat(cs%ilw_direction(:, 1, :, s)))
   endif
   endfunction same_ghost

   logical function same_slice(c, off, a, b) result(same)
   !< The slots off+1 ... off+a%n of the bodies set c equal the single-body set a, every field by bits (body = b,
   !< sub_body = a's body).
   type(load_samples_t), intent(in) :: c, a    !< Bodies set, single-body set.
   integer(I4P),         intent(in) :: off, b  !< Offset, body.
   integer(I4P) :: n                           !< Samples.

   n = a%n
   same = .false.
   if (c%n < off + n) return
   same = eqr(flat(c%points(:, off + 1:off + n)), flat(a%points)) .and. eqr(c%weight(off + 1:off + n), a%weight) .and. &
          eqr(flat(c%normal(:, off + 1:off + n)), flat(a%normal)) .and. &
          all(c%element(off + 1:off + n) == a%element) .and. all(c%patch(off + 1:off + n) == a%patch) .and. &
          all(c%body(off + 1:off + n) == b) .and. all(c%sub_body(off + 1:off + n) == a%body) .and. &
          eqr(reshape(c%probe(:, :, off + 1:off + n), [size(a%probe)]), reshape(a%probe, [size(a%probe)])) .and. &
          all(c%visible(:, off + 1:off + n) .eqv. a%visible) .and. all(c%probe_cell(:, :, off + 1:off + n) == a%probe_cell)
   endfunction same_slice

   ! 2. narrow gap

   subroutine gap_3d()
   !< Part 2, 3-D: two spheres at gaps 0.25 .. 4 h.
   type(surface_stl_object)        :: sph, two(2)
   type(facet_object), allocatable :: fc(:)
   type(facet_object), pointer     :: fp
   type(rigid_pose_t)              :: pp(2)
   type(ib_closure_t)              :: cb
   real(R8P), allocatable          :: x(:), y(:)
   real(R8P)                       :: h, sag, nn(3), ctr(3, 2), dl, a
   integer(I4P)                    :: k, f, st, g

   h = 0.1_R8P
   call make_sphere(4_I4P, fc) ; call sph%adopt_facets(facets=fc)
   sag = 0._R8P
   do f = 1, sph%get_facets_number()
      fp => sph%facet_at(f)
      nn = [fp%normal%x, fp%normal%y, fp%normal%z] ; nn = nn / norm2(nn)
      sag = max(sag, 1._R8P - abs(dot_product(nn, [fp%vertex(1)%x, fp%vertex(1)%y, fp%vertex(1)%z])))
   enddo
   two(1) = sph ; two(2) = sph
   x = [(-2.65_R8P + h * g, g = 0, 53)] ; y = [(-1.25_R8P + h * g, g = 0, 25)]
   print '(A)', '2. Narrow gap: shared ghosts, clearance and visibility against the analytic sphere / circle'
   print '(A,ES9.2,A)', '   3-D two spheres R = 1 (5 120 facets, sag s_B = ', sag, '), h = 0.1, reach 2 (sd = 2 h), probes h, 2 h, 3 h'
   call gap_header
   do k = 1, size(GAPS)
      dl = GAPS(k) * h ; a = 1._R8P + 0.5_R8P * dl
      ctr(:, 1) = [-a, 0.0137_R8P, -0.0071_R8P] ; ctr(:, 2) = [a, 0.0137_R8P, -0.0071_R8P]
      pp(1) = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, ctr(:, 1))
      pp(2) = pose_from_axis_angle([0._R8P, 1._R8P, -1._R8P], -0.9_R8P, ctr(:, 2))
      call fd_closure_3d_bodies(two, x, y, y, [REACH, REACH, REACH], cb, status=st, poses=pp, lengths=[h, 2 * h, 3 * h])
      call stat(st)
      call gap_check(cb, ctr, [1._R8P, 1._R8P], sag, dl, h, GAPS(k), '3-D')
   enddo
   endsubroutine gap_3d

   subroutine gap_2d()
   !< Part 2, 2-D: two circles at gaps 0.25 .. 4 h.
   integer(I4P), parameter :: NV = 1024_I4P
   real(R8P),    parameter :: R = 0.5_R8P
   type(curve_object)     :: circ, two(2)
   type(rigid_pose_2d_t)  :: pp(2)
   type(ib_closure_t)     :: cb
   real(R8P), allocatable :: x(:), y(:), u(:), v(:)
   real(R8P)              :: h, sag, ctr(2, 2), dl, a
   integer(I4P)           :: k, st, g

   h = 0.02_R8P
   u = [(R * cos(2 * PI * g / NV), g = 0, NV - 1)] ; v = [(R * sin(2 * PI * g / NV), g = 0, NV - 1)]
   call circ%adopt_loops(u=u, v=v, loop_sizes=[NV])
   sag = R * (1._R8P - cos(PI / NV))
   two(1) = circ ; two(2) = circ
   x = [(-1.21_R8P + h * g, g = 0, 121)] ; y = [(-0.63_R8P + h * g, g = 0, 63)]
   print '(A,ES9.2,A)', '   2-D two circles R = 0.5 (1 024 segments, sag s_B = ', sag, '), h = 0.02, reach 2 (sd = 2 h), probes h, 2 h, 3 h'
   call gap_header
   do k = 1, size(GAPS)
      dl = GAPS(k) * h ; a = R + 0.5_R8P * dl
      ctr(:, 1) = [-a, 0.0037_R8P] ; ctr(:, 2) = [a, 0.0037_R8P]
      pp(1) = pose_2d(0.3_R8P, ctr(:, 1)) ; pp(2) = pose_2d(-1.1_R8P, ctr(:, 2))
      call fd_closure_2d_bodies(two, x, y, [REACH, REACH], cb, status=st, poses=pp, lengths=[h, 2 * h, 3 * h])
      call stat(st)
      call gap_check(cb, ctr, [R, R], sag, dl, h, GAPS(k), '2-D')
   enddo
   endsubroutine gap_2d

   subroutine gap_header()
   !< Header of the narrow-gap tables.
   print '(A)', '   gap/h  ghosts  shared  N(d<sd-s)  N(d<sd)  inconsistent  clearance - analytic [min, max]   '// &
                '(least clearance - gap)/h  blocked  [N(meet R-s), N(meet R)]  against'
   endsubroutine gap_header

   subroutine gap_check(cb, ctr, rr, sag, dl, h, gh, label)
   !< The narrow-gap checks of one table: shared count against the analytic bracket, clearance against the analytic
   !< clearance, least clearance against the gap, blocked probes against the analytic bracket.
   type(ib_closure_t), intent(in) :: cb         !< Bodies table.
   real(R8P),          intent(in) :: ctr(:,:)   !< Centres [d, 2].
   real(R8P),          intent(in) :: rr(2)      !< Radii.
   real(R8P),          intent(in) :: sag        !< Sag of the other body.
   real(R8P),          intent(in) :: dl, h, gh  !< Gap, spacing, gap / h.
   character(*),       intent(in) :: label      !< 2-D / 3-D.
   real(R8P)    :: an, e, emin, emax, sd, dseg
   integer(I4P) :: g, o, j, nsh, nlo, nhi, ninc, nbl, nblo, nbhi, nbad

   sd = cb%shared_distance
   nsh = 0_I4P ; nlo = 0_I4P ; nhi = 0_I4P ; ninc = 0_I4P ; nbl = 0_I4P ; nblo = 0_I4P ; nbhi = 0_I4P ; nbad = 0_I4P
   emin = huge(1._R8P) ; emax = -huge(1._R8P)
   do g = 1, cb%n
      o = 3_I4P - cb%body(g)
      an = norm2(cb%foot(:, g) - ctr(:, o)) - rr(o)
      e = cb%clearance(g) - an ; emin = min(emin, e) ; emax = max(emax, e)
      if (cb%shared(g)) nsh = nsh + 1_I4P
      if (an < sd - sag) nlo = nlo + 1_I4P
      if (an < sd) nhi = nhi + 1_I4P
      if (cb%shared(g) .and. an >= sd) ninc = ninc + 1_I4P
      if (.not. cb%shared(g) .and. an < sd - sag) ninc = ninc + 1_I4P
      if (cb%second_body(g) /= o) ninc = ninc + 1_I4P
      do j = 1, cb%m
         dseg = seg_dist(ctr(:, o), cb%foot(:, g), cb%probe(:, j, g))
         if (.not. cb%visible(j, g)) nbl = nbl + 1_I4P
         if (dseg < rr(o) - sag) nblo = nblo + 1_I4P
         if (dseg <= rr(o)) nbhi = nbhi + 1_I4P
         if (.not. cb%visible(j, g) .and. dseg > rr(o)) nbad = nbad + 1_I4P
         if (cb%visible(j, g) .and. dseg < rr(o) - sag) nbad = nbad + 1_I4P
      enddo
   enddo
   print '(3X,F5.2,I8,I8,I11,I9,I14,2X,"[",ES10.2,",",ES10.2,"]",F27.4,I9,5X,"[",I6,",",I6,"]",I9)', gh, cb%n, nsh, nlo, &
         nhi, ninc, emin, emax, (minval(cb%clearance) - dl) / h, nbl, nblo, nbhi, nbad
   call check('narrow gap '//label//', gap '//trim(ftoa(gh))//' h: shared within the analytic bracket, clearance, '// &
              'least clearance, visibility', cb%n > 0_I4P .and. nlo <= nsh .and. nsh <= nhi .and. ninc == 0_I4P .and. &
              emin >= -TOL .and. emax <= sag + TOL .and. minval(cb%clearance) >= dl - TOL .and. nbad == 0_I4P .and. &
              cb%overlap_points == 0_I4P)
   endsubroutine gap_check

   pure real(R8P) function seg_dist(c, a, b) result(d)
   !< Distance from c to the segment [a, b].
   real(R8P), intent(in) :: c(:), a(:), b(:) !< Points.
   real(R8P) :: t, ab2

   ab2 = sum((b - a)**2)
   t = 0._R8P ; if (ab2 > 0._R8P) t = max(0._R8P, min(1._R8P, dot_product(c - a, b - a) / ab2))
   d = norm2(c - (a + t * (b - a)))
   endfunction seg_dist

   ! 3. loads through the union

   subroutine loads_3d()
   !< Part 3, 3-D: buoyancy per body through the cut of the union and through the load samples.
   type(surface_stl_object)        :: sph, ell, cube, two(2)
   type(facet_object), allocatable :: fc(:)
   type(rigid_pose_t)              :: pp(2)
   real(R8P), parameter            :: H(2) = [0.13_R8P, 0.125_R8P]
   real(R8P), allocatable          :: gx(:)
   integer(I4P)                    :: k, i

   print '(A)', '3. Loads through the union, p = p0 + g . x: per body F = -g V_k, M = V_k (x_c,k - c) x (-g); error / (eps S)'
   print '(A)', '   case                                        body  wall terms  untrusted  contact terms    FV F  FV M  '// &
                'FD q1 F  q2 F  q2 M  q3 F  q3 M'
   call make_sphere(4_I4P, fc) ; call sph%adopt_facets(facets=fc)
   ell = sph ; call ell%resize(factor=vector_R8P(1._R8P, 0.7_R8P, 0.5_R8P)) ; call ell%analyze
   call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
   two(1) = sph ; two(2) = sph
   pp(1) = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [-1.0125_R8P, 0.0137_R8P, -0.0071_R8P])
   pp(2) = pose_from_axis_angle([0._R8P, 1._R8P, -1._R8P], -0.9_R8P, [1.0125_R8P, 0.0137_R8P, -0.0071_R8P])
   gx = [(-2.13_R8P + 0.1_R8P * i, i = 0, 43)]
   call buoyancy_3d(two, pp, gx, [(-1.11_R8P + 0.1_R8P * i, i = 0, 22)], 'two spheres, gap 0.25 h (h = 0.1)', 0._R8P)
   two(2) = ell
   pp(2) = pose_from_axis_angle([0._R8P, 1._R8P, -1._R8P], -0.9_R8P, [1.05_R8P, 0.0137_R8P, -0.0071_R8P])
   call buoyancy_3d(two, pp, gx, [(-1.11_R8P + 0.1_R8P * i, i = 0, 22)], 'sphere + ellipsoid, gap >= 0.5 h', 0._R8P)
   two(1) = cube ; two(2) = cube
   do k = 1, 2
      pp(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P)
      pp(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [1._R8P, 0._R8P, 0._R8P])
      gx = [(-0.27_R8P + H(k) * i, i = 0, nint(2.6_R8P / H(k)))]
      call buoyancy_3d(two, pp, gx, gx, 'cubes in face contact, h = '//trim(ftoa(H(k))), 1._R8P)
      pp(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [1._R8P, 0.37_R8P, 0.21_R8P])
      call buoyancy_3d(two, pp, gx, gx, 'cubes in face contact, shifted, h = '//trim(ftoa(H(k))), 1._R8P)
   enddo
   endsubroutine loads_3d

   subroutine buoyancy_3d(two, pp, gx, gy, name, xcontact)
   !< Per-body buoyancy of two posed bodies through both paths.
   type(surface_stl_object), intent(in) :: two(2)     !< Bodies (body frame).
   type(rigid_pose_t),       intent(in) :: pp(2)      !< Poses.
   real(R8P),                intent(in) :: gx(:), gy(:) !< Cell edges along x; along y and z.
   character(*),             intent(in) :: name       !< Label.
   real(R8P),                intent(in) :: xcontact   !< Plane x of a face contact (0: none).
   type(cut_cells_3d_t)      :: cut
   type(cut_walls_3d_t)      :: walls
   type(piece_bodies_t)      :: pb
   type(surface_stl_object)  :: mg
   type(load_samples_t)      :: ls
   real(R8P), allocatable    :: p(:), pk(:), gk(:,:), gyr(:,:), pf(:,:), pm(:,:)
   real(R8P)                 :: v(2), xc(3, 2), fr(3), mr(3), f(3), m(3), sf, sm, kg(3), efv, emv, ef(3, 2), em(3, 2)
   integer(I4P)              :: b, w, q, i, st, nun, ncon(2), nw(2)
   logical                   :: good

   do b = 1, 2
      call two(b)%mass_properties(volume=v(b), center=xc(:, b), status=st) ; call stat(st)
      xc(:, b) = to_world(pp(b), xc(:, b))
   enddo
   call cut_cells_3d_bodies(two, pp, gx, gy, gy, cut, status=st, walls=walls, piece_bodies=pb, merged=mg) ; call stat(st)
   nun = 0_I4P ; if (allocated(cut%pieces)) nun = int(count(cut%pieces < 0_I4P), I4P)
   allocate(p(walls%n))
   do w = 1, walls%n
      p(w) = P0 + dot_product(G3, walls%centroid(:, w))
   enddo
   call cut_walls_gyration(walls, mg, gx, gy, gy, gyr, status=st) ; call stat(st)
   ! FD path: per body by integrate_loads
   do q = 1, 3
      call load_samples_3d_bodies(two, samples=ls, status=st, poses=pp, order=q) ; call stat(st)
      if (allocated(pk)) deallocate(pk)
      allocate(pk(ls%n))
      do i = 1, ls%n
         pk(i) = P0 + dot_product(G3, ls%points(:, i))
      enddo
      call integrate_loads(ls, f, m, status=st, pressure=pk, about=C3, per_body_force=pf, per_body_moment=pm)
      call stat(st)
      do b = 1, 2
         fr = -G3 * v(b) ; mr = cross(v(b) * (xc(:, b) - C3), -G3)
         sf = v(b) * norm2(G3) ; sm = v(b) * norm2(G3) * norm2(xc(:, b) - C3)
         do i = 1, ls%n
            if (ls%body(i) /= b) cycle
            sf = sf + ls%weight(i) * abs(pk(i)) ; sm = sm + ls%weight(i) * abs(pk(i)) * norm2(ls%points(:, i) - C3)
         enddo
         ef(q, b) = maxval(abs(pf(:, b) - fr)) / (EPS * sf) ; em(q, b) = maxval(abs(pm(:, b) - mr)) / (EPS * sm)
      enddo
   enddo
   do b = 1, 2
      fr = -G3 * v(b) ; mr = cross(v(b) * (xc(:, b) - C3), -G3)
      if (allocated(pk)) deallocate(pk)
      allocate(pk(walls%n), gk(3, walls%n))
      sf = v(b) * norm2(G3) ; sm = v(b) * norm2(G3) * norm2(xc(:, b) - C3) ; ncon(b) = 0_I4P ; nw(b) = 0_I4P
      do w = 1, walls%n
         if (pb%wall_body(w) == b) then
            pk(w) = p(w) ; gk(:, w) = G3 ; nw(b) = nw(b) + 1_I4P
            kg = [gyr(1, w) * G3(1) + gyr(4, w) * G3(2) + gyr(6, w) * G3(3), &
                  gyr(4, w) * G3(1) + gyr(2, w) * G3(2) + gyr(5, w) * G3(3), &
                  gyr(6, w) * G3(1) + gyr(5, w) * G3(2) + gyr(3, w) * G3(3)]
            sf = sf + abs(p(w)) * walls%wet_area(w)
            sm = sm + abs(p(w)) * walls%wet_area(w) * norm2(walls%centroid(:, w) - C3) + norm2(kg) * walls%wet_area(w)
            if (xcontact > 0._R8P .and. abs(walls%centroid(1, w) - xcontact) < 1.e-12_R8P) ncon(b) = ncon(b) + 1_I4P
         else
            pk(w) = 0._R8P ; gk(:, w) = 0._R8P
         endif
      enddo
      call cut_loads_3d(walls, f, m, status=st, pressure=pk, gradient=gk, gyration=gyr, about=C3) ; call stat(st)
      efv = maxval(abs(f - fr)) / (EPS * sf) ; emv = maxval(abs(m - mr)) / (EPS * sm)
      deallocate(pk, gk)
      print '(3X,A,I4,I12,I11,I15,2F6.2,F9.2,4F6.2)', pad44(name), b, nw(b), nun, ncon(b), efv, emv, ef(1, b), ef(2, b), &
            em(2, b), ef(3, b), em(3, b)
      good = efv <= 1._R8P .and. emv <= 1._R8P .and. all(ef(:, b) <= 1._R8P) .and. em(2, b) <= 1._R8P .and. &
             em(3, b) <= 1._R8P .and. nun == 0_I4P .and. nw(b) > 0_I4P
      if (xcontact > 0._R8P) good = good .and. ncon(b) > 0_I4P
      call check(name//', body '//trim(itoa(b))//': F and M per body within 1 eps S (FV, FD q = 1-3), every cell trusted', &
                 good)
   enddo
   if (walls%n /= sum(nw)) call check(name//': every wall term labelled with a body', .false.)
   endsubroutine buoyancy_3d

   subroutine loads_2d()
   !< Part 3, 2-D: buoyancy per body through the cut of the union and through the load samples.
   integer(I4P), parameter :: NV = 1024_I4P
   real(R8P),    parameter :: R = 0.5_R8P, H(2) = [0.13_R8P, 0.125_R8P]
   type(curve_object)      :: circ, naca, ell, sq, two(2)
   type(rigid_pose_2d_t)   :: pp(2)
   real(R8P), allocatable  :: u(:), v(:), gx(:)
   real(R8P)               :: lo
   integer(I4P)            :: g, k

   print '(A)', '   2-D (per unit depth; M the scalar moment)'
   u = [(R * cos(2 * PI * g / NV), g = 0, NV - 1)] ; v = [(R * sin(2 * PI * g / NV), g = 0, NV - 1)]
   call circ%adopt_loops(u=u, v=v, loop_sizes=[NV])
   two(1) = circ ; two(2) = circ
   pp(1) = pose_2d(0.3_R8P, [-0.5025_R8P, 0.0037_R8P]) ; pp(2) = pose_2d(-1.1_R8P, [0.5025_R8P, 0.0037_R8P])
   gx = [(-1.113_R8P + 0.02_R8P * g, g = 0, 111)]
   call buoyancy_2d(two, pp, gx, [(-0.617_R8P + 0.02_R8P * g, g = 0, 61)], 'two circles, gap 0.25 h (h = 0.02)', 0._R8P)
   call naca0012(400_I4P, u, v)
   call naca%adopt_loops(u=u, v=v, loop_sizes=[size(u, kind=I4P)])
   u = [(0.3_R8P * cos(2 * PI * g / 512), g = 0, 511)] ; v = [(0.12_R8P * sin(2 * PI * g / 512), g = 0, 511)]
   call ell%adopt_loops(u=u, v=v, loop_sizes=[512_I4P])
   two(1) = naca ; two(2) = ell
   lo = sqrt(0.09_R8P * sin(0.5_R8P)**2 + 0.0144_R8P * cos(0.5_R8P)**2) ! lowest point of the turned ellipse below its centre
   pp(1) = pose_2d(0._R8P, [0._R8P, 0._R8P]) ; pp(2) = pose_2d(-0.5_R8P, [0.31_R8P, 0.0601_R8P + lo + 0.005_R8P])
   gx = [(-0.113_R8P + 0.01_R8P * g, g = 0, 125)]
   call buoyancy_2d(two, pp, gx, [(-0.117_R8P + 0.01_R8P * g, g = 0, 57)], 'NACA 0012 + ellipse, gap >= 0.5 h', 0._R8P)
   call sq%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
   two(1) = sq ; two(2) = sq
   do k = 1, 2
      gx = [(-0.27_R8P + H(k) * g, g = 0, nint(2.6_R8P / H(k)))]
      pp(1) = pose_2d(0._R8P, [0._R8P, 0._R8P]) ; pp(2) = pose_2d(0._R8P, [1._R8P, 0._R8P])
      call buoyancy_2d(two, pp, gx, gx, 'squares in edge contact, h = '//trim(ftoa(H(k))), 1._R8P)
      pp(2) = pose_2d(0._R8P, [1._R8P, 0.37_R8P])
      call buoyancy_2d(two, pp, gx, gx, 'squares in edge contact, shifted, h = '//trim(ftoa(H(k))), 1._R8P)
   enddo
   endsubroutine loads_2d

   subroutine buoyancy_2d(two, pp, gx, gy, name, xcontact)
   !< Per-body buoyancy of two posed curves through both paths.
   type(curve_object),    intent(in) :: two(2)       !< Bodies (body frame).
   type(rigid_pose_2d_t), intent(in) :: pp(2)        !< Poses.
   real(R8P),             intent(in) :: gx(:), gy(:) !< Cell edges (u, v).
   character(*),          intent(in) :: name         !< Label.
   real(R8P),             intent(in) :: xcontact     !< Line u of an edge contact (0: none).
   type(cut_cells_2d_t)      :: cut
   type(cut_walls_2d_t)      :: walls
   type(piece_bodies_t)      :: pb
   type(load_samples_t)      :: ls
   real(R8P), allocatable    :: p(:), pk(:), gk(:,:), pf(:,:), pm(:)
   real(R8P)                 :: aa(2), xc(2, 2), fr(2), mr, f(2), m, sf, sm, d(2), efv, emv, ef(3, 2), em(3, 2)
   integer(I4P)              :: b, w, q, i, st, nun, ncon(2), nw(2)
   logical                   :: good

   do b = 1, 2
      call two(b)%mass_properties(area=aa(b), center=xc(:, b), status=st) ; call stat(st)
      xc(:, b) = to_world(pp(b), xc(:, b))
   enddo
   call cut_cells_2d_bodies(two, pp, gx, gy, cut, status=st, piece_bodies=pb, walls=walls) ; call stat(st)
   nun = 0_I4P ; if (allocated(cut%pieces)) nun = int(count(cut%pieces < 0_I4P), I4P)
   allocate(p(walls%n))
   do w = 1, walls%n
      p(w) = P0 + dot_product(G2, walls%centroid(:, w))
   enddo
   do q = 1, 3
      call load_samples_2d_bodies(two, samples=ls, status=st, poses=pp, order=q) ; call stat(st)
      if (allocated(pk)) deallocate(pk)
      allocate(pk(ls%n))
      do i = 1, ls%n
         pk(i) = P0 + dot_product(G2, ls%points(:, i))
      enddo
      call integrate_loads(ls, f, m, status=st, pressure=pk, about=C2, per_body_force=pf, per_body_moment=pm)
      call stat(st)
      do b = 1, 2
         fr = -G2 * aa(b) ; mr = aa(b) * cross2(xc(:, b) - C2, -G2)
         sf = aa(b) * norm2(G2) ; sm = aa(b) * norm2(G2) * norm2(xc(:, b) - C2)
         do i = 1, ls%n
            if (ls%body(i) /= b) cycle
            sf = sf + ls%weight(i) * abs(pk(i)) ; sm = sm + ls%weight(i) * abs(pk(i)) * norm2(ls%points(:, i) - C2)
         enddo
         ef(q, b) = maxval(abs(pf(:, b) - fr)) / (EPS * sf) ; em(q, b) = abs(pm(b) - mr) / (EPS * sm)
      enddo
   enddo
   do b = 1, 2
      fr = -G2 * aa(b) ; mr = aa(b) * cross2(xc(:, b) - C2, -G2)
      if (allocated(pk)) deallocate(pk)
      allocate(pk(walls%n), gk(2, walls%n))
      sf = aa(b) * norm2(G2) ; sm = aa(b) * norm2(G2) * norm2(xc(:, b) - C2) ; ncon(b) = 0_I4P ; nw(b) = 0_I4P
      do w = 1, walls%n
         if (pb%wall_body(w) == b) then
            pk(w) = p(w) ; gk(:, w) = G2 ; nw(b) = nw(b) + 1_I4P
            d = walls%ends(3:4, w) - walls%ends(1:2, w)
            sf = sf + abs(p(w)) * walls%length(w)
            sm = sm + abs(p(w)) * walls%length(w) * norm2(walls%centroid(:, w) - C2) + &
                 norm2(d) * abs(dot_product(d, G2)) / 12._R8P * walls%length(w)
            if (xcontact > 0._R8P .and. abs(walls%centroid(1, w) - xcontact) < 1.e-12_R8P) ncon(b) = ncon(b) + 1_I4P
         else
            pk(w) = 0._R8P ; gk(:, w) = 0._R8P
         endif
      enddo
      call cut_loads_2d(walls, f, m, status=st, pressure=pk, gradient=gk, about=C2) ; call stat(st)
      efv = maxval(abs(f - fr)) / (EPS * sf) ; emv = abs(m - mr) / (EPS * sm)
      deallocate(pk, gk)
      print '(3X,A,I4,I12,I11,I15,2F6.2,F9.2,4F6.2)', pad44(name), b, nw(b), nun, ncon(b), efv, emv, ef(1, b), ef(2, b), &
            em(2, b), ef(3, b), em(3, b)
      good = efv <= 1._R8P .and. emv <= 1._R8P .and. all(ef(:, b) <= 1._R8P) .and. em(2, b) <= 1._R8P .and. &
             em(3, b) <= 1._R8P .and. nun == 0_I4P .and. nw(b) > 0_I4P
      if (xcontact > 0._R8P) good = good .and. ncon(b) > 0_I4P
      call check(name//', body '//trim(itoa(b))//': F and M per body within 1 eps S (FV, FD q = 1-3), every cell trusted', &
                 good)
   enddo
   if (walls%n /= sum(nw)) call check(name//': every wall segment labelled with a body', .false.)
   endsubroutine buoyancy_2d

   ! helpers

   subroutine naca0012(mh, xx, yy)
   !< NACA 0012 of unit chord, closed trailing edge, 2 mh vertices, cosine spacing, counter-clockwise.
   integer(I4P),           intent(in)  :: mh           !< Points per side.
   real(R8P), allocatable, intent(out) :: xx(:), yy(:) !< Vertices.
   real(R8P)    :: xi
   integer(I4P) :: q, j

   allocate(xx(2 * mh), yy(2 * mh)) ; j = 0
   do q = mh, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = yt(xi)
   enddo
   do q = 1, mh - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = -yt(xi)
   enddo
   endsubroutine naca0012

   pure real(R8P) function yt(xx)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: xx !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.1260_R8P * xx - 0.3516_R8P * xx ** 2 + 0.2843_R8P * xx ** 3 - &
                   0.1036_R8P * xx ** 4)
   endfunction yt

   function pad44(s) result(r)
   !< s left-justified in 44 characters.
   character(*), intent(in) :: s !< String.
   character(44)            :: r !< Padded.

   r = s
   endfunction pad44

   pure function flat(a) result(r)
   !< A rank-2 array as rank 1.
   real(R8P), intent(in) :: a(:,:)          !< Array.
   real(R8P)             :: r(size(a))      !< Flattened.

   r = reshape(a, [size(a)])
   endfunction flat

   logical function eqr(a, b)
   !< Equal bits (and sizes).
   real(R8P), intent(in) :: a(:), b(:) !< Reals.

   eqr = size(a) == size(b)
   if (eqr) eqr = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
   endfunction eqr

   pure function cross(u, v) result(w)
   !< Cross product.
   real(R8P), intent(in) :: u(3), v(3) !< Vectors.
   real(R8P)             :: w(3)

   w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
   endfunction cross

   pure real(R8P) function cross2(u, v)
   !< z component of u x v.
   real(R8P), intent(in) :: u(2), v(2) !< Vectors.

   cross2 = u(1) * v(2) - u(2) * v(1)
   endfunction cross2

   function itoa(i) result(s)
   !< Integer to string.
   integer(I4P), intent(in) :: i !< Integer.
   character(16)            :: s !< String.

   write(s, '(I0)') i
   endfunction itoa

   function ftoa(x) result(s)
   !< Short real to string.
   real(R8P), intent(in) :: x !< Real.
   character(16)          :: s !< String.

   write(s, '(G0.4)') x
   s = adjustl(s)
   endfunction ftoa

   subroutine stat(st)
   !< Count a refused call.
   integer(I4P), intent(in) :: st !< Status.

   if (st /= STATUS_OK) nstat = nstat + 1_I4P
   endsubroutine stat

   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (cond) then
      print '(A)', 'ok: '//label
   else
      print '(A)', 'FAILED: '//label ; ok = .false.
   endif
   endsubroutine check
endprogram fossil_val_bodies
