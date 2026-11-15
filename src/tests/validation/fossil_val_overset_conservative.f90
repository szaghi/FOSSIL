!< FOSSIL validation: a passive scalar through a conservative overset interface (issue #106 step 1, solver gate).
program fossil_val_overset_conservative
!< FOSSIL validation: a passive scalar through a conservative overset interface (issue #106 step 1, solver gate).
!<
!< A minor ring block (a cylinder O-grid, r in [1, 2]; or a square ring on grid lines, with degenerate pieces) in a
!< Cartesian major grid on [-3, 3]^2. A scalar blob is advected by the potential flow around the cylinder (uniform flow
!< for the square) with first-order upwind finite volumes and forward Euler, two ways:
!< A. Wang's patched interface (`overset_patch_2d`): the major grid is cut by the minor boundary (holes, cut pieces with
!<    their apertures W_f), the minor cells keep their faces, and the two grids exchange the upwind fluxes through the
!<    smallest faces, each flux added to one side and subtracted from the other. Small major pieces (area below half a
!<    cell, `merge_targets`) and degenerate ones (`merge_into`) are merged into their targets: a group shares one value
!<    and the sum of its members' fluxes. The total mass plus the mass that left through the box changes by round-off.
!< B. The same assembly coupled by interpolation (Chimera): the major cells whose centre lies at r < 1.4 are holes, at
!<    1.4 <= r < 1.6 receptors (bilinear from the minor cell centres, located by FOSSIL in the block of the cell centres,
!<    `cell_centre_block_2d` + `locate_many`), the minor's outer row of cells receives from the major cell centres
!<    (bilinear); each grid advances its own field cells. The composite mass (the major cells with centre at r >= 2, all
!<    minor cells) plus the mass that left the box drifts by the interpolation.
!< The inner minor boundary is a wall (no flux). Exit status 0 if A conserves to 1e-12 relative and stays in [0, 1]
!< (upwind within its CFL limit), and B's defect exceeds 1e-6.
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d, &
                   OVERSET_PATCH_OK, cut_links_t, cut_links_2d, merge_t, merge_targets, MERGE_OK, MERGE_CHAIN, &
                   cell_centre_block_2d, LOC_OUTSIDE
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P) !< Pi.
real(R8P), parameter :: ANG = 0.2_R8P      !< Free-stream angle.
type(curvilinear_block_2d_t) :: blk(1)     !< Minor block.
real(R8P), allocatable       :: x(:), pu(:), pv(:) !< Major grid, square.
real(R8P)                    :: da, db, mx, mn  !< Defects, extremes.
integer(I4P)                 :: m, k, n     !< Counters.
logical                      :: ok          !< Checks.

ok = .true.
do m = 1, 2
   call ogrid(m, blk(1))
   x = grid_1d(-3._R8P, 3._R8P, 60_I4P * m)
   call patched(blk, x, 1._R8P, 0.008_R8P / m, 400_I4P * m, da, mx, mn)
   call chimera(blk(1), x, 1._R8P, 0.008_R8P / m, 400_I4P * m, db)
   print '(A,I0,A,I0,A,ES10.2,A,ES10.2,A,2ES10.2)', '  cylinder, m = ', m, ' (', 60 * m, &
         '^2): mass defect A (patched) ', da, ', B (interpolated) ', db, '; A min, max ', mn, mx
   ok = ok .and. da <= 1e-12_R8P .and. db > 1e-6_R8P .and. mx <= 1._R8P + 1e-12_R8P .and. mn >= -1e-12_R8P
enddo
! a square ring on grid lines: degenerate pieces along its sides, folded by merge_into
x = grid_1d(-3._R8P, 3._R8P, 60_I4P)
allocate(pu(120), pv(120)) ; n = 0
do k = 16, 45
   n = n + 1 ; pu(n) = x(k) ; pv(n) = x(16)
enddo
do k = 16, 45
   n = n + 1 ; pu(n) = x(46) ; pv(n) = x(k)
enddo
do k = 46, 17, -1
   n = n + 1 ; pu(n) = x(k) ; pv(n) = x(46)
enddo
do k = 46, 17, -1
   n = n + 1 ; pu(n) = x(16) ; pv(n) = x(k)
enddo
call ring_block(pu, pv, 3_I4P, blk(1))
call patched(blk, x, 0._R8P, 0.008_R8P, 400_I4P, da, mx, mn)
print '(A,ES10.2,A,2ES10.2)', '  square ring on grid lines (uniform flow): mass defect A ', da, '; min, max ', mn, mx
ok = ok .and. da <= 1e-12_R8P .and. mx <= 1._R8P + 1e-12_R8P .and. mn >= -1e-12_R8P
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine patched(blk, x, rc, dt, nsteps, defect, vmax, vmin)
   !< A: the scalar through Wang's patched interface; defect = max over the steps of |M + out - M0| / M0.
   type(curvilinear_block_2d_t), intent(in)  :: blk(:)  !< Minor block.
   real(R8P),                    intent(in)  :: x(:)    !< Major grid edges (both directions).
   real(R8P),                    intent(in)  :: rc      !< Cylinder radius of the flow (0: uniform).
   real(R8P),                    intent(in)  :: dt      !< Time step.
   integer(I4P),                 intent(in)  :: nsteps  !< Steps.
   real(R8P),                    intent(out) :: defect  !< Mass defect.
   real(R8P),                    intent(out) :: vmax, vmin !< Extreme values reached.
   type(overset_patch_2d_t)  :: pt                     !< Interface.
   type(cut_links_t)         :: links                  !< Piece links.
   type(merge_t)             :: mt                     !< Merge targets.
   integer(I4P), allocatable :: g(:)                   !< Group (root piece) of each piece (0: an empty piece).
   real(R8P),    allocatable :: vg(:), phi(:), r(:), mg(:) !< Group volumes, values, residuals, masses.
   real(R8P),    allocatable :: vm(:,:), pm(:,:), rm(:,:) !< Minor volumes, values, residuals.
   real(R8P)    :: out, m0, mtot, f, fl, mid(2), nrm(2), len !< Work.
   integer(I4P) :: st, np, p, q, k, i, j, ni, nj, step, nx, t, nsm, ndg, nal, ia, ib !< Work.

   call overset_patch_2d(blk, [1_I4P], x, x, pt, status=st)
   if (st /= OVERSET_PATCH_OK) error stop 'overset_patch_2d'
   np = pt%cut%np ; nx = pt%cut%nx
   call cut_links_2d(pt%cut, links)
   call merge_targets(links, pt%cut%piece_cell, pt%cut%piece_area, pt%cut%piece_wet_normal, x, x, mt)
   ! groups: degenerate pieces into merge_into (then its final target), small pieces into their final target
   allocate(g(np), vg(np), phi(np), r(np), mg(np)) ; nsm = 0 ; ndg = 0 ; nal = 0
   do p = 1, np
      g(p) = p
      if (pt%degenerate(p)) then
         t = pt%merge_into(p)
         g(p) = 0
         if (t > 0) then
            ndg = ndg + 1 ; g(p) = t
            if ((mt%status(t) == MERGE_OK .or. mt%status(t) == MERGE_CHAIN) .and. mt%final_target(t) > 0) &
               g(p) = mt%final_target(t)
         endif
      elseif (mt%status(p) == MERGE_OK .or. mt%status(p) == MERGE_CHAIN) then
         if (mt%final_target(p) > 0) then
            nsm = nsm + 1 ; g(p) = mt%final_target(p)
         else
            nal = nal + 1
         endif
      elseif (mt%status(p) /= 0) then
         nal = nal + 1
      endif
   enddo
   vg = 0._R8P ; mg = 0._R8P ; phi = 0._R8P
   do p = 1, np
      if (g(p) == 0) cycle
      vg(g(p)) = vg(g(p)) + pt%cut%piece_area(p)
      mg(g(p)) = mg(g(p)) + pt%cut%piece_area(p) * blob(pt%cut%piece_centroid(:, p))
   enddo
   do p = 1, np
      if (g(p) > 0) phi(p) = mg(g(p)) / vg(g(p))
   enddo
   ni = blk(1)%ni ; nj = blk(1)%nj
   allocate(vm(ni - 1, nj - 1), pm(ni - 1, nj - 1), rm(ni - 1, nj - 1))
   do j = 1, nj - 1
      do i = 1, ni - 1
         vm(i, j) = abs(quad(blk(1), i, j)) ; pm(i, j) = blob(centre(blk(1), i, j))
      enddo
   enddo
   m0 = mass_a(g, phi, vg, pm, vm) ; out = 0._R8P ; defect = 0._R8P ; vmax = 0._R8P ; vmin = 0._R8P
   do step = 1, nsteps
      r = 0._R8P ; rm = 0._R8P
      ! major: across the cell edges (each linked pair once), the shared fluid length of the edge
      do p = 1, np
         do k = links%start(p), links%start(p + 1) - 1
            q = links%piece(k)
            if (q <= p) cycle
            call edge(pt%cut%piece_cell(:, p), links%face(k), x, mid, nrm, len)
            f = dot_product(vel(mid, rc), nrm) * links%area(k)
            fl = f * merge(phi(p), phi(q), f > 0._R8P)
            r(p) = r(p) + fl ; r(q) = r(q) - fl
         enddo
      enddo
      ! major: the box boundary (inflow 0, outflow upwind)
      do p = 1, np
         i = pt%cut%piece_cell(1, p) ; j = pt%cut%piece_cell(2, p)
         do t = 1, 4
            if (.not. ((t == 1 .and. i == 1) .or. (t == 2 .and. i == nx) .or. (t == 3 .and. j == 1) .or. &
                       (t == 4 .and. j == nx))) cycle
            call edge([i, j], t, x, mid, nrm, len)
            f = dot_product(vel(mid, rc), nrm) * pt%cut%piece_aperture(t, p) * len
            if (f <= 0._R8P) cycle
            r(p) = r(p) + f * phi(p) ; out = out + dt * f * phi(p)
         enddo
      enddo
      ! the smallest faces: major piece <-> minor cell (sf_area points into the minor cell)
      do k = 1, pt%ns
         p = pt%sf_piece(k) ; i = pt%sf_minor(2, k) ; j = pt%sf_minor(3, k)
         f = dot_product(vel(pt%sf_centroid(:, k), rc), pt%sf_area(:, k))
         fl = f * merge(phi(p), pm(i, j), f > 0._R8P)
         r(p) = r(p) + fl ; rm(i, j) = rm(i, j) - fl
      enddo
      ! minor: its interior faces (periodic in i); the inner face is a wall, the outer one is the smallest faces
      do j = 1, nj - 1
         do i = 1, ni - 1
            ia = i ; ib = mod(i, ni - 1_I4P) + 1_I4P ! across the vertex line i + 1
            call minor_face(blk(1), [i + 1_I4P, j], [i + 1_I4P, j + 1_I4P], centre(blk(1), ia, j), &
                            centre(blk(1), ib, j), mid, nrm)
            f = dot_product(vel(mid, rc), nrm)
            fl = f * merge(pm(ia, j), pm(ib, j), f > 0._R8P)
            rm(ia, j) = rm(ia, j) + fl ; rm(ib, j) = rm(ib, j) - fl
            if (j == nj - 1) cycle
            call minor_face(blk(1), [i, j + 1_I4P], [i + 1_I4P, j + 1_I4P], centre(blk(1), i, j), &
                            centre(blk(1), i, j + 1_I4P), mid, nrm)
            f = dot_product(vel(mid, rc), nrm)
            fl = f * merge(pm(i, j), pm(i, j + 1), f > 0._R8P)
            rm(i, j) = rm(i, j) + fl ; rm(i, j + 1) = rm(i, j + 1) - fl
         enddo
      enddo
      ! update: the groups by their summed residuals, the minor cells by theirs
      mg = 0._R8P
      do p = 1, np
         if (g(p) > 0) mg(g(p)) = mg(g(p)) + r(p)
      enddo
      do p = 1, np
         if (g(p) == p) mg(p) = phi(p) * vg(p) - dt * mg(p)
      enddo
      do p = 1, np
         if (g(p) > 0) phi(p) = mg(g(p)) / vg(g(p))
      enddo
      pm = pm - dt * rm / vm
      mtot = mass_a(g, phi, vg, pm, vm)
      defect = max(defect, abs(mtot + out - m0) / m0)
      vmax = max(vmax, maxval(phi), maxval(pm)) ; vmin = min(vmin, minval(phi), minval(pm))
   enddo
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.5,A,ES12.5)', '    A: pieces ', np, ', smallest faces ', pt%ns, &
         ', small merged ', nsm, ', degenerate folded ', ndg, ' (alone ', nal, '), steps ', nsteps, '; mass ', mtot, &
         ', left the box ', out
   endsubroutine patched

   real(R8P) function mass_a(g, phi, vg, pm, vm) result(tot)
   !< Total mass of A (Kahan): the groups and the minor cells.
   integer(I4P), intent(in) :: g(:)                 !< Groups.
   real(R8P),    intent(in) :: phi(:), vg(:)        !< Values, group volumes.
   real(R8P),    intent(in) :: pm(:,:), vm(:,:)     !< Minor values, volumes.
   real(R8P)                :: comp                 !< Compensation.
   integer(I4P)             :: p, i, j              !< Counters.

   tot = 0._R8P ; comp = 0._R8P
   do p = 1, size(g, kind=I4P)
      if (g(p) == p) call kahan(tot, comp, phi(p) * vg(p))
   enddo
   do j = 1, size(pm, 2, kind=I4P)
      do i = 1, size(pm, 1, kind=I4P)
         call kahan(tot, comp, pm(i, j) * vm(i, j))
      enddo
   enddo
   endfunction mass_a

   subroutine chimera(blk, x, rc, dt, nsteps, defect)
   !< B: the same assembly coupled by interpolation; defect as in A on the composite mass.
   type(curvilinear_block_2d_t), intent(in)  :: blk    !< Minor block.
   real(R8P),                    intent(in)  :: x(:)   !< Major grid edges (both directions).
   real(R8P),                    intent(in)  :: rc     !< Cylinder radius of the flow.
   real(R8P),                    intent(in)  :: dt     !< Time step.
   integer(I4P),                 intent(in)  :: nsteps !< Steps.
   real(R8P),                    intent(out) :: defect !< Mass defect.
   type(curvilinear_block_2d_t) :: dual                !< Block of the minor cell centres.
   integer(I4P), allocatable :: kd(:,:), cell(:,:), lst(:), rcv(:,:) !< Major kinds, donor cells, status, receptors.
   real(R8P),    allocatable :: phi(:,:), r(:,:), pm(:,:), rm(:,:), vm(:,:), xi(:,:), qx(:), qy(:) !< Work.
   real(R8P)    :: h, out, m0, mtot, f, fl, mid(2), nrm(2), c(2), rr, s(2), w(2), len !< Work.
   integer(I4P) :: nx, ni, nj, i, j, k, step, nr, t, ib, i0, j0, nbad !< Work.

   nx = size(x, kind=I4P) - 1 ; h = x(2) - x(1) ; ni = blk%ni ; nj = blk%nj
   allocate(kd(nx, nx), phi(nx, nx), r(nx, nx), pm(ni - 1, nj - 1), rm(ni - 1, nj - 1), vm(ni - 1, nj - 1))
   ! major kinds: 0 hole, 1 receptor, 2 field
   nr = 0
   do j = 1, nx
      do i = 1, nx
         c = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (x(j) + x(j + 1))] ; rr = norm2(c)
         kd(i, j) = 2 ; if (rr < 1.6_R8P) kd(i, j) = 1 ; if (rr < 1.4_R8P) kd(i, j) = 0
         if (kd(i, j) == 1) nr = nr + 1
         phi(i, j) = blob(c)
      enddo
   enddo
   nbad = 0 ! a field cell next to a hole: the fringe is too thin
   do j = 1, nx
      do i = 1, nx
         if (kd(i, j) /= 2) cycle
         if (i > 1) then
            if (kd(i - 1, j) == 0) nbad = nbad + 1
         endif
         if (i < nx) then
            if (kd(i + 1, j) == 0) nbad = nbad + 1
         endif
         if (j > 1) then
            if (kd(i, j - 1) == 0) nbad = nbad + 1
         endif
         if (j < nx) then
            if (kd(i, j + 1) == 0) nbad = nbad + 1
         endif
      enddo
   enddo
   do j = 1, nj - 1
      do i = 1, ni - 1
         vm(i, j) = abs(quad(blk, i, j)) ; pm(i, j) = blob(centre(blk, i, j))
      enddo
   enddo
   ! FOSSIL donors of the major receptors in the block of the minor cell centres
   call cell_centre_block_2d(blk, dual)
   call dual%build_search
   allocate(qx(nr), qy(nr), rcv(2, nr)) ; k = 0
   do j = 1, nx
      do i = 1, nx
         if (kd(i, j) /= 1) cycle
         k = k + 1 ; rcv(:, k) = [i, j] ; qx(k) = 0.5_R8P * (x(i) + x(i + 1)) ; qy(k) = 0.5_R8P * (x(j) + x(j + 1))
      enddo
   enddo
   call dual%locate_many(qx, qy, cell, xi, lst)
   if (any(lst == LOC_OUTSIDE)) error stop 'chimera: a major receptor outside the minor block'
   m0 = mass_b(x, phi, pm, vm) ; out = 0._R8P ; defect = 0._R8P
   do step = 1, nsteps
      ! interpolation: major receptors from the minor cell centres, the minor outer row from the major cell centres
      do k = 1, nr
         i0 = cell(1, k) ; j0 = cell(2, k)
         phi(rcv(1, k), rcv(2, k)) = (1 - xi(1, k)) * (1 - xi(2, k)) * dual_value(pm, i0, j0) + &
                                     xi(1, k) * (1 - xi(2, k)) * dual_value(pm, i0 + 1_I4P, j0) + &
                                     xi(1, k) * xi(2, k) * dual_value(pm, i0 + 1_I4P, j0 + 1_I4P) + &
                                     (1 - xi(1, k)) * xi(2, k) * dual_value(pm, i0, j0 + 1_I4P)
      enddo
      do i = 1, ni - 1
         c = centre(blk, i, nj - 1_I4P)
         s = (c - (x(1) + 0.5_R8P * h)) / h ; i0 = int(s(1)) + 1 ; j0 = int(s(2)) + 1 ; w = s - [i0 - 1, j0 - 1]
         if (any(kd(i0:i0 + 1, j0:j0 + 1) /= 2)) error stop 'chimera: a minor receptor with a major non-field donor'
         pm(i, nj - 1) = (1 - w(1)) * (1 - w(2)) * phi(i0, j0) + w(1) * (1 - w(2)) * phi(i0 + 1, j0) + &
                         w(1) * w(2) * phi(i0 + 1, j0 + 1) + (1 - w(1)) * w(2) * phi(i0, j0 + 1)
      enddo
      ! major field cells: the faces with a field cell on one side at least
      r = 0._R8P
      do j = 1, nx
         do i = 1, nx
            do t = 2, 4, 2
               ib = merge(i + 1, i, t == 2)
               if (t == 4) ib = j + 1
               if (ib > nx) cycle
               if (t == 2) then
                  if (kd(i, j) /= 2 .and. kd(i + 1, j) /= 2) cycle
                  call edge([i, j], t, x, mid, nrm, len)
                  f = dot_product(vel(mid, rc), nrm) * len
                  fl = f * merge(phi(i, j), phi(i + 1, j), f > 0._R8P)
                  r(i, j) = r(i, j) + fl ; r(i + 1, j) = r(i + 1, j) - fl
               else
                  if (kd(i, j) /= 2 .and. kd(i, j + 1) /= 2) cycle
                  call edge([i, j], t, x, mid, nrm, len)
                  f = dot_product(vel(mid, rc), nrm) * len
                  fl = f * merge(phi(i, j), phi(i, j + 1), f > 0._R8P)
                  r(i, j) = r(i, j) + fl ; r(i, j + 1) = r(i, j + 1) - fl
               endif
            enddo
            do t = 1, 4 ! the box boundary
               if (.not. ((t == 1 .and. i == 1) .or. (t == 2 .and. i == nx) .or. (t == 3 .and. j == 1) .or. &
                          (t == 4 .and. j == nx))) cycle
               call edge([i, j], t, x, mid, nrm, len)
               f = dot_product(vel(mid, rc), nrm) * len
               if (f <= 0._R8P) cycle
               r(i, j) = r(i, j) + f * phi(i, j) ; out = out + dt * f * phi(i, j)
            enddo
         enddo
      enddo
      ! minor field cells (j < nj - 1)
      rm = 0._R8P
      do j = 1, nj - 1
         do i = 1, ni - 1
            if (j < nj - 1) then
               ib = mod(i, ni - 1_I4P) + 1_I4P
               call minor_face(blk, [i + 1_I4P, j], [i + 1_I4P, j + 1_I4P], centre(blk, i, j), centre(blk, ib, j), mid, &
                               nrm)
               f = dot_product(vel(mid, rc), nrm)
               fl = f * merge(pm(i, j), pm(ib, j), f > 0._R8P)
               rm(i, j) = rm(i, j) + fl ; rm(ib, j) = rm(ib, j) - fl
               call minor_face(blk, [i, j + 1_I4P], [i + 1_I4P, j + 1_I4P], centre(blk, i, j), &
                               centre(blk, i, j + 1_I4P), mid, nrm)
               f = dot_product(vel(mid, rc), nrm)
               fl = f * merge(pm(i, j), pm(i, j + 1), f > 0._R8P)
               rm(i, j) = rm(i, j) + fl ; rm(i, j + 1) = rm(i, j + 1) - fl
            endif
         enddo
      enddo
      where (kd == 2) phi = phi - dt * r / h ** 2
      pm(:, 1:nj - 2) = pm(:, 1:nj - 2) - dt * rm(:, 1:nj - 2) / vm(:, 1:nj - 2)
      mtot = mass_b(x, phi, pm, vm)
      defect = max(defect, abs(mtot + out - m0) / m0)
   enddo
   print '(A,I0,A,I0,A,ES12.5,A,ES12.5)', '    B: major receptors ', nr, ', field cells next to holes ', nbad, &
         '; mass ', mtot, ', left the box ', out
   if (nbad > 0) error stop 'chimera: fringe too thin'
   endsubroutine chimera

   real(R8P) function mass_b(x, phi, pm, vm) result(tot)
   !< Composite mass of B (Kahan): the major cells with centre at r >= 2, all minor cells.
   real(R8P), intent(in) :: x(:)             !< Major grid.
   real(R8P), intent(in) :: phi(:,:)         !< Major values.
   real(R8P), intent(in) :: pm(:,:), vm(:,:) !< Minor values, volumes.
   real(R8P)             :: comp, c(2)       !< Work.
   integer(I4P)          :: i, j             !< Counters.

   tot = 0._R8P ; comp = 0._R8P
   do j = 1, size(phi, 2, kind=I4P)
      do i = 1, size(phi, 1, kind=I4P)
         c = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (x(j) + x(j + 1))]
         if (norm2(c) >= 2._R8P) call kahan(tot, comp, phi(i, j) * (x(2) - x(1)) ** 2)
      enddo
   enddo
   do j = 1, size(pm, 2, kind=I4P)
      do i = 1, size(pm, 1, kind=I4P)
         call kahan(tot, comp, pm(i, j) * vm(i, j))
      enddo
   enddo
   endfunction mass_b

   pure real(R8P) function dual_value(pm, i, j) result(val)
   !< The value at vertex (i, j) of the block of the cell centres: cell (i - 1, j - 1) inside; on its boundary layer the
   !< mean across the periodic seam (i = 1 or ni + 1), the adjacent cell at the wall and at the outer face.
   real(R8P),    intent(in) :: pm(:,:) !< Minor values [ni - 1, nj - 1].
   integer(I4P), intent(in) :: i, j    !< Dual vertex.
   integer(I4P)             :: jc, nc  !< Cell row, cells around.

   nc = size(pm, 1, kind=I4P)
   jc = min(max(j - 1_I4P, 1_I4P), size(pm, 2, kind=I4P))
   if (i == 1 .or. i == nc + 2) then
      val = 0.5_R8P * (pm(1, jc) + pm(nc, jc))
   else
      val = pm(i - 1, jc)
   endif
   endfunction dual_value

   pure subroutine kahan(tot, comp, a)
   !< tot += a, compensated.
   real(R8P), intent(inout) :: tot, comp !< Sum, compensation.
   real(R8P), intent(in)    :: a         !< Term.
   real(R8P)                :: y, t      !< Work.

   y = a - comp ; t = tot + y ; comp = (t - tot) - y ; tot = t
   endsubroutine kahan

   pure function vel(p, rc) result(u)
   !< The potential flow past the cylinder of radius rc at the origin (free stream 1 at angle ANG); uniform if rc = 0.
   real(R8P), intent(in) :: p(2) !< Point.
   real(R8P), intent(in) :: rc   !< Radius.
   real(R8P)             :: u(2) !< Velocity.
   complex(R8P)          :: z, dw !< Work.

   if (rc == 0._R8P) then
      u = [cos(ANG), sin(ANG)] ; return
   endif
   z = cmplx(p(1), p(2), R8P)
   dw = exp(cmplx(0._R8P, -ANG, R8P)) - rc ** 2 * exp(cmplx(0._R8P, ANG, R8P)) / z ** 2
   u = [real(dw, R8P), -aimag(dw)]
   endfunction vel

   pure real(R8P) function blob(p)
   !< The initial scalar: a Gaussian upstream of the minor block.
   real(R8P), intent(in) :: p(2) !< Point.

   blob = exp(-((p(1) + 2.3_R8P) ** 2 + (p(2) - 0.4_R8P) ** 2) / 0.25_R8P ** 2)
   endfunction blob

   pure subroutine edge(c, t, x, mid, nrm, len)
   !< Face t (1 -u, 2 +u, 3 -v, 4 +v) of the major cell c: midpoint, outward unit normal, length.
   integer(I4P), intent(in)  :: c(2), t !< Cell, face.
   real(R8P),    intent(in)  :: x(:)    !< Grid edges.
   real(R8P),    intent(out) :: mid(2), nrm(2), len !< Face.

   select case(t)
   case(1)
      mid = [x(c(1)), 0.5_R8P * (x(c(2)) + x(c(2) + 1))] ; nrm = [-1._R8P, 0._R8P] ; len = x(c(2) + 1) - x(c(2))
   case(2)
      mid = [x(c(1) + 1), 0.5_R8P * (x(c(2)) + x(c(2) + 1))] ; nrm = [1._R8P, 0._R8P] ; len = x(c(2) + 1) - x(c(2))
   case(3)
      mid = [0.5_R8P * (x(c(1)) + x(c(1) + 1)), x(c(2))] ; nrm = [0._R8P, -1._R8P] ; len = x(c(1) + 1) - x(c(1))
   case default
      mid = [0.5_R8P * (x(c(1)) + x(c(1) + 1)), x(c(2) + 1)] ; nrm = [0._R8P, 1._R8P] ; len = x(c(1) + 1) - x(c(1))
   endselect
   endsubroutine edge

   pure subroutine minor_face(blk, v0, v1, ca, cb, mid, nrm)
   !< The face of the minor block from vertex v0 to v1: midpoint and length x normal pointing from ca to cb.
   type(curvilinear_block_2d_t), intent(in)  :: blk       !< Block.
   integer(I4P),                 intent(in)  :: v0(2), v1(2) !< Vertices.
   real(R8P),                    intent(in)  :: ca(2), cb(2) !< Centres of the cells on either side.
   real(R8P),                    intent(out) :: mid(2), nrm(2) !< Face.
   real(R8P)                                 :: a(2), b(2) !< Ends.

   a = [blk%u(v0(1), v0(2)), blk%v(v0(1), v0(2))] ; b = [blk%u(v1(1), v1(2)), blk%v(v1(1), v1(2))]
   mid = 0.5_R8P * (a + b) ; nrm = [b(2) - a(2), -(b(1) - a(1))]
   if (dot_product(nrm, cb - ca) < 0._R8P) nrm = -nrm
   endsubroutine minor_face

   pure function centre(blk, i, j) result(c)
   !< Vertex mean of cell (i, j).
   type(curvilinear_block_2d_t), intent(in) :: blk  !< Block.
   integer(I4P),                 intent(in) :: i, j !< Cell.
   real(R8P)                                :: c(2) !< Centre.

   c = 0.25_R8P * [blk%u(i, j) + blk%u(i + 1, j) + blk%u(i + 1, j + 1) + blk%u(i, j + 1), &
                   blk%v(i, j) + blk%v(i + 1, j) + blk%v(i + 1, j + 1) + blk%v(i, j + 1)]
   endfunction centre

   pure real(R8P) function quad(blk, i, j) result(a)
   !< Signed area of cell (i, j).
   type(curvilinear_block_2d_t), intent(in) :: blk  !< Block.
   integer(I4P),                 intent(in) :: i, j !< Cell.

   a = 0.5_R8P * ((blk%u(i + 1, j + 1) - blk%u(i, j)) * (blk%v(i, j + 1) - blk%v(i + 1, j)) - &
                  (blk%u(i, j + 1) - blk%u(i + 1, j)) * (blk%v(i + 1, j + 1) - blk%v(i, j)))
   endfunction quad

   subroutine ogrid(m, blk)
   !< The cylinder O-grid: r in [1, 2], 64 m x 8 m cells, its outer face (j = nj) named 'overlap'.
   integer(I4P),                 intent(in)    :: m   !< Refinement.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   type(topo_record_t)    :: rec   !< Record.
   real(R8P), allocatable :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)           :: ni, nj, i, j !< Counters.
   real(R8P)              :: r     !< Radius.

   ni = 64 * m + 1 ; nj = 8 * m + 1
   allocate(u(ni, nj), v(ni, nj))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = r * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni, :) = u(1, :) ; v(ni, :) = v(1, :)
   call blk%destroy
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine ogrid

   subroutine ring_block(pu, pv, nl, blk)
   !< A block of nl layers between the polygon (outer, named 'overlap') and its copy shrunk by half towards its vertex
   !< mean (inner).
   real(R8P),                    intent(in)    :: pu(:), pv(:) !< Polygon (counter-clockwise).
   integer(I4P),                 intent(in)    :: nl           !< Layers.
   type(curvilinear_block_2d_t), intent(inout) :: blk          !< Block.
   type(topo_record_t)    :: rec            !< Record.
   real(R8P), allocatable :: u(:,:), v(:,:) !< Vertices.
   real(R8P)              :: c(2), s        !< Vertex mean, scale.
   integer(I4P)           :: n, j           !< Counters.

   n = size(pu, kind=I4P)
   allocate(u(n + 1, nl + 1), v(n + 1, nl + 1))
   c = [sum(pu), sum(pv)] / n
   u(1:n, nl + 1) = pu ; v(1:n, nl + 1) = pv ; u(n + 1, nl + 1) = pu(1) ; v(n + 1, nl + 1) = pv(1)
   do j = 1, nl
      s = 0.5_R8P + 0.5_R8P * (j - 1) / nl
      u(:, j) = c(1) + s * (u(:, nl + 1) - c(1)) ; v(:, j) = c(2) + s * (v(:, nl + 1) - c(2))
   enddo
   call blk%destroy
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, nl + 1, 1] ; rec%hi = [n + 1, nl + 1, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine ring_block

   pure function grid_1d(a, b, n) result(x)
   !< n + 1 equally spaced edges on [a, b].
   real(R8P),    intent(in) :: a, b     !< Ends.
   integer(I4P), intent(in) :: n        !< Cells.
   real(R8P)                :: x(n + 1) !< Edges.
   integer(I4P)             :: i        !< Counter.

   do i = 1, n + 1
      x(i) = a + (b - a) * (i - 1) / n
   enddo
   x(n + 1) = b
   endfunction grid_1d
endprogram fossil_val_overset_conservative
