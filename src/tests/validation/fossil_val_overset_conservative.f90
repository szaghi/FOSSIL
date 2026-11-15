!< FOSSIL validation: a passive scalar through a conservative overset interface (issue #106 steps 1-4, solver gate).
program fossil_val_overset_conservative
!< FOSSIL validation: a passive scalar through a conservative overset interface (issue #106 steps 1-4, solver gate).
!<
!< A minor ring block (a cylinder O-grid, r in [1, 2]; or a square ring on grid lines, with degenerate pieces) in a
!< Cartesian major grid on [-3, 3]^2. A scalar blob is advected by the potential flow around the cylinder (uniform flow
!< for the square) with first-order upwind finite volumes and forward Euler, three ways:
!< A. Wang's patched interface (`overset_patch_2d`): the major grid is cut by the minor boundary (holes, cut pieces with
!<    their apertures W_f), the minor cells keep their faces, and the two grids exchange the upwind fluxes through the
!<    smallest faces, each flux added to one side and subtracted from the other. Small major pieces (area below half a
!<    cell, `merge_targets`) and degenerate ones (`merge_into`) are merged into their targets: a group shares one value
!<    and the sum of its members' fluxes. The total mass plus the mass that left through the box changes by round-off.
!< I. The same assembly coupled by interpolation (Chimera): the major cells whose centre lies at r < 1.4 are holes, at
!<    1.4 <= r < 1.6 receptors (bilinear from the minor cell centres, located by FOSSIL in the block of the cell centres,
!<    `cell_centre_block_2d` + `locate_many`), the minor's outer row of cells receives from the major cell centres
!<    (bilinear); each grid advances its own field cells. The composite mass (the major cells with centre at r >= 2, all
!<    minor cells) plus the mass that left the box drifts by the interpolation.
!< B. The assembly of I with Picard's intersection-weighted transfer (`overlap_weights_2d`, issue #106 step 2) in place
!<    of the interpolation: a receptor takes q_c = sum w (q_c' + grad q_c' . (x_cc' - x_c')) over the cells c' it
!<    intersects, w = |V_c n V_c'| / |V_c| and x_cc' the centroid of V_c n V_c' (the minor outer row from the major
!<    cells; the major receptors from the minor cells, by the same intersections read the other way). Gradients: central
!<    differences on the major grid (one-sided next to a hole), Green-Gauss on the minor cells. It transfers states,
!<    not fluxes: its composite mass (as I's) drifts too.
!< The inner minor boundary is a wall (no flux). Three measures (issue #106 step 4) separate what each number says:
!< - the composite defect: max over the steps of |M + out - M0| / M0, M the composite mass (A: its groups and the minor
!<   cells, an exact partition of the box; I and B: the major cells with centre at r >= 2, counted whole or not at all,
!<   and every minor cell, a bookkeeping that miscounts the cells straddling r = 2);
!< - I and B on the exact partition: the major cells weighted by their fluid fraction outside the minor polygon (the
!<   cut of `overset_patch_2d`), the minor cells: the coupling alone, without the bookkeeping;
!< - one transfer of a smooth field (`transfer_once`, no time step): the mass the transfer puts in the receptors against
!<   the donor-side mass over the same region (the donor's linear reconstruction integrated over the intersections),
!<   and against the exact integral: B conserves the donor's mass per transfer to round-off, I does not; both are
!<   second order. B's composite defect is then its coupling (states overwritten, not fluxes exchanged), not its
!<   weights.
!< Exit status 0 if A conserves to 1e-12 relative and stays in [0, 1] (upwind within its CFL limit), the defects of I
!< and B exceed 1e-6 on both masses, B's transfer conserves the donor's mass to 1e-14 and I's misses it by more than
!< 1e-6, and both transfers are within 1e-2 of the exact integral.
!<
!< 3-D (issue #106 step 3): a sphere shell (r in [1, 2], poles collapsed) in a Cartesian [-3, 3]^3, the scalar advected
!< by the potential flow past the unit sphere (inner face a wall), time step from the outflows at CFL 0.9:
!< A. `overset_patch_3d`: first the geometry (closure per piece to 1e-12 h^2, face_covered = face_vector to 1e-13, the
!<    pieces' volume = box - polyhedron to 1e-13, sum sf_area = - piece_wet_normal), then the scalar as in 2-D (links
!<    of the cut, merge targets, smallest faces): mass to 1e-12;
!< I. interpolation: holes at centre radius < 1.2, receptors below 1.6 (trilinear in the shell's index space of cell
!<    centres), the shell's outer layer from the major cell centres (trilinear): its defect must exceed 1e-6.
!< B. `overlap_weights_3d` with every shell cell a receptor: weights sum to 1 (2 eps), a constant field exact, a linear
!<    one with Picard's gradient term at the receptor polyhedron's centroid (1e-13), and per donor cell the sum of its
!<    intersections over the receptors = |V_c' n polyhedral shell| from `overset_patch_3d` of the shell's inner and outer
!<    faces (1e-13 h^3): the receptors' splits agree with each other and with the boundary's.
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d, &
                   OVERSET_PATCH_OK, cut_links_t, cut_links_2d, merge_t, merge_targets, MERGE_OK, MERGE_CHAIN, &
                   cell_centre_block_2d, LOC_OUTSIDE, overlap_weights_t, overlap_weights_2d, OVERLAP_WEIGHTS_OK, &
                   curvilinear_block_t, overset_patch_3d_t, overset_patch_3d, overlap_weights_3d
use penf,   only : I8P
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P) !< Pi.
real(R8P), parameter :: ANG = 0.2_R8P      !< Free-stream angle.
real(R8P), parameter :: U3(3) = [cos(0.2_R8P) * cos(0.15_R8P), sin(0.2_R8P) * cos(0.15_R8P), sin(0.15_R8P)] !< 3-D free stream.
type :: minor_3d_t
   !< The cells and interior faces of the 3-D minor shell (the scalar's finite volumes).
   real(R8P),    allocatable :: vol(:,:,:)    !< Cell volumes.
   real(R8P),    allocatable :: cen(:,:,:,:)  !< Cell centres (vertex means) [3, ...].
   real(R8P),    allocatable :: val0(:,:,:)   !< Initial values.
   integer(I4P)              :: nf = 0_I4P    !< Interior faces.
   integer(I4P), allocatable :: fa(:,:), fb(:,:) !< Their cells (from, to) [3, nf].
   real(R8P),    allocatable :: fvec(:,:)     !< Face vector from fa to fb [3, nf].
   real(R8P),    allocatable :: fmid(:,:)     !< Vertex mean [3, nf].
endtype minor_3d_t
type(curvilinear_block_t)    :: blk3(1)    !< 3-D minor block.
type(curvilinear_block_2d_t) :: blk(1)     !< Minor block.
real(R8P), allocatable       :: x(:), pu(:), pv(:) !< Major grid, square.
real(R8P)                    :: da, db, di, mx, mn  !< Defects, extremes.
real(R8P)                    :: dbp, dip, cb, ci, ab, ai !< Defects on the exact partition; one transfer (B, I).
integer(I4P)                 :: m, k, n     !< Counters.
logical                      :: ok          !< Checks.

ok = .true.
do m = 1, 2
   call ogrid(m, blk(1))
   x = grid_1d(-3._R8P, 3._R8P, 60_I4P * m)
   call patched(blk, x, 1._R8P, 0.008_R8P / m, 400_I4P * m, da, mx, mn)
   call chimera(blk(1), x, 1._R8P, 0.008_R8P / m, 400_I4P * m, .false., di, dip)
   call chimera(blk(1), x, 1._R8P, 0.008_R8P / m, 400_I4P * m, .true., db, dbp)
   call transfer_once(blk(1), x, cb, ci, ab, ai)
   print '(A,I0,A,I0,A,ES10.2,A,ES10.2,A,ES10.2,A,2ES10.2)', '  cylinder, m = ', m, ' (', 60 * m, &
         '^2): mass defect A (patched) ', da, ', B (Picard weights) ', db, ', I (interpolated) ', di, '; A min, max ', &
         mn, mx
   print '(A,ES10.2,A,ES10.2,A)', '    on the exact partition (no bookkeeping): B ', dbp, ', I ', dip, &
         ' (A: as above, its pieces and minor cells are the partition)'
   print '(A,ES10.2,A,ES10.2,A,ES10.2,A,ES10.2)', '    one transfer of a smooth field: mass vs donor side B ', cb, &
         ', I ', ci, '; vs exact integral B ', ab, ', I ', ai
   ok = ok .and. da <= 1e-12_R8P .and. di > 1e-6_R8P .and. db > 1e-6_R8P .and. mx <= 1._R8P + 1e-12_R8P .and. &
        mn >= -1e-12_R8P .and. dip > 1e-6_R8P .and. dbp > 1e-6_R8P .and. cb <= 1e-14_R8P .and. ci > 1e-6_R8P .and. &
        ab <= 1e-2_R8P .and. ai <= 1e-2_R8P
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
! 3-D: a sphere shell (r in [1, 2], 32 x 16 x 8 cells) in [-3, 3]^3 of 40^3 cells, potential flow past the unit sphere
call shell(32_I4P, 8_I4P, blk3(1))
x = grid_1d(-3._R8P, 3._R8P, 40_I4P)
call patched_3d(blk3, x, da, mx, mn)
call chimera_3d(blk3(1), x, di)
print '(A,ES10.2,A,ES10.2,A,2ES10.2)', '  sphere shell (40^3): mass defect A (patched) ', da, ', I (interpolated) ', di, &
      '; A min, max ', mn, mx
ok = ok .and. da <= 1e-12_R8P .and. di > 1e-6_R8P .and. mx <= 1._R8P + 1e-12_R8P .and. mn >= -1e-12_R8P
! 3-D B: the whole shell (24 x 12 x 4 cells) as receptors over [-2.6, 2.6]^3 of 24^3 cells
call shell(24_I4P, 4_I4P, blk3(1))
x = grid_1d(-2.6_R8P, 2.6_R8P, 24_I4P)
call weights_3d(blk3(1), x, ok)
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

   subroutine chimera(blk, x, rc, dt, nsteps, picard, defect, defect_p)
   !< I (interpolation) or B (Picard's weights, `picard`): the same overlapping assembly; defect as in A on the composite
   !< mass (major cells with centre at r >= 2 and every minor cell: a bookkeeping that counts a cell straddling the
   !< interface whole or not at all), defect_p on an exact partition of the box (each major cell weighted by its fluid
   !< fraction outside the minor polygon, the cut of `overset_patch_2d`, and every minor cell): the coupling alone.
   type(curvilinear_block_2d_t), intent(in)  :: blk    !< Minor block.
   real(R8P),                    intent(in)  :: x(:)   !< Major grid edges (both directions).
   real(R8P),                    intent(in)  :: rc     !< Cylinder radius of the flow.
   real(R8P),                    intent(in)  :: dt     !< Time step.
   integer(I4P),                 intent(in)  :: nsteps !< Steps.
   logical,                      intent(in)  :: picard !< B: the intersection weights instead of the interpolation.
   real(R8P),                    intent(out) :: defect !< Mass defect (composite bookkeeping).
   real(R8P),                    intent(out) :: defect_p !< Mass defect (exact partition).
   type(overlap_weights_t)   :: ow                     !< Intersections of every minor cell with the major cells.
   type(overset_patch_2d_t)  :: pt                     !< The cut of the major grid by the minor polygon (alpha).
   integer(I4P), allocatable :: tstart(:), tcell(:)    !< Major receptor -> its intersections in ow (CSR by cell).
   real(R8P),    allocatable :: gx(:,:,:), gm(:,:,:), xm(:,:,:) !< Major, minor gradients; minor centroids.
   type(curvilinear_block_2d_t) :: dual                !< Block of the minor cell centres.
   integer(I4P), allocatable :: kd(:,:), cell(:,:), lst(:), rcv(:,:) !< Major kinds, donor cells, status, receptors.
   real(R8P),    allocatable :: phi(:,:), r(:,:), pm(:,:), rm(:,:), vm(:,:), xi(:,:), qx(:), qy(:) !< Work.
   real(R8P)    :: h, out, m0, mtot, f, fl, mid(2), nrm(2), c(2), rr, len, m0p, mtp !< Work.
   integer(I4P) :: nx, ni, nj, i, j, k, step, nr, t, ib, nbad, st !< Work.

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
   if (picard) call picard_setup(blk, x, kd, ow, tstart, tcell, xm)
   if (picard) allocate(gx(2, nx, nx), gm(2, ni - 1, nj - 1))
   call overset_patch_2d([blk], [1_I4P], x, x, pt, status=st)
   if (st /= OVERSET_PATCH_OK) error stop 'chimera: overset_patch_2d'
   m0 = mass_b(x, phi, pm, vm) ; m0p = mass_p(x, pt%cut%alpha, phi, pm, vm) ; out = 0._R8P ; defect = 0._R8P
   defect_p = 0._R8P
   do step = 1, nsteps
      if (picard) then
         call picard_transfer(blk, x, kd, ow, tstart, tcell, xm, gx, gm, phi, pm)
      else
         call interp_transfer(blk, x, kd, rcv, cell, xi, phi, pm)
      endif
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
      mtot = mass_b(x, phi, pm, vm) ; mtp = mass_p(x, pt%cut%alpha, phi, pm, vm)
      defect = max(defect, abs(mtot + out - m0) / m0) ; defect_p = max(defect_p, abs(mtp + out - m0p) / m0p)
   enddo
   print '(A,A,I0,A,I0,A,ES12.5,A,ES12.5,A,I0,A,ES10.2)', merge('    B', '    I', picard), ': major receptors ', nr, &
         ', field cells next to holes ', nbad, '; mass ', mtot, ', left the box ', out, &
         '; intersections ', merge(ow%nw, 0_I4P, picard), '; defect on the exact partition ', defect_p
   if (nbad > 0) error stop 'chimera: fringe too thin'
   endsubroutine chimera

   subroutine picard_setup(blk, x, kd, ow, tstart, tcell, xm)
   !< B: the intersections of every minor cell with the major cells (`overlap_weights_2d`), checked (the minor outer row
   !< takes field cells only, the major receptors are fully covered by minor field cells), read the other way for the
   !< major receptors (tcell: the intersections of major cell c in ow, tstart by linear cell index), and the minor cell
   !< centroids.
   type(curvilinear_block_2d_t), intent(in)  :: blk        !< Minor block.
   real(R8P),                    intent(in)  :: x(:)       !< Major grid.
   integer(I4P),                 intent(in)  :: kd(:,:)    !< Major kinds.
   type(overlap_weights_t),      intent(out) :: ow         !< Intersections.
   integer(I4P), allocatable,    intent(out) :: tstart(:), tcell(:) !< Major cell -> intersections.
   real(R8P),    allocatable,    intent(out) :: xm(:,:,:)  !< Minor centroids [2, ni - 1, nj - 1].
   integer(I4P), allocatable :: cells(:,:), pos(:)         !< Receptor cells, fill positions.
   real(R8P)    :: ar, cv                                  !< Work.
   integer(I4P) :: ni, nj, nx, i, j, k, r, c, st           !< Work.

   ni = blk%ni ; nj = blk%nj ; nx = size(kd, 1, kind=I4P)
   allocate(cells(2, (ni - 1) * (nj - 1)), xm(2, ni - 1, nj - 1))
   do j = 1, nj - 1
      do i = 1, ni - 1
         cells(:, i + (j - 1) * (ni - 1)) = [i, j]
         call polygon_centroid(blk, i, j, ar, xm(:, i, j))
      enddo
   enddo
   call overlap_weights_2d(blk, cells, x, x, ow, status=st)
   if (st /= OVERLAP_WEIGHTS_OK .or. any(ow%flag /= OVERLAP_WEIGHTS_OK)) error stop 'picard: overlap_weights_2d'
   if (any(abs(ow%covered - 1._R8P) > 0._R8P)) error stop 'picard: a minor cell leaves the major box'
   allocate(tstart(nx * nx + 1), pos(nx * nx)) ; tstart = 0
   do r = 1, ow%nr
      do k = ow%start(r), ow%start(r + 1) - 1
         c = ow%donor(1, k) + (ow%donor(2, k) - 1) * nx
         if (cells(2, r) == nj - 1 .and. ow%weight(k) > 0._R8P .and. kd(ow%donor(1, k), ow%donor(2, k)) /= 2) &
            error stop 'picard: a minor receptor with a major non-field donor'
         tstart(c + 1) = tstart(c + 1) + 1
      enddo
   enddo
   tstart(1) = 1
   do c = 1, nx * nx
      tstart(c + 1) = tstart(c + 1) + tstart(c)
   enddo
   allocate(tcell(tstart(nx * nx + 1) - 1)) ; pos = tstart(1:nx * nx)
   do r = 1, ow%nr
      do k = ow%start(r), ow%start(r + 1) - 1
         c = ow%donor(1, k) + (ow%donor(2, k) - 1) * nx
         tcell(pos(c)) = k ; pos(c) = pos(c) + 1
      enddo
   enddo
   do j = 1, nx ! a major receptor: covered by minor field cells
      do i = 1, nx
         if (kd(i, j) /= 1) cycle
         c = i + (j - 1) * nx ; cv = 0._R8P
         do k = tstart(c), tstart(c + 1) - 1
            cv = cv + ow%volume(tcell(k))
         enddo
         if (abs(cv - (x(i + 1) - x(i)) * (x(j + 1) - x(j))) > 1e-13_R8P * (x(2) - x(1)) ** 2) &
            error stop 'picard: a major receptor not covered by the minor block'
      enddo
   enddo
   endsubroutine picard_setup

   subroutine picard_transfer(blk, x, kd, ow, tstart, tcell, xm, gx, gm, phi, pm)
   !< B: the major receptors from the minor cells, then the minor outer row from the major cells (the order of I),
   !< q_c = sum w (q_c' + grad q_c' . (x_cc' - x_c')).
   type(curvilinear_block_2d_t), intent(in)    :: blk       !< Minor block.
   real(R8P),                    intent(in)    :: x(:)      !< Major grid.
   integer(I4P),                 intent(in)    :: kd(:,:)   !< Major kinds.
   type(overlap_weights_t),      intent(in)    :: ow        !< Intersections.
   integer(I4P),                 intent(in)    :: tstart(:), tcell(:) !< Major cell -> intersections.
   real(R8P),                    intent(in)    :: xm(:,:,:) !< Minor centroids.
   real(R8P),                    intent(inout) :: gx(:,:,:), gm(:,:,:) !< Gradients (work).
   real(R8P),                    intent(inout) :: phi(:,:), pm(:,:) !< Major, minor values.
   real(R8P)    :: q, v, xc(2)                    !< Work.
   integer(I4P) :: nx, ni, nj, i, j, k, l, r, c, a, b !< Work.

   nx = size(phi, 1, kind=I4P) ; ni = blk%ni ; nj = blk%nj
   call minor_gradients(blk, xm, pm, gm)
   do j = 1, nx
      do i = 1, nx
         if (kd(i, j) /= 1) cycle
         c = i + (j - 1) * nx ; q = 0._R8P ; v = 0._R8P
         do l = tstart(c), tstart(c + 1) - 1
            k = tcell(l)
            call receptor_of(ow, k, r) ; a = mod(r - 1, ni - 1) + 1 ; b = (r - 1) / (ni - 1) + 1
            q = q + ow%volume(k) * (pm(a, b) + dot_product(gm(:, a, b), ow%centroid(:, k) - xm(:, a, b)))
            v = v + ow%volume(k)
         enddo
         phi(i, j) = q / v
      enddo
   enddo
   call major_gradients(x, kd, phi, gx)
   b = nj - 1
   do a = 1, ni - 1
      r = a + (b - 1) * (ni - 1) ; q = 0._R8P
      do k = ow%start(r), ow%start(r + 1) - 1
         i = ow%donor(1, k) ; j = ow%donor(2, k)
         xc = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (x(j) + x(j + 1))]
         q = q + ow%weight(k) * (phi(i, j) + dot_product(gx(:, i, j), ow%centroid(:, k) - xc))
      enddo
      pm(a, b) = q
   enddo
   endsubroutine picard_transfer

   pure subroutine receptor_of(ow, k, r)
   !< The receptor of intersection k (bisection on ow%start).
   type(overlap_weights_t), intent(in)  :: ow !< Intersections.
   integer(I4P),            intent(in)  :: k  !< Intersection.
   integer(I4P),            intent(out) :: r  !< Receptor.
   integer(I4P)                         :: lo, hi, mid !< Bisection.

   lo = 1 ; hi = ow%nr
   do while (lo < hi) ! start(lo) <= k < start(hi + 1)
      mid = (lo + hi + 1) / 2
      if (ow%start(mid) <= k) then
         lo = mid
      else
         hi = mid - 1
      endif
   enddo
   r = lo
   endsubroutine receptor_of

   pure subroutine major_gradients(x, kd, phi, gx)
   !< Central differences on the major cells that are not holes (one-sided next to a hole or the box, 0 if both sides
   !< are missing).
   real(R8P),    intent(in)  :: x(:)       !< Major grid.
   integer(I4P), intent(in)  :: kd(:,:)    !< Major kinds.
   real(R8P),    intent(in)  :: phi(:,:)   !< Values.
   real(R8P),    intent(out) :: gx(:,:,:)  !< Gradients [2, nx, nx].
   real(R8P)    :: h                       !< Cell size.
   integer(I4P) :: nx, i, j                !< Work.
   logical      :: lo, hi                  !< Neighbours available.

   nx = size(phi, 1, kind=I4P) ; h = x(2) - x(1) ; gx = 0._R8P
   do j = 1, nx
      do i = 1, nx
         if (kd(i, j) == 0) cycle
         lo = i > 1 ; if (lo) lo = kd(i - 1, j) /= 0
         hi = i < nx ; if (hi) hi = kd(i + 1, j) /= 0
         if (lo .and. hi) then
            gx(1, i, j) = (phi(i + 1, j) - phi(i - 1, j)) / (2 * h)
         elseif (hi) then
            gx(1, i, j) = (phi(i + 1, j) - phi(i, j)) / h
         elseif (lo) then
            gx(1, i, j) = (phi(i, j) - phi(i - 1, j)) / h
         endif
         lo = j > 1 ; if (lo) lo = kd(i, j - 1) /= 0
         hi = j < nx ; if (hi) hi = kd(i, j + 1) /= 0
         if (lo .and. hi) then
            gx(2, i, j) = (phi(i, j + 1) - phi(i, j - 1)) / (2 * h)
         elseif (hi) then
            gx(2, i, j) = (phi(i, j + 1) - phi(i, j)) / h
         elseif (lo) then
            gx(2, i, j) = (phi(i, j) - phi(i, j - 1)) / h
         endif
      enddo
   enddo
   endsubroutine major_gradients

   pure subroutine minor_gradients(blk, xm, pm, gm)
   !< Green-Gauss on the minor cells: face values the mean of the two cells (periodic in i; the cell's own value on the
   !< wall and the outer face).
   type(curvilinear_block_2d_t), intent(in)  :: blk       !< Minor block.
   real(R8P),                    intent(in)  :: xm(:,:,:) !< Centroids.
   real(R8P),                    intent(in)  :: pm(:,:)   !< Values.
   real(R8P),                    intent(out) :: gm(:,:,:) !< Gradients [2, ni - 1, nj - 1].
   real(R8P)    :: mid(2), fm(2), nrm(2), ar, g(2) !< Work.
   integer(I4P) :: ni, nj, i, j, t, ia, ja   !< Work.
   integer(I4P) :: v0(2), v1(2)              !< Face vertices.

   ni = blk%ni ; nj = blk%nj
   do j = 1, nj - 1
      do i = 1, ni - 1
         g = 0._R8P
         do t = 1, 4
            select case(t)
            case(1)
               v0 = [i, j] ; v1 = [i, j + 1] ; ia = mod(i - 2 + ni - 1, ni - 1) + 1 ; ja = j
            case(2)
               v0 = [i + 1, j] ; v1 = [i + 1, j + 1] ; ia = mod(i, ni - 1) + 1 ; ja = j
            case(3)
               v0 = [i, j] ; v1 = [i + 1, j] ; ia = i ; ja = max(j - 1, 1)
            case default
               v0 = [i, j + 1] ; v1 = [i + 1, j + 1] ; ia = i ; ja = min(j + 1, nj - 1)
            endselect
            fm = 0.5_R8P * [blk%u(v0(1), v0(2)) + blk%u(v1(1), v1(2)), blk%v(v0(1), v0(2)) + blk%v(v1(1), v1(2))]
            call minor_face(blk, v0, v1, xm(:, i, j), fm, mid, nrm) ! outward: from the centroid to the midpoint
            g = g + 0.5_R8P * (pm(i, j) + pm(ia, ja)) * nrm
         enddo
         ar = abs(quad(blk, i, j))
         gm(:, i, j) = g / ar
      enddo
   enddo
   endsubroutine minor_gradients

   pure subroutine polygon_centroid(blk, i, j, ar, c)
   !< Area and centroid of minor cell (i, j) (about its first vertex).
   type(curvilinear_block_2d_t), intent(in)  :: blk   !< Block.
   integer(I4P),                 intent(in)  :: i, j  !< Cell.
   real(R8P),                    intent(out) :: ar, c(2) !< Signed area, centroid.
   real(R8P)    :: pu(4), pv(4), a(2), b(2), cr, mom(2) !< Work.
   integer(I4P) :: k                                   !< Counter.

   pu = [blk%u(i, j), blk%u(i + 1, j), blk%u(i + 1, j + 1), blk%u(i, j + 1)]
   pv = [blk%v(i, j), blk%v(i + 1, j), blk%v(i + 1, j + 1), blk%v(i, j + 1)]
   ar = 0._R8P ; mom = 0._R8P
   do k = 2, 3
      a = [pu(k) - pu(1), pv(k) - pv(1)] ; b = [pu(k + 1) - pu(1), pv(k + 1) - pv(1)]
      cr = a(1) * b(2) - a(2) * b(1)
      ar = ar + cr ; mom = mom + cr * (a + b)
   enddo
   c = mom / (3._R8P * ar) + [pu(1), pv(1)]
   ar = 0.5_R8P * ar
   endsubroutine polygon_centroid

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

   subroutine interp_transfer(blk, x, kd, rcv, cell, xi, phi, pm)
   !< I: the major receptors bilinear from the minor cell centres (located in the block of the cell centres), then the
   !< minor outer row bilinear from the major cell centres.
   type(curvilinear_block_2d_t), intent(in)    :: blk       !< Minor block.
   real(R8P),                    intent(in)    :: x(:)      !< Major grid.
   integer(I4P),                 intent(in)    :: kd(:,:)   !< Major kinds.
   integer(I4P),                 intent(in)    :: rcv(:,:)  !< Major receptors [2, nr].
   integer(I4P),                 intent(in)    :: cell(:,:) !< Their cells in the block of the minor cell centres.
   real(R8P),                    intent(in)    :: xi(:,:)   !< Their local coordinates.
   real(R8P),                    intent(inout) :: phi(:,:), pm(:,:) !< Major, minor values.
   real(R8P)    :: h, c(2), s(2), w(2)                      !< Work.
   integer(I4P) :: ni, nj, i, k, i0, j0                     !< Work.

   h = x(2) - x(1) ; ni = blk%ni ; nj = blk%nj
   do k = 1, size(rcv, 2, kind=I4P)
      i0 = cell(1, k) ; j0 = cell(2, k)
      phi(rcv(1, k), rcv(2, k)) = (1 - xi(1, k)) * (1 - xi(2, k)) * dual_value(pm, i0, j0) + &
                                  xi(1, k) * (1 - xi(2, k)) * dual_value(pm, i0 + 1_I4P, j0) + &
                                  xi(1, k) * xi(2, k) * dual_value(pm, i0 + 1_I4P, j0 + 1_I4P) + &
                                  (1 - xi(1, k)) * xi(2, k) * dual_value(pm, i0, j0 + 1_I4P)
   enddo
   do i = 1, ni - 1
      c = centre(blk, i, nj - 1_I4P)
      s = (c - (x(1) + 0.5_R8P * h)) / h ; i0 = int(s(1)) + 1 ; j0 = int(s(2)) + 1 ; w = s - [i0 - 1, j0 - 1]
      if (any(kd(i0:i0 + 1, j0:j0 + 1) /= 2)) error stop 'interp_transfer: a minor receptor with a major non-field donor'
      pm(i, nj - 1) = (1 - w(1)) * (1 - w(2)) * phi(i0, j0) + w(1) * (1 - w(2)) * phi(i0 + 1, j0) + &
                      w(1) * w(2) * phi(i0 + 1, j0 + 1) + (1 - w(1)) * w(2) * phi(i0, j0 + 1)
   enddo
   endsubroutine interp_transfer

   subroutine transfer_once(blk, x, cons_b, cons_i, acc_b, acc_i)
   !< One transfer of a smooth field (`smooth`) through the assembly of I and B, both directions (the major receptors
   !< from the minor cells, then the minor outer row from the major cells), with no time step: the transfer alone, free
   !< of the composite-mass bookkeeping. Major cells hold `smooth` at their centres, minor cells at their centroids.
   !< cons: |M - D| / D summed over the two directions, M the mass the transfer puts in the receptors (sum of
   !< |V_c| q_c), D the donor-side mass over the same region (sum over the intersections of |V_c n V_c'| times the
   !< donor's linear reconstruction, the one B uses, integrated over V_c n V_c': exact at the intersection centroid).
   !< B conserves it by construction (to round-off); I, which interpolates point values, does not.
   !< acc: |M - E| / E, E the exact integral of `smooth` over the receptors (7-point degree-5 rule on triangles): the
   !< accuracy of each transfer, second order for both.
   type(curvilinear_block_2d_t), intent(in)  :: blk                    !< Minor block.
   real(R8P),                    intent(in)  :: x(:)                   !< Major grid.
   real(R8P),                    intent(out) :: cons_b, cons_i, acc_b, acc_i !< Measures.
   type(overlap_weights_t)      :: ow                                  !< Intersections.
   type(curvilinear_block_2d_t) :: dual                                !< Block of the minor cell centres.
   integer(I4P), allocatable    :: tstart(:), tcell(:), kd(:,:), cell(:,:), lst(:), rcv(:,:) !< Work.
   real(R8P),    allocatable    :: xm(:,:,:), gx(:,:,:), gm(:,:,:), phi(:,:), pm(:,:), phi0(:,:), pm0(:,:), vm(:,:)
   real(R8P),    allocatable    :: xi(:,:), qx(:), qy(:)               !< Work.
   real(R8P)    :: h, c(2), rr, d(2), mm(2), e(2), ar, xc(2)           !< Work.
   integer(I4P) :: nx, ni, nj, i, j, k, l, r, a, b, nr, pass          !< Work.

   nx = size(x, kind=I4P) - 1 ; h = x(2) - x(1) ; ni = blk%ni ; nj = blk%nj
   allocate(kd(nx, nx), phi0(nx, nx), pm0(ni - 1, nj - 1), vm(ni - 1, nj - 1))
   nr = 0
   do j = 1, nx
      do i = 1, nx
         c = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (x(j) + x(j + 1))] ; rr = norm2(c)
         kd(i, j) = 2 ; if (rr < 1.6_R8P) kd(i, j) = 1 ; if (rr < 1.4_R8P) kd(i, j) = 0
         if (kd(i, j) == 1) nr = nr + 1
         phi0(i, j) = smooth(c)
      enddo
   enddo
   call picard_setup(blk, x, kd, ow, tstart, tcell, xm)
   do j = 1, nj - 1
      do i = 1, ni - 1
         vm(i, j) = abs(quad(blk, i, j)) ; pm0(i, j) = smooth(xm(:, i, j))
      enddo
   enddo
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
   if (any(lst == LOC_OUTSIDE)) error stop 'transfer_once: a major receptor outside the minor block'
   allocate(gx(2, nx, nx), gm(2, ni - 1, nj - 1))
   ! the exact integrals over the receptors: the major receptors (squares), the minor outer row (polygons)
   e = 0._R8P
   do k = 1, nr
      i = rcv(1, k) ; j = rcv(2, k)
      e(1) = e(1) + tri_integral([x(i), x(j)], [x(i + 1), x(j)], [x(i + 1), x(j + 1)]) + &
                    tri_integral([x(i), x(j)], [x(i + 1), x(j + 1)], [x(i), x(j + 1)])
   enddo
   b = nj - 1
   do a = 1, ni - 1
      e(2) = e(2) + tri_integral(vx2(blk, a, b), vx2(blk, a + 1_I4P, b), vx2(blk, a + 1_I4P, b + 1_I4P)) + &
                    tri_integral(vx2(blk, a, b), vx2(blk, a + 1_I4P, b + 1_I4P), vx2(blk, a, b + 1_I4P))
   enddo
   do pass = 1, 2 ! 1: B, 2: I
      phi = phi0 ; pm = pm0
      ! D(1): the minor donors' mass over the major receptors, before the transfer
      call minor_gradients(blk, xm, pm, gm)
      d = 0._R8P
      do k = 1, nr
         l = rcv(1, k) + (rcv(2, k) - 1) * nx
         do r = tstart(l), tstart(l + 1) - 1
            call receptor_of(ow, tcell(r), a) ; i = mod(a - 1, ni - 1) + 1 ; j = (a - 1) / (ni - 1) + 1
            d(1) = d(1) + ow%volume(tcell(r)) * (pm(i, j) + dot_product(gm(:, i, j), ow%centroid(:, tcell(r)) - &
                                                                        xm(:, i, j)))
         enddo
      enddo
      if (pass == 1) then
         call picard_transfer(blk, x, kd, ow, tstart, tcell, xm, gx, gm, phi, pm)
      else
         call interp_transfer(blk, x, kd, rcv, cell, xi, phi, pm)
      endif
      ! D(2): the major donors' mass over the minor outer row, with the major values the transfer used
      call major_gradients(x, kd, phi, gx)
      mm = 0._R8P
      do k = 1, nr
         mm(1) = mm(1) + h * h * phi(rcv(1, k), rcv(2, k))
      enddo
      do a = 1, ni - 1
         r = a + (b - 1) * (ni - 1)
         mm(2) = mm(2) + vm(a, b) * pm(a, b)
         do l = ow%start(r), ow%start(r + 1) - 1
            i = ow%donor(1, l) ; j = ow%donor(2, l)
            xc = [0.5_R8P * (x(i) + x(i + 1)), 0.5_R8P * (x(j) + x(j + 1))]
            d(2) = d(2) + ow%volume(l) * (phi(i, j) + dot_product(gx(:, i, j), ow%centroid(:, l) - xc))
         enddo
      enddo
      ar = (abs(mm(1) - d(1)) + abs(mm(2) - d(2))) / (d(1) + d(2))
      if (pass == 1) then
         cons_b = ar ; acc_b = (abs(mm(1) - e(1)) + abs(mm(2) - e(2))) / (e(1) + e(2))
      else
         cons_i = ar ; acc_i = (abs(mm(1) - e(1)) + abs(mm(2) - e(2))) / (e(1) + e(2))
      endif
   enddo
   endsubroutine transfer_once

   pure function vx2(blk, i, j) result(p)
   !< Vertex (i, j) of a 2-D block.
   type(curvilinear_block_2d_t), intent(in) :: blk  !< Block.
   integer(I4P),                 intent(in) :: i, j !< Vertex.
   real(R8P)                                :: p(2) !< Coordinates.

   p = [blk%u(i, j), blk%v(i, j)]
   endfunction vx2

   pure real(R8P) function tri_integral(p1, p2, p3) result(s)
   !< The integral of `smooth` over the triangle p1 p2 p3 (7-point rule of degree 5).
   real(R8P), intent(in) :: p1(2), p2(2), p3(2) !< Vertices.
   real(R8P), parameter :: A1 = (6._R8P - sqrt(15._R8P)) / 21._R8P, A2 = (6._R8P + sqrt(15._R8P)) / 21._R8P
   real(R8P), parameter :: W0 = 9._R8P / 40._R8P, W1 = (155._R8P - sqrt(15._R8P)) / 1200._R8P, &
                           W2 = (155._R8P + sqrt(15._R8P)) / 1200._R8P
   real(R8P) :: ar !< Area.

   ar = 0.5_R8P * abs((p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1)))
   s = W0 * smooth((p1 + p2 + p3) / 3._R8P) + &
       W1 * (smooth(A1 * p1 + A1 * p2 + (1 - 2 * A1) * p3) + smooth(A1 * p1 + (1 - 2 * A1) * p2 + A1 * p3) + &
             smooth((1 - 2 * A1) * p1 + A1 * p2 + A1 * p3)) + &
       W2 * (smooth(A2 * p1 + A2 * p2 + (1 - 2 * A2) * p3) + smooth(A2 * p1 + (1 - 2 * A2) * p2 + A2 * p3) + &
             smooth((1 - 2 * A2) * p1 + A2 * p2 + A2 * p3))
   s = s * ar
   endfunction tri_integral

   pure real(R8P) function smooth(p)
   !< A smooth positive field for the one-transfer measure.
   real(R8P), intent(in) :: p(2) !< Point.

   smooth = 1._R8P + 0.5_R8P * sin(1.7_R8P * p(1) + 0.3_R8P) * cos(1.3_R8P * p(2) - 0.2_R8P)
   endfunction smooth

   real(R8P) function mass_p(x, alpha, phi, pm, vm) result(tot)
   !< Composite mass on the exact partition (Kahan): the major cells by their fluid fraction outside the minor polygon,
   !< all minor cells.
   real(R8P), intent(in) :: x(:)             !< Major grid.
   real(R8P), intent(in) :: alpha(:,:)       !< Major fluid fractions.
   real(R8P), intent(in) :: phi(:,:)         !< Major values.
   real(R8P), intent(in) :: pm(:,:), vm(:,:) !< Minor values, volumes.
   real(R8P)             :: comp             !< Compensation.
   integer(I4P)          :: i, j             !< Counters.

   tot = 0._R8P ; comp = 0._R8P
   do j = 1, size(phi, 2, kind=I4P)
      do i = 1, size(phi, 1, kind=I4P)
         if (alpha(i, j) > 0._R8P) call kahan(tot, comp, alpha(i, j) * phi(i, j) * (x(i + 1) - x(i)) * (x(j + 1) - x(j)))
      enddo
   enddo
   do j = 1, size(pm, 2, kind=I4P)
      do i = 1, size(pm, 1, kind=I4P)
         call kahan(tot, comp, pm(i, j) * vm(i, j))
      enddo
   enddo
   endfunction mass_p

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
   ! ------------------------------------------------------------------------------------------------------------- 3-D
   subroutine patched_3d(blk, x, defect, vmax, vmin)
   !< 3-D A: the interface geometry (closure per piece, face sums, total volume, sum sf = - wet_normal), then the scalar
   !< through Wang's patched interface; defect = max over the steps of |M + out - M0| / M0.
   type(curvilinear_block_t), intent(in)  :: blk(:)     !< Minor block.
   real(R8P),                 intent(in)  :: x(:)       !< Major grid edges (all directions).
   real(R8P),                 intent(out) :: defect     !< Mass defect.
   real(R8P),                 intent(out) :: vmax, vmin !< Extreme values reached.
   type(overset_patch_3d_t)  :: pt                      !< Interface.
   type(merge_t)             :: mt                      !< Merge targets.
   type(minor_3d_t)          :: mg                      !< Minor geometry.
   integer(I4P), allocatable :: g(:)                    !< Group (root piece) of each piece (0: none).
   real(R8P),    allocatable :: vg(:), phi(:), r(:), mgr(:), outp(:) !< Group volumes, values, residuals, masses, outflow.
   real(R8P),    allocatable :: pm(:,:,:), rm(:,:,:), outm(:,:,:) !< Minor values, residuals, outflow.
   real(R8P)    :: h, out, m0, mtot, fl, sw(3), pf(6), cmax, wmax, fmax, tot, comp, dt, cfl, tm, ex !< Work.
   integer(I4P) :: st, np, p, k, step, nx, t, nsm, ndg, nal, nbit, nrnd, nsteps, ni, nj, nk !< Work.
   integer(I8P) :: c0, c1, cr                                  !< Clock.

   call system_clock(c0, cr)
   call overset_patch_3d(blk, [1_I4P], x, x, x, pt, status=st)
   call system_clock(c1)
   if (st /= OVERSET_PATCH_OK) error stop 'overset_patch_3d'
   np = pt%cut%np ; nx = pt%cut%nx ; h = x(2) - x(1)
   ! the geometry: closure per piece, sum of a piece's smallest faces against its wet normal, faces, total volume
   cmax = 0._R8P ; wmax = 0._R8P ; nbit = 0 ; nrnd = 0 ; tot = 0._R8P ; comp = 0._R8P
   do p = 1, np
      sw = 0._R8P
      do k = pt%piece_sf_start(p), pt%piece_sf_start(p + 1) - 1
         sw = sw + pt%sf_area(:, k)
      enddo
      pf = pt%cut%piece_face(:, p)
      cmax = max(cmax, norm2(h * h * [pf(2) - pf(1), pf(4) - pf(3), pf(6) - pf(5)] + sw) / h ** 2)
      if (all(sw == -pt%cut%piece_wet_normal(:, p))) then
         nbit = nbit + 1
      else
         nrnd = nrnd + 1 ; wmax = max(wmax, norm2(sw + pt%cut%piece_wet_normal(:, p)) / h ** 2)
      endif
      call kahan(tot, comp, pt%cut%piece_volume(p))
   enddo
   fmax = 0._R8P
   do k = 1, pt%nf
      if (norm2(pt%face_vector(:, k)) > 0._R8P) &
         fmax = max(fmax, norm2(pt%face_covered(:, k) - pt%face_vector(:, k)) / norm2(pt%face_vector(:, k)))
   enddo
   ex = (x(nx + 1) - x(1)) ** 3 - pt%boundary%get_volume()
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F7.3,A)', '    A: ', nx, '^3, minor faces ', pt%nf, ': holes ', count(pt%hole), &
         ', pieces ', np, ', smallest faces ', pt%ns, ', degenerate ', count(pt%degenerate), '; overset_patch_3d ', &
         real(c1 - c0, R8P) / cr, ' s'
   print '(A,ES9.2,A,ES9.2,A,ES9.2,A,I0,A,I0,A,ES9.2,A)', '       closure / h^2 ', cmax, '; covered vs face ', fmax, &
         '; volume vs box - polyhedron ', abs(tot - ex) / ex, '; sum sf = -wet_normal: bit ', nbit, ', round-off ', nrnd, &
         ' (max / h^2 ', wmax, ')'
   if (cmax > 1e-12_R8P .or. fmax > 1e-13_R8P .or. abs(tot - ex) > 1e-13_R8P * ex .or. wmax > 1e-13_R8P .or. pt%nuc > 0) &
      error stop 'patched_3d: geometry'
   ! groups: degenerate pieces into merge_into (then its final target), small pieces into their final target
   call merge_targets(pt%links, pt%cut%piece_cell, pt%cut%piece_volume, pt%cut%piece_wet_normal, x, x, mt, z=x)
   allocate(g(np), vg(np), phi(np), r(np), mgr(np), outp(np)) ; nsm = 0 ; ndg = 0 ; nal = 0
   do p = 1, np
      g(p) = p
      if (pt%degenerate(p)) then
         t = pt%merge_into(p) ; g(p) = 0
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
   vg = 0._R8P ; mgr = 0._R8P ; phi = 0._R8P
   do p = 1, np
      if (g(p) == 0) cycle
      vg(g(p)) = vg(g(p)) + pt%cut%piece_volume(p)
      mgr(g(p)) = mgr(g(p)) + pt%cut%piece_volume(p) * blob3(pt%cut%piece_centroid(:, p))
   enddo
   do p = 1, np
      if (g(p) > 0) phi(p) = mgr(g(p)) / vg(g(p))
   enddo
   call minor_geometry(blk(1), mg)
   ni = size(mg%vol, 1, kind=I4P) ; nj = size(mg%vol, 2, kind=I4P) ; nk = size(mg%vol, 3, kind=I4P)
   allocate(pm(ni, nj, nk), rm(ni, nj, nk), outm(ni, nj, nk))
   pm = mg%val0
   ! the time step: the outflow of every group and minor cell (steady velocity) under CFL 0.9
   outp = 0._R8P ; outm = 0._R8P
   call residuals_3d(pt, mg, x, g, phi, pm, r, rm, out, outp, outm, .true.)
   mgr = 0._R8P
   do p = 1, np
      if (g(p) > 0) mgr(g(p)) = mgr(g(p)) + outp(p)
   enddo
   cfl = 0._R8P
   do p = 1, np
      if (g(p) == p .and. vg(p) > 0._R8P) cfl = max(cfl, mgr(p) / vg(p))
   enddo
   cfl = max(cfl, maxval(outm / mg%vol))
   tm = 1.6_R8P ; dt = 0.9_R8P / cfl ; nsteps = ceiling(tm / dt) ; dt = tm / nsteps
   m0 = mass_3d(g, phi, vg, pm, mg%vol) ; out = 0._R8P ; defect = 0._R8P ; vmax = 0._R8P ; vmin = 0._R8P
   do step = 1, nsteps
      call residuals_3d(pt, mg, x, g, phi, pm, r, rm, fl, outp, outm, .false.)
      out = out + dt * fl
      mgr = 0._R8P
      do p = 1, np
         if (g(p) > 0) mgr(g(p)) = mgr(g(p)) + r(p)
      enddo
      do p = 1, np
         if (g(p) == p) mgr(p) = phi(p) * vg(p) - dt * mgr(p)
      enddo
      do p = 1, np
         if (g(p) > 0) phi(p) = mgr(g(p)) / vg(g(p))
      enddo
      pm = pm - dt * rm / mg%vol
      mtot = mass_3d(g, phi, vg, pm, mg%vol)
      defect = max(defect, abs(mtot + out - m0) / m0)
      vmax = max(vmax, maxval(phi), maxval(pm)) ; vmin = min(vmin, minval(phi), minval(pm))
   enddo
   print '(A,I0,A,I0,A,I0,A,I0,A,ES10.3,A,ES12.5,A,ES12.5)', '       small merged ', nsm, ', degenerate folded ', ndg, &
         ' (alone ', nal, '), steps ', nsteps, ' of ', dt, '; mass ', mtot, ', left the box ', out
   endsubroutine patched_3d

   subroutine residuals_3d(pt, mg, x, g, phi, pm, r, rm, out, outp, outm, cflpass)
   !< The upwind residuals of A in 3-D (r: major pieces, rm: minor cells) and the outflow rate through the box; with
   !< cflpass the outflow of every piece and minor cell (outp, outm) instead (values unused).
   type(overset_patch_3d_t), intent(in)    :: pt          !< Interface.
   type(minor_3d_t),         intent(in)    :: mg          !< Minor geometry.
   real(R8P),                intent(in)    :: x(:)        !< Major grid.
   integer(I4P),             intent(in)    :: g(:)        !< Groups.
   real(R8P),                intent(in)    :: phi(:), pm(:,:,:) !< Values.
   real(R8P),                intent(out)   :: r(:), rm(:,:,:) !< Residuals.
   real(R8P),                intent(out)   :: out         !< Outflow rate through the box.
   real(R8P),                intent(inout) :: outp(:), outm(:,:,:) !< Outflows (cflpass).
   logical,                  intent(in)    :: cflpass     !< Outflows only.
   real(R8P)    :: f, fl, mid(3), nrm(3), area !< Work.
   integer(I4P) :: np, nx, p, q, k, t, c(3), m(3) !< Work.

   np = pt%cut%np ; nx = pt%cut%nx ; r = 0._R8P ; rm = 0._R8P ; out = 0._R8P
   ! major: across the cell faces (each linked pair once), the shared fluid area
   do p = 1, np
      do k = pt%links%start(p), pt%links%start(p + 1) - 1
         q = pt%links%piece(k)
         if (q <= p) cycle
         call face_3d(pt%cut%piece_cell(:, p), pt%links%face(k), x, mid, nrm, area)
         f = dot_product(vel3(mid), nrm) * pt%links%area(k)
         if (cflpass) then
            outp(p) = outp(p) + max(f, 0._R8P) ; outp(q) = outp(q) + max(-f, 0._R8P) ; cycle
         endif
         fl = f * merge(phi(p), phi(q), f > 0._R8P)
         r(p) = r(p) + fl ; r(q) = r(q) - fl
      enddo
   enddo
   ! major: the box boundary (inflow 0, outflow upwind)
   do p = 1, np
      c = pt%cut%piece_cell(:, p)
      do t = 1, 6
         if (.not. (c((t + 1) / 2) == merge(1_I4P, nx, mod(t, 2_I4P) == 1))) cycle
         call face_3d(c, t, x, mid, nrm, area)
         f = dot_product(vel3(mid), nrm) * pt%cut%piece_face(t, p) * area
         if (f <= 0._R8P) cycle
         if (cflpass) then
            outp(p) = outp(p) + f ; cycle
         endif
         r(p) = r(p) + f * phi(p) ; out = out + f * phi(p)
      enddo
   enddo
   ! the smallest faces: major piece <-> minor cell (sf_area points into the minor cell)
   do k = 1, pt%ns
      p = pt%sf_piece(k) ; m = pt%sf_minor(2:4, k)
      f = dot_product(vel3(pt%sf_centroid(:, k)), pt%sf_area(:, k))
      if (g(p) == 0 .and. f /= 0._R8P) error stop 'residuals_3d: a flux through a piece without a group'
      if (cflpass) then
         outp(p) = outp(p) + max(f, 0._R8P) ; outm(m(1), m(2), m(3)) = outm(m(1), m(2), m(3)) + max(-f, 0._R8P) ; cycle
      endif
      fl = f * merge(phi(p), pm(m(1), m(2), m(3)), f > 0._R8P)
      r(p) = r(p) + fl ; rm(m(1), m(2), m(3)) = rm(m(1), m(2), m(3)) - fl
   enddo
   ! minor: its interior faces (periodic in i, the poles closed, the inner face a wall, the outer the smallest faces)
   do k = 1, mg%nf
      f = dot_product(vel3(mg%fmid(:, k)), mg%fvec(:, k))
      c = mg%fa(:, k) ; m = mg%fb(:, k)
      if (cflpass) then
         outm(c(1), c(2), c(3)) = outm(c(1), c(2), c(3)) + max(f, 0._R8P)
         outm(m(1), m(2), m(3)) = outm(m(1), m(2), m(3)) + max(-f, 0._R8P) ; cycle
      endif
      fl = f * merge(pm(c(1), c(2), c(3)), pm(m(1), m(2), m(3)), f > 0._R8P)
      rm(c(1), c(2), c(3)) = rm(c(1), c(2), c(3)) + fl ; rm(m(1), m(2), m(3)) = rm(m(1), m(2), m(3)) - fl
   enddo
   endsubroutine residuals_3d

   subroutine minor_geometry(blk, mg)
   !< The minor shell's cells (volume from the face vectors, vertex-mean centre, the initial values) and its interior faces
   !< (a -> b, the face vector 1/2 (c3 - c1) x (c4 - c2) from a to b, the vertex mean): i periodic, j between rows (the
   !< pole faces are collapsed), k between layers (the inner face is a wall, the outer the smallest faces).
   type(curvilinear_block_t), intent(in)    :: blk !< Shell block.
   type(minor_3d_t),          intent(inout) :: mg  !< Geometry.
   integer(I4P) :: ni, nj, nk, i, j, k, d, n, a(3), b(3), f, l !< Work.
   real(R8P)    :: c(3, 4), s(3), cm(3), v                     !< Work.

   ni = blk%ni - 1 ; nj = blk%nj - 1 ; nk = blk%nk - 1
   allocate(mg%vol(ni, nj, nk), mg%cen(3, ni, nj, nk), mg%val0(ni, nj, nk))
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            v = 0._R8P
            do f = 1, 6
               call cell_face(blk, [i, j, k], f, c)
               s = 0.5_R8P * cross3(c(:, 3) - c(:, 1), c(:, 4) - c(:, 2)) ! along +e of the face's axis
               v = v + merge(-1._R8P, 1._R8P, mod(f, 2_I4P) == 1) * dot_product(0.25_R8P * sum(c, dim=2), s) / 3._R8P
            enddo
            mg%vol(i, j, k) = abs(v)
            cm = 0._R8P
            do l = 0, 7
               cm = cm + vtx(blk, [i + mod(l, 2_I4P), j + mod(l / 2_I4P, 2_I4P), k + l / 4_I4P])
            enddo
            mg%cen(:, i, j, k) = cm / 8._R8P ; mg%val0(i, j, k) = blob3(mg%cen(:, i, j, k))
         enddo
      enddo
   enddo
   mg%nf = ni * nj * nk + ni * (nj - 1) * nk + ni * nj * (nk - 1)
   allocate(mg%fa(3, mg%nf), mg%fb(3, mg%nf), mg%fvec(3, mg%nf), mg%fmid(3, mg%nf)) ; n = 0
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            do d = 1, 3
               a = [i, j, k] ; b = a ; b(d) = b(d) + 1
               if (d == 1) b(1) = mod(i, ni) + 1
               if (d == 2 .and. j == nj) cycle
               if (d == 3 .and. k == nk) cycle
               call cell_face(blk, a, 2_I4P * d, c)
               n = n + 1 ; mg%fa(:, n) = a ; mg%fb(:, n) = b ; mg%fmid(:, n) = 0.25_R8P * sum(c, dim=2)
               mg%fvec(:, n) = 0.5_R8P * cross3(c(:, 3) - c(:, 1), c(:, 4) - c(:, 2))
               if (dot_product(mg%fvec(:, n), mg%cen(:, b(1), b(2), b(3)) - mg%cen(:, i, j, k)) < 0._R8P) &
                  mg%fvec(:, n) = -mg%fvec(:, n)
            enddo
         enddo
      enddo
   enddo
   mg%nf = n
   endsubroutine minor_geometry

   subroutine cell_face(blk, c, f, q)
   !< The corners of face f (1-6: -i, +i, -j, +j, -k, +k) of cell c, cyclic.
   type(curvilinear_block_t), intent(in)  :: blk     !< Block.
   integer(I4P),              intent(in)  :: c(3), f !< Cell, face.
   real(R8P),                 intent(out) :: q(3, 4) !< Corners.
   integer(I4P) :: nd, d1, d2, p(3), e1(3), e2(3)    !< Work.

   nd = (f + 1) / 2 ; d1 = mod(nd, 3_I4P) + 1 ; d2 = mod(nd + 1_I4P, 3_I4P) + 1
   p = c ; if (mod(f, 2_I4P) == 0) p(nd) = p(nd) + 1
   e1 = 0 ; e1(d1) = 1 ; e2 = 0 ; e2(d2) = 1
   q(:, 1) = vtx(blk, p) ; q(:, 2) = vtx(blk, p + e1) ; q(:, 3) = vtx(blk, p + e1 + e2) ; q(:, 4) = vtx(blk, p + e2)
   endsubroutine cell_face

   pure function vtx(blk, p) result(v)
   !< Vertex p of a block.
   type(curvilinear_block_t), intent(in) :: blk  !< Block.
   integer(I4P),              intent(in) :: p(3) !< Vertex.
   real(R8P)                             :: v(3) !< Coordinates.

   v = [blk%x(p(1), p(2), p(3)), blk%y(p(1), p(2), p(3)), blk%z(p(1), p(2), p(3))]
   endfunction vtx

   subroutine chimera_3d(blk, x, defect)
   !< 3-D I: the shell coupled by interpolation (module notes); defect as in A on the composite mass (the major cells with
   !< centre at r >= 2, all minor cells), same time step rule.
   type(curvilinear_block_t), intent(in)  :: blk    !< Minor block.
   real(R8P),                 intent(in)  :: x(:)   !< Major grid.
   real(R8P),                 intent(out) :: defect !< Mass defect.
   type(minor_3d_t)          :: mg                  !< Minor geometry.
   integer(I4P), allocatable :: kd(:,:,:)           !< Major kinds: 0 hole, 1 receptor, 2 field.
   real(R8P),    allocatable :: phi(:,:,:), r(:,:,:), pm(:,:,:), rm(:,:,:) !< Values, residuals.
   real(R8P)    :: h, out, m0, mtot, f, fl, c(3), dt, cfl, tm, s(3), w(3) !< Work.
   integer(I4P) :: nx, ni, nj, nk, i, j, k, t, step, nsteps, nr, nbad, b(3), a(3), l, i0(3) !< Work.
   real(R8P),    allocatable :: outc(:,:,:), outm(:,:,:) !< Outflows.

   nx = size(x, kind=I4P) - 1 ; h = x(2) - x(1)
   allocate(kd(nx, nx, nx), phi(nx, nx, nx), r(nx, nx, nx), outc(nx, nx, nx)) ; nr = 0
   do k = 1, nx
      do j = 1, nx
         do i = 1, nx
            c = 0.5_R8P * [x(i) + x(i + 1), x(j) + x(j + 1), x(k) + x(k + 1)]
            kd(i, j, k) = 2 ; if (norm2(c) < 1.6_R8P) kd(i, j, k) = 1 ; if (norm2(c) < 1.2_R8P) kd(i, j, k) = 0
            if (kd(i, j, k) == 1) nr = nr + 1
            phi(i, j, k) = blob3(c)
         enddo
      enddo
   enddo
   nbad = 0 ! a field cell next to a hole: the fringe is too thin
   do k = 1, nx
      do j = 1, nx
         do i = 1, nx
            if (kd(i, j, k) /= 2) cycle
            do t = 1, 6
               b = [i, j, k] ; b((t + 1) / 2) = b((t + 1) / 2) + merge(-1, 1, mod(t, 2_I4P) == 1)
               if (any(b < 1) .or. any(b > nx)) cycle
               if (kd(b(1), b(2), b(3)) == 0) nbad = nbad + 1
            enddo
         enddo
      enddo
   enddo
   if (nbad > 0) error stop 'chimera_3d: fringe too thin'
   call minor_geometry(blk, mg)
   ni = size(mg%vol, 1, kind=I4P) ; nj = size(mg%vol, 2, kind=I4P) ; nk = size(mg%vol, 3, kind=I4P)
   allocate(pm(ni, nj, nk), rm(ni, nj, nk), outm(ni, nj, nk)) ; pm = mg%val0
   do j = 1, nj ! the minor outer layer takes field cells only
      do i = 1, ni
         s = (mg%cen(:, i, j, nk) - (x(1) + 0.5_R8P * h)) / h ; i0 = int(s) + 1
         if (any(kd(i0(1):i0(1) + 1, i0(2):i0(2) + 1, i0(3):i0(3) + 1) /= 2)) &
            error stop 'chimera_3d: a minor receptor with a major non-field donor'
      enddo
   enddo
   ! the time step from the outflows (field cells and minor field cells)
   outc = 0._R8P ; outm = 0._R8P
   call chimera_residuals(mg, x, kd, phi, pm, r, rm, out, outc, outm, .true.)
   cfl = max(maxval(outc / h ** 3, mask=kd == 2), maxval(outm(:, :, 1:nk - 1) / mg%vol(:, :, 1:nk - 1)))
   tm = 1.6_R8P ; dt = 0.9_R8P / cfl ; nsteps = ceiling(tm / dt) ; dt = tm / nsteps
   m0 = mass_3d_b(x, phi, pm, mg%vol) ; out = 0._R8P ; defect = 0._R8P
   do step = 1, nsteps
      ! interpolation: major receptors from the minor cell centres (trilinear in the shell's index space), the minor
      ! outer layer from the major cell centres (trilinear)
      do k = 1, nx
         do j = 1, nx
            do i = 1, nx
               if (kd(i, j, k) /= 1) cycle
               phi(i, j, k) = shell_interp(pm, 0.5_R8P * [x(i) + x(i + 1), x(j) + x(j + 1), x(k) + x(k + 1)])
            enddo
         enddo
      enddo
      do j = 1, nj
         do i = 1, ni
            s = (mg%cen(:, i, j, nk) - (x(1) + 0.5_R8P * h)) / h ; i0 = int(s) + 1 ; w = s - (i0 - 1)
            f = 0._R8P
            do l = 0, 7
               a = [mod(l, 2_I4P), mod(l / 2_I4P, 2_I4P), l / 4_I4P]
               f = f + product(merge(w, 1._R8P - w, a == 1)) * phi(i0(1) + a(1), i0(2) + a(2), i0(3) + a(3))
            enddo
            pm(i, j, nk) = f
         enddo
      enddo
      call chimera_residuals(mg, x, kd, phi, pm, r, rm, fl, outc, outm, .false.)
      out = out + dt * fl
      where (kd == 2) phi = phi - dt * r / h ** 3
      pm(:, :, 1:nk - 1) = pm(:, :, 1:nk - 1) - dt * rm(:, :, 1:nk - 1) / mg%vol(:, :, 1:nk - 1)
      mtot = mass_3d_b(x, phi, pm, mg%vol)
      defect = max(defect, abs(mtot + out - m0) / m0)
   enddo
   print '(A,I0,A,I0,A,ES10.3,A,ES12.5,A,ES12.5)', '    I: major receptors ', nr, ', steps ', nsteps, ' of ', dt, &
         '; mass ', mtot, ', left the box ', out
   endsubroutine chimera_3d

   subroutine chimera_residuals(mg, x, kd, phi, pm, r, rm, out, outc, outm, cflpass)
   !< The upwind residuals of I in 3-D: major faces with a field cell on one side at least, the box, the minor faces of
   !< its field cells (layers below the outer one); with cflpass the outflows instead.
   type(minor_3d_t), intent(in)    :: mg          !< Minor geometry.
   real(R8P),        intent(in)    :: x(:)        !< Major grid.
   integer(I4P),     intent(in)    :: kd(:,:,:)   !< Major kinds.
   real(R8P),        intent(in)    :: phi(:,:,:), pm(:,:,:) !< Values.
   real(R8P),        intent(out)   :: r(:,:,:), rm(:,:,:) !< Residuals.
   real(R8P),        intent(out)   :: out         !< Outflow rate through the box.
   real(R8P),        intent(inout) :: outc(:,:,:), outm(:,:,:) !< Outflows (cflpass).
   logical,          intent(in)    :: cflpass     !< Outflows only.
   real(R8P)    :: f, fl, mid(3), nrm(3), area    !< Work.
   integer(I4P) :: nx, nk, i, j, k, t, b(3), c(3), m(3) !< Work.

   nx = size(kd, 1, kind=I4P) ; nk = size(pm, 3, kind=I4P) ; r = 0._R8P ; rm = 0._R8P ; out = 0._R8P
   do k = 1, nx
      do j = 1, nx
         do i = 1, nx
            do t = 1, 6
               b = [i, j, k] ; b((t + 1) / 2) = b((t + 1) / 2) + merge(-1, 1, mod(t, 2_I4P) == 1)
               call face_3d([i, j, k], t, x, mid, nrm, area)
               f = dot_product(vel3(mid), nrm) * area
               if (any(b < 1) .or. any(b > nx)) then ! the box
                  if (f <= 0._R8P .or. kd(i, j, k) /= 2) cycle
                  if (cflpass) then
                     outc(i, j, k) = outc(i, j, k) + f ; cycle
                  endif
                  r(i, j, k) = r(i, j, k) + f * phi(i, j, k) ; out = out + f * phi(i, j, k) ; cycle
               endif
               if (mod(t, 2_I4P) == 1) cycle ! each interior face once, from its low cell
               if (kd(i, j, k) /= 2 .and. kd(b(1), b(2), b(3)) /= 2) cycle
               if (cflpass) then
                  outc(i, j, k) = outc(i, j, k) + max(f, 0._R8P)
                  outc(b(1), b(2), b(3)) = outc(b(1), b(2), b(3)) + max(-f, 0._R8P) ; cycle
               endif
               fl = f * merge(phi(i, j, k), phi(b(1), b(2), b(3)), f > 0._R8P)
               r(i, j, k) = r(i, j, k) + fl ; r(b(1), b(2), b(3)) = r(b(1), b(2), b(3)) - fl
            enddo
         enddo
      enddo
   enddo
   do k = 1, mg%nf
      c = mg%fa(:, k) ; m = mg%fb(:, k)
      if (c(3) == nk .and. m(3) == nk) cycle ! faces between two receptors
      f = dot_product(vel3(mg%fmid(:, k)), mg%fvec(:, k))
      if (cflpass) then
         outm(c(1), c(2), c(3)) = outm(c(1), c(2), c(3)) + max(f, 0._R8P)
         outm(m(1), m(2), m(3)) = outm(m(1), m(2), m(3)) + max(-f, 0._R8P) ; cycle
      endif
      fl = f * merge(pm(c(1), c(2), c(3)), pm(m(1), m(2), m(3)), f > 0._R8P)
      rm(c(1), c(2), c(3)) = rm(c(1), c(2), c(3)) + fl ; rm(m(1), m(2), m(3)) = rm(m(1), m(2), m(3)) - fl
   enddo
   endsubroutine chimera_residuals

   pure real(R8P) function shell_interp(pm, p) result(val)
   !< The minor values at p, trilinear in the shell's index space of cell centres (theta periodic, phi clamped at the
   !< poles, r over the layers below the outer one).
   real(R8P), intent(in) :: pm(:,:,:) !< Minor values [na, np, nl].
   real(R8P), intent(in) :: p(3)      !< Point.
   real(R8P)    :: rr, th, ph, s(3), w(3) !< Work.
   integer(I4P) :: na, np, nl, i0, j0, k0, l, a(3), ii !< Work.

   na = size(pm, 1, kind=I4P) ; np = size(pm, 2, kind=I4P) ; nl = size(pm, 3, kind=I4P)
   rr = norm2(p) ; ph = acos(max(-1._R8P, min(1._R8P, p(3) / rr))) ; th = atan2(p(2), p(1))
   if (th < 0._R8P) th = th + 2 * PI
   s = [th * na / (2 * PI) + 0.5_R8P, ph * np / PI + 0.5_R8P, (rr - 1._R8P) * nl + 0.5_R8P]
   i0 = floor(s(1)) ; j0 = min(max(floor(s(2)), 1), np - 1) ; k0 = min(max(floor(s(3)), 1), nl - 2)
   w = [s(1) - i0, s(2) - j0, s(3) - k0] ; w = min(max(w, 0._R8P), 1._R8P)
   val = 0._R8P
   do l = 0, 7
      a = [mod(l, 2_I4P), mod(l / 2_I4P, 2_I4P), l / 4_I4P]
      ii = modulo(i0 + a(1) - 1, na) + 1
      val = val + product(merge(w, 1._R8P - w, a == 1)) * pm(ii, j0 + a(2), k0 + a(3))
   enddo
   endfunction shell_interp

   pure subroutine face_3d(c, t, x, mid, nrm, area)
   !< Face t (1-6: -x, +x, -y, +y, -z, +z) of the major cell c: centre, outward unit normal, area.
   integer(I4P), intent(in)  :: c(3), t        !< Cell, face.
   real(R8P),    intent(in)  :: x(:)           !< Grid edges (all directions).
   real(R8P),    intent(out) :: mid(3), nrm(3), area !< Face.
   integer(I4P)              :: a, d1, d2      !< Axes.

   a = (t + 1) / 2 ; d1 = mod(a, 3_I4P) + 1 ; d2 = mod(a + 1_I4P, 3_I4P) + 1
   mid = 0.5_R8P * [x(c(1)) + x(c(1) + 1), x(c(2)) + x(c(2) + 1), x(c(3)) + x(c(3) + 1)]
   mid(a) = x(c(a) + merge(0, 1, mod(t, 2_I4P) == 1))
   nrm = 0._R8P ; nrm(a) = merge(-1._R8P, 1._R8P, mod(t, 2_I4P) == 1)
   area = (x(c(d1) + 1) - x(c(d1))) * (x(c(d2) + 1) - x(c(d2)))
   endsubroutine face_3d

   real(R8P) function mass_3d(g, phi, vg, pm, vm) result(tot)
   !< Total mass of 3-D A (Kahan): the groups and the minor cells.
   integer(I4P), intent(in) :: g(:)               !< Groups.
   real(R8P),    intent(in) :: phi(:), vg(:)      !< Values, group volumes.
   real(R8P),    intent(in) :: pm(:,:,:), vm(:,:,:) !< Minor values, volumes.
   real(R8P)                :: comp               !< Compensation.
   integer(I4P)             :: p, i, j, k         !< Counters.

   tot = 0._R8P ; comp = 0._R8P
   do p = 1, size(g, kind=I4P)
      if (g(p) == p) call kahan(tot, comp, phi(p) * vg(p))
   enddo
   do k = 1, size(pm, 3, kind=I4P)
      do j = 1, size(pm, 2, kind=I4P)
         do i = 1, size(pm, 1, kind=I4P)
            call kahan(tot, comp, pm(i, j, k) * vm(i, j, k))
         enddo
      enddo
   enddo
   endfunction mass_3d

   real(R8P) function mass_3d_b(x, phi, pm, vm) result(tot)
   !< Composite mass of 3-D I (Kahan): the major cells with centre at r >= 2, all minor cells.
   real(R8P), intent(in) :: x(:)                 !< Major grid.
   real(R8P), intent(in) :: phi(:,:,:)           !< Major values.
   real(R8P), intent(in) :: pm(:,:,:), vm(:,:,:) !< Minor values, volumes.
   real(R8P)             :: comp, c(3)           !< Work.
   integer(I4P)          :: i, j, k              !< Counters.

   tot = 0._R8P ; comp = 0._R8P
   do k = 1, size(phi, 3, kind=I4P)
      do j = 1, size(phi, 2, kind=I4P)
         do i = 1, size(phi, 1, kind=I4P)
            c = 0.5_R8P * [x(i) + x(i + 1), x(j) + x(j + 1), x(k) + x(k + 1)]
            if (norm2(c) >= 2._R8P) call kahan(tot, comp, phi(i, j, k) * (x(2) - x(1)) ** 3)
         enddo
      enddo
   enddo
   do k = 1, size(pm, 3, kind=I4P)
      do j = 1, size(pm, 2, kind=I4P)
         do i = 1, size(pm, 1, kind=I4P)
            call kahan(tot, comp, pm(i, j, k) * vm(i, j, k))
         enddo
      enddo
   enddo
   endfunction mass_3d_b

   subroutine weights_3d(blk, x, ok)
   !< 3-D B: every cell of the shell a receptor over the donor grid x^3 (`overlap_weights_3d`): weights sum to 1, constant
   !< and linear fields (Picard's gradient term) exact, and per donor cell the sum of its intersections over all receptors
   !< = |V_c' n the polyhedral shell| (the cut of the donor grid by the shell's inner and outer faces, `overset_patch_3d`
   !< with the explicit faces: the same splits, so the receptors tile that polyhedron).
   type(curvilinear_block_t), intent(in)    :: blk  !< Shell.
   real(R8P),                 intent(in)    :: x(:) !< Donor grid (all directions).
   logical,                   intent(inout) :: ok   !< Checks.
   real(R8P), parameter      :: C0 = 3.7_R8P, B(3) = [0.7_R8P, -2.1_R8P, 1.3_R8P] !< Fields.
   type(overlap_weights_t)   :: ow                  !< Weights.
   type(overset_patch_3d_t)  :: pt                  !< The shell as a minor region.
   integer(I4P), allocatable :: cells(:,:), faces(:,:) !< Receptors, the shell's faces.
   real(R8P),    allocatable :: acc(:,:,:), cmp(:,:,:) !< Per donor cell sums.
   real(R8P)    :: ws, wc, cs, cc, qs, qc, xc(3), e(4), h, ref, dv, tot, comp, vsh, vol, diam !< Work.
   integer(I4P) :: ni, nj, nk, i, j, k, n, r, l, st, nx, d(3), nf !< Work.
   integer(I8P) :: t0, t1, tr                       !< Clock.

   ni = blk%ni - 1 ; nj = blk%nj - 1 ; nk = blk%nk - 1 ; nx = size(x, kind=I4P) - 1 ; h = x(2) - x(1)
   allocate(cells(3, ni * nj * nk), faces(5, 2 * ni * nj)) ; n = 0 ; nf = 0
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            n = n + 1 ; cells(:, n) = [i, j, k]
            if (k == 1) then
               nf = nf + 1 ; faces(:, nf) = [1_I4P, i, j, k, 5_I4P]
            endif
            if (k == nk) then
               nf = nf + 1 ; faces(:, nf) = [1_I4P, i, j, k, 6_I4P]
            endif
         enddo
      enddo
   enddo
   call system_clock(t0, tr)
   call overlap_weights_3d(blk, cells, x, x, x, ow, status=st)
   call system_clock(t1)
   if (st /= OVERLAP_WEIGHTS_OK .or. any(ow%flag /= OVERLAP_WEIGHTS_OK) .or. any(ow%covered /= 1._R8P)) &
      error stop 'weights_3d: overlap_weights_3d'
   e = 0._R8P
   do r = 1, ow%nr
      ws = 0._R8P ; wc = 0._R8P ; cs = 0._R8P ; cc = 0._R8P ; qs = 0._R8P ; qc = 0._R8P
      do l = ow%start(r), ow%start(r + 1) - 1
         d = ow%donor(:, l)
         xc = 0.5_R8P * [x(d(1)) + x(d(1) + 1), x(d(2)) + x(d(2) + 1), x(d(3)) + x(d(3) + 1)]
         call kahan(ws, wc, ow%weight(l)) ; call kahan(cs, cc, ow%weight(l) * C0)
         call kahan(qs, qc, ow%weight(l) * (C0 + dot_product(B, xc) + dot_product(B, ow%centroid(:, l) - xc)))
      enddo
      e(1) = max(e(1), abs(ws - 1._R8P)) ; e(2) = max(e(2), abs(cs - C0) / C0)
      ! the linear field at the receptor's centroid (its polyhedron rebuilt here)
      call receptor_polyhedron(blk, cells(:, r), vol, xc, diam)
      ref = C0 + dot_product(B, xc)
      e(3) = max(e(3), abs(qs - ref) / (norm2(B) * diam + abs(ref)))
   enddo
   ! per donor cell: the sum over the receptors against the cut of the polyhedral shell
   call overset_patch_3d([blk], [1_I4P], x, x, x, pt, faces=faces(:, 1:nf), status=st)
   if (st /= OVERSET_PATCH_OK) error stop 'weights_3d: overset_patch_3d of the shell'
   allocate(acc(nx, nx, nx), cmp(nx, nx, nx)) ; acc = 0._R8P ; cmp = 0._R8P
   do l = 1, ow%nw
      d = ow%donor(:, l) ; call kahan(acc(d(1), d(2), d(3)), cmp(d(1), d(2), d(3)), ow%volume(l))
   enddo
   tot = 0._R8P ; comp = 0._R8P
   do k = 1, nx
      do j = 1, nx
         do i = 1, nx
            dv = (1._R8P - pt%cut%alpha(i, j, k)) * h ** 3
            e(4) = max(e(4), abs(acc(i, j, k) - dv) / h ** 3)
            call kahan(tot, comp, acc(i, j, k))
         enddo
      enddo
   enddo
   vsh = pt%boundary%get_volume()
   print '(A,I0,A,I0,A,I0,A,F6.2,A,F7.3,A)', '    B: ', nx, '^3, ', ow%nr, ' receptors (the whole shell), ', ow%nw, &
         ' intersections (', real(ow%nw, R8P) / ow%nr, ' each); overlap_weights_3d ', real(t1 - t0, R8P) / tr, ' s'
   print '(A,4ES10.2,A,ES10.2)', '       |sum w - 1|, constant, linear, per donor |sum - |V_c'' n shell|| / h^3: ', e, &
         '; total vs shell volume ', abs(tot - vsh) / vsh
   ok = ok .and. e(1) <= 2 * epsilon(1._R8P) .and. e(2) <= 4 * epsilon(1._R8P) .and. e(3) <= 1e-13_R8P .and. &
        e(4) <= 1e-13_R8P .and. abs(tot - vsh) <= 1e-13_R8P * vsh
   endsubroutine weights_3d

   subroutine receptor_polyhedron(blk, c, vol, cen, diam)
   !< Volume, centroid and diameter of receptor cell c: 12 triangles, faces out of the cell by the block's handedness,
   !< shorter diagonal, ties by the least global vertex index (built here, independently of the library).
   type(curvilinear_block_t), intent(in)  :: blk  !< Block.
   integer(I4P),              intent(in)  :: c(3) !< Cell.
   real(R8P),                 intent(out) :: vol, cen(3), diam !< Volume, centroid, diameter.
   integer(I4P), parameter :: FQ(4, 6) = reshape([1, 5, 7, 3,  2, 4, 8, 6,  1, 2, 6, 5,  3, 7, 8, 4,  1, 3, 4, 2, &
                                                  5, 6, 8, 7], [4, 6]) !< Faces out of a right-handed cell.
   real(R8P)    :: v(3, 8), q(3, 4), t(3, 3, 2), d13, d24, a(3), b(3), e(3), dv, o(3) !< Work.
   integer(I4P) :: g(8), gq(4), m, f, ii, jj, kk               !< Work.
   logical      :: use13                                       !< Diagonal.

   do m = 1, 8
      ii = c(1) + mod(m - 1, 2) ; jj = c(2) + mod((m - 1) / 2, 2) ; kk = c(3) + (m - 1) / 4
      v(:, m) = vtx(blk, [ii, jj, kk]) ; g(m) = ii + blk%ni * ((jj - 1) + blk%nj * (kk - 1))
   enddo
   o = v(:, 1) ; vol = 0._R8P ; cen = 0._R8P ; diam = norm2(maxval(v, dim=2) - minval(v, dim=2))
   do f = 1, 6
      q = v(:, FQ(:, f)) ; gq = g(FQ(:, f))
      if (blk%orientation() < 0._R8P) then
         q = q(:, [1, 4, 3, 2]) ; gq = gq([1, 4, 3, 2])
      endif
      d13 = norm2(q(:, 3) - q(:, 1)) ; d24 = norm2(q(:, 4) - q(:, 2))
      use13 = d13 < d24 .or. (d13 == d24 .and. min(gq(1), gq(3)) <= min(gq(2), gq(4)))
      if (use13) then
         t(:, :, 1) = q(:, [1, 2, 3]) ; t(:, :, 2) = q(:, [1, 3, 4])
      else
         t(:, :, 1) = q(:, [1, 2, 4]) ; t(:, :, 2) = q(:, [2, 3, 4])
      endif
      do m = 1, 2
         a = t(:, 1, m) - o ; b = t(:, 2, m) - o ; e = t(:, 3, m) - o
         dv = dot_product(a, cross3(b, e)) / 6._R8P
         vol = vol + dv ; cen = cen + dv * (a + b + e) / 4._R8P
      enddo
   enddo
   cen = o + cen / vol
   endsubroutine receptor_polyhedron

   pure real(R8P) function blob3(p)
   !< The initial scalar in 3-D: a Gaussian upstream of the minor shell.
   real(R8P), intent(in) :: p(3) !< Point.

   blob3 = exp(-((p(1) + 2.2_R8P) ** 2 + (p(2) - 0.4_R8P) ** 2 + (p(3) - 0.2_R8P) ** 2) / 0.35_R8P ** 2)
   endfunction blob3

   pure function vel3(p) result(u)
   !< The potential flow past the unit sphere at the origin, free stream U3 (unit).
   real(R8P), intent(in) :: p(3) !< Point.
   real(R8P)             :: u(3) !< Velocity.
   real(R8P)             :: r    !< Radius.

   r = norm2(p)
   u = U3 + 0.5_R8P * (U3 / r ** 3 - 3._R8P * dot_product(U3, p) * p / r ** 5)
   endfunction vel3

   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   subroutine shell(na, nl, blk)
   !< A sphere shell, r in [1, 2] over nl layers, na x na / 2 cells on each sphere (i around, j pole to pole, k outward:
   !< left-handed), poles collapsed; its outer face (+k) named 'overlap'.
   integer(I4P),              intent(in)    :: na, nl !< Cells around, radially.
   type(curvilinear_block_t), intent(inout) :: blk    !< Block.
   type(topo_record_t)    :: rec                      !< Record.
   real(R8P), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:) !< Vertices.
   real(R8P)              :: th, ph, rr               !< Work.
   integer(I4P)           :: np, i, j, k              !< Counters.

   np = na / 2
   allocate(xx(na + 1, np + 1, nl + 1), yy(na + 1, np + 1, nl + 1), zz(na + 1, np + 1, nl + 1))
   do k = 1, nl + 1
      rr = 1._R8P + (k - 1._R8P) / nl
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            xx(i, j, k) = rr * sin(ph) * cos(th) ; yy(i, j, k) = rr * sin(ph) * sin(th) ; zz(i, j, k) = rr * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            xx(:, j, k) = 0._R8P ; yy(:, j, k) = 0._R8P ; zz(:, j, k) = merge(rr, -rr, j == 1)
         endif
      enddo
      xx(na + 1, :, k) = xx(1, :, k) ; yy(na + 1, :, k) = yy(1, :, k) ; zz(na + 1, :, k) = zz(1, :, k)
   enddo
   call blk%destroy
   call blk%initialize(xx, yy, zz)
   rec%kind = TOPO_BOUNDARY ; rec%face = 6 ; rec%lo = [1, 1, nl + 1] ; rec%hi = [na + 1, np + 1, nl + 1]
   rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine shell
endprogram fossil_val_overset_conservative
