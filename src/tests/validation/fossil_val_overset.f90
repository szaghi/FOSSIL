!< FOSSIL validation: overset assembly at size (issue #72, step 10).
program fossil_val_overset
!< FOSSIL validation: overset assembly at size (issue #72, step 10).
!<
!< 1. Twilight zone through an assembly: a cylinder O-grid (1 <= r <= 2, overlap boundary outside) in a Cartesian
!<    background cut by the cylinder, refined three times (O-grid 32 x 12 to 128 x 48 cells, background h = 0.15 to
!<    0.0375); explicit assembly with stencils of width q = 2, 3, 4 (order-matched xi); every receptor interpolates
!<    f = sin(2x + 1) cos(3y) and the velocity of an isentropic vortex centred in the overlap from its donor stencil:
!<    the slope of the max error over the refinements is q within 0.3 (fixed overlap: q >= s + 1 suffices, Chesshire &
!<    Henshaw); no orphan at any level.
!< 2. 3-D: a spherical-shell O-grid (64 x 32 x 16 cells, poles collapsed) in a 49^3 background: holes, donor
!<    candidates, classification and fringes, timed; no orphan.
!< Exit status 0 if every check holds.

use fossil, only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, &
                   surface_stl_object, wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, &
                   donor_set_t, donor_candidates, donor_candidates_2d, overset_mask_t, classify, classify_2d, &
                   fringe_levels, SELECT_EXPLICIT
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
real(R8P)            :: e(3, 2), rate(2)
integer(I4P)         :: q, l, norph
logical              :: ok

ok = .true.
do q = 2, 4
   do l = 1, 3
      call twilight(2 ** (l - 1), q, e(l, :), norph)
      ok = ok .and. norph == 0
   enddo
   rate = log(e(1, :) / e(3, :)) / log(4._R8P)
   print '(A,I0,A,3ES10.2,A,F6.2,A,3ES10.2,A,F6.2)', 'q = ', q, ': f errors ', e(:, 1), ' order ', rate(1), &
         '; vortex errors ', e(:, 2), ' order ', rate(2)
   ok = ok .and. all(rate >= q - 0.3_R8P)
enddo
call sphere_assembly(ok)
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   pure real(R8P) function f1(x, y) result(f)
   !< Smooth test function.
   real(R8P), intent(in) :: x, y !< Point.

   f = sin(2 * x + 1) * cos(3 * y)
   endfunction f1

   pure real(R8P) function f2(x, y) result(f)
   !< u velocity of an isentropic vortex of strength 0.8 centred at (1.5, 0.4), core radius 0.4.
   real(R8P), intent(in) :: x, y !< Point.
   real(R8P)             :: r2   !< Squared distance over the core.

   r2 = ((x - 1.5_R8P) ** 2 + (y - 0.4_R8P) ** 2) / 0.16_R8P
   f = -0.8_R8P / (2 * PI) * (y - 0.4_R8P) / 0.4_R8P * exp(0.5_R8P * (1 - r2))
   endfunction f2

   subroutine twilight(m, q, err, norph)
   !< Max interpolation errors (f1, f2) at the receptors of the m-times refined assembly.
   integer(I4P), intent(in)  :: m, q     !< Refinement, width.
   real(R8P),    intent(out) :: err(2)   !< Errors.
   integer(I4P), intent(out) :: norph    !< Orphans.
   type(curvilinear_block_2d_t)       :: blocks(2) !< Blocks.
   type(topo_record_t)                :: rec       !< Records.
   type(curve_object)                 :: wall      !< Cutter.
   type(block_logical_t), allocatable :: holes(:)  !< Holes.
   type(donor_set_t),     allocatable :: donors(:) !< Candidates.
   type(overset_mask_t),  allocatable :: masks(:)  !< Classes.
   real(R8P),    allocatable :: u(:,:), v(:,:), x(:,:), y(:,:), xi(:,:), w(:,:), ru(:), rv(:) !< Data.
   integer(I4P), allocatable :: cell(:,:), nodes(:,:,:), sst(:) !< Stencils.
   integer(I4P)              :: ni, nj, nb, i, j, a, k, c, n, kk !< Counters.
   real(R8P)                 :: r, s1, s2 !< Work.

   ni = 32 * m + 1 ; nj = 12 * m + 1 ; nb = 40 * m + 1
   allocate(u(ni, nj), v(ni, nj), x(nb, nb), y(nb, nb))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         u(i, j) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = r * sin(2 * PI * (i - 1) / (ni - 1))
      enddo
   enddo
   u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   call blocks(1)%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   do j = 1, nb
      do i = 1, nb
         x(i, j) = -3._R8P + 6._R8P * (i - 1) / (nb - 1) ; y(i, j) = -3._R8P + 6._R8P * (j - 1) / (nb - 1)
      enddo
   enddo
   call blocks(2)%initialize(x, y)
   call blocks(1)%build_search ; call blocks(2)%build_search
   call wall_cutter_2d(blocks(1:1), wall)
   call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
   call donor_candidates_2d(blocks, donors)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, q, holes=holes, layers=(q + 1) / 2)
   norph = count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0)
   err = 0._R8P
   do a = 1, 2
      n = count(masks(a)%iota < 0)
      if (n == 0) cycle
      allocate(cell(2, n), xi(2, n), ru(n), rv(n)) ; k = 0
      do i = 1, size(masks(a)%iota)
         if (masks(a)%iota(i) >= 0) cycle
         k = k + 1 ; c = masks(a)%donor(i) ; cell(:, k) = donors(a)%cell(:, c) ; xi(:, k) = donors(a)%xi(:, c)
         if (a == 1) then
            ru(k) = u(mod(i - 1, ni) + 1, (i - 1) / ni + 1) ; rv(k) = v(mod(i - 1, ni) + 1, (i - 1) / ni + 1)
         else
            ru(k) = x(mod(i - 1, nb) + 1, (i - 1) / nb + 1) ; rv(k) = y(mod(i - 1, nb) + 1, (i - 1) / nb + 1)
         endif
      enddo
      call blocks(3 - a)%stencils(cell, xi, q, nodes, w, sst)
      do k = 1, n
         s1 = 0._R8P ; s2 = 0._R8P
         do kk = 1, size(w, 1)
            if (w(kk, k) == 0._R8P) cycle
            associate(db => blocks(3 - a))
               s1 = s1 + w(kk, k) * f1(db%u(nodes(1, kk, k), nodes(2, kk, k)), db%v(nodes(1, kk, k), nodes(2, kk, k)))
               s2 = s2 + w(kk, k) * f2(db%u(nodes(1, kk, k), nodes(2, kk, k)), db%v(nodes(1, kk, k), nodes(2, kk, k)))
            endassociate
         enddo
         err = max(err, [abs(s1 - f1(ru(k), rv(k))), abs(s2 - f2(ru(k), rv(k)))])
      enddo
      deallocate(cell, xi, ru, rv)
   enddo
   endsubroutine twilight

   subroutine sphere_assembly(ok)
   !< The 3-D assembly of a spherical shell in a background, timed.
   logical, intent(inout) :: ok !< Checks.
   type(curvilinear_block_t)          :: blocks(2) !< Blocks.
   type(topo_record_t)                :: rec       !< Records.
   type(surface_stl_object)           :: wall      !< Cutter.
   type(block_logical_t), allocatable :: holes(:)  !< Holes.
   type(donor_set_t),     allocatable :: donors(:) !< Candidates.
   type(overset_mask_t),  allocatable :: masks(:)  !< Classes.
   real(R8P), allocatable    :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)              :: i, j, k, na, np, nr, nb !< Counters.
   integer(I8P)              :: t0, t1, t2, t3, cr !< Clock.
   real(R8P)                 :: r, th, ph   !< Coordinates.

   na = 64 ; np = 32 ; nr = 16 ; nb = 49
   allocate(x(na + 1, np + 1, nr + 1), y(na + 1, np + 1, nr + 1), z(na + 1, np + 1, nr + 1))
   do k = 1, nr + 1
      r = 1._R8P + (k - 1._R8P) / nr
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            x(:, j, k) = 0._R8P ; y(:, j, k) = 0._R8P ; z(:, j, k) = merge(r, -r, j == 1)
         endif
      enddo
   enddo
   x(na + 1,:,:) = x(1,:,:) ; y(na + 1,:,:) = y(1,:,:) ; z(na + 1,:,:) = z(1,:,:)
   call blocks(1)%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [na + 1, np + 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 6 ; rec%lo = [1, 1, nr + 1] ; rec%hi = [na + 1, np + 1, nr + 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   deallocate(x, y, z) ; allocate(x(nb, nb, nb), y(nb, nb, nb), z(nb, nb, nb))
   do k = 1, nb
      do j = 1, nb
         do i = 1, nb
            x(i, j, k) = -2.6_R8P + 5.2_R8P * (i - 1) / (nb - 1) ; y(i, j, k) = -2.6_R8P + 5.2_R8P * (j - 1) / (nb - 1)
            z(i, j, k) = -2.6_R8P + 5.2_R8P * (k - 1) / (nb - 1)
         enddo
      enddo
   enddo
   call blocks(2)%initialize(x, y, z)
   call system_clock(t0, cr)
   call blocks(1)%build_search ; call blocks(2)%build_search
   call wall_cutter(blocks(1:1), wall)
   call hole_flags(blocks, wall, holes, cuts=[.false., .true.])
   call system_clock(t1)
   call donor_candidates(blocks, donors)
   call system_clock(t2)
   call classify(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, layers=2_I4P)
   call fringe_levels(blocks, masks)
   call system_clock(t3)
   print '(A,I0,A,I0,A,I0,A,I0)', 'sphere: holes ', count(masks(2)%iota == 0), ', receptors ', &
         count(masks(1)%iota < 0) + count(masks(2)%iota < 0), ', orphans ', &
         count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0), ', candidates ', &
         size(donors(1)%block) + size(donors(2)%block)
   print '(A,3F9.3,A)', 'sphere: search + holes, donor candidates, classification [s]: ', real(t1 - t0, R8P) / cr, &
         real(t2 - t1, R8P) / cr, real(t3 - t2, R8P) / cr
   ok = ok .and. count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0) == 0 .and. count(masks(2)%iota == 0) > 0
   endsubroutine sphere_assembly
endprogram fossil_val_overset
