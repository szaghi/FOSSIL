!< FOSSIL test: flux-redistribution neighbourhoods of the cut-cell pieces (issue #70, step 6).
program fossil_test_cut_frd
!< FOSSIL test: flux-redistribution neighbourhoods of the cut-cell pieces (issue #70, step 6).
!<
!< `frd_neighbourhoods` on 2-D (a circle of 64 segments, a thin plate) and 3-D cut cells (the 5-degree wedge tilted off
!< every axis, the thin plate):
!< 1. the set of every piece equals a brute-force fixpoint (a pull formulation, independent of the push search): a piece
!<    of the 3^d block with volume joins when it is linked to a member whose offset it extends by one step away from
!<    the centre on one axis;
!< 2. the weights are kappa_j / sum kappa and sum to 1 (1e-14; 8.9e-16 measured); a piece filling its cell has itself
!<    only;
!< 3. no set crosses a plate, for the pieces whose 3^d block lies within the plate's extent (near its ends a monotone
!<    path goes around the end).

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, frd_t, &
                                     frd_neighbourhoods, pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c
type(cut_cells_2d_t)      :: c2
type(cut_cells_3d_t)      :: c3
type(cut_links_t)         :: lk
type(frd_t)               :: frd
real(R8P), allocatable    :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                 :: emax
integer(I4P)              :: i, n, nbad, ncut, nmulti
logical                   :: ok

emax = 0._R8P ; nbad = 0 ; ncut = 0 ; nmulti = 0
call bench_bodies(b)
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call run2(x, y, 0._R8P)
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P], v=[0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P], loop_sizes=[4_I4P])
x = [(0.05_R8P + 0.1_R8P * (i - 1), i = 1, 10)]
call run2(x, x, 0.5_R8P)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call run3(x, y, z, 0._R8P)
call to_surface(b(5)%soup, s) ! plate [0, 1]^2 x [0, 1e-3]
x = [(-0.15_R8P + 0.1_R8P * i, i = 0, 13)] ; z = [-0.3_R8P, -0.1_R8P, 0.1_R8P, 0.3_R8P]
call run3(x, x, z, 5.e-4_R8P)
print '(A,I0,A,I0,A,I0)', 'cut pieces ', ncut, '  sets with 2+ members ', nmulti, '  bad ', nbad
print '(A,ES10.2)', 'max |sum of weights - 1|: ', emax
ok = nbad == 0 .and. ncut > 100 .and. nmulti > 100 .and. emax <= 1.e-14_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine run2(xx, yy, plate)
   !< 2-D cut cells, links, sets and checks.
   real(R8P), intent(in) :: xx(:), yy(:), plate

   call c%cut_cells(x=xx, y=yy, cut=c2)
   call cut_links_2d(c2, lk)
   call frd_neighbourhoods(lk, c2%piece_cell, c2%piece_area, xx, yy, frd)
   call checks(c2%piece_cell, c2%piece_area, c2%piece_centroid, xx, yy, plate)
   endsubroutine run2

   subroutine run3(xx, yy, zz, plate)
   !< 3-D cut cells, links, sets and checks.
   real(R8P), intent(in) :: xx(:), yy(:), zz(:), plate

   call cut_cells_3d(s, xx, yy, zz, c3, links=lk)
   call frd_neighbourhoods(lk, c3%piece_cell, c3%piece_volume, xx, yy, frd, z=zz)
   call checks(c3%piece_cell, c3%piece_volume, c3%piece_centroid, xx, yy, plate, zz)
   endsubroutine run3

   subroutine checks(pc, pv, pcen, xx, yy, plate, zz)
   !< Checks 1-3.
   integer(I4P), intent(in)           :: pc(:,:)
   real(R8P),    intent(in)           :: pv(:), pcen(:,:), xx(:), yy(:), plate
   real(R8P),    intent(in), optional :: zz(:)
   logical,      allocatable          :: inr(:), inf(:)
   real(R8P),    allocatable          :: kap(:)
   integer(I4P), allocatable          :: cs(:), cl(:), blk(:)
   integer(I4P)                       :: p, q, k, d, np, nb(3), nc, ic, ib, i1, i2, i3, jc(3), m, nbk
   real(R8P)                          :: cv, ks
   logical                            :: changed

   d = size(pc, 1) ; np = size(pv)
   allocate(inr(np), inf(np), kap(np))
   nb = 1 ; nb(1:d) = maxval(pc, 2) ; nc = product(nb) ! pieces by cell (CSR)
   allocate(cs(nc + 1), cl(np), blk(27 * 8)) ; cs = 0
   do p = 1, np
      ic = cell_id(pc(:, p), nb) ; cs(ic + 1) = cs(ic + 1) + 1
   enddo
   cs(1) = 1
   do ic = 1, nc
      cs(ic + 1) = cs(ic + 1) + cs(ic)
   enddo
   block
      integer(I4P) :: fill(nc)
      fill = 0
      do p = 1, np
         ic = cell_id(pc(:, p), nb) ; cl(cs(ic) + fill(ic)) = p ; fill(ic) = fill(ic) + 1
      enddo
   endblock
   do p = 1, np
      cv = (xx(pc(1, p) + 1) - xx(pc(1, p))) * (yy(pc(2, p) + 1) - yy(pc(2, p)))
      if (present(zz)) cv = cv * (zz(pc(3, p) + 1) - zz(pc(3, p)))
      kap(p) = pv(p) / cv
   enddo
   do p = 1, np
      if (pv(p) <= 0._R8P) then
         if (frd%start(p + 1) /= frd%start(p)) nbad = nbad + 1
         cycle
      endif
      if (frd%member(frd%start(p)) /= p) nbad = nbad + 1
      inf = .false.
      do k = frd%start(p), frd%start(p + 1) - 1
         inf(frd%member(k)) = .true.
      enddo
      if (kap(p) >= 1._R8P) then
         if (frd%start(p + 1) - frd%start(p) /= 1) nbad = nbad + 1
         cycle
      endif
      ncut = ncut + 1
      if (frd%start(p + 1) - frd%start(p) > 1) nmulti = nmulti + 1
      nbk = 0 ! the pieces of the 3^d block
      do i3 = -1, merge(1, 0, d == 3)
         do i2 = -1, 1
            do i1 = -1, 1
               jc = 1 ; jc(1:d) = pc(:, p) ; jc(1) = jc(1) + i1 ; jc(2) = jc(2) + i2 ; if (d == 3) jc(3) = jc(3) + i3
               if (any(jc < 1) .or. any(jc > nb)) cycle
               ic = cell_id(jc(1:d), nb)
               do m = cs(ic), cs(ic + 1) - 1
                  if (nbk == size(blk)) blk = [blk, blk]
                  nbk = nbk + 1 ; blk(nbk) = cl(m)
               enddo
            enddo
         enddo
      enddo
      inr = .false. ; inr(p) = .true. ; changed = .true. ! brute force: pull until nothing changes
      do while (changed)
         changed = .false.
         do ib = 1, nbk
            q = blk(ib)
            if (inr(q) .or. pv(q) <= 0._R8P) cycle
            do k = lk%start(q), lk%start(q + 1) - 1
               if (.not. inr(lk%piece(k))) cycle
               if (extends(pc(:, q) - pc(:, p), pc(:, lk%piece(k)) - pc(:, p))) then
                  inr(q) = .true. ; changed = .true. ; exit
               endif
            enddo
         enddo
      enddo
      if (any(inr .neqv. inf)) nbad = nbad + 1
      ks = sum(kap, mask=inr)
      do k = frd%start(p), frd%start(p + 1) - 1
         if (abs(frd%weight(k) - kap(frd%member(k)) / ks) > 1.e-15_R8P) nbad = nbad + 1
         if (plate > 0._R8P .and. pv(frd%member(k)) > 0._R8P .and. inner(pc(:, p))) then
            if ((pcen(d, frd%member(k)) - plate) * (pcen(d, p) - plate) < 0._R8P) nbad = nbad + 1
         endif
      enddo
      emax = max(emax, abs(sum(frd%weight(frd%start(p):frd%start(p + 1) - 1)) - 1._R8P))
   enddo
   endsubroutine checks

   logical function inner(cc)
   !< The 3^d block of cell cc lies within the plate's extent along the axes other than the last ([0.2, 0.8] in 2-D,
   !< [0, 1]^2 in 3-D).
   integer(I4P), intent(in) :: cc(:)
   integer(I4P)             :: a
   real(R8P)                :: lo, hi

   inner = .true.
   do a = 1, size(cc) - 1
      if (size(cc) == 2) then
         lo = x(max(cc(a) - 1, 1)) ; hi = x(min(cc(a) + 2, size(x))) ; inner = inner .and. lo >= 0.2_R8P .and. hi <= 0.8_R8P
      else
         lo = x(max(cc(a) - 1, 1)) ; hi = x(min(cc(a) + 2, size(x))) ; inner = inner .and. lo >= 0._R8P .and. hi <= 1._R8P
      endif
   enddo
   endfunction inner

   pure integer(I4P) function cell_id(cc, nb)
   !< Linear index of cell cc in a block of nb cells.
   integer(I4P), intent(in) :: cc(:), nb(3)

   cell_id = cc(1) + nb(1) * (cc(2) - 1)
   if (size(cc) == 3) cell_id = cell_id + nb(1) * nb(2) * (cc(3) - 1)
   endfunction cell_id

   pure logical function extends(oq, orr) result(yes)
   !< Offset oq is offset orr plus one step away from the centre on one axis.
   integer(I4P), intent(in) :: oq(:), orr(:)
   integer(I4P)             :: a, nstep

   yes = .true. ; nstep = 0
   do a = 1, size(oq)
      if (oq(a) == orr(a)) cycle
      nstep = nstep + 1
      if (orr(a) /= 0 .or. abs(oq(a)) /= 1) yes = .false.
   enddo
   yes = yes .and. nstep == 1
   endfunction extends
endprogram fossil_test_cut_frd
