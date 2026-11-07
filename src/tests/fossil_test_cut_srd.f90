!< FOSSIL test: state-redistribution neighbourhoods of the cut-cell pieces (issue #70, step 5).
program fossil_test_cut_srd
!< FOSSIL test: state-redistribution neighbourhoods of the cut-cell pieces (issue #70, step 5).
!<
!< `srd_neighbourhoods` on 2-D (a circle of 64 segments, a thin plate, a half-plane) and 3-D cut cells (the 5-degree
!< wedge tilted off every axis, the cone, the thin plate), tau = 1/2:
!< 1. conservation (Giuliani et al. 2022 §3.1, unit column sums): every piece with volume is in exactly N^ov_j
!<    neighbourhoods, so sum_p sum_{j in N_p} 1/N^ov_j = 1 per j (exact), and sum_p vhat_p is the total fluid volume
!<    (1e-13; 1e-14 measured);
!< 2. every complete neighbourhood reaches tau |cell|; none is incomplete on these bodies; members are distinct, have
!<    volume, lie in the 5^d tile and the piece comes first;
!< 3. xhat_p is the weighted centroid of the members, recomputed (1e-14 of the cell size);
!< 4. weighted factors: 0 <= eta^nb <= 1 for a complete neighbourhood, 0 < eta^self <= 1;
!< 5. the half-plane x < 0.58 on cells of 0.1: each small cell (fraction 0.2) merges with the cell across its +x face
!<    only (the wall normal); the plates: no neighbourhood crosses the plate.

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, &
                                     srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, SRD_NO_VOLUME, &
                                     pose_from_axis_angle
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
type(srd_t)               :: srd
real(R8P), allocatable    :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                 :: emax(3)
integer(I4P)              :: i, n, nbad, nsmall
logical                   :: ok

emax = 0._R8P ; nbad = 0 ; nsmall = 0
call bench_bodies(b)
! 2-D
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call run2(x, y, 0._R8P)
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P], v=[0.4995_R8P, 0.4995_R8P, 0.5005_R8P, 0.5005_R8P], loop_sizes=[4_I4P])
x = [(0.05_R8P + 0.1_R8P * (i - 1), i = 1, 10)]
call run2(x, x, 0.5_R8P)
call c%adopt_loops(u=[-1._R8P, 0.58_R8P, 0.58_R8P, -1._R8P], v=[-1._R8P, -1._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
x = [(0.1_R8P * (i - 1), i = 1, 11)]
call run2(x, x, 0._R8P)
call half_plane
! 3-D
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call run3(x, y, z, 0._R8P)
call to_surface(b(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call run3(x, y, z, 0._R8P)
call to_surface(b(5)%soup, s) ! plate [0, 1]^2 x [0, 1e-3]
x = [(-0.15_R8P + 0.1_R8P * i, i = 0, 13)] ; z = [-0.3_R8P, -0.1_R8P, 0.1_R8P, 0.3_R8P]
call run3(x, x, z, 5.e-4_R8P)
print '(A,I0,A,I0)', 'small pieces ', nsmall, '  bad ', nbad
print '(A,3ES10.2)', 'max conservation / total volume / centroid errors: ', emax
ok = nbad == 0 .and. nsmall > 50 .and. emax(1) <= 1.e-14_R8P .and. emax(2) <= 1.e-13_R8P .and. emax(3) <= 1.e-14_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine run2(xx, yy, plate)
   !< 2-D cut cells, links, neighbourhoods and checks 1-4 (and the plate check if plate > 0).
   real(R8P), intent(in) :: xx(:), yy(:), plate

   call c%cut_cells(x=xx, y=yy, cut=c2)
   call cut_links_2d(c2, lk)
   call srd_neighbourhoods(lk, c2%piece_cell, c2%piece_area, c2%piece_centroid, c2%piece_wet_normal, xx, yy, srd)
   call checks(c2%piece_cell, c2%piece_area, c2%piece_centroid, xx, yy, plate)
   endsubroutine run2

   subroutine run3(xx, yy, zz, plate)
   !< 3-D cut cells, links, neighbourhoods and checks 1-4 (and the plate check if plate > 0).
   real(R8P), intent(in) :: xx(:), yy(:), zz(:), plate

   call cut_cells_3d(s, xx, yy, zz, c3, links=lk)
   call srd_neighbourhoods(lk, c3%piece_cell, c3%piece_volume, c3%piece_centroid, c3%piece_wet_normal, xx, yy, srd, z=zz)
   call checks(c3%piece_cell, c3%piece_volume, c3%piece_centroid, xx, yy, plate, zz)
   endsubroutine run3

   subroutine checks(pc, pv, pcen, xx, yy, plate, zz)
   !< Checks 1-4 on srd.
   integer(I4P), intent(in)           :: pc(:,:)
   real(R8P),    intent(in)           :: pv(:), pcen(:,:), xx(:), yy(:), plate
   real(R8P),    intent(in), optional :: zz(:)
   integer(I4P), allocatable          :: cnt(:)
   real(R8P),    allocatable          :: colsum(:)
   integer(I4P)                       :: p, k, j, k2, d
   real(R8P)                          :: cv, vm, xc(3), w, tot, h

   d = size(pc, 1)
   allocate(cnt(size(pv)), colsum(size(pv))) ; cnt = 0 ; colsum = 0._R8P
   do p = 1, srd%np
      if (srd%status(p) == SRD_NO_VOLUME) then
         if (srd%start(p + 1) /= srd%start(p)) nbad = nbad + 1
         cycle
      endif
      if (srd%member(srd%start(p)) /= p) nbad = nbad + 1
      cv = (xx(pc(1, p) + 1) - xx(pc(1, p))) * (yy(pc(2, p) + 1) - yy(pc(2, p)))
      h = xx(pc(1, p) + 1) - xx(pc(1, p))
      if (present(zz)) cv = cv * (zz(pc(3, p) + 1) - zz(pc(3, p)))
      if (srd%status(p) /= SRD_REGULAR) nsmall = nsmall + 1
      if (srd%status(p) == SRD_INCOMPLETE) nbad = nbad + 1
      vm = 0._R8P ; xc = 0._R8P ; w = 0._R8P
      do k = srd%start(p), srd%start(p + 1) - 1
         j = srd%member(k)
         cnt(j) = cnt(j) + 1
         if (pv(j) <= 0._R8P .or. maxval(abs(pc(:, j) - pc(:, p))) > 2) nbad = nbad + 1
         do k2 = srd%start(p), k - 1
            if (srd%member(k2) == j) nbad = nbad + 1
         enddo
         vm = vm + pv(j)
         xc(1:d) = xc(1:d) + pv(j) / srd%nov(j) * pcen(:, j) ; w = w + pv(j) / srd%nov(j)
         if (plate > 0._R8P .and. pv(p) > 0._R8P) then
            if ((pcen(d, j) - plate) * (pcen(d, p) - plate) < 0._R8P) nbad = nbad + 1
         endif
      enddo
      if (srd%status(p) == SRD_COMPLETE .and. vm < srd%tau * cv) nbad = nbad + 1
      if (w > 0._R8P) emax(3) = max(emax(3), maxval(abs(xc(1:d) / w - srd%xhat(:, p))) / h)
      if (srd%status(p) == SRD_COMPLETE .and. (srd%eta_nb(p) < 0._R8P .or. srd%eta_nb(p) > 1._R8P)) nbad = nbad + 1
      if (srd%eta_self(p) <= 0._R8P .or. srd%eta_self(p) > 1._R8P) nbad = nbad + 1
   enddo
   tot = 0._R8P
   do p = 1, srd%np
      if (pv(p) <= 0._R8P) cycle
      if (cnt(p) /= srd%nov(p)) nbad = nbad + 1
      tot = tot + pv(p)
   enddo
   do p = 1, srd%np ! unit column sums
      do k = srd%start(p), srd%start(p + 1) - 1
         j = srd%member(k) ; colsum(j) = colsum(j) + 1._R8P / srd%nov(j)
      enddo
   enddo
   do p = 1, srd%np
      if (pv(p) > 0._R8P) emax(1) = max(emax(1), abs(colsum(p) - 1._R8P))
   enddo
   emax(2) = max(emax(2), abs(sum(srd%vhat, mask=pv > 0._R8P) - tot) / tot)
   endsubroutine checks

   subroutine half_plane
   !< Check 5 on the half-plane x < 0.58: every small piece is {p, the piece across its +x face}.
   integer(I4P) :: p, q

   call c%cut_cells(x=x, y=x, cut=c2)
   call cut_links_2d(c2, lk)
   call srd_neighbourhoods(lk, c2%piece_cell, c2%piece_area, c2%piece_centroid, c2%piece_wet_normal, x, x, srd)
   n = 0
   do p = 1, srd%np
      if (srd%status(p) /= SRD_COMPLETE) cycle
      n = n + 1
      if (srd%start(p + 1) - srd%start(p) /= 2) then
         nbad = nbad + 1 ; cycle
      endif
      q = srd%member(srd%start(p) + 1)
      if (any(c2%piece_cell(:, q) - c2%piece_cell(:, p) /= [1, 0])) nbad = nbad + 1
   enddo
   if (n /= 10) nbad = nbad + 1 ! the column x in [0.5, 0.6], 10 cells
   endsubroutine half_plane
endprogram fossil_test_cut_srd
