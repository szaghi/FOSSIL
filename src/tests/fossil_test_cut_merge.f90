!< FOSSIL test: merge targets of the cut-cell pieces (issue #70, step 7).
program fossil_test_cut_merge
!< FOSSIL test: merge targets of the cut-cell pieces (issue #70, step 7).
!<
!< `merge_targets` on 2-D (the half-plane x < 0.58, a circle of 64 segments, a translating circle) and 3-D cut cells (the
!< 5-degree wedge tilted off every axis), tau = 1/2:
!< 1. the half-plane on cells of 0.1: each small cell (fraction 0.2) has the cell across its +x face as its single
!<    target and as its only normal-directed target, weight 1;
!< 2. every served piece has a target that is linked to it and has volume; it is suitable (not small, not excluded)
!<    unless the status says MERGE_CHAIN, and then final_target is; no MERGE_CYCLE or MERGE_ALONE on these bodies;
!< 3. normal-directed targets: at most 3 (2-D) or 7 (3-D), in distinct cells at offsets o_a in {0, sign n_a}, weights
!<    summing to 1 (1e-15);
!< 4. a moving body: the circle translated by (0.07, 0.03) over cells of 0.05; the fresh pieces (cells cleared by the
!<    step, at the end pose) and the dead ones (cells covered, at the start pose, cells covered excluded as targets)
!<    all get targets, the dead ones in cells that stay wet.

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, merge_t, &
                                     merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, swept_cells_2d_t, swept_cells_2d, &
                                     pose_2d, pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: c, c1
type(cut_cells_2d_t)      :: c2
type(cut_cells_3d_t)      :: c3
type(cut_links_t)         :: lk
type(merge_t)             :: mt
type(swept_cells_2d_t)    :: sw
real(R8P), allocatable    :: x(:), y(:), z(:), u(:), v(:)
logical,   allocatable    :: msk(:), exc(:)
real(R8P)                 :: emax
integer(I4P)              :: i, n, nbad, nserved, nfresh, ndead, p, q
logical                   :: ok

emax = 0._R8P ; nbad = 0 ; nserved = 0 ; nfresh = 0 ; ndead = 0
call bench_bodies(b)
! 1. half-plane
call c%adopt_loops(u=[-1._R8P, 0.58_R8P, 0.58_R8P, -1._R8P], v=[-1._R8P, -1._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
x = [(0.1_R8P * (i - 1), i = 1, 11)]
call c%cut_cells(x=x, y=x, cut=c2)
call cut_links_2d(c2, lk)
call merge_targets(lk, c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x, mt)
n = 0
do p = 1, mt%np
   if (mt%status(p) == MERGE_NONE) cycle
   n = n + 1 ; q = mt%target(p)
   if (q == 0) then
      nbad = nbad + 1 ; cycle
   endif
   if (any(c2%piece_cell(:, q) - c2%piece_cell(:, p) /= [1, 0])) nbad = nbad + 1
   if (mt%mix_start(p + 1) - mt%mix_start(p) /= 1) then
      nbad = nbad + 1
   elseif (mt%mix_piece(mt%mix_start(p)) /= q .or. abs(mt%mix_weight(mt%mix_start(p)) - 1._R8P) > 0._R8P) then
      nbad = nbad + 1
   endif
enddo
if (n /= 10) nbad = nbad + 1
! 2-3. circle, wedge
n = 64
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.013_R8P + 0.07_R8P * (i - 1), i = 1, 15)] ; y = [(-0.021_R8P + 0.071_R8P * (i - 1), i = 1, 16)]
call c%cut_cells(x=x, y=y, cut=c2)
call cut_links_2d(c2, lk)
call merge_targets(lk, c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, y, mt)
call checks(c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, y)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call cut_cells_3d(s, x, y, z, c3, links=lk)
call merge_targets(lk, c3%piece_cell, c3%piece_volume, c3%piece_wet_normal, x, y, mt, z=z)
call checks(c3%piece_cell, c3%piece_volume, c3%piece_wet_normal, x, y, z)
! 4. a moving circle
x = [(0.05_R8P * (i - 1), i = 1, 21)]
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P, [0.07_R8P, 0.03_R8P]), x, x, sw)
c1 = c ; call c1%move(pose_2d(0._R8P, [0.07_R8P, 0.03_R8P]))
call c1%cut_cells(x=x, y=x, cut=c2) ! fresh: at the end pose
call cut_links_2d(c2, lk)
allocate(msk(c2%np)) ; msk = .false.
do p = 1, c2%np
   msk(p) = sw%cleared(c2%piece_cell(1, p), c2%piece_cell(2, p))
enddo
nfresh = count(msk)
call merge_targets(lk, c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x, mt, mask=msk)
call checks(c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x)
call c%cut_cells(x=x, y=x, cut=c2) ! dead: at the start pose, the covered cells excluded
call cut_links_2d(c2, lk)
deallocate(msk) ; allocate(msk(c2%np), exc(c2%np))
do p = 1, c2%np
   msk(p) = sw%covered(c2%piece_cell(1, p), c2%piece_cell(2, p)) ; exc(p) = msk(p)
enddo
ndead = count(msk)
call merge_targets(lk, c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x, mt, mask=msk, exclude=exc)
call checks(c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x, excl=exc)
print '(A,I0,A,I0,A,I0,A,I0)', 'served ', nserved, '  fresh ', nfresh, '  dead ', ndead, '  bad ', nbad
print '(A,ES10.2)', 'max |sum of mixing weights - 1|: ', emax
ok = nbad == 0 .and. nserved > 100 .and. nfresh > 0 .and. ndead > 0 .and. emax <= 1.e-15_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine checks(pc, pv, wn, xx, yy, zz, excl)
   !< Checks 2-3 on mt.
   integer(I4P), intent(in)           :: pc(:,:)
   real(R8P),    intent(in)           :: pv(:), wn(:,:), xx(:), yy(:)
   real(R8P),    intent(in), optional :: zz(:)
   logical,      intent(in), optional :: excl(:)
   integer(I4P)                       :: pp, k, t, d, j, k2
   logical                            :: linked

   d = size(pc, 1)
   do pp = 1, mt%np
      if (mt%status(pp) == MERGE_NONE) cycle
      nserved = nserved + 1
      t = mt%target(pp)
      if (mt%status(pp) /= MERGE_OK .and. mt%status(pp) /= MERGE_CHAIN) then
         nbad = nbad + 1 ; cycle
      endif
      linked = .false.
      do k = lk%start(pp), lk%start(pp + 1) - 1
         if (lk%piece(k) == t) linked = .true.
      enddo
      if (.not. linked .or. pv(t) <= 0._R8P) nbad = nbad + 1
      if (mt%status(pp) == MERGE_OK .and. .not. suitable(t, pc, pv, xx, yy, zz, excl)) nbad = nbad + 1
      if (mt%status(pp) == MERGE_CHAIN) then
         if (suitable(t, pc, pv, xx, yy, zz, excl) .or. mt%final_target(pp) == 0) then
            nbad = nbad + 1
         elseif (.not. suitable(mt%final_target(pp), pc, pv, xx, yy, zz, excl)) then
            nbad = nbad + 1
         endif
      endif
      if (mt%mix_start(pp + 1) - mt%mix_start(pp) > 2 ** d - 1) nbad = nbad + 1
      do k = mt%mix_start(pp), mt%mix_start(pp + 1) - 1
         j = mt%mix_piece(k)
         do i = 1, d
            if (pc(i, j) - pc(i, pp) /= 0 .and. pc(i, j) - pc(i, pp) /= nint(sign(1._R8P, wn(i, pp)))) nbad = nbad + 1
            if (pc(i, j) - pc(i, pp) /= 0 .and. wn(i, pp) == 0._R8P) nbad = nbad + 1
         enddo
         do k2 = mt%mix_start(pp), k - 1
            if (all(pc(:, mt%mix_piece(k2)) == pc(:, j))) nbad = nbad + 1
         enddo
      enddo
      if (mt%mix_start(pp + 1) > mt%mix_start(pp)) &
         emax = max(emax, abs(sum(mt%mix_weight(mt%mix_start(pp):mt%mix_start(pp + 1) - 1)) - 1._R8P))
   enddo
   endsubroutine checks

   logical function suitable(tt, pc, pv, xx, yy, zz, excl)
   !< Piece tt can end a chain: volume, not small, not excluded.
   integer(I4P), intent(in)           :: tt, pc(:,:)
   real(R8P),    intent(in)           :: pv(:), xx(:), yy(:)
   real(R8P),    intent(in), optional :: zz(:)
   logical,      intent(in), optional :: excl(:)
   real(R8P)                          :: cv

   cv = (xx(pc(1, tt) + 1) - xx(pc(1, tt))) * (yy(pc(2, tt) + 1) - yy(pc(2, tt)))
   if (present(zz)) cv = cv * (zz(pc(3, tt) + 1) - zz(pc(3, tt)))
   suitable = pv(tt) >= 0.5_R8P * cv
   if (present(excl)) suitable = suitable .and. .not. excl(tt)
   endfunction suitable
endprogram fossil_test_cut_merge
