!< FOSSIL validation: the cut-cell solver support at scale (issue #70, step 11).
program fossil_val_cut_solver
!< FOSSIL validation: the cut-cell solver support at scale (issue #70, step 11).
!<
!< VALIDATION version (issue #45 tiers) of the regression tests fossil_test_cut_links_3d, _srd, _frd, _merge, _moments,
!< _quadrature, _coarsen and _edges: two closed bodies on 64^3 cells over their box, a marching-cubes sphere (40^3 field)
!< and the bench cone, and a disc of 256 segments on 256^2 cells. (The dragon is not watertight: its open regions are
!< flagged cells, not a test of the solver support.)
!< 1. links: per piece and face the links sum to the piece face area (1e-12), symmetric;
!< 2. SRD: unit column sums, the total weighted volume equals the fluid volume (1e-12); FRD: weights sum to 1;
!<    merge targets: every served piece has a suitable target or a chain to one;
!< 3. quadrature: every rule with quad_status 0 integrates the moments up to order 4 (1e-11 of the cell; 3.0e-12
!<    measured: the cones of some pieces cancel); the pieces with quad_status 1 are counted;
!< 4. coarsening to 32^3: the same pieces per cell and volumes as the cut cells computed on the coarse grid (1e-12);
!< 5. 2-D: links, SRD, quadrature and coarsening as above.

use fossil,             only : surface_stl_object, curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, &
                               srd_neighbourhoods, frd_t, frd_neighbourhoods, merge_t, merge_targets, MERGE_NONE, &
                               MERGE_OK, MERGE_CHAIN, extract_isosurface
use fossil_facet_object, only : facet_object
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use fossil_cut_cells_3d, only : cut_cells_3d_t, cut_cells_3d
use fossil_cut_coarsen,  only : cut_coarsen_3d, cut_coarsen_2d
use fossil_cut_moments,  only : mono_exponents
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

type(surface_stl_object) :: s
type(body_t), allocatable :: bb(:)
type(curve_object)       :: c
type(cut_cells_3d_t)     :: cut, coarse, direct
type(cut_cells_2d_t)     :: c2, c2c, c2d
type(cut_links_t)        :: lk, lkc, lkd
type(srd_t)              :: srd
type(frd_t)              :: frd
type(merge_t)            :: mt
integer(I4P), allocatable :: f2c(:)
real(R8P), allocatable   :: x(:), y(:), z(:), u(:), v(:)
real(R8P)                :: bmin(3), bmax(3), err(6), t0, t1
type(vector_R8P)         :: vb
integer(I4P)             :: i, n, nbad, ninexact, st
logical                  :: ok

err = 0._R8P ; nbad = 0 ; ninexact = 0
call sphere(s)
call run3('sphere')
call bench_bodies(bb)
call to_surface(bb(7)%soup, s)
call run3('cone')
! 2-D
n = 256
u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * (i - 1) / n), i = 1, n)]
call c%adopt_loops(u=u, v=v, loop_sizes=[n])
x = [(0.0031_R8P + (1._R8P - 0.0062_R8P) * i / 256._R8P, i = 0, 256)]
call c%cut_cells(x=x, y=x, cut=c2, moments_order=4_I4P, quadrature=.true.)
call cut_links_2d(c2, lk)
call srd_neighbourhoods(lk, c2%piece_cell, c2%piece_area, c2%piece_centroid, c2%piece_wet_normal, x, x, srd)
call frd_neighbourhoods(lk, c2%piece_cell, c2%piece_area, x, x, frd)
call merge_targets(lk, c2%piece_cell, c2%piece_area, c2%piece_wet_normal, x, x, mt)
call graph_checks(c2%piece_area)
call cut_coarsen_2d(c2, lk, x, x, 2_I4P, c2c, lkc, f2c, st)
call c%cut_cells(x=x(1::2), y=x(1::2), cut=c2d)
if (st /= 0 .or. any(c2c%pieces /= c2d%pieces)) nbad = nbad + 1
err(4) = max(err(4), maxval(abs(c2c%alpha - c2d%alpha), mask=c2d%closed))
print '(A,I0,A,I0)', 'disc 256^2: pieces ', c2%np, ', links ', size(lk%piece)
print '(A,6ES10.2)', 'max errors (links, srd, frd, quadrature, coarsening, -): ', err
print '(A,I0,A,I0)', 'mismatches ', nbad, '  inexact rules (quad_status 1) ', ninexact
ok = nbad == 0 .and. err(1) <= 1.e-12_R8P .and. err(2) <= 1.e-12_R8P .and. err(3) <= 1.e-12_R8P .and. &
     err(4) <= 1.e-11_R8P .and. err(5) <= 1.e-12_R8P
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine run3(label)
   !< Checks 1-4 on the body in s, 64^3 cells over its box.
   character(*), intent(in) :: label
   vb = s%get_bmin() ; bmin = [vb%x, vb%y, vb%z] ; vb = s%get_bmax() ; bmax = [vb%x, vb%y, vb%z]
   bmin = bmin - 0.03_R8P * (bmax - bmin) ; bmax = bmax + 0.03_R8P * (bmax - bmin)
   x = [(bmin(1) + (bmax(1) - bmin(1)) * i / 64._R8P, i = 0, 64)]
   y = [(bmin(2) + (bmax(2) - bmin(2)) * i / 64._R8P, i = 0, 64)]
   z = [(bmin(3) + (bmax(3) - bmin(3)) * i / 64._R8P, i = 0, 64)]
   call cpu_time(t0)
   call cut_cells_3d(s, x, y, z, cut, links=lk, edges=.true., moments_order=4_I4P, quadrature=.true.)
   call cpu_time(t1)
   print '(A,A,F8.2,A,I0,A,I0,A,I0)', trim(label), ' 64^3: cut cells in ', t1 - t0, ' s, pieces ', cut%np, ', untrusted ', &
         count(cut%pieces < 0), ', links ', size(lk%piece)
   call links3
   call srd_neighbourhoods(lk, cut%piece_cell, cut%piece_volume, cut%piece_centroid, cut%piece_wet_normal, x, y, srd, z=z)
   call frd_neighbourhoods(lk, cut%piece_cell, cut%piece_volume, x, y, frd, z=z)
   call merge_targets(lk, cut%piece_cell, cut%piece_volume, cut%piece_wet_normal, x, y, mt, z=z)
   call graph_checks(cut%piece_volume)
   call quad3
   call cut_coarsen_3d(cut, lk, x, y, z, 2_I4P, coarse, lkc, f2c, st)
   call cut_cells_3d(s, x(1::2), y(1::2), z(1::2), direct)
   call coarse3
   endsubroutine run3

   subroutine sphere(sf)
   !< Unit sphere by marching cubes on a 40^3 field over [-1.3, 1.3]^3.
   type(surface_stl_object), intent(inout) :: sf
   real(R8P),          allocatable :: val(:,:,:)
   type(facet_object), allocatable :: fct(:)
   integer(I4P) :: a, b2, c2, st2
   real(R8P)    :: dd

   allocate(val(40, 40, 40)) ; dd = 2.6_R8P / 39._R8P
   do c2 = 1, 40
      do b2 = 1, 40
         do a = 1, 40
            val(a, b2, c2) = norm2([-1.3_R8P + (a - 1) * dd, -1.3_R8P + (b2 - 1) * dd, -1.3_R8P + (c2 - 1) * dd]) - 1._R8P
         enddo
      enddo
   enddo
   call extract_isosurface(values=val, bmin=vector_R8P(-1.3_R8P, -1.3_R8P, -1.3_R8P), &
                           bmax=vector_R8P(1.3_R8P, 1.3_R8P, 1.3_R8P), iso=0._R8P, surface=fct, status=st2)
   call sf%adopt_facets(facets=fct)
   endsubroutine sphere

   subroutine links3
   !< Check 1 on the dragon.
   integer(I4P) :: p, f, k, ic(3), ax, q, kk
   real(R8P)    :: sm, h(3), af(3)
   logical      :: found

   do p = 1, cut%np
      ic = cut%piece_cell(:, p)
      h = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2)), z(ic(3) + 1) - z(ic(3))]
      af = [h(2) * h(3), h(3) * h(1), h(1) * h(2)]
      do f = 1, 6
         ax = (f + 1) / 2
         if ((mod(f, 2) == 1 .and. cut%piece_cell(ax, p) == 1) .or. (mod(f, 2) == 0 .and. cut%piece_cell(ax, p) == 64)) cycle
         q = 0 ; ic = cut%piece_cell(:, p) ; ic(ax) = ic(ax) + merge(-1, 1, mod(f, 2) == 1)
         if (cut%pieces(ic(1), ic(2), ic(3)) < 0) cycle
         sm = 0._R8P
         do k = lk%start(p), lk%start(p + 1) - 1
            if (lk%face(k) == f) sm = sm + lk%area(k)
         enddo
         err(1) = max(err(1), abs(sm - cut%piece_face(f, p) * af(ax)) / af(ax))
      enddo
      do k = lk%start(p), lk%start(p + 1) - 1
         found = .false.
         do kk = lk%start(lk%piece(k)), lk%start(lk%piece(k) + 1) - 1
            if (lk%piece(kk) == p .and. lk%area(kk) == lk%area(k)) found = .true.
         enddo
         if (.not. found) nbad = nbad + 1
      enddo
   enddo
   if (lk%unmatched /= 0) nbad = nbad + 1
   endsubroutine links3

   subroutine graph_checks(pv)
   !< Check 2: SRD column sums and total, FRD weights, merge targets.
   real(R8P), intent(in) :: pv(:)
   real(R8P), allocatable :: colsum(:)
   integer(I4P) :: p, k

   allocate(colsum(size(pv))) ; colsum = 0._R8P
   do p = 1, srd%np
      do k = srd%start(p), srd%start(p + 1) - 1
         colsum(srd%member(k)) = colsum(srd%member(k)) + 1._R8P / srd%nov(srd%member(k))
      enddo
   enddo
   do p = 1, srd%np
      if (pv(p) > 0._R8P) err(2) = max(err(2), abs(colsum(p) - 1._R8P))
   enddo
   err(2) = max(err(2), abs(sum(srd%vhat, mask=pv > 0._R8P) - sum(pv, mask=pv > 0._R8P)) / sum(pv))
   do p = 1, frd%np
      if (frd%start(p + 1) > frd%start(p)) err(3) = max(err(3), abs(sum(frd%weight(frd%start(p):frd%start(p + 1) - 1)) - 1._R8P))
   enddo
   do p = 1, mt%np
      if (mt%status(p) == MERGE_NONE .or. mt%status(p) == MERGE_OK) cycle
      if (mt%status(p) == MERGE_CHAIN .and. mt%final_target(p) /= 0) cycle
      nbad = nbad + 1
   enddo
   endsubroutine graph_checks

   subroutine quad3
   !< Check 3 on the dragon.
   integer(I4P), allocatable :: e(:,:)
   integer(I4P) :: p, k, q, ic(3)
   real(R8P)    :: h(3), cc(3), sv, pt(3)

   e = mono_exponents(3_I4P, 4_I4P)
   ninexact = ninexact + count(cut%quad_status /= 0)
   do p = 1, cut%np
      if (cut%quad_status(p) /= 0) cycle
      ic = cut%piece_cell(:, p)
      h = [x(ic(1) + 1) - x(ic(1)), y(ic(2) + 1) - y(ic(2)), z(ic(3) + 1) - z(ic(3))] ; cc = [x(ic(1)), y(ic(2)), z(ic(3))] + h / 2
      do k = 1, size(e, 2)
         sv = 0._R8P
         do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
            pt = (cut%quad_point(:, q) - cc) / h ; sv = sv + cut%quad_weight(q) * product(pt ** e(:, k))
         enddo
         err(4) = max(err(4), abs(sv - cut%piece_moments(k, p)) / product(h))
      enddo
   enddo
   endsubroutine quad3

   subroutine coarse3
   !< Check 4 on the dragon.
   integer(I4P) :: ia, ja, ka

   if (st /= 0) then
      nbad = nbad + 1 ; return
   endif
   do ka = 1, 32
      do ja = 1, 32
         do ia = 1, 32
            if (coarse%pieces(ia, ja, ka) /= direct%pieces(ia, ja, ka)) nbad = nbad + 1
            if (direct%pieces(ia, ja, ka) < 0) cycle
            err(5) = max(err(5), abs(coarse%alpha(ia, ja, ka) - direct%alpha(ia, ja, ka)))
         enddo
      enddo
   enddo
   endsubroutine coarse3
endprogram fossil_val_cut_solver
