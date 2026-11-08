!run figref-overset-projection fig_ref_overset_projection
!fig ref-overset-projection
!< Figure of the reference page "Overset assembly kernels": two blocks of an annulus overlap on its concave outer wall
!< (r = 2, fluid inside); the fine block has 41 vertices over 0.9 pi, the coarse one 9 over 0.4 pi, each wall made of
!< its own chords. The fine wall vertices bulge out of the coarse chords: without projection they have no donor (ember);
!< wall_projection_2d moves them onto the coarse wall for the donor search (teal, the shift in gold). The annulus is
!< unrolled (arc length across, 40 times the distance from the circle r = 2 up), over 0.67 pi to 0.83 pi.
program fig_ref_overset_projection
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, donor_set_t, donor_candidates_2d, &
                       wall_projection_2d, block_real_t
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t)       :: q(2)
type(donor_set_t), allocatable     :: dn(:)
type(block_real_t)                 :: sh(2, 2)
type(vector_R8P),  allocatable     :: a(:), b(:)
logical                            :: miss(41), found(41), view(41)
integer(I4P)                       :: i, j, m, v, w

call arc(q(1), 0._R8P, 0.9_R8P * PI, 41_I4P) ; call arc(q(2), 0.62_R8P * PI, 1.02_R8P * PI, 9_I4P)
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn)
do i = 1, 41 ! the fine wall vertices (j = 9) inside the coarse block's angular range
   v = i + 41 * 8 ; miss(i) = 0.9_R8P * PI * (i - 1) / 40 > 0.62_R8P * PI + 1e-9_R8P .and. dn(1)%start(v + 1) == dn(1)%start(v)
enddo
call wall_projection_2d(q, 1_I4P, 2_I4P, sh(1, 2)%v, max_distance=0.3_R8P)
call donor_candidates_2d(q, dn, shifts=sh)
do i = 1, 41
   v = i + 41 * 8 ; found(i) = 0.9_R8P * PI * (i - 1) / 40 > 0.62_R8P * PI + 1e-9_R8P .and. dn(1)%start(v + 1) > dn(1)%start(v)
enddo
print '(A,I0,A,I0)', 'fine wall vertices over the coarse block without a donor ', count(miss), ', with projection ', &
                     count(.not. found .and. 0.9_R8P * PI * [(i - 1, i = 1, 41)] / 40 > 0.62_R8P * PI + 1e-9_R8P)
print '(A,ES10.3)', 'their largest shift ', &
      maxval([(merge(hypot(sh(1, 2)%v(2 * (i + 328) - 1), sh(1, 2)%v(2 * (i + 328))), 0._R8P, miss(i)), i = 1, 41)])

call fig_open(w, 'ref-overset-projection', 'concave wall, unrolled: no donor (ember), projected (teal)', &
              -90._R8P, 90._R8P, height=300_I4P)
allocate(a(4000), b(4000))
view = [(abs(0.9_R8P * (i - 1) / 40 - 0.75_R8P) <= 0.08_R8P, i = 1, 41)]
m = 0 ! the coarse block's two outer layers, each straight segment drawn in 10 pieces
do i = 1, 8
   call piece(2_I4P, i, 9_I4P, i + 1_I4P, 9_I4P) ; call piece(2_I4P, i, 8_I4P, i + 1_I4P, 8_I4P)
enddo
do i = 1, 9
   call piece(2_I4P, i, 8_I4P, i, 9_I4P)
enddo
call fig_segments(w, 'box', a(1:m), b(1:m))
m = 0
do i = 1, 8
   call piece(2_I4P, i, 9_I4P, i + 1_I4P, 9_I4P)
enddo
call fig_segments(w, 'accent', a(1:m), b(1:m))
m = 0
do i = 1, 40
   call piece(1_I4P, i, 9_I4P, i + 1_I4P, 9_I4P)
enddo
call fig_segments(w, 'body', a(1:m), b(1:m))
call fig_points(w, 'hole', pack([(unroll(q(1)%u(i, 9), q(1)%v(i, 9)), i = 1, 41)], miss .and. view))
m = 0
do i = 1, 41
   if (.not. (miss(i) .and. view(i))) cycle
   v = i + 41 * 8 ; m = m + 1 ; a(m) = unroll(q(1)%u(i, 9), q(1)%v(i, 9))
   b(m) = unroll(q(1)%u(i, 9) + sh(1, 2)%v(2 * v - 1), q(1)%v(i, 9) + sh(1, 2)%v(2 * v))
enddo
call fig_segments(w, 'foot', a(1:m), b(1:m))
call fig_points(w, 'probe', b(1:m))
call fig_label(w, 'accent', unroll(1.995_R8P * cos(0.665_R8P * PI), 1.995_R8P * sin(0.665_R8P * PI)), 'coarse wall chords')
call fig_label(w, 'body', unroll(2.0005_R8P * cos(0.665_R8P * PI), 2.0005_R8P * sin(0.665_R8P * PI)), 'fine wall')
call fig_close(w)
contains
   subroutine arc(blk, t0, t1, ni)
   !< Annulus 1 <= r <= 2 over [t0, t1], 9 radial vertices clustered at the outer wall (face 4).
   type(curvilinear_block_2d_t), intent(inout) :: blk    !< Block.
   real(R8P),                    intent(in)    :: t0, t1 !< Arc.
   integer(I4P),                 intent(in)    :: ni     !< Vertices.
   real(R8P)                                   :: u(ni, 9), v(ni, 9), r, t
   integer(I4P)                                :: i, j
   type(topo_record_t)                         :: rec

   do j = 1, 9
      r = 2._R8P - ((9 - j) / 8._R8P) ** 2
      do i = 1, ni
         t = t0 + (t1 - t0) * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%name = 'wall' ; rec%lo = [1, 9, 1] ; rec%hi = [ni, 9, 1]
   call blk%add_record(rec)
   endsubroutine arc

   subroutine piece(k, i0, j0, i1, j1)
   !< The straight segment between vertices (i0, j0) and (i1, j1) of block k, unrolled in 10 pieces (if in view).
   integer(I4P), intent(in) :: k, i0, j0, i1, j1 !< Block, ends.
   real(R8P)                :: x0(2), x1(2)
   integer(I4P)             :: h

   x0 = [q(k)%u(i0, j0), q(k)%v(i0, j0)] ; x1 = [q(k)%u(i1, j1), q(k)%v(i1, j1)]
   if (min(angle(x0), angle(x1)) < 0.67_R8P * PI .or. max(angle(x0), angle(x1)) > 0.83_R8P * PI) return
   do h = 1, 10
      m = m + 1
      a(m) = unroll(x0(1) + (x1(1) - x0(1)) * (h - 1) / 10._R8P, x0(2) + (x1(2) - x0(2)) * (h - 1) / 10._R8P)
      b(m) = unroll(x0(1) + (x1(1) - x0(1)) * h / 10._R8P, x0(2) + (x1(2) - x0(2)) * h / 10._R8P)
   enddo
   endsubroutine piece

   function unroll(x, y) result(r)
   !< The annulus unrolled: arc length along the wall, 40 times the distance from the circle r = 2.
   real(R8P), intent(in) :: x, y !< Point.
   type(vector_R8P)      :: r    !< Drawn point.

   r = vector_R8P(-2._R8P * angle([x, y]), 40._R8P * (hypot(x, y) - 2._R8P), 0._R8P)
   endfunction unroll

   pure real(R8P) function angle(x)
   !< Polar angle in [0, 2 pi).
   real(R8P), intent(in) :: x(2) !< Point.

   angle = modulo(atan2(x(2), x(1)), 2 * PI)
   endfunction angle
endprogram fig_ref_overset_projection
