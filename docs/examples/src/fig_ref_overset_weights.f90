!run figref-overset-weights fig_ref_overset_weights
!fig ref-overset-weights
!< Figure of the reference page "Overset assembly kernels", section "Intersection weights": the receptor cells of the
!< outer ring of a coarse O-grid (16 x 4 cells, r in [1, 2]) against the Cartesian donor grid of spacing 0.15
!< (overlap_weights_2d). One receptor is drawn with its intersections V_c n V_c' filled by weight (cream to sepia) and
!< labelled, their centroids in teal; the others of the ring in outline. The program clips each receptor quad against
!< its donor cells on its own (Sutherland-Hodgman, convex cells) and compares the areas with the library's.
program fig_ref_overset_weights
use fossil,     only : curvilinear_block_2d_t, overlap_weights_t, overlap_weights_2d
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_triangle, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter      :: PI = acos(-1._R8P), H = 0.15_R8P, LIFT = 1.e-3_R8P
integer(I4P), parameter      :: SHOW = 3 !< The receptor drawn in full.
type(curvilinear_block_2d_t) :: block
type(overlap_weights_t)      :: w
type(vector_R8P), allocatable :: a(:), b(:), c(:), z(:)
real(R8P)                    :: u(17, 5), v(17, 5), x(41), r, px(16), py(16), dmax, wmax
integer(I4P)                 :: u_, i, j, k, l, n, st, cells(2, 16), nz
character(150)               :: caption, txt

do j = 1, 5
   r = 1._R8P + (j - 1) / 4._R8P
   do i = 1, 17
      u(i, j) = r * cos(2 * PI * (i - 1) / 16) ; v(i, j) = r * sin(2 * PI * (i - 1) / 16)
   enddo
enddo
u(17,:) = u(1,:) ; v(17,:) = v(1,:)
call block%initialize(u, v)
x = [(-3._R8P + H * (i - 1), i = 1, 41)]
cells(1, :) = [(i, i = 1, 16)] ; cells(2, :) = 4
call overlap_weights_2d(block, cells, x, x, w, status=st)
dmax = 0._R8P
do k = 1, w%nr ! the areas against an independent clip
   do l = w%start(k), w%start(k + 1) - 1
      call clip(k, w%donor(:, l), px, py, n)
      dmax = max(dmax, abs(abs(shoelace(px(1:n), py(1:n))) - w%volume(l)))
   enddo
enddo
print '(A,I0,A,I0,A,I0,A,ES8.1)', 'status ', st, ', receptors ', w%nr, ', intersections ', w%nw, &
                                  ', max |area - independent clip| ', dmax
print '(A,I0,A,I0,A,F8.5)', 'receptor ', SHOW, ': donors ', w%start(SHOW + 1) - w%start(SHOW), ', area ', &
                            w%receptor_volume(SHOW)
wmax = maxval(w%weight(w%start(SHOW):w%start(SHOW + 1) - 1))
write(caption, '(A,I0,A,I0,A)') 'receptor ', SHOW, ': ', w%start(SHOW + 1) - w%start(SHOW), &
                                ' donors, intersections by weight, centroids teal; one of zero area (ember)'
call fig_open(u_, 'ref-overset-weights', trim(caption), -90._R8P, 90._R8P, width=760_I4P, height=600_I4P)
call fig_zoom(u_, pt(0.25_R8P * [sum(u(SHOW:SHOW + 1, 4:5)), sum(v(SHOW:SHOW + 1, 4:5))], 0._R8P), 0.42_R8P)
allocate(c(0), z(0)) ; nz = 0
do l = w%start(SHOW), w%start(SHOW + 1) - 1
   if (w%volume(l) == 0._R8P) then ! a donor of zero area (a side through a grid node): its centroid is the cell centre
      nz = nz + 1 ; z = [z, pt(w%centroid(:, l), 3 * LIFT)] ; cycle
   endif
   call clip(SHOW, w%donor(:, l), px, py, n)
   do i = 2, n - 1 ! fan
      call fig_triangle(u_, 'frac', pt([px(1), py(1)], 0._R8P), pt([px(i), py(i)], 0._R8P), &
                        pt([px(i + 1), py(i + 1)], 0._R8P), w%weight(l) / wmax)
   enddo
   c = [c, pt(w%centroid(:, l), 3 * LIFT)]
   if (w%weight(l) < 0.01_R8P) cycle
   write(txt, '(F5.3)') w%weight(l)
   call fig_label(u_, 'foot', pt(w%centroid(:, l) + [-0.035_R8P, 0.004_R8P], 5 * LIFT), trim(txt))
enddo
print '(A,I0)', 'zero-area donors of that receptor: ', nz
allocate(a(0), b(0))
do i = 1, 41 ! donor grid lines
   a = [a, pt([x(i), x(1)], LIFT), pt([x(1), x(i)], LIFT)] ; b = [b, pt([x(i), x(41)], LIFT), pt([x(41), x(i)], LIFT)]
enddo
call fig_segments(u_, 'box', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do j = 1, 5 ! receptor block lines
   do i = 1, 16
      a = [a, pt([u(i, j), v(i, j)], 2 * LIFT)] ; b = [b, pt([u(i + 1, j), v(i + 1, j)], 2 * LIFT)]
   enddo
enddo
do i = 1, 16
   a = [a, pt([u(i, 1), v(i, 1)], 2 * LIFT)] ; b = [b, pt([u(i, 5), v(i, 5)], 2 * LIFT)]
enddo
call fig_segments(u_, 'patch', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do l = 0, 3 ! the receptor drawn in full
   i = SHOW + merge(1, 0, l == 1 .or. l == 2) ; j = 4 + merge(1, 0, l >= 2)
   k = SHOW + merge(1, 0, l == 0 .or. l == 1) ; n = 4 + merge(1, 0, l == 1 .or. l == 2)
   a = [a, pt([u(i, j), v(i, j)], 3 * LIFT)] ; b = [b, pt([u(k, n), v(k, n)], 3 * LIFT)]
enddo
call fig_segments(u_, 'accent', a, b)
call fig_points(u_, 'probe', c)
call fig_points(u_, 'patch', z)
call fig_label(u_, 'patch', z(1) + 0.01_R8P * ex_R8P, 'zero area')
call fig_label(u_, 'accent', pt([u(SHOW, 5), v(SHOW, 5)], 5 * LIFT), 'receptor')
call fig_close(u_)
contains
   function pt(xy, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: xy(2) !< In-plane coordinates.
   real(R8P), intent(in) :: dz    !< Lift.
   type(vector_R8P)      :: r     !< Point.

   r = xy(1) * ex_R8P + xy(2) * ey_R8P + dz * ez_R8P
   endfunction pt

   subroutine clip(k, d, px, py, n)
   !< Receptor k (a convex quad) clipped to donor cell d: Sutherland-Hodgman on the 4 sides.
   integer(I4P), intent(in)  :: k, d(2)        !< Receptor, donor cell.
   real(R8P),    intent(out) :: px(:), py(:)   !< Polygon.
   integer(I4P), intent(out) :: n             !< Vertices.
   real(R8P)                 :: qx(16), qy(16), lo(2), hi(2), s0, s1, t
   integer(I4P)              :: side, m, l0, l1, ax, ic
   logical                   :: keep_hi

   ic = cells(1, k)
   px(1:4) = [u(ic, 4), u(ic + 1, 4), u(ic + 1, 5), u(ic, 5)] ; py(1:4) = [v(ic, 4), v(ic + 1, 4), v(ic + 1, 5), v(ic, 5)]
   n = 4 ; lo = [x(d(1)), x(d(2))] ; hi = [x(d(1) + 1), x(d(2) + 1)]
   do side = 1, 4
      ax = (side + 1) / 2 ; keep_hi = mod(side, 2) == 1 ! keep >= lo, then <= hi
      m = 0
      do l0 = 1, n
         l1 = mod(l0, n) + 1
         s0 = merge(px(l0), py(l0), ax == 1) - merge(lo(ax), hi(ax), keep_hi)
         s1 = merge(px(l1), py(l1), ax == 1) - merge(lo(ax), hi(ax), keep_hi)
         if (.not. keep_hi) then
            s0 = -s0 ; s1 = -s1
         endif
         if (s0 >= 0._R8P) then
            m = m + 1 ; qx(m) = px(l0) ; qy(m) = py(l0)
         endif
         if ((s0 >= 0._R8P) .neqv. (s1 >= 0._R8P)) then
            t = s0 / (s0 - s1) ; m = m + 1
            qx(m) = px(l0) + t * (px(l1) - px(l0)) ; qy(m) = py(l0) + t * (py(l1) - py(l0))
         endif
      enddo
      n = m ; px(1:n) = qx(1:n) ; py(1:n) = qy(1:n)
   enddo
   endsubroutine clip

   pure real(R8P) function shoelace(px, py) result(a)
   !< Area of a polygon.
   real(R8P), intent(in) :: px(:), py(:) !< Vertices.
   integer(I4P)          :: l             !< Counter.

   a = 0._R8P
   do l = 1, size(px)
      a = a + px(l) * py(mod(l, size(px)) + 1) - px(mod(l, size(px)) + 1) * py(l)
   enddo
   a = 0.5_R8P * a
   endfunction shoelace
endprogram fig_ref_overset_weights
