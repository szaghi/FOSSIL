!< Figure of the reference page "Immersed-boundary queries", section "Wall distance and field diagnostic": the curve of
!< two bodies (a unit square and a 24-gon of radius 0.4), a grid of spacing 1/16 around it.
!< - Left: `wall_distance_many` with cutoff 0.4 at the cell centres, shaded from cream (on the wall) to sepia (the
!<   cutoff), the cells beyond the cutoff (body 0) not drawn; teal edges separate cells whose nearest walls belong to
!<   different bodies.
!< - Right: the signed distance field of the same curve on the same grid (`distance_field_object%compute`) and its
!<   `gradient_diagnostic` q = |1 - |grad phi||: the nodes above the threshold 0.05 in ember (teal: inside, cream:
!<   outside). They trace the medial axis (the square's diagonals, the 24-gon's centre, the line midway between the
!<   bodies), where phi has a kink, and four patches off the square's corners, in the far field the sweep filled.
program fig_ref_ib_queries_wall
use fossil,     only : curve_object, distance_field_object, STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), H = 1._R8P / 16, CUTOFF = 0.4_R8P, X0 = -0.6_R8P, Y0 = -0.6_R8P
real(R8P),    parameter       :: SHIFT = 3.5_R8P
integer(I4P), parameter       :: NX = 52, NY = 35, NC = 24
type(curve_object)            :: c
type(distance_field_object)   :: field
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: u(4 + NC), v(4 + NC), p(2, NX * NY), d(NX * NY), xs(NX), ys(NY), q(NX, NY, 1)
integer(I4P)                  :: bod(NX * NY), i, j, k, m, uu, flagged, status
character(140)                :: caption

u(1:4) = [0._R8P, 1._R8P, 1._R8P, 0._R8P] ; v(1:4) = [0._R8P, 0._R8P, 1._R8P, 1._R8P]
do k = 0, NC - 1
   u(5 + k) = 2._R8P + 0.4_R8P * cos(2 * PI * k / NC) ; v(5 + k) = 0.5_R8P + 0.4_R8P * sin(2 * PI * k / NC)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, NC])
xs = [(X0 + H * (i - 0.5_R8P), i = 1, NX)] ; ys = [(Y0 + H * (j - 0.5_R8P), j = 1, NY)]
do j = 1, NY
   do i = 1, NX
      p(:, i + NX * (j - 1)) = [xs(i), ys(j)]
   enddo
enddo
call c%wall_distance_many(points=p, distances=d, bodies=bod, cutoff=CUTOFF)
print '(A,I0,A,I0,A,I0,A,I0)', 'wall distance, cutoff 0.4: ', NX * NY, ' points, body 1 ', count(bod == 1), &
                               ', body 2 ', count(bod == 2), ', beyond the cutoff ', count(bod == 0)
call field%initialize(x=xs, y=ys)
call field%compute(curve=c, status=status)
if (status /= STATUS_OK) error stop 'field: compute failed'
call field%gradient_diagnostic(q=q, flagged=flagged)
print '(A,I0,A,I0,A,F6.3)', 'gradient diagnostic: ', flagged, ' of ', NX * NY, ' nodes above 0.05, max q ', maxval(q)
print '(A,I0,A,I0)', '  flagged in the band ', count(q(:, :, 1) > 0.05_R8P .and. field%in_band(:, :, 1)), &
                     ', in the swept far field ', count(q(:, :, 1) > 0.05_R8P .and. .not. field%in_band(:, :, 1))
write(caption, '(A,I0,A,I0,A)') 'left: wall distance up to the cutoff 0.4 (', count(bod == 0), &
      ' beyond); right: |1 - |grad phi|| > 0.05 at ', flagged, ' nodes (ember)'
call fig_open(uu, 'ref-ib-queries-wall', trim(caption), -90._R8P, 90._R8P, width=900_I4P, height=360_I4P)
do j = 1, NY
   do i = 1, NX
      k = i + NX * (j - 1)
      if (bod(k) > 0) call fig_quad(uu, 'frac', pt(xs(i) - H / 2, ys(j) - H / 2), pt(xs(i) + H / 2, ys(j) - H / 2), &
                                    pt(xs(i) + H / 2, ys(j) + H / 2), pt(xs(i) - H / 2, ys(j) + H / 2), value=d(k) / CUTOFF)
      call fig_quad(uu, 'field', pt(SHIFT + xs(i) - H / 2, ys(j) - H / 2), pt(SHIFT + xs(i) + H / 2, ys(j) - H / 2), &
                    pt(SHIFT + xs(i) + H / 2, ys(j) + H / 2), pt(SHIFT + xs(i) - H / 2, ys(j) + H / 2), &
                    value=merge(0.9_R8P, merge(-0.45_R8P, 0._R8P, field%phi(i, j, 1) < 0._R8P), q(i, j, 1) > 0.05_R8P))
   enddo
enddo
! the edges between cells whose nearest walls belong to different bodies
allocate(a(2 * NX * NY), b(2 * NX * NY)) ; m = 0
do j = 1, NY
   do i = 1, NX
      k = i + NX * (j - 1)
      if (i < NX) then
         if (bod(k) * bod(k + 1) > 0 .and. bod(k) /= bod(k + 1)) then
            m = m + 1 ; a(m) = pt(xs(i) + H / 2, ys(j) - H / 2, 1.e-3_R8P) ; b(m) = pt(xs(i) + H / 2, ys(j) + H / 2, 1.e-3_R8P)
         endif
      endif
      if (j < NY) then
         if (bod(k) * bod(k + NX) > 0 .and. bod(k) /= bod(k + NX)) then
            m = m + 1 ; a(m) = pt(xs(i) - H / 2, ys(j) + H / 2, 1.e-3_R8P) ; b(m) = pt(xs(i) + H / 2, ys(j) + H / 2, 1.e-3_R8P)
         endif
      endif
   enddo
enddo
call fig_segments(uu, 'probe', a(1:m), b(1:m))
deallocate(a, b) ; allocate(a(2 * (4 + NC)), b(2 * (4 + NC)))
do k = 1, 4
   a(k) = pt(u(k), v(k), 2.e-3_R8P) ; b(k) = pt(u(mod(k, 4) + 1), v(mod(k, 4) + 1), 2.e-3_R8P)
enddo
do k = 1, NC
   a(4 + k) = pt(u(4 + k), v(4 + k), 2.e-3_R8P) ; b(4 + k) = pt(u(4 + mod(k, NC) + 1), v(4 + mod(k, NC) + 1), 2.e-3_R8P)
enddo
a(5 + NC:) = [(pt(a(k)%x + SHIFT, a(k)%y, 2.e-3_R8P), k = 1, 4 + NC)]
b(5 + NC:) = [(pt(b(k)%x + SHIFT, b(k)%y, 2.e-3_R8P), k = 1, 4 + NC)]
call fig_segments(uu, 'accent', a, b)
call fig_close(uu)
contains
   pure function pt(x, y, z) result(r)
   !< A point of the plane, lifted by z (optional).
   real(R8P), intent(in)           :: x, y
   real(R8P), intent(in), optional :: z
   type(vector_R8P)                :: r

   r%x = x ; r%y = y ; r%z = 0._R8P
   if (present(z)) r%z = z
   endfunction pt
endprogram fig_ref_ib_queries_wall
