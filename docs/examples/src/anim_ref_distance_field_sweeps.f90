!run anim-ref-distance-field-sweeps anim_ref_distance_field_sweeps
!fig anim-ref-distance-field-sweeps
!< Animation of the reference page "Signed-distance field", section "Pipeline": the three stages on a 2-D field of an
!< L-shaped body and a disc, on a stretched grid. compute_band, compute_signs, then sweep stopped after 1, 2, 3 and 4
!< directional sweeps (each frame is a copy of the signed field swept with max_sweeps = k), then swept to convergence.
program anim_ref_distance_field_sweeps
use fossil,     only : curve_object, distance_field_object, FAR_VALUE, STATUS_OK
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_quad, fig_points, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter     :: PI = acos(-1._R8P)
integer(I4P), parameter     :: NDISC = 40
type(curve_object)          :: body
type(distance_field_object) :: field, signed, swept, final
real(R8P), allocatable      :: x(:), y(:), u(:), v(:)
real(R8P)                   :: neg, pos
integer(I4P)                :: i, k, sweeps, fallback, status, uf
character(120)              :: caption

! an L (outer loop, counter-clockwise) and a disc of radius 0.45 centred at (0.9, 0.85)
u = [-1.2_R8P, 0.2_R8P, 0.2_R8P, -0.4_R8P, -0.4_R8P, -1.2_R8P]
v = [-1.1_R8P, -1.1_R8P, -0.4_R8P, -0.4_R8P, 1.0_R8P, 1.0_R8P]
u = [u, (0.9_R8P + 0.45_R8P * cos(2 * PI * i / NDISC), i = 0, NDISC - 1)]
v = [v, (0.85_R8P + 0.45_R8P * sin(2 * PI * i / NDISC), i = 0, NDISC - 1)]
call body%adopt_loops(u=u, v=v, loop_sizes=[6_I4P, NDISC], status=status)
if (status /= STATUS_OK) error stop 'adopt_loops failed'
! each axis: 27 nodes with spacing 0.125 on [-1.625, 1.625], then 6 per side growing by 1.3
x = [(-1.625_R8P + 0.125_R8P * i, i = 0, 26)]
call stretch(x) ; y = x
call field%initialize(x=x, y=y, status=status)

!region stages
call field%compute_band(curve=body, status=status)                ! exact distances within the band
if (status /= STATUS_OK) error stop 'compute_band failed'
signed = field                                                    ! a copy: each stage kept for the frames
call signed%compute_signs(curve=body, fallback_cells=fallback, status=status)  ! flood fill of the signs
final = signed
call final%sweep(sweeps=sweeps, status=status)                    ! fast sweeping to convergence
!endregion stages
swept = final
neg = minval(swept%phi) ; pos = maxval(swept%phi)
print '(A,I0,A,I0,A,I0,A,F6.3)', 'grid ', size(x), ' x ', size(y), ' nodes, band nodes ', field%get_band_cells(), &
      ', band width ', field%get_band_width()
print '(A,I0)', 'signs: nodes signed by their own point query ', fallback

call fig_open(uf, 'anim-ref-distance-field-sweeps', '', -90._R8P, 90._R8P, width=640_I4P, height=600_I4P)
call fig_anim(uf, 1.4_R8P)
write(caption, '(A,I0,A)') 'compute_band: ', field%get_band_cells(), ' exact nodes, the rest unknown (blank)'
call frame(field, caption, .false., 1.8_R8P)
write(caption, '(A,I0,A)') 'compute_signs: flood fill, ', count(signed%phi < -0.5_R8P * FAR_VALUE), &
                           ' unknown nodes inside (dots)'
call frame(signed, caption, .true., 1.8_R8P)
do k = 1, 4
   !region partial
   swept = signed
   call swept%sweep(max_sweeps=k, status=status)                  ! stop after k directional sweeps
   !endregion partial
   if (any(abs(swept%phi) >= FAR_VALUE)) then
      write(caption, '(A,I0,A,A,A,I0,A)') 'sweep ', k, ' (', trim(direction(k)), '): ', &
                                          count(abs(swept%phi) >= FAR_VALUE), ' nodes still unknown'
   else
      write(caption, '(A,I0,A,A,A,ES8.2)') 'sweep ', k, ' (', trim(direction(k)), '): all known, off the converged field by ', &
                                          maxval(abs(swept%phi - final%phi))
   endif
   call frame(swept, caption, .true., 1.4_R8P)
enddo
write(caption, '(A,I0,A)') 'converged after ', sweeps, ' sweeps: no change above tol in a cycle of 4'
call frame(final, caption, .true., 2.6_R8P)
call fig_close(uf)
contains
   subroutine stretch(c)
   !< 6 nodes per side growing by 1.3 from the uniform core.
   real(R8P), allocatable, intent(inout) :: c(:)
   real(R8P)                             :: h
   integer(I4P)                          :: m

   h = c(2) - c(1)
   do m = 1, 6
      h = 1.3_R8P * h
      c = [c(1) - h, c, c(size(c)) + h]
   enddo
   endsubroutine stretch

   function direction(kk) result(d)
   !< The sweep direction of directional sweep kk (bit 0: x reversed, bit 1: y reversed).
   integer(I4P), intent(in) :: kk
   character(16)            :: d

   d = trim(merge('x down', 'x up  ', btest(kk - 1, 0))) // ', ' // trim(merge('y down', 'y up  ', btest(kk - 1, 1)))
   endfunction direction

   subroutine frame(f, cap, signs, seconds)
   !< One frame: known nodes coloured by phi (teal inside, ember outside, sqrt scale), the
   !< unknown nodes blank, the unknown inside nodes dotted when the signs are known; the body outline.
   type(distance_field_object), intent(in) :: f
   character(*),                intent(in) :: cap
   logical,                     intent(in) :: signs
   real(R8P),                   intent(in) :: seconds
   type(vector_R8P), allocatable           :: dots(:)
   real(R8P)                               :: val, a(2), b(2)
   integer(I4P)                            :: ii, jj, s

   print '(A)', trim(cap)
   call fig_frame(uf, caption=trim(cap), seconds=seconds)
   allocate(dots(0))
   do jj = 1, size(y)
      do ii = 1, size(x)
         if (abs(f%phi(ii, jj, 1)) >= FAR_VALUE) then
            if (signs .and. f%phi(ii, jj, 1) < 0._R8P) dots = [dots, p(x(ii), y(jj), 1.e-3_R8P)]
         else
            if (f%phi(ii, jj, 1) < 0._R8P) then
               val = -sqrt(f%phi(ii, jj, 1) / neg)
            else
               val = sqrt(f%phi(ii, jj, 1) / pos)
            endif
            call fig_quad(uf, 'field', p(mid(x, ii, -1), mid(y, jj, -1)), p(mid(x, ii, 1), mid(y, jj, -1)), &
                          p(mid(x, ii, 1), mid(y, jj, 1)), p(mid(x, ii, -1), mid(y, jj, 1)), val)
         endif
      enddo
   enddo
   if (size(dots) > 0) call fig_points(uf, 'probe', dots)
   do s = 1, body%get_vertices_number()
      a = body%get_vertex(s) ; b = body%get_vertex(body%get_next(s))
      call fig_segment(uf, 'accent', p(a(1), a(2), 2.e-3_R8P), p(b(1), b(2), 2.e-3_R8P))
   enddo
   endsubroutine frame

   pure function mid(c, ii, s) result(r)
   !< The dual-cell side of node ii towards s = -1 or +1 (the node itself at the ends).
   real(R8P),    intent(in) :: c(:)
   integer(I4P), intent(in) :: ii, s
   real(R8P)                :: r

   if (ii + s < 1 .or. ii + s > size(c)) then
      r = c(ii)
   else
      r = 0.5_R8P * (c(ii) + c(ii + s))
   endif
   endfunction mid

   function p(px, py, dz) result(q)
   !< The point (px, py, dz) of the figure, dz = 0 by default.
   real(R8P), intent(in)           :: px, py
   real(R8P), intent(in), optional :: dz
   type(vector_R8P)                :: q

   q = vector_R8P(px, py, 0._R8P)
   if (present(dz)) q%z = dz
   endfunction p
endprogram anim_ref_distance_field_sweeps
