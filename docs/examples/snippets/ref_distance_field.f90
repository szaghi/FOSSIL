!< Reference page "Signed-distance field": the unit cube on a stretched grid against its exact distance, band and swept
!< far field apart; the stored closest points of the band; the eikonal diagnostic. Figures: the error and q on the
!< grid plane z = 0.5.
program ref_distance_field
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segment, fig_close
use fossil,     only : surface_stl_object, distance_field_object, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)    :: cube
type(distance_field_object) :: field
real(R8P), allocatable      :: x(:), err(:,:,:), q(:,:,:)
real(R8P)                   :: h, e(3), exact, band_err, far_err, foot_err
integer(I4P)                :: i, j, k, m, n, sweeps, fallback, flagged, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! [0, 1]^3
! each axis: 29 nodes with spacing 0.05 on [-0.2, 1.2], then 12 per side growing by 1.15
x = [(-0.2_R8P + 0.05_R8P * i, i = 0, 28)]
h = 0.05_R8P
do m = 1, 12
   h = 1.15_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
n = size(x)
call field%initialize(x=x, y=x, z=x, status=status)
call field%compute(surface=cube, sweeps=sweeps, fallback_cells=fallback, store_closest=.true., status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
print '(A,I0,A,F7.3,A,F7.3)', 'grid ', n, '^3 on ', x(1), ' to ', x(n)
print '(A,I0,A,F6.3,A,I0,A,I0)', 'band cells ', field%get_band_cells(), ' (half-width ', field%get_band_width(), &
                                 '), point-query cells ', fallback, ', sweeps ', sweeps
allocate(err(n, n, n))
band_err = 0._R8P ; far_err = 0._R8P ; foot_err = 0._R8P
do k = 1, n ; do j = 1, n ; do i = 1, n
   e = abs([x(i), x(j), x(k)] - 0.5_R8P) - 0.5_R8P                        ! exact distance of the box
   exact = norm2(max(e, 0._R8P)) + min(maxval(e), 0._R8P)
   err(i, j, k) = abs(field%phi(i, j, k) - exact)
   if (field%in_band(i, j, k)) then
      band_err = max(band_err, err(i, j, k))
      foot_err = max(foot_err, abs(norm2([x(i), x(j), x(k)] - field%closest(:, i, j, k)) - abs(field%phi(i, j, k))))
   else
      far_err = max(far_err, err(i, j, k))
   endif
enddo ; enddo ; enddo
print '(A,ES9.2,A,F6.3,A,F6.3)', 'max |phi - exact|: band ', band_err, ', far field ', far_err, &
                                 ', largest spacing ', maxval(x(2:) - x(:n - 1))
print '(A,ES9.2)', 'band: max ||x - closest| - |phi|| ', foot_err
allocate(q(n, n, n))
call field%gradient_diagnostic(q=q, flagged=flagged)
print '(A,I0,A,I0,A,I0)', 'eikonal q > 0.05 at ', flagged, ' of ', n**3, ' nodes; unknown neighbour ', count(q < 0._R8P)
print '(A,I0,A,I0)', '  of them in the band ', count(q > 0.05_R8P .and. field%in_band), ', inside the cube ', &
                     count(q > 0.05_R8P .and. field%phi < 0._R8P)
call draw
contains
   subroutine draw
   !< The error and q on the grid plane z = 0.5, each node drawn as its dual cell, and the cube's section.
   real(R8P), allocatable :: lo(:), hi(:)
   integer(I4P)           :: u, kz
   character(100)         :: caption

   kz = minloc(abs(x - 0.5_R8P), dim=1)
   lo = [x(1), (0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1)]
   hi = [(0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1), x(n)]
   write(caption, '(A,F5.3,A)') '|phi - exact| on z = 0.5, 0 (cream) to ', maxval(err(:, :, kz)), &
                                ' (sepia); band outlined'
   call fig_open(u, 'ref-distance-field-error', trim(caption), -90._R8P, 90._R8P, height=520_I4P)
   do j = 1, n ; do i = 1, n
      call fig_quad(u, 'frac', p(lo(i), lo(j)), p(hi(i), lo(j)), p(hi(i), hi(j)), p(lo(i), hi(j)), &
                    err(i, j, kz) / maxval(err(:, :, kz)))
   enddo ; enddo
   do j = 1, n ; do i = 1, n                                              ! the band's outer border
      if (.not. field%in_band(i, j, kz)) cycle
      if (i < n) then
         if (.not. field%in_band(i + 1, j, kz)) call fig_segment(u, 'probe', p(hi(i), lo(j), 1.e-3_R8P), &
                                                                 p(hi(i), hi(j), 1.e-3_R8P))
      endif
      if (i > 1) then
         if (.not. field%in_band(i - 1, j, kz)) call fig_segment(u, 'probe', p(lo(i), lo(j), 1.e-3_R8P), &
                                                                 p(lo(i), hi(j), 1.e-3_R8P))
      endif
      if (j < n) then
         if (.not. field%in_band(i, j + 1, kz)) call fig_segment(u, 'probe', p(lo(i), hi(j), 1.e-3_R8P), &
                                                                 p(hi(i), hi(j), 1.e-3_R8P))
      endif
      if (j > 1) then
         if (.not. field%in_band(i, j - 1, kz)) call fig_segment(u, 'probe', p(lo(i), lo(j), 1.e-3_R8P), &
                                                                 p(hi(i), lo(j), 1.e-3_R8P))
      endif
   enddo ; enddo
   call square(u)
   call fig_close(u)

   write(caption, '(A)') 'eikonal q = |1 - |grad phi|| on z = 0.5, 0 (cream) to 0.5 and above (sepia)'
   call fig_open(u, 'ref-distance-field-eikonal', trim(caption), -90._R8P, 90._R8P, height=520_I4P)
   call fig_zoom(u, p(0.5_R8P, 0.5_R8P), 1.1_R8P)
   do j = 1, n ; do i = 1, n
      call fig_quad(u, 'frac', p(lo(i), lo(j)), p(hi(i), lo(j)), p(hi(i), hi(j)), p(lo(i), hi(j)), &
                    min(1._R8P, max(0._R8P, q(i, j, kz)) / 0.5_R8P))
   enddo ; enddo
   call square(u)
   call fig_close(u)
   endsubroutine draw

   subroutine square(u)
   !< The section of the cube, just above the plane.
   integer(I4P), intent(in) :: u

   call fig_segment(u, 'accent', p(0._R8P, 0._R8P, 2.e-3_R8P), p(1._R8P, 0._R8P, 2.e-3_R8P))
   call fig_segment(u, 'accent', p(1._R8P, 0._R8P, 2.e-3_R8P), p(1._R8P, 1._R8P, 2.e-3_R8P))
   call fig_segment(u, 'accent', p(1._R8P, 1._R8P, 2.e-3_R8P), p(0._R8P, 1._R8P, 2.e-3_R8P))
   call fig_segment(u, 'accent', p(0._R8P, 1._R8P, 2.e-3_R8P), p(0._R8P, 0._R8P, 2.e-3_R8P))
   endsubroutine square

   function p(a, b, c) result(r)
   !< The point (a, b, c) of the figure, c = 0 by default.
   real(R8P), intent(in)           :: a, b
   real(R8P), intent(in), optional :: c
   type(vector_R8P)                :: r

   r = vector_R8P(a, b, 0._R8P)
   if (present(c)) r%z = c
   endfunction p
endprogram ref_distance_field
