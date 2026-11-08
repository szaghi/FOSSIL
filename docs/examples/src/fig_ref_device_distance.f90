!run figref-device-distance fig_ref_device_distance
!fig ref-device-distance-marker
!fig ref-device-distance-signs
!< Figures of the reference page "GPU distance queries", on the 24^3 marker block of its example (the open dragon's box
!< grown by 10 %), in planes y = const of cell centres seen from -y, the dragon's section by the plane in front:
!< - the middle plane, the marker: inside (teal), cut (gold: |d| within half a cell diagonal) and outside (ember)
!<   cells, from distance_many_device with the winding-number sign; the caption counts the cells that differ from
!<   distance_many;
!< - the plane where the two device signs disagree most, shaded by the winding-number signed distance, the cells where
!<   the pseudo-normal gives the other sign outlined.
!< Built with gfortran, the device loop runs on the host: the device code path, not a GPU run.
program fig_ref_device_distance
use fossil,     only : surface_stl_object, curve_object, slice_surface, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, &
                       STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_box, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NC = 24
real(R8P),    parameter       :: LIFT = 1.e-2_R8P            ! the section drawn just in front of the cells
type(surface_stl_object)      :: body
type(curve_object)            :: section
type(vector_R8P)              :: bmin, bmax, h, e(4)
type(vector_R8P), allocatable :: centre(:), a(:), b(:)
real(R8P),        allocatable :: d(:), d_host(:), d_pn(:)
real(R8P)                     :: half_diag, dmax, uv(2)
logical,          allocatable :: flip(:)
character(90)                 :: caption
integer(I4P),     allocatable :: jc(:)
integer(I4P)                  :: u, i, j, k, c, kc, m, status

call body%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = body%get_bmin() ; bmax = body%get_bmax() ; h = bmax - bmin
bmin = bmin - 0.1_R8P * h ; bmax = bmax + 0.1_R8P * h
h = (bmax - bmin) / real(NC, R8P)
half_diag = 0.5_R8P * sqrt(h%x**2 + h%y**2 + h%z**2)
allocate(centre(NC**3), d(NC**3), d_host(NC**3), d_pn(NC**3))
c = 0
do k = 1, NC ; do j = 1, NC ; do i = 1, NC
   c = c + 1
   centre(c) = vector_R8P(bmin%x + (i - 0.5_R8P) * h%x, bmin%y + (j - 0.5_R8P) * h%y, bmin%z + (k - 0.5_R8P) * h%z)
enddo ; enddo ; enddo
call body%enter_device
call body%distance_many_device(points=centre, distances=d, is_signed=.true., is_square_root=.true., &
                               sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'distance_many_device failed'
call body%distance_many_device(points=centre, distances=d_pn, is_signed=.true., is_square_root=.true., &
                               sign_algorithm=SIGN_PSEUDO_NORMAL, status=status)
call body%exit_device
call body%distance_many(points=centre, distances=d_host, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
flip = (d_pn < 0._R8P) .neqv. (d < 0._R8P)
allocate(jc(NC**3)) ; jc = [(((j, i = 1, NC), j = 1, NC), k = 1, NC)]   ! the y index of each cell
kc = 1                                                       ! the plane y = const with the most sign disagreements
do j = 2, NC
   if (count(flip .and. jc == j) > count(flip .and. jc == kc)) kc = j
enddo

call section_at(NC / 2)                                      ! the marker: the middle plane
write(caption, '(A,I0,A,I0,A,I0)') 'marker (gold: cut), plane y ', NC / 2, ' of ', NC, ': device vs host mismatches ', &
                                   count(jc == NC / 2 .and. d /= d_host)
call fig_open(u, 'ref-device-distance-marker', trim(caption), -90._R8P, 0._R8P)
do c = 1, NC**3
   if (jc(c) /= NC / 2) cycle
   call corners(centre(c))
   if (abs(d(c)) <= half_diag) then
      call fig_quad(u, 'foot', e(1), e(2), e(3), e(4))
   elseif (d(c) < 0._R8P) then
      call fig_quad(u, 'inside', e(1), e(2), e(3), e(4))
   else
      call fig_quad(u, 'outside', e(1), e(2), e(3), e(4))
   endif
enddo
call fig_segments(u, 'box', a(1:m), b(1:m))
call fig_close(u)

call section_at(kc)                                          ! the signs: the plane of most disagreements
dmax = maxval(abs(d), jc == kc)
write(caption, '(A,I0,A,I0,A)') 'winding-number distance, plane y ', kc, ': ', count(flip .and. jc == kc), &
                                 ' pseudo-normal signs differ (outlined)'
call fig_open(u, 'ref-device-distance-signs', trim(caption), -90._R8P, 0._R8P)
do c = 1, NC**3
   if (jc(c) /= kc) cycle
   call corners(centre(c))
   call fig_quad(u, 'field', e(1), e(2), e(3), e(4), value=max(-1._R8P, min(1._R8P, d(c) / (0.5_R8P * dmax))))
   if (flip(c)) call fig_box(u, 'accent', e(1), e(3))
enddo
call fig_segments(u, 'box', a(1:m), b(1:m))
call fig_close(u)
contains
   subroutine section_at(l)
   !< The dragon's section by the plane of the cell centres y index l, as segments a(1:m)-b(1:m) just in front of it.
   integer(I4P), intent(in) :: l
   real(R8P)                :: yc

   yc = centre((l - 1) * NC + 1)%y
   call slice_surface(body, 2_I4P, yc, section, status=status)  ! (u, v) = (z, x)
   if (allocated(a)) deallocate(a, b)
   allocate(a(section%get_vertices_number()), b(section%get_vertices_number())) ; m = 0
   do i = 1, section%get_vertices_number()
      if (section%get_next(i) == 0_I4P) cycle
      m = m + 1
      uv = section%get_vertex(i) ; a(m) = vector_R8P(uv(2), yc - LIFT, uv(1))
      uv = section%get_vertex(section%get_next(i)) ; b(m) = vector_R8P(uv(2), yc - LIFT, uv(1))
   enddo
   endsubroutine section_at

   subroutine corners(q)
   !< The corners of the cell about the centre q, in its plane y = const, into e.
   type(vector_R8P), intent(in) :: q

   e = q
   e(1)%x = q%x - 0.5_R8P * h%x ; e(1)%z = q%z - 0.5_R8P * h%z
   e(2)%x = q%x + 0.5_R8P * h%x ; e(2)%z = q%z - 0.5_R8P * h%z
   e(3)%x = q%x + 0.5_R8P * h%x ; e(3)%z = q%z + 0.5_R8P * h%z
   e(4)%x = q%x - 0.5_R8P * h%x ; e(4)%z = q%z + 0.5_R8P * h%z
   endsubroutine corners
endprogram fig_ref_device_distance
