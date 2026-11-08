!< Figure of the reference page "Signed-distance field", section "What it does": the torus and the signed distance on
!< the mid plane of a stretched grid, a narrow band of exact values around the surface and the swept far field.
program fig_ref_distance_field_torus
use fossil,     only : surface_stl_object, distance_field_object, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig, only : fig_open, fig_triangle, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter     :: NC = 24, NS = 6
type(surface_stl_object)    :: torus
type(facet_object), allocatable :: f(:)
type(facet_object), pointer :: fp
type(vector_R8P)            :: c
real(R8P), parameter        :: PI = acos(-1._R8P)
type(distance_field_object) :: field
type(vector_R8P)            :: bmin, bmax
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P), allocatable      :: x(:), y(:), z(:)
real(R8P)                   :: phimax, zc, v
integer(I4P)                :: i, j, k, m, n(3), sweeps, fallback, status, u
character(160)              :: caption

call make_torus(f)
call torus%adopt_facets(facets=f) ; call torus%analyze
if (.not. torus%is_watertight()) error stop 'torus not watertight'
print '(A,I0,A)', 'torus: R 1, r 0.35, ', torus%get_facets_number(), ' facets'
bmin = torus%get_bmin() ; bmax = torus%get_bmax()
call axis(bmin%x, bmax%x, x) ; call axis(bmin%y, bmax%y, y) ; call axis(bmin%z, bmax%z, z)
call field%initialize(x=x, y=y, z=z, status=status)
call field%compute(surface=torus, sweeps=sweeps, fallback_cells=fallback, status=status)
if (status /= STATUS_OK) error stop 'compute failed'
n = [size(x), size(y), size(z)]
print '(A,I0,A,I0,A,I0,A,I0,A)', 'grid ', n(1), ' x ', n(2), ' x ', n(3), ' nodes, ', product(n), ' in all'
print '(A,I0,A,ES10.3)', 'band nodes ', field%get_band_cells(), ', band width ', field%get_band_width()
print '(A,I0,A,I0)', 'sweeps ', sweeps, ', nodes signed by their own query ', fallback
k = minloc(abs(z), 1) ; zc = z(k)
phimax = maxval(field%phi(:, :, k))
print '(A,F8.4,A,F8.4,A,F8.4)', 'plane z = ', zc, ': phi from ', minval(field%phi(:, :, k)), ' to ', &
                                maxval(field%phi(:, :, k))

write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'phi on z = 0 of a stretched ', n(1), ' x ', n(2), ' x ', n(3), &
      ' grid: inside teal, outside ember, band (', count(field%in_band(:, :, k)), ' nodes) in gold'
call fig_open(u, 'ref-distance-field-torus', trim(caption), -55._R8P, 32._R8P, width=820_I4P, height=520_I4P)
do i = 1, torus%get_facets_number()                     ! the half of the torus away from the viewer: a cutaway
   fp => torus%facet_at(i)
   c = (fp%vertex(1) + fp%vertex(2) + fp%vertex(3)) / 3._R8P
   if (c%x * cos(-55._R8P * PI / 180) + c%y * sin(-55._R8P * PI / 180) < 0._R8P) &
      call fig_triangle(u, 'body', fp%vertex(1), fp%vertex(2), fp%vertex(3))
enddo
do j = 1, n(2)
   do i = 1, n(1)
      if (field%phi(i, j, k) < 0._R8P) then
         v = -sqrt(field%phi(i, j, k) / minval(field%phi(:, :, k)))
      else
         v = sqrt(field%phi(i, j, k) / phimax)
      endif
      call fig_quad(u, 'field', p(mid(x, i, -1), mid(y, j, -1)), p(mid(x, i, 1), mid(y, j, -1)), &
                    p(mid(x, i, 1), mid(y, j, 1)), p(mid(x, i, -1), mid(y, j, 1)), v)
   enddo
enddo
! the outline of the band nodes on the plane: the dual-cell sides between a band node and a non-band node
allocate(a(4 * n(1) * n(2)), b(4 * n(1) * n(2))) ; m = 0
do j = 1, n(2)
   do i = 1, n(1) - 1
      if (field%in_band(i, j, k) .neqv. field%in_band(i + 1, j, k)) then
         m = m + 1 ; a(m) = p(mid(x, i, 1), mid(y, j, -1), 1.e-3_R8P) ; b(m) = p(mid(x, i, 1), mid(y, j, 1), 1.e-3_R8P)
      endif
   enddo
enddo
do j = 1, n(2) - 1
   do i = 1, n(1)
      if (field%in_band(i, j, k) .neqv. field%in_band(i, j + 1, k)) then
         m = m + 1 ; a(m) = p(mid(x, i, -1), mid(y, j, 1), 1.e-3_R8P) ; b(m) = p(mid(x, i, 1), mid(y, j, 1), 1.e-3_R8P)
      endif
   enddo
enddo
call fig_segments(u, 'accent', a(1:m), b(1:m))
call fig_close(u)
contains
   subroutine make_torus(ff)
   !< A torus about the z axis, tube centre radius 1, tube radius 0.35: 48 x 20 quads, two facets each, outward.
   type(facet_object), allocatable, intent(out) :: ff(:)
   integer(I4P), parameter                      :: NU = 48, NV = 20
   real(R8P),    parameter                      :: PI = acos(-1._R8P)
   integer(I4P)                                 :: iu, iv, mm

   allocate(ff(2 * NU * NV)) ; mm = 0
   do iu = 0, NU - 1
      do iv = 0, NV - 1
         mm = mm + 1 ; call tri(ff(mm), tp(iu, iv), tp(iu + 1, iv), tp(iu + 1, iv + 1), iu)
         mm = mm + 1 ; call tri(ff(mm), tp(iu, iv), tp(iu + 1, iv + 1), tp(iu, iv + 1), iu)
      enddo
   enddo
   endsubroutine make_torus

   pure function tp(iu, iv) result(q)
   !< Vertex (iu, iv) of the torus.
   integer(I4P), intent(in) :: iu, iv
   real(R8P), parameter     :: PI = acos(-1._R8P)
   real(R8P)                :: q(3), t, s

   t = 2 * PI * mod(iu, 48) / 48._R8P ; s = 2 * PI * mod(iv, 20) / 20._R8P
   q = [(1._R8P + 0.35_R8P * cos(s)) * cos(t), (1._R8P + 0.35_R8P * cos(s)) * sin(t), 0.35_R8P * sin(s)]
   endfunction tp

   subroutine tri(ft, p1, p2, p3, iu)
   !< A facet oriented away from the tube centre at the angle of quad column iu.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   integer(I4P),       intent(in)    :: iu
   real(R8P), parameter              :: PI = acos(-1._R8P)
   real(R8P)                         :: e1(3), e2(3), nn(3), t, c(3)

   t = 2 * PI * (iu + 0.5_R8P) / 48._R8P ; c = [cos(t), sin(t), 0._R8P]
   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - c) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri

   subroutine axis(lo, hi, c)
   !< NC + 1 uniform nodes over [lo, hi] widened by 4 %, then NS nodes per side growing by 1.25.
   real(R8P),              intent(in)  :: lo, hi
   real(R8P), allocatable, intent(out) :: c(:)
   real(R8P)                           :: h, l0, l1
   integer(I4P)                        :: ii

   l0 = lo - 0.04_R8P * (hi - lo) ; l1 = hi + 0.04_R8P * (hi - lo) ; h = (l1 - l0) / NC
   c = [(l0 + h * ii, ii = 0, NC)]
   do ii = 1, NS
      h = 1.25_R8P * h
      c = [c(1) - h, c, c(size(c)) + h]
   enddo
   endsubroutine axis

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
   !< The point (px, py) of the plane z = zc, lifted by dz (optional).
   real(R8P), intent(in)           :: px, py
   real(R8P), intent(in), optional :: dz
   type(vector_R8P)                :: q

   q = vector_R8P(px, py, zc)
   if (present(dz)) q%z = zc + dz
   endfunction p
endprogram fig_ref_distance_field_torus
