!run figref-thickness-2d fig_ref_thickness_2d
!fig ref-thickness-2d
!< Figure of the reference page "Thickness and gaps", 2-D: a NACA 0012 section (100 points per side), curve%thickness
!< with the medial method, the medial circle of every 8th segment (centre at the midpoint minus r n, gold), and
!< thin_regions for h = chord / 32: the segments thinner than h in ember, the trailing-edge strip.
program fig_ref_thickness_2d
use fossil,     only : curve_object, STATUS_OK, THICKNESS_MEDIAL
use fossil_fig, only : fig_open, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter          :: PI = acos(-1._R8P)
integer(I4P), parameter       :: NP = 100
type(curve_object)            :: c
type(vector_R8P), allocatable :: a(:), b(:), ta(:), tb(:)
real(R8P), allocatable        :: u(:), v(:), th(:), rl(:)
integer(I4P), allocatable     :: region(:)
real(R8P)                     :: xs(0:NP), ys(0:NP), n(2), m(2), le, r
integer(I4P)                  :: i, j, k, nn, status, uu
character(96)                 :: caption

do i = 0, NP
   xs(i) = 0.5_R8P * (1._R8P - cos(PI * i / NP))
   ys(i) = 0.6_R8P * (0.2969_R8P * sqrt(xs(i)) - 0.1260_R8P * xs(i) - 0.3516_R8P * xs(i) ** 2 + 0.2843_R8P * xs(i) ** 3 - &
                      0.1036_R8P * xs(i) ** 4)
enddo
ys(NP) = 0._R8P
allocate(u(2 * NP), v(2 * NP)) ; nn = 0
do i = NP, 1, -1
   nn = nn + 1 ; u(nn) = xs(i) ; v(nn) = ys(i)
enddo
do i = 0, NP - 1
   nn = nn + 1 ; u(nn) = xs(i) ; v(nn) = -ys(i)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[nn])
call c%thickness(thickness=th, method=THICKNESS_MEDIAL, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call c%thin_regions(thickness=th, h=1._R8P / 32, factor=1._R8P, region=region, region_length=rl)
allocate(a(nn), b(nn))
do i = 1, nn
   j = mod(i, nn) + 1
   a(i) = vector_R8P(u(i), v(i), 0._R8P) ; b(i) = vector_R8P(u(j), v(j), 0._R8P)
enddo
ta = pack(a, region > 0) ; tb = pack(b, region > 0)
write(caption, '(A,I0,A,F6.4,A)') 'medial circles; under h = chord/32: ', count(region > 0), ' segments, ', sum(rl), &
                                  ' long'
print '(A)', trim(caption)
call fig_open(uu, 'ref-thickness-2d', trim(caption), -90._R8P, 90._R8P, width=900_I4P, height=320_I4P)
call fig_segments(uu, 'body', a, b)
do i = 4, nn, 8                                          ! a medial circle every 8 segments
   j = mod(i, nn) + 1
   le = norm2([u(j) - u(i), v(j) - v(i)])
   n = [v(j) - v(i), -(u(j) - u(i))] / le                ! outward normal (solid on the left)
   m = 0.5_R8P * [u(i) + u(j), v(i) + v(j)] ; r = 0.5_R8P * th(i)
   call circle(m - r * n, r)
enddo
call fig_segments(uu, 'hole', ta, tb)
call fig_close(uu)
contains
   subroutine circle(cc, rr)
   !< A circle of 48 segments.
   real(R8P), intent(in) :: cc(2), rr
   type(vector_R8P)      :: ca(48), cb(48)

   do k = 1, 48
      ca(k) = vector_R8P(cc(1) + rr * cos(2 * PI * (k - 1) / 48), cc(2) + rr * sin(2 * PI * (k - 1) / 48), 0._R8P)
      cb(k) = vector_R8P(cc(1) + rr * cos(2 * PI * k / 48), cc(2) + rr * sin(2 * PI * k / 48), 0._R8P)
   enddo
   call fig_segments(uu, 'accent', ca, cb)
   endsubroutine circle
endprogram fig_ref_thickness_2d
