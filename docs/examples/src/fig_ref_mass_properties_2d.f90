!run figref-mass-properties-2d fig_ref_mass_properties_2d
!fig ref-mass-properties-2d
!< Figure of the reference page "Mass properties", 2-D section: a NACA 0012 section of unit chord pitched by 12
!< degrees (200 points, cosine spacing, closed trailing edge), with its centroid and its principal axes from
!< curve%mass_properties, drawn with the length of the radius of gyration sqrt(I / A), twice (along the chord) and four
!< times (across it). The long axis follows the chord: the pitch angle is recovered from principal_angle.
program fig_ref_mass_properties_2d
use fossil,     only : curve_object, STATUS_OK
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: NP = 100                      ! points per side
real(R8P),    parameter       :: PI = acos(-1._R8P), PITCH = -12._R8P * PI / 180._R8P, LIFT = 1.e-3_R8P
type(curve_object)            :: section
type(vector_R8P), allocatable :: a(:), b(:)
type(vector_R8P)              :: c, e1, e2
real(R8P)                     :: x(NP + 1), t(NP + 1), uu(2 * NP), vv(2 * NP), area, center(2), sm(3), pmom(2), angle
real(R8P)                     :: g(2), w(2 * NP)
character(96)                 :: caption, text
integer(I4P)                  :: u, i, n, status

x = [(0.5_R8P * (1._R8P - cos(PI * i / NP)), i = 0, NP)]       ! cosine spacing, x from 0 to 1
t = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x ** 2 + 0.2843_R8P * x ** 3 - 0.1036_R8P * x ** 4)
n = 0
do i = NP + 1, 2, -1                                            ! upper side from the trailing edge to the nose
   n = n + 1 ; uu(n) = x(i) ; vv(n) = t(i)
enddo
do i = 1, NP                                                    ! lower side from the nose to the trailing edge
   n = n + 1 ; uu(n) = x(i) ; vv(n) = -t(i)
enddo
uu = uu - 0.25_R8P                                              ! pitch about the quarter chord
w = cos(PITCH) * uu - sin(PITCH) * vv
vv = sin(PITCH) * uu + cos(PITCH) * vv
uu = w
! Upper side leftwards, lower side rightwards: counter-clockwise, the solid on the left.
call section%adopt_loops(u=uu, v=vv, loop_sizes=[2 * NP])
call section%mass_properties(area=area, center=center, second_moments=sm, principal_moments=pmom, &
                             principal_angle=angle, status=status)
if (status /= STATUS_OK) error stop 'the section is not closed'
print '(A,F9.6,A,2F9.5)', 'area ', area, ', centroid', center
print '(A,3ES12.4)', 'second moments I_uu, I_vv, I_uv:', sm
print '(A,2ES12.4,A,F7.2,A)', 'principal values', pmom, ', long axis at ', angle * 180._R8P / PI, ' degrees'

write(caption, '(A,F6.2,A)') 'NACA 0012 pitched by -12 deg: long principal axis at ', angle * 180._R8P / PI, ' deg'
call fig_open(u, 'ref-mass-properties-2d', trim(caption), -90._R8P, 90._R8P, height=360_I4P)
allocate(a(2 * NP), b(2 * NP))
do i = 1, 2 * NP
   a(i) = p(uu(i), vv(i), 0._R8P)
enddo
b = cshift(a, 1)
call fig_segments(u, 'body', a, b)
c = p(center(1), center(2), LIFT)
g = [4._R8P, 2._R8P] * sqrt(pmom / area)                        ! radii of gyration, magnified 4 and 2 times
e1 = cos(angle) * ex_R8P + sin(angle) * ey_R8P                  ! the axis of the larger value (the chord)
e2 = -sin(angle) * ex_R8P + cos(angle) * ey_R8P
call fig_segment(u, 'accent', c - g(2) * e1, c + g(2) * e1)
call fig_segment(u, 'accent', c - g(1) * e2, c + g(1) * e2)
write(text, '(A,ES9.2)') 'principal_moments(2) = ', pmom(2)
call fig_label(u, 'accent', c - g(2) * e1 + 0.03_R8P * e2, trim(text))
write(text, '(A,ES9.2)') 'principal_moments(1) = ', pmom(1)
call fig_label(u, 'accent', c + g(1) * e2, trim(text))
call fig_points(u, 'patch', [c])
call fig_label(u, 'patch', c + 0.02_R8P * e1 - 0.06_R8P * e2, 'centroid')
call fig_close(u)
contains
   function p(pu, pv, dz) result(r)
   !< A point of the plane of the curve, lifted by dz towards the viewer.
   real(R8P), intent(in) :: pu, pv, dz
   type(vector_R8P)      :: r

   r = pu * ex_R8P + pv * ey_R8P + dz * ez_R8P
   endfunction p
endprogram fig_ref_mass_properties_2d
