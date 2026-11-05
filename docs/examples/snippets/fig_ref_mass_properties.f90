!< Figure of the reference page "Mass properties": a closed body with no symmetry (marching cubes of the union of a
!< tilted ellipsoid and an off-centre ball), its centre of mass, its principal axes and the equivalent inertia
!< ellipsoid (the uniform ellipsoid with the same mass and inertia tensor), drawn as its three principal sections. The
!< longest semi-axis lies along the axis of the smallest moment.
program fig_ref_mass_properties
use fossil,              only : surface_stl_object, facet_object, extract_isosurface, STATUS_OK
use fossil_fig,          only : fig_open, fig_surface, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: N = 30, NS = 72
real(R8P),    parameter         :: PI = acos(-1._R8P)
type(surface_stl_object)        :: body
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: values(:,:,:)
type(vector_R8P)                :: c, a(NS), b(NS), p
real(R8P)                       :: lo(3), hi(3), h(3), x(3), r(3,3), q(3)
real(R8P)                       :: mass, center(3), pm(3), pa(3,3), semi(3), t
character(96)                   :: caption, text
integer(I4P)                    :: u, i, j, k, axis, i1, i2, status

! A tilted ellipsoid (semi-axes 1, 0.55, 0.4, rotated 30 degrees about z then 20 about x) and a ball of radius 0.45.
r = reshape([cos(PI / 6), sin(PI / 6), 0._R8P, -sin(PI / 6), cos(PI / 6), 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3])
r = matmul(reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, cos(PI / 9), sin(PI / 9), 0._R8P, -sin(PI / 9), cos(PI / 9)], &
                   [3, 3]), r)
lo = [-1.3_R8P, -1.0_R8P, -0.9_R8P] ; hi = [1.6_R8P, 1.3_R8P, 1.0_R8P]
h = (hi - lo) / (N - 1)
allocate(values(N, N, N))
do k = 1, N
   do j = 1, N
      do i = 1, N
         x = lo + h * [i - 1, j - 1, k - 1]
         q = matmul(transpose(r), x) / [1._R8P, 0.55_R8P, 0.4_R8P]
         values(i,j,k) = min(norm2(q) - 1._R8P, (norm2(x - [0.95_R8P, 0.6_R8P, 0.25_R8P]) - 0.45_R8P) / 0.45_R8P)
      enddo
   enddo
enddo
call extract_isosurface(values=values, bmin=vector_R8P(lo(1), lo(2), lo(3)), bmax=vector_R8P(hi(1), hi(2), hi(3)), &
                        iso=0._R8P, surface=facets, status=status)
call body%adopt_facets(facets=facets)

call body%mass_properties(mass=mass, center=center, principal_moments=pm, principal_axes=pa, status=status)
if (status /= STATUS_OK) error stop 'mass_properties failed'
! Equivalent uniform ellipsoid: I_k = m (s_j^2 + s_l^2) / 5, so s_k^2 = 5 (I_j + I_l - I_k) / (2 m).
semi = sqrt(5._R8P * ([pm(2) + pm(3) - pm(1), pm(1) + pm(3) - pm(2), pm(1) + pm(2) - pm(3)]) / (2._R8P * mass))

write(caption, '(A,I0,A,3F7.3)') 'mass_properties: ', body%get_facets_number(), ' facets; principal moments', pm
call fig_open(u, 'ref-mass-properties', trim(caption), -20._R8P, 50._R8P)
call fig_surface(u, 'ghost', body)
c = vector_R8P(center(1), center(2), center(3))
do axis = 1, 3                                   ! principal sections of the equivalent ellipsoid
   i1 = mod(axis, 3) + 1 ; i2 = mod(axis + 1, 3) + 1
   do k = 1, NS
      t = 2 * PI * (k - 1) / NS
      x = center + semi(i1) * cos(t) * pa(:, i1) + semi(i2) * sin(t) * pa(:, i2)
      a(k) = vector_R8P(x(1), x(2), x(3))
   enddo
   b = cshift(a, 1)
   call fig_segments(u, 'box', a, b)
enddo
do axis = 1, 3                                   ! principal axes, through the centre, to the ellipsoid
   x = center + semi(axis) * pa(:, axis)
   p = vector_R8P(x(1), x(2), x(3))
   x = center - semi(axis) * pa(:, axis)
   call fig_segment(u, 'accent', vector_R8P(x(1), x(2), x(3)), p)
   write(text, '(A,I0,A,F6.3)') 'I', axis, ' = ', pm(axis)
   if (axis == 2) then                           ! label the axes at alternating ends, so the labels do not meet
      call fig_label(u, 'accent', p, trim(text))
   else
      call fig_label(u, 'accent', vector_R8P(x(1), x(2), x(3)), trim(text))
   endif
enddo
call fig_points(u, 'patch', [c])
call fig_label(u, 'patch', c, 'centre of mass')
call fig_close(u)
print '(A,ES12.5)',    'mass (density 1): ', mass
print '(A,3F9.5)',     'centre of mass:   ', center
print '(A,3F9.5)',     'principal moments:', pm
print '(A,3F9.5)',     'ellipsoid axes:   ', semi
endprogram fig_ref_mass_properties
