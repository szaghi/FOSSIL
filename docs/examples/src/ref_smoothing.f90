!run ref_smoothing ref_smoothing
!fig ref-smoothing-sections
!< Run example of the reference page "Mesh smoothing" and its figure: a marching-cubes unit sphere with deterministic
!< noise (sigma 1 % of the box diagonal, a function of the vertex coordinates), smoothed by 10 explicit Laplacian steps
!< and by 5 Taubin pairs (the same 10 Laplacian applications); the figure is the section z = 0 of the three surfaces,
!< seen from above, its radial error magnified 10 times, with the unit circle (zero error, gold): the explicit section
!< lies inside it all around (shrinkage), the Taubin one on it.
program ref_smoothing
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, SMOOTH_METHOD_EXPLICIT, &
                       SMOOTH_METHOD_TAUBIN, SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT
use fossil_fig, only : fig_open, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
type(surface_stl_object)        :: noisy, explicit, taubin
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
real(R8P)                       :: phi(41, 41, 41)
integer(I4P)                    :: i, j, k, u, status
!region noisy
do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41                     ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call noisy%adopt_facets(facets=facets)
call noisy%sanitize(fill_holes=.true.)
allocate(facets(noisy%get_facets_number()))                    ! every vertex moved by Gaussian noise, sigma 0.035
do i = 1, size(facets)
   f => noisy%facet_at(i)
   facets(i) = f
   do j = 1, 3
      facets(i)%vertex(j) = facets(i)%vertex(j) + 0.035_R8P * noise(facets(i)%vertex(j))
   enddo
   call facets(i)%compute_normal
enddo
call noisy%adopt_facets(facets=facets)
call report('noisy           ', noisy)
!endregion noisy
!region smooth
explicit = noisy
call explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, iterations=10_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'explicit smoothing failed'
call report('explicit, 10    ', explicit)
taubin = noisy
call taubin%smooth(status=status)                              ! Taubin, 5 pairs of lambda 0.5, mu -0.53
if (status /= SMOOTH_STATUS_OK) error stop 'Taubin smoothing failed'
call report('Taubin, 5 pairs ', taubin)
!endregion smooth
!region reject
call taubin%smooth(method=SMOOTH_METHOD_TAUBIN, lambda=0.5_R8P, mu=-0.4_R8P, status=status)   ! mu > -lambda
print '(A,L1)', 'mu = -0.4 > -lambda: SMOOTH_STATUS_BAD_INPUT ', status == SMOOTH_STATUS_BAD_INPUT
!endregion reject
call fig_open(u, 'ref-smoothing-sections', 'section z = 0, radial error x 10: noisy (grey), explicit (ember), Taubin (teal)', &
              0._R8P, 90._R8P)
do i = 1, 128
   call fig_segment(u, 'accent', vector_R8P(cos(2 * PI * (i - 1) / 128), sin(2 * PI * (i - 1) / 128), 0._R8P), &
                    vector_R8P(cos(2 * PI * i / 128), sin(2 * PI * i / 128), 0._R8P))
enddo
call section(noisy, 'box')
call section(explicit, 'hole')
call section(taubin, 'probe')
call fig_close(u)
contains
   subroutine report(label, s)
   !< Radial error and volume of a surface that should be the unit sphere.
   character(*),             intent(in) :: label
   type(surface_stl_object), intent(in) :: s
   real(R8P)                            :: r, e
   integer(I4P)                         :: n, v
   r = 0._R8P ; e = 0._R8P
   do n = 1, s%get_facets_number()
      f => s%facet_at(n)
      do v = 1, 3
         r = r + f%vertex(v)%normL2() ; e = e + (f%vertex(v)%normL2() - 1._R8P)**2
      enddo
   enddo
   print '(A,A,F7.4,A,F7.4,A,F7.4,A,F6.2,A)', label, ': mean |v| ', r / (3 * s%get_facets_number()), &
         ', rms(|v| - 1) ', sqrt(e / (3 * s%get_facets_number())), ', volume ', s%get_volume(), ' (', &
         100 * (s%get_volume() / noisy%get_volume() - 1), ' %)'
   endsubroutine report

   function noise(p) result(g)
   !< Three standard Gaussian deviates, a function of the point (Box-Muller on a coordinate hash): byte-stable.
   type(vector_R8P), intent(in) :: p
   type(vector_R8P)             :: g
   real(R8P)                    :: r(6)
   integer(I4P)                 :: s
   do s = 1, 6
      r(s) = sin(p%x * (12.9898_R8P + s) + p%y * (78.233_R8P + 3 * s) + p%z * (37.719_R8P + 7 * s)) * 43758.5453_R8P
      r(s) = r(s) - floor(r(s))
   enddo
   r(1:3) = max(r(1:3), 1.e-12_R8P)
   g%x = sqrt(-2._R8P * log(r(1))) * cos(2 * PI * r(4))
   g%y = sqrt(-2._R8P * log(r(2))) * cos(2 * PI * r(5))
   g%z = sqrt(-2._R8P * log(r(3))) * cos(2 * PI * r(6))
   endfunction noise

   subroutine section(s, style)
   !< The section of s by the plane z = 0: one segment per facet crossing it.
   type(surface_stl_object), intent(in) :: s
   character(*),             intent(in) :: style
   type(vector_R8P)                     :: q(2), a, b
   integer(I4P)                         :: n, e, m
   do n = 1, s%get_facets_number()
      f => s%facet_at(n)
      m = 0
      do e = 1, 3
         a = f%vertex(e) ; b = f%vertex(mod(e, 3) + 1)
         if ((a%z < 0._R8P .neqv. b%z < 0._R8P) .and. m < 2) then
            m = m + 1 ; q(m) = a + (a%z / (a%z - b%z)) * (b - a)
         endif
      enddo
      if (m == 2) call fig_segment(u, style, magnified(q(1)), magnified(q(2)))
   enddo
   endsubroutine section

   function magnified(p) result(q)
   !< The point at radius 1 + 10 (|p| - 1) in the direction of p: the radial error magnified 10 times.
   type(vector_R8P), intent(in) :: p
   type(vector_R8P)             :: q
   q = ((1._R8P + 10._R8P * (p%normL2() - 1._R8P)) / p%normL2()) * p
   endfunction magnified
endprogram ref_smoothing
