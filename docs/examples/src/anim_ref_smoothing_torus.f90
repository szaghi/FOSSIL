!run anim-ref-smoothing-torus anim_ref_smoothing_torus
!fig ref-smoothing-torus
!fig anim-ref-smoothing-torus
!< Figures of the reference page "Mesh smoothing" on a torus (radii 0.7 and 0.3, marching cubes with cells of 0.071)
!< whose vertices are moved by deterministic Gaussian noise of sigma 0.03 (the helper of ref_smoothing). "What it does":
!< the noisy torus, then 10 explicit steps (lambda 0.5), then 5 Taubin pairs (lambda 0.5, mu -0.53): the same 10
!< applications of the Laplacian, with the volume of each. The animation: explicit (left) against Taubin (right) after
!< 0, 10, 20, 40, 80 and 160 applications: the explicit torus thins and shrinks, the Taubin one keeps its volume.
program anim_ref_smoothing_torus
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, SMOOTH_METHOD_EXPLICIT, &
                       SMOOTH_METHOD_TAUBIN, SMOOTH_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_label, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 32, NZ = 13, STEPS(6) = [0, 10, 20, 40, 80, 160]
real(R8P),    parameter         :: H = 2.2_R8P / (NG - 1)
type(surface_stl_object)        :: noisy, explicit, taubin, shown
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: right
real(R8P)                       :: phi(NG, NG, NZ), p(3), v0
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status

do k = 1, NZ ; do j = 1, NG ; do i = 1, NG
   p = [-1.1_R8P + H * (i - 1), -1.1_R8P + H * (j - 1), -H * (NZ - 1) / 2 + H * (k - 1)]
   phi(i, j, k) = hypot(hypot(p(1), p(2)) - 0.7_R8P, p(3)) - 0.3_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.1_R8P, -1.1_R8P, -H * (NZ - 1) / 2), &
                        bmax=vector_R8P(1.1_R8P, 1.1_R8P, H * (NZ - 1) / 2), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call noisy%adopt_facets(facets=facets)
call noisy%sanitize(fill_holes=.true.)
allocate(facets(noisy%get_facets_number()))
do i = 1, size(facets)
   f => noisy%facet_at(i)
   facets(i) = f
   do j = 1, 3
      facets(i)%vertex(j) = facets(i)%vertex(j) + 0.03_R8P * noise(facets(i)%vertex(j))
   enddo
   call facets(i)%compute_normal
enddo
call noisy%adopt_facets(facets=facets)
v0 = noisy%get_volume()
right = vector_R8P(2.0_R8P * 0.866_R8P, 2.0_R8P * 0.5_R8P, 0._R8P)

explicit = noisy
call explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, iterations=10_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'explicit smoothing failed'
taubin = noisy
call taubin%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=5_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'Taubin smoothing failed'
write(caption, '(A,F6.4,A,F6.4,A,F6.4)') '10 Laplacian applications; volume: noisy ', v0, ', explicit ', &
      explicit%get_volume(), ', Taubin ', taubin%get_volume()
print '(A)', trim(caption)
call fig_open(u, 'ref-smoothing-torus', trim(caption), -60._R8P, 35._R8P, 760_I4P, 300_I4P)
call place(noisy, 0, 'noisy')
call place(explicit, 1, 'explicit, 10 steps')
call place(taubin, 2, 'Taubin, 5 pairs')
call fig_close(u)

explicit = noisy ; taubin = noisy
call fig_open(u, 'anim-ref-smoothing-torus', '', -60._R8P, 35._R8P, 720_I4P, 340_I4P)
call fig_anim(u, 1.4_R8P)
do k = 1, size(STEPS)
   if (k > 1) then
      call explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, iterations=STEPS(k) - STEPS(k - 1), &
                           status=status)
      if (status /= SMOOTH_STATUS_OK) error stop 'explicit smoothing failed'
      call taubin%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=(STEPS(k) - STEPS(k - 1)) / 2, status=status)
      if (status /= SMOOTH_STATUS_OK) error stop 'Taubin smoothing failed'
   endif
   write(caption, '(I0,A,F6.2,A,F6.2,A)') STEPS(k), ' applications; volume change: explicit ', &
         100 * (explicit%get_volume() / v0 - 1), ' %, Taubin ', 100 * (taubin%get_volume() / v0 - 1), ' %'
   print '(A)', trim(caption)
   if (k == size(STEPS)) then
      call fig_frame(u, caption=trim(caption), seconds=3._R8P)
   else
      call fig_frame(u, caption=trim(caption))
   endif
   call place(explicit, 0, 'explicit')
   call place(taubin, 1, 'Taubin')
enddo
call fig_close(u)
contains
   subroutine place(s, slot, label)
   !< Draw s moved to the slot-th place along the screen's right, with a label under it.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   integer(I4P),             intent(in) :: slot  !< Place.
   character(*),             intent(in) :: label !< Label.

   shown = s
   call shown%translate(delta=real(slot, R8P) * right)
   call fig_surface(u, 'body', shown)
   call fig_label(u, 'accent', vector_R8P(0._R8P, -1.05_R8P, -0.45_R8P) + real(slot, R8P) * right, label)
   endsubroutine place

   function noise(q) result(gd)
   !< Three standard Gaussian deviates (Box-Muller), a deterministic function of the point `q` (as ref_smoothing).
   type(vector_R8P), intent(in) :: q  !< Point.
   type(vector_R8P)             :: gd !< Deviates.
   real(R8P)                    :: r(6)
   integer(I4P)                 :: s

   do s = 1, 6
      r(s) = sin(q%x * (12.9898_R8P + s) + q%y * (78.233_R8P + 3 * s) + q%z * (37.719_R8P + 7 * s)) * 43758.5453_R8P
      r(s) = r(s) - floor(r(s))
   enddo
   r(1:3) = max(r(1:3), 1.e-12_R8P)
   gd%x = sqrt(-2._R8P * log(r(1))) * cos(2._R8P * acos(-1._R8P) * r(4))
   gd%y = sqrt(-2._R8P * log(r(2))) * cos(2._R8P * acos(-1._R8P) * r(5))
   gd%z = sqrt(-2._R8P * log(r(3))) * cos(2._R8P * acos(-1._R8P) * r(6))
   endfunction noise
endprogram anim_ref_smoothing_torus
