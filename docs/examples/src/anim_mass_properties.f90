!run anim-mass-properties anim_mass_properties
!fig anim-mass-properties
!< Animation of the reference page "Mass properties": the tennis-racket flip (Dzhanibekov effect). A T-handle, a closed
!< prism with three distinct principal moments, spins about its intermediate axis with a small perturbation (0.05 rad/s about the other two axes); the spin is
!< unstable and the handle flips over periodically. The program takes the inertia tensor from mass_properties,
!< integrates Euler's equations I w' = (I w) x w in the principal frame (RK4, fixed step) together with the
!< orientation R' = R [w]x, and draws the body at each pose (principal_frame, move) with its angular velocity (ember,
!< from the centre of mass) and the line of its angular momentum, fixed in the world frame (gold). The printed drifts of the kinetic energy and of the
!< angular momentum are the integrator's error.
program anim_mass_properties
use fossil,              only : surface_stl_object, facet_object, rigid_pose_t, principal_frame, pose_compose, &
                                pose_inverse, STATUS_OK
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_surface, fig_segment, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: NSTEP = 4800, EVERY = 50         ! RK4 steps, steps per frame (16 per turn)
real(R8P),    parameter         :: DT = 1.25e-3_R8P                 ! time step [s]
type(surface_stl_object)        :: handle, moved
type(rigid_pose_t)              :: frame, pose
real(R8P)                       :: mass, center(3), pm(3), pa(3,3), w(3), r(3,3), lw(3), lw0(3), e0, de, dl
real(R8P)                       :: k1(3,4), k2(3,4), k3(3,4), k4(3,4), y(3,4)
character(96)                   :: caption
integer(I4P)                    :: u, n, status, flips
real(R8P)                       :: w2_prev

call make_handle(handle)
call handle%mass_properties(mass=mass, center=center, principal_moments=pm, principal_axes=pa, status=status)
if (status /= STATUS_OK) error stop 'the handle is not closed'
print '(A,3F10.6)', 'principal moments: ', pm
frame = principal_frame(center, pa)                      ! principal coordinates -> body coordinates

w = [0.05_R8P, 2._R8P * acos(-1._R8P), 0.05_R8P]        ! one turn per second about the intermediate axis
r = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3])
e0 = 0.5_R8P * sum(pm * w * w)
lw0 = matmul(r, pm * w)
de = 0._R8P ; dl = 0._R8P ; flips = 0 ; w2_prev = w(2)

call fig_open(u, 'anim-mass-properties', '', -30._R8P, 20._R8P)
call fig_anim(u, 0.09_R8P)
do n = 0, NSTEP
   if (mod(n, EVERY) == 0) then
      pose%r = r ; pose%t = 0._R8P                         ! world <- principal: rotate about the centre of mass
      pose = pose_compose(pose, pose_inverse(frame))       ! world <- body
      moved = handle
      call moved%move(pose)
      write(caption, '(A,F5.2,A,3F6.2,A,I0)') 't = ', n * DT, ' s  w = (', w, ') rad/s  flips: ', flips
      call fig_frame(u, caption=trim(caption))
      call fig_surface(u, 'body', moved)
      call fig_segment(u, 'accent', vec(-0.15_R8P * lw0 / pm(2)), vec(0.15_R8P * lw0 / pm(2)))  ! L: fixed in space
      call fig_segment(u, 'patch', vector_R8P(0._R8P, 0._R8P, 0._R8P), vec(0.15_R8P * matmul(r, w)))   ! w
   endif
   if (n == NSTEP) exit
   y(:,1) = w ; y(:,2:4) = r
   k1 = rhs(y) ; k2 = rhs(y + 0.5_R8P * DT * k1) ; k3 = rhs(y + 0.5_R8P * DT * k2) ; k4 = rhs(y + DT * k3)
   y = y + DT / 6._R8P * (k1 + 2._R8P * k2 + 2._R8P * k3 + k4)
   w = y(:,1) ; r = y(:,2:4)
   if (w(2) * w2_prev < 0._R8P) flips = flips + 1
   w2_prev = w(2)
   lw = matmul(r, pm * w)
   de = max(de, abs(0.5_R8P * sum(pm * w * w) - e0) / e0)
   dl = max(dl, norm2(lw - lw0) / norm2(lw0))
enddo
call fig_close(u)
print '(A,I0,A,F5.2,A)', 'flips of the intermediate axis: ', flips, ' in ', NSTEP * DT, ' s'
print '(A,ES9.2,A,ES9.2)', 'max relative drift: energy ', de, ', angular momentum ', dl

contains
   function rhs(s) result(ds)
   !< Euler's equations in the principal frame and the orientation: w' = I^-1 ((I w) x w), R' = R [w]x.
   real(R8P), intent(in) :: s(3,4)
   real(R8P)             :: ds(3,4)
   real(R8P)             :: om(3), lb(3), wx(3,3)

   om = s(:,1) ; lb = pm * om
   ds(:,1) = [lb(2) * om(3) - lb(3) * om(2), lb(3) * om(1) - lb(1) * om(3), lb(1) * om(2) - lb(2) * om(1)] / pm
   wx = reshape([0._R8P, om(3), -om(2), -om(3), 0._R8P, om(1), om(2), -om(1), 0._R8P], [3, 3])
   ds(:,2:4) = matmul(s(:,2:4), wx)
   endfunction rhs

   function vec(x) result(v)
   !< A vector_R8P from an array.
   real(R8P), intent(in) :: x(3)
   type(vector_R8P)      :: v

   v = vector_R8P(x(1), x(2), x(3))
   endfunction vec

   subroutine make_handle(surf)
   !< The T-handle: the T polygon below, counter-clockwise, extruded over z in [-0.08, 0.08], one closed shell.
   type(surface_stl_object), intent(inout) :: surf
   real(R8P), parameter            :: P(2,8) = reshape([-0.1_R8P, -0.6_R8P,  0.1_R8P, -0.6_R8P,  0.1_R8P, 0.3_R8P, &
                                                         0.5_R8P,  0.3_R8P,  0.5_R8P,  0.5_R8P, -0.5_R8P, 0.5_R8P, &
                                                        -0.5_R8P,  0.3_R8P, -0.1_R8P,  0.3_R8P], [2, 8])
   integer,   parameter            :: T(3,6) = reshape([1, 2, 3,  1, 3, 8,  6, 7, 8,  6, 8, 3,  6, 3, 4,  6, 4, 5], [3, 6])
   real(R8P), parameter            :: Z0 = -0.08_R8P, Z1 = 0.08_R8P
   type(facet_object), allocatable :: f(:)
   integer                         :: i, j, m

   allocate(f(2 * size(T, 2) + 2 * size(P, 2)))
   m = 0
   do i = 1, size(T, 2)                                    ! caps: top as wound, bottom reversed
      m = m + 1 ; call tri(f(m), [P(:,T(1,i)), Z1], [P(:,T(2,i)), Z1], [P(:,T(3,i)), Z1])
      m = m + 1 ; call tri(f(m), [P(:,T(1,i)), Z0], [P(:,T(3,i)), Z0], [P(:,T(2,i)), Z0])
   enddo
   do i = 1, size(P, 2)                                    ! sides, outward
      j = mod(i, size(P, 2)) + 1
      m = m + 1 ; call tri(f(m), [P(:,i), Z0], [P(:,j), Z0], [P(:,j), Z1])
      m = m + 1 ; call tri(f(m), [P(:,i), Z0], [P(:,j), Z1], [P(:,i), Z1])
   enddo
   call surf%adopt_facets(facets=f)
   endsubroutine make_handle

   subroutine tri(f, a, b, c)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: a(3), b(3), c(3)

   f%vertex(1) = vec(a) ; f%vertex(2) = vec(b) ; f%vertex(3) = vec(c)
   call f%compute_metrix
   endsubroutine tri
endprogram anim_mass_properties
