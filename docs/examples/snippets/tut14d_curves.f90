!< Tutorial 14, run 4: the 2-D twins around the NACA 0012 section (distance_to, time_of_impact, first_impact_2d), then a
!< particle sliding along a wall: the plain bound creeps, the directional bound and the default (adaptive) cross the step.
program tut14d_curves
use fossil, only : surface_stl_object, curve_object, rigid_pose_2d_t, pose_2d, first_impact_2d, slice_surface, &
                   CCD_CONTACT, CCD_SEPARATE, CCD_NOT_CONVERGED
use penf,   only : I4P, R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: naca
type(curve_object)       :: wing, wall, disc, curves(4)
type(rigid_pose_2d_t)    :: poses0(4), poses1(4), q0, q1
real(R8P)                :: d, toi, qa(2), gap
integer(I4P)             :: k, steps, status, pair(2), tested, sa, sb

call naca%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(naca, 3_I4P, 0.5_R8P, wing)                                  ! the section, chord 1, in (x, y)
call make_disc(0.03_R8P, disc)                                                   ! a particle of radius 0.03
curves(1) = wing ; curves(2:4) = disc
poses0(1) = pose_2d(angle=0._R8P) ; poses1(1) = poses0(1)                        ! the wing at rest
poses0(2) = pose_2d(angle=0._R8P, translation=[-0.2_R8P, 0.25_R8P])              ! falls onto the upper side
poses1(2) = pose_2d(angle=PI, translation=[0.4_R8P, -0.25_R8P])
poses0(3) = pose_2d(angle=0._R8P, translation=[-0.4_R8P, -0.12_R8P])             ! passes below
poses1(3) = pose_2d(angle=PI, translation=[0.9_R8P, -0.12_R8P])
poses0(4) = pose_2d(angle=0._R8P, translation=[-0.45_R8P, 0.2_R8P])              ! passes ahead
poses1(4) = pose_2d(angle=PI, translation=[-0.45_R8P, -0.2_R8P])
call wing%distance_to(disc, distance=d, point_self=qa, segment_self=sa, segment_other=sb, pose_other=poses0(2))
print '(A,F8.5,A,I0,A,2F8.4)', 'gap at s = 0 ', d, ', wing segment ', sa, ' at', qa
call wing%time_of_impact(disc, poses0(1), poses1(1), poses0(2), poses1(2), toi, status, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'particle 2: contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
call first_impact_2d(curves, poses0, poses1, toi, pair, status, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', 'first contact ', status == CCD_CONTACT, ' at s = ', toi, ', curves', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
call wall%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[-0.1_R8P, -0.1_R8P, 0._R8P, 0._R8P], loop_sizes=[4_I4P])
do k = 1, 2                                                                      ! sliding 0.8 along the wall
   gap = merge(0.02_R8P, 4.e-6_R8P, k == 1)                                      ! at a gap of 0.02, then 4 tolerances
   q0 = pose_2d(angle=0._R8P, translation=[0.1_R8P, 0.03_R8P + gap])
   q1 = pose_2d(angle=0._R8P, translation=[0.9_R8P, 0.03_R8P + gap])
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps, directional=.false.)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', plain:       ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', default:     ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps, directional=.true.)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', directional: ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
enddo
contains
   function status_name(st) result(name)
   !< The name of a CCD_* status.
   integer(I4P), intent(in) :: st
   character(17)            :: name

   name = 'CCD_CONTACT'
   if (st == CCD_SEPARATE) name = 'CCD_SEPARATE'
   if (st == CCD_NOT_CONVERGED) name = 'CCD_NOT_CONVERGED'
   endfunction status_name

   subroutine make_disc(radius, c)
   !< A regular 32-gon of the given radius about the origin, counter-clockwise (solid inside).
   real(R8P),          intent(in)    :: radius
   type(curve_object), intent(inout) :: c
   real(R8P)                         :: a(32)
   integer(I4P)                      :: m

   a = [(2._R8P * PI * m / 32._R8P, m = 0, 31)]
   call c%adopt_loops(u=radius * cos(a), v=radius * sin(a), loop_sizes=[32_I4P])
   endsubroutine make_disc
endprogram tut14d_curves
