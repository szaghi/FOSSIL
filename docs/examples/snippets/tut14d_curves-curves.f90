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
