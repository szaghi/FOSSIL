call naca%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(naca, 3_I4P, 0.5_R8P, wing)                                  ! the section, chord 1, in (x, y)
call make_disc(0.03_R8P, disc)                                                   ! a particle of radius 0.03
curves(1) = wing ; curves(2:4) = disc
curve_poses0(1) = pose_2d(angle=0._R8P) ; curve_poses1(1) = curve_poses0(1)       ! the wing at rest
curve_poses0(2) = pose_2d(angle=0._R8P, translation=[-0.2_R8P, 0.25_R8P])          ! falls onto the upper side
curve_poses1(2) = pose_2d(angle=PI, translation=[0.4_R8P, -0.25_R8P])
curve_poses0(3) = pose_2d(angle=0._R8P, translation=[-0.4_R8P, -0.12_R8P])         ! passes below
curve_poses1(3) = pose_2d(angle=PI, translation=[0.9_R8P, -0.12_R8P])
curve_poses0(4) = pose_2d(angle=0._R8P, translation=[-0.45_R8P, 0.2_R8P])          ! passes ahead
curve_poses1(4) = pose_2d(angle=PI, translation=[-0.45_R8P, -0.2_R8P])
call wing%distance_to(disc, distance=d, point_self=qa, segment_self=sa, segment_other=sb, &
                      pose_other=curve_poses0(2))
print '(A,F8.5,A,I0,A,2F8.4)', '2-D gap at s = 0 ', d, ', wing segment ', sa, ' at', qa
call wing%time_of_impact(disc, curve_poses0(1), curve_poses1(1), curve_poses0(2), curve_poses1(2), toi, status, &
                         iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', '2-D contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
call first_impact_2d(curves, curve_poses0, curve_poses1, toi, pair, status, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', '2-D first contact ', status == CCD_CONTACT, ' at s = ', toi, ', curves', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
