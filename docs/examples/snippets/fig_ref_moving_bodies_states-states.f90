p0 = pose_2d(angle=0._R8P)
p1 = pose_2d(angle=3._R8P * PI / 180._R8P, translation=[-0.02_R8P, 0.01_R8P])
call foil%state_changes(pose0=p0, pose1=p1, points=p, fresh=fresh, dead=dead, checked=checked)
