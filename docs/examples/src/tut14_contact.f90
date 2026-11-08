!run tut14 tut14_contact
!< Tutorial 14: contact between moving bodies: a pebble thrown through the dragon within one step, the time of first
!< contact (plain and directional bounds), the first contact among several bodies, and the 2-D twins.
program tut14_contact
!region setup
use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, pose_2d, &
                   pose_interpolate, first_impact, first_impact_2d, slice_surface, CCD_CONTACT, CCD_SEPARATE, &
                   CCD_NOT_CONVERGED
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter          :: PI = acos(-1._R8P)
type(surface_stl_object)      :: dragon, pebble, naca, bodies(4)
type(curve_object)            :: wing, wall, disc, curves(4)
type(rigid_pose_t)            :: p0, p1, poses0(4), poses1(4)
type(rigid_pose_2d_t)         :: q0, q1, curve_poses0(4), curve_poses1(4)
type(vector_R8P)              :: bmin, bmax
real(R8P)                     :: scale, d, toi, toi_all, pa(3), pb(3), qa(2), s
integer(I4P)                  :: i, j, k, fa, fb, steps, pairs, status, pair(2), tested, sa, sb
logical                       :: hit

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! the dragon, placed as in chapter 3
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze

call pebble%load_from_file(file_name='cube-01.stl', guess_format=.true.)         ! a pebble: a cube of side 0.04
call pebble%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)                       ! centred on its own origin
call pebble%resize(x=0.04_R8P, y=0.04_R8P, z=0.04_R8P)
call pebble%analyze
!endregion setup
!region sampling
! one step: thrown along x through the dragon's body, turning a quarter turn about y
p0 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
do k = 0, 4                                                                      ! the gap at five instants
   s = 0.25_R8P * k
   call dragon%distance_to(pebble, distance=d, intersecting=hit, pose_other=pose_interpolate(p0, p1, s))
   print '(A,F5.2,A,F8.5,A,L1)', 's = ', s, ': gap ', d, ', touching ', hit
enddo
!endregion sampling
!region toi
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi, status, distance=d, point_self=pa, &
                           point_other=pb, facet_self=fa, facet_other=fb, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
print '(A,ES10.3,A,I0,A,I0)', '  gap there ', d, ', dragon facet ', fa, ', pebble facet ', fb
print '(A,3F8.4)', '  on the dragon', pa
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi_all, status, iterations=steps, &
                           directional=.true., pairs=pairs)
print '(A,F9.6,A,I0,A,I0,A,L1)', 'directional: s = ', toi_all, ' (', steps, ' steps, ', pairs, &
                                 ' exact pair evaluations), same contact to 1e-6: ', abs(toi_all - toi) < 1.e-6_R8P
!endregion toi
!region several
bodies(1) = dragon ; poses0(1) = pose_identity() ; poses1(1) = poses0(1)       ! the dragon at rest
bodies(2:4) = pebble                                                             ! three pebbles, three throws
poses0(2) = p0 ; poses1(2) = p1
poses0(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0._R8P, 0.3_R8P, 0.6_R8P])
poses1(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI, translation=[0._R8P, 0.3_R8P, -0.6_R8P])
poses0(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-0.15_R8P, 0.3_R8P, 0.45_R8P])
poses1(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI, translation=[1.05_R8P, 0.3_R8P, 0.45_R8P])
call first_impact(bodies, poses0, poses1, toi, pair, status, tolerance=1.e-6_R8P, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', 'first contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
toi_all = 1._R8P
do i = 1, 3
   do j = i + 1, 4
      call bodies(i)%time_of_impact(bodies(j), poses0(i), poses1(i), poses0(j), poses1(j), s, status, &
                                     tolerance=1.e-6_R8P)
      if (status == CCD_CONTACT) then
         toi_all = min(toi_all, s)
         print '(A,I0,A,I0,A,F9.6)', '  pair ', i, '-', j, ': contact at s = ', s
      else
         print '(A,I0,A,I0,A)', '  pair ', i, '-', j, ': separate'
      endif
   enddo
enddo
print '(A,L1)', 'equal to the loop over all pairs, bit for bit: ', toi_all == toi
!endregion several
!region curves
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
!endregion curves
!region sliding
! a particle sliding along a flat wall at a gap of 4 tolerances
call wall%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[-0.1_R8P, -0.1_R8P, 0._R8P, 0._R8P], loop_sizes=[4_I4P])
q0 = pose_2d(angle=0._R8P, translation=[0.1_R8P, 0.03_R8P + 4.e-6_R8P])
q1 = pose_2d(angle=0._R8P, translation=[0.9_R8P, 0.03_R8P + 4.e-6_R8P])
call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                         tolerance=1.e-6_R8P, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'sliding, plain: not converged ', status == CCD_NOT_CONVERGED, ', s >= ', toi, ' (', &
                              steps, ' steps)'
call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                         tolerance=1.e-6_R8P, iterations=steps, directional=.true.)
print '(A,L1,A,F9.6,A,I0,A)', 'sliding, directional: separate ', status == CCD_SEPARATE, ', s = ', toi, ' (', &
                              steps, ' steps)'
!endregion sliding
contains
   subroutine make_disc(radius, c)
   !< A regular 32-gon of the given radius about the origin, counter-clockwise (solid inside).
   real(R8P),          intent(in)    :: radius
   type(curve_object), intent(inout) :: c
   real(R8P)                         :: a(32)
   integer(I4P)                      :: m

   a = [(2._R8P * PI * m / 32._R8P, m = 0, 31)]
   call c%adopt_loops(u=radius * cos(a), v=radius * sin(a), loop_sizes=[32_I4P])
   endsubroutine make_disc
endprogram tut14_contact
