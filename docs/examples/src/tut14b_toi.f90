!run tut14b tut14b_toi
!< Tutorial 14, run 2: the throw of run 1 by conservative advancement: the pebble's speed bound, the steps one by one,
!< then time_of_impact with its witnesses, plain and directional.
program tut14b_toi
use fossil, only : surface_stl_object, facet_object, rigid_pose_t, pose_identity, pose_from_axis_angle, &
                   motion_bound_3d, CCD_CONTACT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter        :: PI = acos(-1._R8P)
type(surface_stl_object)    :: dragon, pebble
type(rigid_pose_t)          :: p0, p1
type(facet_object), pointer :: f
type(vector_R8P)            :: bmin, bmax, da, db
real(R8P)                   :: scale, mu, tol, d, s, toi, toi_dir, pa(3), pb(3)
real(R8P), allocatable      :: x(:), y(:), z(:)
integer(I4P)                :: i, k, n, fa, fb, steps, pairs, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! the dragon and the pebble of run 1
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze
call pebble%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call pebble%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
call pebble%resize(x=0.04_R8P, y=0.04_R8P, z=0.04_R8P)
call pebble%analyze
p0 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
!region bound
n = 3 * pebble%get_facets_number()
allocate(x(n), y(n), z(n))
do i = 1, pebble%get_facets_number()                                              ! the pebble's vertices (body frame)
   f => pebble%facet_at(i)
   x(3 * i - 2:3 * i) = f%vertex%x ; y(3 * i - 2:3 * i) = f%vertex%y ; z(3 * i - 2:3 * i) = f%vertex%z
enddo
mu = motion_bound_3d(p0, p1, x, y, z)                                             ! |dt| + |theta| max rho
da = dragon%get_bmax() - dragon%get_bmin() ; db = pebble%get_bmax() - pebble%get_bmin()
tol = 1.e-6_R8P * max(da%normL2(), db%normL2())                                   ! the default tolerance
print '(A,F7.4,A,F7.4,A,ES10.3)', 'pebble speed bound mu = ', mu, ' (translation 1.2, spin ', mu - 1.2_R8P, &
                                  '), tolerance ', tol
!endregion bound
!region steps
call dragon%distance_to(pebble, distance=d, pose_other=p0)                        ! s = 0
print '(A,F9.6,A,ES10.3,A,ES10.3)', 'step 0: s = ', 0._R8P, ', gap ', d, ', next step (d - tol/2) / mu = ', &
                                    (d - 0.5_R8P * tol) / mu
do k = 1, 5                                                                       ! stop after k steps: s_k, d(s_k)
   call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, s, status, distance=d, &
                              max_iterations=k)
   if (d > tol) then
      print '(A,I0,A,F9.6,A,ES10.3,A,ES10.3)', 'step ', k, ': s = ', s, ', gap ', d, ', next step (d - tol/2) / mu = ', &
                                               (d - 0.5_R8P * tol) / mu
   else
      print '(A,I0,A,F9.6,A,ES10.3,A)', 'step ', k, ': s = ', s, ', gap ', d, ' <= tolerance: contact'
   endif
enddo
!endregion steps
!region toi
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi, status, distance=d, point_self=pa, &
                           point_other=pb, facet_self=fa, facet_other=fb, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
print '(A,ES10.3,A,I0,A,I0)', '  gap there ', d, ', dragon facet ', fa, ', pebble facet ', fb
print '(A,3F8.4,A,3F8.4)', '  on the dragon', pa, ', on the pebble', pb
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi_dir, status, iterations=steps, &
                           directional=.true., pairs=pairs)
print '(A,F9.6,A,I0,A,I0,A,L1)', 'directional: s = ', toi_dir, ' (', steps, ' steps, ', pairs, &
                                 ' exact pair evaluations), same contact to 1e-6: ', abs(toi_dir - toi) < 1.e-6_R8P
!endregion toi
endprogram tut14b_toi
