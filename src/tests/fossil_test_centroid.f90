!< FOSSIL test: the centroid of a closed surface is the centre of its enclosed volume (issue #59).
!<
!< The former formula weighted each facet by its unit normal where the unnormalised one is needed: the result scaled
!< with the facet size (2, 8, 32, 128 instead of 0.5 on cube-01..04).
!< 1. cube-01..04 (the unit cube at 48 to 3072 facets): centroid (0.5, 0.5, 0.5) to 1e-14.
!< 2. cube-04 translated by (10, -20, 30): the centroid moves with the body.
!< 3. resize(factor=2, respect_centroid=.true.) keeps the centroid of cube-02 in place and scales the volume by 8.
program fossil_test_centroid

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

type(surface_stl_object) :: s
type(vector_R8P)         :: c, c0, d
character(7), parameter  :: CUBES(4) = ['cube-01', 'cube-02', 'cube-03', 'cube-04']
real(R8P)                :: err, v0
integer(I4P)             :: k
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. the unit cube at four resolutions
do k = 1, 4
   call s%load_from_file(file_name='src/tests/'//CUBES(k)//'.stl', guess_format=.true.)
   c = s%get_centroid()
   err = maxval(abs([c%x, c%y, c%z] - 0.5_R8P))
   ok = err <= 1e-14_R8P
   print '(A,A,A,3ES12.4,A,L1)', '1. ', CUBES(k), ' centroid ', c%x, c%y, c%z, '  ok ', ok
   all_ok = all_ok .and. ok
enddo

! ---- 2. moved with the body
d = 10._R8P * ex_R8P - 20._R8P * ey_R8P + 30._R8P * ez_R8P
call s%translate(delta=d)
call s%analyze
c = s%get_centroid()
err = maxval(abs([c%x - 10.5_R8P, c%y + 19.5_R8P, c%z - 30.5_R8P]))
ok = err <= 1e-12_R8P
print '(A,ES10.3,A,L1)', '2. translated cube-04: centroid err ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. resize about the centroid
call s%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
c0 = s%get_centroid() ; v0 = s%get_volume()
call s%resize(factor=2._R8P * (ex_R8P + ey_R8P + ez_R8P), respect_centroid=.true.)
call s%analyze
c = s%get_centroid()
err = maxval(abs([c%x - c0%x, c%y - c0%y, c%z - c0%z]))
ok = err <= 1e-14_R8P .and. abs(s%get_volume() - 8._R8P * v0) <= 1e-14_R8P * 8._R8P * v0
print '(A,ES10.3,A,L1)', '3. resize about the centroid: centroid shift ', err, '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok
endprogram fossil_test_centroid
