!< FOSSIL test: targets of the constrained projection (issue #74, step 2).
!<
!< 1. The unit cube with one patch per face: 6 non-empty patch trees (the default patch is left empty), 12 constraint
!<    curves (one segment each, length 1), 8 corners;
!<    with one patch: the same 12 curves and 8 corners from the sharp edges alone (the patch boundaries add nothing).
!< 2. A cylinder (32 sides) with patches side / bottom / top: 3 patch trees, 2 closed rim curves of 32 segments, no
!<    corner (each rim turns by 11.25 degrees per vertex); with one patch the same (the rims are sharp, 90 degrees).
!< 3. A UV sphere (16 x 32) with patches south / north: one closed curve, the equator (32 segments), no corner; with one
!<    patch: no curve (every dihedral below 30 degrees).
!< 4. 2-D: a square with a vertex at each side's midpoint, patches changing at two midpoints: 4 corners plus the 2
!<    junctions; a circle of 64 segments: none; an open polyline: its 2 ends.
program fossil_test_projection_targets

use fossil,                   only : surface_stl_object, curve_object, projection_target_t
use fossil_projection_bodies, only : cube_patches, cylinder, uv_sphere
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: s
type(curve_object)        :: c
type(projection_target_t) :: t
real(R8P)                 :: cu(64), cv(64)
integer(I4P)              :: k
logical                   :: all_ok

all_ok = .true.

! 1. cube
call cube_patches(s, .true.) ; call s%projection_target(t)
call check(count(t%patch_root > 0) == 6 .and. t%nc == 12 .and. t%nk == 8 .and. all(seg_lengths(t) == 1._R8P), &
           'cube, 6 patches: 6 trees, 12 curves of one unit segment, 8 corners (curves)', real(t%nc, R8P))
call cube_patches(s, .false.) ; call s%projection_target(t)
call check(t%np == 1 .and. t%nc == 12 .and. t%nk == 8, 'cube, 1 patch: 12 curves and 8 corners from the sharp edges', &
           real(t%nk, R8P))

! 2. cylinder
call cylinder(s, 32_I4P, 1._R8P, 2._R8P, .true.) ; call s%projection_target(t)
call check(count(t%patch_root > 0) == 3 .and. t%nc == 2 .and. t%nk == 0 .and. all(curve_sizes(t) == 32), &
           'cylinder, 3 patches: 2 closed rims of 32 segments, no corner (curves)', real(t%nc, R8P))
call cylinder(s, 32_I4P, 1._R8P, 2._R8P, .false.) ; call s%projection_target(t)
call check(t%np == 1 .and. t%nc == 2 .and. t%nk == 0, 'cylinder, 1 patch: the same rims (sharp)', real(t%nc, R8P))

! 3. sphere
call uv_sphere(s, 16_I4P, 32_I4P, 1._R8P, .true.) ; call s%projection_target(t)
call check(count(t%patch_root > 0) == 2 .and. t%nc == 1 .and. t%nk == 0 .and. all(curve_sizes(t) == 32) .and. &
           maxval(abs(t%seg([3, 6], :))) <= 1e-15_R8P, 'sphere, 2 patches: the equator, closed, 32 segments, no corner', &
           real(t%nc, R8P))
call uv_sphere(s, 16_I4P, 32_I4P, 1._R8P, .false.) ; call s%projection_target(t)
call check(t%nc == 0 .and. t%nk == 0, 'sphere, 1 patch: no curve, no corner', real(t%nc, R8P))

! 4. 2-D
call c%adopt_loops(u=[0._R8P, 0.5_R8P, 1._R8P, 1._R8P, 1._R8P, 0.5_R8P, 0._R8P, 0._R8P], &
                   v=[0._R8P, 0._R8P, 0._R8P, 0.5_R8P, 1._R8P, 1._R8P, 1._R8P, 0.5_R8P], loop_sizes=[8_I4P])
call c%set_patch([2_I4P, 3_I4P, 4_I4P, 5_I4P], 'b')
call c%projection_target(t)
call check(count(t%patch_root > 0) == 2 .and. t%nk == 6 .and. t%ne == 8, '2-D square, 2 patches: 4 corners + 2 junctions', &
           real(t%nk, R8P))
do k = 1, 64
   cu(k) = cos(2 * PI * (k - 1) / 64) ; cv(k) = sin(2 * PI * (k - 1) / 64)
enddo
call c%adopt_loops(u=cu, v=cv, loop_sizes=[64_I4P]) ; call c%projection_target(t)
call check(t%nk == 0 .and. t%ne == 64, '2-D circle: no corner', real(t%nk, R8P))
call c%adopt_loops(u=cu(1:10), v=cv(1:10), loop_sizes=[10_I4P], closed=[.false.]) ; call c%projection_target(t)
call check(t%nk == 2 .and. t%ne == 9, '2-D open polyline: its 2 ends', real(t%nk, R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES23.15)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   function seg_lengths(tt) result(l)
   !< Lengths of the constraint segments.
   type(projection_target_t), intent(in) :: tt   !< Target.
   real(R8P), allocatable                :: l(:) !< Lengths.
   integer(I4P)                          :: i

   allocate(l(size(tt%seg, 2)))
   do i = 1, size(l)
      l(i) = norm2(tt%seg(4:6, i) - tt%seg(1:3, i))
   enddo
   endfunction seg_lengths

   function curve_sizes(tt) result(n)
   !< Segments per curve.
   type(projection_target_t), intent(in) :: tt   !< Target.
   integer(I4P), allocatable             :: n(:) !< Sizes.

   n = tt%curve_first(2:tt%nc + 1) - tt%curve_first(1:tt%nc)
   endfunction curve_sizes
endprogram fossil_test_projection_targets
