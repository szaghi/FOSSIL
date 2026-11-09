!< FOSSIL test: slicing at a pose and placing a curve for line crossings equal moving a copy, bit for bit (issue #94).
module fossil_test_slice_pose_lib
!< Comparison helpers of `fossil_test_slice_pose`.
use fossil, only : curve_object
use penf,   only : I4P, R8P

implicit none
private
public :: same_curve, same_crossings

contains
   function same_curve(a, b) result(same)
   !< True if the curves have the same loops, vertices, successors and patches (exact equality).
   type(curve_object), intent(in) :: a, b !< Curves.
   logical                        :: same !< Result.
   integer(I4P)                   :: k, f1, l1, f2, l2 !< Counter, loop ranges.

   same = a%get_loops_number() == b%get_loops_number() .and. a%get_vertices_number() == b%get_vertices_number() .and. &
          a%get_normal_axis() == b%get_normal_axis()
   if (.not. same) return
   do k = 1_I4P, a%get_loops_number()
      call a%get_loop_range(k, f1, l1) ; call b%get_loop_range(k, f2, l2)
      same = same .and. f1 == f2 .and. l1 == l2
   enddo
   do k = 1_I4P, a%get_vertices_number()
      same = same .and. all(a%get_vertex(k) == b%get_vertex(k)) .and. a%get_next(k) == b%get_next(k)
      if (a%get_next(k) > 0_I4P) same = same .and. a%segment_patch(k) == b%segment_patch(k)
   enddo
   endfunction same_curve

   function same_crossings(a, b) result(same)
   !< True if the boxes and the crossings of 2 x 41 lines across the box (both axes) are identical.
   type(curve_object), intent(in) :: a, b   !< Curves.
   logical                        :: same   !< Result.
   real(R8P),    allocatable      :: ta(:), tb(:) !< Crossings.
   integer(I4P), allocatable      :: ga(:), gb(:) !< Segments.
   logical,      allocatable      :: ea(:), eb(:) !< Entering flags.
   real(R8P)                      :: lo(2), hi(2), pos !< Box, line.
   integer(I4P)                   :: ax, k, na, nb !< Work.

   lo = a%get_bmin() ; hi = a%get_bmax()
   same = all(lo == b%get_bmin()) .and. all(hi == b%get_bmax())
   do ax = 1_I4P, 2_I4P ! lines along u at v = pos (ax = 1), along v at u = pos (ax = 2)
      do k = 0_I4P, 40_I4P
         pos = lo(3 - ax) + (hi(3 - ax) - lo(3 - ax)) * real(k, R8P) / 40._R8P
         call a%line_crossings(ax, pos, lo(ax) - 1._R8P, hi(ax) + 1._R8P, ta, ga, ea, na)
         call b%line_crossings(ax, pos, lo(ax) - 1._R8P, hi(ax) + 1._R8P, tb, gb, eb, nb)
         same = same .and. na == nb
         if (na /= nb) cycle
         same = same .and. all(ta(1:na) == tb(1:nb)) .and. all(ga(1:na) == gb(1:nb)) .and. &
                all(ea(1:na) .eqv. eb(1:nb))
      enddo
   enddo
   endfunction same_crossings
endmodule fossil_test_slice_pose_lib

program fossil_test_slice_pose
!< FOSSIL test: slicing at a pose and placing a curve for line crossings equal moving a copy, bit for bit (issue #94).
!<
!< The swept-cell kernels no longer copy and move the body per Gauss node:
!< - `slice_surface(..., pose=p)` of the body-frame surface must give the curve `slice_surface` gives on a copy moved by
!<   `move(p)`: every vertex, successor, loop and patch identical (compared with ==). The dragon at two generic poses
!<   (rotation about a skew axis plus a translation; three axes, five levels, both `on_plane_low` rules) and cube-02
!<   translated by quarters (vertices exactly on the planes, where the on-plane rule decides).
!< - `place_for_crossings(copy, c, p)` must give the bounding box and the line crossings (parameters, segments,
!<   entering flags) of `c` moved by `move(p)`, also when the copy was placed at another pose before.
use fossil,                     only : surface_stl_object, curve_object, slice_surface, rigid_pose_t, rigid_pose_2d_t, &
                                       pose_from_axis_angle, pose_2d
use fossil_curve_object,        only : place_for_crossings
use fossil_test_slice_pose_lib, only : same_curve, same_crossings
use penf,                       only : I4P, R8P
use vecfor,                     only : vector_R8P

implicit none

type(surface_stl_object) :: dragon, cube, moved
type(curve_object)       :: c, cm, cp
type(rigid_pose_t)       :: poses(3)
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: lo(3), hi(3), lev
integer(I4P)             :: ip, ax, k, il, nsec, nseg
logical                  :: ok(3), low

call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
poses(1) = pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.8_R8P], 0.7_R8P, [0.11_R8P, -0.23_R8P, 0.05_R8P])
poses(2) = pose_from_axis_angle([-1._R8P, 0.2_R8P, 0.1_R8P], 2.9_R8P, [-0.4_R8P, 0.0_R8P, 1.3_R8P])

! the dragon at two generic poses
ok(1) = .true. ; nsec = 0_I4P ; nseg = 0_I4P
do ip = 1, 2
   moved = dragon ; call moved%move(poses(ip))
   bmin = moved%get_bmin() ; bmax = moved%get_bmax()
   lo = [bmin%x, bmin%y, bmin%z] ; hi = [bmax%x, bmax%y, bmax%z]
   do ax = 1_I4P, 3_I4P
      do k = 1, 5
         lev = lo(ax) + (hi(ax) - lo(ax)) * real(k, R8P) / 6._R8P
         do il = 0, 1
            low = il == 1
            call slice_surface(surface=moved, normal_axis=ax, level=lev, curve=cm, on_plane_low=low)
            call slice_surface(surface=dragon, normal_axis=ax, level=lev, curve=cp, on_plane_low=low, pose=poses(ip))
            ok(1) = ok(1) .and. same_curve(cm, cp)
            nsec = nsec + 1_I4P ; nseg = nseg + cm%get_vertices_number()
         enddo
      enddo
   enddo
enddo
ok(1) = ok(1) .and. nseg > 100_I4P ! the sections are not empty
print '(A,L1,A,I0,A,I0,A)', 'dragon, slice at a pose == slice of the moved copy: ', ok(1), ' (', nsec, ' sections, ', &
   nseg, ' vertices)'

! cube-02 translated by quarters: vertices exactly on the planes
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
poses(3) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [0.25_R8P, -0.5_R8P, 0.75_R8P])
moved = cube ; call moved%move(poses(3))
ok(2) = .true.
do ax = 1_I4P, 3_I4P
   do k = -4, 8
      lev = 0.25_R8P * real(k, R8P)
      do il = 0, 1
         low = il == 1
         call slice_surface(surface=moved, normal_axis=ax, level=lev, curve=cm, on_plane_low=low)
         call slice_surface(surface=cube, normal_axis=ax, level=lev, curve=cp, on_plane_low=low, pose=poses(3))
         ok(2) = ok(2) .and. same_curve(cm, cp)
      enddo
   enddo
enddo
print '(A,L1)', 'cube-02 by quarters (vertices on the planes), slice at a pose == slice of the moved copy: ', ok(2)

! a curve placed for the line crossings == the curve moved
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call slice_surface(surface=dragon, normal_axis=3_I4P, level=0.5_R8P * (bmin%z + bmax%z), curve=c)
ok(3) = c%get_vertices_number() > 20_I4P
cp = c
call place_for_crossings(cp, c, pose_2d(-2.1_R8P, [3._R8P, 1._R8P])) ! somewhere else first
do ip = 1, 3
   cm = c ; call cm%move(pose_2d(0.37_R8P * real(ip, R8P), [0.1_R8P, -0.2_R8P] * real(ip, R8P)))
   call place_for_crossings(cp, c, pose_2d(0.37_R8P * real(ip, R8P), [0.1_R8P, -0.2_R8P] * real(ip, R8P)))
   ok(3) = ok(3) .and. same_crossings(cm, cp)
enddo
print '(A,L1)', 'curve placed for crossings == curve moved (box, crossings of 82 lines, 3 poses): ', ok(3)

print '(A,L1)', 'Are all tests passed? ', all(ok)
endprogram fossil_test_slice_pose
