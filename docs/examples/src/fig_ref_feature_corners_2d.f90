!run figref-feature-corners-2d fig_ref_feature_corners_2d
!fig ref-feature-corners-2d
!< Figure of the reference page "Feature edges", 2-D section: an L-shaped section (counter-clockwise), its corners from
!< curve%feature_corners (gold), and curve%tag_block(feature_band=2, feature_tags=) on 32 x 32 cells: the cells holding a
!< corner and the two rings around them shaded, darkest at the corners.
program fig_ref_feature_corners_2d
use fossil,     only : curve_object, STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: NC = 32
real(R8P),    parameter       :: LU(6) = [0._R8P, 2._R8P, 2._R8P, 1._R8P, 1._R8P, 0._R8P]
real(R8P),    parameter       :: LV(6) = [0._R8P, 0._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P]
type(curve_object)            :: c
type(vector_R8P)              :: a(6), b(6)
type(vector_R8P), allocatable :: pc(:)
integer(I4P), allocatable     :: corners(:)
integer(I4P)                  :: tags(NC, NC), ft(NC, NC), u, i, j, status
real(R8P)                     :: x(NC + 1), y(NC + 1)

call c%adopt_loops(u=LU, v=LV, loop_sizes=[6_I4P])
call c%feature_corners(corners, status=status)
x = [(-0.3_R8P + 2.6_R8P * (i - 1) / NC, i = 1, NC + 1)] ; y = x
call c%tag_block(x, y, 1_I4P, tags, status=status, feature_band=2_I4P, feature_tags=ft)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A,I0,A,I0,A,I0)', 'corners ', size(corners), ', corner cells ', count(ft == 0), ', within 2 cells ', count(ft <= 2)
call fig_open(u, 'ref-feature-corners-2d', 'feature_corners and feature_tags (band 2) of an L-shape, 32 x 32 cells', &
              -90._R8P, 90._R8P, height=520_I4P)
do j = 1, NC
   do i = 1, NC
      if (ft(i, j) > 2) cycle
      call fig_quad(u, 'frac', p(x(i), y(j)), p(x(i+1), y(j)), p(x(i+1), y(j+1)), p(x(i), y(j+1)), &
                    value=1._R8P - ft(i, j) / 3._R8P)
   enddo
enddo
do i = 1, 6
   a(i) = p(LU(i), LV(i), 1.e-3_R8P)
enddo
b = cshift(a, 1)
call fig_segments(u, 'accent', a, b)
allocate(pc(size(corners)))
do i = 1, size(corners)
   pc(i) = p(LU(corners(i)), LV(corners(i)), 2.e-3_R8P)
enddo
call fig_points(u, 'accent', pc)
call fig_close(u)
contains
   function p(pu, pv, dz) result(r)
   !< A point of the plane, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: pu, pv
   real(R8P), intent(in), optional :: dz
   type(vector_R8P)                :: r

   r = pu * ex_R8P + pv * ey_R8P
   if (present(dz)) r = r + dz * ez_R8P
   endfunction p
endprogram fig_ref_feature_corners_2d
