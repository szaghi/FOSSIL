!< Tutorial 9, step 4: two more refinement fields on the same 32^3 block: the cells near a feature edge
!< (build_features, feature_band, feature_tags) and the cells near a wall thinner than a cell (thickness, build_thin,
!< thin_factor, thin_band, thin_tags). The figure: the cut cells as voxels, ember on a feature edge, teal on a thin
!< wall only, cream neither.
program tut09_d_features
use fossil_fig, only : fig_open, fig_quad, fig_close
use vecfor,     only : vector_R8P
use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nc = 32                     ! cells per direction of the tagged block
type(surface_stl_object) :: surface
real(R8P),   allocatable :: thickness(:), gap(:)
real(R8P)                :: x(nc + 1)
integer(I4P)             :: tags(nc, nc, nc), ftags(nc, nc, nc), ttags(nc, nc, nc), i, j, k, status, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]
call surface%build_features(angle=60._R8P, status=status) ! creases above 60 degrees, boundary, non-manifold
if (status /= STATUS_OK) error stop 'build_features failed'
call surface%thickness(thickness=thickness, gap=gap, rays=9) ! per facet: shape diameter and gap, 9 rays
call surface%build_thin(thickness=thickness, gap=gap, status=status)
if (status /= STATUS_OK) error stop 'build_thin failed'
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags, feature_band=2, feature_tags=ftags, &
                       thin_factor=1._R8P, thin_band=2, thin_tags=ttags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A,I0,A,I0)', 'facets ', surface%get_facets_number(), ', thinner than a cell side (0.0375): ', &
                     count(min(thickness, gap) < x(2) - x(1))
print '(A)', 'cells      at 0      1      2  beyond'
print '(A,4I7)', 'cut     ', count(tags == 0), count(abs(tags) == 1), count(abs(tags) == 2), count(abs(tags) > 2)
print '(A,4I7)', 'feature ', (count(ftags == i), i = 0, 3)
print '(A,4I7)', 'thin    ', (count(ttags == i), i = 0, 3)
print '(A,I0,A,I0,A,I0)', 'cut cells: on a feature edge ', count(tags == 0 .and. ftags == 0), &
      ', on a thin wall only ', count(tags == 0 .and. ftags /= 0 .and. ttags == 0), &
      ', neither ', count(tags == 0 .and. ftags /= 0 .and. ttags /= 0)
call fig_open(u, 'tut09-features', '32^3 cut cells: feature edge (ember), thin wall only (teal), neither (cream)', &
              -60._R8P, 22._R8P)
do k = 1, nc ; do j = 1, nc ; do i = 1, nc
   if (tags(i, j, k) /= 0) cycle
   if (ftags(i, j, k) == 0) then
      call voxel(0.8_R8P, i, j, k)                     ! on a feature edge: ember
   elseif (ttags(i, j, k) == 0) then
      call voxel(-0.8_R8P, i, j, k)                    ! on a thin wall only: teal
   else
      call voxel(0._R8P, i, j, k)                      ! neither: cream
   endif
enddo ; enddo ; enddo
call fig_close(u)
contains
   subroutine voxel(value, ci, cj, ck)
   !< The faces of cut cell (ci, cj, ck) not shared with another cut cell, coloured by value.
   real(R8P),    intent(in) :: value      !< Value of the field style.
   integer(I4P), intent(in) :: ci, cj, ck !< Cell.
   integer(I4P)             :: d, s, c(3), nb(3), ay, az
   real(R8P)                :: p0(3), p1(3), p2(3), p3(3)

   c = [ci, cj, ck]
   do d = 1, 3
      do s = 0, 1
         nb = c ; nb(d) = nb(d) + 2 * s - 1
         if (all(nb >= 1 .and. nb <= nc)) then
            if (tags(nb(1), nb(2), nb(3)) == 0) cycle
         endif
         ay = mod(d, 3) + 1 ; az = mod(d + 1, 3) + 1
         p0 = x(c) ; p0(d) = x(c(d) + s)
         p1 = p0 ; p1(ay) = x(c(ay) + 1)
         p2 = p1 ; p2(az) = x(c(az) + 1)
         p3 = p0 ; p3(az) = x(c(az) + 1)
         call fig_quad(u, 'field', v(p0), v(p1), v(p2), v(p3), value=value)
      enddo
   enddo
   endsubroutine voxel

   pure function v(q) result(w)
   !< The point q.
   real(R8P), intent(in) :: q(3) !< Coordinates.
   type(vector_R8P)      :: w    !< Point.

   w%x = q(1) ; w%y = q(2) ; w%z = q(3)
   endfunction v
endprogram tut09_d_features
