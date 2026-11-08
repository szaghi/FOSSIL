!< Tutorial 16, step 3: box classification on the device (chapter 9's tagging) and the 2-D twins on the section of
!< chapter 7: box_classify_many_device for the surface and for the curve, distance_many_device for the curve, each
!< against its host twin. The figure: the 32 x 32 boxes the curve's device classification returns, coloured by code,
!< with the section on top. Built with gfortran the device kernels run on the host.
program tut16_c_boxes
use fossil,     only : surface_stl_object, curve_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter   :: nb = 16, nc = 32, np = 64     ! 3-D boxes, 2-D boxes, 2-D points per direction
type(surface_stl_object)  :: surface
type(curve_object)        :: section
real(R8P)                 :: lo(3, nb**3), hi(3, nb**3), lo2(2, nc**2), hi2(2, nc**2), pt(2, np**2), d(np**2, 2), h
integer(I4P)              :: codes(nb**3, 2), codes2(nc**2, 2), i, j, k, b

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%enter_device
h = 1.2_R8P / nb
b = 0
do k = 1, nb ; do j = 1, nb ; do i = 1, nb
   b = b + 1
   lo(:, b) = -0.6_R8P + [i - 1, j - 1, k - 1] * h
   hi(:, b) = lo(:, b) + h
enddo ; enddo ; enddo
call surface%box_classify_many_device(lo=lo, hi=hi, codes=codes(:, 1), radius=0.05_R8P)  ! winding-number sign
call surface%box_classify_many(lo=lo, hi=hi, codes=codes(:, 2), radius=0.05_R8P)
print '(A,I0,A,I0)', 'surface boxes: device vs host mismatches ', count(codes(:, 1) /= codes(:, 2)), ' of ', nb**3
print '(A,I0,A,I0,A,I0,A,I0)', '  cut ', count(codes(:, 1) == BOX_CUT), ', near ', count(codes(:, 1) == BOX_NEAR), &
      ', inside ', count(codes(:, 1) == BOX_INSIDE), ', outside ', count(codes(:, 1) == BOX_OUTSIDE)
call surface%exit_device
call section%load_from_file(file_name='dragon-section.dat')         ! chapter 7: the section z = 0
call section%enter_device                                           ! the segment BVH and its caches
h = 1.2_R8P / nc
b = 0
do j = 1, nc ; do i = 1, nc
   b = b + 1
   lo2(:, b) = -0.6_R8P + [i - 1, j - 1] * h
   hi2(:, b) = lo2(:, b) + h
enddo ; enddo
call section%box_classify_many_device(lo=lo2, hi=hi2, codes=codes2(:, 1), radius=0.02_R8P)
call section%box_classify_many(lo=lo2, hi=hi2, codes=codes2(:, 2), radius=0.02_R8P)
print '(A,I0,A,I0)', 'curve boxes: device vs host mismatches ', count(codes2(:, 1) /= codes2(:, 2)), ' of ', nc**2
b = 0
do j = 1, np ; do i = 1, np
   b = b + 1
   pt(:, b) = -0.6_R8P + 1.2_R8P * [i - 1, j - 1] / (np - 1._R8P)
enddo ; enddo
call section%distance_many_device(points=pt, distances=d(:, 1))     ! signed, negative inside
call section%distance_many(points=pt, distances=d(:, 2))
print '(A,I0,A,I0,A,I0)', 'curve distances: device vs host mismatches ', count(d(:, 1) /= d(:, 2)), ' of ', np**2, &
      ', inside ', count(d(:, 1) < 0._R8P)
call section%exit_device
call draw_figure
contains
   subroutine draw_figure
   !< The 2-D boxes shaded by their device code: cut cream, near light, inside teal, outside ember; the section on top.
   use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_close
   type(vector_R8P), allocatable :: a(:), e(:)
   real(R8P)                     :: uv(2), v
   character(100)                :: caption
   integer(I4P)                  :: u, m, s

   write(caption, '(A,I0,A,I0,A,I0,A,I0)') 'curve box_classify_many_device: cut ', count(codes2(:, 1) == BOX_CUT), &
         ', near ', count(codes2(:, 1) == BOX_NEAR), ', inside ', count(codes2(:, 1) == BOX_INSIDE), &
         ', outside ', count(codes2(:, 1) == BOX_OUTSIDE)
   call fig_open(u, 'tut16-boxes', trim(caption), -90._R8P, 90._R8P)
   call fig_zoom(u, 0._R8P * ex_R8P, 0.5_R8P)                     ! the boxes about the section
   do b = 1, nc**2
      select case (codes2(b, 1))
      case (BOX_CUT)    ; v = 0._R8P
      case (BOX_NEAR)   ; v = 0.45_R8P
      case (BOX_INSIDE) ; v = -0.85_R8P
      case default      ; v = 0.85_R8P
      endselect
      call fig_quad(u, 'field', lo2(1, b) * ex_R8P + lo2(2, b) * ey_R8P, hi2(1, b) * ex_R8P + lo2(2, b) * ey_R8P, &
                    hi2(1, b) * ex_R8P + hi2(2, b) * ey_R8P, lo2(1, b) * ex_R8P + hi2(2, b) * ey_R8P, value=v)
   enddo
   allocate(a(section%get_vertices_number()), e(section%get_vertices_number())) ; m = 0
   do s = 1, section%get_vertices_number()
      if (section%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = section%get_vertex(s) ; a(m) = uv(1) * ex_R8P + uv(2) * ey_R8P + 1.e-3_R8P * ez_R8P
      uv = section%get_vertex(section%get_next(s)) ; e(m) = uv(1) * ex_R8P + uv(2) * ey_R8P + 1.e-3_R8P * ez_R8P
   enddo
   call fig_segments(u, 'accent', a(1:m), e(1:m))
   call fig_close(u)
   endsubroutine draw_figure
endprogram tut16_c_boxes
