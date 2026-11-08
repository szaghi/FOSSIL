!< Tutorial 7, second run: the placed dragon sliced by the planes z = -0.30 ... 0.30, step 0.05. Each frame of the
!< animation draws one plane (its square over the domain), the section curves on it (ember) and the dragon
!< translucent (a copy decimated to 600 facets); the captions give the loops and the area of each section.
program tut07_b_planes
use fossil,     only : surface_stl_object, curve_object, slice_surface
use fossil_fig, only : fig_open, fig_anim, fig_frame, fig_surface, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface, light
type(curve_object)       :: section
type(vector_R8P)         :: a, b
real(R8P)                :: z, uv(2)
integer(I4P)             :: k, s, open_loops, u, status
character(80)            :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
light = surface
call light%decimate(target_facets=600_I4P, status=status)
call fig_open(u, 'anim-tut07-planes', '', -60._R8P, 30._R8P)
call fig_anim(u, 0.7_R8P)
call fig_surface(u, 'ghost', light)
do k = -6, 6
   z = 0.05_R8P * k
   call slice_surface(surface, 3_I4P, z, section, merge_collinear=.true., open_loops=open_loops)
   write(caption, '(A,F5.2,A,I0,A,F7.4,A,I0,A)') 'z = ', z, ': loops ', section%get_loops_number(), ', area ', &
                                                section%get_area(), ', ', open_loops, ' open'
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call square(z)
   do s = 1, section%get_vertices_number()
      if (section%get_next(s) == 0_I4P) cycle
      uv = section%get_vertex(s) ; a = uv(1) * ex_R8P + uv(2) * ey_R8P + z * ez_R8P
      uv = section%get_vertex(section%get_next(s)) ; b = uv(1) * ex_R8P + uv(2) * ey_R8P + z * ez_R8P
      call fig_segment(u, 'hole', a, b)
   enddo
enddo
call fig_close(u)
contains
   subroutine square(z)
   !< The plane z over the placed dragon's box in x and y.
   real(R8P), intent(in) :: z !< Level.
   real(R8P)             :: x0, x1, y0, y1

   x0 = -0.3_R8P ; x1 = 0.3_R8P ; y0 = -0.55_R8P ; y1 = 0.55_R8P
   call fig_segment(u, 'box', x0 * ex_R8P + y0 * ey_R8P + z * ez_R8P, x1 * ex_R8P + y0 * ey_R8P + z * ez_R8P)
   call fig_segment(u, 'box', x1 * ex_R8P + y0 * ey_R8P + z * ez_R8P, x1 * ex_R8P + y1 * ey_R8P + z * ez_R8P)
   call fig_segment(u, 'box', x1 * ex_R8P + y1 * ey_R8P + z * ez_R8P, x0 * ex_R8P + y1 * ey_R8P + z * ez_R8P)
   call fig_segment(u, 'box', x0 * ex_R8P + y1 * ey_R8P + z * ez_R8P, x0 * ex_R8P + y0 * ey_R8P + z * ez_R8P)
   endsubroutine square
endprogram tut07_b_planes
