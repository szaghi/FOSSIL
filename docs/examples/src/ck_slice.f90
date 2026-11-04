!run ck_slice ck_slice
program ck_slice
use fossil, only : surface_stl_object, curve_object, slice_surface
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: wing
type(curve_object)       :: section
integer(I4P)             :: open_loops

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section, open_loops=open_loops) ! the plane z = 0.5
print '(A,I0,A,I0,A,I0)', 'loops ', section%get_loops_number(), ', vertices ', section%get_vertices_number(), &
      ', open ', open_loops
print '(A,F8.5,A,F8.5)', 'area ', section%get_area(), ', perimeter ', section%get_perimeter()
call section%save_into_file(file_name='section.dat')
endprogram ck_slice
