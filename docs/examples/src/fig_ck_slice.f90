!run figck-slice fig_ck_slice
!fig ck-slice
!< Figure of the cookbook recipe "Cut a 2-D section of an STL": the NACA 0012 wing cut by the plane z = 0.5 with
!< slice_surface, as in ck_slice. The wing is a ghost (its facets run root to tip, so clip cannot split it at the
!< plane), the plane is the rectangle of the wing's box at z = 0.5, and the section (one closed loop) is drawn in it.
program fig_ck_slice
use fossil,     only : surface_stl_object, curve_object, slice_surface
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter          :: zs = 0.5_R8P                     ! the plane of the recipe
type(surface_stl_object)      :: wing
type(curve_object)            :: section
type(vector_R8P)              :: bmin, bmax, m
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: uv(2)
integer(I4P)                  :: open_loops, i, n, u
character(80)                 :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)  ! the recipe's call
call slice_surface(wing, 3_I4P, zs, section, open_loops=open_loops)
n = section%get_vertices_number()
allocate(a(n), b(n))
do i = 1, n                                                       ! segment i: vertex i to its next
   uv = section%get_vertex(i) ; a(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + zs * ez_R8P
   uv = section%get_vertex(section%get_next(i)) ; b(i) = uv(1) * ex_R8P + uv(2) * ey_R8P + zs * ez_R8P
enddo

bmin = wing%get_bmin() ; bmax = wing%get_bmax()
m = 0.08_R8P * (bmax - bmin) ; m%z = 0._R8P                       ! the plane drawn a little wider than the wing
write(caption, '(A,I0,A,I0,A,F7.5)') 'slice_surface(z = 0.5): ', section%get_loops_number(), ' closed loop, ', n, &
                                     ' vertices, area ', section%get_area()
call fig_open(u, 'ck-slice', trim(caption), -60._R8P, 45._R8P, height=460_I4P)
call fig_surface(u, 'ghost', wing)
call fig_box(u, 'box', (bmin%x - m%x) * ex_R8P + (bmin%y - m%y) * ey_R8P + zs * ez_R8P, &
                       (bmax%x + m%x) * ex_R8P + (bmax%y + m%y) * ey_R8P + zs * ez_R8P)
call fig_segments(u, 'hole', a, b)
call fig_close(u)
endprogram fig_ck_slice
