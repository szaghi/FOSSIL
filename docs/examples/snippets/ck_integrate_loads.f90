program ck_integrate_loads
use fossil, only : surface_stl_object, load_samples_t, load_samples_3d, integrate_loads
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
type(load_samples_t)     :: ls
real(R8P), allocatable   :: p(:), tau(:,:)
real(R8P)                :: force(3), moment(3), area

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube [0, 1]^3
call load_samples_3d(cube, samples=ls, order=2_I4P)        ! three points per facet
p = 1._R8P - ls%points(3, :)                               ! the solver's pressure: a fluid at rest, p = 1 - z
call integrate_loads(ls, force, moment, pressure=p, area=area)           ! moments about the origin
print '(A,3F9.5,A,3F9.5,A,F8.4)', 'pressure: force', r5(force), '  moment', r5(moment), '  area', area
allocate(tau(3, ls%n)) ; tau(1, :) = 0.1_R8P ; tau(2:3, :) = 0._R8P      ! a uniform shear along x
call integrate_loads(ls, force, moment, pressure=p, shear=tau, about=[0.5_R8P, 0.5_R8P, 0.5_R8P])
print '(A,3F9.5,A,3F9.5)', 'with shear: force', r5(force), '  moment about the centre', r5(moment)

contains
   pure function r5(x) result(y)
   !< Rounded to 5 decimals, no negative zero.
   real(R8P), intent(in) :: x(:)       !< Values.
   real(R8P)             :: y(size(x)) !< Rounded.

   y = anint(x * 1.e5_R8P) / 1.e5_R8P + 0._R8P
   endfunction r5
endprogram ck_integrate_loads
