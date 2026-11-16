program ck_load_samples
use fossil, only : surface_stl_object, curve_object, load_samples_t, load_samples_2d, load_samples_3d
use penf,   only : I4P, R8P
implicit none
type(curve_object)       :: body
type(surface_stl_object) :: cube
type(load_samples_t)     :: ls
real(R8P)                :: u(64), v(64), xc(20), f2(2), f3(3)
integer(I4P)             :: i, s

u = [(0.51_R8P + 0.35_R8P * cos(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]   ! an ellipse, counter-clockwise
v = [(0.52_R8P + 0.20_R8P * sin(2._R8P * acos(-1._R8P) * i / 64), i = 0, 63)]
call body%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xc = [((i - 0.5_R8P) * 0.05_R8P, i = 1, 20)]               ! the solver's cell centres
call load_samples_2d(body, xc, xc, ls, order=2_I4P, lengths=[0.05_R8P])
f2 = 0._R8P
do s = 1, ls%n                                             ! the solver's sum: here p = 1 - v, a fluid at rest
   f2 = f2 - (1._R8P - ls%points(2, s)) * ls%normal(:, s) * ls%weight(s)
enddo
print '(A,I0,A,I0,A,2F9.5,A,F9.5)', 'ellipse: samples ', ls%n, ', probe cells found ', count(ls%probe_cell(1, 1, :) > 0), &
      ', force', anint(f2 * 1.e5_R8P) / 1.e5_R8P + 0._R8P, '  area', body%get_area()
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube
call load_samples_3d(cube, samples=ls, order=1_I4P)        ! no lines: the samples only
f3 = 0._R8P
do s = 1, ls%n                                             ! p = 1 - z
   f3 = f3 - (1._R8P - ls%points(3, s)) * ls%normal(:, s) * ls%weight(s)
enddo
print '(A,I0,A,3F9.5)', 'cube: samples ', ls%n, ', force', anint(f3 * 1.e5_R8P) / 1.e5_R8P + 0._R8P
endprogram ck_load_samples
