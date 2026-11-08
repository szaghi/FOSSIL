use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: part, boss, bore
real(R8P)                :: a_boss, a_bore, v_exact
integer(I4P)             :: status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')      ! a 2 x 1 x 0.5 base
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 32_I4P, 'boss') ! a round boss on it
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 32_I4P, 'bore') ! a bore, 0.06 from x = 0
call report('base', part)
call report('boss', boss)
call report('bore', bore)
