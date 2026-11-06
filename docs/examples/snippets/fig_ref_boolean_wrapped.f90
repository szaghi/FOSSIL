!< Figure of the reference page "Boolean operations" (issue #63): the Stanford bunny and the dragon, closed by
!< alpha_wrap and scaled to a unit diagonal, united, then cut in half by a second boolean (the union intersected with
!< a box holding the half y < 0), so the cut face is part of the result, seen from +y. The caption gives the
!< volume identity V(A u B) + V(A n B) = V(A) + V(B) relative to V(A) + V(B).
program fig_ref_boolean_wrapped
use fossil,     only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: a, b, u, i, half, box
character(120)           :: caption
real(R8P)                :: res
integer(I4P)             :: f, st

call wrapped('bunny.stl', a)
call wrapped('dragon.stl', b)
call b%rotate(axis=vector_R8P(0._R8P, 0._R8P, 1._R8P), angle=1.2_R8P) ; call b%translate(x=0.3_R8P) ; call b%analyze
u = a ; call u%boolean(other=b, op=BOOL_UNION, status=st) ; if (st /= BOOL_STATUS_OK) error stop 'union failed'
i = a ; call i%boolean(other=b, op=BOOL_INTERSECT, status=st) ; if (st /= BOOL_STATUS_OK) error stop 'intersect failed'
res = abs(u%get_volume() + i%get_volume() - a%get_volume() - b%get_volume()) / (a%get_volume() + b%get_volume())
! the cut: a box holding the half y < 0 of the union (cube-01 is the unit cube [0, 1]^3)
call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call box%resize(factor=vector_R8P(2._R8P, 1._R8P, 2._R8P))
call box%translate(x=-1._R8P, y=-1._R8P, z=-1._R8P) ; call box%analyze
half = u ; call half%boolean(other=box, op=BOOL_INTERSECT, status=st) ; if (st /= BOOL_STATUS_OK) error stop 'cut failed'
write(caption, '(A,I0,A,ES7.1)') 'bunny u dragon, cut by a second boolean: ', half%get_facets_number(), &
                                  ' facets, identity residual ', res
print '(A,L1)', trim(caption)//', cut half watertight ', half%is_watertight()
call fig_open(f, 'ref-boolean-wrapped', trim(caption), 75._R8P, 20._R8P)
call fig_surface(f, 'body', half)
call fig_close(f)

contains
   subroutine wrapped(file_name, s)
   !< The alpha wrap of a mesh (alpha 1/20, offset 1/200 of its diagonal), scaled to a unit diagonal, centred.
   character(*),             intent(in)    :: file_name !< STL file.
   type(surface_stl_object), intent(inout) :: s         !< Wrapped surface.
   type(surface_stl_object)                :: raw       !< Input mesh.
   type(vector_R8P)                        :: lo, hi    !< Box.
   real(R8P)                               :: d         !< Diagonal.
   integer(I4P)                            :: st        !< Status.

   call raw%load_from_file(file_name=file_name, guess_format=.true.)
   lo = raw%get_bmin() ; hi = raw%get_bmax() ; d = norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call raw%alpha_wrap(alpha=d / 20._R8P, offset=d / 200._R8P, wrapped=s, status=st)
   lo = s%get_bmin() ; hi = s%get_bmax()
   call s%translate(x=-0.5_R8P * (lo%x + hi%x), y=-0.5_R8P * (lo%y + hi%y), z=-0.5_R8P * (lo%z + hi%z))
   call s%resize(factor=vector_R8P(1._R8P / d, 1._R8P / d, 1._R8P / d))
   call s%rotate(axis=vector_R8P(1._R8P, 0._R8P, 0._R8P), angle=2._R8P * atan(1._R8P)) ! the scans are y-up: z-up
   call s%analyze
   endsubroutine wrapped
endprogram fig_ref_boolean_wrapped
