subroutine cylinder(s, c, r, z0, z1, m, name)
!< A regular m-gon of radius r about c, extruded over [z0, z1], caps fanned, outward, all in one patch.
type(surface_stl_object), intent(inout) :: s
real(R8P),                intent(in)    :: c(2), r, z0, z1
integer(I4P),             intent(in)    :: m
character(*),             intent(in)    :: name
type(facet_object), allocatable         :: f(:)
real(R8P)                               :: x(m), y(m)
integer(I4P)                            :: k, j, n

x = [(c(1) + r * cos(2 * PI * (k - 1) / m), k = 1, m)] ; y = [(c(2) + r * sin(2 * PI * (k - 1) / m), k = 1, m)]
allocate(f(2 * (m - 2) + 2 * m))
n = 0
do k = 2, m - 1
   n = n + 1 ; call tri(f(n), [x(1), y(1), z1], [x(k), y(k), z1], [x(k+1), y(k+1), z1])
   n = n + 1 ; call tri(f(n), [x(1), y(1), z0], [x(k+1), y(k+1), z0], [x(k), y(k), z0])
enddo
do k = 1, m
   j = mod(k, m) + 1
   n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z0], [x(j), y(j), z1])
   n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z1], [x(k), y(k), z1])
enddo
call s%destroy
call s%adopt_facets(facets=f)
call s%set_patch(facets=[(k, k = 1, s%get_facets_number())], name=name)
endsubroutine cylinder
