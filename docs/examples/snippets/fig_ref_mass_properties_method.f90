!< Figure of the reference page "Mass properties", Method: a torus (R = 1, r = 0.4, 32 x 16 vertices) and the
!< reference point o, its bounding-box centre, which lies in the hole, outside the body. Every facet with o forms a
!< signed tetrahedron, D = a . (b x c) with a, b, c relative to o: the facets coloured ember where D > 0 (the outer
!< side, facing away from o) and teal where D < 0 (the inner side, facing o); the tetrahedra of one outer and one
!< inner facet drawn as wireframes. Only the half z < 0 is drawn, so the tetrahedra show through the cut. The negative
!< volumes cancel the part of the positive ones that lies in the hole:
!< the sum of D / 6 is the volume that mass_properties returns.
program fig_ref_mass_properties_method
use fossil,     only : surface_stl_object, facet_object, STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P), RR = 1._R8P, RT = 0.4_R8P
integer(I4P), parameter         :: NU = 32, NV = 16
type(surface_stl_object)        :: torus
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(vector_R8P)                :: o, bmin, bmax, ea(12), eb(12)
real(R8P)                       :: dd, vpos, vneg, vol
integer(I4P)                    :: i, j, m, u, status, show(2)
character(100)                  :: caption

allocate(f(2 * NU * NV))
m = 0
do j = 0, NV - 1
   do i = 0, NU - 1
      call quad(torus_point(i, j), torus_point(i + 1, j), torus_point(i + 1, j + 1), torus_point(i, j + 1))
   enddo
enddo
call torus%adopt_facets(facets=f) ; call torus%analyze
call torus%mass_properties(volume=vol, status=status)
if (status /= STATUS_OK) error stop 'mass_properties failed'
bmin = torus%get_bmin() ; bmax = torus%get_bmax()
o = 0.5_R8P * (bmin + bmax)                                         ! the reference point: the bounding-box centre
vpos = 0._R8P ; vneg = 0._R8P
do i = 1, torus%get_facets_number()
   dd = signed6(i)
   if (dd > 0._R8P) then
      vpos = vpos + dd / 6._R8P
   else
      vneg = vneg + dd / 6._R8P
   endif
enddo
print '(A,F8.5,A,F8.5,A,F8.5)', 'sum of D/6: positive ', vpos, ', negative ', vneg, ', total ', vpos + vneg
print '(A,F8.5,A,F8.5)', 'mass_properties volume ', vol, ', smooth torus 2 pi^2 R r^2 ', 2._R8P * PI**2 * RR * RT**2
write(caption, '(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') 'sum D/6: ember ', vpos, ', teal ', vneg, ', total ', vpos + vneg, &
                                                ' = volume ', vol
show = [2 * ((NV - 2) * NU + 9) + 1, 2 * ((NV / 2 + 1) * NU + 3) + 1] ! an outer facet, an inner one, both at z < 0
call fig_open(u, 'ref-mass-properties-tetrahedra', trim(caption), -80._R8P, 55._R8P)
do i = 1, torus%get_facets_number()                                 ! the half z < 0, cut open towards the viewer
   fp => torus%facet_at(i)
   if (fp%vertex(1)%z + fp%vertex(2)%z + fp%vertex(3)%z > 0._R8P) cycle
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=merge(0.7_R8P, -0.7_R8P, signed6(i) > 0))
enddo
do j = 1, 2                                                         ! the two tetrahedra: 6 edges each
   fp => torus%facet_at(show(j))
   do i = 1, 3
      ea(6 * (j - 1) + i) = o ; eb(6 * (j - 1) + i) = fp%vertex(i)
      ea(6 * (j - 1) + 3 + i) = fp%vertex(i) ; eb(6 * (j - 1) + 3 + i) = fp%vertex(mod(i, 3) + 1)
   enddo
enddo
call fig_segments(u, 'accent', ea, eb)
call fig_points(u, 'foot', [o])
call fig_label(u, 'foot', o, 'o')
call fig_close(u)
contains
   real(R8P) function signed6(k)
   !< D = a . (b x c) of facet k, its vertices relative to o: six times the signed volume of the tetrahedron.
   integer(I4P), intent(in)    :: k
   type(facet_object), pointer :: ft
   real(R8P)                   :: a(3), b(3), c(3)

   ft => torus%facet_at(k)
   a = [ft%vertex(1)%x - o%x, ft%vertex(1)%y - o%y, ft%vertex(1)%z - o%z]
   b = [ft%vertex(2)%x - o%x, ft%vertex(2)%y - o%y, ft%vertex(2)%z - o%z]
   c = [ft%vertex(3)%x - o%x, ft%vertex(3)%y - o%y, ft%vertex(3)%z - o%z]
   signed6 = a(1) * (b(2) * c(3) - b(3) * c(2)) + a(2) * (b(3) * c(1) - b(1) * c(3)) + a(3) * (b(1) * c(2) - b(2) * c(1))
   endfunction signed6

   function torus_point(i, j) result(p)
   !< Vertex (i, j) of the torus, indices taken modulo the grid so the seams close bit for bit.
   integer(I4P), intent(in) :: i, j
   type(vector_R8P)         :: p
   real(R8P)                :: a, b

   a = 2._R8P * PI * mod(i, NU) / NU ; b = 2._R8P * PI * mod(j, NV) / NV
   p = vector_R8P((RR + RT * cos(b)) * cos(a), (RR + RT * cos(b)) * sin(a), RT * sin(b))
   endfunction torus_point

   subroutine quad(p1, p2, p3, p4)
   !< Two facets of the quad p1 p2 p3 p4 (counter-clockwise seen from outside), appended to f.
   type(vector_R8P), intent(in) :: p1, p2, p3, p4

   m = m + 1 ; f(m)%vertex(1) = p1 ; f(m)%vertex(2) = p2 ; f(m)%vertex(3) = p3 ; call f(m)%compute_metrix
   m = m + 1 ; f(m)%vertex(1) = p1 ; f(m)%vertex(2) = p3 ; f(m)%vertex(3) = p4 ; call f(m)%compute_metrix
   endsubroutine quad
endprogram fig_ref_mass_properties_method
