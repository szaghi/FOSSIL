!run figref-alpha-wrap-stages fig_ref_alpha_wrap_stages
!fig ref-alpha-wrap-stages
!< Figure of the reference page "Alpha wrapping", section "Pipeline": the rules of stages 1-3 evaluated on the plane
!< y = 0.8333 through the hole of the page's input (the unit cube minus its first facet), alpha 0.15, offset 0.033,
!< cells of h = alpha/4. Stage 1: the unsigned distance u to the input at every cell centre (surface%distance). Stage 2:
!< free cells, u > alpha/2 + offset, flood-filled from the grid corner (4-neighbours in the plane): the reached centres
!< (teal). Stage 3: a cell within u(r) - offset of a reached centre r is carved (pale teal), unless it touches the
!< input (ember: the barrier, u below the half-diagonal of a cell). The rest is solid (cream). Over the cells, the
!< section of the input (gold) and of the wrap alpha_wrap returns (red): it follows the boundary of the carved region
!< at the offset, outside and inside the walls, entering through the hole.
program fig_ref_alpha_wrap_stages
use fossil,     only : surface_stl_object, facet_object, curve_object, slice_surface, AWRAP_STATUS_OK
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter                 :: ALPHA = 0.15_R8P, OFFSET = 0.033_R8P, LO = -0.3_R8P, HI = 1.3_R8P
type(surface_stl_object)                :: source, holed, wrapped
type(facet_object), pointer             :: facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
type(curve_object)                      :: section
real(R8P),    allocatable               :: d(:,:), xc(:)
integer(I4P), allocatable               :: state(:,:), queue(:,:)
real(R8P)                               :: h, cut, value
character(120)                          :: caption
integer(I4P)                            :: u, n, i, k, ri, rk, head, tail, status, nr

call source%load_from_file(file_name='cube-01.stl', guess_format=.true.)
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1))
do i = 1, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1)
enddo
cut = facets(1)%centroid%y
call holed%adopt_facets(facets=holed_facets)
h = ALPHA / 4._R8P
n = nint((HI - LO) / h)
allocate(xc(n), d(n, n), state(n, n), queue(2, n * n))
do i = 1, n
   xc(i) = LO + (i - 0.5_R8P) * h
enddo
do k = 1, n ; do i = 1, n                                          ! stage 1: u at the cell centres (x, cut, z)
   d(i, k) = holed%distance(point=vector_R8P(xc(i), cut, xc(k)), is_square_root=.true.)
enddo ; enddo
state = 0                                                          ! 0 solid, 1 reached, 2 carved, 3 barrier
where (d < h * sqrt(2._R8P) / 2._R8P) state = 3
head = 1 ; tail = 1 ; queue(:, 1) = [1, 1] ; state(1, 1) = 1     ! stage 2: flood fill of the free cells
do while (head <= tail)
   ri = queue(1, head) ; rk = queue(2, head) ; head = head + 1
   call visit(ri - 1, rk) ; call visit(ri + 1, rk) ; call visit(ri, rk - 1) ; call visit(ri, rk + 1)
enddo
nr = tail
do k = 1, n ; do i = 1, n                                          ! stage 3: carving by the reached balls
   if (state(i, k) /= 0) cycle
   do head = 1, nr
      ri = queue(1, head) ; rk = queue(2, head)
      if (hypot(xc(i) - xc(ri), xc(k) - xc(rk)) <= d(ri, rk) - OFFSET) then
         state(i, k) = 2 ; exit
      endif
   enddo
enddo ; enddo
call holed%alpha_wrap(alpha=ALPHA, offset=OFFSET, wrapped=wrapped, status=status)
if (status /= AWRAP_STATUS_OK) error stop 'alpha_wrap failed'
write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'y = 0.8333: ', count(state == 1), ' reached, ', &
      count(state == 2), ' carved, ', count(state == 3), ' barrier, ', count(state == 0), ' solid cells'
print '(A)', trim(caption)
call fig_open(u, 'ref-alpha-wrap-stages', trim(caption), -90._R8P, 0._R8P, 640_I4P, 520_I4P)
do k = 1, n ; do i = 1, n
   select case (state(i, k))
   case (1) ; value = -1._R8P
   case (2) ; value = -0.35_R8P
   case (3) ; value = 0.55_R8P
   case default ; value = 0._R8P
   endselect
   call fig_quad(u, 'field', p3(xc(i) - h / 2, xc(k) - h / 2), p3(xc(i) + h / 2, xc(k) - h / 2), &
                 p3(xc(i) + h / 2, xc(k) + h / 2), p3(xc(i) - h / 2, xc(k) + h / 2), value)
enddo ; enddo
call draw_section(holed, 'accent')
call draw_section(wrapped, 'hole')
call fig_label(u, 'accent', p3(0._R8P, 0.25_R8P), 'hole')
call fig_label(u, 'probe', p3(-0.25_R8P, 1.2_R8P), 'reached from the corner')
call fig_label(u, 'probe', p3(0.3_R8P, 0.6_R8P), 'reached through the hole')
call fig_label(u, 'accent', p3(0.5_R8P, 1.13_R8P), 'carved')
call fig_label(u, 'patch', p3(0.5_R8P, 1.03_R8P), 'barrier: wall cells')
call fig_label(u, 'hole', p3(0.5_R8P, -0.12_R8P), 'wrap section')
call fig_close(u)
contains
   subroutine visit(a, b)
   !< Queue cell (a, b) if it is free (u > alpha/2 + offset) and not reached yet.
   integer(I4P), intent(in) :: a, b !< Cell.

   if (a < 1 .or. b < 1 .or. a > n .or. b > n) return
   if (state(a, b) /= 0) return
   if (d(a, b) <= ALPHA / 2._R8P + OFFSET) return
   state(a, b) = 1 ; tail = tail + 1 ; queue(:, tail) = [a, b]
   endsubroutine visit

   function p3(x, z) result(p)
   !< The point (x, z) of the plane, slightly in front of it (towards the camera, -y).
   real(R8P), intent(in) :: x, z !< In-plane coordinates.
   type(vector_R8P)      :: p    !< Point.

   p = vector_R8P(x, cut + 0.5_R8P, z)
   endfunction p3

   subroutine draw_section(s, style)
   !< The section of a surface by the plane y = cut, drawn in front of the cells.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: style !< Style.
   type(vector_R8P), allocatable        :: a(:), b(:)
   real(R8P)                            :: p(2), q(2)
   integer(I4P)                         :: j, m

   call slice_surface(surface=s, normal_axis=2_I4P, level=cut, curve=section)
   allocate(a(0), b(0))
   do j = 1, section%get_vertices_number()
      m = section%get_next(j)
      if (m == 0) cycle
      p = section%get_vertex(j) ; q = section%get_vertex(m)       ! (u, v) = (z, x) on a y plane
      a = [a, vector_R8P(p(2), cut - 0.5_R8P, p(1))] ; b = [b, vector_R8P(q(2), cut - 0.5_R8P, q(1))]
   enddo
   call fig_segments(u, style, a, b)
   endsubroutine draw_section
endprogram fig_ref_alpha_wrap_stages
