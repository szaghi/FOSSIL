!run figref-2d-geometry fig_ref_2d_geometry
!fig ref-2d-geometry-winding
!fig ref-2d-geometry-signs
!< Figures of the reference page "2-D geometry (curves)": a section of the dragon with an open loop.
!< The facets of dragon.stl whose centroid lies within 0.5 of a point of the section z = 9.62 are removed, so the
!< slice of what remains has an open chain (its two ends in ember). First figure: the generalized winding number on a
!< grid of cells (teal 1, cream 0.5, ember 0), 0 or 1 away from the gap and fractional across it. Second figure: the
!< inside test with the fast winding number sign (teal inside, ember outside), and in gold the cells where the default
!< pseudo-normal sign disagrees with it.
program fig_ref_2d_geometry
use fossil,              only : surface_stl_object, curve_object, slice_surface, SIGN_FAST_WINDING_NUMBER
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter         :: Z0 = 9.62_R8P, R = 0.5_R8P      ! plane, radius of the cut
integer(I4P), parameter         :: NU = 150                        ! cells along u
type(surface_stl_object)        :: dragon, holed
type(curve_object)              :: section
type(facet_object), allocatable :: kept(:), sub(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: centre
type(vector_R8P), allocatable   :: a(:), b(:)
type(vector_R8P)                :: ends(2)
real(R8P)                       :: uv(2), c(3), lo(2), hi(2), h, zc(2)
integer(I4P)                    :: i, j, l, m, n, nv, n_open, first, last, disagree, fractional, ncells, u

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call slice_surface(dragon, 3_I4P, Z0, section, merge_collinear=.true.)
uv = section%get_vertex(1)                                          ! the cut: around the first vertex of the section
centre = uv(1) * ex_R8P + uv(2) * ey_R8P + Z0 * ez_R8P
allocate(kept(dragon%get_facets_number()))
m = 0
do i = 1, dragon%get_facets_number()
   f => dragon%facet_at(i)
   c = [f%vertex(1)%x + f%vertex(2)%x + f%vertex(3)%x, f%vertex(1)%y + f%vertex(2)%y + f%vertex(3)%y, &
        f%vertex(1)%z + f%vertex(2)%z + f%vertex(3)%z] / 3._R8P
   if (norm2(c - [centre%x, centre%y, centre%z]) > R) then
      m = m + 1 ; kept(m) = f
   endif
enddo
allocate(sub, source=kept(1:m))
call holed%adopt_facets(facets=sub)
call slice_surface(holed, 3_I4P, Z0, section, merge_collinear=.true., open_loops=n_open)
n = 0
do l = 1, section%get_loops_number()                                  ! the segments, closed loops or not
   call section%get_loop_range(l, first, last)
   if (.not. section%is_loop_closed(l)) then
      uv = section%get_vertex(first) ; ends(1) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.02_R8P * ez_R8P
      uv = section%get_vertex(last)  ; ends(2) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.02_R8P * ez_R8P
   endif
   n = n + last - first + merge(1, 0, section%is_loop_closed(l))
enddo
allocate(a(n), b(n))
nv = n ; n = 0
do l = 1, section%get_loops_number()
   call section%get_loop_range(l, first, last)
   do i = first, merge(last, last - 1, section%is_loop_closed(l))
      n = n + 1
      uv = section%get_vertex(i) ; a(n) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.01_R8P * ez_R8P
      uv = section%get_vertex(section%get_next(i)) ; b(n) = uv(1) * ex_R8P + uv(2) * ey_R8P + 0.01_R8P * ez_R8P
   enddo
enddo
ncells = 0
print '(A,I0,A,I0,A,I0,A,F6.3)', 'section of the cut dragon: ', section%get_loops_number(), ' loops, ', n_open, &
                                 ' open, ', nv, ' segments, gap ', norm2([ends(1)%x - ends(2)%x, ends(1)%y - ends(2)%y])

lo = section%get_bmin() ; hi = section%get_bmax()
h = (hi(1) - lo(1)) / (NU - 30) ; lo = lo - 15 * h ; hi = hi + 15 * h
fractional = 0 ; disagree = 0
call fig_open(u, 'ref-2d-geometry-winding', 'winding_number: teal 1, cream 0.5, ember 0; the open ends in ember', &
              -90._R8P, 90._R8P)
call cells(1_I4P, lo, hi, h)
call fig_segments(u, 'accent', a, b) ; call fig_points(u, 'patch', ends)
call fig_close(u)
call fig_open(u, 'ref-2d-geometry-signs', 'near the gap: fast winding number inside (teal); gold: pseudo-normal differs', &
              -90._R8P, 90._R8P)
zc = 0.5_R8P * [ends(1)%x + ends(2)%x, ends(1)%y + ends(2)%y] + [0.5_R8P, 0._R8P]
call fig_zoom(u, zc(1) * ex_R8P + zc(2) * ey_R8P, 1.6_R8P)
call cells(2_I4P, zc - [4._R8P, 2._R8P], zc + [4._R8P, 2._R8P], 0.05_R8P)
call fig_segments(u, 'accent', a, b) ; call fig_points(u, 'patch', ends)
call fig_close(u)
print '(A,I0,A,I0,A,I0,A)', 'cells: ', fractional, ' of ', ncells, ' with 0.05 < w < 0.95; near the gap ', disagree, &
                             ' where the pseudo-normal and the fast winding number signs disagree'
contains
   subroutine cells(kind, lo, hi, h)
   !< The cells of size h over [lo, hi]: kind 1 coloured by the winding number, kind 2 by the inside tests.
   integer(I4P), intent(in) :: kind
   real(R8P),    intent(in) :: lo(2), hi(2), h
   real(R8P)                :: x, y, w
   logical                  :: in_fwn, in_pn
   type(vector_R8P)         :: p1, p2, p3, p4

   do j = 0, nint((hi(2) - lo(2)) / h) - 1
      do i = 0, nint((hi(1) - lo(1)) / h) - 1
         x = lo(1) + (i + 0.5_R8P) * h ; y = lo(2) + (j + 0.5_R8P) * h
         p1 = (x - h / 2) * ex_R8P + (y - h / 2) * ey_R8P ; p2 = (x + h / 2) * ex_R8P + (y - h / 2) * ey_R8P
         p3 = (x + h / 2) * ex_R8P + (y + h / 2) * ey_R8P ; p4 = (x - h / 2) * ex_R8P + (y + h / 2) * ey_R8P
         if (kind == 1) then
            w = section%winding_number([x, y]) ; ncells = ncells + 1
            if (w > 0.05_R8P .and. w < 0.95_R8P) fractional = fractional + 1
            call fig_quad(u, 'field', p1, p2, p3, p4, value=max(-1._R8P, min(1._R8P, 1._R8P - 2._R8P * w)))
         else
            in_fwn = section%is_point_inside([x, y], sign_algorithm=SIGN_FAST_WINDING_NUMBER)
            in_pn = section%is_point_inside([x, y])
            if (in_fwn .neqv. in_pn) then
               disagree = disagree + 1
               call fig_quad(u, 'foot', p1, p2, p3, p4)
            else
               call fig_quad(u, 'field', p1, p2, p3, p4, value=merge(-0.8_R8P, 0.8_R8P, in_fwn))
            endif
         endif
      enddo
   enddo
   endsubroutine cells
endprogram fig_ref_2d_geometry
