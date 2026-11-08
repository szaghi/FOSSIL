!run figref-amr-tagging-features fig_ref_amr_tagging_features
!fig ref-amr-tagging-features
!< Figure of the reference page "AMR tagging", section "Refining along feature edges": a cylinder (48 sides) in a block
!< of 16^3 cells. `build_features` finds the two rims (the side turns by 7.5 degrees per facet, below the 30-degree
!< threshold); `tag_block(feature_band=1)` tags the cells they touch, and `box_feature_many` on the same cells agrees.
program fig_ref_amr_tagging_features
use fossil,              only : surface_stl_object, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_quad, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P)
integer(I4P), parameter         :: NS = 48, NC = 16
type(surface_stl_object)        :: cylinder
type(facet_object), allocatable :: f(:)
real(R8P)                       :: x(NC + 1), y(NC + 1), z(NC + 1)
real(R8P),    allocatable       :: lo(:,:), hi(:,:)
logical,      allocatable       :: touched(:)
integer(I4P)                    :: tags(NC, NC, NC), ftags(NC, NC, NC), i, j, k, m, status, u
character(120)                  :: caption

call make_cylinder(f)
call cylinder%adopt_facets(facets=f) ; call cylinder%analyze
if (.not. cylinder%is_watertight()) error stop 'cylinder not watertight'
call cylinder%build_features(status=status)
if (status /= STATUS_OK) error stop 'build_features failed'
x = [(-0.71_R8P + 1.42_R8P * i / NC, i = 0, NC)] ; y = x ; z = [(-0.21_R8P + 1.42_R8P * i / NC, i = 0, NC)]
call cylinder%tag_block(x=x, y=y, z=z, band=2, tags=tags, feature_band=1, feature_tags=ftags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A,I0,A,I0,A,I0)', 'cells ', NC**3, ': cut ', count(tags == 0), ', touching a feature edge ', count(ftags == 0)
print '(A,I0)', 'cells one cell from a feature edge ', count(ftags == 1)
allocate(lo(3, NC**3), hi(3, NC**3), touched(NC**3))
m = 0
do k = 1, NC ; do j = 1, NC ; do i = 1, NC
   m = m + 1 ; lo(:, m) = [x(i), y(j), z(k)] ; hi(:, m) = [x(i + 1), y(j + 1), z(k + 1)]
enddo ; enddo ; enddo
call cylinder%box_feature_many(lo=lo, hi=hi, touched=touched, status=status)
print '(A,I0,A,L1)', 'box_feature_many on the same cells: ', count(touched), ' touched, same cells as the tags ', &
      all(touched .eqv. reshape(ftags == 0, [NC**3]))

write(caption, '(A,I0,A,I0,A)') 'cylinder in a 16^3 block: ', count(tags == 0), ' cut cells, of which the ', &
      count(ftags == 0), ' on the two rims tagged (ember)'
call fig_open(u, 'ref-amr-tagging-features', trim(caption), -60._R8P, 28._R8P, width=820_I4P, height=560_I4P)
call fig_surface(u, 'body', cylinder)
do k = 1, NC ; do j = 1, NC ; do i = 1, NC
   if (ftags(i, j, k) == 0) call cube(x(i), x(i + 1), y(j), y(j + 1), z(k), z(k + 1))
enddo ; enddo ; enddo
call fig_close(u)
contains
   subroutine cube(x0, x1, y0, y1, z0, z1)
   !< The six faces of a cell, a little smaller than the cell so neighbours stay apart.
   real(R8P), intent(in) :: x0, x1, y0, y1, z0, z1
   real(R8P)             :: a(3), b(3), s

   s = 0.06_R8P * (x1 - x0)
   a = [x0 + s, y0 + s, z0 + s] ; b = [x1 - s, y1 - s, z1 - s]
   call fig_quad(u, 'patch', v(a(1), a(2), a(3)), v(b(1), a(2), a(3)), v(b(1), b(2), a(3)), v(a(1), b(2), a(3)))
   call fig_quad(u, 'patch', v(a(1), a(2), b(3)), v(b(1), a(2), b(3)), v(b(1), b(2), b(3)), v(a(1), b(2), b(3)))
   call fig_quad(u, 'patch', v(a(1), a(2), a(3)), v(b(1), a(2), a(3)), v(b(1), a(2), b(3)), v(a(1), a(2), b(3)))
   call fig_quad(u, 'patch', v(a(1), b(2), a(3)), v(b(1), b(2), a(3)), v(b(1), b(2), b(3)), v(a(1), b(2), b(3)))
   call fig_quad(u, 'patch', v(a(1), a(2), a(3)), v(a(1), b(2), a(3)), v(a(1), b(2), b(3)), v(a(1), a(2), b(3)))
   call fig_quad(u, 'patch', v(b(1), a(2), a(3)), v(b(1), b(2), a(3)), v(b(1), b(2), b(3)), v(b(1), a(2), b(3)))
   endsubroutine cube

   function v(px, py, pz) result(q)
   !< A point.
   real(R8P), intent(in) :: px, py, pz
   type(vector_R8P)      :: q

   q = vector_R8P(px, py, pz)
   endfunction v

   subroutine make_cylinder(ff)
   !< Radius 0.5, z in [0, 1], NS sides; caps fanned from their centres; facets oriented away from (0, 0, 0.5).
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: c0(3), c1(3), d0(3), d1(3)
   integer(I4P)                                 :: s, mm

   allocate(ff(4 * NS)) ; mm = 0
   do s = 0, NS - 1
      c0 = [0.5_R8P * cos(2 * PI * s / NS), 0.5_R8P * sin(2 * PI * s / NS), 0._R8P]
      c1 = [0.5_R8P * cos(2 * PI * mod(s + 1, NS) / NS), 0.5_R8P * sin(2 * PI * mod(s + 1, NS) / NS), 0._R8P]
      d0 = c0 + [0._R8P, 0._R8P, 1._R8P] ; d1 = c1 + [0._R8P, 0._R8P, 1._R8P]
      mm = mm + 1 ; call tri(ff(mm), c0, c1, d1)
      mm = mm + 1 ; call tri(ff(mm), c0, d1, d0)
      mm = mm + 1 ; call tri(ff(mm), [0._R8P, 0._R8P, 0._R8P], c0, c1)
      mm = mm + 1 ; call tri(ff(mm), [0._R8P, 0._R8P, 1._R8P], d0, d1)
   enddo
   endsubroutine make_cylinder

   subroutine tri(ft, p1, p2, p3)
   !< A facet oriented away from the inner point (0, 0, 0.5).
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), nn(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - [0._R8P, 0._R8P, 0.5_R8P]) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri
endprogram fig_ref_amr_tagging_features
