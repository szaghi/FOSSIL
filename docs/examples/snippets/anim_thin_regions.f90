!< Headline figure and animation of the reference page "Thickness and gaps": a NACA 0012 wing tapering from chord 1 at
!< the root to 0.5 at the tip (span 1.5, 60 points per side, 12 layers), the medial thickness of every facet.
!< - Figure: the wing coloured by its thickness (teal thin, ember 0.12 and above), the facets thinner than h = 1/32
!<   outlined in gold: the trailing-edge strip, wider towards the tip where the chord is shorter.
!< - Animation: h halving from 1/4 to 1/128, the facets thinner than h in ember: the under-resolved region shrinks to
!<   the trailing edge, whose thickness tends to zero.
program anim_thin_regions
use fossil,              only : surface_stl_object, STATUS_OK, THICKNESS_MEDIAL
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_triangle, fig_segments, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P), SPAN = 1.5_R8P
integer(I4P), parameter         :: NP = 60, NZ = 12
type(surface_stl_object)        :: wing
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(vector_R8P), allocatable   :: a(:), b(:)
real(R8P), allocatable          :: th(:), area(:)
integer(I4P), allocatable       :: region(:)
real(R8P)                       :: h, su(2 * NP), sv(2 * NP), z(0:NZ), c(0:NZ)
integer(I4P)                    :: u, i, l, status, m
character(96)                   :: caption

call make_wing(f)
call wing%adopt_facets(facets=f) ; call wing%analyze
call wing%thickness(thickness=th, method=THICKNESS_MEDIAL, status=status)
if (status /= STATUS_OK .or. .not. wing%is_watertight()) error stop 'wing: thickness failed'
print '(A,I0,A,F7.4,A,F7.4)', 'wing: ', wing%get_facets_number(), ' facets, medial thickness ', minval(th), ' to ', &
                              maxval(th)

! ---- figure
call wing%thin_regions(thickness=th, h=1._R8P / 32, factor=1._R8P, region=region, region_area=area)
write(caption, '(A,I0,A,F6.4)') 'medial thickness; under h = 1/32: ', count(region > 0), ' facets, area ', sum(area)
print '(A)', trim(caption)
allocate(a(3 * count(region > 0)), b(3 * count(region > 0)))
m = 0
do i = 1, wing%get_facets_number()
   if (region(i) == 0) cycle
   fp => wing%facet_at(i)
   a(m + 1) = fp%vertex(1) ; b(m + 1) = fp%vertex(2) ; a(m + 2) = fp%vertex(2) ; b(m + 2) = fp%vertex(3)
   a(m + 3) = fp%vertex(3) ; b(m + 3) = fp%vertex(1) ; m = m + 3
enddo
call fig_open(u, 'ref-thickness-wing', trim(caption), -60._R8P, 40._R8P, width=900_I4P, height=460_I4P)
do i = 1, wing%get_facets_number()
   fp => wing%facet_at(i)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=2._R8P * min(1._R8P, th(i) / 0.12_R8P) - 1._R8P)
enddo
call fig_segments(u, 'accent', a, b)
call fig_close(u)

! ---- animation
call fig_open(u, 'anim-thin-regions', '', -60._R8P, 40._R8P, width=900_I4P, height=460_I4P)
call fig_anim(u, 1.2_R8P)
h = 0.25_R8P
do l = 1, 6
   call wing%thin_regions(thickness=th, h=h, factor=1._R8P, region=region, region_area=area)
   write(caption, '(A,I0,A,I0,A,F6.4)') 'h = 1/', nint(1._R8P / h), ': ', count(region > 0), ' facets thinner than h, area ', &
                                        sum(area)
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption), seconds=merge(2.5_R8P, 1.2_R8P, l == 6))
   do i = 1, wing%get_facets_number()
      fp => wing%facet_at(i)
      if (region(i) > 0) then
         call fig_triangle(u, 'hole', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      else
         call fig_triangle(u, 'body', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      endif
   enddo
   h = 0.5_R8P * h
enddo
call fig_close(u)
contains
   subroutine make_wing(ff)
   !< The tapered wing, span along y and thickness along z: section k (y = SPAN k / NZ) is the unit NACA 0012 scaled by
   !< c = 1 - 0.5 y / SPAN about the leading edge; the end caps fanned from the chord midpoint; facets oriented away from an inner point (convex body).
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: xs(0:NP), ys(0:NP)
   integer(I4P)                                 :: ii, k, n, mm, i2

   do ii = 0, NP
      xs(ii) = 0.5_R8P * (1._R8P - cos(PI * ii / NP))
      ys(ii) = 0.6_R8P * (0.2969_R8P * sqrt(xs(ii)) - 0.1260_R8P * xs(ii) - 0.3516_R8P * xs(ii) ** 2 + &
                          0.2843_R8P * xs(ii) ** 3 - 0.1036_R8P * xs(ii) ** 4)
   enddo
   ys(NP) = 0._R8P
   n = 0
   do ii = NP, 1, -1
      n = n + 1 ; su(n) = xs(ii) ; sv(n) = ys(ii)
   enddo
   do ii = 0, NP - 1
      n = n + 1 ; su(n) = xs(ii) ; sv(n) = -ys(ii)
   enddo
   z = [(SPAN * k / NZ, k = 0, NZ)] ; c = 1._R8P - 0.5_R8P * z / SPAN
   allocate(ff(2 * n * NZ + 2 * n))
   mm = 0
   do k = 0, NZ - 1
      do ii = 1, n
         i2 = mod(ii, n) + 1
         mm = mm + 1 ; call tri(ff(mm), pt(ii, k), pt(i2, k), pt(i2, k + 1))
         mm = mm + 1 ; call tri(ff(mm), pt(ii, k), pt(i2, k + 1), pt(ii, k + 1))
      enddo
   enddo
   do ii = 1, n
      i2 = mod(ii, n) + 1
      mm = mm + 1 ; call tri(ff(mm), [0.5_R8P * c(0), z(0), 0._R8P], pt(ii, 0), pt(i2, 0))
      mm = mm + 1 ; call tri(ff(mm), [0.5_R8P * c(NZ), z(NZ), 0._R8P], pt(ii, NZ), pt(i2, NZ))
   enddo
   endsubroutine make_wing

   pure function pt(jj, kk) result(p)
   !< Vertex jj of section kk of the wing.
   integer(I4P), intent(in) :: jj, kk
   real(R8P)                :: p(3)

   p = [c(kk) * su(jj), z(kk), c(kk) * sv(jj)]
   endfunction pt

   subroutine tri(ft, p1, p2, p3)
   !< A facet oriented away from the inner point (0.35, SPAN / 2, 0).
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), nn(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - [0.35_R8P, 0.5_R8P * SPAN, 0._R8P]) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri
endprogram anim_thin_regions
