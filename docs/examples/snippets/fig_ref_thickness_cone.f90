!< Figures of the reference page "Thickness and gaps", Ray method and AMR criterion:
!< - the ray set of one facet: a plate 4 x 1 x 0.2, a facet of its top face near the end x = 0; the 31 rays of the
!<   default 60-degree cone (cone_directions about -n) cast by intersect_ray_first(exclude=) up to their first hits, the
!<   ray of the lower median, the facet's thickness, in ember; the oblique rays are longer than t (the cone's bias);
!< - the AMR criterion: the plate tapering from 0.05 to 0.45 along x, coloured by its ray thickness, with a block of
!<   cells of side h = 0.25; build_thin then box_thin_many(factor = 1) marks the cells that meet a facet thinner than h
!<   (gold wireframes): the thin half of the plate, x < 2.
program fig_ref_thickness_cone
use fossil,           only : surface_stl_object, facet_object, ray_hit_t, STATUS_OK
use fossil_thickness, only : cone_directions, lower_median, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS
use fossil_fig,       only : fig_open, fig_zoom, fig_surface, fig_segments, fig_segment, fig_points, fig_box, &
                             fig_triangle, fig_close
use penf,             only : I4P, R8P
use vecfor,           only : vector_R8P

implicit none
integer(I4P), parameter         :: NX = 24, NY = 6, MR = THICKNESS_DEFAULT_RAYS, CX = 16, CY = 4, CZ = 2
real(R8P),    parameter         :: HC = 0.25_R8P
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(ray_hit_t)                 :: hit
type(vector_R8P)                :: o, a(MR), b(MR), d
real(R8P),        allocatable   :: th(:)
real(R8P)                       :: cen(3), dirs(3, MR), len(MR), med, lo(3, CX * CY * CZ), hi(3, CX * CY * CZ)
logical                         :: counted(MR), has_hit, touched(CX * CY * CZ)
integer(I4P)                    :: m, i, j, k, n, u, pick, status, invalid
character(100)                  :: caption

! ---- the ray set of one facet of a flat plate
call plate(0.2_R8P, 0._R8P)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
pick = 2 * (NX * 2 + 0) + 1                                         ! top face, cell (0, 2): next to the end x = 0
fp => s%facet_at(pick)
o = fp%centroid
call cone_directions([-fp%normal%x, -fp%normal%y, -fp%normal%z], THICKNESS_DEFAULT_CONE, MR, dirs)
n = 0
do k = 1, MR
   d = vector_R8P(dirs(1, k), dirs(2, k), dirs(3, k))
   call s%intersect_ray_first(ray_origin=o, ray_direction=d, hit=hit, has_hit=has_hit, exclude=pick)
   fp => s%facet_at(hit%facet_id)
   counted(k) = has_hit .and. fp%normal%x * d%x + fp%normal%y * d%y + fp%normal%z * d%z > 0._R8P   ! leaves the solid
   a(k) = o ; b(k) = hit%point ; len(k) = hit%t
   if (counted(k)) n = n + 1
enddo
med = lower_median(pack(len, counted))
print '(A,I0,A,I0,A,F7.4,A,F7.4,A,F7.4)', 'facet ', pick, ': ', n, ' rays counted, lengths ', minval(len), ' to ', &
                                          maxval(len), ', lower median ', med
print '(A,F7.4,A,L1)', 'thickness(facet) ', th(pick), ', equal to the median: ', th(pick) == med
write(caption, '(I0,A,F5.3,A,F5.3,A,F5.3)') n, ' rays counted, lengths ', minval(len), ' to ', maxval(len), &
                                            '; lower median = thickness ', th(pick)
call fig_open(u, 'ref-thickness-cone', trim(caption), -100._R8P, 8._R8P)
call fig_zoom(u, o - vector_R8P(0._R8P, 0._R8P, 0.1_R8P), 0.2_R8P)
call fig_surface(u, 'ghost', s)
call fig_segments(u, 'accent', pack(a, len /= med), pack(b, len /= med))
k = findloc(len, med, dim=1)
call fig_segment(u, 'patch', a(k), b(k))
call fig_points(u, 'foot', b)
call fig_points(u, 'probe', [o])
call fig_close(u)

! ---- the AMR criterion on the tapered plate
call plate(0.05_R8P, 0.10_R8P)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, status=status)
call s%build_thin(thickness=th, status=status)
if (status /= STATUS_OK) error stop 'build_thin failed'
n = 0
do k = 1, CZ ; do j = 1, CY ; do i = 1, CX
   n = n + 1
   lo(:, n) = [(i - 1) * HC, (j - 1) * HC, -0.25_R8P + (k - 1) * HC] ; hi(:, n) = lo(:, n) + HC
enddo ; enddo ; enddo
call s%box_thin_many(lo=lo, hi=hi, factor=1._R8P, touched=touched, status=status)
print '(A,I0,A,I0,A,F7.4,A)', 'cells ', n, ', thin-touched ', count(touched), ' (h = ', HC, ', factor 1)'
print '(A,F7.4)', 'largest x of a touched cell ', maxval(pack(hi(1, :), touched))
write(caption, '(I0,A,I0,A)') count(touched), ' of ', n, ' cells (h = 0.25) meet a facet thinner than h; teal 0, ember 0.5'
call fig_open(u, 'ref-thickness-amr', trim(caption), -60._R8P, 30._R8P)
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=2._R8P * min(1._R8P, th(i) / 0.5_R8P) - 1._R8P)
enddo
do i = 1, n
   if (touched(i)) call fig_box(u, 'accent', vector_R8P(lo(1, i), lo(2, i), lo(3, i)), vector_R8P(hi(1, i), hi(2, i), hi(3, i)))
enddo
call fig_close(u)
contains
   subroutine plate(t0, t1)
   !< The plate [0, 4] x [0, 1], z in [-t/2, t/2] with t = t0 + t1 x; NX x NY cells on the top and bottom faces (top
   !< first, cell (i, j) at facets 2 (NX j + i) + 1, 2), the side walls split to match, facets oriented outward.
   real(R8P), intent(in) :: t0, t1
   real(R8P)             :: x(0:NX), y(0:NY), zb(0:NX), zt(0:NX)
   integer(I4P)          :: i, j

   x = [(4._R8P * i / NX, i = 0, NX)] ; y = [(1._R8P * j / NY, j = 0, NY)]
   zb = -0.5_R8P * (t0 + t1 * x) ; zt = 0.5_R8P * (t0 + t1 * x)
   cen = [2._R8P, 0.5_R8P, 0._R8P]
   if (allocated(f)) deallocate(f)
   allocate(f(4 * NX * NY + 4 * NX + 4 * NY))
   m = 0
   do j = 0, NY - 1
      do i = 0, NX - 1
         call quad([x(i), y(j), zt(i)], [x(i + 1), y(j), zt(i + 1)], [x(i + 1), y(j + 1), zt(i + 1)], [x(i), y(j + 1), zt(i)])
      enddo
   enddo
   do j = 0, NY - 1
      do i = 0, NX - 1
         call quad([x(i), y(j), zb(i)], [x(i + 1), y(j), zb(i + 1)], [x(i + 1), y(j + 1), zb(i + 1)], [x(i), y(j + 1), zb(i)])
      enddo
   enddo
   do i = 0, NX - 1
      call quad([x(i), y(0), zb(i)], [x(i + 1), y(0), zb(i + 1)], [x(i + 1), y(0), zt(i + 1)], [x(i), y(0), zt(i)])
      call quad([x(i), y(NY), zb(i)], [x(i + 1), y(NY), zb(i + 1)], [x(i + 1), y(NY), zt(i + 1)], [x(i), y(NY), zt(i)])
   enddo
   do j = 0, NY - 1
      call quad([x(0), y(j), zb(0)], [x(0), y(j + 1), zb(0)], [x(0), y(j + 1), zt(0)], [x(0), y(j), zt(0)])
      call quad([x(NX), y(j), zb(NX)], [x(NX), y(j + 1), zb(NX)], [x(NX), y(j + 1), zt(NX)], [x(NX), y(j), zt(NX)])
   enddo
   endsubroutine plate

   subroutine quad(p1, p2, p3, p4)
   !< Two facets of the quad p1 p2 p3 p4, appended to f.
   real(R8P), intent(in) :: p1(3), p2(3), p3(3), p4(3)

   call tri(p1, p2, p3) ; call tri(p1, p3, p4)
   endsubroutine quad

   subroutine tri(p1, p2, p3)
   !< A facet oriented away from the plate's centre `cen` (the plate is convex).
   real(R8P), intent(in) :: p1(3), p2(3), p3(3)
   real(R8P)             :: e1(3), e2(3), nn(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   m = m + 1
   f(m)%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - cen) >= 0._R8P) then
      f(m)%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; f(m)%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      f(m)%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; f(m)%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call f(m)%compute_metrix
   endsubroutine tri
endprogram fig_ref_thickness_cone
