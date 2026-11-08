!run figref-cotangent-laplacian-stencil fig_ref_cotangent_laplacian_stencil
!fig ref-laplacian-onering
!fig ref-laplacian-cotangents
!< Figures of the reference page "Cotangent Laplacian", What it does and Pipeline:
!< - an irregular one-ring, a centre vertex lifted above six neighbours, seen from above: each spoke labelled with its
!<   entry L(c, j) of the centre's row, the centre's barycentric area M(c, c) (a third of every incident triangle: the
!<   polygon through the edge midpoints and the triangle centroids) shaded;
!< - two flat triangles sharing the edge (i, j), the opposite angle alpha acute and beta obtuse, their cotangents and
!<   the entry L(i, j) = (cot alpha + cot beta) / 2 read back from the operator.
program fig_ref_cotangent_laplacian_stencil
use fossil,     only : surface_stl_object, facet_object, csr_matrix_t, LAPL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_segments, fig_points, fig_label, fig_triangle, &
                       fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P)
integer(I4P), parameter         :: NR = 6
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(csr_matrix_t)              :: lap, mass
type(vector_R8P)                :: c, ring(NR), mid, cen, arc_a(16), arc_b(16)
real(R8P)                       :: rad(NR), ang(NR), lsum, w, al, be
integer(I4P)                    :: u, i, j, k, status, ic, ids(NR), ii, jj
character(96)                   :: caption
character(16)                   :: text

! ---- the one-ring: centre lifted by 0.3, six neighbours at uneven radii and angles
rad = [1.0_R8P, 0.7_R8P, 1.1_R8P, 0.8_R8P, 1.2_R8P, 0.9_R8P]
ang = [0._R8P, 70._R8P, 125._R8P, 180._R8P, 230._R8P, 300._R8P] * PI / 180._R8P
c = vector_R8P(0._R8P, 0._R8P, 0.3_R8P)
do i = 1, NR
   ring(i) = vector_R8P(rad(i) * cos(ang(i)), rad(i) * sin(ang(i)), 0._R8P)
enddo
allocate(f(NR))
do i = 1, NR
   call tri(f(i), c, ring(i), ring(mod(i, NR) + 1))
enddo
call s%adopt_facets(facets=f) ; call s%analyze
call s%cotangent_laplacian(L=lap, M=mass, status=status)
if (status /= LAPL_STATUS_OK) error stop 'cotangent_laplacian failed'
fp => s%facet_at(1) ; ic = fp%vertex_id(1)                          ! the centre's pool id, and its neighbours'
do i = 1, NR
   fp => s%facet_at(i) ; ids(i) = fp%vertex_id(2)
enddo
lsum = 0._R8P
do i = 1, NR
   lsum = lsum + lap%get_value(ic, ids(i))
enddo
print '(A,6F7.3)', 'row of the centre, L(c, j): ', [(lap%get_value(ic, ids(i)), i = 1, NR)]
print '(A,F7.3,A,F7.3,A,F7.4)', 'L(c, c) = ', lap%get_value(ic, ic), ', -sum L(c, j) = ', -lsum, ', M(c, c) = ', &
                                mass%get_value(ic, ic)
write(caption, '(A,F6.3,A,F6.4,A)') 'row c of L on the spokes, L(c, c) = ', lap%get_value(ic, ic), &
                                    '; M(c, c) = ', mass%get_value(ic, ic), ' (shaded)'
call fig_open(u, 'ref-laplacian-onering', trim(caption), -90._R8P, 70._R8P)
call fig_surface(u, 'body', s)
do i = 1, NR                                                          ! the barycentric area: a third of each triangle
   j = mod(i, NR) + 1
   cen = (c + ring(i) + ring(j)) / 3._R8P
   call fig_triangle(u, 'patch', c, 0.5_R8P * (c + ring(i)), cen)
   call fig_triangle(u, 'patch', c, cen, 0.5_R8P * (c + ring(j)))
enddo
do i = 1, NR
   call fig_segment(u, 'accent', c, ring(i))
   mid = 0.55_R8P * c + 0.45_R8P * ring(i)
   write(text, '(F6.3)') lap%get_value(ic, ids(i))
   call fig_label(u, 'accent', mid, trim(adjustl(text)))
enddo
call fig_points(u, 'foot', [c])
call fig_label(u, 'foot', c, 'c')
call fig_close(u)

! ---- two triangles on the edge (i, j): alpha at k acute, beta at l obtuse
if (allocated(f)) deallocate(f)
allocate(f(2))
ring(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)                          ! i
ring(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)                          ! j
ring(3) = vector_R8P(0.8_R8P, 1.6_R8P, 0._R8P)                        ! k
ring(4) = vector_R8P(1.25_R8P, -0.45_R8P, 0._R8P)                     ! l
call tri(f(1), ring(1), ring(2), ring(3))
call tri(f(2), ring(2), ring(1), ring(4))
call s%adopt_facets(facets=f) ; call s%analyze
call s%cotangent_laplacian(L=lap, M=mass, status=status)
fp => s%facet_at(1) ; ii = fp%vertex_id(1) ; jj = fp%vertex_id(2)
al = angle_at(ring(3), ring(1), ring(2)) ; be = angle_at(ring(4), ring(1), ring(2))
print '(A,F7.2,A,F7.4)', 'alpha ', al * 180._R8P / PI, ' deg, cot ', 1._R8P / tan(al)
print '(A,F7.2,A,F7.4)', 'beta  ', be * 180._R8P / PI, ' deg, cot ', 1._R8P / tan(be)
w = lap%get_value(ii, jj)
print '(A,F7.4,A,F7.4)', 'L(i, j) = ', w, ', (cot alpha + cot beta) / 2 = ', 0.5_R8P * (1._R8P / tan(al) + 1._R8P / tan(be))
write(caption, '(A,F7.4)') 'L(i, j) = (cot alpha + cot beta) / 2 = ', w
call fig_open(u, 'ref-laplacian-cotangents', trim(caption), -90._R8P, 90._R8P)
call fig_surface(u, 'body', s)
call fig_segment(u, 'accent', ring(1), ring(2))
write(text, '(F7.4)') w
call fig_label(u, 'accent', 0.5_R8P * (ring(1) + ring(2)), 'L(i, j) = '//trim(adjustl(text)))
call arc(ring(3), ring(1), ring(2), arc_a) ; call arc(ring(4), ring(1), ring(2), arc_b)
call fig_segments(u, 'patch', arc_a(1:15), arc_a(2:16))
call fig_segments(u, 'probe', arc_b(1:15), arc_b(2:16))
call fig_points(u, 'foot', ring(1:4))
call fig_label(u, 'foot', ring(1), 'i')
call fig_label(u, 'foot', ring(2), 'j')
write(text, '(F6.3)') 1._R8P / tan(al)
call fig_label(u, 'patch', ring(3), 'k: alpha = '//trim(angle_text(al))//', cot = '//trim(adjustl(text)))
write(text, '(F6.3)') 1._R8P / tan(be)
call fig_label(u, 'probe', ring(4), 'l: beta = '//trim(angle_text(be))//', cot = '//trim(adjustl(text)))
call fig_close(u)
contains
   subroutine tri(ft, a, b, d)
   !< The facet a b d, as given (counter-clockwise seen from +z).
   type(facet_object), intent(inout) :: ft
   type(vector_R8P),   intent(in)    :: a, b, d

   ft%vertex(1) = a ; ft%vertex(2) = b ; ft%vertex(3) = d
   call ft%compute_metrix
   endsubroutine tri

   real(R8P) function angle_at(p, a, b)
   !< The angle at p of the triangle p a b.
   type(vector_R8P), intent(in) :: p, a, b
   real(R8P)                    :: e1(3), e2(3)

   e1 = [a%x - p%x, a%y - p%y, a%z - p%z] ; e2 = [b%x - p%x, b%y - p%y, b%z - p%z]
   angle_at = acos(dot_product(e1, e2) / (norm2(e1) * norm2(e2)))
   endfunction angle_at

   subroutine arc(p, a, b, pts)
   !< An arc of radius 0.25 about p from the direction of a to that of b, just above the plane z = 0.
   type(vector_R8P), intent(in)  :: p, a, b
   type(vector_R8P), intent(out) :: pts(:)
   real(R8P)                     :: t0, t1, t
   integer(I4P)                  :: m

   t0 = atan2(a%y - p%y, a%x - p%x) ; t1 = atan2(b%y - p%y, b%x - p%x)
   if (t1 - t0 > PI) t1 = t1 - 2._R8P * PI
   if (t0 - t1 > PI) t1 = t1 + 2._R8P * PI
   do m = 1, size(pts)
      t = t0 + (t1 - t0) * real(m - 1, R8P) / real(size(pts) - 1, R8P)
      pts(m) = vector_R8P(p%x + 0.25_R8P * cos(t), p%y + 0.25_R8P * sin(t), 0.01_R8P)   ! over the facets
   enddo
   endsubroutine arc

   function angle_text(a) result(txt)
   !< An angle in degrees, as text.
   real(R8P), intent(in) :: a
   character(12)         :: txt

   write(txt, '(F6.1,A)') a * 180._R8P / PI, ' deg'
   txt = adjustl(txt)
   endfunction angle_text
endprogram fig_ref_cotangent_laplacian_stencil
