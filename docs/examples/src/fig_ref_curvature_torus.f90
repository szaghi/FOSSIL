!run figref-curvature-torus fig_ref_curvature_torus
!fig ref-curvature-torus
!fig ref-curvature-defect
!< Figures of the reference page "Curvature", What it does and Pipeline:
!< - a parametric torus (R = 1, r = 0.45, 48 x 24 vertices), each facet coloured by the mean of its vertices' Gaussian
!<   curvature from gaussian_curvature (teal K < 0 on the inner saddle side, ember K > 0 on the outer convex side,
!<   saturating at |K| = 1.5); the caption gives the range and the Gauss-Bonnet sum, 0 for a torus;
!< - the angle defect at the apex of a six-sided cone (an open fan): the angles at the apex (arcs), the apex's
!<   barycentric area (shaded) and its mean-curvature normal H n from mean_curvature (a segment along the axis).
program fig_ref_curvature_torus
use fossil,     only : surface_stl_object, facet_object, csr_matrix_t, CURV_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_triangle, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P), RR = 1._R8P, RT = 0.45_R8P
integer(I4P), parameter         :: NU = 48, NV = 24, NF = 6
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(csr_matrix_t)              :: lap, mass
type(vector_R8P)                :: apex, ring(NF), arcs(12, NF), cen, lift
real(R8P),          allocatable :: K(:), H(:)
real(R8P)                       :: theta, ksum, a_apex
integer(I4P)                    :: u, i, j, m, status, ia
character(96)                   :: caption
character(64)                   :: text

! ---- the torus
allocate(f(2 * NU * NV))
m = 0
do j = 0, NV - 1
   do i = 0, NU - 1
      call quad(torus_point(i, j), torus_point(i + 1, j), torus_point(i + 1, j + 1), torus_point(i, j + 1))
   enddo
enddo
call s%adopt_facets(facets=f) ; call s%analyze
call s%gaussian_curvature(K=K, status=status)
if (status /= CURV_STATUS_OK) error stop 'gaussian_curvature failed'
call s%cotangent_laplacian(L=lap, M=mass)                          ! the barycentric areas A_i, for Gauss-Bonnet
ksum = 0._R8P
do i = 1, size(K)
   ksum = ksum + K(i) * mass%get_value(i, i)
enddo
print '(A,I0,A,F7.3,A,F7.3)', 'torus: ', size(K), ' vertices, K from ', minval(K), ' to ', maxval(K)
print '(A,F7.3,A,F7.3)', 'exact: inner equator ', -1._R8P / (RT * (RR - RT)), ', outer equator ', 1._R8P / (RT * (RR + RT))
print '(A,F9.6)', '|sum K_i A_i| / (2 pi) = |Euler characteristic|: ', abs(ksum) / (2._R8P * PI)
write(caption, '(A,F6.2,A,F5.2,A,F5.2)') 'gaussian_curvature: K from ', minval(K), ' to ', maxval(K), &
                                         '; sum K A / 2 pi = ', abs(ksum) / (2._R8P * PI)
call fig_open(u, 'ref-curvature-torus', trim(caption), -60._R8P, 40._R8P)
do i = 1, s%get_facets_number()
   fp => s%facet_at(i)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), &
                     value=max(-1._R8P, min(1._R8P, sum(K(fp%vertex_id)) / 4.5_R8P)))
enddo
call fig_close(u)

! ---- the angle defect at the apex of an open six-sided cone
if (allocated(f)) deallocate(f)
allocate(f(NF))
apex = vector_R8P(0._R8P, 0._R8P, 0.6_R8P)
lift = vector_R8P(0._R8P, 0._R8P, 0.01_R8P)
do i = 1, NF
   ring(i) = vector_R8P(cos(2._R8P * PI * (i - 1) / NF), sin(2._R8P * PI * (i - 1) / NF), 0._R8P)
enddo
do i = 1, NF
   f(i)%vertex(1) = apex ; f(i)%vertex(2) = ring(i) ; f(i)%vertex(3) = ring(mod(i, NF) + 1)
   call f(i)%compute_metrix
enddo
call s%adopt_facets(facets=f) ; call s%analyze
call s%gaussian_curvature(K=K)
call s%mean_curvature(H=H)
call s%cotangent_laplacian(L=lap, M=mass)
fp => s%facet_at(1) ; ia = fp%vertex_id(1)
theta = NF * acos(dot_product(dir(apex, ring(1)), dir(apex, ring(2))))
a_apex = mass%get_value(ia, ia)
print '(A,F7.4,A,F7.4,A,F7.4)', 'apex: angle sum ', theta, ', defect 2 pi - sum ', 2._R8P * PI - theta, ', A_i ', a_apex
print '(A,F7.4,A,F7.4,A,F7.4)', 'apex: K = defect / A_i = ', (2._R8P * PI - theta) / a_apex, ', gaussian_curvature ', &
                                K(ia), ', mean_curvature H ', H(ia)
write(caption, '(A,F5.3,A,F5.3,A,F5.3,A,F5.3)') 'angle sum ', theta, ', defect ', 2._R8P * PI - theta, &
                                               ', K = defect / A = ', K(ia), ', H = ', H(ia)
call fig_open(u, 'ref-curvature-defect', trim(caption), -70._R8P, 50._R8P)
call fig_surface(u, 'body', s)
do i = 1, NF                                                       ! the barycentric area of the apex
   j = mod(i, NF) + 1
   cen = (apex + ring(i) + ring(j)) / 3._R8P + lift                ! drawn just above the facets
   call fig_triangle(u, 'patch', apex + lift, 0.5_R8P * (apex + ring(i)) + lift, cen)
   call fig_triangle(u, 'patch', apex + lift, cen, 0.5_R8P * (apex + ring(j)) + lift)
   do m = 1, 12                                                    ! the angle at the apex
      arcs(m, i) = apex + 4._R8P * lift + 0.3_R8P * slerp(dir(apex, ring(i)), dir(apex, ring(j)), real(m - 1, R8P) / 11._R8P)
   enddo
   call fig_segments(u, 'accent', arcs(1:11, i), arcs(2:12, i))
enddo
call fig_segment(u, 'probe', apex, apex + vector_R8P(0._R8P, 0._R8P, H(ia)))   ! H n: the pseudo-normal is +z here
call fig_points(u, 'foot', [apex])
write(text, '(A,F5.3)') 'H n, |H| = ', H(ia)
call fig_label(u, 'probe', apex + vector_R8P(0._R8P, 0._R8P, H(ia)), trim(text))
call fig_close(u)
contains
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

   function dir(a, b) result(d)
   !< Unit direction from a to b, as an array.
   type(vector_R8P), intent(in) :: a, b
   real(R8P)                    :: d(3)

   d = [b%x - a%x, b%y - a%y, b%z - a%z] ; d = d / norm2(d)
   endfunction dir

   function slerp(a, b, t) result(v)
   !< The unit vector at fraction t of the arc from a to b.
   real(R8P), intent(in) :: a(3), b(3), t
   type(vector_R8P)      :: v
   real(R8P)             :: w(3)

   w = (1._R8P - t) * a + t * b ; w = w / norm2(w)
   v = vector_R8P(w(1), w(2), w(3))
   endfunction slerp
endprogram fig_ref_curvature_torus
