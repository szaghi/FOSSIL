!< Reference page "Sign algorithms at sharp features": the benchmark's method on one body, a 30-degree wedge built here.
!< Points x + delta d with x on a face or on the apex edge, d in the inside or the outside sector (a sixth of the opening
!< away from every face), delta from 1e-1 to 1e-12: the side is known by construction. Three variants: clean, one
!< facet removed (a hole), one facet flipped. Errors of the four sign algorithms and of the naive closest-facet normal.
program ref_sign_benchmark
use fossil,              only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                SIGN_FAST_WINDING_NUMBER
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = 4._R8P * atan(1._R8P), HALF = 15._R8P * PI / 180._R8P  ! half opening
integer(I4P), parameter         :: ND = 12                                           ! decades 1e-1 .. 1e-12
integer(I4P), parameter         :: ALGO(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                              SIGN_FAST_WINDING_NUMBER]
character(9), parameter         :: VARIANT(3) = ['clean    ', 'hole     ', 'flipped  ']
type(surface_stl_object)        :: wedge
type(facet_object), allocatable :: f(:)
type(vector_R8P),   allocatable :: p(:)
logical,            allocatable :: truth(:)
integer(I4P)                    :: v, a, i, wrong(5), ne, per_decade(ND)

call make_points(p, truth, ne)
print '(A,I0,A,I0,A)', 'wedge of 30 degrees, ', size(p), ' points (', ne, ' at the apex edge)'
print '(A)', 'variant    pseudo-normal  ray  solid angle  winding  naive normal'
do v = 1, 3
   call make_wedge(f, v)
   call wedge%adopt_facets(facets=f)
   call wedge%analyze
   do a = 1, 4
      wrong(a) = 0
      do i = 1, size(p)
         if (wedge%is_point_inside(point=p(i), sign_algorithm=ALGO(a)) .neqv. truth(i)) wrong(a) = wrong(a) + 1
      enddo
   enddo
   wrong(5) = 0
   do i = 1, size(p)
      if (naive_inside(p(i)) .neqv. truth(i)) wrong(5) = wrong(5) + 1
   enddo
   print '(A,I13,I5,I13,I9,I14)', VARIANT(v), wrong
   if (v == 1) then                                ! the naive errors by decade, delta = 1e-1 first
      do i = 1, ND
         per_decade(i) = count([(naive_inside(p(a)) .neqv. truth(a), a = 17 * (i - 1) + 1, 17 * i)])
      enddo
   endif
enddo
print '(A,12I3)', 'naive errors on the clean wedge by decade (1e-1 to 1e-12):', per_decade
contains
   logical function naive_inside(q)
   !< The naive sign: the normal of the closest facet, not the pseudo-normal of the closest feature.
   type(vector_R8P), intent(in) :: q
   type(facet_object), pointer  :: fp
   type(vector_R8P)             :: c
   real(R8P)                    :: d
   integer(I4P)                 :: k

   call wedge%compute_distance(point=q, distance=d, facet_index=k, closest_point=c)
   fp => wedge%facet_at(k)
   naive_inside = fp%normal%x * (q%x - c%x) + fp%normal%y * (q%y - c%y) + fp%normal%z * (q%z - c%z) < 0._R8P
   endfunction naive_inside

   subroutine make_points(q, inside, n_edge)
   !< Face points (x on the upper face at distance 0.5 from the apex, y = 0.25, d = +-n) and apex-edge points
   !< (x = (0, 0.5, 0), 12 outside directions, 3 inside), each at the 12 decades.
   type(vector_R8P), allocatable, intent(out) :: q(:)
   logical,          allocatable, intent(out) :: inside(:)
   integer(I4P),                  intent(out) :: n_edge
   real(R8P)                                  :: phi, delta, cu, su
   integer(I4P)                               :: k, j, n

   allocate(q(17 * ND), inside(17 * ND))
   cu = cos(HALF) ; su = sin(HALF)
   n = 0
   do j = 1, ND
      delta = 10._R8P**(-j)
      n = n + 1 ; q(n) = vector_R8P(0.5_R8P * cu - delta * su, 0.25_R8P, 0.5_R8P * su + delta * cu) ; inside(n) = .false.
      n = n + 1 ; q(n) = vector_R8P(0.5_R8P * cu + delta * su, 0.25_R8P, 0.5_R8P * su - delta * cu) ; inside(n) = .true.
      do k = 0, 11                                    ! outside sector: 15 + 55 to 345 - 55 degrees
         phi = (70._R8P + 220._R8P * k / 11._R8P) * PI / 180._R8P
         n = n + 1 ; q(n) = vector_R8P(delta * cos(phi), 0.5_R8P, delta * sin(phi)) ; inside(n) = .false.
      enddo
      do k = -1, 1                                    ! inside sector: -15 + 5 to 15 - 5 degrees
         phi = 10._R8P * k * PI / 180._R8P
         n = n + 1 ; q(n) = vector_R8P(delta * cos(phi), 0.5_R8P, delta * sin(phi)) ; inside(n) = .true.
      enddo
   enddo
   n_edge = 15 * ND
   endsubroutine make_points

   subroutine make_wedge(ff, variant)
   !< The prism of the triangle A = (0, 0), B = (1, tan 15), C = (1, -tan 15) in (x, z), y in [0, 1]: 8 facets,
   !< outward. Variant 2 drops facet 1 (half of the upper face), variant 3 flips it.
   type(facet_object), allocatable, intent(out) :: ff(:)
   integer(I4P),                    intent(in)  :: variant
   real(R8P)                                    :: t, a0(3), b0(3), c0(3), a1(3), b1(3), c1(3)
   integer(I4P)                                 :: k

   t = tan(HALF)
   a0 = [0._R8P, 0._R8P, 0._R8P] ; b0 = [1._R8P, 0._R8P, t] ; c0 = [1._R8P, 0._R8P, -t]
   a1 = a0 ; a1(2) = 1._R8P ; b1 = b0 ; b1(2) = 1._R8P ; c1 = c0 ; c1(2) = 1._R8P
   allocate(ff(merge(7, 8, variant == 2)))
   k = 0
   if (variant /= 2) call tri(ff, k, a0, b0, a1, variant == 3)       ! upper face: the defect facet, on the apex edge
   call tri(ff, k, b0, b1, a1)
   call tri(ff, k, a0, c0, c1) ; call tri(ff, k, a0, c1, a1)          ! lower face
   call tri(ff, k, b0, b1, c1) ; call tri(ff, k, b0, c1, c0)          ! back face x = 1
   call tri(ff, k, a0, b0, c0) ; call tri(ff, k, a1, c1, b1)          ! caps y = 0, y = 1
   endsubroutine make_wedge

   subroutine tri(ff, k, p1, p2, p3, flip)
   !< The facet p1 p2 p3 at ff(k + 1), oriented away from the wedge's centroid (inward if `flip`).
   type(facet_object), intent(inout)        :: ff(:)
   integer(I4P),       intent(inout)        :: k
   real(R8P),          intent(in)           :: p1(3), p2(3), p3(3)
   logical,            intent(in), optional :: flip
   real(R8P)                                :: e1(3), e2(3), n(3), o

   e1 = p2 - p1 ; e2 = p3 - p1
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   o = dot_product(n, (p1 + p2 + p3) / 3._R8P - [2._R8P / 3._R8P, 0.5_R8P, 0._R8P])
   if (present(flip)) then
      if (flip) o = -o
   endif
   k = k + 1
   ff(k)%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (o >= 0._R8P) then
      ff(k)%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ff(k)%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ff(k)%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ff(k)%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ff(k)%compute_metrix
   endsubroutine tri
endprogram ref_sign_benchmark
