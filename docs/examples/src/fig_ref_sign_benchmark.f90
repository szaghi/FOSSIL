!run figref-sign-benchmark fig_ref_sign_benchmark
!fig ref-sign-benchmark-naive
!fig ref-sign-benchmark-hole
!< Figures of the reference page "Sign algorithms at sharp features", on the 30-degree wedge of its example:
!< - the section y = 0.5 at the apex edge, seen along y, 72 points on a circle of radius 0.05 about the edge, coloured
!<   by the side each sign puts them on (teal inside, ember outside): frame 1 the naive closest-facet normal, frame 2 the
!<   pseudo-normal;
!< - the wedge with the upper-face facet on the apex edge removed, seen along y: its section y = 0.5 (the removed part of
!<   the upper face in ember) and a grid of points in that plane, coloured the same way, one frame per sign algorithm.
program fig_ref_sign_benchmark
use fossil,              only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                SIGN_FAST_WINDING_NUMBER
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_zoom, fig_segments, fig_points, fig_anim, &
                                fig_frame, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),     parameter        :: PI = 4._R8P * atan(1._R8P), HALF = 15._R8P * PI / 180._R8P
integer(I4P),  parameter        :: NR = 72, GX = 40, GZ = 24
integer(I4P),  parameter        :: ALGO(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                              SIGN_FAST_WINDING_NUMBER]
character(24), parameter        :: NAME(4) = ['SIGN_PSEUDO_NORMAL      ', 'SIGN_RAY_INTERSECTIONS  ', &
                                              'SIGN_SOLID_ANGLE        ', 'SIGN_FAST_WINDING_NUMBER']
type(surface_stl_object)        :: wedge
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: ring(NR), grid(GX * GZ), sa(3), sb(3)
logical                         :: inside(GX * GZ), truth(GX * GZ), naive(NR), pn(NR), ring_truth(NR)
real(R8P)                       :: t, phi
character(90)                   :: caption
integer(I4P)                    :: u, i, j, n, a

t = tan(HALF)
call make_wedge(f, .false.)
call wedge%adopt_facets(facets=f) ; call wedge%analyze
do i = 1, NR                                                      ! the ring about the apex edge
   phi = 2._R8P * PI * (i - 0.5_R8P) / NR
   ring(i) = vector_R8P(0.05_R8P * cos(phi), 0.5_R8P, 0.05_R8P * sin(phi))
   ring_truth(i) = abs(ring(i)%z) < t * ring(i)%x
   naive(i) = naive_inside(ring(i))
   pn(i) = wedge%is_point_inside(point=ring(i), sign_algorithm=SIGN_PSEUDO_NORMAL)
enddo
sa(1) = vector_R8P(0._R8P, 0.5_R8P, 0._R8P) ; sb(1) = vector_R8P(1._R8P, 0.5_R8P, t)   ! the section's outline
sa(2) = vector_R8P(0._R8P, 0.5_R8P, 0._R8P) ; sb(2) = vector_R8P(1._R8P, 0.5_R8P, -t)
sa(3) = sb(1) ; sb(3) = sb(2)
call fig_open(u, 'ref-sign-benchmark-naive', 'section y = 0.5 at the apex edge', -90._R8P, 0._R8P)
call fig_zoom(u, vector_R8P(0._R8P, 0.5_R8P, 0._R8P), 0.08_R8P)
call fig_anim(u, 2._R8P)
call fig_segments(u, 'box', sa, sb)
write(caption, '(A,I0,A,I0,A)') 'naive closest-facet normal: ', count(naive .neqv. ring_truth), ' of ', NR, ' wrong'
call fig_frame(u, trim(caption))
call fig_points(u, 'inside', pack(ring, naive)) ; call fig_points(u, 'outside', pack(ring, .not. naive))
write(caption, '(A,I0,A,I0,A)') 'SIGN_PSEUDO_NORMAL: ', count(pn .neqv. ring_truth), ' of ', NR, ' wrong'
call fig_frame(u, trim(caption))
call fig_points(u, 'inside', pack(ring, pn)) ; call fig_points(u, 'outside', pack(ring, .not. pn))
call fig_close(u)

call make_wedge(f, .true.)                                        ! the holed wedge
call wedge%adopt_facets(facets=f) ; call wedge%analyze
n = 0
do j = 1, GZ ; do i = 1, GX
   n = n + 1
   grid(n) = vector_R8P(-0.1_R8P + 1.2_R8P * (i - 0.5_R8P) / GX, 0.5_R8P, -0.36_R8P + 0.72_R8P * (j - 0.5_R8P) / GZ)
   truth(n) = abs(grid(n)%z) < t * grid(n)%x .and. grid(n)%x < 1._R8P
enddo ; enddo
call fig_open(u, 'ref-sign-benchmark-hole', 'wedge with a hole', -90._R8P, 0._R8P)
call fig_anim(u, 1.6_R8P)
call fig_segments(u, 'box', sa(2:3), sb(2:3))                     ! the section: lower and back faces,
call fig_segments(u, 'box', [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P * t)], [sb(1)])   ! the upper face kept,
call fig_segments(u, 'hole', [sa(1)], [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P * t)])  ! its part removed
do a = 1, 4
   do n = 1, GX * GZ
      inside(n) = wedge%is_point_inside(point=grid(n), sign_algorithm=ALGO(a))
   enddo
   write(caption, '(A,A,I0,A,I0,A)') trim(NAME(a)), ': ', count(inside .neqv. truth), ' of ', GX * GZ, ' wrong'
   call fig_frame(u, trim(caption))
   call fig_points(u, 'inside', pack(grid, inside)) ; call fig_points(u, 'outside', pack(grid, .not. inside))
enddo
call fig_close(u)
contains
   logical function naive_inside(q)
   !< The naive sign: the normal of the closest facet.
   type(vector_R8P), intent(in) :: q
   type(facet_object), pointer  :: fp
   type(vector_R8P)             :: c
   real(R8P)                    :: d
   integer(I4P)                 :: k

   call wedge%compute_distance(point=q, distance=d, facet_index=k, closest_point=c)
   fp => wedge%facet_at(k)
   naive_inside = fp%normal%x * (q%x - c%x) + fp%normal%y * (q%y - c%y) + fp%normal%z * (q%z - c%z) < 0._R8P
   endfunction naive_inside

   subroutine make_wedge(ff, holed)
   !< The wedge of the example (8 facets, outward); `holed`: without the upper-face facet on the apex edge.
   type(facet_object), allocatable, intent(out) :: ff(:)
   logical,                         intent(in)  :: holed
   real(R8P)                                    :: a0(3), b0(3), c0(3), a1(3), b1(3), c1(3)
   integer(I4P)                                 :: k

   a0 = [0._R8P, 0._R8P, 0._R8P] ; b0 = [1._R8P, 0._R8P, t] ; c0 = [1._R8P, 0._R8P, -t]
   a1 = a0 ; a1(2) = 1._R8P ; b1 = b0 ; b1(2) = 1._R8P ; c1 = c0 ; c1(2) = 1._R8P
   allocate(ff(merge(7, 8, holed)))
   k = 0
   if (.not. holed) call tri(ff, k, a0, b0, a1)
   call tri(ff, k, b0, b1, a1)
   call tri(ff, k, a0, c0, c1) ; call tri(ff, k, a0, c1, a1)
   call tri(ff, k, b0, b1, c1) ; call tri(ff, k, b0, c1, c0)
   call tri(ff, k, a0, b0, c0) ; call tri(ff, k, a1, c1, b1)
   endsubroutine make_wedge

   subroutine tri(ff, k, p1, p2, p3)
   !< The facet p1 p2 p3 at ff(k + 1), oriented away from the wedge's centroid.
   type(facet_object), intent(inout) :: ff(:)
   integer(I4P),       intent(inout) :: k
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), n(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   k = k + 1
   ff(k)%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(n, (p1 + p2 + p3) / 3._R8P - [2._R8P / 3._R8P, 0.5_R8P, 0._R8P]) >= 0._R8P) then
      ff(k)%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ff(k)%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ff(k)%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ff(k)%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ff(k)%compute_metrix
   endsubroutine tri
endprogram fig_ref_sign_benchmark
