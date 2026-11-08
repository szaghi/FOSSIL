!< Figure of the reference page "Self-intersection", section "Pipeline (detection)": six triangle pairs in one surface,
!< 2 apart along x, each pair a triangle T1 in the plane z = 0 (stone) and a second triangle T2 (teal): T2 pierces T1;
!< T2 touches T1's interior with one vertex; T2 overlaps T1 in the same plane; T2 shares an edge of T1 and folds back
!< over it in the plane; T2 shares an edge of T1 at an angle (an ordinary neighbour); T2 lies above T1. Under each pair,
!< the kind find_self_intersections(include_touching = .true.) reports for it, or "none"; the reported segment
!< p-q in red (the two centroids for a coplanar overlap).
program fig_ref_self_intersection_kinds
use fossil,     only : surface_stl_object, facet_object, intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, &
                       SELFINT_COPLANAR, SELFINT_EDGE_FOLD, SELFINT_SHEETS
use fossil_fig, only : fig_open, fig_triangle, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter                :: NP = 6
type(surface_stl_object)               :: soup
type(facet_object), allocatable        :: facets(:)
type(facet_object), pointer            :: f
type(intersection_pair_t), allocatable :: pairs(:)
type(vector_R8P)                       :: t2(3, NP), shift
character(16)                          :: kind(NP)
character(120)                         :: caption
integer(I4P)                           :: u, i, k, status

t2(:, 1) = [v(0.2_R8P, 0.2_R8P, -0.5_R8P), v(0.7_R8P, 0.2_R8P, 0.5_R8P), v(0.2_R8P, 0.7_R8P, 0.5_R8P)]   ! pierces
t2(:, 2) = [v(0.3_R8P, 0.3_R8P, 0._R8P), v(0.8_R8P, 0.3_R8P, 0.8_R8P), v(0.3_R8P, 0.8_R8P, 0.8_R8P)]    ! touches
t2(:, 3) = [v(0.3_R8P, 0.3_R8P, 0._R8P), v(1.2_R8P, 0.3_R8P, 0._R8P), v(0.3_R8P, 1.2_R8P, 0._R8P)]      ! overlaps
t2(:, 4) = [v(1._R8P, 0._R8P, 0._R8P), v(0._R8P, 0._R8P, 0._R8P), v(0.5_R8P, 0.6_R8P, 0._R8P)]          ! folds
t2(:, 5) = [v(1._R8P, 0._R8P, 0._R8P), v(0._R8P, 0._R8P, 0._R8P), v(0.5_R8P, -0.4_R8P, 0.6_R8P)]        ! neighbour
t2(:, 6) = [v(0._R8P, 0._R8P, 0.5_R8P), v(1._R8P, 0._R8P, 0.5_R8P), v(0._R8P, 1._R8P, 0.5_R8P)]         ! apart
allocate(facets(2 * NP))
do k = 1, NP
   shift = v(2._R8P * (k - 1), 0._R8P, 0._R8P)
   facets(2 * k - 1)%vertex = [v(0._R8P, 0._R8P, 0._R8P), v(1._R8P, 0._R8P, 0._R8P), v(0._R8P, 1._R8P, 0._R8P)]
   facets(2 * k)%vertex = t2(:, k)
   do i = 2 * k - 1, 2 * k
      facets(i)%vertex = facets(i)%vertex + shift
      call facets(i)%compute_normal
   enddo
enddo
call soup%adopt_facets(facets=facets)
call soup%find_self_intersections(pairs=pairs, status=status, include_touching=.true.)
kind = 'none'
call fig_open(u, 'ref-self-intersection-kinds', '', -90._R8P, 40._R8P, 760_I4P, 260_I4P)
do i = 1, size(pairs)
   k = (pairs(i)%a + 1) / 2
   select case (pairs(i)%kind)
   case (SELFINT_CROSSING)  ; kind(k) = 'CROSSING'
   case (SELFINT_TOUCHING)  ; kind(k) = 'TOUCHING'
   case (SELFINT_COPLANAR)  ; kind(k) = 'COPLANAR'
   case (SELFINT_EDGE_FOLD) ; kind(k) = 'EDGE_FOLD'
   case (SELFINT_SHEETS)    ; kind(k) = 'SHEETS'
   endselect
   call fig_segment(u, 'hole', pairs(i)%p, pairs(i)%q)
enddo
do k = 1, NP
   f => soup%facet_at(2 * k - 1)
   call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
   f => soup%facet_at(2 * k)
   call fig_triangle(u, 'probe', f%vertex(1), f%vertex(2), f%vertex(3))
   call fig_label(u, 'accent', v(2._R8P * (k - 1) - 0.1_R8P, -0.5_R8P, -0.5_R8P), trim(kind(k)))
   print '(A,I0,A,A)', 'pair ', k, ': ', trim(kind(k))
enddo
write(caption, '(A,I0,A)') 'find_self_intersections(include_touching=.true.): ', size(pairs), ' pairs of 6'
write(u, '(A)') 'caption '//trim(caption)
call fig_close(u)
contains
   pure function v(x, y, z) result(p)
   !< A point.
   real(R8P), intent(in) :: x, y, z !< Coordinates.
   type(vector_R8P)      :: p       !< Point.

   p = vector_R8P(x, y, z)
   endfunction v
endprogram fig_ref_self_intersection_kinds
