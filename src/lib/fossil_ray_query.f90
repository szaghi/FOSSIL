!< Ray-mesh intersection queries (issue #18 §2.5).
!<
!< Step-1 implementation: flat O(n) facet loop. The tree-accelerated traversal
!< lands in step 2; this module's public API stays the same when the body is
!< swapped, so callers don't need to change.

module fossil_ray_query
!< Ray-mesh intersection queries — driver that consumes a `facet_object` array
!< and returns hit lists for first / all / any queries.

use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private

public :: ray_hit_t
public :: ray_intersect_all_flat
public :: RAY_STATUS_OK, RAY_STATUS_BAD_INPUT

integer(I4P), parameter :: RAY_STATUS_OK        = 0_I4P
integer(I4P), parameter :: RAY_STATUS_BAD_INPUT = 1_I4P  !< Empty facet array, or |dir| < tolerance.

real(R8P), parameter :: DIR_TOL = 1.0e-30_R8P  !< |dir|^2 below this is rejected as degenerate.

type :: ray_hit_t
   !< One ray-facet intersection record.
   integer(I4P)     :: facet_id = 0_I4P     !< 1-based index into the facet array passed in.
   real(R8P)        :: t        = 0._R8P    !< Parametric distance along the ray.
   type(vector_R8P) :: point                !< Hit point in world coordinates.
endtype ray_hit_t

contains

   subroutine ray_intersect_all_flat(facet, ray_origin, ray_direction, hits, status)
   !< Return all ray-facet intersections via flat O(n) loop, sorted ascending by `t`.
   !<
   !< Hits with `t < 0` (behind the ray origin) are excluded.
   type(facet_object),           intent(in)               :: facet(:)        !< Facet array.
   type(vector_R8P),             intent(in)               :: ray_origin      !< Ray origin.
   type(vector_R8P),             intent(in)               :: ray_direction   !< Ray direction.
   type(ray_hit_t), allocatable, intent(out)              :: hits(:)         !< Sorted hit records.
   integer(I4P),                 intent(out),    optional :: status          !< Status code.
   real(R8P)                                              :: t, u, v         !< Möller-Trumbore outputs.
   logical                                                :: hit             !< Per-facet hit flag.
   real(R8P)                                              :: dir_norm2       !< |ray_direction|^2.
   integer(I4P)                                           :: nf, n_hits, f   !< Counters.
   type(ray_hit_t), allocatable                           :: tmp(:)          !< Pre-sort buffer.

   if (present(status)) status = RAY_STATUS_OK
   nf = size(facet, kind=I4P)
   if (nf == 0_I4P) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   dir_norm2 = ray_direction%dotproduct(rhs=ray_direction)
   if (dir_norm2 < DIR_TOL) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif

   allocate(tmp(nf))   ! upper bound: every facet hit
   n_hits = 0_I4P
   do f = 1_I4P, nf
      call facet(f)%intersect_ray(ray_origin=ray_origin, ray_direction=ray_direction, &
                                  t=t, u=u, v=v, hit=hit)
      if (.not. hit) cycle
      if (t < 0._R8P) cycle
      n_hits = n_hits + 1_I4P
      tmp(n_hits)%facet_id = f
      tmp(n_hits)%t        = t
      tmp(n_hits)%point    = ray_origin + t * ray_direction
   enddo

   allocate(hits(n_hits))
   hits(1:n_hits) = tmp(1:n_hits)
   call sort_hits_by_t(hits)
   endsubroutine ray_intersect_all_flat

   pure subroutine sort_hits_by_t(hits)
   !< Insertion sort on `t` ascending. n is small (typically 2..10), so the
   !< O(n^2) constant beats the O(n log n) call overhead.
   type(ray_hit_t), intent(inout) :: hits(:)
   type(ray_hit_t)                :: key
   integer(I4P)                   :: i, j

   do i = 2_I4P, size(hits, kind=I4P)
      key = hits(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (hits(j)%t <= key%t) exit
         hits(j + 1_I4P) = hits(j)
         j = j - 1_I4P
      enddo
      hits(j + 1_I4P) = key
   enddo
   endsubroutine sort_hits_by_t

endmodule fossil_ray_query
