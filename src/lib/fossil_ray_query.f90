!< Ray-mesh intersection queries (issue #18 §2.5): the hit record, the exact ray-facet kernel (issue #82) and the flat
!< O(n) reference driver.

module fossil_ray_query
!< Ray-mesh intersection queries — driver that consumes a `facet_object` array
!< and returns hit lists for first / all / any queries.
!<
!< Exact hits (issue #82). Every decision is exact, so a ray crosses a closed surface once per crossing: a ray through
!< a shared edge or vertex hits one of the facets there, a ray touching the surface at an edge or vertex without
!< crossing it, or lying in a face plane, hits none or two, and the parity of the hit count is the state of the
!< origin. Before #82 Moller-Trumbore with closed barycentric tests hit every facet at a shared edge (12 hits for 2
!< crossings through the face centres of cube-01) and could miss both at round-off.
!<
!< The ray o + t d is the oriented line through o and the far point q = o + s d (rounded once, `ray_prepare`), with s
!< such that q is beyond the box of the surface: every decision is then about two points in doubles.
!<
!< - Axis-aligned rays (two zero components of d) are decided by `tri_line_hit_oac`, the kernel of `line_crossings`
!<   and of the cut cells: the line translated by the global perturbation (eps1, eps2, eps3) of issue #69.
!< - Other rays: the line crosses triangle ABC iff the three Plucker side tests orient3d(o, q, A, B), (o, q, B, C),
!<   (o, q, C, A) have one sign (Woop et al. 2013, made exact). A zero side test (the line meets the edge line) takes
!<   one symbolic perturbation of the ray, o + delta and q + delta' with delta = (e1, e2, e4), delta' = (e8, e16, e32)
!<   (e^k: powers of one infinitesimal, so every product has its own order): orient3d(o + delta, q + delta', X, Y) =
!<   orient3d(o, q, X, Y) + g . delta + g' . delta' + (X - Y) . (delta x delta'), g = (q - Y) x (X - Y), g' = (X - Y) x
!<   (o - Y), whose sign is the sign of the first nonzero term in the order g_x, g_y, g_z, g'_x, -(X - Y)_z, (X - Y)_y,
!<   g'_y, -(X - Y)_x (orders e^1, 2, 4, 8, 10, 12, 16, 20). Each term is an exact orient2d or a comparison, the
!<   expansion is antisymmetric in (X, Y) as orient3d is, and X /= Y ends it: the perturbed line meets no edge, and
!<   the two facets of an edge see it on opposite sides.
!<
!< Both rays: the crossing is ahead iff o and q lie on opposite sides of the facet's plane (orient3d, exact; q is
!< beyond every facet). An origin exactly on the facet's plane is a hit at t = 0 (the crossing is at o): hits at t = 0
!< are kept, so a ray cast from the surface hits the facet it starts on, unless that facet is passed as `exclude`.
!< A facet with collinear vertices is never crossed. The topology is exact, t is rounded: the crossing point is the
!< barycentric combination by the side tests (by the line kernel's weights on axis rays), t its projection on d.
!<
!< The kernels themselves (`ray_prepare_oac`, `ray_crossing_oac`) live in `fossil_ray_kernel` on plain arrays, device
!< routines shared with the segment visibility of the probes (WP1 step 01).
!<
!< The tree traversals prune boxes with the slab test on boxes padded by `ray%pad` (1e-9 of the ray's extent over
!< the surface box, plus the octree's facet reach): a crossing on a box face is never lost to round-off.

use fossil_facet_object, only : facet_object
use fossil_ray_kernel,   only : ray_prepare_oac, ray_crossing_oac
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private

public :: ray_hit_t
public :: ray_exact_t, ray_prepare, ray_facet_crossing, sort_hits_by_t
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

type :: ray_exact_t
   !< A ray prepared for the exact kernel (issue #82): origin, far point beyond the surface box, direction.
   real(R8P)    :: o(3) = 0._R8P   !< Origin.
   real(R8P)    :: q(3) = 0._R8P   !< Far point o + s d, rounded, beyond the surface box.
   real(R8P)    :: d(3) = 0._R8P   !< Direction (as given: t is parametric).
   integer(I4P) :: axis = 0_I4P    !< 1, 2, 3 for an axis-aligned ray (two zero components of d), else 0.
   real(R8P)    :: pad  = 0._R8P   !< Pad of the boxes in the slab test.
endtype ray_exact_t

contains

   pure function ray_prepare(ray_origin, ray_direction, bmin, bmax, reach) result(ray)
   !< Prepare a ray for the exact kernel (issue #82): the far point q = o + s d beyond the box [bmin, bmax] (the surface
   !< box), s = 2 R / |d| with R >= the distance from o to any point of the box (s doubled until q /= o), and the pad
   !< of the slab test: 1e-9 R plus `reach` (the octree's facet reach; 0 for the BVH). Runs `ray_prepare_oac`, the
   !< routine of the device segment visibility (WP1 step 01).
   type(vector_R8P), intent(in)           :: ray_origin    !< Ray origin.
   type(vector_R8P), intent(in)           :: ray_direction !< Ray direction (nonzero).
   type(vector_R8P), intent(in)           :: bmin, bmax    !< Box of the surface.
   real(R8P),        intent(in), optional :: reach         !< Facet reach beyond its tree box (octree).
   type(ray_exact_t)                      :: ray           !< Prepared ray.
   real(R8P)                              :: lo(3), hi(3), r_ !< Box, reach.

   ray%o = [ray_origin%x, ray_origin%y, ray_origin%z]
   ray%d = [ray_direction%x, ray_direction%y, ray_direction%z]
   lo = [bmin%x, bmin%y, bmin%z] ; hi = [bmax%x, bmax%y, bmax%z]
   r_ = 0._R8P ; if (present(reach)) r_ = reach
   call ray_prepare_oac(ray%o, ray%d, lo, hi, r_, ray%q, ray%axis, ray%pad)
   endfunction ray_prepare

   subroutine ray_facet_crossing(ray, facet, hit, t)
   !< Exact crossing of the prepared ray with a facet (issue #82, see the module header): hit and the parameter t >= 0.
   !< Runs `ray_crossing_oac` on the facet's vertices, the routine of the device segment visibility (WP1 step 01).
   type(ray_exact_t),  intent(in)  :: ray    !< Prepared ray.
   type(facet_object), intent(in)  :: facet  !< Facet.
   logical,            intent(out) :: hit    !< The ray crosses the facet (t >= 0).
   real(R8P),          intent(out) :: t      !< Parameter of the crossing (0 if the origin is on the facet's plane).
   real(R8P)                       :: tv(9)  !< Vertices.

   tv(1) = facet%vertex(1)%x ; tv(2) = facet%vertex(1)%y ; tv(3) = facet%vertex(1)%z
   tv(4) = facet%vertex(2)%x ; tv(5) = facet%vertex(2)%y ; tv(6) = facet%vertex(2)%z
   tv(7) = facet%vertex(3)%x ; tv(8) = facet%vertex(3)%y ; tv(9) = facet%vertex(3)%z
   call ray_crossing_oac(ray%o, ray%q, ray%d, ray%axis, tv, hit, t)
   endsubroutine ray_facet_crossing

   subroutine ray_intersect_all_flat(facet, ray_origin, ray_direction, hits, status, exclude)
   !< Return all ray-facet intersections via flat O(n) loop, sorted ascending by `t`: one hit per crossing, exactly
   !< (issue #82, the module header), the facet `exclude` skipped.
   !<
   !< Hits with `t < 0` (behind the ray origin) are excluded; hits at t = 0 are kept.
   type(facet_object),           intent(in)               :: facet(:)        !< Facet array.
   type(vector_R8P),             intent(in)               :: ray_origin      !< Ray origin.
   type(vector_R8P),             intent(in)               :: ray_direction   !< Ray direction.
   type(ray_hit_t), allocatable, intent(out)              :: hits(:)         !< Sorted hit records.
   integer(I4P),                 intent(out),    optional :: status          !< Status code.
   integer(I4P),                 intent(in),     optional :: exclude         !< Facet to skip (issue #82).
   real(R8P)                                              :: t               !< Crossing parameter.
   logical                                                :: hit             !< Per-facet hit flag.
   real(R8P)                                              :: dir_norm2       !< |ray_direction|^2.
   integer(I4P)                                           :: nf, n_hits, f   !< Counters.
   type(ray_hit_t), allocatable                           :: tmp(:)          !< Pre-sort buffer.
   type(ray_exact_t)                                      :: ray             !< Prepared ray.
   type(vector_R8P)                                       :: bmin, bmax      !< Box of the facets.

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

   call facets_box(facet, bmin, bmax)
   ray = ray_prepare(ray_origin, ray_direction, bmin, bmax)
   allocate(tmp(nf))   ! upper bound: every facet hit
   n_hits = 0_I4P
   do f = 1_I4P, nf
      if (present(exclude)) then
         if (f == exclude) cycle
      endif
      call ray_facet_crossing(ray, facet(f), hit, t)
      if (.not. hit) cycle
      n_hits = n_hits + 1_I4P
      tmp(n_hits)%facet_id = f
      tmp(n_hits)%t        = t
      tmp(n_hits)%point    = ray_origin + t * ray_direction
   enddo

   allocate(hits(n_hits))
   hits(1:n_hits) = tmp(1:n_hits)
   call sort_hits_by_t(hits)
   endsubroutine ray_intersect_all_flat

   pure subroutine facets_box(facet, bmin, bmax)
   !< Box of the vertices of the facets.
   type(facet_object), intent(in)  :: facet(:)   !< Facets.
   type(vector_R8P),   intent(out) :: bmin, bmax !< Box.
   integer(I4P)                    :: f, k       !< Counters.

   bmin = vector_R8P(huge(1._R8P), huge(1._R8P), huge(1._R8P)) ; bmax = -bmin
   do f = 1, size(facet, kind=I4P)
      do k = 1, 3
         bmin%x = min(bmin%x, facet(f)%vertex(k)%x) ; bmax%x = max(bmax%x, facet(f)%vertex(k)%x)
         bmin%y = min(bmin%y, facet(f)%vertex(k)%y) ; bmax%y = max(bmax%y, facet(f)%vertex(k)%y)
         bmin%z = min(bmin%z, facet(f)%vertex(k)%z) ; bmax%z = max(bmax%z, facet(f)%vertex(k)%z)
      enddo
   enddo
   endsubroutine facets_box

   pure subroutine sort_hits_by_t(hits)
   !< Insertion sort on `t` ascending, then on the facet id (issue #82). n is small (typically 2..10), so the
   !< O(n^2) constant beats the O(n log n) call overhead.
   type(ray_hit_t), intent(inout) :: hits(:)
   type(ray_hit_t)                :: key
   integer(I4P)                   :: i, j

   do i = 2_I4P, size(hits, kind=I4P)
      key = hits(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (hits(j)%t < key%t .or. (hits(j)%t == key%t .and. hits(j)%facet_id <= key%facet_id)) exit ! t, then id
         hits(j + 1_I4P) = hits(j)
         j = j - 1_I4P
      enddo
      hits(j + 1_I4P) = key
   enddo
   endsubroutine sort_hits_by_t

endmodule fossil_ray_query
