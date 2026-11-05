!< FOSSIL, surface deviation: the one-sided Hausdorff distance between two surfaces (issue #60, roadmap #58 WP2).
module fossil_deviation
!< FOSSIL, surface deviation: the one-sided Hausdorff distance between two surfaces (issue #60, roadmap #58 WP2).
!<
!< h(A, B) = sup over x in A of d(x, B). A is given as a triangle mesh (vertices and facets), B only through a callback
!< that returns d(x, B) for a batch of points (the caller's BVH query), so this module depends on no surface type.
!<
!< Lower bound: the largest distance over deterministic samples of A: every vertex, every facet's edge midpoints, and the
!< centroids of the n_f^2 sub-triangles of facet f (n_f = ceil(sqrt(area_f * samples_per_area))). The same centroids,
!< weighted by area_f / n_f^2 (the centroid rule on the sub-triangles: exact for a linear function), give the mean and
!< the rms of d over A: estimates, not bounds.
use penf, only : I4P, R8P

implicit none
private
public :: deviation_eval_interface
public :: deviation_sampled
public :: DEV_STATUS_OK, DEV_STATUS_BAD_INPUT

integer(I4P), parameter :: DEV_STATUS_OK        = 0_I4P !< Success.
integer(I4P), parameter :: DEV_STATUS_BAD_INPUT = 1_I4P !< Empty surface, or a size mismatch.

abstract interface
   subroutine deviation_eval_interface(points, distances)
   !< Distances d(points(:, i), B) of a batch of points.
   import :: R8P
   real(R8P), intent(in)  :: points(:,:)  !< Points [3, n].
   real(R8P), intent(out) :: distances(:) !< Unsigned distances [n].
   endsubroutine deviation_eval_interface
endinterface

contains
   subroutine deviation_sampled(vertices, facets, eval, samples_per_area, lower, mean, rms, per_vertex, samples, status)
   !< Sampled one-sided deviation of the mesh (vertices, facets) from the surface behind `eval`.
   real(R8P),    intent(in)            :: vertices(:,:)     !< Vertices of A [3, nv].
   integer(I4P), intent(in)            :: facets(:,:)       !< Facets of A, vertex indices [3, nf].
   procedure(deviation_eval_interface) :: eval              !< Batched distance to B.
   real(R8P),    intent(in)            :: samples_per_area  !< Interior samples per unit area (> 0).
   real(R8P),    intent(out)           :: lower             !< Largest sampled distance (a lower bound of h).
   real(R8P),    intent(out)           :: mean              !< Area-weighted mean of d over A (estimate).
   real(R8P),    intent(out)           :: rms               !< Area-weighted rms of d over A (estimate).
   real(R8P),    intent(out)           :: per_vertex(:)     !< d at every vertex [nv].
   integer(I4P), intent(out)           :: samples           !< Number of points evaluated.
   integer(I4P), intent(out)           :: status            !< DEV_STATUS_OK or DEV_STATUS_BAD_INPUT.
   real(R8P),    allocatable           :: p(:,:)            !< Sample points.
   real(R8P),    allocatable           :: w(:)              !< Quadrature weights (0 for vertices and midpoints).
   real(R8P),    allocatable           :: d(:)              !< Distances of the samples.
   integer(I4P), allocatable           :: nsub(:)           !< Subdivisions per facet.
   real(R8P)                           :: a(3), b(3), c(3)  !< Facet corners.
   real(R8P)                           :: area              !< Facet area.
   real(R8P)                           :: wsum              !< Total weight.
   integer(I4P)                        :: nv, nf, np        !< Counters.
   integer(I4P)                        :: f, i, j, k        !< Counters.

   lower = 0._R8P ; mean = 0._R8P ; rms = 0._R8P ; samples = 0_I4P
   nv = size(vertices, 2) ; nf = size(facets, 2)
   status = DEV_STATUS_BAD_INPUT
   if (nv == 0 .or. nf == 0 .or. size(per_vertex) /= nv .or. .not. (samples_per_area > 0._R8P)) return
   if (any(facets < 1) .or. any(facets > nv)) return
   per_vertex = 0._R8P
   allocate(nsub(nf))
   np = 0
   do f = 1, nf
      call corners(f, a, b, c)
      area = triangle_area(a, b, c)
      nsub(f) = max(1_I4P, int(ceiling(sqrt(area * samples_per_area)), I4P))
      np = np + 3 + nsub(f) ** 2
   enddo
   np = np + nv
   allocate(p(3, np), w(np), d(np))
   w = 0._R8P
   p(:, 1:nv) = vertices
   k = nv
   do f = 1, nf
      call corners(f, a, b, c)
      area = triangle_area(a, b, c)
      p(:, k + 1) = 0.5_R8P * (a + b) ; p(:, k + 2) = 0.5_R8P * (b + c) ; p(:, k + 3) = 0.5_R8P * (c + a)
      k = k + 3
      ! The n^2 sub-triangles of the regular subdivision: (i, j) upward ones and (i, j) downward ones.
      do i = 0, nsub(f) - 1
         do j = 0, nsub(f) - 1 - i
            k = k + 1                                                       ! upward: (i, j), (i+1, j), (i, j+1)
            p(:, k) = lattice(a, b, c, nsub(f), 3._R8P * i + 1._R8P, 3._R8P * j + 1._R8P)
            w(k) = area / nsub(f) ** 2
            if (i + j <= nsub(f) - 2) then                                  ! downward: (i+1, j), (i, j+1), (i+1, j+1)
               k = k + 1
               p(:, k) = lattice(a, b, c, nsub(f), 3._R8P * i + 2._R8P, 3._R8P * j + 2._R8P)
               w(k) = area / nsub(f) ** 2
            endif
         enddo
      enddo
   enddo
   call eval(p(:, 1:k), d(1:k))
   per_vertex = d(1:nv)
   lower = maxval(d(1:k))
   wsum = sum(w(1:k))
   if (wsum > 0._R8P) then
      mean = sum(w(1:k) * d(1:k)) / wsum
      rms = sqrt(sum(w(1:k) * d(1:k) ** 2) / wsum)
   endif
   samples = k
   status = DEV_STATUS_OK
   contains
      subroutine corners(ff, a_, b_, c_)
      !< Corners of facet ff.
      integer(I4P), intent(in)  :: ff
      real(R8P),    intent(out) :: a_(3), b_(3), c_(3)

      a_ = vertices(:, facets(1, ff)) ; b_ = vertices(:, facets(2, ff)) ; c_ = vertices(:, facets(3, ff))
      endsubroutine corners
   endsubroutine deviation_sampled

   pure function lattice(a, b, c, n, s, t) result(x)
   !< The point a + (s (b - a) + t (c - a)) / (3 n): s, t in thirds of the sub-triangle lattice of step 1/n.
   real(R8P),    intent(in) :: a(3), b(3), c(3) !< Corners.
   integer(I4P), intent(in) :: n                !< Subdivisions.
   real(R8P),    intent(in) :: s, t             !< Lattice coordinates, in thirds.
   real(R8P)                :: x(3)             !< Point.

   x = a + (s * (b - a) + t * (c - a)) / (3._R8P * n)
   endfunction lattice

   pure function triangle_area(a, b, c) result(area)
   !< Area of the triangle (a, b, c).
   real(R8P), intent(in) :: a(3), b(3), c(3) !< Corners.
   real(R8P)             :: area             !< Area.
   real(R8P)             :: u(3), v(3)       !< Edges.

   u = b - a ; v = c - a
   area = 0.5_R8P * norm2([u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)])
   endfunction triangle_area
endmodule fossil_deviation
