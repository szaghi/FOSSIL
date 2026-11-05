!< FOSSIL, surface deviation: the one-sided Hausdorff distance between two surfaces (issue #60, roadmap #58 WP2).
module fossil_deviation
!< FOSSIL, surface deviation: the one-sided Hausdorff distance between two surfaces (issue #60, roadmap #58 WP2).
!<
!< h(A, B) = sup over x in A of d(x, B). A is given as a triangle mesh (vertices and facets), B only through an
!< evaluator, an extension of `deviation_evaluator_t` whose `eval` returns d(x, B) and the closest points for a batch of
!< points (the caller's BVH query), so this module depends on no surface type. (An abstract type rather than a
!< procedure argument: passing an internal procedure would make gfortran build trampolines on an executable stack.)
!<
!< Lower bound: the largest distance over deterministic samples of A: every vertex, every facet's edge midpoints, and the
!< centroids of the n_f^2 sub-triangles of facet f (n_f = ceil(sqrt(area_f * samples_per_area))). The same centroids,
!< weighted by area_f / n_f^2 (the centroid rule on the sub-triangles: exact for a linear function), give the mean and
!< the rms of d over A: estimates, not bounds.
use penf, only : I4P, R8P

implicit none
private
public :: deviation_evaluator_t
public :: deviation_sampled, deviation_certified, triangle_rho
public :: DEV_STATUS_OK, DEV_STATUS_BAD_INPUT, DEV_STATUS_NOT_CONVERGED

integer(I4P), parameter :: DEV_STATUS_OK        = 0_I4P !< Success.
integer(I4P), parameter :: DEV_STATUS_BAD_INPUT = 1_I4P !< Empty surface, or a size mismatch.
integer(I4P), parameter :: DEV_STATUS_NOT_CONVERGED = 2_I4P !< Work cap reached: the bounds are valid but wider.

type, abstract :: deviation_evaluator_t
   !< The surface B, seen through its batched distance query.
   contains
      procedure(eval_interface), deferred, pass(self) :: eval !< Distances and closest points of a batch.
endtype deviation_evaluator_t

abstract interface
   subroutine eval_interface(self, points, distances, closest, facets)
   !< Distances d(points(:, i), B) of a batch of points, the closest points of B and the convex pieces of B (facets)
   !< they lie on (0 if unknown).
   import :: deviation_evaluator_t, I4P, R8P
   class(deviation_evaluator_t), intent(in)  :: self         !< Evaluator.
   real(R8P),                    intent(in)  :: points(:,:)  !< Points [3, n].
   real(R8P),                    intent(out) :: distances(:) !< Unsigned distances [n].
   real(R8P),                    intent(out) :: closest(:,:) !< Closest points of B [3, n].
   integer(I4P),                 intent(out) :: facets(:)    !< Facet of B of each closest point [n].
   endsubroutine eval_interface
endinterface

contains
   subroutine deviation_sampled(vertices, facets, eval, samples_per_area, lower, mean, rms, per_vertex, samples, status, &
                                vertex_feet, vertex_facets)
   !< Sampled one-sided deviation of the mesh (vertices, facets) from the surface behind `eval`.
   real(R8P),    intent(in)            :: vertices(:,:)     !< Vertices of A [3, nv].
   integer(I4P), intent(in)            :: facets(:,:)       !< Facets of A, vertex indices [3, nf].
   class(deviation_evaluator_t), intent(in) :: eval    !< Batched distance to B.
   real(R8P),    intent(in)            :: samples_per_area  !< Interior samples per unit area (> 0).
   real(R8P),    intent(out)           :: lower             !< Largest sampled distance (a lower bound of h).
   real(R8P),    intent(out)           :: mean              !< Area-weighted mean of d over A (estimate).
   real(R8P),    intent(out)           :: rms               !< Area-weighted rms of d over A (estimate).
   real(R8P),    intent(out)           :: per_vertex(:)     !< d at every vertex [nv].
   integer(I4P), intent(out)           :: samples           !< Number of points evaluated.
   integer(I4P), intent(out)           :: status            !< DEV_STATUS_OK or DEV_STATUS_BAD_INPUT.
   real(R8P),    intent(out), optional :: vertex_feet(:,:)  !< Closest points of B to the vertices [3, nv].
   integer(I4P), intent(out), optional :: vertex_facets(:)  !< Facets of B of those closest points [nv].
   real(R8P),    allocatable           :: p(:,:)            !< Sample points.
   real(R8P),    allocatable           :: q(:,:)            !< Closest points of the samples.
   integer(I4P), allocatable           :: g(:)              !< Facets of B of the closest points.
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
   allocate(p(3, np), w(np), d(np), q(3, np), g(np))
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
   call eval%eval(p(:, 1:k), d(1:k), q(:, 1:k), g(1:k))
   per_vertex = d(1:nv)
   if (present(vertex_feet)) vertex_feet = q(:, 1:nv)
   if (present(vertex_facets)) vertex_facets = g(1:nv)
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

   subroutine deviation_certified(vertices, facets, per_vertex, vertex_feet, vertex_facets, eval, tolerance, &
                                  max_triangles, lower, upper, triangles, status)
   !< Certified bounds lower <= h(A, B) <= upper by branch and bound over the facets of A (issue #60 step 2).
   !<
   !< For a triangle T with corners c_k, their distances d_k = d(c_k, B) and closest points q_k on B, every x in T has
   !<   d(x, B) <= d_k + |x - c_k| (d is 1-Lipschitz), so  sup_T d <= max_k d_k + rho_T  (rho_T: `triangle_rho`), and
   !<   d(x, B) <= |x - q_k|, convex in x, so            sup_T d <= min_k max_j |c_j - q_k|, and
   !<   if q_1, q_2, q_3 lie on one facet F of B (convex): x = sum l_k c_k has y = sum l_k q_k in F, so
   !<   d(x, B) <= |x - y| <= sum l_k |c_k - q_k|, and   sup_T d <= max_k d_k   (no refinement needed).
   !< The smallest is the triangle's bound U_T; the second converges quadratically where d is smooth, the third settles
   !< at once every triangle facing a single facet of B (coincident or parallel surfaces). A triangle
   !< with U_T <= lower + tolerance is settled; the others are split into four at their edge midpoints, whose distances
   !< raise `lower`. `upper` is the largest bound among the settled triangles and those left when `max_triangles` (the
   !< total of triangles created) is reached, `status` then DEV_STATUS_NOT_CONVERGED. Rounds run in a fixed order: the
   !< result does not depend on the thread count.
   real(R8P),    intent(in)            :: vertices(:,:)     !< Vertices of A [3, nv].
   integer(I4P), intent(in)            :: facets(:,:)       !< Facets of A [3, nf].
   real(R8P),    intent(in)            :: per_vertex(:)     !< d at the vertices [nv].
   real(R8P),    intent(in)            :: vertex_feet(:,:)  !< Closest points of B to the vertices [3, nv].
   integer(I4P), intent(in)            :: vertex_facets(:)  !< Facets of B of those closest points [nv].
   class(deviation_evaluator_t), intent(in) :: eval    !< Batched distance to B.
   real(R8P),    intent(in)            :: tolerance         !< Target upper - lower (> 0).
   integer(I4P), intent(in)            :: max_triangles     !< Work cap.
   real(R8P),    intent(inout)         :: lower             !< In: a lower bound (the sampled one); out: raised.
   real(R8P),    intent(out)           :: upper             !< Upper bound of h.
   integer(I4P), intent(out)           :: triangles         !< Triangles created.
   integer(I4P), intent(out)           :: status            !< DEV_STATUS_OK or DEV_STATUS_NOT_CONVERGED.
   real(R8P),    allocatable           :: c(:,:,:), q(:,:,:), d(:,:) !< Active triangles: corners, feet, distances.
   real(R8P),    allocatable           :: c2(:,:,:), q2(:,:,:), d2(:,:) !< Next round.
   real(R8P),    allocatable           :: u(:)              !< Bounds of the active triangles.
   logical,      allocatable           :: act(:)            !< Triangles split in this round.
   real(R8P),    allocatable           :: mp(:,:), md(:), mq(:,:) !< Midpoints, their distances and feet.
   integer(I4P), allocatable           :: g(:,:), g2(:,:), mg(:)  !< Facets of B of the feet.
   integer(I4P)                        :: n, na, t, k, j, m !< Counters.
   integer(I4P), parameter             :: E(2,3) = reshape([1, 2, 2, 3, 3, 1], [2, 3]) !< Edges (corner pairs).

   status = DEV_STATUS_OK
   n = size(facets, 2)
   allocate(c(3, 3, n), q(3, 3, n), d(3, n), g(3, n))
   do t = 1, n
      do k = 1, 3
         c(:, k, t) = vertices(:, facets(k, t)) ; q(:, k, t) = vertex_feet(:, facets(k, t)) ; d(k, t) = per_vertex(facets(k, t))
         g(k, t) = vertex_facets(facets(k, t))
      enddo
   enddo
   triangles = n
   upper = lower
   do
      allocate(u(n), act(n))
      do t = 1, n
         u(t) = triangle_bound(c(:,:,t), d(:,t), q(:,:,t), g(:,t))
      enddo
      act = u > lower + tolerance
      na = count(act)
      if (na < n) upper = max(upper, maxval(u, mask=.not. act))
      if (na == 0) exit
      if (triangles + 4 * na > max_triangles) then
         upper = max(upper, maxval(u))
         status = DEV_STATUS_NOT_CONVERGED
         exit
      endif
      allocate(mp(3, 3 * na), md(3 * na), mq(3, 3 * na), mg(3 * na))
      m = 0
      do t = 1, n
         if (.not. act(t)) cycle
         do j = 1, 3
            m = m + 1 ; mp(:, m) = 0.5_R8P * (c(:, E(1,j), t) + c(:, E(2,j), t))
         enddo
      enddo
      call eval%eval(mp, md, mq, mg)
      lower = max(lower, maxval(md))
      allocate(c2(3, 3, 4 * na), q2(3, 3, 4 * na), d2(3, 4 * na), g2(3, 4 * na))
      m = 0 ; k = 0
      do t = 1, n
         if (.not. act(t)) cycle
         k = k + 1
         call split(c(:,:,t), d(:,t), q(:,:,t), g(:,t), mp(:, 3*k-2:3*k), md(3*k-2:3*k), mq(:, 3*k-2:3*k), &
                    mg(3*k-2:3*k), c2(:,:,4*k-3:4*k), d2(:,4*k-3:4*k), q2(:,:,4*k-3:4*k), g2(:,4*k-3:4*k))
      enddo
      triangles = triangles + 4 * na
      call move_alloc(c2, c) ; call move_alloc(q2, q) ; call move_alloc(d2, d) ; call move_alloc(g2, g)
      deallocate(u, act, mp, md, mq, mg)
      n = 4 * na
   enddo
   upper = max(upper, lower)
   endsubroutine deviation_certified

   pure function triangle_bound(c, d, q, g) result(u)
   !< Upper bound of d(x, B) over the triangle with corners c, their distances d, closest points q and the facets g of B
   !< those lie on (see deviation_certified): max_k d_k if the three feet share a facet, else
   !< min(max_k d_k + rho_T, min_k max_j |c_j - q_k|).
   real(R8P),    intent(in) :: c(3,3) !< Corners (columns).
   real(R8P),    intent(in) :: d(3)   !< Corner distances.
   real(R8P),    intent(in) :: q(3,3) !< Corner closest points (columns).
   integer(I4P), intent(in) :: g(3)   !< Facets of B of the closest points (0: unknown).
   real(R8P)                :: u      !< Bound.
   integer                  :: k      !< Counter.

   if (g(1) > 0 .and. g(1) == g(2) .and. g(1) == g(3)) then
      u = maxval(d)
      return
   endif
   u = maxval(d) + triangle_rho(c(:,1), c(:,2), c(:,3))
   do k = 1, 3
      u = min(u, max(norm2(c(:,1) - q(:,k)), norm2(c(:,2) - q(:,k)), norm2(c(:,3) - q(:,k))))
   enddo
   endfunction triangle_bound

   pure function triangle_rho(a, b, c) result(rho)
   !< rho_T = max over x in T of the distance from x to its nearest corner. With the angles at a, b, c:
   !< - no angle obtuse: the circumcentre lies in T and rho is the circumradius |ab| |bc| |ca| / (4 area);
   !< - obtuse at c (likewise at a or b): the maximum lies on the opposite side ab, where the perpendicular bisectors of
   !<   ca and cb cross it: rho = max(|ca| / (2 cos A), |cb| / (2 cos B)).
   !< Capped by the longest side (any point of T is within it of every corner), which also covers degenerate triangles.
   real(R8P), intent(in) :: a(3), b(3), c(3) !< Corners.
   real(R8P)             :: rho              !< Bound.
   real(R8P)             :: la, lb, lc       !< Side lengths opposite a, b, c.
   real(R8P)             :: ca, cb, cc       !< Dot products at the corners.
   real(R8P)             :: area             !< Area.

   la = norm2(b - c) ; lb = norm2(c - a) ; lc = norm2(a - b)
   ca = dot_product(b - a, c - a) ; cb = dot_product(a - b, c - b) ; cc = dot_product(a - c, b - c)
   area = triangle_area(a, b, c)
   rho = max(la, lb, lc)
   if (cc < 0._R8P) then                       ! obtuse at c: the bisectors of ca, cb meet ab
      if (ca > 0._R8P .and. cb > 0._R8P) rho = min(rho, max(lb * lb * lc / (2._R8P * ca), la * la * lc / (2._R8P * cb)))
   elseif (ca < 0._R8P) then                   ! obtuse at a: the bisectors of ab, ac meet bc
      if (cb > 0._R8P .and. cc > 0._R8P) rho = min(rho, max(lc * lc * la / (2._R8P * cb), lb * lb * la / (2._R8P * cc)))
   elseif (cb < 0._R8P) then                   ! obtuse at b: the bisectors of ba, bc meet ca
      if (ca > 0._R8P .and. cc > 0._R8P) rho = min(rho, max(lc * lc * lb / (2._R8P * ca), la * la * lb / (2._R8P * cc)))
   elseif (area > 0._R8P) then
      rho = min(rho, la * lb * lc / (4._R8P * area))
   endif
   endfunction triangle_rho

   subroutine split(c, d, q, g, mc, md, mq, mg, cc, dd, qq, gg)
   !< The four children of a triangle at its edge midpoints m12, m23, m31 (in the order of `E`).
   real(R8P),    intent(in)  :: c(3,3), d(3), q(3,3)          !< Parent corners, distances, feet.
   integer(I4P), intent(in)  :: g(3)                          !< Parent feet's facets.
   real(R8P),    intent(in)  :: mc(3,3), md(3), mq(3,3)       !< Midpoints m12, m23, m31: points, distances, feet.
   integer(I4P), intent(in)  :: mg(3)                         !< Midpoints' feet's facets.
   real(R8P),    intent(out) :: cc(3,3,4), dd(3,4), qq(3,3,4) !< Children.
   integer(I4P), intent(out) :: gg(3,4)                       !< Children's feet's facets.

   call set(1, c(:,1), d(1), q(:,1), g(1), mc(:,1), md(1), mq(:,1), mg(1), mc(:,3), md(3), mq(:,3), mg(3))
   call set(2, mc(:,1), md(1), mq(:,1), mg(1), c(:,2), d(2), q(:,2), g(2), mc(:,2), md(2), mq(:,2), mg(2))
   call set(3, mc(:,3), md(3), mq(:,3), mg(3), mc(:,2), md(2), mq(:,2), mg(2), c(:,3), d(3), q(:,3), g(3))
   call set(4, mc(:,1), md(1), mq(:,1), mg(1), mc(:,2), md(2), mq(:,2), mg(2), mc(:,3), md(3), mq(:,3), mg(3))
   contains
      subroutine set(i, p1, e1, f1, h1, p2, e2, f2, h2, p3, e3, f3, h3)
      !< Child i from three corners with their distances, feet and feet's facets.
      integer,      intent(in) :: i
      real(R8P),    intent(in) :: p1(3), e1, f1(3), p2(3), e2, f2(3), p3(3), e3, f3(3)
      integer(I4P), intent(in) :: h1, h2, h3

      cc(:,1,i) = p1 ; cc(:,2,i) = p2 ; cc(:,3,i) = p3
      dd(:,i) = [e1, e2, e3]
      qq(:,1,i) = f1 ; qq(:,2,i) = f2 ; qq(:,3,i) = f3
      gg(:,i) = [h1, h2, h3]
      endsubroutine set
   endsubroutine split

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
