!< FOSSIL test: one ray hit per crossing of the surface (issue #82).
!<
!< `intersect_ray_all` decides every hit exactly, with one symbolic perturbation of the ray: a ray through a shared
!< edge or vertex crosses exactly one of the facets there, a ray grazing an edge or lying in a face plane crosses none
!< or two, so on a closed surface the parity of the hit count is the inside/outside state of the origin.
!<
!<   1. The +x ray through the centres of two opposite faces of cube-01 (each centre on a facet diagonal, the corner
!<      of four facets): 2 hits, t = 1 and 2 (12 before #82).
!<   2. An origin on the face x = 0: the hit at t = 0 is kept; `exclude=` drops the start facet from `_all` as from
!<      `_first`, the next hit is the far face at t = 1.
!<   3. Parity: rays from lattice points (axis rays along the lines of the marching-cubes grid, so through its
!<      vertices and edges; face diagonals; rays aimed at every mesh vertex; rays in a face plane; generic rays) on
!<      cube-01 and a marching-cubes sphere, both tree kinds: odd hit count <=> `is_point_inside`. Also `_first` =
!<      `_all(1)` and `_any` = (hits > 0) on every ray.
!<   4. The octree's traversal finds every hit: on the dragon, 300 rays, octree = BVH. An octree node holds
!<      the facets whose centroid it contains, which stick out of its box: the slab test pads the boxes by the
!<      facets' reach (without the pad, 49 of 3000 such rays lost hits against the flat O(n) loop).
program fossil_test_ray_crossings
use fossil, only : surface_stl_object, ray_hit_t, RAY_STATUS_OK, extract_isosurface, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH, &
                   SIGN_FAST_WINDING_NUMBER
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)        :: dragon, dragon_bvh
type(surface_stl_object)        :: cube, cube_oct, sphere, sphere_oct
type(facet_object), allocatable :: sphere_facets(:), copy(:)
type(ray_hit_t),    allocatable :: hits(:), flat(:)
type(ray_hit_t)                 :: first
type(vector_R8P)                :: o, d, lo, hi
real(R8P),          allocatable :: values(:,:,:)
real(R8P)                       :: x, y, z
integer(I4P)                    :: status, i, j, k, start, nbad
logical                         :: has_hit, passed(5)

passed = .false.
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call cube_oct%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)

! 1. through the face centres
call cube%intersect_ray_all(ray_origin=vector_R8P(-1._R8P, 0.5_R8P, 0.5_R8P), ray_direction=ex_R8P, hits=hits, status=status)
print '(A,I0,A,*(F6.3))', 'face centres: hits ', size(hits), ', t', hits%t
passed(1) = status == RAY_STATUS_OK .and. size(hits) == 2
if (passed(1)) passed(1) = hits(1)%t == 1._R8P .and. hits(2)%t == 2._R8P

! 2. origin on the surface, exclude
o = vector_R8P(0._R8P, 0.3_R8P, 0.6_R8P)
call cube%intersect_ray_all(ray_origin=o, ray_direction=ex_R8P, hits=hits)
print '(A,I0,A,*(F6.3))', 'origin on x = 0: hits ', size(hits), ', t', hits%t
if (size(hits) == 2) then
   passed(2) = hits(1)%t == 0._R8P .and. hits(2)%t == 1._R8P
   start = hits(1)%facet_id
   call cube%intersect_ray_all(ray_origin=o, ray_direction=ex_R8P, hits=hits, exclude=start)
   call cube%intersect_ray_first(ray_origin=o, ray_direction=ex_R8P, hit=first, has_hit=has_hit, exclude=start)
   print '(A,I0,A,*(F6.3))', '  exclude the start facet: hits ', size(hits), ', t', hits%t
   passed(2) = passed(2) .and. size(hits) == 1 .and. has_hit
   if (passed(2)) passed(2) = hits(1)%t == 1._R8P .and. hits(1)%facet_id /= start .and. first%facet_id == hits(1)%facet_id
endif

! 3. parity on the cube and on a marching-cubes sphere (radius 0.9 on a 0.25 grid: no grid node on the sphere)
passed(3) = parity_ok(cube, 'cube-01 BVH', 0.5_R8P, -0.5_R8P, 1.5_R8P, 11_I4P)
passed(3) = parity_ok(cube_oct, 'cube-01 octree', 1._R8P, -0.5_R8P, 1.5_R8P, 11_I4P) .and. passed(3)
allocate(values(13, 13, 13))
do k = 1, 13 ; do j = 1, 13 ; do i = 1, 13
   x = -1.5_R8P + (i - 1) * 0.25_R8P ; y = -1.5_R8P + (j - 1) * 0.25_R8P ; z = -1.5_R8P + (k - 1) * 0.25_R8P
   values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 0.9_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=values, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=sphere_facets, status=status)
copy = sphere_facets
call sphere%adopt_facets(facets=sphere_facets)
call sphere_oct%adopt_facets(facets=copy, aabb_tree_kind=AABB_TREE_OCTREE)
passed(4) = parity_ok(sphere, 'sphere BVH', 0.75_R8P, -1.25_R8P, 1._R8P, 401_I4P)
passed(4) = parity_ok(sphere_oct, 'sphere octree', 1.5_R8P, -1.25_R8P, 1.75_R8P, 401_I4P) .and. passed(4)

! 4. octree = BVH on the dragon
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call dragon_bvh%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
lo = dragon%get_bmin() ; hi = dragon%get_bmax()
nbad = 0
do k = 1, 300
   o = vector_R8P(lo%x + frac(0.37_R8P * k) * (hi%x - lo%x), lo%y + frac(0.61_R8P * k + 0.1_R8P) * (hi%y - lo%y), &
                  lo%z + frac(0.83_R8P * k + 0.2_R8P) * (hi%z - lo%z))
   d = vector_R8P(sin(1.3_R8P * k), cos(2.9_R8P * k), sin(0.7_R8P * k + 1._R8P))
   call dragon%intersect_ray_all(ray_origin=o, ray_direction=d, hits=hits)
   call dragon_bvh%intersect_ray_all(ray_origin=o, ray_direction=d, hits=flat)
   if (size(hits) /= size(flat)) then
      nbad = nbad + 1
   elseif (any(hits%facet_id /= flat%facet_id)) then
      nbad = nbad + 1
   endif
enddo
print '(A,I0)', 'dragon octree vs BVH, 300 rays: mismatches ', nbad
passed(5) = nbad == 0

print '(A,*(L2))', 'per-case results:', passed
print '(A,L1)', 'Are all tests passed? ', all(passed)
if (.not. all(passed)) error stop 1

contains
   pure function frac(x) result(f)
   !< Fractional part.
   real(R8P), intent(in) :: x !< Argument.
   real(R8P)             :: f !< x - floor(x).

   f = x - floor(x)
   endfunction frac

   function parity_ok(s, name, h, lo, hi, stride) result(ok)
   !< Rays from the lattice points lo:h:hi (those not on the surface) in fixed directions and aimed at every `stride`-th
   !< vertex: odd hit count <=> inside, `_first` = `_all(1)`, `_any` = (hits > 0).
   type(surface_stl_object), intent(in), target :: s !< Closed surface.
   character(*),             intent(in) :: name     !< Label.
   real(R8P),                intent(in) :: h, lo, hi !< Lattice.
   integer(I4P),             intent(in) :: stride   !< Aim at every stride-th vertex.
   logical                              :: ok       !< Result.
   type(vector_R8P), allocatable        :: dirs(:)  !< Directions.
   type(vector_R8P)                     :: p, d     !< Origin, direction.
   type(ray_hit_t), allocatable         :: hs(:)    !< Hits.
   type(ray_hit_t)                      :: f1       !< First hit.
   logical                              :: inside, has, any_hit
   integer(I4P)                         :: n, i, j, k, m, nrays, nbad, nfirst, nany, nv, v
   real(R8P)                            :: r(3)

   n = nint((hi - lo) / h)
   nv = 3 * s%get_facets_number()
   allocate(dirs(26 + 4))
   m = 0
   do i = -1, 1 ; do j = -1, 1 ; do k = -1, 1           ! axes, face diagonals, body diagonals
      if (i == 0 .and. j == 0 .and. k == 0) cycle
      m = m + 1 ; dirs(m) = vector_R8P(real(i, R8P), real(j, R8P), real(k, R8P))
   enddo ; enddo ; enddo
   do m = 27, 30                                         ! generic
      r = [sin(1.3_R8P * m), cos(2.9_R8P * m), sin(0.7_R8P * m + 1._R8P)]
      dirs(m) = vector_R8P(r(1), r(2), r(3))
   enddo
   nrays = 0 ; nbad = 0 ; nfirst = 0 ; nany = 0
   do k = 0, n ; do j = 0, n ; do i = 0, n
      p = vector_R8P(lo + i * h, lo + j * h, lo + k * h)
      if (s%distance(point=p) == 0._R8P) cycle          ! on the surface: neither inside nor outside
      inside = s%is_point_inside(point=p, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
      do m = 1, size(dirs) + nv / stride                ! the fixed directions, then aimed at vertices
         if (m <= size(dirs)) then
            d = dirs(m)
         else
            v = (m - size(dirs) - 1) * stride
            d = vertex_of(s, v / 3 + 1, mod(v, 3) + 1) - p
            if (d%x == 0._R8P .and. d%y == 0._R8P .and. d%z == 0._R8P) cycle
         endif
         call s%intersect_ray_all(ray_origin=p, ray_direction=d, hits=hs)
         nrays = nrays + 1
         if ((mod(size(hs), 2) == 1) .neqv. inside) then
            nbad = nbad + 1
            if (nbad <= 3) print '(2X,A,3F7.3,A,3F7.3,A,I0,A,L1)', 'parity: origin', p%x, p%y, p%z, ' dir', d%x, d%y, d%z, &
                                                                ' hits ', size(hs), ' inside ', inside
         endif
         call s%intersect_ray_first(ray_origin=p, ray_direction=d, hit=f1, has_hit=has)
         if (has .neqv. size(hs) > 0) then
            nfirst = nfirst + 1
         elseif (has) then
            if (f1%facet_id /= hs(1)%facet_id) nfirst = nfirst + 1
         endif
         call s%intersect_ray_any(ray_origin=p, ray_direction=d, found=any_hit)
         if (any_hit .neqv. size(hs) > 0) nany = nany + 1
      enddo
   enddo ; enddo ; enddo
   print '(A,A,I0,A,I0,A,I0,A,I0)', name, ': rays ', nrays, ', parity mismatches ', nbad, ', first /= all(1) ', nfirst, &
                                     ', any /= all ', nany
   ok = nbad == 0 .and. nfirst == 0 .and. nany == 0
   endfunction parity_ok

   function vertex_of(s, f, v) result(p)
   !< Vertex v of facet f.
   type(surface_stl_object), intent(in), target :: s !< Surface.
   integer(I4P),             intent(in) :: f, v !< Facet, vertex.
   type(vector_R8P)                     :: p    !< Vertex.
   type(facet_object), pointer          :: fc   !< Facet.

   fc => s%facet_at(f)
   p = fc%vertex(v)
   endfunction vertex_of
endprogram fossil_test_ray_crossings
