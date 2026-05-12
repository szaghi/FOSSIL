!< FOSSIL, vertex pool object.
!<
!< Authoritative store of unique-vertex coordinates plus the bidirectional
!< mapping (facet, local_v) <-> pool_id. Built from a facet array via strict-EPS
!< coincidence; consumed by `build_connectivity` (canonical-id edge classifier),
!< `remove_duplicate_facets` (orientation-agnostic dedup keys), and
!< `compute_pseudo_normals_via_pool` (incident-facet enumeration around each
!< unique vertex).
!<
!< Facets carry `vertex_id(3)` as the structural source of truth (issue #5
!< stage 3); they also keep `vertex(3)` as a coordinate cache for hot kernels.
!<
!< API:
!<   call pool%initialize_from_facets(facet, status)
!<   n   = pool%vertex_count()
!<   nf  = pool%facets_count()
!<   v   = pool%coord(vid)
!<   vid = pool%facet_vid(facet_id, local_v)
!<   nk  = pool%facets_at_count(vid)
!<   call pool%facets_at(vid, k, facet_id, local_v)
!<   call pool%destroy
!<
!< Build: O((3N)^2) pairwise tests via union-by-rank + path-compression union-find,
!< plus an O(N) CSR inverted-index pass. The pairwise cost matches the pre-pool
!< coincidence work it subsumes. Stage 5 will replace the pairwise step with a
!< spatial hash for sub-linear coincidence checks.

module fossil_vertex_pool_object

use fossil_facet_object, only : facet_object
use fossil_utils,        only : EPS
use penf,                only : I4P, I8P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: vertex_pool_object

type :: vertex_pool_object
   !< Unique-vertex coordinate pool with (facet,local_v) -> pool_id mapping
   !< and the inverted index pool_id -> [(facet_id, local_v), ...] (stage 3b).
   !<
   !< Invariants after `initialize_from_facets`:
   !<   - `coord(1:n_vertices)` contains pairwise-distinct vertex coordinates
   !<     (pairwise distance > EPS on at least one component).
   !<   - `facet_to_pool(1:3, 1:n_facets)` maps every (local_v, facet_id) onto
   !<     a valid pool index in [1, n_vertices].
   !<   - For every (f, v): `facet_to_pool(v, f) == facet(f)%vertex_id(v)` and
   !<     `coord(facet_to_pool(v, f))` bit-equals `facet(f)%vertex(v)` (the cached
   !<     coordinate is materialized from one of the EPS-coincident inputs at build).
   !<   - The inverted index is stored CSR-style:
   !<       at_offset(vid)..at_offset(vid+1)-1 are the rows of at_pairs(:,:)
   !<       holding (facet_id, local_v) pairs that reference pool vertex vid.
   private
   type(vector_R8P), allocatable :: coord_(:)         !< Unique vertex coordinates.
   integer(I4P),     allocatable :: facet_to_pool(:,:) !< (local_v=1..3, facet_id=1..n_facets) -> pool id.
   integer(I4P),     allocatable :: at_offset(:)      !< CSR offsets, size n_vertices+1.
   integer(I4P),     allocatable :: at_pairs(:,:)     !< CSR payload, shape (2, 3*n_facets); row 1 = facet_id, row 2 = local_v.
   integer(I4P)                  :: n_vertices = 0    !< Pool size.
   integer(I4P)                  :: n_facets   = 0    !< Source facet count.
   logical                       :: is_initialized = .false.
   contains
      procedure, pass(self) :: initialize_from_facets
      procedure, pass(self) :: vertex_count
      procedure, pass(self) :: facets_count
      procedure, pass(self) :: coord
      procedure, pass(self) :: facet_vid
      procedure, pass(self) :: facets_at_count !< Number of (facet, local_v) pairs referencing pool vertex vid.
      procedure, pass(self) :: facets_at       !< k-th (facet_id, local_v) pair referencing pool vertex vid.
      procedure, pass(self) :: get_is_initialized
      procedure, pass(self) :: destroy
      final :: finalize
endtype vertex_pool_object

contains

   subroutine initialize_from_facets(self, facet, status, use_union_find)
   !< Build the pool from a facet array via EPS coincidence.
   !<
   !< Default path: spatial-hash builder (issue #5 stage 5). Single linear pass
   !< over (3 * nf) candidate vertices using an open-addressing hash table keyed
   !< on a 3D bucket index. Per insert, probe the 27-neighbour cell cluster and
   !< EPS-test against existing entries; reuse pool id on match, allocate new on
   !< miss. Average O(N), worst case O(N * k) where k is the per-bucket load.
   !<
   !< Fallback path: pairwise O((3N)^2) union-by-rank with path compression.
   !< Selected via `use_union_find=.true.`; retained to allow A/B verification
   !< of the hash path against a reference implementation.
   class(vertex_pool_object),   intent(inout)         :: self           !< Vertex pool.
   type(facet_object),          intent(in)            :: facet(:)       !< Source facets.
   integer(I4P),                intent(out), optional :: status         !< 0=ok, 1=alloc fail.
   logical,                     intent(in),  optional :: use_union_find !< Force the O((3N)^2) reference path.
   integer(I4P)                                       :: istat
   logical                                            :: use_uf_

   call self%destroy
   if (present(status)) status = 0_I4P

   self%n_facets = size(facet)
   if (self%n_facets <= 0) then
      self%is_initialized = .true.
      return
   endif

   use_uf_ = .false.
   if (present(use_union_find)) use_uf_ = use_union_find

   if (use_uf_) then
      call build_via_union_find(self, facet, istat)
   else
      call build_via_spatial_hash(self, facet, istat)
   endif
   if (istat /= 0) then
      if (present(status)) status = 1_I4P
      return
   endif

   ! Stage 3b: build the inverted index pool_id -> [(facet_id, local_v), ...] in CSR.
   call build_inverted_index(self, istat)
   if (istat /= 0 .and. present(status)) status = 1_I4P

   self%is_initialized = .true.
   endsubroutine initialize_from_facets

   subroutine build_via_union_find(self, facet, istat)
   !< Reference O((3N)^2) pool builder via union-by-rank + path compression.
   !< Kept for A/B verification of the spatial-hash builder; not the default.
   type(vertex_pool_object), intent(inout) :: self
   type(facet_object),       intent(in)    :: facet(:)
   integer(I4P),             intent(out)   :: istat
   integer(I4P), allocatable               :: parent(:), rank_(:), root_to_pid(:)
   integer(I4P)                            :: nf, nv_total
   integer(I4P)                            :: f1, f2, v1, v2, gid1, root, pid

   istat = 0_I4P
   nf = self%n_facets
   nv_total = 3 * nf
   allocate(parent(nv_total), rank_(nv_total), stat=istat)
   if (istat /= 0) return
   do gid1 = 1, nv_total
      parent(gid1) = gid1
      rank_(gid1)  = 0_I4P
   enddo

   ! Pairwise EPS union. Inner loop walks f2 > f1 only, but for each (f1, v1) we
   ! also union with (f1, v2 > v1) -- a degenerate triangle can have coincident
   ! local vertices and we want to detect that.
   do f1 = 1, nf
      do v1 = 1, 3
         gid1 = (f1 - 1) * 3 + v1
         do v2 = v1 + 1, 3
            if (vertices_eps_match(facet(f1)%vertex(v1), facet(f1)%vertex(v2))) &
               call uf_union(gid1, (f1 - 1) * 3 + v2, parent, rank_)
         enddo
         do f2 = f1 + 1, nf
            do v2 = 1, 3
               if (vertices_eps_match(facet(f1)%vertex(v1), facet(f2)%vertex(v2))) &
                  call uf_union(gid1, (f2 - 1) * 3 + v2, parent, rank_)
            enddo
         enddo
      enddo
   enddo

   allocate(root_to_pid(nv_total), self%facet_to_pool(3, nf), stat=istat)
   if (istat /= 0) then
      deallocate(parent, rank_)
      return
   endif
   root_to_pid = 0_I4P
   pid = 0_I4P
   do f1 = 1, nf
      do v1 = 1, 3
         gid1 = (f1 - 1) * 3 + v1
         root = uf_find(gid1, parent)
         if (root_to_pid(root) == 0_I4P) then
            pid = pid + 1_I4P
            root_to_pid(root) = pid
         endif
         self%facet_to_pool(v1, f1) = root_to_pid(root)
      enddo
   enddo
   self%n_vertices = pid

   allocate(self%coord_(self%n_vertices), stat=istat)
   if (istat /= 0) then
      deallocate(parent, rank_, root_to_pid)
      return
   endif
   do f1 = 1, nf
      do v1 = 1, 3
         pid = self%facet_to_pool(v1, f1)
         self%coord_(pid) = facet(f1)%vertex(v1)
      enddo
   enddo

   deallocate(parent, rank_, root_to_pid)
   endsubroutine build_via_union_find

   subroutine build_via_spatial_hash(self, facet, istat)
   !< Default O(N)-expected pool builder via open-addressing spatial hash.
   !<
   !< Algorithm:
   !<   1. Compute the facet-array bbox; choose `cell` adaptively so each occupied
   !<      bucket holds ~1 vertex on average. Floor `cell` at `8 * EPS` so the
   !<      27-neighbour probe is correctness-preserving (any two EPS-coincident
   !<      coordinates fall in the same bucket or in immediate neighbours).
   !<   2. Open-addressing hash table sized as the next power of two >= 4 * (3 * nf).
   !<      Each table slot stores a packed I8P bucket key plus the head of a linked
   !<      list of pool ids hashed to that bucket; `next_in_bucket(pid)` chains the
   !<      list. Linear probing on collision.
   !<   3. For each (f, v) in order: compute the 3D bucket index, probe the
   !<      3x3x3 cell cluster, EPS-test against existing entries, reuse pool id on
   !<      match or insert a new entry. Update `self%facet_to_pool(v, f)`.
   !<   4. Materialize `self%coord_` from `facet` using the assigned ids.
   type(vertex_pool_object), intent(inout) :: self
   type(facet_object),       intent(in)    :: facet(:)
   integer(I4P),             intent(out)   :: istat
   integer(I4P)                            :: nf, nv_total, nbuckets, mask
   integer(I4P)                            :: f, v, pid
   integer(I4P)                            :: ix, iy, iz
   real(R8P)                               :: bmin_x, bmin_y, bmin_z
   real(R8P)                               :: bmax_x, bmax_y, bmax_z
   real(R8P)                               :: ext_x, ext_y, ext_z, cell, inv_cell
   real(R8P)                               :: vx, vy, vz
   integer(I8P), allocatable               :: tab_key(:)       !< Hash-table slot key (0 = empty).
   integer(I4P), allocatable               :: tab_head(:)      !< Hash-table slot head of chain (pool id; 0 = empty).
   integer(I4P), allocatable               :: next_in_chain(:) !< next_in_chain(pid) = next pool id in same bucket (0 = end).
   type(vector_R8P), allocatable           :: coord_buf(:)     !< Growable coord buffer (pid-indexed).
   integer(I4P)                            :: cap              !< Current capacity of coord_buf and next_in_chain.

   istat = 0_I4P
   nf = self%n_facets
   nv_total = 3 * nf

   ! Pass 1: bbox extents.
   bmin_x = facet(1)%vertex(1)%x ; bmin_y = facet(1)%vertex(1)%y ; bmin_z = facet(1)%vertex(1)%z
   bmax_x = bmin_x               ; bmax_y = bmin_y               ; bmax_z = bmin_z
   do f = 1, nf
      do v = 1, 3
         vx = facet(f)%vertex(v)%x
         vy = facet(f)%vertex(v)%y
         vz = facet(f)%vertex(v)%z
         if (vx < bmin_x) bmin_x = vx
         if (vy < bmin_y) bmin_y = vy
         if (vz < bmin_z) bmin_z = vz
         if (vx > bmax_x) bmax_x = vx
         if (vy > bmax_y) bmax_y = vy
         if (vz > bmax_z) bmax_z = vz
      enddo
   enddo
   ext_x = max(bmax_x - bmin_x, EPS)
   ext_y = max(bmax_y - bmin_y, EPS)
   ext_z = max(bmax_z - bmin_z, EPS)

   ! Adaptive cell size: target ~1 vertex per bucket along the longest axis.
   ! Floor at 8 * EPS so EPS-coincident points cannot straddle more than one cell
   ! boundary on any axis (the 27-neighbour probe handles up to one boundary).
   cell = max(max(ext_x, max(ext_y, ext_z)) / real(max(int(nv_total**(1.0_R8P/3.0_R8P), I4P), 1_I4P), R8P), 8._R8P * EPS)
   inv_cell = 1._R8P / cell

   ! Hash table sized to the next power of two >= 4 * nv_total, with a floor of 16.
   nbuckets = 16_I4P
   do while (nbuckets < 4 * nv_total)
      nbuckets = nbuckets * 2_I4P
   enddo
   mask = nbuckets - 1_I4P

   allocate(tab_key(0:nbuckets - 1), tab_head(0:nbuckets - 1), stat=istat)
   if (istat /= 0) return
   tab_key  = 0_I8P
   tab_head = 0_I4P

   ! Conservative initial capacity for the coord buffer; grow if exceeded.
   cap = max(nv_total / 2_I4P, 64_I4P)
   allocate(coord_buf(cap), next_in_chain(cap), self%facet_to_pool(3, nf), stat=istat)
   if (istat /= 0) then
      if (allocated(tab_key))  deallocate(tab_key)
      if (allocated(tab_head)) deallocate(tab_head)
      return
   endif
   next_in_chain = 0_I4P
   pid = 0_I4P

   do f = 1, nf
      do v = 1, 3
         vx = facet(f)%vertex(v)%x
         vy = facet(f)%vertex(v)%y
         vz = facet(f)%vertex(v)%z
         ix = floor((vx - bmin_x) * inv_cell, I4P)
         iy = floor((vy - bmin_y) * inv_cell, I4P)
         iz = floor((vz - bmin_z) * inv_cell, I4P)
         ! Look up in the 27-neighbour cluster.
         self%facet_to_pool(v, f) = find_or_insert(vx, vy, vz, ix, iy, iz, pid, cap, &
                                                   coord_buf, next_in_chain,        &
                                                   tab_key, tab_head, mask, istat)
         if (istat /= 0) then
            if (allocated(tab_key))       deallocate(tab_key)
            if (allocated(tab_head))      deallocate(tab_head)
            if (allocated(coord_buf))     deallocate(coord_buf)
            if (allocated(next_in_chain)) deallocate(next_in_chain)
            return
         endif
      enddo
   enddo

   self%n_vertices = pid
   allocate(self%coord_(self%n_vertices), stat=istat)
   if (istat /= 0) then
      deallocate(tab_key, tab_head, coord_buf, next_in_chain)
      return
   endif
   self%coord_(1:self%n_vertices) = coord_buf(1:self%n_vertices)

   deallocate(tab_key, tab_head, coord_buf, next_in_chain)
   endsubroutine build_via_spatial_hash

   function find_or_insert(vx, vy, vz, ix, iy, iz, pid, cap,    &
                           coord_buf, next_in_chain,            &
                           tab_key, tab_head, mask, istat) result(out_pid)
   !< Spatial-hash insertion helper. Searches the 27-neighbour cell cluster of
   !< (ix, iy, iz) for an EPS-coincident vertex; on hit returns that pool id, on
   !< miss allocates a new pool id (extending `coord_buf` and `next_in_chain` if
   !< necessary) and inserts it into the home bucket's chain.
   real(R8P),                     intent(in)    :: vx, vy, vz
   integer(I4P),                  intent(in)    :: ix, iy, iz
   integer(I4P),                  intent(inout) :: pid                  !< Running pool-id counter; advanced on insert.
   integer(I4P),                  intent(inout) :: cap                  !< Current capacity of growable buffers.
   type(vector_R8P), allocatable, intent(inout) :: coord_buf(:)
   integer(I4P),     allocatable, intent(inout) :: next_in_chain(:)
   integer(I8P),                  intent(inout) :: tab_key(0:)
   integer(I4P),                  intent(inout) :: tab_head(0:)
   integer(I4P),                  intent(in)    :: mask
   integer(I4P),                  intent(out)   :: istat
   integer(I4P)                                 :: out_pid
   integer(I4P)                                 :: dx, dy, dz, slot, probe, walker
   integer(I8P)                                 :: key, hkey
   type(vector_R8P), allocatable                :: tmp_coord(:)
   integer(I4P),     allocatable                :: tmp_next(:)

   istat = 0_I4P

   ! Probe all 27 cells.
   do dz = -1, 1
      do dy = -1, 1
         do dx = -1, 1
            key = pack_bucket_key(ix + dx, iy + dy, iz + dz)
            slot = bucket_index(key, mask)
            do probe = 0, mask
               if (tab_key(slot) == 0_I8P) exit  ! empty slot in open addressing -> bucket not in table
               if (tab_key(slot) == key) then
                  walker = tab_head(slot)
                  do while (walker /= 0_I4P)
                     if (eps_match_scalar(vx, vy, vz, coord_buf(walker))) then
                        out_pid = walker
                        return
                     endif
                     walker = next_in_chain(walker)
                  enddo
                  exit  ! this bucket walked; do not advance to next slot (different cell)
               endif
               slot = iand(slot + 1_I4P, mask)
            enddo
         enddo
      enddo
   enddo

   ! Miss: allocate a new pool id.
   pid = pid + 1_I4P
   if (pid > cap) then
      ! Grow buffers (doubling).
      allocate(tmp_coord(cap * 2_I4P), tmp_next(cap * 2_I4P), stat=istat)
      if (istat /= 0) then
         out_pid = 0_I4P
         return
      endif
      tmp_coord(1:cap) = coord_buf(1:cap)
      tmp_next (1:cap) = next_in_chain(1:cap)
      tmp_next(cap + 1 : cap * 2) = 0_I4P
      call move_alloc(from=tmp_coord, to=coord_buf)
      call move_alloc(from=tmp_next,  to=next_in_chain)
      cap = cap * 2_I4P
   endif
   coord_buf(pid) = vector_R8P(vx, vy, vz)
   out_pid = pid

   ! Insert into the home bucket's chain.
   key = pack_bucket_key(ix, iy, iz)
   slot = bucket_index(key, mask)
   do probe = 0, mask
      if (tab_key(slot) == 0_I8P) then
         tab_key(slot)        = key
         next_in_chain(pid)   = 0_I4P
         tab_head(slot)       = pid
         return
      endif
      if (tab_key(slot) == key) then
         next_in_chain(pid)   = tab_head(slot)
         tab_head(slot)       = pid
         return
      endif
      slot = iand(slot + 1_I4P, mask)
   enddo
   ! Table full -- this should never happen with the 4x sizing.
   istat = 2_I4P
   endfunction find_or_insert

   pure function pack_bucket_key(ix, iy, iz) result(key)
   !< Pack three 21-bit signed integers (range ~+/-1 million) into a single I8P.
   !< Bias each by 2^20 so negatives encode positively, then pack.
   integer(I4P), intent(in) :: ix, iy, iz
   integer(I8P)             :: key
   integer(I8P), parameter  :: BIAS = 1048576_I8P   ! 2^20
   integer(I8P), parameter  :: MASK21 = 2097151_I8P ! 2^21 - 1
   key = iand(int(ix, I8P) + BIAS, MASK21)              + &
         ishft(iand(int(iy, I8P) + BIAS, MASK21), 21)   + &
         ishft(iand(int(iz, I8P) + BIAS, MASK21), 42)
   ! Reserve key=0 as "empty slot" by offsetting by 1.
   key = key + 1_I8P
   endfunction pack_bucket_key

   pure function bucket_index(key, mask) result(slot)
   !< Mix the 63-bit packed bucket key into a hash and mask to a slot index.
   !< Multiplier from Knuth's MMIX (a strong 64-bit mixer for arbitrary keys).
   integer(I8P), intent(in) :: key
   integer(I4P), intent(in) :: mask
   integer(I4P)             :: slot
   integer(I8P), parameter  :: MULT = 6364136223846793005_I8P
   slot = int(iand(ishft(key * MULT, -32), int(mask, I8P)), I4P)
   endfunction bucket_index

   pure function eps_match_scalar(vx, vy, vz, b) result(yes)
   !< Strict EPS coincidence between scalar (vx, vy, vz) and an existing pool coord.
   real(R8P),        intent(in) :: vx, vy, vz
   type(vector_R8P), intent(in) :: b
   logical                      :: yes
   yes = (abs(vx - b%x) <= EPS) .and. &
         (abs(vy - b%y) <= EPS) .and. &
         (abs(vz - b%z) <= EPS)
   endfunction eps_match_scalar

   subroutine build_inverted_index(self, istat)
   !< CSR inverted-index build: pool_id -> all (facet, local_v) referencing it.
   !<
   !< Three-pass O(N): (1) count refs per pool id, (2) prefix-sum into at_offset,
   !< (3) scatter pairs. Total payload size is exactly 3 * n_facets (each facet has
   !< three vertex slots, each contributing one pair).
   type(vertex_pool_object), intent(inout) :: self
   integer(I4P),             intent(out)   :: istat
   integer(I4P)                            :: f, v, vid, k
   integer(I4P), allocatable               :: cursor(:)

   istat = 0_I4P
   if (self%n_vertices <= 0) return

   allocate(self%at_offset(self%n_vertices + 1), stat=istat)
   if (istat /= 0) return
   self%at_offset = 0_I4P

   ! Pass 1: count.
   do f = 1, self%n_facets
      do v = 1, 3
         vid = self%facet_to_pool(v, f)
         self%at_offset(vid + 1) = self%at_offset(vid + 1) + 1_I4P
      enddo
   enddo

   ! Pass 2: prefix sum -> offsets. at_offset(1) = 1 (1-based), at_offset(end) = total + 1.
   self%at_offset(1) = 1_I4P
   do k = 2, self%n_vertices + 1
      self%at_offset(k) = self%at_offset(k - 1) + self%at_offset(k)
   enddo

   ! Pass 3: scatter pairs.
   allocate(self%at_pairs(2, 3 * self%n_facets), cursor(self%n_vertices), stat=istat)
   if (istat /= 0) return
   cursor = self%at_offset(1:self%n_vertices)
   do f = 1, self%n_facets
      do v = 1, 3
         vid = self%facet_to_pool(v, f)
         k = cursor(vid)
         self%at_pairs(1, k) = f
         self%at_pairs(2, k) = v
         cursor(vid) = k + 1_I4P
      enddo
   enddo
   deallocate(cursor)
   endsubroutine build_inverted_index

   pure function vertex_count(self) result(n)
   !< Return number of unique vertices in the pool.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P)                          :: n
   n = self%n_vertices
   endfunction vertex_count

   pure function facets_count(self) result(n)
   !< Return the source facet count the pool was built from.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P)                          :: n
   n = self%n_facets
   endfunction facets_count

   pure function coord(self, vid) result(v)
   !< Return the coordinate of pool vertex `vid`. Returns origin if out of range.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P),              intent(in) :: vid
   type(vector_R8P)                      :: v
   if (vid >= 1 .and. vid <= self%n_vertices) then
      v = self%coord_(vid)
   else
      v = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   endif
   endfunction coord

   pure function facet_vid(self, facet_id, local_v) result(vid)
   !< Return the pool vertex id for (facet_id, local_v). Returns 0 if out of range.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P),              intent(in) :: facet_id
   integer(I4P),              intent(in) :: local_v
   integer(I4P)                          :: vid
   if (facet_id >= 1 .and. facet_id <= self%n_facets .and. &
       local_v  >= 1 .and. local_v  <= 3) then
      vid = self%facet_to_pool(local_v, facet_id)
   else
      vid = 0_I4P
   endif
   endfunction facet_vid

   pure function facets_at_count(self, vid) result(n)
   !< Number of (facet, local_v) pairs referencing pool vertex `vid`. 0 if out of range.
   class(vertex_pool_object), intent(in) :: self
   integer(I4P),              intent(in) :: vid
   integer(I4P)                          :: n
   if (vid < 1 .or. vid > self%n_vertices) then
      n = 0_I4P
   else if (.not. allocated(self%at_offset)) then
      n = 0_I4P
   else
      n = self%at_offset(vid + 1) - self%at_offset(vid)
   endif
   endfunction facets_at_count

   pure subroutine facets_at(self, vid, k, facet_id, local_v)
   !< Return the `k`-th (1-based) (facet_id, local_v) pair referencing pool vertex `vid`.
   !< Out-of-range queries set facet_id = 0, local_v = 0.
   class(vertex_pool_object), intent(in)  :: self
   integer(I4P),              intent(in)  :: vid
   integer(I4P),              intent(in)  :: k
   integer(I4P),              intent(out) :: facet_id, local_v
   integer(I4P)                           :: row

   facet_id = 0_I4P
   local_v  = 0_I4P
   if (vid < 1 .or. vid > self%n_vertices) return
   if (.not. allocated(self%at_offset)) return
   if (k < 1 .or. k > self%at_offset(vid + 1) - self%at_offset(vid)) return
   row = self%at_offset(vid) + k - 1_I4P
   facet_id = self%at_pairs(1, row)
   local_v  = self%at_pairs(2, row)
   endsubroutine facets_at

   pure function get_is_initialized(self) result(yes)
   !< Has the pool been built?
   class(vertex_pool_object), intent(in) :: self
   logical                               :: yes
   yes = self%is_initialized
   endfunction get_is_initialized

   pure subroutine destroy(self)
   !< Release all storage and reset to uninitialized.
   class(vertex_pool_object), intent(inout) :: self
   if (allocated(self%coord_))        deallocate(self%coord_)
   if (allocated(self%facet_to_pool)) deallocate(self%facet_to_pool)
   if (allocated(self%at_offset))     deallocate(self%at_offset)
   if (allocated(self%at_pairs))      deallocate(self%at_pairs)
   self%n_vertices     = 0_I4P
   self%n_facets       = 0_I4P
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine finalize(self)
   !< Finaliser. `class` not allowed on `final`; takes `type`.
   type(vertex_pool_object), intent(inout) :: self
   if (allocated(self%coord_))        deallocate(self%coord_)
   if (allocated(self%facet_to_pool)) deallocate(self%facet_to_pool)
   if (allocated(self%at_offset))     deallocate(self%at_offset)
   if (allocated(self%at_pairs))      deallocate(self%at_pairs)
   endsubroutine finalize

   ! ---- internal helpers ----

   pure function vertices_eps_match(a, b) result(yes)
   !< Strict EPS coincidence -- same rule as `vertices_match` inside
   !< compute_connectivity. Kept private to the pool so the coincidence
   !< definition lives in one place.
   type(vector_R8P), intent(in) :: a, b
   logical                      :: yes
   yes = (abs(a%x - b%x) <= EPS) .and. &
         (abs(a%y - b%y) <= EPS) .and. &
         (abs(a%z - b%z) <= EPS)
   endfunction vertices_eps_match

   pure function uf_find(x, parent) result(root)
   !< Union-find root with path halving. Pure: parent is intent(in), so the
   !< usual path-compression write is replaced by re-traversal -- O(log* N)
   !< amortized in practice for our densities, acceptable in stage 1.
   integer(I4P), intent(in) :: x
   integer(I4P), intent(in) :: parent(:)
   integer(I4P)             :: root
   root = x
   do while (parent(root) /= root)
      root = parent(root)
   enddo
   endfunction uf_find

   subroutine uf_union(a, b, parent, rank_)
   !< Union-by-rank on `parent(:)`. Path compression on lookup.
   integer(I4P), intent(in)    :: a, b
   integer(I4P), intent(inout) :: parent(:), rank_(:)
   integer(I4P)                :: ra, rb

   ra = uf_find(a, parent)
   rb = uf_find(b, parent)
   if (ra == rb) return
   if (rank_(ra) < rank_(rb)) then
      parent(ra) = rb
   else if (rank_(ra) > rank_(rb)) then
      parent(rb) = ra
   else
      parent(rb) = ra
      rank_(ra)  = rank_(ra) + 1_I4P
   endif
   endsubroutine uf_union

endmodule fossil_vertex_pool_object
