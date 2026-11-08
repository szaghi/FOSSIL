!< FOSSIL test: read-only access to the acceleration trees (issue #89).
!<
!< `surface%tree_boxes` on cube-01 and the dragon, both tree kinds, and `curve%tree_boxes` on a 64-gon:
!<   - the root first (depth 0, parent 0), every other box after its parent, one level deeper, inside its parent's box;
!<   - the leaves' facets (segments) are every facet once; a BVH leaf box contains its facets' vertices, an octree leaf
!<     its facets' centroids (a cell of the space partition); the BVH root box is the surface box;
!<   - `max_depth = k` returns the prefix of the whole tree's list with depth <= k; max_depth < 0 is invalid input;
!<   - `save_aabb_into_file` writes one STL per non-empty leaf, for the BVH too (it wrote nothing before #89).
program fossil_test_tree_boxes
use fossil, only : surface_stl_object, curve_object, AABB_TREE_OCTREE, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object) :: cube, cube_oct, dragon, dragon_oct, empty
type(curve_object)       :: curve
type(vector_R8P), allocatable :: bmin(:), bmax(:), lo_v(:), hi_v(:)
real(R8P),    allocatable :: u(:), v(:), lo(:,:), hi(:,:)
integer(I4P), allocatable :: depth(:), parent(:), start(:), ids(:)
logical,      allocatable :: leaf(:)
integer(I4P)              :: st, k, nfiles
logical                   :: passed(6), there
character(64)             :: name

passed = .false.
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call cube_oct%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon_oct%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
passed(1) = tree_ok(cube, 'cube-01 BVH', .true.) .and. tree_ok(cube_oct, 'cube-01 octree', .false.)
passed(2) = tree_ok(dragon, 'dragon BVH', .true.) .and. tree_ok(dragon_oct, 'dragon octree', .false.)

! invalid depth, empty surface
call dragon%tree_boxes(bmin=bmin, bmax=bmax, max_depth=-1_I4P, status=st)
passed(3) = st == STATUS_INVALID_INPUT .and. size(bmin) == 0
call empty%tree_boxes(bmin=bmin, bmax=bmax, status=st)
passed(3) = passed(3) .and. st == STATUS_OK .and. size(bmin) == 0
print '(A,L1)', 'max_depth = -1 invalid, empty surface no boxes: ', passed(3)

! curve twin: a 64-gon
allocate(u(64), v(64))
do k = 1, 64
   u(k) = cos(6.283185307179586_R8P * (k - 1) / 64) ; v(k) = sin(6.283185307179586_R8P * (k - 1) / 64)
enddo
call curve%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
call curve%tree_boxes(bmin=lo, bmax=hi, depth=depth, leaf=leaf, parent=parent, segment_start=start, segment_ids=ids, &
                      status=st)
passed(4) = st == STATUS_OK .and. size(depth) > 1 .and. parent(1) == 0 .and. depth(1) == 0
do k = 2, size(depth)
   passed(4) = passed(4) .and. parent(k) < k .and. depth(k) == depth(parent(k)) + 1 .and. &
               all(lo(:, k) >= lo(:, parent(k))) .and. all(hi(:, k) <= hi(:, parent(k)))
enddo
passed(4) = passed(4) .and. once(ids, 64_I4P) .and. count(leaf) == count(start(2:) > start(:size(start) - 1))
do k = 1, size(depth)
   if (leaf(k)) passed(4) = passed(4) .and. all(u(ids(start(k):start(k+1)-1)) >= lo(1, k)) .and. &
                            all(u(ids(start(k):start(k+1)-1)) <= hi(1, k))
enddo
print '(A,I0,A,I0,A,L1)', 'curve: boxes ', size(depth), ', leaves ', count(leaf), ', ok ', passed(4)
call curve%tree_boxes(bmin=lo, bmax=hi, max_depth=2_I4P, depth=depth)
passed(4) = passed(4) .and. maxval(depth) == 2 .and. size(depth) == 7

! save_aabb_into_file on the BVH
call cube%save_aabb_into_file(base_file_name='fossil_test_tree_boxes_', is_ascii=.true.)
call cube%tree_boxes(bmin=bmin, bmax=bmax, depth=depth, leaf=leaf)
nfiles = 0
do k = 1, size(bmin)
   if (.not. leaf(k)) cycle
   write(name, '(A,I0,A,I0,A)') 'fossil_test_tree_boxes_aabb-l_', depth(k), '-b_', k, '.stl'
   inquire(file=trim(name), exist=there)
   if (there) then
      nfiles = nfiles + 1
      open(newunit=st, file=trim(name)) ; close(st, status='delete')
   endif
enddo
passed(5) = nfiles == count(leaf) .and. nfiles > 0
print '(A,I0,A,I0)', 'save_aabb_into_file (BVH): leaves ', count(leaf), ', files ', nfiles

! a depth cap is the prefix of the whole list
call dragon%tree_boxes(bmin=bmin, bmax=bmax, depth=depth)
k = count(depth <= 4)
call dragon%tree_boxes(bmin=lo_v, bmax=hi_v, max_depth=4_I4P, depth=parent)
passed(6) = size(lo_v) == k .and. all(parent == depth(1:k)) .and. all(lo_v%x == bmin(1:k)%x) .and. &
            all(hi_v%z == bmax(1:k)%z)
print '(A,I0,A,L1)', 'dragon BVH, max_depth 4: boxes ', k, ', prefix of the whole list ', passed(6)

print '(A,*(L2))', 'per-case results:', passed
print '(A,L1)', 'Are all tests passed? ', all(passed)
if (.not. all(passed)) error stop 1

contains
   function tree_ok(s, label, is_bvh) result(ok)
   !< Order, nesting and payload of the whole tree of s.
   type(surface_stl_object), intent(in), target :: s      !< Surface.
   character(*),             intent(in)         :: label  !< Label.
   logical,                  intent(in)         :: is_bvh !< BVH (else octree).
   logical                                      :: ok     !< Result.
   type(vector_R8P), allocatable                :: b0(:), b1(:) !< Boxes.
   integer(I4P),     allocatable                :: dp(:), pa(:), fs(:), fi(:) !< Depths, parents, leaf facets.
   logical,          allocatable                :: lf(:)  !< Leaves.
   type(facet_object), pointer                  :: f      !< Facet.
   type(vector_R8P)                             :: p      !< Point.
   integer(I4P)                                 :: i, j, m, st !< Counters.
   logical                                      :: nested, cover, inside

   call s%tree_boxes(bmin=b0, bmax=b1, depth=dp, leaf=lf, parent=pa, facet_start=fs, facet_ids=fi, status=st)
   ok = st == STATUS_OK .and. size(b0) > 0
   if (.not. ok) return
   nested = pa(1) == 0 .and. dp(1) == 0
   do i = 2, size(b0)
      j = pa(i)
      nested = nested .and. j >= 1 .and. j < i .and. dp(i) == dp(j) + 1 .and. inbox(b0(i), b0(j), b1(j)) .and. &
               inbox(b1(i), b0(j), b1(j))
   enddo
   cover = once(fi, s%get_facets_number())
   inside = .true.
   do i = 1, size(b0)
      if (.not. lf(i)) then
         inside = inside .and. fs(i + 1) == fs(i)
         cycle
      endif
      do m = fs(i), fs(i + 1) - 1
         f => s%facet_at(fi(m))
         if (is_bvh) then
            do j = 1, 3
               inside = inside .and. inbox(f%vertex(j), b0(i), b1(i))
            enddo
         else
            p = f%centroid
            inside = inside .and. inbox(p, b0(i), b1(i))
         endif
      enddo
   enddo
   if (is_bvh) then
      p = s%get_bmin() ; inside = inside .and. p%x == b0(1)%x .and. p%y == b0(1)%y .and. p%z == b0(1)%z
      p = s%get_bmax() ; inside = inside .and. p%x == b1(1)%x .and. p%y == b1(1)%y .and. p%z == b1(1)%z
   endif
   ok = nested .and. cover .and. inside
   print '(A,A,I0,A,I0,A,I0,3(A,L1))', label, ': boxes ', size(b0), ', leaves ', count(lf), ', depth ', maxval(dp), &
                                       ', nested ', nested, ', facets once ', cover, ', in their leaf ', inside
   endfunction tree_ok

   pure function inbox(p, lo, hi) result(yes)
   !< p in the closed box [lo, hi].
   type(vector_R8P), intent(in) :: p, lo, hi !< Point, box.
   logical                      :: yes       !< Result.

   yes = p%x >= lo%x .and. p%x <= hi%x .and. p%y >= lo%y .and. p%y <= hi%y .and. p%z >= lo%z .and. p%z <= hi%z
   endfunction inbox

   pure function once(ids, n) result(yes)
   !< ids is a permutation of 1..n.
   integer(I4P), intent(in) :: ids(:) !< Ids.
   integer(I4P), intent(in) :: n      !< Count.
   logical                  :: yes    !< Result.
   integer(I4P)             :: seen(n), i

   yes = size(ids) == n
   if (.not. yes) return
   seen = 0
   do i = 1, n
      if (ids(i) < 1 .or. ids(i) > n) then
         yes = .false. ; return
      endif
      seen(ids(i)) = seen(ids(i)) + 1
   enddo
   yes = all(seen == 1)
   endfunction once
endprogram fossil_test_tree_boxes
