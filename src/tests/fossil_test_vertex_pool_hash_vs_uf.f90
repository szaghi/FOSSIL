!< FOSSIL test: spatial-hash pool builder produces the same partition as the
!< O((3N)^2) union-find reference (issue #5 stage 5).
!<
!< Two pools are built from the dragon facet array: one via the spatial hash
!< (default), one via the union-find reference (use_union_find=.true.). Numeric
!< pool ids may differ, but the partition of (facet, local_v) slots into
!< equivalence classes MUST be identical.
!<
!< Partition equality check (O(N)):
!<   For each path, build a canonical-form vector by relabelling every id with
!<   the smallest global index (facet*3 + local_v) that shares it. If the two
!<   canonical vectors agree element-wise, the partitions are identical.

program fossil_test_vertex_pool_hash_vs_uf

use fossil,                    only : surface_stl_object
use fossil_facet_object,       only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P

implicit none

type(surface_stl_object)        :: surface
type(vertex_pool_object)        :: pool_hash, pool_uf
type(facet_object), allocatable :: facet(:)
integer(I4P)                    :: nf, n_slots, f, v
integer(I4P), allocatable       :: canon_hash(:), canon_uf(:)
logical                         :: are_tests_passed(3)

are_tests_passed = .false.

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
nf = surface%get_facets_number()
n_slots = 3 * nf

! Snapshot the facet array (the surface owns it as `private`, so we go through
! the public `facet_at` accessor and copy).
allocate(facet(nf))
do f = 1, nf
   associate (fp => surface%facet_at(f))
      facet(f) = fp
   end associate
enddo

call pool_hash%initialize_from_facets(facet=facet)                       ! default: hash
call pool_uf  %initialize_from_facets(facet=facet, use_union_find=.true.)

! Vertex counts must agree exactly -- same EPS coincidence rule, same input.
are_tests_passed(1) = (pool_hash%vertex_count() == pool_uf%vertex_count())

! Build canonical-form vectors. canon(k) = min { k' : pool_path%facet_vid(f,v) == pool_path%facet_vid(f',v') }
! where k = (f-1)*3 + v. Two paths producing the same partition will yield
! identical canon(:) arrays.
allocate(canon_hash(n_slots), canon_uf(n_slots))
call canonicalize(pool_hash, nf, n_slots, canon_hash)
call canonicalize(pool_uf,   nf, n_slots, canon_uf)

are_tests_passed(2) = all(canon_hash == canon_uf)

! Sanity: both partitions must have exactly pool_hash%vertex_count() distinct classes.
block
   logical, allocatable :: seen(:)
   integer(I4P)         :: k, n_classes
   allocate(seen(n_slots))
   seen = .false.
   n_classes = 0
   do k = 1, n_slots
      if (.not. seen(canon_hash(k))) then
         seen(canon_hash(k)) = .true.
         n_classes = n_classes + 1
      endif
   enddo
   are_tests_passed(3) = (n_classes == pool_hash%vertex_count())
   deallocate(seen)
end block

print '(A,I0,A,I0)',     'pool_hash unique vertices=', pool_hash%vertex_count(), &
                         '  pool_uf unique vertices=', pool_uf%vertex_count()
print '(A,3L2)',         'per-case results: ', are_tests_passed
print '(A,L1)',          'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains

   subroutine canonicalize(pool, nf, n_slots, canon)
   !< Build a canonical-form array: canon(k) = minimum slot index that maps to
   !< the same pool id as slot k. Two pools producing equivalent partitions
   !< yield equal canon(:) arrays.
   type(vertex_pool_object), intent(in)  :: pool
   integer(I4P),             intent(in)  :: nf, n_slots
   integer(I4P),             intent(out) :: canon(:)
   integer(I4P), allocatable             :: first_slot_for_pid(:)
   integer(I4P)                          :: f, v, k, pid

   allocate(first_slot_for_pid(pool%vertex_count()))
   first_slot_for_pid = 0
   k = 0
   do f = 1, nf
      do v = 1, 3
         k = k + 1
         pid = pool%facet_vid(f, v)
         if (first_slot_for_pid(pid) == 0) first_slot_for_pid(pid) = k
         canon(k) = first_slot_for_pid(pid)
      enddo
   enddo
   deallocate(first_slot_for_pid)
   endsubroutine canonicalize

endprogram fossil_test_vertex_pool_hash_vs_uf
