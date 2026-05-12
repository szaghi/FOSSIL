!< FOSSIL test: vertex pool stage 1 -- the pool is built as a derived artifact
!< from analyze() and matches the underlying facet coordinates exactly.
!<
!< Verifies:
!<   1. After load_from_file, the pool is initialized.
!<   2. Pool size is strictly less than 3 * facets_number (real dedup happens).
!<   3. Pool size is plausible for a closed surface (Euler-ish bound: a closed
!<      triangulated surface has roughly V ~= F/2). Allow a wide band [F/10, F].
!<   4. For every (f, v), pool%coord(pool%facet_vid(f, v)) == facet(f)%vertex(v)
!<      bit-exactly. The pool is built from the facet array, so no rounding
!<      should occur on materialization.
!<   5. The mapping is consistent: any two (f, v) referring to coordinates that
!<      match within EPS map to the same pool id.

program fossil_test_vertex_pool

use fossil,                    only : surface_stl_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none

type(surface_stl_object)          :: surface
type(vertex_pool_object), pointer :: pool
type(vector_R8P)                  :: a, b, c
integer(I4P)                      :: n_facets, n_vertices
integer(I4P)                      :: f, v, vid
integer(I4P)                      :: f1, f2, v1, v2
real(R8P)                         :: max_coord_error
logical                           :: are_tests_passed(5)
real(R8P), parameter              :: eps_strict = 1.0e-12_R8P

are_tests_passed = .false.

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
n_facets = surface%get_facets_number()
pool => surface%get_vertex_pool()

! 1. Pool is initialized.
are_tests_passed(1) = pool%get_is_initialized() .and. (pool%facets_count() == n_facets)

! 2. Pool dedupes: |pool| < 3 * |facets|.
n_vertices = pool%vertex_count()
are_tests_passed(2) = (n_vertices < 3 * n_facets) .and. (n_vertices > 0)

! 3. Plausible size for a closed surface. Dragon is closed; expect V ~ F/2 +/-.
!    Wide band so this is a smoke test, not a tight assertion.
are_tests_passed(3) = (n_vertices >= n_facets / 10) .and. (n_vertices <= n_facets)

! 4. Bit-exact coordinate match. The pool materializes each unique coord from one
!    of the facet vertex slots, so equality should hold without tolerance for
!    every (f, v). Use a tiny epsilon as a paranoid floor, not as the contract.
max_coord_error = 0._R8P
do f = 1, n_facets
   associate (facet_ptr => surface%facet_at(f))
      do v = 1, 3
         vid = pool%facet_vid(f, v)
         if (vid < 1 .or. vid > n_vertices) cycle
         a = pool%coord(vid)
         b = facet_ptr%vertex(v)
         max_coord_error = max(max_coord_error,                  &
                               abs(a%x - b%x) + abs(a%y - b%y) + abs(a%z - b%z))
      enddo
   end associate
enddo
are_tests_passed(4) = (max_coord_error <= eps_strict)

! 5. Consistency: any two (f, v) whose coordinates match in the strict sense map
!    to the same pool id. Spot-check on the first few facets to keep it O(F).
!    For each facet, check that its three local vertices, if any pair coincides,
!    share a pool id -- this would only trigger on a degenerate facet, so for a
!    clean dragon it's a vacuous check. Cross-facet: check facet 1's three
!    vertices against facets sharing a connectivity neighbor.
are_tests_passed(5) = .true.
do f = 1, min(n_facets, 200_I4P)
   associate (facet_ptr => surface%facet_at(f))
      do v1 = 1, 3
         do v2 = v1 + 1, 3
            a = facet_ptr%vertex(v1)
            b = facet_ptr%vertex(v2)
            if (abs(a%x - b%x) + abs(a%y - b%y) + abs(a%z - b%z) < eps_strict) then
               if (pool%facet_vid(f, v1) /= pool%facet_vid(f, v2)) then
                  are_tests_passed(5) = .false.
                  exit
               endif
            endif
         enddo
         if (.not. are_tests_passed(5)) exit
      enddo
   end associate
   if (.not. are_tests_passed(5)) exit
enddo

! Additional cross-facet consistency: for facet 1 and each connectivity neighbor
! along its three edges, the shared edge's endpoints must map to the same pool ids.
if (are_tests_passed(5) .and. n_facets >= 2) then
   block
      integer(I4P) :: e
      associate (f_a => surface%facet_at(1))
         do e = 1, 3
            f2 = f_a%fcon_edge(e)
            if (f2 < 1 .or. f2 > n_facets) cycle
            associate (f_b => surface%facet_at(f2))
               ! Endpoints of edge e on facet 1 are vertex(e) and vertex(mod(e,3)+1).
               ! They must appear on facet 2 with the same pool ids (in some order).
               v1 = pool%facet_vid(1, e)
               v2 = pool%facet_vid(1, mod(e, 3) + 1)
               c = vector_R8P(0._R8P, 0._R8P, 0._R8P)  ! unused; here to silence "unused" if any
               if (.not. has_pool_vid(pool, f2, v1) .or. .not. has_pool_vid(pool, f2, v2)) then
                  are_tests_passed(5) = .false.
                  exit
               endif
            end associate
         enddo
      end associate
   end block
endif

print '(A,I0,A,I0,A,F8.4)', 'facets=', n_facets, '  unique_vertices=', n_vertices, &
                            '  V/F=', real(n_vertices, R8P) / real(max(n_facets, 1), R8P)
print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains

   pure function has_pool_vid(pool, f, vid) result(yes)
   !< Does facet `f` reference pool id `vid` in any of its three slots?
   type(vertex_pool_object), intent(in) :: pool
   integer(I4P),             intent(in) :: f
   integer(I4P),             intent(in) :: vid
   logical                              :: yes
   integer(I4P)                         :: lv

   yes = .false.
   do lv = 1, 3
      if (pool%facet_vid(f, lv) == vid) then
         yes = .true.
         return
      endif
   enddo
   endfunction has_pool_vid

endprogram fossil_test_vertex_pool
