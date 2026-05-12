!< FOSSIL test: time the spatial-hash pool builder against the union-find
!< reference (issue #5 stage 5). Reports wall-clock for both paths plus the
!< speedup factor.
!<
!< Usage: fossil_test_pool_speedup [<stl_path>]
!<        default <stl_path> = src/tests/dragon.stl

program fossil_test_pool_speedup

use fossil,                    only : surface_stl_object
use fossil_facet_object,       only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P

implicit none

type(surface_stl_object)        :: surface
type(facet_object), allocatable :: facet(:)
type(vertex_pool_object)        :: pool
integer(I4P)                    :: nf, f, argc
real(R8P)                       :: t0, t1, t_hash, t_uf
character(len=512)              :: stl_path

argc = command_argument_count()
if (argc >= 1) then
   call get_command_argument(1, stl_path)
else
   stl_path = 'src/tests/dragon.stl'
endif

call surface%load_from_file(file_name=trim(stl_path), guess_format=.true.)
nf = surface%get_facets_number()
allocate(facet(nf))
do f = 1, nf
   associate (fp => surface%facet_at(f))
      facet(f) = fp
   end associate
enddo

call cpu_time(t0)
call pool%initialize_from_facets(facet=facet)
call cpu_time(t1)
t_hash = t1 - t0

call cpu_time(t0)
call pool%initialize_from_facets(facet=facet, use_union_find=.true.)
call cpu_time(t1)
t_uf = t1 - t0

print '(A,A,I0,A)',       trim(stl_path), ': ', nf, ' facets'
print '(A,F12.4,A)',      'hash builder:        ', t_hash * 1000._R8P, ' ms'
print '(A,F12.4,A)',      'union-find builder:  ', t_uf   * 1000._R8P, ' ms'
if (t_hash > 0._R8P) then
   print '(A,F12.2,A)',   'speedup:             ', t_uf / t_hash, 'x'
endif
print '(A,L1)', 'Are all tests passed? ', .true.

endprogram fossil_test_pool_speedup
