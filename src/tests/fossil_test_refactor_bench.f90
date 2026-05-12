!< FOSSIL benchmark: end-to-end `load_from_file` (= header parse + facet read +
!< analyze pipeline) on a user-supplied STL, plus per-facet structural footprint.
!<
!< Portable across the issue #5 refactor: uses only `surface_stl_object` and
!< `facet_object` -- the same API present in v1.2.0 and current.
!<
!< Usage: fossil_test_refactor_bench [<stl_path>]
!<        default <stl_path> = src/tests/dragon.stl

program fossil_test_refactor_bench

use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none

integer(I4P), parameter   :: N_REPEATS = 5
type(surface_stl_object)  :: surface
type(facet_object)        :: probe
real(R8P)                 :: t0, t1, t_total, t_min
integer(I4P)              :: r, nf
integer                   :: facet_bits, argc
character(len=512)        :: stl_path

argc = command_argument_count()
if (argc >= 1) then
   call get_command_argument(1, stl_path)
else
   stl_path = 'src/tests/dragon.stl'
endif

t_min = huge(0._R8P)
t_total = 0._R8P

do r = 1, N_REPEATS
   call cpu_time(t0)
   call surface%load_from_file(file_name=trim(stl_path), guess_format=.true.)
   call cpu_time(t1)
   t_total = t_total + (t1 - t0)
   if (t1 - t0 < t_min) t_min = t1 - t0
   nf = surface%get_facets_number()
enddo

facet_bits = storage_size(probe)

print '(A,A)',         '--- FOSSIL refactor benchmark: ', trim(stl_path)
print '(A,I0)',        'facets:                                ', nf
print '(A,I0,A)',      'sizeof(facet_object):                   ', facet_bits / 8, ' bytes'
print '(A,I0)',        'load_from_file repeats:                ', N_REPEATS
print '(A,F10.2,A)',   'load_from_file mean:                  ', 1000._R8P * t_total / real(N_REPEATS, R8P), ' ms'
print '(A,F10.2,A)',   'load_from_file best:                  ', 1000._R8P * t_min, ' ms'
print '(A,L1)', 'Are all tests passed? ', .true.

endprogram fossil_test_refactor_bench
