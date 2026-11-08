!run ref_cut_cells_lines ref_cut_cells_lines
!< Reference page "Cut-cell geometry": batched grid-line crossings. 24 x 24 lines along x through the dragon's box,
!< into fixed buffers of 16 crossings per line (surface%line_crossings_many), checked line by line against
!< surface%line_crossings, and the entering and leaving crossings of every line counted.
program ref_cut_cells_lines
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
integer(I4P), parameter   :: N = 24, NMAX = 16
type(surface_stl_object)  :: dragon
type(vector_R8P)          :: lo, hi
real(R8P)                 :: a(N * N), b(N * N), t(NMAX, N * N)
integer(I4P)              :: facets(NMAX, N * N), counts(N * N), i, j, m, nodd, ndiff
logical                   :: entering(NMAX, N * N), overflow(N * N)
real(R8P),    allocatable :: t1(:)
integer(I4P), allocatable :: f1(:)
logical,      allocatable :: e1(:)

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
lo = dragon%get_bmin() ; hi = dragon%get_bmax()
do j = 1, N
   do i = 1, N                                         ! line (i, j) at y = a, z = b, cell-centred in the box
      a(i + (j - 1) * N) = lo%y + (hi%y - lo%y) * (i - 0.5_R8P) / N
      b(i + (j - 1) * N) = lo%z + (hi%z - lo%z) * (j - 0.5_R8P) / N
   enddo
enddo
call dragon%line_crossings_many(axis=1_I4P, a=a, b=b, t0=lo%x - 1._R8P, t1=hi%x + 1._R8P, t=t, facets=facets, &
                                entering=entering, counts=counts, overflow=overflow)
print '(A,I0,A,I0,A,I0,A,I0)', 'lines ', N * N, ', crossing the dragon ', count(counts > 0), ', crossings ', sum(counts), &
                               ', overflowed ', count(overflow)
nodd = 0 ; ndiff = 0
do m = 1, N * N
   if (count(entering(1:counts(m), m)) /= counts(m) - count(entering(1:counts(m), m))) nodd = nodd + 1
   call dragon%line_crossings(axis=1_I4P, a=a(m), b=b(m), t0=lo%x - 1._R8P, t1=hi%x + 1._R8P, t=t1, facets=f1, &
                              entering=e1, n=i)
   if (i /= counts(m)) then
      ndiff = ndiff + 1
   elseif (any(t1(1:i) /= t(1:i, m)) .or. any(f1(1:i) /= facets(1:i, m))) then
      ndiff = ndiff + 1
   endif
enddo
print '(A,L1)', 'dragon watertight ', dragon%is_watertight()
print '(A,I0,A,I0)', 'lines with unequal entering/leaving counts ', nodd, ', lines differing from line_crossings ', ndiff
m = maxloc(counts, 1)
print '(A,I0,A,F7.3,A,F7.3,A,I0,A)', 'most crossed line: ', counts(m), ' crossings, sorted from x = ', t(1, m), ' to ', &
                                    t(counts(m), m), ', ', count(entering(1:counts(m), m)), ' entering'
endprogram ref_cut_cells_lines
