!< FOSSIL test support: bitwise comparison of two piece maps (issue #105 step 1a).
module fossil_piece_map_compare
!< FOSSIL test support: bitwise comparison of two piece maps (issue #105 step 1a).
!<
!< `piece_map_2d_mismatches(a, b)` and `piece_map_3d_mismatches(a, b)` count the differing elements of every field of
!< `piece_map_2d_t` / `piece_map_3d_t`: the scalars (nx, ny, nz, nt, np0, np1, exact, untrusted_cells,
!< multi_track_cells, evaluated_cells, evaluations, skipped_nodes, contact_cells, narrow_cells) and the counters
!< n_status(-1:7) one by one, every integer array elementwise (tracks, first_track, track_cell, status, old_start,
!< old_piece, new_start, new_piece, old_track, new_track) and every real array by its bits (`transfer` to a 64-bit
!< integer, so a NaN equals itself and -0 differs from +0: volume0, volume1, faces, wall, residual, outflow). An array
!< allocated in one map only, or of another shape, counts as max(1, its size). With `verbose=.true.` every field that
!< differs is printed with its count. Used by `fossil_test_piece_map_compare` (CI) and by the later steps of #105 (the
!< refactored host build and the device maps against the reference build).
use fossil, only : piece_map_2d_t, piece_map_3d_t
use penf,   only : I4P, I8P, R8P

implicit none
private
public :: piece_map_2d_mismatches, piece_map_3d_mismatches

interface cmp
   module procedure cmp_i0, cmp_l0, cmp_i1, cmp_i2, cmp_i3, cmp_r1, cmp_r2
endinterface cmp

contains
   function piece_map_2d_mismatches(a, b, verbose) result(n)
   !< Differing elements of two 2-D piece maps (see the module notes).
   type(piece_map_2d_t), intent(in)           :: a, b    !< Maps.
   logical,              intent(in), optional :: verbose !< Print every differing field.
   integer(I4P)                               :: n       !< Mismatches.
   logical                                    :: v       !< Verbose.
   integer(I4P)                               :: k       !< Counter.

   v = .false. ; if (present(verbose)) v = verbose
   n = 0_I4P
   call cmp(n, v, 'nx', a%nx, b%nx)
   call cmp(n, v, 'ny', a%ny, b%ny)
   call cmp(n, v, 'nt', a%nt, b%nt)
   call cmp(n, v, 'np0', a%np0, b%np0)
   call cmp(n, v, 'np1', a%np1, b%np1)
   call cmp(n, v, 'tracks', a%tracks, b%tracks)
   call cmp(n, v, 'first_track', a%first_track, b%first_track)
   call cmp(n, v, 'track_cell', a%track_cell, b%track_cell)
   call cmp(n, v, 'status', a%status, b%status)
   call cmp(n, v, 'old_start', a%old_start, b%old_start)
   call cmp(n, v, 'old_piece', a%old_piece, b%old_piece)
   call cmp(n, v, 'new_start', a%new_start, b%new_start)
   call cmp(n, v, 'new_piece', a%new_piece, b%new_piece)
   call cmp(n, v, 'old_track', a%old_track, b%old_track)
   call cmp(n, v, 'new_track', a%new_track, b%new_track)
   call cmp(n, v, 'volume0', a%volume0, b%volume0)
   call cmp(n, v, 'volume1', a%volume1, b%volume1)
   call cmp(n, v, 'faces', a%faces, b%faces)
   call cmp(n, v, 'wall', a%wall, b%wall)
   call cmp(n, v, 'residual', a%residual, b%residual)
   call cmp(n, v, 'outflow', a%outflow, b%outflow)
   call cmp(n, v, 'exact', a%exact, b%exact)
   do k = -1_I4P, 7_I4P
      call cmp(n, v, 'n_status', a%n_status(k), b%n_status(k))
   enddo
   call cmp(n, v, 'untrusted_cells', a%untrusted_cells, b%untrusted_cells)
   call cmp(n, v, 'multi_track_cells', a%multi_track_cells, b%multi_track_cells)
   call cmp(n, v, 'evaluated_cells', a%evaluated_cells, b%evaluated_cells)
   call cmp(n, v, 'evaluations', a%evaluations, b%evaluations)
   call cmp(n, v, 'skipped_nodes', a%skipped_nodes, b%skipped_nodes)
   call cmp(n, v, 'contact_cells', a%contact_cells, b%contact_cells)
   call cmp(n, v, 'narrow_cells', a%narrow_cells, b%narrow_cells)
   endfunction piece_map_2d_mismatches

   function piece_map_3d_mismatches(a, b, verbose) result(n)
   !< Differing elements of two 3-D piece maps (see the module notes).
   type(piece_map_3d_t), intent(in)           :: a, b    !< Maps.
   logical,              intent(in), optional :: verbose !< Print every differing field.
   integer(I4P)                               :: n       !< Mismatches.
   logical                                    :: v       !< Verbose.
   integer(I4P)                               :: k       !< Counter.

   v = .false. ; if (present(verbose)) v = verbose
   n = 0_I4P
   call cmp(n, v, 'nx', a%nx, b%nx)
   call cmp(n, v, 'ny', a%ny, b%ny)
   call cmp(n, v, 'nz', a%nz, b%nz)
   call cmp(n, v, 'nt', a%nt, b%nt)
   call cmp(n, v, 'np0', a%np0, b%np0)
   call cmp(n, v, 'np1', a%np1, b%np1)
   call cmp(n, v, 'tracks', a%tracks, b%tracks)
   call cmp(n, v, 'first_track', a%first_track, b%first_track)
   call cmp(n, v, 'track_cell', a%track_cell, b%track_cell)
   call cmp(n, v, 'status', a%status, b%status)
   call cmp(n, v, 'old_start', a%old_start, b%old_start)
   call cmp(n, v, 'old_piece', a%old_piece, b%old_piece)
   call cmp(n, v, 'new_start', a%new_start, b%new_start)
   call cmp(n, v, 'new_piece', a%new_piece, b%new_piece)
   call cmp(n, v, 'old_track', a%old_track, b%old_track)
   call cmp(n, v, 'new_track', a%new_track, b%new_track)
   call cmp(n, v, 'volume0', a%volume0, b%volume0)
   call cmp(n, v, 'volume1', a%volume1, b%volume1)
   call cmp(n, v, 'faces', a%faces, b%faces)
   call cmp(n, v, 'outflow', a%outflow, b%outflow)
   call cmp(n, v, 'wall', a%wall, b%wall)
   call cmp(n, v, 'residual', a%residual, b%residual)
   call cmp(n, v, 'exact', a%exact, b%exact)
   do k = -1_I4P, 7_I4P
      call cmp(n, v, 'n_status', a%n_status(k), b%n_status(k))
   enddo
   call cmp(n, v, 'untrusted_cells', a%untrusted_cells, b%untrusted_cells)
   call cmp(n, v, 'multi_track_cells', a%multi_track_cells, b%multi_track_cells)
   call cmp(n, v, 'evaluated_cells', a%evaluated_cells, b%evaluated_cells)
   call cmp(n, v, 'evaluations', a%evaluations, b%evaluations)
   call cmp(n, v, 'skipped_nodes', a%skipped_nodes, b%skipped_nodes)
   call cmp(n, v, 'contact_cells', a%contact_cells, b%contact_cells)
   call cmp(n, v, 'narrow_cells', a%narrow_cells, b%narrow_cells)
   endfunction piece_map_3d_mismatches

   subroutine tally(n, verbose, name, m)
   !< Add m mismatches of a field; print it if verbose.
   integer(I4P), intent(inout) :: n       !< Mismatches so far.
   logical,      intent(in)    :: verbose !< Print.
   character(*), intent(in)    :: name    !< Field.
   integer(I4P), intent(in)    :: m       !< Its mismatches.

   if (m == 0_I4P) return
   n = n + m
   if (verbose) print '(A,A,A,I0)', '    piece map field ', name, ' differs: ', m
   endsubroutine tally

   subroutine cmp_i0(n, verbose, name, a, b)
   !< Integer scalar.
   integer(I4P), intent(inout) :: n       !< Mismatches so far.
   logical,      intent(in)    :: verbose !< Print.
   character(*), intent(in)    :: name    !< Field.
   integer(I4P), intent(in)    :: a, b    !< Values.

   call tally(n, verbose, name, merge(1_I4P, 0_I4P, a /= b))
   endsubroutine cmp_i0

   subroutine cmp_l0(n, verbose, name, a, b)
   !< Logical scalar.
   integer(I4P), intent(inout) :: n       !< Mismatches so far.
   logical,      intent(in)    :: verbose !< Print.
   character(*), intent(in)    :: name    !< Field.
   logical,      intent(in)    :: a, b    !< Values.

   call tally(n, verbose, name, merge(1_I4P, 0_I4P, a .neqv. b))
   endsubroutine cmp_l0

   pure function alloc_shape(aa, ab, sa, sb) result(m)
   !< Mismatches from allocation status and shape alone (-1: both allocated with one shape, compare the elements).
   logical,      intent(in) :: aa, ab !< Allocated.
   integer(I4P), intent(in) :: sa(:)  !< Shape of a (meaningful if allocated).
   integer(I4P), intent(in) :: sb(:)  !< Shape of b.
   integer(I4P)             :: m      !< Mismatches, or -1.

   m = -1_I4P
   if (aa .neqv. ab) then
      if (aa) m = max(1_I4P, product(sa))
      if (ab) m = max(1_I4P, product(sb))
      return
   endif
   if (.not. aa) then
      m = 0_I4P ; return
   endif
   if (any(sa /= sb)) m = max(1_I4P, product(sa))
   endfunction alloc_shape

   subroutine cmp_i1(n, verbose, name, a, b)
   !< Integer rank-1 array, elementwise.
   integer(I4P),              intent(inout) :: n       !< Mismatches so far.
   logical,                   intent(in)    :: verbose !< Print.
   character(*),              intent(in)    :: name    !< Field.
   integer(I4P), allocatable, intent(in)    :: a(:), b(:) !< Arrays.
   integer(I4P)                             :: m       !< Mismatches.
   integer(I4P)                             :: sa(1), sb(1) !< Shapes (0 if not allocated).

   sa = 0_I4P ; sb = 0_I4P
   if (allocated(a)) sa = shape(a, kind=I4P)
   if (allocated(b)) sb = shape(b, kind=I4P)
   m = alloc_shape(allocated(a), allocated(b), sa, sb)
   if (m < 0_I4P) m = count(a /= b, kind=I4P)
   call tally(n, verbose, name, m)
   endsubroutine cmp_i1

   subroutine cmp_i2(n, verbose, name, a, b)
   !< Integer rank-2 array, elementwise.
   integer(I4P),              intent(inout) :: n       !< Mismatches so far.
   logical,                   intent(in)    :: verbose !< Print.
   character(*),              intent(in)    :: name    !< Field.
   integer(I4P), allocatable, intent(in)    :: a(:,:), b(:,:) !< Arrays.
   integer(I4P)                             :: m       !< Mismatches.
   integer(I4P)                             :: sa(2), sb(2) !< Shapes (0 if not allocated).

   sa = 0_I4P ; sb = 0_I4P
   if (allocated(a)) sa = shape(a, kind=I4P)
   if (allocated(b)) sb = shape(b, kind=I4P)
   m = alloc_shape(allocated(a), allocated(b), sa, sb)
   if (m < 0_I4P) m = count(a /= b, kind=I4P)
   call tally(n, verbose, name, m)
   endsubroutine cmp_i2

   subroutine cmp_i3(n, verbose, name, a, b)
   !< Integer rank-3 array, elementwise.
   integer(I4P),              intent(inout) :: n       !< Mismatches so far.
   logical,                   intent(in)    :: verbose !< Print.
   character(*),              intent(in)    :: name    !< Field.
   integer(I4P), allocatable, intent(in)    :: a(:,:,:), b(:,:,:) !< Arrays.
   integer(I4P)                             :: m       !< Mismatches.
   integer(I4P)                             :: sa(3), sb(3) !< Shapes (0 if not allocated).

   sa = 0_I4P ; sb = 0_I4P
   if (allocated(a)) sa = shape(a, kind=I4P)
   if (allocated(b)) sb = shape(b, kind=I4P)
   m = alloc_shape(allocated(a), allocated(b), sa, sb)
   if (m < 0_I4P) m = count(a /= b, kind=I4P)
   call tally(n, verbose, name, m)
   endsubroutine cmp_i3

   subroutine cmp_r1(n, verbose, name, a, b)
   !< Real rank-1 array, by bits.
   integer(I4P),           intent(inout) :: n          !< Mismatches so far.
   logical,                intent(in)    :: verbose    !< Print.
   character(*),           intent(in)    :: name       !< Field.
   real(R8P), allocatable, intent(in)    :: a(:), b(:) !< Arrays.
   integer(I4P)                          :: m          !< Mismatches.
   integer(I4P)                          :: sa(1), sb(1) !< Shapes (0 if not allocated).

   sa = 0_I4P ; sb = 0_I4P
   if (allocated(a)) sa = shape(a, kind=I4P)
   if (allocated(b)) sb = shape(b, kind=I4P)
   m = alloc_shape(allocated(a), allocated(b), sa, sb)
   if (m < 0_I4P) m = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   call tally(n, verbose, name, m)
   endsubroutine cmp_r1

   subroutine cmp_r2(n, verbose, name, a, b)
   !< Real rank-2 array, by bits.
   integer(I4P),           intent(inout) :: n              !< Mismatches so far.
   logical,                intent(in)    :: verbose        !< Print.
   character(*),           intent(in)    :: name           !< Field.
   real(R8P), allocatable, intent(in)    :: a(:,:), b(:,:) !< Arrays.
   integer(I4P)                          :: m              !< Mismatches.
   integer(I4P)                          :: sa(2), sb(2)   !< Shapes (0 if not allocated).

   sa = 0_I4P ; sb = 0_I4P
   if (allocated(a)) sa = shape(a, kind=I4P)
   if (allocated(b)) sb = shape(b, kind=I4P)
   m = alloc_shape(allocated(a), allocated(b), sa, sb)
   if (m < 0_I4P) m = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   call tally(n, verbose, name, m)
   endsubroutine cmp_r2
endmodule fossil_piece_map_compare
