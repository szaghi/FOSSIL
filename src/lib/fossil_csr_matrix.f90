!< Compressed-sparse-row (CSR) matrix type (issue #18 §2.1 foundation).
!<
!< Minimal CSR storage with the small set of operations the Tier-2 DDG
!< primitives need: build by row-coordinate triplet append + finalize,
!< matrix-vector multiply, row sum, symmetry check, value lookup. No
!< solver, no factorization — this module is *storage and accessors*
!< only. Callers bring their own solver (SuiteSparse, MUMPS, PETSc, ...)
!< when they need one; the same `csr_matrix_t` is the interchange format.
!<
!< Storage convention (standard CSR):
!<   - `row_ptr(1:n_rows + 1)` — cumulative count of nonzeros per row;
!<     `row_ptr(i):row_ptr(i+1)-1` is the slice of `col_idx` / `val`
!<     belonging to row `i` (1-based on the storage indices themselves).
!<   - `col_idx(1:nnz)` — column indices, sorted ascending within each row.
!<   - `val(1:nnz)`     — corresponding nonzero values.
!<
!< Building a matrix:
!<   1. `call M%initialize(n_rows, n_cols, est_nnz)` reserves capacity.
!<   2. Call `M%append(i, j, v)` for every triplet. Duplicates accumulate
!<      (so two appends to (i, j) sum their values — this is what we want
!<      for the cotangent Laplacian, where each interior edge contributes
!<      a cotangent term from each of its two opposite vertices).
!<   3. `call M%finalize()` sorts each row's column indices ascending and
!<      coalesces duplicate column entries within a row.

module fossil_csr_matrix
!< CSR sparse-matrix type used as the storage format for Tier-2 DDG.

use penf, only : I4P, R8P

implicit none
private

public :: csr_matrix_t
public :: CSR_STATUS_OK, CSR_STATUS_BAD_INPUT, CSR_STATUS_OUT_OF_RANGE

integer(I4P), parameter :: CSR_STATUS_OK            = 0_I4P
integer(I4P), parameter :: CSR_STATUS_BAD_INPUT     = 1_I4P
integer(I4P), parameter :: CSR_STATUS_OUT_OF_RANGE  = 2_I4P  !< Row / column index outside [1, n_rows] / [1, n_cols].

integer(I4P), parameter :: CSR_INITIAL_NNZ_PER_ROW  = 8_I4P  !< Per-row capacity guess; rows grow geometrically beyond it.

type :: csr_matrix_t
   !< CSR-stored sparse matrix.
   integer(I4P)              :: n_rows  = 0_I4P     !< Number of rows.
   integer(I4P)              :: n_cols  = 0_I4P     !< Number of columns.
   integer(I4P)              :: nnz     = 0_I4P     !< Number of nonzeros (set by finalize; 0 before).
   integer(I4P), allocatable :: row_ptr(:)          !< Length n_rows + 1.
   integer(I4P), allocatable :: col_idx(:)          !< Length nnz; per-row sorted ascending after finalize.
   real(R8P),    allocatable :: val(:)              !< Length nnz; aligned with col_idx.
   logical                   :: is_finalized = .false. !< True after `finalize` runs (no further `append` allowed).
   ! Build-time scratch — discarded by finalize.
   integer(I4P), allocatable :: build_row_counts(:) !< Per-row count of appended (col_idx, val) entries.
   integer(I4P)              :: build_capacity = 0_I4P !< Leading dimension of build_col_idx / build_val. Shared across all rows (2D layout).
   integer(I4P), allocatable :: build_col_idx(:,:)  !< (build_capacity, n_rows) growable; trimmed at finalize.
   real(R8P),    allocatable :: build_val(:,:)      !< (build_capacity, n_rows) growable.
   contains
      procedure, pass(self) :: initialize         => csr_initialize
      procedure, pass(self) :: destroy            => csr_destroy
      procedure, pass(self) :: append             => csr_append
      procedure, pass(self) :: finalize           => csr_finalize
      procedure, pass(self) :: get_nrows          => csr_get_nrows
      procedure, pass(self) :: get_ncols          => csr_get_ncols
      procedure, pass(self) :: get_nnz            => csr_get_nnz
      procedure, pass(self) :: multiply_vector    => csr_multiply_vector
      procedure, pass(self) :: row_sum            => csr_row_sum
      procedure, pass(self) :: get_value          => csr_get_value
      procedure, pass(self) :: is_symmetric       => csr_is_symmetric
endtype csr_matrix_t

contains

   subroutine csr_initialize(self, n_rows, n_cols, status)
   !< Reserve build-time capacity. Per-row growable arrays start at
   !< `CSR_INITIAL_NNZ_PER_ROW` and double on demand.
   class(csr_matrix_t), intent(inout)         :: self
   integer(I4P),        intent(in)            :: n_rows
   integer(I4P),        intent(in)            :: n_cols
   integer(I4P),        intent(out), optional :: status

   if (present(status)) status = CSR_STATUS_OK
   call self%destroy
   if (n_rows < 0_I4P .or. n_cols < 0_I4P) then
      if (present(status)) status = CSR_STATUS_BAD_INPUT
      return
   endif
   self%n_rows = n_rows
   self%n_cols = n_cols
   self%nnz    = 0_I4P
   self%is_finalized = .false.
   if (n_rows == 0_I4P) return
   allocate(self%build_row_counts(n_rows))
   self%build_row_counts = 0_I4P
   self%build_capacity   = CSR_INITIAL_NNZ_PER_ROW
   allocate(self%build_col_idx(CSR_INITIAL_NNZ_PER_ROW, n_rows))
   allocate(self%build_val(CSR_INITIAL_NNZ_PER_ROW, n_rows))
   endsubroutine csr_initialize

   subroutine csr_destroy(self)
   !< Release all storage. Safe to call repeatedly.
   class(csr_matrix_t), intent(inout) :: self
   if (allocated(self%row_ptr))          deallocate(self%row_ptr)
   if (allocated(self%col_idx))          deallocate(self%col_idx)
   if (allocated(self%val))              deallocate(self%val)
   if (allocated(self%build_row_counts)) deallocate(self%build_row_counts)
   if (allocated(self%build_col_idx))    deallocate(self%build_col_idx)
   if (allocated(self%build_val))        deallocate(self%build_val)
   self%n_rows = 0_I4P
   self%n_cols = 0_I4P
   self%nnz    = 0_I4P
   self%build_capacity = 0_I4P
   self%is_finalized = .false.
   endsubroutine csr_destroy

   subroutine csr_append(self, row, col, value, status)
   !< Append a single triplet `(row, col, value)` to the build buffers.
   !< Duplicates (same row + col on multiple calls) are NOT coalesced
   !< here — they are summed by `finalize`. This matters for the
   !< cotangent Laplacian builder, which contributes a cotangent term to
   !< the same `(i, j)` entry from each of the two triangles sharing
   !< edge `(i, j)`.
   class(csr_matrix_t), intent(inout)         :: self
   integer(I4P),        intent(in)            :: row
   integer(I4P),        intent(in)            :: col
   real(R8P),           intent(in)            :: value
   integer(I4P),        intent(out), optional :: status
   integer(I4P)                               :: slot

   if (present(status)) status = CSR_STATUS_OK
   if (self%is_finalized) then
      if (present(status)) status = CSR_STATUS_BAD_INPUT
      return
   endif
   if (row < 1_I4P .or. row > self%n_rows .or. col < 1_I4P .or. col > self%n_cols) then
      if (present(status)) status = CSR_STATUS_OUT_OF_RANGE
      return
   endif
   slot = self%build_row_counts(row) + 1_I4P
   if (slot > self%build_capacity) call grow_capacity(self=self)
   self%build_col_idx(slot, row) = col
   self%build_val(slot, row)     = value
   self%build_row_counts(row)    = slot
   endsubroutine csr_append

   subroutine grow_capacity(self)
   !< Double the leading dimension (per-row capacity) of the build arrays.
   !< The 2D layout means all rows share the same physical capacity —
   !< we widen them together.
   !<
   !< (Earlier version of this routine tracked a per-row capacity array
   !< and updated only one entry on growth; that lost data because other
   !< rows still saw the stale "old_cap" but the underlying 2D buffer
   !< was already wider. Single shared `build_capacity` scalar fixes it.)
   type(csr_matrix_t), intent(inout) :: self
   integer(I4P)                      :: new_cap
   integer(I4P), allocatable         :: tmp_col(:,:)
   real(R8P),    allocatable         :: tmp_val(:,:)
   integer(I4P)                      :: i, count_i

   new_cap = 2_I4P * self%build_capacity
   allocate(tmp_col(new_cap, self%n_rows))
   allocate(tmp_val(new_cap, self%n_rows))
   do i = 1_I4P, self%n_rows
      count_i = self%build_row_counts(i)
      if (count_i > 0_I4P) then
         tmp_col(1:count_i, i) = self%build_col_idx(1:count_i, i)
         tmp_val(1:count_i, i) = self%build_val(1:count_i, i)
      endif
   enddo
   call move_alloc(from=tmp_col, to=self%build_col_idx)
   call move_alloc(from=tmp_val, to=self%build_val)
   self%build_capacity = new_cap
   endsubroutine grow_capacity

   subroutine csr_finalize(self, status)
   !< Sort each row's column indices ascending, coalesce duplicates by
   !< summing their values, pack into the canonical `(row_ptr, col_idx,
   !< val)` arrays, and release the build scratch.
   !<
   !< After this call, `append` returns BAD_INPUT (the matrix is read-only)
   !< and the row_ptr/col_idx/val arrays are the canonical CSR storage.
   class(csr_matrix_t), intent(inout)         :: self
   integer(I4P),        intent(out), optional :: status
   integer(I4P)                               :: i, k, count_i, nnz_total
   integer(I4P), allocatable                  :: row_col(:)
   real(R8P),    allocatable                  :: row_val(:)
   integer(I4P)                               :: unique_count, j

   if (present(status)) status = CSR_STATUS_OK
   if (self%is_finalized) return
   if (self%n_rows == 0_I4P) then
      allocate(self%row_ptr(1))
      self%row_ptr(1) = 1_I4P
      self%nnz = 0_I4P
      allocate(self%col_idx(0))
      allocate(self%val(0))
      self%is_finalized = .true.
      return
   endif

   ! First pass: per-row sort + duplicate coalesce, in place on the build buffers.
   nnz_total = 0_I4P
   do i = 1_I4P, self%n_rows
      count_i = self%build_row_counts(i)
      if (count_i == 0_I4P) cycle
      allocate(row_col(count_i))
      allocate(row_val(count_i))
      row_col(1:count_i) = self%build_col_idx(1:count_i, i)
      row_val(1:count_i) = self%build_val(1:count_i, i)
      call sort_and_coalesce(col=row_col, val=row_val, n=count_i, unique=unique_count)
      self%build_col_idx(1:unique_count, i) = row_col(1:unique_count)
      self%build_val(1:unique_count, i)     = row_val(1:unique_count)
      self%build_row_counts(i) = unique_count
      nnz_total = nnz_total + unique_count
      deallocate(row_col, row_val)
   enddo

   ! Second pass: pack into canonical CSR.
   allocate(self%row_ptr(self%n_rows + 1_I4P))
   allocate(self%col_idx(nnz_total))
   allocate(self%val(nnz_total))
   k = 1_I4P
   self%row_ptr(1) = 1_I4P
   do i = 1_I4P, self%n_rows
      count_i = self%build_row_counts(i)
      do j = 1_I4P, count_i
         self%col_idx(k) = self%build_col_idx(j, i)
         self%val(k)     = self%build_val(j, i)
         k = k + 1_I4P
      enddo
      self%row_ptr(i + 1_I4P) = k
   enddo
   self%nnz = nnz_total

   ! Release build scratch.
   if (allocated(self%build_row_counts)) deallocate(self%build_row_counts)
   if (allocated(self%build_col_idx))    deallocate(self%build_col_idx)
   if (allocated(self%build_val))        deallocate(self%build_val)
   self%build_capacity = 0_I4P
   self%is_finalized = .true.
   endsubroutine csr_finalize

   pure subroutine sort_and_coalesce(col, val, n, unique)
   !< In-place insertion sort by column, then merge runs of equal columns
   !< (summing their values). On exit, `col(1:unique)` is sorted ascending
   !< with no duplicates and `val(1:unique)` carries the summed values.
   integer(I4P), intent(inout) :: col(:)
   real(R8P),    intent(inout) :: val(:)
   integer(I4P), intent(in)    :: n
   integer(I4P), intent(out)   :: unique
   integer(I4P)                :: i, j, key_col
   real(R8P)                   :: key_val
   integer(I4P)                :: out_idx

   ! Insertion sort by col, carrying val along.
   do i = 2_I4P, n
      key_col = col(i)
      key_val = val(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (col(j) <= key_col) exit
         col(j + 1_I4P) = col(j)
         val(j + 1_I4P) = val(j)
         j = j - 1_I4P
      enddo
      col(j + 1_I4P) = key_col
      val(j + 1_I4P) = key_val
   enddo

   ! Coalesce in place: scan and merge equal-column runs.
   if (n == 0_I4P) then
      unique = 0_I4P
      return
   endif
   out_idx = 1_I4P
   do i = 2_I4P, n
      if (col(i) == col(out_idx)) then
         val(out_idx) = val(out_idx) + val(i)
      else
         out_idx = out_idx + 1_I4P
         col(out_idx) = col(i)
         val(out_idx) = val(i)
      endif
   enddo
   unique = out_idx
   endsubroutine sort_and_coalesce

   pure function csr_get_nrows(self) result(n)
   class(csr_matrix_t), intent(in) :: self
   integer(I4P)                    :: n
   n = self%n_rows
   endfunction csr_get_nrows

   pure function csr_get_ncols(self) result(n)
   class(csr_matrix_t), intent(in) :: self
   integer(I4P)                    :: n
   n = self%n_cols
   endfunction csr_get_ncols

   pure function csr_get_nnz(self) result(n)
   class(csr_matrix_t), intent(in) :: self
   integer(I4P)                    :: n
   n = self%nnz
   endfunction csr_get_nnz

   subroutine csr_multiply_vector(self, x, y, status)
   !< Compute `y = A * x`. Requires the matrix to be finalized.
   class(csr_matrix_t), intent(in)            :: self
   real(R8P),           intent(in)            :: x(:)
   real(R8P),           intent(out)           :: y(:)
   integer(I4P),        intent(out), optional :: status
   integer(I4P)                               :: i, k

   if (present(status)) status = CSR_STATUS_OK
   if (.not. self%is_finalized) then
      y = 0._R8P
      if (present(status)) status = CSR_STATUS_BAD_INPUT
      return
   endif
   if (size(x, kind=I4P) /= self%n_cols .or. size(y, kind=I4P) /= self%n_rows) then
      y = 0._R8P
      if (present(status)) status = CSR_STATUS_OUT_OF_RANGE
      return
   endif
   y = 0._R8P
   do i = 1_I4P, self%n_rows
      do k = self%row_ptr(i), self%row_ptr(i + 1_I4P) - 1_I4P
         y(i) = y(i) + self%val(k) * x(self%col_idx(k))
      enddo
   enddo
   endsubroutine csr_multiply_vector

   pure function csr_row_sum(self, row) result(s)
   !< Sum of all nonzero values in `row`. Out of range → 0.
   class(csr_matrix_t), intent(in) :: self
   integer(I4P),        intent(in) :: row
   real(R8P)                       :: s
   integer(I4P)                    :: k

   s = 0._R8P
   if (.not. self%is_finalized) return
   if (row < 1_I4P .or. row > self%n_rows) return
   do k = self%row_ptr(row), self%row_ptr(row + 1_I4P) - 1_I4P
      s = s + self%val(k)
   enddo
   endfunction csr_row_sum

   pure function csr_get_value(self, row, col) result(v)
   !< Return the value at `(row, col)`. Linear scan within the row — fine
   !< for the small row sizes typical of mesh-derived matrices.
   class(csr_matrix_t), intent(in) :: self
   integer(I4P),        intent(in) :: row, col
   real(R8P)                       :: v
   integer(I4P)                    :: k

   v = 0._R8P
   if (.not. self%is_finalized) return
   if (row < 1_I4P .or. row > self%n_rows) return
   do k = self%row_ptr(row), self%row_ptr(row + 1_I4P) - 1_I4P
      if (self%col_idx(k) == col) then
         v = self%val(k)
         return
      endif
   enddo
   endfunction csr_get_value

   function csr_is_symmetric(self, tol) result(yes)
   !< True iff `|A(i,j) - A(j,i)| <= tol` for every stored nonzero. The
   !< check is O(nnz × log row_density) due to the linear-scan lookup of
   !< `A(j, i)` for each stored `A(i, j)`. Fine for the typical mesh-
   !< Laplacian sizes (~10k vertices × ~6 incident edges each).
   class(csr_matrix_t), intent(in) :: self
   real(R8P),           intent(in) :: tol
   logical                         :: yes
   integer(I4P)                    :: i, k, j
   real(R8P)                       :: aij, aji

   yes = .true.
   if (.not. self%is_finalized) then
      yes = .false.
      return
   endif
   if (self%n_rows /= self%n_cols) then
      yes = .false.
      return
   endif
   do i = 1_I4P, self%n_rows
      do k = self%row_ptr(i), self%row_ptr(i + 1_I4P) - 1_I4P
         j = self%col_idx(k)
         aij = self%val(k)
         aji = self%get_value(row=j, col=i)
         if (abs(aij - aji) > tol) then
            yes = .false.
            return
         endif
      enddo
   enddo
   endfunction csr_is_symmetric

endmodule fossil_csr_matrix
