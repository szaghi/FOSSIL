!< FOSSIL test support: bitwise comparison of two 2-D cut-cell results (issue #100 step 2).
module fossil_cut_cells_compare
!< FOSSIL test support: bitwise comparison of two 2-D cut-cell results (issue #100 step 2).
!<
!< `cut_cells_2d_mismatches(a, b)` counts the differing elements of every output of `curve%cut_cells_device` (moments
!< order 0, no boundaries, no quadrature): the scalars nx, ny, np, moments_order; the cell arrays alpha, aperture,
!< centroid, wet_length, wet_centroid, wet_normal, pieces, first_piece, closed, node_in; the piece arrays piece_cell,
!< piece_area, piece_centroid, piece_wet_length, piece_wet_normal, piece_aperture, piece_iv_start, piece_iv_edge,
!< piece_iv_lo, piece_iv_hi. Reals compare by their bits (`transfer` to a 64-bit integer: -0 /= +0, NaN = NaN of the
!< same payload). An array allocated in one result only, or of another shape, counts as max(1, its size) mismatches;
!< the host-only outputs (moments, boundaries, quadrature) must be unallocated in both. `cut_cells_twin_mismatches(c,
!< x, y)` runs `c%cut_cells` and `c%cut_cells_device` on the block x, y, with and without the intervals, and sums
!< them. Used by `fossil_test_cut_cells_device_twin` (CI, host) and by the curve bench (device vs serial).
use fossil, only : curve_object, cut_cells_2d_t
use penf,   only : I4P, I8P, R8P

implicit none
private
public :: cut_cells_2d_mismatches, cut_cells_twin_mismatches

contains
   function cut_cells_twin_mismatches(c, x, y, status) result(n)
   !< Mismatches of `c%cut_cells_device` against `c%cut_cells` on the block x, y, with the intervals and without (sum);
   !< status: the sum of the statuses (0: all ok). Each twin call is made twice (issue #100 step 6a): copied back, and
   !< resident (`resident=.true.`, then `update_host`, compared, then the next call into the same `d`, which releases it;
   !< `exit_device` at the end). A resident result must say `on_device` (unless the edges were refused); a copy of it
   !< must not be (its `exit_device` leaves the original mapped: the following `update_host` would fail otherwise).
   type(curve_object), intent(in)            :: c          !< Curve.
   real(R8P),          intent(in)            :: x(:), y(:) !< Cell edges.
   integer(I4P),       intent(out), optional :: status     !< Sum of the statuses.
   integer(I4P)                              :: n          !< Mismatches.
   type(cut_cells_2d_t)                      :: h, d, b    !< Host, twin, copy of the twin.
   integer(I4P)                              :: s(6)       !< Statuses.
   logical                                   :: wi         !< With the intervals.
   integer(I4P)                              :: k          !< Pass.

   n = 0_I4P
   do k = 1, 2
      wi = k == 1
      call c%cut_cells(x=x, y=y, cut=h, status=s(3 * k - 2), intervals=wi)
      call c%cut_cells_device(x=x, y=y, cut=d, status=s(3 * k - 1), intervals=wi)
      n = n + cut_cells_2d_mismatches(h, d)
      call poison(d) ! the next call reuses the host arrays of d (same shapes): every element must be rewritten
      call c%cut_cells_device(x=x, y=y, cut=d, status=s(3 * k), intervals=wi, resident=.true.)
      if (d%on_device .neqv. allocated(d%alpha)) n = n + 1_I4P
      if (k == 1) then ! a copy is not resident: releasing it leaves d mapped
         b = d
         call b%exit_device
         if (b%on_device) n = n + 1_I4P
      endif
      call d%update_host
      n = n + cut_cells_2d_mismatches(h, d)
   enddo
   call d%exit_device
   if (d%on_device) n = n + 1_I4P
   if (present(status)) status = sum(s)
   endfunction cut_cells_twin_mismatches

   subroutine poison(d)
   !< Overwrite every output of d with values no result holds (issue #100 step 6a: the device path reuses them).
   type(cut_cells_2d_t), intent(inout) :: d !< Result.
   real(R8P), parameter :: R = -7.77e300_R8P !< Real poison.
   integer(I4P), parameter :: I = -777777_I4P !< Integer poison.

   if (.not. allocated(d%alpha)) return
   d%alpha = R ; d%aperture = R ; d%centroid = R ; d%wet_length = R ; d%wet_centroid = R ; d%wet_normal = R
   d%pieces = I ; d%first_piece = I ; d%closed = .not. d%closed ; d%node_in = .not. d%node_in
   d%piece_cell = I ; d%piece_area = R ; d%piece_centroid = R ; d%piece_wet_length = R ; d%piece_wet_normal = R
   d%piece_aperture = R ; d%piece_iv_start = I ; d%piece_iv_edge = I ; d%piece_iv_lo = R ; d%piece_iv_hi = R
   endsubroutine poison

   function cut_cells_2d_mismatches(a, b) result(n)
   !< Differing elements of the device-twin outputs of a and b (0: bit-identical).
   type(cut_cells_2d_t), intent(in) :: a, b !< Results.
   integer(I4P)                     :: n    !< Mismatches.

   n = 0_I4P
   if (a%nx /= b%nx) n = n + 1_I4P
   if (a%ny /= b%ny) n = n + 1_I4P
   if (a%np /= b%np) n = n + 1_I4P
   if (a%moments_order /= b%moments_order) n = n + 1_I4P
   ! cell arrays
   n = n + real_mismatches(allocated(a%alpha), allocated(b%alpha), shape_of(a%alpha), shape_of(b%alpha), &
                           flat_r2(a%alpha), flat_r2(b%alpha))
   n = n + real_mismatches(allocated(a%aperture), allocated(b%aperture), shape_of3(a%aperture), shape_of3(b%aperture), &
                           flat_r3(a%aperture), flat_r3(b%aperture))
   n = n + real_mismatches(allocated(a%centroid), allocated(b%centroid), shape_of3(a%centroid), shape_of3(b%centroid), &
                           flat_r3(a%centroid), flat_r3(b%centroid))
   n = n + real_mismatches(allocated(a%wet_length), allocated(b%wet_length), shape_of(a%wet_length), &
                           shape_of(b%wet_length), flat_r2(a%wet_length), flat_r2(b%wet_length))
   n = n + real_mismatches(allocated(a%wet_centroid), allocated(b%wet_centroid), shape_of3(a%wet_centroid), &
                           shape_of3(b%wet_centroid), flat_r3(a%wet_centroid), flat_r3(b%wet_centroid))
   n = n + real_mismatches(allocated(a%wet_normal), allocated(b%wet_normal), shape_of3(a%wet_normal), &
                           shape_of3(b%wet_normal), flat_r3(a%wet_normal), flat_r3(b%wet_normal))
   n = n + int_mismatches(allocated(a%pieces), allocated(b%pieces), shape_of_i(a%pieces), shape_of_i(b%pieces), &
                          flat_i2(a%pieces), flat_i2(b%pieces))
   n = n + int_mismatches(allocated(a%first_piece), allocated(b%first_piece), shape_of_i(a%first_piece), &
                          shape_of_i(b%first_piece), flat_i2(a%first_piece), flat_i2(b%first_piece))
   n = n + log_mismatches(allocated(a%closed), allocated(b%closed), shape_of_l(a%closed), shape_of_l(b%closed), &
                          flat_l2(a%closed), flat_l2(b%closed))
   n = n + log_mismatches(allocated(a%node_in), allocated(b%node_in), shape_of_l(a%node_in), shape_of_l(b%node_in), &
                          flat_l2(a%node_in), flat_l2(b%node_in))
   ! piece arrays
   n = n + int_mismatches(allocated(a%piece_cell), allocated(b%piece_cell), shape_of_i(a%piece_cell), &
                          shape_of_i(b%piece_cell), flat_i2(a%piece_cell), flat_i2(b%piece_cell))
   n = n + real_mismatches(allocated(a%piece_area), allocated(b%piece_area), shape_of1(a%piece_area), &
                           shape_of1(b%piece_area), flat_r1(a%piece_area), flat_r1(b%piece_area))
   n = n + real_mismatches(allocated(a%piece_centroid), allocated(b%piece_centroid), shape_of(a%piece_centroid), &
                           shape_of(b%piece_centroid), flat_r2(a%piece_centroid), flat_r2(b%piece_centroid))
   n = n + real_mismatches(allocated(a%piece_wet_length), allocated(b%piece_wet_length), shape_of1(a%piece_wet_length), &
                           shape_of1(b%piece_wet_length), flat_r1(a%piece_wet_length), flat_r1(b%piece_wet_length))
   n = n + real_mismatches(allocated(a%piece_wet_normal), allocated(b%piece_wet_normal), shape_of(a%piece_wet_normal), &
                           shape_of(b%piece_wet_normal), flat_r2(a%piece_wet_normal), flat_r2(b%piece_wet_normal))
   n = n + real_mismatches(allocated(a%piece_aperture), allocated(b%piece_aperture), shape_of(a%piece_aperture), &
                           shape_of(b%piece_aperture), flat_r2(a%piece_aperture), flat_r2(b%piece_aperture))
   n = n + int_mismatches(allocated(a%piece_iv_start), allocated(b%piece_iv_start), shape_of_i1(a%piece_iv_start), &
                          shape_of_i1(b%piece_iv_start), flat_i1(a%piece_iv_start), flat_i1(b%piece_iv_start))
   n = n + int_mismatches(allocated(a%piece_iv_edge), allocated(b%piece_iv_edge), shape_of_i1(a%piece_iv_edge), &
                          shape_of_i1(b%piece_iv_edge), flat_i1(a%piece_iv_edge), flat_i1(b%piece_iv_edge))
   n = n + real_mismatches(allocated(a%piece_iv_lo), allocated(b%piece_iv_lo), shape_of1(a%piece_iv_lo), &
                           shape_of1(b%piece_iv_lo), flat_r1(a%piece_iv_lo), flat_r1(b%piece_iv_lo))
   n = n + real_mismatches(allocated(a%piece_iv_hi), allocated(b%piece_iv_hi), shape_of1(a%piece_iv_hi), &
                           shape_of1(b%piece_iv_hi), flat_r1(a%piece_iv_hi), flat_r1(b%piece_iv_hi))
   ! host-only outputs: unallocated in both
   if (allocated(a%piece_moments) .or. allocated(b%piece_moments)) n = n + 1_I4P
   if (allocated(a%piece_wall_moments) .or. allocated(b%piece_wall_moments)) n = n + 1_I4P
   if (allocated(a%piece_seg_start) .or. allocated(b%piece_seg_start)) n = n + 1_I4P
   if (allocated(a%piece_seg) .or. allocated(b%piece_seg)) n = n + 1_I4P
   if (allocated(a%piece_seg_wall) .or. allocated(b%piece_seg_wall)) n = n + 1_I4P
   if (allocated(a%quad_start) .or. allocated(b%quad_start)) n = n + 1_I4P
   if (allocated(a%wall_quad_start) .or. allocated(b%wall_quad_start)) n = n + 1_I4P
   endfunction cut_cells_2d_mismatches

   function real_mismatches(la, lb, sa, sb, fa, fb) result(n)
   !< Bitwise mismatches of two flattened real arrays.
   logical,      intent(in) :: la, lb       !< Allocated.
   integer(I4P), intent(in) :: sa(3), sb(3) !< Shapes (padded with 1).
   real(R8P),    intent(in) :: fa(:), fb(:) !< Elements.
   integer(I4P)             :: n            !< Mismatches.
   integer(I4P)             :: k            !< Counter.

   n = 0_I4P
   if (.not. (la .or. lb)) return
   if ((la .neqv. lb) .or. any(sa /= sb)) then
      n = max(1_I4P, size(fa, kind=I4P), size(fb, kind=I4P)) ; return
   endif
   do k = 1_I4P, size(fa, kind=I4P)
      if (transfer(fa(k), 0_I8P) /= transfer(fb(k), 0_I8P)) n = n + 1_I4P
   enddo
   endfunction real_mismatches

   function int_mismatches(la, lb, sa, sb, fa, fb) result(n)
   !< Mismatches of two flattened integer arrays.
   logical,      intent(in) :: la, lb       !< Allocated.
   integer(I4P), intent(in) :: sa(3), sb(3) !< Shapes.
   integer(I4P), intent(in) :: fa(:), fb(:) !< Elements.
   integer(I4P)             :: n            !< Mismatches.

   n = 0_I4P
   if (.not. (la .or. lb)) return
   if ((la .neqv. lb) .or. any(sa /= sb)) then
      n = max(1_I4P, size(fa, kind=I4P), size(fb, kind=I4P)) ; return
   endif
   n = count(fa /= fb, kind=I4P)
   endfunction int_mismatches

   function log_mismatches(la, lb, sa, sb, fa, fb) result(n)
   !< Mismatches of two flattened logical arrays.
   logical,      intent(in) :: la, lb       !< Allocated.
   integer(I4P), intent(in) :: sa(3), sb(3) !< Shapes.
   logical,      intent(in) :: fa(:), fb(:) !< Elements.
   integer(I4P)             :: n            !< Mismatches.

   n = 0_I4P
   if (.not. (la .or. lb)) return
   if ((la .neqv. lb) .or. any(sa /= sb)) then
      n = max(1_I4P, size(fa, kind=I4P), size(fb, kind=I4P)) ; return
   endif
   n = count(fa .neqv. fb, kind=I4P)
   endfunction log_mismatches

   function shape_of1(a) result(s)
   !< Padded shape of a rank-1 real array (zeros if unallocated).
   real(R8P), allocatable, intent(in) :: a(:) !< Array.
   integer(I4P)                       :: s(3) !< Shape.

   s = [0_I4P, 1_I4P, 1_I4P] ; if (allocated(a)) s(1) = size(a, kind=I4P)
   endfunction shape_of1

   function shape_of(a) result(s)
   !< Padded shape of a rank-2 real array.
   real(R8P), allocatable, intent(in) :: a(:,:) !< Array.
   integer(I4P)                       :: s(3)   !< Shape.

   s = [0_I4P, 0_I4P, 1_I4P] ; if (allocated(a)) s(1:2) = shape(a, kind=I4P)
   endfunction shape_of

   function shape_of3(a) result(s)
   !< Shape of a rank-3 real array.
   real(R8P), allocatable, intent(in) :: a(:,:,:) !< Array.
   integer(I4P)                       :: s(3)     !< Shape.

   s = 0_I4P ; if (allocated(a)) s = shape(a, kind=I4P)
   endfunction shape_of3

   function shape_of_i1(a) result(s)
   !< Padded shape of a rank-1 integer array.
   integer(I4P), allocatable, intent(in) :: a(:) !< Array.
   integer(I4P)                          :: s(3) !< Shape.

   s = [0_I4P, 1_I4P, 1_I4P] ; if (allocated(a)) s(1) = size(a, kind=I4P)
   endfunction shape_of_i1

   function shape_of_i(a) result(s)
   !< Padded shape of a rank-2 integer array.
   integer(I4P), allocatable, intent(in) :: a(:,:) !< Array.
   integer(I4P)                          :: s(3)   !< Shape.

   s = [0_I4P, 0_I4P, 1_I4P] ; if (allocated(a)) s(1:2) = shape(a, kind=I4P)
   endfunction shape_of_i

   function shape_of_l(a) result(s)
   !< Padded shape of a rank-2 logical array.
   logical, allocatable, intent(in) :: a(:,:) !< Array.
   integer(I4P)                     :: s(3)   !< Shape.

   s = [0_I4P, 0_I4P, 1_I4P] ; if (allocated(a)) s(1:2) = shape(a, kind=I4P)
   endfunction shape_of_l

   function flat_r1(a) result(f)
   !< Elements of a rank-1 real array (none if unallocated).
   real(R8P), allocatable, intent(in) :: a(:) !< Array.
   real(R8P), allocatable             :: f(:) !< Elements.

   if (allocated(a)) then
      f = a
   else
      allocate(f(0))
   endif
   endfunction flat_r1

   function flat_r2(a) result(f)
   !< Elements of a rank-2 real array in array element order.
   real(R8P), allocatable, intent(in) :: a(:,:) !< Array.
   real(R8P), allocatable             :: f(:)   !< Elements.

   if (allocated(a)) then
      f = reshape(a, [size(a)])
   else
      allocate(f(0))
   endif
   endfunction flat_r2

   function flat_r3(a) result(f)
   !< Elements of a rank-3 real array in array element order.
   real(R8P), allocatable, intent(in) :: a(:,:,:) !< Array.
   real(R8P), allocatable             :: f(:)     !< Elements.

   if (allocated(a)) then
      f = reshape(a, [size(a)])
   else
      allocate(f(0))
   endif
   endfunction flat_r3

   function flat_i1(a) result(f)
   !< Elements of a rank-1 integer array.
   integer(I4P), allocatable, intent(in) :: a(:) !< Array.
   integer(I4P), allocatable             :: f(:) !< Elements.

   if (allocated(a)) then
      f = a
   else
      allocate(f(0))
   endif
   endfunction flat_i1

   function flat_i2(a) result(f)
   !< Elements of a rank-2 integer array.
   integer(I4P), allocatable, intent(in) :: a(:,:) !< Array.
   integer(I4P), allocatable             :: f(:)   !< Elements.

   if (allocated(a)) then
      f = reshape(a, [size(a)])
   else
      allocate(f(0))
   endif
   endfunction flat_i2

   function flat_l2(a) result(f)
   !< Elements of a rank-2 logical array.
   logical, allocatable, intent(in) :: a(:,:) !< Array.
   logical, allocatable             :: f(:)   !< Elements.

   if (allocated(a)) then
      f = reshape(a, [size(a)])
   else
      allocate(f(0))
   endif
   endfunction flat_l2
endmodule fossil_cut_cells_compare
