!< FOSSIL, signed-distance field on a rectilinear (possibly stretched) grid, issue #32.
module fossil_distance_field
!< FOSSIL, signed-distance field on a rectilinear (possibly stretched) grid, issue #32.
!<
!< `distance_field_object` holds the signed distance `phi(i,j,k)` from the points
!< `(x(i), y(j), z(k))` of a rectilinear grid to a surface. The grid is given by three
!< strictly increasing 1D coordinate arrays (typically cell centres); the spacing may vary
!< along each axis (stretched grids).
!<
!< The field is built in stages:
!<  1. **Band** (`compute_band`): the cells within `band` of the surface get the exact
!<     signed distance, computed in one `surface%distance_many` batch. Candidates come from
!<     rasterizing each facet's bounding box, inflated by `band`, onto the index space
!<     (binary search on the 1D coordinates, so stretched grids cost the same).
!<
!< Halo: the outer `ng` layers of cells on each side are halo cells. They are never
!< computed; a caller with a block decomposition fills them from the neighbouring blocks.
!< `phi` is public for that reason. FOSSIL itself stays MPI-free.
!<
!< Sign convention: negative inside, positive outside, as `surface%distance`.

use fossil, only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL, STATUS_OK, STATUS_INVALID_INPUT
use penf, only : I4P, I8P, R8P, MaxR8P
use vecfor, only : vector_R8P

implicit none
private
public :: distance_field_object
public :: FAR_VALUE

real(R8P), parameter :: FAR_VALUE = MaxR8P !< Magnitude of a cell whose distance is not known yet.

type :: distance_field_object
   !< Signed-distance field on a rectilinear grid (issue #32).
   private
   integer(I4P)                       :: ni = 0_I4P          !< Points along x (halo included).
   integer(I4P)                       :: nj = 0_I4P          !< Points along y (halo included).
   integer(I4P)                       :: nk = 0_I4P          !< Points along z (halo included).
   integer(I4P)                       :: ng = 0_I4P          !< Halo layers on each side.
   real(R8P),             allocatable :: x(:)                !< x coordinates, strictly increasing.
   real(R8P),             allocatable :: y(:)                !< y coordinates, strictly increasing.
   real(R8P),             allocatable :: z(:)                !< z coordinates, strictly increasing.
   real(R8P)                          :: band_width = 0._R8P !< Band half-width used by the last compute_band.
   integer(I4P)                       :: band_cells = 0_I4P  !< Number of band cells.
   logical                            :: is_signed = .true.  !< Signed field (else unsigned distance).
   real(R8P),    allocatable, public  :: phi(:,:,:)          !< Signed distance; FAR_VALUE where unknown. PUBLIC for halo exchange.
   logical,      allocatable, public  :: in_band(:,:,:)      !< True on band cells (exact distance). PUBLIC, read-only.
   contains
      procedure, pass(self) :: initialize     !< Set the grid; reset the field.
      procedure, pass(self) :: destroy        !< Free the field.
      procedure, pass(self) :: compute_band   !< Exact signed distances in the narrow band.
      procedure, pass(self) :: get_dims       !< Points per axis and halo width.
      procedure, pass(self) :: get_band_width !< Band half-width of the last compute_band.
      procedure, pass(self) :: get_band_cells !< Number of band cells.
      procedure, pass(self) :: default_band   !< Default band for a surface: 2 x the local max spacing.
endtype distance_field_object

contains
   subroutine initialize(self, x, y, z, ng, status)
   !< Set the grid and reset the field (every interior and halo cell at FAR_VALUE, no band).
   !<
   !< `status = STATUS_INVALID_INPUT` (field left empty) if an axis has no point, is not
   !< strictly increasing, or leaves no interior point after removing `2 * ng` halo layers.
   class(distance_field_object), intent(inout)         :: self   !< Distance field.
   real(R8P),                    intent(in)            :: x(:)   !< x coordinates (halo included).
   real(R8P),                    intent(in)            :: y(:)   !< y coordinates (halo included).
   real(R8P),                    intent(in)            :: z(:)   !< z coordinates (halo included).
   integer(I4P),                 intent(in),  optional :: ng     !< Halo layers per side (default 0).
   integer(I4P),                 intent(out), optional :: status !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                                        :: ng_    !< Resolved halo width.

   call self%destroy
   if (present(status)) status = STATUS_OK
   ng_ = 0_I4P ; if (present(ng)) ng_ = ng
   if (ng_ < 0_I4P .or. .not. (is_valid_axis(x, ng_) .and. is_valid_axis(y, ng_) .and. is_valid_axis(z, ng_))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   self%ni = size(x, kind=I4P)
   self%nj = size(y, kind=I4P)
   self%nk = size(z, kind=I4P)
   self%ng = ng_
   allocate(self%x, source=x)
   allocate(self%y, source=y)
   allocate(self%z, source=z)
   allocate(self%phi(self%ni, self%nj, self%nk))
   allocate(self%in_band(self%ni, self%nj, self%nk))
   self%phi     = FAR_VALUE
   self%in_band = .false.
   endsubroutine initialize

   pure subroutine destroy(self)
   !< Free the field.
   class(distance_field_object), intent(inout) :: self  !< Distance field.

   if (allocated(self%x))       deallocate(self%x)
   if (allocated(self%y))       deallocate(self%y)
   if (allocated(self%z))       deallocate(self%z)
   if (allocated(self%phi))     deallocate(self%phi)
   if (allocated(self%in_band)) deallocate(self%in_band)
   self%ni = 0_I4P ; self%nj = 0_I4P ; self%nk = 0_I4P ; self%ng = 0_I4P
   self%band_width = 0._R8P
   self%band_cells = 0_I4P
   self%is_signed  = .true.
   endsubroutine destroy

   pure subroutine get_dims(self, ni, nj, nk, ng)
   !< Points per axis (halo included) and halo width.
   class(distance_field_object), intent(in)            :: self       !< Distance field.
   integer(I4P),                 intent(out)           :: ni, nj, nk !< Points per axis.
   integer(I4P),                 intent(out), optional :: ng         !< Halo layers per side.

   ni = self%ni ; nj = self%nj ; nk = self%nk
   if (present(ng)) ng = self%ng
   endsubroutine get_dims

   pure function get_band_width(self) result(band)
   !< Band half-width used by the last `compute_band` (0 before it).
   class(distance_field_object), intent(in) :: self !< Distance field.
   real(R8P)                                :: band !< Band half-width.

   band = self%band_width
   endfunction get_band_width

   pure function get_band_cells(self) result(n)
   !< Number of band cells after the last `compute_band`.
   class(distance_field_object), intent(in) :: self !< Distance field.
   integer(I4P)                             :: n    !< Band cells.

   n = self%band_cells
   endfunction get_band_cells

   function default_band(self, surface) result(band)
   !< Default band half-width: twice the largest grid spacing over the index box that covers
   !< the surface bounding box (one extra point per side). Two spacings guarantee that the two
   !< ends of any grid edge crossing the surface are band cells, which is what the sign flood
   !< fill needs. Outside the surface the grid may be stretched freely without widening the
   !< band. If the grid does not overlap the surface, the largest spacing of the whole grid.
   class(distance_field_object), intent(in) :: self    !< Distance field.
   type(surface_stl_object),     intent(in) :: surface !< Surface.
   real(R8P)                                :: band    !< Band half-width.
   type(vector_R8P)                         :: bmin    !< Surface bbox minimum.
   type(vector_R8P)                         :: bmax    !< Surface bbox maximum.

   bmin = surface%get_bmin()
   bmax = surface%get_bmax()
   band = 2._R8P * max(max_spacing(self%x, bmin%x, bmax%x), &
                       max_spacing(self%y, bmin%y, bmax%y), &
                       max_spacing(self%z, bmin%z, bmax%z))
   endfunction default_band

   subroutine compute_band(self, surface, band, is_signed, sign_algorithm, status)
   !< Exact (signed) distances on the narrow band: the interior cells with |d| <= `band`.
   !<
   !< 1. Candidates: each facet's bounding box, inflated by `band`, is mapped to an index box
   !<    by binary search on the coordinates; the union of those boxes is the candidate set.
   !<    Every cell within `band` of the surface is a candidate (it lies within `band` of some
   !<    facet, hence inside that facet's inflated box).
   !< 2. One `surface%distance_many` batch over the candidates (OpenMP-parallel when built
   !<    with OpenMP), with the requested sign algorithm.
   !< 3. Candidates with |d| <= band become band cells: `phi = d`, `in_band = .true.`.
   !<    Every other cell is reset to FAR_VALUE.
   !<
   !< The band values are bit-exact against a per-cell `surface%distance` call with the same
   !< options. The surface must be analysed (`load_from_file`, `sanitize`, `adopt_facets`).
   !<
   !< `status = STATUS_INVALID_INPUT` if the field is not initialized, `band <= 0`, or
   !< `distance_many` refuses (it never does for a consistent call).
   class(distance_field_object), intent(inout)         :: self           !< Distance field.
   type(surface_stl_object),     intent(inout), target :: surface        !< Surface (inout/target: device path).
   real(R8P),                    intent(in),  optional :: band           !< Band half-width (default `default_band`).
   logical,                      intent(in),  optional :: is_signed      !< Signed distances (default .true.).
   integer(I4P),                 intent(in),  optional :: sign_algorithm !< SIGN_* code (default SIGN_PSEUDO_NORMAL).
   integer(I4P),                 intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                           :: band_          !< Resolved band.
   integer(I4P)                                        :: algo           !< Resolved sign algorithm.
   logical,                      allocatable           :: candidate(:,:,:) !< Candidate mask.
   type(vector_R8P),             allocatable           :: points(:)      !< Candidate points.
   real(R8P),                    allocatable           :: d(:)           !< Candidate distances.
   integer(I4P),                 allocatable           :: cell(:,:)      !< Candidate indices (3, n).
   type(facet_object),           pointer               :: f              !< Facet cursor.
   integer(I4P)                                        :: i, j, k, m, n, st
   integer(I4P)                                        :: ilo, ihi, jlo, jhi, klo, khi
   real(R8P)                                           :: lo(3), hi(3)   !< Inflated facet bbox.

   if (present(status)) status = STATUS_OK
   if (.not. allocated(self%phi)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   band_ = 0._R8P
   if (present(band)) then
      band_ = band
   else
      band_ = self%default_band(surface)
   endif
   if (.not. (band_ > 0._R8P)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   self%is_signed = .true. ; if (present(is_signed)) self%is_signed = is_signed
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm

   self%phi        = FAR_VALUE
   self%in_band    = .false.
   self%band_width = band_
   self%band_cells = 0_I4P

   ! 1. candidates, restricted to the interior (halo cells are the caller's)
   allocate(candidate(self%ni, self%nj, self%nk))
   candidate = .false.
   do m = 1_I4P, surface%get_facets_number()
      f => surface%facet_at(m)
      lo(1) = min(f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x) - band_
      lo(2) = min(f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y) - band_
      lo(3) = min(f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z) - band_
      hi(1) = max(f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x) + band_
      hi(2) = max(f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y) + band_
      hi(3) = max(f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z) + band_
      call index_range(self%x, self%ng, lo(1), hi(1), ilo, ihi)
      call index_range(self%y, self%ng, lo(2), hi(2), jlo, jhi)
      call index_range(self%z, self%ng, lo(3), hi(3), klo, khi)
      if (ilo > ihi .or. jlo > jhi .or. klo > khi) cycle
      candidate(ilo:ihi, jlo:jhi, klo:khi) = .true.
   enddo

   ! 2. one batch over the candidates
   n = count(candidate, kind=I4P)
   if (n == 0_I4P) return
   allocate(points(n), d(n), cell(3, n))
   m = 0_I4P
   do k = 1_I4P, self%nk
      do j = 1_I4P, self%nj
         do i = 1_I4P, self%ni
            if (.not. candidate(i, j, k)) cycle
            m = m + 1_I4P
            points(m)%x = self%x(i)
            points(m)%y = self%y(j)
            points(m)%z = self%z(k)
            cell(1, m) = i ; cell(2, m) = j ; cell(3, m) = k
         enddo
      enddo
   enddo
   call surface%distance_many(points=points, distances=d, is_signed=self%is_signed, sign_algorithm=algo, &
                              is_square_root=.true., status=st)
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif

   ! 3. keep |d| <= band
   do m = 1_I4P, n
      if (abs(d(m)) <= band_) then
         self%phi(cell(1, m), cell(2, m), cell(3, m))     = d(m)
         self%in_band(cell(1, m), cell(2, m), cell(3, m)) = .true.
         self%band_cells = self%band_cells + 1_I4P
      endif
   enddo
   endsubroutine compute_band

   ! private helpers
   pure function is_valid_axis(c, ng) result(ok)
   !< True if the axis has an interior point after `ng` halo layers per side and is strictly
   !< increasing.
   real(R8P),    intent(in) :: c(:) !< Coordinates.
   integer(I4P), intent(in) :: ng   !< Halo layers per side.
   logical                  :: ok   !< Validity.
   integer(I4P)             :: i    !< Counter.

   ok = size(c, kind=I4P) > 2_I4P * ng
   if (.not. ok) return
   do i = 2_I4P, size(c, kind=I4P)
      if (.not. (c(i) > c(i - 1_I4P))) then
         ok = .false.
         return
      endif
   enddo
   endfunction is_valid_axis

   pure subroutine index_range(c, ng, lo, hi, ilo, ihi)
   !< Interior index range [ilo, ihi] of the coordinates inside [lo, hi] (empty: ilo > ihi).
   !< Binary search, so the cost does not depend on the spacing.
   real(R8P),    intent(in)  :: c(:)     !< Strictly increasing coordinates.
   integer(I4P), intent(in)  :: ng       !< Halo layers per side.
   real(R8P),    intent(in)  :: lo, hi   !< Interval.
   integer(I4P), intent(out) :: ilo, ihi !< Index range.
   integer(I4P)              :: n        !< Points.

   n   = size(c, kind=I4P)
   ilo = max(first_at_least(c, lo), ng + 1_I4P)
   ihi = min(last_at_most(c, hi), n - ng)
   endsubroutine index_range

   pure function first_at_least(c, v) result(i)
   !< Smallest i with c(i) >= v (size(c) + 1 if none).
   real(R8P), intent(in) :: c(:)   !< Strictly increasing coordinates.
   real(R8P), intent(in) :: v      !< Value.
   integer(I4P)          :: i      !< Index.
   integer(I4P)          :: a, b, m

   a = 1_I4P ; b = size(c, kind=I4P) + 1_I4P
   do while (a < b)
      m = (a + b) / 2_I4P
      if (c(m) >= v) then
         b = m
      else
         a = m + 1_I4P
      endif
   enddo
   i = a
   endfunction first_at_least

   pure function last_at_most(c, v) result(i)
   !< Largest i with c(i) <= v (0 if none).
   real(R8P), intent(in) :: c(:)   !< Strictly increasing coordinates.
   real(R8P), intent(in) :: v      !< Value.
   integer(I4P)          :: i      !< Index.
   integer(I4P)          :: a, b, m

   a = 0_I4P ; b = size(c, kind=I4P)
   do while (a < b)
      m = (a + b + 1_I4P) / 2_I4P
      if (c(m) <= v) then
         a = m
      else
         b = m - 1_I4P
      endif
   enddo
   i = a
   endfunction last_at_most

   pure function max_spacing(c, lo, hi) result(h)
   !< Largest spacing c(i+1) - c(i) over the index range covering [lo, hi] plus one point per
   !< side; the largest spacing of the whole axis if the range is empty. 0 for one point.
   real(R8P), intent(in) :: c(:)   !< Strictly increasing coordinates.
   real(R8P), intent(in) :: lo, hi !< Interval.
   real(R8P)             :: h      !< Largest spacing.
   integer(I4P)          :: a, b, n, i

   n = size(c, kind=I4P)
   h = 0._R8P
   if (n < 2_I4P) return
   a = max(first_at_least(c, lo) - 1_I4P, 1_I4P)
   b = min(last_at_most(c, hi) + 1_I4P, n)
   if (a >= b) then
      a = 1_I4P ; b = n
   endif
   do i = a, b - 1_I4P
      h = max(h, c(i + 1_I4P) - c(i))
   enddo
   endfunction max_spacing
endmodule fossil_distance_field
