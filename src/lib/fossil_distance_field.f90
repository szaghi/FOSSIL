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
!<     signed distance, computed in one `surface%distance_many` batch, or on the GPU through
!<     `surface%distance_many_device` with `use_device=.true.`. Candidates come from
!<     rasterizing each facet's bounding box, inflated by `band`, onto the index space
!<     (binary search on the 1D coordinates, so stretched grids cost the same).
!<  2. **Signs** (`compute_signs`): the cells outside the band are grouped into 6-connected
!<     components; the band separates inside from outside, so each component takes the sign
!<     of the band cells it touches. Point queries cover what the band cannot decide.
!<  3. **Far field** (`sweep`): the Fast Sweeping Method (Zhao 2005) solves |grad phi| = 1
!<     outside the band, with the band values fixed. Godunov upwind update with the
!<     non-uniform one-sided spacings of each axis; sweeps in hyperplane order (Detrixhe,
!<     Gibou & Min 2013), so the cells of one plane are independent and are updated in
!<     parallel (OpenMP) with the serial result bit for bit. Far-field error O(h |log h|).
!<  `compute` runs the three stages.
!<
!< Halo: the outer `ng` layers of cells on each side are halo cells. They are never
!< computed; a caller with a block decomposition fills them from the neighbouring blocks.
!< `phi` is public for that reason. FOSSIL itself stays MPI-free.
!<
!< Sign convention: negative inside, positive outside, as `surface%distance`.

! The surface and facet modules, not `fossil`: `fossil` re-exports this module, so it is built
! into the library.
use fossil_surface_stl_object, only : surface_stl_object, SIGN_PSEUDO_NORMAL, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object, only : facet_object
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
   real(R8P)                          :: local_spacing = 0._R8P !< Largest spacing near the surface (default_band / 2).
   integer(I4P)                       :: sign_algorithm = SIGN_PSEUDO_NORMAL !< Sign algorithm of the last compute_band.
   logical                            :: use_device = .false. !< Device flag of the last compute_band (point queries follow it).
   integer(I4P)                       :: band_cells = 0_I4P  !< Number of band cells.
   logical                            :: is_signed = .true.  !< Signed field (else unsigned distance).
   real(R8P),    allocatable, public  :: phi(:,:,:)          !< Signed distance; FAR_VALUE where unknown. PUBLIC for halo exchange.
   logical,      allocatable, public  :: in_band(:,:,:)      !< True on band cells (exact distance). PUBLIC, read-only.
   contains
      procedure, pass(self) :: initialize     !< Set the grid; reset the field.
      procedure, pass(self) :: destroy        !< Free the field.
      procedure, pass(self) :: compute_band   !< Exact signed distances in the narrow band.
      procedure, pass(self) :: compute_signs  !< Signs outside the band by flood fill.
      procedure, pass(self) :: sweep          !< Far-field magnitudes by fast sweeping.
      procedure, pass(self) :: compute        !< compute_band + compute_signs + sweep.
      procedure, pass(self) :: get_dims       !< Points per axis and halo width.
      procedure, pass(self) :: get_band_width !< Band half-width of the last compute_band.
      procedure, pass(self) :: get_band_cells !< Number of band cells.
      procedure, pass(self) :: default_band   !< Default band for a surface: 2 x the local max spacing.
      procedure, pass(self), private :: linear   !< Linear index of a cell.
      procedure, pass(self), private :: sweep_once !< One directional sweep.
      procedure, pass(self), private :: godunov  !< Godunov update of one cell.
      procedure, pass(self), private :: unlinear !< Cell of a linear index.
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
   self%local_spacing = 0._R8P
   self%sign_algorithm = SIGN_PSEUDO_NORMAL
   self%use_device = .false.
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

   subroutine compute_band(self, surface, band, is_signed, sign_algorithm, use_device, status)
   !< Exact (signed) distances on the narrow band: the interior cells with |d| <= `band`.
   !<
   !< 1. Candidates: each facet's bounding box, inflated by `band`, is mapped to an index box
   !<    by binary search on the coordinates; the union of those boxes is the candidate set.
   !<    Every cell within `band` of the surface is a candidate (it lies within `band` of some
   !<    facet, hence inside that facet's inflated box).
   !< 2. One `surface%distance_many` batch over the candidates (OpenMP-parallel when built
   !<    with OpenMP), with the requested sign algorithm. With `use_device=.true.` the batch
   !<    runs through `surface%distance_many_device` instead (bit-exact with the host batch
   !<    where validated); the surface is made resident for the call and released after it,
   !<    unless it was already resident, in which case residency is left to the caller.
   !<    The device supports SIGN_PSEUDO_NORMAL and SIGN_FAST_WINDING_NUMBER only, on the SAH
   !<    BVH; anything else is refused (STATUS_INVALID_INPUT, stderr message, band empty).
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
   logical,                      intent(in),  optional :: use_device     !< Band batch on the GPU (default .false.).
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
   self%sign_algorithm = algo
   self%use_device = .false. ; if (present(use_device)) self%use_device = use_device
   self%local_spacing  = 0.5_R8P * self%default_band(surface)

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
   call batch_distances(surface=surface, points=points, d=d, is_signed=self%is_signed, sign_algorithm=algo, &
                        use_device=self%use_device, status=st)
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

   subroutine compute_signs(self, surface, fallback_cells, status)
   !< Signs of the cells outside the band (call after `compute_band`, before `sweep`).
   !<
   !< The non-band interior cells are split into 6-connected components (halo cells take no
   !< part). When the band is at least as thick as the largest grid spacing near the surface,
   !< both ends of every grid edge that crosses the surface are band cells, so no component
   !< contains both inside and outside cells. Each component then takes:
   !<  - the sign of the band cells it touches, if they all agree (cells exactly on the
   !<    surface, d = 0, carry no sign);
   !<  - one point query, at its first cell, if it touches no band cell (a block that does not
   !<    see the surface, or a cavity thinner than the band);
   !<  - a point query per cell if it touches band cells of both signs, which can only happen
   !<    with an inconsistent band (defensive).
   !< If the band is thinner than the largest spacing near the surface, a component can leak
   !< across the surface while touching band cells of one sign only, so every non-band cell
   !< gets a point query instead: correct, at the cost of one distance query per cell.
   !<
   !< Point queries use `surface%distance_many` with the sign algorithm of `compute_band`, or
   !< `surface%distance_many_device` if `compute_band` ran with `use_device=.true.` (bit-exact
   !< with the host queries where validated).
   !< Magnitudes are kept (`phi = sign(|phi|, s)`). No-op on an unsigned field.
   class(distance_field_object), intent(inout)         :: self           !< Distance field.
   type(surface_stl_object),     intent(inout), target :: surface        !< Surface.
   integer(I4P),                 intent(out), optional :: fallback_cells !< Cells signed by a point query.
   integer(I4P),                 intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),                 allocatable           :: label(:,:,:)   !< Component id; 0 unvisited, -1 band/halo.
   integer(I4P),                 allocatable           :: queue(:)       !< BFS queue (linear indices).
   integer(I4P),                 allocatable           :: comp_sign(:)   !< Per component: +1, -1, 0 (query each cell), 2 (query one).
   integer(I4P),                 allocatable           :: comp_first(:)  !< Per component: first cell (linear index).
   type(vector_R8P),             allocatable           :: points(:)      !< Query points.
   real(R8P),                    allocatable           :: d(:)           !< Query results.
   integer(I4P),                 allocatable           :: qcell(:)       !< Query cells (linear indices) or components.
   integer(I4P)                                        :: i, j, k, l, m, c, n_comp, head, tail, nb
   integer(I4P)                                        :: ii, jj, kk, n_query, st
   logical                                             :: touch_pos, touch_neg, thin
   integer(I4P), parameter                             :: OFF(3, 6) = reshape([-1, 0, 0,  1, 0, 0,  0, -1, 0, &
                                                                               0, 1, 0,  0, 0, -1,  0, 0, 1], [3, 6])

   if (present(status)) status = STATUS_OK
   if (present(fallback_cells)) fallback_cells = 0_I4P
   if (.not. allocated(self%phi)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (.not. self%is_signed) return
   thin = self%band_width < self%local_spacing

   allocate(label(self%ni, self%nj, self%nk))
   label = -1_I4P
   label(self%ng + 1:self%ni - self%ng, self%ng + 1:self%nj - self%ng, self%ng + 1:self%nk - self%ng) = 0_I4P
   where (self%in_band) label = -1_I4P
   allocate(queue(count(label == 0_I4P, kind=I4P)))
   allocate(comp_sign(size(queue)), comp_first(size(queue)))

   ! label the components, recording the band signs each one touches
   n_comp = 0_I4P
   do k = 1_I4P, self%nk
      do j = 1_I4P, self%nj
         do i = 1_I4P, self%ni
            if (label(i, j, k) /= 0_I4P) cycle
            n_comp = n_comp + 1_I4P
            label(i, j, k) = n_comp
            comp_first(n_comp) = self%linear(i, j, k)
            touch_pos = .false. ; touch_neg = .false.
            head = 1_I4P ; tail = 1_I4P ; queue(1) = comp_first(n_comp)
            do while (head <= tail)
               call self%unlinear(queue(head), ii, jj, kk)
               head = head + 1_I4P
               do nb = 1, 6
                  l = ii + OFF(1, nb) ; m = jj + OFF(2, nb) ; c = kk + OFF(3, nb)
                  if (l < 1_I4P .or. l > self%ni .or. m < 1_I4P .or. m > self%nj .or. c < 1_I4P .or. c > self%nk) cycle
                  if (self%in_band(l, m, c)) then
                     if (self%phi(l, m, c) > 0._R8P) touch_pos = .true.
                     if (self%phi(l, m, c) < 0._R8P) touch_neg = .true.
                  elseif (label(l, m, c) == 0_I4P) then
                     label(l, m, c) = n_comp
                     tail = tail + 1_I4P
                     queue(tail) = self%linear(l, m, c)
                  endif
               enddo
            enddo
            if (thin .or. (touch_pos .and. touch_neg)) then
               comp_sign(n_comp) = 0_I4P
            elseif (touch_pos) then
               comp_sign(n_comp) = 1_I4P
            elseif (touch_neg) then
               comp_sign(n_comp) = -1_I4P
            else
               comp_sign(n_comp) = 2_I4P
            endif
         enddo
      enddo
   enddo

   ! point queries: every cell of a 0-component, the first cell of a 2-component
   n_query = count(comp_sign(1:n_comp) == 2_I4P, kind=I4P)
   do k = 1_I4P, self%nk
      do j = 1_I4P, self%nj
         do i = 1_I4P, self%ni
            if (label(i, j, k) > 0_I4P) then
               if (comp_sign(label(i, j, k)) == 0_I4P) n_query = n_query + 1_I4P
            endif
         enddo
      enddo
   enddo
   allocate(points(n_query), d(n_query), qcell(n_query))
   m = 0_I4P
   do c = 1_I4P, n_comp
      if (comp_sign(c) /= 2_I4P) cycle
      m = m + 1_I4P
      qcell(m) = comp_first(c)
      call self%unlinear(comp_first(c), ii, jj, kk)
      points(m)%x = self%x(ii) ; points(m)%y = self%y(jj) ; points(m)%z = self%z(kk)
   enddo
   do k = 1_I4P, self%nk
      do j = 1_I4P, self%nj
         do i = 1_I4P, self%ni
            if (label(i, j, k) <= 0_I4P) cycle
            if (comp_sign(label(i, j, k)) /= 0_I4P) cycle
            m = m + 1_I4P
            qcell(m) = self%linear(i, j, k)
            points(m)%x = self%x(i) ; points(m)%y = self%y(j) ; points(m)%z = self%z(k)
         enddo
      enddo
   enddo
   if (n_query > 0_I4P) then
      call batch_distances(surface=surface, points=points, d=d, is_signed=.true., sign_algorithm=self%sign_algorithm, &
                           use_device=self%use_device, status=st)
      if (st /= STATUS_OK) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   if (present(fallback_cells)) fallback_cells = n_query

   ! component-wide signs from the one-point queries
   do m = 1_I4P, n_query
      call self%unlinear(qcell(m), ii, jj, kk)
      c = label(ii, jj, kk)
      if (comp_sign(c) == 2_I4P) then
         comp_sign(c) = 1_I4P ; if (d(m) < 0._R8P) comp_sign(c) = -1_I4P
      else
         self%phi(ii, jj, kk) = sign(abs(self%phi(ii, jj, kk)), d(m))
      endif
   enddo
   ! apply
   do k = 1_I4P, self%nk
      do j = 1_I4P, self%nj
         do i = 1_I4P, self%ni
            if (label(i, j, k) <= 0_I4P) cycle
            c = comp_sign(label(i, j, k))
            if (c == 1_I4P)  self%phi(i, j, k) =  abs(self%phi(i, j, k))
            if (c == -1_I4P) self%phi(i, j, k) = -abs(self%phi(i, j, k))
         enddo
      enddo
   enddo
   endsubroutine compute_signs

   subroutine sweep(self, tol, max_sweeps, sweeps, max_change, status)
   !< Far-field magnitudes by the Fast Sweeping Method (call after `compute_band` and, for a
   !< signed field, `compute_signs`).
   !<
   !< Solves |grad phi| = 1 on the non-band interior cells with the band cells and the halo
   !< cells fixed. Each directional sweep visits the hyperplanes i' + j' + k' = s in order,
   !< where (i', j', k') are the interior indices counted from the corner the sweep starts at;
   !< a cell's six neighbours lie on the planes s - 1 and s + 1, so the cells of one plane are
   !< independent and the result does not depend on the order inside a plane. The 8 sweep
   !< directions alternate; after each full cycle of 8 the largest change of |phi| is
   !< compared with `tol` times the diagonal of the interior box.
   !<
   !< The update of a cell takes, per axis, the neighbour that gives the larger one-sided slope
   !< (u - a) / h at the solution (h is that side's spacing), and solves the Godunov quadratic
   !< over the axes whose neighbour is below u. Magnitudes only decrease; signs are kept.
   !<
   !< Accuracy: first order away from the band (the band values are exact). A caller with a
   !< block decomposition alternates halo exchanges with calls of `sweep(max_sweeps=8)` until
   !< the global `max_change` is below its tolerance.
   class(distance_field_object), intent(inout)         :: self       !< Distance field.
   real(R8P),                    intent(in),  optional :: tol        !< Relative tolerance (default 1e-12).
   integer(I4P),                 intent(in),  optional :: max_sweeps !< Maximum directional sweeps (default 400).
   integer(I4P),                 intent(out), optional :: sweeps     !< Directional sweeps done.
   real(R8P),                    intent(out), optional :: max_change !< Largest |phi| change of the last cycle.
   integer(I4P),                 intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                           :: tol_       !< Absolute tolerance.
   integer(I4P)                                        :: max_sweeps_, n_sweeps, dir
   real(R8P)                                           :: cycle_change, change

   if (present(status)) status = STATUS_OK
   if (present(sweeps)) sweeps = 0_I4P
   if (present(max_change)) max_change = 0._R8P
   if (.not. allocated(self%phi)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   tol_ = 1.e-12_R8P ; if (present(tol)) tol_ = tol
   tol_ = tol_ * sqrt((self%x(self%ni - self%ng) - self%x(self%ng + 1_I4P))**2 + &
                      (self%y(self%nj - self%ng) - self%y(self%ng + 1_I4P))**2 + &
                      (self%z(self%nk - self%ng) - self%z(self%ng + 1_I4P))**2)
   max_sweeps_ = 400_I4P ; if (present(max_sweeps)) max_sweeps_ = max_sweeps

   n_sweeps     = 0_I4P
   cycle_change = 0._R8P
   do while (n_sweeps < max_sweeps_)
      dir = mod(n_sweeps, 8_I4P)
      call self%sweep_once(dir=dir, change=change)
      n_sweeps = n_sweeps + 1_I4P
      cycle_change = max(cycle_change, change)
      if (dir == 7_I4P) then
         if (cycle_change <= tol_) exit
         if (n_sweeps < max_sweeps_) cycle_change = 0._R8P
      endif
   enddo
   if (present(sweeps)) sweeps = n_sweeps
   if (present(max_change)) max_change = cycle_change
   endsubroutine sweep

   subroutine sweep_once(self, dir, change)
   !< One directional sweep in hyperplane order. `dir` in 0..7 selects the signs of the three
   !< axes (bit 0: x reversed, bit 1: y, bit 2: z). `change` is the largest |phi| decrease
   !< (FAR_VALUE counts as a change of FAR_VALUE).
   class(distance_field_object), intent(inout) :: self   !< Distance field.
   integer(I4P),                 intent(in)    :: dir    !< Sweep direction.
   real(R8P),                    intent(out)   :: change !< Largest change.
   integer(I4P)                                :: n1, n2, n3, s, ip, jp, kp, i, j, k
   logical                                     :: rx, ry, rz
   real(R8P)                                   :: u, old

   n1 = self%ni - 2_I4P * self%ng
   n2 = self%nj - 2_I4P * self%ng
   n3 = self%nk - 2_I4P * self%ng
   rx = btest(dir, 0) ; ry = btest(dir, 1) ; rz = btest(dir, 2)
   change = 0._R8P
   ! One parallel region per sweep, one worksharing loop per plane (its implicit barrier
   ! orders the planes). The cells of a plane are independent, so any thread split gives the
   ! serial result bit for bit; max is an exact reduction.
   !$omp parallel default(none) shared(self, n1, n2, n3, rx, ry, rz) private(s, kp, jp, ip, i, j, k, old, u) &
   !$omp   reduction(max:change)
   do s = 0_I4P, n1 + n2 + n3 - 3_I4P
      !$omp do schedule(static)
      do kp = max(0_I4P, s - (n1 - 1_I4P) - (n2 - 1_I4P)), min(n3 - 1_I4P, s)
         do jp = max(0_I4P, s - kp - (n1 - 1_I4P)), min(n2 - 1_I4P, s - kp)
            ip = s - kp - jp
            i = self%ng + 1_I4P + ip ; if (rx) i = self%ni - self%ng - ip
            j = self%ng + 1_I4P + jp ; if (ry) j = self%nj - self%ng - jp
            k = self%ng + 1_I4P + kp ; if (rz) k = self%nk - self%ng - kp
            if (self%in_band(i, j, k)) cycle
            old = abs(self%phi(i, j, k))
            u   = self%godunov(i, j, k)
            if (u < old) then
               self%phi(i, j, k) = sign(u, self%phi(i, j, k))
               change = max(change, old - u)
            endif
         enddo
      enddo
      !$omp end do
   enddo
   !$omp end parallel
   endsubroutine sweep_once

   pure function godunov(self, i, j, k) result(u)
   !< Godunov upwind solution of |grad u| = 1 at cell (i, j, k) from the magnitudes of its six
   !< neighbours, with non-uniform one-sided spacings. Missing neighbours (grid boundary) and
   !< FAR_VALUE neighbours do not contribute. Returns FAR_VALUE when no neighbour is known.
   class(distance_field_object), intent(in) :: self    !< Distance field.
   integer(I4P),                 intent(in) :: i, j, k !< Cell.
   real(R8P)                                :: u       !< Solution.
   real(R8P)                                :: a(2, 3) !< Neighbour magnitudes (side, axis).
   real(R8P)                                :: h(2, 3) !< One-sided spacings (side, axis).
   real(R8P)                                :: aa(3), hh(3) !< Chosen side per axis.
   integer(I4P)                             :: side(3), side_new(3), d, it
   real(R8P)                                :: ub      !< One-axis upper bound.

   a = FAR_VALUE ; h = 1._R8P
   if (i > 1_I4P)       then ; a(1, 1) = abs(self%phi(i - 1_I4P, j, k)) ; h(1, 1) = self%x(i) - self%x(i - 1_I4P) ; endif
   if (i < self%ni)     then ; a(2, 1) = abs(self%phi(i + 1_I4P, j, k)) ; h(2, 1) = self%x(i + 1_I4P) - self%x(i) ; endif
   if (j > 1_I4P)       then ; a(1, 2) = abs(self%phi(i, j - 1_I4P, k)) ; h(1, 2) = self%y(j) - self%y(j - 1_I4P) ; endif
   if (j < self%nj)     then ; a(2, 2) = abs(self%phi(i, j + 1_I4P, k)) ; h(2, 2) = self%y(j + 1_I4P) - self%y(j) ; endif
   if (k > 1_I4P)       then ; a(1, 3) = abs(self%phi(i, j, k - 1_I4P)) ; h(1, 3) = self%z(k) - self%z(k - 1_I4P) ; endif
   if (k < self%nk)     then ; a(2, 3) = abs(self%phi(i, j, k + 1_I4P)) ; h(2, 3) = self%z(k + 1_I4P) - self%z(k) ; endif

   ! upper bound: the best one-axis solution
   u = FAR_VALUE
   do d = 1, 3
      if (a(1, d) < FAR_VALUE) u = min(u, a(1, d) + h(1, d))
      if (a(2, d) < FAR_VALUE) u = min(u, a(2, d) + h(2, d))
   enddo
   if (u >= FAR_VALUE) return
   ub = u
   ! choose per axis the side with the larger slope (u - a) / h at the current u, solve, repeat
   side = 0_I4P
   do it = 1, 4
      do d = 1, 3
         ! a FAR_VALUE side never wins; no division by h of a FAR_VALUE difference (overflow)
         if (a(2, d) >= FAR_VALUE) then
            side_new(d) = 1_I4P
         elseif (a(1, d) >= FAR_VALUE) then
            side_new(d) = 2_I4P
         elseif ((u - a(1, d)) / h(1, d) >= (u - a(2, d)) / h(2, d)) then
            side_new(d) = 1_I4P
         else
            side_new(d) = 2_I4P
         endif
         aa(d) = a(side_new(d), d) ; hh(d) = h(side_new(d), d)
      enddo
      if (all(side_new == side)) exit
      side = side_new
      u = solve_quadratic(aa, hh)
   enddo
   u = min(u, ub)   ! the multi-axis solution never exceeds the best one-axis one
   endfunction godunov

   pure function linear(self, i, j, k) result(l)
   !< Linear index of cell (i, j, k).
   class(distance_field_object), intent(in) :: self    !< Distance field.
   integer(I4P),                 intent(in) :: i, j, k !< Cell.
   integer(I4P)                             :: l       !< Linear index.

   l = i + self%ni * ((j - 1_I4P) + self%nj * (k - 1_I4P))
   endfunction linear

   pure subroutine unlinear(self, l, i, j, k)
   !< Cell (i, j, k) of linear index l.
   class(distance_field_object), intent(in)  :: self    !< Distance field.
   integer(I4P),                 intent(in)  :: l       !< Linear index.
   integer(I4P),                 intent(out) :: i, j, k !< Cell.
   integer(I4P)                              :: r       !< Remainder.

   r = l - 1_I4P
   i = mod(r, self%ni) + 1_I4P
   r = r / self%ni
   j = mod(r, self%nj) + 1_I4P
   k = r / self%nj + 1_I4P
   endsubroutine unlinear

   subroutine compute(self, surface, band, is_signed, sign_algorithm, use_device, tol, max_sweeps, sweeps, fallback_cells, &
                      status)
   !< The whole field: `compute_band`, `compute_signs` (signed fields), `sweep`. See each one
   !< for the arguments; `sweeps` and `fallback_cells` report their counters.
   class(distance_field_object), intent(inout)         :: self           !< Distance field.
   type(surface_stl_object),     intent(inout), target :: surface        !< Surface.
   real(R8P),                    intent(in),  optional :: band           !< Band half-width (default `default_band`).
   logical,                      intent(in),  optional :: is_signed      !< Signed distances (default .true.).
   integer(I4P),                 intent(in),  optional :: sign_algorithm !< SIGN_* code (default SIGN_PSEUDO_NORMAL).
   logical,                      intent(in),  optional :: use_device     !< Band batch on the GPU (default .false.).
   real(R8P),                    intent(in),  optional :: tol            !< Sweep tolerance (default 1e-12, relative).
   integer(I4P),                 intent(in),  optional :: max_sweeps     !< Maximum directional sweeps (default 400).
   integer(I4P),                 intent(out), optional :: sweeps         !< Directional sweeps done.
   integer(I4P),                 intent(out), optional :: fallback_cells !< Cells signed by a point query.
   integer(I4P),                 intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                                        :: st             !< Stage status.

   if (present(status)) status = STATUS_OK
   call self%compute_band(surface=surface, band=band, is_signed=is_signed, sign_algorithm=sign_algorithm, &
                          use_device=use_device, status=st)
   if (st == STATUS_OK) call self%compute_signs(surface=surface, fallback_cells=fallback_cells, status=st)
   if (st == STATUS_OK) call self%sweep(tol=tol, max_sweeps=max_sweeps, sweeps=sweeps, status=st)
   if (present(status)) status = st
   endsubroutine compute

   ! private helpers
   subroutine batch_distances(surface, points, d, is_signed, sign_algorithm, use_device, status)
   !< One batch of distance queries, on the host (`distance_many`) or on the device
   !< (`distance_many_device`). The device path makes the surface resident for the call and
   !< releases it after, unless it was already resident.
   type(surface_stl_object), intent(inout), target :: surface        !< Surface.
   type(vector_R8P),         intent(in)            :: points(:)      !< Query points.
   real(R8P),                intent(out)           :: d(:)           !< Distances.
   logical,                  intent(in)            :: is_signed      !< Signed distances.
   integer(I4P),             intent(in)            :: sign_algorithm !< SIGN_* code.
   logical,                  intent(in)            :: use_device     !< Device batch.
   integer(I4P),             intent(out)           :: status         !< distance_many(_device) status.
   logical                                         :: was_resident   !< Surface already on the device.

   if (use_device) then
      was_resident = surface%is_on_device
      if (.not. was_resident) call surface%enter_device
      call surface%distance_many_device(points=points, distances=d, is_signed=is_signed, sign_algorithm=sign_algorithm, &
                                        is_square_root=.true., status=status)
      if (.not. was_resident) call surface%exit_device
   else
      call surface%distance_many(points=points, distances=d, is_signed=is_signed, sign_algorithm=sign_algorithm, &
                                 is_square_root=.true., status=status)
   endif
   endsubroutine batch_distances

   pure function solve_quadratic(a, h) result(u)
   !< Solution of sum over the active axes of ((u - a(d)) / h(d))**2 = 1, where an axis is
   !< active when u > a(d): axes are added in increasing a(d) and the first consistent solution
   !< is returned. The one-axis solution min(a + h) is always consistent.
   real(R8P), intent(in) :: a(3)  !< Neighbour magnitudes, one per axis.
   real(R8P), intent(in) :: h(3)  !< Matching spacings.
   real(R8P)             :: u     !< Solution.
   real(R8P)             :: as(3), hs(3), t, qa, qb, qc, disc
   integer(I4P)          :: m, p, q

   as = a ; hs = h
   do p = 1, 2   ! sort by a (3 elements)
      do q = 1, 3 - p
         if (as(q) > as(q + 1)) then
            t = as(q) ; as(q) = as(q + 1) ; as(q + 1) = t
            t = hs(q) ; hs(q) = hs(q + 1) ; hs(q + 1) = t
         endif
      enddo
   enddo
   u = as(1) + hs(1)
   qa = 0._R8P ; qb = 0._R8P ; qc = -1._R8P
   do m = 1, 3
      if (as(m) >= FAR_VALUE) exit
      qa = qa + 1._R8P / hs(m)**2
      qb = qb + as(m) / hs(m)**2
      qc = qc + (as(m) / hs(m))**2
      disc = qb * qb - qa * qc
      if (disc < 0._R8P) exit
      t = (qb + sqrt(disc)) / qa
      if (t < as(m)) exit
      u = t
      if (m < 3) then
         if (t <= as(m + 1)) exit
      endif
   enddo
   endfunction solve_quadratic

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
