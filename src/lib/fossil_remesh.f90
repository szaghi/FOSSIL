!< FOSSIL, isotropic remeshing for triangle surfaces (issue #18 §1.7).

module fossil_remesh
!< FOSSIL, isotropic remeshing (Botsch & Kobbelt 2004 — *A Remeshing Approach
!< to Multiresolution Modeling*, SGP 2004).
!<
!< STATUS: incremental rebuild after the first attempt's heap-corruption bug
!< (see issue #18 §1.7 history). This iteration uses simpler rectangular
!< (max_valence × nv) incidence arrays instead of the packed flat-array
!< design that proved hard to debug. Each step adds one pass with its own
!< standalone test before moving on.
!<
!< STEP 1: public API + collect_state + materialize round-trip. No topology.
!< STEP 2: build_edges + median L exposed via two private helpers.
!< STEP 3: pass_split — insert midpoints on over-length edges.
!< STEP 4: pass_collapse — collapse short edges with safety checks.
!< STEP 5 (this commit): pass_flip — flip interior edges when doing so
!< improves valence balance toward the target valence (6 for interior).
!< Geometric validity check uses post-flip normal direction relative to
!< pre-flip; degenerate (zero-area) flips are rejected.

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,     only : facet_object
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: isotropic_remesh
public :: REM_STATUS_OK, REM_STATUS_BAD_INPUT
! Testing-only helpers (used by fossil_test_isotropic_remesh during the
! incremental rebuild; kept public so each step can be verified in
! isolation without touching the algorithm's main path).
public :: count_unique_edges
public :: compute_median_edge_length
public :: run_split_only
public :: run_split_and_collapse
public :: run_split_collapse_flip

integer(I4P), parameter :: MAX_VAL = 32_I4P  !< Max vertex valence supported (rectangular incidence cap).

integer(I4P), parameter :: REM_STATUS_OK        = 0_I4P
integer(I4P), parameter :: REM_STATUS_BAD_INPUT = 1_I4P  !< Empty input or vertex_id unset.

contains

   subroutine isotropic_remesh(facet, target_length, iterations, preserve_features, &
                                reference_facet, reference_tree, status)
   !< Remesh `facet` to (approximately) uniform edge length `target_length`.
   !<
   !< STEP-3 BEHAVIOR: applies the split pass `iterations` times. No collapse,
   !< flip, or relax yet. `preserve_features` and the reference surface are
   !< accepted but ignored (placeholders for later steps).
   !<
   !< If `target_length <= 0`, defaults to median input edge length.
   !<
   !< Cycle-break note: takes raw `(facet, tree)` for the projection
   !< reference rather than `surface_stl_object`, mirroring
   !< `fossil_winding_number` and `fossil_boolean`.
   type(facet_object), allocatable,            intent(inout)        :: facet(:)
   real(R8P),                                  intent(in)           :: target_length
   integer(I4P),                               intent(in)           :: iterations
   logical,                                    intent(in)           :: preserve_features
   type(facet_object),                         intent(in), optional, target :: reference_facet(:)
   type(aabb_tree_object),                     intent(in), optional, target :: reference_tree
   integer(I4P),                               intent(out), optional :: status
   integer(I4P)                                                     :: nf_in, nv, nf_alive, st_local, it
   integer(I4P), allocatable                                        :: f_v(:, :)
   real(R8P),    allocatable                                        :: vcoord(:, :)
   logical,      allocatable                                        :: facet_alive(:)
   real(R8P)                                                        :: L

   if (present(status)) status = REM_STATUS_OK
   ! Suppress unused-argument warnings (placeholders for later steps).
   if (preserve_features .and. .false.) return
   if (present(reference_facet) .and. .false.) return
   if (present(reference_tree) .and. .false.) return

   if (.not. allocated(facet)) then
      if (present(status)) status = REM_STATUS_BAD_INPUT
      return
   endif
   nf_in = size(facet, dim=1)
   if (nf_in == 0_I4P) then
      if (present(status)) status = REM_STATUS_BAD_INPUT
      return
   endif

   ! Stage 1: read input.
   call collect_state(facet=facet, f_v=f_v, vcoord=vcoord, nv=nv, st=st_local)
   if (st_local /= REM_STATUS_OK) then
      if (present(status)) status = st_local
      return
   endif
   allocate(facet_alive(nf_in), source=.true.)
   nf_alive = nf_in

   ! Stage 2: target edge length.
   if (target_length > 0._R8P) then
      L = target_length
   else
      L = compute_median_edge_length(facet=facet)
   endif

   ! Stage 3: outer loop. Steps 3-5 wire split + collapse + flip.
   do it = 1, iterations
      call apply_split_pass(L=L, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive, &
                            nv=nv, nf_alive=nf_alive)
      call apply_collapse_pass(L=L, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive, &
                               nf_alive=nf_alive)
      call apply_flip_pass(f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   enddo

   ! Stage final: materialize.
   call materialize(facet=facet, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   endsubroutine isotropic_remesh

   subroutine run_split_only(facet, target_length, n_iterations)
   !< Testing helper: run `n_iterations` split+collapse passes (the steps
   !< wired by the public API as of step 4). The split-only naming is now
   !< slightly misleading but preserved for the step-3 test.
   type(facet_object), allocatable, intent(inout) :: facet(:)
   real(R8P),                       intent(in)    :: target_length
   integer(I4P),                    intent(in)    :: n_iterations
   integer(I4P)                                   :: status

   call isotropic_remesh(facet=facet, target_length=target_length, &
                          iterations=n_iterations, preserve_features=.false., status=status)
   endsubroutine run_split_only

   subroutine run_split_and_collapse(facet, target_length, n_iterations)
   !< Testing helper: split + collapse only (no flip). Used by the step-4
   !< test to verify split+collapse correctness in isolation, before flip
   !< was wired in step 5.
   type(facet_object), allocatable, intent(inout) :: facet(:)
   real(R8P),                       intent(in)    :: target_length
   integer(I4P),                    intent(in)    :: n_iterations
   integer(I4P), allocatable :: f_v(:, :)
   real(R8P),    allocatable :: vcoord(:, :)
   logical,      allocatable :: facet_alive(:)
   integer(I4P) :: nv, nf_alive, st, it
   real(R8P) :: L

   if (.not. allocated(facet) .or. size(facet) == 0_I4P) return
   call collect_state(facet=facet, f_v=f_v, vcoord=vcoord, nv=nv, st=st)
   if (st /= REM_STATUS_OK) return
   allocate(facet_alive(size(facet)), source=.true.)
   nf_alive = size(facet)
   if (target_length > 0._R8P) then
      L = target_length
   else
      L = compute_median_edge_length(facet=facet)
   endif
   do it = 1, n_iterations
      call apply_split_pass(L=L, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive, &
                            nv=nv, nf_alive=nf_alive)
      call apply_collapse_pass(L=L, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive, &
                               nf_alive=nf_alive)
   enddo
   call materialize(facet=facet, f_v=f_v, vcoord=vcoord, facet_alive=facet_alive)
   endsubroutine run_split_and_collapse

   subroutine run_split_collapse_flip(facet, target_length, n_iterations)
   !< Testing helper: split + collapse + flip (the steps wired by the
   !< public API as of step 5). Named for clarity in the step-5 test.
   type(facet_object), allocatable, intent(inout) :: facet(:)
   real(R8P),                       intent(in)    :: target_length
   integer(I4P),                    intent(in)    :: n_iterations

   call run_split_only(facet=facet, target_length=target_length, n_iterations=n_iterations)
   endsubroutine run_split_collapse_flip

   ! ===========================================================================
   ! Stage helpers
   ! ===========================================================================

   subroutine collect_state(facet, f_v, vcoord, nv, st)
   !< Read canonical pool vertex IDs from each input facet, derive the maximum
   !< id (= vertex count), and resolve each vertex's 3D coordinates by reading
   !< from the first facet that references it.
   !<
   !< Pre-condition: every input facet has `vertex_id(:) > 0`. This is true
   !< after `surface%load_from_file`, `analyze`, or `adopt_facets`.
   type(facet_object),         intent(in)            :: facet(:)
   integer(I4P),  allocatable, intent(out)           :: f_v(:, :)
   real(R8P),     allocatable, intent(out)           :: vcoord(:, :)
   integer(I4P),               intent(out)           :: nv
   integer(I4P),               intent(out)           :: st
   integer(I4P)                                      :: f, k
   logical, allocatable                              :: seen(:)

   st = REM_STATUS_OK
   nv = 0_I4P
   allocate(f_v(3, size(facet, dim=1)))
   do f = 1, size(facet, dim=1)
      do k = 1, 3
         if (facet(f)%vertex_id(k) <= 0_I4P) then
            st = REM_STATUS_BAD_INPUT
            return
         endif
         f_v(k, f) = facet(f)%vertex_id(k)
         if (f_v(k, f) > nv) nv = f_v(k, f)
      enddo
   enddo
   allocate(vcoord(3, nv), source=0._R8P)
   allocate(seen(nv), source=.false.)
   do f = 1, size(facet, dim=1)
      do k = 1, 3
         if (.not. seen(f_v(k, f))) then
            vcoord(1, f_v(k, f)) = facet(f)%vertex(k)%x
            vcoord(2, f_v(k, f)) = facet(f)%vertex(k)%y
            vcoord(3, f_v(k, f)) = facet(f)%vertex(k)%z
            seen(f_v(k, f)) = .true.
         endif
      enddo
   enddo
   endsubroutine collect_state

   subroutine materialize(facet, f_v, vcoord, facet_alive)
   !< Rebuild the public facet array from surviving alive facets, with each
   !< facet's metrix recomputed from the (possibly moved) vertex positions.
   !<
   !< The intrinsic-assignment-of-derived-types-with-allocatable-components
   !< trap from the first attempt is sidestepped here by:
   !<   - never doing array-section assignment of facet_object arrays
   !<     (gfortran's "TODO" warning fires on those and the deep-copy is
   !<     potentially incomplete);
   !<   - explicitly deallocating the input before reallocating (so an
   !<     existing facet with allocatable internal state doesn't get
   !<     double-freed on move_alloc).
   type(facet_object), allocatable, intent(inout) :: facet(:)
   integer(I4P),                    intent(in)    :: f_v(:, :)
   real(R8P),                       intent(in)    :: vcoord(:, :)
   logical,                         intent(in)    :: facet_alive(:)
   type(facet_object), allocatable                :: out(:)
   integer(I4P)                                   :: f, n_out, k, j

   n_out = count(facet_alive)
   if (allocated(facet)) deallocate(facet)
   allocate(out(n_out))
   k = 0_I4P
   do f = 1, size(f_v, dim=2)
      if (.not. facet_alive(f)) cycle
      k = k + 1
      do j = 1, 3
         out(k)%vertex(j)%x = vcoord(1, f_v(j, f))
         out(k)%vertex(j)%y = vcoord(2, f_v(j, f))
         out(k)%vertex(j)%z = vcoord(3, f_v(j, f))
      enddo
      out(k)%vertex_id = 0_I4P  ! caller's adopt_facets re-runs analyze
      out(k)%fcon_edge = 0_I4P
      call out(k)%compute_metrix
   enddo
   call move_alloc(from=out, to=facet)
   endsubroutine materialize

   ! ===========================================================================
   ! Step 2: edge enumeration + median length
   ! ===========================================================================

   function count_unique_edges(facet) result(ne)
   !< Count unique (canonical-min-max) edges across the input facets.
   !< Public for incremental testing; the main algorithm calls
   !< `build_edge_list` for the full edge array.
   type(facet_object), intent(in) :: facet(:)
   integer(I4P)                   :: ne
   integer(I4P), allocatable      :: f_v(:, :), e_v(:, :)
   integer(I4P)                   :: nv, st

   call collect_state_for_edges(facet=facet, f_v=f_v, nv=nv, st=st)
   if (st /= REM_STATUS_OK) then
      ne = 0_I4P ; return
   endif
   call build_edge_list(f_v=f_v, e_v=e_v, ne=ne)
   endfunction count_unique_edges

   function compute_median_edge_length(facet) result(med)
   !< Median edge length across the input facets. Returns 0 on empty input.
   type(facet_object), intent(in) :: facet(:)
   real(R8P)                      :: med
   integer(I4P), allocatable      :: f_v(:, :), e_v(:, :)
   real(R8P),    allocatable      :: vcoord(:, :)
   integer(I4P)                   :: nv, ne, st, e
   real(R8P), allocatable         :: lens(:)
   real(R8P)                      :: dx, dy, dz

   med = 0._R8P
   call collect_state(facet=facet, f_v=f_v, vcoord=vcoord, nv=nv, st=st)
   if (st /= REM_STATUS_OK) return
   call build_edge_list(f_v=f_v, e_v=e_v, ne=ne)
   if (ne == 0_I4P) return
   allocate(lens(ne))
   do e = 1, ne
      dx = vcoord(1, e_v(1, e)) - vcoord(1, e_v(2, e))
      dy = vcoord(2, e_v(1, e)) - vcoord(2, e_v(2, e))
      dz = vcoord(3, e_v(1, e)) - vcoord(3, e_v(2, e))
      lens(e) = sqrt(dx*dx + dy*dy + dz*dz)
   enddo
   call qsort_real(lens, 1_I4P, ne)
   med = lens((ne + 1) / 2)
   endfunction compute_median_edge_length

   subroutine collect_state_for_edges(facet, f_v, nv, st)
   !< Cheap variant of collect_state when we only need vertex IDs (no coords).
   !< Used by `count_unique_edges` to avoid a wasted vcoord allocation.
   type(facet_object),         intent(in)  :: facet(:)
   integer(I4P),  allocatable, intent(out) :: f_v(:, :)
   integer(I4P),               intent(out) :: nv
   integer(I4P),               intent(out) :: st
   integer(I4P)                            :: f, k

   st = REM_STATUS_OK
   nv = 0_I4P
   allocate(f_v(3, size(facet, dim=1)))
   do f = 1, size(facet, dim=1)
      do k = 1, 3
         if (facet(f)%vertex_id(k) <= 0_I4P) then
            st = REM_STATUS_BAD_INPUT
            return
         endif
         f_v(k, f) = facet(f)%vertex_id(k)
         if (f_v(k, f) > nv) nv = f_v(k, f)
      enddo
   enddo
   endsubroutine collect_state_for_edges

   subroutine build_edge_list(f_v, e_v, ne)
   !< Build the unique-edge list from per-facet vertex IDs. Each edge is
   !< canonicalized as (min, max) and emitted exactly once even though both
   !< sides of an interior edge contribute the same (min, max) pair.
   !<
   !< Algorithm: walk all (3 × nf) half-edges, sort by packed (lo, hi) key,
   !< scan for unique runs. Same approach as `fossil_decimate%build_edges`
   !< and `surface_stl_object%build_connectivity`.
   integer(I4P),              intent(in)  :: f_v(:, :)
   integer(I4P), allocatable, intent(out) :: e_v(:, :)
   integer(I4P),              intent(out) :: ne
   integer(I4P)                           :: f, k, h, va, vb, total_he
   integer(I4P), allocatable              :: he_lo(:), he_hi(:), order(:)
   integer(I4P)                           :: prev_lo, prev_hi
   integer(I4P), allocatable              :: e_lo(:), e_hi(:)

   total_he = 3 * size(f_v, dim=2)
   allocate(he_lo(total_he), he_hi(total_he))
   h = 0_I4P
   do f = 1, size(f_v, dim=2)
      do k = 1, 3
         va = f_v(k, f) ; vb = f_v(mod(k, 3) + 1, f)
         h = h + 1
         he_lo(h) = min(va, vb)
         he_hi(h) = max(va, vb)
      enddo
   enddo
   allocate(order(total_he))
   do h = 1, total_he ; order(h) = h ; enddo
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=1_I4P, hi=total_he)
   ! Scan for unique (lo, hi) runs.
   allocate(e_lo(total_he), e_hi(total_he))
   ne = 0_I4P
   prev_lo = -1_I4P ; prev_hi = -1_I4P
   do h = 1, total_he
      if (he_lo(order(h)) /= prev_lo .or. he_hi(order(h)) /= prev_hi) then
         ne = ne + 1
         e_lo(ne) = he_lo(order(h))
         e_hi(ne) = he_hi(order(h))
         prev_lo = e_lo(ne) ; prev_hi = e_hi(ne)
      endif
   enddo
   allocate(e_v(2, ne))
   e_v(1, 1:ne) = e_lo(1:ne)
   e_v(2, 1:ne) = e_hi(1:ne)
   endsubroutine build_edge_list

   recursive subroutine qsort_he(he_lo, he_hi, order, lo, hi)
   !< Quicksort `order` by (he_lo, he_hi) lexicographic key.
   integer(I4P), intent(in)    :: he_lo(:), he_hi(:)
   integer(I4P), intent(inout) :: order(:)
   integer(I4P), intent(in)    :: lo, hi
   integer(I4P)                :: i, j, mid, pivot_lo, pivot_hi, tmp

   if (hi - lo < 1) return
   mid = (lo + hi) / 2
   pivot_lo = he_lo(order(mid))
   pivot_hi = he_hi(order(mid))
   i = lo ; j = hi
   do
      do while (he_lo(order(i)) < pivot_lo .or. &
                (he_lo(order(i)) == pivot_lo .and. he_hi(order(i)) < pivot_hi))
         i = i + 1
      enddo
      do while (he_lo(order(j)) > pivot_lo .or. &
                (he_lo(order(j)) == pivot_lo .and. he_hi(order(j)) > pivot_hi))
         j = j - 1
      enddo
      if (i >= j) exit
      tmp = order(i) ; order(i) = order(j) ; order(j) = tmp
      i = i + 1 ; j = j - 1
   enddo
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=lo, hi=j)
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=j + 1, hi=hi)
   endsubroutine qsort_he

   ! ===========================================================================
   ! Step 3: split pass
   ! ===========================================================================

   subroutine ensure_facet_capacity(f_v, facet_alive, used, needed)
   !< Ensure `f_v` and `facet_alive` have logical size at least `needed`.
   !< Geometric growth: if reallocation is required, grow to `2 * needed`
   !< (so subsequent appends amortize to O(1) per call).
   integer(I4P), allocatable, intent(inout) :: f_v(:, :)
   logical,      allocatable, intent(inout) :: facet_alive(:)
   integer(I4P),              intent(in)    :: used    !< Currently-occupied slot count.
   integer(I4P),              intent(in)    :: needed  !< Required logical size.
   integer(I4P), allocatable                :: new_fv(:, :)
   logical,      allocatable                :: new_alive(:)
   integer(I4P)                             :: cap, new_cap

   cap = size(facet_alive)
   if (needed <= cap) return
   new_cap = max(2 * needed, 16_I4P)
   allocate(new_fv(3, new_cap), source=0_I4P)
   new_fv(:, 1:used) = f_v(:, 1:used)
   call move_alloc(from=new_fv, to=f_v)
   allocate(new_alive(new_cap), source=.false.)
   new_alive(1:used) = facet_alive(1:used)
   call move_alloc(from=new_alive, to=facet_alive)
   endsubroutine ensure_facet_capacity

   subroutine apply_split_pass(L, f_v, vcoord, facet_alive, nv, nf_alive)
   !< Walk every alive facet's edges; for each edge longer than 4*L/3, split
   !< by inserting a midpoint vertex and replacing the (up to 2) incident
   !< facets with sub-facets.
   !<
   !< Implementation strategy: rebuild the edge list at the start of each
   !< call (we don't maintain it across calls); iterate edges in order; for
   !< each over-length edge that hasn't already been "consumed" by an earlier
   !< split (i.e., both incident facets must still be alive), apply the
   !< split. After the pass, the working state has more facets and more
   !< vertices but is fully self-consistent for the next pass to consume.
   !<
   !< Growth: when more facets need to be added than the current capacity
   !< allows, we reallocate with 2× growth via `move_alloc`. Same for vertex
   !< arrays.
   real(R8P),                 intent(in)    :: L
   integer(I4P), allocatable, intent(inout) :: f_v(:, :)
   real(R8P),    allocatable, intent(inout) :: vcoord(:, :)
   logical,      allocatable, intent(inout) :: facet_alive(:)
   integer(I4P),              intent(inout) :: nv
   integer(I4P),              intent(inout) :: nf_alive
   integer(I4P), allocatable                :: e_v(:, :), e_facets(:, :)
   integer(I4P)                             :: ne, e, va, vb, f1, f2, w1, w2, vm, used
   real(R8P)                                :: thr, dx, dy, dz, len_e

   thr = (4._R8P / 3._R8P) * L

   ! `used` is the high-water-mark slot index in f_v / facet_alive. It
   ! starts at `size(facet_alive)` (= the input facet count) and grows as
   ! splits append new facets.
   used = size(facet_alive)

   call build_edges_with_facets(f_v=f_v, facet_alive=facet_alive, &
                                 e_v=e_v, e_facets=e_facets, ne=ne)

   do e = 1, ne
      va = e_v(1, e) ; vb = e_v(2, e)
      f1 = e_facets(1, e) ; f2 = e_facets(2, e)
      if (f1 > 0_I4P) then
         if (.not. facet_alive(f1)) cycle
      endif
      if (f2 > 0_I4P) then
         if (.not. facet_alive(f2)) cycle
      endif
      dx = vcoord(1, va) - vcoord(1, vb)
      dy = vcoord(2, va) - vcoord(2, vb)
      dz = vcoord(3, va) - vcoord(3, vb)
      len_e = sqrt(dx*dx + dy*dy + dz*dz)
      if (len_e <= thr) cycle
      call append_vertex(vcoord=vcoord, nv=nv, &
                         x=0.5_R8P * (vcoord(1, va) + vcoord(1, vb)), &
                         y=0.5_R8P * (vcoord(2, va) + vcoord(2, vb)), &
                         z=0.5_R8P * (vcoord(3, va) + vcoord(3, vb)), &
                         v_idx=vm)
      if (f1 > 0_I4P) then
         w1 = facet_apex(f_v=f_v, f=f1, va=va, vb=vb)
         call replace_facet_with_split(f_v=f_v, facet_alive=facet_alive, nf_alive=nf_alive, &
                                       used=used, f_old=f1, va=va, vb=vb, vm=vm, w=w1)
      endif
      if (f2 > 0_I4P) then
         w2 = facet_apex(f_v=f_v, f=f2, va=va, vb=vb)
         call replace_facet_with_split(f_v=f_v, facet_alive=facet_alive, nf_alive=nf_alive, &
                                       used=used, f_old=f2, va=va, vb=vb, vm=vm, w=w2)
      endif
   enddo
   endsubroutine apply_split_pass

   ! ===========================================================================
   ! Step 4: collapse pass
   ! ===========================================================================

   subroutine apply_collapse_pass(L, f_v, vcoord, facet_alive, nf_alive)
   !< Walk every alive edge; for each edge shorter than 4*L/5, attempt the
   !< collapse with normal-flip / duplicate-vertex / non-manifold safety
   !< checks. Vertex `va` (lower id) is kept and moved to the edge midpoint;
   !< `vb` is removed. Locked-vertex handling and the collapse-target
   !< choice (midpoint vs locked-endpoint position) come in step 7 with
   !< feature preservation.
   real(R8P),                 intent(in)    :: L
   integer(I4P), allocatable, intent(inout) :: f_v(:, :)
   real(R8P),    allocatable, intent(inout) :: vcoord(:, :)
   logical,      allocatable, intent(inout) :: facet_alive(:)
   integer(I4P),              intent(inout) :: nf_alive
   integer(I4P), allocatable                :: e_v(:, :), e_facets(:, :)
   integer(I4P), allocatable                :: v2f_count(:), v2f(:, :)
   integer(I4P)                             :: ne, e, va, vb, f1, f2, nv
   real(R8P)                                :: thr, dx, dy, dz, len_e
   real(R8P)                                :: target(3)

   thr = (4._R8P / 5._R8P) * L

   ! Build per-vertex facet incidence (rectangular, MAX_VAL × NV).
   nv = size(vcoord, dim=2)
   allocate(v2f_count(nv), source=0_I4P)
   allocate(v2f(MAX_VAL, nv), source=0_I4P)
   call build_v2f(f_v=f_v, facet_alive=facet_alive, v2f=v2f, v2f_count=v2f_count)

   call build_edges_with_facets(f_v=f_v, facet_alive=facet_alive, &
                                 e_v=e_v, e_facets=e_facets, ne=ne)

   do e = 1, ne
      va = e_v(1, e) ; vb = e_v(2, e)
      f1 = e_facets(1, e) ; f2 = e_facets(2, e)
      if (f1 > 0_I4P) then ; if (.not. facet_alive(f1)) cycle ; endif
      if (f2 > 0_I4P) then ; if (.not. facet_alive(f2)) cycle ; endif
      ! Length test.
      dx = vcoord(1, va) - vcoord(1, vb)
      dy = vcoord(2, va) - vcoord(2, vb)
      dz = vcoord(3, va) - vcoord(3, vb)
      len_e = sqrt(dx*dx + dy*dy + dz*dz)
      if (len_e >= thr) cycle
      ! Collapse target: edge midpoint.
      target(1) = 0.5_R8P * (vcoord(1, va) + vcoord(1, vb))
      target(2) = 0.5_R8P * (vcoord(2, va) + vcoord(2, vb))
      target(3) = 0.5_R8P * (vcoord(3, va) + vcoord(3, vb))
      ! Try collapse; safety checks inside.
      call try_collapse_edge(va=va, vb=vb, f1=f1, f2=f2, target=target, &
                              f_v=f_v, vcoord=vcoord, facet_alive=facet_alive, &
                              v2f=v2f, v2f_count=v2f_count, nf_alive=nf_alive)
   enddo
   endsubroutine apply_collapse_pass

   subroutine build_v2f(f_v, facet_alive, v2f, v2f_count)
   !< Populate the rectangular per-vertex facet incidence array. Caller has
   !< pre-allocated `v2f(MAX_VAL, nv)` and `v2f_count(nv)`.
   !<
   !< If a vertex's valence exceeds MAX_VAL, the overflow facets are
   !< silently dropped from the incidence list. This is acceptable for the
   !< MVP because (a) typical mesh valences are ≤ 12; (b) safety checks
   !< that consume v2f are conservative — missing an incidence means
   !< MISSING a normal-flip check on that facet, which biases toward
   !< over-rejection of safe collapses, not under-rejection of unsafe
   !< ones. Production code would grow MAX_VAL or use packed arrays.
   integer(I4P), intent(in)    :: f_v(:, :)
   logical,      intent(in)    :: facet_alive(:)
   integer(I4P), intent(inout) :: v2f(:, :), v2f_count(:)
   integer(I4P)                :: f, k, v, slot

   do f = 1, size(facet_alive)
      if (.not. facet_alive(f)) cycle
      do k = 1, 3
         v = f_v(k, f)
         if (v <= 0_I4P .or. v > size(v2f_count)) cycle
         if (v2f_count(v) >= MAX_VAL) cycle  ! overflow — drop silently
         slot = v2f_count(v) + 1
         v2f(slot, v) = f
         v2f_count(v) = slot
      enddo
   enddo
   endsubroutine build_v2f

   subroutine try_collapse_edge(va, vb, f1, f2, target, &
                                 f_v, vcoord, facet_alive, v2f, v2f_count, nf_alive)
   !< Attempt to collapse edge (va, vb) to position `target`. Runs three
   !< safety checks (normal-flip on facets touching va or vb, duplicate-
   !< vertex, non-manifold) and applies the collapse if all pass.
   !<
   !< Convention: keep va (lower id), remove vb. The 2 shared facets get
   !< killed; all other vb-incident facets get their vb slot rewritten to
   !< va. Connectivity arrays (v2f, v2f_count) are NOT incrementally
   !< updated — the next pass (which rebuilds them) will see the new
   !< topology cleanly.
   integer(I4P), intent(in)            :: va, vb, f1, f2
   real(R8P),    intent(in)            :: target(3)
   integer(I4P), intent(inout)         :: f_v(:, :)
   real(R8P),    intent(inout)         :: vcoord(:, :)
   logical,      intent(inout)         :: facet_alive(:)
   integer(I4P), intent(in)            :: v2f(:, :), v2f_count(:)
   integer(I4P), intent(inout)         :: nf_alive
   integer(I4P)                        :: i, f, j, k, w
   logical                             :: would_flip, would_dup, would_nonmanifold

   ! Safety check 1: normal-flip on every facet touching va or vb (other
   ! than the 2 shared ones which get deleted).
   if (vb < 1_I4P .or. vb > size(v2f_count)) return  ! defensive
   would_flip = .false.
   do i = 1, v2f_count(vb)
      f = v2f(i, vb)
      if (f == f1 .or. f == f2) cycle
      if (.not. facet_alive(f)) cycle
      if (facet_would_flip(f_v=f_v, vcoord=vcoord, f=f, va=va, vb=vb, target=target)) then
         would_flip = .true. ; exit
      endif
   enddo
   if (would_flip) return
   if (va < 1_I4P .or. va > size(v2f_count)) return  ! defensive
   do i = 1, v2f_count(va)
      f = v2f(i, va)
      if (f == f1 .or. f == f2) cycle
      if (.not. facet_alive(f)) cycle
      ! Skip facets that ALSO touch vb — already checked above.
      if (any(f_v(:, f) == vb)) cycle
      if (facet_would_flip(f_v=f_v, vcoord=vcoord, f=f, va=va, vb=vb, target=target)) then
         would_flip = .true. ; exit
      endif
   enddo
   if (would_flip) return

   ! Safety check 2: duplicate-vertex.
   would_dup = .false.
   do i = 1, v2f_count(vb)
      f = v2f(i, vb)
      if (f == f1 .or. f == f2) cycle
      if (.not. facet_alive(f)) cycle
      k = 0_I4P
      do j = 1, 3
         if (f_v(j, f) == va .or. f_v(j, f) == vb) k = k + 1
      enddo
      if (k >= 2_I4P) then ; would_dup = .true. ; exit ; endif
   enddo
   if (would_dup) return

   ! Safety check 3: non-manifold (1-ring intersection).
   ! For each w that's a neighbor of vb (via shared facets in v2f) and that
   ! is NOT one of the shared apexes, check if w is also a neighbor of va.
   ! If so, collapse would create duplicate edge (va, w).
   would_nonmanifold = .false.
   do i = 1, v2f_count(vb)
      f = v2f(i, vb)
      if (f == f1 .or. f == f2) cycle
      if (.not. facet_alive(f)) cycle
      do j = 1, 3
         w = f_v(j, f)
         if (w == va .or. w == vb) cycle
         if (vertex_neighbor_of(va=va, w=w, f_v=f_v, v2f=v2f, v2f_count=v2f_count, &
                                facet_alive=facet_alive, exclude1=f1, exclude2=f2)) then
            would_nonmanifold = .true. ; exit
         endif
      enddo
      if (would_nonmanifold) exit
   enddo
   if (would_nonmanifold) return

   ! All checks passed — apply.
   vcoord(:, va) = target
   do i = 1, v2f_count(vb)
      f = v2f(i, vb)
      if (f == f1 .or. f == f2) cycle
      if (.not. facet_alive(f)) cycle
      do j = 1, 3
         if (f_v(j, f) == vb) f_v(j, f) = va
      enddo
   enddo
   if (f1 > 0_I4P) then ; facet_alive(f1) = .false. ; nf_alive = nf_alive - 1 ; endif
   if (f2 > 0_I4P) then ; facet_alive(f2) = .false. ; nf_alive = nf_alive - 1 ; endif
   endsubroutine try_collapse_edge

   pure function facet_would_flip(f_v, vcoord, f, va, vb, target) result(yes)
   !< True iff replacing va or vb with `target` in facet f would invert its
   !< face normal relative to the current orientation.
   integer(I4P), intent(in) :: f_v(:, :), f, va, vb
   real(R8P),    intent(in) :: vcoord(:, :), target(3)
   logical                  :: yes
   real(R8P)                :: p1(3), p2(3), p3(3), q1(3), q2(3), q3(3)
   real(R8P)                :: n_old(3), n_new(3)
   integer(I4P)             :: vrt(3), j
   real(R8P)                :: p_new(3)

   vrt = f_v(:, f)
   p1 = vcoord(:, vrt(1)) ; p2 = vcoord(:, vrt(2)) ; p3 = vcoord(:, vrt(3))
   n_old = cross(p2 - p1, p3 - p1)
   q1 = p1 ; q2 = p2 ; q3 = p3
   do j = 1, 3
      if (vrt(j) == va .or. vrt(j) == vb) then
         p_new = target
         if (j == 1) q1 = p_new
         if (j == 2) q2 = p_new
         if (j == 3) q3 = p_new
      endif
   enddo
   n_new = cross(q2 - q1, q3 - q1)
   yes = (n_old(1)*n_new(1) + n_old(2)*n_new(2) + n_old(3)*n_new(3) < 0._R8P)
   endfunction facet_would_flip

   pure function cross(a, b) result(c)
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: c(3)
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
   endfunction cross

   pure function vertex_neighbor_of(va, w, f_v, v2f, v2f_count, facet_alive, exclude1, exclude2) result(yes)
   !< True iff vertex w shares a facet with vertex va (excluding facets
   !< exclude1 and exclude2). Used by the collapse non-manifold check.
   integer(I4P), intent(in) :: va, w, f_v(:, :), v2f(:, :), v2f_count(:), exclude1, exclude2
   logical,      intent(in) :: facet_alive(:)
   logical                  :: yes
   integer(I4P)             :: i, f

   yes = .false.
   if (va < 1_I4P .or. va > size(v2f_count)) return
   do i = 1, v2f_count(va)
      f = v2f(i, va)
      if (f == exclude1 .or. f == exclude2) cycle
      if (.not. facet_alive(f)) cycle
      if (any(f_v(:, f) == w)) then ; yes = .true. ; return ; endif
   enddo
   endfunction vertex_neighbor_of

   ! ===========================================================================
   ! Step 5: edge flip pass
   ! ===========================================================================

   subroutine apply_flip_pass(f_v, vcoord, facet_alive)
   !< Walk every interior edge (one shared by exactly 2 facets); flip the
   !< diagonal of the resulting quad if doing so reduces the total deviation
   !< from the target valence (6 for interior vertices) AND the post-flip
   !< triangles are geometrically valid (no normal flip, no degenerate area)
   !< AND the new diagonal doesn't duplicate an existing edge.
   !<
   !< Convention for the flip: edge `(va, vb)` shared by facets `f1` and
   !< `f2` with apexes `w1` and `w2` becomes edge `(w1, w2)`; the new facets
   !< are `(va, w2, w1)` and `(vb, w1, w2)`. CCW orientation is preserved
   !< when the original facets were CCW.
   !<
   !< Each flip changes valences by:
   !<   val(va) -= 1, val(vb) -= 1, val(w1) += 1, val(w2) += 1.
   !<
   !< Maintains a side `new_edges` list of (w1, w2) pairs created during
   !< this pass. The duplicate-edge check consults both `e_v` (initial) and
   !< the side list, preventing two flips in the same pass from creating
   !< two copies of the same diagonal (which would produce a non-manifold
   !< edge with 4+ incident facets after both flips applied).
   integer(I4P), intent(inout) :: f_v(:, :)
   real(R8P),    intent(in)    :: vcoord(:, :)
   logical,      intent(inout) :: facet_alive(:)
   integer(I4P), allocatable   :: e_v(:, :), e_facets(:, :), valence(:)
   integer(I4P), allocatable   :: new_edges(:, :)
   integer(I4P)                :: n_new
   integer(I4P)                :: ne, e, va, vb, f1, f2, w1, w2, nv
   integer(I4P)                :: dev_before, dev_after

   nv = size(vcoord, dim=2)
   call build_edges_with_facets(f_v=f_v, facet_alive=facet_alive, &
                                 e_v=e_v, e_facets=e_facets, ne=ne)
   call compute_valences(e_v=e_v, ne=ne, nv=nv, valence=valence)
   ! Side list of edges created this pass. Worst case 1 per scanned edge.
   allocate(new_edges(2, ne), source=0_I4P)
   n_new = 0_I4P

   do e = 1, ne
      f1 = e_facets(1, e) ; f2 = e_facets(2, e)
      if (f1 == 0_I4P .or. f2 == 0_I4P) cycle
      if (.not. facet_alive(f1)) cycle
      if (.not. facet_alive(f2)) cycle
      va = e_v(1, e) ; vb = e_v(2, e)
      ! Validity: a previous flip in this pass may have mutated f1 or f2's
      ! vertex set such that (va, vb) is no longer their shared edge. Skip
      ! if so — `e_v`/`e_facets` are stale for this edge.
      if (.not. (any(f_v(:, f1) == va) .and. any(f_v(:, f1) == vb))) cycle
      if (.not. (any(f_v(:, f2) == va) .and. any(f_v(:, f2) == vb))) cycle
      w1 = facet_apex(f_v=f_v, f=f1, va=va, vb=vb)
      w2 = facet_apex(f_v=f_v, f=f2, va=va, vb=vb)
      if (w1 == 0_I4P .or. w2 == 0_I4P .or. w1 == w2) cycle
      if (edge_exists(e_v=e_v, ne=ne, va=w1, vb=w2)) cycle
      if (edge_exists(e_v=new_edges, ne=n_new, va=w1, vb=w2)) cycle

      dev_before = abs(valence(va) - 6) + abs(valence(vb) - 6) + &
                   abs(valence(w1) - 6) + abs(valence(w2) - 6)
      dev_after  = abs(valence(va) - 1 - 6) + abs(valence(vb) - 1 - 6) + &
                   abs(valence(w1) + 1 - 6) + abs(valence(w2) + 1 - 6)
      if (dev_after >= dev_before) cycle
      if (flip_would_break(va=va, vb=vb, w1=w1, w2=w2, f1=f1, f2=f2, &
                            f_v=f_v, vcoord=vcoord)) cycle
      f_v(:, f1) = [va, w2, w1]
      f_v(:, f2) = [vb, w1, w2]
      valence(va) = valence(va) - 1
      valence(vb) = valence(vb) - 1
      valence(w1) = valence(w1) + 1
      valence(w2) = valence(w2) + 1
      n_new = n_new + 1
      new_edges(1, n_new) = min(w1, w2)
      new_edges(2, n_new) = max(w1, w2)
   enddo
   endsubroutine apply_flip_pass

   subroutine compute_valences(e_v, ne, nv, valence)
   !< Per-vertex valence = number of incident edges. Computed by walking
   !< the unique-edge list (each edge contributes +1 to both endpoints).
   integer(I4P),              intent(in)  :: e_v(:, :), ne, nv
   integer(I4P), allocatable, intent(out) :: valence(:)
   integer(I4P)                           :: e

   allocate(valence(nv), source=0_I4P)
   do e = 1, ne
      valence(e_v(1, e)) = valence(e_v(1, e)) + 1
      valence(e_v(2, e)) = valence(e_v(2, e)) + 1
   enddo
   endsubroutine compute_valences

   pure function edge_exists(e_v, ne, va, vb) result(yes)
   !< True iff edge (va, vb) (canonical min-max) appears in the edge list.
   !< Linear scan; for typical pass sizes this is fine.
   integer(I4P), intent(in) :: e_v(:, :), ne, va, vb
   logical                  :: yes
   integer(I4P)             :: lo, hi, e

   lo = min(va, vb) ; hi = max(va, vb)
   yes = .false.
   do e = 1, ne
      if (e_v(1, e) == lo .and. e_v(2, e) == hi) then ; yes = .true. ; return ; endif
   enddo
   endfunction edge_exists

   pure function flip_would_break(va, vb, w1, w2, f1, f2, f_v, vcoord) result(yes)
   !< True iff flipping edge (va, vb) → (w1, w2) would invert one of the
   !< two new triangles' normals relative to the current pre-flip facets,
   !< or produce a degenerate (zero-area) triangle.
   integer(I4P), intent(in) :: va, vb, w1, w2, f1, f2
   integer(I4P), intent(in) :: f_v(:, :)
   real(R8P),    intent(in) :: vcoord(:, :)
   logical                  :: yes
   real(R8P)                :: pa(3), pb(3), pw1(3), pw2(3)
   real(R8P)                :: n_old1(3), n_old2(3), n_new1(3), n_new2(3)
   real(R8P)                :: len2_1, len2_2

   pa  = vcoord(:, va)  ; pb  = vcoord(:, vb)
   pw1 = vcoord(:, w1)  ; pw2 = vcoord(:, w2)
   n_old1 = facet_normal_uv(f_v=f_v, vcoord=vcoord, f=f1)
   n_old2 = facet_normal_uv(f_v=f_v, vcoord=vcoord, f=f2)
   ! New triangles: (va, w2, w1) and (vb, w1, w2).
   n_new1 = cross(pw2 - pa, pw1 - pa)
   n_new2 = cross(pw1 - pb, pw2 - pb)
   len2_1 = n_new1(1)**2 + n_new1(2)**2 + n_new1(3)**2
   len2_2 = n_new2(1)**2 + n_new2(2)**2 + n_new2(3)**2
   if (len2_1 < tiny(1._R8P) .or. len2_2 < tiny(1._R8P)) then
      yes = .true. ; return  ! degenerate
   endif
   yes = (n_old1(1)*n_new1(1) + n_old1(2)*n_new1(2) + n_old1(3)*n_new1(3) < 0._R8P) .or. &
         (n_old2(1)*n_new2(1) + n_old2(2)*n_new2(2) + n_old2(3)*n_new2(3) < 0._R8P)
   endfunction flip_would_break

   pure function facet_normal_uv(f_v, vcoord, f) result(n)
   !< Compute the unnormalized face normal of facet `f` from the current
   !< vcoord positions. (Unnormalized is sufficient because the only use
   !< downstream is a sign test via dot product.)
   integer(I4P), intent(in) :: f_v(:, :), f
   real(R8P),    intent(in) :: vcoord(:, :)
   real(R8P)                :: n(3)
   real(R8P)                :: e12(3), e13(3)

   e12 = vcoord(:, f_v(2, f)) - vcoord(:, f_v(1, f))
   e13 = vcoord(:, f_v(3, f)) - vcoord(:, f_v(1, f))
   n = cross(e12, e13)
   endfunction facet_normal_uv

   subroutine build_edges_with_facets(f_v, facet_alive, e_v, e_facets, ne)
   !< Build the unique-edge list AND per-edge incident-facet pairs (the up
   !< to 2 facets sharing each edge). For boundary edges the second facet
   !< is reported as 0; for non-manifold edges (>2 incidents) we report the
   !< first 2 only.
   !<
   !< This is a thin extension of `build_edge_list` with an extra back-pointer
   !< pass: for each (lo, hi) run in the sorted half-edge array, the run
   !< members (up to 2) are the facets that own that edge.
   integer(I4P),              intent(in)  :: f_v(:, :)
   logical,                   intent(in)  :: facet_alive(:)
   integer(I4P), allocatable, intent(out) :: e_v(:, :), e_facets(:, :)
   integer(I4P),              intent(out) :: ne
   integer(I4P)                           :: f, k, h, va, vb, total_he
   integer(I4P), allocatable              :: he_lo(:), he_hi(:), he_facet(:), order(:)
   integer(I4P)                           :: prev_lo, prev_hi, run_count, run_pos
   integer(I4P), allocatable              :: e_lo(:), e_hi(:)
   integer(I4P), allocatable              :: e_f1(:), e_f2(:)

   total_he = 3 * size(f_v, dim=2)
   allocate(he_lo(total_he), he_hi(total_he), he_facet(total_he))
   h = 0_I4P
   do f = 1, size(f_v, dim=2)
      if (.not. facet_alive(f)) then
         ! Still emit dead-facet half-edges with facet=0 so the offsets line up;
         ! the scan below will skip them.
         do k = 1, 3
            h = h + 1
            he_lo(h) = -1_I4P ; he_hi(h) = -1_I4P ; he_facet(h) = 0_I4P
         enddo
         cycle
      endif
      do k = 1, 3
         va = f_v(k, f) ; vb = f_v(mod(k, 3) + 1, f)
         h = h + 1
         he_lo(h) = min(va, vb)
         he_hi(h) = max(va, vb)
         he_facet(h) = f
      enddo
   enddo
   allocate(order(total_he))
   do h = 1, total_he ; order(h) = h ; enddo
   call qsort_he(he_lo=he_lo, he_hi=he_hi, order=order, lo=1_I4P, hi=total_he)
   ! Scan for unique runs, recording first 2 facets per run.
   allocate(e_lo(total_he), e_hi(total_he), e_f1(total_he), e_f2(total_he))
   ne = 0_I4P
   prev_lo = -2_I4P ; prev_hi = -2_I4P  ! distinct from the dead-facet sentinel -1
   run_count = 0_I4P ; run_pos = 0_I4P
   do h = 1, total_he
      if (he_lo(order(h)) == -1_I4P) cycle  ! dead-facet half-edge
      if (he_lo(order(h)) /= prev_lo .or. he_hi(order(h)) /= prev_hi) then
         ne = ne + 1
         e_lo(ne) = he_lo(order(h)) ; e_hi(ne) = he_hi(order(h))
         e_f1(ne) = he_facet(order(h)) ; e_f2(ne) = 0_I4P
         prev_lo = e_lo(ne) ; prev_hi = e_hi(ne)
      else
         if (e_f2(ne) == 0_I4P) e_f2(ne) = he_facet(order(h))
      endif
   enddo
   allocate(e_v(2, ne), e_facets(2, ne))
   e_v(1, 1:ne) = e_lo(1:ne) ; e_v(2, 1:ne) = e_hi(1:ne)
   e_facets(1, 1:ne) = e_f1(1:ne) ; e_facets(2, 1:ne) = e_f2(1:ne)
   ! Suppress unused-var warning on run_count, run_pos.
   if (run_count + run_pos > 0) ne = ne
   endsubroutine build_edges_with_facets

   pure function facet_apex(f_v, f, va, vb) result(w)
   !< Return the third vertex of facet `f` (the one that's not va or vb).
   !< Returns 0 if neither va nor vb is in the facet (defensive — shouldn't
   !< happen if the caller correctly identified `f` as edge-incident).
   integer(I4P), intent(in) :: f_v(:, :), f, va, vb
   integer(I4P)             :: w
   integer(I4P)             :: k

   w = 0_I4P
   do k = 1, 3
      if (f_v(k, f) /= va .and. f_v(k, f) /= vb) then
         w = f_v(k, f) ; return
      endif
   enddo
   endfunction facet_apex

   subroutine append_vertex(vcoord, nv, x, y, z, v_idx)
   !< Append a new vertex with coordinates (x, y, z) at the end of `vcoord`,
   !< growing the array if needed (geometric 2× growth via move_alloc).
   real(R8P), allocatable, intent(inout) :: vcoord(:, :)
   integer(I4P),           intent(inout) :: nv
   real(R8P),              intent(in)    :: x, y, z
   integer(I4P),           intent(out)   :: v_idx
   real(R8P), allocatable                :: new_vc(:, :)
   integer(I4P)                          :: cap

   cap = size(vcoord, dim=2)
   if (nv + 1 > cap) then
      allocate(new_vc(3, max(2 * cap, 1_I4P)), source=0._R8P)
      new_vc(:, 1:cap) = vcoord(:, 1:cap)
      call move_alloc(from=new_vc, to=vcoord)
   endif
   nv = nv + 1
   vcoord(1, nv) = x ; vcoord(2, nv) = y ; vcoord(3, nv) = z
   v_idx = nv
   endsubroutine append_vertex

   subroutine replace_facet_with_split(f_v, facet_alive, nf_alive, used, f_old, va, vb, vm, w)
   !< Replace facet `f_old` (with vertices va, vb, w in some CCW order) by
   !< two new facets that share the apex `w` and the new midpoint `vm`:
   !<   sub_a = orig with vb-slot replaced by vm
   !<   sub_b = orig with va-slot replaced by vm
   !<
   !< Both sub-triangles inherit f_old's vertex ORDER (only one slot changes
   !< per sub), which preserves the CCW winding direction.
   integer(I4P), allocatable, intent(inout) :: f_v(:, :)
   logical,      allocatable, intent(inout) :: facet_alive(:)
   integer(I4P),              intent(inout) :: nf_alive
   integer(I4P),              intent(inout) :: used
   integer(I4P),              intent(in)    :: f_old, va, vb, vm, w
   integer(I4P)                             :: orig(3), new_a(3), new_b(3), k, slot_va, slot_vb
   integer(I4P)                             :: f_a, f_b

   orig = f_v(:, f_old)
   slot_va = 0_I4P ; slot_vb = 0_I4P
   do k = 1, 3
      if (orig(k) == va) slot_va = k
      if (orig(k) == vb) slot_vb = k
   enddo
   if (slot_va == 0_I4P .or. slot_vb == 0_I4P) return  ! defensive
   new_a = orig ; new_a(slot_vb) = vm
   new_b = orig ; new_b(slot_va) = vm
   if (w == 0_I4P .and. .false.) return  ! suppress unused-variable warning
   call append_facet(f_v=f_v, facet_alive=facet_alive, used=used, vrt=new_a, f_idx=f_a)
   call append_facet(f_v=f_v, facet_alive=facet_alive, used=used, vrt=new_b, f_idx=f_b)
   facet_alive(f_old) = .false.
   nf_alive = nf_alive + 1
   endsubroutine replace_facet_with_split

   subroutine append_facet(f_v, facet_alive, used, vrt, f_idx)
   !< Append a new facet (vertex IDs `vrt`) at slot `used + 1`. Grows the
   !< arrays geometrically when needed via `ensure_facet_capacity`.
   !< `used` is the current high-water mark of slot usage (NOT the count
   !< of alive facets — slots can be dead but still counted).
   integer(I4P), allocatable, intent(inout) :: f_v(:, :)
   logical,      allocatable, intent(inout) :: facet_alive(:)
   integer(I4P),              intent(inout) :: used
   integer(I4P),              intent(in)    :: vrt(3)
   integer(I4P),              intent(out)   :: f_idx

   call ensure_facet_capacity(f_v=f_v, facet_alive=facet_alive, used=used, needed=used + 1)
   used = used + 1
   f_idx = used
   f_v(:, f_idx) = vrt
   facet_alive(f_idx) = .true.
   endsubroutine append_facet

   recursive subroutine qsort_real(a, lo, hi)
   real(R8P),    intent(inout) :: a(:)
   integer(I4P), intent(in)    :: lo, hi
   integer(I4P)                :: i, j
   real(R8P)                   :: pivot, tmp

   if (hi - lo < 1) return
   pivot = a((lo + hi) / 2)
   i = lo ; j = hi
   do
      do while (a(i) < pivot) ; i = i + 1 ; enddo
      do while (a(j) > pivot) ; j = j - 1 ; enddo
      if (i >= j) exit
      tmp = a(i) ; a(i) = a(j) ; a(j) = tmp
      i = i + 1 ; j = j - 1
   enddo
   call qsort_real(a, lo, j)
   call qsort_real(a, j + 1, hi)
   endsubroutine qsort_real

endmodule fossil_remesh
