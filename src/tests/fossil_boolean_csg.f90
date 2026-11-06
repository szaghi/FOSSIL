!< FOSSIL test-side module: CSG chains for the boolean validation (issue #63, roadmap #58 WP5, step 4).
module fossil_boolean_csg
!< FOSSIL test-side module: CSG chains for the boolean validation (issue #63, roadmap #58 WP5, step 4).
!<
!< Chains of booleans whose exact volume and second moments follow from independent booleans, so that every step is
!< checked and drift along the chain shows:
!< - `drill_chain(k, n, ...)`: the unit cube minus k^2 disjoint vertical capped cylinders (marching cubes, n^3) on a
!<   k x k grid, each through the block. Step j: V(R_j) = V(R_j-1) - V(cube ∩ C_j), the intersection computed apart
!<   against the original cube (the cylinders are disjoint, so C_j ∩ R_j-1 = C_j ∩ cube);
!< - `lattice_chain(k, n, ...)`: the union of k^2 spheres (n^3) on a k x k grid in the plane z = 0, radius 0.6 of the
!<   spacing: neighbours along x and y overlap, diagonal ones do not, so no three overlap. Step j:
!<   V(R_j) = V(R_j-1) + V(S_j) - sum over the earlier neighbours i of V(S_i ∩ S_j).
!< The second moments (raw, about the origin) follow the same sums. After each step: the residuals of both, relative
!< to the first operand, watertight and manifold, the facet count and the time. The status of every boolean must be
!< OK (a non-OK status ends the chain and is reported).
use fossil,                 only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_boolean_bodies,  only : BODY_SPHERE, BODY_CYLINDER, mc_body
use fossil_boolean_harness, only : moments
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none
private
public :: chain_report_t, drill_chain, lattice_chain, unit_cube

type :: chain_report_t
   !< Outcome of a chain.
   integer(I4P)              :: steps = 0_I4P      !< Steps done.
   integer(I4P)              :: bad_status = 0_I4P !< Non-OK statuses (the chain stops at the first).
   real(R8P)                 :: max_vol = 0._R8P   !< Worst volume residual of a step.
   real(R8P)                 :: max_mom = 0._R8P   !< Worst moment residual of a step.
   real(R8P)                 :: final_vol = 0._R8P !< Volume residual of the final body against the sums.
   real(R8P)                 :: final_mom = 0._R8P !< Moment residual of the final body against the sums.
   logical                   :: valid = .true.     !< Every step watertight and manifold.
   integer(I4P)              :: facets = 0_I4P     !< Facets of the final body.
   real(R8P)                 :: seconds = 0._R8P   !< Time of the chain booleans.
endtype chain_report_t

contains
   subroutine unit_cube(s)
   !< The unit cube centred at the origin.
   type(surface_stl_object), intent(inout) :: s     !< Cube.
   type(surface_stl_object)                :: fresh !< Reset.

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
   call s%analyze
   endsubroutine unit_cube

   subroutine drill_chain(k, n, rep, verbose)
   !< The unit cube drilled by k^2 cylinders.
   integer(I4P),         intent(in)  :: k       !< Cylinders per row.
   integer(I4P),         intent(in)  :: n       !< Marching-cubes samples per direction.
   type(chain_report_t), intent(out) :: rep     !< Outcome.
   logical,              intent(in)  :: verbose !< One line per step.
   type(surface_stl_object)          :: cube, r, c, ic
   real(R8P)                         :: s, f, x, y, vr, cr(6), vi, ci(6), ve, ce(6), vref, cref
   integer(I4P)                      :: i, j, st

   call unit_cube(cube) ; r = cube
   s = 0.8_R8P / max(k - 1, 1) ; f = 0.5_R8P * s   ! radius 0.15 s: disjoint, inside the cube
   call moments(cube, [0._R8P, 0._R8P, 0._R8P], ve, ce)
   vref = ve ; cref = sum(abs(ce))
   do j = 1, k ; do i = 1, k
      x = -0.4_R8P + s * (i - 1) ; y = -0.4_R8P + s * (j - 1)
      if (k == 1) then
         x = 0._R8P ; y = 0._R8P
      endif
      call mc_body(BODY_CYLINDER, n, c)
      call c%resize(factor=vector_R8P(f, f, 1.2_R8P), respect_centroid=.true.)
      call c%translate(x=x, y=y) ; call c%analyze
      ic = cube ; call ic%boolean(other=c, op=BOOL_INTERSECT, status=st)
      if (st /= BOOL_STATUS_OK) then
         rep%bad_status = rep%bad_status + 1 ; return
      endif
      call moments(ic, [0._R8P, 0._R8P, 0._R8P], vi, ci)
      ve = ve - vi ; ce = ce - ci
      call step(r, c, BOOL_DIFFERENCE, ve, ce, vref, cref, rep, verbose)
      if (rep%bad_status > 0) return
   enddo ; enddo
   call moments(r, [0._R8P, 0._R8P, 0._R8P], vr, cr)
   rep%final_vol = abs(vr - ve) / vref ; rep%final_mom = sum(abs(cr - ce)) / cref
   rep%facets = r%get_facets_number()
   endsubroutine drill_chain

   subroutine lattice_chain(k, n, rep, verbose)
   !< The union of k^2 spheres on a grid.
   integer(I4P),         intent(in)  :: k       !< Spheres per row.
   integer(I4P),         intent(in)  :: n       !< Marching-cubes samples per direction.
   type(chain_report_t), intent(out) :: rep     !< Outcome.
   logical,              intent(in)  :: verbose !< One line per step.
   type(surface_stl_object)          :: r, p
   type(surface_stl_object), allocatable :: sph(:)
   real(R8P)                         :: s, f, vr, cr(6), vi, ci(6), ve, ce(6), vref, cref
   integer(I4P)                      :: i, j, m, st

   allocate(sph(k * k))
   s = 0.8_R8P / max(k - 1, 1) ; f = 1.2_R8P * s   ! radius 0.6 s
   do j = 1, k ; do i = 1, k
      m = i + k * (j - 1)
      call mc_body(BODY_SPHERE, n, sph(m))
      call sph(m)%resize(factor=vector_R8P(f, f, f), respect_centroid=.true.)
      call sph(m)%translate(x=-0.4_R8P + s * (i - 1), y=-0.4_R8P + s * (j - 1)) ; call sph(m)%analyze
   enddo ; enddo
   r = sph(1)
   call moments(r, [0._R8P, 0._R8P, 0._R8P], ve, ce)
   vref = ve * k * k ; cref = sum(abs(ce)) * k * k
   do m = 2, k * k
      i = mod(m - 1, k) + 1 ; j = (m - 1) / k + 1
      call moments(sph(m), [0._R8P, 0._R8P, 0._R8P], vi, ci)
      ve = ve + vi ; ce = ce + ci
      ! earlier neighbours: left (i - 1, j) and below (i, j - 1)
      if (i > 1) call overlap(sph(m - 1), sph(m))
      if (j > 1) call overlap(sph(m - k), sph(m))
      if (rep%bad_status > 0) return
      call step(r, sph(m), BOOL_UNION, ve, ce, vref, cref, rep, verbose)
      if (rep%bad_status > 0) return
   enddo
   call moments(r, [0._R8P, 0._R8P, 0._R8P], vr, cr)
   rep%final_vol = abs(vr - ve) / vref ; rep%final_mom = sum(abs(cr - ce)) / cref
   rep%facets = r%get_facets_number()

   contains
      subroutine overlap(a, b)
      !< Subtract V(a ∩ b) and its moments from the expected sums.
      type(surface_stl_object), intent(in) :: a, b

      p = a ; call p%boolean(other=b, op=BOOL_INTERSECT, status=st)
      if (st /= BOOL_STATUS_OK) then
         rep%bad_status = rep%bad_status + 1 ; return
      endif
      call moments(p, [0._R8P, 0._R8P, 0._R8P], vi, ci)
      ve = ve - vi ; ce = ce - ci
      endsubroutine overlap
   endsubroutine lattice_chain

   subroutine step(r, b, op, ve, ce, vref, cref, rep, verbose)
   !< r = op(r, b); compare with the expected volume and moments; record.
   type(surface_stl_object), intent(inout) :: r          !< Running body.
   type(surface_stl_object), intent(in)    :: b          !< Operand.
   integer(I4P),             intent(in)    :: op         !< BOOL_*.
   real(R8P),                intent(in)    :: ve, ce(6)  !< Expected volume and moments.
   real(R8P),                intent(in)    :: vref, cref !< Scales.
   type(chain_report_t),     intent(inout) :: rep        !< Outcome.
   logical,                  intent(in)    :: verbose    !< One line per step.
   real(R8P)                               :: vr, cr(6), ev, em
   integer(I8P)                            :: t0, t1, rate
   integer(I4P)                            :: st
   logical                                 :: ok

   call system_clock(t0, rate)
   call r%boolean(other=b, op=op, status=st)
   call system_clock(t1)
   rep%seconds = rep%seconds + real(t1 - t0, R8P) / rate
   rep%steps = rep%steps + 1
   if (st /= BOOL_STATUS_OK) then
      rep%bad_status = rep%bad_status + 1
      if (verbose) print '(A,I4,A,I0)', '      step ', rep%steps, ': status ', st
      return
   endif
   call moments(r, [0._R8P, 0._R8P, 0._R8P], vr, cr)
   ev = abs(vr - ve) / vref ; em = sum(abs(cr - ce)) / cref
   rep%max_vol = max(rep%max_vol, ev) ; rep%max_mom = max(rep%max_mom, em)
   ok = r%is_watertight() .and. r%is_manifold()
   rep%valid = rep%valid .and. ok
   if (verbose) print '(A,I4,A,I7,A,ES8.1,A,ES8.1,A,L1,A,F8.3,A)', '      step ', rep%steps, ': facets ', &
                      r%get_facets_number(), ', vol ', ev, ', mom ', em, ', valid ', ok, ', ', real(t1 - t0, R8P) / rate, ' s'
   endsubroutine step
endmodule fossil_boolean_csg
