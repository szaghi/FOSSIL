!< FOSSIL test support: the checks of a 3-D piece map (issue #102 step 2).
module fossil_piece_map_3d_checks
!< FOSSIL test support: the checks of a 3-D piece map (issue #102 step 2), shared by the regression tests and the
!< validation program: every piece at pose0 and at pose1 in exactly one track of its cell, numbered as the cut cells
!< of a copy moved to the pose (the track volumes equal to their pieces' volumes bit for bit), per-cell sums over the
!< tracks equal to the swept geometry (faces, wall, residual, swept = V1 - V0, outflows = residual - swept, to 1e-12 of
!< the cell on the exact paths; on the Gauss path, where a multi-track cell integrates per event interval and `sw` by
!< the step's nodes, faces and wall to 4 error_estimate, residual and outflows to that times the largest speed),
!< residual = V1 - V0 + sum of the outflows per track, every per-track residual at round-off (1e-12) on the exact
!< paths.
use fossil, only : surface_stl_object, rigid_pose_t, swept_cells_3d_t, piece_map_3d_t, cut_cells_3d_t, cut_cells_3d
use fossil_rigid_pose, only : pose_log_3d
use penf,   only : I4P, R8P

implicit none
private
public :: check_map

contains
   logical function check_map(s, sw, map, q0, q1, g, exact) result(pass)
   !< The map checks (see the program notes) against an explicit `cut_cells_3d` of copies moved to the end poses.
   type(surface_stl_object), intent(in) :: s      !< Surface (body frame).
   type(swept_cells_3d_t),   intent(in) :: sw     !< Swept geometry of the step.
   type(piece_map_3d_t),     intent(in) :: map    !< Its piece map.
   type(rigid_pose_t), intent(in) :: q0, q1 !< End poses.
   real(R8P),          intent(in) :: g(:)   !< Cell edges.
   logical,            intent(in) :: exact  !< Exact paths (per-track residuals at round-off).
   type(surface_stl_object)       :: c0, c1 !< Copies.
   type(cut_cells_3d_t)           :: k0, k1 !< Cuts.
   real(R8P)                      :: ef, ew, er, ev, eo, ei, rt, vs, sf(6), cf(6), so, w(3), h(3), area(6) !< Errors, sums.
   integer(I4P)                   :: i1, j1, l1, t, p, nn, f, bad !< Work.
   logical                        :: numbering !< Pieces numbered as the explicit cuts.
   real(R8P)                      :: tg, tv, kax(3), ang, em !< Tolerances of the per-cell sums, rotation of the step.

   c0 = s ; call c0%move(q0) ; call cut_cells_3d(c0, g, g, g, k0)
   c1 = s ; call c1%move(q1) ; call cut_cells_3d(c1, g, g, g, k1)
   numbering = k0%np == map%np0 .and. k1%np == map%np1
   bad = 0_I4P
   if (numbering) then
      do t = 1_I4P, map%nt ! volumes: the explicit pieces' volumes in piece order, bit for bit
         vs = 0._R8P
         do p = map%old_start(t), map%old_start(t + 1) - 1_I4P
            vs = vs + k0%piece_volume(map%old_piece(p))
            if (map%old_track(map%old_piece(p)) /= t) bad = bad + 1_I4P
            if (any(k0%piece_cell(:, map%old_piece(p)) /= map%track_cell(:, t))) bad = bad + 1_I4P
         enddo
         if (vs /= map%volume0(t)) bad = bad + 1_I4P
         vs = 0._R8P
         do p = map%new_start(t), map%new_start(t + 1) - 1_I4P
            vs = vs + k1%piece_volume(map%new_piece(p))
            if (map%new_track(map%new_piece(p)) /= t) bad = bad + 1_I4P
            if (any(k1%piece_cell(:, map%new_piece(p)) /= map%track_cell(:, t))) bad = bad + 1_I4P
         enddo
         if (vs /= map%volume1(t)) bad = bad + 1_I4P
      enddo
      numbering = bad == 0_I4P .and. map%old_start(map%nt + 1) == k0%np + 1 .and. map%new_start(map%nt + 1) == k1%np + 1
   endif
   tg = 1e-12_R8P ; tv = 1e-12_R8P
   if (.not. exact) then ! the per-interval integration against the step's Gauss nodes: the quadrature error
      call pose_log_3d(q0, q1, kax, ang)
      tg = max(tg, 4._R8P * sw%error_estimate)
      tv = max(tv, tg * (norm2(q1%t - q0%t) + ang * (maxval(abs(g)) + norm2(q0%t))))
   endif
   ef = 0._R8P ; ew = 0._R8P ; er = 0._R8P ; ev = 0._R8P ; eo = 0._R8P ; ei = 0._R8P ; rt = 0._R8P ; em = 0._R8P
   do l1 = 1, size(g) - 1
      do j1 = 1, size(g) - 1
         do i1 = 1, size(g) - 1
            h = [g(i1 + 1) - g(i1), g(j1 + 1) - g(j1), g(l1 + 1) - g(l1)]
            area = [h(2) * h(3), h(2) * h(3), h(3) * h(1), h(3) * h(1), h(1) * h(2), h(1) * h(2)]
            sf = 0._R8P ; w = 0._R8P ; so = 0._R8P ; vs = 0._R8P ; rt = 0._R8P
            do t = map%first_track(i1, j1, l1), map%first_track(i1, j1, l1) + map%tracks(i1, j1, l1) - 1_I4P
               sf = sf + map%faces(:, t) ; w = w + map%wall(:, t) ; rt = rt + map%residual(t)
               vs = vs + map%volume1(t) - map%volume0(t)
               so = so + sum(map%outflow(:, t))
               ei = max(ei, abs(map%residual(t) - (map%volume1(t) - map%volume0(t) + sum(map%outflow(:, t)))))
               if (exact .and. map%status(t) /= -1_I4P) eo = max(eo, abs(map%residual(t)))
               if (map%tracks(i1, j1, l1) > 1_I4P) em = max(em, abs(map%residual(t)))
            enddo
            if (map%tracks(i1, j1, l1) == 0_I4P) cycle
            cf = [sw%face_x(i1, j1, l1), sw%face_x(i1 + 1, j1, l1), sw%face_y(i1, j1, l1), sw%face_y(i1, j1 + 1, l1), &
                  sw%face_z(i1, j1, l1), sw%face_z(i1, j1, l1 + 1)]
            do f = 1, 6
               ef = max(ef, abs(sf(f) - cf(f)) / merge(area(f), 1._R8P, exact))
            enddo
            ew = max(ew, maxval(abs(w - sw%wall(:, i1, j1, l1))))
            er = max(er, abs(rt - sw%residual(i1, j1, l1)), abs(so - (sw%residual(i1, j1, l1) - sw%swept(i1, j1, l1))))
            ev = max(ev, abs(vs - sw%swept(i1, j1, l1)))
         enddo
      enddo
   enddo
   nn = map%nt
   pass = numbering .and. ef <= tg .and. ew <= 2._R8P * tg .and. er <= 6._R8P * tv .and. ev <= 1e-12_R8P .and. &
          ei <= 1e-12_R8P
   if (exact) pass = pass .and. eo <= 1e-12_R8P
   print '(A,L1,A,I0,A,6ES9.1)', '    numbering as explicit cuts ', numbering, '  tracks ', nn, &
         '  per-cell sums - sw: faces, wall, residual/outflow, swept; residual identity; |per-track residual| ', ef, ew, er, &
         ev, ei, eo
   print '(A,ES9.1,A,ES9.1)', '    max |residual| of the multi-track cells'' tracks ', em, '  per-cell tolerance ', tg
   endfunction check_map
endmodule fossil_piece_map_3d_checks
