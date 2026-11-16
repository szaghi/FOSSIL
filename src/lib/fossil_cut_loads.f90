!< FOSSIL, cut-cell loads: force and moment on a body from the wall terms of its cut cells (#109 WP2 step 05).
module fossil_cut_loads
!< FOSSIL, cut-cell loads: force and moment on a body from the wall terms of its cut cells (#109 WP2 step 05).
!<
!< A finite-volume solver on cut cells knows the pressure (and the viscous traction) per cut cell, or per fluid piece; the
!< force on the body is then the sum over the wall terms of the cells. The terms come from the cut itself:
!<
!< - 3-D: `cut_cells_3d(..., walls=)` (`cut_walls_3d_t`): one term per (cell, facet) whose clipped polygon has 3 vertices
!<   or more, with its vector area A_w (solid to fluid: the body's outward normal times the area), its area |A_w| and its
!<   centroid x_w, labelled by piece and facet.
!< - 2-D: `cut_walls_2d(cut, walls)` (`cut_walls_2d_t`) lists the wall segments of the pieces' oriented boundaries of
!<   `curve%cut_cells(..., boundaries=.true.)` (`piece_seg`, the fluid on the left; cell-edge segments left out), in piece
!<   order then boundary order, each with its piece, its curve segment (`piece_seg_src`), its end points, its vector area
!<   A_w = (-(b_v - a_v), b_u - a_u) (the left normal times the length: out of the solid, into the fluid), length and
!<   midpoint.
!<
!< `cut_loads_3d` / `cut_loads_2d` sum, as `integrate_loads` sums the load samples (the conventions L1, L4 of #109 WP2):
!<
!<   f_w = -p_w A_w + tau_w |A_w|,   F = sum_w f_w,   M = sum_w (x_w - c) x f_w  [- (K_w g_w) x A_w],
!<
!< p_w the pressure at the wall's centroid (`pressure(n)`, or `pressure_per_piece(np)` spread to the walls of each piece),
!< tau_w the viscous traction (`shear(d, n)`), c = `about` (default the origin). In 2-D both are per unit depth and the
!< moment is the scalar r_u f_v - r_v f_u. The sums are compensated (Neumaier) in wall order, per patch and per body (the
!< facet's or segment's, given the surface or curve), deterministic.
!<
!< Exactness. The walls are flat, so the centroid value integrates a linear pressure exactly for the force: in p = p0 + g . x
!< the force on a body inside the grid is -g V to round-off. The moment's integrand (x - c) p is then quadratic: the
!< centroid value misses int (x - x_w)(x - x_w)^T dA g = |A_w| K_w g. With the pressure's gradient per wall (`gradient`)
!< the moment adds the term -(K_w g_w) x A_w and is exact for a pressure linear on each wall: in 2-D K_w = d d^T / 12 from
!< the segment d = b - a, in 3-D from `cut_walls_gyration` (the wall polygons clipped again by the cut's own clip,
!< `cut_wall_polygon_3d`, and their exact second moments about the stored centroid). Without `gradient` the pressure is
!< taken constant on each wall (a cell value), which the moment integrates exactly too.
!<
!< CPU only; free routines, re-exported by `use fossil`.
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_cut_cells_3d,       only : cut_walls_3d_t, cut_wall_polygon_3d
use fossil_facet_object,       only : facet_object
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use iso_fortran_env,           only : stderr => error_unit
use penf,                      only : I4P, R8P

implicit none
private
public :: cut_walls_2d_t, cut_walls_2d, cut_walls_gyration, cut_loads_3d, cut_loads_2d, cut_loads_linear_pressure

type :: cut_walls_2d_t
   !< The wall segments of a 2-D cut (see the module header), in piece order, each piece's in the order of its oriented
   !< boundary (`piece_seg`).
   integer(I4P)              :: n = 0_I4P     !< Wall segments.
   integer(I4P), allocatable :: piece(:)      !< Piece of each segment [n].
   integer(I4P), allocatable :: segment(:)    !< Curve segment it lies on (its first vertex, `piece_seg_src`) [n].
   real(R8P),    allocatable :: ends(:,:)     !< End points (a, b), the fluid on the left [4, n].
   real(R8P),    allocatable :: area(:,:)     !< Vector area: (-(b_v - a_v), b_u - a_u), out of the solid [2, n].
   real(R8P),    allocatable :: length(:)     !< Length [n].
   real(R8P),    allocatable :: centroid(:,:) !< Midpoint [2, n].
endtype cut_walls_2d_t

contains
   subroutine cut_walls_2d(cut, walls, status)
   !< The wall segments of a 2-D cut made with `boundaries=.true.` (see the module header). `status =
   !< STATUS_INVALID_INPUT` (stderr message, no segment) if the cut has no oriented boundaries.
   type(cut_cells_2d_t), intent(in)            :: cut    !< 2-D cut cells (with boundaries).
   type(cut_walls_2d_t), intent(out)           :: walls  !< Wall segments.
   integer(I4P),         intent(out), optional :: status !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: p, k, n, w                            !< Work.
   real(R8P)    :: du, dv                                !< Segment.

   if (present(status)) status = STATUS_OK
   if (.not. (allocated(cut%piece_seg_start) .and. allocated(cut%piece_seg) .and. allocated(cut%piece_seg_wall) .and. &
              allocated(cut%piece_seg_src))) then
      write(stderr, '(A)') 'fossil cut_walls_2d: refused -- the cut has no boundaries (cut_cells(..., boundaries=.true.))'
      if (present(status)) status = STATUS_INVALID_INPUT
      call walls_2d_alloc(walls, 0_I4P)
      return
   endif
   n = int(count(cut%piece_seg_wall), I4P)
   call walls_2d_alloc(walls, n)
   w = 0_I4P
   do p = 1_I4P, cut%np
      do k = cut%piece_seg_start(p), cut%piece_seg_start(p + 1) - 1_I4P
         if (.not. cut%piece_seg_wall(k)) cycle
         w = w + 1_I4P
         walls%piece(w) = p ; walls%segment(w) = cut%piece_seg_src(k) ; walls%ends(:, w) = cut%piece_seg(1:4, k)
         du = cut%piece_seg(3, k) - cut%piece_seg(1, k) ; dv = cut%piece_seg(4, k) - cut%piece_seg(2, k)
         walls%area(1, w) = -dv ; walls%area(2, w) = du
         walls%length(w) = sqrt(du * du + dv * dv)
         walls%centroid(1, w) = 0.5_R8P * (cut%piece_seg(1, k) + cut%piece_seg(3, k))
         walls%centroid(2, w) = 0.5_R8P * (cut%piece_seg(2, k) + cut%piece_seg(4, k))
      enddo
   enddo
   endsubroutine cut_walls_2d

   subroutine cut_walls_gyration(walls, surface, x, y, z, gyration, status)
   !< The normalised second moments K_w = int (x - x_w)(x - x_w)^T dA / |A_w| of the 3-D wall terms, about their stored
   !< centroids x_w, as gyration(1:6, w) = (K_xx, K_yy, K_zz, K_xy, K_yz, K_zx): each term's polygon clipped again from its
   !< facet and cell by the cut's own clip (`cut_wall_polygon_3d`), then its fans' exact second moments (int r_i r_j dA =
   !< A_t / 12 (sum_k r_ki r_kj + s_i s_j), s the sum of the fan's three vertices). The surface and the cell edges must be
   !< those of the cut. `status = STATUS_INVALID_INPUT` (stderr message, gyration of size [6, 0]) if the edges do not match
   !< the terms' cells, a facet is out of range or a polygon has fewer than 3 vertices (another surface or grid).
   type(cut_walls_3d_t),     intent(in)            :: walls       !< Wall terms of `cut_cells_3d(..., walls=)` (host).
   type(surface_stl_object), intent(in)            :: surface     !< The surface of the cut.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< The cell edges of the cut.
   real(R8P), allocatable,   intent(out)           :: gyration(:,:) !< K_w [6, n].
   integer(I4P),             intent(out), optional :: status      !< STATUS_OK or STATUS_INVALID_INPUT.
   type(facet_object), pointer :: fp                              !< Facet.
   real(R8P)    :: poly(3, 16), t(3, 3), c0(3), c1(3), xb(3), a(3), b(3), cc(3), s(3), e1(3), e2(3), nn(3), j6(6), at, ta
   integer(I4P) :: nx, ny, nz, c, i, j, k, w, f, np, m, nf, st   !< Work.

   st = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   nf = surface%get_facets_number()
   if (nx < 1_I4P .or. ny < 1_I4P .or. nz < 1_I4P .or. .not. allocated(walls%cell_start)) then
      st = STATUS_INVALID_INPUT
   elseif (size(walls%cell_start, kind=I4P) /= nx * ny * nz + 1_I4P) then
      st = STATUS_INVALID_INPUT
   endif
   allocate(gyration(6, max(0_I4P, walls%n))) ; gyration = 0._R8P
   if (st == STATUS_OK) then
      cells: do c = 1_I4P, nx * ny * nz
         if (walls%cell_start(c + 1) <= walls%cell_start(c)) cycle
         i = mod(c - 1_I4P, nx) + 1_I4P ; j = mod((c - 1_I4P) / nx, ny) + 1_I4P ; k = (c - 1_I4P) / (nx * ny) + 1_I4P
         c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
         do w = walls%cell_start(c), walls%cell_start(c + 1) - 1_I4P
            f = walls%facet(w)
            if (f < 1_I4P .or. f > nf) then
               st = STATUS_INVALID_INPUT
               exit cells
            endif
            fp => surface%facet_at(f)
            do m = 1, 3
               t(:, m) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
            enddo
            call cut_wall_polygon_3d(t, c0, c1, poly, np)
            if (np < 3_I4P) then
               st = STATUS_INVALID_INPUT
               exit cells
            endif
            xb = walls%centroid(:, w) ; j6 = 0._R8P ; ta = 0._R8P
            do m = 2, np - 1
               a = poly(:, 1) - xb ; b = poly(:, m) - xb ; cc = poly(:, m + 1) - xb
               e1 = b - a ; e2 = cc - a
               nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
               at = 0.5_R8P * sqrt(nn(1) * nn(1) + nn(2) * nn(2) + nn(3) * nn(3))
               s = a + b + cc ; ta = ta + at
               j6(1) = j6(1) + at * (a(1) * a(1) + b(1) * b(1) + cc(1) * cc(1) + s(1) * s(1)) / 12._R8P
               j6(2) = j6(2) + at * (a(2) * a(2) + b(2) * b(2) + cc(2) * cc(2) + s(2) * s(2)) / 12._R8P
               j6(3) = j6(3) + at * (a(3) * a(3) + b(3) * b(3) + cc(3) * cc(3) + s(3) * s(3)) / 12._R8P
               j6(4) = j6(4) + at * (a(1) * a(2) + b(1) * b(2) + cc(1) * cc(2) + s(1) * s(2)) / 12._R8P
               j6(5) = j6(5) + at * (a(2) * a(3) + b(2) * b(3) + cc(2) * cc(3) + s(2) * s(3)) / 12._R8P
               j6(6) = j6(6) + at * (a(3) * a(1) + b(3) * b(1) + cc(3) * cc(1) + s(3) * s(1)) / 12._R8P
            enddo
            if (ta > 0._R8P) gyration(:, w) = j6 / ta
         enddo
      enddo cells
   endif
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil cut_walls_gyration: refused -- the terms do not match the surface and the cell edges'
      deallocate(gyration) ; allocate(gyration(6, 0))
   endif
   if (present(status)) status = st
   endsubroutine cut_walls_gyration

   subroutine cut_loads_3d(walls, force, moment, status, pressure, pressure_per_piece, shear, gradient, gyration, about, &
                           surface, per_patch_force, per_patch_moment, per_body_force, per_body_moment, area, unassigned)
   !< Force and moment on the body from the wall terms of its 3-D cut (see the module header): f_w = -p_w A_w + tau_w |A_w|,
   !< F = sum_w f_w, M = sum_w (x_w - about) x f_w, minus (K_w g_w) x A_w with `gradient` (g_w, the pressure's gradient on
   !< wall w) and `gyration` (K_w from `cut_walls_gyration`), in wall order with Neumaier compensation.
   !<
   !< The pressure at the walls' centroids: `pressure(n)`, or `pressure_per_piece(np)` (the value of each fluid piece on its
   !< walls; np >= the largest piece of the terms); at most one of them, and a pressure and/or `shear(3, n)`. With
   !< `pressure_per_piece` a term of no piece (piece 0: a cell without trusted pieces) takes no pressure: their count is
   !< `unassigned` (and one stderr line). Per patch and per body (the facet's, `surface` required): allocatable [3, np],
   !< [3, nb], np and nb the largest ids of the terms' facets. `area` = sum_w |A_w|. `status = STATUS_INVALID_INPUT`
   !< (stderr message, zero outputs, per-patch / per-body arrays of size 0) on a missing or doubled pressure, a size that
   !< differs, `gradient` without `gyration` or without a pressure, or per patch / per body without the surface.
   type(cut_walls_3d_t),     intent(in)                      :: walls                 !< Wall terms (host).
   real(R8P),                intent(out)                     :: force(:)              !< Force [3].
   real(R8P),                intent(out)                     :: moment(:)             !< Moment about `about` [3].
   integer(I4P),             intent(out),           optional :: status                !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                intent(in),            optional :: pressure(:)           !< Pressure at the walls [n].
   real(R8P),                intent(in),            optional :: pressure_per_piece(:) !< Pressure per piece [np].
   real(R8P),                intent(in),            optional :: shear(:,:)            !< Viscous traction [3, n].
   real(R8P),                intent(in),            optional :: gradient(:,:)         !< Pressure gradient per wall [3, n].
   real(R8P),                intent(in),            optional :: gyration(:,:)         !< K_w [6, n].
   real(R8P),                intent(in),            optional :: about(:)              !< Point of the moments [3].
   type(surface_stl_object), intent(in),            optional :: surface               !< The surface (patches, bodies).
   real(R8P), allocatable,   intent(out),           optional :: per_patch_force(:,:)  !< Force per patch [3, np].
   real(R8P), allocatable,   intent(out),           optional :: per_patch_moment(:,:) !< Moment per patch [3, np].
   real(R8P), allocatable,   intent(out),           optional :: per_body_force(:,:)   !< Force per body [3, nb].
   real(R8P), allocatable,   intent(out),           optional :: per_body_moment(:,:)  !< Moment per body [3, nb].
   real(R8P),                intent(out),           optional :: area                  !< Wetted area sum_w |A_w|.
   integer(I4P),             intent(out),           optional :: unassigned            !< Terms of no piece (per piece p).
   real(R8P),    allocatable :: ctr(:,:), tot(:), pp(:,:), pb(:,:)                    !< Contributions, totals.
   integer(I4P), allocatable :: idp(:), idb(:)                                         !< Patch, body ids.
   real(R8P)    :: ab(3), pw, f(3), m(3), kg(3), r(3), av(3)                         !< Work.
   integer(I4P) :: st, n, w, nun                                                       !< Work.
   logical      :: wp, wpp, ws, wg                                                     !< Inputs.

   n = max(0_I4P, walls%n) ; ab = 0._R8P ; nun = 0_I4P
   wp = present(pressure) ; wpp = present(pressure_per_piece) ; ws = present(shear) ; wg = present(gradient)
   st = STATUS_OK
   if (wp .and. wpp) st = STATUS_INVALID_INPUT
   if (.not. (wp .or. wpp .or. ws)) st = STATUS_INVALID_INPUT
   if (size(force) /= 3 .or. size(moment) /= 3) st = STATUS_INVALID_INPUT
   if (n > 0_I4P .and. .not. (allocated(walls%area) .and. allocated(walls%centroid) .and. allocated(walls%wet_area) .and. &
                              allocated(walls%piece) .and. allocated(walls%facet))) st = STATUS_INVALID_INPUT
   if (wp) then
      if (size(pressure, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (wpp .and. n > 0_I4P .and. st == STATUS_OK) then
      if (size(pressure_per_piece, kind=I4P) < maxval(walls%piece(1:n))) st = STATUS_INVALID_INPUT
   endif
   if (ws) then
      if (size(shear, 1) /= 3 .or. size(shear, 2, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (wg) then
      if (.not. (wp .or. wpp) .or. .not. present(gyration)) then
         st = STATUS_INVALID_INPUT
      elseif (size(gradient, 1) /= 3 .or. size(gradient, 2, kind=I4P) /= n .or. size(gyration, 1) /= 6 .or. &
              size(gyration, 2, kind=I4P) /= n) then
         st = STATUS_INVALID_INPUT
      endif
   endif
   if (present(about)) then
      if (size(about) /= 3) then
         st = STATUS_INVALID_INPUT
      else
         ab = about
      endif
   endif
   if ((present(per_patch_force) .or. present(per_patch_moment) .or. present(per_body_force) .or. &
        present(per_body_moment)) .and. .not. present(surface)) st = STATUS_INVALID_INPUT
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil cut_loads_3d: refused -- one pressure (per wall or per piece) and/or shear of the '// &
                           'terms'' sizes; gradient with gyration; per patch / per body with the surface'
      allocate(tot(7), pp(7, 0), pb(7, 0)) ; tot = 0._R8P
   else
      allocate(ctr(7, n), idp(n), idb(n)) ; idp = 0_I4P ; idb = 0_I4P
      do w = 1_I4P, n
         av = walls%area(:, w) ; f = 0._R8P ; m = 0._R8P
         pw = 0._R8P
         if (wp) then
            pw = pressure(w)
         elseif (wpp) then
            if (walls%piece(w) >= 1_I4P) then
               pw = pressure_per_piece(walls%piece(w))
            else
               nun = nun + 1_I4P
            endif
         endif
         if (wp .or. wpp) f = -(pw * av)
         if (ws) f = f + shear(:, w) * walls%wet_area(w)
         r = walls%centroid(:, w) - ab
         m = cross3(r, f)
         if (wg) then
            kg = sym_times(gyration(:, w), gradient(:, w))
            m = m - cross3(kg, av)
         endif
         ctr(1:3, w) = f ; ctr(4:6, w) = m ; ctr(7, w) = walls%wet_area(w)
         if (present(surface)) then
            idp(w) = surface%facet_patch(walls%facet(w)) ; idb(w) = surface%facet_body(walls%facet(w))
         endif
      enddo
      call reduce(7_I4P, n, ctr, idp, idb, tot, pp, pb)
      if (nun > 0_I4P) write(stderr, '(A,I0,A)') 'fossil cut_loads_3d: ', nun, ' wall terms of no piece took no pressure'
   endif
   force = 0._R8P ; moment = 0._R8P
   if (size(force) == 3) force = tot(1:3)
   if (size(moment) == 3) moment = tot(4:6)
   if (present(per_patch_force)) per_patch_force = pp(1:3, :)
   if (present(per_patch_moment)) per_patch_moment = pp(4:6, :)
   if (present(per_body_force)) per_body_force = pb(1:3, :)
   if (present(per_body_moment)) per_body_moment = pb(4:6, :)
   if (present(area)) area = tot(7)
   if (present(unassigned)) unassigned = nun
   if (present(status)) status = st
   endsubroutine cut_loads_3d

   subroutine cut_loads_2d(walls, force, moment, status, pressure, pressure_per_piece, shear, gradient, about, curve, &
                           per_patch_force, per_patch_moment, per_body_force, per_body_moment, area)
   !< 2-D twin of `cut_loads_3d` per unit depth on the wall segments of `cut_walls_2d`: shear(2, n), force(2), about(2),
   !< the scalar moment M = sum_w (r_u f_v - r_v f_u), r = x_w - about; with `gradient(2, n)` the term -(K_w g_w) x A_w,
   !< K_w = d d^T / 12 from the segment d = b - a (no gyration input). Per patch and per body (the curve segment's,
   !< `curve` required): per_patch_force(2, np), per_patch_moment(np), per_body_force(2, nb), per_body_moment(nb).
   !< `area` = the wetted length sum_w |A_w|. Every segment has a piece. Refusals as in 3-D.
   type(cut_walls_2d_t),   intent(in)            :: walls                 !< Wall segments.
   real(R8P),              intent(out)           :: force(:)              !< Force [2].
   real(R8P),              intent(out)           :: moment                !< Moment about `about`.
   integer(I4P),           intent(out), optional :: status                !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),              intent(in),  optional :: pressure(:)           !< Pressure at the walls [n].
   real(R8P),              intent(in),  optional :: pressure_per_piece(:) !< Pressure per piece [np].
   real(R8P),              intent(in),  optional :: shear(:,:)            !< Viscous traction [2, n].
   real(R8P),              intent(in),  optional :: gradient(:,:)         !< Pressure gradient per wall [2, n].
   real(R8P),              intent(in),  optional :: about(:)              !< Point of the moments [2].
   type(curve_object),     intent(in),  optional :: curve                 !< The curve (patches, bodies).
   real(R8P), allocatable, intent(out), optional :: per_patch_force(:,:)  !< Force per patch [2, np].
   real(R8P), allocatable, intent(out), optional :: per_patch_moment(:)   !< Moment per patch [np].
   real(R8P), allocatable, intent(out), optional :: per_body_force(:,:)   !< Force per body [2, nb].
   real(R8P), allocatable, intent(out), optional :: per_body_moment(:)    !< Moment per body [nb].
   real(R8P),              intent(out), optional :: area                  !< Wetted length sum_w |A_w|.
   real(R8P),    allocatable :: ctr(:,:), tot(:), pp(:,:), pb(:,:)        !< Contributions, totals.
   integer(I4P), allocatable :: idp(:), idb(:), sbody(:)                   !< Patch, body ids; segment bodies.
   real(R8P)    :: ab(2), pw, f(2), d(2), kg(2), av(2), r(2), m        !< Work.
   integer(I4P) :: st, n, w, l, i0, i1, sg                                 !< Work.
   logical      :: wp, wpp, ws, wg                                         !< Inputs.

   n = max(0_I4P, walls%n) ; ab = 0._R8P
   wp = present(pressure) ; wpp = present(pressure_per_piece) ; ws = present(shear) ; wg = present(gradient)
   st = STATUS_OK
   if (wp .and. wpp) st = STATUS_INVALID_INPUT
   if (.not. (wp .or. wpp .or. ws)) st = STATUS_INVALID_INPUT
   if (size(force) /= 2) st = STATUS_INVALID_INPUT
   if (n > 0_I4P .and. .not. (allocated(walls%area) .and. allocated(walls%centroid) .and. allocated(walls%length) .and. &
                              allocated(walls%piece) .and. allocated(walls%ends))) st = STATUS_INVALID_INPUT
   if (wp) then
      if (size(pressure, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (wpp .and. n > 0_I4P .and. st == STATUS_OK) then
      if (size(pressure_per_piece, kind=I4P) < maxval(walls%piece(1:n))) st = STATUS_INVALID_INPUT
   endif
   if (ws) then
      if (size(shear, 1) /= 2 .or. size(shear, 2, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (wg) then
      if (.not. (wp .or. wpp)) then
         st = STATUS_INVALID_INPUT
      elseif (size(gradient, 1) /= 2 .or. size(gradient, 2, kind=I4P) /= n) then
         st = STATUS_INVALID_INPUT
      endif
   endif
   if (present(about)) then
      if (size(about) /= 2) then
         st = STATUS_INVALID_INPUT
      else
         ab = about
      endif
   endif
   if ((present(per_patch_force) .or. present(per_patch_moment) .or. present(per_body_force) .or. &
        present(per_body_moment)) .and. .not. present(curve)) st = STATUS_INVALID_INPUT
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil cut_loads_2d: refused -- one pressure (per wall or per piece) and/or shear of the '// &
                           'segments'' sizes; gradient with a pressure; per patch / per body with the curve'
      allocate(tot(4), pp(4, 0), pb(4, 0)) ; tot = 0._R8P
   else
      allocate(ctr(4, n), idp(n), idb(n)) ; idp = 0_I4P ; idb = 0_I4P
      if (present(curve)) then
         allocate(sbody(max(1_I4P, curve%get_vertices_number()))) ; sbody = 0_I4P
         do l = 1_I4P, curve%get_loops_number()
            call curve%get_loop_range(l, i0, i1)
            if (i1 >= i0) sbody(i0:i1) = curve%loop_body(l)
         enddo
      endif
      do w = 1_I4P, n
         av = walls%area(:, w) ; f = 0._R8P
         pw = 0._R8P
         if (wp) pw = pressure(w)
         if (wpp) pw = pressure_per_piece(walls%piece(w))
         if (wp .or. wpp) f = -(pw * av)
         if (ws) f = f + shear(:, w) * walls%length(w)
         r = walls%centroid(:, w) - ab
         m = r(1) * f(2) - r(2) * f(1)
         if (wg) then
            d = walls%ends(3:4, w) - walls%ends(1:2, w)
            kg = d * ((d(1) * gradient(1, w) + d(2) * gradient(2, w)) / 12._R8P)
            m = m - (kg(1) * av(2) - kg(2) * av(1))
         endif
         ctr(1:2, w) = f ; ctr(3, w) = m ; ctr(4, w) = walls%length(w)
         if (present(curve)) then
            sg = walls%segment(w)
            if (sg >= 1_I4P .and. sg <= size(sbody, kind=I4P)) then
               idp(w) = curve%segment_patch(sg) ; idb(w) = sbody(sg)
            endif
         endif
      enddo
      call reduce(4_I4P, n, ctr, idp, idb, tot, pp, pb)
   endif
   force = 0._R8P
   if (size(force) == 2) force = tot(1:2)
   moment = tot(3)
   if (present(per_patch_force)) per_patch_force = pp(1:2, :)
   if (present(per_patch_moment)) per_patch_moment = pp(3, :)
   if (present(per_body_force)) per_body_force = pb(1:2, :)
   if (present(per_body_moment)) per_body_moment = pb(3, :)
   if (present(area)) area = tot(4)
   if (present(status)) status = st
   endsubroutine cut_loads_2d

   pure function cut_loads_linear_pressure(points, p0, gradient) result(p)
   !< The linear pressure p = p0 + g . x at points(d, n) (the walls' centroids), the products added in axis order: the
   !< field of the buoyancy checks (a fluid at rest, g = rho times minus the gravity), with no callback.
   real(R8P), intent(in) :: points(:,:)           !< Points [d, n].
   real(R8P), intent(in) :: p0                    !< Pressure at the origin.
   real(R8P), intent(in) :: gradient(:)           !< Gradient [d].
   real(R8P)             :: p(size(points, 2))    !< Pressure [n].
   real(R8P)             :: s                     !< Sum.
   integer(I4P)          :: w, a                  !< Counters.

   do w = 1, size(points, 2)
      s = 0._R8P
      do a = 1, min(size(points, 1), size(gradient))
         s = s + gradient(a) * points(a, w)
      enddo
      p(w) = p0 + s
   enddo
   endfunction cut_loads_linear_pressure

   ! private

   subroutine walls_2d_alloc(walls, n)
   !< Allocate the wall segments' arrays.
   type(cut_walls_2d_t), intent(inout) :: walls !< Wall segments.
   integer(I4P),         intent(in)    :: n     !< Segments.

   walls%n = n
   allocate(walls%piece(n), walls%segment(n), walls%ends(4, n), walls%area(2, n), walls%length(n), walls%centroid(2, n))
   endsubroutine walls_2d_alloc

   subroutine reduce(nc, n, ctr, idp, idb, tot, pp, pb)
   !< The compensated (Neumaier) sums of the contributions in wall order: all, per patch, per body (ids below 1 enter the
   !< totals only).
   integer(I4P),           intent(in)  :: nc, n          !< Components, walls.
   real(R8P),              intent(in)  :: ctr(nc, n)     !< Contributions.
   integer(I4P),           intent(in)  :: idp(n), idb(n) !< Patch and body ids.
   real(R8P), allocatable, intent(out) :: tot(:), pp(:,:), pb(:,:) !< Totals.
   real(R8P), allocatable :: ct(:), cp(:,:), cb(:,:)     !< Compensations.
   integer(I4P)           :: np, nb, s, c, ip, ib        !< Work.

   np = 0_I4P ; nb = 0_I4P
   if (n > 0_I4P) then
      np = max(0_I4P, maxval(idp)) ; nb = max(0_I4P, maxval(idb))
   endif
   allocate(tot(nc), pp(nc, np), pb(nc, nb), ct(nc), cp(nc, np), cb(nc, nb))
   tot = 0._R8P ; pp = 0._R8P ; pb = 0._R8P ; ct = 0._R8P ; cp = 0._R8P ; cb = 0._R8P
   do s = 1_I4P, n
      ip = idp(s) ; ib = idb(s)
      do c = 1_I4P, nc
         call neumaier_add(tot(c), ct(c), ctr(c, s))
         if (ip >= 1_I4P) call neumaier_add(pp(c, ip), cp(c, ip), ctr(c, s))
         if (ib >= 1_I4P) call neumaier_add(pb(c, ib), cb(c, ib), ctr(c, s))
      enddo
   enddo
   tot = tot + ct ; pp = pp + cp ; pb = pb + cb
   endsubroutine reduce

   pure subroutine neumaier_add(sm, cm, x)
   !< Add x to the running sum sm with the compensation cm (Neumaier): the total is sm + cm.
   real(R8P), intent(inout) :: sm, cm !< Sum, compensation.
   real(R8P), intent(in)    :: x      !< Term.
   real(R8P) :: t                     !< New sum.

   t = sm + x
   if (abs(sm) >= abs(x)) then
      cm = cm + ((sm - t) + x)
   else
      cm = cm + ((x - t) + sm)
   endif
   sm = t
   endsubroutine neumaier_add

   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c(1) = a(2) * b(3) - a(3) * b(2) ; c(2) = a(3) * b(1) - a(1) * b(3) ; c(3) = a(1) * b(2) - a(2) * b(1)
   endfunction cross3

   pure function sym_times(k6, g) result(kg)
   !< The symmetric matrix k6 = (xx, yy, zz, xy, yz, zx) times g.
   real(R8P), intent(in) :: k6(6), g(3) !< Matrix, vector.
   real(R8P)             :: kg(3)       !< Product.

   kg(1) = k6(1) * g(1) + k6(4) * g(2) + k6(6) * g(3)
   kg(2) = k6(4) * g(1) + k6(2) * g(2) + k6(5) * g(3)
   kg(3) = k6(6) * g(1) + k6(5) * g(2) + k6(3) * g(3)
   endfunction sym_times
endmodule fossil_cut_loads
