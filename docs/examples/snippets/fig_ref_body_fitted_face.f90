!< Figure of the reference page "Body-fitted grid support": a 9 x 9 squircle face a little below the top cap of a
!< cylinder and a little wider than it, put on the body twice: by project_block_face (boundary on the rim, interior on
!< the cap) and with tolerance 0 (every vertex on the patch of its nearest facet: the boundary slides down the side).
program fig_ref_body_fitted_face
use fossil,              only : surface_stl_object, projection_target_t, curvilinear_block_t, project_block_face, &
                                PROJ_CURVE, PROJ_PATCH
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_zoom, fig_surface, fig_segments, fig_label, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
integer(I4P), parameter         :: NS = 32
type(surface_stl_object)        :: cyl
type(projection_target_t)       :: target
type(curvilinear_block_t)       :: constrained, free
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: b0, b1, t0, t1
real(R8P)                       :: x(9, 9, 2), y(9, 9, 2), z(9, 9, 2), s, t
integer(I4P), allocatable       :: kinds(:,:)
integer(I4P)                    :: i, j, m, u

allocate(f(26 * NS)) ; m = 0 ! a cylinder of radius 1 and height 2, patches side, bottom, top (the top in 12 rings)
do i = 0, NS - 1
   b0 = vector_R8P(cos(2 * PI * i / NS), sin(2 * PI * i / NS), 0._R8P)
   b1 = vector_R8P(cos(2 * PI * mod(i + 1, NS) / NS), sin(2 * PI * mod(i + 1, NS) / NS), 0._R8P)
   t0 = b0 ; t0%z = 2._R8P ; t1 = b1 ; t1%z = 2._R8P
   call add(b0, b1, t1) ; call add(b0, t1, t0) ; call add(vector_R8P(0._R8P, 0._R8P, 0._R8P), b1, b0)
enddo
do i = 0, NS - 1
   do j = 1, 12
      t0 = vector_R8P(cos(2 * PI * i / NS), sin(2 * PI * i / NS), 0._R8P)
      t1 = vector_R8P(cos(2 * PI * mod(i + 1, NS) / NS), sin(2 * PI * mod(i + 1, NS) / NS), 0._R8P)
      b0 = t0 * ((j - 1) / 12._R8P) ; b1 = t1 * ((j - 1) / 12._R8P) ; t0 = t0 * (j / 12._R8P) ; t1 = t1 * (j / 12._R8P)
      b0%z = 2._R8P ; b1%z = 2._R8P ; t0%z = 2._R8P ; t1%z = 2._R8P
      call add(b0, t0, t1)
      if (j > 1) call add(b0, t1, b1)
   enddo
enddo
call cyl%adopt_facets(facets=f)
call cyl%set_patch([(3 * i + 1, 3 * i + 2, i = 0, NS - 1)], 'side')
call cyl%set_patch([(3 * i + 3, i = 0, NS - 1)], 'bottom')
call cyl%set_patch([(i, i = 3 * NS + 1, m)], 'top')
call cyl%analyze
call cyl%projection_target(target)
do j = 1, 9 ! the squircle: the disc map of the square, radius 1.04, 0.05 below the cap
   do i = 1, 9
      s = -1._R8P + (i - 1) / 4._R8P ; t = -1._R8P + (j - 1) / 4._R8P
      x(i, j, :) = 1.04_R8P * s * sqrt(1._R8P - t**2 / 2) ; y(i, j, :) = 1.04_R8P * t * sqrt(1._R8P - s**2 / 2)
      z(i, j, :) = [1.95_R8P, 1.5_R8P]
   enddo
enddo
call constrained%initialize(x, y, z) ; free = constrained
call project_block_face(constrained, face=5_I4P, target=target, tolerance=0.07_R8P, kinds=kinds)
print '(A,I0,A,I0,A,F7.4)', 'constrained: on the rim ', count(kinds == PROJ_CURVE), ', on the cap ', &
                            count(kinds == PROJ_PATCH), ', lowest z ', minval(constrained%z(:, :, 1))
call project_block_face(free, face=5_I4P, target=target, tolerance=0._R8P, kinds=kinds)
print '(A,I0,A,F7.4)', 'tolerance 0: on a patch ', count(kinds == PROJ_PATCH), ', lowest z ', minval(free%z(:, :, 1))
call fig_open(u, 'ref-body-fitted-face', 'a face put on the cap: constrained (gold), tolerance 0 (ember)', &
              -60._R8P, 35._R8P)
call fig_zoom(u, vector_R8P(0.35_R8P, -0.45_R8P, 1.9_R8P), 0.65_R8P)
call fig_surface(u, 'body', cyl)
call grid(free, 'hole', .true.)
call grid(constrained, 'foot', .false.)
call fig_label(u, 'hole', vector_R8P(0.72_R8P, -0.72_R8P, 1.7_R8P), 'boundary on the side')
call fig_label(u, 'accent', vector_R8P(0.3_R8P, -0.2_R8P, 2.02_R8P), 'interior on the cap')
call fig_close(u)
contains
   subroutine add(a, b, c)
   !< Append the triangle (a, b, c).
   type(vector_R8P), intent(in) :: a, b, c !< Vertices.

   m = m + 1 ; f(m)%vertex(1) = a ; f(m)%vertex(2) = b ; f(m)%vertex(3) = c ; f(m)%id = m
   endsubroutine add

   subroutine grid(b, style, rim)
   !< The lines of the face k = 1 of a block (rim: only the segments that touch its boundary), each segment drawn in 8
   !< pieces so that the depth sort keeps them over the long facets of the cap.
   type(curvilinear_block_t), intent(in) :: b     !< Block.
   character(*),              intent(in) :: style !< Style.
   logical,                   intent(in) :: rim   !< Boundary segments only.
   type(vector_R8P)                      :: p(1152), q(1152), e(2, 2)
   integer(I4P)                          :: k, l, n, h, o

   n = 0
   do l = 1, 9
      do k = 1, 8
         if (rim .and. .not. (l == 1 .or. l == 9 .or. k == 1 .or. k == 8)) cycle
         e(:, 1) = [vertex(b, k, l), vertex(b, k + 1_I4P, l)]
         e(:, 2) = [vertex(b, l, k), vertex(b, l, k + 1_I4P)]
         do o = 1, 2
            do h = 1, 8
               n = n + 1
               p(n) = e(1, o) * ((8 - h + 1) / 8._R8P) + e(2, o) * ((h - 1) / 8._R8P)
               q(n) = e(1, o) * ((8 - h) / 8._R8P) + e(2, o) * (h / 8._R8P)
            enddo
         enddo
      enddo
   enddo
   call fig_segments(u, style, p(1:n), q(1:n))
   endsubroutine grid

   function vertex(b, i, j) result(p)
   !< Vertex (i, j, 1) of a block.
   type(curvilinear_block_t), intent(in) :: b    !< Block.
   integer(I4P),              intent(in) :: i, j !< Indices.
   type(vector_R8P)                      :: p    !< Vertex.

   p = vector_R8P(b%x(i, j, 1), b%y(i, j, 1), b%z(i, j, 1))
   endfunction vertex
endprogram fig_ref_body_fitted_face
