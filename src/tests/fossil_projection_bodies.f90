!< FOSSIL test-side module: bodies with patches for the projection tests (issue #74).
module fossil_projection_bodies
!< FOSSIL test-side module: bodies with patches for the projection tests (issue #74).
!<
!< - cube_patches: src/tests/cube.stl ([0, 1]^3), one patch per face (named by the outward axis), or one patch;
!< - cylinder: radius r, height h along z from z = 0, n sides, caps as fans about their centres; patches `side`,
!<   `bottom`, `top` (or one);
!< - uv_sphere: radius r about the origin, nt latitude bands (even: the equator is a mesh line), np longitudes, poles as
!<   fans; patches `south` and `north` split at the equator (or one).
use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: cube_patches, cylinder, uv_sphere

real(R8P), parameter :: PI = acos(-1._R8P)

contains
   subroutine cube_patches(s, six)
   !< The unit cube, one patch per face if six.
   type(surface_stl_object), intent(inout) :: s    !< Cube.
   logical,                  intent(in)    :: six  !< One patch per face.
   type(surface_stl_object)                :: fresh
   type(facet_object), pointer             :: f
   character(2), parameter                 :: NAMES(6) = ['-x', '+x', '-y', '+y', '-z', '+z']
   integer(I4P)                            :: i, k, a

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   if (.not. six) return
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      a = 0
      do k = 1, 3
         if (all([f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x] == 0._R8P)) a = 1
         if (all([f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x] == 1._R8P)) a = 2
         if (all([f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y] == 0._R8P)) a = 3
         if (all([f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y] == 1._R8P)) a = 4
         if (all([f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z] == 0._R8P)) a = 5
         if (all([f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z] == 1._R8P)) a = 6
      enddo
      call s%set_patch([i], NAMES(a))
   enddo
   call s%analyze
   endsubroutine cube_patches

   subroutine cylinder(s, n, r, h, patches)
   !< Closed cylinder, outward triangles.
   type(surface_stl_object), intent(inout) :: s       !< Cylinder.
   integer(I4P),             intent(in)    :: n       !< Sides.
   real(R8P),                intent(in)    :: r, h    !< Radius, height.
   logical,                  intent(in)    :: patches !< Patches side / bottom / top.
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   type(vector_R8P)                        :: b0, b1, t0, t1, cb, ct
   integer(I4P)                            :: i, m

   s = fresh
   allocate(f(4 * n)) ; m = 0
   cb = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; ct = vector_R8P(0._R8P, 0._R8P, h)
   do i = 0, n - 1
      b0 = vector_R8P(r * cos(2 * PI * i / n), r * sin(2 * PI * i / n), 0._R8P)
      b1 = vector_R8P(r * cos(2 * PI * mod(i + 1, n) / n), r * sin(2 * PI * mod(i + 1, n) / n), 0._R8P)
      t0 = b0 ; t0%z = h ; t1 = b1 ; t1%z = h
      call add_tri(f, m, b0, b1, t1) ; call add_tri(f, m, b0, t1, t0) ! side
      call add_tri(f, m, cb, b1, b0)                                  ! bottom (outward -z)
      call add_tri(f, m, ct, t0, t1)                                  ! top
   enddo
   call s%adopt_facets(facets=f)
   if (patches) then
      call s%set_patch([(4 * i + 1, 4 * i + 2, i = 0, n - 1)], 'side')
      call s%set_patch([(4 * i + 3, i = 0, n - 1)], 'bottom')
      call s%set_patch([(4 * i + 4, i = 0, n - 1)], 'top')
      call s%analyze
   endif
   endsubroutine cylinder

   subroutine uv_sphere(s, nt, np, r, patches)
   !< UV sphere, outward triangles; nt even.
   type(surface_stl_object), intent(inout) :: s       !< Sphere.
   integer(I4P),             intent(in)    :: nt, np  !< Latitude bands, longitudes.
   real(R8P),                intent(in)    :: r       !< Radius.
   logical,                  intent(in)    :: patches !< Patches south / north.
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   integer(I4P), allocatable               :: north(:), south(:)
   integer(I4P)                            :: i, j, m, nn, ns

   s = fresh
   allocate(f(2 * np * nt), north(2 * np * nt), south(2 * np * nt)) ; m = 0 ; nn = 0 ; ns = 0
   do j = 0, nt - 1
      do i = 0, np - 1
         if (j < nt - 1) then
            call add_tri(f, m, sphere_point(i, j, nt, np, r), sphere_point(i + 1, j + 1, nt, np, r), &
                         sphere_point(i, j + 1, nt, np, r))
            if (j < nt / 2) then
               ns = ns + 1 ; south(ns) = m
            else
               nn = nn + 1 ; north(nn) = m
            endif
         endif
         if (j > 0) then
            call add_tri(f, m, sphere_point(i, j, nt, np, r), sphere_point(i + 1, j, nt, np, r), &
                         sphere_point(i + 1, j + 1, nt, np, r))
            if (j < nt / 2) then
               ns = ns + 1 ; south(ns) = m
            else
               nn = nn + 1 ; north(nn) = m
            endif
         endif
      enddo
   enddo
   f = f(1:m)
   call s%adopt_facets(facets=f)
   if (patches) then
      call s%set_patch(south(1:ns), 'south') ; call s%set_patch(north(1:nn), 'north')
      call s%analyze
   endif
   endsubroutine uv_sphere

   subroutine add_tri(f, m, a, b, c)
   !< Append the triangle (a, b, c).
   type(facet_object), intent(inout) :: f(:)    !< Facets.
   integer(I4P),       intent(inout) :: m       !< Count.
   type(vector_R8P),   intent(in)    :: a, b, c !< Vertices.

   m = m + 1 ; f(m)%vertex(1) = a ; f(m)%vertex(2) = b ; f(m)%vertex(3) = c ; f(m)%id = m
   endsubroutine add_tri

   function sphere_point(i, j, nt, np, r) result(x)
   !< Vertex (longitude i, latitude j) of the UV sphere (poles exact).
   integer(I4P), intent(in) :: i, j, nt, np !< Indices, counts.
   real(R8P),    intent(in) :: r            !< Radius.
   type(vector_R8P)         :: x            !< Vertex.
   real(R8P)                :: th, ph

   th = PI * j / nt ; ph = 2 * PI * mod(i, np) / np
   x = vector_R8P(r * sin(th) * cos(ph), r * sin(th) * sin(ph), -r * cos(th))
   if (j == 0) x = vector_R8P(0._R8P, 0._R8P, -r)
   if (j == nt) x = vector_R8P(0._R8P, 0._R8P, r)
   endfunction sphere_point
endmodule fossil_projection_bodies
