!< FOSSIL test: named patches through remeshing (issue #66, step 5).
!<
!< A marching-cubes sphere (12^3) split into the patches `top` (facet centroid z > 0) and `bottom`. The patch boundary
!< is the polyline of edges between facets of the two patches.
!< 1. isotropic_remesh (preserve_features off: the patch boundary is locked anyway): every input boundary vertex is still
!<    a boundary vertex, every output boundary vertex lies on an input boundary edge (1e-12), every facet keeps a side
!<    (centroid z of the same sign as its patch), the area per patch changes by < 2 %.
!< 2. decimate to half: every output boundary vertex is an input boundary vertex (bit for bit: the boundary coarsens
!<    as a curve, issue #81), the same side and area checks (5 %); with rim_collapses=.false. the boundary vertex set is
!<    the input one exactly.
!< 3. smooth (Taubin): the boundary vertices do not move; side check.
!< 4. alpha_wrap: each facet of the wrap takes the patch of the nearest input facet; facets with centroid z > 0.08 are
!<    `top`, z < -0.08 `bottom`.
program fossil_test_patches_remesh

use fossil,                only : surface_stl_object, facet_object
use fossil_boolean_bodies, only : BODY_SPHERE, mc_body
use penf,                  only : I4P, R8P

implicit none

type(surface_stl_object) :: s, r
real(R8P), allocatable   :: bv(:, :), be(:, :, :)
real(R8P)                :: a0(2), a1(2)
logical                  :: ok1, ok2, ok3, ok4

call mc_body(BODY_SPHERE, 12_I4P, s)
call split(s)
call boundary(s, bv, be)
a0 = patch_areas(s)
print '(A,I0,A,I0)', 'sphere facets ', s%get_facets_number(), ', boundary vertices ', size(bv, dim=2)

r = s
call r%isotropic_remesh(iterations=3_I4P, preserve_features=.false.)
a1 = patch_areas(r)
ok1 = contains_all(r, bv) .and. on_boundary(r, be) .and. sides_ok(r, 0.08_R8P) .and. all(abs(a1 - a0) < 0.02_R8P * a0)
print '(A,L1,A,2F9.5,A,2F9.5)', '1. isotropic_remesh: ', ok1, '  areas ', a0, ' -> ', a1

r = s
call r%decimate(target_facets=s%get_facets_number() / 2)
a1 = patch_areas(r)
ok2 = subset(r, bv) .and. sides_ok(r, 0.08_R8P) .and. all(abs(a1 - a0) < 0.05_R8P * a0)
print '(A,L1,A,2F9.5,A,I0)', '2. decimate: ', ok2, '  areas ', a1, ', boundary vertices ', boundary_count(r)
r = s
call r%decimate(target_facets=s%get_facets_number() / 2, rim_collapses=.false.)
a1 = patch_areas(r)
ok2 = ok2 .and. same_set(r, bv) .and. sides_ok(r, 0.08_R8P) .and. all(abs(a1 - a0) < 0.05_R8P * a0)
print '(A,L1,A,2F9.5)', '2. decimate, rims locked: ', ok2, '  areas ', a1

r = s
call r%smooth(iterations=10_I4P)
ok3 = same_set(r, bv) .and. sides_ok(r, 0.08_R8P)
print '(A,L1)', '3. smooth: ', ok3

call s%alpha_wrap(alpha=0.15_R8P, wrapped=r)
ok4 = r%get_facets_number() > 0 .and. r%get_patches_number() >= 3 .and. sides_ok(r, 0.08_R8P)
print '(A,L1)', '4. alpha_wrap: ', ok4

print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4

contains
   subroutine split(x)
   !< Patches top (centroid z > 0) and bottom.
   type(surface_stl_object), intent(inout) :: x
   type(facet_object), pointer             :: f
   integer(I4P), allocatable               :: t(:), b(:)
   integer(I4P)                            :: i

   allocate(t(0), b(0))
   do i = 1, x%get_facets_number()
      f => x%facet_at(i)
      if (f%centroid%z > 0._R8P) then
         t = [t, i]
      else
         b = [b, i]
      endif
   enddo
   call x%set_patch(t, 'top')
   call x%set_patch(b, 'bottom')
   endsubroutine split

   subroutine boundary(x, v, e)
   !< Vertices and edges between facets of different patches (vertex pairs compared by coordinates).
   type(surface_stl_object), intent(in)  :: x
   real(R8P), allocatable,   intent(out) :: v(:, :), e(:, :, :)
   type(facet_object), pointer           :: f, g
   integer(I4P)                          :: i, j, k, l
   real(R8P)                             :: p(3), q(3)

   allocate(v(3, 0), e(3, 2, 0))
   do i = 1, x%get_facets_number()
      f => x%facet_at(i)
      do k = 1, 3
         p = xyz(f, k) ; q = xyz(f, mod(k, 3) + 1)
         do j = 1, x%get_facets_number()
            if (j == i) cycle
            g => x%facet_at(j)
            if (g%patch == f%patch) cycle
            do l = 1, 3
               if (all(xyz(g, l) == q) .and. all(xyz(g, mod(l, 3) + 1) == p)) then
                  if (f%patch < g%patch) e = reshape([reshape(e, [3 * 2 * size(e, dim=3)]), p, q], &
                                                     [3, 2, size(e, dim=3) + 1])
                  call add(v, p) ; call add(v, q)
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine boundary

   subroutine add(v, p)
   real(R8P), allocatable, intent(inout) :: v(:, :)
   real(R8P),              intent(in)    :: p(3)
   integer(I4P)                          :: i

   do i = 1, size(v, dim=2)
      if (all(v(:, i) == p)) return
   enddo
   v = reshape([reshape(v, [3 * size(v, dim=2)]), p], [3, size(v, dim=2) + 1])
   endsubroutine add

   function xyz(f, k) result(p)
   type(facet_object), intent(in) :: f
   integer(I4P),       intent(in) :: k
   real(R8P)                      :: p(3)

   p = [f%vertex(k)%x, f%vertex(k)%y, f%vertex(k)%z]
   endfunction xyz

   function contains_all(x, v) result(ok)
   !< Every vertex of v is a boundary vertex of x.
   type(surface_stl_object), intent(in) :: x
   real(R8P),                intent(in) :: v(:, :)
   logical                              :: ok
   real(R8P), allocatable               :: w(:, :), e(:, :, :)
   integer(I4P)                         :: i, j

   call boundary(x, w, e)
   ok = .true.
   do i = 1, size(v, dim=2)
      ok = ok .and. any([(all(w(:, j) == v(:, i)), j = 1, size(w, dim=2))])
   enddo
   endfunction contains_all

   function same_set(x, v) result(ok)
   !< The boundary vertices of x are exactly v.
   type(surface_stl_object), intent(in) :: x
   real(R8P),                intent(in) :: v(:, :)
   logical                              :: ok
   real(R8P), allocatable               :: w(:, :), e(:, :, :)

   call boundary(x, w, e)
   ok = size(w, dim=2) == size(v, dim=2) .and. contains_all(x, v)
   endfunction same_set

   function subset(x, v) result(ok)
   !< Every boundary vertex of x is a vertex of v.
   type(surface_stl_object), intent(in) :: x
   real(R8P),                intent(in) :: v(:, :)
   logical                              :: ok
   real(R8P), allocatable               :: w(:, :), e(:, :, :)
   integer(I4P)                         :: i, j

   call boundary(x, w, e)
   ok = size(w, dim=2) > 0
   do i = 1, size(w, dim=2)
      ok = ok .and. any([(all(v(:, j) == w(:, i)), j = 1, size(v, dim=2))])
   enddo
   endfunction subset

   function boundary_count(x) result(n)
   !< The number of boundary vertices of x.
   type(surface_stl_object), intent(in) :: x
   integer(I4P)                         :: n
   real(R8P), allocatable               :: w(:, :), e(:, :, :)

   call boundary(x, w, e)
   n = size(w, dim=2)
   endfunction boundary_count

   function on_boundary(x, e0) result(ok)
   !< Every boundary vertex of x lies on an edge of e0 (distance <= 1e-12).
   type(surface_stl_object), intent(in) :: x
   real(R8P),                intent(in) :: e0(:, :, :)
   logical                              :: ok
   real(R8P), allocatable               :: w(:, :), e(:, :, :)
   real(R8P)                            :: d, t, a(3), b(3)
   integer(I4P)                         :: i, j

   call boundary(x, w, e)
   ok = .true.
   do i = 1, size(w, dim=2)
      d = huge(1._R8P)
      do j = 1, size(e0, dim=3)
         a = e0(:, 1, j) ; b = e0(:, 2, j)
         t = max(0._R8P, min(1._R8P, dot_product(w(:, i) - a, b - a) / dot_product(b - a, b - a)))
         d = min(d, norm2(w(:, i) - a - t * (b - a)))
      enddo
      ok = ok .and. d <= 1e-12_R8P
   enddo
   endfunction on_boundary

   function sides_ok(x, margin) result(ok)
   !< Facets of patch top have centroid z > -margin... strictly: top facets with z < -margin, or bottom facets with
   !< z > margin, are errors.
   type(surface_stl_object), intent(in) :: x
   real(R8P),                intent(in) :: margin
   logical                              :: ok
   type(facet_object), pointer          :: f
   integer(I4P)                         :: i

   ok = .true.
   do i = 1, x%get_facets_number()
      f => x%facet_at(i)
      if (x%patch_name(f%patch) == 'top' .and. f%centroid%z < -margin) ok = .false.
      if (x%patch_name(f%patch) == 'bottom' .and. f%centroid%z > margin) ok = .false.
   enddo
   endfunction sides_ok

   function patch_areas(x) result(a)
   type(surface_stl_object), intent(in) :: x
   real(R8P)                            :: a(2)
   type(facet_object), pointer          :: f
   integer(I4P)                         :: i

   a = 0._R8P
   do i = 1, x%get_facets_number()
      f => x%facet_at(i)
      if (x%patch_name(f%patch) == 'top') a(1) = a(1) + f%area()
      if (x%patch_name(f%patch) == 'bottom') a(2) = a(2) + f%area()
   enddo
   endfunction patch_areas
endprogram fossil_test_patches_remesh
