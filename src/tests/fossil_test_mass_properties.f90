!< FOSSIL test: mass properties of closed surfaces (issue #59, roadmap #58 WP1, step 1).
!<
!< 1. Box a x b x c: mass, centre and the diagonal tensor m (b^2 + c^2)/12, ... (1e-14 relative); about a corner by
!<    the parallel-axis theorem.
!< 2. L-shaped prism (one closed shell): the two boxes' tensors moved to the common centre.
!< 3. Hollow box (outer shell, inner shell reversed): outer minus inner.
!< 4. The box under a proper rotation R: R I R^T; principal moments the box's, axes R's columns up to sign.
!< 5. The box translated by 1e6: the tensor about the centre unchanged to 1e-12 relative (conditioning).
!< 6. cube-04 (3072 facets, closed): volume = get_volume() = 1 and centre = get_centroid() = 0.5 to 1e-13 (dragon.stl is
!<    not watertight: it is rejected, check 8 covers that path).
!< 7. Marching-cubes sphere at three resolutions: I -> 2/5 m r^2 with observed order ~2.
!< 8. Rejections: an open box and an inside-out box give STATUS_INVALID_INPUT.
!< 9. symmetric3_eigen: random symmetric matrices (fixed seed) give A v = w v to 1e-13, det(v) = +1; repeated
!<    eigenvalues of a diagonal matrix give the identity axes.
program fossil_test_mass_properties

use fossil, only : surface_stl_object, facet_object, extract_isosurface, pose_from_axis_angle, rigid_pose_t, &
                   STATUS_OK, STATUS_INVALID_INPUT
use fossil_mass_properties, only : symmetric3_eigen
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: s
type(rigid_pose_t)       :: pose
type(vector_R8P)         :: pv
real(R8P)                :: mass, center(3), inertia(3,3), expected(3,3), pm(3), pa(3,3), vol
real(R8P)                :: m1, m2, c1(3), c2(3), cc(3), t1(3,3), t2(3,3), ref(3,3), err, e(3), order
real(R8P)                :: a, b, c
integer(I4P)             :: status, k
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. box a x b x c, density 2
a = 1._R8P ; b = 2._R8P ; c = 3._R8P
call make_prism(s, rect(0._R8P, 0._R8P, a, b), 0._R8P, c)
call s%mass_properties(mass=mass, center=center, inertia=inertia, density=2._R8P, status=status)
expected = box_tensor(2._R8P * a * b * c, a, b, c)
err = rel_err(inertia, expected)
ok = status == STATUS_OK .and. abs(mass - 12._R8P) <= 1e-14_R8P * 12._R8P .and. &
     maxval(abs(center - [0.5_R8P, 1._R8P, 1.5_R8P])) <= 1e-14_R8P * 3._R8P .and. err <= 1e-14_R8P
print '(A,ES10.3,A,L1)', '1. box: tensor rel err ', err, '  ok ', ok
all_ok = all_ok .and. ok
call s%mass_properties(inertia=inertia, density=2._R8P, about=[0._R8P, 0._R8P, 0._R8P], status=status)
err = rel_err(inertia, expected + shift_tensor(12._R8P, [0.5_R8P, 1._R8P, 1.5_R8P]))
ok = status == STATUS_OK .and. err <= 1e-14_R8P
print '(A,ES10.3,A,L1)', '1. box about a corner: rel err ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. L-shaped prism, depth 1: rectangles [0,2]x[0,1] and [0,1]x[1,2]
call make_prism(s, reshape([0._R8P, 0._R8P, 2._R8P, 0._R8P, 2._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P, 2._R8P, 0._R8P, 2._R8P], &
                           [2, 6]), 0._R8P, 1._R8P, &
                tri=reshape([4, 5, 6, 4, 6, 1, 4, 1, 2, 4, 2, 3], [3, 4]))
call s%mass_properties(mass=mass, center=center, inertia=inertia, status=status)
m1 = 2._R8P ; c1 = [1._R8P, 0.5_R8P, 0.5_R8P] ; t1 = box_tensor(m1, 2._R8P, 1._R8P, 1._R8P)
m2 = 1._R8P ; c2 = [0.5_R8P, 1.5_R8P, 0.5_R8P] ; t2 = box_tensor(m2, 1._R8P, 1._R8P, 1._R8P)
cc = (m1 * c1 + m2 * c2) / (m1 + m2)
ref = t1 + shift_tensor(m1, c1 - cc) + t2 + shift_tensor(m2, c2 - cc)
err = rel_err(inertia, ref)
ok = status == STATUS_OK .and. abs(mass - 3._R8P) <= 1e-14_R8P * 3._R8P .and. &
     maxval(abs(center - cc)) <= 1e-14_R8P * 2._R8P .and. err <= 1e-14_R8P
print '(A,ES10.3,A,L1)', '2. L prism: tensor rel err ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. hollow box: [0,4]^3 minus [1,2]x[1,3]x[1,2] (inner shell reversed)
call make_hollow(s)
call s%mass_properties(mass=mass, center=center, inertia=inertia, status=status)
m1 = 64._R8P ; c1 = [2._R8P, 2._R8P, 2._R8P] ; t1 = box_tensor(m1, 4._R8P, 4._R8P, 4._R8P)
m2 = 2._R8P  ; c2 = [1.5_R8P, 2._R8P, 1.5_R8P] ; t2 = box_tensor(m2, 1._R8P, 2._R8P, 1._R8P)
cc = (m1 * c1 - m2 * c2) / (m1 - m2)
ref = t1 + shift_tensor(m1, c1 - cc) - t2 - shift_tensor(m2, c2 - cc)
err = rel_err(inertia, ref)
ok = status == STATUS_OK .and. abs(mass - 62._R8P) <= 1e-14_R8P * 62._R8P .and. &
     maxval(abs(center - cc)) <= 1e-14_R8P * 4._R8P .and. err <= 1e-14_R8P
print '(A,ES10.3,A,L1)', '3. hollow box: tensor rel err ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 4. rotated box
pose = pose_from_axis_angle(axis=[1._R8P, 2._R8P, 3._R8P], angle=0.7_R8P)
call make_prism(s, rect(0._R8P, 0._R8P, a, b), 0._R8P, c, r=pose%r)
call s%mass_properties(inertia=inertia, principal_moments=pm, principal_axes=pa, status=status)
expected = box_tensor(a * b * c, a, b, c)
err = rel_err(inertia, matmul(pose%r, matmul(expected, transpose(pose%r))))
e(1) = maxval(abs(pm - [expected(3,3), expected(2,2), expected(1,1)])) / expected(1,1)
e(2) = 0._R8P
do k = 1, 3 ! principal axis k (ascending moments: z, y, x of the box) against R's column 4-k, up to sign
   e(2) = max(e(2), 1._R8P - abs(dot_product(pa(:,k), pose%r(:, 4 - k))))
enddo
ok = status == STATUS_OK .and. err <= 1e-14_R8P .and. e(1) <= 1e-14_R8P .and. e(2) <= 1e-14_R8P .and. &
     abs(det3(pa) - 1._R8P) <= 1e-14_R8P
print '(A,3ES10.3,A,L1)', '4. rotated box: tensor, moments, axes errs ', err, e(1), e(2), '  ok ', ok
all_ok = all_ok .and. ok

! ---- 5. the box translated by 1e6
call make_prism(s, rect(1.e6_R8P, 1.e6_R8P, a, b), 1.e6_R8P, 1.e6_R8P + c)
call s%mass_properties(center=center, inertia=inertia, status=status)
err = rel_err(inertia, box_tensor(a * b * c, a, b, c))
ok = status == STATUS_OK .and. err <= 1e-12_R8P .and. &
     maxval(abs(center - (1.e6_R8P + [0.5_R8P, 1._R8P, 1.5_R8P]))) <= 1e-9_R8P
print '(A,ES10.3,A,L1)', '5. box at 1e6: tensor rel err ', err, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 6. cube-04: consistency with get_volume and get_centroid, and the exact values
call s%load_from_file(file_name='src/tests/cube-04.stl', guess_format=.true.)
call s%mass_properties(center=center, volume=vol, inertia=inertia, status=status)
pv = s%get_centroid() ; cc = [pv%x, pv%y, pv%z]
e(1) = max(abs(vol - s%get_volume()), abs(vol - 1._R8P))
e(2) = max(maxval(abs(center - cc)), maxval(abs(center - 0.5_R8P)))
e(3) = rel_err(inertia, box_tensor(1._R8P, 1._R8P, 1._R8P, 1._R8P))
ok = status == STATUS_OK .and. e(1) <= 1e-13_R8P .and. e(2) <= 1e-13_R8P .and. e(3) <= 1e-13_R8P
print '(A,3ES10.3,A,L1)', '6. cube-04: volume, centre, tensor errs ', e, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 7. marching-cubes sphere, r = 1: I -> 2/5 m r^2 at second order
do k = 1, 3
   call make_sphere(s, 8 * 2 ** k + 1)
   call s%mass_properties(mass=mass, inertia=inertia, status=status)
   e(k) = abs((inertia(1,1) + inertia(2,2) + inertia(3,3)) / 3._R8P - 0.4_R8P * mass) / (0.4_R8P * mass)
   e(k) = e(k) + abs(mass - 4._R8P / 3._R8P * PI) / (4._R8P / 3._R8P * PI)
enddo
order = log(e(2) / e(3)) / log(2._R8P)
ok = status == STATUS_OK .and. order > 1.5_R8P .and. order < 2.5_R8P .and. e(3) < e(2) .and. e(2) < e(1)
print '(A,3ES10.3,A,F5.2,A,L1)', '7. MC sphere: errs ', e, '  order ', order, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 8. rejections
call make_prism(s, rect(0._R8P, 0._R8P, a, b), 0._R8P, c, drop_last=.true.)
call s%mass_properties(mass=mass, status=status)
ok = status == STATUS_INVALID_INPUT .and. mass == 0._R8P
call make_prism(s, rect(0._R8P, 0._R8P, a, b), 0._R8P, c, reversed=.true.)
call s%mass_properties(mass=mass, status=status)
ok = ok .and. status == STATUS_INVALID_INPUT .and. mass == 0._R8P
print '(A,L1)', '8. open and inside-out boxes rejected: ', ok
all_ok = all_ok .and. ok

! ---- 9. symmetric3_eigen
ok = check_eigen()
print '(A,L1)', '9. symmetric3_eigen: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   pure function rect(x0, y0, dx, dy) result(p)
   !< Counter-clockwise rectangle.
   real(R8P), intent(in) :: x0, y0, dx, dy
   real(R8P)             :: p(2,4)

   p = reshape([x0, y0, x0 + dx, y0, x0 + dx, y0 + dy, x0, y0 + dy], [2, 4])
   endfunction rect

   subroutine make_prism(surf, poly, z0, z1, tri, r, reversed, drop_last)
   !< Closed prism over the counter-clockwise polygon `poly`, from z0 to z1; caps triangulated by `tri` (default: a fan
   !< from vertex 1). Optional rotation `r` of every vertex; `reversed` flips every facet; `drop_last` leaves a hole.
   type(surface_stl_object), intent(inout)        :: surf
   real(R8P),                intent(in)           :: poly(:,:), z0, z1
   integer,                  intent(in), optional :: tri(:,:)
   real(R8P),                intent(in), optional :: r(3,3)
   logical,                  intent(in), optional :: reversed, drop_last
   real(R8P), allocatable                         :: v(:,:,:)
   integer,   allocatable                         :: t(:,:)
   integer                                        :: n, nt, i, j, f

   n = size(poly, 2)
   if (present(tri)) then
      t = tri
   else
      allocate(t(3, n - 2))
      do i = 1, n - 2
         t(:, i) = [1, i + 1, i + 2]
      enddo
   endif
   nt = size(t, 2)
   allocate(v(3, 3, 2 * nt + 2 * n))
   f = 0
   do i = 1, nt
      f = f + 1 ; v(:,1,f) = [poly(:,t(1,i)), z1] ; v(:,2,f) = [poly(:,t(2,i)), z1] ; v(:,3,f) = [poly(:,t(3,i)), z1]
      f = f + 1 ; v(:,1,f) = [poly(:,t(1,i)), z0] ; v(:,2,f) = [poly(:,t(3,i)), z0] ; v(:,3,f) = [poly(:,t(2,i)), z0]
   enddo
   do i = 1, n
      j = mod(i, n) + 1
      f = f + 1 ; v(:,1,f) = [poly(:,i), z0] ; v(:,2,f) = [poly(:,j), z0] ; v(:,3,f) = [poly(:,j), z1]
      f = f + 1 ; v(:,1,f) = [poly(:,i), z0] ; v(:,2,f) = [poly(:,j), z1] ; v(:,3,f) = [poly(:,i), z1]
   enddo
   if (present(r)) then
      do f = 1, size(v, 3)
         v(:,:,f) = matmul(r, v(:,:,f))
      enddo
   endif
   if (present(reversed)) then
      if (reversed) v = v(:, [1, 3, 2], :)
   endif
   nt = size(v, 3)
   if (present(drop_last)) then
      if (drop_last) nt = nt - 1
   endif
   call adopt(surf, v(:,:,1:nt))
   endsubroutine make_prism

   subroutine make_hollow(surf)
   !< [0,4]^3 with the void [1,2]x[1,3]x[1,2] (inner shell wound inward).
   type(surface_stl_object), intent(inout) :: surf
   type(surface_stl_object)                :: outer, inner
   real(R8P), allocatable                  :: v(:,:,:), w(:,:,:)

   call make_prism(outer, rect(0._R8P, 0._R8P, 4._R8P, 4._R8P), 0._R8P, 4._R8P)
   call make_prism(inner, rect(1._R8P, 1._R8P, 1._R8P, 2._R8P), 1._R8P, 2._R8P, reversed=.true.)
   v = vertices_of(outer) ; w = vertices_of(inner)
   call adopt(surf, reshape([v, w], [3, 3, size(v, 3) + size(w, 3)]))
   endsubroutine make_hollow

   function vertices_of(surf) result(v)
   !< The facet vertices of a surface.
   type(surface_stl_object), intent(in) :: surf
   real(R8P), allocatable               :: v(:,:,:)
   type(facet_object), pointer          :: fct
   integer                              :: f, k

   allocate(v(3, 3, surf%get_facets_number()))
   do f = 1, surf%get_facets_number()
      fct => surf%facet_at(f)
      do k = 1, 3
         v(:,k,f) = [fct%vertex(k)%x, fct%vertex(k)%y, fct%vertex(k)%z]
      enddo
   enddo
   endfunction vertices_of

   subroutine adopt(surf, v)
   !< Build the surface from facet vertices.
   type(surface_stl_object), intent(inout) :: surf
   real(R8P),                intent(in)    :: v(:,:,:)
   type(facet_object), allocatable         :: facets(:)
   integer                                 :: f, k

   allocate(facets(size(v, 3)))
   do f = 1, size(v, 3)
      do k = 1, 3
         facets(f)%vertex(k) = vector_R8P(v(1,k,f), v(2,k,f), v(3,k,f))
      enddo
      call facets(f)%compute_metrix
   enddo
   call surf%destroy
   call surf%adopt_facets(facets=facets)
   endsubroutine adopt

   subroutine make_sphere(surf, n)
   !< Marching-cubes unit sphere on an n^3 grid over [-1.5, 1.5]^3.
   type(surface_stl_object), intent(inout) :: surf
   integer,                  intent(in)    :: n
   real(R8P), allocatable                  :: values(:,:,:)
   type(facet_object), allocatable         :: facets(:)
   real(R8P)                               :: h
   integer                                 :: i, j, k
   integer(I4P)                            :: st

   allocate(values(n, n, n))
   h = 3._R8P / (n - 1)
   do k = 1, n
      do j = 1, n
         do i = 1, n
            values(i,j,k) = norm2([-1.5_R8P + (i - 1) * h, -1.5_R8P + (j - 1) * h, -1.5_R8P + (k - 1) * h]) - 1._R8P
         enddo
      enddo
   enddo
   call extract_isosurface(values=values, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                           bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=st)
   call surf%destroy
   call surf%adopt_facets(facets=facets)
   endsubroutine make_sphere

   pure function box_tensor(m, dx, dy, dz) result(t)
   !< Inertia tensor of a box about its centre.
   real(R8P), intent(in) :: m, dx, dy, dz
   real(R8P)             :: t(3,3)

   t = 0._R8P
   t(1,1) = m * (dy ** 2 + dz ** 2) / 12._R8P
   t(2,2) = m * (dx ** 2 + dz ** 2) / 12._R8P
   t(3,3) = m * (dx ** 2 + dy ** 2) / 12._R8P
   endfunction box_tensor

   pure function shift_tensor(m, d) result(t)
   !< Parallel-axis term m (|d|^2 1 - d d^T).
   real(R8P), intent(in) :: m, d(3)
   real(R8P)             :: t(3,3)
   integer               :: i

   t = -m * spread(d, 2, 3) * spread(d, 1, 3)
   do i = 1, 3
      t(i,i) = t(i,i) + m * dot_product(d, d)
   enddo
   endfunction shift_tensor

   pure function rel_err(x, y) result(e)
   !< Max-norm error relative to the largest entry of y.
   real(R8P), intent(in) :: x(3,3), y(3,3)
   real(R8P)             :: e

   e = maxval(abs(x - y)) / maxval(abs(y))
   endfunction rel_err

   pure function det3(m) result(d)
   !< Determinant.
   real(R8P), intent(in) :: m(3,3)
   real(R8P)             :: d

   d = m(1,1) * (m(2,2) * m(3,3) - m(3,2) * m(2,3)) - m(1,2) * (m(2,1) * m(3,3) - m(3,1) * m(2,3)) + &
       m(1,3) * (m(2,1) * m(3,2) - m(3,1) * m(2,2))
   endfunction det3

   function check_eigen() result(good)
   !< Random symmetric matrices (Park-Miller, fixed seed) and a repeated-eigenvalue case.
   logical      :: good
   real(R8P)    :: m(3,3), w(3), v(3,3), resid
   integer      :: trial, i, j
   integer(8)   :: seed

   good = .true.
   seed = 12345_8
   do trial = 1, 200
      do j = 1, 3
         do i = j, 3
            seed = mod(16807_8 * seed, 2147483647_8)
            m(i,j) = 2._R8P * real(seed, R8P) / 2147483647._R8P - 1._R8P
            m(j,i) = m(i,j)
         enddo
      enddo
      call symmetric3_eigen(m, w, v)
      resid = maxval(abs(matmul(m, v) - v * spread(w, 1, 3)))
      good = good .and. resid <= 1e-13_R8P .and. abs(det3(v) - 1._R8P) <= 1e-13_R8P .and. w(1) <= w(2) .and. w(2) <= w(3)
      good = good .and. maxval(abs(matmul(transpose(v), v) - reshape([1,0,0,0,1,0,0,0,1], [3,3]))) <= 1e-13_R8P
   enddo
   m = 0._R8P ; m(1,1) = 2._R8P ; m(2,2) = 2._R8P ; m(3,3) = 5._R8P
   call symmetric3_eigen(m, w, v)
   good = good .and. all(w == [2._R8P, 2._R8P, 5._R8P]) .and. all(v == reshape([1,0,0,0,1,0,0,0,1], [3,3]))
   endfunction check_eigen
endprogram fossil_test_mass_properties
