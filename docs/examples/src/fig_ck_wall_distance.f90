!run figck-wall-distance fig_ck_wall_distance
!fig ck-wall-distance
!< Figure of the cookbook recipe "Compute the wall distance up to a cutoff": wall_distance_many with the recipe's
!< cutoff 0.25 on the cells of the plane z = 0.5 through cube-01, seen from above, coloured by d / cutoff (cream at
!< the wall, sepia at the cutoff: the plateau is every cell with no wall closer than the cutoff). The recipe's first
!< three points lie on this plane and are labelled with their distances; the fourth, (3, 3, 3), is off the frame.
program fig_ck_wall_distance
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_quad, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: nc = 60                              ! cells per side of the plane drawn
real(R8P),    parameter  :: x0 = -0.5_R8P, x1 = 1.5_R8P, cutoff = 0.25_R8P
type(surface_stl_object) :: cube
type(vector_R8P)         :: points(4), centres(nc * nc), a(4), b(4), o, q, off(3)
real(R8P)                :: d(4), dc(nc * nc), h, cd
integer(I4P)             :: body(4), bc(nc * nc), i, j, n, u
character(80)            :: text

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the recipe's points and call
points = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.1_R8P, 0.5_R8P, 0.5_R8P), &
          vector_R8P(1.3_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(3.0_R8P, 3.0_R8P, 3.0_R8P)]
call cube%wall_distance_many(points=points, distances=d, bodies=body, cutoff=cutoff)

h = (x1 - x0) / nc                                                ! the same call on the cell centres of z = 0.5
n = 0
do j = 1, nc ; do i = 1, nc
   n = n + 1
   centres(n) = (x0 + (i - 0.5_R8P) * h) * ex_R8P + (x0 + (j - 0.5_R8P) * h) * ey_R8P + 0.5_R8P * ez_R8P
enddo ; enddo
call cube%wall_distance_many(points=centres, distances=dc, bodies=bc, cutoff=cutoff)

write(text, '(A,I0,A)') 'wall_distance_many(cutoff=0.25) on z = 0.5: ', count(bc == 0), ' cells at the cutoff'
call fig_open(u, 'ck-wall-distance', trim(text), -90._R8P, 90._R8P, height=520_I4P)
n = 0
do j = 1, nc ; do i = 1, nc
   n = n + 1
   call fig_quad(u, 'frac', p(x0 + (i - 1) * h, x0 + (j - 1) * h), p(x0 + i * h, x0 + (j - 1) * h), &
                 p(x0 + i * h, x0 + j * h), p(x0 + (i - 1) * h, x0 + j * h), value=dc(n) / cutoff)
enddo ; enddo
a = [p(0._R8P, 0._R8P, 1.e-3_R8P), p(1._R8P, 0._R8P, 1.e-3_R8P), p(1._R8P, 1._R8P, 1.e-3_R8P), &
     p(0._R8P, 1._R8P, 1.e-3_R8P)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'hole', a, b)                                ! the section of the cube
call cube%compute_distance(point=points(2), distance=cd, closest_point=q)   ! the wall point of the one hit
q%z = q%z + 2.e-3_R8P
call fig_segment(u, 'probe', points(2) + 2.e-3_R8P * ez_R8P, q)
call fig_points(u, 'probe', [(points(i) + 2.e-3_R8P * ez_R8P, i = 1, 3)])
call fig_points(u, 'foot', [q])
off = [-0.1_R8P * ex_R8P + 0.06_R8P * ey_R8P, -0.12_R8P * ex_R8P - 0.1_R8P * ey_R8P, &  ! label offsets
       -0.02_R8P * ex_R8P + 0.06_R8P * ey_R8P]
do i = 1, 3
   if (body(i) == 0_I4P) then
      write(text, '(F4.2,A)') d(i), ' = cutoff'
   else
      write(text, '(F4.2,A,I0)') d(i), ', body ', body(i)
   endif
   call fig_label(u, 'probe', points(i) + off(i) + 3.e-3_R8P * ez_R8P, trim(text))
enddo
o = p(x0 + 0.05_R8P, x0 + 0.05_R8P, 3.e-3_R8P)                     ! axes glyph at a corner
call fig_segment(u, 'box', o, o + 0.2_R8P * ex_R8P) ; call fig_label(u, 'box', o + 0.2_R8P * ex_R8P, 'x')
call fig_segment(u, 'box', o, o + 0.2_R8P * ey_R8P) ; call fig_label(u, 'box', o + 0.2_R8P * ey_R8P, 'y')
call fig_close(u)
contains
   function p(xx, yy, dz) result(r)
   !< A point of the plane z = 0.5, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: xx, yy !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: r      !< Point.

   r = xx * ex_R8P + yy * ey_R8P + 0.5_R8P * ez_R8P
   if (present(dz)) r%z = r%z + dz
   endfunction p
endprogram fig_ck_wall_distance
