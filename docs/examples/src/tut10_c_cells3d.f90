!run tut10c tut10_c_cells3d
!fig tut10-cells3d
!fig anim-tut10-sweep
!< Tutorial 10, step 3: the 3-D cut cells of the dragon (cut_cells_3d): volume fractions, face fractions, pieces,
!< flagged cells; the faces of the plane z = 0 against the 2-D cut cells of the section of chapter 7.
!< The figures (outside the regions): the cells below z = 0 that are not all fluid as cubes, and an animation of the
!< plane sweeping the block: per layer of cells, its fluid fractions and the sections at its bottom and top planes.
program tut10_c_cells3d
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_quad, fig_segments, fig_anim, fig_frame, fig_close
use vecfor,     only : vector_R8P
use fossil,     only : slice_surface
!region cells
use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(cut_cells_3d_t)     :: cut
type(cut_cells_2d_t)     :: cut2
real(R8P)                :: x(25), h
integer(I4P)             :: i, c(3), p, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)    ! chapter 3
h = 0.05_R8P
x = [(-0.6_R8P + h * i, i = 0, 24)]                  ! 24^3 cells, h = 0.05; z = 0 is the plane x(13)
call cut_cells_3d(surface=surface, x=x, y=x, z=x, cut=cut)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells: fluid ', count(cut%alpha == 1._R8P), ', solid ', count(cut%alpha == 0._R8P), &
                                ', cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P), &
                                ', flagged ', count(cut%pieces < 0)
!endregion cells
!region faces
call section%load_from_file(file_name='dragon-section.dat', status=status)        ! chapter 7: the section z = 0
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
call section%cut_cells(x=x, y=x, cut=cut2)                                         ! its 2-D cut cells, same grid
print '(A,I0,A,ES9.2)', 'faces z = 0: cut ', count(cut%face_z(:, :, 13) > 0._R8P .and. cut%face_z(:, :, 13) < 1._R8P), &
                        ', largest |face_z - 2-D alpha| ', maxval(abs(cut%face_z(:, :, 13) - cut2%alpha))
!endregion faces
!region pieces
print '(A,I0,A,I0)', 'pieces ', cut%np, ', cells with two or more ', count(cut%pieces > 1)
c = findloc(cut%pieces, 2)                                                         ! the first cell with two pieces
p = cut%first_piece(c(1), c(2), c(3))
print '(A,I0,A,I0,A,I0,A,F8.5,A,2F8.5)', 'cell (', c(1), ',', c(2), ',', c(3), ') alpha ', &
   cut%alpha(c(1), c(2), c(3)), ' = ', cut%piece_volume(p:p + 1) / h**3
print '(A,6F7.3)', '  faces of piece 1 (-x +x -y +y -z +z)', cut%piece_face(:, p)
print '(A,6F7.3)', '  faces of piece 2                    ', cut%piece_face(:, p + 1)
print '(A,6F7.3)', '  faces of the cell                   ', cut%face_x(c(1):c(1) + 1, c(2), c(3)), &
                   cut%face_y(c(1), c(2):c(2) + 1, c(3)), cut%face_z(c(1), c(2), c(3):c(3) + 1)
!endregion pieces
call figures
contains
   subroutine figures
   !< The cubes below z = 0, then the animation of the layers.
   type(surface_stl_object)      :: top
   type(curve_object)            :: cs
   type(vector_R8P), allocatable :: a(:), b(:), a1(:), b1(:)
   character(96)                 :: cap
   integer(I4P)                  :: u, ii, jj, kk, k

   top = surface                                                                   ! the dragon above z = 0
   call top%clip(bmin=vec(-1._R8P, -1._R8P, 0._R8P), bmax=vec(1._R8P, 1._R8P, 1._R8P))
   call fig_open(u, 'tut10-cells3d', 'cut_cells_3d below z = 0: fluid fraction, solid cells cream, flagged in ember', &
                 -60._R8P, 50._R8P)
   call fig_surface(u, 'ghost', top)
   do kk = 1, 12 ; do jj = 1, 24 ; do ii = 1, 24
      if (.not. is_drawn(ii, jj, kk)) cycle
      if (.not. is_drawn(ii - 1, jj, kk)) call cube_face(u, ii, jj, kk, 1)
      if (.not. is_drawn(ii + 1, jj, kk)) call cube_face(u, ii, jj, kk, 2)
      if (.not. is_drawn(ii, jj - 1, kk)) call cube_face(u, ii, jj, kk, 3)
      if (.not. is_drawn(ii, jj + 1, kk)) call cube_face(u, ii, jj, kk, 4)
      if (.not. is_drawn(ii, jj, kk - 1)) call cube_face(u, ii, jj, kk, 5)
      if (.not. is_drawn(ii, jj, kk + 1)) call cube_face(u, ii, jj, kk, 6)
   enddo ; enddo ; enddo
   call fig_close(u)

   call fig_open(u, 'anim-tut10-sweep', '', 0._R8P, 90._R8P)
   call fig_zoom(u, vec(-0.025_R8P, 0._R8P, 0._R8P), 0.36_R8P)
   call fig_anim(u, 0.6_R8P)
   do k = 1, 24
      if (all(cut%alpha(:, :, k) == 1._R8P)) cycle                                ! the layers the dragon reaches
      write(cap, '(A,F6.3,A,F6.3,A,I0,A,I0,A)') 'z in [', x(k), ',', x(k + 1), ']: ', &
            count(cut%alpha(:, :, k) > 0._R8P .and. cut%alpha(:, :, k) < 1._R8P), ' cut, ', &
            count(cut%pieces(:, :, k) < 0), ' flagged; section at bottom teal, top red'
      call fig_frame(u, caption=trim(cap))
      do jj = 1, 24 ; do ii = 1, 24
         if (cut%pieces(ii, jj, k) < 0) then
            call fig_quad(u, 'patch', vec(x(ii), x(jj), 0._R8P), vec(x(ii + 1), x(jj), 0._R8P), &
                          vec(x(ii + 1), x(jj + 1), 0._R8P), vec(x(ii), x(jj + 1), 0._R8P))
         elseif (cut%alpha(ii, jj, k) < 1._R8P) then
            call fig_quad(u, 'frac', vec(x(ii), x(jj), 0._R8P), vec(x(ii + 1), x(jj), 0._R8P), &
                          vec(x(ii + 1), x(jj + 1), 0._R8P), vec(x(ii), x(jj + 1), 0._R8P), value=cut%alpha(ii, jj, k))
         endif
      enddo ; enddo
      call slice_surface(surface, 3_I4P, x(k), cs, status=status)                  ! the bottom plane of the layer
      call curve_segments(cs, 1.e-3_R8P, a, b)
      call slice_surface(surface, 3_I4P, x(k + 1), cs, status=status)              ! and its top plane
      call curve_segments(cs, 2.e-3_R8P, a1, b1)
      if (size(a) > 0) call fig_segments(u, 'probe', a, b)
      if (size(a1) > 0) call fig_segments(u, 'hole', a1, b1)
   enddo
   call fig_close(u)
   endsubroutine figures

   subroutine curve_segments(cv, lift, a, b)
   !< The segments of a section (normal z), at height lift.
   type(curve_object),            intent(in)  :: cv         !< Section.
   real(R8P),                     intent(in)  :: lift       !< Height.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Ends.
   real(R8P)                                  :: uv(2)
   integer(I4P)                               :: s, m

   allocate(a(cv%get_vertices_number()), b(cv%get_vertices_number())) ; m = 0
   do s = 1, cv%get_vertices_number()
      if (cv%get_next(s) == 0_I4P) cycle
      m = m + 1
      uv = cv%get_vertex(s) ; a(m) = vec(uv(1), uv(2), lift)
      uv = cv%get_vertex(cv%get_next(s)) ; b(m) = vec(uv(1), uv(2), lift)
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine curve_segments

   pure function is_drawn(ii, jj, kk) result(yes)
   !< The cell is below z = 0 and not all fluid: cut, solid or flagged.
   integer(I4P), intent(in) :: ii, jj, kk !< Cell.
   logical                  :: yes        !< Drawn.

   yes = .false.
   if (min(ii, jj, kk) < 1 .or. max(ii, jj) > 24 .or. kk > 12) return
   yes = cut%alpha(ii, jj, kk) < 1._R8P
   endfunction is_drawn

   subroutine cube_face(u, ii, jj, kk, f)
   !< Face f (-x, +x, -y, +y, -z, +z) of a cell, by alpha, or in ember if the cell is flagged.
   integer(I4P), intent(in) :: u, ii, jj, kk, f !< Scene, cell, face.
   real(R8P)                :: lo(3), hi(3), q(3, 4)
   integer(I4P)             :: m, ia, ib

   lo = [x(ii), x(jj), x(kk)] ; hi = [x(ii + 1), x(jj + 1), x(kk + 1)]
   m = (f + 1) / 2 ; ia = mod(m, 3) + 1 ; ib = mod(m + 1, 3) + 1
   q(m, :) = merge(lo(m), hi(m), mod(f, 2) == 1)
   q(ia, :) = [lo(ia), hi(ia), hi(ia), lo(ia)] ; q(ib, :) = [lo(ib), lo(ib), hi(ib), hi(ib)]
   if (cut%pieces(ii, jj, kk) < 0) then
      call fig_quad(u, 'patch', vec(q(1, 1), q(2, 1), q(3, 1)), vec(q(1, 2), q(2, 2), q(3, 2)), &
                    vec(q(1, 3), q(2, 3), q(3, 3)), vec(q(1, 4), q(2, 4), q(3, 4)))
   else
      call fig_quad(u, 'frac', vec(q(1, 1), q(2, 1), q(3, 1)), vec(q(1, 2), q(2, 2), q(3, 2)), &
                    vec(q(1, 3), q(2, 3), q(3, 3)), vec(q(1, 4), q(2, 4), q(3, 4)), value=cut%alpha(ii, jj, kk))
   endif
   endsubroutine cube_face

   function vec(p, q, z) result(v)
   !< The point (p, q, z).
   real(R8P), intent(in) :: p, q, z !< Coordinates.
   type(vector_R8P)      :: v       !< Point.

   v%x = p ; v%y = q ; v%z = z
   endfunction vec
endprogram tut10_c_cells3d
