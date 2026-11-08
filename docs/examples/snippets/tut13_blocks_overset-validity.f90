call blocks(1)%geometry(area=area)
print '(A,F10.6,A,F10.6)', 'O-grid area ', sum(area), ', annulus between the ellipses ', &
      PI * c**2 * (cosh(mu) * sinh(mu) - cosh(mu0) * sinh(mu0))
call blocks(1)%validity(valid)
print '(A,I0,A,I0)', 'valid cells ', count(valid == CELL_VALID), ' of ', size(valid)
call blocks(1)%quality(skewness=skew, orthogonality=ortho, aspect=aspect)
print '(A,F6.3,A,F6.3,A,F6.2)', 'max skewness ', maxval(skew), ', min orthogonality ', minval(ortho), &
      ', max aspect ', maxval(aspect)
folded = blocks(1)
folded%u(13, 5) = folded%u(13, 7) ; folded%v(13, 5) = folded%v(13, 7)  ! a vertex dragged two layers out
call folded%validity(valid)
print '(A,I0,A)', 'folded copy: ', count(valid == CELL_INVALID), ' invalid cells'
do j = 1, NJ - 1
   do i = 1, NI - 1
      if (valid(i, j) == CELL_INVALID) print '(A,2I3)', '  invalid cell', i, j
   enddo
enddo
