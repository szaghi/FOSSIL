x = [(-0.6_R8P + 0.015_R8P * i, i = 0, 80)]
call field%initialize(x=x, y=x)                                   ! no z: a 2-D field
call field%compute(curve=section, status=status)
print '(A,I0,A,I0)', '2-D field: band cells ', field%get_band_cells(), ', inside ', count(field%phi < 0._R8P)
