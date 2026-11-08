call c%query_many(points=p, distances=d, closest_points=foot, vertices=vtx, normals=nrm, bodies=bod)
near = abs(d) < BAND
err = 0._R8P
do k = 1, n
   if (near(k) .and. vtx(k) == 0) err = max(err, norm2(p(:, k) - foot(:, k) - d(k) * nrm(:, k)))
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'bodies ', c%get_bodies_number(), ', points within 0.3: ', count(near), &
                               ', body 1: ', count(near .and. bod == 1), ', body 2: ', count(near .and. bod == 2)
print '(A,I0,A,I0,A,I0)', 'feet inside a segment ', count(near .and. vtx == 0), ', at a vertex ', &
                          count(near .and. vtx /= 0), ', of them at the square corners ', &
                          count(near .and. vtx /= 0 .and. bod == 1)
print '(A,ES9.2)', 'feet inside a segment: max |p - foot - d n| ', err
