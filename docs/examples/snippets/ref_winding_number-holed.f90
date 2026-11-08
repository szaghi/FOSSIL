facets => cube%facets_ref()                                   ! drop the 8 facets of the top face z = 1
allocate(kept(count(facets%normal%z < 0.5_R8P)))
j = 0
do i = 1, size(facets)
   if (facets(i)%normal%z > 0.5_R8P) cycle
   j = j + 1 ; kept(j) = facets(i)
enddo
call holed%adopt_facets(facets=kept)
print '(A,L1,A,F6.4,A,F6.4)', 'open box: watertight ', holed%is_watertight(), ', w at the centre ', &
      holed%winding_number(point=q(1)), ', below the opening ', holed%winding_number(point=vector_R8P(0.5_R8P, &
      0.5_R8P, 0.9_R8P))
