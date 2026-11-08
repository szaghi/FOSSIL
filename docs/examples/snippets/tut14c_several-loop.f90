toi_all = 1._R8P                                                                 ! the same by every pair
do i = 1, 3
   do j = i + 1, 4
      call bodies(i)%time_of_impact(bodies(j), poses0(i), poses1(i), poses0(j), poses1(j), s, status, &
                                     tolerance=1.e-6_R8P)
      if (status == CCD_CONTACT) then
         toi_all = min(toi_all, s)
         print '(A,I0,A,I0,A,F9.6)', '  pair ', i, '-', j, ': contact at s = ', s
      else
         print '(A,I0,A,I0,A)', '  pair ', i, '-', j, ': separate'
      endif
   enddo
enddo
print '(A,L1)', 'equal to the loop over all pairs, bit for bit: ', toi_all == toi
