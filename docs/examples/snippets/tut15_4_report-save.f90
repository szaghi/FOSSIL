call save_block_report(report, 'tut15-report.txt')
open(newunit=unit, file='tut15-report.txt', action='read')
do i = 1, 7
   read(unit, '(A)') line
   if (i /= 5) print '(A)', trim(line)                                  ! line 5, the cell quality, is long
enddo
close(unit)
