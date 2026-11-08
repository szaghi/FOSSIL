call swept_cells_2d(section, q0, q1, x, x, sw, check_input=.true.)    ! the section crosses itself (chapter 7)
print '(A,L1)', '  check_input: self-intersecting ', sw%self_intersecting
